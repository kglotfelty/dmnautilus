/*                                                                
**  Copyright (C) 2004-2008  Smithsonian Astrophysical Observatory 
*/                                                                

/*                                                                          */
/*  This program is free software; you can redistribute it and/or modify    */
/*  it under the terms of the GNU General Public License as published by    */
/*  the Free Software Foundation; either version 3 of the License, or       */
/*  (at your option) any later version.                                     */
/*                                                                          */
/*  This program is distributed in the hope that it will be useful,         */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/*  GNU General Public License for more details.                            */
/*                                                                          */
/*  You should have received a copy of the GNU General Public License along */
/*  with this program; if not, write to the Free Software Foundation, Inc., */
/*  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.             */
/*                                                                          */




#include <dslib.h>
#include <ascdm.h>
#include <stdlib.h>
#include <math.h>
#include <histlib.h>

#define FLOOR(x)  floor((x))
#define CEIL(x)   ceil((x))

#include <cxcregion.h>

regRegion *maskRegion;
double Crpix[2], Cdelt[2], Crval[2];
dmDescriptor *PhysDes;




/* Load an image into a float array */
void get_data( dmDescriptor *inDes, long npix, float *data )
{
  long ii;
  /* Read data image */
  switch ( dmGetDataType(inDes) ) {

  case dmBYTE: {
    unsigned char *sd = (unsigned char*)calloc(npix,sizeof(unsigned char));
    dmGetArray_ub( inDes, sd, npix );
    for (ii=0;ii<npix;ii++) { data[ii] = sd[ii]; }
    free(sd);
    break;
  }

  case dmSHORT: {
    short *sd = (short*)calloc(npix,sizeof(short));
    dmGetArray_s( inDes, sd, npix );
    for (ii=0;ii<npix;ii++) { data[ii] = sd[ii]; }
    free(sd);
    break;
  }
    
  case dmUSHORT: {
    unsigned short *sd = (unsigned short*)calloc(npix,sizeof(unsigned short));
    dmGetArray_us( inDes, sd, npix );
    for (ii=0;ii<npix;ii++) { data[ii] = sd[ii]; }
    free(sd);
    break;
  }
  case dmLONG: {
    long *sd = (long*)calloc(npix,sizeof(long));
    dmGetArray_l( inDes, sd, npix );
    for (ii=0;ii<npix;ii++) { data[ii] = sd[ii]; }
    free(sd);
    break;
  }
    
  case dmULONG: {
    unsigned long *sd = (unsigned long*)calloc(npix,sizeof(unsigned long));
    dmGetArray_ul( inDes, sd, npix );
    for (ii=0;ii<npix;ii++) { data[ii] = sd[ii]; }
    free(sd);
    break;
  }

  case dmFLOAT: {
    dmGetArray_f( inDes, data, npix );
    break;
  }
  case dmDOUBLE: {
    double *sd = (double*)calloc(npix,sizeof(double));
    dmGetArray_d( inDes, sd, npix );
    for (ii=0;ii<npix;ii++) { data[ii] = sd[ii]; }
    free(sd);
    break;
  }


  default:
    err_msg("ERROR: Unknown datatype\n");
    exit(-1); /* Should never get to here */
    break;

  }
  
}




double get_snr( 
        float *data,   /* i: data array */
        float *derr,   /* i: error array */
        long   xlen,   /* i: length of x-axis (full img) */
        long   ylen,   /* i: length of y-axis (full img) */
        long   xs,     /* i: start of x-axis (sub img) */
        long   ys,     /* i: start of y-axis (sub img) */
        long   xl,     /* i: length of x-axis (sub img) */
        long   yl,      /* i: length of y-axis (sub img) */
        float  *oval
  )
{
  
  float val;
  float noise;
  float locsnr;
  long ii,jj;

  val = 0.0;
  noise = 0.0;

  /* Determine SNR for current sub-image */
  for (ii=xs; ii<(xl+xs); ii++ ) {

    for (jj=ys; jj<(yl+ys); jj++) {
      long pix;
      if ( ii >= xlen ) {
    continue; 
      }
      if ( jj >= ylen ) {
    continue; 
      }
    pix = ii+(jj*xlen);

      val += data[pix];
      noise += ( derr[pix] * derr[pix] );
    }
  }
  locsnr = val / sqrt(noise);
  *oval = val;

  return locsnr;
}




/* Recursive binning routine */

void abin_rec ( float *data,   /* i: data array */
        float *derr,   /* i: error array */
        float *outd,   /* o: output array */
        float *outa,   /* o: output area array */
        unsigned long *mask, /* o: output mask */
        float *osnr,   /* o: output SNR */
        long   xlen,   /* i: length of x-axis (full img) */
        long   ylen,   /* i: length of y-axis (full img) */
        long   xs,     /* i: start of x-axis (sub img) */
        long   ys,     /* i: start of y-axis (sub img) */
        long   xl,     /* i: length of x-axis (sub img) */
        long   yl,     /* i: length of y-axis (sub img) */
        float  snr )   /* i: SNR threshold */
{
  
  static unsigned long mask_no; /* keep as static to avoid yet another
                   param to function; gets updated for
                   each recursive call */


  /* Determine SNR for current sub-image */

  float ll, lr, ul, ur, oval;
  ll = get_snr( data, derr, xlen, ylen, 
          xs, ys, FLOOR(xl/2.0), FLOOR(yl/2.0), &oval ); /* low-left */
  lr = get_snr( data, derr, xlen, ylen, 
          xs+FLOOR(xl/2.0), ys, CEIL(xl/2.0), FLOOR(yl/2.0), &oval ); /* low-rite*/
  ul = get_snr( data, derr, xlen, ylen, 
          xs, ys+FLOOR(yl/2.0), FLOOR(xl/2.0), CEIL(yl/2.0), &oval); /* up-left */
  ur = get_snr( data, derr, xlen, ylen, 
          xs+FLOOR(xl/2.0), ys+FLOOR(yl/2.0), CEIL(xl/2.0), CEIL(yl/2.0), &oval ); /* up-rite */

  /* Algorithm is basically if the SNR of all subimages is > limit, then
     split into 4 and repeat. */


  short check = ( (ll > snr ) && (lr > snr) && (ul > snr) && (ur > snr ));         

  if (( check )&&(xl>1)&&(yl>1)) {
    /* Enter recursion */

    /* need to use floor() and ceil() because input image may
       not be square, or 2^n.  This will bias the left-upper image w/ 1 more 
       pixel per bin; but that's a limit of not using square images.
       The alternative is only use square smallest sub-image or pad image to 2**N.*/
    
    abin_rec( data, derr, outd, outa, mask, osnr, xlen, ylen, 
          xs, ys, FLOOR(xl/2.0), FLOOR(yl/2.0), snr ); /* low-left */
    abin_rec( data, derr, outd, outa, mask, osnr, xlen, ylen, 
          xs+FLOOR(xl/2.0), ys, CEIL(xl/2.0), FLOOR(yl/2.0), snr ); /* low-rite*/
    abin_rec( data, derr, outd, outa, mask, osnr, xlen, ylen, 
          xs, ys+FLOOR(yl/2.0), FLOOR(xl/2.0), CEIL(yl/2.0), snr ); /* up-left */
    abin_rec( data, derr, outd, outa, mask, osnr, xlen, ylen, 
          xs+FLOOR(xl/2.0), ys+FLOOR(yl/2.0), CEIL(xl/2.0), CEIL(yl/2.0), 
          snr ); /* up-rite */
    
    
  } else { /* else do bin */
    double regx,regy;
    double regr[2];
    double logcoor[2];
    double physical[2];
    float locsnr;
    long ii,jj;
    float val;

    locsnr = get_snr( data, derr, xlen, ylen, xs, ys, xl, yl, &val );
    
    double area = ( xl * yl );
    val /= area;

    mask_no += 1; /* statically increases per bin */


    /* store output values */
    for (ii=xs; ii<xs+xl; ii++ ) {
      for (jj=ys; jj<ys+yl; jj++) {
    long pix;
    
    if ( ii >= xlen ) { /* shouldn't be needed anymore */
      continue; 
    }
    if ( jj >= ylen ) { /* shouldn't be needed anymore */
      continue; 
    }
    
    pix = ii+(jj*xlen);
    outd[pix] = val;
    outa[pix] = area;
    mask[pix] = mask_no;
    osnr[pix] = locsnr;
      }
      
    }


    regx = xs+(xl/2.0) + 1; /* -> to phys */
    regy = ys+(yl/2.0) + 1; /* -> to phys */

    logcoor[0]=regx;
    logcoor[1]=regy;
    dmCoordCalc_d( PhysDes, logcoor, physical );
    regx=physical[0];
    regy=physical[1];

    regr[0] = xl * fabs( Cdelt[0]);       /* -> to phys */
    regr[1] = yl * fabs( Cdelt[1]);       /* -> to phys */

    regAppendShape( maskRegion, "Box", 1, 1, &regx, &regy,
            1, regr, NULL, 0, 0 );


  } /* end else */

  return;
}





/* Main routine; does all the work of a quad-tree adaptive binning routine*/
int abin ()
{

  char infile[DS_SZ_FNAME];
  char errimg[DS_SZ_FNAME];

  char outfile[DS_SZ_FNAME];
  char areafile[DS_SZ_FNAME];
  char maskfile[DS_SZ_FNAME];
  char snrfile[DS_SZ_FNAME];
  float snr;
  short clobber;

  long nAxes;
  long *lAxes;
  float *data;
  float *derr;
  float *outd;
  float *outa;
  float *osnr;
  unsigned long *mask;
  long npix;

  dmBlock *inBlock;
  dmBlock *outBlock;
  dmDescriptor *inDes;
  dmDescriptor *outDes;

  /* Read in all the data */
  clgetstr( "infile", infile, DS_SZ_FNAME );
  clgetstr( "outfile", outfile, DS_SZ_FNAME );
  snr = clgetd( "snr" );
  clgetstr( "inerrfile",   errimg, DS_SZ_FNAME );
  clgetstr( "outmaskfile", maskfile, DS_SZ_FNAME );
  clgetstr( "outsnrfile",  snrfile,  DS_SZ_FNAME );
  clgetstr( "outareafile", areafile, DS_SZ_FNAME );
  clobber = clgetb( "clobber" );


  /* Go ahead and take care of the autonaming */
  ds_autoname( infile, outfile, "abinimg", DS_SZ_FNAME );
  ds_autoname( outfile, maskfile, "maskimg", DS_SZ_FNAME );
  ds_autoname( outfile, snrfile, "snrimg", DS_SZ_FNAME );
  ds_autoname( outfile, areafile, "areaimg", DS_SZ_FNAME );

  /* Read the data */
  inBlock = dmImageOpen( infile );
  if ( !inBlock ) {
    err_msg("ERROR: Could not open infile='%s'\n", infile );
  }
  inDes = dmImageGetDataDescriptor(inBlock);
  nAxes = dmGetArrayDimensions( inDes, &lAxes );
  if ( nAxes != 2 ) {
    err_msg("ERROR: Image must have 2 axes\n" );
    return(-1);
  }
  npix = ( lAxes[0]*lAxes[1]);
  if (npix ==0 ) {
    err_msg("ERROR: Image is empty (one axis is 0 length)\n");
    return(-1);
  }

  PhysDes = dmArrayGetAxisGroup( inDes, 1 );
  dmCoordGetTransform_d( PhysDes, Crpix, Crval, Cdelt, 2 );


  /* Allocate memory for the products */
  data = (float*)calloc(npix,sizeof(float));
  derr = (float*)calloc(npix,sizeof(float));
  outd = (float*)calloc(npix,sizeof(float));
  outa = (float*)calloc(npix,sizeof(float));
  osnr = (float*)calloc(npix,sizeof(float));
  mask = (unsigned long*)calloc(npix,sizeof(unsigned long));


  get_data( inDes, npix, data );
  

  /* Read Error Image */
  if ( ( strlen(errimg) == 0 ) ||
       ( ds_strcmp_cis(errimg,"none" ) == 0 ) ) {
    long jj;
    
    for (jj=0;jj<npix;jj++) { 
      derr[jj] = sqrt( data[jj] );
    }

  } else {
    long enAxes;
    long *elAxes;
    dmDescriptor *errDs;
    dmBlock *erBlock;

    erBlock = dmImageOpen( errimg );
    if ( erBlock == NULL ) {
      err_msg("ERROR: Could not open file '%s'\n", errimg);
      return(-1);
    }
    errDs = dmImageGetDataDescriptor(erBlock );
    enAxes = dmGetArrayDimensions( errDs, &elAxes );
    if ( (enAxes != 2 ) || 
     ( elAxes[0] != lAxes[0] ) ||
     ( elAxes[1] != lAxes[1] )    ) {
      err_msg("ERROR: Error image must be 2D image with non-zero axes\n");
      return(-1);
    }
    get_data( errDs, npix, derr );
    dmImageClose( erBlock );
  }


  /* Start Algorithm */

  maskRegion = regCreateEmptyRegion();

  abin_rec( data, derr, outd, outa, mask, osnr, lAxes[0], lAxes[1], 
    0, 0, lAxes[0], lAxes[1], snr );

  /* Write out file */

  if ( ds_clobber( outfile, clobber, NULL) == 0 ) {
    outBlock = dmImageCreate(outfile, dmFLOAT, lAxes, nAxes );
    if ( outBlock == NULL ) {
      err_msg("ERROR: Could not create output '%s'\n", outfile);
    }
    outDes = dmImageGetDataDescriptor( outBlock );
    dmBlockCopy( inBlock, outBlock, "HEADER"); 
    ds_copy_full_header( inBlock, outBlock, "dmnautilus", 0 );
    put_param_hist_info( outBlock, "dmnautilus", NULL, 0 );
    dmBlockCopyWCS( inBlock, outBlock);
    dmSetArray_f( outDes, outd, npix );
    dmImageClose( outBlock );
  } else {
    return(-1);
  }

  if ( (strlen(areafile)>0) && (ds_strcmp_cis(areafile,"none")!=0) ) {
    if ( ds_clobber( areafile, clobber, NULL) == 0 ) {
      outBlock = dmImageCreate(areafile, dmFLOAT, lAxes, nAxes );
      if ( outBlock == NULL ) {
    err_msg("ERROR: Could not create output '%s'\n", areafile);
      }
      outDes = dmImageGetDataDescriptor( outBlock );
      dmBlockCopy( inBlock, outBlock, "HEADER");
      ds_copy_full_header( inBlock, outBlock, "dmnautilus", 0 );
      put_param_hist_info( outBlock, "dmnautilus", NULL, 0 );
      dmBlockCopyWCS( inBlock, outBlock);
      dmSetArray_f( outDes, outa, npix );
      dmImageClose( outBlock );
    } else {
      return(-1);
    }
  }

  if ( (strlen(snrfile)>0) && (ds_strcmp_cis(snrfile,"none")!=0) ) {
    if ( ds_clobber( snrfile, clobber, NULL) == 0 ) {
      outBlock = dmImageCreate(snrfile, dmFLOAT, lAxes, nAxes );
      if ( outBlock == NULL ) {
    err_msg("ERROR: Could not create output '%s'\n", snrfile);
      }
      outDes = dmImageGetDataDescriptor( outBlock );
      dmBlockCopy( inBlock, outBlock, "HEADER");
      ds_copy_full_header( inBlock, outBlock, "dmnautilus", 0 );
      put_param_hist_info( outBlock, "dmnautilus", NULL, 0 );
      dmBlockCopyWCS( inBlock, outBlock);
      dmSetArray_f( outDes, osnr, npix );
      dmImageClose( outBlock );
    } else {
      return(-1);
    }
  }


  if ( (strlen(maskfile)>0) && (ds_strcmp_cis(maskfile,"none")!=0) ) {
    if ( ds_clobber( maskfile, clobber, NULL) == 0 ) {
      outBlock = dmImageCreate(maskfile, dmULONG, lAxes, nAxes );
      if ( outBlock == NULL ) {
    err_msg("ERROR: Could not create output '%s'\n", maskfile);
      }
      outDes = dmImageGetDataDescriptor( outBlock );
      dmBlockCopy( inBlock, outBlock, "HEADER");
      ds_copy_full_header( inBlock, outBlock, "dmnautilus", 0 );
      put_param_hist_info( outBlock, "dmnautilus", NULL, 0 );
      dmBlockCopyWCS( inBlock, outBlock);
      dmSetArray_ul( outDes, mask, npix );

      dmBlockClose( dmTableWriteRegion( dmBlockGetDataset( outBlock ),
              "REGION", NULL, maskRegion ));

      dmImageClose( outBlock );
    } else {
      return(-1);
    }
  }
  
  dmImageClose( inBlock ); /* Must keep open to do all the wcs/hdr
                  copies */
  
  free(data);
  free(derr);
  free(outd);
  free(outa);
  free(osnr);
  free(mask);


  return(0);

}
