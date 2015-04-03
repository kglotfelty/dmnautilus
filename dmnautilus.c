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

#include <cxcregion.h>
#include <dsnan.h>

#define FLOOR(x)  floor((x))
#define CEIL(x)   ceil((x))

regRegion *maskRegion;
dmDescriptor *PhysDes;


/* Okay, global variables are a bad idea, excpet when doing recursion and
 * we are passing the same data to each level ... more data gets pushed on the
 * stack and causes crashes.  Also using global variables can make it
 * run much faster so we'll bite the bullet here. */
float *GlobalData;     /* i: data array */
float *GlobalDErr;     /* i: error array */
float *GlobalOutData;  /* o: output array */
float *GlobalOutArea;  /* o: output area array */
float *GlobalOutSNR;   /* o: output SNR */
unsigned long *GlobalMask; /* o: output mask */
long GlobalXLen;        /* i: length of x-axis (full img) */
long GlobalYLen;        /* i: length of y-axis (full img) */
long GlobalSNRThresh;   /* i: SNR threshold */
enum { ZERO_ABOVE=0, ONE_ABOVE, TWO_ABOVE, THREE_ABOVE, ALL_ABOVE} GlobalSplitCriteria;


/* TODO:  Replace i/o with dmimgio routines incl null/nan/indef checking 
#include "dmimgio.h"

*/




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
        long   xs,     /* i: start of x-axis (sub img) */
        long   ys,     /* i: start of y-axis (sub img) */
        long   xl,     /* i: length of x-axis (sub img) */
        long   yl,      /* i: length of y-axis (sub img) */
        float  *oval,   /* o: sum of pixel values */
        long   *area    /* o: number of pixels */
  )
{
  
  float val;
  float noise;
  float locsnr;
  long ii,jj;

  val = 0.0;
  noise = 0.0;
  *area = 0;

  /* Determine SNR for current sub-image */
  for (ii=xs; ii<(xl+xs); ii++ ) {
      if ( ii >= GlobalXLen ) {
        continue; 
      }
    for (jj=ys; jj<(yl+ys); jj++) {
      long pix;
      if ( jj >= GlobalYLen ) {
        continue; 
      }
      pix = ii+(jj*GlobalXLen);


      /* TODO: check for valid pixels */

      val += GlobalData[pix];
      noise += ( GlobalDErr[pix] * GlobalDErr[pix] );
      *area += 1;
    }
  }
  locsnr = val / sqrt(noise);
  *oval = val;

  return locsnr;
}




/* Recursive binning routine */
void abin_rec ( 
        long   xs,     /* i: start of x-axis (sub img) */
        long   ys,     /* i: start of y-axis (sub img) */
        long   xl,     /* i: length of x-axis (sub img) */
        long   yl      /* i: length of y-axis (sub img) */
       )   
{
  
  static unsigned long mask_no; /* keep as static to avoid yet another
                   param to function; gets updated for
                   each recursive call */

  
  short check = 0;


  if ( ZERO_ABOVE == GlobalSplitCriteria ) {
      /* This is the original method -- if the current block is above SNR
       * then split it. */
      float at, oval;
      long npix;  
      at = get_snr( xs, ys, xl, yl , &oval, &npix );
      check = ( at > GlobalSNRThresh );

  } else {
      /* Determine SNR for current sub-image */

      /* This new method will only split the 2x2 if ALL 4 of the subimages
       * exceed the SNR threshold 
       */
      
      float ll, lr, ul, ur;
      float oval_ll, oval_lr, oval_ul, oval_ur;
      long npix;  
      short ill, ilr, iul, iur;

      ll = get_snr( xs, ys, FLOOR(xl/2.0), FLOOR(yl/2.0), &oval_ll, &npix ); /* low-left */
      lr = get_snr( xs+FLOOR(xl/2.0), ys, CEIL(xl/2.0), FLOOR(yl/2.0), &oval_lr, &npix ); /* low-rite*/
      ul = get_snr( xs, ys+FLOOR(yl/2.0), FLOOR(xl/2.0), CEIL(yl/2.0), &oval_ul, &npix); /* up-left */
      ur = get_snr( xs+FLOOR(xl/2.0), ys+FLOOR(yl/2.0), CEIL(xl/2.0), CEIL(yl/2.0), &oval_ur, &npix ); /* up-rite */

      ill = ( ll >= GlobalSNRThresh);
      ilr = ( lr >= GlobalSNRThresh);
      iul = ( ul >= GlobalSNRThresh);
      iur = ( ur >= GlobalSNRThresh);

      if ( ONE_ABOVE == GlobalSplitCriteria ) {
          /* If any one of the sub images is above snr, then split */
          check = ( ill + ilr + iul + iur  );            
      } else if ( TWO_ABOVE == GlobalSplitCriteria ) {
          /* If two, then they have to be side-by side, not diagonal */
         check = ( (ill && ilr ) || 
                   (ilr && iur ) ||
                   (iur && iul ) ||
                   (iul && ill )
                  ); 
      } else if ( THREE_ABOVE == GlobalSplitCriteria ) {
          check = ( ill + ilr + iul + iur ) >= 3 ? 1 : 0;

      } else if ( ALL_ABOVE == GlobalSplitCriteria ) {
          check = ( ill + ilr + iul + iur ) == 4 ? 1 : 0;
      } else {
          err_msg("This should not have happened, something's amiss");
          return;
      }
      
  } // end else 

  
  if (( check )&&(xl>1)&&(yl>1)) {
    /* Enter recursion */

    /* need to use floor() and ceil() because input image may
       not be square, or 2^n.  This will bias the left-upper image w/ 1 more 
       pixel per bin; but that's a limit of not using square images.
       The alternative is only use square smallest sub-image or pad image to 2**N.*/
        abin_rec( xs, ys, FLOOR(xl/2.0), FLOOR(yl/2.0) ); /* low-left */
        abin_rec( xs+FLOOR(xl/2.0), ys, CEIL(xl/2.0), FLOOR(yl/2.0) ); /* low-rite*/
        abin_rec( xs, ys+FLOOR(yl/2.0), FLOOR(xl/2.0), CEIL(yl/2.0) ); /* up-left */
        abin_rec( xs+FLOOR(xl/2.0), ys+FLOOR(yl/2.0), CEIL(xl/2.0), CEIL(yl/2.0)); /* up-rite */
        return; 
    }
      
    float locsnr;
    long ii,jj;
    float val;
    long area;

    locsnr = get_snr( xs, ys, xl, yl, &val, &area );

    val /= area;

    mask_no += 1; /* statically increases per bin */


    /* store output values */
    for (ii=xs; ii<xs+xl; ii++ ) {
      if ( ii >= GlobalXLen ) { /* shouldn't be needed anymore */
      continue; 
      }
      for (jj=ys; jj<ys+yl; jj++) {
      long pix;

      if ( jj >= GlobalYLen ) { /* shouldn't be needed anymore */
        continue; 
      }

    pix = ii+(jj*GlobalXLen);
    GlobalOutData[pix] = val;
    GlobalOutArea[pix] = area;
    GlobalMask[pix] = mask_no;
    GlobalOutSNR[pix] = locsnr;
      }
      
    }


    /*   
     * The original hacky way to use the Cdelt[]'s doesn't work. There are 
     * cases where the regions overlap which is bad.
     * 
     * So instead of using a box, we use a rectangle since the
     * edge points are explicitly specified.
     * 
    */

    double regx[2], regy[2];
    double logcorr[2];
    double phycorr[2];
    
    logcorr[0] = xs;
    logcorr[1] = ys;
    dmCoordCalc_d( PhysDes, logcorr, phycorr);
    regx[0] = phycorr[0];
    regy[0] = phycorr[1];
    
    logcorr[0] = xs+xl;
    logcorr[1] = ys+yl;
    dmCoordCalc_d( PhysDes, logcorr, phycorr);
    regx[1] = phycorr[0];
    regy[1] = phycorr[1];
    
    regAppendShape( maskRegion, "Rectangle", 1, 1, regx, regy,
            1, NULL, NULL, 0, 0 );

    return;
}



int load_error_image( char *errimg ) {
    
  /* Read Error Image */
  unsigned long npix = GlobalXLen*GlobalYLen;

  if ( ( strlen(errimg) == 0 ) ||
       ( ds_strcmp_cis(errimg,"none" ) == 0 ) ) {
    long jj;    

    for (jj=0;jj<npix;jj++) { 
      GlobalDErr[jj] =  sqrt(GlobalData[jj]);  // assumes Gaussian stats
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
     ( elAxes[0] != GlobalXLen ) ||
     ( elAxes[1] != GlobalYLen )    ) {
      err_msg("ERROR: Error image must be 2D image with non-zero axes\n");
      return(-1);
    }
    get_data( errDs, npix, GlobalDErr );
    dmImageClose( erBlock );
  }

  return 0;
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
  short method;
  short clobber;

  long nAxes;
  long *lAxes;
  long npix;

  dmBlock *inBlock;
  dmBlock *outBlock;
  dmDescriptor *inDes;
  dmDescriptor *outDes;

  /* Read in all the data */
  clgetstr( "infile", infile, DS_SZ_FNAME );
  clgetstr( "outfile", outfile, DS_SZ_FNAME );
  GlobalSNRThresh = clgetd( "snr" );
  method = clgeti("method");
  clgetstr( "inerrfile",   errimg, DS_SZ_FNAME );
  clgetstr( "outmaskfile", maskfile, DS_SZ_FNAME );
  clgetstr( "outsnrfile",  snrfile,  DS_SZ_FNAME );
  clgetstr( "outareafile", areafile, DS_SZ_FNAME );
  clobber = clgetb( "clobber" );

  switch (method) 
  {
    case 0: GlobalSplitCriteria = ZERO_ABOVE; break;
    case 1: GlobalSplitCriteria = ONE_ABOVE; break;
    case 2: GlobalSplitCriteria = TWO_ABOVE; break;
    case 3: GlobalSplitCriteria = THREE_ABOVE; break;
    case 4: GlobalSplitCriteria = ALL_ABOVE; break;
    default:
      err_msg("Invalid method parameter value");
      break;
  };


  /* Go ahead and take care of the autonaming */
  ds_autoname( infile, outfile, "abinimg", DS_SZ_FNAME );
  ds_autoname( outfile, maskfile, "maskimg", DS_SZ_FNAME );
  ds_autoname( outfile, snrfile, "snrimg", DS_SZ_FNAME );
  ds_autoname( outfile, areafile, "areaimg", DS_SZ_FNAME );


  /* Read the data */
  /* TODO: Replace with dmimgio routines */

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
  GlobalXLen = lAxes[0];
  GlobalYLen = lAxes[1];


  /* Allocate memory for the products */
  GlobalData = (float*)calloc(npix,sizeof(float));
  GlobalDErr = (float*)calloc(npix,sizeof(float));
  GlobalOutData = (float*)calloc(npix,sizeof(float));
  GlobalOutArea = (float*)calloc(npix,sizeof(float));
  GlobalOutSNR = (float*)calloc(npix,sizeof(float));
  GlobalMask = (unsigned long*)calloc(npix,sizeof(unsigned long));
  maskRegion = regCreateEmptyRegion();


  get_data( inDes, npix, GlobalData );


  if ( 0 != load_error_image( errimg ) ) {
        return -1;
  }
  
  /* Start Algorithm */

  abin_rec( 0, 0, GlobalXLen, GlobalYLen);


  /* Write out files -- NB: mask file has different datatypes and different extensions */

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
    dmSetArray_f( outDes, GlobalOutData, npix );
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
      dmSetArray_f( outDes, GlobalOutArea, npix );
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
      dmSetArray_f( outDes, GlobalOutSNR, npix );
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
      dmSetArray_ul( outDes, GlobalMask, npix );

      dmBlockClose( dmTableWriteRegion( dmBlockGetDataset( outBlock ),
              "REGION", NULL, maskRegion ));

      dmImageClose( outBlock );
    } else {
      return(-1);
    }
  }
  
  dmImageClose( inBlock ); /* Must keep open to do all the wcs/hdr
                  copies */
  
  /* make valgrind happy */
  free(GlobalData);
  free(GlobalDErr);
  free(GlobalOutData);
  free(GlobalOutArea);
  free(GlobalOutSNR);
  free(GlobalMask);


  return(0);

}
