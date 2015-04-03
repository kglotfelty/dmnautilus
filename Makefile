##############################################################################

#MK_TOP = ../../../..
MK_TOP  = /export/ciao_from_source/ciao-4.7/src/

include $(MK_TOP)/Makefile.master
include $(MK_TOP)/include/Makefile.scidev

EXEC              = dmnautilus++
LIB_FILES = 
PAR_FILES         = dmnautilus.par
INC_FILES         = 
XML_FILES         = dmnautilus.xml

SRCS	=           dmnautilus.c t_dmnautilus.c

LOCAL_LIBS = -L$(MK_TOP)/da/analysis/dmtools/dmimgio/ -ldmimgio
LOCAL_INC  = -I$(MK_TOP)/da/analysis/dmtools/dmimgio/

OBJS = $(SRCS:.c=.o)


MAKETEST_SCRIPT   = dmnautilus.t


include $(MK_TOP)/Makefile.all

#-----------------------------------------------------------------------
# 			MAKEFILE DEPENDENCIES	
#-----------------------------------------------------------------------



$(EXEC): $(OBJS)
	$(LINK)
	@echo


announce1:
	@echo "   /----------------------------------------------------------\ "
	@echo "   |                Building dmnautilus                       | "
	@echo "   \----------------------------------------------------------/ "

