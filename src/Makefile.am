tool = dmnautilus

binexedir = $(prefix)/binexe
paramdir = $(prefix)/param
ahelpdir = $(prefix)/share/doc/xml

binexe_PROGRAMS = $(tool)

dmnautilus_SOURCES = dmnautilus.c t_dmnautilus.c
dmnautilus_CPPFLAGS = $(CIAO_CFLAGS)
dmnautilus_LDADD = $(CIAO_LIBS) -ldmimgio
dmnautilus_LINK = $(CXX) -o $@ -Wl,-rpath,$(prefix)/lib -Wl,-rpath,$(prefix)/ots/lib 

# problems with libstdc++ and cxcparam on Linux :(
if LINUX
dmnautilus_LDADD += -L$(prefix)/ots/lib -lstdc++
endif


dist_param_DATA = $(tool).par
dist_ahelp_DATA = $(tool).xml

install-data-hook:
	chmod a-w $(paramdir)/$(dist_param_DATA)
	chmod a-w $(ahelpdir)/$(dist_ahelp_DATA)

install-exec-local:
	$(RM) $(binexedir)/$(tool)
	$(RM) $(paramdir)/$(dist_param_DATA)
	$(RM) $(prefix)/bin/$(tool)
	$(RM) $(ahelpdir)/$(dist_ahelp_DATA)

