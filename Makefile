# Makefile to build class 'pan~' for Pure Data.
# Needs Makefile.pdlibbuilder as helper makefile for platform-dependent build
# settings and rules.

# library name
lib.name = dsf

# input source file (class name == source file basename)
class.sources = dsf.c
dsf.class.sources = dsf_pd.c

# all extra files to be included in binary distribution of the library
datafiles = 

# include Makefile.pdlibbuilder from submodule directory 'pd-lib-builder'
PDLIBBUILDER_DIR=.
include $(PDLIBBUILDER_DIR)/Makefile.pdlibbuilder
