#============================================================================================
# definitions
#============================================================================================
DEV_DIR = ../../CosmoRec.v2.0b_prep/Development

MAKE_COMMAND = make -f Makefile.in DEV_DIR=$(DEV_DIR)

#============================================================================================
# rules
#============================================================================================
all: 
	$(MAKE_COMMAND) all

lib: 
	$(MAKE_COMMAND) lib

bin: 
	$(MAKE_COMMAND) bin

pub:
	make -f Makefile.pub all

pubtar:
	make -f Makefile.pub tarball

plug:
	make -f Makefile.plug all

plugtar:
	make -f Makefile.plug tarball

#============================================================================================
# rules to clean up
#============================================================================================
clean:
	rm -f *.o

cleanall:
	rm -f *.o *~ Helium libHelium.a

cleanallDEV:
	make cleanall
	rm -f $(DEV_DIR)/Definitions/*.o $(DEV_DIR)/Definitions/*~
	rm -f $(DEV_DIR)/Simple_routines/*.o $(DEV_DIR)/Simple_routines/*~
	rm -f $(DEV_DIR)/Hydrogenic/*.o $(DEV_DIR)/Hydrogenic/*~
	rm -f $(DEV_DIR)/Line_profiles/*.o $(DEV_DIR)/Line_profiles/*~
	rm -f $(DEV_DIR)/Integration/*.o $(DEV_DIR)/Integration/*~
	rm -f $(DEV_DIR)/Hydrogenic/*.o $(DEV_DIR)/Hydrogenic/*~

tidy:
	make cleanallDEV

wipeDS:
	find . -type f -name \.DS_Store -print | xargs rm

#============================================================================================
#============================================================================================
