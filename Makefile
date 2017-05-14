MAKE = make		#change this line if you are using a different GNU make software

dirFASTAtoRF = ./src/FASTAtoRF
dirRCOpt = ./src/RCOpt
dirRC = ./src/RC

all: MK_dir CC_FASTAtoRF CC_RCOpt CC_RC RM_objectFiles

MK_dir:
	mkdir -p ./bin

CC_FASTAtoRF: $(dirFASTAtoRF)/Makefile
	$(MAKE) -C $(dirFASTAtoRF)
	
CC_RCOpt: $(dirRCOpt)/Makefile
	$(MAKE) -C $(dirRCOpt)

CC_RC: $(dirRC)/Makefile
	$(MAKE) -C $(dirRC)

RM_objectFiles:
	rm -f $(dirFASTAtoRF)/*.o $(dirRCOpt)/*.o $(dirRC)/*.o 

clean:
	rm -f $(dirFASTAtoRF)/*.o $(dirRCOpt)/*.o $(dirRC)/*.o ./bin/*
