MAKE = make		#change this line if you are using a different GNU make software

dirFASTAtoRF = ./src/FASTAtoRF
dirRCOpt = ./src/RCOpt

all: MK_dir CC_FASTAtoRF CC_RCOpt RM_objectFiles

MK_dir:
	mkdir -p ./bin

CC_FASTAtoRF: $(dirFASTAtoRF)/Makefile
	$(MAKE) -C $(dirFASTAtoRF)
	
CC_RCOpt: $(dirRCOpt)/Makefile
	$(MAKE) -C $(dirRCOpt)

RM_objectFiles:
	rm -f $(dirFASTAtoRF)/*.o $(dirRCOpt)/*.o 

clean:
	rm -f $(dirFASTAtoRF)/*.o $(dirRCOpt)/*.o ./bin/*
