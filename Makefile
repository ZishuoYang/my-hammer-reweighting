BINPATH	:=	bin
VPATH	:=	$(BINPATH)

# Compiler settings
COMPILER	:=	$(shell root-config --cxx)
CXXFLAGS	:=	$(shell root-config --cflags) -Iinclude
LINKFLAGS	:=	$(shell root-config --libs)
ADDLINKFLAGS	:=	-lHammerTools -lHammerBase -lHammerCore -lFormFactors -lAmplitudes -lRates


###########
# General #
###########

tools: Bc2JpsiMuNu Bc2JpsiMuNu_SaveWeights

.PHONY: clean
clean:
	@rm -rf ./bin/*


####################
# Generic patterns #
####################

.SECONDARY:

# Reweighters with HAMMER
%: %.cc
	$(COMPILER) $(CXXFLAGS) -o $(BINPATH)/$@ $< $(LINKFLAGS) $(ADDLINKFLAGS)
