Module 		:=RICH

include	lib$(Module).pkg

Sources   	:=$(SRCS)

OutputDir	:=$(KIR_DIR)/$(Module)
SoLib		:=$(KIR_DIR)/lib$(Module).so

CxxFlags	:=-O -g -W -Wall -fPIC -pipe -fmessage-length=0 -I$(ROOTSYS)/include -I$(ALICE_ROOT)/include

ifdef ALIVERBOSE
	Mute		:=
else
	Mute		:=@        
endif        

Headers          := $(Sources:.cxx=.h)  $(Module)LinkDef.h

DictSource    := $(OutputDir)/G__$(Module).cxx
DictHeader    := $(DictSource:.cxx=.h)
DictObject    := $(DictSource:.cxx=.o)


Objects         := $(patsubst %.cxx,$(OutputDir)/%.o,$(Sources)) $(DictObject)


DepFile		:= $(OutputDir)/$(Module).d

##### TARGETS #####

$(SoLib):	 $(Objects)
	@echo "Creating $@"
	$(Mute)g++ -Wl,-soname,lib$(Module).so -shared -O -g $^ -o $@	

$(DepFile):	$(Headers)
	@[ -d $(OutputDir) ] || mkdir -p $(OutputDir)
	@touch $(DepFile)
	@echo "Generating dependency $@"
	$(Mute)rmkdepend -f$(DepFile) -p$(OutputDir)/ -- $(CxxFlags) -- $(Sources) 2>/dev/null

$(DictSource):  $(Headers)
	@echo "Generating $@"
	$(Mute)rootcint -f $@ -c $(filter -I%,$(CxxFlags)) $^

$(DictObject) : $(DictSource)
	@echo "Compiling $^"        
	$(Mute)g++ $(CxxFlags) -I. -c $^ -o $@

show:
	@echo "Headers: $(Headers)"
	@echo "Sources: $(Sources)"
	@echo "Depend : $(DepFile)"
	@echo "Objects: $(Objects)"
	@echo "Library: $(SoLib)"

spec:	$(Sources)
	@echo $^
	@echo $@
                
clean:
	@echo "Cleaning..."
	$(Mute)rm -rf $(Objects) $(DictSource) $(DictHeader) $(SoLib) $(DepFile)		  

############################ cxx rule #########################################
$(OutputDir)/%.o : %.cxx
	@echo $*.cxx
	$(Mute)g++ $(CxxFlags) -c $*.cxx -o $(OutputDir)/$*.o
############################ Dependencies #####################################

-include $(DepFile) 
