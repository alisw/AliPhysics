Module 		:=RICH

include	lib$(Module).pkg

Sources   	:=$(SRCS)

RootTarget      :=$(shell root-config --arch)
OutputDir	:=$(KIR_DIR)/tgt_$(RootTarget)/$(Module)
SoLib		:=$(KIR_DIR)/tgt_$(RootTarget)/lib$(Module).so

ifeq ($(RootTarget),linuxicc)
 Compiler	:=icc
 CompilerOptions:=-O0 -fpstkchk -I$(shell root-config --incdir) -I$(ALICE_ROOT)/include
 SoOptions      := -Wl,-soname,lib$(Module).so -shared -O -g	
else
 Compiler       :=g++
 CompilerOptions:=-g -W -Wall -Werror -fPIC -pipe -fmessage-length=0 -Wno-long-long -pedantic-errors -ansi -I$(shell root-config --incdir) -I$(ALICE_ROOT)/include
 SoOptions      :=-O -g -shared -Wl,-soname,lib$(Module).so  
endif

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
	$(Mute)$(Compiler) $(SoOptions) $^ -o $@	

$(DepFile):	$(Headers)
	@[ -d $(OutputDir) ] || mkdir -p $(OutputDir)
	@touch $(DepFile)
	@echo "Generating dependency $@"
	$(Mute)rmkdepend -f$(DepFile) -p$(OutputDir)/ -- $(CxxFlags) -- $(Sources) 2>/dev/null

$(DictSource):  $(Headers)
	@echo "Generating $@"
	$(Mute)rootcint -f $@ -c $(filter -I%,$(CompilerOptions)) $^

$(DictObject) : $(DictSource)
	@echo "Compiling $^"        
	$(Mute)$(Compiler) $(CompilerOptions) -I. -c $^ -o $@

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
	$(Mute)$(Compiler) $(CompilerOptions) -c $*.cxx -o $(OutputDir)/$*.o
############################ Dependencies #####################################

-include $(DepFile) 
