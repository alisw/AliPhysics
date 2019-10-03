TASK = TOFSpectraPbPb
JDL = $(TASK).jdl
EXE = $(TASK).sh
TGZ = $(TASK).tgz
DIR = /alice/cern.ch/user/r/rpreghen/$(TASK)
SE = ALICE::CNAF::SE

FILES = AliAnalysisTask$(TASK).h \
	AliAnalysisTask$(TASK).cxx \
	AliAnalysisEvent.h \
	AliAnalysisEvent.cxx \
	AliAnalysisTrack.h \
	AliAnalysisTrack.cxx \
	AliAnalysisParticle.h \
	AliAnalysisParticle.cxx \
	AddAnalysisTask$(TASK).C \
	SteerAnalysisTask$(TASK).C

all: $(TGZ)

$(TGZ): $(FILES)
	@echo "creating archive: "$(TGZ)
	@tar zcvfh $(TGZ) $(FILES)

install: $(TGZ)
	@echo "installing archive: "$(TGZ)
	@alien_cp -n $(TGZ) alien://$(DIR)/$(TGZ)@$(SE)
	@echo "installing JDL: "$(JDL)
	@alien_cp -n $(JDL) alien://$(DIR)/$(JDL)@$(SE)
	@echo "installing executable: "$(EXE)
	@alien_cp -n $(EXE) alien://$(DIR)/$(EXE)@$(SE)

uninstall:
	@alien_erase $(DIR)/$(TGZ)

clean:
	@rm -rf *~ $(TGZ)
