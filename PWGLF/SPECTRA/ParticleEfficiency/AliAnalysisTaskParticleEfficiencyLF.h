#ifndef ALIANALYSISTASKPARTICLEEFFICIENCYLF
#define ALIANALYSISTASKPARTICLEEFFICIENCYLF

#include "AliAnalysisTaskSE.h"

class AliESDtrackCuts;
class TList;
class TH1F;
class TH2F;

class AliAnalysisTaskParticleEfficiencyLF :
public AliAnalysisTaskSE
{

 public:

  AliAnalysisTaskParticleEfficiencyLF(const Char_t *partName); // default constructor
  ~AliAnalysisTaskParticleEfficiencyLF(); // default destructor

  void UserCreateOutputObjects(); // user create output objects
  void UserExec(Option_t *option); // user exec
  
 protected:

  AliAnalysisTaskParticleEfficiencyLF(const AliAnalysisTaskParticleEfficiencyLF &); // copy constructor
  AliAnalysisTaskParticleEfficiencyLF &operator=(const AliAnalysisTaskParticleEfficiencyLF &); // operator=

  Int_t fParticlePdgCode; // particle PDG code
  AliESDtrackCuts *fTrackCuts; // ESD track cuts

  TList *fHistoList; // histo list
  TH2F *fHistoEvents; // histo events
  TH2F *fHistoGenerated; // histo generated
  TH2F *fHistoReconstructed; // histo reconstructed
  TH2F *fHistoGeneratedDaughter[2]; // histo generated daughter
  TH2F *fHistoReconstructedDaughter[2]; // histo reconstructed daughter

  ClassDef(AliAnalysisTaskParticleEfficiencyLF, 1);

};

#endif /* ALIANALYSISTASKPARTICLEEFFICIENCYLF */
