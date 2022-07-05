#ifndef AliPMDpPbAnalysisTaskSim_cxx
#define AliPMDpPbAnalysisTaskSim_cxx

//PMD pA data analysis
//Author: Abhi Modak (abhi.modak@cern.ch)
//Co-Authors: Sudipan De (Sudipan.De@cern.ch) and Sidharth Kumar Prasad (sidharth.kumar.prasad@cern.ch)
//Taken the help from the code: alisw/AliPhysics/PWGLF/FORWARD/photons/AliAnalysisTaskPMD.h

class TH1F;
class TH2F;
class TString;
class TNtuple;
class TArrayF;
class AliESDEvent;
class AliESDPmdTrack;
class AliHeader;
class TParticle;
class AliMultSelection;
class AliStack;
class AliPPVsMultUtils;
class AliGenEventHeader;

#include "AliAnalysisTaskSE.h"

class AliPMDpPbAnalysisTaskSim : public AliAnalysisTaskSE {
 public:
 AliPMDpPbAnalysisTaskSim() : AliAnalysisTaskSE(),
    fESD(0), 
    fEsdV0(0),
    fOutputList(0),
    fCentEstimator(0),
    fTrigSel(0),
    fHistClsXYPre(0),
    fHistClsXYCpv(0),
    fHistADCPre(0),
    fHistADCCpv(0),
    fHistVtxZ(0),
    fHistTotEvent(0),
    fHistEtaCut1(0),
    fHistEtaCut2(0),
    fHistNcellCut1(0),
    fHistNcellCut2(0),
    fHistEtaChTrue(0),
    fHistEtaPhTrue(0),
    ntMeas1(0),
    ntMeas2(0),
    ntCent(0),
    ntCorr(0),
    ntTrue(0),
    ntDet1(0),
    ntDet2(0) {
    for(Int_t i=0; i<7; i++){
      fHistEtaChCentClass[i] = 0;
      fHistEtaGammaCentClass[i] = 0;
    }
  }

  AliPMDpPbAnalysisTaskSim(const char *name);
  virtual ~AliPMDpPbAnalysisTaskSim(){}

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);
  void SetTrigger(const char* trigger){
    fTrigSel = trigger;
  };
  void SetCentralityEstimator(const char* Estimator){
    fCentEstimator = Estimator;
  };
  
 private:
  AliESDEvent *fESD;
  AliESDVZERO *fEsdV0;
  TList       *fOutputList;
  TNtuple     *ntMeas1;
  TNtuple     *ntMeas2;
  TNtuple     *ntCent;
  TNtuple     *ntCorr;
  TNtuple     *ntTrue;
  TNtuple     *ntDet1;
  TNtuple     *ntDet2;
  TH1F        *fHistTotEvent;
  TH2F        *fHistClsXYPre;
  TH2F        *fHistClsXYCpv;
  TH1F        *fHistADCPre;
  TH1F        *fHistADCCpv;
  TH1F        *fHistVtxZ;
  TH1F        *fHistEtaCut1;
  TH1F        *fHistEtaCut2;
  TH1F        *fHistNcellCut1;
  TH1F        *fHistNcellCut2;
  TH1F        *fHistEtaChTrue;
  TH1F        *fHistEtaPhTrue;
  TH1F        *fHistEtaChCentClass[7];//For charged-particle
  TH1F        *fHistEtaGammaCentClass[7];//For photons
  TString     fCentEstimator;
  TString     fTrigSel;
  
    
  AliPMDpPbAnalysisTaskSim(const AliPMDpPbAnalysisTaskSim&);
  AliPMDpPbAnalysisTaskSim& operator = (const AliPMDpPbAnalysisTaskSim&);

  ClassDef(AliPMDpPbAnalysisTaskSim, 1);
};

#endif 
