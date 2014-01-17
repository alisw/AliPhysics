#ifndef AliAnalysisTaskPMDSim_cxx
#define AliAnalysisTaskPMDSim_cxx

// AnalysisTask For PMD
// Authors: Sudipan De, Subhash Singha

class TH1F;
class TH2F;
class AliESDEvent;
class AliESDPmdTrack;
class AliESDVertex;
class AliStack;
class AliHeader;
class AliGenEventHeader;
class TParticle;

#include <AliAnalysisTaskSE.h>

class AliAnalysisTaskPMDSim : public AliAnalysisTaskSE {
 public:
 AliAnalysisTaskPMDSim() : AliAnalysisTaskSE(), 
    fESD(0), 
    fOutputList(0), 
    fHistTotEvent(0),
    fHistTotEventAfterPhySel(0),
    fHistTotEventAfterVtx(0),
    fVtxZ(0),
    fHistXYPre(0),
    fHistEtaPhM(0),
    fHistEtaPhM1(0),
    fHistEtaT(0),
    fMultMeasured(0),
    fMultMeasured1(0),
    fMultTrue(0),
    fMultCorr(0),
    fMultCorr1(0) {
    for(Int_t i=0; i<10; i++){
      fHistMultMeasEtaBinA[i] = 0;
      fHistMultMeasEtaBinA1[i] = 0;
      fHistMultTrueEtaBinA[i] = 0;
      fHistMultCorrEtaBinA[i] = 0;
      fHistMultCorrEtaBinA1[i] = 0;
    }
  }
  AliAnalysisTaskPMDSim(const char *name);
  virtual ~AliAnalysisTaskPMDSim() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
 
 private:
  AliESDEvent *fESD;    //! ESD object
  TList *fOutputList; //! Output list
  TH1F *fHistTotEvent; //total event
  TH1F *fHistTotEventAfterPhySel; //total event after physel
  TH1F *fHistTotEventAfterVtx; //# event after vertex cut
  TH1F *fVtxZ;//Vertex Z
  TH2F *fHistXYPre;//2d scatter plot pre
  TH1F *fHistEtaPhM;
  TH1F *fHistEtaPhM1;
  TH1F *fHistEtaT;
  TH1F *fMultMeasured;
  TH1F *fMultMeasured1;
  TH1F *fMultTrue;
  TH2F *fMultCorr;
  TH2F *fMultCorr1;
  
  TH2F *fHistMultCorrEtaBinA[10];//mult. corr. for diff. eta bin
  TH2F *fHistMultCorrEtaBinA1[10];//mult. corr. for diff. eta bin
  TH1F *fHistMultTrueEtaBinA[10];//multTrue
  TH1F *fHistMultMeasEtaBinA[10];//meas. mult. dist. for diff. eta bins
  TH1F *fHistMultMeasEtaBinA1[10];//meas. mult. dist. for diff. eta bins
  
  AliAnalysisTaskPMDSim(const AliAnalysisTaskPMDSim&); // not implemented
  AliAnalysisTaskPMDSim& operator=(const AliAnalysisTaskPMDSim&); // not implemented
  
  ClassDef(AliAnalysisTaskPMDSim, 1); // example of analysis
};

#endif
