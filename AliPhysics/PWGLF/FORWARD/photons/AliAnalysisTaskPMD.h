#ifndef AliAnalysisTaskPMD_cxx
#define AliAnalysisTaskPMD_cxx

// AnalysisTask For PMD 
// Authors: Sudipan De, Subhash Singha

class TH1F;
class TH2F;
class AliESDEvent;
class AliESDPmdTrack;
class AliESDVertex;

#include <AliAnalysisTaskSE.h>

class AliAnalysisTaskPMD : public AliAnalysisTaskSE {
 public:
 AliAnalysisTaskPMD() : AliAnalysisTaskSE(), 
    fESD(0), 
    fOutputList(0), 
    fHistTotEvent(0),
    fHistTotEventAfterPhySel(0),
    fHistTotEventAfterVtx(0),
    fHistVtxZ(0),
    fHistXYPre(0),
    fHistEta(0),
    fHistEta1(0),
    fHistMultMeas(0),
    fHistMultMeas1(0)
      {
	for(Int_t i=0; i<10; i++){
	  fHistMultMeasEtaBinA[i] = 0;
	  fHistMultMeasEtaBinA1[i] = 0;
	}
      }
  AliAnalysisTaskPMD(const char *name);
  virtual ~AliAnalysisTaskPMD() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
 private:
  AliESDEvent *fESD;    //! ESD object
  TList *fOutputList; //! Output list
  TH1F *fHistTotEvent; //total event
  TH1F *fHistTotEventAfterPhySel; //# event after physics selection
  TH1F *fHistTotEventAfterVtx; //# event after vertex cut
  TH1F *fHistVtxZ;//z vertex distribution
  TH2F *fHistXYPre;//2d scatter plot pre
  TH1F *fHistEta; // eta distribution in PMD coverage
  TH1F *fHistEta1; // eta distribution in PMD coverage
  TH1F *fHistMultMeas;//measured multiplicity (2.3-3.9)
  TH1F *fHistMultMeas1;//measured multiplicity (2.3-3.9)
  
  TH1F *fHistMultMeasEtaBinA[10];//meas. mult. dist. for diff. eta bins
  TH1F *fHistMultMeasEtaBinA1[10];//meas. mult. dist. for diff. eta bins
  
  AliAnalysisTaskPMD(const AliAnalysisTaskPMD&); // not implemented
  AliAnalysisTaskPMD& operator=(const AliAnalysisTaskPMD&); // not implemented
  
  ClassDef(AliAnalysisTaskPMD, 1); // example of analysis
};

#endif
