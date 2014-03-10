#ifndef AliEbyEMultFluctuationTask_cxx
#define AliEbyEMultFluctuationTask_cxx

class TString;
class TH1F;
class TH2D;
class TH1I;
class TNtuple;
class AliESDTrack;
class AliAODEvent;
class AliAODHeader;
class AliVEvent;
class AliAODVertex;
class AliAODVZERO;
class AliAODTrack;
class AliAODTracklets;
#include "AliAnalysisTaskSE.h"

class AliEbyEMultFluctuationTask : public AliAnalysisTaskSE {
 public:
 AliEbyEMultFluctuationTask() : AliAnalysisTaskSE(), fAOD(0), fAODVertex(0),fHistNchPt(0),fHistNchEta(0),fHistNchEtaCent(0),fHistNchPhi(0),fHistDCAxy(0),fHistDCAz(0),fHistnclus(0),fHistchi2ndf(0),fHistchi2ndfvscs(0),fHistVz(0),fHistMultV0A(0),fHistMultV0C(0),fHistMultV0total(0),My_ntuple(0),fOutputList(0),fCentralityEstimator("V0M"),fCentralityBins20(kFALSE),fCentralityCounter(0),fEventCounter(0),histcounter(0){
		for(Int_t ibin=0;ibin<91;ibin++)
		{
			fMult[ibin]=NULL;
		}
		for(Int_t jbin=0;jbin<46;jbin++)
		{
		fMultTwo[jbin]=NULL;
		}
		for(Int_t kbin=0;kbin<15;kbin++)
		{
		fMultFive[kbin]=NULL;
		}
	
	  
	

	}
  AliEbyEMultFluctuationTask(const char *name);
  virtual ~AliEbyEMultFluctuationTask() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
	

  
  void SetCentralityEstimator(const char* centralityEstimator) {
    fCentralityEstimator = centralityEstimator;}
  void SetCentralityBins20() {fCentralityBins20 = kTRUE;}
private:
//functions
	
	Bool_t SelectEvent(AliAODVertex *vertex);
	Int_t SelectTrack(AliAODTrack *track);
	
  
 private://objects
 AliAODEvent *fAOD;    //! AOD object
AliAODVertex *fAODVertex;

	
TH1D        *fHistNchPt; //! 
TH1D *fHistNchEta;//!
 TH1D *fHistNchEtaCent;
TH1D *fHistNchPhi;
 TH1D *fHistDCAxy;
 TH1D *fHistDCAz;
 TH1D *fHistnclus;
 TH1D *fHistchi2ndf;  
 TH2D *fHistchi2ndfvscs;
 

TH1F *fMult[91];
TH1F *fMultTwo[46];
TH1F *fMultFive[15];

 TH1D *fHistVz;
 	TH1F	*fHistMultV0A;
	TH1F	*fHistMultV0C;
		TH1F	*fHistMultV0total;
		TNtuple *My_ntuple;
TList       *fOutputList; //! Output list
	

 TString fCentralityEstimator;//"V0M","TRK","TKL","ZDC","FMD"
 Bool_t fCentralityBins20;//centrality bins of 5% width

 TH1D *fCentralityCounter;
 TH1D *fEventCounter;
TH1D *histcounter;

  AliEbyEMultFluctuationTask(const AliEbyEMultFluctuationTask&); 
  AliEbyEMultFluctuationTask& operator=(const AliEbyEMultFluctuationTask&);

  ClassDef(AliEbyEMultFluctuationTask, 1); 
};

#endif
