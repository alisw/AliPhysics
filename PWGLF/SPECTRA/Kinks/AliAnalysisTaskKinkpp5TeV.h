#ifndef AliAnalysisTaskKinkpp5TeV_cxx 
#define AliAnalysisTaskKinkpp5TeV_cxx

//Author :: Nur Hussain, Gauhati University 
// Thanks to Martha Spyropoulou-Stassinaki for her suggestions for the modification


class TH1F;
class TH2F;
class TH3F;
class AliPIDResponse;
class AliESDtrackCuts;
//class fITSPIDResponse;
//class AliPIDResponse;

///class TFile;
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskKinkpp5TeV : public AliAnalysisTaskSE {
 public:
  	AliAnalysisTaskKinkpp5TeV();
	AliAnalysisTaskKinkpp5TeV(const char *name, Float_t lRadiusKUp,  Float_t lRadiusKLow, Int_t lNCluster, Float_t lLowQtValue, Float_t yRange);
  	virtual ~AliAnalysisTaskKinkpp5TeV() {}
  
  	virtual void   UserCreateOutputObjects();
  	virtual void   UserExec(Option_t *option);
  	virtual void   Terminate(Option_t *);
	
//	virtual Float_t GetVertex(AliESDEvent* fESD) const;
	const AliESDVertex *GetEventVertex(const AliESDEvent* esd) const;  
	void SetYRange(Float_t  yRange){fRapiK=yRange;}
        void SetKinkRadius(Float_t lRadiusKLow, Float_t lRadiusKUp)  { fKinkRadLow=lRadiusKLow; fKinkRadUp=lRadiusKUp;}
        void SetNCluster(Int_t lNCluster){fLowCluster=lNCluster;}
        void SetLowQtValue(Float_t   lLowQtValue){fLowQt=lLowQtValue;}

 private:
	void   Process();
  	TList       *fOutputList; //! Output list
  	TH1F        *fHistPt; //! Pt spectrum
	TH1F *fMultiplicity;
	TH1F *fIncompletEvent;
	TH1F *fMultpileup;
	TH1F *fMultV0trigger;
	TH1F *fZvertex;
	TH1F *fEventVertex;
	TH1F * fRatioCrossedRows;
	TH2F *fZvXv;
	TH2F *fZvYv;
	TH2F *fXvYv;
	TH1D *fRpr;
	TH1D *fdcaToVertexXY;
	TH1D *fdcaToVertexXYafterCut;
	TH1F *fptAllKink;
	TH1F *fRatioCrossedRowsKink;
	TH2F *fPosiKink;
	TH1F *fQtAll;
	TH1F *fptKink;
	TH2F *fQtMothP;
	TH1F *fqT1;
	TH1F *fEta;
	TH1F *fqT2;
	TH1F * fKinkKaonBackg;
	TF1 *f1;
   	TF1 *f2;
	TH1F * fPtCut1;
	TH2F *fAngMotherPi;
	TH2F *fQtInvM;
	TH1F *fInvMuNuAll;
	TH2F *fInvMassMuNuPtAll;
	TH2F * fRadiusPt;
	TH2F *fAngMotherKC;
	TH1F *fkaonInvaiant;	
	TH2F *fRadiusNcl;
	TH1F *fPtKPDG;
	TH2F *fAngMotherKKinks;
	TH1F *fPtCut2;
	TH1F *fPtCut3;
	TH2F *fTPCSignlMotherK;
	TH1F *fPtKaon;
	TH1F *fPtKaonP;
	TH1F *fPtKaonN;
	TH2F *fTPCSignalP;

	TH2F *fRadiusNclCln;
	TH2F *fRadiusPtcln;
	TH2F *fInvMassMuNuPt;
	TH2F *fTPCSignlPtpc;	
	TH2F *fMothKinkMomSignl;
	TH2F *fTPCSignlKinkDau;
	TH2F *fTPCMomNSigmaAllKaon;
	TH1F *fnSigmaTPC;
	TH1F *fradiurKink;
	TH1F *fLenthKink;
	TH1F *fEtaK;
	TH1F *frapiKESD;
	TH2F *fzVertexPositionKinKvsKinkRad;
	TH2F * fSignPtNcl;
	TH2F *fSignPtrapiK;
	TH2F *frapiKNcl;
	TH1F *fSignPt;
	TH2F *fChi2NclTPC;
	TH1F *fRatioChi2Ncl;	
	TH1F *flifetime;
	TH1F *fPtKinkKaon;
	TH1F *fDCAkink;
	TH1F *fPtKink;
	TH1F *fPtKinkPos;
	TH1F *fPtKinkNeg;
	TH1F *fPtKinkK0;
	TH1F *fPtKinkK0P;
	TH1F *fPtKinkK0N;
	TH1F *fPtKinkGyu;
	TH1F *fPtKinkGyuP;
	TH1F *fPtKinkGyuN;
	TH1F *fKinKRbn;
	TH3F *fradPtRpDt;
	TH2F *fAngMomK;
	TH2F *fPosiKinkK;
	TH2F *fPosiKinKXZ;	
	TH2F *fPosiKinKYZ;	
	TH1F *fNumberOfEvent;	
	TH1F *fNumberOfEvent_cent;
	TH1F  *fbgCleaningHigh;
	TH2F *fTPCSignalPt;
        TH2F *fqTvsPt;
        TH2F *fInvMassPt;
	TH1F *fCent;
	TH2F *fEventVsCentrality;

	Float_t fKinkRadUp;
	Float_t fKinkRadLow;
	Int_t fLowCluster;
	Float_t fLowQt;
	Float_t fRapiK;



	Float_t fVtxCut;	
	AliPIDResponse *fPIDResponse;
	AliESDtrackCuts * fESDtrackCuts;

 	AliAnalysisTaskKinkpp5TeV(const AliAnalysisTaskKinkpp5TeV&); // not implemented
  	AliAnalysisTaskKinkpp5TeV& operator=(const AliAnalysisTaskKinkpp5TeV&); // not implemented
  
  	ClassDef(AliAnalysisTaskKinkpp5TeV, 1); // example of analysis
};

#endif
