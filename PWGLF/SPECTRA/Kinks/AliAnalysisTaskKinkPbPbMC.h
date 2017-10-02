#ifndef AliAnalysisTaskKinkPbPbMC_cxx 
#define AliAnalysisTaskKinkPbPbMC_cxx

class TH1F;
class TH2F;
class TH3F;
class AliPIDResponse;
class AliESDtrackCuts;

///class TFile;
#include "AliAnalysisTaskSE.h"
#include "AliITSPIDResponse.h"

class AliAnalysisTaskKinkPbPbMC : public AliAnalysisTaskSE {
 public:
  	AliAnalysisTaskKinkPbPbMC();
	AliAnalysisTaskKinkPbPbMC(const char *name, Float_t lRadiusKUp,  Float_t lRadiusKLow, Int_t lNCluster, Float_t lLowQtValue, Float_t yRange);
  	virtual ~AliAnalysisTaskKinkPbPbMC() {}
  
  	virtual void   UserCreateOutputObjects();
  	virtual void   UserExec(Option_t *option);
  	virtual void   Terminate(Option_t *);
	
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
	TH2F *fPtKPDG;
	TH2F *fAngMotherKKinks;
	TH2F *fPtCut2;
	TH2F *fPtCut3;
	TH2F *fTPCSignlMotherK;
	TH2F *fPtKaon;
	TH2F *fPtKaonP;
	TH2F *fPtKaonN;
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
	TH2F *fPtKinkKaon;
	TH1F *fDCAkink;
	TH2F *fPtKink;
	TH2F *fPtKinkPos;
	TH2F *fPtKinkNeg;
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
	TH1F *fMultMC_wo_any_cut;
	TH1F *fMultMC_incompleteDAQ;
	TH1F *fMultMC_AfterPileUp;
	TH1F *fMultTriggerMCAfterV0;
	TH1F *frapidKMC;
	TH2F *fptKMC;
	TH1F *fSignPtGen;
	TH2F *fPtKPlMC;
	TH2F *fPtKMnMC;
	TH1F *flengthTrackRef;
	TH1F *flifetimeMC;
	TH1F *flifeSmallMC;
	TH1F *flifeInt;
	TH1F *flifeYuri;
	TH1F *flifetiMCK;
	TH1F *flenYuri;
	TH1F *flengthMCK;
	TH1F *flifetime_MCprocess4;
	TH3F *fradPtRapMC;
	TH1F *flifetime_kaonpionPDG;
	TH1F *flifetime_kaonmuonPDG;
	TH2F *fmaxAngMomKmuMC;
	TH3F *fradPtRapDC;
	TH1F *fradMC;
	TH1F *fQtKMuMC;
	TH2F *fgenPtEtR;
	TH2F *fgenPtEtRP;
	TH1F *fMCEtaKaon;
	TH2F *fSignPtEtaMC;
	TH1F *fSignPtMC;
	TH2F *fgenPtEtRN;
	TH1F *fQtKElMC;
	TH1F *fQtKPiMC;
	TH1F * fESDMult;
	TH2F *fFakepipi;
	TH2F *fFakeKPi;
	TH2F *fRadiusPtPion;
	TH2F *fRadiusPtKaon;

	TH1F *fQtKMu;
	TH1F *fQtKPi;
	TH1F *fQtKEl;
	TH1F *fQtK3PiP;
	TH1F *fQtK3PiM;
	TH2F *fHistPtKPDG;
 	TH2F *fHiPtKPDGP;
 	TH2F *fHiPtKPDGN;
 	TH1F *fHistEta;
 	TH1F *frapidESDK;
 	TH1F *fHistQt2;
	TH2F *fAngMomPi;	
	TH1F *fMinvPi;
	TH1F *fMinvKa;
	TH2F *fcodeH;
	TH2F *fZkinkZDau;
	TH2F *fRadiusPtFake;
	TH2F *fTPCMomNSgnl;
	TH2F *fPtPrKink;
	TH1F *flifTiESDK;
	TH2F *fKinkKaon;
	TH2F *fkinkKaonP;
	TH2F *fkinkKaonN;
	TH2F *fcode2;
	TH2F *fTPCSgnlPtpc;
	TH2F *fMothKinkMomSgnlD;
	TH2F * fKinkKaonBg;
	TH2F *fMothKinkMomSgnl;
	TH2F *fcodeDau2;
	TH2F *fTPCSgnlKinkDau;
	TH1F *fMinvPr;
	TH1F *fDCAkinkBG;
	TH2F *fPosiKinKBgXY;
	TH2F *fPosiKinKBgZY;
	TH2F *fPosiKinKBgZX;
	TH2F *fKinKBGP;
	TH2F *fKinKBGN;
	TH2F *fdcodeH;
	TH2F *fcode4;
	TH1F *fNumberOfEvent_cent;
	TH2F *fEventVsCentrality;

	Float_t fKinkRadUp;
	Float_t fKinkRadLow;
	Int_t fLowCluster;
	Float_t fLowQt;
	Float_t fRapiK;


	UInt_t fTrigSel;
	Float_t fVtxCut;	
	AliPIDResponse *fPIDResponse;
	AliESDtrackCuts * fESDtrackCuts;

 	AliAnalysisTaskKinkPbPbMC(const AliAnalysisTaskKinkPbPbMC&); // not implemented
  	AliAnalysisTaskKinkPbPbMC& operator=(const AliAnalysisTaskKinkPbPbMC&); // not implemented
  
  	ClassDef(AliAnalysisTaskKinkPbPbMC, 1); // example of analysis
};

#endif
