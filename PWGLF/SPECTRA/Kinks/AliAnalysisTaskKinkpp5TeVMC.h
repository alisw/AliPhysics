#ifndef AliAnalysisTaskKinkpp5TeVMC_cxx 
#define AliAnalysisTaskKinkpp5TeVMC_cxx

// example of an analysis task creating a p_t spectrum

class TH1F;
class TH2F;
class TH3F;
class AliPIDResponse;
class AliESDtrackCuts;
//class fITSPIDResponse;
//class AliPIDResponse;

///class TFile;
#include "AliAnalysisTaskSE.h"
#include "AliITSPIDResponse.h"
#include "THnSparse.h"

class AliAnalysisTaskKinkpp5TeVMC : public AliAnalysisTaskSE {
 public:
 AliAnalysisTaskKinkpp5TeVMC() : AliAnalysisTaskSE(), fOutputList(0), fHistPt(0),fVtxCut(10.),fMultiplicity(0),fIncompletEvent(0),fMultpileup(0),fMultV0trigger(0),fZvertex(0),
	fEventVertex(0),fRatioCrossedRows(0),fZvXv(0), fZvYv(0), fXvYv(0),fRpr(0),fdcaToVertexXY(0),fdcaToVertexXYafterCut(0),fptAllKink(0),fRatioCrossedRowsKink(0),fPosiKink(0),
	fQtAll(0),fptKink(0),fQtMothP(0),fqT1(0),fEta(0), fqT2(0),fKinkKaonBackg(0),f1(0),f2(0),fPtCut1(0),fAngMotherPi(0),
	fQtInvM(0),fInvMuNuAll(0),fInvMassMuNuPtAll(0),fRadiusPt(0), fKinkRadUp(200.), fKinkRadLow(130.), fLowCluster(30), fLowQt(.12), fRapiK(0.5),
	fAngMotherKC(0),fkaonInvaiant(0),fRadiusNcl(0),fPtKPDG(0),fAngMotherKKinks(0),fPtCut2(0),fPtCut3(0),fTPCSignlMotherK(0),fPtKaon(0), fPtKaonP(0), fPtKaonN(0),
	fTPCSignalP(0),fRadiusNclCln(0),fRadiusPtcln(0),fInvMassMuNuPt(0),fTPCSignlPtpc(0), fMothKinkMomSignl(0),fTPCSignlKinkDau(0),fTPCMomNSigmaAllKaon(0),
	fnSigmaTPC(0),fradiurKink(0),fLenthKink(0),fEtaK(0),frapiKESD(0),fzVertexPositionKinKvsKinkRad(0),fSignPtNcl(0),fSignPtrapiK(0),frapiKNcl(0),fSignPt(0),
	fChi2NclTPC(0),fRatioChi2Ncl(0),flifetime(0),fPtKinkKaon(0),fDCAkink(0),fPtKink(0),fPtKinkPos(0),fPtKinkNeg(0),fPtKinkK0(0),fPtKinkK0P(0),fPtKinkK0N(0),
	fPtKinkGyu(0),fPtKinkGyuP(0),fPtKinkGyuN(0),fKinKRbn(0),fradPtRpDt(0),fAngMomK(0),fPosiKinkK(0),fPosiKinKXZ(0), fPosiKinKYZ(0),fPIDResponse(0),fNumberOfEvent(0),
	fMultMC_wo_any_cut(0),fMultMC_incompleteDAQ(0),fMultMC_AfterPileUp(0),fMultTriggerMCAfterV0(0),fTrigSel(AliVEvent::kINT7),
	frapidKMC(0), fptKMC(0),fSignPtGen(0),fPtKPlMC(0),fPtKMnMC(0),flengthTrackRef(0),flifetimeMC(0),flifeSmallMC(0),flifeInt(0),flifeYuri(0),flifetiMCK(0),
	flenYuri(0),flengthMCK(0),flifetime_MCprocess4(0),fradPtRapMC(0),flifetime_kaonpionPDG(0),flifetime_kaonmuonPDG(0),fmaxAngMomKmuMC(0),fradPtRapDC(0),
	fradMC(0),fQtKMuMC(0),fgenPtEtR(0),fgenPtEtRP(0),fMCEtaKaon(0),fSignPtEtaMC(0),fSignPtMC(0),fgenPtEtRN(0),fQtKElMC(0),fQtKPiMC(0),fESDMult(0),
	fFakepipi(0),fFakeKPi(0),fRadiusPtPion(0),fRadiusPtKaon(0),fQtKMu(0),fQtKPi(0),fQtKEl(0),fQtK3PiP(0),fQtK3PiM(0),fHistPtKPDG(0), fHiPtKPDGP(0), fHiPtKPDGN(0),
        fHistEta(0), frapidESDK(0), fHistQt2(0),fAngMomPi(0),fMinvPi(0),fMinvKa(0),fcodeH(0), fZkinkZDau(0),fRadiusPtFake(0),fTPCMomNSgnl(0),fPtPrKink(0),flifTiESDK(0),
	fKinkKaon(0),fkinkKaonP(0),fkinkKaonN(0),fcode2(0),fTPCSgnlPtpc(0),fMothKinkMomSgnlD(0),fKinkKaonBg(0),fMothKinkMomSgnl(0),fcodeDau2(0),fTPCSgnlKinkDau(0),
	fMinvPr(0),fDCAkinkBG(0),fPosiKinKBgXY(0),fPosiKinKBgZY(0),fPosiKinKBgZX(0),fKinKBGP(0),fKinKBGN(0),fdcodeH(0),fcode4(0),
	fNumberOfEvent_cent(0),fESDtrackCuts(0),fEventVsCentrality(0)



  {}
  	AliAnalysisTaskKinkpp5TeVMC(const char *name);
  	virtual ~AliAnalysisTaskKinkpp5TeVMC() {}
  
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
	TH1F *fMultMC_wo_any_cut;
	TH1F *fMultMC_incompleteDAQ;
	TH1F *fMultMC_AfterPileUp;
	TH1F *fMultTriggerMCAfterV0;
	TH1F *frapidKMC;
	TH1F *fptKMC;
	TH1F *fSignPtGen;
	TH1F *fPtKPlMC;
	TH1F *fPtKMnMC;
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
	TH1F *fgenPtEtR;
	TH1F *fgenPtEtRP;
	TH1F *fMCEtaKaon;
	TH2F *fSignPtEtaMC;
	TH1F *fSignPtMC;
	TH1F *fgenPtEtRN;
	TH1F *fQtKElMC;
	TH1F *fQtKPiMC;
	TH1F * fESDMult;
	TH1F *fFakepipi;
	TH1F *fFakeKPi;
	TH2F *fRadiusPtPion;
	TH2F *fRadiusPtKaon;

	TH1F *fQtKMu;
	TH1F *fQtKPi;
	TH1F *fQtKEl;
	TH1F *fQtK3PiP;
	TH1F *fQtK3PiM;
	TH1F *fHistPtKPDG;
 	TH1F *fHiPtKPDGP;
 	TH1F *fHiPtKPDGN;
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
	TH1F *fPtPrKink;
	TH1F *flifTiESDK;
	TH1F *fKinkKaon;
	TH1F *fkinkKaonP;
	TH1F *fkinkKaonN;
	TH2F *fcode2;
	TH2F *fTPCSgnlPtpc;
	TH2F *fMothKinkMomSgnlD;
	TH1F * fKinkKaonBg;
	TH2F *fMothKinkMomSgnl;
	TH2F *fcodeDau2;
	TH2F *fTPCSgnlKinkDau;
	TH1F *fMinvPr;
	TH1F *fDCAkinkBG;
	TH2F *fPosiKinKBgXY;
	TH2F *fPosiKinKBgZY;
	TH2F *fPosiKinKBgZX;
	TH1F *fKinKBGP;
	TH1F *fKinKBGN;
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

 	AliAnalysisTaskKinkpp5TeVMC(const AliAnalysisTaskKinkpp5TeVMC&); // not implemented
  	AliAnalysisTaskKinkpp5TeVMC& operator=(const AliAnalysisTaskKinkpp5TeVMC&); // not implemented
  
  	ClassDef(AliAnalysisTaskKinkpp5TeVMC, 1); // example of analysis
};

#endif
