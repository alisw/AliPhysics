#ifndef AliAnalysisPionKinksMCESD_h
#define AliAnalysisPionKinksMCESD_h


/*  See cxx source for full Copyright notice */

//-----------------------------------------------------------------
//                 AliAnalysisPionKinksMCESD class
//         This task is an example of an analysis task
//                  for kink topology Study
//          Authors: Eftychios Cheiladakis, under the supervision
//           of Martha's Spyropoulou-Stassinaki   
//          Physics Department of Athens University
//------------------------------------------------------------------
class AliESDEvent;
class TF1;
class TH1;
class TH2;
class TH3;
class TParticle;
class TParticle;
class AliESDtrackCuts;
class AliPhysicsSelection;
class AliPIDResponse;

#include "AliAnalysisTaskSE.h"

class AliAnalysisPionKinksMCESD : public AliAnalysisTaskSE {
	public:
		AliAnalysisPionKinksMCESD(const char *name = "AliAnalysisPionKinksMCESD");
		virtual ~AliAnalysisPionKinksMCESD() {}
  
		virtual void   UserCreateOutputObjects();
		virtual void   UserExec(Option_t *option);
		virtual void   Terminate(Option_t *);

  		const AliESDVertex *GetEventVertex(AliESDEvent* esd);
		Double_t Energy(AliESDtrack* track) const;
		Double_t GetMCRapidity(TParticle* particle) const;
		Double_t fuRapidity(AliESDtrack* track) const;
		Double_t MCPQt(AliMCEvent* mcEvent, Int_t iMC, TParticle* MCdaughter) const;
		Double_t fuMCKinkAngle(AliMCEvent* mcEvent, Int_t iMC, TParticle* MCdaughter, Bool_t degrees) const;
		Bool_t IsGoodTrack(AliESDtrack* ESDTrack) const;
		Bool_t IsPrimaryTrack(AliESDtrack* ESDTrack) const;
		Bool_t IsRealPionKink(AliESDkink* kink, AliStack* MCstack) const;

		// Set limits & cuts
		void SetMulCut(Int_t low, Int_t up){fLowMulcut=low; fUpMulcut=up;}	  
		void SetPtCut(Double_t PtCut){cLowPt=PtCut;}
		void SetRapidityLimits(Double_t RapidityLim){cRapidityLim=RapidityLim;}
		void SetRadiusRange(Double_t LowR, Double_t UpR){cLowR=LowR; cUpR=UpR;}
		void SetZRange(Double_t LowZ, Double_t UpZ){cLowZ=LowZ; cUpZ=UpZ;}		
		void SetKinkAngleCut(Double_t LowAngle){cLowKinkAngle=LowAngle;}
		void SetQtRange(Double_t LowQt, Double_t UpQt){cLowQt=LowQt; cUpQt=UpQt;}
		void SetInvMassRange(Double_t LowInvMass, Double_t UpInvMass){cLowInvMass=LowInvMass; cUpInvMass=UpInvMass;}
		void SetSigmaCut(Double_t SigmaCut){cSigmaCut=SigmaCut;}
		void SetPdgCodes(Int_t kaon, Int_t pion, Int_t muon, Int_t electron){cPdgKaon=kaon; cPdgPion=pion; cPdgMuon=muon; cPdgElectron=electron;}
		void SetMasses(Double_t KaonMass, Double_t PionMass, Double_t MuonMass, Double_t ElectronMass){cKaonMass=KaonMass; cPionMass=PionMass; cMuonMass=MuonMass; cElectronMass=ElectronMass;}

		// Set histos limits
		void SetMultHistos(Int_t BinsMult, Int_t LowMult, Int_t UpMult){nBinsMult=BinsMult; hLowPt=LowMult; hUpMult=UpMult;}
		void SetPtHistos(Int_t BinsPt, Double_t LowPt, Double_t UpPt){nBinsPt=BinsPt; hLowPt=LowPt; hUpPt=UpPt;}   
		void SetEtaHistos(Int_t BinsEta, Double_t LowEta, Double_t UpEta){nBinsEta=BinsEta; hLowEta=LowEta; hUpEta=UpEta;}
		void SetQtHistos(Int_t BinsQt, Double_t LowQt, Double_t UpQt){nBinsQt=BinsQt; hLowQt=LowQt; hUpQt=UpQt;}
		void SetPdgHistos(Int_t BinsPdg, Double_t LowPdg, Double_t UpPdg){nBinsPdg=BinsPdg; hLowPdg=LowPdg; hUpPdg=UpPdg;}
		void SetPdg2Histos(Int_t BinsPdg2, Double_t LowPdg2, Double_t UpPdg2){nBinsPdg2=BinsPdg2; hLowPdg2=LowPdg2; hUpPdg2=UpPdg2;}
		void SetUIDHistos(Int_t BinsUID, Double_t LowUID, Double_t UpUID){nBinsUID=BinsUID; hLowUID=LowUID; hUpUID=UpUID;}
		void SetRHistos(Int_t BinsR, Double_t LowR, Double_t UpR){nBinsR=BinsR; hLowR=LowR; hUpR=UpR;}
		void SetZHistos(Int_t BinsZ, Double_t LowZ, Double_t UpZ){nBinsZ=BinsZ; hLowZ=LowZ; hUpZ=UpZ;}
		void SetXYHistos(Int_t BinsXY, Double_t LowXY, Double_t UpXY){nBinsXY=BinsXY; hLowXY=LowXY; hUpXY=UpXY;}
		void SetAngleHistos(Int_t BinsAngle, Double_t LowAngle, Double_t UpAngle){nBinsAngle=BinsAngle; hLowAngle=LowAngle; hUpAngle=UpAngle;}
		void SetZVHistos(Int_t BinsZV, Double_t LowZV, Double_t UpZV){nBinsZV=BinsZV; hLowZV=LowZV; hUpZV=UpZV;}
		void SetXYVHistos(Int_t BinsXYV, Double_t LowXYV, Double_t UpXYV){nBinsXYV=BinsXYV; hLowXYV=LowXYV; hUpXYV=UpXYV;}
		void SetInvMassHistos(Int_t BinsInvMass, Double_t LowInvMass, Double_t UpInvMass){nBinsInvMass=BinsInvMass; hLowInvMass=LowInvMass; hUpInvMass=UpInvMass;}
		void SetdEdxHistos(Int_t BinsdEdx, Double_t LowdEdx, Double_t UpdEdx){nBinsdEdx=BinsdEdx; hLowdEdx=LowdEdx; hUpdEdx=UpdEdx;}
	private:
		TF1* fMaxKinkAngKmu; 
		TF1* fMaxKinkAngPimu;

		TH1F* hMCMult;
		TH1F* hMCMultPrim;
		TH1F* hMCPtAll;
		TH1F* hMCEtaAll;
		TH1F* hMCPtPrim; 
		TH1F* hMCEtaPrim; 
		TH1F* hMCPt; 
		TH1F* hMCEta;  
		TH1F* hMCPdg; 
		TH1F* hMCMultPiPlus; 
		TH1F* hMCPtPiPlus; 
		TH1F* hMCEtaPiPlus; 
		TH1F* hMCRapidityPiPlus;
		TH1F* hMCNDaughtersPlus;
		TH1F* hMCRadPiDauPlus; 
		TH1F* hMCKinkPosZPlus;
		TH1F* hMCUIDPiDauPlus; 
		TH1F* hMCPdgPiNonDecayedPlus; 
		TH1F* hMCPdgPiDauPlus;
		TH1F* hMCQtPlus; 
		TH1F* hMCKinkAnglePlus;
		TH2F* hMCPKinkAngPlus;
		TH2F* hMCPdgCodemdPlus; 
		TH2F* hMCPtmdPlus; 
		TH1F* hMCPtPimuonPlus;
		TH1F* hMCEtaPimuonPlus;
		TH1F* hMCRapidityPimuonPlus;
		TH1F* hMCQtPimuonPlus; 
		TH1F* hMCPKinkAngPimuonPlus;
		TH1F* hMCPtPiotherPlus;
		TH1F* hMCEtaPiotherPlus;
		TH1F* hMCRapidityPiotherPlus;
		TH1F* hMCQtPiotherPlus; 
		TH1F* hMCPKinkAngPiotherPlus;
		TH1F* hMCMultPiMinus; 
		TH1F* hMCPtPiMinus; 
		TH1F* hMCEtaPiMinus; 
		TH1F* hMCRapidityPiMinus;
		TH1F* hMCNDaughtersMinus;
		TH1F* hMCRadPiDauMinus; 
		TH1F* hMCKinkPosZMinus;
		TH1F* hMCUIDPiDauMinus; 
		TH1F* hMCPdgPiNonDecayedMinus; 
		TH1F* hMCPdgPiDauMinus;
		TH1F* hMCQtMinus; 
		TH1F* hMCKinkAngleMinus;
		TH2F* hMCPKinkAngMinus;
		TH2F* hMCPdgCodemdMinus; 
		TH2F* hMCPtmdMinus; 
		TH1F* hMCPtPimuonMinus;
		TH1F* hMCEtaPimuonMinus;
		TH1F* hMCRapidityPimuonMinus;
		TH1F* hMCQtPimuonMinus; 
		TH1F* hMCPKinkAngPimuonMinus;
		TH1F* hMCPtPiotherMinus;
		TH1F* hMCEtaPiotherMinus;
		TH1F* hMCRapidityPiotherMinus;
		TH1F* hMCQtPiotherMinus; 
		TH1F* hMCPKinkAngPiotherMinus;

		TH1F* hMult;
		TH1F* hAcceptedMult;
		TH1F* hMultPS;
		TH3F* hvtx;
		TH2F* hvtxy;
		TH2F* hvtyz;
		TH2F* hvtxz;
		TH1F* hMultPSV;
		TH1F* hPtAll;
		TH1F* hEtaAll;
		TH3F* hTrackPos;
		TH2F* hTrackPosxy;
		TH2F* hTrackPosyz;
		TH2F* hTrackPosxz;
		//TH1F* hTPCchi2clusters;
		//TH1F* hdcaToVertexXY;
		//TH1F* hdcaToVertexZ;
		TH1F* hMultPrim;
		TH1F* hPtPrim;
		TH1F* hEtaPrim;
		TH3F* hPrimTrackPos;
		TH2F* hPrimTrackPosxy;
		TH2F* hPrimTrackPosyz;
		TH2F* hPrimTrackPosxz;
		TH1F* hPt;
		TH1F* hEta;
		//TH1F* hRapidity;
		TH1F* hPtKink;
		TH1F* hEtaKink;
		TH1F* hRapidityKink;
		TH2F* hPmP;
		TH1F* hKinkPosRTPCclusters1;
		TH2F* hKinkPosRTPCclusters2;
		TH1F* hQt;
		TH1F* hKinkAngle;
		TH1F* hDCAkink;
		TH2F* hPmKinkAng;
		TH2F* hKinkPosXY;
		TH2F* hKinkPosZY;
		TH2F* hKinkPosZR;
		TH1F* hKinkPosR;
		TH1F* hKinkPosZ;
		TH2F* hKinkPosZMCKinkPosZ;
		TH2F* hPdgCodemd;
		TH2F* hPmd;
		TH1F* hMinvPimu;
		TH1F* hUIDKinkDau;
		TH2F* hdEdx;
		TH1F* hPtKinkFake;
		TH1F* hEtaKinkFake;
		TH1F* hRapidityKinkFake;
		TH2F* hPmPFake;
		TH1F* hKinkPosRTPCclusters1Fake;
		TH2F* hKinkPosRTPCclusters2Fake;
		TH1F* hQtFake;
		TH1F* hKinkAngleFake;
		TH1F* hDCAkinkFake;
		TH2F* hPmKinkAngFake;
		TH2F* hKinkPosXYFake;
		TH2F* hKinkPosZYFake;
		TH2F* hKinkPosZRFake;
		TH1F* hKinkPosRFake;
		TH1F* hKinkPosZFake;
		TH2F* hKinkPosZMCKinkPosZFake;
		TH2F* hPdgCodemdFake;
		TH2F* hPmdFake;
		TH1F* hMinvPimuFake;
		TH1F* hUIDKinkDauFake;
		TH2F* hdEdxFake;
		TH1F* hPtPosRSelected;
		TH2F* hPdgCodemdZRejected;
		TH1F* hPtZSelected;
		TH2F* hPdgCodemdAngRejected;
		TH1F* hPtAngSelected;
		TH2F* hPdgCodemdPmRejected;
		TH1F* hPtPmSelected;
		TH2F* hPdgCodemdQtLowRejected;
		TH1F* hPtGoodKink; 
		TH1F* hEtaGoodKink; 
		TH1F* hRapidityGoodKink; 
		TH1F* hQtGoodKink; 
		TH2F* hPmGoodKinkAng;  
		TH2F* hPdgCodemdGoodKink; 
		TH2F* hPmdGoodKink;
		TH1F* hUIDGoodKinkDau;
		TH2F* hdEdxGoodKink;
		TH1F* hUIDPiDauPlus;
		TH1F* hMultPiPlus;
		TH1F* hPtPiPlus;
		TH1F* hEtaPiPlus;
		TH1F* hRapidityPiPlus;
		TH1F* hQtPiPlus;
		TH1F* hKinkAnglePiPlus;
		TH2F* hPmKinkAngPiPlus;
		TH2F* hKinkPosXYPiPlus;
		TH2F* hKinkPosZRPiPlus;
		TH1F* hKinkPosRPiPlus;
		TH1F* hDCAkinkPiPlus;
		TH2F* hPdgCodemdPiPlus;
		TH2F* hPmdPiPlus;
		TH2F* hdEdxPiPlus;
		TH1F* hQtPimuPlus;
		TH2F* hPmKinkAngPimuPlus;
		TH1F* hQtPiotherPlus;
		TH2F* hPmKinkAngPiotherPlus; 
		TH2F* hPdgCodemdPiotherPlus;
		TH1F* hUIDPiDauMinus;
		TH1F* hMultPiMinus;
		TH1F* hPtPiMinus;
		TH1F* hEtaPiMinus;
		TH1F* hRapidityPiMinus;
		TH1F* hQtPiMinus;
		TH1F* hKinkAnglePiMinus;
		TH2F* hPmKinkAngPiMinus;
		TH2F* hKinkPosXYPiMinus;
		TH2F* hKinkPosZRPiMinus;
		TH1F* hKinkPosRPiMinus;
		TH1F* hDCAkinkPiMinus;
		TH2F* hPdgCodemdPiMinus;
		TH2F* hPmdPiMinus;
		TH2F* hdEdxPiMinus;
		TH1F* hQtPimuMinus;
		TH2F* hPmKinkAngPimuMinus;
		TH1F* hQtPiotherMinus;
		TH2F* hPmKinkAngPiotherMinus; 
		TH2F* hPdgCodemdPiotherMinus;
		TH2F* hPdgCodemdQtRejected;
		TH1F* hPtQtSelected;
		TH2F* hPdgCodemdMaxAngRejected;
		TH1F* hPtMaxAngSelected;
		TH2F* hPdgCodemdRTPCclustersRejected;
		TH1F* hPtRTPCclustersSelected;
		TH2F* hRTPCclustersRTPCclustersSelected;
		TH2F* hPdgCodemdMinvRejected;
		TH1F* hPtSelected; 
		TH1F* hEtaSelected; 
		TH1F* hRapiditySelected; 
		TH1F* hQtSelected; 
		TH1F* hKinkAngleSelected;
		TH1F* hDCAkinkSelected;
		TH2F* hPmKinkAngSelected; 
		TH2F* hKinkPosXYSelected; 
		TH2F* hKinkPosZRSelected; 
		TH1F* hKinkPosRSelected; 
		TH2F* hPdgCodemdSelected; 
		TH2F* hPmdSelected;
		TH1F* hMinvPimuSelected;  
		TH1F* hUIDKinkDauSelected; 
		TH2F* hdEdxSelected; 
		TH1F* hPtSelectedFake; 
		TH1F* hEtaSelectedFake; 
		TH1F* hRapiditySelectedFake; 
		TH1F* hQtSelectedFake; 
		TH1F* hKinkAngleSelectedFake;
		TH1F* hDCAkinkSelectedFake;
		TH2F* hPmKinkAngSelectedFake;
		TH2F* hKinkPosXYSelectedFake;
		TH2F* hKinkPosZRSelectedFake;
		TH1F* hKinkPosRSelectedFake;
		TH2F* hPmdSelectedFake;
		TH1F* hMinvPimuSelectedFake; 
		TH2F* hdEdxSelectedFake;
		TH2F* hPdgCodemddEdxRejected;
		TH1F* hPtPiSelected; 
		TH1F* hEtaPiSelected; 
		TH1F* hRapidityPiSelected; 
		TH1F* hQtPiSelected; 
		TH1F* hKinkAnglePiSelected;
		TH1F* hDCAkinkPiSelected;
		TH2F* hPmKinkAngPiSelected; 
		TH1F* hKinkPosRTPCclusters1PiSelected;
		TH2F* hKinkPosRTPCclusters2PiSelected;
		TH2F* hKinkPosXYPiSelected; 
		TH2F* hKinkPosZRPiSelected; 
		TH1F* hKinkPosRPiSelected; 
		TH1F* hKinkPosZPiSelected; 
		TH2F* hPmPPiSelected; 
		TH2F* hPdgCodemdPiSelected; 
		TH2F* hPmdPiSelected;
		TH1F* hMinvPimuPiSelected;  
		TH1F* hUIDKinkDauPiSelected; 
		TH2F* hdEdxPiSelected;

		TH1F* hPtPiSelectedPlus; //plus
		TH1F* hEtaPiSelectedPlus; 
		TH1F* hRapidityPiSelectedPlus; 
		TH1F* hQtPiSelectedPlus; 
		TH1F* hKinkAnglePiSelectedPlus;
		TH1F* hDCAkinkPiSelectedPlus;
		TH2F* hPmKinkAngPiSelectedPlus; 
		TH2F* hKinkPosXYPiSelectedPlus; 
		TH2F* hKinkPosZRPiSelectedPlus;  
		TH2F* hPdgCodemdPiSelectedPlus; 
		TH2F* hPmdPiSelectedPlus;
		TH1F* hMinvPimuPiSelectedPlus;  
		TH1F* hUIDPiDaumuSelectedPlus; 
		TH2F* hdEdxPiSelectedPlus;
		TH1F* hPtPrimPiKinksPlus; 
		TH1F* hEtaPrimPiKinksPlus; 
		TH1F* hRapidityPrimPiKinksPlus;
		TH1F* hPtSecondPiKinksPlus; 
		TH1F* hEtaSecondPiKinksPlus; 
		TH1F* hRapiditySecondPiKinksPlus;
		TH1F* hPtNonPiKinksPlus; 
		TH1F* hEtaNonPiKinksPlus; 
		TH1F* hRapidityNonPiKinksPlus; 
		TH2F* hPdgCodemdNonPiKinksPlus; 
		TH1F* hPtPiSelectedMinus;//minus
		TH1F* hEtaPiSelectedMinus; 
		TH1F* hRapidityPiSelectedMinus; 
		TH1F* hQtPiSelectedMinus; 
		TH1F* hKinkAnglePiSelectedMinus;
		TH1F* hDCAkinkPiSelectedMinus;
		TH2F* hPmKinkAngPiSelectedMinus; 
		TH2F* hKinkPosXYPiSelectedMinus; 
		TH2F* hKinkPosZRPiSelectedMinus;  
		TH2F* hPdgCodemdPiSelectedMinus; 
		TH2F* hPmdPiSelectedMinus;
		TH1F* hMinvPimuPiSelectedMinus;  
		TH1F* hUIDPiDaumuSelectedMinus; 
		TH2F* hdEdxPiSelectedMinus;
		TH1F* hPtPrimPiKinksMinus; 
		TH1F* hEtaPrimPiKinksMinus; 
		TH1F* hRapidityPrimPiKinksMinus;
		TH1F* hPtSecondPiKinksMinus; 
		TH1F* hEtaSecondPiKinksMinus; 
		TH1F* hRapiditySecondPiKinksMinus;
		TH1F* hPtNonPiKinksMinus; 
		TH1F* hEtaNonPiKinksMinus; 
		TH1F* hRapidityNonPiKinksMinus; 
		TH2F* hPdgCodemdNonPiKinksMinus;// reconstruction histograms
		
		TList* fListOfHistos;

		// Limits and cuts
		Int_t fLowMulcut;
		Int_t fUpMulcut;
		
		Double_t cLowPt;
		Double_t cRapidityLim;
		Double_t cLowR, cUpR;
		Double_t cLowZ, cUpZ;
		Double_t cLowKinkAngle;
		Double_t cLowQt, cUpQt;
		Double_t cLowInvMass, cUpInvMass;
		Double_t cSigmaCut;
		Int_t cPdgKaon, cPdgPion, cPdgMuon, cPdgElectron;
		Double_t cKaonMass, cPionMass, cMuonMass, cElectronMass;

		// Histos limits
		Int_t nBinsMult, hLowMult, hUpMult;
		Int_t nBinsPt;
		Double_t hLowPt, hUpPt;
		Int_t nBinsEta;
		Double_t hLowEta, hUpEta;
		Int_t nBinsQt;
		Double_t hLowQt, hUpQt;
		Int_t nBinsPdg;
		Double_t hLowPdg, hUpPdg;
		Int_t nBinsPdg2;
		Double_t hLowPdg2, hUpPdg2;
		Int_t nBinsUID;
		Double_t hLowUID, hUpUID;
		Int_t nBinsR;
		Double_t hLowR, hUpR;
		Int_t nBinsZ;
		Double_t hLowZ, hUpZ;
		Int_t nBinsXY;
		Double_t hLowXY, hUpXY;
		Int_t nBinsAngle;
		Double_t hLowAngle, hUpAngle;
		Int_t nBinsZV;
		Double_t hLowZV, hUpZV;
		Int_t nBinsXYV;
		Double_t hLowXYV, hUpXYV;
		Int_t nBinsInvMass;
		Double_t hLowInvMass, hUpInvMass;
		Int_t nBinsdEdx;
		Double_t hLowdEdx, hUpdEdx;

		AliPIDResponse *fPIDResponse;     //! PID response object
		AliESDtrackCuts* fMaxDCAtoVtxCut;
		AliESDtrackCuts* fTrackCuts;

		
		AliAnalysisPionKinksMCESD(const AliAnalysisPionKinksMCESD&); 
		AliAnalysisPionKinksMCESD& operator=(const AliAnalysisPionKinksMCESD&);
  
	ClassDef(AliAnalysisPionKinksMCESD, 1); 
};
#endif	


