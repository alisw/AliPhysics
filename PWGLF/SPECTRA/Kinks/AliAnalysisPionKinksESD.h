#ifndef AliAnalysisPionKinksESD_h
#define AliAnalysisPionKinksESD_h


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

class AliAnalysisPionKinksESD : public AliAnalysisTaskSE {
	public:
		AliAnalysisPionKinksESD(const char *name = "AliAnalysisPionKinksESD");
		virtual ~AliAnalysisPionKinksESD() {}
  
		virtual void   UserCreateOutputObjects();
		virtual void   UserExec(Option_t *option);
		virtual void   Terminate(Option_t *);

  		const AliESDVertex *GetEventVertex(AliESDEvent* esd);
		Double_t Energy(AliESDtrack* track) const;
		Double_t fuRapidity(AliESDtrack* track) const;
		Bool_t IsGoodTrack(AliESDtrack* ESDTrack) const;
		Bool_t IsPrimaryTrack(AliESDtrack* ESDTrack) const;

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
		TH2F* hPmd;
		TH1F* hMinvPimu;
		TH2F* hdEdx;
		TH1F* hPtPosRSelected;
		TH1F* hPtZSelected;
		TH1F* hPtAngSelected;
		TH1F* hPtPmSelected;
		TH1F* hPtGoodKink; 
		TH1F* hEtaGoodKink; 
		TH1F* hRapidityGoodKink; 
		TH1F* hQtGoodKink; 
		TH2F* hPmGoodKinkAng;  
		TH2F* hPmdGoodKink;
		TH2F* hdEdxGoodKink;
		TH1F* hPtQtSelected;
		TH1F* hPtMaxAngSelected;
		TH1F* hPtRTPCclustersSelected;
		TH2F* hRTPCclustersRTPCclustersSelected;
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
		TH2F* hPmdSelected;
		TH1F* hMinvPimuSelected;  
		TH2F* hdEdxSelected; 
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
		TH2F* hPmdPiSelected;
		TH1F* hMinvPimuPiSelected;  
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
		TH2F* hPmdPiSelectedPlus;
		TH1F* hMinvPimuPiSelectedPlus;  
		TH2F* hdEdxPiSelectedPlus;
		TH1F* hPtPiSelectedMinus;//minus
		TH1F* hEtaPiSelectedMinus; 
		TH1F* hRapidityPiSelectedMinus; 
		TH1F* hQtPiSelectedMinus; 
		TH1F* hKinkAnglePiSelectedMinus;
		TH1F* hDCAkinkPiSelectedMinus;
		TH2F* hPmKinkAngPiSelectedMinus; 
		TH2F* hKinkPosXYPiSelectedMinus; 
		TH2F* hKinkPosZRPiSelectedMinus;  
		TH2F* hPmdPiSelectedMinus;
		TH1F* hMinvPimuPiSelectedMinus;  
		TH2F* hdEdxPiSelectedMinus; // reconstruction histograms
		
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

		
		AliAnalysisPionKinksESD(const AliAnalysisPionKinksESD&); 
		AliAnalysisPionKinksESD& operator=(const AliAnalysisPionKinksESD&);
  
	ClassDef(AliAnalysisPionKinksESD, 1); 
};
#endif	


