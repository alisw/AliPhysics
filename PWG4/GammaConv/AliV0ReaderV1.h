#ifndef ALIV0READERV1_H
#define ALIV0READERV1_H

#include "AliAnalysisTaskSE.h"
#include "TClonesArray.h"
#include "AliKFParticle.h"
#include "AliAODv0.h"
#include "AliESDv0.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliKFConversionPhoton.h"
#include "AliConversionPhotonBase.h"
#include "AliConversionAODBGHandlerRP.h"
#include "TVector.h"
#include "AliKFVertex.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliVTrack.h"
#include "AliStack.h"
#include "AliEventplane.h"

#include "TH1.h"
#include "TH2.h"

class AliAODConversionPhoton;

using namespace std;

class AliV0ReaderV1 : public AliAnalysisTaskSE {
	
 public: 
	
  AliV0ReaderV1(const char *name="V0ReaderV1");
  virtual ~AliV0ReaderV1();                            //virtual destructor
  void UserCreateOutputObjects();

  virtual void   UserExec(Option_t *option);
  virtual void Terminate(Option_t *);

  // Reconstruct Gammas
  void ProcessV0();

  Bool_t CheckV0Status();

  void ProcessMC(AliKFConversionPhoton*);
  void ProcessMCGammasForEfficiency();

  void SetIsHeavyIon(Bool_t heavyion){fIsHeavyIon=heavyion;if(!fIsHeavyIon){fNCentralityBins=1;}}

  // Return Reconstructed Gammas
  TClonesArray *GetReconstructedGammas(){return fConversionGammas;}
  Int_t GetNReconstructedGammas(){if(fConversionGammas){return fConversionGammas->GetEntriesFast();}else{return 0;}}

  // Process Gamma Candidates

  Bool_t IsGammaCandidate(AliConversionPhotonBase *fPhotonCandidate);
  Bool_t IsMCConversionGammaInAcceptance(TParticle *particle);

  void PrintCuts();

  Bool_t CentralitySelection();
  Bool_t VertexZCut();
  Double_t GetCentrality();
  Bool_t SetEventPlane();
  AliEventplane *GetEventPlane();
  AliVTrack *GetTrack(Int_t label);

  // Getter Functions

  const AliExternalTrackParam *GetExternalTrackParam(Int_t charge);
  const AliExternalTrackParam *GetExternalTrackParamP(){return GetExternalTrackParam(1);};
  const AliExternalTrackParam *GetExternalTrackParamN(){return GetExternalTrackParam(-1);};


  Bool_t GetConversionPoint(const AliExternalTrackParam *pparam,const AliExternalTrackParam *nparam,Double_t convpos[3]);
  Bool_t GetHelixCenter(const AliExternalTrackParam *track, Double_t b,Int_t charge, Double_t center[2]);

  Double_t GetMagneticField() const{if(fESDEvent){return fESDEvent->GetMagneticField();}if(fAODEvent){return fAODEvent->GetMagneticField();}return 0;}


  const AliVertex *GetPrimaryVertex();

  AliKFParticle *GetPositiveKFParticle(AliAODv0 *fCurrentV0,Int_t fTrackLabel[2]);
  AliKFParticle *GetNegativeKFParticle(AliAODv0 *fCurrentV0,Int_t fTrackLabel[2]);
  AliKFParticle *GetPositiveKFParticle(AliESDv0 *fCurrentV0,Int_t fTrackLabel[2]);
  AliKFParticle *GetNegativeKFParticle(AliESDv0 *fCurrentV0,Int_t fTrackLabel[2]);

  Int_t GetNumberOfContributorsVtx();

  // Cut Functions

  Bool_t EventIsSelected(){return fEventIsSelected;}
  Bool_t TrackCuts(AliConversionPhotonBase *fPhotonCandidate);
  Bool_t AcceptanceCuts(AliConversionPhotonBase *fPhotonCandidate);
  Bool_t AcceptanceCut(TParticle *particle, TParticle * ePos,TParticle* eNeg);
  Bool_t dEdxCuts(AliConversionPhotonBase *fPhotonCandidate);
  Bool_t ArmenterosQtCut(AliConversionPhotonBase *fPhotonCandidate);
  Bool_t AsymmetryCut(AliConversionPhotonBase *fPhotonCandidate);
  Bool_t PIDProbabilityCut(AliConversionPhotonBase *fPhotonCandidate);
  Bool_t EventCuts();

  // Set Mass Zero for Gammas
  void SetGammaMassZero(){fCurrentMotherKFCandidate->E()=fCurrentMotherKFCandidate->GetP();}


  // Set Output

  void SetUseAODConversionPhoton(Bool_t b){if(b){cout<<"Setting Outputformat to AliAODConversionPhoton "<<endl;}else{cout<<"Setting Outputformat to AliKFConversionPhoton "<<endl;};kUseAODConversionPhoton=b;}

  void SetCreateAODs(Bool_t k){fCreateAOD=k;}
  void SetDeltaAODFilename(TString s){fDeltaAODFilename=s;}




  // Set Cuts

  void SetMaxVertexZ(Double_t maxVertexZ){fMaxVertexZ=maxVertexZ;}

  void SetUseImprovedVertex(Bool_t useImprovedVertex){fUseImprovedVertex=useImprovedVertex;}

   void SetOnFlyFlag(Bool_t flag){fUseOnFlyV0Finder = flag;}

  void SetMaxRCut(Double_t maxR){fMaxR=maxR;}
  void SetMinRCut(Double_t minR){fMinR=minR;}

  void SetEtaCut(Double_t etaCut){fEtaCut=etaCut;}

  void SetEtaCutMin(Double_t etaCutMin){fEtaCutMin=etaCutMin;}

  void SetPtCut(Double_t ptCut){fPtCut=ptCut;}
	
  void SetSinglePtCut(Double_t singleptCut){fSinglePtCut=singleptCut;}
	
  void SetMaxZCut(Double_t maxZ){fMaxZ=maxZ;}
	
  void SetMinClsTPCCut(Double_t minClsTPC){fMinClsTPC=minClsTPC;}

  void SetMinClsTPCCutToF(Double_t minClsTPCToF){fMinClsTPCToF=minClsTPCToF;}
	
  void SetLineCutZRSlope(Double_t LineCutZRSlope){fLineCutZRSlope=LineCutZRSlope;}
  void SetLineCutZValue(Double_t LineCutZValue){fLineCutZValue=LineCutZValue;}
	
  void SetLineCutZRSlopeMin(Double_t LineCutZRSlopeMin){fLineCutZRSlopeMin=LineCutZRSlopeMin;}
  void SetLineCutZValueMin(Double_t LineCutZValueMin){fLineCutZValueMin=LineCutZValueMin;}
		
  void SetChi2CutConversion(Double_t chi2){fChi2CutConversion=chi2;}

  void SetPIDProbability(Double_t pidProb){fPIDProbabilityCutPositiveParticle=pidProb; fPIDProbabilityCutNegativeParticle=pidProb;}
	
  void SetPIDProbabilityNegativeParticle(Double_t pidProb){fPIDProbabilityCutNegativeParticle=pidProb;}
	
  void SetPIDProbabilityPositiveParticle(Double_t pidProb){fPIDProbabilityCutPositiveParticle=pidProb;}

  void SetPIDnSigmaAboveElectronLine(Double_t nSigmaAbove){fPIDnSigmaAboveElectronLine=nSigmaAbove;}
  void SetTofPIDnSigmaAboveElectronLine(Double_t nTofSigmaAbove){fTofPIDnSigmaAboveElectronLine=nTofSigmaAbove;} // RRnewTOF
	
  void SetPIDnSigmaBelowElectronLine(Double_t nSigmaBelow){fPIDnSigmaBelowElectronLine=nSigmaBelow;}
  void SetTofPIDnSigmaBelowElectronLine(Double_t nTofSigmaBelow){fTofPIDnSigmaBelowElectronLine=nTofSigmaBelow;} // RRnewTOF
	
  /*
   * Sets the PIDnSigmaAbovePion cut value for the tracks.
   */
  void SetPIDnSigmaAbovePionLine(Double_t nSigmaAbovePion){fPIDnSigmaAbovePionLine=nSigmaAbovePion;}

  /*
   * Sets the PIDnSigmaAbovePion cut value for the tracks.
   */
  void SetPIDnSigmaAbovePionLineHighPt(Double_t nSigmaAbovePionHighPt){fPIDnSigmaAbovePionLineHighPt=nSigmaAbovePionHighPt;}

  /*
   * Sets the PIDMinPnSigmaAbovePion cut value for the tracks.
   */
  void SetPIDMinPnSigmaAbovePionLine(Double_t MinPnSigmaAbovePion){fPIDMinPnSigmaAbovePionLine=MinPnSigmaAbovePion;}

 /*
   * Sets the PIDMinPnSigmaAbovePion cut value for the tracks.
   */
  void SetPIDMaxPnSigmaAbovePionLine(Double_t MaxPnSigmaAbovePion){fPIDMaxPnSigmaAbovePionLine=MaxPnSigmaAbovePion;}

  /*
   * Sets the SigmaMassCut value.
   */
  void SetSigmaMass(Double_t sigmaMass){fNSigmaMass=sigmaMass;}
	
  void SetTRDEfficiency(Double_t eff){if(eff<=1&&eff>=0){fPIDTRDEfficiency=eff;}}
  void SetDoTRDPID(Bool_t doTRD){fDoTRDPID=doTRD;}

  void SetDodEdxSigmaCut( Bool_t dodEdxSigmaCut){fDodEdxSigmaCut=dodEdxSigmaCut;}
  void SetDoTOFsigmaCut( Bool_t doTOFsigmaCut){fDoTOFsigmaCut=doTOFsigmaCut;} //RRnewTOF
  void SetDoPhotonAsymmetryCut( Bool_t doPhotonAsymmetryCut){fDoPhotonAsymmetryCut=doPhotonAsymmetryCut;}

  void SetMinPPhotonAsymmetryCut(Double_t minPPhotonAsymmetryCut){fMinPPhotonAsymmetryCut=minPPhotonAsymmetryCut;}
  void SetMinPhotonAsymmetry(Double_t minPhotonAsymmetry){fMinPhotonAsymmetry=minPhotonAsymmetry;}
  /*
   * Sets the flag to enable/disable the cut dedx N sigma for Kaon Rejection at low p 
   */
  void SetDoKaonRejectionLowP( Bool_t doKaonRejectionLowP){fDoKaonRejectionLowP=doKaonRejectionLowP;}
  /*
   * Sets the flag to enable/disable the cut dedx N sigma for Proton Rejection at low p 
   */
  void SetDoProtonRejectionLowP( Bool_t doProtonRejectionLowP){fDoProtonRejectionLowP=doProtonRejectionLowP;}

  /*
   * Sets the flag to enable/disable the cut dedx N sigma for Pion Rejection at low p 
   */
  void SetDoPionRejectionLowP( Bool_t doPionRejectionLowP){fDoPionRejectionLowP=doPionRejectionLowP;}

  /*
   * Sets the PIDMinPnSigmaAroundKaon cut value for the tracks.
   */
  void SetPIDnSigmaAtLowPAroundKaonLine(Double_t nSigmaAtLowPAroundKaon){fPIDnSigmaAtLowPAroundKaonLine =nSigmaAtLowPAroundKaon;}

  /*
   * Sets the PIDMinPnSigmaAroundProton cut value for the tracks.
   */
  void SetPIDnSigmaAtLowPAroundProtonLine(Double_t nSigmaAtLowPAroundProton){fPIDnSigmaAtLowPAroundProtonLine =nSigmaAtLowPAroundProton;}

  /*
   * Sets the PIDMinPnSigmaAroundPion cut value for the tracks.
   */
  void SetPIDnSigmaAtLowPAroundPionLine(Double_t nSigmaAtLowPAroundPion){fPIDnSigmaAtLowPAroundPionLine =nSigmaAtLowPAroundPion;}

  /*
   * Sets the PIDMinPnSigmaAbovePion cut value for the tracks.
   */
  void SetPIDMinPKaonRejectionLowP(Double_t PIDMinPKaonRejectionLowP ){fPIDMinPKaonRejectionLowP=PIDMinPKaonRejectionLowP;}

  /*
   * Sets the PIDMinPnSigmaAbovePion cut value for the tracks.
   */
  void SetPIDMinPProtonRejectionLowP(Double_t PIDMinPProtonRejectionLowP ){fPIDMinPProtonRejectionLowP=PIDMinPProtonRejectionLowP;}
  /*
   * Sets the PIDMinPnSigmaAbovePion cut value for the tracks.
   */
  void SetPIDMinPPionRejectionLowP(Double_t PIDMinPPionRejectionLowP ){fPIDMinPPionRejectionLowP=PIDMinPPionRejectionLowP;}

  /*
   *Set if we want to use Gamma Selection based on Qt from Armenteros
   */
  void SetDoQtGammaSelection(Bool_t doQtGammaSelection){fDoQtGammaSelection=doQtGammaSelection;}
  void SetDoHighPtQtGammaSelection(Bool_t doHighPtQtGammaSelection){fDoHighPtQtGammaSelection=doHighPtQtGammaSelection;} // RRnew

  void SetQtMax(Double_t qtMax){fQtMax=qtMax;}
  void SetHighPtQtMax(Double_t qtMaxHighPt){fHighPtQtMax=qtMaxHighPt;} // RRnew
  void SetPtBorderForQt(Double_t ptBorderForQt){fPtBorderForQt=ptBorderForQt;} // RRnew

  void SetUseOwnXYZCalculation(Bool_t flag){fUseOwnXYZCalculation=flag;}

  void SetUseConstructGamma(Bool_t flag){fUseConstructGamma=flag;}

protected:
    TClonesArray *fConversionGammas;
    AliESDEvent *fESDEvent;
    AliAODEvent *fAODEvent;
    AliStack *fMCStack;
    TList *fOutputList;
    AliKFConversionPhoton *fCurrentMotherKFCandidate;
    AliKFParticle *fCurrentPositiveKFParticle;
    AliKFParticle *fCurrentNegativeKFParticle;
    Int_t *fCurrentTrackLabels;
    Int_t fCurrentV0Index;
    Int_t fNCentralityBins;
    Int_t fCentralityBin;
    Bool_t fEventIsSelected;
    AliConversionAODBGHandlerRP *fBGHandler;
    Double_t fVertexZ;
    Double_t fCentrality;
    Double_t fEPAngle;

private:
    //Event Cuts
    Double_t fMaxVertexZ; // max z vertex cut
  //cuts
  Double_t fMaxR; //r cut
  Double_t fMinR; //r cut
  Double_t fEtaCut; //eta cut
  Double_t fEtaCutMin; //eta cut
  Double_t fPtCut; // pt cut
  Double_t fSinglePtCut; // pt cut for electron/positron
  Double_t fMaxZ; //z cut
  Double_t fMinClsTPC; // minimum clusters in the TPC
  Double_t fMinClsTPCToF; // minimum clusters to findable clusters
  Double_t fLineCutZRSlope; //linecut
  Double_t fLineCutZValue; //linecut
  Double_t fLineCutZRSlopeMin; //linecut
  Double_t fLineCutZValueMin; //linecut
  Double_t fChi2CutConversion; //chi2cut
  Double_t fPIDProbabilityCutNegativeParticle;
  Double_t fPIDProbabilityCutPositiveParticle;
  Bool_t   fDodEdxSigmaCut; // flag to use the dEdxCut based on sigmas
  Bool_t   fDoTOFsigmaCut; // flag to use TOF pid cut RRnewTOF
  Double_t fPIDTRDEfficiency;
  Bool_t fDoTRDPID;
  Double_t fPIDnSigmaAboveElectronLine; // sigma cut
  Double_t fPIDnSigmaBelowElectronLine; // sigma cut
  Double_t fTofPIDnSigmaAboveElectronLine; // sigma cut RRnewTOF
  Double_t fTofPIDnSigmaBelowElectronLine; // sigma cut RRnewTOF 
  Double_t fPIDnSigmaAbovePionLine;     // sigma cut
  Double_t fPIDnSigmaAbovePionLineHighPt;     // sigma cut
  Double_t fPIDMinPnSigmaAbovePionLine; // sigma cut
  Double_t fPIDMaxPnSigmaAbovePionLine; // sigma cut
  Double_t fDoKaonRejectionLowP;   // Kaon rejection at low p
  Double_t fDoProtonRejectionLowP; // Proton rejection at low p
  Double_t fDoPionRejectionLowP;   // Pion rejection at low p
  Double_t fPIDnSigmaAtLowPAroundKaonLine; // sigma cut
  Double_t fPIDnSigmaAtLowPAroundProtonLine; // sigma cut
  Double_t fPIDnSigmaAtLowPAroundPionLine; // sigma cut
  Double_t fPIDMinPKaonRejectionLowP; // Momentum limit to apply kaon rejection
  Double_t fPIDMinPProtonRejectionLowP; // Momentum limit to apply proton rejection
  Double_t fPIDMinPPionRejectionLowP; // Momentum limit to apply proton rejection
  Bool_t   fDoQtGammaSelection; // Select gammas using qtMax
  Bool_t   fDoHighPtQtGammaSelection; // RRnew Select gammas using qtMax for high pT
  Double_t fQtMax; // Maximum Qt from Armenteros to select Gammas
  Double_t fHighPtQtMax; // RRnew Maximum Qt for High pT from Armenteros to select Gammas
  Double_t fPtBorderForQt; // RRnew 
  Double_t fXVertexCut; //vertex cut
  Double_t fYVertexCut; //vertex cut
  Double_t fZVertexCut; // vertexcut
	
  Double_t fNSigmaMass; //nsigma cut
	
  Bool_t fUseImprovedVertex; //flag

  Bool_t fUseOwnXYZCalculation; //flag that determines if we use our own calculation of xyz (markus)

  Bool_t fUseConstructGamma; //flag that determines if we use ConstructGamma method from AliKF

  Bool_t fUseEtaMinCut; //flag

  Bool_t fUseOnFlyV0Finder; //flag

  Bool_t   fDoPhotonAsymmetryCut; // flag to use the PhotonAsymetryCut
  Double_t fMinPPhotonAsymmetryCut;
  Double_t fMinPhotonAsymmetry;

  Bool_t kUseAODConversionPhoton;

  Bool_t fIsHeavyIon;
  Bool_t fCreateAOD;
  TString fDeltaAODFilename;



  // Histograms
  
  TH1F *hV0CurrentFinder;
  TH2F *hV0AllArmenteros;
  TH1F *hV0Good;
  TH2F *hV0GoodArmenteros;

  TH1F *hV0CutLikeSign;
  TH1F *hV0CutRefit;
  TH1F *hV0CutKinks;
  TH1F *hV0CutMinNclsTPCToF;

  TH1F *hV0CutNContributors;
  TH1F *hV0CutVertexZ;

  TH1F *hV0CutdEdxElectron;
  TH1F *hV0CutdEdxPion;
  TH1F *hV0CutdEdxKaonLowP;
  TH1F *hV0CutdEdxProtonLowP;
  TH1F *hV0CutdEdxPionLowP;
  TH1F *hV0CutdEdxTOFElectron;
  TH1F * hV0CutdEdxTRD;

  TH1F *hV0CutQt;

  TH1F *hV0CutR;
  TH1F *hV0CutMinR;
  TH1F *hV0CutLine;
  TH1F *hV0CutZ;
  TH1F *hV0CutEta;
  TH1F *hV0CutSinglePt;
  TH1F *hV0CutPt;
  TH1F *hV0CutNDF;
  TH1F *hV0CutChi2;

  TH1F *hV0CutPIDProb;
  TH1F* hV0CutAsymmetry;

  TH2F *hMCPtResolution;
  TH1F **hGammaPt;
  TH1F **hMCPtRECOTRUE;
  TH1F **hMCPtTRUE;
  TH2F *hMCPtResolutionPhi;
  TH2F *hMCRResolutionvsR;
  TH2F *hMCZResolutionvsZ;
  TH1F *hGammaPhi;
  TH2F *hGammadEdxbefore;
  TH2F *hGammadEdxafter;
  TH2F *hGammaConversionMapXY;
  TH2F *hGammaConversionMapZR;

  // Event

  TH1F *hV0EventCuts;
  TH1F *hNEvents;
  TH1F *hCentrality;
  TH1F *hVertexZ;


  ClassDef(AliV0ReaderV1,1)
};


#endif





