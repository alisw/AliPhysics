#ifndef ALICONVERSIONCUTS_H
#define ALICONVERSIONCUTS_H

// Class handling all kinds of selection cuts for Gamma Conversion analysis
// Authors: Svein Lindal, Daniel Lohner												*

#include "AliAODpidUtil.h"
#include "AliConversionPhotonBase.h"
#include "AliAODConversionMother.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliVTrack.h"
#include "AliAODTrack.h"
#include "AliStack.h"
#include "AliAnalysisCuts.h"
#include "TH1F.h"

class AliESDEvent;
class AliAODEvent;
class AliConversionPhotonBase;
class AliKFVertex;
class TH1F;
class TH2F;
class AliPIDResponse;
class AliAnalysisCuts;
class iostream;
class TList;
class AliAnalysisManager;


using namespace std;

class AliConversionCuts : public AliAnalysisCuts {
	
 public: 


  enum cutIds {
	kgoodId=0, 
	kv0FinderType,                
	keProbCut,                    
	kededxSigmaCut,               
	kpidedxSigmaCut,              
	kpiMomdedxSigmaCut,           
	kchi2GammaCut,                
	ksinglePtCut,                 
	kclsTPCCut,                   
	ketaCut,                      
	kchi2MesonCut,                
	kLowPRejectionSigmaCut,       
	kQtMaxCut,                    
	kpiMaxMomdedxSigmaCut,        
	kalphaMesonCut,               
	kminRCut,                     
	kRapidityMesonCut,            
	kBackgroundScheme,            
	kDegreesForRotationMethod,    
	kNumberOfRotations,           
	kremovePileUp,                
	kselectV0AND,                 
	kmultiplicityBin,             
	kisHeavyIon,                  
	kuseCentrality,               
	kcentralityBin,               
	kTOFelectronPID,              
	kuseMCPSmearing,              
	kdoPhotonAsymmetryCut,
	kPsiPair, 
	kCosPAngle,
     //   kHBTmultiplicityBin,
     //   kprimaryCutNumber,
     //   kuseBayesPID,
     //   kdalitzelectronsPID,
     //   kpsiCutNumber,
     //   kdalitzBackgroundType,
     //   keleclsTPCCut,
	kNCuts
  };

  enum photonCuts {
      kPhotonIn=0,
      kOnFly,
      kNoTracks,
      kTrackCuts,
      kdEdxCuts,
      kPhotonCuts,
      kPhotonOut
  };


  Bool_t SetCutIds(TString cutString); 
  Int_t fCuts[kNCuts];
  Bool_t SetCut(cutIds cutID, Int_t cut);
  Bool_t UpdateCutString(cutIds cutID, Int_t value);


  static const char * fgkCutNames[kNCuts];

  Double_t GetPsiPair(const AliESDv0* v0, const AliESDEvent * event) const;
  Double_t GetCosineOfPointingAngle(const AliConversionPhotonBase * photon, const AliVEvent * event) const; 


  Bool_t InitializeCutsFromCutString(const TString analysisCutSelection);

  AliConversionCuts(const char *name="V0Cuts", const char * title="V0 Cuts");
  virtual ~AliConversionCuts();                            //virtual destructor

  virtual Bool_t IsSelected(TObject* /*obj*/){return kTRUE;}
  virtual Bool_t IsSelected(TList* /*list*/) {return kTRUE;}

  TString GetCutNumber();

  // Cut Selection
  Bool_t EventIsSelected(AliVEvent *fInputEvent);
  Bool_t PhotonIsSelected(AliConversionPhotonBase * photon, AliVEvent  * event);
  Bool_t PhotonIsSelectedMC(TParticle *particle,AliStack *fMCStack);
  Bool_t TracksAreSelected(AliVTrack * negTrack, AliVTrack * posTrack);
  Bool_t MesonIsSelected(AliAODConversionMother *pi0,Bool_t IsSignal=kTRUE);
  Bool_t MesonIsSelectedMC(TParticle *fMCMother,AliStack *fMCStack,Bool_t bMCDaughtersInAcceptance=kFALSE);

  void InitAODpidUtil(Int_t type);
  static AliConversionCuts * GetStandardCuts2010PbPb();
  Bool_t InitPIDResponse();
  
  void SetPIDResponse(AliPIDResponse * pidResponse) {fPIDResponse = pidResponse;}
  AliPIDResponse * GetPIDResponse() { return fPIDResponse;}
  
  void PrintCuts();

  void InitCutHistograms();
  void SetFillCutHistograms(){if(!fHistograms){InitCutHistograms();};}
  TList *GetCutHistograms(){return fHistograms;}
  void FillPhotonCutIndex(Int_t photoncut){if(hCutIndex)hCutIndex->Fill(photoncut);}

  AliVTrack * GetTrack(AliVEvent * event, Int_t label) const;

  ///Cut functions
  Bool_t SpecificTrackCuts(AliAODTrack * negTrack, AliAODTrack * posTrack,Int_t &cutIndex);
  Bool_t SpecificTrackCuts(AliESDtrack * negTrack, AliESDtrack * posTrack,Int_t &cutIndex);
  Bool_t AcceptanceCuts(AliConversionPhotonBase *photon);
  Bool_t AcceptanceCut(TParticle *particle, TParticle * ePos,TParticle* eNeg);
  Bool_t dEdxCuts(AliVTrack * track);
  Bool_t ArmenterosQtCut(AliConversionPhotonBase *photon);
  Bool_t AsymmetryCut(AliConversionPhotonBase *photon,AliVEvent *event);
  Bool_t PIDProbabilityCut(AliConversionPhotonBase *photon, AliVEvent * event);
  Bool_t SelectV0Finder(Bool_t onfly){return onfly&&fUseOnFlyV0Finder;}
  Bool_t PhotonCuts(AliConversionPhotonBase *photon,AliVEvent *event);
  Bool_t CorrectedTPCClusterCut(AliConversionPhotonBase *photon, AliVEvent * event);
  Bool_t PsiPairCut(const AliConversionPhotonBase * photon, const AliVEvent * event) const;
  Bool_t CosinePAngleCut(const AliConversionPhotonBase * photon, const AliVEvent * event) const;
  // Event Cuts
  Bool_t IsCentralitySelected(AliVEvent *fInputEvent);
  Double_t GetCentrality(AliVEvent *event);
  Int_t GetNumberOfContributorsVtx(AliVEvent *event);
  Bool_t VertexZCut(AliVEvent *fInputEvent);

  // Set Individual Cuts
  Bool_t SetRCut(Int_t RCut);
  Bool_t SetV0Finder(Int_t v0FinderType);
  Bool_t SetElectronProbCut(Int_t eProbCut);
  Bool_t SetChi2GammaCut(Int_t chi2GammaCut);
  Bool_t SetTPCdEdxCutPionLine(Int_t pidedxSigmaCut);
  Bool_t SetTPCdEdxCutElectronLine(Int_t ededxSigmaCut);
  Bool_t SetSinglePtCut(Int_t singlePtCut);
  Bool_t SetTPCClusterCut(Int_t clsTPCCut);
  Bool_t SetEtaCut(Int_t etaCut);
  Bool_t SetChi2MesonCut(Int_t chi2MesonCut);
  Bool_t SetMinMomPiondEdxCut(Int_t piMinMomdedxSigmaCut);
  Bool_t SetMaxMomPiondEdxCut(Int_t piMaxMomdedxSigmaCut);
  Bool_t SetLowPRejectionCuts(Int_t LowPRejectionSigmaCut);
  Bool_t SetQtMaxCut(Int_t QtMaxCut);
  Bool_t SetAlphaMesonCut(Int_t alphaMesonCut);
  Bool_t SetTOFElectronPIDCut(Int_t TOFelectronPID);
  Bool_t SetTRDElectronCut(Int_t TRDElectronCut);
  Bool_t SetRapidityMesonCut(Int_t RapidityMesonCut);
  Bool_t SetUseCentrality(Int_t useCentrality);
  Bool_t SetIsHeavyIon(Int_t isHeavyIon);
  Bool_t SetCentralityBin(Int_t centralityBin);
  Bool_t SetPhotonAsymmetryCut(Int_t doPhotonAsymmetryCut);
  Bool_t SetRemovePileUp(Int_t removePileUp);
  Bool_t SetBackgroundScheme(Int_t BackgroundScheme);
  Bool_t SetNDegreesForRotationMethod(Int_t DegreesForRotationMethod);
  Bool_t SetNumberOfRotations(Int_t NumberOfRotations);
  Bool_t SetMCPSmearing(Int_t useMCPSmearing);
  Bool_t SetMultiplicityBin(Int_t multiplicityBin);
  Bool_t SetSelectV0AND(Int_t selectV0AND);
  Bool_t SetCosPAngleCut(Int_t cosCut);
  Bool_t SetPsiPairCut(Int_t psiCut);

  // Request Flags

  Bool_t UseRotationMethod(){return fUseRotationMethodInBG;}
  Int_t NumberOfRotationEvents(){return fnumberOfRotationEventsForBG;}
  Int_t NDegreesRotation(){return fnDegreeRotationPMForBG;}
  Bool_t IsHeavyIon(){return fIsHeavyIon;}
  Bool_t UseTrackMultiplicity(){return fUseTrackMultiplicityForBG;}


  Int_t GetFirstTPCRow(Double_t radius);
  

  protected:
  TList *fHistograms;
  AliPIDResponse *fPIDResponse;

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
  Double_t fChi2CutMeson; //chicut meson
  Double_t fPIDProbabilityCutNegativeParticle;
  Double_t fPIDProbabilityCutPositiveParticle;
  Bool_t   fDodEdxSigmaCut; // flag to use the dEdxCut based on sigmas
  Bool_t   fDoTOFsigmaCut; // flag to use TOF pid cut RRnewTOF
  Double_t fPIDTRDEfficiency; // required electron efficiency for TRD PID
  Bool_t   fDoTRDPID; // flag to use TRD pid
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
  Bool_t fUseEtaMinCut; //flag
  Bool_t fUseOnFlyV0Finder; //flag
  Bool_t   fDoPhotonAsymmetryCut; // flag to use the PhotonAsymetryCut
  Double_t fMinPPhotonAsymmetryCut; // Min Momentum for Asymmetry Cut
  Double_t fMinPhotonAsymmetry;  // Asymmetry Cut
  Bool_t fIsHeavyIon;               // flag for heavy ion
  Double_t fMaxVertexZ;    // max z offset of vertex
  Int_t fUseCentrality;  // centrality selection
  Int_t fUseCentralityBin; // centrality selection with individual bins
  Bool_t fUseCorrectedTPCClsInfo; // flag to use corrected tpc cl info
  Bool_t fUseTOFpid; // flag to use tof pid
  Double_t fAlphaMinCutMeson; // min value for meson alpha cut
  Double_t fAlphaCutMeson; // max value for meson alpha cut
  Double_t fRapidityCutMeson; // max value for meson rapidity
  Bool_t fUseRotationMethodInBG; // flag to apply rotation method for meson bg estimation
  Bool_t fdoBGProbability; // flag to use probability method for meson bg estimation
  Bool_t fUseTrackMultiplicityForBG; // flag to use track multiplicity for meson bg estimation (else V0 mult)
  Int_t fnDegreeRotationPMForBG; //
  Int_t fnumberOfRotationEventsForBG; //
  Bool_t fUseMCPSmearing; // flag
  Double_t fPBremSmearing;//
  Double_t fPSigSmearing; //
  Double_t fPSigSmearingCte; //
  Bool_t fUseMultiplicity; // flag
  Int_t fUseMultiplicityBin; // selected multiplicity bin
  Bool_t fSelectV0AND; // flag
  Bool_t fRemovePileUp; //flag
  Float_t fOpeningAngle; // min opening angle for meson
  Float_t fPsiPairCut;
  Float_t fCosPAngleCut;

  // Histograms
  TObjString *fCutString; // cut number used for analysis
  TH1F *hdEdxCuts;  // bookkeeping for dEdx cuts
  TH2F *hTPCdEdxbefore; // TPC dEdx before cuts
  TH2F *hTPCdEdxafter; // TPC dEdx after cuts
  TH1F *hTrackCuts; // bookkeeping for track cuts
  TH1F *hPhotonCuts; // bookkeeping for photon specific cuts
  TH1F *hInvMassbefore; // e+e- inv mass distribution before cuts
  TH2F *hArmenterosbefore; // armenteros podolanski plot before cuts
  TH1F *hInvMassafter; // e+e- inv mass distribution after cuts
  TH2F *hArmenterosafter;  // armenteros podolanski plot after cuts
  TH1F *hAcceptanceCuts; // bookkeeping for acceptance cuts
  TH1F *hCutIndex; // bookkeeping for cuts
  TH1F *hV0EventCuts; // bookkeeping for event selection cuts
  TH1F *hCentrality; // centrality distribution for selected events
  TH1F *hVertexZ; // vertex z distribution for selected events
  TH1F *hMesonCuts; // bookkeeping for meson cuts
  TH1F *hMesonBGCuts; // bookkeeping for meson bg cuts


private:

  AliConversionCuts(const AliConversionCuts&); // not implemented
  AliConversionCuts& operator=(const AliConversionCuts&); // not implemented


  ClassDef(AliConversionCuts,2)
};


inline void AliConversionCuts::InitAODpidUtil(Int_t type) {
  if (!fPIDResponse) fPIDResponse = new AliAODpidUtil();
  Double_t alephParameters[5];
  // simulation
  alephParameters[0] = 2.15898e+00/50.;
  alephParameters[1] = 1.75295e+01;
  alephParameters[2] = 3.40030e-09;
  alephParameters[3] = 1.96178e+00;
  alephParameters[4] = 3.91720e+00;
  fPIDResponse->GetTOFResponse().SetTimeResolution(80.);
  
  // data
  if (type==1){
    alephParameters[0] = 0.0283086/0.97;
    alephParameters[1] = 2.63394e+01;
    alephParameters[2] = 5.04114e-11;
    alephParameters[3] = 2.12543e+00;
    alephParameters[4] = 4.88663e+00;
    fPIDResponse->GetTOFResponse().SetTimeResolution(130.);
    fPIDResponse->GetTPCResponse().SetMip(50.);
  }
  
  fPIDResponse->GetTPCResponse().SetBetheBlochParameters(
    alephParameters[0],alephParameters[1],alephParameters[2],
    alephParameters[3],alephParameters[4]);
  
  fPIDResponse->GetTPCResponse().SetSigma(3.79301e-03, 2.21280e+04);
}


#endif
