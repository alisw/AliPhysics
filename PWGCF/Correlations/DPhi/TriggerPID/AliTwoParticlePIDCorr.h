#ifndef AliTwoParticlePIDCorr_H
#define AliTwoParticlePIDCorr_H

#include "THn.h" // in cxx file causes .../THn.h:257: error: conflicting declaration ‘typedef class THnT<float> THnF’


class TH1F;
class TH2F;
class TH3F;
class THnSparse;
class TString;
class TList;
//class AliESDtrackCuts;
class TSeqCollection;
class AliPIDResponse;
class AliPIDCombined;  
class AliAODEvent;
class AliVEvent;
class AliAODTrack;
class AliVTrack;
class AliAODv0;
class AliAODVertex;
class AliEventPoolManager;
class TFormula;
class AliAnalysisUtils;
class AliPPVsMultUtils;
class LRCParticlePID;
class AliVParticle;
class AliCFContainer;
class AliCFGridSparse;
class THnBase;
//class AliTHn;//Showing problem as a redefinition of ALiTHn, after commenting running fine
class TProfile;


#include <TObject.h> //LRCParticlePID is a derived class from"TObject"
#include "TMath.h"
#include "TNamed.h"
#include "AliUEHist.h"
#include "AliPID.h"
#include "AliAnalysisTask.h"
#include "AliUEHist.h"
#include "TString.h"
#include "AliVParticle.h"
#include "TParticle.h"
#include "AliLog.h"
#include "AliTHn.h"
#include "TBits.h"


#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

namespace AliPIDNameSpace {
  
  enum PIDTypes
  {
    NSigmaTPC = 0,
    NSigmaTOF,
    NSigmaTPCTOF,  // squared sum
    Bayes    
  };
  const Int_t NSigmaPIDType=NSigmaTPCTOF;//number of Nsigma PID types
    
  enum AliDetectorType
  {
    fITS = 0,
    fTPC,
    fTOF,
    fNDetectors
  };
  
  enum AliParticleSpecies
  {
    SpPion = 0,
    SpKaon,
    SpProton,
    unidentified,
    SpKs0,
    SpKs0_LS_Bckg,
    SpKs0_RS_Bckg,
    SpLam,
    SpLam_LS_Bckg,
    SpLam_RS_Bckg,
    NSpecies=unidentified,//for pion, kaon and proton part only not for v0s
    SpUndefined=999
  }; // Particle species used in plotting
  
  
  enum AliCharge
  {
    Posch = 0,
    Negch,
    NCharge
  };
}


using namespace AliPIDNameSpace;

class AliTwoParticlePIDCorr : public AliAnalysisTaskSE {
 public:
    AliTwoParticlePIDCorr();
    AliTwoParticlePIDCorr(const char *name);
    virtual ~AliTwoParticlePIDCorr();
    
    virtual void     UserCreateOutputObjects();
    virtual void     UserExec(Option_t *option);
    virtual void    doAODevent();
    virtual void    doMCAODevent();
    virtual void     Terminate(Option_t *);
  void	   SetSharedClusterCut(Double_t value) { fSharedClusterCut = value; }
  void	   SetSharedTPCmapCut(Double_t value1) { fSharedTPCmapCut = value1; }
  void	   SetSharedfraction_Pair_cut(Double_t value2) { fSharedfraction_Pair_cut = value2; }


  void SettwoTrackEfficiencyCutDataReco(Bool_t twoTrackEfficiencyCutDataReco,Float_t twoTrackEfficiencyCutValue1,Float_t TwoTrackCutMinRadius,Float_t TwoTrackCutMaxRadius)
  {
    ftwoTrackEfficiencyCutDataReco=twoTrackEfficiencyCutDataReco;
    twoTrackEfficiencyCutValue=twoTrackEfficiencyCutValue1;
    fTwoTrackCutMinRadius=TwoTrackCutMinRadius;
    fTwoTrackCutMaxRadius=TwoTrackCutMaxRadius;
  }
  void SetVertextype(Int_t Vertextype){fVertextype=Vertextype;}                                                 //Check it every time
    void SetZvtxcut(Double_t zvtxcut) {fzvtxcut=zvtxcut;}
    void SetZvtxcut_MC(Double_t VxMax_MC,Double_t VyMax_MC,Double_t VzMax_MC) {
fVxMax_MC=VxMax_MC;
fVyMax_MC=VyMax_MC;
fVzMax_MC=VzMax_MC;
}

    void SetCustomBinning(TString receivedCustomBinning) { fCustomBinning = receivedCustomBinning; }
    void SetMaxNofMixingTracks(Int_t MaxNofMixingTracks) {fMaxNofMixingTracks=MaxNofMixingTracks;}               //Check it every time
  void SetCentralityEstimator(TString CentralityMethod) { fCentralityMethod = CentralityMethod;}
  void SetPPVsMult(Bool_t val, Bool_t PileUp_zvtx_INEL_evsel){
    fPPVsMult = val;
    fPileUp_zvtx_INEL_evsel=PileUp_zvtx_INEL_evsel;
  }
  void SetSampleType(TString SampleType) {fSampleType=SampleType;}
  void SetRequestEventPlane(Bool_t RequestEventPlane,Bool_t RequestEventPlanemixing,Bool_t V2,Bool_t V3,TString EPdetector,Bool_t IsAfter2011){
fRequestEventPlane=RequestEventPlane;
fRequestEventPlanemixing=RequestEventPlanemixing;
fV2=V2;
fV3=V3;
fEPdet=EPdetector;
fIsAfter2011=IsAfter2011;
}
  void SetAnalysisType(TString AnalysisType){fAnalysisType=AnalysisType;}
  void SetFilterBit(Int_t FilterBit) {fFilterBit=FilterBit;}
  void SetTrackStatus(UInt_t status) { fTrackStatus = status; }
  void SetfilltrigassoUNID(Bool_t filltrigassoUNID){ffilltrigassoUNID=filltrigassoUNID;}
  void SetfilltrigUNIDassoID(Bool_t filltrigUNIDassoID){ffilltrigUNIDassoID=filltrigUNIDassoID;}
  void SetfilltrigIDassoUNID(Bool_t filltrigIDassoUNID){ffilltrigIDassoUNID=filltrigIDassoUNID;}
  void SetfilltrigIDassoID(Bool_t filltrigIDassoID){ ffilltrigIDassoID=filltrigIDassoID;}
  void SetfilltrigIDassoIDMCTRUTH(Bool_t filltrigIDassoIDMCTRUTH){ffilltrigIDassoIDMCTRUTH=filltrigIDassoIDMCTRUTH;}
  void SetSelectHighestPtTrig(Bool_t SelectHighestPtTrig){fSelectHighestPtTrig=SelectHighestPtTrig;}
  void SetTriggerSpeciesSelection(Bool_t TriggerSpeciesSelection,Int_t TriggerSpecies,Bool_t containPIDtrig){
    fTriggerSpeciesSelection=TriggerSpeciesSelection;//if it is KTRUE then Set containPIDtrig=kFALSE
    fTriggerSpecies=TriggerSpecies;
    fcontainPIDtrig=containPIDtrig;
  }
  void SetAssociatedSpeciesSelection(Bool_t AssociatedSpeciesSelection,Int_t AssociatedSpecies, Bool_t containPIDasso)
  {
    fAssociatedSpeciesSelection=AssociatedSpeciesSelection;//if it is KTRUE then Set containPIDasso=kFALSE
    fAssociatedSpecies=AssociatedSpecies;
    fcontainPIDasso=containPIDasso;
  }

  void SettingChargeCounting(Int_t val) {SetChargeAxis=val;}

 void SetFIllPIDQAHistos(Bool_t FIllPIDQAHistos){fFIllPIDQAHistos=FIllPIDQAHistos;}
  void SetRejectPileUp(Bool_t rejectPileUp) {frejectPileUp=rejectPileUp;}
  void SetCheckFirstEventInChunk(Bool_t CheckFirstEventInChunk) {fCheckFirstEventInChunk=CheckFirstEventInChunk;}

  void SetKinematicCuts(Float_t minPt, Float_t maxPt,Float_t mineta,Float_t maxeta)
  {
    fminPt=minPt;
    fmaxPt=maxPt;
    fmineta=mineta;
    fmaxeta=maxeta;
  }

  void SetDCACut(TFormula* value) { fDCAXYCut = value; }

  void SetfillHistQA(Bool_t fillhistQAReco,Bool_t fillhistQATruth)
  {
    ffillhistQAReco=fillhistQAReco;
    ffillhistQATruth=fillhistQATruth;
  }
  void SetPtordering(Bool_t PtOrderDataReco,Bool_t PtOrderMCTruth)
  {
    fPtOrderDataReco=PtOrderDataReco;
    fPtOrderMCTruth=PtOrderMCTruth;
  }
  void SetWeightPerEvent(Bool_t flag) {  fWeightPerEvent = flag;}
  void Setselectprimarydatareco(Bool_t onlyprimarydatareco) {fonlyprimarydatareco=onlyprimarydatareco;}
  void SetselectprimaryTruth(Bool_t selectprimaryTruth) {fselectprimaryTruth=selectprimaryTruth;}
  void SetCombinedNSigmaCut(Double_t NSigmaPID) {fNSigmaPID=NSigmaPID;}
  void SetHighPtKaonNSigmaPID(Float_t HighPtKaonSigma,Float_t HighPtKaonNSigmaPID)
  {
    fHighPtKaonSigma=HighPtKaonSigma;
    fHighPtKaonNSigmaPID=HighPtKaonNSigmaPID;
  }
  void IgnoreoverlappedTracks(Bool_t UseExclusiveNSigma){fUseExclusiveNSigma=UseExclusiveNSigma;}
  void SetRemoveTracksT0Fill( Bool_t RemoveTracksT0Fill){fRemoveTracksT0Fill=RemoveTracksT0Fill;}
  void SetPairSelectCharge(Int_t SelectCharge){fSelectCharge=SelectCharge;}
  void SetTrigAssoSelectcharge( Int_t TriggerSelectCharge,Int_t AssociatedSelectCharge)
  {
    fTriggerSelectCharge=TriggerSelectCharge;
    fAssociatedSelectCharge=AssociatedSelectCharge;
  }
  void SetEtaOrdering(Bool_t EtaOrdering){fEtaOrdering=EtaOrdering;}
  void SetCutConversionsResonances( Bool_t CutConversions,Bool_t CutResonances)
  {
     fCutConversions=CutConversions;
     fCutResonances=CutResonances;
  }
  void SetCleanUp(Bool_t InjectedSignals,Bool_t RemoveWeakDecays,Bool_t RemoveDuplicates)
  {
    fInjectedSignals=InjectedSignals;
    fRemoveWeakDecays=RemoveWeakDecays;
    fRemoveDuplicates=RemoveDuplicates;
  }
  void SetEfficiency(Bool_t fillefficiency,Bool_t applyTrigefficiency,Bool_t applyAssoefficiency)
    {
      ffillefficiency=fillefficiency;
      fapplyTrigefficiency=applyTrigefficiency;
      fapplyAssoefficiency=applyAssoefficiency;

    }
  void SetComboeffPtRange(Double_t minPtComboeff,Double_t maxPtComboeff) {
    fminPtComboeff=minPtComboeff;
    fmaxPtComboeff=maxPtComboeff;}
  //only one can be kTRUE at a time(for the next two Setters)
  void Setmesoneffrequired(Bool_t mesoneffrequired) {fmesoneffrequired=mesoneffrequired;}
  void Setkaonprotoneffrequired(Bool_t kaonprotoneffrequired){fkaonprotoneffrequired=kaonprotoneffrequired;}
   void SetOnlyOneEtaSide(Int_t OnlyOneEtaSide){fOnlyOneEtaSide=OnlyOneEtaSide;}
   void SetRejectResonanceDaughters(Int_t RejectResonanceDaughters){fRejectResonanceDaughters=RejectResonanceDaughters;}

void SetTOFPIDVal(Bool_t RequestTOFPID,Float_t PtTOFPIDmin,Float_t PtTOFPIDmax)
   {
fRequestTOFPID=RequestTOFPID;
fPtTOFPIDmin=PtTOFPIDmin;
fPtTOFPIDmax=PtTOFPIDmax;
}
 void SetRunShufflingbalance(Bool_t RunShufflingbalance){
   fRunShufflingbalance=RunShufflingbalance;
 }

 void setElectronRejectionTRUTH(Bool_t ElectronRejectionTRUTH)
 {
   fElectronRejectionTRUTH=ElectronRejectionTRUTH;
 }

 void setElectronRejection(Bool_t ElectronRejection, Bool_t ElectronOnlyRejection, Double_t ElectronRejectionMinPt, Double_t ElectronRejectionMaxPt,Double_t ElectronRejectionNSigma){
   fElectronRejection=ElectronRejection;
   fElectronOnlyRejection=ElectronOnlyRejection;
   fElectronRejectionMinPt=ElectronRejectionMinPt;
   fElectronRejectionMaxPt=ElectronRejectionMaxPt;
   fElectronRejectionNSigma=ElectronRejectionNSigma;   
 }

 void setdeltacutPion(Double_t Piondeltacutmin,Double_t Piondeltacutmax){
   fPiondeltacutmin=Piondeltacutmin;
   fPiondeltacutmax=Piondeltacutmax;
 }
 void setdeltacutProton(Double_t Protondeltacutmin, Double_t Protondeltacutmax){
   fProtondeltacutmin=Protondeltacutmin;
   fProtondeltacutmax=Protondeltacutmax;
 }

 void SetEffcorectionfilePathName(TString efffilename) {fefffilename=efffilename;}
 void SetV0MeanSigmafilePathName(TString efffilename) {ffilenamesigmaV0=efffilename;}


  //PID Type
  void SetPIDTypes(PIDTypes PIDmethod) { fPIDType = PIDmethod; }
  PIDTypes GetPIDTypes() {return fPIDType; }
  //NSigma cut
  //set cut on beyesian probability
  void SetBayesCut(Double_t cut){fBayesCut=cut;}
  void SetdiffPtRanges_PIDcut(Double_t Pt1, Double_t Pt2, Double_t Pt3)//set Pt ranges for diff. pid cut values < fPtTOFPIDMax
  {
    fPt1=Pt1;
    fPt2=Pt2;
    fPt3=Pt3;
  }
  void SetdiffPIDcutvalues(Bool_t diffPIDcutvalues,Double_t PIDCutval1, Double_t PIDCutval2, Double_t PIDCutval3,Double_t PIDCutval4){
    fdiffPIDcutvalues=diffPIDcutvalues;
    fPIDCutval1=PIDCutval1;
    fPIDCutval2=PIDCutval2;
    fPIDCutval3=PIDCutval3;
    fPIDCutval4=PIDCutval4;
  }
 void  SetRandomizeReactionPlane(Bool_t RandomizeReactionPlane){fRandomizeReactionPlane=RandomizeReactionPlane;}
 
   //****************************************************************************************EP related part
  void OpenInfoCalbration(Int_t run);
  void SetTPCclusterN(Int_t ncl){fNcluster=ncl;};
   //****************************************************************************************EP related part
//--------------------------------------------------------------------------//
//v0 daughters

  void SetV0TrigCorr(Bool_t V0TrigCorr){
    fV0TrigCorr=V0TrigCorr;
  }
  void SetfillofflineV0(Bool_t fillofflineV0){
    ffillofflineV0=fillofflineV0;
  }
void SetUsev0DaughterPID(Bool_t Usev0DaughterPID){fUsev0DaughterPID=Usev0DaughterPID;}

 void SetCutsForV0AndDaughters(Double_t MinPtDaughter,Double_t MaxPtDaughter ,Double_t DCAtoPrimVtx, Double_t MaxDCADaughter,Double_t MinCPA,Double_t MaxBoundary,Double_t DaughNClsTPC,Float_t FracSharedTPCcls,Bool_t CutDaughterPtV0)
{
  //fEtaLimitDaughter=EtaLimit;//0.8
fMinPtDaughter=MinPtDaughter;//1.0 GeV/c for our AliHelper
fMaxPtDaughter=MaxPtDaughter;//4.0 GeV/c
fDCAToPrimVtx=DCAtoPrimVtx;//0.1 cm
fMaxDCADaughter=MaxDCADaughter;//1.0 cm
fMinCPA=MinCPA;//0.998
lMax=MaxBoundary;//100 cm
fDaugNClsTPC=DaughNClsTPC;//70
fFracTPCcls=FracSharedTPCcls;//0.4
fCutDaughterPtV0=CutDaughterPtV0;//switch to cut on the daughter of the V0 particles to constrain them within a Pt range where ttrack by track PID can be applied;kFALSE by defaul
}

 void SetCtauCut3D(Bool_t CtauCut3D) {
   fCtauCut3D=CtauCut3D;
 }

 private:                                                                                      
 //histograms
    TList *fOutput;        //! Output list
    TList *fOutputList;        //! Output list
    TList *fList;              //! List for output objects


    TString    fCentralityMethod;     // Method to determine centrality
    Bool_t fPPVsMult;//switch to ON quantile information for pp 7 TeV case
    Bool_t fPileUp_zvtx_INEL_evsel;
    TString    fSampleType;     // pp,p-Pb,Pb-Pb
    Bool_t fRequestEventPlane; //only for PbPb
    Bool_t fRequestEventPlanemixing; //only for PbPb
    Int_t    fnTracksVertex;        // QA tracks pointing to principal vertex
    AliAODVertex* trkVtx;//!
    Float_t zvtx;
    Int_t    fFilterBit;         // track selection cuts
     UInt_t         fTrackStatus;       // if non-0, the bits set in this variable are required for each track
    Double_t       fSharedClusterCut;  // cut on shared clusters (only for AOD, give the actual cut value)
    Double_t fSharedTPCmapCut;//cut on TPC shared map(set any non negative value to implement this cut automatically, no meaning of the value itself)
    Double_t fSharedfraction_Pair_cut;//cut on pairs at the correlation level to check whether the correlating pair has large shared clusters(set fraction percentage to be set as cut off)
    Int_t fVertextype;
    Int_t skipParticlesAbove;
    Double_t fzvtxcut;
    Double_t fVxMax_MC;
    Double_t fVyMax_MC;
    Double_t fVzMax_MC;
    Bool_t ffilltrigassoUNID;
    Bool_t ffilltrigUNIDassoID;
    Bool_t ffilltrigIDassoUNID;
    Bool_t ffilltrigIDassoID;
    Bool_t ffilltrigIDassoIDMCTRUTH;
    Int_t fMaxNofMixingTracks;
    Bool_t fPtOrderMCTruth;
    Bool_t fPtOrderDataReco;
    Bool_t fWeightPerEvent;
    Bool_t fTriggerSpeciesSelection;
    Bool_t fAssociatedSpeciesSelection;
    Bool_t fRandomizeReactionPlane;
    Int_t fTriggerSpecies;
    Int_t fAssociatedSpecies;
    TString fCustomBinning;//for setting customized binning
    TString fBinningString;//final binning string
    Bool_t fSelectHighestPtTrig;
    Bool_t fcontainPIDtrig;
    Bool_t fcontainPIDasso;
    Int_t SetChargeAxis;
     Bool_t frejectPileUp;
     Bool_t fCheckFirstEventInChunk;
    Float_t fminPt;
    Float_t fmaxPt;
    Float_t fmineta;
    Float_t fmaxeta;
    Bool_t fselectprimaryTruth;
    Bool_t fonlyprimarydatareco;
    Double_t fdcacutvalue;
    Bool_t ffillhistQAReco;
    Bool_t ffillhistQATruth;
    Int_t kTrackVariablesPair ;
    Double_t fminPtTrig;
    Double_t fmaxPtTrig;
    Double_t fminPtComboeff;
    Double_t fmaxPtComboeff;
    Double_t fminPtAsso;
    Double_t fmaxPtAsso;
    Double_t fmincentmult;
    Double_t fmaxcentmult;
    TH1F *fPriHistShare;//!
    TH1F *fhistcentrality;//!
    TH1F *fhistImpactParm;//!
    TH2F *fhistImpactParmvsMult;//!
    TH2F *fNchNpartCorr;//!
    TH1F *fEventCounter; //!
    TH2F *fEtaSpectrasso;//!
    TH2F *fphiSpectraasso;//!
    TH1F *MCtruthpt;//! 
    TH1F *MCtrutheta;//! 
    TH1F *MCtruthphi;//!
    TH1F *MCtruthpionpt;//!
    TH1F *MCtruthpioneta;//!
    TH1F *MCtruthpionphi;//!
    TH1F *MCtruthkaonpt;//!
    TH1F *MCtruthkaoneta;//!
    TH1F *MCtruthkaonphi;//!
    TH1F *MCtruthprotonpt;//!
    TH1F *MCtruthprotoneta;//!
    TH1F *MCtruthprotonphi;//! 
    TH2F *fPioncont;//!
    TH2F *fKaoncont;//!
    TH2F *fProtoncont;//!
    TH2F *fUNIDcont;//!
    TH2F *fEventno;//!
    TH2F *fEventnobaryon;//!
    TH2F *fEventnomeson;//!
    TH2F *fhistJetTrigestimate;//!
    TH3F* fTwoTrackDistancePtdip;//!
    TH3F* fTwoTrackDistancePtdipmix;//!
    TH3F* fTwoTrackDistancePt[2];    //! control histograms for two-track efficiency study: dphi*_min vs deta (0 = before cut, 1 = after cut)
    TH3F* fTwoTrackDistancePtmix[2];    //! control histograms for two-track efficiency study: dphi*_min vs deta (0 = before cut, 1 = after cut)

    TH2D* fCentralityCorrelation;  //! centrality vs Tracks multiplicity
    TH2D* fCentralityCorrelationMC;  //! centrality vs Tracks multiplicity in MC

 //VZERO calibration
  TH1F *fHistVZEROAGainEqualizationMap;//VZERO calibration map
  TH1F *fHistVZEROCGainEqualizationMap;//VZERO calibration map
  TH2F *fHistVZEROChannelGainEqualizationMap; //VZERO calibration map
  TH1* fCentralityWeights;		     // for centrality flattening

    TH2F *fHistCentStats; //!centrality stats
    TH2F *fHistRefmult;//!
    TH2F *fHistEQVZEROvsTPCmultiplicity;//!
    TH2F *fHistEQVZEROAvsTPCmultiplicity;//!
    TH2F *fHistEQVZEROCvsTPCmultiplicity;//!
    TH2F *fHistVZEROCvsEQVZEROCmultiplicity;//!
    TH2F *fHistVZEROAvsEQVZEROAmultiplicity;//!
    TH2F *fHistVZEROCvsVZEROAmultiplicity;//!
    TH2F *fHistEQVZEROCvsEQVZEROAmultiplicity;//!
    TH2F *fHistVZEROSignal;//!
    TH2F *fHistEventPlaneTruth;//!
    TH2D *fHistPsiMinusPhi;//! psi - phi QA histogram
    TH3F *fEventPlanePID;//!
   //****************************************************************************************EP related part

    Float_t evplaneMC,fgPsi2v0a,fgPsi2v0c,fgPsi2tpc; // current Psi2
   Float_t fgPsi3v0a,fgPsi3v0c,fgPsi3tpc; // current Psi3
   Float_t fgPsi2v0aMC,fgPsi2v0cMC,fgPsi2tpcMC; // current Psi2
   Float_t fgPsi3v0aMC,fgPsi3v0cMC,fgPsi3tpcMC,gReactionPlane; // current Psi3
  Bool_t fV2; // switch to set the harmonics
  Bool_t fV3; // switch to set the harmonics
  Bool_t fIsAfter2011; // switch for 2011 and later runs

  //    Int_t nCentrBin = 9;          //  cenrality bins


  
  //
  // Cuts and options
  //
 //centrality binning for lambda , kshort a0,a1 factor//*********************************hardcoded values in array, be careful(now in custom binning)

     Int_t kNCent_ds;
     // Double_t *kBinCent_ds;
     
     /*    Double_t *kK0s_a0;
     Double_t *kK0s_a1;
     Double_t *kK0s_a2;
     Double_t *kLambda_a0;
     Double_t *kLambda_a1;
     Double_t *kLambda_a2;
     Double_t *kAntiLambda_a0;
     Double_t *kAntiLambda_a1;
     Double_t *kAntiLambda_a2;*/


  Int_t fRun;                       // current run checked to load VZERO calibrations

  Int_t fNcluster;           // Numer of TPC cluster required

  TString fEPdet; //Set the name of the event plane to be used to reconstruct the event plane
    
  // Output objects
  TProfile *fMultV0;                //! object containing VZERO calibration information
  Float_t fV0Cpol;          //! loaded by OADB
  Float_t fV0Apol;          //! loaded by OADB
  Float_t fMeanQ[9][2][2];           // and recentering
  Float_t fWidthQ[9][2][2];          // ...
  Float_t fMeanQv3[9][2][2];         // also for v3
  Float_t fWidthQv3[9][2][2];        // ...

  TProfile *fHResTPCv0A2;   //! TProfile for subevent resolution (output)
  TProfile *fHResTPCv0C2;   //! TProfile for subevent resolution (output)
  TProfile *fHResv0Cv0A2;   //! TProfile for subevent resolution (output)
  TProfile *fHResTPCv0A3;    //! also for v3
  TProfile *fHResTPCv0C3;   //! also for v3
  TProfile *fHResv0Cv0A3;   //! also for v3

 TProfile *fHResMA2;   //! TProfile for subevent resolution (output)
  TProfile *fHResMC2;   //! TProfile for subevent resolution (output)
  TProfile *fHResAC2;   //! TProfile for subevent resolution (output)
  TProfile *fHResMA3;    //! also for v3
  TProfile *fHResMC3;   //! also for v3
  TProfile *fHResAC3;   //! also for v3

  TH2F *fPhiRPTPC;          //! EP distribution vs. centrality (v2)
  TH2F *fPhiRPTPCv3;          //! EP distribution vs. centrality (v2)
  TH2F *fPhiRPv0A;          //! EP distribution vs. centrality (v2)
  TH2F *fPhiRPv0C;          //! EP distribution vs. centrality (v2)
  TH2F *fPhiRPv0Av3;      //! EP distribution vs. centrality (v3)
  TH2F *fPhiRPv0Cv3;      //! EP distribution vs. centrality (v3)
    //****************************************************************************************EP related part

    TH2F* fControlConvResoncances; //! control histograms for cuts on conversions and resonances

    TH2F *fHistoTPCdEdx;//!
    TH2F *fHistoTOFbeta;//!
    TH3F *fTPCTOFPion3d;//!
    TH3F *fTPCTOFKaon3d;//!
    TH3F *fTPCTOFProton3d;//!
    TH1F *fPionPt;//!
    TH1F *fPionEta;//!
    TH1F *fPionPhi;//!
    TH1F *fKaonPt;//!
    TH1F *fKaonEta;//!
    TH1F *fKaonPhi;//!
    TH1F *fProtonPt;//!
    TH1F *fProtonEta;//!
    TH1F *fProtonPhi;//!

    TH2F *kShortSigmahisto;//!
    TH2F *LambdaSigmahisto;//!
    
  TH2D *fHistdEdxVsPTPCbeforePIDelectron; //!
  TH2D *fHistNSigmaTPCvsPtbeforePIDelectron; //!
  TH2D *fHistdEdxVsPTPCafterPIDelectron; //!
  TH2D *fHistNSigmaTPCvsPtafterPIDelectron; //!
    // TH3F *fHistocentNSigmaTPC;//! nsigma TPC
    // TH3F *fHistocentNSigmaTOF;//! nsigma TOF 
 
    AliTHn *fCorrelatonTruthPrimary;//!
    AliTHn *fCorrelatonTruthPrimarymix;//!
    AliTHn *fTHnCorrUNID;//!
    AliTHn *fTHnCorrUNIDmix;//!
    AliTHn *fTHnCorrID;//!
    AliTHn *fTHnCorrIDmix;//!
    AliTHn *fTHnCorrIDUNID;//!
    AliTHn *fTHnCorrIDUNIDmix;//!
    AliTHn *fTHnTrigcount;//!
    AliTHn *fTHnTrigcountMCTruthPrim;//!
    AliTHn* fTrackHistEfficiency[6]; //! container for tracking efficiency and contamination (all particles filled including leading one): axes: eta, pT, particle species:::::::::0 pion, 1 kaon,2 proton,3 mesons,4 kaons+protons,5 all

    
    TH1F *fHistQA[16]; //!                  
     
   
    THnSparse *effcorection[6];//!

   
    
    // THnF *effmap[6];  

    Int_t ClassifyTrack(AliAODTrack* track,AliAODVertex* vertex,Float_t magfield,Bool_t fill);
  Double_t* GetBinning(const char* configuration, const char* tag, Int_t& nBins);


  void Fillcorrelation(Float_t ReactionPlane,TObjArray *trackstrig,TObjArray *tracksasso,Double_t cent,Float_t vtx,Float_t weight,Bool_t firstTime,Float_t bSign,Bool_t fPtOrder,Bool_t twoTrackEfficiencyCut,Bool_t mixcase,TString fillup);//mixcase=kTRUE in case of mixing; 
 Bool_t CalculateSharedFraction(const TBits *triggerPadMap,const TBits *assocPadMap,const TBits *triggerShareMap,const TBits *assocShareMap);
 Float_t GetTrackbyTrackeffvalue(AliAODTrack* track,Double_t cent,Float_t evzvtx, Int_t parpid);
 Float_t GetV0_MeanSigma_CentPt(Double_t cent, Float_t V0Pt, Int_t parpid);
 
 //Fill PID and Event planes
 void FillPIDEventPlane(Double_t centrality,Int_t par,Float_t trigphi,Float_t fReactionPlane);
 //V0-h correlation related functions
 // Int_t GetCentBin(Double_t cent);
 TObjArray* GetV0Particles(AliVEvent* event,Double_t cent);
 Bool_t CheckStatusv0Daughter(AliAODTrack *t1 ,AliAODTrack *t2);
 Float_t  GetFractionTPCSharedCls( AliAODTrack *track);
 Bool_t  CheckStatusv0(AliAODv0 *v1);
 Bool_t IsTrackFromV0(AliAODEvent* fAOD,AliAODTrack* track);





//Mixing functions
 // void DefineEventPool();
  AliEventPoolManager    *fPoolMgr;//! 
  TClonesArray          *fArrayMC;//!
  TString          fAnalysisType;          // "MCAOD", "MC", "AOD"
  TString fefffilename;
  TString ffilenamesigmaV0;

    //PID part histograms

  //PID functions
    Bool_t HasTPCPID(AliAODTrack *track) const; // has TPC PID
    Bool_t HasTOFPID(AliAODTrack *track) const; // has TOF PID
    Double_t GetBeta(AliAODTrack *track);
    void CalculateNSigmas(AliAODTrack *track, Bool_t FIllQAHistos);
    Int_t FindMinNSigma(AliAODTrack *track, Bool_t FIllQAHistos);
    Bool_t* GetDoubleCounting(AliAODTrack * trk, Bool_t FIllQAHistos);
    Int_t GetParticle(AliAODTrack * trk, Bool_t FIllQAHistos);

    Bool_t TPCCutMIGeo(AliAODTrack* track,AliAODEvent* evt);//TPCSectoredge cut for TPC only PID at high Pt
 
     TH2F* GetHistogram2D(const char * name);//!return histogram "name" from fOutputList

     Bool_t ftwoTrackEfficiencyCutDataReco; 
    Float_t fTwoTrackCutMinRadius;
    Float_t fTwoTrackCutMaxRadius;	   
   Float_t twoTrackEfficiencyCutValue;
  //Pid objects
  AliPIDResponse *fPID; //! PID
  AliPIDCombined   *fPIDCombined;     //! PIDCombined

  //set PID Combined
  void SetPIDCombined(AliPIDCombined *obj){fPIDCombined=obj;}
  AliPIDCombined *GetPIDCombined(){return fPIDCombined;}
 
  Double_t GetBayesCut(){return fBayesCut;}
  Int_t GetIDBayes(AliAODTrack *trk, Bool_t FIllQAHistos);//calculate the PID according to bayesian PID
  UInt_t CalcPIDCombined(AliAODTrack *track,Int_t detMask, Double_t* prob) const;
  Bool_t* GetAllCompatibleIdentitiesNSigma(AliAODTrack * trk, Bool_t FIllQAHistos);//All the identities are true


  Int_t eventno;
  Float_t fPtTOFPIDmin; //lower pt bound for the TOCTOF combined circular pid
  Float_t fPtTOFPIDmax; //uper pt bound for the TOCTOF combined circular pid
  Bool_t fRequestTOFPID;//if true returns kSpUndefined if the TOF signal is missing
  PIDTypes fPIDType; // PID type  Double_t fNSigmaPID; // number of sigma for PID cut
  Bool_t fFIllPIDQAHistos; //Switch for filling the nSigma histos
  Double_t fNSigmaPID; // number of sigma for PID cut
  Double_t fBayesCut; // Cut on Bayesian probability
 Bool_t fdiffPIDcutvalues;
 Double_t fPt1;//set Pt ranges for diff. pid cut values < fPtTOFPIDMax
 Double_t fPt2;
 Double_t fPt3;
 Double_t fPIDCutval1;
 Double_t fPIDCutval2;
 Double_t fPIDCutval3;
 Double_t fPIDCutval4;

  Float_t fHighPtKaonNSigmaPID;// number of sigma for PID cut for Kaons above fHighPtKaonSigma(-1 default, no cut applied)
  Float_t fHighPtKaonSigma;//lower pt bound for the fHighPtKaonNSigmaPID to be set >0(i.e. to make it applicable)
  Bool_t fUseExclusiveNSigma;//if true returns the identity only if no double counting(i.e not in the overlap area)
  Bool_t fRemoveTracksT0Fill;//if true remove tracks for which only StartTime from To-Fill is available (worst resolution)
 Int_t fSelectCharge;           // (un)like sign selection when building correlations: 0: no selection; 1: unlike sign; 2: like sign
    Int_t fTriggerSelectCharge;    // select charge of trigger particle: 1: positive; -1 negative
    Int_t fAssociatedSelectCharge; // select charge of associated particle: 1: positive; -1 negative
    Float_t fTriggerRestrictEta;   // restrict eta range for trigger particle (default: -1 [off])
    Bool_t fEtaOrdering;           // eta ordering, see AliUEHistograms.h for documentation
    Bool_t fCutConversions;        // cut on conversions (inv mass)
    Bool_t fCutResonances;         // cut on resonances (inv mass)
    Int_t fRejectResonanceDaughters; // reject all daughters of all resonance candidates (1: test method (cut at m_inv=0.9); 2: k0; 3: lambda)
    Int_t fOnlyOneEtaSide;       // decides that only trigger particle from one eta side are considered (0 = all; -1 = negative, 1 = positive)
    Bool_t fInjectedSignals;	  // check header to skip injected signals in MC
    Bool_t fRemoveWeakDecays;	   // remove secondaries from weak decays from tracks and particles
    Bool_t fRemoveDuplicates;// remove particles with the same label (double reconstruction)
    Bool_t fapplyTrigefficiency;//if kTRUE then eff correction calculation starts
    Bool_t fapplyAssoefficiency;//if kTRUE then eff correction calculation starts
    Bool_t ffillefficiency;  //if kTRUE then THNsparses used for eff. calculation are filled up
    Bool_t fmesoneffrequired;
    Bool_t fkaonprotoneffrequired;
    Bool_t fRunShufflingbalance;

  Bool_t fElectronRejectionTRUTH;//flag to use electron rejection by PDG code in TRUTH
  Bool_t   fElectronRejection;//flag to use electron rejection by Nsigma method in data and MC
  Bool_t   fElectronOnlyRejection;//flag to use electron rejection with exclusive electron PID (no other particle in nsigma range)
  Double_t fElectronRejectionNSigma;//nsigma cut for electron rejection
  Double_t fElectronRejectionMinPt;//minimum pt for electron rejection (default = 0.)
  Double_t fElectronRejectionMaxPt;//maximum pt for electron rejection (default = 1000.)

    
   AliAnalysisUtils*     fAnalysisUtils;      // points to class with common analysis utilities
   AliPPVsMultUtils*  fPPVsMultUtils;

  TFormula*      fDCAXYCut;          // additional pt dependent cut on DCA XY (only for AOD)
  //*****************************************************************************V0 related objects are here
  Bool_t fV0TrigCorr;
  Bool_t ffillofflineV0;
  Bool_t fUsev0DaughterPID;
 Double_t fMinPtDaughter ;// to be decided to make it compatible with AliHelperPID so far we keep it 1GeV/C
  Double_t fMaxPtDaughter; //same statement as above
  Double_t fDCAToPrimVtx ;//put standard cuts
  Double_t fMaxDCADaughter;//put standard cuts
  Double_t fMinCPA; //cosine of pointing angle
  Double_t lMax;
TH3F*  fHistRawPtCentInvK0s;//!
TH3F*  fHistRawPtCentInvLambda;//!
TH3F*  fHistRawPtCentInvAntiLambda;//!
TH3F*  fHistFinalPtCentInvK0s;//!
TH3F*  fHistFinalPtCentInvLambda;//!
TH3F*  fHistFinalPtCentInvAntiLambda;//!
 Bool_t fCtauCut3D;
  Double_t NCtau;
  Double_t fCutctauK0s; //ctau cut for kShort
  Double_t fCutctauLambda;
  Double_t fCutctauAntiLambda;
  Double_t fRapCutK0s;
  Double_t fRapCutLambda; 
Int_t fDaugNClsTPC;
Float_t fFracTPCcls;
 Bool_t fCutDaughterPtV0;

 Bool_t TPCSectoredgecut;
 Int_t fNclsusedfordEdXdtr;
 Double_t fPiondeltacutmin;
 Double_t fPiondeltacutmax;  
 Double_t fProtondeltacutmin;
 Double_t fProtondeltacutmax;
 Double_t deltapion_val;//global deltapion variable;no need toset any value for it


  Float_t fnsigmas[NSpecies][NSigmaPIDType+1]; //nsigma values
  Bool_t fHasDoubleCounting[NSpecies];//array with compatible identities

  //Int_t fPIDMethod; // PID method

 //functions
  Bool_t CheckTrack(AliAODTrack * part);
  Float_t PhiRange(Float_t DPhi);
  Float_t GetInvMassSquared(Float_t pt1, Float_t eta1, Float_t phi1, Float_t pt2, Float_t eta2, Float_t phi2, Float_t m0_1, Float_t m0_2);
Float_t GetInvMassSquaredCheap(Float_t pt1, Float_t eta1, Float_t phi1, Float_t pt2, Float_t eta2, Float_t phi2, Float_t m0_1, Float_t m0_2);
  Float_t GetDPhiStar(Float_t phi1, Float_t pt1, Float_t charge1, Float_t phi2, Float_t pt2, Float_t charge2, Float_t radius, Float_t bSign);
  TObjArray* GetShuffledTracks(TObjArray *tracks);
  TObjArray* CloneAndReduceTrackList(TObjArray* tracks);


  void ShiftTracks(TObjArray* tracks, Double_t angle);
  Bool_t AcceptEventCentralityWeight(Double_t centrality);

  //get event plane
  Float_t GetEventPlane(AliVEvent *event,Bool_t truth,Double_t v0Centr);
  Double_t GetAcceptedEventMultiplicity(AliVEvent *aod,Bool_t truth);//returns centrality after event(mainly vertex) selection IsEventAccepted  GetAcceptedEventMultiplicity
  
  //get vzero equalization
  Double_t GetEqualizationFactor(Int_t run, const char* side);
  Double_t GetChannelEqualizationFactor(Int_t run,Int_t channel);
  void SetVZEROCalibrationFile(const char* filename,const char* lhcPeriod);
  void SetCentralityWeights(TH1* hist) { fCentralityWeights = hist; }

  Double_t GetRefMultiOrCentrality(AliVEvent *event, Bool_t truth);
  Double_t GetReferenceMultiplicityVZEROFromAOD(AliVEvent *event);//mainly important for pp 7 TeV

    
    AliTwoParticlePIDCorr(const AliTwoParticlePIDCorr&); // not implemented
    AliTwoParticlePIDCorr& operator=(const AliTwoParticlePIDCorr&); // not implemented
    
    ClassDef(AliTwoParticlePIDCorr, 1); // example of analysis
};
class LRCParticlePID : public TObject {
public:
 LRCParticlePID(Int_t par,Double_t Invmass,Short_t icharge,Float_t pt,Float_t eta, Float_t phi,Float_t effcorrectionval,const TBits *clustermap,const TBits *sharemap)
   :fparticle(par),fInvmass(Invmass),fcharge(icharge),fPt(pt), fEta(eta), fPhi(phi),feffcorrectionval(effcorrectionval),fTPCClusterMap(clustermap),fTPCHitShareMap(sharemap) {}
  virtual ~LRCParticlePID() {}
  
    virtual Float_t Eta()        const { return fEta; }
    virtual Float_t Phi()        const { return fPhi; }
    virtual Float_t Pt() const { return fPt; }
    Int_t getparticle() const {return fparticle;}
    Double_t GetInvMass() const {return fInvmass;}
    virtual Short_t Charge()      const { return fcharge; }
    Float_t geteffcorrectionval() const {return feffcorrectionval;}
    virtual Bool_t IsEqual(const TObject* obj) const { return (obj->GetUniqueID() == GetUniqueID()); }
    virtual void SetPhi(Double_t phiv) { fPhi = phiv; }
    virtual const TBits * GetTPCPadMap() {return fTPCClusterMap; }
    virtual const TBits * GetTPCSharedMap() {return fTPCHitShareMap; }

private:
  LRCParticlePID(const LRCParticlePID&);  // not implemented
   LRCParticlePID& operator=(const LRCParticlePID&);  // not implemented
  
  Int_t fparticle;
  Double_t fInvmass;
  Short_t fcharge;
  Float_t fPt;
  Float_t fEta;
  Float_t fPhi;
  Float_t feffcorrectionval;
   const TBits   *fTPCClusterMap;
   const TBits   *fTPCHitShareMap;
  ClassDef(LRCParticlePID, 1);
} ;

#endif

//(fSampleType=="pp_2_76" || fCentralityMethod.EndsWith("_MANUAL"))
//(fSampleType=="pp_2_76" || fCentralityMethod.EndsWith("_MANUAL") || (fSampleType=="pp_7" && fPPVsMultUtils==kFALSE))
//(fCentralityMethod.EndsWith("_MANUAL"))
/*
fMinPtDaughter
fMaxPtDaughter
fDCAToPrimVtx
fMaxDCADaughter
fMinCPA
lMax
fCentralityMethod == "MC_b"
fCentralityCorrelation
*/
//fV0TrigCorr
//ParticlePID_InvMass

//particlepidtrig
//fRequestEventPlanemixing
/*
AliAnalysisTask *AddTaskPIDResponse(Bool_t isMC=kFALSE, Bool_t autoMCesd=kTRUE,
                                    Bool_t tuneOnData=kFALSE, Int_t recoPass=2,
                                    Bool_t cachePID=kFALSE, TString detResponse="",
                                    Bool_t useTPCEtaCorrection = kTRUE,//Please use default value! Otherwise splines can be off
                                    Bool_t useTPCMultiplicityCorrection = kTRUE,//Please use default value! Otherwise splines can be off
                                    Int_t  recoDataPass = -1)*/
