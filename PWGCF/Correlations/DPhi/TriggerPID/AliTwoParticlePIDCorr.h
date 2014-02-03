#ifndef ALITWOPARTICLEPIDCORR_H
#define ALITWOPARTICLEPIDCORR_H

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
class AliAODEvent;
class AliAODTrack;
class AliAODVertex;
class AliEventPoolManager;
class TFormula;
//class AliAnalysisUtils;
class LRCParticlePID;
class AliVParticle;


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


#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

namespace AliPIDNameSpace {
  
  enum PIDType
  {
    NSigmaTPC = 0,
    NSigmaTOF,
    NSigmaTPCTOF,  // squared sum
    NSigmaPIDType=NSigmaTPCTOF
  };
    
  enum AliDetectorType
  {
    TPC = 0,
    TOF,
    NDetectors
  };
  
  
  enum AliParticleSpecies
  {
    SpPion = 0,
    SpKaon,
    SpProton,
    unidentified,
    NSpecies=unidentified,
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
  void		 SetSharedClusterCut(Double_t value) { fSharedClusterCut = value; }

    void SetVertextype(Int_t Vertextype){fVertextype=Vertextype;}
    void SetZvtxcut(Double_t zvtxcut) {fzvtxcut=zvtxcut;}
    void SetCustomBinning(TString receivedCustomBinning) { fCustomBinning = receivedCustomBinning; }
    void SetAsymmetricBin(THnSparse *h,Int_t axisno,Double_t *arraybin,Int_t arraybinsize,TString axisTitle); 
    void SetMaxNofMixingTracks(Int_t MaxNofMixingTracks) {fMaxNofMixingTracks=MaxNofMixingTracks;}
  void SetCentralityEstimator(TString CentralityMethod) { fCentralityMethod = CentralityMethod;}
  void SetSampleType(TString SampleType) {fSampleType=SampleType;}
  void SetAnalysisType(TString AnalysisType){fAnalysisType=AnalysisType;}
  void SetFilterBit(Int_t FilterBit) {fFilterBit=FilterBit;}
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
  // void SetRejectPileUp(Bool_t rejectPileUp) {frejectPileUp=rejectPileUp;}
  void SetKinematicCuts(Float_t minPt, Float_t maxPt,Float_t mineta,Float_t maxeta)
  {
    fminPt=minPt;
    fmaxPt=maxPt;
    fmineta=mineta;
    fmaxeta=maxeta;
  }
  void SetAsymmetricnSigmaCut( Float_t minprotonsigmacut,Float_t maxprotonsigmacut,Float_t minpionsigmacut,Float_t maxpionsigmacut)
  {
    fminprotonsigmacut=minprotonsigmacut;
    fmaxprotonsigmacut=maxprotonsigmacut;
    fminpionsigmacut=minpionsigmacut;
    fmaxpionsigmacut=maxpionsigmacut;
  }
  void SetDcaCut(Bool_t dcacut,Double_t dcacutvalue)
  {
    fdcacut=dcacut;
    fdcacutvalue=dcacutvalue;
  }
  void SetfillHistQA(Bool_t fillhistQAReco,Bool_t fillhistQATruth)
  {
    ffillhistQAReco=fillhistQAReco;
    ffillhistQATruth=fillhistQATruth;
  }
  void Setselectprimarydatareco(Bool_t onlyprimarydatareco) {fonlyprimarydatareco=onlyprimarydatareco;}
  void SetselectprimaryTruth(Bool_t selectprimaryTruth) {fselectprimaryTruth=selectprimaryTruth;}
  void SetCombinedNSigmaCut(Double_t NSigmaPID) {fNSigmaPID=NSigmaPID;}
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

 void SetEffcorectionfilePathName(TString efffilename) {fefffilename=efffilename;}
 
 private:
 //histograms
    TList *fOutput;        //! Output list
    TString    fCentralityMethod;     // Method to determine centrality
    TString    fSampleType;     // pp,p-Pb,Pb-Pb
    Int_t    fnTracksVertex;        // QA tracks pointing to principal vertex
    AliAODVertex* trkVtx;//!
    Float_t zvtx;
    Int_t    fFilterBit;         // track selection cuts
    Double_t       fSharedClusterCut;  // cut on shared clusters (only for AOD)
    Int_t fVertextype;
    Double_t fzvtxcut;
    Bool_t ffilltrigassoUNID;
    Bool_t ffilltrigUNIDassoID;
    Bool_t ffilltrigIDassoUNID;
    Bool_t ffilltrigIDassoID;
    Bool_t ffilltrigIDassoIDMCTRUTH;
    Int_t fMaxNofMixingTracks;
    Bool_t fPtOrderMCTruth;
    Bool_t fTriggerSpeciesSelection;
    Bool_t fAssociatedSpeciesSelection;
    Int_t fTriggerSpecies;
    Int_t fAssociatedSpecies;
    TString fCustomBinning;//for setting customized binning
    TString fBinningString;//final binning string
    Bool_t fSelectHighestPtTrig;
    Bool_t fcontainPIDtrig;
    Bool_t fcontainPIDasso;
    // Bool_t frejectPileUp;
    Float_t fminPt;
    Float_t fmaxPt;
    Float_t fmineta;
    Float_t fmaxeta;
    Float_t fminprotonsigmacut;
    Float_t fmaxprotonsigmacut;
    Float_t fminpionsigmacut;
    Float_t fmaxpionsigmacut;
    Bool_t fselectprimaryTruth;
    Bool_t fonlyprimarydatareco;
    Bool_t fdcacut;
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
    TH1F *fhistcentrality;//!
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

    TH2D* fCentralityCorrelation;  //! centrality vs multiplicity

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
    // TH3F *fHistocentNSigmaTPC;//! nsigma TPC
    // TH3F *fHistocentNSigmaTOF;//! nsigma TOF 
    
    THnSparse *fCorrelatonTruthPrimary;//!
    THnSparse *fCorrelatonTruthPrimarymix;//!
    THnSparse *fTHnCorrUNID;//!
    THnSparse *fTHnCorrUNIDmix;//!
    THnSparse *fTHnCorrID;//!
    THnSparse *fTHnCorrIDmix;//!
    THnSparse *fTHnCorrIDUNID;//!
    THnSparse *fTHnCorrIDUNIDmix;//!
    THnSparse *fTHnTrigcount;//!
    THnSparse *fTHnTrigcountMCTruthPrim;//!

    
    TH1F *fHistQA[16]; //!
     
    THnSparse *fTHnrecomatchedallPid[6];//!                 //0 pion, 1 kaon,2 proton,3 mesons,4 kaons+protons,5 all
    THnSparse *fTHngenprimPidTruth[6];//!
    THnSparse *effcorection[6];//!
    // THnF *effmap[6];  
    //TH2F* fControlConvResoncances; //! control histograms for cuts on conversions and resonances

    Int_t ClassifyTrack(AliAODTrack* track,AliAODVertex* vertex,Float_t magfield);
  Double_t* GetBinning(const char* configuration, const char* tag, Int_t& nBins);


  void Fillcorrelation(TObjArray *trackstrig,TObjArray *tracksasso,Double_t cent,Float_t vtx,Float_t weight,Float_t bSign,Bool_t fPtOrder,Bool_t twoTrackEfficiencyCut,Bool_t mixcase,TString fillup);//mixcase=kTRUE in case of mixing; 
 Float_t GetTrackbyTrackeffvalue(AliAODTrack* track,Double_t cent,Float_t evzvtx, Int_t parpid);

//Mixing functions
  void DefineEventPool();
  // AliAnalysisUtils *fUtils;
  AliEventPoolManager    *fPoolMgr;//! 
  TClonesArray          *fArrayMC;//!
  TString          fAnalysisType;          // "MC", "ESD", "AOD"
  TString fefffilename;

    //PID part histograms

  //PID functions
    Bool_t HasTPCPID(AliAODTrack *track) const; // has TPC PID
    Bool_t HasTOFPID(AliAODTrack *track) const; // has TOF PID
    Float_t GetBeta(AliAODTrack *track);
    void CalculateNSigmas(AliAODTrack *track);
    Int_t FindMinNSigma(AliAODTrack *track);
    Bool_t* GetDoubleCounting(AliAODTrack * trk);
    Int_t GetParticle(AliAODTrack * trk);  
   
   
     	   
   Float_t twoTrackEfficiencyCutValue;
  //Pid objects
  AliPIDResponse *fPID; //! PID
  Int_t eventno;
  Float_t fPtTOFPIDmin; //lower pt bound for the TOCTOF combined circular pid
  Float_t fPtTOFPIDmax; //uper pt bound for the TOCTOF combined circular pid
  Bool_t fRequestTOFPID;//if true returns kSpUndefined if the TOF signal is missing
  PIDType fPIDType; // PID type  Double_t fNSigmaPID; // number of sigma for PID cut
  Double_t fNSigmaPID; // number of sigma for PID cut
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
    Bool_t ffillefficiency;//if kTRUE then THNsparses used for eff. calculation are filled up
    Bool_t fmesoneffrequired;
    Bool_t fkaonprotoneffrequired;
    //  AliAnalysisUtils*     fAnalysisUtils;      // points to class with common analysis utilities
  TFormula*      fDCAXYCut;          // additional pt dependent cut on DCA XY (only for AOD)


  Float_t fnsigmas[NSpecies][NSigmaPIDType+1]; //nsigma values
  Bool_t fHasDoubleCounting[NSpecies];//array with compatible identities

  //Int_t fPIDMethod; // PID method

 //functions
  Float_t PhiRange(Float_t DPhi);
  Float_t GetInvMassSquared(Float_t pt1, Float_t eta1, Float_t phi1, Float_t pt2, Float_t eta2, Float_t phi2, Float_t m0_1, Float_t m0_2);
Float_t GetInvMassSquaredCheap(Float_t pt1, Float_t eta1, Float_t phi1, Float_t pt2, Float_t eta2, Float_t phi2, Float_t m0_1, Float_t m0_2);
  Float_t GetDPhiStar(Float_t phi1, Float_t pt1, Float_t charge1, Float_t phi2, Float_t pt2, Float_t charge2, Float_t radius, Float_t bSign);
  TObjArray* CloneAndReduceTrackList(TObjArray* tracks);
	   
    
    AliTwoParticlePIDCorr(const AliTwoParticlePIDCorr&); // not implemented
    AliTwoParticlePIDCorr& operator=(const AliTwoParticlePIDCorr&); // not implemented
    
    ClassDef(AliTwoParticlePIDCorr, 1); // example of analysis
};
class LRCParticlePID : public TObject {
public:
 LRCParticlePID(Int_t par,Short_t icharge,Float_t pt,Float_t eta, Float_t phi,Float_t effcorrectionval)
   :fparticle(par),fcharge(icharge),fPt(pt), fEta(eta), fPhi(phi),feffcorrectionval(effcorrectionval)  {}
  virtual ~LRCParticlePID() {}

  
    virtual Float_t Eta()        const { return fEta; }
    virtual Float_t Phi()        const { return fPhi; }
    virtual Float_t Pt() const { return fPt; }
    Int_t getparticle() const {return fparticle;}
    virtual Short_t Charge()      const { return fcharge; }
    Float_t geteffcorrectionval() const {return feffcorrectionval;}
    virtual Bool_t IsEqual(const TObject* obj) const { return (obj->GetUniqueID() == GetUniqueID()); }


private:
  LRCParticlePID(const LRCParticlePID&);  // not implemented
   LRCParticlePID& operator=(const LRCParticlePID&);  // not implemented
  
  Int_t fparticle;
  Short_t fcharge;
  Float_t fPt;
  Float_t fEta;
  Float_t fPhi;
  Float_t feffcorrectionval;
  ClassDef(LRCParticlePID, 1);
} ;

#endif

