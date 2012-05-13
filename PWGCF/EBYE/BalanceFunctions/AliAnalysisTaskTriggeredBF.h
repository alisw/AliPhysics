#ifndef ALIANALYSISTASKTRIGGEREDBF_CXX
#define ALIANALYSISTASKTRIGGEREDBF_CXX

// Analysis task for the TriggeredBF code
// Authors: Panos Cristakoglou@cern.ch, m.weber@cern.ch

#include "AliLog.h"
#include "AliAnalysisTaskSE.h"
#include "AliBalanceTriggered.h"

class TList;
class TH1F;
class TH2F;

class AliBalanceTriggered;
class AliEventPoolManager;


class AliAnalysisTaskTriggeredBF : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskTriggeredBF(const char *name = "AliAnalysisTaskTriggeredBF");
  virtual ~AliAnalysisTaskTriggeredBF(); 
  
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   UserExecMix(Option_t*);
  virtual void   FinishTaskOutput();
  virtual void   Terminate(Option_t *);

  void SetAnalysisObject(AliBalanceTriggered *const analysis) {
    fBalance         = analysis;
    }
  void SetShufflingObject(AliBalanceTriggered *const analysisShuffled) {
    fRunShuffling = kTRUE;
    fShuffledBalance = analysisShuffled;
  }
  void SetMixingObject(AliBalanceTriggered *const analysisMixed) {
    fRunMixing = kTRUE;
    fMixedBalance = analysisMixed;
  }
  void SetMixingTracks(Int_t tracks) { fMixingTracks = tracks; }

 
  void SetVertexDiamond(Double_t vx, Double_t vy, Double_t vz) {
    fVxMax = vx;
    fVyMax = vy;
    fVzMax = vz;
  }

  //==============AOD analysis==============//
  void SetAODtrackCutBit(Int_t bit){
    nAODtrackCutBit = bit;
  }

  void SetKinematicsCutsAOD(Double_t ptmin, Double_t ptmax, Double_t etamin, Double_t etamax){
    fPtMin  = ptmin;
    fPtMax  = ptmax;
    fEtaMin = etamin;
    fEtaMax = etamax;

  }

  void SetExtraDCACutsAOD(Double_t DCAxy, Double_t DCAz){
    fDCAxyCut  = DCAxy;
    fDCAzCut = DCAz;
  }

   void SetExtraTPCCutsAOD(Double_t maxTPCchi2, Int_t minNClustersTPC){
    fTPCchi2Cut      = maxTPCchi2;
    fNClustersTPCCut = minNClustersTPC;
  }

  //Centrality
  void SetCentralityEstimator(const char* centralityEstimator) {fCentralityEstimator = centralityEstimator;}
  const char* GetCentralityEstimator(void)                     {return fCentralityEstimator;}
  void SetCentralityPercentileRange(Double_t min, Double_t max) { 
    fUseCentrality = kTRUE;
    fCentralityPercentileMin=min;
    fCentralityPercentileMax=max;
  }
  void SetImpactParameterRange(Double_t min, Double_t max) { 
    fUseCentrality = kTRUE;
    fImpactParameterMin=min;
    fImpactParameterMax=max;
  }

  //multiplicity
  void SetMultiplicityRange(Int_t min, Int_t max) {
    fUseMultiplicity = kTRUE;
    fNumberOfAcceptedTracksMin = min;
    fNumberOfAcceptedTracksMax = max;}
  
  void UseOfflineTrigger() {fUseOfflineTrigger = kTRUE;}
  

 private:
  Float_t    IsEventAccepted(AliVEvent* event);
  TObjArray* GetAcceptedTracks(AliVEvent* event);
  TObjArray* GetShuffledTracks(TObjArray* tracks);

  AliBalanceTriggered *fBalance; //TriggeredBF object
  Bool_t fRunShuffling;//run shuffling or not
  AliBalanceTriggered *fShuffledBalance; //TriggeredBF object (shuffled)
  Bool_t fRunMixing;//run mixing or not
  Int_t  fMixingTracks;
  AliBalanceTriggered *fMixedBalance; //TriggeredBF object (mixed)
  AliEventPoolManager*     fPoolMgr;         //! event pool manager
    

  TList *fList; //fList object
  TList *fListTriggeredBF; //fList object
  TList *fListTriggeredBFS; //fList object (shuffling)
  TList *fListTriggeredBFM; //fList object (mixing)
  TList *fHistListPIDQA;  //! list of histograms

  TH1F *fHistEventStats; //event stats
  TH2F *fHistCentStats; //centrality stats
  TH1F *fHistTriggerStats; //trigger stats
  TH1F *fHistTrackStats; //Track filter bit stats
  TH1F *fHistVx; //x coordinate of the primary vertex
  TH1F *fHistVy; //y coordinate of the primary vertex
  TH1F *fHistVz; //z coordinate of the primary vertex

  TH2F *fHistClus;//
  TH2F *fHistDCA;//
  TH1F *fHistChi2;//
  TH1F *fHistPt;//
  TH1F *fHistEta;//
  TH1F *fHistPhi;//
  TH1F *fHistPhiBefore;//
  TH1F *fHistPhiAfter;//
  TH2F *fHistV0M;//
  TH2F *fHistRefTracks;//

 
  TString fCentralityEstimator;      //"V0M","TRK","TKL","ZDC","FMD"
  Bool_t fUseCentrality;//use the centrality (PbPb) or not (pp)
  Double_t fCentralityPercentileMin;//centrality percentile min
  Double_t fCentralityPercentileMax;//centrality percentile max
  Double_t fImpactParameterMin;//impact parameter min (used for MC)
  Double_t fImpactParameterMax;//impact parameter max (used for MC)

  Bool_t fUseMultiplicity;//use the multiplicity cuts
  Int_t fNumberOfAcceptedTracksMin;//min. number of number of accepted tracks (used for the multiplicity dependence study - pp)
  Int_t fNumberOfAcceptedTracksMax;//max. number of number of accepted tracks (used for the multiplicity dependence study - pp)
  TH1F *fHistNumberOfAcceptedTracks;//hisot to store the number of accepted tracks

  Bool_t fUseOfflineTrigger;//Usage of the offline trigger selection

  Double_t fVxMax;//vxmax
  Double_t fVyMax;//vymax
  Double_t fVzMax;//vzmax

  Int_t nAODtrackCutBit;//track cut bit from track selection (only used for AODs)

  Double_t fPtMin;//only used for AODs
  Double_t fPtMax;//only used for AODs
  Double_t fEtaMin;//only used for AODs
  Double_t fEtaMax;//only used for AODs

  Double_t fDCAxyCut;//only used for AODs
  Double_t fDCAzCut;//only used for AODs

  Double_t fTPCchi2Cut;//only used for AODs
  Int_t fNClustersTPCCut;//only used for AODs

 
  AliAnalysisTaskTriggeredBF(const AliAnalysisTaskTriggeredBF&); // not implemented
  AliAnalysisTaskTriggeredBF& operator=(const AliAnalysisTaskTriggeredBF&); // not implemented
  
  ClassDef(AliAnalysisTaskTriggeredBF, 1); // example of analysis
};

// class used for io with AliBalance (taken from AliAnalysisTaskPhiCorrelations)
class AliBFBasicParticle : public AliVParticle
{ 
  public:
    AliBFBasicParticle(Float_t eta, Float_t phi, Float_t pt, Short_t charge)
      : fEta(eta), fPhi(phi), fpT(pt), fCharge(charge)
    {
    }
    ~AliBFBasicParticle() {}
    
    // kinematics
    virtual Double_t Px() const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Py() const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Pz() const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Pt() const { return fpT; }
    virtual Double_t P() const { AliFatal("Not implemented"); return 0; }
    virtual Bool_t   PxPyPz(Double_t[3]) const { AliFatal("Not implemented"); return 0; }

    virtual Double_t Xv() const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Yv() const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Zv() const { AliFatal("Not implemented"); return 0; }
    virtual Bool_t   XvYvZv(Double_t[3]) const { AliFatal("Not implemented"); return 0; }

    virtual Double_t OneOverPt()  const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Phi()        const { return fPhi; }
    virtual Double_t Theta()      const { AliFatal("Not implemented"); return 0; }


    virtual Double_t E()          const { AliFatal("Not implemented"); return 0; }
    virtual Double_t M()          const { AliFatal("Not implemented"); return 0; }
    
    virtual Double_t Eta()        const { return fEta; }
    virtual Double_t Y()          const { AliFatal("Not implemented"); return 0; }
    
    virtual Short_t Charge()      const { return fCharge; }
    virtual Short_t SetCharge(Short_t charge)   { fCharge = charge; }
    virtual Int_t   GetLabel()    const { AliFatal("Not implemented"); return 0; }
    // PID
    virtual Int_t   PdgCode()     const { AliFatal("Not implemented"); return 0; }      
    virtual const Double_t *PID() const { AliFatal("Not implemented"); return 0; }
    
  private:
    Float_t fEta;      // eta
    Float_t fPhi;      // phi
    Float_t fpT;       // pT
    Short_t fCharge;   // charge
    
    ClassDef( AliBFBasicParticle, 1); // class which contains only quantities requires for this analysis to reduce memory consumption for event mixing
};

#endif
