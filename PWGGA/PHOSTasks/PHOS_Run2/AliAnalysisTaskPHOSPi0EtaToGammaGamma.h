#ifndef AliAnalysisTaskPHOSPi0EtaToGammaGamma_cxx
#define AliAnalysisTaskPHOSPi0EtaToGammaGamma_cxx

// Author: Daiki Sekihata (Hiroshima University)

class TList;
class AliVEvent;
class AliESDEvent;
class AliAODEvent;
class AliPHOSGeometry;
class TClonesArray;
class AliESDtrackCuts;
class AliMultSelection;
class AliStack;
class AliQnCorrectionsManager;

#include "AliAnalysisTaskSE.h"

#include "AliPHOSTriggerHelper.h"
#include "AliPHOSEventCuts.h"
#include "AliPHOSClusterCuts.h"
#include "AliPHOSJetJetMC.h"

class AliAnalysisTaskPHOSPi0EtaToGammaGamma : public AliAnalysisTaskSE {
  public:

    AliAnalysisTaskPHOSPi0EtaToGammaGamma(const char *name = "Pi0EtaToGammaGamma");
    virtual ~AliAnalysisTaskPHOSPi0EtaToGammaGamma(); 
    void SetNonLinearityStudyFlag(Bool_t flag) {fIsNonLinStudyNeeded = flag;}

    void SetESDtrackCutsForGlobal(AliESDtrackCuts* cuts){fESDtrackCutsGlobal = cuts;}
    void SetESDtrackCutsForGlobalConstrained(AliESDtrackCuts* cuts){fESDtrackCutsGlobalConstrained = cuts;}

    void SetTenderFlag(Bool_t tender) {fUsePHOSTender = tender;}
    void SetMCFlag(Bool_t mc) {fIsMC = mc;}
    void SetCoreEnergyFlag(Bool_t iscore) {fUseCoreEnergy = iscore;}
    void SetBunchSpace(Double_t bs) {fBunchSpace = bs;}
    void SetCollisionSystem(Int_t id) {fCollisionSystem = id;}
    void SetQnVectorTask(Bool_t flag) {fIsFlowTask = flag;}
    void SetJetJetMC(Bool_t flag){ fIsJJMC = flag; }
    void SetMCType(TString type){fMCType = type;}
    void SetTOFCutEfficiencyFunction(TF1 *f1) {fTOFEfficiency = f1;}

    void SetEventCuts(Bool_t isMC){
      fPHOSEventCuts = new AliPHOSEventCuts("PHOSEventCuts");
      fPHOSEventCuts->SetMCFlag(isMC);
      fPHOSEventCuts->SetMaxAbsZvtx(10.);
      fPHOSEventCuts->SetRejectPileup(kTRUE);
      fPHOSEventCuts->SetRejectDAQIncompleteEvent(kTRUE);
    }

    void SetClusterCuts(Bool_t useCoreDisp, Double_t NsigmaCPV, Double_t NsigmaDisp){
      fPHOSClusterCuts = new AliPHOSClusterCuts("PHOSClusterCuts");
      fPHOSClusterCuts->SetUseCoreDispersion(useCoreDisp);
      fPHOSClusterCuts->SetNsigmaCPV(NsigmaCPV);
      fPHOSClusterCuts->SetNsigmaDisp(NsigmaDisp);
    }

    void SetAdditionalPi0PtWeightFunction(TArrayD *centarray, TObjArray *funcarray) {
      Int_t Ncen = centarray->GetSize();
      fCentArrayPi0 = centarray;
      for(Int_t icen=0;icen<Ncen-1;icen++){
        fAdditionalPi0PtWeight[icen] = (TF1*)funcarray->At(icen);
      }
    }

    void SetAdditionalK0SPtWeightFunction(TArrayD *centarray, TObjArray *funcarray) {
      Int_t Ncen = centarray->GetSize();
      fCentArrayK0S = centarray;
      for(Int_t icen=0;icen<Ncen-1;icen++){
        fAdditionalK0SPtWeight[icen] = (TF1*)funcarray->At(icen);
      }
    }//adjust charged K/pi ratio

    void SetCentralityMin(Float_t min) {fCentralityMin = min;}
    void SetCentralityMax(Float_t max) {fCentralityMax = max;}
    void SetDepthNMixed(Int_t Nmix)    {fNMixed        = Nmix;}
    void SetCentralityEstimator(TString estimator) {fEstimator = estimator;}
    void SetPHOSTriggerAnalysis(TString selection){
      fIsPHOSTriggerAnalysis = kTRUE;
      fPHOSTriggerHelper  = new AliPHOSTriggerHelper(selection);
    }

  protected:
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *option);
    virtual void ProcessMC();
    void TrackQA();
    void CellQA();
    void TriggerQA();
    virtual void ClusterQA();
    virtual void FillPhoton();
    virtual void FillMgg();
    virtual void FillMixMgg();
    virtual void EstimatePIDCutEfficiency();
    void EstimateTOFCutEfficiency();
    void EstimateTriggerEfficiency();
    void SelectTriggeredCluster();
    void FillRejectionFactorMB();
    void FillEpRatio();

    virtual void SetMCWeight();//set weight related to M.C. (pT slope of mother pi0/eta/K0S/gamma)

    Bool_t IsFrom(Int_t label, Double_t &TruePt, const Int_t target_pdg);
    Bool_t IsPhoton(Int_t label);

    Double_t R(AliAODMCParticle *p);//in cylindrical system
    Double_t Rho(AliAODMCParticle *p);//in sperical system
    Double_t DeltaPhiIn0Pi(Double_t dphi);//this returns dphi in 0-pi range.

    virtual Int_t FindCommonParent(Int_t iPart, Int_t jPart);

    virtual void DoNonLinearityStudy();

    void FillHistogramTH1(TList *list, const Char_t *name, Double_t x, Double_t w=1., Option_t *opt = "") const ;
    void FillHistogramTH2(TList *list, const Char_t *name, Double_t x, Double_t y, Double_t w=1., Option_t *opt = "") const ;
    void FillHistogramTH3(TList *list, const Char_t *name, Double_t x, Double_t y, Double_t z, Double_t w=1., Option_t *opt = "") const ;
    void FillProfile(TList *list, const Char_t *name, Double_t x, Double_t y) const ;
    void FillSparse(TList *list, const Char_t *name, Double_t *x, Double_t w=1.) const;

    TF1 *GetTOFCutEfficiencyFunction() {return fTOFEfficiency;}

    TF1 *GetAdditionalPi0PtWeightFunction(Float_t centrality){
      if(fCentArrayPi0){
        if(fCollisionSystem==0) centrality -= 0.5;//for pp, in case that centrality is defind by multiplicity 
        Int_t lastBinUpperIndex = fCentArrayPi0->GetSize()-1;
        Int_t index = TMath::BinarySearch<Double_t>( lastBinUpperIndex, fCentArrayPi0->GetArray(), centrality);
        return fAdditionalPi0PtWeight[index];
      }
      else
        return fAdditionalPi0PtWeight[0]; 
    }

    TF1 *GetAdditionalK0SPtWeightFunction(Float_t centrality){
      if(fCentArrayK0S){
        if(fCollisionSystem==0) centrality -= 0.5;//for pp, in case that centrality is defind by multiplicity 
        Int_t lastBinUpperIndex = fCentArrayK0S->GetSize()-1;
        Int_t index = TMath::BinarySearch<Double_t>( lastBinUpperIndex, fCentArrayK0S->GetArray(), centrality);
        return fAdditionalK0SPtWeight[index];
      }
      else 
        return fAdditionalK0SPtWeight[0];
    }

    AliStack *GetMCInfoESD();
    TClonesArray *GetMCInfoAOD();

    AliGenPythiaEventHeader* GetPythiaEventHeader(AliVEvent *event);
    Bool_t ComparePtHardWithJet(AliVEvent *event);
    Bool_t ComparePtHardWithSingleParticle(AliVEvent *event);
    Int_t FindPrimaryMotherESD(Int_t label);
    AliPHOSGeometry *GetPHOSGeometry();

  protected:
    Bool_t fIsMC;
    Bool_t fIsJJMC;//jet jet MC
    Int_t fPtHardBin;//your selected pT hard bin for jet jet MC
    TString fMCType;//MBMC, JJMC
    Bool_t fUsePHOSTender;
    Bool_t fUseCoreEnergy;
    AliPHOSEventCuts *fPHOSEventCuts;
    AliPHOSClusterCuts *fPHOSClusterCuts;
    Double_t fBunchSpace;// in unit of ns.
    Double_t fMinDistBCM;//minimum distance to the closest bad channel
    Int_t fCollisionSystem;//colliions system : pp=0, PbPb=1, pPb (Pbp)=2;
    TF1 *fTOFEfficiency;//TOF cut efficiency as a function of cluster energy;
    AliESDtrackCuts *fESDtrackCutsGlobal;//good global track
    AliESDtrackCuts *fESDtrackCutsGlobalConstrained;//global track but constrained to IP because of SPD dead area
    TF1 *fAdditionalPi0PtWeight[10];//weight function for pT distribution
    TF1 *fAdditionalK0SPtWeight[10];//weight function for pT distribution
    TArrayD *fCentArrayPi0;
    TArrayD *fCentArrayK0S;

    THashList *fOutputContainer;
    AliVEvent *fEvent;
    AliESDEvent *fESDEvent;
    AliAODEvent *fAODEvent;
    AliStack *fMCArrayESD;     //MC particles array in ESD
    TClonesArray *fMCArrayAOD; //MC particles array in AOD
    AliPHOSJetJetMC *fJJMCHandler;
    Int_t fRunNumber;
    AliPHOSGeometry *fPHOSGeo;
    TList *fPHOSEvents[10][12];
    TClonesArray *fPHOSClusterArray;
    TString fEstimator;
    AliMultSelection *fMultSelection;
    Float_t fCentralityMain;
    Float_t fCentralityMin;
    Float_t fCentralityMax;
    Int_t fNMixed;
    Double_t fVertex[3];
    Int_t fZvtx;
    Bool_t fIsFlowTask;
    AliQnCorrectionsManager *fFlowQnVectorMgr;
    TString fTPCEPName[3]; 
    TString fV0EPName[3]; 
    Double_t fEPV0A;
    Double_t fEPV0C;
    Double_t fEPTPC;
    Int_t fEPBin;
    Int_t fNHybridTrack;
    Bool_t fIsNonLinStudyNeeded;
    TF1 *fNonLin[7][7];
    Bool_t fIsPHOSTriggerAnalysis;
    AliPHOSTriggerHelper *fPHOSTriggerHelper;//for real PHOS triggered data analysis
    AliPHOSTriggerHelper *fPHOSTriggerHelperL0; //only for rejection factor in MB
    AliPHOSTriggerHelper *fPHOSTriggerHelperL1H;//only for rejection factor in MB
    AliPHOSTriggerHelper *fPHOSTriggerHelperL1M;//only for rejection factor in MB
    AliPHOSTriggerHelper *fPHOSTriggerHelperL1L;//only for rejection factor in MB
    AliPIDResponse *fPIDResponse;

  private:
    AliAnalysisTaskPHOSPi0EtaToGammaGamma(const AliAnalysisTaskPHOSPi0EtaToGammaGamma&);
    AliAnalysisTaskPHOSPi0EtaToGammaGamma& operator=(const AliAnalysisTaskPHOSPi0EtaToGammaGamma&);

    ClassDef(AliAnalysisTaskPHOSPi0EtaToGammaGamma, 22);
};

#endif
