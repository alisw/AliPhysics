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

#include "TVector2.h"
#include "AliAnalysisTaskSE.h"
#include "AliQnCorrectionsQnVector.h"

#include "AliPHOSTriggerHelper.h"
#include "AliPHOSEventCuts.h"
#include "AliPHOSClusterCuts.h"
#include "AliPHOSJetJetMC.h"

class AliAnalysisTaskPHOSPi0EtaToGammaGamma : public AliAnalysisTaskSE {
  public:

    enum FlowMethod{
      kOFF  = -1,
      kEP   = 0,
      kSP   = 1
    };

    enum QnDetector{
      kNone      = -1,
      kFullTPC   = 0,
      kTPCNegEta = 1,
      kTPCPosEta = 2,
      kFullV0    = 3,
      kV0A       = 4,
      kV0C       = 5
    };

    AliAnalysisTaskPHOSPi0EtaToGammaGamma(const char *name = "Pi0EtaToGammaGamma");
    virtual ~AliAnalysisTaskPHOSPi0EtaToGammaGamma(); 

    void SetESDtrackCutsForGlobal(AliESDtrackCuts* cuts){fESDtrackCutsGlobal = cuts;}
    void SetESDtrackCutsForGlobalConstrained(AliESDtrackCuts* cuts){fESDtrackCutsGlobalConstrained = cuts;}

    void SetTenderFlag(Bool_t tender) {fUsePHOSTender = tender;}
    void SetMCFlag(Bool_t mc) {fIsMC = mc;}
    void SetCoreEnergyFlag(Bool_t iscore) {fUseCoreEnergy = iscore;}
    void SetBunchSpace(Double_t bs) {fBunchSpace = bs;}
    void SetCollisionSystem(Int_t id) {fCollisionSystem = id;}
    void SetQnVectorTask(Bool_t flag) {fIsFlowTask = flag;}
    void SetHarmonics(Int_t harmonics) {fHarmonics = harmonics;}
    void SetJetJetMC(Bool_t flag){ fIsJJMC = flag; }
    void SetMCType(TString type){fMCType = type;}
    void SetTOFCutEfficiencyFunction(TF1 *f1) {fTOFEfficiency = f1;}
    void SetNonLinearityStudy(Bool_t flag) {fIsNonLinStudy = flag;}

    void SetTriggerEfficiency(TF1 *f1) {fTriggerEfficiency = f1;}

    void SetEventCuts(Bool_t isMC){
      fPHOSEventCuts = new AliPHOSEventCuts("PHOSEventCuts");
      fPHOSEventCuts->SetMCFlag(isMC);
      fPHOSEventCuts->SetMaxAbsZvtx(10.);
      fPHOSEventCuts->SetRejectPileup(kTRUE);
      fPHOSEventCuts->SetRejectDAQIncompleteEvent(kTRUE);
    }

    void SetClusterCuts(Bool_t useCoreDisp, Double_t NsigmaCPV, Double_t NsigmaDisp, Double_t distBC){
      fPHOSClusterCuts = new AliPHOSClusterCuts("PHOSClusterCuts");
      fPHOSClusterCuts->SetUseCoreDispersion(useCoreDisp);
      fPHOSClusterCuts->SetNsigmaCPV(NsigmaCPV);
      fPHOSClusterCuts->SetNsigmaDisp(NsigmaDisp);
      fPHOSClusterCuts->SetMinDistanceFromBC(distBC);
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

    void SetAdditionalEtaPtWeightFunction(TArrayD *centarray, TObjArray *funcarray) {
      Int_t Ncen = centarray->GetSize();
      fCentArrayEta = centarray;
      for(Int_t icen=0;icen<Ncen-1;icen++){
        fAdditionalEtaPtWeight[icen] = (TF1*)funcarray->At(icen);
      }
    }

    void SetAdditionalGammaPtWeightFunction(TArrayD *centarray, TObjArray *funcarray) {
      Int_t Ncen = centarray->GetSize();
      fCentArrayGamma = centarray;
      for(Int_t icen=0;icen<Ncen-1;icen++){
        fAdditionalGammaPtWeight[icen] = (TF1*)funcarray->At(icen);
      }
    }

    void SetCentralityMin(Float_t min) {fCentralityMin = min;}
    void SetCentralityMax(Float_t max) {fCentralityMax = max;}
    void SetDepthNMixed(Int_t Nmix)    {fNMixed        = Nmix;}
    void SetCentralityEstimator(TString estimator) {fEstimator = estimator;}
    void SetQNomalization(TString norm) {fQNormalization = norm;}
    void SetFlowMethod(Int_t fm) {fFM = fm;}
    void SetQnDetector(Int_t det){
      fQnDetectorMain = det;

      if(fQnDetectorMain == AliAnalysisTaskPHOSPi0EtaToGammaGamma::kFullTPC){
        fQnDetectorSub1 = AliAnalysisTaskPHOSPi0EtaToGammaGamma::kV0A;
        fQnDetectorSub2 = AliAnalysisTaskPHOSPi0EtaToGammaGamma::kV0C;
      }
      else if(fQnDetectorMain == AliAnalysisTaskPHOSPi0EtaToGammaGamma::kTPCNegEta){
        fQnDetectorSub1 = AliAnalysisTaskPHOSPi0EtaToGammaGamma::kV0A;
        fQnDetectorSub2 = AliAnalysisTaskPHOSPi0EtaToGammaGamma::kV0C;
      }
      else if(fQnDetectorMain == AliAnalysisTaskPHOSPi0EtaToGammaGamma::kTPCPosEta){
        fQnDetectorSub1 = AliAnalysisTaskPHOSPi0EtaToGammaGamma::kV0A;
        fQnDetectorSub2 = AliAnalysisTaskPHOSPi0EtaToGammaGamma::kV0C;
      }
      else if(fQnDetectorMain == AliAnalysisTaskPHOSPi0EtaToGammaGamma::kFullV0){
        fQnDetectorSub1 = AliAnalysisTaskPHOSPi0EtaToGammaGamma::kTPCNegEta;
        fQnDetectorSub2 = AliAnalysisTaskPHOSPi0EtaToGammaGamma::kTPCPosEta;
      }
      else if(fQnDetectorMain == AliAnalysisTaskPHOSPi0EtaToGammaGamma::kV0A){
        fQnDetectorSub1 = AliAnalysisTaskPHOSPi0EtaToGammaGamma::kV0C;
        fQnDetectorSub2 = AliAnalysisTaskPHOSPi0EtaToGammaGamma::kFullTPC;
      }
      else if(fQnDetectorMain == AliAnalysisTaskPHOSPi0EtaToGammaGamma::kV0C){
        fQnDetectorSub1 = AliAnalysisTaskPHOSPi0EtaToGammaGamma::kV0A;
        fQnDetectorSub2 = AliAnalysisTaskPHOSPi0EtaToGammaGamma::kFullTPC;
      }
      else{ 
        AliInfo("Your choice is ignored. QnVector is measured by FullV0-TPCNegEta-TPCPosEta");
        fQnDetectorMain = AliAnalysisTaskPHOSPi0EtaToGammaGamma::kFullV0;
        fQnDetectorSub1 = AliAnalysisTaskPHOSPi0EtaToGammaGamma::kTPCNegEta;
        fQnDetectorSub2 = AliAnalysisTaskPHOSPi0EtaToGammaGamma::kTPCPosEta;
      }
    }

    void SetPHOSTriggerAnalysis(TString selection, Bool_t isMC){
      //obsolete
      fIsPHOSTriggerAnalysis = kTRUE;
      fPHOSTriggerHelper  = new AliPHOSTriggerHelper(selection,isMC);
    }

    void SetPHOSTriggerAnalysis(Int_t L1input, Int_t L0input, Double_t Ethre, Bool_t isMC){
      fIsPHOSTriggerAnalysis = kTRUE;
      fEnergyThreshold = Ethre;
      fPHOSTriggerHelper  = new AliPHOSTriggerHelper(L1input,L0input,isMC);
    }

    void SetTriggerMatchingDeltaR(Double_t DeltaR){
      fPHOSTriggerHelper->SetMatchingDeltaR(DeltaR);
    }
    void SetTriggerMatchingDXZ(Int_t xmin, Int_t zmin, Int_t xmax, Int_t zmax){
      //this is a default setting
      fPHOSTriggerHelper->SetMatchingDistance(xmin,zmin,xmax,zmax);
    }

    void SetForceActiveTRU(Int_t L1input, Int_t L0input, Double_t Ethre, Bool_t isMC){
      //this function should not be called together with SetPHOSTriggerAnalysis
      AliInfo("Force active TRU region!");
      if(fIsPHOSTriggerAnalysis) AliInfo("fIsPHOSTriggerAnalysis = kTRUE. Are you sure you want to use active TRU separately?");

      fForceActiveTRU = kTRUE;
      fIsPHOSTriggerAnalysis = kFALSE;
      fEnergyThreshold = Ethre;
      fPHOSTriggerHelper  = new AliPHOSTriggerHelper(L1input,L0input,isMC);
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
    const AliQnCorrectionsQnVector *GetQnVectorFromList(const TList *qnlist, const char* subdetector, const char *expcorr, const char *altcorr);

    virtual void SetMCWeight();//set weight related to M.C. (pT slope of mother pi0/eta/K0S/gamma)

    Bool_t IsFrom(Int_t label, Double_t &TruePt, const Int_t target_pdg);
    Bool_t IsPhoton(Int_t label);

    Double_t R(AliAODMCParticle *p);//in cylindrical system
    Double_t Rho(AliAODMCParticle *p);//in sperical system
    Double_t DeltaPhiIn0Pi(Double_t dphi);//this returns dphi in 0-pi range.

    virtual Int_t FindCommonParent(Int_t iPart, Int_t jPart);

    void FillHistogramTH1(TList *list, const Char_t *name, Double_t x, Double_t w=1., Option_t *opt = "") const ;
    void FillHistogramTH2(TList *list, const Char_t *name, Double_t x, Double_t y, Double_t w=1., Option_t *opt = "") const ;
    void FillHistogramTH3(TList *list, const Char_t *name, Double_t x, Double_t y, Double_t z, Double_t w=1., Option_t *opt = "") const ;
    void FillProfile(TList *list, const Char_t *name, Double_t x, Double_t y) const ;
    void FillSparse(TList *list, const Char_t *name, Double_t *x, Double_t w=1.) const;

    TF1 *GetTOFCutEfficiencyFunction() {return fTOFEfficiency;}
    TF1 *GetTriggerEfficiencyFunction() {return fTriggerEfficiency;}

    TF1 *GetAdditionalPi0PtWeightFunction(Float_t centrality){
      if(fCentArrayPi0){
        Int_t lastBinUpperIndex = fCentArrayPi0->GetSize()-1;
        Int_t index = TMath::BinarySearch<Double_t>( lastBinUpperIndex, fCentArrayPi0->GetArray(), centrality);
        return fAdditionalPi0PtWeight[index];
      }
      else
        return fAdditionalPi0PtWeight[0]; 
    }

    TF1 *GetAdditionalK0SPtWeightFunction(Float_t centrality){
      if(fCentArrayK0S){
        Int_t lastBinUpperIndex = fCentArrayK0S->GetSize()-1;
        Int_t index = TMath::BinarySearch<Double_t>( lastBinUpperIndex, fCentArrayK0S->GetArray(), centrality);
        return fAdditionalK0SPtWeight[index];
      }
      else 
        return fAdditionalK0SPtWeight[0];
    }

    TF1 *GetAdditionalEtaPtWeightFunction(Float_t centrality){
      if(fCentArrayEta){
        Int_t lastBinUpperIndex = fCentArrayEta->GetSize()-1;
        Int_t index = TMath::BinarySearch<Double_t>( lastBinUpperIndex, fCentArrayEta->GetArray(), centrality);
        return fAdditionalEtaPtWeight[index];
      }
      else
        return fAdditionalEtaPtWeight[0]; 
    }

    TF1 *GetAdditionalGammaPtWeightFunction(Float_t centrality){
      if(fCentArrayGamma){
        Int_t lastBinUpperIndex = fCentArrayGamma->GetSize()-1;
        Int_t index = TMath::BinarySearch<Double_t>( lastBinUpperIndex, fCentArrayGamma->GetArray(), centrality);
        return fAdditionalGammaPtWeight[index];
      }
      else
        return fAdditionalGammaPtWeight[0]; 
    }

    Bool_t ExtractQnVector();

    AliStack *GetMCInfoESD();
    TClonesArray *GetMCInfoAOD();

    AliGenPythiaEventHeader* GetPythiaEventHeader(AliVEvent *event);
    Bool_t ComparePtHardWithJet(AliVEvent *event);
    Bool_t ComparePtHardWithSingleParticle(AliVEvent *event);
    Int_t FindPrimaryMotherESD(Int_t label);
    AliPHOSGeometry *GetPHOSGeometry();

    virtual void DoNonLinearityStudy();

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
    Int_t fCollisionSystem;//colliions system : pp=0, PbPb=1, pPb (Pbp)=2;
    TF1 *fTOFEfficiency;//TOF cut efficiency as a function of cluster energy;
    TF1 *fTriggerEfficiency;//TOF cut efficiency as a function of cluster energy;
    AliESDtrackCuts *fESDtrackCutsGlobal;//good global track
    AliESDtrackCuts *fESDtrackCutsGlobalConstrained;//global track but constrained to IP because of SPD dead area
    TF1 *fAdditionalPi0PtWeight[10];//weight function for pT distribution
    TF1 *fAdditionalEtaPtWeight[10];//weight function for pT distribution
    TF1 *fAdditionalGammaPtWeight[10];//weight function for pT distribution
    TF1 *fAdditionalK0SPtWeight[10];//weight function for pT distribution. note that this weight is aiming to reproduce K/pi ratio.
    TArrayD *fCentArrayPi0;
    TArrayD *fCentArrayEta;
    TArrayD *fCentArrayGamma;
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
    TString fEstimator;//V0[M|A|C], ZN[A|C], CL[0|1], HybridTrack, SPDTracklet
    AliMultSelection *fMultSelection;
    Float_t fCentralityMain;
    Float_t fCentralityMin;
    Float_t fCentralityMax;
    Int_t fNMixed;
    Double_t fVertex[3];
    Int_t fZvtx;
    Int_t fEPBin;
    Bool_t fIsFlowTask;
    Int_t fHarmonics;
    TString fQNormalization;
    Int_t fFM;//kEP or kSP
    Int_t fQnDetectorMain;
    Int_t fQnDetectorSub1;
    Int_t fQnDetectorSub2;
    AliQnCorrectionsManager *fFlowQnVectorMgr;
    TString fTPCEPName[3]; 
    TString fV0EPName[3]; 
    Double_t fEventPlane;
    TVector2 fQVector1;//x,y
    Int_t fNHybridTrack;
    Bool_t fIsPHOSTriggerAnalysis;
    Double_t fEnergyThreshold;
    AliPHOSTriggerHelper *fPHOSTriggerHelper;//for real PHOS triggered data analysis
    AliPHOSTriggerHelper *fPHOSTriggerHelperL0; //only for rejection factor in MB
    AliPHOSTriggerHelper *fPHOSTriggerHelperL1H;//only for rejection factor in MB
    AliPHOSTriggerHelper *fPHOSTriggerHelperL1M;//only for rejection factor in MB
    AliPHOSTriggerHelper *fPHOSTriggerHelperL1L;//only for rejection factor in MB
    Bool_t fForceActiveTRU;
    AliPIDResponse *fPIDResponse;
    Bool_t fIsNonLinStudy;
    TF1 *fNonLin[7][7];

  private:
    AliAnalysisTaskPHOSPi0EtaToGammaGamma(const AliAnalysisTaskPHOSPi0EtaToGammaGamma&);
    AliAnalysisTaskPHOSPi0EtaToGammaGamma& operator=(const AliAnalysisTaskPHOSPi0EtaToGammaGamma&);

    ClassDef(AliAnalysisTaskPHOSPi0EtaToGammaGamma, 40);
};

#endif
