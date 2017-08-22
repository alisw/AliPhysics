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
//class AliEventCuts;

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
    void SetEventCuts(AliPHOSEventCuts *cuts) {fPHOSEventCuts = cuts;}
    void SetClusterCuts(AliPHOSClusterCuts *cuts) {fPHOSClusterCuts = cuts;}
    void SetBunchSpace(Double_t bs) {fBunchSpace = bs;}
    void SetCollisionSystem(Int_t id) {fCollisionSystem = id;}

    void SetPtHardBin(Int_t pThardbin) {
      fPtHardBin = pThardbin;
      fIsJJMC = kTRUE;
      fJJMCHandler = new AliPHOSJetJetMC(pThardbin);
    }

    void SetTOFCutEfficiencyFunction(TF1 *f1) {fTOFEfficiency = f1;}
    void SetAdditionalPi0PtWeightFunction(TF1 *f1) {fAdditionalPi0PtWeight = f1;}
    void SetAdditionalK0SPtWeightFunction(TF1 *f1) {fAdditionalK0SPtWeight = f1;}

    void SetCentralityBinningPbPb(const TArrayD& edges, const TArrayI& nMixed)
    {
      // Define centrality bins by their edges
      for(int i=0; i<edges.GetSize()-1; ++i)
        if(edges.At(i) > edges.At(i+1)) AliFatal("edges are not sorted");
      if(edges.GetSize() != nMixed.GetSize()+1) AliFatal("edges and nMixed don't have appropriate relative sizes");//Pb-Pb,pPb

      fCentEdges = edges;
      fCentNMixed = nMixed;
    }

    void SetCentralityBinningPP(const TArrayD& edges, const TArrayI& nMixed)
    {
      // Define centrality bins by their edges
      for(int i=0; i<edges.GetSize()-1; ++i)
        if(edges.At(i) > edges.At(i+1)) AliFatal("edges are not sorted");
      if(edges.GetSize() != nMixed.GetSize()) AliFatal("edges and nMixed don't have appropriate relative sizes");//pp

      fCentEdges = edges;
      fCentNMixed = nMixed;
    }

  void SetCentralityEstimator(TString estimator) {fEstimator = estimator;}

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

    virtual void EstimateNcellCutEfficiency();
    virtual void EstimateM20CutEfficiency();
    virtual void EstimateM02CutEfficiency();
    virtual void EstimateSTDCutEfficiency();

    void EstimateTOFCutEfficiency();
    void EstimateTriggerEfficiency();
    void SelectTriggeredCluster();

    Bool_t IsFromLambda0(Int_t label1, Int_t label2, Double_t &TrueL0Pt);
    Bool_t IsFromK0S(Int_t label1, Int_t label2, Double_t &TrueK0SPt);
    Bool_t IsFromPi0(Int_t label1, Double_t &TruePi0Pt);
    Bool_t IsPhoton(Int_t label);

    Double_t R(AliAODMCParticle *p);//in cylindrical system
    Double_t Rho(AliAODMCParticle *p);//in sperical system

    virtual Int_t FindCommonParent(Int_t iPart, Int_t jPart);

    virtual void DoNonLinearityStudy();

    void FillHistogramTH1(TList *list, const Char_t *name, Double_t x, Double_t w=1., Option_t *opt = "") const ;
    void FillHistogramTH2(TList *list, const Char_t *name, Double_t x, Double_t y, Double_t w=1., Option_t *opt = "") const ;
    void FillHistogramTH3(TList *list, const Char_t *name, Double_t x, Double_t y, Double_t z, Double_t w=1., Option_t *opt = "") const ;
    void FillProfile(TList *list, const Char_t *name, Double_t x, Double_t y) const ;

    TF1 *GetTOFCutEfficiencyFunction() {return fTOFEfficiency;}
    TF1 *GetAdditionalPi0PtWeightFunction() {return fAdditionalPi0PtWeight;}
    TF1 *GetAdditionalK0SPtWeightFunction() {return fAdditionalK0SPtWeight;}

    Int_t GetCentralityBinPbPb(Float_t centrality)
    {
      if(centrality < fCentEdges[0]) return -1;
      Int_t lastBinUpperIndex = fCentEdges.GetSize()-1;

      //if(centrality >= fCentEdges[lastBinUpperIndex]) return lastBinUpperIndex-1;
      //if(centrality >= fCentEdges[lastBinUpperIndex]) return lastBinUpperIndex;
      if(centrality >= fCentEdges[lastBinUpperIndex]) return -1;
      else return TMath::BinarySearch<Double_t>( lastBinUpperIndex, fCentEdges.GetArray(), centrality);
    }

    Int_t GetCentralityBinPP(Int_t multiplicity)
    {
      if((Double_t)multiplicity < fCentEdges[0]) return -1;//this should not happen. only for protection
      Int_t lastBinUpperIndex = fCentEdges.GetSize()-1;

      if((Double_t)multiplicity >= fCentEdges[lastBinUpperIndex]) return lastBinUpperIndex;
      else return TMath::BinarySearch<Double_t>(lastBinUpperIndex, fCentEdges.GetArray(), (Double_t)multiplicity);
    }

    AliStack *GetMCInfoESD();
    TClonesArray *GetMCInfoAOD();
    AliGenPythiaEventHeader* GetPythiaEventHeader(AliVEvent *event);
    Bool_t ComparePtHardWithJet(AliVEvent *event);
    Bool_t ComparePtHardWithSingleParticle(AliVEvent *event);

    AliPHOSGeometry *GetPHOSGeometry();

  protected:
    Bool_t fIsMC;
    Bool_t fIsJJMC;//jet jet MC
    Int_t fPtHardBin;//your selected pT hard bin for jet jet MC
    Bool_t fUsePHOSTender;
    Bool_t fUseCoreEnergy;
    AliPHOSEventCuts *fPHOSEventCuts;
    AliPHOSClusterCuts *fPHOSClusterCuts;
    Double_t fBunchSpace;// in unit of ns.
    TArrayD fCentEdges;  // Centrality Bin Lower edges
    TArrayI fCentNMixed; // Number of mixed events for each centrality bin
    Int_t fCollisionSystem;//colliions system : pp=0, PbPb=1, pPb (Pbp)=2;
    TF1 *fTOFEfficiency;//TOF cut efficiency as a function of cluster energy;
    AliESDtrackCuts *fESDtrackCutsGlobal;//good global track
    AliESDtrackCuts *fESDtrackCutsGlobalConstrained;//global track but constrained to IP because of SPD dead area
    TF1 *fAdditionalPi0PtWeight;//weight function for pT distribution
    TF1 *fAdditionalK0SPtWeight;//weight function for pT distribution

    THashList *fOutputContainer;
    AliVEvent *fEvent;
    AliESDEvent *fESDEvent;
    AliAODEvent *fAODEvent;
    AliStack *fMCArrayESD;     //MC particles array in ESD
    TClonesArray *fMCArrayAOD; //MC particles array in AOD
    AliPHOSJetJetMC *fJJMCHandler;
    Int_t fRunNumber;
    AliPHOSGeometry *fPHOSGeo;
    TList *fPHOSEvents[10][10];
    TClonesArray *fPHOSClusterArray;
    TString fEstimator;
    AliMultSelection *fMultSelection;
    Float_t fCentralityV0M;
    Float_t fCentralityCL0;
    Float_t fCentralityCL1;
    Float_t fCentralityV0A;
    Float_t fCentralityV0C;
    Float_t fCentralityZNA;
    Float_t fCentralityZNC;
    Float_t fCentralityMain;
    Double_t fVertex[3];
    Int_t fZvtx;
    Int_t fCenBin;
    Int_t fNHybridTrack;
    Bool_t fIsNonLinStudyNeeded;
    TF1 *fNonLin[7][7][7];
//    AliEventCuts fEventCuts;

  private:
    AliAnalysisTaskPHOSPi0EtaToGammaGamma(const AliAnalysisTaskPHOSPi0EtaToGammaGamma&); // not implemented
    AliAnalysisTaskPHOSPi0EtaToGammaGamma& operator=(const AliAnalysisTaskPHOSPi0EtaToGammaGamma&); // not implemented

    ClassDef(AliAnalysisTaskPHOSPi0EtaToGammaGamma, 13); // example of analysis
};

#endif
