#ifndef AliAnalysisTaskPHOSObjectCreator_cxx
#define AliAnalysisTaskPHOSObjectCreator_cxx

// Author: Daiki Sekihata (Hiroshima University)
class TF1;
class TObjArray;
class TClonesArray;
class AliPHOSGeometry;
class TParticle;
class AliStack;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskPHOSObjectCreator : public AliAnalysisTaskSE {
  public:

    AliAnalysisTaskPHOSObjectCreator(const char *name = "PHOSObjectCreator");
    virtual ~AliAnalysisTaskPHOSObjectCreator(); 

    //AnalysisTaskSE method is overloaded.
    void UserCreateOutputObjects();
    void UserExec(Option_t *option);
    void Terminate(Option_t *option);
    //at least, these 3 functions above are necessary.

    void SetTenderFlag(Bool_t tender) {fUsePHOSTender = tender;}
    void SetMCFlag(Bool_t mc) {fIsMC = mc;}
    void SetBunchSpace(Double_t bs) {fBunchSpace = bs;}

    void SetPHOSBadMap(Int_t mod,TH2I *h)
    {
      if(fPHOSBadMap[mod]) delete fPHOSBadMap[mod];
      fPHOSBadMap[mod] = new TH2I(*h);
      AliInfo(Form("Setting Bad Map Histogram  %s",fPHOSBadMap[mod]->GetName()));
    }

    void SetUserDefinedNonlinearity(TF1 *f1nonlin) {fUserNonLinCorr = f1nonlin;}
    void ExcludeM4(Bool_t flag) {fIsM4Excluded = flag;}

    void SetSingleSim(Bool_t flag) {fIsSingleSim = flag;}
    void SetEmbedding(Bool_t flag) {fIsEmbedding = flag;}

  protected:
    void Init();
    void InitBadMap();//this is for internal

    void SetGeometry();
    Bool_t IsGoodChannel(const char * det, Int_t mod,Int_t ix, Int_t iz);
    Double_t TestCPVRun2(Double_t dx, Double_t dz, Double_t pt, Int_t charge);
    Double_t TestLambda(Double_t e,Double_t l1,Double_t l2);
    Double_t TestCoreLambda(Double_t e,Double_t l1,Double_t l2);
    void EvalCoreLambdas(AliVCluster * clu, AliVCaloCells * cells, Double_t &m02, Double_t &m20);
    void EvalLambdas(AliVCluster * clu, Double_t &m02, Double_t &m20);
    Double_t CoreEnergy(AliVCluster *clu, AliVCaloCells *cells);
    Double_t CoreEnergy(AliVCluster *clu);
    Int_t FindTrackMatching(Int_t mod,TVector3 *locpos,Double_t &dx, Double_t &dz, Double_t &pttrack, Int_t &charge);
    void DistanceToBadChannel(Int_t mod, TVector3 * locPos, Double_t &minDist) ;
    Int_t FindPrimary(AliCaloPhoton *ph,  Bool_t&sure);

    void EstimateSTDCutEfficiency(TClonesArray *array);
    Bool_t PassSTDCut(AliVCluster *cluster);

  protected:
    THashList *fOutputContainer;
    TH2F *fHistoMggvsEProbe;
    TH2F *fHistoMggvsEPassingProbe;
    AliVEvent *fEvent;
    AliAODEvent *fAODEvent;
    AliESDEvent *fESDEvent;
    Double_t fVertex[3];
    TClonesArray *fPHOSObjectArray;
    AliPHOSGeometry *fPHOSGeo;
    TH2I *fPHOSBadMap[6];
    TF1 *fNonLinCorr;
    TF1 *fUserNonLinCorr;
    Int_t fRunNumber;
    Bool_t fUsePHOSTender;
    Bool_t fIsMC;
    Double_t fBunchSpace;
    AliStack *fMCArrayESD;
    TClonesArray *fMCArrayAOD;
    Bool_t fIsM4Excluded;
    Bool_t fIsSingleSim;
    Bool_t fIsEmbedding;

  private:
    AliAnalysisTaskPHOSObjectCreator(const AliAnalysisTaskPHOSObjectCreator&);
    AliAnalysisTaskPHOSObjectCreator& operator=(const AliAnalysisTaskPHOSObjectCreator&);

    ClassDef(AliAnalysisTaskPHOSObjectCreator, 17);
};

#endif

