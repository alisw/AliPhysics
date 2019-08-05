#ifndef AliAnalysisTaskReducedTreeHypertriton_cxx
#define AliAnalysisTaskReducedTreeHypertriton_cxx


#include "AliAnalysisTaskSE.h"
#include "AliMultSelection.h"
#include "AliAnalysisUtils.h"
#include "AliAnalysisTask.h"
#include "AliPIDResponse.h"
#include "AliAODVertex.h"
#include "AliEventCuts.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "TObjArray.h"
#include "TVector2.h"
#include "TVector3.h"
#include "AliAODv0.h"
#include "TList.h"
#include "TTree.h"
#include "TH1F.h"

class AliAnalysisTaskReducedTreeHypertriton : public AliAnalysisTaskSE {
    
public:
    AliAnalysisTaskReducedTreeHypertriton();
    AliAnalysisTaskReducedTreeHypertriton(const char *name);
    virtual ~AliAnalysisTaskReducedTreeHypertriton();
    
    virtual void   UserCreateOutputObjects();
    virtual void   UserExec (Option_t *option);
    
    void SetCentrality (Double_t centralityMin, Double_t centralityMax)  {
        fcentralityMin = centralityMin;
        fcentralityMax = centralityMax;
    }
    
    Bool_t   GetInputEvent ();
    Bool_t   PassedBasicTrackQualityCuts (AliAODTrack *track);
    Bool_t   PassedV0QualityCuts         (AliAODv0 *V0);
    Bool_t   IsHyperTritonCandidate      (AliAODv0 *V0);
    Bool_t   Is3HeCandidate              (AliAODTrack *track);
    Bool_t   IsPionCandidate             (AliAODTrack *track);

    virtual void   Terminate(Option_t *);
    
private:
    AliAODEvent      *fAODevent;//!
    AliPIDResponse   *fPIDResponse;//!
    AliEventCuts      fAODeventCuts;// 
    AliAnalysisUtils *fUtils;//!
    TList            *fOutputList;//!
    TList            *fQAList;//!
    Double_t          fcentralityMin;//
    Double_t          fcentralityMax;//

    //Event Selection Histogram
    TH1F *hEvents;//!
    
    //Reduced Tree
    TTree *reducedTree_HyperTriton;//!
    
    //Global Variables
    Double_t centrality;//

    //Variables for HyperTriton - First Daughter
    Double_t px_Daughter1;//
    Double_t py_Daughter1;//
    Double_t pz_Daughter1;//
    Int_t    q_Daughter1;//
    Double_t dcaxy_Daughter1;//
    Int_t    nTPC_Clusters_Daughter1;//
    Int_t    nTPC_FindableClusters_Daughter1;//
    Int_t    nTPC_CrossedRows_Daughter1;//
    Int_t    nTPC_Clusters_dEdx_Daughter1;//
    Double_t chi2_TPC_Daughter1;//
    Double_t nSigmaTPC_He3_Daughter1;//
    Double_t nSigmaTPC_Pion_Daughter1;//

    //Variables for HyperTriton - Second Daughter
    Double_t px_Daughter2;//
    Double_t py_Daughter2;//
    Double_t pz_Daughter2;//
    Int_t    q_Daughter2;//
    Double_t dcaxy_Daughter2;//
    Int_t    nTPC_Clusters_Daughter2;//
    Int_t    nTPC_FindableClusters_Daughter2;//
    Int_t    nTPC_CrossedRows_Daughter2;//
    Int_t    nTPC_Clusters_dEdx_Daughter2;//
    Double_t chi2_TPC_Daughter2;//
    Double_t nSigmaTPC_He3_Daughter2;//
    Double_t nSigmaTPC_Pion_Daughter2;//
    
    //Pair Variables
    Int_t    isOnTheFlyV0;//
    Double_t cosPointingAngle;//
    Double_t dcaV0Daughters;//
    Double_t dcaV0ToVertex;//
    Double_t radius;//
    Double_t decayLength;//

    
    AliAnalysisTaskReducedTreeHypertriton(const AliAnalysisTaskReducedTreeHypertriton&);
    AliAnalysisTaskReducedTreeHypertriton& operator=(const AliAnalysisTaskReducedTreeHypertriton&);
    
    ClassDef(AliAnalysisTaskReducedTreeHypertriton, 1);
};
#endif
