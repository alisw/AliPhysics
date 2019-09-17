#ifndef AliAnalysisTaskReducedTreeHypertritonBindingEnergy_cxx
#define AliAnalysisTaskReducedTreeHypertritonBindingEnergy_cxx


#include "AliAnalysisTaskSE.h"
#include "AliMultSelection.h"
#include "AliAnalysisUtils.h"
#include "AliESDtrackCuts.h"
#include "AliAnalysisTask.h"
#include "AliTimeRangeCut.h"
#include "AliPIDResponse.h"
#include "AliESDVertex.h"
#include "AliEventCuts.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "TVector2.h"
#include "TVector3.h"
#include "AliESDv0.h"
#include "TList.h"
#include "TTree.h"
#include "TH1F.h"


class AliAnalysisTaskReducedTreeHypertritonBindingEnergy : public AliAnalysisTaskSE {
    
public:
    AliAnalysisTaskReducedTreeHypertritonBindingEnergy();
    AliAnalysisTaskReducedTreeHypertritonBindingEnergy(const char *name);
    virtual ~AliAnalysisTaskReducedTreeHypertritonBindingEnergy();
    virtual void   UserCreateOutputObjects();
    virtual void   UserExec (Option_t *option);
    
    void     SetCentrality (Double_t centrMin, Double_t centrMax)  { fCentralityMin = centrMin; fCentralityMax = centrMax; }
    Bool_t   GetInputEvent ();
    Bool_t   PassedBasicTrackQualityCuts_Pos (AliESDtrack *track);
    Bool_t   PassedBasicTrackQualityCuts_Neg (AliESDtrack *track);
    Double_t GetTransverseDCA                (AliESDtrack *track);
    Bool_t   PassedMinimalQualityCutsV0      (AliESDv0 *V0);
    Bool_t   IsHyperTritonCandidate          (AliESDv0 *V0);
    Double_t GetDecayLengthV0                (AliESDv0 *V0);
    Bool_t   Is3HeCandidate                  (AliESDtrack *track);
    Bool_t   IsPionCandidate                 (AliESDtrack *track);
    Double_t InvariantMassHypertriton        (TVector3 P1, TVector3 P2);
    
    virtual void   Terminate(Option_t *);
    
private:
    AliESDEvent      *fESDevent;//!
    AliPIDResponse   *fPIDResponse;//!
    AliESDtrackCuts  *fESDtrackCuts_Pos;//!
    AliESDtrackCuts  *fESDtrackCuts_Neg;//!
    AliEventCuts      fESDeventCuts;//
    AliTimeRangeCut   fTimeRangeCut;//
    AliAnalysisUtils *fUtils;//!
    TList            *fOutputList;//!
    TList            *fQAList;//!
    Double_t          fCentralityMin;//
    Double_t          fCentralityMax;//
    
    
    
    
    //Event Selection Histogram
    TH1F *hEvents;//!
    
    
    //Reduced Tree
    TTree *reducedTree_HyperTriton;//!
    
    //Global Variables
    Int_t    iEvent;//
    Double_t zVertex;//
    Double_t centrality;//
    
    //Variables for HyperTriton - First Daughter
    Double_t px_Daughter1;//
    Double_t py_Daughter1;//
    Double_t pz_Daughter1;//
    Int_t    q_Daughter1;//
    Double_t dcaxy_Daughter1;//
    Int_t    nTPC_Clusters_Daughter1;//
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
    Int_t    nTPC_Clusters_dEdx_Daughter2;//
    Double_t chi2_TPC_Daughter2;//
    Double_t nSigmaTPC_He3_Daughter2;//
    Double_t nSigmaTPC_Pion_Daughter2;//
    
    //Pair Variables
    Int_t    isOnTheFlyV0;//
    Double_t cosPointingAngle;//
    Double_t dcaV0Daughters;//
    Double_t radius;//
    Double_t chi2V0;//
    Double_t decayLength;//
    
    
    AliAnalysisTaskReducedTreeHypertritonBindingEnergy(const AliAnalysisTaskReducedTreeHypertritonBindingEnergy&);
    AliAnalysisTaskReducedTreeHypertritonBindingEnergy& operator=(const AliAnalysisTaskReducedTreeHypertritonBindingEnergy&);
    
    ClassDef(AliAnalysisTaskReducedTreeHypertritonBindingEnergy, 1);
};
#endif
