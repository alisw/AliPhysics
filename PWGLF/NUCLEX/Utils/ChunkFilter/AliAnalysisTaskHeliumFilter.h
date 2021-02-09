#ifndef AliAnalysisTaskHeliumFilter_cxx
#define AliAnalysisTaskHeliumFilter_cxx


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
#include "TString.h"
#include "TList.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"

//___________________________________________________________________________________________________________________________________
class AliAnalysisTaskHeliumFilter : public AliAnalysisTaskSE {
    
public:
    AliAnalysisTaskHeliumFilter();
    AliAnalysisTaskHeliumFilter(const char *name);
    virtual ~AliAnalysisTaskHeliumFilter();
    virtual void   UserCreateOutputObjects();
    virtual void   UserExec (Option_t *option);
    
    Bool_t   GetInputEvent ();
    Bool_t   PassedBasicTrackQualityCuts (AliESDtrack *track);
    Bool_t   IsHeliumCandidate           (AliESDtrack *track);
    
    virtual void   Terminate(Option_t *);
    
private:
    AliESDEvent      *fESDevent;//!
    AliPIDResponse   *fPIDResponse;//!
    AliESDtrackCuts  *fESDtrackCuts;//!
    AliEventCuts      fESDeventCuts;//
    AliTimeRangeCut   fTimeRangeCut;//
    AliAnalysisUtils *fUtils;//!
    TList            *fOutputList;//!

    
    //Number of Events 
    TH1F *hEvents;//!
    
    //Histograms
    TH2F *hdEdx_vs_p;//!
    TH2F *hdEdx_vs_p_Helium;//!
    
    //Reduced Tree
    TTree *tree_ListOfFiles;//!
    Int_t   fEventIdFile;//! Event id in file
    TString fFileName;//! Chunk File Name

    //Bethe-Bloch parameterizations
    Double_t fParamDeuteron[5];
    Double_t fParamTriton[5];
    Double_t fParamHe3[5];
    Double_t fParamDeuteronMC[5];
    Double_t fParamTritonMC[5];
    Double_t fParamHe3MC[5];
    
    
    AliAnalysisTaskHeliumFilter(const AliAnalysisTaskHeliumFilter&);
    AliAnalysisTaskHeliumFilter& operator=(const AliAnalysisTaskHeliumFilter&);
    
    ClassDef(AliAnalysisTaskHeliumFilter, 1);
};
//___________________________________________________________________________________________________________________________________

#endif
