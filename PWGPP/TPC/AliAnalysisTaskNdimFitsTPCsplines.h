#ifndef  AliAnalysisTaskNdimFitsTPCsplines_H
#define  AliAnalysisTaskNdimFitsTPCsplines_H


#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliESDtrackCuts.h"
#include "AliAnalysisTask.h"
#include "AliTimeRangeCut.h"
#include "AliPIDResponse.h"
#include "AliESDVertex.h"
#include "AliEventCuts.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "THnSparse.h"
#include "AliESDv0.h"
#include "TVector3.h"
#include "TList.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"


//___________________________________________________________________________________________________________________________________
class AliAnalysisTaskNdimFitsTPCsplines : public AliAnalysisTaskSE {
    
public:
    AliAnalysisTaskNdimFitsTPCsplines();
    AliAnalysisTaskNdimFitsTPCsplines(const char *name);
    virtual ~AliAnalysisTaskNdimFitsTPCsplines();
    virtual void  UserCreateOutputObjects();
    virtual void  UserExec (Option_t *option);
    virtual void  Terminate(Option_t *);
    void BinLogAxis (THnSparseF *h, Int_t axisNumber);

    //User Functions
    Bool_t    GetInputEvent ();
    Bool_t    PassedTrackSelectionV0daugh    (AliESDtrack *track);
    Bool_t    PassedPrimaryTrackSelection       (AliESDtrack *track);
    Bool_t    PassedNucleiTrackSelection     (AliESDtrack *track);
    Bool_t    PassedV0Selection              (AliESDv0 *V0);
    Bool_t    PassedLambdaSelection          (AliESDv0 *V0);
    Bool_t    PassedAntiLambdaSelection      (AliESDv0 *V0);
    Bool_t    PassedK0shortSelection         (AliESDv0 *V0);
    Bool_t    PassedGammaConversionSelection (AliESDv0 *V0);
    Double_t  GetPhiV                        (AliESDtrack *track1, AliESDtrack *track2);
    Double_t  MassLambda     (TVector3 Ppion,  TVector3 Pprot);
    Double_t  MassDielectron (TVector3 Pelec1, TVector3 Pelec2);
    Double_t  MassK0short    (TVector3 Ppion1, TVector3 Ppion2);

    //Standard Event Selection
    AliEventCuts fESDeventCuts;//

private:
    AliESDEvent      *fESDevent;//!
    AliPIDResponse   *fPIDResponse;//!
    AliESDtrackCuts  *fESDtrackCuts_V0daugh;//!
    AliESDtrackCuts  *fESDtrackCuts_Primary;//!
    AliESDtrackCuts  *fESDtrackCuts_Nuclei;//!
    AliAnalysisUtils *fUtils;//!
    TList            *fOutputList;//!
    TList            *fQAList;//!
    AliTimeRangeCut   fTimeRangeCut;//

    
    //Number of Events
    TH1F *hNumberOfEvents;//!
    
    //nDimensional Histograms: Raw TPC dE/dx vs. p,eta,centrality
    THnSparseF *hTPCdEdx_Electrons;//!
    THnSparseF *hTPCdEdx_Pions;//!
    THnSparseF *hTPCdEdx_Kaons;//!
    THnSparseF *hTPCdEdx_Protons;//!
    THnSparseF *hTPCdEdx_Tritons;//!
    THnSparseF *hTPCdEdx_Helium3;//!

    
    AliAnalysisTaskNdimFitsTPCsplines(const AliAnalysisTaskNdimFitsTPCsplines&);
    AliAnalysisTaskNdimFitsTPCsplines& operator=(const AliAnalysisTaskNdimFitsTPCsplines&);
    
    ClassDef(AliAnalysisTaskNdimFitsTPCsplines, 1);
};
//___________________________________________________________________________________________________________________________________

#endif
