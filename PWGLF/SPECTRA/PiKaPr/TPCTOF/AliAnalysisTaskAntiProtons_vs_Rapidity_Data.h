#ifndef AliAnalysisTaskAntiProtons_vs_Rapidity_Data_cxx
#define AliAnalysisTaskAntiProtons_vs_Rapidity_Data_cxx

//======================== Antiprotons vs. Rapidity ========================//
//                                                                          //
//    Antiproton production vs. transverse momentum and rapidity            //
//    for the calculation of the coalescence parameters B_{2} and B_{3}.    //
//                                                                          //
//==========================================================================//

#include "AliAnalysisTaskSE.h"
#include "AliMultSelection.h"
#include "AliAnalysisUtils.h"
#include "AliESDtrackCuts.h"
#include "AliAnalysisTask.h"
#include "AliPIDResponse.h"
#include "AliESDVertex.h"
#include "AliEventCuts.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "THnSparse.h"
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"

//__________________________________________________________________________________________________________
class AliAnalysisTaskAntiProtons_vs_Rapidity_Data : public AliAnalysisTaskSE {
        
public:
    AliAnalysisTaskAntiProtons_vs_Rapidity_Data();
    AliAnalysisTaskAntiProtons_vs_Rapidity_Data(const char *name);
    virtual ~AliAnalysisTaskAntiProtons_vs_Rapidity_Data();
        
    //ITS Re-calibration
    void SetRunningMode (Bool_t isITSrecalib) {fIsITSrecalib = isITSrecalib;}
    
    //Set ITS Recalibration maps
    void SetITSRecalibrationMaps (TH2F *hITSnsigma_Mean, TH2F *hITSnsigma_Width)  {
        
        hMean  = hITSnsigma_Mean;
        hWidth = hITSnsigma_Width;
    }
    
    //General Functions
    virtual void UserCreateOutputObjects();
    virtual void UserExec  (Option_t *option);
    virtual void Terminate (Option_t *);
    
    //User Functions
    Bool_t   GetEvent ();
    Bool_t   PassedTrackSelection     (AliESDtrack *track, Int_t isyst);
    Bool_t   IsHighPurityProton       (AliESDtrack *track);
    Bool_t   IsProtonCandidate        (AliESDtrack *track);
    Double_t GetDCAtoPrimaryVertex    (AliESDtrack *track, Int_t index);
    Double_t GetRapidity              (AliESDtrack *track, Double_t mass);
    Double_t GetRecalibratedITSnsigma (Double_t nsigma, Double_t eta, Double_t p);

    //Standard Event Selection
    AliEventCuts  fESDEventSelection;//

private:
    AliESDEvent     *fESDEvent;//!
    AliPIDResponse  *fPIDResponse;//!
    AliESDtrackCuts *fESDtrackCuts[50];//!
    TList           *fOutputList;//!
    TList           *fQAList;//!
    Bool_t           fIsITSrecalib;//
    TH2F *hMean;//
    TH2F *hWidth;//

    //Event Counter and Centrality Distribution
    TH1F *hNumberOfEvents;//!
    TH1F *hMultiplicity;//!
    
    //n-Dimensional Histograms
    THnSparse *hTPCnsigma;//!
    THnSparse *hTOFnsigma;//!
    THnSparse *hDCAxy;//!
    
    //2D ITS Recalibration Map
    TH3F *hITSnsigma;//!
    
    //n-Dimensional Histograms (y>0 vs. y<0)
    THnSparse *hTPCnsigma_vs_rap;//!
    THnSparse *hTOFnsigma_vs_rap;//!
    
    
    AliAnalysisTaskAntiProtons_vs_Rapidity_Data(const AliAnalysisTaskAntiProtons_vs_Rapidity_Data&);
    AliAnalysisTaskAntiProtons_vs_Rapidity_Data& operator=(const AliAnalysisTaskAntiProtons_vs_Rapidity_Data&);
        
    ClassDef(AliAnalysisTaskAntiProtons_vs_Rapidity_Data, 1);
    
};
//__________________________________________________________________________________________________________

#endif
