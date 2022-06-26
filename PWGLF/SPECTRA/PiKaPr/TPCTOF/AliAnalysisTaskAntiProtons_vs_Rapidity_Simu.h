#ifndef AliAnalysisTaskAntiProtons_vs_Rapidity_Simu_cxx
#define AliAnalysisTaskAntiProtons_vs_Rapidity_Simu_cxx

//======================== Antiprotons vs. Rapidity ========================//
//                                                                          //
//    Antiproton production vs. transverse momentum and rapidity            //
//    for the calculation of the coalescence parameters B_{2} and B_{3}.    //
//    Efficiency and DCA_{xy} templates for antiprotons.                    //
//                                                                          //
//==========================================================================//

#include "AliAnalysisTaskSE.h"
#include "AliMCEventHandler.h"
#include "AliMultSelection.h"
#include "AliAnalysisUtils.h"
#include "AliAnalysisTask.h"
#include "AliESDtrackCuts.h"
#include "AliPIDResponse.h"
#include "AliMCParticle.h"
#include "AliEventCuts.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliMCEvent.h"
#include "THnSparse.h"
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"

//__________________________________________________________________________________________________________
class AliAnalysisTaskAntiProtons_vs_Rapidity_Simu : public AliAnalysisTaskSE {
        
public:
    AliAnalysisTaskAntiProtons_vs_Rapidity_Simu();
    AliAnalysisTaskAntiProtons_vs_Rapidity_Simu(const char *name);
    virtual ~AliAnalysisTaskAntiProtons_vs_Rapidity_Simu();
        
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
    Bool_t   GetSimEvent ();
    Bool_t   GetRecEvent ();
    Bool_t   IsINELgtzero();
    Bool_t   PassedTrackSelection     (AliESDtrack *track, Int_t isyst);
    Bool_t   IsHighPurityProton       (AliESDtrack *track);
    Double_t GetDCAtoPrimaryVertex    (AliESDtrack *track, Int_t index);
    Double_t GetRapidity              (AliESDtrack *track, Double_t mass);
    Double_t GetRecalibratedITSnsigma (Double_t nsigma, Double_t eta, Double_t p);

    //Standard Event Selection
    AliEventCuts  fESDEventSelection;//

private:
    AliESDEvent       *fESDEvent;//!
    AliMCEvent        *fMCEvent;//!
    AliMCEventHandler *fMCEventHandler;//!
    AliESDtrackCuts   *fESDtrackCuts[50];//!
    AliPIDResponse    *fPIDResponse;//!
    TList             *fOutputList;//!
    TList             *fQAList;//!
    
    //ITS Recalibration Maps
    TH2F *hMean;//
    TH2F *hWidth;//

    //Event Counter and Multiplicity Distribution
    TH1F *hNumberOfEvents;//!
    TH1F *hMultiplicity;//!
    
    //n-Dimensional Histograms
    THnSparse *hGen;//!
    THnSparse *hTPCnsigma;//!
    THnSparse *hTOFnsigma;//!
    THnSparse *hDCAxy_prim;//!
    THnSparse *hDCAxy_sec;//!
    
    //2D ITS Recalibration Map
    TH3F *hITSnsigma;//!
    
    //Pile-up Correction
    TH2F *hGen_pileup;//!
    TH2F *hRec_TPC_pileup;//!
    TH2F *hRec_TOF_pileup;//!

    //Checks
    TH1F *hGen1d;//!
    TH1F *hRec_TPC1d;//!
    TH1F *hRec_TOF1d;//!
    TH2F *hGen2d;//!
    TH2F *hRec_TPC2d;//!
    TH2F *hRec_TOF2d;//!
    
    //INEL
    TH1F *hGenINELgtzero;//!
    TH1F *hGenINEL;//!
    
    AliAnalysisTaskAntiProtons_vs_Rapidity_Simu(const AliAnalysisTaskAntiProtons_vs_Rapidity_Simu&);
    AliAnalysisTaskAntiProtons_vs_Rapidity_Simu& operator=(const AliAnalysisTaskAntiProtons_vs_Rapidity_Simu&);
        
    ClassDef(AliAnalysisTaskAntiProtons_vs_Rapidity_Simu, 1);
    
};
//__________________________________________________________________________________________________________

#endif
