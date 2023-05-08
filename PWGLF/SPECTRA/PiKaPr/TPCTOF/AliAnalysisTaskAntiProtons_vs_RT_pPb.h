#ifndef AliAnalysisTaskAntiProtons_vs_RT_pPb_cxx
#define AliAnalysisTaskAntiProtons_vs_RT_pPb_cxx

//========================= Antiprotons vs R_{T} =========================//        
//                                                                        //
//      Antiprotons analysis vs. R_{T} in p-Pb collisions @ 5.02 TeV      //            
//          for the calculation of coalescence parameter B_{2}            //  
//                                                                        //
//========================================================================//

#include "AliMCEventHandler.h"
#include "AliAnalysisTaskSE.h"
#include "AliMultSelection.h"
#include "AliAnalysisUtils.h"
#include "AliESDtrackCuts.h"
#include "AliAnalysisTask.h"
#include "AliPIDResponse.h"
#include "AliMCParticle.h"
#include "AliEventCuts.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliMCEvent.h"
#include "THnSparse.h"
#include "TList.h"
#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"

//______________________________________________________________________________________________________________________________

class AliAnalysisTaskAntiProtons_vs_RT_pPb : public AliAnalysisTaskSE {

public:
    AliAnalysisTaskAntiProtons_vs_RT_pPb();
    AliAnalysisTaskAntiProtons_vs_RT_pPb(const char *name);
    virtual ~AliAnalysisTaskAntiProtons_vs_RT_pPb();

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

    //General Settings
    void SetTriggerType            (UInt_t triggerType)                 { fTriggerType = triggerType;}
    void SetMultiplicityInterval   (Double_t multMin, Double_t multMax) { fMultMin = multMin; fMultMax = multMax; }
    void SetAverageTransverseMult  (Double_t average_Nch_Transv)        { fAverage_Nch_Transv = average_Nch_Transv; }
    void SetInputData              (Bool_t isMC)                        { fIsMC = isMC; }
    void SetMinPtLeadingTrack      (Double_t pt)                        { fPt_min_leading = pt; }
    void SetParticleType           (Bool_t isPion)                      { fIsPion = isPion; }

    //Process Real and Simulated Event
    void ProcessRealEvent ();
    void ProcessSimEvent ();

    //User Functions
    Bool_t   GetESDEvent ();
    Bool_t   GetMCEvent ();
    Int_t    GetLeadingTrack();
    void     FillHistograms_StandardCuts                (Int_t mult_Transverse, Int_t leading_track_ID, AliESDtrack *track);
    void     FillHistograms_Systematics                 (Int_t mult_Transverse, Int_t leading_track_ID, AliESDtrack *track, Int_t isyst);
    void     FillHistograms_StandardCuts_Sim            (AliESDtrack *track);
    void     FillHistograms_Systematics_Sim             (AliESDtrack *track, Int_t isyst);
    Bool_t   PassedTrackQualityCuts_LeadingTrack        (AliESDtrack *track);
    Bool_t   PassedTrackQualityCuts_TransverseMult      (AliESDtrack *track);
    Bool_t   IsTrackInTransverseRegion                  (AliESDtrack *track, Int_t leading_track_ID);
    Bool_t   IsTrackInTowardRegion                      (AliESDtrack *track, Int_t leading_track_ID);
    Bool_t   IsTrackInAwayRegion                        (AliESDtrack *track, Int_t leading_track_ID);
    Bool_t   IsHighPurityProton                         (AliESDtrack *track);
    Bool_t   IsHighPurityPion                           (AliESDtrack* track);
    Bool_t   IsTrackCandidate                           (AliESDtrack *track);
    Bool_t   PassedTrackQualityCuts_Syst                (AliESDtrack *track, Int_t isyst);
    Double_t GetTransverseDCA                           (AliESDtrack *track);
    Double_t GetLongitudinalDCA                         (AliESDtrack *track);
    Double_t GetRapidity                                (AliESDtrack *track);
    Double_t GetRecalibratedITSnsigma                   (Double_t nsigma, Double_t eta, Double_t p);

    //Standard Event Cuts
    AliEventCuts  fESDEventSelection;//

private:
    AliESDEvent       *fESDevent;//!
    AliMCEvent        *fMCevent;//!
    AliMCEventHandler *fMCEventHandler;//!
    AliESDtrackCuts   *fESDtrackCuts_AntiProton;//!
    AliESDtrackCuts   *fESDtrackCuts_Pions;//!
    AliESDtrackCuts   *fESDtrackCuts_LeadingTrack;//!
    AliESDtrackCuts   *fESDtrackCuts_TransverseMult;//!
    AliPIDResponse    *fPIDResponse;//!
    TList             *fOutputList;//!
    TList             *fQAList;//!
    
    //Input from AddTask
    UInt_t    fTriggerType;//
    Double_t  fMultMin;//
    Double_t  fMultMax;//
    Double_t  fAverage_Nch_Transv;//
    Double_t  fPt_min_leading;//
    Bool_t    fIsMC;//
    Bool_t    fIsITSrecalib;//
    Bool_t    fIsPion;//
    TH2F      *hMean;//
    TH2F      *hWidth;//
    
    //Event Counter and Multiplicity Distributions
    TH1F *hNumberOfEvents;//!
    TH1F *hMultPercentile;//!
    TH1I *hMultTransverse;//!
    TH1I *hMultToward;//!
    TH1I *hMultAway;//!
    TH1F *hRtDistribution;//!
    TH1F *hEventsWithLeadingTrack;//!
    TH1I *hNumberOfAntiProtons;//!
    TH1I *hNumberOfPions;//!


    //Correlations between Transverse and Integrated Mult
    TH2F *hNchTransv_NchTot;//!
    TH2F *hRt_NchTot;//!
    TH2F *hNchTransv_MultPercentile;//!
    TH2F *hRt_MultPercentile;//!

    //3-Dimensional Histogram for ITS recalibration
    THnSparseF *hnsigmaITS_antiprotons;//!
    

    //***************************************** Data *****************************************

    //antiprotons histograms

    //3-Dimensional Histograms (low p_{T})
    THnSparseF *hnsigmaTPC_antiprotons_Toward;//!
    THnSparseF *hnsigmaTPC_antiprotons_Away;//!
    THnSparseF *hnsigmaTPC_antiprotons_Transverse;//!

    THnSparseF *hnsigmaTPC_antiprotons_Toward_noITSpresel;//!
    THnSparseF *hnsigmaTPC_antiprotons_Away_noITSpresel;//!
    THnSparseF *hnsigmaTPC_antiprotons_Transverse_noITSpresel;//!

    //3-Dimensional Histograms (high p_{T})
    THnSparseF *hnsigmaTOF_antiprotons_Toward;//!
    THnSparseF *hnsigmaTOF_antiprotons_Away;//!
    THnSparseF *hnsigmaTOF_antiprotons_Transverse;//!
    
    //DCA_{xy} Distributions
    THnSparseF *hDCAxy_antiprotons_Toward;//!
    THnSparseF *hDCAxy_antiprotons_Away;//!
    THnSparseF *hDCAxy_antiprotons_Transverse;//!
    
    //4-Dimensional Histograms for Syst. Uncertainties
    THnSparseF *hnsigmaTPC_antiprotons_Syst;//!
    THnSparseF *hnsigmaTOF_antiprotons_Syst;//!
    THnSparseF *hDCAxy_antiprotons_Syst;//!

    //pos and neg pions histograms

    //3-Dimensional Histograms (high p_{T})
    THnSparseF *hnsigmaTOF_pos_pions_Toward;//!
    THnSparseF *hnsigmaTOF_pos_pions_Away;//!
    THnSparseF *hnsigmaTOF_pos_pions_Transverse;//!

    THnSparseF *hnsigmaTOF_neg_pions_Toward;//!
    THnSparseF *hnsigmaTOF_neg_pions_Away;//!
    THnSparseF *hnsigmaTOF_neg_pions_Transverse;//!
    
    //DCA_{xy} Distributions
    THnSparseF *hDCAxy_pos_pions_Toward;//!
    THnSparseF *hDCAxy_pos_pions_Away;//!
    THnSparseF *hDCAxy_pos_pions_Transverse;//!

    THnSparseF *hDCAxy_neg_pions_Toward;//!
    THnSparseF *hDCAxy_neg_pions_Away;//!
    THnSparseF *hDCAxy_neg_pions_Transverse;//!

    //4-Dimensional Histograms for Syst. Uncertainties
    THnSparseF *hnsigmaTOF_pos_pions_Syst;//!
    THnSparseF *hDCAxy_pos_pions_Syst;//!

    THnSparseF *hnsigmaTOF_neg_pions_Syst;//!
    THnSparseF *hDCAxy_neg_pions_Syst;//!


    //****************************************** MC ******************************************
    
    //antiprotons histograms

    //Generated p_{T} Spectra
    TH1F *h_antiprotons_Gen;//!

    //Reconstructed p_{T} Spectra (low p_{T})
    TH2F *hnsigmaTPC_antiprotons_Rec;//!
        
    //Reconstructed p_{T} Spectra (high p_{T})
    TH2F *hnsigmaTOF_antiprotons_Rec;//!

    //DCA_{xy} distributions
    TH2F *hDCAxy_antiprotons_prim;//!
    TH2F *hDCAxy_antiprotons_sec;//!

    //Histograms for Syst. Uncertainties
    TH2F *hnsigmaTPC_antiprotons_Rec_Syst;//!
    TH2F *hnsigmaTOF_antiprotons_Rec_Syst;//!
    THnSparseF *hDCAxy_antiprotons_prim_Syst;//!
    THnSparseF *hDCAxy_antiprotons_sec_Syst;//! 

    //pos and neg pions histograms

    //Generated p_{T} Spectra
    TH1F *h_pos_pions_Gen;//!
    TH1F *h_neg_pions_Gen;//!

    //Reconstructed p_{T} Spectra
    TH2F *hnsigmaTOF_pos_pions_Rec;//!
    TH2F *hnsigmaTOF_neg_pions_Rec;//!

    //DCA_{xy} distributions
    TH2F *hDCAxy_pos_pions_prim;//!
    TH2F *hDCAxy_pos_pions_sec;//!
    TH2F *hDCAxy_neg_pions_prim;//!
    TH2F *hDCAxy_neg_pions_sec;//!

    //Histograms for Syst. Uncertainties
    TH2F *hnsigmaTOF_pos_pions_Rec_Syst;//!
    TH2F *hnsigmaTOF_neg_pions_Rec_Syst;//!
    THnSparseF *hDCAxy_pos_pions_prim_Syst;//! 
    THnSparseF *hDCAxy_neg_pions_prim_Syst;//! 
    THnSparseF *hDCAxy_pos_pions_sec_Syst;//! 
    THnSparseF *hDCAxy_neg_pions_sec_Syst;//! 

    //****************************************************************************************
    
    AliAnalysisTaskAntiProtons_vs_RT_pPb(const AliAnalysisTaskAntiProtons_vs_RT_pPb&);
    AliAnalysisTaskAntiProtons_vs_RT_pPb& operator=(const AliAnalysisTaskAntiProtons_vs_RT_pPb&);
        
    ClassDef(AliAnalysisTaskAntiProtons_vs_RT_pPb, 1);
};


#endif
