#ifndef AliAnalysisTaskDeuteronsRT_cxx
#define AliAnalysisTaskDeuteronsRT_cxx

//==================== DEUTERON ANALYSIS VS. RT =====================//
//                                                                   //
//      Deuteron analysis in pp collisions at 13 TeV vs. R_{T}.      //
//                                                                   //
//===================================================================//

#include "AliMCEventHandler.h"
#include "AliAnalysisTaskSE.h"
#include "AliESDtrackCuts.h"
#include "AliAnalysisTask.h"
#include "AliPIDResponse.h"
#include "AliMCParticle.h"
#include "AliEventCuts.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliMCEvent.h"
#include "THnSparse.h"
#include "TList.h"
#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"

//_________________________________________________________________________________________________________________________________________________
class AliAnalysisTaskDeuteronsRT : public AliAnalysisTaskSE {
        
public:
    AliAnalysisTaskDeuteronsRT();
    AliAnalysisTaskDeuteronsRT(const char *name);
    virtual ~AliAnalysisTaskDeuteronsRT();

    virtual void UserCreateOutputObjects();
    virtual void UserExec  (Option_t *option);
    virtual void Terminate (Option_t *);
    
    //General Settings
    void SetTriggerType            (UInt_t triggerType)                 { fTriggerType = triggerType;}
    void SetMultiplicityInterval   (Double_t multMin, Double_t multMax) { fMultMin = multMin; fMultMax = multMax; }
    void SetAverageTransverseMult  (Double_t average_Nch_Transv)        { fAverage_Nch_Transv = average_Nch_Transv; }
    void SetAnalysisParametersSyst (TH2F *h2Dmatrix)                    { hAnalysisParameters = h2Dmatrix; }
    void SetInputData              (Bool_t ispPb, Bool_t isMC)          { fIspPb = ispPb; fIsMC = isMC; }
    void SetMinPtLeadingTrack      (Double_t pt)                        { fPt_min_leading = pt; }
    void SetIsUEAnalysis           (Bool_t isUEanalysis)                { fIsUEanalysis = isUEanalysis; }

    //Process Real and Simulated Event
    void ProcessRealEvent ();
    void ProcessRealEventRapidityDependence();
    void ProcessSimEvent ();
    
    //User Functions
    Bool_t   GetESDEvent ();
    Bool_t   GetMCEvent ();
    Int_t    GetLeadingTrack();
    void     FillHistograms_StandardCuts                (Int_t mult_Transverse, Int_t leading_track_ID, AliESDtrack *track);
    void     FillHistograms_Systematics                 (Int_t mult_Transverse, Int_t leading_track_ID, AliESDtrack *track, Int_t isyst);
    void     FillHistograms_RapidityDependence          (AliESDtrack *track);
    void     FillHistograms_StandardCuts_Sim            (AliESDtrack *track);
    void     FillHistograms_Systematics_Sim             (AliESDtrack *track, Int_t isyst);
    void     FillHistograms_Rapidity_Systematics        (AliESDtrack *track, Int_t isyst);
    void     FillHistograms_Rapidity_Systematics_Sim    (AliESDtrack *track, Int_t isyst);
    Bool_t   PassedBasicTrackQualityCuts               (AliESDtrack *track);
    Bool_t   PassedBasicTrackQualityCuts_NoRapidityCut (AliESDtrack *track);
    Bool_t   PassedTrackQualityCuts_LeadingTrack       (AliESDtrack *track);
    Bool_t   PassedTrackQualityCuts_TransverseMult     (AliESDtrack *track);
    Bool_t   IsTrackInTransverseRegion           (AliESDtrack *track, Int_t leading_track_ID);
    Bool_t   IsTrackInTowardRegion               (AliESDtrack *track, Int_t leading_track_ID);
    Bool_t   IsTrackInAwayRegion                 (AliESDtrack *track, Int_t leading_track_ID);
    Bool_t   IsCleanDeuteron                     (AliESDtrack *track);
    Bool_t   IsDeuteronCandidate                 (AliESDtrack *track);
    Bool_t   PassedTrackQualityCuts_Syst         (AliESDtrack *track, Int_t isyst);
    Double_t GetTransverseDCA                    (AliESDtrack *track);
    Double_t GetLongitudinalDCA                  (AliESDtrack *track);
    Double_t GetRapidity (AliESDtrack *track);

    //Standard Event Cuts
    AliEventCuts  fESDeventCuts;//
    
private:
    AliESDEvent       *fESDevent;//!
    AliMCEvent        *fMCevent;//!
    AliMCEventHandler *fMCEventHandler;//!
    AliESDtrackCuts   *fESDtrackCuts_Deuteron;//!
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
    TH2F     *hAnalysisParameters;//
    Bool_t    fIspPb;//
    Bool_t    fIsMC;//
    Double_t  fPt_min_leading;//
    Bool_t    fIsUEanalysis;//
    
    //Event Counter and Multiplicity Distributions
    TH1F *hNumberOfEvents;//!
    TH1F *hMultPercentile;//!
    TH1I *hMultTransverse;//!
    TH1I *hMultToward;//!
    TH1I *hMultAway;//!
    TH1F *hRtDistribution;//!
    TH1F *hEventsWithLeadingTrack;//!
    TH1I *hNumberOfDeuterons;//
    
    //Correlations between Transverse and Integrated Mult
    TH2F *hNchTransv_NchTot;//!
    TH2F *hRt_NchTot;//!
    TH2F *hNchTransv_MultPercentile;//!
    TH2F *hRt_MultPercentile;//!
    

    //***************************************** Data *****************************************
    
    //3-Dimensional Histograms (low p_{T})
    THnSparseF *hnsigmaTPC_deuterons_Toward;//!
    THnSparseF *hnsigmaTPC_deuterons_Away;//!
    THnSparseF *hnsigmaTPC_deuterons_Transverse;//!
    THnSparseF *hnsigmaTPC_antideuterons_Toward;//!
    THnSparseF *hnsigmaTPC_antideuterons_Away;//!
    THnSparseF *hnsigmaTPC_antideuterons_Transverse;//!
    
    //3-Dimensional Histograms (high p_{T})
    THnSparseF *hnsigmaTOF_deuterons_Toward;//!
    THnSparseF *hnsigmaTOF_deuterons_Away;//!
    THnSparseF *hnsigmaTOF_deuterons_Transverse;//!
    THnSparseF *hnsigmaTOF_antideuterons_Toward;//!
    THnSparseF *hnsigmaTOF_antideuterons_Away;//!
    THnSparseF *hnsigmaTOF_antideuterons_Transverse;//!
    
    //DCA_{xy} Distributions
    THnSparseF *hDCAxy_deuterons_Toward;//!
    THnSparseF *hDCAxy_deuterons_Away;//!
    THnSparseF *hDCAxy_deuterons_Transverse;//!
    THnSparseF *hDCAxy_antideuterons_Toward;//!
    THnSparseF *hDCAxy_antideuterons_Away;//!
    THnSparseF *hDCAxy_antideuterons_Transverse;//!
    
    //4-Dimensional Histograms for Syst. Uncertainties
    THnSparseF *hnsigmaTPC_deuterons_Syst;//!
    THnSparseF *hnsigmaTPC_antideuterons_Syst;//!
    THnSparseF *hnsigmaTOF_deuterons_Syst;//!
    THnSparseF *hnsigmaTOF_antideuterons_Syst;//!
    THnSparseF *hDCAxy_deuterons_Syst;//!
    THnSparseF *hDCAxy_antideuterons_Syst;//!

    
    //3-Dimensional Histograms for Deuterons vs. y
    THnSparseF *hnsigmaTPC_deuterons_rap;//!
    THnSparseF *hnsigmaTPC_antideuterons_rap;//!
    THnSparseF *hnsigmaTOF_deuterons_rap;//!
    THnSparseF *hnsigmaTOF_antideuterons_rap;//!
    
    //3-Dimensional Histograms for Deuterons vs. y
    THnSparseF *hnsigmaTPC_deuterons_rap_Syst;//!
    THnSparseF *hnsigmaTPC_antideuterons_rap_Syst;//!
    THnSparseF *hnsigmaTOF_deuterons_rap_Syst;//!
    THnSparseF *hnsigmaTOF_antideuterons_rap_Syst;//!

    //****************************************** MC ******************************************
    
    //Generated p_{T} Spectra
    TH1F *h_deuterons_Gen;//!
    TH1F *h_antideuterons_Gen;//!

    //Reconstructed p_{T} Spectra (low p_{T})
    TH2F *hnsigmaTPC_deuterons_Rec;//!
    TH2F *hnsigmaTPC_antideuterons_Rec;//!
        
    //Reconstructed p_{T} Spectra (high p_{T})
    TH2F *hnsigmaTOF_deuterons_Rec;//!
    TH2F *hnsigmaTOF_antideuterons_Rec;//!

    //DCA_{xy} Distributions
    TH2F *hDCAxy_deuterons_Prim;//!
    TH2F *hDCAxy_deuterons_Sec;//!
    
    //Histograms for Syst. Uncertainties
    TH2F *hnsigmaTPC_deuterons_Rec_Syst;//!
    TH2F *hnsigmaTPC_antideuterons_Rec_Syst;//!
    TH2F *hnsigmaTOF_deuterons_Rec_Syst;//!
    TH2F *hnsigmaTOF_antideuterons_Rec_Syst;//!

    //Efficiency vs. Rapidity
    TH2F *hGeneratedDeuterons_vs_Rapidity;//!
    TH2F *hGeneratedAntiDeuterons_vs_Rapidity;//!
    TH2F *hReconstructedDeuterons_TPC_vs_Rapidity;//!
    TH2F *hReconstructedAntiDeuterons_TPC_vs_Rapidity;//!
    TH2F *hReconstructedDeuterons_TOF_vs_Rapidity;//!
    TH2F *hReconstructedAntiDeuterons_TOF_vs_Rapidity;//!
    
    //Histograms for Syst. Uncertainties in Rapidity analysis
    THnSparseF *hnsigmaTPC_deuterons_Rec_rap_Syst;//!
    THnSparseF *hnsigmaTPC_antideuterons_Rec_rap_Syst;//!
    THnSparseF *hnsigmaTOF_deuterons_Rec_rap_Syst;//!
    THnSparseF *hnsigmaTOF_antideuterons_Rec_rap_Syst;//!

    //****************************************************************************************

    
    AliAnalysisTaskDeuteronsRT(const AliAnalysisTaskDeuteronsRT&);
    AliAnalysisTaskDeuteronsRT& operator=(const AliAnalysisTaskDeuteronsRT&);
        
    ClassDef(AliAnalysisTaskDeuteronsRT, 1);
    
};
//_________________________________________________________________________________________________________________________________________________

#endif
