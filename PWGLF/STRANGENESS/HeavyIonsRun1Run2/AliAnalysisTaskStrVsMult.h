#ifndef AliAnalysisTaskStrVsMult_H
#define AliAnalysisTaskStrVsMult_H

//class AliESDpid;
//class AliESDEvent;

#include "TString.h"
#include "AliPIDResponse.h"
#include "AliAnalysisTaskSE.h"
#include "THistManager.h"
#include "AliEventCuts.h"
#include "AliESDtrackCuts.h"

class AliAnalysisTaskStrVsMult : public AliAnalysisTaskSE {
  public:
    AliAnalysisTaskStrVsMult();
    AliAnalysisTaskStrVsMult(const char *name, TString lExtraOptions = "");
    virtual ~AliAnalysisTaskStrVsMult();

    // enum and names.
    enum cutnumb_V0{kV0_DcaV0Daught, kV0_DcaPosToPV, kV0_DcaNegToPV, kV0_V0CosPA, kV0_V0Rad, kV0_MaxV0Rad, kV0_y, kV0_etaDaugh, kV0_LeastCRaws, kV0_LeastCRawsOvF, kV0_TrackLengthCut, kV0_MaxChi2perCls, kV0_NSigPID, kV0_PropLifetK0s, kV0_PropLifetLam, kV0_ITSTOFtracks, kV0_ArmenterosK0s, kV0cutsnum};
    enum cutnumb_Casc{kCasc_DcaCascDaught, kCasc_CascCosPA, kCasc_CascRad, kCasc_NSigPID, kCasc_LeastCRaws, kCasc_LeastCRawsOvF, kCasc_TrackLengthCut, kCasc_MaxChi2perCls, kCasc_InvMassLam, kCasc_DcaV0Daught, kCasc_V0CosPA, kCasc_DcaV0ToPV, kCasc_DcaBachToPV, kCasc_ITSTOFtracks, kCasc_y, kCasc_etaDaugh, kCasc_PropLifetXi, kCasc_PropLifetOm, kCasc_CompetingXiMass, kCasc_V0Rad, kCasc_MaxV0Rad, kCasc_DcaMesToPV, kCasc_DcaBarToPV, kCasc_BacBarCosPA, kCasccutsnum};
    enum particles{kK0s, kLam, kXi, kOm, knumpart}; 
    enum signedparticles{kk0s, klam, kalam, kxip, kxim, komp, komm, ksignednumpart}; 

    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *);

    //cut values setter
    void SetDefOnly(bool);
    void SetCutVal(bool, bool, int, double);
    void SetParametricBacBarCosPA(int, float*, float*, int);
    void SetParametricTrackLengthCut(int, float*, int*);

    //TrackLength Cut setters
    void SetDeadZoneWidthGeoCut(float DeadZoneWidth){fDeadZoneWidth_GeoCut = DeadZoneWidth;};
    void SetNcrNclLengthGeoCut(float NcrNclLength){fNcrNclLength_GeoCut = NcrNclLength;};
    void SetTPCsignalNCut(int TPCsignalNCut){fTPCsignalNCut = TPCsignalNCut;};

    //binning setters
    void SetCentbinning(int, double*);
    void SetMassbinning(int, int, double, double);
    void SetPtbinning(int, int, double*);

    //set and get which particle specie are analysed
    void SetParticleAnalysisStatus(bool, bool, bool, bool);
    bool GetParticleAnalysisStatus(int);

    //MC-related setters and getters
    void SetIsMC(bool IsMC){fisMC = IsMC;};
    void SetIsMCassoc(bool IsMCassoc){fisMCassoc = IsMCassoc;};
    void SetIsMaterialAnalysis(bool IsMaterialAnalysis){fisMaterialAnalysis = IsMaterialAnalysis;};

    //pile-up rejection setter
    void SetRejectPileUpEvts(bool RejectPileupEvts, int PileupCut=1){if(RejectPileupEvts==kTRUE) fPileupCut = PileupCut;};

    //centrality estimator setter
    void SetCentralityEstimator(TString CentEstimator){fCentEstimator = CentEstimator;};

  private:
    THistManager* fHistos_eve;                                //!
    THistManager* fHistos_K0S;                                //!
    THistManager* fHistos_Lam;                                //!
    THistManager* fHistos_ALam;                               //!
    THistManager* fHistos_XiMin;                              //!
    THistManager* fHistos_XiPlu;                              //!
    THistManager* fHistos_OmMin;                              //!
    THistManager* fHistos_OmPlu;                              //!

    //objects retreived from input handler
    AliPIDResponse *fPIDResponse;                             //!
    UInt_t fTriggerMask;                                      //!

    //AliEventCuts object
    AliEventCuts fEventCuts;                                  //

    //pile-up rejection flag
    int fPileupCut;                                           //

    //Centrality estimator
    TString fCentEstimator;                                   //
    
    //MC-realted variables
    bool fisMC;                                               //
    bool fisMCassoc;                                          //
    bool fisMaterialAnalysis;                                 //

    //Default cut configuration
    bool fDefOnly;                                            //
    double fV0_Cuts[kV0cutsnum];                              //
    double fCasc_Cuts[kCasccutsnum];                          //

    //particles to be analysed
    bool fParticleAnalysisStatus[ksignednumpart];             //

    //geometrical cut usage
    AliESDtrackCuts fESDTrackCuts;                            //

    //variables for V0 analysis
    double fV0_DcaV0Daught;                                   //!
    double fV0_DcaPosToPV;                                    //!
    double fV0_DcaNegToPV;                                    //!
    double fV0_V0CosPA;                                       //!
    double fV0_V0Rad;                                         //!
    double fV0_Pt;                                            //!
    double fV0_yK0S;                                          //!
    double fV0_yLam;                                          //!
    double fV0_etaPos;                                        //!
    double fV0_etaNeg;                                        //!
    double fV0_InvMassK0s;                                    //!
    double fV0_InvMassLam;                                    //!
    double fV0_InvMassALam;                                   //!
    double fV0_LeastCRaws;                                    //!
    double fV0_LeastCRawsOvF;                                 //!
    int fV0_TrackLengthCut;                                   //!
    double fV0_MaxChi2perCls;                                 //!
    double fV0_NSigPosProton;                                 //!
    double fV0_NSigPosPion;                                   //!
    double fV0_NSigNegProton;                                 //!
    double fV0_NSigNegPion;                                   //!
    double fV0_DistOverTotP;                                  //!
    int fV0_ITSTOFtracks;                                     //!
    ULong64_t fV0_NegTrackStatus;                             //!
    ULong64_t fV0_PosTrackStatus;                             //!
    double fV0_kinkidx;                                       //!
    double fV0_AlphaArm;                                      //!
    double fV0_pTArm;                                         //!

    //variables for Cascade analysis
    double fCasc_DcaCascDaught;                               //!
    double fCasc_CascCosPA;                                   //!
    double fCasc_CascRad;                                     //!
    double fCasc_etaPos;                                      //!
    double fCasc_etaNeg;                                      //!
    double fCasc_etaBac;                                      //!
    double fCasc_kinkidx;                                     //!
    double fCasc_NSigPosProton;                               //!
    double fCasc_NSigPosPion;                                 //!
    double fCasc_NSigNegProton;                               //!
    double fCasc_NSigNegPion;                                 //!
    double fCasc_NSigBacPion;                                 //!
    double fCasc_NSigBacKaon;                                 //!
    double fCasc_LeastCRaws;                                  //!
    double fCasc_LeastCRawsOvF;                               //!
    int fCasc_TrackLengthCut;                                 //!
    double fCasc_MaxChi2perCls;                               //!
    double fCasc_InvMassLam;                                  //!
    double fCasc_DcaV0Daught;                                 //!
    double fCasc_V0CosPA;                                     //!
    double fCasc_DcaV0ToPV;                                   //!
    double fCasc_DcaBachToPV;                                 //!
    int fCasc_ITSTOFtracks;                                   //!
    double fCasc_yXi;                                         //!
    double fCasc_yOm;                                         //!
    int fCasc_charge;                                         //!
    double fCasc_Pt;                                          //!
    double fCasc_DistOverTotP;                                //!
    double fCasc_InvMassXiMin;                                //!
    double fCasc_InvMassXiPlu;                                //!
    double fCasc_InvMassOmMin;                                //!
    double fCasc_InvMassOmPlu;                                //!
    double fCasc_V0Rad;                                       //!
    double fCasc_DcaPosToPV;                                  //!
    double fCasc_DcaNegToPV;                                  //!
    ULong64_t fCasc_NegTrackStatus;                           //!
    ULong64_t fCasc_PosTrackStatus;                           //!
    ULong64_t fCasc_BacTrackStatus;                           //!
    double fCasc_BacBarCosPA;                                 //!

    bool fisParametricBacBarCosPA;                            //
    TH1F *fHist_PtBacBarCosPA;                                //
    int fCentLimit_BacBarCosPA;                               //

    bool fisParametricTrackLengthCut;                         //
    TH1I *fHist_CentTrackLengthCut;                           //
    float fDeadZoneWidth_GeoCut;                              //
    float fNcrNclLength_GeoCut;                               //
    int fTPCsignalNCut;                                       //

    //cut values to be set
    double cutval_V0[kV0cutsnum];                             //
    int nvarcut_V0[kV0cutsnum];                               //
    double varlowcut_V0[kV0cutsnum];                          //
    double varhighcut_V0[kV0cutsnum];                         //
    double cutval_Casc[kCasccutsnum];                         //
    int nvarcut_Casc[kCasccutsnum];                           //
    double varlowcut_Casc[kCasccutsnum];                      //
    double varhighcut_Casc[kCasccutsnum];                     //

    //variables to handle binning
    int fncentbins;                                           //
    double fcentbinning[150];                                 //
    int fnmassbins[knumpart];                                 //
    double fmassbinning[knumpart][1000];                      //
    int fnptbins[knumpart];                                   //
    double fptbinning[knumpart][2000];                        //

    //functions to allow flushing part of code out of UserExec
    bool ApplyCuts(int);
    void DataPosting();
    void FillHistCutVariations(bool, double, bool, bool*, double);
    //functions to allow the correct streaming of the cut variation
    void SetDefCutVals();
    void SetCutVariation(bool, int, int, double, double);
    void SetDefCutVariations();

    AliAnalysisTaskStrVsMult(const AliAnalysisTaskStrVsMult&);            // not implemented
    AliAnalysisTaskStrVsMult& operator=(const AliAnalysisTaskStrVsMult&); // not implemented

    ClassDef(AliAnalysisTaskStrVsMult, 21); 
    //version 21: add histogram for secondaries from material
};

#endif
