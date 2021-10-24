#ifndef AliAnalysisTaskStrVsMult_BumpStudies_H
#define AliAnalysisTaskStrVsMult_BumpStudies_H

//class AliESDpid;
//class AliESDEvent;

#include "TString.h"
#include "AliPIDResponse.h"
#include "AliAnalysisTaskSE.h"
#include "THistManager.h"

class AliAnalysisTaskStrVsMult_BumpStudies : public AliAnalysisTaskSE {
  public:
    AliAnalysisTaskStrVsMult_BumpStudies();
    AliAnalysisTaskStrVsMult_BumpStudies(const char *name, TString lExtraOptions = "");
    virtual ~AliAnalysisTaskStrVsMult_BumpStudies();

    // enum and names.
    enum cutnumb_V0{kV0_DcaV0Daught, kV0_DcaPosToPV, kV0_DcaNegToPV, kV0_V0CosPA, kV0_V0Rad, kV0_y, kV0_etaDaugh, kV0_LeastCRaws, kV0_LeastCRawsOvF, kV0_NSigPID, kV0_PropLifetK0s, kV0_PropLifetLam, kV0_TOFBunchCrossing, kV0cutsnum}; 
    enum cutnumb_Casc{kCasc_DcaCascDaught, kCasc_CascCosPA, kCasc_CascRad, kCasc_NSigPID, kCasc_LeastCRaws, kCasc_LeastCRawsOvF, kCasc_InvMassLam, kCasc_DcaV0Daught, kCasc_V0CosPA, kCasc_DcaV0ToPV, kCasc_DcaBachToPV, kCasc_TOFBunchCrossing, kCasc_y, kCasc_etaDaugh, kCasc_PropLifetXi, kCasc_PropLifetOm, kCasc_V0Rad, kCasc_DcaMesToPV, kCasc_DcaBarToPV, kCasc_BacBarCosPA, kCasccutsnum}; // kCasc_etaPos, kCasc_etaNeg, kCasc_etaBac, kCasc_kinkidx,
    enum particles{kK0s, kLam, kXi, kOm, knumpart}; 
    enum signedparticles{kk0s, klam, kalam, kxip, kxim, komp, komm, ksignednumpart}; 

    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *);

    //cut values setter
    void SetDefOnly(bool);
    void SetCutVal(bool, bool, int, double);

    //binning setters
    void SetCentbinning(int, int, double*);
    void SetMassbinning(int, int, double, double);
    void SetPtbinning(int, int, double*);

    //set and get which particle specie are analysed
    void SetParticleAnalysisStatus(bool, bool, bool, bool);
    bool GetParticleAnalysisStatus(int);

    //MC-related setters and getters
    void SetIsMC(bool IsMC){fisMC = IsMC;};
    void SetIsMCantiassoc(bool IsMCantiassoc){fisMCantiassoc = IsMCantiassoc;};

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

    //MC-realted variables
    bool fisMC;                                               //
    bool fisMCantiassoc;                                          //

    //Default cut configuration
    bool fDefOnly;                                            //
    double fV0_Cuts[kV0cutsnum];                              //
    double fCasc_Cuts[kCasccutsnum];                          //

    //particles to be analysed
    bool fParticleAnalysisStatus[ksignednumpart];             //

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
    double fV0_NSigPosProton;                                 //!
    double fV0_NSigPosPion;                                   //!
    double fV0_NSigNegProton;                                 //!
    double fV0_NSigNegPion;                                   //!
    double fV0_DistOverTotP;                                  //!
    double fV0_NegTOFBunchCrossing;                           //!
    double fV0_PosTOFBunchCrossing;                           //!
    ULong64_t fV0_NegTrackStatus;                             //!
    ULong64_t fV0_PosTrackStatus;                             //!
    double fV0_kinkidx;                                       //!

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
    double fCasc_InvMassLam;                                  //!
    double fCasc_DcaV0Daught;                                 //!
    double fCasc_V0CosPA;                                     //!
    double fCasc_DcaV0ToPV;                                   //!
    double fCasc_DcaBachToPV;                                 //!
    double fCasc_PosTOFBunchCrossing;                         //!
    double fCasc_NegTOFBunchCrossing;                         //!
    double fCasc_BacTOFBunchCrossing;                         //!
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
    double fCasc_BacBarCosPA;                                   //!

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
    int fncentbins[knumpart];                                 //
    double fcentbinning[knumpart][50];                        //
    int fnmassbins[knumpart];                                 //
    double fmassbinning[knumpart][1000];                      //
    int fnptbins[knumpart];                                   //
    double fptbinning[knumpart][600];                         //

    //functions to allow flushing part of code out of UserExec
    bool ApplyCuts(int);
    void DataPosting();
    void FillHistCutVariations(bool, double, bool, bool*, double);
    //functions to allow the correct streaming of the cut variation
    void SetDefCutVals();
    void SetCutVariation(bool, int, int, double, double);
    void SetDefCutVariations();

    AliAnalysisTaskStrVsMult_BumpStudies(const AliAnalysisTaskStrVsMult_BumpStudies&);            // not implemented
    AliAnalysisTaskStrVsMult_BumpStudies& operator=(const AliAnalysisTaskStrVsMult_BumpStudies&); // not implemented

    ClassDef(AliAnalysisTaskStrVsMult_BumpStudies, 1); 
    //version 4: introduced MC handeling
};

#endif
