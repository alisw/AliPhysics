/// \class AliAnalysisTaskStrAODqa
/// \brief This task fills histograms required to perform the QA analysis on V0s and cascades.
/// \authors Chiara De Martin <chiara.de.martin@cern.ch>, University and INFN Trieste - Alessandro Balbino <alessandro.balbino@cern.ch>, Politecnico and INFN Torino
/// \date May 8, 2020

#ifndef AliAnalysisTaskStrAODqa_H
#define AliAnalysisTaskStrAODqa_H

#include "TString.h"
#include "AliPIDResponse.h"
#include "AliAnalysisTaskSE.h"
#include "THistManager.h"

class AliAnalysisTaskStrAODqa : public AliAnalysisTaskSE {
public:
    AliAnalysisTaskStrAODqa();
    AliAnalysisTaskStrAODqa(const char *name, TString lExtraOptions = "");
    virtual ~AliAnalysisTaskStrAODqa();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);

    void SetMC(Bool_t isMC){fReadMCTruth = isMC;}

private:
    THistManager* fHistos_eve;   //!
    THistManager* fHistos_V0;   //!
    THistManager* fHistos_Casc;   //!

    //objects retreived from input handler
    AliPIDResponse *fPIDResponse;     //! PID response object

    //variables for MC analysis
    AliMCEvent *            fMCEvent;         //!                                                                                                       
    Bool_t                  fReadMCTruth;

    //variables for V0 analysis
    double fV0_DcaV0Daught;  //!
    double fV0_DcaPosToPV;   //!
    double fV0_DcaNegToPV;   //!
    double fV0_V0CosPA;      //!
    double fV0_V0Rad;        //!
    double fV0_Pt;           //!
    double fV0_yK0S;         //!
    double fV0_yLam;         //!
    double fV0_etaPos;       //!
    double fV0_etaNeg;       //!
    double fV0_InvMassK0s;   //!
    double fV0_InvMassLam;   //!
    double fV0_InvMassALam;  //!
    double fV0_LeastCRaws;   //!
    double fV0_LeastCRawsOvF;//!
    double fV0_NSigPosProton;//!
    double fV0_NSigPosPion;  //!
    double fV0_NSigNegProton;//!
    double fV0_NSigNegPion;  //!
    double fV0_DistOverTotP; //!

    //variables for Cascade analysis
    bool   fCasc_isNotTPCRefit;      //!
    double fCasc_DcaCascDaught;      //!
    double fCasc_CascCosPA;          //!
    double fCasc_CascRad;            //!
    double fCasc_etaPos;             //!
    double fCasc_etaNeg;             //!
    double fCasc_etaBac;             //!
    double fCasc_kinkidx;            //!
    double fCasc_NSigPosProton;      //!
    double fCasc_NSigPosPion;        //!
    double fCasc_NSigNegProton;      //!
    double fCasc_NSigNegPion;        //!
    double fCasc_NSigBacPion;        //!
    double fCasc_NSigBacKaon;        //!
    double fCasc_LeastCRaws;         //!
    double fCasc_LeastCRawsOvF;      //!
    double fCasc_DcaV0Daught;        //!
    double fCasc_V0CosPA;            //!
    double fCasc_V0CosPAToXi;        //!
    double fCasc_DcaV0ToPV;          //!
    double fCasc_DcaBachToPV;        //!
    double fCasc_yXi;                //!
    double fCasc_yOm;                //!
    int    fCasc_charge;             //!
    double fCasc_Pt;                 //!
    double fCasc_Ptot;               //!
    double fCasc_DistOverTotP;       //!
    double fCasc_V0DistOverTotP;     //!
    double fCasc_CascCtauXi;         //!
    double fCasc_CascCtauOmega;      //!
    double fCasc_V0Ctau;             //!
    double fCasc_InvMassXi;          //!
    double fCasc_InvMassOm;          //!
    double fCasc_V0Rad;              //!
    double fCasc_DcaPosToPV;         //!
    double fCasc_DcaNegToPV;         //!
    double fCasc_InvMassLambda;      //!

    //list of analysed particles
    //    static const int kNParticles = 7;
    //    const char *kParticleNames[kNParticles]= {"K0s", "Lambda", "anti-Lambda", "Xi-", "Xi+", "Om-", "Om+"}; 

    //cuts application
    bool ApplyCuts(int, bool, bool);
    bool ApplyCutsNSigmaTPC(int);

    AliAnalysisTaskStrAODqa(const AliAnalysisTaskStrAODqa&);            // not implemented
    AliAnalysisTaskStrAODqa& operator=(const AliAnalysisTaskStrAODqa&); // not implemented

    ClassDef(AliAnalysisTaskStrAODqa, 1);
    //1: first implementation
};

#endif
