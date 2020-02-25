#ifndef ALIANALYSISTASKSEXICZERO2XIPIFROMKFP_H
#define ALIANALYSISTASKSEXICZERO2XIPIFROMKFP_H

/// \class AliAnalysisTaskSEXicZero2XiPifromKFP
/// \brief This is a brief description of my class
///
/// This is a longer description of the class. This longer description is
/// formatted using the Markdown syntax (see below) and can span on multiple
/// lines.
///
/// \author Jianhui Zhu <zjh@mail.ccnu.edu.cn>, Central China Normal University & GSI Helmholtz Centre for Heavy Ion Research
/// \date Apr 26, 2019

/* $Id$ */

#ifndef HomogeneousField
#define HomogeneousField
#endif

#include "AliAnalysisTaskSE.h"
#include "AliAODMCParticle.h"
#include "AliRDHFCutsKFP.h"
#include "AliNormalizationCounter.h"
#include "THnSparse.h"
#include "AliPIDResponse.h"
#include "AliAODInputHandler.h"

// includes added to play with KFParticle
#include <vector>
#include "KFParticleBase.h"
#include "KFParticle.h"
#include "KFPTrack.h"
#include "KFPVertex.h"
#include "KFVertex.h"

class AliAnalysisTaskSEXicZero2XiPifromKFP : public AliAnalysisTaskSE  
{
    public:
                                AliAnalysisTaskSEXicZero2XiPifromKFP();
                                AliAnalysisTaskSEXicZero2XiPifromKFP(const char *name, AliRDHFCutsKFP* cuts);
        virtual                 ~AliAnalysisTaskSEXicZero2XiPifromKFP();

        virtual void            UserCreateOutputObjects();
        virtual void            Init();
        virtual void            LocalInit() {Init();} 
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);

        void                    SetMC(Bool_t IsMC) {fIsMC=IsMC;}
        void                    SelectTrack(AliVEvent *event, Int_t trkEntries, Int_t &nSeleTrks, Bool_t *seleFlags);
        Bool_t                  MakeMCAnalysis(TClonesArray *mcArray);
        void                    MakeAnaXicZeroFromV0(AliAODEvent *AODEvent, TClonesArray *mcArray, KFParticle PV);
        void                    MakeAnaXicZeroFromCasc(AliAODEvent *AODEvent, TClonesArray *mcArray, KFParticle PV);
        Int_t                   MatchToMCXic0(AliAODTrack *trackProton, AliAODTrack *trackPion3, AliAODTrack *trackPion2, AliAODTrack *trackAntiPion1, TClonesArray *mcArray);
        Int_t                   MatchToMCAntiXic0(AliAODTrack *trackAntiProton, AliAODTrack *trackAntiPion3, AliAODTrack *trackAntiPion2, AliAODTrack *trackPion1, TClonesArray *mcArray);
        Int_t                   MatchToMCXiMinus(AliAODTrack *trackProton, AliAODTrack *trackPion3, AliAODTrack *trackPion2, TClonesArray *mcArray);
        Int_t                   MatchToMCXiPlus(AliAODTrack *trackAntiProton, AliAODTrack *trackAntiPion3, AliAODTrack *trackAntiPion2, TClonesArray *mcArray);
        Int_t                   MatchToMCLambda(AliAODTrack *trackProton, AliAODTrack *trackPion3, TClonesArray *mcArray);
        Int_t                   MatchToMCAntiLambda(AliAODTrack *trackAntiProton, AliAODTrack *trackAntiPion3, TClonesArray *mcArray);
        Int_t                   MatchToMCLambdaFromXi(AliAODTrack *trackProton, AliAODTrack *trackPion3, TClonesArray *mcArray);
        Int_t                   MatchToMCAntiLambdaFromXi(AliAODTrack *trackAntiProton, AliAODTrack *trackAntiPion3, TClonesArray *mcArray);
        Int_t                   MatchToMCPion(AliAODTrack *track, TClonesArray *mcArray);
        Double_t                InvMassV0atPV(AliAODTrack *trk1, AliAODTrack *trk2, Int_t pdg1, Int_t pdg2);
        Double_t                CosPointingAngleKF(KFParticle kfp, KFParticle kfpmother);
        Double_t                CosThetaStarKF(Int_t ip, UInt_t pdgvtx, UInt_t pdgprong0, UInt_t pdgprong1, KFParticle kfpvtx, KFParticle kfpprong0, KFParticle kfpprong1);
        Bool_t                  CheckVertexCov(AliAODVertex *vtx);
        Bool_t                  CheckTrackCov(AliAODTrack *track);
        Bool_t                  CheckKFParticleCov(KFParticle kfp);
        KFParticle              CreateKFTrack(Double_t *param, Double_t *cov, Float_t Chi2perNDF, Int_t charge, Int_t pdg);
        KFVertex                CreateKFVertex(Double_t *param, Double_t *cov);
        KFParticle              CreateKFParticleFromAODtrack(AliAODTrack *track, Int_t pdg);
        KFParticle              CreateKFMotherParticle(AliAODTrack *track1, AliAODTrack *track2, Int_t pdg1, Int_t pdg2);
        KFParticle              CreateSecKFParticle(KFParticle kfp1, AliAODTrack *track2, Int_t pdg1, Int_t pdg2);
        Int_t                   MatchToXicZeroMC(TClonesArray *mcArray, Int_t PDGXicZero, const Int_t nDaughters, const Int_t *daughterIndex, const Int_t *daughterPDG);

        /// set MC usage
        void SetWriteXic0MCGenTree(Bool_t a) {fWriteXic0MCGenTree = a;}
        Bool_t GetWriteXic0MCGenTree() const {return fWriteXic0MCGenTree;}

        void SetWriteXic0Tree(Bool_t a) {fWriteXic0Tree = a;}
        Bool_t GetWriteXic0Tree() const {return fWriteXic0Tree;}

        void FillEventROOTObjects();
        void FillTreeGenXic0(AliAODMCParticle *mcpart, AliAODMCParticle *mcpipart, AliAODMCParticle *mccascpart, Int_t decaytype);
        void FillTreeRecXic0FromV0(KFParticle kfpXicZero, AliAODTrack *trackPi, KFParticle kfpBP, KFParticle kfpXiMinus, KFParticle kfpXiMinus_m, AliAODTrack *trackPiFromXi, AliAODv0 *v0, KFParticle kfpK0Short, KFParticle kfpLambda, KFParticle kfpLambda_m, AliAODTrack *trkP, AliAODTrack *trkN, KFParticle PV, TClonesArray *mcArray, Int_t lab_Xic0);
        void FillTreeRecXic0FromCasc(KFParticle kfpXicZero, AliAODTrack *trackPi, KFParticle kfpBP, KFParticle kfpXiMinus, KFParticle kfpXiMinus_m, AliAODTrack *trackPiFromXi, AliAODcascade *casc, KFParticle kfpK0Short, KFParticle kfpLambda, KFParticle kfpLambda_m, AliAODTrack *trkP, AliAODTrack *trkN, KFParticle PV, TClonesArray *mcArray, Int_t lab_Xic0);

    private:
        void                    DefineEvent();
        void                    DefineTreeRecXic0();
        void                    DefineTreeGenXic0();
        void                    DefineAnaHist();
        AliPIDResponse*         fPID;                 ///<
        AliRDHFCutsKFP*         fAnaCuts;             ///< Cuts
        AliAODVertex*           fpVtx;                //!<! primary vertex
        AliMCEvent*             fMCEvent;             //!<! corresponding mc event
        Double_t                fBzkG;                ///< magnetic field value [kG]
        Float_t                 fCentrality;           ///< Centrality
        Int_t                   fRunNumber;            ///< Run Number
        Int_t                   fEvNumberCounter;      ///< EvNumber counter
//        TObjArray               fMapParticle;         ///< Map of particles in the supporting TClonesArray
        vector<Int_t>           fAodTrackInd;         ///< Translation table: aodTrackInd(mcTrackIndex) = aodTrackIndex
        TList*                  fOutputList;          //!<! Output list
        TTree*                  fTree_Event;          //!<! tree of event
        Float_t*                fVar_Event;           //!<! variables of event to be written to the tree
        TTree*                  fTree_Xic0;             //!<! tree of the candidate variables
        Float_t*                fVar_Xic0;         //!<! variables of Xic0 to be written to the tree
        TTree*                  fTree_Xic0MCGen; //!<! tree of the candidate variables after track selection on output slot
        Float_t*                fVar_Xic0MCGen;   //!<! variables to be written to the tree
//        TTree*                  fVarTree_AntiXicZero;             //!<! tree of the candidate variables
//        Float_t*                fVar_AntiXicZero;         //!<! variables of Anti-Xic0 to be written to the tree
        TList*                  fListCuts;           //!<! User output slot 3 // Cuts 

        Bool_t                  fIsMC;                ///< Flag of MC analysis

        AliNormalizationCounter* fCounter; //!<! Counter for normalization
        TH1F*                   fHistMCGen_Lambda_Pt; //!<! Pt distribution of lambda at gen. level
        TH1F*                   fHistMCGen_AntiLambda_Pt; //!<! Pt distribution of lambda at gen. level
        TH1F*                   fHistMCGen_Lambda_Pt_wYcut; //!<! Pt distribution of lambda at gen. level
        TH1F*                   fHistMCGen_AntiLambda_Pt_wYcut; //!<! Pt distribution of lambda at gen. level
        TH1F*                   fCounterGen_Cuts_Lambda; //!<! 1D counter for Lambda (Cuts)
        TH1F*                   fCounterGen_Cuts_AntiLambda; //!<! 1D counter for Anti-Lambda (Cuts)
        TH1F*                   fCounterRecMC_Cuts_Lambda; //!<! 1D counter for Lambda (Cuts)
        TH1F*                   fCounterRecMC_Cuts_AntiLambda; //!<! 1D counter for Anti-Lambda (Cuts)
        TH1F*                   fCounterRec_Cuts_Lambda; //!<! 1D counter for Lambda (Cuts)
        TH1F*                   fCounterRec_Cuts_AntiLambda; //!<! 1D counter for Anti-Lambda (Cuts)
        TH1F*                   fHistMCGen_XiMinus_Pt; //!<! Pt distribution of Xi- at gen. level
        TH1F*                   fHistMCGen_XiPlus_Pt; //!<! Pt distribution of Xi+ at gen. level
        TH1F*                   fHistMCGen_XiMinus_Pt_wYcut; //!<! Pt distribution of Xi- at gen. level
        TH1F*                   fHistMCGen_XiPlus_Pt_wYcut; //!<! Pt distribution of Xi+ at gen. level
        TH1F*                   fCounterGen_Cuts_XiMinus; //!<! 1D counter for Xi- (Cuts)
        TH1F*                   fCounterGen_Cuts_XiPlus; //!<! 1D counter for Xi+ (Cuts)
        TH1F*                   fCounterRecMC_Cuts_XiMinus; //!<! 1D counter for Xi- (Cuts)
        TH1F*                   fCounterRecMC_Cuts_XiPlus; //!<! 1D counter for Xi+ (Cuts)
        TH1F*                   fCounterRec_Cuts_XiMinus; //!<! 1D counter for Xi- (Cuts)
        TH1F*                   fCounterRec_Cuts_XiPlus; //!<! 1D counter for Xi+ (Cuts)
        TH2F*                   f2DCounterRecMC_CutsVsPt_Lambda; //!<! 2D counter for Lambda (Pt. vs. Cuts)
        TH2F*                   f2DCounterRecMC_CutsVsPt_AntiLambda; //!<! 2D counter for Anti-Lambda (Pt. vs. Cuts)
        TH2F*                   f2DCounterRecMC_CutsVsPt_XiMinus; //!<! 2D counter for Xi- (Pt. vs. Cuts)
        TH2F*                   f2DCounterRecMC_CutsVsPt_XiPlus; //!<! 2D counter for Xi+ (Pt. vs. Cuts)
        TH2F*                   f2DHistMassPtLambda;  //!<! 2D (Pt vs. Mass) histogram of Lambda
        TH2F*                   f2DHistMassPtAntiLambda;  //!<! 2D (Pt vs. Mass) histogram of Anti-Lambda
        TH2F*                   f2DHistMassPtXiMinus;      //!<! 2D (Pt vs. Mass) histogram of Xi-
        TH2F*                   f2DHistMassPtXiPlus;      //!<! 2D (Pt vs. Mass) histogram of Xi+
        TH2F*                   f2DHistMassPtXicZero_woQtCut; //!<! 2D (Pt vs. Mass) histogram of XicZero
        TH2F*                   f2DHistMassPtXicZero; //!<! 2D (Pt vs. Mass) histogram of XicZero
        TH2F*                   f2DHistMassPtAntiXicZero; //!<! 2D (Pt vs. Mass) histogram of Anti-XicZero
        TH2F*                   f2DHistMassPtXicZeroTot; //!<! 2D (Pt vs. Mass) histogram of XicZero + Anti-XicZero

        TH2F*                   f2DHistChi2vsNDF_Lambda; //!<! 2D (Chi2 vs. NDF) histogram of Lambda
        TH2F*                   f2DHistChi2vsNDF_Lambda_Match; //!<! 2D (Chi2 vs. NDF) histogram of Lambda
        TH2F*                   f2DHistChi2vsNDF_AntiLambda; //!<! 2D (Chi2 vs. NDF) histogram of AntiAnti--Lambda
        TH2F*                   f2DHistXKFMCvsPt_Proton; //!<! 2D ((x_KF-x_MC) vs Pt) histogram of Proton
        TH2F*                   f2DHistYKFMCvsPt_Proton; //!<! 2D ((y_KF-y_MC) vs Pt) histogram of Proton
        TH2F*                   f2DHistZKFMCvsPt_Proton; //!<! 2D ((z_KF-z_MC) vs Pt) histogram of Proton
        TH2F*                   f2DHistXKFMCvsPt_Pion; //!<! 2D ((x_KF-x_MC) vs Pt) histogram of Pion
        TH2F*                   f2DHistYKFMCvsPt_Pion; //!<! 2D ((y_KF-y_MC) vs Pt) histogram of Pion
        TH2F*                   f2DHistZKFMCvsPt_Pion; //!<! 2D ((z_KF-z_MC) vs Pt) histogram of Pion
        TH2F*                   f2DHistXPULLvsPt_Proton; //!<! 2D (x_PLLL vs Pt) histogram of Proton
        TH2F*                   f2DHistYPULLvsPt_Proton; //!<! 2D (y_PLLL vs Pt) histogram of Proton
        TH2F*                   f2DHistZPULLvsPt_Proton; //!<! 2D (z_PLLL vs Pt) histogram of Proton
        TH2F*                   f2DHistXPULLvsPt_Pion; //!<! 2D (x_PLLL vs Pt) histogram of Pion
        TH2F*                   f2DHistYPULLvsPt_Pion; //!<! 2D (y_PLLL vs Pt) histogram of Pion
        TH2F*                   f2DHistZPULLvsPt_Pion; //!<! 2D (z_PLLL vs Pt) histogram of Pion
        TH2F*                   f2DHistXRecMCvsPt_Lambda; //!<! 2D ((x_KF-x_MC) vs Pt) histogram of Lambda
        TH2F*                   f2DHistXRecMCvsPt_AntiLambda; //!<! 2D (x_KF-x_MC vs Pt) histogram of Anti-Lambda
        TH2F*                   f2DHistXV0MCvsPt_Lambda; //!<! 2D ((x_V0-x_MC) vs Pt) histogram of Lambda
        TH2F*                   f2DHistPtRecMCvsPt_Lambda; //!<! 2D ((Pt_Rec-Pt_MC)/Pt_Rec vs Pt) histogram of Lambda
        TH2F*                   f2DHistPtV0MCvsPt_Lambda; //!<! 2D ((Pt_V0-Pt_MC)/Pt_V0 vs Pt) histogram of Lambda
        TH2F*                   f2DHistPtRecMCvsPt_AntiLambda; //!<! 2D ((Pt_Rec-Pt_MC)/Pt_Rec vs Pt) histogram of Anti-Lambda
        TH2F*                   f2DHistMassRecMCvsPt_Lambda; //!<! 2D (m_Rec-m_MC vs Pt) histogram of Lambda
        TH2F*                   f2DHistMassV0MCvsPt_Lambda; //!<! 2D (m_V0-m_MC vs Pt) histogram of Lambda
        TH2F*                   f2DHistMassRecMCvsPt_AntiLambda; //!<! 2D (m_Rec-m_MC vs Pt) histogram of Anti-Lambda
        TH2F*                   f2DHistXPULLvsPt_Lambda; //!<! 2D (x_PLLL vs Pt) histogram of Lambda
        TH2F*                   f2DHistXPULLvsPt_Lambda_V0; //!<! 2D (x_PLLL vs Pt) histogram of Lambda
        TH2F*                   f2DHistXPULLvsPt_AntiLambda; //!<! 2D (x_PLLL vs Pt) histogram of Anti-Lambda
        TH2F*                   f2DHistPtPULLvsPt_Lambda; //!<! 2D (Pt_PLLL vs Pt) histogram of Lambda
        TH2F*                   f2DHistPtPULLvsPt_Lambda_V0; //!<! 2D (Pt_PLLL vs Pt) histogram of Lambda
        TH2F*                   f2DHistPtPULLvsPt_AntiLambda; //!<! 2D (Pt_PLLL vs Pt) histogram of Anti-Lambda
        TH2F*                   f2DHistMassPULLvsPt_Lambda; //!<! 2D (m_PLLL vs Pt) histogram of Lambda
        TH2F*                   f2DHistMassPULLvsPt_Lambda_V0; //!<! 2D (m_PLLL vs Pt) histogram of Lambda
        TH2F*                   f2DHistMassPULLvsPt_AntiLambda; //!<! 2D (m_PLLL vs Pt) histogram of Anti-Lambda
        TH2F*                   f2DHistMassPULLvsRadius_Lambda; //!<! 2D (m_Rec-m_MC vs Radius) histogram of Lambda
        TH2F*                   f2DHistMassPULLvsRadius_AntiLambda; //!<! 2D (m_Rec-m_MC vs Radius) histogram of Anti-Lambda
        TH2F*                   f2DHistArmenterosPodolanski_FirstDaugPos; //!<! Histogram of Armenteros-Podolanski of V0 with positive first daughter
        TH2F*                   f2DHistArmenterosPodolanski_FirstDaugNeg; //!<! Histogram of Armenteros-Podolanski of V0 with negative first daughter
        TH2F*                   f2DHistArmenterosPodolanski_candidate; //!<! Histogram of Armenteros-Podolanski of V0 with positive first daughter after selection
        TH2F*                   f2DHistArmenterosPodolanski_Lam; //!<! Histogram of Armenteros-Podolanski of V0 with positive first daughter after selection
        TH2F*                   f2DHistArmenterosPodolanski_AntiLam; //!<! Histogram of Armenteros-Podolanski of V0 with positive first daughter after selection
        TH2F*                   f2DHistChargeDaughters; //!<! 2D Histogram of daughters charges
        TH2F*                   f2DHistV0XY_OnFly; //!<! 2D Histogram of y vs. x of V0 (on-the-fly)
        TH2F*                   f2DHistV0XY_Offline; //!<! 2D Histogram of y vs. x of V0 (offline)
        TH2F*                   f2DHistV0XY_FirstDaugPos; //!<! 2D Histogram of y vs. x of V0 with positive first daughter 
        TH2F*                   f2DHistV0XY_FirstDaugNeg; //!<! 2D Histogram of y vs. x of V0 with negative first daughter
        TH2F*                   f2DHistLambdaXY; //!<! 2D Histogram of y vs. x of Lambda at decay vertex
        TH2F*                   f2DHistXiMinusXY_DV; //!<! 2D Histogram of y vs. x of Xi- at decay vertex
        TH2F*                   f2DHistXiMinusXY_PV; //!<! 2D Histogram of y vs. x of Xi- at production vertex
        TH2F*                   f2DHistXiPlusXY_DV; //!<! 2D Histogram of y vs. x of Xi+ at decay vertex
        TH2F*                   f2DHistXiPlusXY_PV; //!<! 2D Histogram of y vs. x of Xi+ at production vertex

        TH1F*                   fHistEvents;          //!<! Histogram of selected events
        TH1F*                   fHTrigger;            //!<! Histogram of trigger
        TH1F*                   fHistOnFlyStatus; //!<! =1, fOnFlyStatus=kTRUE, this V0 is reconstructed "on fly" during the tracking; =-1, fOnFlyStatus=kFALSE
        TH1F*                   fHistOnFlyStatus_FirstDaugPos; //!<! =1, fOnFlyStatus=kTRUE, this V0 with positive first daughter is reconstructed "on fly" during the tracking; =-1, fOnFlyStatus=kFALSE. 
        TH1F*                   fHistOnFlyStatus_FirstDaugNeg; //!<! =1, fOnFlyStatus=kTRUE, this V0 with negative first daughter is reconstructed "on fly" during the tracking; =-1, fOnFlyStatus=kFALSE
        TH1F*                   fHistChargeV0; //!<! Histogram of V0 charge
        TH1F*                   fHistNProngV0; //!<! Histogram of number of V0 prong
        TH1F*                   fHistNDaughterV0; //!<! Histogram of Duaghter number of V0
        TH1F*                   fHistChargeFirstDaughter; //!<! Histogram of first daughter charge
        TH1F*                   fHistChargeSecondDaughter; //!<! Histogram of second daughter charge
        TH1F*                   fHistXtrkP;         //!<! Histogram of x of positive track
        TH1F*                   fHistYtrkP;         //!<! Histogram of y of positive track
        TH1F*                   fHistZtrkP;         //!<! Histogram of z of positive track
        TH1F*                   fHistXtrkP_XYZv;    //!<! Histogram of x of positive track at SecVtx
        TH1F*                   fHistYtrkP_XYZv;    //!<! Histogram of y of positive track at SecVtx
        TH1F*                   fHistZtrkP_XYZv;    //!<! Histogram of z of positive track at SecVtx
        TH1F*                   fHistXtrkP_Rec_MC;  //!<! Histogram of x(Rec-MC) of positive track
        TH1F*                   fHistYtrkP_Rec_MC;  //!<! Histogram of y(Rec-MC) of positive track
        TH1F*                   fHistZtrkP_Rec_MC;  //!<! Histogram of z(Rec-MC) of positive track
        TH1F*                   fHistXtrkP_Rec_MC_XYZv;  //!<! Histogram of x(Rec-MC) of positive track
        TH1F*                   fHistYtrkP_Rec_MC_XYZv;  //!<! Histogram of y(Rec-MC) of positive track
        TH1F*                   fHistZtrkP_Rec_MC_XYZv;  //!<! Histogram of z(Rec-MC) of positive track
        TH1F*                   fHistXtrkN;           //!<! Histogram of x of negative track
        TH1F*                   fHistYtrkN;           //!<! Histogram of Y of negative track
        TH1F*                   fHistZtrkN;           //!<! Histogram of z of negative track
        TH1F*                   fHistXtrkN_XYZv;    //!<! Histogram of x of negative track at SecVtx
        TH1F*                   fHistYtrkN_XYZv;    //!<! Histogram of y of negative track at SecVtx
        TH1F*                   fHistZtrkN_XYZv;    //!<! Histogram of z of negative track at SecVtx
        TH1F*                   fHistXtrkN_Rec_MC;  //!<! Histogram of x(Rec-MC) of negative track
        TH1F*                   fHistYtrkN_Rec_MC;  //!<! Histogram of y(Rec-MC) of negative track
        TH1F*                   fHistZtrkN_Rec_MC;  //!<! Histogram of z(Rec-MC) of negative track
        TH1F*                   fHistXtrkN_Rec_MC_XYZv;  //!<! Histogram of x(Rec-MC) of negative track
        TH1F*                   fHistYtrkN_Rec_MC_XYZv;  //!<! Histogram of y(Rec-MC) of negative track
        TH1F*                   fHistZtrkN_Rec_MC_XYZv;  //!<! Histogram of z(Rec-MC) of negative track
        TH1F*                   fHistLDeltaLRec_Lambda;      //!<! Histogram of l/DeltaL of Lambda
        TH1F*                   fHistLDeltaLRecMC_Lambda;      //!<! Histogram of l/DeltaL of Lambda
        TH1F*                   fHistLDeltaLRecMC_LambdaFromXi;      //!<! Histogram of l/DeltaL of Lambda
        TH1F*                   fHistLDeltaLRec_AntiLambda;      //!<! Histogram of l/DeltaL of Anti-Lambda
        TH1F*                   fHistLDeltaLRecMC_AntiLambda;      //!<! Histogram of l/DeltaL of Anti-Lambda
        TH1F*                   fHistLDeltaLRecMC_AntiLambdaFromXi;      //!<! Histogram of l/DeltaL of Anti-Lambda
        TH1F*                   fHistXLambdaTot;      //!<! Histogram of x of Lambda + Anti-Lambda
        TH1F*                   fHistYLambdaTot;      //!<! Histogram of y of Lambda + Anti-Lambda
        TH1F*                   fHistZLambdaTot;      //!<! Histogram of z of Lambda + Anti-Lambda
        TH1F*                   fHistXLambda_KF_MC;  //!<! Histogram of x(KF-MC) of Lambda
        TH1F*                   fHistYLambda_KF_MC;  //!<! Histogram of y(KF-MC) of Lambda
        TH1F*                   fHistZLambda_KF_MC;  //!<! Histogram of z(KF-MC) of Lambda
        TH1F*                   fHistXProton_KF_MC;  //!<! Histogram of x(KF-MC) of Proton
        TH1F*                   fHistYProton_KF_MC;  //!<! Histogram of y(KF-MC) of Proton
        TH1F*                   fHistZProton_KF_MC;  //!<! Histogram of z(KF-MC) of Proton
        TH1F*                   fHistXPion_KF_MC;  //!<! Histogram of x(KF-MC) of Pion
        TH1F*                   fHistYPion_KF_MC;  //!<! Histogram of y(KF-MC) of Pion
        TH1F*                   fHistZPion_KF_MC;  //!<! Histogram of z(KF-MC) of Pion
        TH1F*                   fHistXLambda_V0_MC;  //!<! Histogram of x(V0-MC) of Lambda
        TH1F*                   fHistXAntiLambda_Rec_MC;  //!<! Histogram of x(Rec-MC) of Anti-Lambda
        TH1F*                   fHistYAntiLambda_Rec_MC;  //!<! Histogram of y(Rec-MC) of Anti-Lambda
        TH1F*                   fHistZAntiLambda_Rec_MC;  //!<! Histogram of z(Rec-MC) of Anti-Lambda
        TH1F*                   fHistXLambda_PULL;  //!<! Histogram of x(PULL) of Lambda
        TH1F*                   fHistYLambda_PULL;  //!<! Histogram of y(PULL) of Lambda
        TH1F*                   fHistZLambda_PULL;  //!<! Histogram of z(PULL) of Lambda
        TH1F*                   fHistXProton_PULL;  //!<! Histogram of x(PULL) of Proton
        TH1F*                   fHistYProton_PULL;  //!<! Histogram of y(PULL) of Proton
        TH1F*                   fHistZProton_PULL;  //!<! Histogram of z(PULL) of Proton
        TH1F*                   fHistXPion_PULL;  //!<! Histogram of x(PULL) of Pion
        TH1F*                   fHistYPion_PULL;  //!<! Histogram of y(PULL) of Pion
        TH1F*                   fHistZPion_PULL;  //!<! Histogram of z(PULL) of Pion
        TH1F*                   fHistXAntiLambda_PULL;  //!<! Histogram of x(PULL) of Anti-Lambda
        TH1F*                   fHistYAntiLambda_PULL;  //!<! Histogram of y(PULL) of Anti-Lambda
        TH1F*                   fHistZAntiLambda_PULL;  //!<! Histogram of z(PULL) of Anti-Lambda
        TH1F*                   fHistXXiTot;           //!<! Histogram of x of XiMinus + XiPlus
        TH1F*                   fHistYXiTot;           //!<! Histogram of y of XiMinus + XiPlus
        TH1F*                   fHistZXiTot;           //!<! Histogram of z of XiMinus + XiPlus
        TH1F*                   fHistXXicZeroTot;     //!<! Histogram of x of XicZero + Anti-XicZero
        TH1F*                   fHistYXicZeroTot;     //!<! Histogram of y of XicZero + Anti-XicZero
        TH1F*                   fHistZXicZeroTot;     //!<! Histogram of z of XicZero + Anti-XicZero
        TH1F*                   fGenHistRapidity_Lambda; //!<! Histogram of rapidity of Lambda at Gen. level
        TH1F*                   fGenHistRapidity_AntiLambda; //!<! Histogram of rapidity of AntiLambda at Gen. level
        TH1F*                   fRecHistRapidity_Lambda_offline; //!<! Histogram of rapidity of offline Lambda
        TH1F*                   fRecHistRapidity_Lambda_wSTD; //!<! Histogram of rapidity of Lambda after STD cut
        TH1F*                   fRecHistRapidity_AntiLambda_wSTD; //!<! Histogram of rapidity of AntiLambda after STD cut
        TH1F*                   fHistPtLambda;        //!<! Pt histogram of Lambda
        TH1F*                   fHistPtAntiLambda;    //!<! Pt histogram of Anti-Lambda
        TH1F*                   fHistPtLambdaTot;     //!<! Pt histogram of Lambda + Anti-Lambda
        TH1F*                   fHistPtXiMinus;       //!<! Pt histogram of XiMinus
        TH1F*                   fHistPtXiPlus;        //!<! Pt histogram of XiPlus
        TH1F*                   fHistPtXiTot;        //!<! Pt histogram of XiMinus + XiPlus
        TH1F*                   fHistPtXicZero;       //!<! Pt histogram of XicZero
        TH1F*                   fHistPtAntiXicZero;   //!<! Pt histogram of Anti-XicZero
        TH1F*                   fHistPtXicZeroTot;    //!<! Pt histogram of XicZero + Anti-XicZero
        TH1F*                   fHistMassK0S;         //!<! Mass of K0S
        TH1F*                   fHistMassLambda_woCut; //!<! Mass histogram of Lambda without any cut
        TH1F*                   fHistMassAntiLambda_woCut; //!<! Mass histogram of Anti-Lambda without any cut
        TH1F*                   fHistMassLambdaCheck;      //!<! Mass histogram of Lambda (check v0->MassLambda())
        TH1F*                   fHistMassAntiLambdaCheck;      //!<! Mass histogram of Anti-Lambda (check v0->MassAntiLambda())
        TH1F*                   fHistMassLambda_wSTDv0Cut; //!<!
        TH1F*                   fHistMassAntiLambda_wSTDv0Cut; //!<!
        TH1F*                   fHistMassLambda_BeforeSecSel; //!<!
        TH1F*                   fHistMassAntiLambda_BeforeSecSel; //!<!
        TH1F*                   fHistMassLambda;      //!<! Mass histogram of Lambda
        TH1F*                   fHistMassLambda_woArmenterosPodolanskiCut; //!<! Mass histogram of Lambda without ArmenterosPodolanski cut
        TH1F*                   fHistMassLambda_woMassCut;      //!<! Mass histogram of Lambda before mass window cut
        TH1F*                   fHistMassAntiLambda;      //!<! Mass histogram of Anti-Lambda
        TH1F*                   fHistMassAntiLambda_woMassCut;  //!<! Mass histogram of Anti-Lambda before mass window cut
        TH1F*                   fHistMassLambdaTot;   //!<! Mass histogram of Lambda + Anti-Lambda
        TH1F*                   fHistMassLambda_Match;  //!<! Mass histogram of Lambda after MatchToMC
        TH1F*                   fHistMassAntiLambda_Match;  //!<! Mass histogram of Anti-Lambda after MatchToMC
        TH1F*                   fHistMassLambdaTot_Match;   //!<! Mass histogram of Lambda + Anti-Lambda after MatchToMC
        TH1F*                   fHistMassLambda_V0;      //!<! Mass histogram of Lambda
        TH1F*                   fHistMassAntiLambda_V0;  //!<! Mass histogram of Anti-Lambda
        TH1F*                   fHistMassLambdaTot_V0;   //!<! Mass histogram of Lambda + Anti-Lambda
        TH1F*                   fHistMassLambda_KF_V0;      //!<! Mass histogram of Lambda
        TH1F*                   fHistMassAntiLambda_KF_V0;  //!<! Mass histogram of Anti-Lambda
        TH1F*                   fHistMassLambda_KF_MC;      //!<! Mass histogram of Lambda
        TH1F*                   fHistMassLambda_V0_MC;      //!<! Mass histogram of Lambda
        TH1F*                   fHistMassAntiLambda_KF_MC;  //!<! Mass histogram of Anti-Lambda
        TH1F*                   fHistMassLambda_PULL_KF;  //!<! Mass histogram of Anti-Lambda
        TH1F*                   fHistMassAntiLambda_PULL_KF;  //!<! Mass histogram of Anti-Lambda
        TH1F*                   fHistMassLambda_M;      //!<! Mass histogram of Lambda with mass constraint
        TH1F*                   fHistMassAntiLambda_M;  //!<! Mass histogram of Anti-Lambda with mass constraint
        TH1F*                   fHistMassLambdaTot_M;   //!<! Mass histogram of Lambda + Anti-Lambda with mass constraint
        TH1F*                   fHistMassLambda_MV;      //!<! Mass histogram of Lambda with mass and vertex constraint
        TH1F*                   fHistMassAntiLambda_MV;  //!<! Mass histogram of Anti-Lambda with mass and vertex constraint
        TH1F*                   fHistMassLambdaTot_MV;   //!<! Mass histogram of Lambda + Anti-Lambda with mass and vertex constraint
        TH1F*                   fHistMassXiMinus;     //!<! Mass histogram of XiMinus
        TH1F*                   fHistMassXiMinus_M;     //!<! Mass histogram of XiMinus after mass constraint
        TH1F*                   fHistMassXiMinus_Match;     //!<! Mass histogram of XiMinus after Matching to MC
        TH1F*                   fHistMassXiPlus_Match;     //!<! Mass histogram of XiPlus after Matching to MC
        TH1F*                   fHistMassXiPlus;      //!<! Mass histogram of XiPlus
        TH1F*                   fHistMassXiPlus_M;      //!<! Mass histogram of XiPlus after mass constraint
        TH1F*                   fHistMassXiTot;      //!<! Mass histogram of XiMinus + XiPlus
        TH1F*                   fHistMassXicZero_woQtCut;     //!<! Mass histogram of XicZero
        TH1F*                   fHistMassXicZero;     //!<! Mass histogram of XicZero
        TH1F*                   fHistMassAntiXicZero; //!<! Mass histogram of Anti-XicZero
        TH1F*                   fHistMassXicZeroTot;   //!<! Mass histogram of XicZero + Anti-XicZero
        TH1F*                   fHistQtDiffPionXiMinus; //!<!
        TH1F*                   fHistQtDiffPionXiPlus; //!<!
        TH1F*                   fHistMassXiMinus2;      //!<! Mass histogram of XiMinus2
        TH1F*                   fHistChi2ndfProton;   //!<! Chi2/NDF of Proton
        TH1F*                   fHistChi2ndfPion;   //!<! Chi2/NDF of Pion
        TH1F*                   fHistChi2ndfLambda;   //!<! Chi2/NDF of Lambda
        TH1F*                   fHistChi2ndfLambda_Match;   //!<! Chi2/NDF of Lambda
        TH1F*                   fHistChi2ndfAntiLambda;   //!<! Chi2/NDF of Anti-Lambda
        TH1F*                   fHistChi2ndfAntiLambda_Match;   //!<! Chi2/NDF of Lambda
        TH1F*                   fHistChi2ndfLambdaTot;   //!<! Chi2/NDF of Lambda + Anti-Lambda
        TH1F*                   fHistProbProton; //!<! Probability of Proton
        TH1F*                   fHistProbPion; //!<! Probability of Pion
        TH1F*                   fHistProbLambda; //!<! Probability of Lambda
        TH1F*                   fHistProbLambda_chi2cut; //!<! Probability of Lambda
        TH1F*                   fHistProbLambda_Match; //!<! Probability of Lambda
        TH1F*                   fHistProbAntiLambda; //!<! Probability of Anti-Lambda
        TH1F*                   fHistProbAntiLambda_chi2cut; //!<! Probability of Anti-Lambda
        TH1F*                   fHistProbAntiLambda_Match; //!<! Probability of Anti-Lambda
        TH1F*                   fHistProbLambdaTot; //!<! Probability of Lambda + Anti-Lambda
        TH1F*                   fHistProbXiMinus; //!<! Probability of XiMinus
        TH1F*                   fHistProbXiMinus_chi2cut; //!<! Probability of XiMinus
        TH1F*                   fHistProbXiPlus; //!<! Probability of XiPlus
        TH1F*                   fHistProbXiPlus_chi2cut; //!<! Probability of XiPlus
        TH1F*                   fHistProbXicZero; //!<! Probability of XicZero
        TH1F*                   fHistProbXicZero_chi2cut; //!<! Probability of XicZero
        TH1F*                   fHistProbAntiXicZero; //!<! Probability of Anti-XicZero
        TH1F*                   fHistProbAntiXicZero_chi2cut; //!<! Probability of Anti-XicZero
        TH1F*                   fHistChi2ndfXiMinus;   //!<! Chi2/NDF of Xi-
        TH1F*                   fHistChi2ndfXiPlus;   //!<! Chi2/NDF of Xi+
        TH1F*                   fHistChi2ndfXiTot;   //!<! Chi2/NDF of Xi- + Xi+
        TH1F*                   fHistChi2ndfXicZero;   //!<! Chi2/NDF of Xic0
        TH1F*                   fHistChi2ndfAntiXicZero;   //!<! Chi2/NDF of Anti-Xic0
        TH1F*                   fHistChi2ndfXicZeroTot;   //!<! Chi2/NDF of Xic0 + Anti-Xic0
        TH1F*                   fHistDecayLLambda; //!<! Decay length histogram of Lambda
        TH1F*                   fHistDecayLAntiLambda; //!<! Decay length histogram of Anti-Lambda
        TH1F*                   fHistDecayLLambdaTot; //!<! Decay length histogram of Lambda + Anti-Lambda
        TH1F*                   fHistDecayLXiMinus; //!<! Decay length histogram of Xi-
        TH1F*                   fHistDecayLXiPlus; //!<! Decay length histogram of Xi+
        TH1F*                   fHistDecayLXiTot; //!<! Decay length histogram of Xi- + Xi+
        TH1F*                   fHistDecayLXicZero; //!<! Decay length histogram of Xic0
        TH1F*                   fHistDecayLAntiXicZero; //!<! Decay length histogram of Anti-Xic0
        TH1F*                   fHistDecayLXicZeroTot; //!<! Decay length histogram of Xic0 + Anti-Xic0
        TH1F*                   fHistCosPA_Lambda;      //!<! Histogram of cos(PA) of Lambda
        TH1F*                   fHistCosPA_AntiLambda; //!<! Histogram of cos(PA) of Anti-Lambda
        TH1F*                   fHistCosPA_LambdaTot;  //!<! Histogram of cos(PA) of Lambda + Anti-Lambda
        TH1F*                   fHistCosPA_XiMinus;    //!<! Histogram of cos(PA) of Xi-
        TH1F*                   fHistCosPA_XiPlus;    //!<! Histogram of cos(PA) of Xi+
        TH1F*                   fHistCosPA_XiTot;    //!<! Histogram of cos(PA) of Xi- + Xi+
        TH1F*                   fHistPVx;              //!<! Histogram of primary vertex in x
        TH1F*                   fHistPVy;              //!<! Histogram of primary vertex in y
        TH1F*                   fHistPVz;              //!<! Histogram of primary vertex in z
        TH1F*                   fHCentrality;          //!<! Histogram of centrality
        TH1F*                   fHistMCXicZeroDecayType; //!<! MC event type of Xic0
        TH1F*                   fHistMCXiDecayType; //!<! MC event type of Xi
        TH1F*                   fHistMCpdg_All;     //!<! PDG of all particle
        TH1F*                   fHistMCpdg_Dau_XicZero;     //!<! PDG of all particle from Xic0 decay
        TH1F*                   fHistMCpdg_Dau_XicPM;     //!<! PDG of all particle from Xic+- decay
        THnSparseF*             fHistMCGen_XicZeroTot;  //!<! mcArray
        THnSparseF*             fHistMCGen_XicZero;     //!<! mcArray
        THnSparseF*             fHistMCGen_AntiXicZero; //!<! mcArray
        THnSparseF*             fHistMCGen_PionTot; //!<! mcArray
        THnSparseF*             fHistMCGen_PionPlus; //!<! mcArray
        THnSparseF*             fHistMCGen_PionMinus; //!<! mcArray
        THnSparseF*             fHistMCGen_XiTot; //!<! mcArray
        THnSparseF*             fHistMCGen_XiMinus; //!<! mcArray
        THnSparseF*             fHistMCGen_XiPlus; //!<! mcArray
        THnSparseF*             fHistMCGen_Lambda; //!<! mcArray
        THnSparseF*             fHistMCGen_AntiLambda; //!<! mcArray
        THnSparseF*             fHistMCGen_PiXiInvMass; //!<! mcArray
        THnSparseF*             fHistMCGen_PiXiMassvsPiPt; //!<! mcArray
        THnSparseF*             fHistMCGen_PiXiMassvsPiPt_PionPlus; //!<! mcArray
        THnSparseF*             fHistMCGen_PiXiMassvsPiPt_PionMinus; //!<! mcArray
        Bool_t                  fWriteXic0MCGenTree; ///< flag to decide whether to write the MC candidate variables on a tree variables
        Bool_t                  fWriteXic0Tree; ///< flag to decide whether to write XicZero tree



        AliAnalysisTaskSEXicZero2XiPifromKFP(const AliAnalysisTaskSEXicZero2XiPifromKFP &source); // not implemented
        AliAnalysisTaskSEXicZero2XiPifromKFP& operator=(const AliAnalysisTaskSEXicZero2XiPifromKFP& source); // not implemented

        ClassDef(AliAnalysisTaskSEXicZero2XiPifromKFP, 2);
};

#endif
