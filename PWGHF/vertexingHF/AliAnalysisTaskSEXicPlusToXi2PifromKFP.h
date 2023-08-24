#ifndef ALIANALYSISTASKSEXICPLUSTOXI2PIFROMKFP_H
#define ALIANALYSISTASKSEXICPLUSTOXI2PIFROMKFP_H

/// \class AliAnalysisTaskSEXicPlusToXi2PifromKFP
/// \brief This is a brief description of my class
///
/// This is a longer description of the class. This longer description is
/// formatted using the Markdown syntax (see below) and can span on multiple
/// lines.
///
/// \author Jianhui Zhu <zjh@ccnu.edu.cn>, Central China Normal University & GSI Helmholtz Centre for Heavy Ion Research
/// \date Jun 1, 2021

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
#include "AliVertexingHFUtils.h"
#include "AliVertexerTracks.h"

// includes added to play with KFParticle
#include <vector>
#include "KFParticleBase.h"
#include "KFParticle.h"
#include "KFPTrack.h"
#include "KFPVertex.h"
#include "KFVertex.h"

class AliAnalysisTaskSEXicPlusToXi2PifromKFP : public AliAnalysisTaskSE
{
    public:
                                AliAnalysisTaskSEXicPlusToXi2PifromKFP();
                                AliAnalysisTaskSEXicPlusToXi2PifromKFP(const char *name, AliRDHFCutsKFP* cuts);
        virtual                 ~AliAnalysisTaskSEXicPlusToXi2PifromKFP();

        virtual void            UserCreateOutputObjects();
        virtual void            Init();
        virtual void            LocalInit() {Init();} 
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);

        void                    SetMC(Bool_t IsMC) {fIsMC=IsMC;}
        void                    SelectTrack(AliVEvent *event, Int_t trkEntries, Int_t &nSeleTrks, Bool_t *seleFlags);
        Bool_t                  MakeMCAnalysis(TClonesArray *mcArray);
        void                    MakeAnaXicPlusFromCasc(AliAODEvent *AODEvent, TClonesArray *mcArray, KFParticle PV);
        Int_t                   MatchToMCXicPlus(AliAODTrack *trackProton, AliAODTrack *trackPionMinus_2, AliAODTrack *trackPionMinus_1, AliAODTrack *trackPionPlus_trk1, AliAODTrack *trackPionPlus_trk2, TClonesArray *mcArray);
        Int_t                   MatchToMCAntiXicPlus(AliAODTrack *trackAntiProton, AliAODTrack *trackPionPlus_2, AliAODTrack *trackPionPlus_1, AliAODTrack *trackPionMinus_trk1, AliAODTrack *trackPionMinus_trk2, TClonesArray *mcArray);
        Int_t                   MatchToMCXiMinus(AliAODTrack *trackProton, AliAODTrack *trackPion3, AliAODTrack *trackPion2, TClonesArray *mcArray);
        Int_t                   MatchToMCXiPlus(AliAODTrack *trackAntiProton, AliAODTrack *trackAntiPion3, AliAODTrack *trackAntiPion2, TClonesArray *mcArray);
        Int_t                   MatchToMCLambda(AliAODTrack *trackProton, AliAODTrack *trackPion3, TClonesArray *mcArray);
        Int_t                   MatchToMCAntiLambda(AliAODTrack *trackAntiProton, AliAODTrack *trackAntiPion3, TClonesArray *mcArray);
        Int_t                   MatchToMCLambdaFromXi(AliAODTrack *trackProton, AliAODTrack *trackPion3, TClonesArray *mcArray);
        Int_t                   MatchToMCAntiLambdaFromXi(AliAODTrack *trackAntiProton, AliAODTrack *trackAntiPion3, TClonesArray *mcArray);
        Int_t                   MatchToMCPion(AliAODTrack *track, TClonesArray *mcArray);
        Double_t                InvMassV0atPV(AliAODTrack *trk1, AliAODTrack *trk2, Int_t pdg1, Int_t pdg2);
        ULong64_t               GetEventIdAsLong(AliVHeader* header);

        /// set MC usage
        void SetWriteXicPlusMCGenTree(Bool_t a) {fWriteXicPlusMCGenTree = a;}
        Bool_t GetWriteXicPlusMCGenTree() const {return fWriteXicPlusMCGenTree;}

        void SetWriteXicPlusTree(Bool_t a) {fWriteXicPlusTree = a;}
        Bool_t GetWriteXicPlusTree() const {return fWriteXicPlusTree;}

        void SetWriteXicPlusQATree(Bool_t a) {fWriteXicPlusQATree = a;}
        Bool_t GetWriteXicPlusQATree() const {return fWriteXicPlusQATree;}

        void FillEventROOTObjects();
        void FillTreeGenXicPlus(AliAODMCParticle *mcpart, Int_t CheckOrigin, Double_t MLoverP);
        void FillTreeRecXicPlusFromCasc(AliAODEvent *AODEvent, AliAODcascade *casc, KFParticle kfpXicPlus, AliAODTrack *trackPiFromXicPlus_trk1, KFParticle kfpBP_trk1, KFParticle kfpXiMinus, KFParticle kfpXiMinus_m, KFParticle kfpPionOrKaon, AliAODTrack *trackPiFromXiOrKaonFromOmega, KFParticle kfpK0Short, KFParticle kfpGamma, KFParticle kfpLambda, KFParticle kfpLambda_m, AliAODTrack *trkProton, AliAODTrack *trkPion, AliAODTrack *trackPiFromXicPlus_trk2, KFParticle kfpBP_trk2, KFParticle kfpProtonFromLam, KFParticle kfpPionFromLam, KFParticle PV, KFParticle PV_KF_Refit, TClonesArray *mcArray, Int_t lab_XicPlus);
        void FillQATreeXicPlusFromCasc_woMassConstForLamAndXi(KFParticle kfpLambda, KFParticle kfpXiMinus_woMassConstForLamAndXi, KFParticle kfpXicPlus_woMassConstForLamAndXi, KFParticle PV, AliAODTrack *trackPiFromXicPlus_HighPt, TClonesArray *mcArray, Int_t lab_XicPlus, AliAODEvent *AODEvent);
        void FillQATreeXicPlusFromCasc_woMassConstForLam_wMassConstForXi(KFParticle kfpLambda, KFParticle kfpXiMinus_woMassConstForLam_wMassConstForXi, KFParticle kfpXicPlus_woMassConstForLam_wMassConstForXi, KFParticle PV, AliAODTrack *trackPiFromXicPlus_HighPt, TClonesArray *mcArray, Int_t lab_XicPlus, AliAODEvent *AODEvent);
        void FillQATreeXicPlusFromCasc_wMassConstForLam_woMassConstForXi(KFParticle kfpLambda_wMassConst, KFParticle kfpXiMinus_wMassConstForLam_woMassConstForXi, KFParticle kfpXicPlus_wMassConstForLam_woMassConstForXi, KFParticle PV, AliAODTrack *trackPiFromXicPlus_HighPt, TClonesArray *mcArray, Int_t lab_XicPlus, AliAODEvent *AODEvent);
        void FillQATreeXicPlusFromCasc_wMassConstForLamAndXi(KFParticle kfpLambda_wMassConst, KFParticle kfpXiMinus_wMassConstForLamAndXi, KFParticle kfpXicPlus_wMassConstForLamAndXi, KFParticle PV, AliAODTrack *trackPiFromXicPlus_HighPt, TClonesArray *mcArray, Int_t lab_XicPlus, AliAODEvent *AODEvent);
        void FillQATreeXicPlusFromCasc_wMassAndTopoConstForLam_wMassConstForXi(KFParticle kfpLambda_wMassConst_To_Xi_wMassConst, KFParticle kfpXiMinus_wMassAndTopoConstForLam_wMassConstForXi, KFParticle kfpXicPlus_wMassAndTopoConstForLam_wMassConstForXi, KFParticle PV, AliAODTrack *trackPiFromXicPlus_HighPt, TClonesArray *mcArray, Int_t lab_XicPlus, AliAODEvent *AODEvent);
        void FillQATreeXicPlusFromCasc_wMassAndTopoConstForLam_wMassAndTopoConstForXi(KFParticle kfpLambda_wMassConst_To_Xi_wMassConst, KFParticle kfpXiMinus_wMassAndTopoConstForLam_wMassAndTopoConstForXi, KFParticle kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi, KFParticle PV, AliAODTrack *trackPiFromXicPlus_HighPt, TClonesArray *mcArray, Int_t lab_XicPlus, AliAODEvent *AODEvent);
        void FillQATreeXicPlusFromCasc_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic(KFParticle kfpLambda_wMassConst_To_Xi_wMassConst, KFParticle kfpXiMinus_wMassAndTopoConstForLam_wMassAndTopoConstForXi, KFParticle kfpXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic, KFParticle PV, AliAODTrack *trackPiFromXicPlus_HighPt, AliAODTrack *trackPiFromXicPlus_LowPt, AliAODTrack *trackPiFromXiOrKaonFromOmega, AliAODTrack *trackPrFromLam, AliAODTrack *trackPiFromLam, KFParticle kfpBP_HighPt, KFParticle kfpBP_LowPt, KFParticle kfpPion_ForXi, KFParticle kfpProton_ForLam, KFParticle kfpPion_ForLam, TClonesArray *mcArray, Int_t lab_XicPlus, AliAODEvent *AODEvent);

        AliAODVertex* PrimaryVertex(const TObjArray *trkArray, AliVEvent *event);
        AliAODVertex* CallPrimaryVertex(AliAODcascade *casc, AliAODTrack *trk1, AliAODTrack *trk2, AliAODEvent *aodEvent);

        unsigned int GetMCEventID();

    private:
        void                    DefineEvent();
        void                    DefineTreeRecXicPlus();
        void                    DefineTreeGenXicPlus();
        void                    DefineAnaHist();
        void                    DefineTreeQAXicPlus();
        void                    DefineTreeQAXicPlus_woMassConstForLamAndXi();
        void                    DefineTreeQAXicPlus_woMassConstForLam_wMassConstForXi();
        void                    DefineTreeQAXicPlus_wMassConstForLam_woMassConstForXi();
        void                    DefineTreeQAXicPlus_wMassConstForLamAndXi();
        void                    DefineTreeQAXicPlus_wMassAndTopoConstForLam_wMassConstForXi();
        void                    DefineTreeQAXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi();
        void                    DefineTreeQAXicPlus_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic();
        AliPIDResponse*         fPID;                 ///<
        AliRDHFCutsKFP*         fAnaCuts;             ///< Cuts
        AliAODVertex*           fpVtx;                //!<! primary vertex
        AliMCEvent*             fMCEvent;             //!<! corresponding mc event
        Double_t                fBzkG;                ///< magnetic field value [kG]
        Float_t                 fCentrality;          ///< Centrality
        vector<Int_t>           fAodTrackInd;         ///< Translation table: aodTrackInd(mcTrackIndex) = aodTrackIndex
        TList*                  fOutputList;          //!<! Output list
        TTree*                  fTree_Event;          //!<! tree of event
        Float_t*                fVar_Event;           //!<! variables of event to be written to the tree
        TTree*                  fTree_XicPlus;        //!<! tree of the candidate variables
        Float_t*                fVar_XicPlus;         //!<! variables of Xic+ to be written to the tree
        TTree*                  fTree_XicPlus_QA;     //!<! QA tree of the candidate variables
        Float_t*                fVar_XicPlus_QA;      //!<! variables of Xic+ to be written to the QA tree
        TTree*                  fTree_XicPlusMCGen;   //!<! tree of the candidate variables after track selection on output slot
        Float_t*                fVar_XicPlusMCGen;    //!<! variables to be written to the tree
        TTree*                  fTree_XicPlus_QA_woMassConstForLamAndXi;     //!<! QA tree without mass constraint for both Lambda and Xi
        Float_t*                fVar_XicPlus_QA_woMassConstForLamAndXi;      //!<! variables of QA tree without mass constraint for both Lambda and Xi
        TTree*                  fTree_XicPlus_QA_wMassConstForLam_woMassConstForXi;     //!<! QA tree with mass constraint for Lambda and without mass constraint for Xi
        Float_t*                fVar_XicPlus_QA_wMassConstForLam_woMassConstForXi;      //!<! variables of QA tree with mass constraint for Lambda and without mass constraint for Xi
        TTree*                  fTree_XicPlus_QA_woMassConstForLam_wMassConstForXi;     //!<! QA tree without mass constraint for Lambda and with mass constraint for Xi
        Float_t*                fVar_XicPlus_QA_woMassConstForLam_wMassConstForXi;      //!<! variables of QA tree without mass constraint for Lambda and with mass constraint for Xi
        TTree*                  fTree_XicPlus_QA_wMassConstForLamAndXi;     //!<! QA tree with mass constraint for both Lambda and Xi
        Float_t*                fVar_XicPlus_QA_wMassConstForLamAndXi;      //!<! variables of QA tree with mass constraint for both Lambda and Xi
        TTree*                  fTree_XicPlus_QA_wMassAndTopoConstForLam_wMassConstForXi;     //!<! QA tree (Lambda with mass constraint and topo to mass constrained Xi)
        Float_t*                fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassConstForXi;      //!<! variables of QA tree
        TTree*                  fTree_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi;     //!<! QA tree (Lambda with mass constraint and topo to mass constrained Xi + Xi with topo constraint to primary)
        Float_t*                fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi;      //!<! variables of QA tree
        TTree*                  fTree_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic;     //!<! QA tree (Lambda with mass constraint and topo to mass constrained Xi + Xi with topo constraint to primary)
        Float_t*                fVar_XicPlus_QA_wMassAndTopoConstForLam_wMassAndTopoConstForXi_wTopoConstForXic;      //!<! variables of QA tree
        TList*                  fListCuts;            //!<! User output slot 3 // Cuts 
        ULong64_t               fVar_XicPlus_EvtID;   //!<! Event ID

        Bool_t                  fIsMC; ///< Flag of MC analysis

        AliNormalizationCounter* fCounter; //!<! Counter for normalization

        TH1F*                   fHistEvents;          //!<! Histogram of selected events
        TH1F*                   fHTrigger;            //!<! Histogram of trigger
        TH1F*                   fHCentrality;          //!<! Histogram of centrality
        TH1F*                   fHCountUsedForPrimVtxFit; //!<! Histogram of frequency of counting AOD track used for primary vertex fit
        TH1F*                   fHNumberOfCasc; //!<! Histogram of frequency of number of cascade
        TH1F*                   fHPrimVtx_woDau_x; //!<! Histogram of PV after removing daughter tracks (x)
        TH1F*                   fHPrimVtx_woDau_y; //!<! Histogram of PV after removing daughter tracks (y)
        TH1F*                   fHPrimVtx_woDau_z; //!<! Histogram of PV after removing daughter tracks (z)
        TH1F*                   fHPrimVtx_woDau_err_x; //!<! Histogram of PV after removing daughter tracks (err_x)
        TH1F*                   fHPrimVtx_woDau_err_y; //!<! Histogram of PV after removing daughter tracks (err_y)
        TH1F*                   fHPrimVtx_woDau_err_z; //!<! Histogram of PV after removing daughter tracks (err_z)
        TH1F*                   fHNumOfCandidatePerEvent_In3sigma; //!<! Histogram of number of Xic+ candidate per event within 3 sigma (assuming sigma=0.01)
        TH1F*                   fHPrimVtx_PV_KF_Refit_Minus_PVrec_x; //!<! Histogram of difference between PV(rec) and PV(rec) from KF (x)
        TH1F*                   fHPrimVtx_PV_KF_Refit_Minus_PVrec_y; //!<! Histogram of difference between PV(rec) and PV(rec) from KF (y)
        TH1F*                   fHPrimVtx_PV_KF_Refit_Minus_PVrec_z; //!<! Histogram of difference between PV(rec) and PV(rec) from KF (z)
        TH1F*                   fHPrimVtx_recalPV_Minus_PVrec_x; //!<! Histogram of difference between recalPV and PV(rec) (x)
        TH1F*                   fHPrimVtx_recalPV_Minus_PVrec_y; //!<! Histogram of difference between recalPV and PV(rec) (y)
        TH1F*                   fHPrimVtx_recalPV_Minus_PVrec_z; //!<! Histogram of difference between recalPV and PV(rec) (z)
        TH1F*                   fHPrimVtx_recalPV_KF_Refit_Minus_PV_KF_Refit_x; //!<! Histogram of difference between recalPV_KF and PV_KF (w/ Refit) (x)
        TH1F*                   fHPrimVtx_recalPV_KF_Refit_Minus_PV_KF_Refit_y; //!<! Histogram of difference between recalPV_KF and PV_KF (w/ Refit) (y)
        TH1F*                   fHPrimVtx_recalPV_KF_Refit_Minus_PV_KF_Refit_z; //!<! Histogram of difference between recalPV_KF and PV_KF (w/ Refit) (z)
        TH1F*                   fHPrimVtx_recalPV_KF_Minus_PV_x; //!<! Histogram of difference between recalPV_KF and PV (w/o Refit) (x)
        TH1F*                   fHPrimVtx_recalPV_KF_Minus_PV_y; //!<! Histogram of difference between recalPV_KF and PV (w/o Refit) (y)
        TH1F*                   fHPrimVtx_recalPV_KF_Minus_PV_z; //!<! Histogram of difference between recalPV_KF and PV (w/o Refit) (z)
        TH1F*                   fHPrimVtx_PVrec_Minus_PVgen_x; //!<! Histogram of difference between PV(rec) and PV(gen) (x)
        TH1F*                   fHPrimVtx_PVrec_Minus_PVgen_y; //!<! Histogram of difference between PV(rec) and PV(gen) (y)
        TH1F*                   fHPrimVtx_PVrec_Minus_PVgen_z; //!<! Histogram of difference between PV(rec) and PV(gen) (z)
        TH1F*                   fHPrimVtx_PV_KF_Refit_Minus_PVgen_x; //!<! Histogram of difference between PV(rec) from KF and PV(gen) (x)
        TH1F*                   fHPrimVtx_PV_KF_Refit_Minus_PVgen_y; //!<! Histogram of difference between PV(rec) from KF and PV(gen) (y)
        TH1F*                   fHPrimVtx_PV_KF_Refit_Minus_PVgen_z; //!<! Histogram of difference between PV(rec) from KF and PV(gen) (z)
        TH1F*                   fHPrimVtx_recalPV_Minus_PVgen_x; //!<! Histogram of difference between recalPV(rec) and PV(gen) (x)
        TH1F*                   fHPrimVtx_recalPV_Minus_PVgen_y; //!<! Histogram of difference between recalPV(rec) and PV(gen) (y)
        TH1F*                   fHPrimVtx_recalPV_Minus_PVgen_z; //!<! Histogram of difference between recalPV(rec) and PV(gen) (z)
        TH1F*                   fHPrimVtx_recalPV_KF_Refit_Minus_PVgen_x; //!<! Histogram of difference between recalPV(rec) from KF after checking track used for PV and PV(gen) (x)
        TH1F*                   fHPrimVtx_recalPV_KF_Refit_Minus_PVgen_y; //!<! Histogram of difference between recalPV(rec) from KF after checking track used for PV and PV(gen) (y)
        TH1F*                   fHPrimVtx_recalPV_KF_Refit_Minus_PVgen_z; //!<! Histogram of difference between recalPV(rec) from KF after checking track used for PV and PV(gen) (z)
        TH1F*                   fHPrimVtx_PV_PULL_x; //!<! Histogram of PULL of PV(rec) (x)
        TH1F*                   fHPrimVtx_PV_PULL_y; //!<! Histogram of PULL of PV(rec) (y)
        TH1F*                   fHPrimVtx_PV_PULL_z; //!<! Histogram of PULL of PV(rec) (z)
        TH1F*                   fHPrimVtx_PV_KF_PULL_x; //!<! Histogram of PULL of PV(rec) from KF (x)
        TH1F*                   fHPrimVtx_PV_KF_PULL_y; //!<! Histogram of PULL of PV(rec) from KF (y)
        TH1F*                   fHPrimVtx_PV_KF_PULL_z; //!<! Histogram of PULL of PV(rec) from KF (z)
        TH1F*                   fHPrimVtx_PV_KF_Refit_PULL_x; //!<! Histogram of PULL of PV(rec) from KF (x) (w/ Refit)
        TH1F*                   fHPrimVtx_PV_KF_Refit_PULL_y; //!<! Histogram of PULL of PV(rec) from KF (y) (w/ Refit)
        TH1F*                   fHPrimVtx_PV_KF_Refit_PULL_z; //!<! Histogram of PULL of PV(rec) from KF (z) (w/ Refit)
        TH1F*                   fHPrimVtx_recalPV_PULL_x; //!<! Histogram of PULL of recalPV(rec) (x)
        TH1F*                   fHPrimVtx_recalPV_PULL_y; //!<! Histogram of PULL of recalPV(rec) (y)
        TH1F*                   fHPrimVtx_recalPV_PULL_z; //!<! Histogram of PULL of recalPV(rec) (z)
        TH1F*                   fHPrimVtx_recalPV_KF_Refit_PULL_x; //!<! Histogram of PULL of recalPV(rec) from KF after checking track used for PV (x) (w/ Refit)
        TH1F*                   fHPrimVtx_recalPV_KF_Refit_PULL_y; //!<! Histogram of PULL of recalPV(rec) from KF after checking track used for PV (y) (w/ Refit)
        TH1F*                   fHPrimVtx_recalPV_KF_Refit_PULL_z; //!<! Histogram of PULL of recalPV(rec) from KF after checking track used for PV (z) (w/ Refit)
        TH1F*                   fHPrimVtx_recalPV_KF_Minus_PVgen_x; //!<! Histogram of difference between recalPV(rec) from KF after checking track used for PV and PV(gen) (x) (w/o Refit)
        TH1F*                   fHPrimVtx_recalPV_KF_Minus_PVgen_y; //!<! Histogram of difference between recalPV(rec) from KF after checking track used for PV and PV(gen) (y) (w/o Refit)
        TH1F*                   fHPrimVtx_recalPV_KF_Minus_PVgen_z; //!<! Histogram of difference between recalPV(rec) from KF after checking track used for PV and PV(gen) (z) (w/o Refit)
        TH1F*                   fHPrimVtx_recalPV_KF_PULL_x; //!<! Histogram of PULL of recalPV(rec) from KF after checking track used for PV (x) (w/o Refit)
        TH1F*                   fHPrimVtx_recalPV_KF_PULL_y; //!<! Histogram of PULL of recalPV(rec) from KF after checking track used for PV (y) (w/o Refit)
        TH1F*                   fHPrimVtx_recalPV_KF_PULL_z; //!<! Histogram of PULL of recalPV(rec) from KF after checking track used for PV (z) (w/o Refit)
        Bool_t                  fWriteXicPlusMCGenTree; ///< flag to decide whether to write the MC candidate variables on a tree variables
        Bool_t                  fWriteXicPlusTree; ///< flag to decide whether to write Xic+ tree
        Bool_t                  fWriteXicPlusQATree; ///< flag to decide whether to write Xic+ QA tree
        Int_t                   fCount_NumOfCandidatePerEvent_In3Sigma; ///< Count number of Xic+ candidate per event within 3 sigma (assuming sigma=0.01)
        TString                 fFileName;
        unsigned int            fEventNumber;
        unsigned int            fDirNumber;

        AliAnalysisTaskSEXicPlusToXi2PifromKFP(const AliAnalysisTaskSEXicPlusToXi2PifromKFP &source); // not implemented
        AliAnalysisTaskSEXicPlusToXi2PifromKFP& operator=(const AliAnalysisTaskSEXicPlusToXi2PifromKFP& source); // not implemented

        ClassDef(AliAnalysisTaskSEXicPlusToXi2PifromKFP, 7);
};

#endif
