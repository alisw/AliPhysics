#ifndef ALIANALYSISTASKSEXICCPPTOXICPLUSPIPLUSFROMKFP_H
#define ALIANALYSISTASKSEXICCPPTOXICPLUSPIPLUSFROMKFP_H

/// \class AliAnalysisTaskSEXiccPPToXicPlusPiPlusfromKFP
/// \brief This is a brief description of my class
///
/// This is a longer description of the class. This longer description is
/// formatted using the Markdown syntax (see below) and can span on multiple
/// lines.
///
/// \author Jianhui Zhu <zjh@ccnu.edu.cn>, Central China Normal University & GSI Helmholtz Centre for Heavy Ion Research
/// \date Jun 22, 2021

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

// includes added to play with KFParticle
#include <vector>
#include "KFParticleBase.h"
#include "KFParticle.h"
#include "KFPTrack.h"
#include "KFPVertex.h"
#include "KFVertex.h"

class AliAnalysisTaskSEXiccPPToXicPlusPiPlusfromKFP : public AliAnalysisTaskSE
{
    public:
                                AliAnalysisTaskSEXiccPPToXicPlusPiPlusfromKFP();
                                AliAnalysisTaskSEXiccPPToXicPlusPiPlusfromKFP(const char *name, AliRDHFCutsKFP* cuts);
        virtual                 ~AliAnalysisTaskSEXiccPPToXicPlusPiPlusfromKFP();

        virtual void            UserCreateOutputObjects();
        virtual void            Init();
        virtual void            LocalInit() {Init();} 
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);

        void                    SetMC(Bool_t IsMC) {fIsMC=IsMC;}
        void                    SelectTrack(AliVEvent *event, Int_t trkEntries, Int_t &nSeleTrks, Bool_t *seleFlags);
        Bool_t                  MakeMCAnalysis(TClonesArray *mcArray);
        void                    MakeAnaXiccPPFromCasc(AliAODEvent *AODEvent, TClonesArray *mcArray, KFParticle PV);
        Int_t                   MatchToMCXiccPP(AliAODTrack *trackProton, AliAODTrack *trackPionMinus_2, AliAODTrack *trackPionMinus_1, AliAODTrack *trackPionPlus_trk1, AliAODTrack *trackPionPlus_trk2, AliAODTrack *trackPionPlus_trk3, TClonesArray *mcArray);
        Int_t                   MatchToMCAntiXiccPP(AliAODTrack *trackAntiProton, AliAODTrack *trackPionPlus_2, AliAODTrack *trackPionPlus_1, AliAODTrack *trackPionMinus_trk1, AliAODTrack *trackPionMinus_trk2, AliAODTrack *trackPionMinus_trk3, TClonesArray *mcArray);
        Int_t                   MatchToMCXiMinus(AliAODTrack *trackProton, AliAODTrack *trackPion3, AliAODTrack *trackPion2, TClonesArray *mcArray);
        Int_t                   MatchToMCXiPlus(AliAODTrack *trackAntiProton, AliAODTrack *trackAntiPion3, AliAODTrack *trackAntiPion2, TClonesArray *mcArray);
        Int_t                   MatchToMCLambda(AliAODTrack *trackProton, AliAODTrack *trackPion3, TClonesArray *mcArray);
        Int_t                   MatchToMCAntiLambda(AliAODTrack *trackAntiProton, AliAODTrack *trackAntiPion3, TClonesArray *mcArray);
        Int_t                   MatchToMCLambdaFromXi(AliAODTrack *trackProton, AliAODTrack *trackPion3, TClonesArray *mcArray);
        Int_t                   MatchToMCAntiLambdaFromXi(AliAODTrack *trackAntiProton, AliAODTrack *trackAntiPion3, TClonesArray *mcArray);
        Int_t                   MatchToMCPion(AliAODTrack *track, TClonesArray *mcArray);
        Double_t                InvMassV0atPV(AliAODTrack *trk1, AliAODTrack *trk2, Int_t pdg1, Int_t pdg2);

        /// set MC usage
        void SetWriteXiccPPMCGenTree(Bool_t a) {fWriteXiccPPMCGenTree = a;}
        Bool_t GetWriteXiccPPMCGenTree() const {return fWriteXiccPPMCGenTree;}

        void SetWriteXiccPPTree(Bool_t a) {fWriteXiccPPTree = a;}
        Bool_t GetWriteXiccPPTree() const {return fWriteXiccPPTree;}

        void FillEventROOTObjects();
        void FillTreeGenXiccPP(AliAODMCParticle *mcpart, Int_t CheckOrigin);
        void FillTreeRecXiccPPFromCasc(KFParticle kfpXiccPP, KFParticle kfpXicPlus, AliAODTrack *trackPiFromXiccPP, KFParticle kfpPiFromXiccPP, AliAODTrack *trackPiFromXicPlus_trk1, KFParticle kfpBP_trk1, KFParticle kfpXiMinus, KFParticle kfpXiMinus_m, KFParticle kfpPionOrKaon, AliAODTrack *trackPiFromXiOrKaonFromOmega, KFParticle kfpK0Short, KFParticle kfpGamma, KFParticle kfpLambda, KFParticle kfpLambda_m, AliAODTrack *trkProton, AliAODTrack *trkPion, AliAODTrack *trackPiFromXicPlus_trk2, KFParticle kfpBP_trk2, KFParticle kfpProtonFromLam, KFParticle kfpPionFromLam, KFParticle PV, TClonesArray *mcArray, Int_t lab_XiccPP);

    private:
        void                    DefineEvent();
        void                    DefineTreeRecXiccPP();
        void                    DefineTreeGenXiccPP();
        void                    DefineAnaHist();
        AliPIDResponse*         fPID;                 ///<
        AliRDHFCutsKFP*         fAnaCuts;             ///< Cuts
        AliAODVertex*           fpVtx;                //!<! primary vertex
        AliMCEvent*             fMCEvent;             //!<! corresponding mc event
        Double_t                fBzkG;                ///< magnetic field value [kG]
        Float_t                 fCentrality;           ///< Centrality
        vector<Int_t>           fAodTrackInd;         ///< Translation table: aodTrackInd(mcTrackIndex) = aodTrackIndex
        TList*                  fOutputList;          //!<! Output list
        TTree*                  fTree_Event;          //!<! tree of event
        Float_t*                fVar_Event;           //!<! variables of event to be written to the tree
        TTree*                  fTree_XiccPP;             //!<! tree of the candidate variables
        Float_t*                fVar_XiccPP;         //!<! variables of XiccPP to be written to the tree
        TTree*                  fTree_XiccPPMCGen; //!<! tree of the candidate variables after track selection on output slot
        Float_t*                fVar_XiccPPMCGen;   //!<! variables to be written to the tree
        TList*                  fListCuts;           //!<! User output slot 3 // Cuts 

        Bool_t                  fIsMC; ///< Flag of MC analysis

        AliNormalizationCounter* fCounter; //!<! Counter for normalization

        TH1F*                   fHistEvents;          //!<! Histogram of selected events
        TH1F*                   fHTrigger;            //!<! Histogram of trigger
        TH1F*                   fHCentrality;          //!<! Histogram of centrality
        Bool_t                  fWriteXiccPPMCGenTree; ///< flag to decide whether to write the MC candidate variables on a tree variables
        Bool_t                  fWriteXiccPPTree; ///< flag to decide whether to write XicZero tree

        AliAnalysisTaskSEXiccPPToXicPlusPiPlusfromKFP(const AliAnalysisTaskSEXiccPPToXicPlusPiPlusfromKFP &source); // not implemented
        AliAnalysisTaskSEXiccPPToXicPlusPiPlusfromKFP& operator=(const AliAnalysisTaskSEXiccPPToXicPlusPiPlusfromKFP& source); // not implemented

        ClassDef(AliAnalysisTaskSEXiccPPToXicPlusPiPlusfromKFP, 1);
};

#endif
