#ifndef ALIANALYSISTASKSELC2PKS0FROMKFP_H
#define ALIANALYSISTASKSELC2PKS0FROMKFP_H

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
/// \class AliAnalysisTaskSELc2pKs0fromKFP
/// \brief This is a brief description of my class
///
/// This is a longer description of the class. This longer description is
/// formatted using the Markdown syntax (see below) and can span on multiple
/// lines.
///
/// \author Jianhui Zhu <zjh@mail.ccnu.edu.cn>, Central China Normal University & GSI Helmholtz Centre for Heavy Ion Research
/// \date Jul 27, 2020
////////////////////////////////////////////////////////////////////////////

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

class AliAnalysisTaskSELc2pKs0fromKFP : public AliAnalysisTaskSE
{
    public:
                                AliAnalysisTaskSELc2pKs0fromKFP();
                                AliAnalysisTaskSELc2pKs0fromKFP(const char *name, AliRDHFCutsKFP* cuts);
        virtual                 ~AliAnalysisTaskSELc2pKs0fromKFP();

        virtual void            UserCreateOutputObjects();
        virtual void            Init();
        virtual void            LocalInit() {Init();} 
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);

        void                    SetMC(Bool_t IsMC) {fIsMC=IsMC;}
        void                    SetAnaLc2Lpi(Bool_t IsAnaLc2Lpi) {fIsAnaLc2Lpi=IsAnaLc2Lpi;}
        void                    SelectTrack(AliVEvent *event, Int_t trkEntries, Int_t &nSeleTrks, Bool_t *seleFlags);
        Bool_t                  MakeMCAnalysis(TClonesArray *mcArray);
        void                    MakeAnaLcFromCascadeHF(TClonesArray *arrayLc2pKs0orLpi, AliAODEvent *aodEvent, TClonesArray *mcArray, KFParticle PV);
        Double_t                InvMassV0atPV(AliAODTrack *trk1, AliAODTrack *trk2, Int_t pdg1, Int_t pdg2);
        Int_t                   MatchToMCKs0(AliAODTrack *v0Pos, AliAODTrack *v0Neg, TClonesArray *mcArray);
        Int_t                   MatchToMCLam(AliAODTrack *v0Pos, AliAODTrack *v0Neg, TClonesArray *mcArray, Bool_t IsParticle);
        Int_t                   MatchToMCLc2pKs0(AliAODTrack *v0Pos, AliAODTrack *v0Neg, AliAODTrack *bachPart, TClonesArray *mcArray);
        Int_t                   MatchToMCLc2Lpi(AliAODTrack *v0Pos, AliAODTrack *v0Neg, AliAODTrack *bachPart, TClonesArray *mcArray, Bool_t IsParticle);

        /// set MC usage
        void SetWriteLcMCGenTree(Bool_t a) {fWriteLcMCGenTree = a;}
        Bool_t GetWriteLcMCGenTree() const {return fWriteLcMCGenTree;}

        void SetWriteLcTree(Bool_t a) {fWriteLcTree = a;}
        Bool_t GetWriteLcTree() const {return fWriteLcTree;}

        void SetWriteLcQATree(Bool_t a) {fWriteLcQATree = a;}
        Bool_t GetWriteLcQATree() const {return fWriteLcQATree;}
        void FillEventROOTObjects();
        void FillTreeGenLc(AliAODMCParticle *mcpart, Int_t CheckOrigin);
        void FillTreeRecLcFromCascadeHF(AliAODRecoCascadeHF *Lc2pKs0orLpi, KFParticle kfpLc, AliAODTrack *trackBach, KFParticle kfpBach, KFParticle kfpV0, KFParticle kfpV0_massConstraint, AliAODTrack *v0Pos, AliAODTrack *v0Neg, KFParticle PV, TClonesArray *mcArray, Int_t lab_V0, Int_t lab_Lc, KFParticle kfpLc_woV0MassConst);
        void SetWeightFunction(TF1* weight) {fWeight=weight;}

    private:
        void                    DefineEvent();
        void                    DefineTreeLc_Rec();
        void                    DefineTreeLc_Rec_QA();
        void                    DefineTreeLc_Gen();
        void                    DefineAnaHist();
        AliPIDResponse*         fPID;                 ///<
        AliPIDCombined*         fPIDCombined;         //!<! combined PID response object
        AliRDHFCutsKFP*         fAnaCuts;             ///< Cuts
        AliAODVertex*           fpVtx;                //!<! primary vertex
        AliMCEvent*             fMCEvent;             //!<! corresponding mc event
        Double_t                fBzkG;                ///< magnetic field value [kG]
        Float_t                 fCentrality;           ///< Centrality
//        TObjArray               fMapParticle;         ///< Map of particles in the supporting TClonesArray
        vector<Int_t>           fAodTrackInd;         ///< Translation table: aodTrackInd(mcTrackIndex) = aodTrackIndex
        TList*                  fOutputList;          //!<! Output list
        TList*                  fOutputWeight;        //!<! Output list after weight
        TTree*                  fTree_Event;          //!<! tree of event
        Float_t*                fVar_Event;           //!<! variables of event to be written to the tree
        TTree*                  fTree_Lc;             //!<! tree of the candidate variables
        Float_t*                fVar_Lc;         //!<! variables of Lc to be written to the tree
        TTree*                  fTree_Lc_QA;             //!<! tree of the candidate variables
        Float_t*                fVar_Lc_QA;         //!<! variables of Lc to be written to the tree
        TTree*                  fTree_LcMCGen; //!<! tree of the candidate variables after track selection on output slot
        Float_t*                fVar_LcMCGen;   //!<! variables to be written to the tree
        TList*                  fListCuts;           //!<! User output slot 3 // Cuts 

        Bool_t                  fIsMC; ///< Flag of MC analysis
        Bool_t                  fIsAnaLc2Lpi; ///< Flag of Lc->Lpi analysis

        AliNormalizationCounter* fCounter; //!<! Counter for normalization

        TH1F*                   fHistEvents;          //!<! Histogram of selected events
        TH1F*                   fHTrigger;            //!<! Histogram of trigger
        Bool_t                  fWriteLcMCGenTree; ///< flag to decide whether to write the MC candidate variables on a tree variables
        Bool_t                  fWriteLcTree; ///< flag to decide whether to write Lc tree
        Bool_t                  fWriteLcQATree; ///< flag to decide whether to write QA output tree
        TF1*                    fWeight; ///< weight of Data/MC_gen
        TH1D*                   fHistMCGen_LcPt_weight; //!<! pt of Lc after weight at gen. level
        TH2D*                   f2DHistMCRec_LcPt_weight; //!<! pt of Lc after weight at rec. level

        AliAnalysisTaskSELc2pKs0fromKFP(const AliAnalysisTaskSELc2pKs0fromKFP &source); // not implemented
        AliAnalysisTaskSELc2pKs0fromKFP& operator=(const AliAnalysisTaskSELc2pKs0fromKFP& source); // not implemented

        ClassDef(AliAnalysisTaskSELc2pKs0fromKFP, 5);
};

#endif
