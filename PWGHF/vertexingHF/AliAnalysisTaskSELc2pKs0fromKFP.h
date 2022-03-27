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
#include "AliVVertex.h"

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

        enum EAnalysisType { /// enum for setting analysis system/year (for loading profile histograms for multiplicity correction)
           kpPb2016,
           kpp2016
           };
  
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
        Bool_t                  MakeMCAnalysis(TClonesArray *mcArray, AliAODMCHeader *mcHeader, AliAODEvent *aodEvent);
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
        void FillEventROOTObjects(AliAODEvent* aodEvent);
        void FillTreeGenLc(TClonesArray *mcArray, AliAODMCParticle *mcpart, Int_t CheckOrigin, AliAODMCHeader *mcHeader, AliAODEvent *aodEvent);
        void FillTreeRecLcFromCascadeHF(AliAODRecoCascadeHF *Lc2pKs0orLpi, KFParticle kfpLc, AliAODTrack *trackBach, KFParticle kfpBach, KFParticle kfpV0, KFParticle kfpV0_massConstraint, AliAODTrack *v0Pos, AliAODTrack *v0Neg, KFParticle PV, TClonesArray *mcArray, Int_t lab_V0, Int_t lab_Lc, KFParticle kfpLc_woV0MassConst, AliAODEvent *aodEvent, AliAODVertex *ownPVtx);
        void SetWeightFunction(TF1* weight) {fWeight=weight;}

        void SetUseWeights(Bool_t opt) { fUseWeights = opt;}
        void SetUseMult(Bool_t opt) { fUseMult = opt;}
        void SetKeepOnlyMCSignal(Bool_t opt) {fKeepOnlyMCSignal = opt;}
        void SetKeepAllVariables (Bool_t opt) {fKeepAllVariables = opt;}
        void SetAnalysisType(Int_t opt) { fAnalysisType = opt;}
        void SetReferenceMultiplicity(Double_t opt) {fRefMult = opt;}
        void SetMultVsZProfileLHC16qt1stBunch(TProfile* hprof){
          if(fMultEstimatorAvg[0]) delete fMultEstimatorAvg[0];
          fMultEstimatorAvg[0]=new TProfile(*hprof);
        }
        void SetMultVsZProfileLHC16qt2ndBunch(TProfile* hprof){
          if(fMultEstimatorAvg[1]) delete fMultEstimatorAvg[1];
          fMultEstimatorAvg[1]=new TProfile(*hprof);
        }
        void SetMultVsZProfileLHC16qt3rdBunch(TProfile* hprof){
          if(fMultEstimatorAvg[2]) delete fMultEstimatorAvg[2];
          fMultEstimatorAvg[2]=new TProfile(*hprof);
        }
        void SetMultVsZProfileLHC16qt4thBunch(TProfile* hprof){
          if(fMultEstimatorAvg[3]) delete fMultEstimatorAvg[3];
          fMultEstimatorAvg[3]=new TProfile(*hprof);
        }

        void SetMultVsZProfileLHC16j(TProfile* hprof){
          if(fMultEstimatorAvg[0]) delete fMultEstimatorAvg[0];
          fMultEstimatorAvg[0]=new TProfile(*hprof);
        }
        void SetMultVsZProfileLHC16k(TProfile* hprof){
          if(fMultEstimatorAvg[1]) delete fMultEstimatorAvg[1];
          fMultEstimatorAvg[1]=new TProfile(*hprof);
        }
        void SetMultVsZProfileLHC16l(TProfile* hprof){
          if(fMultEstimatorAvg[2]) delete fMultEstimatorAvg[2];
          fMultEstimatorAvg[2]=new TProfile(*hprof);
        }
        void SetUseOnTheFlyV0(Bool_t opt) {fUseOnTheFlyV0 = opt;}
        Bool_t GetUseOnTheFlyV0() {return fUseOnTheFlyV0;}

        TProfile* GetEstimatorHistogram(const AliVEvent* event);

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
        AliVVertex*           fpVtxOff;                //!<! primary vertex const off
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
        Bool_t                  fUseWeights; ///< Flag to use pT weight functions
        Bool_t                  fKeepOnlyMCSignal; ///< flag to keep only signal candidates
        Bool_t                  fKeepAllVariables; ///<flag to keep all possible variables that were removed to reduce the tree size
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
        TF1*                    fFuncWeightPythia; ///<flat pT weight vs pythia
        TF1*                    fFuncWeightFONLL5overLHC13d3; ///< pT weight vs FONLL (D meson case)
        TF1*                    fFuncWeightFONLL5overLHC13d3Lc; ///< pT weight vs FONLL (Lc case)
        
        Bool_t                  fUseMult; /// switch for multiplicity in tree 
        TProfile* fMultEstimatorAvg[4]; /// TProfile with mult vs. Z per period
        Double_t                fRefMult;      ///reference multiplicity for ntrk correction
        Int_t                   fAnalysisType; ///< switch for analysis period (for multiplicity corrections)
        Bool_t                  fUseOnTheFlyV0; ///< switch for use of on-the-fly V0s

        AliAnalysisTaskSELc2pKs0fromKFP(const AliAnalysisTaskSELc2pKs0fromKFP &source); // not implemented
        AliAnalysisTaskSELc2pKs0fromKFP& operator=(const AliAnalysisTaskSELc2pKs0fromKFP& source); // not implemented

        ClassDef(AliAnalysisTaskSELc2pKs0fromKFP, 9);
};

#endif
