#ifndef ALIANALYSISTASKSESEMILEPTONICOMEGAC0KFP_H
#define ALIANALYSISTASKSESEMILEPTONICOMEGAC0KFP_H

/// \class ALIANALYSISTASKSESEMILEPTONICOMEGAC0KFP
/// \brief class for cuts on AOD reconstructed on Omegac0->e+Omega
///
//***********************************************************
/// \author Tiantian Cheng <chengtiantian@mails.ccnu.edu.cn>, Central China Normal University & GSI Helmholtz Centre for Heavy Ion Research
/// \date Aug 27, 2021
//***********************************************************

/* $Id$ */

#ifndef HomogeneousField
#define HomogeneousField
#endif

#include "TROOT.h"
#include "TVector.h"
#include "TVector2.h"
#include "TSystem.h"
#include "TProfile.h"
#include "THistManager.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliAnalysisTaskSE.h"
#include "AliAODMCParticle.h"
#include "AliNormalizationCounter.h"
#include "THnSparse.h"
#include "AliPIDResponse.h"
#include "AliAODInputHandler.h"
#include "AliVertexingHFUtils.h"
#include "AliPID.h"
#include "AliRDHFCutsOmegactoeleOmegafromKFP.h"

// includes added to play with KFParticle
#include <vector>
#include "KFParticleBase.h"
#include "KFParticle.h"
#include "KFPTrack.h"
#include "KFPVertex.h"
#include "KFVertex.h"


class AliAnalysisTaskSESemileptonicOmegac0KFP : public AliAnalysisTaskSE
{
    public:
                            AliAnalysisTaskSESemileptonicOmegac0KFP ();
                            AliAnalysisTaskSESemileptonicOmegac0KFP ( const char * name, AliRDHFCutsOmegactoeleOmegafromKFP * cuts);
    virtual                  ~AliAnalysisTaskSESemileptonicOmegac0KFP();
    
    virtual void            UserCreateOutputObjects();
    virtual void            Init();
    virtual void            LocalInit() {Init();}
    virtual void            UserExec(Option_t* option);
    virtual void            Terminate(Option_t* option);
    
   
    void                    SelectTrack(const AliVEvent *event, Int_t trkEntries, Int_t &nSeleTrks, Bool_t *seleFlags, TClonesArray *mcArray);
    void                    SelectCascade(const AliVEvent *event,Int_t nCasc,Int_t &nSeleCasc, Bool_t *seleCascFlags, TClonesArray *mcArray);
    Bool_t                  PrefilterElectronULS(AliAODTrack *trk, AliAODEvent *aodEvent, Double_t &mass);
    Bool_t                  PrefilterElectronLS(AliAODTrack *trk, AliAODEvent *aodEvent,  Double_t &samesign_mass);
    
    Bool_t                  MakeMCAnalysis(TClonesArray *mcArray);
    void                    MakeAnaOmegacZeroFromCasc(AliAODEvent *aodEvent, TClonesArray *mcArray, KFParticle PV);
    Int_t                   MatchToMCOmegac0(AliAODTrack* trackProton, AliAODTrack *trackPionMinus, AliAODTrack *trackKaon, AliAODTrack *trackEletron, TClonesArray *mcArray);
    Int_t                   MatchToMCAntiOmegac0(AliAODTrack* trackAntiProton, AliAODTrack *trackPionPlus, AliAODTrack *trackKaon, AliAODTrack *trackEletron, TClonesArray *mcArray);
    
    //---- set MC usage
    void                    SetMC(Bool_t theMCon) {fUseMCInfo = theMCon;}
    Bool_t                  GetMC() const {return fUseMCInfo;}
    
    void                    SetWriteOmegac0Tree(Bool_t a) {fWriteOmegac0Tree = a;}
    Bool_t                  GetWriteOmegac0Tree() const {return fWriteOmegac0Tree;}
    
    void                    SetWriteOmegac0QATree(Bool_t a){fWriteOmegac0QATree = a;}
    Bool_t                  GetWriteOmegac0QATree() const {return fWriteOmegac0QATree;}
    
    void                    SetWriteOmegac0MCGenTree(Bool_t a) {fWriteOmegac0MCGenTree = a;}
    Bool_t                  GetWriteOmegac0MCGenTree() const {return fWriteOmegac0MCGenTree;}
    
    void                    FillTreeGenOmegac0(AliAODMCParticle *mcpart, Int_t CheckOrigin);
    
    void                    FillEventROOTObjects();

    void                    FillTreeRecOmegac0FromCasc(KFParticle kfpOmegac0, KFParticle kfpOmegac0_woMassConst, AliAODTrack *trackElectronFromOmegac0, KFParticle kfpBE, KFParticle kfpOmegaMinus, KFParticle kfpOmegaMinus_m, KFParticle kfpKaon, AliAODTrack *trackKaonFromOmega, AliAODcascade *casc, KFParticle kfpK0Short,  KFParticle kfpLambda, KFParticle kfpLambda_m, AliAODTrack *trkProton, AliAODTrack *trkPion, KFParticle PV, TClonesArray *mcArray, AliAODEvent *aodEvent, Int_t lab_Omegac0, Int_t decaytype);
    
    
    //--- private
    
private:
    void                     DefineEvent();
    void                     DefineTreeRecoOmegac0();
    void                     DefineTreeRecoOmegac0_QA();
    void                     DefineTreeMCGenOmegac0();
    void                     DefineAnaHist();
  
    AliPIDResponse*          fPID; ///<
    AliRDHFCutsOmegactoeleOmegafromKFP*   fAnalCuts; /// !<! Cuts
    AliAODVertex*           fpVtx;                //!<! primary vertex
    AliMCEvent*             fMCEvent;             //!<! corresponding mc event
    Double_t                fBzkG;                ///< magnetic field value [kG]
    TList*                  fListCuts;            //!<! User output
    TList*                  fOutputList;          //!<! Output list
    TTree*                  fTree_Event;          //!<! tree of event
    Float_t*                fVar_Event;           //!<! variables of event to be written to the tree
    TTree*                  fTree_Omegac0;        //!<! tree of the candidate variables
    Float_t*                fVar_Omegac0;         //!<! variables of Omegac0 to be written to the tree
    TTree*                  fTree_Omegac0_QA;        //!<! QA tree of the candidate variables
    Float_t*                fVar_Omegac0_QA;         //!<! QA check for the candidate variables
    TTree*                  fTree_Omegac0MCGen;        //!<! tree of the candidate variables after the track selection on output slot
    Float_t*                fVar_Omegac0MCGen;         //!<! variables of Omegac0 to be written to the tree
   
    AliNormalizationCounter* fCounter; //!<! Counter for normalization
    Bool_t                  fUseMCInfo; ///< Flag of MC analysis
    Bool_t                  fWriteOmegac0Tree;   ///< flag to decide whether to write Omegac0 tree
    Bool_t                  fWriteOmegac0QATree; ///< flag to decide whether to write Omegac0QA tree
    Bool_t                  fWriteOmegac0MCGenTree;  ///<flag to decide whether to write the MC candidate variables on a tree variables
    
    TH1F*                   fHistEvents;          //!<! Histogram of selected events
    TH1F*                   fHTrigger;            //!<! Histogram of trigger
    TH2F*                   fHistoElectronTPCPID;     //!<! TPC electron PID
    TH2F*                   fHistoElectronTOFPID;     //!<! TOF electron PID
       
    THnSparse*              fHistoElectronTPCPIDSelTOF;     //!<! TPC electron PID
    THnSparse*              fHistoMassConversions;          //!<! electron-pairs mass conversion
    THnSparse*              fHistoElectronTPCTOFSelPID;     //!<! TPC, TOF electron PID
    
    /// \cond CLASSIMP
    ClassDef(AliAnalysisTaskSESemileptonicOmegac0KFP,1);   // class for Omegac0 -> e+Omega KFP
    /// \endcond
};

#endif

















