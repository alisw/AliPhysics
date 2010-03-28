#ifndef ALIANALYSISTASKCHECKV0_H
#define ALIANALYSISTASKCHECKV0_H

/*  See cxx source for full Copyright notice */

//-----------------------------------------------------------------
//                 AliAnalysisTaskCheckV0 class
//            This task is for QAing the V0s from ESD/AOD
//              Origin: B.H. Nov2007, hippolyt@in2p3.fr
//-----------------------------------------------------------------

class TString;
class TList;
class TH1F;
class TH2F;

#include "AliAnalysisTaskSE.h"

class AliESDEvent;
class AliESDVertex;
class AliAODEvent;

class AliAnalysisTaskCheckV0 : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskCheckV0();
  AliAnalysisTaskCheckV0(const char *name);
 ~AliAnalysisTaskCheckV0();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void   SetCollidingSystems(Short_t collidingSystems = 0) {fCollidingSystems = collidingSystems;}
  void   SetUsePhysicsSelection(Bool_t usePhysicsSelection = 0) {fUsePhysicsSelection = usePhysicsSelection;}
  void   SetAnalysisType(const char* analysisType) {fAnalysisType = analysisType;}
  void   SetMaxPrimaryVtxPosZ(const Float_t maxPrimaryVtxPosZ = 100.) {fMaxPrimaryVtxPosZ = maxPrimaryVtxPosZ;}
  void   SetMinV0Pt(const Float_t minV0Pt =   0.) {fMinV0Pt = minV0Pt;}
  void   SetMaxV0Pt(const Float_t maxV0Pt = 100.) {fMaxV0Pt = maxV0Pt;}
  void   SetMaxV0Rapidity(const Float_t maxV0Rapidity = 1.) {fMaxV0Rapidity = maxV0Rapidity;}
  void   SetMinDaughterTpcClusters(const Int_t minDaughterTpcClusters = 80) {fMinDaughterTpcClusters = minDaughterTpcClusters;}
  
 private:
  TString      fAnalysisType;                   //  ESD or AOD
  Short_t      fCollidingSystems;               //  Colliding systems 0/1 for pp/PbPb
  Bool_t       fUsePhysicsSelection;            //  Delegate event selection to AliPhysicsSelectionTask
  Float_t      fMaxPrimaryVtxPosZ;              //  Primary vertex selection in Z
  Float_t      fMinV0Pt;                        //  Minimum pt selection for the V0
  Float_t      fMaxV0Pt;                        //  Maximum pt selection for the V0
  Float_t      fMaxV0Rapidity;                  //  Maximum rapidity selection for the V0
  Int_t        fMinDaughterTpcClusters;         //  Minimum number of TPC clusters for the both daughter tracks of the V0
  TList       *fListHist;                       //! List of histograms
  TH1F        *fHistPrimaryVertexPosX;          //! Primary vertex position in X
  TH1F        *fHistPrimaryVertexPosY;          //! Primary vertex position in Y
  TH1F        *fHistPrimaryVertexPosZ;          //! Primary vertex position in Z
  TH1F        *fHistKeptPrimaryVertexPosX;      //! Primary vertex position in X after event selection
  TH1F        *fHistKeptPrimaryVertexPosY;      //! Primary vertex position in Y after event selection
  TH1F        *fHistKeptPrimaryVertexPosZ;      //! Primary vertex position in Z after event selection
  TH1F        *fHistTrackMultiplicity;          //! Track multiplicity distribution
  TH1F        *fHistV0Multiplicity;             //! V0 multiplicity distribution
  TH1F        *fHistV0OnFlyStatus;              //! V0 on fly status distribution

              // V0 offline distributions
  TH1F        *fHistV0MultiplicityOff;          //! V0 multiplicity distribution offline
  TH1F        *fHistV0Chi2Off;                  //! V0 chi2 distribution
  TH1F        *fHistDcaV0DaughtersOff;          //! Dca between V0 daughters
  TH1F        *fHistV0CosineOfPointingAngleOff; //! Cosine of V0 pointing angle
  TH1F        *fHistV0RadiusOff;                //! V0 radial distance distribution
  TH1F        *fHistDcaV0ToPrimVertexOff;       //! Dca of V0 to primary vertex
  TH1F        *fHistDcaPosToPrimVertexOff;      //! Dca of V0 positive daughter to primary vertex
  TH1F        *fHistDcaNegToPrimVertexOff;      //! Dca of V0 negative daughter to primary vertex

  TH1F        *fHistMassK0sOff;                 //! Invariant mass of K0s
  TH1F        *fHistMassLambdaOff;              //! Invariant mass of Lambda
  TH1F        *fHistMassAntiLambdaOff;          //! Invariant mass of Anti-Lambda
  TH2F        *fHistMassK0sOffVsPt;             //! Invariant mass of K0s
  TH2F        *fHistMassLambdaOffVsPt;          //! Invariant mass of Lambda
  TH2F        *fHistMassAntiLambdaOffVsPt;      //! Invariant mass of Anti-Lambda
  TH2F        *fHistArmenterosPodolanskiOff;    //! Armenteros-Podolanski distribution       

              // V0 on-the-fly distributions
  TH1F        *fHistV0MultiplicityOn;           //! V0 multiplicity distribution on-the-fly
  TH1F        *fHistV0Chi2On;                   //! V0 chi2 distribution
  TH1F        *fHistDcaV0DaughtersOn;           //! Dca between V0 daughters
  TH1F        *fHistV0CosineOfPointingAngleOn;  //! Cosine of V0 pointing angle
  TH1F        *fHistV0RadiusOn;                 //! V0 radial distance distribution
  TH1F        *fHistDcaV0ToPrimVertexOn;        //! Dca of V0 to primary vertex
  TH1F        *fHistDcaPosToPrimVertexOn;       //! Dca of V0 positive daughter to primary vertex
  TH1F        *fHistDcaNegToPrimVertexOn;       //! Dca of V0 negative daughter to primary vertex

  TH1F        *fHistMassK0sOn;                  //! Invariant mass of K0s
  TH1F        *fHistMassLambdaOn;               //! Invariant mass of Lambda
  TH1F        *fHistMassAntiLambdaOn;           //! Invariant mass of Anti-Lambda
  TH2F        *fHistMassK0sOnVsPt;              //! Invariant mass of K0s
  TH2F        *fHistMassLambdaOnVsPt;           //! Invariant mass of Lambda
  TH2F        *fHistMassAntiLambdaOnVsPt;       //! Invariant mass of Anti-Lambda
  TH2F        *fHistArmenterosPodolanskiOn;     //! Armenteros-Podolanski distribution       
   
  AliAnalysisTaskCheckV0(const AliAnalysisTaskCheckV0&);            // not implemented
  AliAnalysisTaskCheckV0& operator=(const AliAnalysisTaskCheckV0&); // not implemented
  
  ClassDef(AliAnalysisTaskCheckV0, 1);
};

#endif
