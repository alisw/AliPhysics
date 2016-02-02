#ifndef ALIANALYSISTASKNOMEGALPK_H
#define ALIANALYSISTASKNOMEGALPK_H

/* $Id$ */ 

#include "TROOT.h"
#include "TSystem.h"

#include "AliAnalysisTaskSE.h"
#include "AliAODEvent.h"
#include "AliPID.h"
#include "AliAODTrack.h"
#include "AliNormalizationCounter.h"

class TH1F;
class TClonesArray;
class AliAODRecoCascade;

class AliAnalysisTaskNOmegaLPK : public AliAnalysisTaskSE 
{
 public:
  AliAnalysisTaskNOmegaLPK();
  AliAnalysisTaskNOmegaLPK(const Char_t* name);
  virtual ~AliAnalysisTaskNOmegaLPK();

  // Implementation of interface methods  
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

//  void FillSpectrum(AliAODcascade *casc, AliAODTrack *part, Int_t &nSelectedAnal, TClonesArray *mcArray);

  void SetFiducialVolMin1(double minfv1) {fFiducialVolMin1 = minfv1;}
  void SetFiducialVolMax1(double maxfx1) {fFiducialVolMax1 = maxfx1;}
  void SetPosDCAToPVMin1(double minposdca1) {fPosDCAToPVMin1 = minposdca1;}
  void SetNegDCAToPVMin1(double minnegdca1) {fNegDCAToPVMin1 = minnegdca1;}
  void SetDCADaughterMax1(double maxdcad1) {fDCADaughterMax1 = maxdcad1;}
  void SetCPAMin1(double mincpa1) {fCPAMin1 = mincpa1;}
  void SetDCAToPVMin1(double mindcapv1) {fDCAToPVMin1 = mindcapv1;}
  void SetCOADaughterMin1(double mincoad1) {fCOADaughterMin1 = mincoad1;}
  void SetDCAZDaughterMax1(double maxdcazd1) {fDCAZDaughterMax1 = maxdcazd1;}

  void SetWindowLambda(double winlam = 0.0045) {fWindowLambda = winlam;}
  void SetCPADibaryon(double cpadi = 0.99875) {fCPADibaryon = cpadi;}

  void MakeAnalysis(TClonesArray *mcArray, Int_t &nSelectedAnal, AliESDEvent *esdEvent);


  // set MC usage
  void SetMC(Bool_t theMCon) {fUseMCInfo = theMCon;}
  Bool_t GetMC() const {return fUseMCInfo;}

 private:

  AliAnalysisTaskNOmegaLPK(const AliAnalysisTaskNOmegaLPK &source);
  AliAnalysisTaskNOmegaLPK& operator=(const AliAnalysisTaskNOmegaLPK& source); 

	AliESDtrackCuts  *fESDtrackCutsV0; // basic cut variables for v0's
//	AliESDv0Cuts     *fESDCutsV0;      // V0 track cuts
	AliESDtrackCuts  *fESDtrackCuts;   // track cuts

	Bool_t PreTrackCut(AliESDtrack *track);
	Bool_t GetEsdV0Momentum(AliESDEvent *esdEvent, Int_t id1,Int_t id2,Double_t cutV0Parameters[9],Double_t ReturnV0[6],Double_t MomV0[3],Double_t PosV0[3],Double_t MomV0Pos[3],Double_t MomV0Neg[3]);
	Bool_t GetSelfV0Momentum(AliESDtrack *posTrk,AliESDtrack *negTrk,Double_t cutV0Parameters[9],Double_t ReturnV0[6],Double_t MomV0[3],Double_t PosV0[3],Double_t MomPos[3],Double_t MomNeg[3]);
	Bool_t GetSelfV0CrossMomentum(AliESDtrack *posTrk,AliESDtrack *negTrk,Double_t cutV0Parameters[9],Double_t ReturnV0[6],Double_t MomV0[3],Double_t PosV0[3],Double_t MomPos[3],Double_t MomNeg[3]);
	Double_t InvMassLambda(Double_t MomPos[3],Double_t MomNeg[3],AliESDtrack *pos,AliESDtrack *neg,Double_t v0Return[1]);
	Double_t InvMassLambdaStar(Double_t MomPos[3],Double_t MomNeg[3],AliESDtrack *pos,AliESDtrack *neg,Double_t v0Return[1]);
	Double_t GetPaFromPxPyPz(Double_t Momentum[3]);
	Double_t GetPtFromPxPyPz(Double_t Momentum[3]);

  void DefineTreeVariables();

  Bool_t fUseMCInfo;             // Use MC info
  AliPIDResponse *fPIDResponse;  // PID response object
  AliPIDCombined *fPIDCombined;  // combined PID response object
  Bool_t fIsEventSelected;       // flag for event selected
  TTree *fVariablesTree1;        // tree of the candidate variables on output slot 1
  TTree *fVariablesTree2;        // tree of the candidate variables on output slot 2
  Float_t *fCandidateVariables1; // variables to be written to the tree1
  Float_t *fCandidateVariables2; // variables to be written to the tree2
	Bool_t fMixedEvent;            // Use mixed event
	AliESDVertex *fVtx1;           // reconstructed vertex from 2 tracks (new v0)
  Double_t fBzkG;                // magnetic field value [kG]
	Int_t fCentrality;
	Double_t fSigmaTPC;
	Int_t fClustersTPC;
	Double_t fSigmaTOF;
	Double_t fPseudoRap;
	Int_t fCountMatch;
	Int_t fCountOnlySelf;
	Int_t fCountOnlyV0;
	Int_t fCountLambda;
	Int_t fCountAntiLambda;
	Int_t fCountLambdaSelf;
	Int_t fCountAntiLambdaSelf;
	Int_t fCountMatchLambda;
	Int_t fCountMatchAntiLambda;
	Int_t fCountEvent;
	Int_t fIsDCA;
	Int_t fNAll;

	Double_t fFiducialVolMin1;  // Min radius of the fiducial volume
	Double_t fFiducialVolMax1;  // Max radius of the fiducial volume
	Double_t fPosDCAToPVMin1;   // Min length of impact parameter for the positive track
	Double_t fNegDCAToPVMin1;   // Min length of impact parameter for the negative track
	Double_t fDCADaughterMax1;  // Max DCA between the daughter tracks
	Double_t fCPAMin1;          // Min cosine of V0's pointing angle
	Double_t fDCAToPVMin1;      // Min DCA V0 to PV
	Double_t fCOADaughterMin1;  // Min cosine between the daughter tracks
	Double_t fDCAZDaughterMax1; // Max DCAZ V0 to PV

	Double_t fWindowLambda;     //Mass window cut for Lambda
	Double_t fCPADibaryon;      //Mass window cut for Lambda

  ClassDef(AliAnalysisTaskNOmegaLPK,5); // class for Lc->p K0
};
#endif

