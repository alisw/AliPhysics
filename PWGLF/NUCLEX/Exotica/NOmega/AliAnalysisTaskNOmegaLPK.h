#ifndef ALIANALYSISTASKNOMEGALPK_H
#define ALIANALYSISTASKNOMEGALPK_H

  //------------------------------------------------------------------------------------------
  // version 1.60 (2016/03/01)
  //------------------------------------------------------------------------------------------

/* $Id$ */ 

#include "TROOT.h"
#include "TSystem.h"

#include "AliAnalysisTaskSE.h"

class TH1F;
class AliESDtrack;
class TClonesArray;

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

	void SetReqSigmaTPC(double tpcsigma) {fReqSigmaTPC = tpcsigma;}
	void SetReqClustersTPC(int tpcclust) {fReqClustersTPC = tpcclust;}
	void SetReqSigmaTOF(double tofsigma) {fReqSigmaTOF = tofsigma;}
	void SetReqPseudoRap(double pseudorap) {fReqPseudoRap = pseudorap;}
	void SetProtonPMax(double maxprotonp) {fProtonPMax = maxprotonp;}
	void SetPionPMax(double maxpionp) {fPionPMax = maxpionp;}
	void SetKaonPMax(double maxkaonp) {fKaonPMax = maxkaonp;}

	void SetFiducialVolMin1(double minfv1) {fFiducialVolMin1 = minfv1;}
	void SetFiducialVolMax1(double maxfx1) {fFiducialVolMax1 = maxfx1;}
	void SetPosDCAToPVMin1(double minposdca1) {fPosDCAToPVMin1 = minposdca1;}
	void SetNegDCAToPVMin1(double minnegdca1) {fNegDCAToPVMin1 = minnegdca1;}
	void SetDCADaughterMax1(double maxdcad1) {fDCADaughterMax1 = maxdcad1;}
	void SetCPAMin1(double mincpa1) {fCPAMin1 = mincpa1;}
	void SetDCAToPVMin1(double mindcapv1) {fDCAToPVMin1 = mindcapv1;}
	void SetCOADaughterMin1(double mincoad1) {fCOADaughterMin1 = mincoad1;}
	void SetDCAZDaughterMax1(double maxdcazd1) {fDCAZDaughterMax1 = maxdcazd1;}

	void SetWindowLambda(double winlam) {fWindowV02 = winlam;}
	void SetCPADibaryon(double cpadi) {fCPADibaryonNO = cpadi;}

	void SetRecoTypeV0(int recotypev0) {fRecoTypeV0 = recotypev0;}
	void SetLikeSignDB(int likesigndb) {fLikeSignDB = likesigndb;}

	void MakeAnalysis(TClonesArray *mcArray, AliESDEvent *esdEvent);

	// set MC usage
	void SetMC(Bool_t theMCon) {fUseMCInfo = theMCon;}
	Bool_t GetMC() const {return fUseMCInfo;}

 private:

	AliAnalysisTaskNOmegaLPK(const AliAnalysisTaskNOmegaLPK &source);
	AliAnalysisTaskNOmegaLPK& operator=(const AliAnalysisTaskNOmegaLPK& source); 

	// for V0 (Vertexing)
	Bool_t GetEsdV0Momentum(Int_t typeV0,Int_t id1,Int_t id2,Double_t ReturnV0[6],Double_t MomV0[3],Double_t PosV0[3],Double_t MomV0Pos[3],Double_t MomV0Neg[3]);
	Bool_t GetSelfV0Momentum(Int_t typeV0,AliESDtrack *posTrk,AliESDtrack *negTrk,Double_t ReturnV0[6],Double_t MomV0[3],Double_t PosV0[3],Double_t MomPos[3],Double_t MomNeg[3]);

	// for V0 (Invariant mass)
	Double_t InvMassLambda(Double_t MomPos[3],Double_t MomNeg[3],AliESDtrack *pos,AliESDtrack *neg,Double_t v0Return[1]);
	Double_t InvMassLambdaStar(Double_t MomPos[3],Double_t MomNeg[3],AliESDtrack *pos,AliESDtrack *neg,Double_t v0Return[1]);

	// others
	Double_t GetPaFromPxPyPz(Double_t Momentum[3]) {return TMath::Sqrt( Momentum[0]*Momentum[0] + Momentum[1]*Momentum[1] + Momentum[2]*Momentum[2] );}
	Double_t GetPtFromPxPyPz(Double_t Momentum[3]) {return TMath::Sqrt( Momentum[0]*Momentum[0] + Momentum[1]*Momentum[1] );}

	void DefineTreeVariables();

	AliESDEvent *fESDEvent;            // ESD event
	AliVEvent *fVEvent;                // ESD event
	AliESDtrackCuts *fESDtrackCutsNeg; // basic cut variables for v0's
//	AliESDv0Cuts     *fESDCutsV        // V0 track cuts
	AliESDtrackCuts *fESDtrackCuts;    // track cuts
	Bool_t fUseMCInfo;                 // Use MC info
	AliPIDResponse *fPIDResponse;      // PID response object
	Bool_t fIsEventSelected;           // flag for event selected
	Bool_t fMixedEvent;                // Use mixed event
	TTree *fParametersTree;            // tree of the cut parameters on output slot 1
	TTree *fVariablesTree;             // tree of the candidate variables on output slot 2
	Float_t *fParameters;              // cut parameters to be written to the tree1
	Float_t *fCandidateVariables;      // variables to be written to the tree2
	AliESDVertex *fVtx1;               // reconstructed vertex from 2 tracks (new v0)
	Double_t fBzkG;                    // magnetic field value [kG]
	Int_t fCentrality;
	Int_t fCountEvent;
	Int_t fIsDCA;
	Int_t fNAll;

	// settings for analysis
	Int_t fRecoTypeV0;          // Reconstruction type of V0 (0:N-Omega, 1:H-dibaryon)
	Int_t fLikeSignDB;          // Like-sign of DB (0:ALL, 1:LL, 2:(Lbar)(Lbar), 3:L(Lbar), 4:(Lbar)L)

	// cut parameters for tracks
	Double_t fReqSigmaTPC;      // TPC PIDcut sigma
	Int_t fReqClustersTPC;      // TPC number of clusters
	Double_t fReqSigmaTOF;      // TOF PIDcut sigma
	Double_t fReqPseudoRap;     // PseudoRapidity
	Double_t fProtonPMax;       // Max momentum of proton
	Double_t fPionPMax;         // Max momentum of pion
	Double_t fKaonPMax;         // Max momentum of kaon

	// cut parameters for V01 (Inner V0)
	Double_t fFiducialVolMin1;  // Min radius of the fiducial volume
	Double_t fFiducialVolMax1;  // Max radius of the fiducial volume
	Double_t fPosDCAToPVMin1;   // Min length of impact parameter for the positive track
	Double_t fNegDCAToPVMin1;   // Min length of impact parameter for the negative track
	Double_t fDCADaughterMax1;  // Max DCA between the daughter tracks
	Double_t fCPAMin1;          // Min cosine of V0's pointing angle to PV
	Double_t fDCAToPVMin1;      // Min DCA V0 to PV
	Double_t fCOADaughterMin1;  // Min cosine between the daughter tracks
	Double_t fDCAZDaughterMax1; // Max DCAZ V0 to PV
	Double_t fWindowV01;        // Mass window cut for Lambda

	// cut paramters for V02 (Outer V0)
	Double_t fFiducialVolMin2;  // Min radius of the fiducial volume
	Double_t fFiducialVolMax2;  // Max radius of the fiducial volume
	Double_t fPosDCAToPVMin2;   // Min length of impact parameter for the positive track
	Double_t fNegDCAToPVMin2;   // Min length of impact parameter for the negative track
	Double_t fDCADaughterMax2;  // Max DCA between the daughter tracks
	Double_t fCPAMin2;          // Min cosine of V0's pointing angle to PV
	Double_t fDCAToPVMin2;      // Min DCA V0 to PV
	Double_t fCOADaughterMin2;  // Min cosine between the daughter tracks
	Double_t fDCAZDaughterMax2; // Max DCAZ V0 to PV
	Double_t fWindowV02;        // Mass window cut for Lambda

	// cut parameters for dibaryon
	Double_t fCPAV01toV02;      // Min cosine of V02's pointing angle to V01
	Double_t fCPADibaryonNO;    // Min cosine of dibaryon's pointing angle


	ClassDef(AliAnalysisTaskNOmegaLPK,2); // analysisclass
};
#endif

