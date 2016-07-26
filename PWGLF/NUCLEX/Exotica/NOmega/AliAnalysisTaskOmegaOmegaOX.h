#ifndef ALIANALYSISTASKOMEGAOMEGAOX_H
#define ALIANALYSISTASKOMEGAOMEGAOX_H

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
class AliAODRecoCascade;

class AliAnalysisTaskOmegaOmegaOX : public AliAnalysisTaskSE 
{
 public:
	AliAnalysisTaskOmegaOmegaOX();
	AliAnalysisTaskOmegaOmegaOX(const Char_t* name);
	virtual ~AliAnalysisTaskOmegaOmegaOX();

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

	void SetCPADibaryon(double cpadi) {fCPADibaryon = cpadi;}
	void SetMassWinCascade(double masswinc) {fMassWinCascade = masswinc;}

	void SetCascChi2max(double cschi2max) {fCsChi2max = cschi2max;}
	void SetCascDV0Min(double csdv0min) {fCsDV0min = csdv0min;}
	void SetCascMassWinLambda(double csmasswinl) {fCsMassWinLambda = csmasswinl;}
	void SetCascDBachMin(double csdbachmin) {fCsDBachMin = csdbachmin;}
	void SetCascDCAmax(double csdcamax) {fCsDCAmax = csdcamax;}
	void SetCascCPAmin(double cscpamin) {fCsCPAmin = cscpamin;}
	void SetCascRmin(double csrmin) {fCsRmin = csrmin;}
	void SetCascRmax(double csrmax) {fCsRmax = csrmax;}

	void SetRecoTypeDB(int recotype) {fRecoTypeDB = recotype;}
	void SetLikeSignDB(int likesigndb) {fLikeSignDB = likesigndb;}

	void MakeAnalysis(TClonesArray *mcArray, AliESDEvent *fESDEvent);


	// set MC usage
	void SetMC(Bool_t theMCon) {fUseMCInfo = theMCon;}
	Bool_t GetMC() const {return fUseMCInfo;}

 private:

	AliAnalysisTaskOmegaOmegaOX(const AliAnalysisTaskOmegaOmegaOX &source);
	AliAnalysisTaskOmegaOmegaOX& operator=(const AliAnalysisTaskOmegaOmegaOX& source); 

	AliESDtrackCuts  *fESDtrackCutsV0; // basic cut variables for v0's
	AliESDtrackCuts  *fESDtrackCuts;   // track cuts

	Bool_t PreTrackCut(AliESDtrack *track);
	Double_t InvMassLambda(Double_t MomPos[3],Double_t MomNeg[3],AliESDtrack *pos,AliESDtrack *neg,Double_t v0Return[1]);
	Double_t InvMassLambdaStar(Double_t MomPos[3],Double_t MomNeg[3],AliESDtrack *pos,AliESDtrack *neg,Double_t v0Return[1]);
	Double_t GetPaFromPxPyPz(Double_t Momentum[3]);
	Double_t GetPtFromPxPyPz(Double_t Momentum[3]);
	Bool_t ReconstructCascade();

  Double_t Det(Double_t a00, Double_t a01, Double_t a10, Double_t a11) const;
  Double_t Det(Double_t a00,Double_t a01,Double_t a02,
         Double_t a10,Double_t a11,Double_t a12,
         Double_t a20,Double_t a21,Double_t a22) const;
  Double_t PropagateToDCA(AliESDv0 *vtx,AliExternalTrackParam *trk,Double_t b);

	void Rotate(Double_t x,Double_t y,Double_t angle);
	void Rotate(Double_t x,Double_t y,Double_t angle,Double_t xCenter,Double_t yCenter);
	Double_t GetAngleFromCosSin(Double_t cos,Double_t sin);

	void DefineTreeVariables();

	AliESDEvent *fESDEvent;       // ESD event
	Bool_t fUseMCInfo;            // Use MC info
	AliPIDResponse *fPIDResponse; // PID response object
	Bool_t fIsEventSelected;      // flag for event selected
	TTree *fParametersTree;       // tree of the cut parameters on output slot 1
	TTree *fVariablesTree;        // tree of the candidate variables on output slot 2
	Float_t *fParameters;         // cut parameters to be written to the tree1
	Float_t *fCandidateVariables; // variables to be written to the tree2
	Bool_t fMixedEvent;           // Use mixed event
	AliESDVertex *fVtx1;          // reconstructed vertex from 2 tracks (new v0)
	Double_t fBzkG;               // magnetic field value [kG]
	Int_t fCentrality;
	Int_t fCountEvent;
	Int_t fIsDCA;
	Int_t fNAll;

	// settings for analysis
	Int_t fRecoTypeDB;   // Reconstruction type of DiBaryon (0:All, 1:OmOm, 2:OmXi, 3:XiOm, 4:XiXi)
	Int_t fLikeSignDB;   // Like-sign of DB (0:ALL, 1:OO, 2:(Obar)(Obar))
	Int_t fRecoSelfCasc; // Cascade reconstruction is made by (0:ESD class, 1:by myself)

	// cut parameters for tracks
	Double_t fReqSigmaTPC;  // TPC PIDcut sigma
	Int_t fReqClustersTPC;  // TPC number of clusters
	Double_t fReqSigmaTOF;  // TOF PIDcut sigma
	Double_t fReqPseudoRap; // PseudoRapidity

	// cut parameters for dibaryon
	Double_t fCPADibaryon; // Min cosine of dibaryon's pointing angle
	Double_t fMassWinCascade; //"window" around the Cascade mass

	// cut parameters for reconstruction of cascade
	Double_t fCsChi2max;       //maximal allowed chi2 
	Double_t fCsDV0min;        //min V0 impact parameter
	Double_t fCsMassWinLambda; //"window" around the Lambda mass
	Double_t fCsDBachMin;      //min bachelor impact parameter
	Double_t fCsDCAmax;        //max DCA between the V0 and the track 
	Double_t fCsCPAmin;        //min cosine of the cascade pointing angle
	Double_t fCsRmin;          //min radius of the fiducial volume
	Double_t fCsRmax;          //max radius of the fiducial volume

	ClassDef(AliAnalysisTaskOmegaOmegaOX,2); // analysisclass
};
#endif

