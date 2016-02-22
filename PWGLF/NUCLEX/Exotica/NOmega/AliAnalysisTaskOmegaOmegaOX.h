#ifndef ALIANALYSISTASKOMEGAOMEGAOX_H
#define ALIANALYSISTASKOMEGAOMEGAOX_H

	//------------------------------------------------------------------------------------------
	// version 1.50 (2016/02/22)
	//------------------------------------------------------------------------------------------

/* $Id$ */ 

#include "TROOT.h"
#include "TSystem.h"

#include "AliAnalysisTaskSE.h"

class TH1F;
class TClonesArray;
class AliAODRecoCascade;

/// Missing forward declarations...
class AliESDtrack;
class AliESDEvent;
class AliESDVertex;
class AliESDtrackCuts;
class TTree;
class AliPIDCombined;
class AliPIDResponse;

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
	void SetReqClustersTPC(double tpcclust) {fReqClustersTPC = tpcclust;}
	void SetReqSigmaTOF(double tofsigma) {fReqSigmaTOF = tofsigma;}
	void SetReqPseudoRap(double pseudorap) {fReqPseudoRap = pseudorap;}

	void SetCPADibaryon(double cpadi) {fCPADibaryon = cpadi;}

	void SetRecoTypeDB(int recotype) {fRecoTypeDB = recotype;}

	void MakeAnalysis(TClonesArray *mcArray, AliESDEvent *fESDEvent);


	// set MC usage
	void SetMC(Bool_t theMCon) {fUseMCInfo = theMCon;}
	Bool_t GetMC() const {return fUseMCInfo;}

 private:

	AliAnalysisTaskOmegaOmegaOX(const AliAnalysisTaskOmegaOmegaOX &source);
	AliAnalysisTaskOmegaOmegaOX& operator=(const AliAnalysisTaskOmegaOmegaOX& source); 

	AliESDtrackCuts  *fESDtrackCutsV0; // basic cut variables for v0's
//	AliESDv0Cuts     *fESDCutsV0;      // V0 track cuts
	AliESDtrackCuts  *fESDtrackCuts;   // track cuts

	Bool_t PreTrackCut(AliESDtrack *track);
	Double_t InvMassLambda(Double_t MomPos[3],Double_t MomNeg[3],AliESDtrack *pos,AliESDtrack *neg,Double_t v0Return[1]);
	Double_t InvMassLambdaStar(Double_t MomPos[3],Double_t MomNeg[3],AliESDtrack *pos,AliESDtrack *neg,Double_t v0Return[1]);
	Double_t GetPaFromPxPyPz(Double_t Momentum[3]);
	Double_t GetPtFromPxPyPz(Double_t Momentum[3]);
	void Rotate(Double_t x,Double_t y,Double_t angle);
	void Rotate(Double_t x,Double_t y,Double_t angle,Double_t xCenter,Double_t yCenter);
	Double_t GetAngleFromCosSin(Double_t cos,Double_t sin);

	void DefineTreeVariables();

	AliESDEvent *fESDEvent;       // ESD event
	Bool_t fUseMCInfo;            // Use MC info
	AliPIDResponse *fPIDResponse; // PID response object
	AliPIDCombined *fPIDCombined; // combined PID response object
	Bool_t fIsEventSelected;      // flag for event selected
	TTree *fParametersTree;       // tree of the cut parameters on output slot 1
	TTree *fVariablesTree;        // tree of the candidate variables on output slot 2
	Float_t *fParameters;         // cut parameters to be written to the tree1
	Float_t *fCandidateVariables; // variables to be written to the tree2
	Bool_t fMixedEvent;           // Use mixed event
	AliESDVertex *fVtx1;          // reconstructed vertex from 2 tracks (new v0)
	Double_t fBzkG;               // magnetic field value [kG]
	Int_t fCentrality;
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

	// settings for analysis
	Int_t fRecoTypeDB;            // Reconstruction type of DiBaryon (0:All, 1:OmOm, 2:OmXi, 3:XiOm, 4:XiXi)

	// cut parameters for tracks
	Double_t fReqSigmaTPC;         // TPC PIDcut sigma
	Int_t fReqClustersTPC;         // TPC number of clusters
	Double_t fReqSigmaTOF;         // TOF PIDcut sigma
	Double_t fReqPseudoRap;        // PseudoRapidity

	// cut parameters for dibaryon
	Double_t fCPADibaryon;      // Min cosine of dibaryon's pointing angle

	ClassDef(AliAnalysisTaskOmegaOmegaOX,2); // analysisclass
};
#endif

