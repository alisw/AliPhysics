#ifndef ALIANALYSISTASKNOmegaLX_H
#define ALIANALYSISTASKNOmegaLX_H

	//------------------------------------------------------------------------------------------
	// version 1.62 (2016/03/03)
	//------------------------------------------------------------------------------------------

/* $Id$ */ 

#include "TROOT.h"
#include "TSystem.h"

#include "AliAnalysisTaskSE.h"

class TH1F;
class AliESDtrack;
class TClonesArray;
class AliAODRecoCascade;

class AliAnalysisTaskNOmegaLX : public AliAnalysisTaskSE 
{
 public:
	AliAnalysisTaskNOmegaLX();
	AliAnalysisTaskNOmegaLX(const Char_t* name);
	virtual ~AliAnalysisTaskNOmegaLX();

	// Implementation of interface methods  
	virtual void UserCreateOutputObjects();
	virtual void Init();
	virtual void LocalInit() {Init();}
	virtual void UserExec(Option_t *option);
	virtual void Terminate(Option_t *option);

	void SetReqSigmaTPC(Double_t  tpcsigma) {fReqSigmaTPC = tpcsigma;}
	void SetReqClustersTPC(Int_t  tpcclust) {fReqClustersTPC = tpcclust;}
	void SetReqSigmaTOF(Double_t tofsigma) {fReqSigmaTOF = tofsigma;}
	void SetReqPseudoRap(Double_t pseudorap) {fReqPseudoRap = pseudorap;}
	void SetProtonPMax(Double_t maxprotonp) {fProtonPMax = maxprotonp;}
	void SetPionPMax(Double_t maxpionp) {fPionPMax = maxpionp;}
	void SetKaonPMax(Double_t maxkaonp) {fKaonPMax = maxkaonp;}
	void SetTrackPMin(Double_t mintrkp) {fTrackPMin = mintrkp;}

	void SetCPADibaryon(Double_t cpadi) {fCPADibaryon = cpadi;}
	void SetDCADibaryon(Double_t dcadi) {fDCADibaryon = dcadi;}
	void SetMassWinCascade(Double_t masswinc) {fMassWinCascade = masswinc;}
	void SetWindowLambda(Double_t winlam) {fWindowV02 = winlam;}

	void SetCascChi2max(Double_t cschi2max) {fCsChi2max = cschi2max;}
	void SetCascDV0Min(Double_t csdv0min) {fCsDV0min = csdv0min;}
	void SetCascMassWinLambda(Double_t csmasswinl) {fCsMassWinLambda = csmasswinl;}
	void SetCascDBachMin(Double_t csdbachmin) {fCsDBachMin = csdbachmin;}
	void SetCascDCAmax(Double_t csdcamax) {fCsDCAmax = csdcamax;}
	void SetCascCPAmin(Double_t cscpamin) {fCsCPAmin = cscpamin;}
	void SetCascRmin(Double_t csrmin) {fCsRmin = csrmin;}
	void SetCascRmax(Double_t csrmax) {fCsRmax = csrmax;}

	void SetFiducialVolMin2(Double_t minfv2) {fFiducialVolMin2 = minfv2;}
	void SetFiducialVolMax2(Double_t maxfx2) {fFiducialVolMax2 = maxfx2;}
	void SetPosDCAToPVMin2(Double_t minposdca2) {fPosDCAToPVMin2 = minposdca2;}
	void SetNegDCAToPVMin2(Double_t minnegdca2) {fNegDCAToPVMin2 = minnegdca2;}
	void SetDCADaughterMax2(Double_t maxdcad2) {fDCADaughterMax2 = maxdcad2;}
	void SetCPAMin2(Double_t mincpa2) {fCPAMin2 = mincpa2;}
	void SetDCAToPVMin2(Double_t mindcapv2) {fDCAToPVMin2 = mindcapv2;}
	void SetCOADaughterMin2(Double_t mincoad2) {fCOADaughterMin2 = mincoad2;}
	void SetDCAZDaughterMax2(Double_t maxdcazd2) {fDCAZDaughterMax2 = maxdcazd2;}

	void SetRecoTypeDB(Int_t recotype) {fRecoTypeDB = recotype;}
	void SetLikeSignDB(Int_t likesigndb) {fLikeSignDB = likesigndb;}

	void MakeAnalysis(TClonesArray *mcArray, AliESDEvent *fESDEvent);


	// set MC usage
	void SetMC(Bool_t theMCon) {fUseMCInfo = theMCon;}
	Bool_t GetMC() const {return fUseMCInfo;}

 private:

	AliAnalysisTaskNOmegaLX(const AliAnalysisTaskNOmegaLX &source);
	AliAnalysisTaskNOmegaLX& operator=(const AliAnalysisTaskNOmegaLX& source); 

	Bool_t PreTrackCut(AliESDtrack *track);
	Bool_t ReconstructCascade();

  Double_t Det(Double_t a00, Double_t a01, Double_t a10, Double_t a11) const;
  Double_t Det(Double_t a00,Double_t a01,Double_t a02,
         Double_t a10,Double_t a11,Double_t a12,
         Double_t a20,Double_t a21,Double_t a22) const;
  Double_t PropagateToDCA(AliESDv0 *vtx,AliExternalTrackParam *trk,Double_t b);

	// for V0 (Invariant mass)
	Double_t InvMassLambda(Double_t MomPos[3],Double_t MomNeg[3],AliESDtrack *pos,AliESDtrack *neg,Double_t v0Return[1]);
	Double_t InvMassLambdaStar(Double_t MomPos[3],Double_t MomNeg[3],AliESDtrack *pos,AliESDtrack *neg,Double_t v0Return[1]);

	// others
	Double_t GetPaFromPxPyPz(Double_t Momentum[3]) {return TMath::Sqrt( Momentum[0]*Momentum[0] + Momentum[1]*Momentum[1] + Momentum[2]*Momentum[2] );}
	Double_t GetPtFromPxPyPz(Double_t Momentum[3]) {return TMath::Sqrt( Momentum[0]*Momentum[0] + Momentum[1]*Momentum[1] );}

	void DefineTreeVariables();

	AliESDEvent *fESDEvent;            // ESD event
	Bool_t fUseMCInfo;                 // Use MC info
	AliESDtrackCuts  *fESDtrackCutsV0; // basic cut variables for v0's
	AliESDtrackCuts  *fESDtrackCuts;   // track cuts
	AliPIDResponse *fPIDResponse;      // PID response object
	Bool_t fIsEventSelected;           // flag for event selected
	TTree *fParametersTree;            // tree of the cut parameters on output slot 1
	TTree *fVariablesTree;             // tree of the candidate variables on output slot 2
	Float_t *fParameters;              // cut parameters to be written to the tree1
	Float_t *fCandidateVariables;      // variables to be written to the tree2
	Bool_t fMixedEvent;                // Use mixed event
	AliESDVertex *fVtx1;               // reconstructed vertex from 2 tracks (new v0)
	Double_t fBzkG;                    // magnetic field value [kG]
	Int_t fCentrality;
	Int_t fCountEvent;
	Int_t fIsDCA;
	Int_t fNAll;

	// settings for analysis
	Int_t fRecoTypeDB;   // Reconstruction type of DiBaryon (0:All, 1:OmOm, 2:OmXi, 3:XiOm, 4:XiXi)
	Int_t fLikeSignDB;   // Like-sign of DB (0:ALL, 1:XL, 2:(Xbar)(Lbar), 3:X(Lbar), 4:(Xbar)L)
	Int_t fRecoSelfCasc; // Cascade reconstruction is made by (0:ESD class, 1:by myself)

	// cut parameters for tracks
	Double_t fReqSigmaTPC;  // TPC PIDcut sigma
	Int_t fReqClustersTPC;  // TPC number of clusters
	Double_t fReqSigmaTOF;  // TOF PIDcut sigma
	Double_t fReqPseudoRap; // PseudoRapidity
	Double_t fProtonPMax;   // Max momentum of proton
	Double_t fPionPMax;     // Max momentum of pion
	Double_t fKaonPMax;     // Max momentum of kaon
	Double_t fTrackPMin;    // Min momentum of track

	// cut parameters for dibaryon
	Double_t fCPADibaryon;    // Min cosine of dibaryon's pointing angle
	Double_t fMassWinCascade; //"window" around the Cascade mass
	Double_t fDCADibaryon;					// Max DCA Between two baryons

	// cut parameters for reconstruction of cascade
	Double_t fCsChi2max;       //maximal allowed chi2 
	Double_t fCsDV0min;        //min V0 impact parameter
	Double_t fCsMassWinLambda; //"window" around the Lambda mass
	Double_t fCsDBachMin;      //min bachelor impact parameter
	Double_t fCsDCAmax;        //max DCA between the V0 and the track 
	Double_t fCsCPAmin;        //min cosine of the cascade pointing angle
	Double_t fCsRmin;          //min radius of the fiducial volume
	Double_t fCsRmax;          //max radius of the fiducial volume

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

	ClassDef(AliAnalysisTaskNOmegaLX,2); // analysisclass
};
#endif

