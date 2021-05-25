#ifndef AliAnalysisTaskLambdaNNRun2_H
#define AliAnalysisTaskLambdaNNRun2_H

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliAODTrack.h"

class AliAnalysisTaskLambdaNNRun2 : public AliAnalysisTaskSE
{
public:
	class AnalysisV0 {
	public:
		enum Type { Unknown, LambdaNN, AntiLambdaNN };

		Type     topology;
		Bool_t   onlineV0; // V0 finder online/offline
		Double_t dca;
		Double_t cosPointing;
		Double_t decayRadius;
		Double_t mass;

		Int_t    triton_Charge;
		Double_t triton_P;
		Double_t triton_PTPC;
		Double_t triton_Px;
		Double_t triton_Py;
		Double_t triton_Pz;
		Double_t triton_Eta;
		UShort_t triton_TPCNcls;
		Double_t triton_Chi2perNDF;
		Double_t triton_TPCsignal;
		Double_t triton_NSigmaTPC;
		Double_t triton_NSigmaTOF;
		Double_t triton_TPCExpSignal;
		Double_t triton_TOFExpSignal;

		Int_t    pion_Charge;
		Double_t pion_P;
		Double_t pion_PTPC;
		Double_t pion_Px;
		Double_t pion_Py;
		Double_t pion_Pz;
		Double_t pion_Eta;
		UShort_t pion_TPCNcls;
		Double_t pion_Chi2perNDF;
		Double_t pion_TPCsignal;
		Double_t pion_NSigmaTPC;
		Double_t pion_NSigmaTOF;
		Double_t pion_TPCExpSignal;
		Double_t pion_TOFExpSignal;

		Bool_t isCowboy;

		void Reset() {
			topology = Type::Unknown;
			dca = -1.;
			cosPointing = -2.;
			decayRadius = -1.;
			mass = 0.;

			triton_Charge       = -2;
			triton_P            = 0.;
			triton_PTPC         = 0.;
			triton_Px           = 0.;
			triton_Py           = 0.;
			triton_Pz           = 0.;
			triton_Eta          = -2.;
			triton_TPCNcls      = 0;
			triton_Chi2perNDF   = -1.;
			triton_TPCsignal    = -1.;
			triton_NSigmaTPC    = -20.;
			triton_NSigmaTOF    = -20.;
			triton_TPCExpSignal = -1.;
			triton_TOFExpSignal = -1.;

			pion_Charge       = -2;
			pion_P            = 0.;
			pion_PTPC         = 0.;
			pion_Px           = 0.;
			pion_Py           = 0.;
			pion_Pz           = 0.;
			pion_Eta          = -2.;
			pion_TPCNcls      = 0;
			pion_Chi2perNDF   = -1.;
			pion_TPCsignal    = -1.;
			pion_NSigmaTPC    = -20.;
			pion_NSigmaTOF    = -20.;
			pion_TPCExpSignal = -1.;
			pion_TOFExpSignal = -1.;

			isCowboy = 0;

		};
	};

	class AnalysisEvent {
	public:
		Int_t                   RunNumber;                  // Run number of the event
		Double_t                MagField;                   // Magentic field
		Float_t                 Centrality;                 // Centrality
		size_t                  NV0s;                       // Number of V0s
		std::vector<AnalysisV0> EventV0s;                   // V0s in the event
		void Init() {
			RunNumber  = -1;
			MagField   = -1.;
			Centrality = -1.;
			NV0s       = 0;
			EventV0s.clear();
		}
		void AddV0(const AnalysisV0& v0) { EventV0s.push_back(v0); }
		size_t Size() const { return EventV0s.size(); }
	};


private:
	enum MassType_t {
		kMassPion,
		kMassProton,
		kMasstriton,
		kMassTriton,
		kMassHelium3,
		kMassMax  // Number of enum entries
	};
	static const Double_t	fgkMass[];              //! masses

public:
	AliAnalysisTaskLambdaNNRun2();
	AliAnalysisTaskLambdaNNRun2(const char *name, Bool_t isMC=kFALSE);
	virtual                 ~AliAnalysisTaskLambdaNNRun2();

	virtual void            UserCreateOutputObjects();
	virtual void            UserExec(Option_t* option);
	virtual void            Terminate(Option_t* option);
	//	void SetMC(Bool_t isMC)                  {fMC   = isMC;};
private:
	bool                    LooseTrackCuts(AliAODTrack*);
	void                    FillEvent(AnalysisV0::Type, AliAODTrack*, AliAODTrack*);


	TTree*                  fOutputTree;            //! output tree
	AnalysisEvent*          fOutputEvent;           //! event class
	Bool_t                  fMC;               			// isMC
	AliAODEvent*            fAOD;                   //! input event
	AliMCEvent*             fMCEvent;               //! corresponding MC event
	AliEventCuts            fEventCut;
	AliPIDResponse*         fPID;                   //! PID response
	TList*                  fOutputList;            //! output list

	AnalysisV0              fAnalysis_V0;
	AliAnalysisTaskLambdaNNRun2(const AliAnalysisTaskLambdaNNRun2&); // not implemented
	AliAnalysisTaskLambdaNNRun2& operator=(const AliAnalysisTaskLambdaNNRun2&); // not implemented

	ClassDef(AliAnalysisTaskLambdaNNRun2, 1);
};

#endif
