#ifndef AliAnalysisTaskLambdaNRun2_H
#define AliAnalysisTaskLambdaNRun2_H

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliAODTrack.h"

class AliAnalysisTaskLambdaNRun2 : public AliAnalysisTaskSE
{
public:
	class AnalysisV0 {
	public:
		enum Type { Unknown, LambdaN, AntiLambdaN };

		Type     topology;
		Bool_t   onlineV0; // V0 finder online/offline
		Double_t dca;
		Double_t cosPointing;
		Double_t decayRadius;
		Double_t mass;

		Int_t    deuteron_Charge;
		Double_t deuteron_P;
		Double_t deuteron_PTPC;
		Double_t deuteron_Px;
		Double_t deuteron_Py;
		Double_t deuteron_Pz;
		Double_t deuteron_Eta;
		UShort_t deuteron_TPCNcls;
		Double_t deuteron_Chi2perNDF;
		Double_t deuteron_TPCsignal;
		Double_t deuteron_NSigmaTPC;
		Double_t deuteron_NSigmaTOF;
		Double_t deuteron_TPCExpSignal;
		Double_t deuteron_TOFExpSignal;

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


		void Reset() {
			topology = Type::Unknown;
			dca = -1.;
			cosPointing = -2.;
			decayRadius = -1.;
			mass = 0.;

			deuteron_Charge       = -2;
			deuteron_P            = 0.;
			deuteron_PTPC         = 0.;
			deuteron_Px           = 0.;
			deuteron_Py           = 0.;
			deuteron_Pz           = 0.;
			deuteron_Eta          = -2.;
			deuteron_TPCNcls      = 0;
			deuteron_Chi2perNDF   = -1.;
			deuteron_TPCsignal    = -1.;
			deuteron_NSigmaTPC    = -20.;
			deuteron_NSigmaTOF    = -20.;
			deuteron_TPCExpSignal = -1.;
			deuteron_TOFExpSignal = -1.;

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
		kMassDeuteron,
		kMassTriton,
		kMassHelium3,
		kMassMax  // Number of enum entries
	};
	static const Double_t	fgkMass[];              //! masses

public:
	AliAnalysisTaskLambdaNRun2();
	AliAnalysisTaskLambdaNRun2(const char *name);
	virtual                 ~AliAnalysisTaskLambdaNRun2();

	virtual void            UserCreateOutputObjects();
	virtual void            UserExec(Option_t* option);
	virtual void            Terminate(Option_t* option);

private:
	bool                    LooseTrackCuts(AliAODTrack*);
	void                    FillEvent(AnalysisV0::Type, AliAODTrack*, AliAODTrack*);

	TTree*                  fOutputTree;            //! output tree
	AnalysisEvent*          fOutputEvent;           //! event class
	AliAODEvent*            fAOD;                   //! input event
	AliEventCuts            fEventCut;
	AliPIDResponse*         fPID;                   //! PID response
	TList*                  fOutputList;            //! output list

	AnalysisV0              fAnalysis_V0;
	AliAnalysisTaskLambdaNRun2(const AliAnalysisTaskLambdaNRun2&); // not implemented
	AliAnalysisTaskLambdaNRun2& operator=(const AliAnalysisTaskLambdaNRun2&); // not implemented

	ClassDef(AliAnalysisTaskLambdaNRun2, 1);
};

#endif
