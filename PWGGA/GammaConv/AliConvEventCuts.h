#ifndef ALICONVEVENTCUTS_H
#define ALICONVEVENTCUTS_H

// Class handling all kinds of selection cuts for Gamma Conversion analysis
// Authors: Svein Lindal, Daniel Lohner                                    *

#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliVEvent.h"
#include "AliMCEvent.h"
#include "AliVTrack.h"
#include "AliStack.h"
#include "AliAnalysisCuts.h"
#include "TH1F.h"
#include "TF1.h"
#include "AliAnalysisUtils.h"
#include "AliAnalysisManager.h"
#include "TRandom3.h"

class AliESDEvent;
class AliAODEvent;
class TH1F;
class TH2F;
class TF1;
class AliAnalysisCuts;
class iostream;
class TList;
class AliAnalysisManager;
class AliAODMCParticle;

using namespace std;

class AliConvEventCuts : public AliAnalysisCuts {
		
	public: 
		enum cutIds {
			kisHeavyIon,                  
			kCentralityMin,               
			kCentralityMax,               
			kSelectSpecialTriggerAlias,                 
			kSelectSubTriggerClass,             
			kremovePileUp,                
			kExtraSignals, 
			kNCuts
		};

		Bool_t SetCutIds(TString cutString); 
		Int_t fCuts[kNCuts];
		Bool_t SetCut(cutIds cutID, Int_t cut);
		Bool_t UpdateCutString();
		static const char * fgkCutNames[kNCuts];

		Bool_t InitializeCutsFromCutString(const TString analysisCutSelection);
		void SelectCollisionCandidates(UInt_t offlineTriggerMask = AliVEvent::kAny) {
			fOfflineTriggerMask = offlineTriggerMask;
			fTriggerSelectedManually = kTRUE;
		}
		void SelectSpecialTrigger(UInt_t offlineTriggerMask = AliVEvent::kAny, TString TriggerClassName = "AliVEvent::kAny" ) {
			fOfflineTriggerMask = offlineTriggerMask;
			fSpecialTriggerName = TriggerClassName;
			cout << fSpecialTriggerName.Data() << endl;
			
		}   

		void SetAcceptedHeader(TList *HeaderList){fHeaderList = HeaderList;}   
		TString *GetFoundHeader(){return fGeneratorNames;}

		Int_t GetEventQuality(){return fEventQuality;}
		Bool_t GetIsFromPileup(){return fRemovePileUp;}
		
		AliConvEventCuts(const char *name="EventCuts", const char * title="Event Cuts");
		AliConvEventCuts(const AliConvEventCuts&);
		AliConvEventCuts& operator=(const AliConvEventCuts&);

		virtual ~AliConvEventCuts();                            //virtual destructor

// 		static AliConvEventCuts * GetStandardCuts2010PbPb();
// 		static AliConvEventCuts * GetStandardCuts2010pp();

		virtual Bool_t IsSelected(TObject* /*obj*/){return kTRUE;}
		virtual Bool_t IsSelected(TList* /*list*/) {return kTRUE;}

		TString GetCutNumber();
		
		void GetCentralityRange(Double_t range[2]){range[0]=10*fCentralityMin;range[1]=10*fCentralityMax;};
		
		// Cut Selection
		Bool_t EventIsSelected(AliVEvent *fInputEvent, AliVEvent *fMCEvent);
		Int_t IsEventAcceptedByCut(AliConvEventCuts *ReaderCuts, AliVEvent *InputEvent, AliMCEvent *MCEvent, Int_t isHeavyIon);
			
		void PrintCuts();
		void PrintCutsWithValues();

		void InitCutHistograms(TString name="",Bool_t preCut = kTRUE);
		void SetFillCutHistograms(TString name="",Bool_t preCut = kTRUE){if(!fHistograms){InitCutHistograms(name,preCut);};}
		TList *GetCutHistograms(){return fHistograms;}
		
		///Cut functions
		Int_t IsParticleFromBGEvent(Int_t index, AliStack *MCStack, AliVEvent *InputEvent = 0x0);
		void GetNotRejectedParticles(Int_t rejection, TList *HeaderList, AliVEvent *MCEvent);
		void SetUseReweightingWithHistogramFromFile( Bool_t pi0reweight=kTRUE, Bool_t etareweight=kFALSE, Bool_t k0sreweight=kFALSE, TString path="$ALICE_ROOT/PWGGA/GammaConv/MCSpectraInput.root", 
														TString histoNamePi0 = "", TString histoNameEta = "", TString histoNameK0s = "",
														TString fitNamePi0 = "", TString fitNameEta = "", TString fitNameK0s ="" ) {
			AliInfo(Form("enabled reweighting for: pi0 : %i, eta: %i, K0s: %i",pi0reweight, etareweight, k0sreweight));
			fDoReweightHistoMCPi0 = pi0reweight; 
			fDoReweightHistoMCEta = etareweight; 
			fDoReweightHistoMCK0s = k0sreweight; 
			fPathTrFReweighting=path;
			fNameHistoReweightingPi0 =histoNamePi0;
			fNameHistoReweightingEta =histoNameEta;
			fNameHistoReweightingK0s =histoNameK0s; 
			fNameFitDataPi0 =fitNamePi0;
			fNameFitDataEta =fitNameEta;
			fNameFitDataK0s =fitNameK0s;      
		} 
		void LoadReweightingHistosMCFromFile ();
		void SetAddedSignalPDGCode(Int_t addedSignalPDGcode) {fAddedSignalPDGCode = addedSignalPDGcode;}

		// Event Cuts
		Bool_t IsCentralitySelected(AliVEvent *fInputEvent, AliVEvent *fMCEvent = NULL);
		Double_t GetCentrality(AliVEvent *event);
		Int_t GetNumberOfContributorsVtx(AliVEvent *event);
		Bool_t VertexZCut(AliVEvent *fInputEvent);
		Bool_t IsTriggerSelected(AliVEvent *fInputEvent);
		Bool_t HasV0AND(){return fHasV0AND;}
		Bool_t IsSDDFired(){return fIsSDDFired;}
		Int_t IsSpecialTrigger(){return fSpecialTrigger;}
		TString GetSpecialTriggerName(){return fSpecialTriggerName;}

		// Request Flags
		Int_t IsHeavyIon(){return fIsHeavyIon;}
		Float_t GetWeightForMeson(TString period, Int_t index, AliStack *MCStack, AliVEvent *InputEvent = 0x0);
		void SetPreSelectionCutFlag(Bool_t preSelFlag){fPreSelCut = preSelFlag;}   

		
		Int_t GetMultiplicityMethod(){return fMultiplicityMethod;}
		Int_t GetSignalRejection(){return fRejectExtraSignals;}
		Int_t GetNAcceptedHeaders(){return fnHeaders; }
		TString * GetAcceptedHeaderNames(){return fGeneratorNames;}
		Int_t * GetAcceptedHeaderStart(){return fNotRejectedStart;}
		Int_t * GetAcceptedHeaderEnd(){return fNotRejectedEnd;}
		TList* GetAcceptedHeader(){return fHeaderList;}
		
		
		// Eta shift Setting
		void SetEtaShift(Double_t etaShift) {
			fEtaShift = etaShift;
		}
		void SetEtaShift(TString pPbOrPbp) {
			Double_t etaShift = 0.0;
			if(!pPbOrPbp.CompareTo("pPb"))      etaShift = -0.465;
			else if(!pPbOrPbp.CompareTo("Pbp")) etaShift =  0.465;
			
			fEtaShift = etaShift;
		}
		Double_t GetEtaShift() {return fEtaShift;}
		Bool_t GetDoEtaShift(){return fDoEtaShift;}
		void DoEtaShift(Bool_t doEtaShift){fDoEtaShift = doEtaShift;}
		void GetCorrectEtaShiftFromPeriod(TString periodName);
	
		Bool_t SetIsHeavyIon(Int_t isHeavyIon);
		Bool_t SetCentralityMax(Int_t centralityBin);
		Bool_t SetCentralityMin(Int_t centralityBin);
		Bool_t SetRemovePileUp(Int_t removePileUp);  
		Bool_t SetMultiplicityMethod(Int_t multiplicityMethod);
		Bool_t SetSelectSpecialTrigger(Int_t selectSpecialTrigger);
		Bool_t SetSelectSubTriggerClass (Int_t selectSpecialSubTriggerClass);
		Bool_t SetRejectExtraSignalsCut(Int_t extraSignal);
		void SetV0ReaderName(TString name) {fV0ReaderName = name;}
		
	protected:
		TList 				*fHistograms;
		TList 				*fHeaderList;

		Int_t 				fEventQuality; 							// EventQuality
		//cuts
		Int_t 				fIsHeavyIon;							// flag for heavy ion
		Int_t 				fDetectorCentrality;					// centrality detecotor V0M or CL1
		Int_t 				fModCentralityClass;					// allows to select smaller centrality classes
		Double_t 			fMaxVertexZ;							// max z offset of vertex
		Int_t 				fCentralityMin;							// centrality selection lower bin value
		Int_t 				fCentralityMax;							// centrality selection upper bin value
		Int_t 				fMultiplicityMethod;					// selected multiplicity method
		Int_t 				fSpecialTrigger;						// flag
		Int_t 				fSpecialSubTrigger;						// flag
		Bool_t 				fRemovePileUp;							// flag
		Int_t 				fRejectExtraSignals;					//
		UInt_t 				fOfflineTriggerMask;   					// Task processes collision candidates only
		Bool_t 				fHasV0AND; 								// V0AND Offline Trigger
		Bool_t 				fIsSDDFired; 							// SDD FIRED to select with SDD events
		TRandom3		 	fRandom; 								//
		Int_t 				fnHeaders; 								// Number of Headers
		Int_t 				*fNotRejectedStart; 					//[fnHeaders]
		Int_t 				*fNotRejectedEnd; 						//[fnHeaders]
		TString 			*fGeneratorNames; 						//[fnHeaders]
		TObjString 			*fCutString; 							// cut number used for analysis
		AliAnalysisUtils 	*fUtils;
		Double_t			fEtaShift;
		Bool_t 				fDoEtaShift;							// Flag for Etashift
		Bool_t 				fDoReweightHistoMCPi0; 					// Flag for reweighting Pi0 input with histogram
		Bool_t				fDoReweightHistoMCEta;					// Flag for reweighting Eta input with histogram
		Bool_t 				fDoReweightHistoMCK0s;					// Flag for reweighting K0s input with histogram
		TString 			fPathTrFReweighting;					// Path for file used in reweighting
		TString 			fNameHistoReweightingPi0;				// Histogram name for reweighting Pi0
		TString 			fNameHistoReweightingEta;				// Histogram name for reweighting Eta
		TString 			fNameHistoReweightingK0s;				// Histogram name for reweighting K0s
		TString 			fNameFitDataPi0;						// Fit name for fit to spectrum of pi0s in Data
		TString 			fNameFitDataEta;						// Fit name for fit to spectrum of etas in Data
		TString 			fNameFitDataK0s;						// Fit name for fit to spectrum of k0s in Data
		// Histograms
		TH1F 				*fHistoEventCuts;						// bookkeeping for event selection cuts
		TH1F 				*hCentrality;							// centrality distribution for selected events
		TH2F 				*hCentralityVsNumberOfPrimaryTracks;	// centrality distribution for selected events
		TH1F 				*hVertexZ; 								// vertex z distribution for selected events
		TH1F 				*hTriggerClass; 						// fired offline trigger class
		TH1F 				*hTriggerClassSelected;					// selected fired offline trigger class
		TH1D 				*hReweightMCHistPi0;					// histogram input for reweighting Pi0
		TH1D 				*hReweightMCHistEta; 					// histogram input for reweighting Eta
		TH1D 				*hReweightMCHistK0s; 					// histogram input for reweighting K0s
		TF1  				*fFitDataPi0;							// fit to pi0 spectrum in Data
		TF1  				*fFitDataEta;							// fit to eta spectrum in Data
		TF1 				*fFitDataK0s;							// fit to K0s spectrum in Data
		Int_t 				fAddedSignalPDGCode;
		Bool_t 				fPreSelCut; 							// Flag for preselection cut used in V0Reader
		Bool_t 				fTriggerSelectedManually; 				// Flag for manual trigger selection
		TString 			fSpecialTriggerName; 					// Name of the Special Triggers
		TString 			fSpecialSubTriggerName; 				// Name of the Special Triggers
		Int_t 				fNSpecialSubTriggerOptions;
		TString 			fV0ReaderName;							// Name of V0Reader
		
		
	private:

		ClassDef(AliConvEventCuts,2)
};


#endif
