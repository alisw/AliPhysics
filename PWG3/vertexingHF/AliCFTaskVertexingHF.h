#ifndef ALICFTASKVERTEXINGHF_H
#define ALICFTASKVERTEXINGHF_H
/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */ 

//-----------------------------------------------------------------------
// Class for HF corrections as a function of many variables and step 
// Author : C. Zampolli, CERN
//			D. Caffarri, Univ & INFN Padova caffarri@pd.infn.it
// Base class for HF Unfolding - agrelli@uu.nl
//-----------------------------------------------------------------------


#include "AliAnalysisTaskSE.h"
#include "AliCFVertexingHF2Prong.h"
#include "AliCFVertexingHF.h"

class TH1I;
class TParticle ;
class TFile ;
class TClonesArray ;
class AliCFManager;
class AliAODRecoDecay;
class AliAODRecoDecayHF2Prong;
class AliAODMCParticle;
class THnSparse;
class AliRDHFCuts;
class AliCFVertexingHF2Prong;
class AliCFVertexingHF3Prong;

class AliCFTaskVertexingHF: public AliAnalysisTaskSE {
public:
	
	enum {
	        kStepGeneratedLimAcc = 0,
		kStepGenerated       = 1,
		kStepAcceptance      = 2,
		kStepVertex          = 3,
		kStepRefit           = 4,
		kStepReconstructed   = 5,
		kStepRecoAcceptance  = 6,
		kStepRecoITSClusters = 7,
		kStepRecoPPR         = 8,
		kStepRecoPID         = 9
	};

	enum {
		kSnail = 0,    // slow configuration, all variables
		kCheetah = 1   // fast configuration, only a subset of variables
	};
	
	AliCFTaskVertexingHF();
	AliCFTaskVertexingHF(const Char_t* name, AliRDHFCuts* cuts);
	AliCFTaskVertexingHF& operator= (const AliCFTaskVertexingHF& c);
	AliCFTaskVertexingHF(const AliCFTaskVertexingHF& c);
 	virtual ~AliCFTaskVertexingHF();
	
	// ANALYSIS FRAMEWORK STUFF to loop on data and fill output objects
	void     UserCreateOutputObjects();
	void     UserExec(Option_t *option);
	void     Init();
	void     LocalInit() {Init();}
	void     Terminate(Option_t *);
	
	// UNFOLDING
	void     SetCorrelationMatrix(THnSparse* h) {fCorrelation=h;}
	void     SetAcceptanceUnf(Bool_t AcceptanceUnf) {fAcceptanceUnf = AcceptanceUnf;}
	Bool_t   GetAcceptanceUnf() const {return fAcceptanceUnf;}
	
	
	// CORRECTION FRAMEWORK RELATED FUNCTIONS
	void           SetCFManager(AliCFManager* io) {fCFManager = io;}   // global correction manager
	AliCFManager * GetCFManager()                 {return fCFManager;} // get corr manager
	
	// Setters (and getters) for the config macro
	void    SetFillFromGenerated(Bool_t flag) {fFillFromGenerated = flag;}
	Bool_t  GetFillFromGenerated() const {return fFillFromGenerated;}
	void    SetDecayChannel (Int_t decayChannel) {fDecayChannel = decayChannel;}
	Int_t   GetDecayChannel () {return fDecayChannel;}
	void     SetUseWeight(Bool_t useWeight){fUseWeight=useWeight;}
	Bool_t   GetUseWeight() const {return fUseWeight;}
	Double_t GetWeight(Float_t pt);
	Double_t dNdptFit(Float_t pt, Double_t* par);
 	
	void   SetDselection(UShort_t originDselection) {fOriginDselection=originDselection;}
	UShort_t GetDselection (){return fOriginDselection;}
	void SetSign(Char_t isSign) {fSign = isSign;}
	Char_t GetSign() {return fSign;}
	 
	void SetCentralitySelection(Bool_t centSelec = kTRUE) {fCentralitySelection = centSelec;}   
	Bool_t GetCentralitySelection() {return fCentralitySelection;} 

	void SetFakeSelection(Int_t fakeSel = 0) {fFakeSelection=fakeSel;}
	Int_t GetFakeSelection(){return fFakeSelection;}

	void SetRejectCandidateIfNotFromQuark(Bool_t opt){fRejectIfNoQuark=opt;}
	Bool_t GetRejectCandidateIfNotFromQuark(){return fRejectIfNoQuark;}

	void SetUseMCVertex(Bool_t opt){fUseMCVertex=opt;}
	Bool_t GetUseMCVertex(){return fUseMCVertex;}
	

	void SetKeepDsViaPhi(){fDsOption=1;}
	void SetKeepDsViaK0star(){fDsOption=2;}
	void SetKeepAllDs(){fDsOption=3;}

	Bool_t ProcessDs(Int_t returnCodeDs) const;

	void SetConfiguration(Int_t configuration) {(configuration == kSnail) ? Printf("Slow configuration chosen, all variables will be used!") : Printf("Fast configuration chosen, all variablesOnly pt, y, phi, ct, fake, z_vtx, centrality and multiplicity will be used!"); fConfiguration = configuration;} 
	Int_t GetConfiguration() const {return fConfiguration;} 
	
protected:
	AliCFManager   *fCFManager;   //  pointer to the CF manager
	TH1I *fHistEventsProcessed;   //! simple histo for monitoring the number of events processed
	THnSparse* fCorrelation;      //  response matrix for unfolding
	Int_t fCountMC;               //  MC particle found
	Int_t fCountAcc;              //  MC particle found that satisfy acceptance cuts
	Int_t fCountVertex;       //  Reco particle found that satisfy vertex constrained
	Int_t fCountRefit;        //  Reco particle found that satisfy kTPCrefit and kITSrefit
	Int_t fCountReco;             //  Reco particle found that satisfy cuts
	Int_t fCountRecoAcc;          //  Reco particle found that satisfy cuts in requested acceptance
	Int_t fCountRecoITSClusters;  //  Reco particle found that satisfy cuts in n. of ITS clusters
	Int_t fCountRecoPPR;          //  Reco particle found that satisfy cuts in PPR
	Int_t fCountRecoPID;          //Reco PID step 
	Int_t fEvents;                //  n. of events
	Int_t fDecayChannel;          // decay channel to configure the task
	Bool_t fFillFromGenerated;    //  flag to indicate whether data container should be filled with generated values also for reconstructed particles
	UShort_t fOriginDselection;      // flag to select D0 origins. 0 Only from charm 1 only from beauty 2 both from charm and beauty
	Bool_t fAcceptanceUnf;        //  flag for unfolding before or after cuts.
	AliRDHFCuts* fCuts;            // cuts
	Bool_t fUseWeight;             //flag to decide whether to use weights != 1 when filling the container or not
	Double_t fWeight;              //weight used to fill the container
	Int_t fNvar;                   // number of variables for the container
	TString fPartName;    // D meson name
	TString fDauNames;    // daughter in fin state
	Char_t fSign;                 // flag to decide wheter to keep D0 only (0), D0bar only (1), or both D0 and D0bar (2)
        Bool_t fCentralitySelection;  //flag to switch off the centrality selection
	Int_t  fFakeSelection;  //selection flag for fakes tracks 
	Bool_t fRejectIfNoQuark;  // flag to remove events not geenrated with PYTHIA
	Bool_t fUseMCVertex;  // flag to use MC vertex (useful when runnign in pp)
	Int_t  fDsOption;     // Ds decay option
	Int_t fConfiguration; // configuration (slow / fast) of the CF --> different variables will be allocated (all / reduced number)

	ClassDef(AliCFTaskVertexingHF,8); // class for HF corrections as a function of many variables
};

#endif
