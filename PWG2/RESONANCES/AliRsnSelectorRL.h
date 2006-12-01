/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 **************************************************************************/

//-------------------------------------------------------------------------
//  Class AliRsnSelectorRL
// ------------------------
// Reader for conversion of ESD output into the internal format
// used for resonance study.
// ---
// author: A. Pulvirenti             (email: alberto.pulvirenti@ct.infn.it)
// ---
// adapted for TSelector compliance
// by    : R. Vernet                 (email: renaud.vernet@cern.ch)
//-------------------------------------------------------------------------

#ifndef ALIRSNSELECTORRL_H
#define ALIRSNSELECTORRL_H

#include "AliPID.h"
#include "AliSelectorRL.h"

class TH3D;
class TH1D;
class TOrdCollection;
class TTree;
class TBranch;
class AliRunLoader;

class AliRsnSelectorRL : public AliSelectorRL
{
public:

	enum EPIDMethod {
    	kNoPID = 0,                  // no PID performed
		kPerfectPID = 1,             // use Kinematics to simulate a 100% PID efficiency
		kESDPID = 2,                 // use experimental ESD weights
		kPerfectPIDwithRecSign = 3   // get particle type from Kine and charge sign from reconstruction
	};
	
	                      AliRsnSelectorRL(TTree *tree = 0);
						  AliRsnSelectorRL(const AliRsnSelectorRL&);
	virtual              ~AliRsnSelectorRL();
	AliRsnSelectorRL&     operator=(const AliRsnSelectorRL&);
	
	// TSelector-inherited member functions
	virtual Int_t   Version() const {return 1;}
	virtual void    Begin(TTree *tree);
	virtual void    SlaveBegin(TTree *tree);
	virtual void    Init(TTree *tree);
//	virtual Bool_t  Notify();
	virtual Bool_t  Process(Long64_t entry);
	virtual void    SetOption(const char *option) {fOption = option;}
	virtual void    SetObject(TObject *obj) {fObject = obj;}
	virtual void    SetInputList(TList *input) {fInput = input;}
	virtual TList  *GetOutputList() const {return fOutput;}
	virtual void    SlaveTerminate();
	virtual void    Terminate();
	
	// Parameter/flag setting
	void            SetDebugFlag(Bool_t flag)     {fDebugFlag=flag;}
	void            SetOutputFile(char* file)     {fOutputPath=new TString(file);}
	void            SetStoreKineInfo(Bool_t flag) {fStoreKineInfo=flag;}
	void            SetCheckITSRefit(Bool_t b)    {fCheckITSRefit=b;}
	void            SetRejectFakes(Bool_t b)      {fRejectFakes=b;}
	void            SetCopyMomentum(Bool_t b)     {fCopyMomentum=b;}
	
	// Other
	void            Clear(const Option_t *option = "");
	Double_t*       GetPIDprobabilities(AliRsnDaughter track);
	void            Identify(AliRsnDaughter &track);
	void            SetMaxRadius(Double_t value) {fMaxRadius=value;}
	void            SetPIDMethod(AliRsnSelectorRL::EPIDMethod pm) {fPIDMethod=pm;}
	void            SetPriorProbabilities(Double_t *prior);
	void            SetPriorProbability(AliPID::EParticleType type, Double_t p);
	void            SetProbabilityThreshold(Double_t p) {fProbThreshold=p;}
	void            SetPtLimit4PID(Double_t value) {fPtLimit4PID=value;}
	
 protected:
 
 	// Parameters/flags
	TString*      fOutputPath;
	Bool_t        fDebugFlag;
	Bool_t        fStoreKineInfo;
	Bool_t        fCheckITSRefit;
	Bool_t        fRejectFakes;
	Bool_t        fCopyMomentum;
	
	// Workaround for AliSelectorRL:
	Bool_t        fIsRunLoaderOpen;
	
	// IO tree
	TTree*        fRsnEventTree;
	AliRsnEvent*  fRsnEvent;
	TBranch*      fRsnEventBranch;
	
	// PID
	EPIDMethod    fPIDMethod;                //  PID method
	Double_t      fPrior[AliPID::kSPECIES];  //  prior probabilities (in case of REAL pid)
	Double_t      fPtLimit4PID;              //  maximum transverse momentum to accept realistic PID
	Double_t      fProbThreshold;            //  minimum required probability to accept realistic PID
	Double_t      fMaxRadius;                //  maximum allowed distance from primary vertex
	
	// Functions
	AliPID::EParticleType FindType(Int_t pdg);
	
	// Rsn event reader implementation
	ClassDef(AliRsnSelectorRL,1)
};

#endif

