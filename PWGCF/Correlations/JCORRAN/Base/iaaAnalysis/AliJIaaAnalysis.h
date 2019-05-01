/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

#ifndef AliJIaaANALYSIS_H
#define AliJIaaANALYSIS_H

#include "../AliJDataManager.h"
#include "../AliJConst.h"
#include <TH1D.h>

// iaaAnalysis main class
// used in local and grid execution

class AliJCard;
class AliJIaaHistograms;
class AliJIaaCorrelations;
class AliJEventHeader;
class AliJEventPool;
class AliJRunHeader;
class AliJEfficiency;
class AliJTrackCounter;
class AliJAcceptanceCorrection;

class TClonesArray;
class TF1;
class AliJRunTable;
class TRandom3;

class AliJIaaAnalysis : public TObject
{
public:
	AliJIaaAnalysis(); // Default contructor
	AliJIaaAnalysis(Bool_t execLocal); // Constructor
	virtual ~AliJIaaAnalysis(); // Destructor
	AliJIaaAnalysis(const AliJIaaAnalysis& obj); // Copy constructor
	AliJIaaAnalysis& operator=(const AliJIaaAnalysis& obj); // Equal sign operator

	void Initialize() const; // Initializer
	void Init(){ Initialize(); } // Initializer
	void UserCreateOutputObjects(); // Output object creation
	void UserExec(); // Event by event functionality
	void Terminate(); // Closing formalities

	AliJIaaHistograms *GetHistos() { return fhistos; } // Getter for histogram container
	AliJIaaCorrelations *GetCorrelations() { return fcorrelations; } // Getter for correlation analysis
	AliJEventPool *GetAssocPool() { return fassocPool; } // Getter for associated particle pool
	AliJCard *GetCard() { return fcard; } // Getter for JCard

	void SetCard( AliJCard *c ) { fcard = c; } // Setter for JCrad
	void SetTrigger( const char* p ) { fjtrigg = GetParticleType(p); } // Setter for trigger particle type
	void SetAssoc( const char* p ) { fjassoc = GetParticleType(p); } // Setter for associated particle type
	void SetInclusiveFile( const char *f ){ fInclusiveFile = f; } // Setter for inclusive histogram file
	void SetInputFile( const char *f ) { finputFile = f; } // Setter for Data Manager configuration file

	void SetTrackList( TClonesArray *a ) { fdmg->SetTrackList( a ); } // Setter for Data Manager track list
	void SetMCTrackList( TClonesArray *a ) { fdmg->SetMCTrackList( a ); } // Setter for Data Manager MC track list
	void SetHeaderList( TClonesArray *a ) { fdmg->SetHeaderList( a ); } // Setter for Data Manager header list
	void SetRunHeader( AliJRunHeader *a ) { frunHeader = a; } // Setter for run header
	void SetRunInfoList( TList *a ) { fdmg->SetRunInfoList( a ); } // Setter for Data Manager run info list

	double DeltaPhi(double phi1, double phi2); // Calculate deltaPhi from two phi values
	particleType  GetParticleType(const char *inchar); // Get particleType from string

private:

	Bool_t fExecLocal; // Execution mode (local vs. grid)
	Bool_t fFirstEvent; // Flag for first event in the analysis

	particleType fjtrigg; // Associated particle type
	particleType fjassoc; // Trigger particle type

	AliJCard *fcard; // JCard containing the binning information etc.
	const char *finputFile; //! Name of the Data Manager initialization file for local analysis
	TString fInclusiveFile; // File for inclusive distributions

	Int_t fevt; // Event counter
	AliJIaaHistograms *fhistos; //! Histogram container
	AliJIaaCorrelations *fcorrelations; //! Correlation analysis details
	AliJAcceptanceCorrection *fAcceptanceCorrection; //! Class for acceptance correction
	AliJEventPool *fassocPool; //! Pool of associated particles for event mixing
	TClonesArray *fchargedHadronList; //! List of charged particles for correlation analysis
	TClonesArray *ftriggList; //! List of trigger particles
	TClonesArray *fassocList; //! List of associated particles
	TClonesArray *finputList; //! List of particles currently used as input

	AliJDataManager *fdmg; //! Pointer to Data Manager
	AliJEventHeader *feventHeader; //! Pointer to Event Header
	AliJRunHeader *frunHeader; //! Pointer to Run Header

	double fcent; //! Centrality percentile
	bool fbTriggCorrel; //! Flag for triggered correlation
	bool fbLPCorrel; //! Flag for leading particle correlation
	double fMinimumPt; //!  Minimum pT value for a particle to be still accepted to analysis
    bool fMCTruthRun; //! false = regular run, true = read particles from MC particle list
	Int_t fEventBC; //! Selector for some BC%4

	AliJEfficiency *fEfficiency; // Efficience class
	AliJRunTable *fRunTable; // Run table
	int fHadronSelectionCut; // Used hadron selection cut

	ClassDef(AliJIaaAnalysis, 1); // ClassDef needed if inheriting from TObject

};

#endif
























