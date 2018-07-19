/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

#ifndef AliJIaaANA_H
#define AliJIaaANA_H

#include "../AliJDataManager.h"
#include "../AliJConst.h"
#include <TH1D.h>

// iaaAnalysis main class
// used in local and grid execution

class AliJCard;
class AliJIaaHistograms;
class AliJIaaCorrelations;
class AliJEventPool;
class AliJEfficiency;
class AliJTrackCounter;
//class AliJAcceptanceCorrection;

class TClonesArray;
class TF1;
class AliJRunTable;
class TRandom3;

class AliJIaaAna : public TObject
{
public:
	AliJIaaAna(); // Default contructor
	AliJIaaAna(Bool_t execLocal); // Constructor
	virtual ~AliJIaaAna(); // Destructor
	AliJIaaAna(const AliJIaaAna& obj); // Copy constructor
	AliJIaaAna& operator=(const AliJIaaAna& obj); // Equal sign operator

	void Initialize() const; // Initializer
	void Init(){ Initialize(); } // Initializer
	void UserCreateOutputObjects(); // Output object creation
	void UserExec(); // Event by event functionality
	void Terminate(); // Closing formalities

	AliJIaaHistograms *GetHistos() { return fhistos; } // Getter for histogram container
	AliJIaaCorrelations *GetCorrelations() { return fcorrelations; } // Getter for correlation analysis
	AliJEventPool *GetAssocPool() { return fassocPool; } // Getter for associated particle pool
	AliJCard *GetCard() { return fcard; } // Getter for JCard

	void SetTrackList( TClonesArray *inlist ) { finputList = inlist; } // Setter for Data Manager track list
	void SetCard( AliJCard *c ) { fcard = c; } // Setter for JCrad
	void SetTrigger( char const *p ) { fjtrigg = GetParticleType(p); } // Setter for trigger particle type
	void SetAssoc( char const *p ) { fjassoc = GetParticleType(p); } // Setter for associated particle type
	void SetInclusiveFile( const char *f ){ fInclusiveFile = f; } // Setter for inclusive histogram file
	void SetInputFile( char *f ) { finputFile = f; } // Setter for Data Manager configuration file

	//Event information to be set from Analysis task AliJJtTask
	void SetRunNumber(int runN) { fRunNumber= runN;}
	void SetCentrality(float cent) { fcent = cent;}
	void SetZVertex(double zvtx) { fZvert = zvtx;}
	void SetEventPlane(double ep2) { fPsi2 = ep2;}
	double GetPhiS2(double phit, double ep);
	void SetEPmin(double min) { fEPmin = min;}
	void SetEPmax(double max) { fEPmax = max;}
	Bool_t AccecptEPBins(double phis, double min, double max);
	void SetEnableEP(Bool_t enable) { fenableEP = enable;}
	void RunCorrelations(TClonesArray *triggList, TClonesArray *assoList, int noAllTriggTracks, int cbin, int zbin); 

	double DeltaPhi(double phi1, double phi2); // Calculate deltaPhi from two phi values
	particleType  GetParticleType(char const *inchar); // Get particleType from string
	TClonesArray * GetInputList() const{return finputList;}

private:
	Bool_t fFirstEvent; // Flag for first event in the analysis
	double fRunNumber;
	double fcent;
	double fZvert;
	Bool_t fenableEP;
	double fPsi2;
	double fEPmin;
	double fEPmax;

	particleType fjtrigg; // Associated particle type
	particleType fjassoc; // Trigger particle type

	AliJCard *fcard; // JCard containing the binning information etc.
	char *finputFile; //! Name of the Data Manager initialization file for local analysis
	TString fInclusiveFile; // File for inclusive distributions

	Int_t fevt; // Event counter
	AliJIaaHistograms *fhistos; //! Histogram container
	AliJIaaCorrelations *fcorrelations; //! Correlation analysis details
	//AliJAcceptanceCorrection *fAcceptanceCorrection; //! Class for acceptance correction
	AliJEventPool *fassocPool; //! Pool of associated particles for event mixing
	TClonesArray *ftriggList; //! List of trigger particles
	TClonesArray *fassocList; //! List of associated particles
	TClonesArray *finputList; //! List of particles currently used as input

	bool fbTriggCorrel; //! Flag for triggered correlation
	bool fbLPCorrel; //! Flag for leading particle correlation
	double fMinimumPt; //!  Minimum pT value for a particle to be still accepted to analysis
	bool fMCTruthRun; //! false = regular run, true = read particles from MC particle list
	Int_t fEventBC; //! Selector for some BC%4

	AliJEfficiency *fEfficiency; // Efficience class
	AliJRunTable *fRunTable; // Run table
	int fHadronSelectionCut; // Used hadron selection cut

	ClassDef(AliJIaaAna, 1); // ClassDef needed if inheriting from TObject

};

#endif
























