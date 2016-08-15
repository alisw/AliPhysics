/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

#ifndef ALIJJTANALYSIS_H
#define ALIJJTANALYSIS_H

#include "../AliJDataManager.h"
#include "../AliJConst.h"
#include <TH1D.h>

// jtAnalysis main class
// used in local and grid execution

class AliJCard;
class AliJJtHistograms;
class AliJJtCorrelations;
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

class AliJJtAnalysis : public TObject
{
public:
  AliJJtAnalysis(); // Default contructor
  AliJJtAnalysis(Bool_t execLocal); // Constructor
  virtual ~AliJJtAnalysis(); // Destructor
  AliJJtAnalysis(const AliJJtAnalysis& obj); // Copy constructor
  AliJJtAnalysis& operator=(const AliJJtAnalysis& obj); // Equal sign operator
  
  void Initialize() const; // Initializer
  void Init(){ Initialize(); } // Initializer
  void UserCreateOutputObjects(); // Output object creation
  void UserExec(); // Event by event functionality
  void Terminate(); // Closing formalities
  
  AliJJtHistograms *GetHistos() { return fhistos; } // Getter for histogram container
  AliJJtCorrelations *GetCorrelations() { return fcorrelations; } // Getter for correlation analysis
  AliJEventPool *GetAssocPool() { return fassocPool; } // Getter for associated particle pool
  AliJCard *GetCard() { return fcard; } // Getter for JCard
  
  void SetCard( AliJCard *c ) { fcard = c; } // Setter for JCrad
  void SetTrigger( char* p ) { fjtrigg = GetParticleType(p); } // Setter for trigger particle type
  void SetAssoc( char* p ) { fjassoc = GetParticleType(p); } // Setter for associated particle type
  void SetInclusiveFile( const char *f ){ fInclusiveFile = f; } // Setter for inclusive histogram file
  void SetInputFile( char *f ) { finputFile = f; } // Setter for Data Manager configuration file
  
  void SetTrackList( TClonesArray *a ) { fdmg->SetTrackList( a ); } // Setter for Data Manager track list
  void SetPhotonList( TClonesArray *a ) { fdmg->SetPhotonList( a ); } // Setter for Data Manager photon list
  void SetCaloCellList( TClonesArray *a ) { fdmg->SetCaloCellList( a ); } // Setter for Calo Cell track list
  void SetMCTrackList( TClonesArray *a ) { fdmg->SetMCTrackList( a ); } // Setter for Data Manager MC track list
  void SetHeaderList( TClonesArray *a ) { fdmg->SetHeaderList( a ); } // Setter for Data Manager header list
  void SetRunHeader( AliJRunHeader *a ) { frunHeader = a; } // Setter for run header
  void SetRunInfoList( TList *a ) { fdmg->SetRunInfoList( a ); } // Setter for Data Manager run info list
  
  double DeltaPhi(double phi1, double phi2); // Calculate deltaPhi from two phi values
  particleType  GetParticleType(char *inchar); // Get particleType from string
  
private:
  
  Bool_t fExecLocal; // Execution mode (local vs. grid)
  Bool_t fFirstEvent; // Flag for first event in the analysis
  
  particleType fjtrigg; // Associated particle type
  particleType fjassoc; // Trigger particle type
  
  AliJCard *fcard; // JCard containing the binning information etc.
  char *finputFile; //! Name of the Data Manager initialization file for local analysis
  TString fInclusiveFile; // File for inclusive distributions
  
  Int_t fevt; // Event counter
  AliJJtHistograms *fhistos; //! Histogram container
  AliJJtCorrelations *fcorrelations; //! Correlation analysis details
  AliJAcceptanceCorrection *fAcceptanceCorrection; //! Class for acceptance correction
  AliJEventPool *fassocPool; //! Pool of associated particles for event mixing
  TClonesArray *fphotonList; //! List of photons for photon analysis
  TClonesArray *fchargedHadronList; //! List of charged particles for correlation analysis
  TClonesArray *fpizeroList; //! List of pi zeros for pi zero analysis
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
  
  Int_t fEventBC; //! Selector for some BC%4
  
  AliJEfficiency *fEfficiency; // Efficience class
  AliJRunTable *fRunTable; // Run table
  int fHadronSelectionCut; // Used hadron selection cut
  
  ClassDef(AliJJtAnalysis, 1); // ClassDef needed if inheriting from TObject
  
};

#endif
























