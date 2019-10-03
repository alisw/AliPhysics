/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

#ifndef ALIJJTANA_H
#define ALIJJTANA_H

#include "../AliJConst.h"
#include <TH1D.h>
#include <TRandom3.h>

// jtAnalysis main class
// used in local and grid execution

class AliJCard;
class AliJJtHistograms;
class AliJJtCorrelations;
class AliJEventHeader;
class AliJEventPool;
class AliJEfficiency;
class AliJTrackCounter;
class AliJAcceptanceCorrection;

class TClonesArray;
class TF1;
class AliJRunTable;
class TRandom3;

using namespace std;

class AliJJtAna : public TObject
{
public:
  AliJJtAna(); // Default contructor
  AliJJtAna(Bool_t execLocal); // Constructor
  virtual ~AliJJtAna(); // Destructor
  AliJJtAna(const AliJJtAna& obj); // Copy constructor
  AliJJtAna& operator=(const AliJJtAna& obj); // Equal sign operator
  
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
  void SetTrackList( TClonesArray *inlist ) { finputList = inlist; } // Setter for Data Manager track list
  void SetTrigger( char const *p ) { fjtrigg = GetParticleType(p); } // Setter for trigger particle type
  void SetAssoc( char const *p ) { fjassoc = GetParticleType(p); } // Setter for associated particle type
  void SetInclusiveFile( const char *f ){ fInclusiveFile = f; } // Setter for inclusive histogram file
  void SetInputFile( char *f ) { finputFile = f; } // Setter for Data Manager configuration file

  //Event information to be set from Analysis task AliJJtTask
  void SetRunNumber(int runN) { fRunNumber= runN;}
  void SetCentrality(float cent) { fcent = cent;}
  void SetZVertex(double zvtx) { fZvert = zvtx;}

  double DeltaPhi(double phi1, double phi2); // Calculate deltaPhi from two phi values
  particleType  GetParticleType(char const *inchar); // Get particleType from string
  TClonesArray * GetInputList() const{return finputList;}

private:

  Bool_t fFirstEvent; // Flag for first event in the analysis
  double fRunNumber;
  double fcent;
  double fZvert;

  particleType fjtrigg; // Associated particle type
  particleType fjassoc; // Trigger particle type

  AliJCard *fcard; // JCard containing the binning information etc.
  char *finputFile; //! Name of the Data Manager initialization file for local analysis
  TString fInclusiveFile; // File for inclusive distributions

  TRandom3 *frandom; // Random number generator

  Int_t fevt; // Event counter
  AliJJtHistograms *fhistos; //! Histogram container
  AliJJtCorrelations *fcorrelations; //! Correlation analysis details
  AliJAcceptanceCorrection *fAcceptanceCorrection; //! Class for acceptance correction
  AliJEventPool *fassocPool; //! Pool of associated particles for event mixing
  TClonesArray *ftriggList; //! List of trigger particles
  TClonesArray *fassocList; //! List of associated particles
  TClonesArray *finputList; // List of particles currently used as input

  bool fbTriggCorrel; //! Flag for triggered correlation
  bool fbLPCorrel; //! Flag for leading particle correlation
  double fMinimumPt; //!  Minimum pT value for a particle to be still accepted to analysis
  bool fbLPSystematics; //! false = regular run, true = do systematic error estimate for missed leading particles
  bool fMCTruthRun; //! false = regular run, true = read particles from MC particle list

  AliJEfficiency *fEfficiency; // Efficience class
  AliJRunTable *fRunTable; // Run table
  int fHadronSelectionCut; // Used hadron selection cut

  ClassDef(AliJJtAna, 1); // ClassDef needed if inheriting from TObject

};

#endif

