#ifndef AliEventPoolMuon_H
#define AliEventPoolMuon_H

/* $Id$ */ 

#include "AliVEventPool.h"
#include "AliRunTagCuts.h"
#include "AliLHCTagCuts.h"
#include "AliDetectorTagCuts.h"
#include "AliEventTagCuts.h"
#include "AliTagAnalysis.h"

// Realisation of an AliVEventPool via
// on the flight generation of the bin using AliTagAnalysis.
// Created expanding AliEventPoolOTF class functionalities
//
// Authors Alessandro De Falco and Antonio Uras, INFN Cagliari
// alessandro.de.falco@ca.infn.it  antonio.uras@ca.infn.it

//======================================================================================================

class AliEventPoolMuon : public AliVEventPool {

 public:
  AliEventPoolMuon();
  AliEventPoolMuon(const Char_t *name, const Char_t *title = "AOD");
  
  virtual ~AliEventPoolMuon() {;}

  virtual TChain* GetNextChain();
  virtual void  GetCurrentBin(Float_t* /*bin*/);
  virtual Int_t GetDimension();
  virtual void  Init();
  virtual void  SetMultiplicityRange(Int_t min, Int_t max, Int_t step)
  { fMultiplicityMin = min; fMultiplicityMax = max; fMultiplicityStep = step; }
  virtual void  SetNFWMuonRange(Int_t min, Int_t max, Int_t step)
  { fNFWMuonMin = min; fNFWMuonMax = max; fNFWMuonStep = step; }
  virtual void  SetPrimaryVertexZRange(Int_t min, Int_t max, Int_t step)
  { fPrimaryVertexZMin = min; fPrimaryVertexZMax = max; fPrimaryVertexZStep = step; }
  virtual Double_t GetMeanPrimaryVertexZ() { return fPrimaryVertexZ + 0.5*fPrimaryVertexZStep; }
  
  void SetTagDirectory(const Char_t *dirname) {fTagDirectory = dirname;};
  virtual Int_t BinNumber() const {return fBinNumber;}
  
 private:
  AliEventPoolMuon(const AliEventPoolMuon& obj);
  AliEventPoolMuon& operator=(const AliEventPoolMuon& other);

 protected:

  AliTagAnalysis      *fTagAnalysis;  // Pointer to tag analysis
  AliRunTagCuts       *fRunCuts;      // Run      cuts
  AliLHCTagCuts       *fLHCCuts;      // LHC      cuts
  AliDetectorTagCuts  *fDetectorCuts; // Detector cuts
  AliEventTagCuts     *fEventCuts;    // Event    cuts

  const Char_t        *fTagDirectory; // Directory with local tag files

  Int_t  fMultiplicityMin;       // Minimum multiplicity
  Int_t  fMultiplicityMax;       // Maximum multiplicity
  Int_t  fMultiplicityStep;      // Multiplicity step-size 
  Int_t  fMultiplicity;          // Minimum multiplicity for the current bin
  	
  Int_t  fNFWMuonMin;            // Minimum NFWMuon
  Int_t  fNFWMuonMax;            // Maximum NFWMuon
  Int_t  fNFWMuonStep;           // NFWMuon step-size 
  Int_t  fNFWMuon;               // Minimum NFWMuon for the current bin
  	
  Double_t  fPrimaryVertexZMin;  // Minimum PrimaryVertexZ
  Double_t  fPrimaryVertexZMax;  // Maximum PrimaryVertexZ
  Double_t  fPrimaryVertexZStep; // PrimaryVertexZ step-size
  Double_t  fPrimaryVertexZ;     // Minimum PrimaryVertexZ for the current bin
  
  Int_t  fBinNumber;             // Current bin number
  
  ClassDef(AliEventPoolMuon, 0); 

};

//======================================================================================================
 
#endif
