/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

// A class for calculating acceptance correction


#ifndef ALIJACCEPTANCECORRECTION_H
#define ALIJACCEPTANCECORRECTION_H

#ifndef ROOT_TObject
#include <TObject.h>
#endif

#include "AliJHistManager.h"
#include "AliJCard.h"
#include "AliJConst.h"

class AliJAcceptanceCorrection{

public:


  AliJAcceptanceCorrection();                     // default constructor
  AliJAcceptanceCorrection(AliJCard *inputCard);  // constructor with card definition
  AliJAcceptanceCorrection(const AliJAcceptanceCorrection& a); // copy constructor
  ~AliJAcceptanceCorrection(){;}    //destructor

  void ReadInclusiveHistos(const char *inclusFileName);     // Read acceptance histograms from input file
  double GetAcceptanceCorrectionTriangle(double deltaEta);  // Calculate acceptance correction from deltaEta triangle
  double GetAcceptanceCorrectionInclusive(double deltaEta, int centralityBin, int triggerBin, int assocBin, int assocType);  // Calculate the acceptance correction from a one dimensional input histogram

  AliJAcceptanceCorrection& operator=(const AliJAcceptanceCorrection& a);


private:
  void NormalizeAcceptanceTraditional(AliJTH1D &acceptanceHisto, corrType assocType);
  
  AliJCard *fCard;                    // Card containing analysis details
  AliJTH1D fDEtaNearAcceptance;       // DeltaEta acceptance histogram for the near side
  AliJTH1D fDEta3DNearAcceptance;     // DeltaEta acceptance histogram for the 3D near side
  AliJTH2D fDEtaDPhiNearAcceptance;   // DeltaEta DeltaPhi acceptance histogram for the near side
  AliJTH2D fDEtaDPhi3DNearAcceptance; // DeltaEta DeltaPhi acceptance histogram for the 3D near side
  
  bool fInclusiveHistogramsNotFound;  // Flag if inclusive acceptance histograms are in the input file

  ClassDef(AliJAcceptanceCorrection,1)
};

#endif
