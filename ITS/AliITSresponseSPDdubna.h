#ifndef ALIITSRESPONSESPDDUBNA_H
#define ALIITSRESPONSESPDDUBNA_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TString.h>
#include <TRandom.h>

#include "AliITSresponse.h"

//----------------------------------------------
// ITS response class for SPD
class AliITSresponseSPDdubna : public AliITSresponse {
 public:
  AliITSresponseSPDdubna();
  virtual ~AliITSresponseSPDdubna(){}// destructror
  //
  // set noise and baseline
  virtual void SetNoiseParam(Float_t n=200., Float_t b=0.) {
    fNoise=n; fBaseline=b;}
  // Returns just noise value
  virtual Float_t GetNoise() const {return fNoise;} 
  //Returns just baseline value
  virtual Float_t GetBaseline() const {return fBaseline;}
  virtual void GetNoiseParam(Float_t &n, Float_t &b) {
    n=fNoise; b=fBaseline;} // get noise and baseline
  // Applies a randome noise and addes the baseline
  virtual Float_t ApplyBaselineAndNoise() const {return fBaseline+
                                                 fNoise*gRandom->Gaus();}
  // Zero-suppression option threshold
  virtual void SetThreshold(Int_t th=2000) {fThreshold=th;}
  // Get zero-suppression threshold
  virtual Float_t GetThreshold() const {return fThreshold;}
  // Sets the coupling probabilities for columns and rows
  virtual void SetCouplings(Float_t ccol=0.,Float_t crow=0.){fCouplCol=ccol;
                                                       fCouplRow=crow;}
  // Gets the coupling probabilities for columns and rows
  virtual void GetCouplings(Float_t &ccol,Float_t &crow){ccol=fCouplCol;
                                                         crow=fCouplRow;}
  // Retruns the fraction of dead pixels
  virtual Float_t GetFractionDeadPixels() const {return fDeadPixels;}
  // Sets the fraction of dead pixels
  virtual void SetFractionDeadPixels(Float_t f=0.01){fDeadPixels = f;}
  // Returns a logical if a pixels is dead or not
  virtual Bool_t IsPixelDead(Int_t mod,Int_t ix,Int_t iz) const ;
  // Sets the Type of data -real or simulated
  virtual void SetDataType(const char *data="simulated") {fDataType=data;}
  // Get data type, real or simulated.
  virtual const char* DataType() const {return fDataType.Data();}
  // Prints out the content of this class in ASCII format.
  virtual void Print(ostream *os);
  // Reads in the content of this class in the format of Print
  virtual void Read(istream *is);

  ClassDef(AliITSresponseSPDdubna,2) // SPD response
    ;
 protected:
  Float_t fNoise;           // Noise value
  Float_t fBaseline;        // Baseline value
  Float_t fCouplCol;        // Coupling Probability along columns
  Float_t fCouplRow;        // Coupling Probability along rows
  Float_t fThreshold;       // Zero-Suppression threshold
  Float_t fDeadPixels;      // Fraction of dead pixels.
  TString fDataType;        // Type of data - real or simulated
};
// Input and output function for standard C++ input/output.
ostream& operator<<(ostream &os,AliITSresponseSPDdubna &source);
istream& operator>>(istream &os,AliITSresponseSPDdubna &source);

#endif
