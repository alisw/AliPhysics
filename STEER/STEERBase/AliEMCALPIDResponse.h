#ifndef AliEMCALPIDResponse_h
#define AliEMCALPIDResponse_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliEMCALPIDResponse                                                  //
//                                                                      //
// EMCAL class to perfom PID                                            //
// This is a prototype and still under development                      //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "AliPID.h"
class TF1;

class AliEMCALPIDResponse: public TObject 
{
public : 
    AliEMCALPIDResponse();    //ctor
    AliEMCALPIDResponse( const AliEMCALPIDResponse& other);                //copy ructor
    AliEMCALPIDResponse &operator=( const AliEMCALPIDResponse& other);     //assignment operator

    virtual ~AliEMCALPIDResponse();     //dtor
  

    // Getters
    Int_t     GetPtBin(Float_t pt) const;
    Int_t     GetNPtBins() const {return fNptBins;};

    Float_t   GetLowEoP() const {return fLowEoP;};
    Float_t   GetHighEoP() const {return fHighEoP;};

    Double_t  GetExpectedSignal( Float_t pt, AliPID::EParticleType n,  Int_t charge) const;
    Double_t  GetExpectedSigma ( Float_t pt, AliPID::EParticleType n,  Int_t charge) const;
    Double_t  GetNumberOfSigmas( Float_t pt,  Float_t eop, AliPID::EParticleType n,  Int_t charge) const;
    Double_t  GetExpectedNorm  ( Float_t pt, AliPID::EParticleType n,  Int_t charge) const;
    Double_t  GetLowProb       ( Float_t pt, AliPID::EParticleType n,  Int_t charge) const;
    Double_t  GetHighProb      ( Float_t pt, AliPID::EParticleType n,  Int_t charge) const;

    //Setters
    void   SetPtBoundary();
    void   SetParametrizations();

    // EMCAL probability -> should go to another place?
    Double_t ComputeEMCALProbability( Float_t pt, Float_t eop, Int_t charge, Double_t *pEMCAL) const;

protected:
  
private:

  TF1 *fNorm;                            // Gauss function for normalizing NON electron probabilities 

  static const Int_t fNptBins   = 6;     // number of momentum bins
  static const Float_t fLowEoP  = 0.5;   // lower E/p threshold for NON electrons
  static const Float_t fHighEoP = 1.5;   // upper E/p threshold for NON electrons

  Float_t fPtCutMin[fNptBins+1];                       // min values for pt bins
  Float_t fMeanEoP[2*AliPID::kSPECIES][fNptBins];      // mean value of E/p distribution (charge dependent)
  Float_t fSigmaEoP[2*AliPID::kSPECIES][fNptBins];     // mean value of E/p distribution (charge dependent)
  Float_t fProbLow[2*AliPID::kSPECIES][fNptBins];      // probability below E/p threshold for NON electrons (charge dependent)
  Float_t fProbHigh[2*AliPID::kSPECIES][fNptBins];     // probability above E/p threshold for NON electrons (charge dependent)


  ClassDef(AliEMCALPIDResponse, 1)
};

#endif // #ifdef AliEMCALPIDResponse_cxx

