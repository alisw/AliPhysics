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
// Author: Michael Weber (m.weber@cern.ch)                              //
//////////////////////////////////////////////////////////////////////////

#include "AliPID.h"
#include <TVectorD.h>

class TF1;

class AliEMCALPIDResponse: public TObject 
{
public : 
    AliEMCALPIDResponse();    //ctor
    AliEMCALPIDResponse( const AliEMCALPIDResponse& other);                //copy ructor
    AliEMCALPIDResponse &operator=( const AliEMCALPIDResponse& other);     //assignment operator

    virtual ~AliEMCALPIDResponse();     //dtor
  

    // Getters
    Double_t  GetNumberOfSigmas( Float_t pt,  Float_t eop, AliPID::EParticleType n,  Int_t charge) const;
    Double_t  GetExpectedNorm  ( Float_t pt, AliPID::EParticleType n,  Int_t charge) const;
  
    //Setters
    void   SetPIDParams(const TObjArray * params) { fkPIDParams = params; }
    

    // EMCAL probability -> should go to another place?
    Double_t ComputeEMCALProbability( Float_t pt, Float_t eop, Int_t charge, Double_t *pEMCAL) const;

protected:
  
private:

  TF1 *fNorm;                            // Gauss function for normalizing NON electron probabilities 

  const TObjArray *fkPIDParams;               // PID Params

  const TVectorD* GetParams(Int_t nParticle, Float_t fPt) const; 

  ClassDef(AliEMCALPIDResponse, 1)
};

#endif // #ifdef AliEMCALPIDResponse_cxx

