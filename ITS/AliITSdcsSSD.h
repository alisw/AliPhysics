#ifndef ALIITSDCSSSD_H
#define ALIITSDCSSSD_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>
#include <TArrayS.h>
#include <TRandom.h>

//____________________________________________________________________
//
//  Class AliITSdcsSSD 
//  describes Detector Control System parameters for one SSD module.
//   
//  This class stores parametrers such as gain, threshold
//  capacitive coupling.
//  
//  Class takes care of invalid strip menagement during 
//  simulation and runtime
//
//
//  created at: Warsaw University of Technology
//  ver. 1.0    WARSAW, 23.12.1999  
// 
//___________________________________________________________________
//


class AliITSsegmentation;
class AliITSresponse;


class AliITSdcsSSD: public TObject {

public:       
    
    //________________________________________________________________
    //
    // Constructors and deconstructor
    //________________________________________________________________
    //
    
    AliITSdcsSSD(AliITSsegmentation *s, AliITSresponse *r);
    ~AliITSdcsSSD();
    AliITSdcsSSD(const AliITSdcsSSD &source); // copy constructor
    AliITSdcsSSD& operator=(const AliITSdcsSSD &source); // assignment operator
    
    //________________________________________________________________
    //	
    // Invalid strips management methods
    //________________________________________________________________
    //
    
    // Parameters for invalid strips MonteCarlo
    
    void SetInvalidParam(Float_t mean, Float_t sigma);
    void GetInvalidParam(Float_t &mean, Float_t &sigma); 
    
    
    // Methods for creating invalid strips
    
    void SetInvalidMC(Float_t mean, Float_t sigma);
    void SetInvalidMC();
    
    
    // Testing if strip is valid
    
    Bool_t  IsValidN(Int_t strip);      //True if strip works properly
    Bool_t  IsValidP(Int_t strip);      //True if strip works properly
    
    
    // Access to invalid strips
    
    TArrayS *GetInvalidP();             //Array of invalid P strips
    TArrayS *GetInvalidN();             //Array of invalid N strips
    Int_t    GetNInvalidP();            //Number of invalid P strips
    Int_t    GetNInvalidN();            //Number of invalid N strips
    
    
    // Creating invalid strips
    
    void    SetInvalidP(Int_t strip, Bool_t side){
      //Set invalid if true   
    }
    
    void    SetInvalidN(Int_t strip, Bool_t side){
      //Set invalid if true
    }

    Float_t  GetCouplingPR() {
      // couplings
      return fCouplingPR;
    }

    Float_t  GetCouplingPL() {
      // couplings
      return fCouplingPL;
    }

    Float_t  GetCouplingNR() {
      // couplings
      return fCouplingNR;
    }

    Float_t  GetCouplingNL() {
      // couplings
      return fCouplingNL;
    }
    
 protected:   
    
    //_____________________________________________
    //
    // Parameters for invalid strips simulatation 
    //_____________________________________________
    
    Float_t  fCouplingPR;  // couplings
    Float_t  fCouplingPL;  // couplings
    Float_t  fCouplingNR;  // couplings
    Float_t  fCouplingNL;  // couplings   

    Float_t   fNstrips ;          //Number of strips
    Float_t   fNInvalid;          //Mean number of invalid strips (for simulation) 
    Float_t   fISigma;            //RMS of invalid strips (Gaussian)

    TArrayS  *fInvalidP;          //Array of invalid strips on P-side
    TArrayS  *fInvalidN;          //Array of invalid strips on N-side
   
    TRandom  *fRandom;            //!Random numbers generator

    ClassDef(AliITSdcsSSD, 1)     // ITS SSD DCS specific class

};


#endif
