#ifndef ALIITSDCSSSD_H
#define ALIITSDCSSSD_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///////////////////////////////////////////////////////////////////////////
//                                                                       //
//  Class AliITSdcsSSD                                                   //
//  describes Detector Control System parameters for one SSD module.     //
//                                                                       //
//  This class stores parametrers such as gain, threshold                //
//  capacitive coupling.                                                 //
//                                                                       //
//  Class takes care of invalid strip menagement during                  //
//  simulation and runtime                                               //
//                                                                       //
//                                                                       //
//  created at: Warsaw University of Technology                          //
//  ver. 1.0    WARSAW, 23.12.1999                                       //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include <TObject.h>

class TArrayS;
class TRandom;
class AliITSsegmentation;
class AliITSCalibration;

class AliITSdcsSSD: public TObject {

 public:
    AliITSdcsSSD(); // Default constructor
    // Standard constructor
    AliITSdcsSSD(AliITSsegmentation *s, AliITSCalibration *r); 
    virtual ~AliITSdcsSSD(); // Destructor
    AliITSdcsSSD(const AliITSdcsSSD &source); // copy constructor
    AliITSdcsSSD& operator=(const AliITSdcsSSD &source); // assignment operator
    //________________________________________________________________
    //	
    // Invalid strips management methods
    //________________________________________________________________
    // Parameters for invalid strips MonteCarlo
    void SetInvalidParam(Float_t mean, Float_t sigma);
    void GetInvalidParam(Float_t &mean, Float_t &sigma) const;
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
    void    SetInvalidP(Int_t,Bool_t){}//Set invalid if true
    void    SetInvalidN(Int_t,Bool_t){}//Set invalid if true
    Float_t  GetCouplingPR() const {// couplings
      return fCouplingPR;}
    Float_t  GetCouplingPL() const {// couplings
      return fCouplingPL;}
    Float_t  GetCouplingNR() const {// couplings
      return fCouplingNR;}
    Float_t  GetCouplingNL() const {// couplings
      return fCouplingNL;}
    
 protected:
    //_____________________________________________
    //
    // Parameters for invalid strips simulatation 
    //_____________________________________________    
    Float_t  fCouplingPR;  // couplings
    Float_t  fCouplingPL;  // couplings
    Float_t  fCouplingNR;  // couplings
    Float_t  fCouplingNL;  // couplings   

    Float_t   fNstrips;    //Number of strips
    Float_t   fNInvalid;   //Mean number of invalid strips (for simulation) 
    Float_t   fISigma;     //RMS of invalid strips (Gaussian)

    TArrayS  *fInvalidP;   //Array of invalid strips on P-side
    TArrayS  *fInvalidN;   //Array of invalid strips on N-side

    ClassDef(AliITSdcsSSD, 1)     // ITS SSD DCS specific class

};
#endif
