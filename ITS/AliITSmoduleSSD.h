#ifndef ALIITSMODLUESSD_H
#define ALIITSMODLUESSD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "TObject.h"
#include "TArrayS.h"
#include "TClonesArray.h"
#include "TRandom.h"
//#include "AliConst.h"
#include "TMath.h"

#include "AliITSdigitSSD.h"
#include "AliITS.h"
#include "AliITSmodule.h"



//____________________________________________________________________
//
//  Class AliITSmoduleSSD
//  describes one SSD module
//  The main function of modules is to simulate DIGITS from  
//  GEANT HITS and produce POINTS from DIGITS
//
//  This class keep track of Detector Control System parameters
//  of one SSD module: invalid (dead) strips, gain, capacitive
//  coupling.
//
//  In this release it make only simulation, and is not completly
//  tuned - up. Simulation will be improved after tests of SSD 
//  modules on SPS. Improved simulation and reconstruction will
//  appear in next releas.
//  
//  created by: A.Boucham, W.Peryt, S.Radomski, P.Skowronski
//  ver. 1.0    CERN, 16.09.1999  
// 
//___________________________________________________________________
//



class AliITSmoduleSSD: public AliITSmodule {

    
public:       
    
    //________________________________________________________________
    //
    // Constructors and deconstructor
    //________________________________________________________________
    //
    
    AliITSmoduleSSD();
    AliITSmoduleSSD(Int_t index);
    ~AliITSmoduleSSD();
  
    //________________________________________________________________
    //
    // Data process methods
    //________________________________________________________________
    //
    
    void AddDigit(Int_t, Int_t,  Bool_t);
    void HitToDigit();           // Process all hits in module
    void HitToDigit(Int_t);      // Proces one hit
    void DigitToPoint() {};      // Not impemented yet
    void HitToPoint() {};        // Not impemented y
   
    //________________________________________________________________
    //	
    //Invalid strips menagement methods
    //________________________________________________________________
    //
    
    // Parameters for invalid strips MonteCarlo
//    void SetInvalidParam(Float_t mean, Float_t sigma);
//    void GetInvalParam(Float_t &mean, Float_t &sigma); 
    
    // Methods for creating invalid strips
    void SetInvalidMC(Float_t mean, Float_t sigma);
    void SetInvalidMC();
    
    // Testing if strip is valid
    Bool_t  IsValidN(Int_t strip);      //True if strip work properly
    Bool_t  IsValidP(Int_t strip);
//    TArrayI GetInvalidP();              //Array of invalid strips
//    TArrayI GetInvalidN();
//    Int_t   GetNInvalidP();             //Number of invalid srtips
//    Int_t   GrtNInvalidN();
    
    // Creating invalid strips
    void      SetInvalidP(Int_t strip, Bool_t side);   //Set invalid if true 
//    void      SetInvalidN(Int_t strip, Bool_t side);          
       
  
protected:
         
    //________________________________________________________________
    //
    //Private methods for geometry
    //________________________________________________________________
    //
    
    Int_t   GetStripN(Float_t x, Float_t z);    // Nearest strip number P-side
    Int_t   GetStripP(Float_t x, Float_t z);    // Nearest strip number N-side
    
    Float_t Get2StripN(Float_t, Float_t);       // Ditance do the nearest strip P 
    Float_t Get2StripP(Float_t, Float_t);       // Ditance do the nearest strip N
     
    Bool_t  GetCrossing(Float_t&, Float_t&);      //x, y of strips crossing 
    
    //________________________________________________________________
    //
    //Private methods for simulation  
    //________________________________________________________________
    //
    
    void ApplyNoise();
    void ApplyCoupling();

    Float_t   F(Float_t x, Float_t s) {return (TMath::Erf(x*kPitch/s)+1) /2;}              
    
    // Proceding part should be in SSDgeo ----->
    
    //Technical parameters of detector
    static const Float_t   kStereo = 0.0175;  //Stereo Angle 17.5 mrad
    static const Float_t   kTan = 0.0175;  
    static const Int_t     kNStrips = 768;    //Number of strips on each side
    static const Float_t   kPitch = 0.095;    //Distance strip - strip (mm)
    static const Float_t   kX = 72.96;        //X size (mm)
    static const Float_t   kY = 0.3;          //Y size (mm)
    static const Float_t   kZ = 40;           //Thickness (mm)
    
    // <------------------------------
  
    //______________________________________________________________
    //  
    // Parameters for simulation
    //______________________________________________________________
      
    static const Float_t   kSigmaP = 0.003;     //Gaussian sigm
    static const Float_t   kSigmaN = 0.002;
    static const Int_t     kSteps  = 10;        //Number of steps 
    static const Int_t     kTresholdP = 1500;    
    static const Int_t     kTresholdN = 2500; 
   
    //________________________________________________________________
    //      
    // DCS Parameters
    //________________________________________________________________
    //
    
    Float_t   fSNRatioP;      //Signal - Noise Ratio P-side
    Float_t   fSNRatioN;      //Signal - Noise RatioNP-side
    
    Float_t   fGainP;         //Charge - ADC conversion parameter P
    Float_t   fGainN;         //Charge - ADC conversi parameter N
    
    Int_t     fNInvalidP;     //Number of invalid strips P
    TArrayS  *fInvalidP;      //Invalid strips P-side
    Int_t     fNInvalidN;     //Number of invalid strips N
    TArrayS  *fInvalidN;      //Invalid strips N-side 


    //________________________________________________________________
    //
    // Capacitive coupling parameters
    //________________________________________________________________
    //
    
    Float_t   fCouplingPR;
    Float_t   fCouplingPL;
    Float_t   fCouplingNR;
    Float_t   fCouplingNL;     

    //________________________________________________________________
    //
    // Parameters for invalid strips simulatation 
    //________________________________________________________________
    
    Float_t     fNInvalid;             //Meam number of invalid strips 
    Float_t     fISigma;               //RMS of invalid strips (Gaussian)

    //________________________________________________________________
    //
    // temp for simulation
    //________________________________________________________________
    //
    
    TArrayI *fN;         // for signal
    TArrayI *fP;        
    
    TArrayI *fNtrack1;   // for tracks, signal orgin N-side
    TArrayI *fNtrack2;
    TArrayI *fNtrack3;
       
    TArrayI *fPtrack1;   // for tracks, signal orgin P-side
    TArrayI *fPtrack2;
    TArrayI *fPtrack3;

public:
    ClassDef(AliITSmoduleSSD, 1)
};

#endif
