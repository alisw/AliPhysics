#ifndef ROOT_TMevSimPartTypeParams
#define ROOT_TMevSimPartTypeParams
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TMevSimPartTypeParams                                                //
//                                                                      //
// This class is a helper class for TMevSim interface. It stores        //
// (and writes out on demand) the properties of a single particle       //
// type.                                                                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#ifndef ROOT_TObject
#include <TObject.h>
#endif

#include "MevSimCommon.h"

class TMevSimPartTypeParams : public TObject {

 protected:
  
  Int_t       fGPid;                     
  Int_t       fMultMean;                
  Int_t       fMultVarianceControl;     
  Float_t     fTempMean, fTempStDev;
  Float_t     fSigmaMean, fSigmaStDev;
  Float_t     fExpVelMean, fExpVelStDev;
  Float_t     fVnMean[NFLOWTERMS][4];
  Float_t     fVnStDev[NFLOWTERMS][4];
 
  // GPid : particle ID ala geant3
  // MultMean, MultVarianceControl: mean multiplicy and variance control
  //           MultVarianceControl 0; for no variance in multiplicity
  //           MultVarianceControl 1; to allow Poisson distribution for particle multiplicities 
  // TempMean, TempStDev: Temperature parameter (in GeV)
  //           mean and standard deviation (Gaussian distribution assumed) 
  // SigmaMean, SigmaStDev: Rapidity distribution width (sigma)  
  //           mean and standard deviation (Gaussian distribution assumed) 
  // ExpVelMean, ExpVelStDev: Expansion velocity ala Scott Pratt (in units of c)  
  //           mean and standard deviation (Gaussian distribution assumed) 
  // VnMean VnStDev: Anisotropic flow parameters for Fourier components NFLOWTERMS=1,6
  //                 mean and standard deviation
  //                 include all 6 sets of parameters, set them to 0 if not used                                

 public:
  
   // Constructors and destructors
   
   TMevSimPartTypeParams();
   TMevSimPartTypeParams(Int_t agpid, Int_t amultmean, Int_t amultvc, 
			 Float_t atempmean, Float_t atempstdev, Float_t asigmamean,
			 Float_t asigmastdev, Float_t aexpvelmean, Float_t aexpvelstdev);
   virtual ~TMevSimPartTypeParams();
   
   // Copy and assignment operators;
   
   TMevSimPartTypeParams (const TMevSimPartTypeParams& pars);                    // copy constructor
   virtual TMevSimPartTypeParams& operator=(const TMevSimPartTypeParams& pars);  // assignment operator 
   
   // Parameters of the particle type
   
   virtual void        SetGPid(Int_t gpid);
   virtual Int_t       GetGPid() const;
   
   virtual void        SetMultMean(Int_t multmean);
   virtual Int_t       GetMultMean() const;
   
   virtual void        SetMultVarianceControl(Int_t multvc);
   virtual Int_t       GetMultVarianceControl() const;
   
   virtual void        SetTempParams(Float_t tempmean, Float_t tempstdev);
   virtual Float_t     GetTempMean() const;
   virtual Float_t     GetTempStDev() const;
   
   virtual void        SetSigmaPrams(Float_t sigmamean, Float_t sigmastdev);
   virtual Float_t     GetSigmaMean() const;
   virtual Float_t     GetSigmaStDev() const;
   
   virtual void        SetExpVelParams(Float_t expvelmean, Float_t expvelstdev);
   virtual Float_t     GetExpVelMean() const;
   virtual Float_t     GetExpVelStDev() const;
   
   virtual void        SetVnMeanComponent(Int_t nComponent, Float_t mean1, Float_t mean2, 
					  Float_t mean3, Float_t mean4);
   virtual void        SetVnStDevComponent(Int_t nComponent, Float_t stdev1, Float_t stdev2, 
					   Float_t stdev3, Float_t stdev4);
					   
   virtual void        SetVnMeanComponent (Int_t nComponent, Float_t mean[4]); 
   virtual void        SetVnStDevComponent(Int_t nComponent, Float_t stdev[4]);			   
					   
   virtual Float_t     GetVnMeanComponent(Int_t nComponent, Int_t nMean) const;
   virtual Float_t     GetVnStDevComponent(Int_t nComponent, Int_t nStDev) const;
   

   ClassDef(TMevSimPartTypeParams,1)            //Parameters of the type of particle for MevSim event generator

     
};
#endif












