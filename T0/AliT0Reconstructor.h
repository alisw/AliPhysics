#ifndef ALIT0RECONSTRUCTOR_H
#define ALIT0RECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/*  Alla Maevskaya INR RAS alla@inr.ru */
/* $Id$ */


#include "AliReconstructor.h"
#include "AliT0Parameters.h"
#include "AliT0.h"

class AliT0Reconstructor: public AliReconstructor {
 public:
  AliT0Reconstructor();
  virtual ~AliT0Reconstructor() {};
  AliT0Reconstructor( const AliT0Reconstructor&r );
  AliT0Reconstructor& operator=(const AliT0Reconstructor&r); 


  virtual  void   Reconstruct(TTree* fdigits, TTree * frecpoints) const;
  virtual  void   Reconstruct(AliRawReader*rawReader , TTree* recTree) const;
  
  virtual void     FillESD( AliRawReader*/*rawReader*/,  TTree*clustersTree, AliESDEvent*esd ) const
  {FillESD((TTree*)NULL,clustersTree,esd);}
  virtual void     FillESD( TTree* digitsTree,  TTree*clustersTree, AliESDEvent*esd ) const;

  virtual Bool_t   HasDigitConversion() const {return kFALSE;}
   
 protected:
  Float_t             fdZonA;             // Zideal - Zreal side A 
  Float_t             fdZonC;             // Zideal - Zreal side C
  Float_t             fZposition;        // vertex position
  Float_t             fTime0vertex[24];  // time position if Zvertex=0
  AliT0Parameters     *fParam;           //pointer to T0 parameters class     
  TObjArray           fAmpLEDrec;        // amp LED-CFD 

  ClassDef(AliT0Reconstructor, 1)   // class for the T0 reconstruction

};

typedef AliT0Reconstructor AliSTARTReconstructor; // for backward compatibility

#endif
