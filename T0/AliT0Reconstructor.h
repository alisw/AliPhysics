#ifndef ALIT0RECONSTRUCTOR_H
#define ALIT0RECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliReconstructor.h"
#include "AliT0digit.h"
#include "AliT0.h"

class AliT0Reconstructor: public AliReconstructor {
 public:
  AliT0Reconstructor();
  virtual ~AliT0Reconstructor() {};
  AliT0Reconstructor( const AliT0Reconstructor& );
  AliT0Reconstructor& operator=(const AliT0Reconstructor&); 


  virtual  void   Reconstruct(TTree* fdigits, TTree * frecpoints) const;
  virtual  void   Reconstruct(AliRawReader*rawReader , TTree* recTree) const;
  
  virtual void         FillESD( AliRawReader*,  TTree*clustersTree, AliESDEvent*esd ) const
  {FillESD((TTree*)NULL,clustersTree,esd);}
  virtual void         FillESD( TTree*,  TTree*, AliESDEvent* ) const;

  virtual Bool_t       HasDigitConversion() const {return kFALSE;}
 public:
 
  Float_t              fZposition; // vertex position
  
 protected:
  AliT0Parameters     *fParam;           //pointer to T0 parameters class     
  TObjArray           fAmpLEDrec;        // amp LED-CFD 
  Float_t             fTime0vertex[24];  // time position if Zvertex=0
  Float_t             fdZ_A;   // Zideal - Zreal side A 
  Float_t             fdZ_C; // Zideal - Zreal side C


  ClassDef(AliT0Reconstructor, 0)   // class for the T0 reconstruction

};

typedef AliT0Reconstructor AliSTARTReconstructor; // for backward compatibility

#endif
