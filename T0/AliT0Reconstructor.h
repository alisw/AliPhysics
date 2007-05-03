#ifndef ALIT0RECONSTRUCTOR_H
#define ALIT0RECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliReconstructor.h"
#include "AliT0digit.h"
#include "AliT0.h"
class AliRunLoader;

class AliT0Reconstructor: public AliReconstructor {
 public:
  AliT0Reconstructor();
  virtual ~AliT0Reconstructor() {};
  AliT0Reconstructor( const AliT0Reconstructor& );
  AliT0Reconstructor& operator=(const AliT0Reconstructor&); 


   virtual  void   Reconstruct(TTree* fdigits, TTree * frecpoints) const;
  //   virtual  void   Reconstruct(TTree* , TTree * ) const {};
//  virtual  void   Reconstruct(TTree* , TTree * ) const {};
 
  virtual  void   Reconstruct(AliRunLoader* , AliRawReader*   ) const {};
  virtual  void   Reconstruct(AliRawReader* ) const {};
  virtual  void   Reconstruct(AliRawReader*rawReader , TTree* recTree) const;
  virtual  void   Reconstruct(AliRunLoader* ) const  {};
  
  virtual void         FillESD(AliRunLoader* runLoader, AliESD* esd) const;
  virtual void         FillESD(AliRunLoader* , AliRawReader*, AliESD* ) const  {};
  virtual void         FillESD(  AliRawReader*,  TTree*, AliESD* ) const  {};
  virtual void         FillESD( TTree*,  TTree*, AliESD* ) const  {};
  virtual Bool_t       HasLocalReconstruction() const {return kTRUE;}
  virtual Bool_t       HasDigitConversion() const {return kFALSE;}
 public:
 
  Float_t fZposition; // vertex position

 protected:

  ClassDef(AliT0Reconstructor, 0)   // class for the T0 reconstruction

};

typedef AliT0Reconstructor AliSTARTReconstructor; // for backward compatibility

#endif
