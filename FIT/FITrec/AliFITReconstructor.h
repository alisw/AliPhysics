#ifndef ALIFITRECONSTRUCTOR_H
#define ALIFITRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/******************************************************************** 
 * header class FIT reconstruction 
 * Alla Maevskaya INR RAS alla@inr.ru      *
 * Alla.Maevskaya@cern.ch
 *******************************************************************/

#include "AliReconstructor.h"
#include "AliESDFIT.h"

class TClonesArray;

class AliFITReconstructor: public AliReconstructor {
 public:
  AliFITReconstructor();
  virtual ~AliFITReconstructor() {};

  virtual Bool_t HasDigitConversion() const { return kTRUE; }
  virtual void   ConvertDigits(AliRawReader* rawReader,
			       TTree* digitsTree) const;

  virtual void Init();
  virtual  void   Reconstruct(TTree* /*fdigits*/, TTree * /*frecpoints*/) const {};
  virtual  void   Reconstruct(AliRawReader*/*rawReader*/ , TTree* /*recTree*/) const {};
  
  virtual void     FillESD( AliRawReader*/*rawReader*/,  TTree* /*clustersTree*/, AliESDEvent* /*esd*/ ) const
  {}
  virtual void     FillESD( TTree* digitsTree,  TTree* /*clustersTree*/, AliESDEvent* esd ) const;
 
 // static const AliFITRecoParam* GetRecoParam()
   // { return dynamic_cast<const AliFITRecoParam*>(AliReconstructor::GetRecoParam(11)); } // getting RecoParam obj  TObjArray           fQTC;        // QTC vs #MIPs
 
 protected:

   AliESDEvent*             fESD;       // ESD object
   AliESDFIT*        fESDFIT;       // ESD output object  
   mutable TClonesArray * fDigits;
  
 private:
  AliFITReconstructor( const AliFITReconstructor&r ); //Not implemented
  AliFITReconstructor& operator=(const AliFITReconstructor&r); //Not implemented

  ClassDef(AliFITReconstructor, 1)   // class for the T0 reconstruction

};


#endif
