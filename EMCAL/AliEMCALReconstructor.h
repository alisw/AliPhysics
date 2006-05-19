#ifndef ALIEMCALRECONSTRUCTOR_H
#define ALIEMCALRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  Wrapping class for reconstruction
//*--
//*-- Author: Yves Schutz (SUBATECH) 
//*--         Dmitri Peressounko (SUBATECH & Kurchatov Institute)


// --- ROOT system ---

#include "AliReconstructor.h" 
class AliEMCALDigitizer ;
class AliEMCALClusterizer ;
class AliEMCALSDigitizer ;
class AliESD ;
class AliRawReader ;

// --- Standard library ---

// --- AliRoot header files ---

class AliEMCALReconstructor : public AliReconstructor {

public:

  AliEMCALReconstructor() ; //ctor            
  AliEMCALReconstructor(const AliEMCALReconstructor & rec) : AliReconstructor(rec) {
    // cpy ctor: 
    // requested by the Coding Convention
    Fatal("cpy ctor", "not implemented") ;
  }
   
  virtual ~AliEMCALReconstructor() ; //dtor

  Bool_t       Debug() const { return fDebug ; }

  using AliReconstructor::FillESD;
  virtual void FillESD(AliRunLoader* runLoader, AliESD* esd) const ;

  using AliReconstructor::Reconstruct;
  virtual void Reconstruct(AliRunLoader* runLoader) const ;
  virtual void Reconstruct(AliRunLoader* runLoader, AliRawReader* rawreader) const ;
  
  
  AliEMCALReconstructor & operator = (const AliEMCALReconstructor & /*rvalue*/)  {
    // assignement operator requested by coding convention but not needed
    Fatal("operator =", "not implemented") ;
    return *this ; 
  }
  

private:
  
  Bool_t fDebug; //! verbosity controller
 
  ClassDef(AliEMCALReconstructor,1)  // Reconstruction algorithm class (Base Class)

}; 

#endif // ALIEMCALRECONSTRUCTOR_H
