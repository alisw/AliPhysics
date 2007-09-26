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
#include "AliEMCALTracker.h" 
class AliEMCALDigitizer ;
class AliEMCALClusterizer ;
class AliEMCALSDigitizer ;
class AliEMCALRecParam;
class AliESDEvent ;
class AliRawReader ;

// --- Standard library ---

// --- AliRoot header files ---

class AliEMCALReconstructor : public AliReconstructor {

public:

  AliEMCALReconstructor() ; //ctor            
  AliEMCALReconstructor(const AliEMCALReconstructor & rec);
   
  virtual ~AliEMCALReconstructor() ; //dtor

  Bool_t       Debug() const { return fDebug ; }

  using AliReconstructor::FillESD;
  virtual void FillESD(TTree* digitsTree, TTree* clustersTree, 
		       AliESDEvent* esd) const;
  AliTracker*  CreateTracker () const 
  {return new AliEMCALTracker;} 
  using AliReconstructor::Reconstruct;
  virtual void Reconstruct(TTree* digitsTree, TTree* clustersTree) const;

  virtual Bool_t             HasDigitConversion() const {return kTRUE;};
  virtual void               ConvertDigits(AliRawReader* rawReader, TTree* digitsTree) const;
  
  
  AliEMCALReconstructor & operator = (const AliEMCALReconstructor & /*rvalue*/)  {
    // assignement operator requested by coding convention but not needed
    Fatal("operator =", "not implemented") ;
    return *this ; 
  }
  
  void SetRecParam(AliEMCALRecParam * recParam){ fgkRecParam = recParam;}

  static const AliEMCALRecParam* GetRecParam(){ return fgkRecParam;}

private:
  
  Bool_t fDebug; //! verbosity controller
  static AliEMCALRecParam*   fgkRecParam; // reconstruction parameters for EMCAL

  ClassDef(AliEMCALReconstructor,2)  // Reconstruction algorithm class (Base Class)

}; 

#endif // ALIEMCALRECONSTRUCTOR_H
