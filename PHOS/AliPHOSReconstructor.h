#ifndef ALIPHOSRECONSTRUCTOR_H
#define ALIPHOSRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  Wrapping class for reconstruction
//*--
//*-- Author: Yves Schutz (SUBATECH) 
// Reconstruction class. Redesigned from the old AliReconstructionner class and 
// derived from STEER/AliReconstructor. 
//_________________________________________________________________________

// --- ROOT system ---

#include "AliReconstructor.h" 
class AliPHOSDigitizer ;
class AliPHOSClusterizer ;
class AliPHOSTrackSegmentMaker ;
class AliPHOSPID ;
class AliPHOSSDigitizer ;
class AliESD ;

// --- Standard library ---

// --- AliRoot header files ---

class AliPHOSReconstructor : public AliReconstructor {

public:

  AliPHOSReconstructor() ; //ctor            
  AliPHOSReconstructor(const AliPHOSReconstructor & rec) : AliReconstructor(rec) {
    // cpy ctor: 
    // requested by the Coding Convention
    Fatal("cpy ctor", "not implemented") ;
  }
   
  virtual ~AliPHOSReconstructor() {} ;

  Bool_t                     Debug() const { return fDebug ; }
  virtual void               FillESD(AliRunLoader* runLoader, AliESD* esd) const ;
  virtual void               Reconstruct(AliRunLoader* runLoader) const ;

  AliPHOSReconstructor & operator = (const AliPHOSReconstructor & /*rvalue*/)  {
    // assignement operator requested by coding convention but not needed
    Fatal("operator =", "not implemented") ;
    return *this ; 
  }
  
private:
  
  Bool_t fDebug; //! verbosity controller

  ClassDef(AliPHOSReconstructor,2)  // Reconstruction algorithm class (Base Class)

}; 

#endif // ALIPHOSRECONSTRUCTOR_H
