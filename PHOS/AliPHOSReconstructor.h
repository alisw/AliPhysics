#ifndef ALIPHOSRECONSTRUCTOR_H
#define ALIPHOSRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.10  2007/07/11 13:43:30  hristov
 * New class AliESDEvent, backward compatibility with the old AliESD (Christian)
 *
 * Revision 1.9  2006/11/15 16:05:03  kharlov
 * New FillESD() for raw data is added
 *
 * Revision 1.8  2005/05/28 14:19:04  schutz
 * Compilation warnings fixed by T.P.
 *
 */

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
class AliESDEvent ;
class AliRawReader; 
class AliPHOSRecoParam;

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
  ~AliPHOSReconstructor() ; //dtor            
  static void                SetDebug()   { fgDebug = kTRUE ; }
  static void                ResetDebug() { fgDebug = kFALSE ; }
  static Bool_t              Debug() { return fgDebug ; }
  AliTracker *CreateTracker(AliRunLoader* runLoader) const;
  using AliReconstructor::FillESD;
  virtual void               FillESD(AliRunLoader* runLoader, AliESDEvent* esd) const ;
  virtual void FillESD(AliRunLoader* runLoader,AliRawReader* rawReader,AliESDEvent* esd) const;
  using AliReconstructor::Reconstruct;
  virtual void               Reconstruct(AliRunLoader* runLoader) const ;
  virtual void               Reconstruct(AliRunLoader* runLoader, AliRawReader * rawreader) const ;

  AliPHOSReconstructor & operator = (const AliPHOSReconstructor & /*rvalue*/)  {
    // assignement operator requested by coding convention but not needed
    Fatal("operator =", "not implemented") ;
    return *this ; 
  }
  
  void SetRecoParamEmc(AliPHOSRecoParam * param){ fgkRecoParamEmc = param;}
  void SetRecoParamCpv(AliPHOSRecoParam * param){ fgkRecoParamCpv = param;}

  static const AliPHOSRecoParam* GetRecoParamEmc(){ return fgkRecoParamEmc;}
  static const AliPHOSRecoParam* GetRecoParamCpv(){ return fgkRecoParamCpv;}

private:
  
  static Bool_t fgDebug ; //! verbosity controller
  static AliPHOSRecoParam*   fgkRecoParamEmc; // reconstruction parameters for EMC
  static AliPHOSRecoParam*   fgkRecoParamCpv; // reconstruction parameters for EMC

  ClassDef(AliPHOSReconstructor,2)  // Reconstruction algorithm class (Base Class)

}; 

#endif // ALIPHOSRECONSTRUCTOR_H
