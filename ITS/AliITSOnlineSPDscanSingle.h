#ifndef ALIITSONLINESPDSCANSINGLE_H
#define ALIITSONLINESPDSCANSINGLE_H  
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                 //
// Interface class to the containers of an online scan    //
// with only one step.                                    //
////////////////////////////////////////////////////////////

#include "AliITSOnlineSPDscan.h"

class AliITSOnlineSPDscanSingle :  public AliITSOnlineSPDscan {

 public:
  AliITSOnlineSPDscanSingle() {}
  AliITSOnlineSPDscanSingle(const Char_t *fileName, Bool_t readFromGridFile=kFALSE);
  virtual ~AliITSOnlineSPDscanSingle();
  
  // SET METHODS ***********************************
  void     SetHits(UInt_t hs, UInt_t chipi, UInt_t coli, UInt_t rowi, UInt_t val) ;
  void     IncrementTriggers() ;
  void     SetTriggers(UInt_t val) ;
  void     IncrementHits(UInt_t hs, UInt_t chipi, UInt_t coli, UInt_t rowi) ;
  void     SetHitEvents(UInt_t hs, UInt_t chipi, UInt_t val) ;
  void     SetHitEventsTot(UInt_t hs, UInt_t val) ;
  void     IncrementHitEvents(UInt_t hs, UInt_t chipi) ;
  void     IncrementHitEventsTot(UInt_t hs) ;
  // GET METHODS ***********************************
  UInt_t   GetTriggers() ;
  UInt_t   GetHits(UInt_t hs, UInt_t chipi, UInt_t coli, UInt_t rowi) ;
  Float_t  GetHitsEfficiency(UInt_t hs, UInt_t chipi, UInt_t coli, UInt_t rowi) ;
  Float_t  GetHitsEfficiencyError(UInt_t hs, UInt_t chipi, UInt_t coli, UInt_t rowi) ;
  UInt_t   GetHitEvents(UInt_t hs, UInt_t chipi) ;
  UInt_t   GetHitEventsTot(UInt_t hs) ;
  Float_t  GetHitEventsEfficiency(UInt_t hs, UInt_t chipi) ;
  Float_t  GetHitEventsTotEfficiency(UInt_t hs) ;
  Float_t  GetHitEventsEfficiencyError(UInt_t hs, UInt_t chipi) ;
  Float_t  GetHitEventsTotEfficiencyError(UInt_t hs) ;
  Float_t  GetAverageMultiplicity(UInt_t hs, UInt_t chipi) ;
  Float_t  GetAverageMultiplicityTot(UInt_t hs) ;


};

#endif
