#ifndef ALICALOTRACKESDREADER_H
#define ALICALOTRACKESDREADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id:  $ */

//_________________________________________________________________________
// Class for reading data (ESDs) in order to do prompt gamma 
//  or other particle identification and correlations
// 
//
//
//
//*-- Author: Gustavo Conesa (INFN-LNF)

// --- ROOT system --- 

// --- AliRoot system ---
#include "AliCaloTrackReader.h" 
#include "AliESDEvent.h"

class AliCaloTrackESDReader : public AliCaloTrackReader {
  
  public: 
  
  AliCaloTrackESDReader() ; // ctor
  //AliCaloTrackESDReader(const AliCaloTrackESDReader & g) ; // cpy ctor
  //AliCaloTrackESDReader & operator = (const AliCaloTrackESDReader & g) ;//cpy assignment
  virtual ~AliCaloTrackESDReader() {;} //virtual dtor

  Double_t GetBField() const;
  void SetInputOutputMCEvent(AliVEvent* esd, AliAODEvent* aod, AliMCEvent* mc) ; 
	
  TString GetFiredTriggerClasses() {return ((AliESDEvent*)GetInputEvent())->GetFiredTriggerClasses();}

  ClassDef(AliCaloTrackESDReader,1)
} ;


#endif //ALICALOTRACKESDREADER_H



