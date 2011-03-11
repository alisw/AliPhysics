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

#include "AliESDEvent.h"
#include "AliCaloTrackReader.h" 

class AliCaloTrackESDReader : public AliCaloTrackReader {
  
  public: 
  
  AliCaloTrackESDReader() ; // ctor
  virtual ~AliCaloTrackESDReader() {;} //virtual dtor

  AliCentrality* GetCentrality() const {
    AliESDEvent* event = dynamic_cast<AliESDEvent*> (fInputEvent);	 
    if(event) return event->GetCentrality() ; 
    else      return 0x0                    ; }  
  
  void SetInputOutputMCEvent(AliVEvent* esd, AliAODEvent* aod, AliMCEvent* mc) ; 
	
  TString GetFiredTriggerClasses() {return ((AliESDEvent*)GetInputEvent())->GetFiredTriggerClasses();}
    
  ClassDef(AliCaloTrackESDReader,1)
} ;


#endif //ALICALOTRACKESDREADER_H



