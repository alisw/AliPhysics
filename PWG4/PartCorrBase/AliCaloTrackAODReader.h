#ifndef ALICALOTRACKAODREADER_H
#define ALICALOTRACKAODREADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id: $ */

//_________________________________________________________________________
// Class for reading data (AODs) in order to do prompt gamma or other particle
// identification and correlations.
// This part is commented: Mixing analysis can be done, input AOD with events
// is opened in the AliCaloTrackReader::Init()
//
//
// -- Author: Gustavo Conesa (INFN-LNF)

#include "AliAnalysisTaskSE.h"

#include "AliAODEvent.h"
#include "AliCaloTrackReader.h" 

class AliCaloTrackAODReader : public AliCaloTrackReader {
	
public: 
	
  AliCaloTrackAODReader() ; // ctor
  virtual ~AliCaloTrackAODReader() {;} //virtual dtor
     
  AliCentrality* GetCentrality() const ;  
  void SetInputOutputMCEvent(AliVEvent* esd, AliAODEvent* aod, AliMCEvent* mc) ; 
  
  AliVEvent* GetOriginalInputEvent() const { return fOrgInputEvent; }
  
  TString GetFiredTriggerClasses() {return ((AliAODEvent*)GetInputEvent())->GetFiredTriggerClasses();}
  
private:
  
  AliVEvent *fOrgInputEvent; //! Original input event, not from filtering
  
  AliCaloTrackAODReader(const AliCaloTrackAODReader & ) ; // cpy ctor
  AliCaloTrackAODReader & operator = (const AliCaloTrackAODReader & g) ;//cpy assignment
  
  ClassDef(AliCaloTrackAODReader,6)
} ;

#endif //ALICALOTRACKAODREADER_H



