#ifndef ALICALOTRACKAODREADER_H
#define ALICALOTRACKAODREADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
// Class for reading data (AODs) in order to do prompt gamma or other particle
// identification and correlations.
// This part is commented: Mixing analysis can be done, input AOD with events
// is opened in the AliCaloTrackReader::Init()
//
//
// -- Author: Gustavo Conesa (INFN-LNF)

class AliAODEvent;

#include "AliCaloTrackReader.h" 

class AliCaloTrackAODReader : public AliCaloTrackReader {
	
public: 
	
  AliCaloTrackAODReader() ;            // ctor
  
  virtual ~AliCaloTrackAODReader() {;} // virtual dtor
     
  AliCentrality* GetCentrality() const ;  
  
  void SetInputOutputMCEvent(AliVEvent* esd, AliAODEvent* aod, AliMCEvent* mc) ; 
  
  AliVEvent* GetOriginalInputEvent() const { return fOrgInputEvent; }
    
private:
  
  AliVEvent *fOrgInputEvent; //! Original input event, not from filtering
  
  AliCaloTrackAODReader(              const AliCaloTrackAODReader & r) ; // cpy ctor
  AliCaloTrackAODReader & operator = (const AliCaloTrackAODReader & r) ; // cpy assignment
  
  ClassDef(AliCaloTrackAODReader,6)
  
} ;

#endif //ALICALOTRACKAODREADER_H



