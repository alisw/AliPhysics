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

#include "AliCaloTrackReader.h" 

class AliCaloTrackAODReader : public AliCaloTrackReader {
	
public: 
	
  AliCaloTrackAODReader() ; // ctor
  virtual ~AliCaloTrackAODReader() {;} //virtual dtor
  
  void SetInputOutputMCEvent(AliVEvent* esd, AliAODEvent* aod, AliMCEvent* mc) ; 
  
  ClassDef(AliCaloTrackAODReader,5)
} ;

#endif //ALICALOTRACKAODREADER_H



