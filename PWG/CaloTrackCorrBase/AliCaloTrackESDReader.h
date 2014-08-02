#ifndef ALICALOTRACKESDREADER_H
#define ALICALOTRACKESDREADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
// Class for reading data (ESDs) in order to do prompt gamma 
//  or other particle identification and correlations
// 
//
//
//
//*-- Author: Gustavo Conesa (INFN-LNF)

class AliESDEvent;

#include "AliCaloTrackReader.h" 

class AliCaloTrackESDReader : public AliCaloTrackReader {
  
public:
  
                   AliCaloTrackESDReader() ; // ctor
  
  virtual         ~AliCaloTrackESDReader() ; // virtual dtor

  Bool_t           CheckForPrimaryVertex() const ;
  
  void             Init();
  
  Bool_t           SelectTrack(AliVTrack* track, Double_t* pTrack);
  
  AliESDtrackCuts* GetTrackCuts()                    const { return fESDtrackCuts     ; }
  void             SetTrackCuts(AliESDtrackCuts * cuts) ;
  
  AliESDtrackCuts* GetTrackComplementaryCuts()       const { return fESDtrackComplementaryCuts ; }
  void             SetTrackComplementaryCuts(AliESDtrackCuts * cuts)  ;

  void             SwitchOnConstrainTrackToVertex()        { fConstrainTrack = kTRUE  ; }
  void             SwitchOffConstrainTrackToVertex()       { fConstrainTrack = kFALSE ; }
  
  void             SetInputOutputMCEvent(AliVEvent* esd, AliAODEvent* aod, AliMCEvent* mc) ;
	 
private:
  
  Bool_t           fConstrainTrack;            // Constrain Track to vertex
  AliESDtrackCuts* fESDtrackCuts ;             // Track cut
  AliESDtrackCuts* fESDtrackComplementaryCuts; // Track cut, complementary cuts for hybrids
  
  ClassDef(AliCaloTrackESDReader,2)
  
} ;


#endif //ALICALOTRACKESDREADER_H



