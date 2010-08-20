/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

// AliStarEventReader, a class for reading the star ntuples
// origin: Jim Thomas,        jhthomas@lbl.gov
//         Mikolaj Krzewicki, mikolaj.krzewicki@cern.ch

#ifndef ALISTAREVENTREADER_H
#define ALISTAREVENTREADER_H

#include <TObject.h>

class TList   ;
class TNtuple ;
class AliStarTrack;
class AliStarEvent;

class AliStarEventReader : public TObject {

 private:
  TList *fFileList ;  //file list 

  AliStarEventReader& operator=( const AliStarEventReader& event ); //not implemented
  AliStarEventReader(const AliStarEventReader& event); //not implemented

 public:
  TNtuple *fEventHeader;  //star event header
  TNtuple *fTracks;       //track data
  AliStarEvent* fEvent;   //encapsulated star event

  AliStarEventReader();
  AliStarEventReader( const char* inputFileDirectory );      
  virtual ~AliStarEventReader();

  virtual Bool_t  AcceptEvent( AliStarEvent* event ); 
  virtual Bool_t  AcceptTrack( AliStarTrack* track ); 
  virtual Bool_t  GetNextEvent();
  virtual Int_t   Centrality( Int_t referenceMultiplicity );   
  virtual Bool_t  MakeFileList( const char* inputFileDirectory );
  virtual Bool_t  MakeFileListFromDir( const char* inputFileDirectory );
  virtual Bool_t  MakeFileListFromFile( const char* inputFileName );
  virtual Int_t   ParticleID( AliStarTrack* track ) ;

  virtual void    PrintEventHeader( ) ;
  virtual void    PrintTrack( Int_t j ) ;
  const AliStarEvent* GetEvent() const {return fEvent;}

  ClassDef(AliStarEventReader,1)         // Base class
};
#endif

