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

 public:

  AliStarEventReader();
  AliStarEventReader( const char* inputFileDirectory );      
  virtual ~AliStarEventReader();

  virtual Bool_t  GetNextEvent();
  virtual Bool_t  MakeFileList( const char* inputFileDirectory );
  virtual Bool_t  MakeFileListFromDir( const char* inputFileDirectory );
  virtual Bool_t  MakeFileListFromFile( const char* inputFileName );

  const AliStarEvent* GetEvent() const {return fEvent;}

 private:
  TList *fFileList ;  //file list 
  AliStarEvent* fEvent;   //encapsulated star event

  AliStarEventReader& operator=( const AliStarEventReader& event ); //not implemented
  AliStarEventReader(const AliStarEventReader& event); //not implemented

  ClassDef(AliStarEventReader,1)         // Base class
};
#endif

