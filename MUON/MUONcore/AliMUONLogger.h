#ifndef ALIMUONLOGGER_H
#define ALIMUONLOGGER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup core
/// \class AliMUONLogger
/// \brief A logger that keeps track of the number of times a message appeared
/// 
//  Author Laurent Aphecetche

#include <Riostream.h>

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

#ifndef ROOT_TString
#  include "TString.h"
#endif

class AliMUONStringIntMap;

using std::ofstream;

class AliMUONLogger : public TObject
{
public:
  AliMUONLogger(Int_t maxNumberOfEntries=-1, const char* name="AliMUONLogger");
  virtual ~AliMUONLogger();
  
  Int_t  Log(const char* message);
  
  void   Print(Option_t* opt="") const;
  
  void   Print(TString& key, ofstream& out) const;
  
  void   Clear(Option_t* /*option*/ ="");
  
  Bool_t Next(TString& msg, Int_t& occurance);
  
  void   ResetItr();
  
  Int_t NumberOfEntries() const;
  
  Long64_t Merge(TCollection* list);

  const char* GetName() const { return fName.Data(); }

  ULong_t Hash() const { return fName.Hash(); }

private:
  /// Not implemented
  AliMUONLogger(const AliMUONLogger& rhs); // not implemented
  /// Not implemented
  AliMUONLogger& operator=(const AliMUONLogger& rhs); // not implemented
  
private:
  
  Int_t fMaxNumberOfEntries; //!<! after this number, print and reset
  AliMUONStringIntMap* fLog; //!<! map from message to number of times the message was issued
  TString fName; //!<! object name
  
  ClassDef(AliMUONLogger,2) // A logger that keeps track of the number of times a message appeared
};

#endif
