#ifndef ALIMUONCHECK_H
#define ALIMUONCHECK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup base
/// \class AliMUONCheck
/// \brief A helper class to dump data from AliRoot-generated root files.
/// 
//  Author Laurent Aphecetche

#ifndef ROOT_TObject
#  include "TObject.h"
#endif
#ifndef ROOT_TString
#  include "TString.h"
#endif

class AliMUONData;
class AliRunLoader;

class AliMUONCheck : public TObject
{
public:
  AliMUONCheck(const char* galiceFile, Int_t firstEvent=0, Int_t lastEvent=-1);
  virtual ~AliMUONCheck();
  
  /// Return true if contains valid data
  Bool_t IsValid() const { return (fData!=0); }
  
  void DumpDigits(Option_t* opt="") const;
  
private:
  AliMUONCheck(const AliMUONCheck& rhs);
  AliMUONCheck& operator=(const AliMUONCheck& rhs);
  
private:
  TString fFileName;   //!< File (galice.root) to read from
  Int_t   fFirstEvent; //!< First event to consider
  Int_t   fLastEvent;  //!< Last event to consider
  AliRunLoader* fRunLoader; //!< AliRunLoader pointer
  AliMUONData*  fData; //!< AliMUONData pointer (to access containers)

  ClassDef(AliMUONCheck,0) // Dumper of MUON related data
};

#endif
