#ifndef AliOADBFillingScheme_H
#define AliOADBFillingScheme_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
//     OADB container for filling scheme information (BX ids, name ...)
//     Author: Michele Floris, CERN
//-------------------------------------------------------------------------

#include <TNamed.h>
#include "TMap.h"
#include "TObjString.h"


class AliOADBFillingScheme : public TNamed {

 public :
  AliOADBFillingScheme();
  AliOADBFillingScheme(char* name);
  virtual ~AliOADBFillingScheme();
  AliOADBFillingScheme(const AliOADBFillingScheme& cont); 
  AliOADBFillingScheme& operator=(const AliOADBFillingScheme& cont);
  void Init();
  
  // Getters
  const char * GetBXIDs(const char * beamSide) const; 
  const char * GetFillingSchemeName() const { return fFSName; } 
  // Setters
  void SetBXIDs(const char * beamSide, const char * bxids) { fBXIds->Add(new TObjString(beamSide), new TObjString(bxids)); }
  void SetFillingSchemeName(const char * name) { fFSName = name; }
  // Browse
  virtual Bool_t	IsFolder() const { return kTRUE; }
  void Browse(TBrowser *b);
  // Print
  virtual void	Print(Option_t* option = "") const;

 private :
  
  TString fFSName               ; // Name of the filling scheme 
  TMap * fBXIds              ; // Map from the beam side bunch crossing number. Beam side is "B", "A", "C", "E".

  ClassDef(AliOADBFillingScheme, 1);
};

#endif
