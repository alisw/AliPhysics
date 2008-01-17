#ifndef ALIMUONPAINTERENV_H
#define ALIMUONPAINTERENV_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup graphics
/// \class AliMUONPainterEnv
/// \brief Resource file handling
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class TEnv;

class AliMUONPainterEnv : public TObject
{
public:
  AliMUONPainterEnv(const char* resourceFile=".mchviewrc");
  virtual ~AliMUONPainterEnv();
  
  const char* String(const char* resourceName, const char* defaultValue="");
  
  Int_t Integer(const char* resourceName, Int_t defaultValue=0);
  
  Double_t Double(const char* resourceName, Double_t defaultValue=0.0);
  
  void Save();
  
  void Set(const char* resourceName, Int_t value);

  void Set(const char* resourceName, const char* value);

  void Set(const char* resourceName, Double_t value);

private:
  /// Not implemented
  AliMUONPainterEnv(const AliMUONPainterEnv& rhs);
  /// Not implemented
  AliMUONPainterEnv& operator=(const AliMUONPainterEnv& rhs);
  
  TEnv* fEnv; ///< the worker class
  
  ClassDef(AliMUONPainterEnv,1) // Painter display resource file
};

#endif
