//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTESDEVENT_H
#define ALIHLTESDEVENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTESDEvent.h
/// @author Matthias Richter
/// @date   2010-10-29
/// @brief  A streamlined container class for AliESDEvent.
/// @note   

#include "AliESDEvent.h"

/**
 * @class AliHLTESDEvent
 * @brief A streamlined container class for AliESDEvent.
 *
 * The class inherits from AliESDEvent and can be used like that, it only
 * implements customized streamers to treat some of the objects in the list
 * specifically.
 */
class AliHLTESDEvent : public AliESDEvent {
 public:
  /// standard constructor
  AliHLTESDEvent();
  /// copy constructor
  AliHLTESDEvent(const AliHLTESDEvent& src);
  /// assignement operator
  AliHLTESDEvent& operator=(const AliHLTESDEvent& src);
  /// destructor
  virtual ~AliHLTESDEvent();

  AliHLTESDEvent& operator=(const AliESDEvent& esd) {
    if (this!=&esd) AliESDEvent::operator=(esd);
    return *this;
  }

  /// overloaded from TObject, print info
  virtual void        Print(const char* options) const;

  /// overloaded from TObject, more crude data dump
  virtual void        Dump() const;

  /// overloaded from TObject, clear object
  virtual void        Clear(Option_t * option="");

  /// overloaded from TObject, clone object
  virtual TObject    *Clone(const char *newname="") const;

  /// overloaded from TObject, copy object
  virtual void        Copy(TObject &object) const;

  /// overloaded from TObject, execute custum function
  /// implemented: LoadTemplate
  virtual void        Execute(const char *method,  const char *params, Int_t *error=0);

  /// load a template from OCDB or create the default template
  int LoadTemplate(const char* cdbPath=NULL);

private:
  /// the template instance
  AliESDEvent* fTemplateEsd; //!

  ClassDef(AliHLTESDEvent, 1); // AliESDEvent instance optimized for HLT
};

#endif
