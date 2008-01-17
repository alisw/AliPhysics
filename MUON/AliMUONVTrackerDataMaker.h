#ifndef ALIMUONVTRACKERDATAMAKER_H
#define ALIMUONVTRACKERDATAMAKER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup graphics
/// \class AliMUONVTrackerDataMaker
/// \brief Producer of some AliMUONVTrackerData
/// 
//  Author Laurent Aphecetche, Subatech

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class AliMUONVTrackerData;

class AliMUONVTrackerDataMaker : public TObject
{
public:
  AliMUONVTrackerDataMaker();
  virtual ~AliMUONVTrackerDataMaker();
  
  virtual Bool_t IsValid() const = 0;
  
  virtual AliMUONVTrackerData* Data() const = 0;

  virtual Bool_t IsRunnable() const = 0;

  virtual Bool_t IsRunning() const = 0;
  
  virtual void SetRunning(Bool_t flag) = 0;
  
  virtual Bool_t NextEvent() = 0;
  
  virtual void Rewind() = 0;
  
  /// Whether we're owner of our data
  virtual void SetOwner(Bool_t flag) = 0; 
  
  virtual void SetSource(const char* source) = 0;
  
  virtual TString Source() const = 0;
  
  ClassDef(AliMUONVTrackerDataMaker,1) // Producer of AliMUONVTrackerData
};

#endif
