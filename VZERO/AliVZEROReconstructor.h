#ifndef ALIVZERORECONSTRUCTOR_H
#define ALIVZERORECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
///                                                                          //
/// class for VZERO reconstruction                                           //
///                                                                          //
///////////////////////////////////////////////////////////////////////////////

#include "AliReconstructor.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliVZEROCalibData.h"
#include "AliCDBEntry.h"

class AliLoader;

class AliVZEROReconstructor: public AliReconstructor {
public:
  AliVZEROReconstructor();
  virtual ~AliVZEROReconstructor();

  AliCDBStorage     *SetStorage(const char* uri);
  AliVZEROCalibData *GetCalibData() const; 

private:

  AliVZEROReconstructor(const AliVZEROReconstructor& reconstructor);
  AliVZEROReconstructor& operator = (const AliVZEROReconstructor& reconstructor);
  
  AliVZEROCalibData *fCalibData;       //! calibration data
 
  ClassDef(AliVZEROReconstructor, 0)   // class for the VZERO reconstruction
};

#endif
