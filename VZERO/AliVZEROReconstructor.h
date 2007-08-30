#ifndef ALIVZERORECONSTRUCTOR_H
#define ALIVZERORECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.*/
/* See cxx source for full Copyright notice                              */
/* $Id$  */

///////////////////////////////////////////////////////////////////////////
///                                                                      //
/// class for VZERO reconstruction                                       //
///                                                                      //
///////////////////////////////////////////////////////////////////////////

#include "AliReconstructor.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"

#include "AliVZEROCalibData.h"
#include "AliLog.h"
#include "AliESDVZERO.h"

class AliESDEvent;

class AliVZEROReconstructor: public AliReconstructor {
public:
  AliVZEROReconstructor();
  virtual ~AliVZEROReconstructor();
  virtual void   Init(AliRunLoader* /*runLoader*/);
  virtual void   Reconstruct(AliRunLoader* /*runLoader*/) const {
    AliError("Method not implemented"); return;};
  
  virtual void   Reconstruct(AliRawReader* /*rawReader*/, 
		             TTree* /*clustersTree*/) const {
    AliError("Method not implemented"); return;};
  virtual void   Reconstruct(AliRunLoader* /*runLoader*/, 
                             AliRawReader* /*rawReader*/) const {
    AliError("Method not implemented"); return;};
  virtual void   Reconstruct(TTree*, TTree*) const {return;};
  
  virtual void   FillESD(AliRunLoader* /*runLoader*/, AliESDEvent* /*esd*/) const {
    AliInfo("Method is not used"); return;};
  
  virtual void   FillESD(TTree* digitsTree, TTree* /*clustersTree*/, 
			 AliESDEvent* esd) const;

  virtual void   FillESD(AliRawReader* /*rawReader*/, TTree* /*clustersTree*/, 
			 AliESDEvent* /*esd*/) const {
    AliError("Method not implemented"); return;};
  
  virtual void   FillESD(AliRunLoader* /*runLoader*/, 
			 AliRawReader* /*rawReader*/, AliESDEvent* /*esd*/) const {
    AliInfo("Method is not used"); return;};
  
  virtual Bool_t HasDigitConversion() const { return kTRUE; }
  virtual void   ConvertDigits(AliRawReader* rawReader,
			       TTree* digitsTree) const;
  virtual Bool_t HasLocalReconstruction() const { return kTRUE; }

  AliCDBStorage     *SetStorage(const char* uri);
  AliVZEROCalibData *GetCalibData() const; 

protected:
  AliESDVZERO*        fESDVZERO;      // ESD output object  
  AliESDEvent*                  fESD;      // ESD object
  
private:
  AliVZEROReconstructor(const AliVZEROReconstructor& reconstructor);
  AliVZEROReconstructor& operator = (const AliVZEROReconstructor& reconstructor);
  
private:
  AliVZEROCalibData* fCalibData;      //! calibration data
 
  ClassDef(AliVZEROReconstructor, 0)  // class for the VZERO reconstruction
};

#endif
