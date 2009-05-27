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
class AliESDVZEROfriend;

class AliVZEROReconstructor: public AliReconstructor {
public:
  AliVZEROReconstructor();
  virtual ~AliVZEROReconstructor();
  virtual void   Init();
  
  virtual void   Reconstruct(AliRawReader* /*rawReader*/, 
		             TTree* /*clustersTree*/) const {
                 AliError("Method not implemented"); return;};
  virtual void   Reconstruct(TTree*, TTree*) const {return;};
  
  virtual void   FillESD(TTree* digitsTree, TTree* /*clustersTree*/, 
			 AliESDEvent* esd) const;

  virtual void   FillESD(AliRawReader* /*rawReader*/, TTree* /*clustersTree*/, 
			 AliESDEvent* /*esd*/) const {
                 AliError("Method not implemented"); return;};
  
  virtual Bool_t HasDigitConversion() const { return kTRUE; }
  virtual void   ConvertDigits(AliRawReader* rawReader,
			       TTree* digitsTree) const;
                 		 
  AliCDBStorage     *SetStorage(const char* uri);
  void GetCollisionMode();
  
  AliVZEROCalibData *GetCalibData() const; 

protected:
  AliESDVZERO*        fESDVZERO;       // ESD output object  
  AliESDEvent*             fESD;       // ESD object
  AliESDVZEROfriend*  fESDVZEROfriend; // ESD friend object (complete raw data)
  
private:
  AliVZEROReconstructor(const AliVZEROReconstructor& reconstructor);
  AliVZEROReconstructor& operator = (const AliVZEROReconstructor& reconstructor);
  
  AliVZEROCalibData* fCalibData;      //! calibration data

  Int_t              fCollisionMode;  // =0->p-p, =1->A-A
  Float_t            fBeamEnergy;     // beam energy
  
  ClassDef(AliVZEROReconstructor, 0)  // class for the VZERO reconstruction
};

#endif
