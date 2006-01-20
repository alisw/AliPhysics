#ifndef ALIZDCRECONSTRUCTOR_H
#define ALIZDCRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for ZDC reconstruction                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliReconstructor.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliZDCCalibData.h"

class TF1;
class AliLoader;

class AliZDCReconstructor: public AliReconstructor {
public:
  AliZDCReconstructor();
  virtual ~AliZDCReconstructor();

  virtual void   Reconstruct(AliRunLoader* runLoader) const;
  virtual void   Reconstruct(AliRunLoader* runLoader,
  		  AliRawReader* rawReader) const;
  virtual void   Reconstruct(TTree* digitsTree, TTree* clustersTree) const 
  		  {AliReconstructor::Reconstruct(digitsTree,clustersTree);}
  virtual void   Reconstruct(AliRawReader* rawReader, TTree* clustersTree) const 
  		  {AliReconstructor::Reconstruct(rawReader,clustersTree);}
  virtual void   FillESD(AliRunLoader* runLoader, AliESD* esd) const;
  virtual void   FillESD(TTree* digitsTree, TTree* clustersTree, AliESD* esd) const 
  		  {AliReconstructor::FillESD(digitsTree,clustersTree,esd);}
  virtual void   FillESD(AliRawReader* rawReader, TTree* clustersTree, AliESD* esd) const 
  		  {AliReconstructor::FillESD(rawReader,clustersTree,esd);}
  virtual void   FillESD(AliRunLoader* runLoader, AliRawReader* rawReader, AliESD* esd) const 
  		  {AliReconstructor::FillESD(runLoader,rawReader,esd);}
  
  void   GetStorage(const char* uri) {fStorage = AliCDBManager::Instance()->GetStorage(uri);}
  AliCDBStorage   *SetStorage(const char* uri);
  AliZDCCalibData *GetCalibData(int runNumber) const; 

private:
  AliZDCReconstructor(const AliZDCReconstructor& reconstructor);
  AliZDCReconstructor& operator = (const AliZDCReconstructor& reconstructor);

  void   ReconstructEvent(AliLoader* loader, Float_t zncorr, Float_t zpcorr, Float_t zemcorr) const;

  TF1*   fZNCen;     //! Nspectator n true vs. EZN
  TF1*   fZNPer;     //! Nspectator n true vs. EZN
  TF1*   fZPCen;     //! Nspectator p true vs. EZP
  TF1*   fZPPer;     //! Nspectator p true vs. EZP
  TF1*   fZDCCen;    //! Nspectators true vs. EZDC
  TF1*   fZDCPer;    //! Nspectators true vs. EZDC
  TF1*   fbCen;      //! b vs. EZDC
  TF1*   fbPer;      //! b vs. EZDC
  TF1*   fZEMn;      //! Nspectators n from ZEM energy
  TF1*   fZEMp;      //! Nspectators p from ZEM energy
  TF1*   fZEMsp;     //! Nspectators from ZEM energy
  TF1*   fZEMb;      //! b from ZEM energy
  
  AliCDBStorage *fStorage; 	//! storage
  AliZDCCalibData *fCalibData; 	//! calibration data

  ClassDef(AliZDCReconstructor, 1)   // class for the ZDC reconstruction
};

#endif
