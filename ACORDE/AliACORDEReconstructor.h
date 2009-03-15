#ifndef ALIACORDERECONSTRUCTOR_H
#define ALIACORDERECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.*/
/* See cxx source for full Copyright notice                              */
/* $Id: AliACORDEReconstructor.h 20956 2007-09-26 14:22:18Z cvetan $  */

///////////////////////////////////////////////////////////////////////////
///                                                                      //
/// class for ACORDE reconstruction                                       //
///                                                                      //
///////////////////////////////////////////////////////////////////////////

#include "AliReconstructor.h"
#include "AliLog.h"
#include "AliACORDERecoParam.h"

class AliACORDECalibData;
class AliESDACORDE;
class AliESDEvent;

class AliACORDEReconstructor: public AliReconstructor {
public:
  AliACORDEReconstructor();
  virtual ~AliACORDEReconstructor();
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
  virtual void ConvertDigits(AliRawReader* rawReader, TTree* digitsTree) const;

  AliACORDECalibData *GetCalibData() const; 

  AliACORDERecoParam *GetRecoParam() const;


protected:

  AliESDACORDE*        fESDACORDE;      // ACORDE ESD object  
  AliACORDERecoParam* fAcordeRecoParam; // Pointer to the ACORDE's RecoParam

private:
  AliACORDEReconstructor(const AliACORDEReconstructor& reconstructor);
  AliACORDEReconstructor& operator = (const AliACORDEReconstructor& reconstructor);
  
  AliACORDECalibData* fCalibData;      //! calibration data
 
  ClassDef(AliACORDEReconstructor, 1)  // class for the ACORDE reconstruction
};

#endif
