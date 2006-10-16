#ifndef ALIMUONRECONSTRUCTOR_H
#define ALIMUONRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

/// \ingroup rec
/// \class AliMUONReconstructor
/// \brief Class for the MUON reconstruction

#include "AliReconstructor.h"

class AliMUONCalibrationData;
class AliMUONData;
class TTask;
class AliMUONDigitMaker;
class AliMUONTriggerCrateStore;
class TClonesArray;
class AliMpSegFactory;
class AliMUONGeometryTransformer;

class AliMUONReconstructor: public AliReconstructor 
{
  public:
    AliMUONReconstructor();
    virtual ~AliMUONReconstructor();

    virtual void         Init(AliRunLoader* runLoader);

    virtual void         Reconstruct(TTree* /*digitsTree*/, 
				     TTree* /*clustersTree*/) const {return;}
    virtual void         Reconstruct(AliRawReader* /*rawReader*/, 
				     TTree* /*clustersTree*/) const {return;}
    virtual void         Reconstruct(AliRunLoader* runLoader) const;
    virtual void         Reconstruct(AliRunLoader* runLoader, 
                                   AliRawReader* rawReader) const;

    virtual void         FillESD(TTree* /*digitsTree*/, TTree* /*clustersTree*/, 
				 AliESD* /*esd*/) const {return;}
    virtual void         FillESD(AliRawReader* /*rawReader*/, TTree* /*clustersTree*/, 
				 AliESD* /*esd*/) const {return;}
    virtual void         FillESD(AliRunLoader* runLoader, AliESD* esd) const;
    virtual void         FillESD(AliRunLoader* runLoader, 
				 AliRawReader* /*rawReader*/, AliESD* esd) const;
     
    enum {kNone, kOriginal, kKalman, kCombi};

private:

    TTask* GetCalibrationTask(AliMUONData* data) const;
    AliMUONReconstructor(const AliMUONReconstructor& right);
    AliMUONReconstructor&  operator = (const AliMUONReconstructor& right);

private:
    AliRunLoader* fRunLoader;       //!< pointer to runloader
    AliMUONDigitMaker* fDigitMaker; //!< pointer to the digit maker class

    mutable AliMUONCalibrationData* fCalibrationData; //!< pointer to calibration data
    
    AliMUONTriggerCrateStore* fCrateManager;     //!< Crate array

    TClonesArray* fTriggerCircuit;   //!< trigger circuit
 
    AliMpSegFactory* fSegFactory;        //!< Mapping segmentation factory

    AliMUONGeometryTransformer* fTransformer; //!< pointer to transformation

  ClassDef(AliMUONReconstructor, 0)   // class for the MUON reconstruction
};

#endif
