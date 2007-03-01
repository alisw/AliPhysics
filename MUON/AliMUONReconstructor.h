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
class AliMUONDigitMaker;
class AliMUONTriggerCrateStore;
class AliMUONGeometryTransformer;
class AliTracker;
class AliMUONClusterReconstructor;
class AliMUONSegmentation;

class TTask;
class TClonesArray;

class AliMUONReconstructor: public AliReconstructor 
{
  public:
    AliMUONReconstructor();
    virtual ~AliMUONReconstructor();

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

    AliTracker*          CreateTracker(AliRunLoader* runLoader) const;
 
    enum {kNone, kOriginal, kKalman, kCombi};

private:

    TTask* GetCalibrationTask() const;
    AliMUONClusterReconstructor* CreateClusterReconstructor() const;
    
    AliMUONReconstructor(const AliMUONReconstructor& right);
    AliMUONReconstructor&  operator = (const AliMUONReconstructor& right);

private:
    AliMUONDigitMaker* fDigitMaker; //!< pointer to the digit maker class

    mutable AliMUONCalibrationData* fCalibrationData; //!< pointer to calibration data
    
    AliMUONTriggerCrateStore* fCrateManager;     //!< Crate array

    TClonesArray* fTriggerCircuit;   //!< trigger circuit
 
    AliMUONGeometryTransformer* fTransformer; //!< pointer to transformation
    AliMUONSegmentation*        fSegmentation; //!< pointer to segmentation

    AliMUONData* fMUONData;                    //!< pointer to container

  ClassDef(AliMUONReconstructor, 0)   // class for the MUON reconstruction
};

#endif
