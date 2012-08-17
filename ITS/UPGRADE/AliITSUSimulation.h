#ifndef ALIITSUSIMULATION_H
#define ALIITSUSIMULATION_H
/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include <TObject.h>
#include "AliITSUSensMap.h"
#include "AliITSsegmentation.h"

class AliITSCalibration;
class AliITSUSimuParam;
class AliITSUModule;
class TRandom;
class TSegCollection;

// This is the base class for ITS detector signal simulations. Data members
// include are a pointer to the detectors specific response and segmentation
// classes. See the detector specific implementations for the propper code.
// The detector-specific implementations are responsible for processing detector type, 
// which still may have different segmentations!
// Note: the detector specific objects (segmentation, calib, etc.) are not owned
// by the object but set externaly

class AliITSUSimulation : public TObject 
{
 public:
  AliITSUSimulation();
  AliITSUSimulation(AliITSUSimuParam* sim, AliITSUSensMap* map);
  virtual ~AliITSUSimulation() {} 
  AliITSUSimulation(const AliITSUSimulation &source);
  AliITSUSimulation& operator=(const AliITSUSimulation &source);
  virtual void Init() = 0;
  //
  void UpdateMapSignal(UInt_t dim0,UInt_t dim1, Int_t trk,Int_t ht,Double_t signal);
  void UpdateMapNoise(UInt_t dim0,UInt_t dim1, Double_t noise);
  virtual void InitSimulationModule(Int_t, Int_t, AliITSsegmentation*);
  //
  // Hits -> SDigits
  virtual void SDigitiseModule(AliITSUModule* mod, Int_t mask, Int_t event, AliITSsegmentation* seg) = 0;
  virtual void FinishSDigitiseModule() = 0;
  virtual Bool_t AddSDigitsToModule( TSeqCollection *pItemArray, Int_t mask );
  //
  // Hits -> Digits
  virtual void DigitiseModule(AliITSUModule* mod, Int_t mask, Int_t event, AliITSsegmentation* seg) = 0;
  virtual void CreateFastRecPoints(AliITSUModule *,Int_t,TRandom *,TClonesArray* /*recp*/) {}
  //
  AliITSCalibration*  GetCalibDead()                   const {return fCalibDead;}
  AliITSCalibration*  GetCalibNoisy()                  const {return fCalibNoisy;}
  AliITSsegmentation* GetSegmentation()                const {return fSeg;}
  AliITSUSimuParam* GetSimuParam()                   const {return fSimuParam;}
  AliITSUSensMap*   GetMap()                         const {return fSensMap;}
  Int_t               GetModule()                      const {return fModule;}
  Int_t               GetEvent()                       const {return fEvent;}
  Bool_t              GetDebug(Int_t level=1)          const {return fDebug>=level;}

  //
  void SetCalibDead(AliITSCalibration *calib)              {fCalibDead = calib;}
  void SetCalibNoisy(AliITSCalibration *calib)             {fCalibNoisy = calib;}
  void SetSegmentation(AliITSsegmentation *seg)            {fSeg = seg; if (seg&&fSensMap) fSensMap->SetDimensions(seg->Npz(),seg->Npx());}
  void SetSimuParam(AliITSUSimuParam *sp)                {fSimuParam = sp;}
  void SetMap(AliITSUSensMap *p)                         {fSensMap = p;}
  void SetModule(Int_t mod)                                {fModule=mod;} 
  void SetEvent(Int_t evnt)                                {fEvent=evnt;} 
  void SetDebug(Int_t level=5)                             {fDebug=level;}
  void SetNoDebug()                                        {fDebug=0;}
  void ClearMap()                                          {fSensMap->Clear();}
  //
  static  Int_t GenOrderedSample(UInt_t nmax,UInt_t ngen,TArrayI &vals,TArrayI &ind);
  //
 protected:
  AliITSsegmentation  *fSeg;            //! segmentation
  AliITSCalibration   *fCalibDead;      //! dead channels
  AliITSCalibration   *fCalibNoisy;     //! noisy channels
  AliITSUSensMap    *fSensMap;        //! sensor map for hits manipulations
  AliITSUSimuParam  *fSimuParam;      //! simulation parameters
  Int_t                fModule;         //! module number being processed
  Int_t                fEvent;          //! event number being processed
  Int_t                fDebug;          //!  debug flag

  ClassDef(AliITSUSimulation,1)       // Simulation base class 
    
};

#endif
