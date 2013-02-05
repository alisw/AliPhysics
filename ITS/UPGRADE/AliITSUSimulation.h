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
#include "AliMathBase.h"

class AliITSCalibration;
class AliITSUSimuParam;
class AliITSUModule;
class TRandom;
class TSegCollection;
class AliParamList;

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
  virtual void InitSimulationModule(AliITSUModule* mod, Int_t ev, AliITSsegmentation* seg, AliParamList* resp);
  //
  // Hits -> SDigits
  virtual void SDigitiseModule() = 0;
  virtual void FinishSDigitiseModule() = 0;
  virtual Bool_t AddSDigitsToModule( TSeqCollection *pItemArray, Int_t mask );
  //
  // Hits -> Digits
  virtual void DigitiseModule() = 0;
  virtual void CreateFastRecPoints(AliITSUModule *,Int_t,TRandom *,TClonesArray* /*recp*/) {}
  //
  // readout phase (strobe, timing etc) generation
  virtual void GenerateStrobePhase() {}
  //
  AliITSCalibration*  GetCalibDead()                   const {return fCalibDead;}
  AliITSCalibration*  GetCalibNoisy()                  const {return fCalibNoisy;}
  AliITSsegmentation* GetSegmentation()                const {return fSeg;}
  AliITSUSimuParam*   GetSimuParam()                   const {return fSimuParam;}
  AliITSUSensMap*     GetMap()                         const {return fSensMap;}
  AliITSUModule*      GetModule()                      const {return fModule;}
  AliParamList*       GetResponseParam()               const {return fResponseParam;}
  Int_t               GetEvent()                       const {return fEvent;}
  Bool_t              GetDebug(Int_t level=1)          const {return fDebug>=level;}
  
  //
  void SetCalibDead(AliITSCalibration *calib)              {fCalibDead = calib;}
  void SetCalibNoisy(AliITSCalibration *calib)             {fCalibNoisy = calib;}
  void SetSegmentation(AliITSsegmentation *seg)            {fSeg = seg; if (seg&&fSensMap) fSensMap->SetDimensions(seg->Npz(),seg->Npx());}
  void SetSimuParam(AliITSUSimuParam *sp)                  {fSimuParam = sp;}
  virtual void SetResponseParam(AliParamList* resp)        {fResponseParam = resp;}
  void SetMap(AliITSUSensMap *p)                           {fSensMap = p;}
  void SetModule(AliITSUModule* mod)                       {fModule=mod;} 
  void SetEvent(Int_t evnt)                                {fEvent=evnt;} 
  void SetDebug(Int_t level=5)                             {fDebug=level;}
  void SetNoDebug()                                        {fDebug=0;}
  void ClearMap()                                          {fSensMap->Clear();}
  //
  static  Int_t GenOrderedSample(UInt_t nmax,UInt_t ngen,TArrayI &vals,TArrayI &ind);
  //
  static  Double_t GausInt1D(Double_t sig,Double_t a,Double_t b);
  static  Double_t GausInt2D(Double_t sig0,Double_t a0,Double_t b0,
			     Double_t sig1,Double_t a1,Double_t b1);
  //
 protected:
  AliITSsegmentation  *fSeg;            //! segmentation
  AliITSCalibration   *fCalibDead;      //! dead channels
  AliITSCalibration   *fCalibNoisy;     //! noisy channels
  AliITSUSensMap      *fSensMap;        //! sensor map for hits manipulations
  AliITSUSimuParam    *fSimuParam;      //! simulation parameters
  AliParamList        *fResponseParam;  //! response parameterization data
  AliITSUModule       *fModule;         //! module being processed
  Int_t                fEvent;          //! event number being processed
  Int_t                fDebug;          //!  debug flag

  ClassDef(AliITSUSimulation,1)       // Simulation base class 
    
};

//_____________________________________________________________________________
inline Double_t AliITSUSimulation::GausInt1D(Double_t sig,Double_t a,Double_t b)
{
  // calculate gaussian integral from a to b (with respecto to mean)
  const Double_t kRoot2 = 1.414213562; // Sqrt(2).
  double sp = 1.0/(sig*kRoot2);
  return 0.5*TMath::Abs( AliMathBase::ErfcFast(sp*a) - AliMathBase::ErfcFast(sp*b) );
}

//_____________________________________________________________________________
inline Double_t AliITSUSimulation::GausInt2D(Double_t sig0,Double_t a0,Double_t b0,
					     Double_t sig1,Double_t a1,Double_t b1)
{
  // calculate gaussian 2D integral from x0 to x1 (with respect to mean)
  const Double_t kRoot2 = 1.414213562; // Sqrt(2).
  double sp0 = 1.0/(sig0*kRoot2);
  double sp1 = 1.0/(sig1*kRoot2);
  return 0.25*TMath::Abs( (AliMathBase::ErfcFast(sp0*a0) - AliMathBase::ErfcFast(sp0*b0)) *
			  (AliMathBase::ErfcFast(sp1*a1) - AliMathBase::ErfcFast(sp1*b1)));
}

#endif
