#ifndef ALIITSMFTSIMULATION_H
#define ALIITSMFTSIMULATION_H
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
#include "AliITSMFTSensMap.h"
#include "AliITSMFTSegmentationPix.h"
#include "AliMathBase.h"
#include <TArrayS.h>

class AliITSMFTCalibration;
class AliITSMFTSimuParam;
class AliITSMFTChip;
class TArrayI;
class TRandom;
class TSegCollection;
class AliITSMFTParamList;

// This is the base class for ITS detector signal simulations. Data members
// include are a pointer to the detectors specific response and segmentation
// classes. See the detector specific implementations for the propper code.
// The detector-specific implementations are responsible for processing detector type, 
// which still may have different segmentations!
// Note: the detector specific objects (segmentation, calib, etc.) are not owned
// by the object but set externaly

class AliITSMFTSimulation : public TObject 
{
 public:
  enum {kMaxROCycleAccept=126};             // flag for read-out cycle to discard
  //
  AliITSMFTSimulation();
  AliITSMFTSimulation(AliITSMFTSimuParam* sim, AliITSMFTSensMap* map);
  virtual ~AliITSMFTSimulation() {} 
  AliITSMFTSimulation(const AliITSMFTSimulation &source);
  AliITSMFTSimulation& operator=(const AliITSMFTSimulation &source);
  virtual void Init() = 0;
  //
  void UpdateMapSignal(UInt_t col,UInt_t row, Int_t trk,Int_t ht,Double_t signal, Int_t roCycle=0);
  void UpdateMapNoise(UInt_t col,UInt_t row, Double_t noise, Int_t roCycle=0);
  virtual void InitSimulationChip(AliITSMFTChip* mod, Int_t ev, AliITSMFTSegmentationPix* seg, AliITSMFTParamList* resp);
  //
  // Hits -> SDigits
  virtual void SDigitiseChip(TClonesArray *) = 0;
  virtual void FinishSDigitiseChip(TObjArray *) = 0;
  virtual Bool_t AddSDigitsToChip( TSeqCollection *pItemArray, Int_t mask );
  //
  // Hits -> Digits
  virtual void DigitiseChip(TObjArray *) = 0;
  virtual void CreateFastRecPoints(AliITSMFTChip *,Int_t,TRandom *,TClonesArray* /*recp*/) {}
  //
  // readout phase (strobe, timing etc) generation
  virtual Double_t GenerateReadOutCycleOffset();
  //
  AliITSMFTCalibration*  GetCalibDead()                   const {return fCalibDead;}
  AliITSMFTCalibration*  GetCalibNoisy()                  const {return fCalibNoisy;}
  AliITSMFTSegmentationPix* GetSegmentation()                const {return fSeg;}
  AliITSMFTSimuParam*   GetSimuParam()                   const {return fSimuParam;}
  AliITSMFTSensMap*     GetMap()                         const {return fSensMap;}
  AliITSMFTChip*      GetChip()                      const {return fChip;}
  AliITSMFTParamList*   GetResponseParam()               const {return fResponseParam;}
  Int_t               GetEvent()                       const {return fEvent;}
  Bool_t              GetDebug(Int_t level=1)          const {return fDebug>=level;}
  
  //
  void SetCalibDead(AliITSMFTCalibration *calib)              {fCalibDead = calib;}
  void SetCalibNoisy(AliITSMFTCalibration *calib)             {fCalibNoisy = calib;}
  void SetSegmentation(AliITSMFTSegmentationPix *seg)            {fSeg = seg; if (seg&&fSensMap) fSensMap->SetDimensions(seg->Npz(),seg->Npx(),kMaxROCycleAccept);}
  void SetSimuParam(AliITSMFTSimuParam *sp)                  {fSimuParam = sp;}
  virtual void SetResponseParam(AliITSMFTParamList* resp)    {fResponseParam = resp;}
  void SetMap(AliITSMFTSensMap *p)                           {fSensMap = p;}
  void SetChip(AliITSMFTChip* mod)                       {fChip=mod;} 
  void SetEvent(Int_t evnt)                                {fEvent=evnt;} 
  void SetDebug(Int_t level=5)                             {fDebug=level;}
  void SetNoDebug()                                        {fDebug=0;}
  void ClearMap()                                          {fSensMap->Clear();}
  //
  Double_t GetReadOutCycleOffset()              const      {return fReadOutCycleOffset;}
  void     SetReadOutCycleOffset(Double_t v=0)             {fReadOutCycleOffset = v;}
  //
  Double_t GetReadOutCycleLength()              const      {return fReadOutCycleLength;}
  void     SetReadOutCycleLength(Double_t v=25e-6)         {fReadOutCycleLength = v;}
  //
  static  Int_t GenOrderedSample(UInt_t nmax,UInt_t ngen,TArrayI &vals,TArrayI &ind);
  //
  static  Double_t GausInt1D(Double_t sig,Double_t a,Double_t b);
  static  Double_t GausInt2D(Double_t sig0,Double_t a0,Double_t b0,
			     Double_t sig1,Double_t a1,Double_t b1);
  //
 protected:
  AliITSMFTSegmentationPix  *fSeg;            //! segmentation
  AliITSMFTCalibration   *fCalibDead;      //! dead channels
  AliITSMFTCalibration   *fCalibNoisy;     //! noisy channels
  AliITSMFTSensMap      *fSensMap;        //! sensor map for hits manipulations
  AliITSMFTSimuParam    *fSimuParam;      //! simulation parameters
  AliITSMFTParamList    *fResponseParam;  //! response parameterization data
  AliITSMFTChip       *fChip;         //! chip being processed
  Float_t              fReadOutCycleOffset; //! The phase of the RO with respect to the trigger
  Float_t              fReadOutCycleLength; //! readout cycle lenght in s

  Int_t                fEvent;          //! event number being processed
  Int_t                fDebug;          //!  debug flag
  Bool_t               fCyclesID[2*kMaxROCycleAccept+1]; //! status of RO cycles

  ClassDef(AliITSMFTSimulation,3)       // Simulation base class 
    
};

//_____________________________________________________________________________
inline Double_t AliITSMFTSimulation::GausInt1D(Double_t sig,Double_t a,Double_t b)
{
  // calculate gaussian integral from a to b (with respecto to mean)
  const Double_t kRoot2 = 1.414213562; // Sqrt(2).
  double sp = 1.0/(sig*kRoot2);
  return 0.5*TMath::Abs( AliMathBase::ErfcFast(sp*a) - AliMathBase::ErfcFast(sp*b) );
}

//_____________________________________________________________________________
inline Double_t AliITSMFTSimulation::GausInt2D(Double_t sig0,Double_t a0,Double_t b0,
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
