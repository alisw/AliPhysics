#ifndef AliMFT_H
#define AliMFT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//====================================================================================================================================================
//
//      Geometry of the Muon Forward Tracker based on TGeo
//
//      Contact author: antonio.uras@cern.ch
//
//====================================================================================================================================================

#include "TFile.h"
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoMatrix.h"
#include "TVirtualMC.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TGeoGlobalMagField.h"
#include "AliRun.h"
#include "AliLoader.h"
#include "AliDetector.h"
#include "AliMC.h"
#include "AliMagF.h"
#include "AliMFTHit.h"
#include "AliMFTDigit.h"
#include "AliMFTCluster.h"
#include "AliTrackReference.h"
#include "AliMFTSegmentation.h"
#include "AliMFTDigitizer.h"
#include "AliMFTPlane.h"
#include "TString.h"
#include "TObjArray.h"
#include "AliMFTConstants.h"

//====================================================================================================================================================

class AliMFT : public AliDetector  {
  
public:
  
  AliMFT();
  AliMFT(const Char_t *name, const Char_t *title);
  AliMFT(const Char_t *name, const Char_t *title, Char_t *nameGeomFile);
  
  virtual ~AliMFT();
  
  Int_t IsVersion() const { return fVersion; }
  
  // ------- framework part -----------------------------------------------------------------------------------
  void CreateMaterials();                       // from AliModule invoked from AliMC
  void CreateGeometry();                        // from AliModule invoked from AliMC
  void AddAlignableVolumes();
  void StepManager();                           // from AliModule invoked from AliMC::Stepping()
  void Hits2SDigits();
  void Hits2SDigitsLocal(TClonesArray *hits, const TObjArray *pSDig, Int_t track);
  void MakeBranch(Option_t *option="");
  void SetTreeAddress();
  
  // ------- create containers -----------------------------------------------------------------------------------
  void CreateHits();
  void CreateSDigits();
  void CreateDigits();
  void CreateRecPoints();
  
  TObjArray*    GetSDigitsList()            const { return fSDigitsPerPlane; }     // get sdigits list for all planes
  TClonesArray* GetSDigitsList(Int_t plane) const { return fSDigitsPerPlane ? (TClonesArray*) fSDigitsPerPlane->At(plane):0; }
  
  TObjArray*    GetDigitsList()            const{return fDigitsPerPlane;}          // get digits list for all layers
  TClonesArray* GetDigitsList(Int_t plane) const{return fDigitsPerPlane ? (TClonesArray*) fDigitsPerPlane->At(plane):0; }
  
  TObjArray*    GetRecPointsList()            const{return fRecPointsPerPlane;}    // get cluster list for all layers
  TClonesArray* GetRecPointsList(Int_t plane) const{return fRecPointsPerPlane ? (TClonesArray*) fRecPointsPerPlane->At(plane):0; }
  
  void ResetSDigits()   { if(fSDigitsPerPlane)   for(int iPlane=0; iPlane<AliMFTConstants::kNDisks; iPlane++) ((TClonesArray*) fSDigitsPerPlane  ->At(iPlane))->Clear(); }   // reset sdigits list
  void ResetDigits()    { if(fDigitsPerPlane)    for(int iPlane=0; iPlane<AliMFTConstants::kNDisks; iPlane++) ((TClonesArray*) fDigitsPerPlane   ->At(iPlane))->Clear(); }   // reset digits list
  void ResetRecPoints() { if(fRecPointsPerPlane) for(int iPlane=0; iPlane<AliMFTConstants::kNDisks; iPlane++) ((TClonesArray*) fRecPointsPerPlane->At(iPlane))->Clear(); }   // reset recPoints list
  
  AliDigitizer* CreateDigitizer(AliDigitizationInput *digInp) const { return new AliMFTDigitizer(digInp); }
  
  AliMFTSegmentation* GetSegmentation() const { return fSegmentation; }
  
  enum EMedia{kZero,kAir, kVacuum, kSi, kReadout, kSupport, kCarbon, kBe, kAlu, kWater, kSiO2, kInox, kKapton, kEpoxy, kCarbonFiber, kCarbonEpoxy, kRohacell, kPolyimide, kPEEK, kFR4, kCu};  // media IDs used in CreateMaterials
  
  // Geometry/segmentation creation part
  void AddAlignableVolumes() const { /* not needed */ return; }
  void SetGeometry();
  
  Int_t GetNPlanes() const { return fNPlanes; }
  
  void SetChargeDispersion(Double_t chargeDispersion) { fChargeDispersion = chargeDispersion; }
  Double_t GetChargeDispersion() { return fChargeDispersion; }
  void SetNStepForChargeDispersion(Int_t nStepForChargeDispersion) { fNStepForChargeDispersion = nStepForChargeDispersion; }
  Int_t GetNStepForChargeDispersion() { return fNStepForChargeDispersion; }
  Double_t GetSingleStepForChargeDispersion() { return fSingleStepForChargeDispersion; }
  
  void SetDensitySupportOverSi(Double_t density) { if (density>1e-6) fDensitySupportOverSi=density; else fDensitySupportOverSi=1e-6; }
  
  //--------- for underlying and pile-up events --------------------
  
  void SetFileNameForUnderlyingEvent(TString fileName) { if (fileName.EndsWith("MFT.RecPoints.root")) fFileNameForUnderyingEvent += fileName; }
  void SetFileNameForPileUpEvents(TString fileName)    { if (fileName.EndsWith("MFT.RecPoints.root")) fFileNameForPileUpEvents   += fileName; }
  
  void SetUnderlyingEventID(Short_t eventID) { fUnderlyingEventID = eventID; }
  void SetPileUpEventID(Short_t i, Short_t eventID) { if (i>=0 && i<AliMFTConstants::fNMaxPileUpEvents) fPileUpEventsIDs[i] = eventID; }
  
  const Char_t* GetFileNameForUnderlyingEvent() { return fFileNameForUnderyingEvent; }
  const Char_t* GetFileNameForPileUpEvents()    { return fFileNameForPileUpEvents; }
  Short_t GetUnderlyingEventID() { return fUnderlyingEventID; }
  Short_t GetPileUpEventID(Short_t i) { if (i>=0 && i<AliMFTConstants::fNMaxPileUpEvents) return fPileUpEventsIDs[i]; else return -1; }
  
protected:
  
  static const Int_t fNMaxPlanes = AliMFTConstants::fNMaxPlanes;        // max number of MFT planes
  
  Int_t fVersion;
  
  Int_t fNPlanes;                             // # of MFT planes
  
  TObjArray *fSDigitsPerPlane;                // ![fNPlanes] list of sdigits [per plane]
  TObjArray *fDigitsPerPlane;                 // ![fNPlanes] list of digits [per plane]
  TObjArray *fRecPointsPerPlane;              // ![fNPlanes] list of recPoints [per plane]
  TClonesArray *fSideDigits;                  // ! list of digits fired by the charge dispersion coming from the main hit
  
  AliMFTSegmentation *fSegmentation;
  
  TString fNameGeomFile;
  
  Double_t fChargeDispersion;
  Double_t fSingleStepForChargeDispersion;
  Int_t fNStepForChargeDispersion;
  
  Double_t fDensitySupportOverSi;
  
  TString fFileNameForUnderyingEvent, fFileNameForPileUpEvents;
  Short_t fNPileUpEvents, fUnderlyingEventID, fPileUpEventsIDs[AliMFTConstants::fNMaxPileUpEvents];
  
private:
  
  AliMFT (const AliMFT& mft);             // dummy copy constructor
  AliMFT &operator=(const AliMFT& mft);   // dummy assignment operator
  
  ClassDef(AliMFT,2)
  
};

//====================================================================================================================================================

#endif

