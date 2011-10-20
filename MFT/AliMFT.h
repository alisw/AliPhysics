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
#include "TClonesArray.h"
#include "TGeoGlobalMagField.h"
#include "AliRun.h"
#include "AliLoader.h"
#include "AliDetector.h"
#include "AliMC.h"
#include "AliMagF.h"
#include "AliMFT.h"
#include "AliMFTHit.h"
#include "AliMFTDigit.h"
#include "AliMFTCluster.h"
#include "AliTrackReference.h"
#include "AliMFTSegmentation.h"
#include "AliMFTDigitizer.h"
#include "AliMFTPlane.h"
#include "TString.h"
#include "TObjArray.h"

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

  TObjArray*    GetSDigitsList()            const { return fSDigitsPerPlane; }                                                  // get sdigits list for all planes
  TClonesArray* GetSDigitsList(Int_t plane) const { return fSDigitsPerPlane ? (TClonesArray*) fSDigitsPerPlane->At(plane):0; }  // get sdigits list for a plane

  TObjArray*    GetDigitsList()            const{return fDigitsPerPlane;}                                                   // get digits list for all layers
  TClonesArray* GetDigitsList(Int_t plane) const{return fDigitsPerPlane ? (TClonesArray*) fDigitsPerPlane->At(plane):0; }   // get digits list for a plane

  TObjArray*    GetRecPointsList()            const{return fRecPointsPerPlane;}                                                      // get digits list for all layers
  TClonesArray* GetRecPointsList(Int_t plane) const{return fRecPointsPerPlane ? (TClonesArray*) fRecPointsPerPlane->At(plane):0; }   // get digits list for a plane

  void ResetSDigits()   { if(fSDigitsPerPlane)   for(int iPlane=0; iPlane<fNPlanes; iPlane++) ((TClonesArray*) fSDigitsPerPlane  ->At(iPlane))->Clear(); }   // reset sdigits list  
  void ResetDigits()    { if(fDigitsPerPlane)    for(int iPlane=0; iPlane<fNPlanes; iPlane++) ((TClonesArray*) fDigitsPerPlane   ->At(iPlane))->Clear(); }   // reset digits list
  void ResetRecPoints() { if(fRecPointsPerPlane) for(int iPlane=0; iPlane<fNPlanes; iPlane++) ((TClonesArray*) fRecPointsPerPlane->At(iPlane))->Clear(); }   // reset recPoints list
  
  AliDigitizer* CreateDigitizer(AliRunDigitizer *manager) const { return new AliMFTDigitizer(manager); }   // from AliModule invoked from AliSimulation::RunDigitization()
  
  AliMFTSegmentation* GetSegmentation() const { return fSegmentation; }

  enum EMedia{kAir, kSi, kReadout, kSupport};  // media IDs used in CreateMaterials  
    
  // Geometry/segmentation creation part
  TGeoVolumeAssembly* CreateVol();
  void AddAlignableVolumes() const { /* not needed */ return; }
  void SetGeometry();

  void  SetNSlices(Int_t nSlices) { fNSlices = nSlices; }
  Int_t GetNSlices() const { return fNSlices; }

  Int_t GetNPlanes() const { return fNPlanes; }

  void SetChargeDispersion(Double_t chargeDispersion) { fChargeDispersion = chargeDispersion; }
  Double_t GetChargeDispersion() { return fChargeDispersion; }
  void SetNStepForChargeDispersion(Int_t nStepForChargeDispersion) { fNStepForChargeDispersion = nStepForChargeDispersion; }
  Int_t GetNStepForChargeDispersion() { return fNStepForChargeDispersion; }
  Double_t GetSingleStepForChargeDispersion() { return fSingleStepForChargeDispersion; }

  void SetDensitySiOverSupport(Double_t density) { fDensitySiOverSupport = density; }

protected:

  static const Int_t fNMaxPlanes = 20;        // max number of MFT planes

  Int_t fVersion;

  Int_t fNPlanes;                             // # of MFT planes
  Int_t fNSlices;                             // # of slices per MFT plane
  
  TObjArray *fSDigitsPerPlane;                // ![fNPlanes] list of sdigits [per plane]
  TObjArray *fDigitsPerPlane;                 // ![fNPlanes] list of digits [per plane]
  TObjArray *fRecPointsPerPlane;              // ![fNPlanes] list of recPoints [per plane]
  TClonesArray *fSideDigits;                  // ! list of digits fired by the charge dispersion coming from the main hit

  AliMFTSegmentation *fSegmentation;

  TString fNameGeomFile;

  Double_t fChargeDispersion;
  Double_t fSingleStepForChargeDispersion;
  Int_t fNStepForChargeDispersion;

  Double_t fDensitySiOverSupport;
  
private:

  AliMFT (const AliMFT& mft);             // dummy copy constructor
  AliMFT &operator=(const AliMFT& mft);   // dummy assignment operator
  
  ClassDef(AliMFT,1) 
    
};

//====================================================================================================================================================

#endif

