#ifndef AliMFTPlane_H
#define AliMFTPlane_H 

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//====================================================================================================================================================
//
//      Class for the description of the structure for the planes of the ALICE Muon Forward Tracker
//
//      Contact author: antonio.uras@cern.ch
//
//====================================================================================================================================================

#include "TNamed.h"
#include "THnSparse.h"
#include "TClonesArray.h"
#include "TAxis.h"
#include "TPave.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TEllipse.h"
#include "TMath.h"
#include "AliLog.h"

//====================================================================================================================================================

class AliMFTPlane : public TNamed {

public:

  AliMFTPlane();
  AliMFTPlane(const Char_t *name, const Char_t *title);
  AliMFTPlane(const AliMFTPlane& pt);
  AliMFTPlane& operator=(const AliMFTPlane &source);

  virtual ~AliMFTPlane() {};  // destructor

  Bool_t Init(Int_t    planeNumber,
	      Double_t zCenter, 
	      Double_t rMin, 
	      Double_t rMax, 
	      Double_t pixelSizeX, 
	      Double_t pixelSizeY, 
	      Double_t thicknessActive, 
	      Double_t thicknessSupport, 
	      Double_t thicknessReadout);
  
  Bool_t CreateStructure();

  Int_t GetNActiveElements()  { return fActiveElements->GetEntries();  }
  Int_t GetNReadoutElements() { return fReadoutElements->GetEntries(); }
  Int_t GetNSupportElements() { return fSupportElements->GetEntries(); }

  TClonesArray* GetActiveElements()  { return fActiveElements;  }
  TClonesArray* GetReadoutElements() { return fReadoutElements; }
  TClonesArray* GetSupportElements() { return fSupportElements; }

  THnSparseC* GetActiveElement(Int_t id);
  THnSparseC* GetReadoutElement(Int_t id);
  THnSparseC* GetSupportElement(Int_t id);

  Bool_t IsFront(THnSparseC *element) { return (element->GetAxis(2)->GetXmin() < fZCenter); }

  void DrawPlane(Char_t *opt="");

  Double_t GetRMinSupport() { return fRMinSupport; }
  Double_t GetRMaxSupport() { return fRMaxSupport; }
  Double_t GetThicknessSupport() { return GetSupportElement(0)->GetAxis(2)->GetXmax() - GetSupportElement(0)->GetAxis(2)->GetXmin(); }
  
  Double_t GetZCenter() { return fZCenter; }
  Double_t GetZCenterActiveFront() { return fZCenterActiveFront; }
  Double_t GetZCenterActiveBack()  { return fZCenterActiveBack; }

  void SetEquivalentSilicon(Double_t equivalentSilicon)                       { fEquivalentSilicon            = equivalentSilicon; }
  void SetEquivalentSiliconBeforeFront(Double_t equivalentSiliconBeforeFront) { fEquivalentSiliconBeforeFront = equivalentSiliconBeforeFront; }
  void SetEquivalentSiliconBeforeBack(Double_t equivalentSiliconBeforeBack)   { fEquivalentSiliconBeforeBack  = equivalentSiliconBeforeBack; }
  Double_t GetEquivalentSilicon()            { return fEquivalentSilicon; }
  Double_t GetEquivalentSiliconBeforeFront() { return fEquivalentSiliconBeforeFront; }
  Double_t GetEquivalentSiliconBeforeBack()  { return fEquivalentSiliconBeforeBack; }
  
private:

  // measures in cm

  static const Double_t fRadiusMin = 2.225;           // minimum radial distance of the MFT sensors. To be carefully coordinated with fActiveSuperposition

  static const Double_t fActiveSuperposition = 0.05;  // superposition between the active elements tasselling the MFT planes, for having a 
                                                      // full acceptance coverage even in case of 10 degrees inclined tracks
  static const Double_t fHeightActive = 0.5;          // height of the active elements
  static const Double_t fHeightReadout = 0.3;         // height of the readout elements attached to the active ones

  static const Double_t fSupportExtMargin = 0.3 + 0.3;    // minimum border size between the end of the support plane and the sensors: fHeightReadout + 0.3

  Int_t fPlaneNumber;

  Double_t fZCenter, fRMinSupport, fRMax, fRMaxSupport, fPixelSizeX, fPixelSizeY, fThicknessActive, fThicknessSupport, fThicknessReadout;
  Double_t fZCenterActiveFront, fZCenterActiveBack, fEquivalentSilicon, fEquivalentSiliconBeforeFront, fEquivalentSiliconBeforeBack;

  TClonesArray *fActiveElements, *fReadoutElements, *fSupportElements;

  ClassDef(AliMFTPlane, 1)

};

//====================================================================================================================================================
	
#endif

