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
  
  virtual ~AliMFTPlane();  // destructor
  virtual void Clear(const Option_t* /*opt*/);
  
  Bool_t Init(Int_t    planeNumber,
	      Double_t zCenter, 
	      Double_t rMin, 
	      Double_t rMax, 
	      Double_t pixelSizeX, 
	      Double_t pixelSizeY, 
	      Double_t thicknessActive, 
	      Double_t thicknessSupport, 
	      Double_t thicknessReadout,
	      Bool_t   hasPixelRectangularPatternAlongY);
  
  Bool_t CreateStructure();

  Int_t GetNActiveElements()  const { return fActiveElements->GetEntries();  }
  Int_t GetNReadoutElements() const { return fReadoutElements->GetEntries(); }
  Int_t GetNSupportElements() const { return fSupportElements->GetEntries(); }

  TClonesArray* GetActiveElements()  { return fActiveElements;  }
  TClonesArray* GetReadoutElements() { return fReadoutElements; }
  TClonesArray* GetSupportElements() { return fSupportElements; }

  THnSparseC* GetActiveElement(Int_t id);
  THnSparseC* GetReadoutElement(Int_t id);
  THnSparseC* GetSupportElement(Int_t id);

  Bool_t IsFront(THnSparseC *element) const { return (element->GetAxis(2)->GetXmin() < fZCenter); }

  void DrawPlane(Option_t *opt="");

  Double_t GetRMinSupport() const { return fRMinSupport; }
  Double_t GetRMaxSupport() const { return fRMaxSupport; }
  Double_t GetThicknessSupport() { return GetSupportElement(0)->GetAxis(2)->GetXmax() - GetSupportElement(0)->GetAxis(2)->GetXmin(); }
  
  Double_t GetZCenter()            const { return fZCenter; }
  Double_t GetZCenterActiveFront() const { return fZCenterActiveFront; }
  Double_t GetZCenterActiveBack()  const { return fZCenterActiveBack; }

  void SetEquivalentSilicon(Double_t equivalentSilicon)                       { fEquivalentSilicon            = equivalentSilicon; }
  void SetEquivalentSiliconBeforeFront(Double_t equivalentSiliconBeforeFront) { fEquivalentSiliconBeforeFront = equivalentSiliconBeforeFront; }
  void SetEquivalentSiliconBeforeBack(Double_t equivalentSiliconBeforeBack)   { fEquivalentSiliconBeforeBack  = equivalentSiliconBeforeBack; }
  Double_t GetEquivalentSilicon()            const { return fEquivalentSilicon; }
  Double_t GetEquivalentSiliconBeforeFront() const { return fEquivalentSiliconBeforeFront; }
  Double_t GetEquivalentSiliconBeforeBack()  const { return fEquivalentSiliconBeforeBack; }

  Int_t GetNumberOfChips(Option_t *opt);
  Bool_t HasPixelRectangularPatternAlongY() { return fHasPixelRectangularPatternAlongY; }
  
private:

  // measures in cm

  static const Double_t fRadiusMin;            // minimum radial distance of the MFT sensors. To be carefully coordinated with fActiveSuperposition

  static const Double_t fActiveSuperposition;  // superposition between the active elements tasselling the MFT planes, for having a 
                                               // full acceptance coverage even in case of 10 degrees inclined tracks
  static const Double_t fHeightActive;         // height of the active elements
  static const Double_t fHeightReadout;        // height of the readout elements attached to the active ones

  static const Double_t fSupportExtMargin;     // minimum border size between the end of the support plane and the sensors: fHeightReadout + 0.3

  Int_t fPlaneNumber;

  Double_t fZCenter, fRMinSupport, fRMax, fRMaxSupport, fPixelSizeX, fPixelSizeY, fThicknessActive, fThicknessSupport, fThicknessReadout;
  Double_t fZCenterActiveFront, fZCenterActiveBack, fEquivalentSilicon, fEquivalentSiliconBeforeFront, fEquivalentSiliconBeforeBack;

  TClonesArray *fActiveElements, *fReadoutElements, *fSupportElements;

  Bool_t fHasPixelRectangularPatternAlongY, fPlaneIsOdd;

  ClassDef(AliMFTPlane, 1)

};

//====================================================================================================================================================
	
#endif

