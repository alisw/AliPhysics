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
#include "AliMFTConstants.h"
#include "AliMFTPlane.h"

const Double_t AliMFTPlane::fRadiusMin           = AliMFTConstants::fRadiusMin;
const Double_t AliMFTPlane::fActiveSuperposition = AliMFTConstants::fActiveSuperposition;
const Double_t AliMFTPlane::fHeightActive        = AliMFTConstants::fHeightActive;
const Double_t AliMFTPlane::fHeightReadout       = AliMFTConstants::fHeightReadout;
const Double_t AliMFTPlane::fSupportExtMargin    = AliMFTConstants::fSupportExtMargin;

ClassImp(AliMFTPlane)

//====================================================================================================================================================

AliMFTPlane::AliMFTPlane():
  TNamed(),
  fPlaneNumber(-1),
  fZCenter(0), 
  fRMinSupport(0), 
  fRMax(0),
  fRMaxSupport(0),
  fPixelSizeX(0), 
  fPixelSizeY(0), 
  fThicknessActive(0), 
  fThicknessSupport(0), 
  fThicknessReadout(0),
  fZCenterActiveFront(0),
  fZCenterActiveBack(0),
  fEquivalentSilicon(0),
  fEquivalentSiliconBeforeFront(0),
  fEquivalentSiliconBeforeBack(0),
  fActiveElements(0),
  fReadoutElements(0),
  fSupportElements(0),
  fHasPixelRectangularPatternAlongY(kFALSE),
  fPlaneIsOdd(kFALSE)
{

  // default constructor

}

//====================================================================================================================================================

AliMFTPlane::AliMFTPlane(const Char_t *name, const Char_t *title):
  TNamed(name, title),
  fPlaneNumber(-1),
  fZCenter(0), 
  fRMinSupport(0), 
  fRMax(0),
  fRMaxSupport(0),
  fPixelSizeX(0), 
  fPixelSizeY(0), 
  fThicknessActive(0), 
  fThicknessSupport(0), 
  fThicknessReadout(0),
  fZCenterActiveFront(0),
  fZCenterActiveBack(0),
  fEquivalentSilicon(0),
  fEquivalentSiliconBeforeFront(0),
  fEquivalentSiliconBeforeBack(0),
  fActiveElements(0),
  fReadoutElements(0),
  fSupportElements(0),
  fHasPixelRectangularPatternAlongY(kFALSE),
  fPlaneIsOdd(kFALSE)
{

  // constructor
  fActiveElements  = new TClonesArray("THnSparseC");
  fReadoutElements = new TClonesArray("THnSparseC");
  fSupportElements = new TClonesArray("THnSparseC");
  fActiveElements->SetOwner(kTRUE);
  fReadoutElements->SetOwner(kTRUE);
  fSupportElements->SetOwner(kTRUE);
  
}

//====================================================================================================================================================

AliMFTPlane::AliMFTPlane(const AliMFTPlane& plane):
  TNamed(plane),
  fPlaneNumber(plane.fPlaneNumber),
  fZCenter(plane.fZCenter), 
  fRMinSupport(plane.fRMinSupport), 
  fRMax(plane.fRMax),
  fRMaxSupport(plane.fRMaxSupport),
  fPixelSizeX(plane.fPixelSizeX), 
  fPixelSizeY(plane.fPixelSizeY), 
  fThicknessActive(plane.fThicknessActive), 
  fThicknessSupport(plane.fThicknessSupport), 
  fThicknessReadout(plane.fThicknessReadout),
  fZCenterActiveFront(plane.fZCenterActiveFront),
  fZCenterActiveBack(plane.fZCenterActiveBack),
  fEquivalentSilicon(plane.fEquivalentSilicon),
  fEquivalentSiliconBeforeFront(plane.fEquivalentSiliconBeforeFront),
  fEquivalentSiliconBeforeBack(plane.fEquivalentSiliconBeforeBack),
  fActiveElements(0),
  fReadoutElements(0),
  fSupportElements(0),
  fHasPixelRectangularPatternAlongY(plane.fHasPixelRectangularPatternAlongY),
  fPlaneIsOdd(plane.fPlaneIsOdd)
{

  // copy constructor
  fActiveElements  = new TClonesArray(*(plane.fActiveElements));
  fActiveElements  -> SetOwner(kTRUE);
  fReadoutElements = new TClonesArray(*(plane.fReadoutElements));
  fReadoutElements -> SetOwner(kTRUE);
  fSupportElements = new TClonesArray(*(plane.fSupportElements));
  fSupportElements -> SetOwner(kTRUE);

	
}

//====================================================================================================================================================

AliMFTPlane::~AliMFTPlane() {

  AliInfo("Delete AliMFTPlane");
  if(fActiveElements) fActiveElements->Delete();
  delete fActiveElements; 
  if(fReadoutElements) fReadoutElements->Delete();
  delete fReadoutElements; 
  if(fSupportElements) fSupportElements->Delete();
  delete fSupportElements; 

}

//====================================================================================================================================================

void AliMFTPlane::Clear(const Option_t* /*opt*/) {

  AliInfo("Clear AliMFTPlane");
  if(fActiveElements) fActiveElements->Delete();
  delete fActiveElements; fActiveElements=NULL;
  if(fReadoutElements) fReadoutElements->Delete();
  delete fReadoutElements;  fReadoutElements=NULL; 
  if(fSupportElements) fSupportElements->Delete();
  delete fSupportElements;   fSupportElements=NULL;

}

//====================================================================================================================================================

AliMFTPlane& AliMFTPlane::operator=(const AliMFTPlane& plane) {

  // Assignment operator
  
  // check assignement to self
  if (this != &plane) {
    
    // base class assignement
    TNamed::operator=(plane);
    
    // clear memory
    Clear("");
    
    fPlaneNumber                      = plane.fPlaneNumber;
    fZCenter                          = plane.fZCenter; 
    fRMinSupport                      = plane.fRMinSupport; 
    fRMax                             = plane.fRMax;
    fRMaxSupport                      = plane.fRMaxSupport;
    fPixelSizeX                       = plane.fPixelSizeX;
    fPixelSizeY                       = plane.fPixelSizeY; 
    fThicknessActive                  = plane.fThicknessActive; 
    fThicknessSupport                 = plane.fThicknessSupport; 
    fThicknessReadout                 = plane.fThicknessReadout;
    fZCenterActiveFront               = plane.fZCenterActiveFront;
    fZCenterActiveBack                = plane.fZCenterActiveBack;
    fEquivalentSilicon                = plane.fEquivalentSilicon;
    fEquivalentSiliconBeforeFront     = plane.fEquivalentSiliconBeforeFront;
    fEquivalentSiliconBeforeBack      = plane.fEquivalentSiliconBeforeBack;
    fActiveElements = new TClonesArray(*(plane.fActiveElements));
    fActiveElements -> SetOwner(kTRUE);
    fReadoutElements = new TClonesArray(*(plane.fReadoutElements));
    fReadoutElements -> SetOwner(kTRUE);
    fSupportElements = new TClonesArray(*(plane.fSupportElements));
    fSupportElements -> SetOwner(kTRUE);
    fHasPixelRectangularPatternAlongY = plane.fHasPixelRectangularPatternAlongY;
    fPlaneIsOdd                       = plane.fPlaneIsOdd;

  }
  
  return *this;
  
}

//====================================================================================================================================================

Bool_t AliMFTPlane::Init(Int_t    planeNumber,
			 Double_t zCenter, 
			 Double_t rMin, 
			 Double_t rMax, 
			 Double_t pixelSizeX, 
			 Double_t pixelSizeY, 
			 Double_t thicknessActive, 
			 Double_t thicknessSupport, 
			 Double_t thicknessReadout,
			 Bool_t   hasPixelRectangularPatternAlongY) {

  AliDebug(1, Form("Initializing Plane Structure for Plane %s", GetName()));

  fPlaneNumber      = planeNumber;
  fZCenter          = zCenter;
  fRMinSupport      = rMin;
  fRMax             = rMax;
  fPixelSizeX       = pixelSizeX;
  fPixelSizeY       = pixelSizeY;
  fThicknessActive  = thicknessActive;
  fThicknessSupport = thicknessSupport;
  fThicknessReadout = thicknessReadout;

  fHasPixelRectangularPatternAlongY = hasPixelRectangularPatternAlongY;

  fZCenterActiveFront = fZCenter - 0.5*fThicknessSupport - 0.5*fThicknessActive;
  fZCenterActiveBack  = fZCenter + 0.5*fThicknessSupport + 0.5*fThicknessActive;

//   if (fRMinSupport <= fRadiusMin) fRMinSupport = fRadiusMin;
//   else {
//     fRMinSupport = fRadiusMin + (fHeightActive-fActiveSuperposition) * Int_t((fRMinSupport-fRadiusMin)/(fHeightActive-fActiveSuperposition));
//   }
  
  if (fRMax < fRMinSupport+fHeightActive) fRMax = fRMinSupport + fHeightActive;

  Int_t nLaddersWithinPipe = Int_t(fRMinSupport/(fHeightActive-fActiveSuperposition));
  if (fRMinSupport-nLaddersWithinPipe*(fHeightActive-fActiveSuperposition) > 0.5*(fHeightActive-2*fActiveSuperposition)) fPlaneIsOdd = kTRUE;
  else fPlaneIsOdd = kFALSE;

  fRMax = fRMinSupport + (fHeightActive-fActiveSuperposition) * 
    (Int_t((fRMax-fRMinSupport-fHeightActive)/(fHeightActive-fActiveSuperposition))+1) + fHeightActive;

  fRMaxSupport = TMath::Sqrt(fHeightActive*(2.*rMax-fHeightActive) + fRMax*fRMax) + fSupportExtMargin;
   
  return kTRUE;
 
}

//====================================================================================================================================================

Bool_t AliMFTPlane::CreateStructure() {

  Int_t nBins[3]={0};
  Double_t minPosition[3]={0}, maxPosition[3]={0};
  
  // ------------------- det elements: active + readout ----------------------------------

  // 1st Section : below and above the beam pipe

  Double_t lowEdgeActive = -1.*fRMax;
  Double_t supEdgeActive = lowEdgeActive + fHeightActive;
  Double_t zMinFront = fZCenter - 0.5*fThicknessSupport - fThicknessActive;
  Double_t zMinBack  = fZCenter + 0.5*fThicknessSupport;
  Double_t zMin = 0.;
  Bool_t isFront = kTRUE;
  
  while (lowEdgeActive < 0) {
    
    Double_t extLimitAtLowEdgeActive = TMath::Sqrt((fRMax-TMath::Abs(lowEdgeActive)) * TMath::Abs(2*fRMax - (fRMax-TMath::Abs(lowEdgeActive))));
    Double_t extLimitAtSupEdgeActive = TMath::Sqrt((fRMax-TMath::Abs(supEdgeActive)) * TMath::Abs(2*fRMax - (fRMax-TMath::Abs(supEdgeActive))));

    // creating new det element: active + readout
    
    Double_t extLimitDetElem = TMath::Max(extLimitAtLowEdgeActive, extLimitAtSupEdgeActive);
    
    if (supEdgeActive<-1.*fRMinSupport+0.01 || lowEdgeActive>1.*fRMinSupport-0.01) {     // single element covering the row
      
      nBins[0] = TMath::Nint(2.*extLimitDetElem/fPixelSizeX);
      nBins[1] = TMath::Nint(fHeightActive/fPixelSizeY);
      nBins[2] = 1;

      // element below the pipe
      
      if (isFront) zMin = zMinFront;
      else         zMin = zMinBack;

      minPosition[0] = -1.*extLimitDetElem;
      minPosition[1] = lowEdgeActive;
      minPosition[2] = zMin;
      
      maxPosition[0] = +1.*extLimitDetElem;
      maxPosition[1] = supEdgeActive;
      maxPosition[2] = zMin+fThicknessActive; 
      
      new ((*fActiveElements)[fActiveElements->GetEntries()]) THnSparseC(Form("MFTActiveElemHist_%02d%03d", fPlaneNumber, fActiveElements->GetEntries()), 
									 Form("MFTActiveElemHist_%02d%03d", fPlaneNumber, fActiveElements->GetEntries()), 
									 3, nBins, minPosition, maxPosition);

      minPosition[1] = lowEdgeActive-fHeightReadout;
      maxPosition[1] = lowEdgeActive;
      
      new ((*fReadoutElements)[fReadoutElements->GetEntries()]) THnSparseC(Form("MFTReadoutElemHist_%02d%03d", fPlaneNumber, fReadoutElements->GetEntries()), 
									   Form("MFTReadoutElemHist_%02d%03d", fPlaneNumber, fReadoutElements->GetEntries()), 
									   3, nBins, minPosition, maxPosition);

      // specular element above the pipe

      if (fPlaneIsOdd) {
	if (isFront) zMin = zMinBack;
	else         zMin = zMinFront;
      }

      minPosition[0] = -1.*extLimitDetElem;
      minPosition[1] = -1.*supEdgeActive;
      minPosition[2] = zMin;
      
      maxPosition[0] = +1.*extLimitDetElem;
      maxPosition[1] = -1.*lowEdgeActive;
      maxPosition[2] = zMin+fThicknessActive; 
      
      new ((*fActiveElements)[fActiveElements->GetEntries()]) THnSparseC(Form("MFTActiveElemHist_%02d%03d", fPlaneNumber, fActiveElements->GetEntries()), 
									 Form("MFTActiveElemHist_%02d%03d", fPlaneNumber, fActiveElements->GetEntries()), 
									 3, nBins, minPosition, maxPosition);

      minPosition[1] = -1.*lowEdgeActive;
      maxPosition[1] = -1.*(lowEdgeActive-fHeightReadout);

      new ((*fReadoutElements)[fReadoutElements->GetEntries()]) THnSparseC(Form("MFTReadoutElemHist_%02d%03d", fPlaneNumber, fReadoutElements->GetEntries()), 
									   Form("MFTReadoutElemHist_%02d%03d", fPlaneNumber, fReadoutElements->GetEntries()), 
									   3, nBins, minPosition, maxPosition);

    }
    
    else {     // two elements covering the row
      
      Double_t intLimitAtLowEdge = 0., intLimitAtSupEdge = 0.;
      if (fRMinSupport-TMath::Abs(lowEdgeActive)>0.) intLimitAtLowEdge = TMath::Sqrt((fRMinSupport-TMath::Abs(lowEdgeActive)) * TMath::Abs(2*fRMinSupport - (fRMinSupport-TMath::Abs(lowEdgeActive))));
      if (fRMinSupport-TMath::Abs(supEdgeActive)>0.) intLimitAtSupEdge = TMath::Sqrt((fRMinSupport-TMath::Abs(supEdgeActive)) * TMath::Abs(2*fRMinSupport - (fRMinSupport-TMath::Abs(supEdgeActive))));
      Double_t intLimitDetElem = TMath::Max(intLimitAtLowEdge, intLimitAtSupEdge);
      
      nBins[0] = TMath::Nint((extLimitDetElem-intLimitDetElem)/fPixelSizeX);
      nBins[1] = TMath::Nint(fHeightActive/fPixelSizeY);
      nBins[2] = 1;
      
      // left element: below the beam pipe
      
      if (isFront) zMin = zMinFront;
      else         zMin = zMinBack;

      minPosition[0] = -1.*extLimitDetElem;
      minPosition[1] = lowEdgeActive;
      minPosition[2] = zMin;
      
      maxPosition[0] = -1.*intLimitDetElem;
      maxPosition[1] = supEdgeActive;
      maxPosition[2] = zMin+fThicknessActive; 
      
      new ((*fActiveElements)[fActiveElements->GetEntries()]) THnSparseC(Form("MFTActiveElemHist_%02d%03d", fPlaneNumber, fActiveElements->GetEntries()), 
									 Form("MFTActiveElemHist_%02d%03d", fPlaneNumber, fActiveElements->GetEntries()), 
									 3, nBins, minPosition, maxPosition);	
      
      minPosition[1] = lowEdgeActive-fHeightReadout;
      maxPosition[1] = lowEdgeActive;
      
      new ((*fReadoutElements)[fReadoutElements->GetEntries()]) THnSparseC(Form("MFTReadoutElemHist_%02d%03d", fPlaneNumber, fReadoutElements->GetEntries()), 
									   Form("MFTReadoutElemHist_%02d%03d", fPlaneNumber, fReadoutElements->GetEntries()), 
									   3, nBins, minPosition, maxPosition);

      // left element: above the beam pipe
      
      if (supEdgeActive < 0.5*fHeightActive) {
	
	if (fPlaneIsOdd) {
	  if (isFront) zMin = zMinBack;
	  else         zMin = zMinFront;
	}
	
	minPosition[0] = -1.*extLimitDetElem;
	minPosition[1] = -1.*supEdgeActive;
	minPosition[2] = zMin;
	
	maxPosition[0] = -1.*intLimitDetElem;
	maxPosition[1] = -1.*lowEdgeActive;
	maxPosition[2] = zMin+fThicknessActive; 
	
	new ((*fActiveElements)[fActiveElements->GetEntries()]) THnSparseC(Form("MFTActiveElemHist_%02d%03d", fPlaneNumber, fActiveElements->GetEntries()), 
									   Form("MFTActiveElemHist_%02d%03d", fPlaneNumber, fActiveElements->GetEntries()), 
									   3, nBins, minPosition, maxPosition);	
	
	minPosition[1] = -1.*lowEdgeActive;
	maxPosition[1] = -1.*(lowEdgeActive-fHeightReadout);
	
	new ((*fReadoutElements)[fReadoutElements->GetEntries()]) THnSparseC(Form("MFTReadoutElemHist_%02d%03d", fPlaneNumber, fReadoutElements->GetEntries()), 
									     Form("MFTReadoutElemHist_%02d%03d", fPlaneNumber, fReadoutElements->GetEntries()), 
									     3, nBins, minPosition, maxPosition);
      
      }

      // right element: below the beam pipe
      
      if (isFront) zMin = zMinFront;
      else         zMin = zMinBack;

      minPosition[0] = +1.*intLimitDetElem;
      minPosition[1] = lowEdgeActive;
      minPosition[2] = zMin;
      
      maxPosition[0] = +1.*extLimitDetElem;
      maxPosition[1] = supEdgeActive;
      maxPosition[2] = zMin+fThicknessActive; 
      
      new ((*fActiveElements)[fActiveElements->GetEntries()]) THnSparseC(Form("MFTActiveElemHist_%02d%03d", fPlaneNumber, fActiveElements->GetEntries()), 
									 Form("MFTActiveElemHist_%02d%03d", fPlaneNumber, fActiveElements->GetEntries()), 
									 3, nBins, minPosition, maxPosition);	
      
      minPosition[1] = lowEdgeActive-fHeightReadout;
      maxPosition[1] = lowEdgeActive;

      new ((*fReadoutElements)[fReadoutElements->GetEntries()]) THnSparseC(Form("MFTReadoutElemHist_%02d%03d", fPlaneNumber, fReadoutElements->GetEntries()), 
									   Form("MFTReadoutElemHist_%02d%03d", fPlaneNumber, fReadoutElements->GetEntries()), 
									   3, nBins, minPosition, maxPosition);

      // right element: above the beam pipe
      
      if (supEdgeActive < 0.5*fHeightActive) {

	if (fPlaneIsOdd) {
	  if (isFront) zMin = zMinBack;
	  else         zMin = zMinFront;
	}
	
	minPosition[0] = +1.*intLimitDetElem;
	minPosition[1] = -1.*supEdgeActive;
	minPosition[2] = zMin;
	
	maxPosition[0] = +1.*extLimitDetElem;
	maxPosition[1] = -1.*lowEdgeActive;
	maxPosition[2] = zMin+fThicknessActive; 
	
	new ((*fActiveElements)[fActiveElements->GetEntries()]) THnSparseC(Form("MFTActiveElemHist_%02d%03d", fPlaneNumber, fActiveElements->GetEntries()), 
									   Form("MFTActiveElemHist_%02d%03d", fPlaneNumber, fActiveElements->GetEntries()), 
									   3, nBins, minPosition, maxPosition);	
	
	minPosition[1] = -1.*lowEdgeActive;
	maxPosition[1] = -1.*(lowEdgeActive-fHeightReadout);
	
	new ((*fReadoutElements)[fReadoutElements->GetEntries()]) THnSparseC(Form("MFTReadoutElemHist_%02d%03d", fPlaneNumber, fReadoutElements->GetEntries()), 
									     Form("MFTReadoutElemHist_%02d%03d", fPlaneNumber, fReadoutElements->GetEntries()), 
									     3, nBins, minPosition, maxPosition);

      }
      
    }
    
    lowEdgeActive += fHeightActive - fActiveSuperposition;
    supEdgeActive = lowEdgeActive + fHeightActive;
    isFront = !isFront;
    
  }
  
  // ------------------- support element -------------------------------------------------
  
  nBins[0] = 1;
  nBins[1] = 1;
  nBins[2] = 1;
  
  minPosition[0] = -1.*fRMaxSupport;
  minPosition[1] = -1.*fRMaxSupport;
  minPosition[2] = fZCenter - 0.5*fThicknessSupport;
  
  maxPosition[0] = +1.*fRMaxSupport;
  maxPosition[1] = +1.*fRMaxSupport;
  maxPosition[2] = fZCenter + 0.5*fThicknessSupport;
  
  new ((*fSupportElements)[fSupportElements->GetEntries()]) THnSparseC(Form("MFTSupportElemHist_%02d%03d", fPlaneNumber, fSupportElements->GetEntries()), 
								       Form("MFTSupportElemHist_%02d%03d", fPlaneNumber, fSupportElements->GetEntries()), 
								       3, nBins, minPosition, maxPosition);

  // --------------------------------------------------------------------------------------

  AliDebug(1, Form("Structure completed for MFT plane %s", GetName()));

  return kTRUE;
  
}

//====================================================================================================================================================

THnSparseC* AliMFTPlane::GetActiveElement(Int_t id) {

  if (id<0 || id>=GetNActiveElements()) return NULL;
  else return (THnSparseC*) fActiveElements->At(id);

}

//====================================================================================================================================================

THnSparseC* AliMFTPlane::GetReadoutElement(Int_t id) {

  if (id<0 || id>=GetNReadoutElements()) return NULL;
  else return (THnSparseC*) fReadoutElements->At(id);

}

//====================================================================================================================================================

THnSparseC* AliMFTPlane::GetSupportElement(Int_t id) {

  if (id<0 || id>=GetNSupportElements()) return NULL;
  else return (THnSparseC*) fSupportElements->At(id);

}

//====================================================================================================================================================

void AliMFTPlane::DrawPlane(Option_t *opt) {

  // ------------------- "FRONT" option ------------------

  if (!strcmp(opt, "front")) {

    TCanvas *cnv = new TCanvas("cnv", GetName(), 900, 900);
    cnv->Draw();

    TH2D *h = new TH2D("tmp", GetName(), 
		       1, 1.1*GetSupportElement(0)->GetAxis(0)->GetXmin(), 1.1*GetSupportElement(0)->GetAxis(0)->GetXmax(), 
		       1, 1.1*GetSupportElement(0)->GetAxis(1)->GetXmin(), 1.1*GetSupportElement(0)->GetAxis(1)->GetXmax());
    h->SetXTitle("x [cm]");
    h->SetYTitle("y [cm]");
    h->Draw();

    AliInfo("Created hist");

    TEllipse *supportExt = new TEllipse(0.0, 0.0, fRMaxSupport, fRMaxSupport);
    TEllipse *supportInt = new TEllipse(0.0, 0.0, fRMinSupport, fRMinSupport);
    supportExt->SetFillColor(kCyan-10);
    supportExt -> Draw("same");
    supportInt -> Draw("same");

    for (Int_t iEl=0; iEl<GetNActiveElements(); iEl++) {
      if (!IsFront(GetActiveElement(iEl))) continue;
      TPave *pave = new TPave(GetActiveElement(iEl)->GetAxis(0)->GetXmin(), 
			      GetActiveElement(iEl)->GetAxis(1)->GetXmin(), 
			      GetActiveElement(iEl)->GetAxis(0)->GetXmax(), 
			      GetActiveElement(iEl)->GetAxis(1)->GetXmax(), 1);
      pave -> SetFillColor(kGreen);
      pave -> Draw("same");
    }

    for (Int_t iEl=0; iEl<GetNReadoutElements(); iEl++) {
      if (!IsFront(GetReadoutElement(iEl))) continue;
      TPave *pave = new TPave(GetReadoutElement(iEl)->GetAxis(0)->GetXmin(), 
			      GetReadoutElement(iEl)->GetAxis(1)->GetXmin(), 
			      GetReadoutElement(iEl)->GetAxis(0)->GetXmax(), 
			      GetReadoutElement(iEl)->GetAxis(1)->GetXmax(), 1);
      pave -> SetFillColor(kRed);
      pave -> Draw("same");
    }

  }
    
  // ------------------- "BACK" option ------------------

  else if (!strcmp(opt, "back")) {

    TCanvas *cnv = new TCanvas("cnv", GetName(), 900, 900);
    cnv->Draw();
    
    TH2D *h = new TH2D("tmp", GetName(), 
		       1, 1.1*GetSupportElement(0)->GetAxis(0)->GetXmin(), 1.1*GetSupportElement(0)->GetAxis(0)->GetXmax(), 
		       1, 1.1*GetSupportElement(0)->GetAxis(1)->GetXmin(), 1.1*GetSupportElement(0)->GetAxis(1)->GetXmax());
    h->SetXTitle("x [cm]");
    h->SetYTitle("y [cm]");
    h->Draw();

    TEllipse *supportExt = new TEllipse(0.0, 0.0, fRMaxSupport, fRMaxSupport);
    TEllipse *supportInt = new TEllipse(0.0, 0.0, fRMinSupport, fRMinSupport);
    supportExt -> SetFillColor(kCyan-10);
    supportExt -> Draw("same");
    supportInt -> Draw("same");

    for (Int_t iEl=0; iEl<GetNActiveElements(); iEl++) {
      if (IsFront(GetActiveElement(iEl))) continue;
      TPave *pave = new TPave(GetActiveElement(iEl)->GetAxis(0)->GetXmin(), 
			      GetActiveElement(iEl)->GetAxis(1)->GetXmin(), 
			      GetActiveElement(iEl)->GetAxis(0)->GetXmax(), 
			      GetActiveElement(iEl)->GetAxis(1)->GetXmax(), 1);
      pave -> SetFillColor(kGreen);
      pave -> Draw("same");
    }

    for (Int_t iEl=0; iEl<GetNReadoutElements(); iEl++) {
      if (IsFront(GetReadoutElement(iEl))) continue;
      TPave *pave = new TPave(GetReadoutElement(iEl)->GetAxis(0)->GetXmin(), 
			      GetReadoutElement(iEl)->GetAxis(1)->GetXmin(), 
			      GetReadoutElement(iEl)->GetAxis(0)->GetXmax(), 
			      GetReadoutElement(iEl)->GetAxis(1)->GetXmax(), 1);
      pave -> SetFillColor(kRed);
      pave -> Draw("same");
    }

  }

  // ------------------- "BOTH" option ------------------

  else if (!strcmp(opt, "both")) {

    TCanvas *cnv = new TCanvas("cnv", GetName(), 900, 900);
    cnv->Draw();

    TH2D *h = new TH2D("tmp", GetName(), 
		       1, 1.1*GetSupportElement(0)->GetAxis(0)->GetXmin(), 1.1*GetSupportElement(0)->GetAxis(0)->GetXmax(), 
		       1, 1.1*GetSupportElement(0)->GetAxis(1)->GetXmin(), 1.1*GetSupportElement(0)->GetAxis(1)->GetXmax());
    h->SetXTitle("x [cm]");
    h->SetYTitle("y [cm]");
    h->Draw();

    TEllipse *supportExt = new TEllipse(0.0, 0.0, fRMaxSupport, fRMaxSupport);
    TEllipse *supportInt = new TEllipse(0.0, 0.0, fRMinSupport, fRMinSupport);
    supportExt -> SetFillColor(kCyan-10);
    supportExt -> Draw("same");
    supportInt -> Draw("same");

    for (Int_t iEl=0; iEl<GetNActiveElements(); iEl++) {
      if (IsFront(GetActiveElement(iEl)) && GetActiveElement(iEl)->GetAxis(0)->GetXmin()<0.) {
	TPave *pave = new TPave(GetActiveElement(iEl)->GetAxis(0)->GetXmin(), 
				GetActiveElement(iEl)->GetAxis(1)->GetXmin(), 
				TMath::Min(GetActiveElement(iEl)->GetAxis(0)->GetXmax(), 0.),
				GetActiveElement(iEl)->GetAxis(1)->GetXmax(), 1);
	pave -> SetFillColor(kGreen);
	pave -> Draw("same");
      }
      else if (!IsFront(GetActiveElement(iEl)) && GetActiveElement(iEl)->GetAxis(0)->GetXmax()>0.) {
	TPave *pave = new TPave(TMath::Max(GetActiveElement(iEl)->GetAxis(0)->GetXmin(), 0.), 
				GetActiveElement(iEl)->GetAxis(1)->GetXmin(), 
				GetActiveElement(iEl)->GetAxis(0)->GetXmax(), 
				GetActiveElement(iEl)->GetAxis(1)->GetXmax(), 1);
	pave -> SetFillColor(kGreen);
	pave -> Draw("same");
      }
    }
    
    for (Int_t iEl=0; iEl<GetNReadoutElements(); iEl++) {
      if (IsFront(GetReadoutElement(iEl)) && GetReadoutElement(iEl)->GetAxis(0)->GetXmin()<0.) {
	TPave *pave = new TPave(GetReadoutElement(iEl)->GetAxis(0)->GetXmin(), 
				GetReadoutElement(iEl)->GetAxis(1)->GetXmin(), 
				TMath::Min(GetReadoutElement(iEl)->GetAxis(0)->GetXmax(), 0.), 
				GetReadoutElement(iEl)->GetAxis(1)->GetXmax(), 1);
	pave -> SetFillColor(kRed);
	pave -> Draw("same");
      }
      else if (!IsFront(GetReadoutElement(iEl)) && GetReadoutElement(iEl)->GetAxis(0)->GetXmax()>0.) {
	TPave *pave = new TPave(TMath::Max(GetReadoutElement(iEl)->GetAxis(0)->GetXmin(), 0.),  
				GetReadoutElement(iEl)->GetAxis(1)->GetXmin(), 
				GetReadoutElement(iEl)->GetAxis(0)->GetXmax(), 
				GetReadoutElement(iEl)->GetAxis(1)->GetXmax(), 1);
	pave -> SetFillColor(kRed);
	pave -> Draw("same");
      }
    }
    
  }

  // ------------------- "PROFILE" option ------------------

  else if (!strcmp(opt, "profile")) {

    TCanvas *cnv = new TCanvas("cnv", GetName(), 300, 900);
    cnv->Draw();

    TH2D *h = new TH2D("tmp", GetName(), 
		       1, fZCenter-0.5, fZCenter+0.5, 
		       1, 1.1*GetSupportElement(0)->GetAxis(1)->GetXmin(), 1.1*GetSupportElement(0)->GetAxis(1)->GetXmax());
    h->SetXTitle("z [cm]");
    h->SetYTitle("y [cm]");
    h->Draw();

    TPave *supportExt = new TPave(GetSupportElement(0)->GetAxis(2)->GetXmin(), -fRMaxSupport, 
				  GetSupportElement(0)->GetAxis(2)->GetXmax(),  fRMaxSupport);
    TPave *supportInt = new TPave(GetSupportElement(0)->GetAxis(2)->GetXmin(), -fRMinSupport, 
				  GetSupportElement(0)->GetAxis(2)->GetXmax(),  fRMinSupport);
    supportExt -> SetFillColor(kCyan-10);
    supportInt -> SetFillColor(kCyan-10);
    supportExt -> SetBorderSize(1);
    supportInt -> SetBorderSize(1);
    supportExt -> Draw("same");
    supportInt -> Draw("same");

    for (Int_t iEl=0; iEl<GetNActiveElements(); iEl++) {
      TPave * pave = 0;
      if (IsFront(GetActiveElement(iEl))) {
	pave = new TPave(GetActiveElement(iEl)->GetAxis(2)->GetXmax() - 
			 5*(GetActiveElement(iEl)->GetAxis(2)->GetXmax()-GetActiveElement(iEl)->GetAxis(2)->GetXmin()), 
			 GetActiveElement(iEl)->GetAxis(1)->GetXmin(), 
			 GetActiveElement(iEl)->GetAxis(2)->GetXmax(), 
			 GetActiveElement(iEl)->GetAxis(1)->GetXmax(), 1);
      }
      else {
	pave = new TPave(GetActiveElement(iEl)->GetAxis(2)->GetXmin(), 
			 GetActiveElement(iEl)->GetAxis(1)->GetXmin(), 
			 GetActiveElement(iEl)->GetAxis(2)->GetXmin() + 
			 5*(GetActiveElement(iEl)->GetAxis(2)->GetXmax()-GetActiveElement(iEl)->GetAxis(2)->GetXmin()), 
			 GetActiveElement(iEl)->GetAxis(1)->GetXmax(), 1);
      }	
      pave -> SetFillColor(kGreen);
      pave -> Draw("same");
    }
    
    for (Int_t iEl=0; iEl<GetNReadoutElements(); iEl++) {
      TPave *pave = 0;
      if (IsFront(GetReadoutElement(iEl))) {
	pave = new TPave(GetReadoutElement(iEl)->GetAxis(2)->GetXmax() - 
			 5*(GetReadoutElement(iEl)->GetAxis(2)->GetXmax()-GetReadoutElement(iEl)->GetAxis(2)->GetXmin()), 
			 GetReadoutElement(iEl)->GetAxis(1)->GetXmin(), 
			 GetReadoutElement(iEl)->GetAxis(2)->GetXmax(), 
			 GetReadoutElement(iEl)->GetAxis(1)->GetXmax(), 1);
      }
      else {
	pave = new TPave(GetReadoutElement(iEl)->GetAxis(2)->GetXmin(), 
			 GetReadoutElement(iEl)->GetAxis(1)->GetXmin(), 
			 GetReadoutElement(iEl)->GetAxis(2)->GetXmin() + 
			 5*(GetReadoutElement(iEl)->GetAxis(2)->GetXmax()-GetReadoutElement(iEl)->GetAxis(2)->GetXmin()), 
			 GetReadoutElement(iEl)->GetAxis(1)->GetXmax(), 1);
      }	
      pave -> SetFillColor(kRed);
      pave -> Draw("same");
    }
    
  }

}

//====================================================================================================================================================

Int_t AliMFTPlane::GetNumberOfChips(Option_t *opt) {

  Int_t nChips = 0;

  if (!strcmp(opt, "front")) {
    for (Int_t iEl=0; iEl<GetNActiveElements(); iEl++) {
      if (!IsFront(GetActiveElement(iEl))) continue;
      Double_t length = GetActiveElement(iEl)->GetAxis(0)->GetXmax() - GetActiveElement(iEl)->GetAxis(0)->GetXmin();
      nChips += Int_t (length/AliMFTConstants::fWidthChip) + 1;
    }
  }

  else if (!strcmp(opt, "back")) {
    for (Int_t iEl=0; iEl<GetNActiveElements(); iEl++) {
      if (IsFront(GetActiveElement(iEl))) continue;
      Double_t length = GetActiveElement(iEl)->GetAxis(0)->GetXmax() - GetActiveElement(iEl)->GetAxis(0)->GetXmin();
      nChips += Int_t (length/AliMFTConstants::fWidthChip) + 1;
    }
  }

  return nChips;

}

//====================================================================================================================================================
