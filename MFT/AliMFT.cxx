// **************************************************************************
// * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
// *                                                                        *
// * Author: The ALICE Off-line Project.                                    *
// * Contributors are mentioned in the code where appropriate.              *
// *                                                                        *
// * Permission to use, copy, modify and distribute this software and its   *
// * documentation strictly for non-commercial purposes is hereby granted   *
// * without fee, provided that the above copyright notice appears in all   *
// * copies and that both the copyright notice and this permission notice   *
// * appear in the supporting documentation. The authors make no claims     *
// * about the suitability of this software for any purpose. It is          *
// * provided "as is" without express or implied warranty.                  *
// **************************************************************************

//====================================================================================================================================================
//
//      Main class of the ALICE Muon Forward Tracker
//
//      Contact author: antonio.uras@cern.ch
//
//====================================================================================================================================================

#include "AliLog.h"
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

ClassImp(AliMFT) 

//====================================================================================================================================================

AliMFT::AliMFT():
  AliDetector(),
  fVersion(1),
  fNPlanes(0),
  fNSlices(0),
  fSDigitsPerPlane(0),
  fDigitsPerPlane(0),
  fRecPointsPerPlane(0),
  fSideDigits(0),
  fSegmentation(0),
  fNameGeomFile(0),
  fChargeDispersion(0),
  fSingleStepForChargeDispersion(0),
  fNStepForChargeDispersion(0),
  fDensitySiOverSupport(10)
{

  // default constructor

}

//====================================================================================================================================================

AliMFT::AliMFT(const Char_t *name, const Char_t *title):
  AliDetector(name, title),
  fVersion(1),
  fNPlanes(0),
  fNSlices(0),
  fSDigitsPerPlane(0),
  fDigitsPerPlane(0),
  fRecPointsPerPlane(0),
  fSideDigits(0),
  fSegmentation(0),
  fNameGeomFile(0),
  fChargeDispersion(0),
  fSingleStepForChargeDispersion(0),
  fNStepForChargeDispersion(0),
  fDensitySiOverSupport(10)
{

  fNameGeomFile = "AliMFTGeometry.root";

  SetGeometry();

  Init();

}

//====================================================================================================================================================

AliMFT::AliMFT(const Char_t *name, const Char_t *title, Char_t *nameGeomFile):
  AliDetector(name, title),
  fVersion(1),
  fNPlanes(0),
  fNSlices(0),
  fSDigitsPerPlane(0),
  fDigitsPerPlane(0),
  fRecPointsPerPlane(0),
  fSideDigits(0),
  fSegmentation(0),
  fNameGeomFile(0),
  fChargeDispersion(0),
  fSingleStepForChargeDispersion(0),
  fNStepForChargeDispersion(0),
  fDensitySiOverSupport(10)
{

  fNameGeomFile = nameGeomFile;

  SetGeometry();

  Init();

}

//====================================================================================================================================================

AliMFT::~AliMFT() {

  if (fSDigitsPerPlane)   { fSDigitsPerPlane->Delete();    delete fSDigitsPerPlane;   }
  if (fDigitsPerPlane)    { fDigitsPerPlane->Delete();     delete fDigitsPerPlane;    }
  if (fRecPointsPerPlane) { fRecPointsPerPlane->Delete();  delete fRecPointsPerPlane; }

}

//====================================================================================================================================================

void AliMFT::CreateMaterials() {

  // Definition of MFT materials  - to be updated to the most recent values

  AliInfo("Start MFT materials");

  // data from PDG booklet 2002     density [gr/cm^3] rad len [cm] abs len [cm]    
  Float_t   aAir[4]={12,14,16,36} ,   zAir[4]={6,7,8,18} ,   wAir[4]={0.000124,0.755267,0.231781,0.012827} , dAir=0.00120479; Int_t nAir=4;   // Air mixture
  Float_t   aSi = 28.085 ,            zSi   = 14 ,           dSi   =  2.329    ,   radSi   =  21.82/dSi , absSi   = 108.4/dSi  ;              // Silicon

  Int_t   matId  = 0;                        // tmp material id number
  Int_t   unsens = 0, sens=1;                // sensitive or unsensitive medium
  Int_t   itgfld = 3;			     // type of field intergration 0 no field -1 user in guswim 1 Runge Kutta 2 helix 3 const field along z
  Float_t maxfld = 5.; 		             // max field value

  Float_t tmaxfd = -10.0;                    // max deflection angle due to magnetic field in one step
  Float_t stemax =  0.001;                   // max step allowed [cm]
  Float_t deemax = -0.2;                     // maximum fractional energy loss in one step 0<deemax<=1  
  Float_t epsil  =  0.001;                   // tracking precision [cm]   
  Float_t stmin  = -0.001;                   // minimum step due to continuous processes [cm] (negative value: choose it automatically)

  Float_t tmaxfdSi =  0.1;                    // max deflection angle due to magnetic field in one step
  Float_t stemaxSi =  5.0e-4;                 // maximum step allowed [cm]
  Float_t deemaxSi =  0.1;                    // maximum fractional energy loss in one step 0<deemax<=1
  Float_t epsilSi  =  0.5e-4;                 // tracking precision [cm]
  Float_t stminSi  = -0.001;                 // minimum step due to continuous processes [cm] (negative value: choose it automatically)

  Int_t   isxfld = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Integ();   // from CreateMaterials in STRUCT/AliPIPEv3.cxx
  Float_t sxmgmx = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Max();     // from CreateMaterials in STRUCT/AliPIPEv3.cxx
      
  AliMixture(++matId,"Air", aAir,  zAir,   dAir,   nAir,   wAir); 
  AliMedium(kAir,    "Air", matId, unsens, itgfld, maxfld, tmaxfd, stemax, deemax, epsil, stmin);
  
  AliMaterial(++matId, "Si", aSi,   zSi,  dSi,    radSi,  absSi  );  
  AliMedium(kSi,       "Si", matId, sens, isxfld, sxmgmx, tmaxfdSi, stemaxSi, deemaxSi, epsilSi, stminSi);

  AliMaterial(++matId, "Readout", aSi,   zSi,    dSi,    radSi,  absSi  );  
  AliMedium(kReadout,  "Readout", matId, unsens, isxfld, sxmgmx, tmaxfdSi, stemaxSi, deemaxSi, epsilSi, stminSi);

  AliMaterial(++matId, "Support", aSi,   zSi,    dSi/fDensitySiOverSupport, fDensitySiOverSupport*radSi, fDensitySiOverSupport*absSi);  
  AliMedium(kSupport,  "Support", matId, unsens, isxfld,  sxmgmx,    tmaxfdSi, stemaxSi, deemaxSi, epsilSi, stminSi);
    
  AliInfo("End MFT materials");
          
}

//====================================================================================================================================================

void AliMFT::CreateGeometry() {

  // Creates detailed geometry simulation (currently GEANT volumes tree)        

  AliInfo("Start MFT preliminary version building");
  if(!gMC->IsRootGeometrySupported()) return;                
  TGeoVolumeAssembly *vol = CreateVol();
  AliInfo("TGeoVolumeAssembly created!");
  gGeoManager->GetVolume("ALIC")->AddNode(vol,0);
  AliInfo("Stop MFT preliminary version building");

  if (fNStepForChargeDispersion) fSingleStepForChargeDispersion = fChargeDispersion/Double_t(fNStepForChargeDispersion);

} 

//====================================================================================================================================================

void AliMFT::StepManager() {

  // Full Step Manager

  if (!fSegmentation) AliFatal("No segmentation available");    // DO WE HAVE A SEGMENTATION???

  if (!(this->IsActive())) return;
  if (!(gMC->TrackCharge())) return;

  TString planeNumber   = gMC->CurrentVolName();
  TString detElemNumber = gMC->CurrentVolName();
  if (planeNumber.Contains("support")) return;
  if (planeNumber.Contains("readout")) return;
  planeNumber.Remove(0,9);
  detElemNumber.Remove(0,19);
  planeNumber.Remove(2);
  detElemNumber.Remove(3);
  Int_t detElemID = fSegmentation->GetDetElemID(planeNumber.Atoi(), detElemNumber.Atoi());

  if (gMC->IsTrackExiting()) {
    AddTrackReference(gAlice->GetMCApp()->GetCurrentTrackNumber(), AliTrackReference::kMFT);
  }

  static TLorentzVector position, momentum;
  static AliMFTHit hit;
  
  Int_t  status = 0;
  
  // Track status
  if (gMC->IsTrackInside())      status +=  1;
  if (gMC->IsTrackEntering())    status +=  2;
  if (gMC->IsTrackExiting())     status +=  4;
  if (gMC->IsTrackOut())         status +=  8;
  if (gMC->IsTrackDisappeared()) status += 16;
  if (gMC->IsTrackStop())        status += 32;
  if (gMC->IsTrackAlive())       status += 64;

  // ---------- Fill hit structure

  hit.SetDetElemID(detElemID);
  hit.SetPlane(planeNumber.Atoi());
  hit.SetTrack(gAlice->GetMCApp()->GetCurrentTrackNumber());
    
  gMC->TrackPosition(position);
  gMC->TrackMomentum(momentum);

  AliDebug(1, Form("AliMFT::StepManager()->%s Hit #%06d (z=%f) belongs to track %02d\n", 
		   gMC->CurrentVolName(), fNhits, position.Z(), gAlice->GetMCApp()->GetCurrentTrackNumber())); 

  hit.SetPosition(position);
  hit.SetTOF(gMC->TrackTime());
  hit.SetMomentum(momentum);
  hit.SetStatus(status);
  hit.SetEloss(gMC->Edep());
  //  hit.SetShunt(GetIshunt());
//   if (gMC->IsTrackEntering()) {
//     hit.SetStartPosition(position);
//     hit.SetStartTime(gMC->TrackTime());
//     hit.SetStartStatus(status);
//     return; // don't save entering hit.
//   } 

  // Fill hit structure with this new hit.
  new ((*fHits)[fNhits++]) AliMFTHit(hit);

  // Save old position... for next hit.
//   hit.SetStartPosition(position);
//   hit.SetStartTime(gMC->TrackTime());
//   hit.SetStartStatus(status);

  return;

}

//====================================================================================================================================================

TGeoVolumeAssembly* AliMFT::CreateVol() {

  // method to create the MFT Geometry (silicon circular planes)

  if (!fSegmentation) CreateGeometry();
  
  TGeoVolumeAssembly *vol = new TGeoVolumeAssembly("MFT");
  TGeoMedium *silicon = gGeoManager->GetMedium("MFT_Si");
  TGeoMedium *readout = gGeoManager->GetMedium("MFT_Readout");
  TGeoMedium *support = gGeoManager->GetMedium("MFT_Support");
  for (Int_t iPar=0; iPar<8; iPar++) AliInfo(Form("silicon->GetParam(%d) = %f", iPar, silicon->GetParam(iPar)));
  for (Int_t iPar=0; iPar<8; iPar++) AliInfo(Form("readout->GetParam(%d) = %f", iPar, readout->GetParam(iPar)));
  for (Int_t iPar=0; iPar<8; iPar++) AliInfo(Form("support->GetParam(%d) = %f", iPar, support->GetParam(iPar)));

  Double_t origin[3] = {0};

  for (Int_t iPlane=0; iPlane<fNPlanes; iPlane++) {

    AliDebug(1, Form("Creating volumes for MFT plane %02d",iPlane));

    AliMFTPlane *plane = fSegmentation->GetPlane(iPlane);

    // --------- support element(s)
    
    origin[0] =  0.5*(plane->GetSupportElement(0)->GetAxis(0)->GetXmax() + plane->GetSupportElement(0)->GetAxis(0)->GetXmin());
    origin[1] =  0.5*(plane->GetSupportElement(0)->GetAxis(1)->GetXmax() + plane->GetSupportElement(0)->GetAxis(1)->GetXmin());
    origin[2] = -0.5*(plane->GetSupportElement(0)->GetAxis(2)->GetXmax() + plane->GetSupportElement(0)->GetAxis(2)->GetXmin());
    TGeoVolume *supportElem = gGeoManager->MakeTube(Form("MFT_plane%02d_support", iPlane), support, 
						    plane->GetRMinSupport(), 
						    plane->GetRMaxSupport(),
						    0.5*(plane->GetSupportElement(0)->GetAxis(2)->GetXmax() - 
							 plane->GetSupportElement(0)->GetAxis(2)->GetXmin()) );
    vol -> AddNode(supportElem, 0, new TGeoTranslation(origin[0], origin[1], origin[2]));

    AliDebug(1, "support elements created!");
    
    // --------- active elements

    for (Int_t iActive=0; iActive<plane->GetNActiveElements(); iActive++) {
      
      Double_t dx = 0.5*TMath::Abs(plane->GetActiveElement(iActive)->GetAxis(0)->GetXmax() - plane->GetActiveElement(iActive)->GetAxis(0)->GetXmin());
      Double_t dy = 0.5*TMath::Abs(plane->GetActiveElement(iActive)->GetAxis(1)->GetXmax() - plane->GetActiveElement(iActive)->GetAxis(1)->GetXmin());
      Double_t dz = 0.5*TMath::Abs(plane->GetActiveElement(iActive)->GetAxis(2)->GetXmax() - plane->GetActiveElement(iActive)->GetAxis(2)->GetXmin());
      dz /= Double_t(fNSlices);

      origin[0] =  0.5*(plane->GetActiveElement(iActive)->GetAxis(0)->GetXmax() + plane->GetActiveElement(iActive)->GetAxis(0)->GetXmin());
      origin[1] =  0.5*(plane->GetActiveElement(iActive)->GetAxis(1)->GetXmax() + plane->GetActiveElement(iActive)->GetAxis(1)->GetXmin());

      for (Int_t iSlice=0; iSlice<fNSlices; iSlice++) {
	origin[2] = -0.5*(plane->GetActiveElement(iActive)->GetAxis(2)->GetXmin() + 2*dz*(iSlice+1) + plane->GetActiveElement(iActive)->GetAxis(2)->GetXmin() + 2*dz*(iSlice) );
	TGeoVolume *activeElem = gGeoManager->MakeBox(Form("MFT_plane%02d_active%03d_slice%02d", iPlane, iActive, iSlice), silicon, dx, dy, dz);
	vol -> AddNode(activeElem, 0, new TGeoTranslation(origin[0], origin[1], origin[2]));
      }

    }

    AliDebug(1, "active elements created!");

    // --------- readout elements

    for (Int_t iReadout=0; iReadout<plane->GetNReadoutElements(); iReadout++) {
      
      Double_t dx = 0.5*TMath::Abs(plane->GetReadoutElement(iReadout)->GetAxis(0)->GetXmax() - plane->GetReadoutElement(iReadout)->GetAxis(0)->GetXmin());
      Double_t dy = 0.5*TMath::Abs(plane->GetReadoutElement(iReadout)->GetAxis(1)->GetXmax() - plane->GetReadoutElement(iReadout)->GetAxis(1)->GetXmin());
      Double_t dz = 0.5*TMath::Abs(plane->GetReadoutElement(iReadout)->GetAxis(2)->GetXmax() - plane->GetReadoutElement(iReadout)->GetAxis(2)->GetXmin());

      origin[0] =  0.5*(plane->GetReadoutElement(iReadout)->GetAxis(0)->GetXmax() + plane->GetReadoutElement(iReadout)->GetAxis(0)->GetXmin());
      origin[1] =  0.5*(plane->GetReadoutElement(iReadout)->GetAxis(1)->GetXmax() + plane->GetReadoutElement(iReadout)->GetAxis(1)->GetXmin());
      origin[2] = -0.5*(plane->GetReadoutElement(iReadout)->GetAxis(2)->GetXmax() + plane->GetReadoutElement(iReadout)->GetAxis(2)->GetXmin());

      TGeoVolume *readoutElem = gGeoManager->MakeBox(Form("MFT_plane%02d_readout%03d", iPlane, iReadout), readout, dx, dy, dz);
      vol -> AddNode(readoutElem, 0, new TGeoTranslation(origin[0], origin[1], origin[2]));
      
    }

    AliDebug(1, "readout elements created!");

  }

  return vol;

}

//====================================================================================================================================================

void AliMFT::Hits2SDigits(){
  
  // Interface method invoked from AliSimulation to create a list of sdigits corresponding to list of hits. Every hit generates one sdigit.

  AliDebug(1,"Start Hits2SDigits.");
  
  if (!fSegmentation) CreateGeometry();
 
  if (!fLoader->TreeH()) fLoader->LoadHits();

  if (!fLoader->TreeS()) {

    for (Int_t iEvt=0;iEvt<fLoader->GetRunLoader()->GetNumberOfEvents(); iEvt++) {

      fLoader->GetRunLoader()->GetEvent(iEvt);
      fLoader->MakeTree("S");
      MakeBranch("S");
      SetTreeAddress();

      AliDebug(1, Form("Event %03d: fLoader->TreeH()->GetEntries() = %2d", iEvt, Int_t(fLoader->TreeH()->GetEntries())));

      for (Int_t iTrack=0; iTrack<fLoader->TreeH()->GetEntries(); iTrack++) {
	fLoader->TreeH()->GetEntry(iTrack);     
	Hits2SDigitsLocal(Hits(), GetSDigitsList(), iTrack);    // convert these hits to a list of sdigits  
      }
      
      fLoader->TreeS()->Fill();
      fLoader->WriteSDigits("OVERWRITE");
      ResetSDigits();

    }
  }

  fLoader->UnloadHits();
  fLoader->UnloadSDigits();

  AliDebug(1,"Stop Hits2SDigits.");
  
}

//====================================================================================================================================================

void AliMFT::Hits2SDigitsLocal(TClonesArray *hits, const TObjArray *pSDig, Int_t track) {

  //  Add sdigits of these hits to the list
  
  AliDebug(1, "Start Hits2SDigitsLocal");
  
  if (!fSegmentation) CreateGeometry();
  
  TClonesArray *pSDigList[fNMaxPlanes];
  for (Int_t iPlane=0; iPlane<fNMaxPlanes; iPlane++) pSDigList[iPlane] = NULL; 
  for (Int_t iPlane=0; iPlane<fNPlanes;    iPlane++) { 
    pSDigList[iPlane] = (TClonesArray*) (*pSDig)[iPlane];
    AliDebug(1,Form("Entries of pSDigList %3d; plane: %02d,",pSDigList[iPlane]->GetEntries(),iPlane));
    if (!track && pSDigList[iPlane]->GetEntries()!=0) AliErrorClass("Some of sdigits lists is not empty");
  }
  
  for (Int_t iHit=0; iHit<hits->GetEntries(); iHit++) {

    AliMFTHit *hit = (AliMFTHit*) hits->At(iHit);

    AliMFTDigit sDigit;
    sDigit.SetEloss(hit->GetEloss());
    sDigit.SetDetElemID(hit->GetDetElemID());
    sDigit.SetPlane(hit->GetPlane());
    sDigit.AddMCLabel(hit->GetTrack()); 

    Int_t xPixel = -1;
    Int_t yPixel = -1;
    if (fSegmentation->Hit2PixelID(hit->X(), hit->Y(), sDigit.GetDetElemID(), xPixel, yPixel)) {
      sDigit.SetPixID(xPixel, yPixel, 0);
      sDigit.SetPixWidth(fSegmentation->GetPixelSizeX(sDigit.GetDetElemID()), 
			 fSegmentation->GetPixelSizeY(sDigit.GetDetElemID()),
			 fSegmentation->GetPixelSizeZ(sDigit.GetDetElemID()));  
      sDigit.SetPixCenter(fSegmentation->GetPixelCenterX(sDigit.GetDetElemID(), xPixel), 
			  fSegmentation->GetPixelCenterY(sDigit.GetDetElemID(), yPixel),
			  fSegmentation->GetPixelCenterZ(sDigit.GetDetElemID(), 0));  
      new ((*pSDigList[sDigit.GetPlane()])[pSDigList[sDigit.GetPlane()]->GetEntries()]) AliMFTDigit(sDigit);
      AliDebug(1, Form("Created new sdigit (%f, %f, %f) from hit (%f, %f, %f)",
		       sDigit.GetPixelCenterX(), sDigit.GetPixelCenterY(), sDigit.GetPixelCenterZ(), hit->X(), hit->Y(), hit->Z()));
//       AliDebug(1, Form("Created new sdigit from hit: residual is (%f, %f, %f)",
// 		       sDigit.GetPixelCenterX()-hit->X(), sDigit.GetPixelCenterY()-hit->Y(), sDigit.GetPixelCenterZ()-hit->Z()));
    }

    // creating "side hits" to simulate the effect of charge dispersion

    Int_t xPixelNew = -1;
    Int_t yPixelNew = -1;
    Double_t x0 = hit->X();
    Double_t y0 = hit->Y();
    Double_t pi4 = TMath::Pi()/4.;
    for (Int_t iStep=0; iStep<fNStepForChargeDispersion; iStep++) {
      Double_t shift = (iStep+1) * fSingleStepForChargeDispersion;
      for (Int_t iAngle=0; iAngle<8; iAngle++) {
	Double_t shiftX = shift*TMath::Cos(iAngle*pi4);
	Double_t shiftY = shift*TMath::Sin(iAngle*pi4);
	if (fSegmentation->Hit2PixelID(x0+shiftX, y0+shiftY, hit->GetDetElemID(), xPixelNew, yPixelNew)) {
	  Bool_t digitExists = kFALSE;
	  if (xPixelNew==xPixel && yPixelNew==yPixel) digitExists = kTRUE;
	  if (!digitExists) {
	    for (Int_t iSideDigit=0; iSideDigit<fSideDigits->GetEntries(); iSideDigit++) {
	      if (xPixelNew==((AliMFTDigit*) fSideDigits->At(iSideDigit))->GetPixelX() && 
		  yPixelNew==((AliMFTDigit*) fSideDigits->At(iSideDigit))->GetPixelY()) {
		digitExists = kTRUE;
		break;
	      }
	    }
	  }
	  if (!digitExists) {
	    AliMFTDigit newSDigit;
	    newSDigit.SetEloss(0.);
	    newSDigit.SetDetElemID(hit->GetDetElemID());
	    newSDigit.SetPlane(hit->GetPlane());
	    newSDigit.AddMCLabel(hit->GetTrack());
	    newSDigit.SetPixID(xPixelNew, yPixelNew, 0);
	    newSDigit.SetPixWidth(fSegmentation->GetPixelSizeX(sDigit.GetDetElemID()), 
				  fSegmentation->GetPixelSizeY(sDigit.GetDetElemID()),
				  fSegmentation->GetPixelSizeZ(sDigit.GetDetElemID()));  
	    newSDigit.SetPixCenter(fSegmentation->GetPixelCenterX(sDigit.GetDetElemID(), xPixelNew), 
				   fSegmentation->GetPixelCenterY(sDigit.GetDetElemID(), yPixelNew),
				   fSegmentation->GetPixelCenterZ(sDigit.GetDetElemID(), 0)); 
	    new ((*fSideDigits)[fSideDigits->GetEntries()]) AliMFTDigit(newSDigit);
	  }
	}
      }
    }

    for (Int_t iSideDigit=0; iSideDigit<fSideDigits->GetEntries(); iSideDigit++) {
      AliMFTDigit *newDigit = (AliMFTDigit*) fSideDigits->At(iSideDigit);
      new ((*pSDigList[sDigit.GetPlane()])[pSDigList[sDigit.GetPlane()]->GetEntries()]) AliMFTDigit(*newDigit);
    }

    fSideDigits->Clear();    

  }
  
  AliDebug(1,"Stop Hits2SDigitsLocal");

}

//====================================================================================================================================================

void AliMFT::MakeBranch(Option_t *option) {

  printf("AliMFT::MakeBranch(...)\n");

  // Create Tree branches 
  AliDebug(1, Form("Start with option= %s.",option));
  
  const Int_t kBufSize = 4000;
  
  const Char_t *cH = strstr(option,"H");
  const Char_t *cD = strstr(option,"D");
  const Char_t *cS = strstr(option,"S");

  if (cH && fLoader->TreeH()) {
    CreateHits();
    MakeBranchInTree(fLoader->TreeH(), "MFT", &fHits, kBufSize, 0);   
  }

  if (cS && fLoader->TreeS()) {
    CreateSDigits();
    for(Int_t iPlane=0; iPlane<fNPlanes; iPlane++) MakeBranchInTree(fLoader->TreeS(), 
								    Form("Plane_%02d",iPlane), 
								    &((*fSDigitsPerPlane)[iPlane]), 
								    kBufSize, 0);
  }
  
  if (cD && fLoader->TreeD()) {
    CreateDigits();
    for(Int_t iPlane=0; iPlane<fNPlanes; iPlane++) MakeBranchInTree(fLoader->TreeD(), 
								    Form("Plane_%02d",iPlane), 
								    &((*fDigitsPerPlane)[iPlane]),
								    kBufSize, 0);
  }

  AliDebug(1,"Stop.");

}

//====================================================================================================================================================

void AliMFT::SetTreeAddress() {

  AliDebug(1, "AliMFT::SetTreeAddress()");

  //Set branch address for the Hits and Digits Tree.
  AliDebug(1, "Start.");

  AliDebug(1, Form("AliMFT::SetTreeAddress Hits  fLoader->TreeH() = %p\n", fLoader->TreeH()));
  if (fLoader->TreeH() && fLoader->TreeH()->GetBranch("MFT")) {
    CreateHits();
    fLoader->TreeH()->SetBranchAddress("MFT", &fHits);
  }
    
  AliDebug(1, Form("AliMFT::SetTreeAddress SDigits  fLoader->TreeS() = %p\n", fLoader->TreeS()));
  if (fLoader->TreeS() && fLoader->TreeS()->GetBranch("Plane_00")) {
    CreateSDigits();
    for(Int_t iPlane=0; iPlane<fNPlanes; iPlane++) {
      fLoader->TreeS()->SetBranchAddress(Form("Plane_%02d",iPlane), &((*fSDigitsPerPlane)[iPlane]));
    }
  }
    
  AliDebug(1, Form("AliMFT::SetTreeAddress Digits  fLoader->TreeD() = %p\n", fLoader->TreeD()));
  if (fLoader->TreeD() && fLoader->TreeD()->GetBranch("Plane_00")) {
    CreateDigits(); 
    for(Int_t iPlane=0; iPlane<fNPlanes; iPlane++) {
      fLoader->TreeD()->SetBranchAddress(Form("Plane_%02d",iPlane), &((*fDigitsPerPlane)[iPlane]));
    }
  }

  AliDebug(1, Form("AliMFT::SetTreeAddress RecPoints  fLoader->TreeR() = %p\n", fLoader->TreeR()));
  if (fLoader->TreeR() && fLoader->TreeR()->GetBranch("Plane_00")) {
    CreateRecPoints(); 
    for(Int_t iPlane=0; iPlane<fNPlanes; iPlane++) {
      fLoader->TreeR()->SetBranchAddress(Form("Plane_%02d",iPlane), &((*fRecPointsPerPlane)[iPlane]));
    }
  }

  AliDebug(1,"Stop.");

}

//====================================================================================================================================================

void AliMFT::SetGeometry() {

  printf("AliMFT::SetGeometry\n");

  fSegmentation = new AliMFTSegmentation(fNameGeomFile.Data());

  fNPlanes = fSegmentation->GetNPlanes();

}

//====================================================================================================================================================

void AliMFT::CreateHits() { 

  // create array of hits

  AliDebug(1, "AliMFT::CreateHits()");

  if (fHits) return;    
  fHits = new TClonesArray("AliMFTHit");  

}

//====================================================================================================================================================

void AliMFT::CreateSDigits() { 
 
  // create sdigits list

  AliDebug(1, "AliMFT::CreateSDigits()");
 
  if (fSDigitsPerPlane) return; 
  fSDigitsPerPlane = new TObjArray(fNPlanes); 
  for (Int_t iPlane=0; iPlane<fNPlanes; iPlane++) fSDigitsPerPlane->AddAt(new TClonesArray("AliMFTDigit"), iPlane);

  fSideDigits = new TClonesArray("AliMFTDigit");

}

//====================================================================================================================================================

void AliMFT::CreateDigits() {

  // create digits list

  AliDebug(1, "AliMFT::CreateDigits()");

  if (fDigitsPerPlane) return; 
  fDigitsPerPlane = new TObjArray(fNPlanes);
  for(Int_t iPlane=0; iPlane<fNPlanes; iPlane++) fDigitsPerPlane->AddAt(new TClonesArray("AliMFTDigit"), iPlane);

}

//====================================================================================================================================================

void AliMFT::CreateRecPoints() {

  // create recPoints list

  AliDebug(1, "AliMFT::CreateRecPoints()");

  if (fRecPointsPerPlane) return; 
  fRecPointsPerPlane = new TObjArray(fNPlanes);
  for(Int_t iPlane=0; iPlane<fNPlanes; iPlane++) fRecPointsPerPlane->AddAt(new TClonesArray("AliMFTCluster"), iPlane);

}

//====================================================================================================================================================
