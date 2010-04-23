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

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Transition Radiation Detector version 1 -- slow simulator             //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <TLorentzVector.h>
#include <TMath.h>
#include <TRandom.h>
#include <TVirtualMC.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoPhysicalNode.h>

#include "AliTrackReference.h"
#include "AliMC.h"
#include "AliRun.h"
#include "AliGeomManager.h"

#include "AliTRDgeometry.h"
#include "AliTRDCommonParam.h"
#include "AliTRDsimTR.h"
#include "AliTRDv1.h"

ClassImp(AliTRDv1)
 
//_____________________________________________________________________________
AliTRDv1::AliTRDv1()
  :AliTRD()
  ,fTRon(kTRUE)
  ,fTR(NULL)
  ,fStepSize(0)
  ,fWion(0)
{
  //
  // Default constructor
  //

}

//_____________________________________________________________________________
AliTRDv1::AliTRDv1(const char *name, const char *title) 
  :AliTRD(name,title) 
  ,fTRon(kTRUE)
  ,fTR(NULL)
  ,fStepSize(0.1)
  ,fWion(0)
{
  //
  // Standard constructor for Transition Radiation Detector version 1
  //

  SetBufferSize(128000);

  if      (AliTRDCommonParam::Instance()->IsXenon()) {
    fWion = 23.53; // Ionization energy XeCO2 (85/15)
  }
  else if (AliTRDCommonParam::Instance()->IsArgon()) {
    fWion = 27.21; // Ionization energy ArCO2 (82/18)
  }
  else {
    AliFatal("Wrong gas mixture");
    exit(1);
  }

}

//_____________________________________________________________________________
AliTRDv1::~AliTRDv1()
{
  //
  // AliTRDv1 destructor
  //

  if (fTR) {
    delete fTR;
    fTR     = 0;
  }

}
 
//_____________________________________________________________________________
void AliTRDv1::AddAlignableVolumes() const
{
  //
  // Create entries for alignable volumes associating the symbolic volume
  // name with the corresponding volume path. Needs to be syncronized with
  // eventual changes in the geometry.
  //

  TString volPath;
  TString symName;

  TString vpStr   = "ALIC_1/B077_1/BSEGMO";
  TString vpApp1  = "_1/BTRD";
  TString vpApp2  = "_1";
  TString vpApp3a = "/UTR1_1/UTS1_1/UTI1_1/UT";
  TString vpApp3b = "/UTR2_1/UTS2_1/UTI2_1/UT";
  TString vpApp3c = "/UTR3_1/UTS3_1/UTI3_1/UT";

  TString snStr   = "TRD/sm";
  TString snApp1  = "/st";
  TString snApp2  = "/pl";

  //
  // The super modules
  // The symbolic names are: TRD/sm00
  //                           ...
  //                         TRD/sm17
  //
  for (Int_t isector = 0; isector < AliTRDgeometry::Nsector(); isector++) {

    volPath  = vpStr;
    volPath += isector;
    volPath += vpApp1;
    volPath += isector;
    volPath += vpApp2;

    symName  = snStr;
    symName += Form("%02d",isector);

    gGeoManager->SetAlignableEntry(symName.Data(),volPath.Data());

  }

  //
  // The readout chambers
  // The symbolic names are: TRD/sm00/st0/pl0
  //                           ...
  //                         TRD/sm17/st4/pl5
  //
  AliGeomManager::ELayerID idTRD1 = AliGeomManager::kTRD1;
  Int_t layer, modUID;
  
  for (Int_t isector = 0; isector < AliTRDgeometry::Nsector(); isector++) {

    if (fGeometry->GetSMstatus(isector) == 0) continue;

    for (Int_t istack = 0; istack < AliTRDgeometry::Nstack(); istack++) {
      for (Int_t ilayer = 0; ilayer < AliTRDgeometry::Nlayer(); ilayer++) {

	layer = idTRD1 + ilayer;
	modUID = AliGeomManager::LayerToVolUIDSafe(layer,isector*5+istack);

        Int_t idet = AliTRDgeometry::GetDetectorSec(ilayer,istack);

        volPath  = vpStr;
        volPath += isector;
        volPath += vpApp1;
        volPath += isector;
        volPath += vpApp2;
        switch (isector) {
        case 13:
        case 14:
        case 15:
          if (istack == 2) {
            continue;
	  }
          volPath += vpApp3c;
          break;
        case 11:
        case 12:
          volPath += vpApp3b;
          break;
        default:
          volPath += vpApp3a;
	};
        volPath += Form("%02d",idet);
        volPath += vpApp2;

        symName  = snStr;
        symName += Form("%02d",isector);
        symName += snApp1;
        symName += istack;
        symName += snApp2;
        symName += ilayer;

        TGeoPNEntry *alignableEntry = 
	  gGeoManager->SetAlignableEntry(symName.Data(),volPath.Data(),modUID);

	// Add the tracking to local matrix following the TPC example
	if (alignableEntry) {
	  TGeoHMatrix *globMatrix = alignableEntry->GetGlobalOrig();
	  Double_t sectorAngle = 20.0 * (isector % 18) + 10.0;
	  TGeoHMatrix *t2lMatrix  = new TGeoHMatrix();
	  t2lMatrix->RotateZ(sectorAngle);
	  t2lMatrix->MultiplyLeft(&(globMatrix->Inverse()));
	  alignableEntry->SetMatrix(t2lMatrix);
	}
	else {
	  AliError(Form("Alignable entry %s is not valid!",symName.Data()));
	}

      }
    }
  }

}

//_____________________________________________________________________________
void AliTRDv1::CreateGeometry()
{
  //
  // Create the GEANT geometry for the Transition Radiation Detector - Version 1
  // This version covers the full azimuth. 
  //

  // Check that FRAME is there otherwise we have no place where to put the TRD
  AliModule* frame = gAlice->GetModule("FRAME");
  if (!frame) {
    AliError("TRD needs FRAME to be present\n");
    return;
  }

  // Define the chambers
  AliTRD::CreateGeometry();

}

//_____________________________________________________________________________
void AliTRDv1::CreateMaterials()
{
  //
  // Create materials for the Transition Radiation Detector version 1
  //

  AliTRD::CreateMaterials();

}

//_____________________________________________________________________________
void AliTRDv1::CreateTRhit(Int_t det)
{
  //
  // Creates an electron cluster from a TR photon.
  // The photon is assumed to be created a the end of the radiator. The 
  // distance after which it deposits its energy takes into account the 
  // absorbtion of the entrance window and of the gas mixture in drift
  // volume.
  //

  // Maximum number of TR photons per track
  const Int_t   kNTR         = 50;

  TLorentzVector mom;
  TLorentzVector pos;

  Float_t eTR[kNTR];
  Int_t   nTR;

  // Create TR photons
  gMC->TrackMomentum(mom);
  Float_t pTot = mom.Rho();
  fTR->CreatePhotons(11,pTot,nTR,eTR);
  if (nTR > kNTR) {
    AliFatal(Form("Boundary error: nTR = %d, kNTR = %d",nTR,kNTR));
  }

  // Loop through the TR photons
  for (Int_t iTR = 0; iTR < nTR; iTR++) {

    Float_t energyMeV = eTR[iTR] * 0.001;
    Float_t energyeV  = eTR[iTR] * 1000.0;
    Float_t absLength = 0.0;
    Float_t sigma     = 0.0;

    // Take the absorbtion in the entrance window into account
    Double_t muMy = fTR->GetMuMy(energyMeV);
    sigma         = muMy * fFoilDensity;
    if (sigma > 0.0) {
      absLength = gRandom->Exp(1.0/sigma);
      if (absLength < AliTRDgeometry::MyThick()) {
        continue;
      }
    }
    else {
      continue;
    }

    // The absorbtion cross sections in the drift gas
    // Gas-mixture (Xe/CO2)
    Double_t muNo = 0.0;
    if      (AliTRDCommonParam::Instance()->IsXenon()) {
      muNo = fTR->GetMuXe(energyMeV);
    }
    else if (AliTRDCommonParam::Instance()->IsArgon()) {
      muNo = fTR->GetMuAr(energyMeV);
    }
    Double_t muCO = fTR->GetMuCO(energyMeV);
    sigma = (fGasNobleFraction * muNo + (1.0 - fGasNobleFraction) * muCO) 
          * fGasDensity 
          * fTR->GetTemp();

    // The distance after which the energy of the TR photon
    // is deposited.
    if (sigma > 0.0) {
      absLength = gRandom->Exp(1.0/sigma);
      if (absLength > (AliTRDgeometry::DrThick()
                     + AliTRDgeometry::AmThick())) {
        continue;
      }
    }
    else {
      continue;
    }

    // The position of the absorbtion
    Float_t posHit[3];
    gMC->TrackPosition(pos);
    posHit[0] = pos[0] + mom[0] / pTot * absLength;
    posHit[1] = pos[1] + mom[1] / pTot * absLength;
    posHit[2] = pos[2] + mom[2] / pTot * absLength;

    // Create the charge 
    Int_t q = ((Int_t) (energyeV / fWion));

    // Add the hit to the array. TR photon hits are marked 
    // by negative charge
    AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber()
          ,det
          ,posHit
          ,-q
          ,gMC->TrackTime()*1.0e06
          ,kTRUE);

  }

}

//_____________________________________________________________________________
void AliTRDv1::Init() 
{
  //
  // Initialise Transition Radiation Detector after geometry has been built.
  //

  AliTRD::Init();

  AliDebug(1,"Slow simulator\n");

  // Switch on TR simulation as default
  if (!fTRon) {
    AliInfo("TR simulation off");
  }
  else {
    fTR = new AliTRDsimTR();
  }

  AliDebug(1,"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");

}

//_____________________________________________________________________________
void AliTRDv1::StepManager()
{
  //
  // Slow simulator. Every charged track produces electron cluster as hits 
  // along its path across the drift volume. The step size is fixed in
  // this version of the step manager.
  //
  // Works for Xe/CO2 as well as Ar/CO2
  //

  // PDG code electron
  const Int_t   kPdgElectron = 11;

  Int_t    layer  = 0;
  Int_t    stack  = 0;
  Int_t    sector = 0;
  Int_t    det    = 0;
  Int_t    qTot;

  Float_t  hits[3];
  Double_t eDep;

  Bool_t   drRegion = kFALSE;
  Bool_t   amRegion = kFALSE;

  TString  cIdPath;
  Char_t   cIdSector[3];
           cIdSector[2]  = 0;

  TString  cIdCurrent;
  TString  cIdSensDr = "J";
  TString  cIdSensAm = "K";
  Char_t   cIdChamber[3];
           cIdChamber[2] = 0;

  TLorentzVector pos;
  TLorentzVector mom;

  const Int_t    kNlayer      = AliTRDgeometry::Nlayer();
  const Int_t    kNstack      = AliTRDgeometry::Nstack();
  const Int_t    kNdetsec     = kNlayer * kNstack;

  const Double_t kBig         = 1.0e+12;
  const Float_t  kEkinMinStep = 1.0e-5;  // Minimum energy for the step size adjustment

  // Set the maximum step size to a very large number for all 
  // neutral particles and those outside the driftvolume
  gMC->SetMaxStep(kBig); 

  // If not charged track or already stopped or disappeared, just return.
  if ((!gMC->TrackCharge()) || 
        gMC->IsTrackDisappeared()) {
    return;
  }

  // Inside a sensitive volume?
  cIdCurrent = gMC->CurrentVolName();

  if (cIdSensDr == cIdCurrent[1]) {
    drRegion = kTRUE;
  }
  if (cIdSensAm == cIdCurrent[1]) {
    amRegion = kTRUE;
  }

  if ((!drRegion) && 
      (!amRegion)) {
    return;
  }

  // The hit coordinates and charge
  gMC->TrackPosition(pos);
  hits[0] = pos[0];
  hits[1] = pos[1];
  hits[2] = pos[2];

  // The sector number (0 - 17), according to standard coordinate system
  cIdPath      = gGeoManager->GetPath();
  cIdSector[0] = cIdPath[21];
  cIdSector[1] = cIdPath[22];
  sector = atoi(cIdSector);

  // The plane and chamber number
  cIdChamber[0]   = cIdCurrent[2];
  cIdChamber[1]   = cIdCurrent[3];
  Int_t idChamber = (atoi(cIdChamber) % kNdetsec);
  stack = ((Int_t) idChamber / kNlayer);
  layer = ((Int_t) idChamber % kNlayer);

  // The detector number
  det = fGeometry->GetDetector(layer,stack,sector);

  // 0: InFlight 1:Entering 2:Exiting
  Int_t trkStat = 0;

  // Special hits only in the drift region
  if      ((drRegion) &&
           (gMC->IsTrackEntering())) {

    // Create a track reference at the entrance of each
    // chamber that contains the momentum components of the particle
    gMC->TrackMomentum(mom);
    AddTrackReference(gAlice->GetMCApp()->GetCurrentTrackNumber(), AliTrackReference::kTRD);
    trkStat = 1;

    // Create the hits from TR photons if electron/positron is
    // entering the drift volume
    if ((fTR)   &&
        (fTRon) &&
        (TMath::Abs(gMC->TrackPid()) == kPdgElectron)) {
      CreateTRhit(det);
    }

  }
  else if ((amRegion) && 
           (gMC->IsTrackExiting())) {

    // Create a track reference at the exit of each
    // chamber that contains the momentum components of the particle
    gMC->TrackMomentum(mom);
    AddTrackReference(gAlice->GetMCApp()->GetCurrentTrackNumber(), AliTrackReference::kTRD);
    trkStat = 2;

  }
  
  // Calculate the charge according to GEANT Edep
  // Create a new dEdx hit
  eDep = TMath::Max(gMC->Edep(),0.0) * 1.0e+09;
  qTot = (Int_t) (eDep / fWion);
  if ((qTot) ||
      (trkStat)) {
    AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber()
          ,det
          ,hits
          ,qTot
          ,gMC->TrackTime()*1.0e06
          ,drRegion);
  }

  // Set Maximum Step Size
  // Produce only one hit if Ekin is below cutoff
  if ((gMC->Etot() - gMC->TrackMass()) < kEkinMinStep) {
    return;
  }
  gMC->SetMaxStep(fStepSize);

}
