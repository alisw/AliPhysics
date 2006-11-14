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

#include <stdlib.h> 

#include <TF1.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TRandom.h>
#include <TVector.h>
#include <TVirtualMC.h>
#include <TGeoManager.h>

#include "AliConst.h"
#include "AliLog.h"
#include "AliMC.h"
#include "AliRun.h"

#include "AliTRDgeometry.h"
#include "AliTRDhit.h"
#include "AliTRDsim.h"
#include "AliTRDv1.h"

ClassImp(AliTRDv1)
 
//_____________________________________________________________________________
AliTRDv1::AliTRDv1()
  :AliTRD()
  ,fTRon(kFALSE)
  ,fTR(NULL)
  ,fTypeOfStepManager(0)
  ,fStepSize(0)
  ,fDeltaE(NULL)
  ,fDeltaG(NULL)
  ,fTrackLength0(0)
  ,fPrimaryTrackPid(0)
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
  ,fTypeOfStepManager(2)
  ,fStepSize(0.1)
  ,fDeltaE(NULL)
  ,fDeltaG(NULL)
  ,fTrackLength0(0)
  ,fPrimaryTrackPid(0)
{
  //
  // Standard constructor for Transition Radiation Detector version 1
  //

  SetBufferSize(128000);

}

//_____________________________________________________________________________
AliTRDv1::~AliTRDv1()
{
  //
  // AliTRDv1 destructor
  //

  if (fDeltaE) {
    delete fDeltaE;
    fDeltaE = 0;
  }

  if (fDeltaG) {
    delete fDeltaG;
    fDeltaG = 0;
  }

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

  TString vpStr  = "ALIC_1/B077_1/BSEGMO";
  TString vpApp1 = "_1/BTRD";
  TString vpApp2 = "_1";
  TString vpApp3 = "/UTR1_1/UTS1_1/UTI1_1/UT";

  TString snStr  = "TRD/sm";
  TString snApp1 = "/st";
  TString snApp2 = "/pl";

  //
  // The super modules
  // The symbolic names are: TRD/sm00
  //                           ...
  //                         TRD/sm17
  //
  for (Int_t isect = 0; isect < AliTRDgeometry::Nsect(); isect++) {

    volPath  = vpStr;
    volPath += isect;
    volPath += vpApp1;
    volPath += isect;
    volPath += vpApp2;

    symName  = snStr;
    symName += Form("%02d",isect);

    gGeoManager->SetAlignableEntry(symName.Data(),volPath.Data());

  }

  //
  // The readout chambers
  // The symbolic names are: TRD/sm00/st0/pl0
  //                           ...
  //                         TRD/sm17/st4/pl5
  //
  for (Int_t isect = 0; isect < AliTRDgeometry::Nsect(); isect++) {
    for (Int_t icham = 0; icham < AliTRDgeometry::Ncham(); icham++) {
      for (Int_t iplan = 0; iplan < AliTRDgeometry::Nplan(); iplan++) {

        Int_t idet = AliTRDgeometry::GetDetectorSec(iplan,icham);

        volPath  = vpStr;
        volPath += isect;
        volPath += vpApp1;
        volPath += isect;
        volPath += vpApp2;
        volPath += vpApp3;
        volPath += Form("%02d",idet);
        volPath += vpApp2;

        symName  = snStr;
        symName += Form("%02d",isect);
        symName += snApp1;
        symName += icham;
        symName += snApp2;
        symName += iplan;

        gGeoManager->SetAlignableEntry(symName.Data(),volPath.Data());

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

  // PDG code electron
  const Int_t   kPdgElectron = 11;

  // Ionization energy
  const Float_t kWion        = 23.53;

  // Maximum number of TR photons per track
  const Int_t   kNTR         = 50;

  TLorentzVector mom;
  TLorentzVector pos;

  // Create TR at the entrance of the chamber
  if (gMC->IsTrackEntering()) {

    // Create TR only for electrons 
    Int_t iPdg = gMC->TrackPid();
    if (TMath::Abs(iPdg) != kPdgElectron) {
      return;
    }

    Float_t eTR[kNTR];
    Int_t   nTR;

    // Create TR photons
    gMC->TrackMomentum(mom);
    Float_t pTot = mom.Rho();
    fTR->CreatePhotons(iPdg,pTot,nTR,eTR);
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
      Double_t muXe = fTR->GetMuXe(energyMeV);
      Double_t muCO = fTR->GetMuCO(energyMeV);
      sigma = (0.85 * muXe + 0.15 * muCO) * fGasDensity * fTR->GetTemp();

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
      Int_t q = ((Int_t) (energyeV / kWion));

      // Add the hit to the array. TR photon hits are marked 
      // by negative charge
      AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber()
            ,det
            ,posHit
            ,-q
            ,kTRUE); 

    }

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
    fTR = new AliTRDsim();
  }

  // First ionization potential (eV) for the gas mixture (90% Xe + 10% CO2)
  const Float_t kPoti = 12.1;
  // Maximum energy (50 keV);
  const Float_t kEend = 50000.0;
  // Ermilova distribution for the delta-ray spectrum
  Float_t poti        = TMath::Log(kPoti);
  Float_t eEnd        = TMath::Log(kEend);

  // Ermilova distribution for the delta-ray spectrum
  fDeltaE = new TF1("deltae" ,Ermilova ,poti,eEnd,0);

  // Geant3 distribution for the delta-ray spectrum
  fDeltaG = new TF1("deltag",IntSpecGeant,2.421257,28.536469,0);

  AliDebug(1,"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");

}

//_____________________________________________________________________________
void AliTRDv1::StepManager()
{
  //
  // Slow simulator. Every charged track produces electron cluster as hits
  // along its path across the drift volume. 
  //

  switch (fTypeOfStepManager) {
   case 0: 
    StepManagerErmilova();
    break;  
   case 1: 
    StepManagerGeant();  
    break;  
   case 2: 
    StepManagerFixedStep();
    break;  
   default: 
    AliWarning("Not a valid Step Manager.");
  }

}

//_____________________________________________________________________________
void AliTRDv1::SelectStepManager(Int_t t)
{
  //
  // Selects a step manager type:
  //   0 - Ermilova
  //   1 - Geant3
  //   2 - Fixed step size
  //

  fTypeOfStepManager = t;
  AliInfo(Form("Step Manager type %d was selected",fTypeOfStepManager));

}

//_____________________________________________________________________________
void AliTRDv1::StepManagerGeant()
{
  //
  // Slow simulator. Every charged track produces electron cluster as hits
  // along its path across the drift volume. The step size is set acording
  // to Bethe-Bloch. The energy distribution of the delta electrons follows
  // a spectrum taken from Geant3.
  //
  // Version by A. Bercuci
  //

  Int_t    pla = 0;
  Int_t    cha = 0;
  Int_t    sec = 0;
  Int_t    det = 0;
  Int_t    iPdg;
  Int_t    qTot;

  Float_t  hits[3];
  Float_t  charge;
  Float_t  aMass;

  Double_t pTot     = 0;
  Double_t eDelta;
  Double_t betaGamma;
  Double_t pp;
  Double_t stepSize = 0;

  Bool_t   drRegion = kFALSE;
  Bool_t   amRegion = kFALSE;

  TString  cIdCurrent;
  TString  cIdSensDr = "J";
  TString  cIdSensAm = "K";
  Char_t   cIdChamber[3];
           cIdChamber[2] = 0;

  TLorentzVector pos;
  TLorentzVector mom;

  TArrayI        processes;

  const Int_t    kNplan       = AliTRDgeometry::Nplan();
  const Int_t    kNcham       = AliTRDgeometry::Ncham();
  const Int_t    kNdetsec     = kNplan * kNcham;

  const Double_t kBig         = 1.0e+12; // Infinitely big
  const Float_t  kWion        = 23.53;   // Ionization energy
  const Float_t  kPTotMaxEl   = 0.002;   // Maximum momentum for e+ e- g

  // Minimum energy for the step size adjustment
  const Float_t  kEkinMinStep = 1.0e-5;
  // energy threshold for production of delta electrons
  const Float_t  kECut        = 1.0e4;
  // Parameters entering the parametrized range for delta electrons
  const Float_t  kRa          = 5.37e-4;
  const Float_t  kRb          = 0.9815;
  const Float_t  kRc          = 3.123e-3;
  // Gas density -> To be made user adjustable !
  // [0.85*0.00549+0.15*0.00186 (Xe-CO2 85-15)]
  const Float_t  kRho         = 0.004945 ; 

  // Plateau value of the energy-loss for electron in xenon
  // The averaged value (26/3/99)
  const Float_t  kPlateau     = 1.55;
  // dN1/dx|min for the gas mixture (90% Xe + 10% CO2)
  const Float_t  kPrim        = 19.34;  
  // First ionization potential (eV) for the gas mixture (90% Xe + 10% CO2)
  const Float_t  kPoti        = 12.1;
  // PDG code electron
  const Int_t    kPdgElectron = 11;  

  // Set the maximum step size to a very large number for all
  // neutral particles and those outside the driftvolume
  gMC->SetMaxStep(kBig);

  // Use only charged tracks
  if (( gMC->TrackCharge()       ) &&
      (!gMC->IsTrackDisappeared())) {

    // Inside a sensitive volume?
    drRegion = kFALSE;
    amRegion = kFALSE;
    cIdCurrent = gMC->CurrentVolName();
    if (cIdSensDr == cIdCurrent[1]) {
      drRegion = kTRUE;
    }
    if (cIdSensAm == cIdCurrent[1]) {
      amRegion = kTRUE;
    }
    if (drRegion || amRegion) {

      // The hit coordinates and charge
      gMC->TrackPosition(pos);
      hits[0] = pos[0];
      hits[1] = pos[1];
      hits[2] = pos[2];

      // The sector number (0 - 17)
      // The numbering goes clockwise and starts at y = 0
      Float_t phi = kRaddeg*TMath::ATan2(pos[0],pos[1]);
      if (phi < 90.0) {
        phi = phi + 270.0;
      }
      else {
        phi = phi -  90.0;
      }
      sec = ((Int_t) (phi / 20.0));

      // The plane and chamber number
      cIdChamber[0]   = cIdCurrent[2];
      cIdChamber[1]   = cIdCurrent[3];
      Int_t idChamber = (atoi(cIdChamber) % kNdetsec);
      cha = kNcham - ((Int_t) idChamber / kNplan) - 1;
      pla = ((Int_t) idChamber % kNplan);

      // Check on selected volumes
      Int_t addthishit = 1;

      // Add this hit
      if (addthishit) {

	// The detector number
        det = fGeometry->GetDetector(pla,cha,sec);

	// Special hits only in the drift region
        if (drRegion) {

          // Create a track reference at the entrance and
          // exit of each chamber that contain the
	  // momentum components of the particle
          if (gMC->IsTrackEntering() || 
              gMC->IsTrackExiting()) {
            gMC->TrackMomentum(mom);
            AddTrackReference(gAlice->GetMCApp()->GetCurrentTrackNumber());
          }

	  if (gMC->IsTrackEntering() && 
              !gMC->IsNewTrack()) {
	    // determine if hit belong to primary track 
	    fPrimaryTrackPid = gAlice->GetMCApp()->GetCurrentTrackNumber();
	    // determine track length when entering the detector
	    fTrackLength0    = gMC->TrackLength();
	  }
					
	  // Create the hits from TR photons
          if (fTR) CreateTRhit(det);

        }

	// Calculate the energy of the delta-electrons
	// modified by Alex Bercuci (A.Bercuci@gsi.de) on 26.01.06
	// take into account correlation with the underlying GEANT tracking
	// mechanism. see
        // http://www-linux.gsi.de/~abercuci/Contributions/TRD/index.html
	//
	// determine the most significant process (last on the processes list)
	// which caused this hit
        gMC->StepProcesses(processes);
        Int_t nofprocesses = processes.GetSize();
        Int_t pid;
	if (!nofprocesses) {
          pid = 0;
	}
	else {
          pid =	processes[nofprocesses-1];		
	}		
		
	// generate Edep according to GEANT parametrisation
	eDelta = TMath::Exp(fDeltaG->GetRandom()) - kPoti;
        eDelta = TMath::Max(eDelta,0.0);
	Float_t prRange = 0.0;
	Float_t range   = gMC->TrackLength() - fTrackLength0;
	// merge GEANT tracker information with locally cooked one
	if (gAlice->GetMCApp()->GetCurrentTrackNumber() == fPrimaryTrackPid) {
	  if      (pid == 27) { 
	    if (eDelta >= kECut) {                
	      prRange = kRa * eDelta * 0.001
                      * (1.0 - kRb / (1.0 + kRc * eDelta * 0.001)) / kRho;
              if (prRange >= (3.7 - range)) {
                eDelta *= 0.1;
	      }
	    }
	  } 
          else if (pid ==  1) {	
	    if (eDelta <  kECut) {
              eDelta *= 0.5;
	    }
	    else {                
	      prRange = kRa * eDelta * 0.001
                      * (1.0 - kRb / (1.0 + kRc * eDelta * 0.001)) / kRho;
              if (prRange >= ((AliTRDgeometry::DrThick()
                             + AliTRDgeometry::AmThick()) - range)) {
                eDelta *= 0.05;
	      }
	      else {
                eDelta *= 0.5;
	      }
	    }
	  } 
          else {
            eDelta = 0.0;
	  }	
	} 
        else {
          eDelta = 0.0;
	}

        // Generate the electron cluster size
        if (eDelta == 0.0) {
          qTot = 0;
	}
	else {
          qTot = ((Int_t) (eDelta / kWion) + 1);
	}

	// Create a new dEdx hit
        AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber()
              ,det
              ,hits
              ,qTot
              ,drRegion);
				
        // Calculate the maximum step size for the next tracking step
	// Produce only one hit if Ekin is below cutoff
        aMass = gMC->TrackMass();
        if ((gMC->Etot() - aMass) > kEkinMinStep) {

          // The energy loss according to Bethe Bloch
          iPdg = TMath::Abs(gMC->TrackPid());
          if ((iPdg != kPdgElectron) ||
	      ((iPdg == kPdgElectron) && 
               (pTot  < kPTotMaxEl))) {
            gMC->TrackMomentum(mom);
            pTot      = mom.Rho();
            betaGamma = pTot / aMass;
            pp        = BetheBlochGeant(betaGamma);
	    // Take charge > 1 into account
            charge     = gMC->TrackCharge();
            if (TMath::Abs(charge) > 1) {
              pp = pp * charge*charge;
	    }
          } 
          else { 
            // Electrons above 20 Mev/c are at the plateau
	    pp = kPrim * kPlateau;
          }

	  Int_t nsteps = 0;
	  do {
            nsteps = gRandom->Poisson(pp);
          } while(!nsteps);
          stepSize = 1.0 / nsteps;
	  gMC->SetMaxStep(stepSize);

	}

      }

    }

  }

}

//_____________________________________________________________________________
void AliTRDv1::StepManagerErmilova()
{
  //
  // Slow simulator. Every charged track produces electron cluster as hits 
  // along its path across the drift volume. The step size is set acording
  // to Bethe-Bloch. The energy distribution of the delta electrons follows
  // a spectrum taken from Ermilova et al.
  //

  Int_t    pla = 0;
  Int_t    cha = 0;
  Int_t    sec = 0;
  Int_t    det = 0;
  Int_t    iPdg;
  Int_t    qTot;

  Float_t  hits[3];
  Double_t random[1];
  Float_t  charge;
  Float_t  aMass;

  Double_t pTot     = 0.0;
  Double_t eDelta;
  Double_t betaGamma;
  Double_t pp;
  Double_t stepSize;

  Bool_t   drRegion = kFALSE;
  Bool_t   amRegion = kFALSE;

  TString  cIdCurrent;
  TString  cIdSensDr = "J";
  TString  cIdSensAm = "K";
  Char_t   cIdChamber[3];
           cIdChamber[2] = 0;

  TLorentzVector pos;
  TLorentzVector mom;

  const Int_t    kNplan       = AliTRDgeometry::Nplan();
  const Int_t    kNcham       = AliTRDgeometry::Ncham();
  const Int_t    kNdetsec     = kNplan * kNcham;

  const Double_t kBig         = 1.0e+12; // Infinitely big
  const Float_t  kWion        = 23.53;   // Ionization energy
  const Float_t  kPTotMaxEl   = 0.002;   // Maximum momentum for e+ e- g 

  // Minimum energy for the step size adjustment
  const Float_t  kEkinMinStep = 1.0e-5;

  // Plateau value of the energy-loss for electron in xenon
  // The averaged value (26/3/99)
  const Float_t  kPlateau     = 1.55;
  // dN1/dx|min for the gas mixture (90% Xe + 10% CO2)
  const Float_t  kPrim        = 48.0;  
  // First ionization potential (eV) for the gas mixture (90% Xe + 10% CO2)
  const Float_t  kPoti        = 12.1;
  // PDG code electron
  const Int_t    kPdgElectron = 11;  

  // Set the maximum step size to a very large number for all 
  // neutral particles and those outside the driftvolume
  gMC->SetMaxStep(kBig); 

  // Use only charged tracks 
  if (( gMC->TrackCharge()       ) &&
      (!gMC->IsTrackDisappeared())) {

    // Inside a sensitive volume?
    drRegion = kFALSE;
    amRegion = kFALSE;
    cIdCurrent = gMC->CurrentVolName();
    if (cIdSensDr == cIdCurrent[1]) {
      drRegion = kTRUE;
    }
    if (cIdSensAm == cIdCurrent[1]) {
      amRegion = kTRUE;
    }
    if (drRegion || amRegion) {

      // The hit coordinates and charge
      gMC->TrackPosition(pos);
      hits[0] = pos[0];
      hits[1] = pos[1];
      hits[2] = pos[2];

      // The sector number (0 - 17)
      // The numbering goes clockwise and starts at y = 0
      Float_t phi = kRaddeg*TMath::ATan2(pos[0],pos[1]);
      if (phi < 90.0) { 
        phi = phi + 270.0;
      }
      else {
        phi = phi -  90.0;
      }
      sec = ((Int_t) (phi / 20.0));

      // The plane and chamber number
      cIdChamber[0] = cIdCurrent[2];
      cIdChamber[1] = cIdCurrent[3];
      Int_t idChamber = (atoi(cIdChamber) % kNdetsec);
      cha = kNcham - ((Int_t) idChamber / kNplan) - 1;
      pla = ((Int_t) idChamber % kNplan);

      // Check on selected volumes
      Int_t addthishit = 1;

      // Add this hit
      if (addthishit) {

	// The detector number
        det = fGeometry->GetDetector(pla,cha,sec);

	// Special hits only in the drift region
        if (drRegion) {

          // Create a track reference at the entrance and
          // exit of each chamber that contain the 
	  // momentum components of the particle
          if (gMC->IsTrackEntering() || 
              gMC->IsTrackExiting()) {
            gMC->TrackMomentum(mom);
            AddTrackReference(gAlice->GetMCApp()->GetCurrentTrackNumber());
          }
          // Create the hits from TR photons
          if (fTR) {
            CreateTRhit(det);
	  }

	}

        // Calculate the energy of the delta-electrons
        eDelta = TMath::Exp(fDeltaE->GetRandom()) - kPoti;
        eDelta = TMath::Max(eDelta,0.0);
        // Generate the electron cluster size
        if (eDelta == 0.0) {
          qTot = 0;
	}
	else {
          qTot = ((Int_t) (eDelta / kWion) + 1);
	}

	// Create a new dEdx hit
        if (drRegion) {
          AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber()
                ,det
                ,hits
                ,qTot
                ,kTRUE);
	}
        else {
          AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber()
                ,det
                ,hits
                ,qTot
                ,kFALSE);
	}

        // Calculate the maximum step size for the next tracking step
	// Produce only one hit if Ekin is below cutoff 
        aMass = gMC->TrackMass();
        if ((gMC->Etot() - aMass) > kEkinMinStep) {

          // The energy loss according to Bethe Bloch
          iPdg  = TMath::Abs(gMC->TrackPid());
          if ((iPdg != kPdgElectron) ||
	      ((iPdg == kPdgElectron) && 
               (pTot  < kPTotMaxEl))) {
            gMC->TrackMomentum(mom);
            pTot      = mom.Rho();
            betaGamma = pTot / aMass;
            pp        = kPrim * BetheBloch(betaGamma);
	    // Take charge > 1 into account
            charge = gMC->TrackCharge();
            if (TMath::Abs(charge) > 1) {
              pp = pp * charge*charge;
	    }
          } 
          else { 
            // Electrons above 20 Mev/c are at the plateau
	    pp = kPrim * kPlateau;
          }
      
          if (pp > 0.0) {
            do {
              gMC->GetRandom()->RndmArray(1,random);
	    }
            while ((random[0] == 1.0) || 
                   (random[0] == 0.0));
            stepSize = - TMath::Log(random[0]) / pp; 
            gMC->SetMaxStep(stepSize);
	  }

	}

      }

    }

  }

}

//_____________________________________________________________________________
void AliTRDv1::StepManagerFixedStep()
{
  //
  // Slow simulator. Every charged track produces electron cluster as hits 
  // along its path across the drift volume. The step size is fixed in
  // this version of the step manager.
  //

  Int_t    pla = 0;
  Int_t    cha = 0;
  Int_t    sec = 0;
  Int_t    det = 0;
  Int_t    qTot;

  Float_t  hits[3];
  Double_t eDep;

  Bool_t   drRegion = kFALSE;
  Bool_t   amRegion = kFALSE;

  TString  cIdCurrent;
  TString  cIdSensDr = "J";
  TString  cIdSensAm = "K";
  Char_t   cIdChamber[3];
  cIdChamber[2] = 0;

  TLorentzVector pos;
  TLorentzVector mom;

  const Int_t    kNplan       = AliTRDgeometry::Nplan();
  const Int_t    kNcham       = AliTRDgeometry::Ncham();
  const Int_t    kNdetsec     = kNplan * kNcham;

  const Double_t kBig         = 1.0e+12;

  const Float_t  kWion        = 23.53;   // Ionization energy
  const Float_t  kEkinMinStep = 1.0e-5;  // Minimum energy for the step size adjustment

  // Set the maximum step size to a very large number for all 
  // neutral particles and those outside the driftvolume
  gMC->SetMaxStep(kBig); 

  // If not charged track or already stopped or disappeared, just return.
  if ((!gMC->TrackCharge()) || 
        gMC->IsTrackDisappeared()) return;

  // Inside a sensitive volume?
  cIdCurrent = gMC->CurrentVolName();

  if (cIdSensDr == cIdCurrent[1]) drRegion = kTRUE;
  if (cIdSensAm == cIdCurrent[1]) amRegion = kTRUE;

  if ((!drRegion) && 
      (!amRegion)) {
    return;
  }

  // The hit coordinates and charge
  gMC->TrackPosition(pos);
  hits[0] = pos[0];
  hits[1] = pos[1];
  hits[2] = pos[2];

  // The sector number (0 - 17)
  // The numbering goes clockwise and starts at y = 0
  Float_t phi = kRaddeg*TMath::ATan2(pos[0],pos[1]);
  if (phi < 90.0) {
    phi = phi + 270.0;
  }
  else {          
    phi = phi -  90.0;
  }
  sec = ((Int_t) (phi / 20.0));

  // The plane and chamber number
  cIdChamber[0]   = cIdCurrent[2];
  cIdChamber[1]   = cIdCurrent[3];
  Int_t idChamber = (atoi(cIdChamber) % kNdetsec);
  cha = kNcham - ((Int_t) idChamber / kNplan) - 1;
  pla = ((Int_t) idChamber % kNplan);

  // Check on selected volumes
  Int_t addthishit = 1;

  if (!addthishit) {
    return;
  }

  // The detector number
  det = fGeometry->GetDetector(pla,cha,sec);

  // 0: InFlight 1:Entering 2:Exiting
  Int_t trkStat = 0;

  // Special hits only in the drift region
  if (drRegion) {

    // Create a track reference at the entrance and exit of each
    // chamber that contain the momentum components of the particle

    if (gMC->IsTrackEntering()) {
      gMC->TrackMomentum(mom);
      AddTrackReference(gAlice->GetMCApp()->GetCurrentTrackNumber());
      trkStat = 1;
    }
    if (gMC->IsTrackExiting()) {
      gMC->TrackMomentum(mom);
      AddTrackReference(gAlice->GetMCApp()->GetCurrentTrackNumber());
      trkStat = 2;
    }

    // Create the hits from TR photons
    if (fTR) {
      CreateTRhit(det);    
    }

  }
  
  // Calculate the charge according to GEANT Edep
  // Create a new dEdx hit
  eDep = TMath::Max(gMC->Edep(),0.0) * 1.0e+09;
  qTot = (Int_t) (eDep / kWion);
  AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber()
        ,det
        ,hits
        ,qTot
        ,drRegion);

  // Set Maximum Step Size
  // Produce only one hit if Ekin is below cutoff
  if ((gMC->Etot() - gMC->TrackMass()) < kEkinMinStep) {
    return;
  }
  gMC->SetMaxStep(fStepSize);

}

//_____________________________________________________________________________
Double_t AliTRDv1::BetheBloch(Double_t bg) 
{
  //
  // Parametrization of the Bethe-Bloch-curve
  // The parametrization is the same as for the TPC and is taken from Lehrhaus.
  //

  // This parameters have been adjusted to averaged values from GEANT
  const Double_t kP1    = 7.17960e-02;
  const Double_t kP2    = 8.54196;
  const Double_t kP3    = 1.38065e-06;
  const Double_t kP4    = 5.30972;
  const Double_t kP5    = 2.83798;

  // Lower cutoff of the Bethe-Bloch-curve to limit step sizes
  const Double_t kBgMin = 0.8;
  const Double_t kBBMax = 6.83298;

  if (bg > kBgMin) {
    Double_t yy = bg / TMath::Sqrt(1.0 + bg*bg);
    Double_t aa = TMath::Power(yy,kP4);
    Double_t bb = TMath::Power((1.0/bg),kP5);
             bb = TMath::Log(kP3 + bb);
    return ((kP2 - aa - bb) * kP1 / aa);
  }
  else {
    return kBBMax;
  }

}

//_____________________________________________________________________________
Double_t AliTRDv1::BetheBlochGeant(Double_t bg)
{
  //
  // Return dN/dx (number of primary collisions per centimeter)
  // for given beta*gamma factor.
  //
  // Implemented by K.Oyama according to GEANT 3 parametrization shown in
  // A.Andronic's webpage: http://www-alice.gsi.de/trd/papers/dedx/dedx.html
  // This must be used as a set with IntSpecGeant.
  //

  Int_t    i = 0;

  Double_t arrG[20]  = {     1.100000,     1.200000,     1.300000,     1.500000
                       ,     1.800000,     2.000000,     2.500000,     3.000000
                       ,     4.000000,     7.000000,    10.000000,    20.000000
                       ,    40.000000,    70.000000,   100.000000,   300.000000
                       ,   600.000000,  1000.000000,  3000.000000, 10000.000000 };

  Double_t arrNC[20] = {    75.009056,    45.508083,    35.299252,    27.116327
                       ,    22.734999,    21.411915,    19.934095,    19.449375
                       ,    19.344431,    20.185553,    21.027925,    22.912676
                       ,    24.933352,    26.504053,    27.387468,    29.566597
                       ,    30.353779,    30.787134,    31.129285,    31.157350 };

  // Betagamma to gamma
  Double_t g = TMath::Sqrt(1.0 + bg*bg);

  // Find the index just before the point we need.
  for (i = 0; i < 18; i++) {
    if ((arrG[i]   < g) && 
        (arrG[i+1] > g)) {
      break;
    }
  }

  // Simple interpolation.
  Double_t pp = ((arrNC[i+1] - arrNC[i]) / (arrG[i+1]  - arrG[i])) 
              * (g - arrG[i]) + arrNC[i];

  return pp;

}

//_____________________________________________________________________________
Double_t Ermilova(Double_t *x, Double_t *)
{
  //
  // Calculates the delta-ray energy distribution according to Ermilova.
  // Logarithmic scale !
  //

  Double_t energy;
  Double_t dpos;
  Double_t dnde;

  Int_t    pos1;
  Int_t    pos2;

  const Int_t kNv = 31;

  Float_t vxe[kNv] = {  2.3026,  2.9957,  3.4012,  3.6889,  3.9120  
                     ,  4.0943,  4.2485,  4.3820,  4.4998,  4.6052
                     ,  4.7005,  5.0752,  5.2983,  5.7038,  5.9915
                     ,  6.2146,  6.5221,  6.9078,  7.3132,  7.6009
                     ,  8.0064,  8.5172,  8.6995,  8.9872,  9.2103
                     ,  9.4727,  9.9035, 10.3735, 10.5966, 10.8198
                     , 11.5129 };

  Float_t vye[kNv] = { 80.0,    31.0,    23.3,    21.1,    21.0
                     , 20.9,    20.8,    20.0,    16.0,    11.0
                     ,  8.0,     6.0,     5.2,     4.6,     4.0
                     ,  3.5,     3.0,     1.4,     0.67,    0.44
                     ,  0.3,     0.18,    0.12,    0.08,    0.056
                     ,  0.04,    0.023,   0.015,   0.011,    0.01
		     ,  0.004  };

  energy = x[0];

  // Find the position 
  pos1 = 0;
  pos2 = 0;
  dpos = 0;
  do {
    dpos = energy - vxe[pos2++];
  } 
  while (dpos > 0);
  pos2--; 
  if (pos2 > kNv) {
    pos2 = kNv - 1;
  }
  pos1 = pos2 - 1;

  // Differentiate between the sampling points
  dnde = (vye[pos1] - vye[pos2]) / (vxe[pos2] - vxe[pos1]);

  return dnde;

}

//_____________________________________________________________________________
Double_t IntSpecGeant(Double_t *x, Double_t *)
{
  //
  // Integrated spectrum from Geant3
  //

  const Int_t npts = 83;
  Double_t arre[npts]    = {  2.421257,     2.483278,    2.534301,     2.592230
                           ,  2.672067,     2.813299,    3.015059,     3.216819
                           ,  3.418579,     3.620338,    3.868209,     3.920198
                           ,  3.978284,     4.063923,    4.186264,     4.308605
                           ,  4.430946,     4.553288,    4.724261,     4.837736
                           ,  4.999842,     5.161949,    5.324056,     5.486163
                           ,  5.679688,     5.752998,    5.857728,     5.962457
                           ,  6.067185,     6.171914,    6.315653,     6.393674
                           ,  6.471694,     6.539689,    6.597658,     6.655627
                           ,  6.710957,     6.763648,    6.816338,     6.876198
                           ,  6.943227,     7.010257,    7.106285,     7.252151
                           ,  7.460531,     7.668911,    7.877290,     8.085670
                           ,  8.302979,     8.353585,    8.413120,     8.483500
                           ,  8.541030,     8.592857,    8.668865,     8.820485
                           ,  9.037086,     9.253686,    9.470286,     9.686887
                           ,  9.930838,     9.994655,   10.085822,    10.176990
                           , 10.268158,    10.359325,   10.503614,    10.627565
                           , 10.804637,    10.981709,   11.158781,    11.335854
                           , 11.593397,    11.781165,   12.049404,    12.317644
                           , 12.585884,    12.854123,   14.278421,    16.975889
                           , 20.829416,    24.682943,   28.536469 };

  Double_t arrdnde[npts] = { 10.960000,    10.960000,   10.359500,     9.811340
                           ,  9.1601500,    8.206670,    6.919630,     5.655430
                           ,  4.6221300,    3.777610,    3.019560,     2.591950
                           ,  2.5414600,    2.712920,    3.327460,     4.928240
                           ,  7.6185300,   10.966700,   12.225800,     8.094750
                           ,  3.3586900,    1.553650,    1.209600,     1.263840
                           ,  1.3241100,    1.312140,    1.255130,     1.165770
                           ,  1.0594500,    0.945450,    0.813231,     0.699837
                           ,  0.6235580,    2.260990,    2.968350,     2.240320
                           ,  1.7988300,    1.553300,    1.432070,     1.535520
                           ,  1.4429900,    1.247990,    1.050750,     0.829549
                           ,  0.5900280,    0.395897,    0.268741,     0.185320
                           ,  0.1292120,    0.103545,    0.0949525,    0.101535
                           ,  0.1276380,    0.134216,    0.123816,     0.104557
                           ,  0.0751843,    0.0521745,   0.0373546,    0.0275391
                           ,  0.0204713,    0.0169234,   0.0154552,    0.0139194
                           ,  0.0125592,    0.0113638,   0.0107354,    0.0102137
                           ,  0.00845984,   0.00683338,  0.00556836,   0.00456874
                           ,  0.0036227,    0.00285991,  0.00226664,   0.00172234
                           ,  0.00131226,   0.00100284,  0.000465492,  7.26607e-05
                           ,  3.63304e-06,  0.0000000,   0.0000000   };

  Int_t    i;
  Double_t energy = x[0];

  for (i = 0; i < npts; i++) {
    if (energy < arre[i]) {
      break;
    }
  }

  if (i == 0) {
    AliErrorGeneral("AliTRDv1::IntSpecGeant","Given energy value is too small or zero");
  }

  return arrdnde[i];

}
