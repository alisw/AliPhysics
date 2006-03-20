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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Transition Radiation Detector version 1 -- slow simulator                //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliTRDfullClass.gif">
*/
//End_Html
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <stdlib.h> 

#include <TF1.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TRandom.h>
#include <TVector.h>
#include <TVirtualMC.h>

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
AliTRDv1::AliTRDv1():AliTRD()
{
  //
  // Default constructor
  //

  fSensSelect        =  0;
  fSensPlane         = -1;
  fSensChamber       = -1;
  fSensSector        = -1;
  fSensSectorRange   =  0;

  fDeltaE            = NULL;
  fDeltaG            = NULL;
  fTR                = NULL;

  fStepSize          = 0.1;
  fTypeOfStepManager = 1;

}

//_____________________________________________________________________________
AliTRDv1::AliTRDv1(const char *name, const char *title) 
         :AliTRD(name, title) 
{
  //
  // Standard constructor for Transition Radiation Detector version 1
  //

  fSensSelect        =  0;
  fSensPlane         = -1;
  fSensChamber       = -1;
  fSensSector        = -1;
  fSensSectorRange   =  0;

  fDeltaE            = NULL;
  fDeltaG            = NULL;
  fTR                = NULL;
  fStepSize          = 0.1;
  fTypeOfStepManager = 1;

  SetBufferSize(128000);

}

//_____________________________________________________________________________
AliTRDv1::AliTRDv1(const AliTRDv1 &trd):AliTRD(trd)
{
  //
  // Copy constructor
  //

  ((AliTRDv1 &) trd).Copy(*this);

}

//_____________________________________________________________________________
AliTRDv1::~AliTRDv1()
{
  //
  // AliTRDv1 destructor
  //

  if (fDeltaE) delete fDeltaE;
  if (fDeltaG) delete fDeltaG;
  if (fTR)     delete fTR;

}
 
//_____________________________________________________________________________
AliTRDv1 &AliTRDv1::operator=(const AliTRDv1 &trd)
{
  //
  // Assignment operator
  //

  if (this != &trd) ((AliTRDv1 &) trd).Copy(*this);
  return *this;

}
 
//_____________________________________________________________________________
void AliTRDv1::Copy(TObject &trd) const
{
	printf("void AliTRDv1::Copy(TObject &trd) const\n");
  //
  // Copy function
  //

  ((AliTRDv1 &) trd).fSensSelect        = fSensSelect;
  ((AliTRDv1 &) trd).fSensPlane         = fSensPlane;
  ((AliTRDv1 &) trd).fSensChamber       = fSensChamber;
  ((AliTRDv1 &) trd).fSensSector        = fSensSector;
  ((AliTRDv1 &) trd).fSensSectorRange   = fSensSectorRange;

  ((AliTRDv1 &) trd).fTypeOfStepManager = fTypeOfStepManager;
  ((AliTRDv1 &) trd).fStepSize          = fStepSize;

  fDeltaE->Copy(*((AliTRDv1 &) trd).fDeltaE);
  fDeltaG->Copy(*((AliTRDv1 &) trd).fDeltaG);
  fTR->Copy(*((AliTRDv1 &) trd).fTR);

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
  if (!frame) return;

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

  TLorentzVector mom, pos;

  // Create TR at the entrance of the chamber
  if (gMC->IsTrackEntering()) {

    // Create TR only for electrons 
    Int_t iPdg = gMC->TrackPid();
    if (TMath::Abs(iPdg) != kPdgElectron) return;

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
      Float_t absLength = 0;
      Float_t sigma     = 0;

      // Take the absorbtion in the entrance window into account
      Double_t muMy = fTR->GetMuMy(energyMeV);
      sigma = muMy * fFoilDensity;
      if (sigma > 0.0) {
        absLength = gRandom->Exp(1.0/sigma);
        if (absLength < AliTRDgeometry::MyThick()) continue;
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
     AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber(),det,posHit,-q,kTRUE); 

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
  if (fSensSelect) {
    if (fSensPlane   >= 0)
      AliInfo(Form("Only plane %d is sensitive"));
    if (fSensChamber >= 0)   
      AliInfo(Form("Only chamber %d is sensitive",fSensChamber));
    if (fSensSector  >= 0) {
      Int_t sens1  = fSensSector;
      Int_t sens2  = fSensSector + fSensSectorRange;
            sens2 -= ((Int_t) (sens2 / AliTRDgeometry::Nsect())) 
                   * AliTRDgeometry::Nsect();
	    AliInfo(Form("Only sectors %d - %d are sensitive\n",sens1,sens2-1));
    }
  }
  if (fTR) 
    AliInfo("TR simulation on")
  else
    AliInfo("TR simulation off");

  // First ionization potential (eV) for the gas mixture (90% Xe + 10% CO2)
  const Float_t kPoti = 12.1;
  // Maximum energy (50 keV);
  const Float_t kEend = 50000.0;
  // Ermilova distribution for the delta-ray spectrum
  Float_t poti = TMath::Log(kPoti);
  Float_t eEnd = TMath::Log(kEend);

  // Ermilova distribution for the delta-ray spectrum
  fDeltaE = new TF1("deltae" ,Ermilova ,poti,eEnd,0);

  // Geant3 distribution for the delta-ray spectrum
  fDeltaG = new TF1("deltag",IntSpecGeant,2.421257,28.536469,0);

  AliDebug(1,"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");

}

//_____________________________________________________________________________
AliTRDsim *AliTRDv1::CreateTR()
{
  //
  // Enables the simulation of TR
  //

  fTR = new AliTRDsim();
  return fTR;

}

//_____________________________________________________________________________
void AliTRDv1::SetSensPlane(Int_t iplane)
{
  //
  // Defines the hit-sensitive plane (0-5)
  //

  if ((iplane < 0) || (iplane > 5)) {
    AliWarning(Form("Wrong input value:%d",iplane));
    AliWarning("Use standard setting");
    fSensPlane  = -1;
    fSensSelect =  0;
    return;
  }

  fSensSelect = 1;
  fSensPlane  = iplane;

}

//_____________________________________________________________________________
void AliTRDv1::SetSensChamber(Int_t ichamber)
{
  //
  // Defines the hit-sensitive chamber (0-4)
  //

  if ((ichamber < 0) || (ichamber > 4)) {
    AliWarning(Form("Wrong input value: %d",ichamber));
    AliWarning("Use standard setting");
    fSensChamber = -1;
    fSensSelect  =  0;
    return;
  }

  fSensSelect  = 1;
  fSensChamber = ichamber;

}

//_____________________________________________________________________________
void AliTRDv1::SetSensSector(Int_t isector)
{
  //
  // Defines the hit-sensitive sector (0-17)
  //

  SetSensSector(isector,1);

}

//_____________________________________________________________________________
void AliTRDv1::SetSensSector(Int_t isector, Int_t nsector)
{
  //
  // Defines a range of hit-sensitive sectors. The range is defined by
  // <isector> (0-17) as the starting point and <nsector> as the number 
  // of sectors to be included.
  //

  if ((isector < 0) || (isector > 17)) {
    AliWarning(Form("Wrong input value <isector>: %d",isector));
    AliWarning("Use standard setting");
    fSensSector      = -1;
    fSensSectorRange =  0;
    fSensSelect      =  0;
    return;
  }

  if ((nsector < 1) || (nsector > 18)) {
    AliWarning(Form("Wrong input value <nsector>: %d",nsector));
    AliWarning("Use standard setting");
    fSensSector      = -1;
    fSensSectorRange =  0;
    fSensSelect      =  0;
    return;
  }

  fSensSelect      = 1;
  fSensSector      = isector;
  fSensSectorRange = nsector;

}

//_____________________________________________________________________________
void AliTRDv1::StepManager()
{
  //
  // Slow simulator. Every charged track produces electron cluster as hits
  // along its path across the drift volume. 
  //

  switch (fTypeOfStepManager) {
    case 0  : StepManagerErmilova();  break;  // 0 is Ermilova
    case 1  : StepManagerGeant();     break;  // 1 is Geant
    case 2  : StepManagerFixedStep(); break;  // 2 is fixed step
    default : AliWarning("Not a valid Step Manager.");
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

/*  if (t == 1) {
    AliWarning("Sorry, Geant parametrization step manager is not implemented yet. Please ask K.Oyama for detail.");
  }
*/
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
  Int_t    pla = 0;
  Int_t    cha = 0;
  Int_t    sec = 0;
  Int_t    det = 0;
  Int_t    iPdg;
  Int_t    qTot;

  Float_t  hits[3];
  Float_t  charge;
  Float_t  aMass;

  Double_t pTot = 0;
  Double_t eDelta;
  Double_t betaGamma, pp;
  Double_t stepSize=0;

  Bool_t   drRegion = kFALSE;
  Bool_t   amRegion = kFALSE;

  TString  cIdCurrent;
  TString  cIdSensDr = "J";
  TString  cIdSensAm = "K";
  Char_t   cIdChamber[3];
           cIdChamber[2] = 0;

  TLorentzVector pos, mom;

  const Int_t    kNplan       = AliTRDgeometry::Nplan();
  const Int_t    kNcham       = AliTRDgeometry::Ncham();
  const Int_t    kNdetsec     = kNplan * kNcham;

  const Double_t kBig         = 1.0E+12; // Infinitely big
  const Float_t  kWion        = 23.53;   // Ionization energy
  const Float_t  kPTotMaxEl   = 0.002;   // Maximum momentum for e+ e- g

  // Minimum energy for the step size adjustment
  const Float_t  kEkinMinStep = 1.0e-5;
  // energy threshold for production of delta electrons
	const Float_t  kECut = 1.0e4;
	// Parameters entering the parametrized range for delta electrons
	const float ra=5.37E-4, rb=0.9815, rc=3.123E-3;
	// Gas density -> To be made user adjustable !
	const float rho=0.004945 ; //[0.85*0.00549+0.15*0.00186 (Xe-CO2 85-15)]

  // Plateau value of the energy-loss for electron in xenon
  // taken from: Allison + Comb, Ann. Rev. Nucl. Sci. (1980), 30, 253
  //const Double_t kPlateau = 1.70;
  // the averaged value (26/3/99)
  const Float_t  kPlateau     = 1.55;

  const Float_t  kPrim        = 19.34;  // dN1/dx|min for the gas mixture (90% Xe + 10% CO2)
  // First ionization potential (eV) for the gas mixture (90% Xe + 10% CO2)
  const Float_t  kPoti        = 12.1;

  const Int_t    kPdgElectron = 11;  // PDG code electron

  // Set the maximum step size to a very large number for all
  // neutral particles and those outside the driftvolume
  gMC->SetMaxStep(kBig);

  // Use only charged tracks
  if (( gMC->TrackCharge()       ) &&
      (!gMC->IsTrackStop()       ) &&
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
      if (phi < 90.)
        phi = phi + 270.;
      else
        phi = phi -  90.;
      sec = ((Int_t) (phi / 20));

      // The plane and chamber number
      cIdChamber[0] = cIdCurrent[2];
      cIdChamber[1] = cIdCurrent[3];
      Int_t idChamber = (atoi(cIdChamber) % kNdetsec);
      cha = kNcham - ((Int_t) idChamber / kNplan) - 1;
      pla = ((Int_t) idChamber % kNplan);

      // Check on selected volumes
      Int_t addthishit = 1;
      if (fSensSelect) {
        if ((fSensPlane   >= 0) && (pla != fSensPlane  )) addthishit = 0;
        if ((fSensChamber >= 0) && (cha != fSensChamber)) addthishit = 0;
        if (fSensSector  >= 0) {
          Int_t sens1  = fSensSector;
          Int_t sens2  = fSensSector + fSensSectorRange;
                sens2 -= ((Int_t) (sens2 / AliTRDgeometry::Nsect()))
                       * AliTRDgeometry::Nsect();
          if (sens1 < sens2) {
            if ((sec < sens1) || (sec >= sens2)) addthishit = 0;
	  			}
          else {
            if ((sec < sens1) && (sec >= sens2)) addthishit = 0;
	  			}
				}
      }

      // Add this hit
      if (addthishit) {

	// The detector number
        det = fGeometry->GetDetector(pla,cha,sec);

	// Special hits only in the drift region
        if (drRegion) {
          // Create a track reference at the entrance and
          // exit of each chamber that contain the
	        // momentum components of the particle
          if (gMC->IsTrackEntering() || gMC->IsTrackExiting()) {
            gMC->TrackMomentum(mom);
            AddTrackReference(gAlice->GetMCApp()->GetCurrentTrackNumber());
          }
					
					if (gMC->IsTrackEntering() && !gMC->IsNewTrack()) {
						// determine if hit belong to primary track 
						fPrimaryTrackPid=gAlice->GetMCApp()->GetCurrentTrackNumber();
						//determine track length when entering the detector
						fTrackLength0=gMC->TrackLength();
					}
					
					// Create the hits from TR photons
          if (fTR) CreateTRhit(det);
				}

				// Calculate the energy of the delta-electrons
				// modified by Alex Bercuci (A.Bercuci@gsi.de) on 26.01.06
				// take into account correlation with the underlying GEANT tracking
				// mechanism. see
        // http://www-linux.gsi.de/~abercuci/Contributions/TRD/index.html

				// determine the most significant process (last on the processes list)
				// which caused this hit
        TArrayI processes;
        gMC->StepProcesses(processes);
        int nofprocesses=processes.GetSize(), pid;
				if(!nofprocesses) pid=0;
				else pid=	processes[nofprocesses-1];		
				
				// generate Edep according to GEANT parametrisation
				eDelta =TMath::Exp(fDeltaG->GetRandom()) - kPoti;
        eDelta=TMath::Max(eDelta,0.0);
				float pr_range=0.;
				float range=gMC->TrackLength()-fTrackLength0;
				// merge GEANT tracker information with localy cooked one
				if(gAlice->GetMCApp()->GetCurrentTrackNumber()==fPrimaryTrackPid) {
//					printf("primary pid=%d eDelta=%f\n",pid,eDelta);
					if(pid==27){ 
						if(eDelta>=kECut){                
							pr_range=ra*eDelta*.001*(1.-rb/(1.+rc*eDelta*0.001))/rho;
        			if(pr_range>=(3.7-range)) eDelta*=.1;
						}
					} else if(pid==1){	
						if(eDelta<kECut) eDelta*=.5;
						else {                
							pr_range=ra*eDelta*.001*(1.-rb/(1.+rc*eDelta*0.001))/rho;
        			if(pr_range>=((AliTRDgeometry::DrThick()
                       + AliTRDgeometry::AmThick())-range)) eDelta*=.05;
							else eDelta*=.5;
						}
					} else eDelta=0.;	
				} else eDelta=0.;

        // Generate the electron cluster size
        if(eDelta==0.) qTot=0;
				else qTot = ((Int_t) (eDelta / kWion) + 1);
				// Create a new dEdx hit
        AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber(),det,hits,qTot, drRegion);
				
        // Calculate the maximum step size for the next tracking step
	// Produce only one hit if Ekin is below cutoff
        aMass = gMC->TrackMass();
        if ((gMC->Etot() - aMass) > kEkinMinStep) {

          // The energy loss according to Bethe Bloch
          iPdg  = TMath::Abs(gMC->TrackPid());
          if ( (iPdg != kPdgElectron) ||
	      			((iPdg == kPdgElectron) && (pTot < kPTotMaxEl))) {
            gMC->TrackMomentum(mom);
            pTot      = mom.Rho();
            betaGamma = pTot / aMass;
            pp        = BetheBlochGeant(betaGamma);
			// Take charge > 1 into account
            charge = gMC->TrackCharge();
            if (TMath::Abs(charge) > 1) pp = pp * charge*charge;
          } else { // Electrons above 20 Mev/c are at the plateau
						pp = kPrim * kPlateau;
          }

          stepSize = 1./gRandom->Poisson(pp);
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

  Double_t pTot = 0;
  Double_t eDelta;
  Double_t betaGamma, pp;
  Double_t stepSize;

  Bool_t   drRegion = kFALSE;
  Bool_t   amRegion = kFALSE;

  TString  cIdCurrent;
  TString  cIdSensDr = "J";
  TString  cIdSensAm = "K";
  Char_t   cIdChamber[3];
           cIdChamber[2] = 0;

  TLorentzVector pos, mom;

  const Int_t    kNplan       = AliTRDgeometry::Nplan();
  const Int_t    kNcham       = AliTRDgeometry::Ncham();
  const Int_t    kNdetsec     = kNplan * kNcham;

  const Double_t kBig         = 1.0E+12; // Infinitely big
  const Float_t  kWion        = 23.53;   // Ionization energy
  const Float_t  kPTotMaxEl   = 0.002;   // Maximum momentum for e+ e- g 

  // energy threshold for production of delta electrons
  //const Float_t  kECut = 1.0e4;
  // Parameters entering the parametrized range for delta electrons
  //const float ra=5.37E-4, rb=0.9815, rc=3.123E-3;
  // Gas density -> To be made user adjustable !
  //const float rho=0.004945 ; //[0.85*0.00549+0.15*0.00186 (Xe-CO2 85-15)]

  // Minimum energy for the step size adjustment
  const Float_t  kEkinMinStep = 1.0e-5;

  // Plateau value of the energy-loss for electron in xenon
  // taken from: Allison + Comb, Ann. Rev. Nucl. Sci. (1980), 30, 253
  //const Double_t kPlateau = 1.70;
  // the averaged value (26/3/99)
  const Float_t  kPlateau     = 1.55;

  const Float_t  kPrim        = 48.0;  // dN1/dx|min for the gas mixture (90% Xe + 10% CO2)
  // First ionization potential (eV) for the gas mixture (90% Xe + 10% CO2)
  const Float_t  kPoti        = 12.1;

  const Int_t    kPdgElectron = 11;  // PDG code electron

  // Set the maximum step size to a very large number for all 
  // neutral particles and those outside the driftvolume
  gMC->SetMaxStep(kBig); 

  // Use only charged tracks 
  if (( gMC->TrackCharge()       ) &&
      (!gMC->IsTrackStop()       ) && 
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
      if (phi < 90.) 
        phi = phi + 270.;
      else
        phi = phi -  90.;
      sec = ((Int_t) (phi / 20));

      // The plane and chamber number
      cIdChamber[0] = cIdCurrent[2];
      cIdChamber[1] = cIdCurrent[3];
      Int_t idChamber = (atoi(cIdChamber) % kNdetsec);
      cha = kNcham - ((Int_t) idChamber / kNplan) - 1;
      pla = ((Int_t) idChamber % kNplan);

      // Check on selected volumes
      Int_t addthishit = 1;
      if (fSensSelect) {
        if ((fSensPlane   >= 0) && (pla != fSensPlane  )) addthishit = 0;
        if ((fSensChamber >= 0) && (cha != fSensChamber)) addthishit = 0;
        if (fSensSector  >= 0) {
          Int_t sens1  = fSensSector;
          Int_t sens2  = fSensSector + fSensSectorRange;
                sens2 -= ((Int_t) (sens2 / AliTRDgeometry::Nsect())) 
                       * AliTRDgeometry::Nsect();
          if (sens1 < sens2) {
            if ((sec < sens1) || (sec >= sens2)) addthishit = 0;
	  			}
          else {
            if ((sec < sens1) && (sec >= sens2)) addthishit = 0;
	  			}
				}
      }

      // Add this hit
      if (addthishit) {

	// The detector number
        det = fGeometry->GetDetector(pla,cha,sec);

	// Special hits only in the drift region
        if (drRegion) {

          // Create a track reference at the entrance and
          // exit of each chamber that contain the 
	  // momentum components of the particle
          if (gMC->IsTrackEntering() || gMC->IsTrackExiting()) {
            gMC->TrackMomentum(mom);
            AddTrackReference(gAlice->GetMCApp()->GetCurrentTrackNumber());
          }
          // Create the hits from TR photons
          if (fTR) CreateTRhit(det);
				}


        // Calculate the energy of the delta-electrons
        eDelta = TMath::Exp(fDeltaE->GetRandom()) - kPoti;
        eDelta = TMath::Max(eDelta,0.0);
        // Generate the electron cluster size
        if(eDelta==0.) qTot=0;
				else qTot = ((Int_t) (eDelta / kWion) + 1);

				// Create a new dEdx hit
        if (drRegion) {
          AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber()
                ,det,hits,qTot, kTRUE);
				}
        else {
          AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber()
                ,det,hits,qTot,kFALSE);
				}

        // Calculate the maximum step size for the next tracking step
	// Produce only one hit if Ekin is below cutoff 
        aMass = gMC->TrackMass();
        if ((gMC->Etot() - aMass) > kEkinMinStep) {

          // The energy loss according to Bethe Bloch
          iPdg  = TMath::Abs(gMC->TrackPid());
          if ( (iPdg != kPdgElectron) ||
	      			((iPdg == kPdgElectron) && (pTot < kPTotMaxEl))) {
            gMC->TrackMomentum(mom);
            pTot      = mom.Rho();
            betaGamma = pTot / aMass;
            pp        = kPrim * BetheBloch(betaGamma);
	    // Take charge > 1 into account
            charge = gMC->TrackCharge();
            if (TMath::Abs(charge) > 1) pp = pp * charge*charge;
          } else { // Electrons above 20 Mev/c are at the plateau
						pp = kPrim * kPlateau;
          }
      
          if (pp > 0) {
            do 
            gMC->GetRandom()->RndmArray(1, random);
            while ((random[0] == 1.) || (random[0] == 0.));
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

  TLorentzVector pos, mom;

  const Int_t    kNplan       = AliTRDgeometry::Nplan();
  const Int_t    kNcham       = AliTRDgeometry::Ncham();
  const Int_t    kNdetsec     = kNplan * kNcham;

  const Double_t kBig         = 1.0E+12;

  const Float_t  kWion        = 23.53;   // Ionization energy
  const Float_t  kEkinMinStep = 1.0e-5;  // Minimum energy for the step size adjustment

  // Set the maximum step size to a very large number for all 
  // neutral particles and those outside the driftvolume
  gMC->SetMaxStep(kBig); 

  // If not charged track or already stopped or disappeared, just return.
  if ((!gMC->TrackCharge()) || 
        gMC->IsTrackStop()  || 
        gMC->IsTrackDisappeared()) return;

  // Inside a sensitive volume?
  cIdCurrent = gMC->CurrentVolName();

  if (cIdSensDr == cIdCurrent[1]) drRegion = kTRUE;
  if (cIdSensAm == cIdCurrent[1]) amRegion = kTRUE;

  if ((!drRegion) && (!amRegion)) return;

  // The hit coordinates and charge
  gMC->TrackPosition(pos);
  hits[0] = pos[0];
  hits[1] = pos[1];
  hits[2] = pos[2];

  // The sector number (0 - 17)
  // The numbering goes clockwise and starts at y = 0
  Float_t phi = kRaddeg*TMath::ATan2(pos[0],pos[1]);
  if (phi < 90.) phi += 270.;
  else           phi -=  90.;
  sec = ((Int_t) (phi / 20.));

  // The plane and chamber number
  cIdChamber[0] = cIdCurrent[2];
  cIdChamber[1] = cIdCurrent[3];
  Int_t idChamber = (atoi(cIdChamber) % kNdetsec);
  cha = kNcham - ((Int_t) idChamber / kNplan) - 1;
  pla = ((Int_t) idChamber % kNplan);

  // Check on selected volumes
  Int_t addthishit = 1;
  if(fSensSelect) {
    if ((fSensPlane   >= 0) && (pla != fSensPlane  )) addthishit = 0;
    if ((fSensChamber >= 0) && (cha != fSensChamber)) addthishit = 0;
    if (fSensSector  >= 0) {
      Int_t sens1  = fSensSector;
      Int_t sens2  = fSensSector + fSensSectorRange;
      sens2 -= ((Int_t) (sens2 / AliTRDgeometry::Nsect())) * AliTRDgeometry::Nsect();
      if (sens1 < sens2) {
        if ((sec < sens1) || (sec >= sens2)) addthishit = 0;
      }
      else {
        if ((sec < sens1) && (sec >= sens2)) addthishit = 0;
      }
    }
  }

  if (!addthishit) return;

  det = fGeometry->GetDetector(pla,cha,sec);  // The detector number
  
  Int_t trkStat = 0;  // 0: InFlight 1:Entering 2:Exiting

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
    if (fTR) CreateTRhit(det);    

  }
  
  // Calculate the charge according to GEANT Edep
  // Create a new dEdx hit
  eDep = TMath::Max(gMC->Edep(),0.0) * 1.0e+09;
  qTot = (Int_t) (eDep / kWion);
  AddHit(gAlice->GetMCApp()->GetCurrentTrackNumber()
        ,det,hits,qTot,drRegion);

  // Set Maximum Step Size
  // Produce only one hit if Ekin is below cutoff
  if ((gMC->Etot() - gMC->TrackMass()) < kEkinMinStep) return;
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
  const Double_t kP1 = 7.17960e-02;
  const Double_t kP2 = 8.54196;
  const Double_t kP3 = 1.38065e-06;
  const Double_t kP4 = 5.30972;
  const Double_t kP5 = 2.83798;

  // This parameters have been adjusted to Xe-data found in:
  // Allison & Cobb, Ann. Rev. Nucl. Sci. (1980), 30, 253
  //const Double_t kP1 = 0.76176E-1;
  //const Double_t kP2 = 10.632;
  //const Double_t kP3 = 3.17983E-6;
  //const Double_t kP4 = 1.8631;
  //const Double_t kP5 = 1.9479;

  // Lower cutoff of the Bethe-Bloch-curve to limit step sizes
  const Double_t kBgMin = 0.8;
  const Double_t kBBMax = 6.83298;
  //const Double_t kBgMin = 0.6;
  //const Double_t kBBMax = 17.2809;
  //const Double_t kBgMin = 0.4;
  //const Double_t kBBMax = 82.0;

  if (bg > kBgMin) {
    Double_t yy = bg / TMath::Sqrt(1. + bg*bg);
    Double_t aa = TMath::Power(yy,kP4);
    Double_t bb = TMath::Power((1./bg),kP5);
             bb = TMath::Log(kP3 + bb);
    return ((kP2 - aa - bb)*kP1 / aa);
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

  Double_t arr_g[20] = {
    1.100000,   1.200000,    1.300000,    1.500000,
    1.800000,   2.000000,    2.500000,    3.000000,
    4.000000,   7.000000,    10.000000,   20.000000,
    40.000000,  70.000000,   100.000000,  300.000000,
    600.000000, 1000.000000, 3000.000000, 10000.000000 };

  Double_t arr_nc[20] = {
    75.009056,   45.508083,   35.299252,   27.116327,
    22.734999,   21.411915,   19.934095,   19.449375,
    19.344431,   20.185553,   21.027925,   22.912676,
    24.933352,   26.504053,   27.387468,   29.566597,
    30.353779,   30.787134,   31.129285,   31.157350 };

  // betagamma to gamma
  Double_t g = TMath::Sqrt( 1. + bg*bg );

  // Find the index just before the point we need.
  int i;
  for( i = 0 ; i < 18 ; i++ )
    if( arr_g[i] < g && arr_g[i+1] > g )
      break;

  // Simple interpolation.
  Double_t pp = ((arr_nc[i+1] - arr_nc[i]) / 
		 (arr_g[i+1]-arr_g[i])) * (g-arr_g[i]) + arr_nc[i];

  return pp; //arr_nc[8];

}

//_____________________________________________________________________________
void AliTRDv1::Stepping()
{    
// Stepping info
// ---

   cout << "X(cm)    "
       << "Y(cm)    "
       << "Z(cm)  "
       << "KinE(MeV)   "
       << "dE(MeV) "
       << "Step(cm) "
       << "TrackL(cm) "
       << "Volume  "
       << "Process "  
       << endl;

   // Position
    //
    Double_t x, y, z;
    gMC->TrackPosition(x, y, z);
    cout << setw(8) << setprecision(3) << x << " "
         << setw(8) << setprecision(3) << y << " "
         << setw(8) << setprecision(3) << z << "  ";

    // Kinetic energy
    //
    Double_t px, py, pz, etot;
    gMC->TrackMomentum(px, py, pz, etot);
    Double_t ekin = etot - gMC->TrackMass();
    cout << setw(9) << setprecision(4) << ekin*1e03 << " ";

    // Energy deposit
    //
    cout << setw(9) << setprecision(4) << gMC->Edep()*1e03 << " ";

    // Step length
    //
    cout << setw(8) << setprecision(3) << gMC->TrackStep() << " ";

    // Track length
    //
    cout << setw(8) << setprecision(3) << gMC->TrackLength() << "     ";

    // Volume
    //
    if (gMC->CurrentVolName() != 0)
      cout << setw(4) << gMC->CurrentVolName() << "  ";
    else
      cout << setw(4) << "None"  << "  ";

    // Process
    //
    TArrayI processes;
    Int_t nofProcesses = gMC->StepProcesses(processes);
    for(int ip=0;ip<nofProcesses; ip++)
      cout << TMCProcessName[processes[ip]]<<" / ";

    cout << endl;
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

  Int_t    pos1, pos2;

  const Int_t kNv = 31;

  Float_t vxe[kNv] = { 2.3026, 2.9957, 3.4012, 3.6889, 3.9120  
                     , 4.0943, 4.2485, 4.3820, 4.4998, 4.6052
                     , 4.7005, 5.0752, 5.2983, 5.7038, 5.9915
                     , 6.2146, 6.5221, 6.9078, 7.3132, 7.6009
                     , 8.0064, 8.5172, 8.6995, 8.9872, 9.2103
                     , 9.4727, 9.9035,10.3735,10.5966,10.8198
                     ,11.5129 };

  Float_t vye[kNv] = { 80.0  , 31.0  , 23.3  , 21.1  , 21.0
                     , 20.9  , 20.8  , 20.0  , 16.0  , 11.0
                     ,  8.0  ,  6.0  ,  5.2  ,  4.6  ,  4.0
                     ,  3.5  ,  3.0  ,  1.4  ,  0.67 ,  0.44
                     ,  0.3  ,  0.18 ,  0.12 ,  0.08 ,  0.056
                     ,  0.04 ,  0.023,  0.015,  0.011,  0.01
		     ,  0.004 };

  energy = x[0];

  // Find the position 
  pos1 = pos2 = 0;
  dpos = 0;
  do {
    dpos = energy - vxe[pos2++];
  } 
  while (dpos > 0);
  pos2--; 
  if (pos2 > kNv) pos2 = kNv - 1;
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

  const Int_t n_pts = 83;
  Double_t arr_e[n_pts] = {
    2.421257,  2.483278,  2.534301,  2.592230,
    2.672067,  2.813299,  3.015059,  3.216819,
    3.418579,  3.620338,  3.868209,  3.920198,
    3.978284,  4.063923,  4.186264,  4.308605,
    4.430946,  4.553288,  4.724261,  4.837736,
    4.999842,  5.161949,  5.324056,  5.486163,
    5.679688,  5.752998,  5.857728,  5.962457,
    6.067185,  6.171914,  6.315653,  6.393674,
    6.471694,  6.539689,  6.597658,  6.655627,
    6.710957,  6.763648,  6.816338,  6.876198,
    6.943227,  7.010257,  7.106285,  7.252151,
    7.460531,  7.668911,  7.877290,  8.085670,
    8.302979,  8.353585,  8.413120,  8.483500,
    8.541030,  8.592857,  8.668865,  8.820485,
    9.037086,  9.253686,  9.470286,  9.686887,
    9.930838,  9.994655, 10.085822, 10.176990,
    10.268158, 10.359325, 10.503614, 10.627565,
    10.804637, 10.981709, 11.158781, 11.335854,
    11.593397, 11.781165, 12.049404, 12.317644,
    12.585884, 12.854123, 14.278421, 16.975889,
    20.829416, 24.682943, 28.536469
  };
  Double_t arr_dndx[n_pts] = {
    19.344431, 18.664679, 18.136106, 17.567745,
    16.836426, 15.677382, 14.281277, 13.140237,
    12.207677, 11.445510, 10.697049, 10.562296,
    10.414673, 10.182341,  9.775256,  9.172330,
    8.240271,  6.898587,  4.808303,  3.889751,
    3.345288,  3.093431,  2.897347,  2.692470,
    2.436222,  2.340029,  2.208579,  2.086489,
    1.975535,  1.876519,  1.759626,  1.705024,
    1.656374,  1.502638,  1.330566,  1.200697,
    1.101168,  1.019323,  0.943867,  0.851951,
    0.755229,  0.671576,  0.570675,  0.449672,
    0.326722,  0.244225,  0.188225,  0.149608,
    0.121529,  0.116289,  0.110636,  0.103490,
    0.096147,  0.089191,  0.079780,  0.063927,
    0.047642,  0.036341,  0.028250,  0.022285,
    0.017291,  0.016211,  0.014802,  0.013533,
    0.012388,  0.011352,  0.009803,  0.008537,
    0.007039,  0.005829,  0.004843,  0.004034,
    0.003101,  0.002564,  0.001956,  0.001494,
    0.001142,  0.000873,  0.000210,  0.000014,
    0.000000,  0.000000,  0.000000
  };

  Int_t i;
  Double_t energy = x[0];
  Double_t dnde;

  for( i = 0 ; i < n_pts ; i++ )
    if( energy < arr_e[i] ) break;

  if( i == 0 )
    AliErrorGeneral("AliTRDv1","Given energy value is too small or zero");

  // Differentiate
  dnde = (arr_dndx[i-1] - arr_dndx[i]) / (arr_e[i] - arr_e[i-1]);

  return dnde;

}
