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

/*
$Log$
Revision 1.26  2000/11/30 17:38:08  cblume
Changes to get in line with new STEER and EVGEN

Revision 1.25  2000/11/15 14:30:16  cblume
Fixed bug in calculating detector no. of extra hit

Revision 1.24  2000/11/10 14:58:36  cblume
Introduce additional hit with amplitude 0 at the chamber borders

Revision 1.23  2000/11/01 14:53:21  cblume
Merge with TRD-develop

Revision 1.17.2.5  2000/10/15 23:40:01  cblume
Remove AliTRDconst

Revision 1.17.2.4  2000/10/06 16:49:46  cblume
Made Getters const

Revision 1.17.2.3  2000/10/04 16:34:58  cblume
Replace include files by forward declarations

Revision 1.17.2.2  2000/09/18 13:50:17  cblume
Include TR photon generation and adapt to new AliTRDhit

Revision 1.22  2000/06/27 13:08:50  cblume
Changed to Copy(TObject &A) to appease the HP-compiler

Revision 1.21  2000/06/09 11:10:07  cblume
Compiler warnings and coding conventions, next round

Revision 1.20  2000/06/08 18:32:58  cblume
Make code compliant to coding conventions

Revision 1.19  2000/06/07 16:27:32  cblume
Try to remove compiler warnings on Sun and HP

Revision 1.18  2000/05/08 16:17:27  cblume
Merge TRD-develop

Revision 1.17.2.1  2000/05/08 14:59:16  cblume
Made inline function non-virtual. Bug fix in setting sensitive chamber

Revision 1.17  2000/02/28 19:10:26  cblume
Include the new TRD classes

Revision 1.16.4.1  2000/02/28 18:04:35  cblume
Change to new hit version, introduce geometry class, and move digitization and clustering to AliTRDdigitizer/AliTRDclusterizerV1

Revision 1.16  1999/11/05 22:50:28  fca
Do not use Atan, removed from ROOT too

Revision 1.15  1999/11/02 17:20:19  fca
initialise nbytes before using it

Revision 1.14  1999/11/02 17:15:54  fca
Correct ansi scoping not accepted by HP compilers

Revision 1.13  1999/11/02 17:14:51  fca
Correct ansi scoping not accepted by HP compilers

Revision 1.12  1999/11/02 16:35:56  fca
New version of TRD introduced

Revision 1.11  1999/11/01 20:41:51  fca
Added protections against using the wrong version of FRAME

Revision 1.10  1999/09/29 09:24:35  fca
Introduction of the Copyright and cvs Log

*/

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

#include <TMath.h>
#include <TVector.h>
#include <TRandom.h>
#include <TF1.h>
#include <TLorentzVector.h>

#include "AliRun.h"
#include "AliMC.h"
#include "AliConst.h"

#include "AliTRDv1.h"
#include "AliTRDhit.h"
#include "AliTRDmatrix.h"
#include "AliTRDgeometry.h"
#include "AliTRDsim.h"

ClassImp(AliTRDv1)
 
//_____________________________________________________________________________
AliTRDv1::AliTRDv1():AliTRD()
{
  //
  // Default constructor
  //

  fIdSens          =  0;

  fIdChamber1      =  0;
  fIdChamber2      =  0;
  fIdChamber3      =  0;

  fSensSelect      =  0;
  fSensPlane       = -1;
  fSensChamber     = -1;
  fSensSector      = -1;
  fSensSectorRange =  0;

  fDeltaE          = NULL;
  fTR              = NULL;

}

//_____________________________________________________________________________
AliTRDv1::AliTRDv1(const char *name, const char *title) 
         :AliTRD(name, title) 
{
  //
  // Standard constructor for Transition Radiation Detector version 1
  //

  fIdSens          =  0;

  fIdChamber1      =  0;
  fIdChamber2      =  0;
  fIdChamber3      =  0;

  fSensSelect      =  0;
  fSensPlane       = -1;
  fSensChamber     = -1;
  fSensSector      = -1;
  fSensSectorRange =  0;

  fDeltaE          = NULL;
  fTR              = NULL;

  SetBufferSize(128000);

}

//_____________________________________________________________________________
AliTRDv1::AliTRDv1(const AliTRDv1 &trd)
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
void AliTRDv1::Copy(TObject &trd)
{
  //
  // Copy function
  //

  ((AliTRDv1 &) trd).fIdSens          = fIdSens;

  ((AliTRDv1 &) trd).fIdChamber1      = fIdChamber1;
  ((AliTRDv1 &) trd).fIdChamber2      = fIdChamber2;
  ((AliTRDv1 &) trd).fIdChamber3      = fIdChamber3;

  ((AliTRDv1 &) trd).fSensSelect      = fSensSelect;
  ((AliTRDv1 &) trd).fSensPlane       = fSensPlane;
  ((AliTRDv1 &) trd).fSensChamber     = fSensChamber;
  ((AliTRDv1 &) trd).fSensSector      = fSensSector;
  ((AliTRDv1 &) trd).fSensSectorRange = fSensSectorRange;

  fDeltaE->Copy(*((AliTRDv1 &) trd).fDeltaE);
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
  const Float_t kWion        = 22.04;

  // Maximum number of TR photons per track
  const Int_t   kNTR         = 50;

  TLorentzVector mom, pos;
  TClonesArray  &lhits = *fHits;

  // Create TR only for electrons 
  Int_t iPdg = gMC->TrackPid();
  if (TMath::Abs(iPdg) != kPdgElectron) return;

  // Create TR at the entrance of the chamber
  if (gMC->IsTrackEntering()) {

    Float_t eTR[kNTR];
    Int_t   nTR;

    // Create TR photons
    gMC->TrackMomentum(mom);
    Float_t pTot = mom.Rho();
    fTR->CreatePhotons(iPdg,pTot,nTR,eTR);
    if (nTR > kNTR) {
      printf("AliTRDv1::CreateTRhit -- ");
      printf("Boundary error: nTR = %d, kNTR = %d\n",nTR,kNTR);
      exit(1);
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
      absLength = gRandom->Exp(sigma);
      if (absLength < AliTRDgeometry::MyThick()) continue;

      // The absorbtion cross sections in the drift gas
      if (fGasMix == 1) {
        // Gas-mixture (Xe/CO2)
        Double_t muXe = fTR->GetMuXe(energyMeV);
        Double_t muCO = fTR->GetMuCO(energyMeV);
        sigma = (0.90 * muXe + 0.10 * muCO) * fGasDensity;
      }
      else {
        // Gas-mixture (Xe/Isobutane) 
        Double_t muXe = fTR->GetMuXe(energyMeV);
        Double_t muBu = fTR->GetMuBu(energyMeV);
        sigma = (0.97 * muXe + 0.03 * muBu) * fGasDensity;
      }

      // The distance after which the energy of the TR photon
      // is deposited.
      absLength = gRandom->Exp(sigma);
      if (absLength > AliTRDgeometry::DrThick()) continue;

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
      new(lhits[fNhits++]) AliTRDhit(fIshunt,gAlice->CurrentTrack()
                                    ,det,posHit,-q);

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

  printf("          Slow simulator\n\n");
  if (fSensSelect) {
    if (fSensPlane   >= 0)
      printf("          Only plane %d is sensitive\n",fSensPlane);
    if (fSensChamber >= 0)   
      printf("          Only chamber %d is sensitive\n",fSensChamber);
    if (fSensSector  >= 0) {
      Int_t sens1  = fSensSector;
      Int_t sens2  = fSensSector + fSensSectorRange;
            sens2 -= ((Int_t) (sens2 / AliTRDgeometry::Nsect())) 
                   * AliTRDgeometry::Nsect();
      printf("          Only sectors %d - %d are sensitive\n",sens1,sens2-1);
    }
  }
  if (fTR) 
    printf("          TR simulation on\n");
  else
    printf("          TR simulation off\n");
  printf("\n");

  // First ionization potential (eV) for the gas mixture (90% Xe + 10% CO2)
  const Float_t kPoti = 12.1;
  // Maximum energy (50 keV);
  const Float_t kEend = 50000.0;
  // Ermilova distribution for the delta-ray spectrum
  Float_t poti = TMath::Log(kPoti);
  Float_t eEnd = TMath::Log(kEend);
  fDeltaE = new TF1("deltae",Ermilova,poti,eEnd,0);

  // Identifier of the sensitive volume (drift region)
  fIdSens     = gMC->VolId("UL05");

  // Identifier of the TRD-driftchambers
  fIdChamber1 = gMC->VolId("UCIO");
  fIdChamber2 = gMC->VolId("UCIM");
  fIdChamber3 = gMC->VolId("UCII");

  for (Int_t i = 0; i < 80; i++) printf("*");
  printf("\n");

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
    printf("Wrong input value: %d\n",iplane);
    printf("Use standard setting\n");
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
    printf("Wrong input value: %d\n",ichamber);
    printf("Use standard setting\n");
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
    printf("Wrong input value <isector>: %d\n",isector);
    printf("Use standard setting\n");
    fSensSector      = -1;
    fSensSectorRange =  0;
    fSensSelect      =  0;
    return;
  }

  if ((nsector < 1) || (nsector > 18)) {
    printf("Wrong input value <nsector>: %d\n",nsector);
    printf("Use standard setting\n");
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
  // along its path across the drift volume. The step size is set acording
  // to Bethe-Bloch. The energy distribution of the delta electrons follows
  // a spectrum taken from Ermilova et al.
  //

  Int_t    iIdSens, icSens;
  Int_t    iIdSpace, icSpace;
  Int_t    iIdChamber, icChamber;
  Int_t    pla = 0;
  Int_t    cha = 0;
  Int_t    sec = 0;
  Int_t    det = 0;
  Int_t    iPdg;
  Int_t    qTot;

  Float_t  hits[3];
  Float_t  moms[3];
  Float_t  random[1];
  Float_t  charge;
  Float_t  aMass;

  Double_t pTot;
  Double_t eDelta;
  Double_t betaGamma, pp;

  TLorentzVector pos, mom;
  TClonesArray  &lhits = *fHits;

  const Double_t kBig     = 1.0E+12;

  // Ionization energy
  const Float_t  kWion        = 22.04;
  // Maximum momentum for e+ e- g 
  const Float_t  kPTotMaxEl   = 0.002;
  // Minimum momentum for the step size adjustment
  const Float_t  kPTotMinStep = 1.0e-5;
  // Plateau value of the energy-loss for electron in xenon
  // taken from: Allison + Comb, Ann. Rev. Nucl. Sci. (1980), 30, 253
  //const Double_t kPlateau = 1.70;
  // the averaged value (26/3/99)
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

  // Create some special hits with amplitude 0 at the entrance and
  // exit of each chamber that contain the momentum components of the particle
  if (gMC->TrackCharge() &&
     (gMC->IsTrackEntering() || gMC->IsTrackExiting())) {

    // Inside a sensitive volume?
    iIdSens = gMC->CurrentVolID(icSens);
    if (iIdSens == fIdSens) { 

      iIdSpace   = gMC->CurrentVolOffID(4,icSpace  );
      iIdChamber = gMC->CurrentVolOffID(1,icChamber);

      // The hit coordinates
      gMC->TrackPosition(pos);
      hits[0] = pos[0];
      hits[1] = pos[1];
      hits[2] = pos[2];

      // The track momentum
      gMC->TrackMomentum(mom);
      moms[0] = mom[0];
      moms[1] = mom[1];
      moms[2] = mom[2];

      // The sector number (0 - 17)
      // The numbering goes clockwise and starts at y = 0
      Float_t phi = kRaddeg*TMath::ATan2(pos[0],pos[1]);
      if (phi < 90.) 
        phi = phi + 270.;
      else
        phi = phi -  90.;
      sec = ((Int_t) (phi / 20));

      // The chamber number 
      //   0: outer left
      //   1: middle left
      //   2: inner
      //   3: middle right
      //   4: outer right
      if      (iIdChamber == fIdChamber1)
        cha = (hits[2] < 0 ? 0 : 4);
      else if (iIdChamber == fIdChamber2)       
        cha = (hits[2] < 0 ? 1 : 3);
      else if (iIdChamber == fIdChamber3)       
        cha = 2;

      // The plane number
      // The numbering starts at the innermost plane
      pla = icChamber - TMath::Nint((Float_t) (icChamber / 7)) * 6 - 1;

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
        det = fGeometry->GetDetector(pla,cha,sec);
        new(lhits[fNhits++]) AliTRDhit(fIshunt
                                      ,gAlice->CurrentTrack()
                                      ,det
                                      ,moms
                                      ,0);
      }

    }

  }

  // Use only charged tracks 
  if (( gMC->TrackCharge()       ) &&
      (!gMC->IsTrackStop()       ) && 
      (!gMC->IsTrackDisappeared())) {

    // Inside a sensitive volume?
    iIdSens = gMC->CurrentVolID(icSens);
    if (iIdSens == fIdSens) { 

      iIdSpace   = gMC->CurrentVolOffID(4,icSpace  );
      iIdChamber = gMC->CurrentVolOffID(1,icChamber);

      // Calculate the energy of the delta-electrons
      eDelta = TMath::Exp(fDeltaE->GetRandom()) - kPoti;
      eDelta = TMath::Max(eDelta,0.0);

      // The number of secondary electrons created
      qTot = ((Int_t) (eDelta / kWion) + 1);

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

      // The chamber number 
      //   0: outer left
      //   1: middle left
      //   2: inner
      //   3: middle right
      //   4: outer right
      if      (iIdChamber == fIdChamber1)
        cha = (hits[2] < 0 ? 0 : 4);
      else if (iIdChamber == fIdChamber2)       
        cha = (hits[2] < 0 ? 1 : 3);
      else if (iIdChamber == fIdChamber3)       
        cha = 2;

      // The plane number
      // The numbering starts at the innermost plane
      pla = icChamber - TMath::Nint((Float_t) (icChamber / 7)) * 6 - 1;

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

        det = fGeometry->GetDetector(pla,cha,sec);

        // Create the electron cluster from TR photons
        if (fTR) CreateTRhit(det);

        new(lhits[fNhits++]) AliTRDhit(fIshunt
                                      ,gAlice->CurrentTrack()
                                      ,det
                                      ,hits
                                      ,qTot);

        // The energy loss according to Bethe Bloch
        gMC->TrackMomentum(mom);
        pTot = mom.Rho();
        iPdg = TMath::Abs(gMC->TrackPid());
        if ( (iPdg != kPdgElectron) ||
	    ((iPdg == kPdgElectron) && (pTot < kPTotMaxEl))) {
          aMass     = gMC->TrackMass();
          betaGamma = pTot / aMass;
          pp        = kPrim * BetheBloch(betaGamma);
	  // Take charge > 1 into account
          charge = gMC->TrackCharge();
          if (TMath::Abs(charge) > 1) pp = pp * charge*charge;
        }
        // Electrons above 20 Mev/c are at the plateau
        else {
          pp = kPrim * kPlateau;
        }
      
        // Calculate the maximum step size for the next tracking step
        // introduce a lower momentum cut
        //if ((pp > 0) && (pTot > kPTotMinStep)) {
        if (pp > 0) {
          do 
            gMC->Rndm(random,1);
          while ((random[0] == 1.) || (random[0] == 0.));
          gMC->SetMaxStep( - TMath::Log(random[0]) / pp);
	}

      }
      else {
        // set step size to maximal value
        gMC->SetMaxStep(kBig); 
      }

    }

  }

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

  if (bg > 0) {
    Double_t yy = bg / TMath::Sqrt(1. + bg*bg);
    Double_t aa = TMath::Power(yy,kP4);
    Double_t bb = TMath::Power((1./bg),kP5);
             bb = TMath::Log(kP3 + bb);
    return ((kP2 - aa - bb)*kP1 / aa);
  }
  else
    return 0;

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
  if (pos2 > kNv) pos2 = kNv;
  pos1 = pos2 - 1;

  // Differentiate between the sampling points
  dnde = (vye[pos1] - vye[pos2]) / (vxe[pos2] - vxe[pos1]);

  return dnde;

}
