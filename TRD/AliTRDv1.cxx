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
//  Transition Radiation Detector version 2 -- slow simulator                //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliTRDfullClass.gif">
*/
//End_Html
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TMath.h>
#include <TVector.h>
#include <TRandom.h>

#include "AliRun.h"
#include "AliMC.h"
#include "AliConst.h"

#include "AliTRDv1.h"
#include "AliTRDmatrix.h"
#include "AliTRDgeometry.h"

ClassImp(AliTRDv1)

//_____________________________________________________________________________
AliTRDv1::AliTRDv1(const char *name, const char *title) 
         :AliTRD(name, title) 
{
  //
  // Standard constructor for Transition Radiation Detector version 1
  //

  fIdSens        =  0;

  fIdChamber1    =  0;
  fIdChamber2    =  0;
  fIdChamber3    =  0;

  fSensSelect    =  0;
  fSensPlane     = -1;
  fSensChamber   = -1;
  fSensSector    = -1;

  fDeltaE        = NULL;

  SetBufferSize(128000);

}

//_____________________________________________________________________________
AliTRDv1::~AliTRDv1()
{

  if (fDeltaE) delete fDeltaE;

}
 
//_____________________________________________________________________________
void AliTRDv1::CreateGeometry()
{
  //
  // Create the GEANT geometry for the Transition Radiation Detector - Version 1
  // This version covers the full azimuth. 
  //

  // Check that FRAME is there otherwise we have no place where to put the TRD
  AliModule* FRAME = gAlice->GetModule("FRAME");
  if (!FRAME) return;

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
    if (fSensSector  >= 0)
      printf("          Only sector %d is sensitive\n",fSensSector);
  }
  printf("\n");

  // First ionization potential (eV) for the gas mixture (90% Xe + 10% CO2)
  const Float_t kPoti = 12.1;
  // Maximum energy (50 keV);
  const Float_t kEend = 50000.0;
  // Ermilova distribution for the delta-ray spectrum
  Float_t Poti = TMath::Log(kPoti);
  Float_t Eend = TMath::Log(kEend);
  fDeltaE  = new TF1("deltae",Ermilova,Poti,Eend,0);

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

  if ((isector < 0) || (isector > 17)) {
    printf("Wrong input value: %d\n",isector);
    printf("Use standard setting\n");
    fSensSector = -1;
    fSensSelect =  0;
    return;
  }

  fSensSelect = 1;
  fSensSector = isector;

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
  Int_t    iPdg;

  Float_t  hits[4];
  Float_t  random[1];
  Float_t  charge;
  Float_t  aMass;

  Double_t pTot;
  Double_t qTot;
  Double_t eDelta;
  Double_t betaGamma, pp;

  TLorentzVector pos, mom;
  TClonesArray  &lhits = *fHits;

  const Double_t kBig     = 1.0E+12;

  // Ionization energy
  const Float_t  kWion    = 22.04;
  // Maximum energy for e+ e- g for the step-size calculation
  const Float_t  kPTotMax = 0.002;
  // Plateau value of the energy-loss for electron in xenon
  // taken from: Allison + Comb, Ann. Rev. Nucl. Sci. (1980), 30, 253
  //const Double_t kPlateau = 1.70;
  // the averaged value (26/3/99)
  const Float_t  kPlateau = 1.55;
  // dN1/dx|min for the gas mixture (90% Xe + 10% CO2)
  const Float_t  kPrim    = 48.0;
  // First ionization potential (eV) for the gas mixture (90% Xe + 10% CO2)
  const Float_t  kPoti    = 12.1;

  // PDG code electron
  const Int_t    pdgElectron = 11;

  // Set the maximum step size to a very large number for all 
  // neutral particles and those outside the driftvolume
  gMC->SetMaxStep(kBig); 

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
      qTot = (Double_t) ((Int_t) (eDelta / kWion) + 1);

      // The hit coordinates and charge
      gMC->TrackPosition(pos);
      hits[0] = pos[0];
      hits[1] = pos[1];
      hits[2] = pos[2];
      hits[3] = qTot;

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
        if ((fSensSector  >= 0) && (sec != fSensSector )) addthishit = 0;
      }

      // Add this hit
      if (addthishit) {

        new(lhits[fNhits++]) AliTRDhit(fIshunt
                                      ,gAlice->CurrentTrack()
                                      ,fGeometry->GetDetector(pla,cha,sec)
                                      ,hits);

        // The energy loss according to Bethe Bloch
        gMC->TrackMomentum(mom);
        pTot = mom.Rho();
        iPdg = TMath::Abs(gMC->TrackPid());
        if ( (iPdg != pdgElectron) ||
	    ((iPdg == pdgElectron) && (pTot < kPTotMax))) {
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

  const Int_t nV = 31;

  Float_t vxe[nV] = { 2.3026, 2.9957, 3.4012, 3.6889, 3.9120  
                    , 4.0943, 4.2485, 4.3820, 4.4998, 4.6052
                    , 4.7005, 5.0752, 5.2983, 5.7038, 5.9915
                    , 6.2146, 6.5221, 6.9078, 7.3132, 7.6009
                    , 8.0064, 8.5172, 8.6995, 8.9872, 9.2103
                    , 9.4727, 9.9035,10.3735,10.5966,10.8198
                    ,11.5129 };

  Float_t vye[nV] = { 80.0  , 31.0  , 23.3  , 21.1  , 21.0
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
  if (pos2 > nV) pos2 = nV;
  pos1 = pos2 - 1;

  // Differentiate between the sampling points
  dnde = (vye[pos1] - vye[pos2]) / (vxe[pos2] - vxe[pos1]);

  return dnde;

}
