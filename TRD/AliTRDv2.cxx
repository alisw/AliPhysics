///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Transition Radiation Detector version 2 -- detailed simulation           //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliTRDv2Class.gif">
*/
//End_Html
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TMath.h>
#include <TVector.h>

#include "AliTRDv2.h"
#include "AliRun.h"
#include "AliMC.h"
#include "AliConst.h"

ClassImp(AliTRDv2)

//_____________________________________________________________________________
AliTRDv2::AliTRDv2(const char *name, const char *title) 
         :AliTRD(name, title) 
{
  //
  // Standard constructor for Transition Radiation Detector version 2
  //

  fIdSens      = 0;

  fIdSpace1    = 0;
  fIdSpace2    = 0;
  fIdSpace3    = 0;

  fIdChamber1  = 0;
  fIdChamber2  = 0;
  fIdChamber3  = 0;

  fSensSelect  = 0;
  fSensPlane   = 0;
  fSensChamber = 0;
  fSensSector  = 0;

  fDeltaE      = NULL;

  SetBufferSize(128000);

}

AliTRDv2::~AliTRDv2()
{

  if (fDeltaE)  delete fDeltaE;

}
 
//_____________________________________________________________________________
void AliTRDv2::CreateGeometry()
{
  //
  // Create the GEANT geometry for the Transition Radiation Detector - Version 2
  // This version covers the full azimuth. 
  //
  // Author:  Christoph Blume (C.Blume@gsi.de) 20/07/99 
  //

  Float_t xpos, ypos, zpos;

  // Check that FRAME is there otherwise we have no place where to put the TRD
  AliModule* FRAME = gAlice->GetModule("FRAME");
  if (!FRAME) return;

  // Define the chambers
  AliTRD::CreateGeometry();

  // Position the the TRD-sectors in all TRD-volumes in the spaceframe
  xpos     = 0.;
  ypos     = 0.;
  zpos     = 0.;
  gMC->Gspos("TRD ",1,"BTR1",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("TRD ",2,"BTR2",xpos,ypos,zpos,0,"ONLY");
  gMC->Gspos("TRD ",3,"BTR3",xpos,ypos,zpos,0,"ONLY");

}

//_____________________________________________________________________________
void AliTRDv2::CreateMaterials()
{
  //
  // Create materials for the Transition Radiation Detector version 2
  //

  AliTRD::CreateMaterials();

}

//_____________________________________________________________________________
void AliTRDv2::Init() 
{
  //
  // Initialise Transition Radiation Detector after geometry has been built
  //

  // First ionization potential (eV) for the gas mixture (90% Xe + 10% CO2)
  const Float_t kPoti = 12.1;
  // Maximum energy (50 keV);
  const Float_t kEend = 50000.0;

  AliTRD::Init();

  if (fSensPlane)
    printf("          Only plane %d is sensitive\n",fSensPlane);
  if (fSensChamber)   
    printf("          Only chamber %d is sensitive\n",fSensChamber);
  if (fSensSector)
    printf("          Only sector %d is sensitive\n",fSensSector);

  for(Int_t i = 0; i < 80; i++) printf("*");
  printf("\n");

  // Ermilova distribution for the delta-ray spectrum
  Float_t Poti = TMath::Log(kPoti);
  Float_t Eend = TMath::Log(kEend);
  fDeltaE  = new TF1("deltae",Ermilova,Poti,Eend,0);

  // Identifier of the sensitive volume (drift region)
  fIdSens     = gMC->VolId("UL05");

  // Identifier of the TRD-spaceframe volumina
  fIdSpace1   = gMC->VolId("B028");
  fIdSpace2   = gMC->VolId("B029");
  fIdSpace3   = gMC->VolId("B030");

  // Identifier of the TRD-driftchambers
  fIdChamber1 = gMC->VolId("UCIO");
  fIdChamber2 = gMC->VolId("UCIM");
  fIdChamber3 = gMC->VolId("UCII");

}

//_____________________________________________________________________________
void AliTRDv2::SetSensPlane(Int_t iplane)
{
  //
  // Defines the hit-sensitive plane (1-6)
  //

  if ((iplane < 0) || (iplane > 6)) {
    printf("Wrong input value: %d\n",iplane);
    printf("Use standard setting\n");
    fSensPlane  = 0;
    fSensSelect = 0;
    return;
  }

  fSensSelect = 1;
  fSensPlane  = iplane;

}

//_____________________________________________________________________________
void AliTRDv2::SetSensChamber(Int_t ichamber)
{
  //
  // Defines the hit-sensitive chamber (1-5)
  //

  if ((ichamber < 0) || (ichamber > 5)) {
    printf("Wrong input value: %d\n",ichamber);
    printf("Use standard setting\n");
    fSensChamber = 0;
    fSensSelect  = 0;
    return;
  }

  fSensSelect  = 1;
  fSensChamber = ichamber;

}

//_____________________________________________________________________________
void AliTRDv2::SetSensSector(Int_t isector)
{
  //
  // Defines the hit-sensitive sector (1-18)
  //

  if ((isector < 0) || (isector > 18)) {
    printf("Wrong input value: %d\n",isector);
    printf("Use standard setting\n");
    fSensSector = 0;
    fSensSelect = 0;
    return;
  }

  fSensSelect = 1;
  fSensSector = isector;

}

//_____________________________________________________________________________
void AliTRDv2::StepManager()
{
  //
  // Called at every step in the Transition Radiation Detector version 2.
  // Slow simulator. Every charged track produces electron cluster as hits 
  // along its path across the drift volume. The step size is set acording
  // to Bethe-Bloch. The energy distribution of the delta electrons follows
  // a spectrum taken from Ermilova et al.
  //

  Int_t    iIdSens, icSens;
  Int_t    iIdSpace, icSpace;
  Int_t    iIdChamber, icChamber;
  Int_t    vol[3]; 
  Int_t    iPid;

  Int_t    secMap1[10] = {  3,  7,  8,  9, 10, 11,  2,  1, 18, 17 };
  Int_t    secMap2[ 5] = { 16, 15, 14, 13, 12 };
  Int_t    secMap3[ 3] = {  5,  6,  4 };

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

  const Double_t kBig = 1.0E+12;

  // Ionization energy
  const Float_t kWion    = 22.04;
  // Maximum energy for e+ e- g for the step-size calculation
  const Float_t kPTotMax = 0.002;
  // Plateau value of the energy-loss for electron in xenon
  // taken from: Allison + Comb, Ann. Rev. Nucl. Sci. (1980), 30, 253
  //const Double_t kPlateau = 1.70;
  // the averaged value (26/3/99)
  const Float_t kPlateau = 1.55;
  // dN1/dx|min for the gas mixture (90% Xe + 10% CO2)
  const Float_t kPrim    = 48.0;
  // First ionization potential (eV) for the gas mixture (90% Xe + 10% CO2)
  const Float_t kPoti    = 12.1;

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

      // The sector number
      if      (iIdSpace == fIdSpace1) 
        vol[0] = secMap1[icSpace-1];
      else if (iIdSpace == fIdSpace2) 
        vol[0] = secMap2[icSpace-1];
      else if (iIdSpace == fIdSpace3) 
        vol[0] = secMap3[icSpace-1];

      // The chamber number 
      //   1: outer left
      //   2: middle left
      //   3: inner
      //   4: middle right
      //   5: outer right
      if      (iIdChamber == fIdChamber1)
        vol[1] = (hits[2] < 0 ? 1 : 5);
      else if (iIdChamber == fIdChamber2)       
        vol[1] = (hits[2] < 0 ? 2 : 4);
      else if (iIdChamber == fIdChamber3)       
        vol[1] = 3;

      // The plane number
      vol[2] = icChamber - TMath::Nint((Float_t) (icChamber / 7)) * 6;

      // Check on selected volumes
      Int_t addthishit = 1;
      if (fSensSelect) {
        if ((fSensPlane)   && (vol[2] != fSensPlane  )) addthishit = 0;
        if ((fSensChamber) && (vol[1] != fSensChamber)) addthishit = 0;
        if ((fSensSector)  && (vol[0] != fSensSector )) addthishit = 0;
      }

      // Add this hit
      if (addthishit) {

        new(lhits[fNhits++]) AliTRDhit(fIshunt,gAlice->CurrentTrack(),vol,hits);

        // The energy loss according to Bethe Bloch
        gMC->TrackMomentum(mom);
        pTot = mom.Rho();
        iPid = gMC->TrackPid();
        if ( (iPid >  3) ||
	    ((iPid <= 3) && (pTot < kPTotMax))) {
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
Double_t AliTRDv2::BetheBloch(Double_t bg) 
{
  //
  // Parametrization of the Bethe-Bloch-curve
  // The parametrization is the same as for the TPC and is taken from Lehrhaus.
  //

  // The parameters have been adjusted to Xe-data found in:
  // Allison & Cobb, Ann. Rev. Nucl. Sci. (1980), 30, 253
  //const Double_t kP1 = 0.76176E-1;
  //const Double_t kP2 = 10.632;
  //const Double_t kP3 = 3.17983E-6;
  //const Double_t kP4 = 1.8631;
  //const Double_t kP5 = 1.9479;

  // This parameters have been adjusted to averaged values from GEANT
  const Double_t kP1 = 7.17960e-02;
  const Double_t kP2 = 8.54196;
  const Double_t kP3 = 1.38065e-06;
  const Double_t kP4 = 5.30972;
  const Double_t kP5 = 2.83798;

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
  // Calculates the delta-ray energy distribution according to Ermilova
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
