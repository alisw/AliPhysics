///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Time Projection Chamber version 3 -- detailed TPC and slow simulation    //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliTPCv3Class.gif">
*/
//End_Html
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <TMath.h>

#include "AliTPCv3.h"
#include "AliRun.h"
#include "AliConst.h"

ClassImp(AliTPCv3)
 
//_____________________________________________________________________________
AliTPCv3::AliTPCv3(const char *name, const char *title) :
  AliTPC(name, title) 
{
  //
  // Standard constructor for Time Projection Chamber version 2
  //

  SetBufferSize(128000);
}
 
//_____________________________________________________________________________
void AliTPCv3::CreateGeometry()
{
  //
  // Create the geometry of Time Projection Chamber version 3
  //
  //Begin_Html
  /*
    <img src="picts/AliTPCv3.gif">
  */
  //End_Html
  //Begin_Html
  /*
    <img src="picts/AliTPCv3Tree.gif">
  */
  //End_Html

  Int_t *idtmed = fIdtmed->GetArray()-399;

  Float_t tana, rlsl, wlsl, rssl, rlsu, wssl, wlsu,
    rssu, wssu, alpha, x, y, sec_thick;
  
  Float_t x1, z0, z1, x2, theta1, theta2, theta3, dm[21];
  Int_t il, iu;
  Float_t z_side;
  Int_t idrotm[100];
  
  Float_t x0l, x0u;
  Int_t idr;
  //Float_t thl, thu;
  Float_t opl, opu, phi1, phi2, phi3;
  
  // ---------------------------------------------------- 
  //          FIELD CAGE WITH ENDCAPS - CARBON FIBER 
  //          THIS IS ALSO A TPC MOTHER VOLUME 
  // ---------------------------------------------------- 
  dm[0] = 76.;
  dm[1] = 278.;
  dm[2] = 275.;
  
  gMC->Gsvolu("TPC ", "TUBE", idtmed[407], dm, 3);
  // ------------------------------------------------------- 
  //     drift gas Ne/CO2 (90/10 volume) - sensitive 
  //     field cage thickness = 0.52% X0 
  // ---------------------------------------------------- 
  dm[0] = 76.+0.09776;
  dm[1] = 257.;
  dm[2] = 250.;

  gMC->Gsvolu("TGAS", "TUBE", idtmed[403], dm, 3);
  // ------------------------------------------------------ 
  //     "side" gas volume (the same as drift gas) 
  //     here the readout chambers are positioned 
  // ------------------------------------------------------ 
  dm[2]  = 0.5*(275.-250.);
  z_side = dm[2];
  
  gMC->Gsvolu("TPSG", "TUBE", idtmed[401], dm, 3);
  // ------------------------------------------------------ 
  //      HV midplane - 20 microns of mylar 
  // ----------------------------------------------------- 
  dm[2] = .001;
  
  gMC->Gsvolu("TPHV", "TUBE", idtmed[405], dm, 3);
  
  // ==================================================== 
  //   lower and upper readout chambers 
  // ==================================================== 
  //   sectros opening angles in degrees 
  // --------------------------------------------------- 
  opl = 30.;
  opu = 15.;
  //thl = TMath::Tan(opl * .5 * kDegrad);
  //thu = TMath::Tan(opu * .5 * kDegrad);
  // --------------------------------------------------- 
  //         S and L-sectors radii 
  // --------------------------------------------------- 
  rssl = 88.;
  rssu = 136.;
  rlsl = 142.;
  rlsu = 250.;
  // -------------------------------------------------- 
  //          Sectors widths 
  // -------------------------------------------------- 
  wssl = 46.5;
  wssu = 72.2;
  wlsl = 37.;
  wlsu = 65.4;
  // --------------------------------------------------- 
  //    Sector thickness 25% of X0 (Al) 
  // --------------------------------------------------- 
  sec_thick = 2.225;
  // --------------------------------------------------- 
  //     S-sectors readout chambers (lower sectors) 
  // --------------------------------------------------- 
  dm[0] = wssl * .5;
  dm[1] = wssu * .5;
  dm[2] = sec_thick * .5;
  dm[3] = (rssu - rssl) * .5;
  
  x0l = rssl + dm[3];
  
  gMC->Gsvolu("TRCS", "TRD1", idtmed[399], dm, 4);
  // --------------------------------------------------- 
  //     L-sectors readout chambers (upper sectors) 
  // --------------------------------------------------- 
  dm[0] = wlsl * .5;
  dm[1] = wlsu * .5;
  dm[2] = sec_thick * .5;
  dm[3] = (rlsu - rlsl) * .5;
  
  x0u = rlsl + dm[3];
  
  gMC->Gsvolu("TRCL", "TRD1", idtmed[399], dm, 4);
  // ---------------------------------------------------- 
  //    positioning of the S-sector readout chambers 
  //    rotation matices 1-12 
  // ---------------------------------------------------- 
  z1 = -z_side + sec_thick * .5;

  for (il = 1; il <= 12; ++il) {
    phi1 = (il - 1) * opl + 270.;
    if (phi1 > 360.) {
      phi1 += -360.;
    }
    theta1 = 90.;
    phi2   = 90.;
    theta2 = 180.;
    phi3   = (il - 1) * opl;
    theta3 = 90.;
    
    idr = il;
    
    alpha = (il - 1) * opl * kDegrad;
    x     = x0l * TMath::Cos(alpha);
    y     = x0l * TMath::Sin(alpha);
    
    AliMatrix(idrotm[idr], theta1, phi1, theta2, phi2, theta3, phi3);
    gMC->Gspos("TRCS", il, "TPSG", x, y, z1, idrotm[idr], "ONLY");
  }
  // ---------------------------------------------------- 
  //    positioning of the L-sector readout chambers 
  //    rotation matices 13-36 
  // ---------------------------------------------------- 
  for (iu = 1; iu <= 24; ++iu) {
    phi1 = (iu - 1) * opu + 270.;
    if (phi1 > 360.) {
      phi1 += -360.;
    }
    theta1 = 90.;
    phi2   = 90.;
    theta2 = 180.;
    phi3   = (iu - 1) * opu;
    theta3 = 90.;
    
    idr = iu + 12;
    AliMatrix(idrotm[idr], theta1, phi1, theta2, phi2, theta3, phi3);
    
    alpha = (iu - 1) * opu * kDegrad;
    x     = x0u * TMath::Cos(alpha);
    y     = x0u * TMath::Sin(alpha);
    
    gMC->Gspos("TRCL", iu, "TPSG", x, y, z1, idrotm[idr], "ONLY");
  }
  // -------------------------------------------------------- 
  //             Spoke wheel structures 
  // -------------------------------------------------------- 
  gMC->Gsvolu("TSWS", "TUBE", idtmed[399], dm, 0);
  
  z0 = -z_side + 2.;
  
  dm[0] = 82.;
  dm[1] = 86.;
  dm[2] = 1.;
  
  gMC->Gsposp("TSWS", 1, "TPSG", 0, 0, z0, 0, "ONLY", dm, 3);
  
  dm[0] = 253.;
  dm[1] = 257.;
  
  gMC->Gsposp("TSWS", 2, "TPSG", 0, 0, z0, 0, "ONLY", dm, 3);
  
  dm[0] = 140.9;
  dm[1] = 141.9;
  
  gMC->Gsposp("TSWS", 3, "TPSG", 0, 0, z0, 0, "ONLY", dm, 3);
  
  // ------------------------------------------------------- 
  //    this volumes are to avoid overlaping 
  // ------------------------------------------------------- 
  z0 = 253.;
  
  dm[0] = 76.;
  dm[1] = 76.+0.09776;
  
  gMC->Gsposp("TSWS", 4, "TPC ", 0, 0, z0, 0, "ONLY", dm, 3);
  gMC->Gsposp("TSWS", 5, "TPC ", 0, 0, -z0,0, "ONLY", dm, 3);
  
  z0 += 21.;
  
  gMC->Gsposp("TSWS", 6, "TPC ", 0, 0, z0, 0, "ONLY", dm, 3);
  gMC->Gsposp("TSWS", 7, "TPC ", 0, 0, -z0, 0, "ONLY", dm, 3);
  
  dm[0] = 257.;
  dm[1] = 257.+0.09776;
  dm[2] = 11.5;
  
  z0 = 263.5;
  
  gMC->Gsposp("TSWS", 8, "TPC ", 0, 0, z0, 0, "ONLY", dm, 3);
  gMC->Gsposp("TSWS", 9, "TPC ", 0, 0, -z0, 0, "ONLY", dm, 3);
  // ========================================================== 
  //                  wheels 
  // ========================================================== 
  // ---------------------------------------------------------- 
  //       Large wheel -> positioned in the TPC 
  // ---------------------------------------------------------- 
  dm[0] = 257.+0.09776;
  dm[1] = 278.;
  dm[2] = 11.5;
  gMC->Gsvolu("TPW1", "TUBE", idtmed[399], dm, 3);
  
  dm[0] = 259.;
  dm[1] = 278.;
  dm[2] = 9.5;
  
  gMC->Gsvolu("TPW2", "TUBE", idtmed[498], dm, 3);
  
  gMC->Gspos("TPW2", 1, "TPW1", 0, 0, 0, 0, "ONLY");
  
  gMC->Gspos("TPW1", 1, "TPC ", 0, 0, z0, 0, "ONLY");
  gMC->Gspos("TPW1", 2, "TPC ", 0, 0, -z0, 0, "ONLY");
  // ----------------------------------------------------------- 
  //     Small wheel -> positioned in the TPSG 
  // ----------------------------------------------------------- 
  dm[0] = 76.+0.09776;
  dm[1] = 82.;
  dm[2] = 11.5;
  
  gMC->Gsvolu("TPW3", "TUBE", idtmed[399], dm, 3);
  
  dm[0] = 76.+0.09776;
  dm[1] = 80.;
  dm[2] = 9.5;
  
  gMC->Gsvolu("TPW4", "TUBE", idtmed[401], dm, 3);
  
  gMC->Gspos("TPW4", 1, "TPW3", 0, 0, 0, 0, "ONLY");
  
  z0 = 1.;
  
  gMC->Gspos("TPW3", 1, "TPSG", 0, 0, z0, 0, "ONLY");
  // --------------------------------------------------------- 
  //       spokes, inner and outer, also the inner ring 
  // --------------------------------------------------------- 
  dm[0] = 0.5*(135.9-82.1);
  dm[1] = 3.;
  dm[2] = 2.;
  
  x1 = dm[0] + 82.;
  
  gMC->Gsvolu("TSPI", "BOX ", idtmed[399], dm, 3);
  
  dm[1] = 2.;
  dm[2] = 1.;
  
  gMC->Gsvolu("TSP1", "BOX ", idtmed[498], dm, 3);

  gMC->Gspos("TSP1", 1, "TSPI", 0, 0, 0, 0, "ONLY");
  
  dm[0] = 0.5*(256.9-142.1);
  dm[1] = 3.;
  dm[2] = 2.;
  
  x2 = dm[0] + 142.;
  
  gMC->Gsvolu("TSPO", "BOX ", idtmed[399], dm, 3);
  
  dm[1] = 2.;
  dm[2] = 1.;
  
  gMC->Gsvolu("TSP2", "BOX ", idtmed[498], dm, 3);
  
  gMC->Gspos("TSP2", 1, "TSPO", 0, 0, 0, 0, "ONLY");
  // -------------------------------------------------------- 
  dm[0] = 136.;
  dm[1] = 142.;
  dm[2] = 2.;
  
  gMC->Gsvolu("TSWH", "TUBE", idtmed[399], dm, 3);
  
  dm[0] = 137.;
  dm[1] = 141.;
  dm[2] = 1.;
  
  gMC->Gsvolu("TSW1", "TUBE", idtmed[498], dm, 3);
  
  gMC->Gspos("TSW1", 1, "TSWH", 0, 0, 0, 0, "ONLY");
  
  z0 = z_side - .16168 - 2.;
  // -------------------------------------------------------- 
  gMC->Gspos("TSWH", 1, "TPSG", 0, 0, z0, 0, "ONLY");
  // ------------------------------------------------------- 
  //     posiioning of the inner spokes 
  // ------------------------------------------------------- 
  for (il = 1; il <= 6; ++il) {
    phi1 = opl * .5 + (il - 1) * 2. * opl;
    theta1 = 90.;
    phi2 = opl * .5 + 90. + (il - 1) * 2. * opl;
    if (phi2 > 360.) {
      phi2 += -360.;
    }
    theta2 = 90.;
    phi3   = 0.;
    theta3 = 0.;
    
    alpha = phi1 * kDegrad;
    x     = x1 * TMath::Cos(alpha);
    y     = x1 * TMath::Sin(alpha);
    
    idr = il + 36;
    
    AliMatrix(idrotm[idr], theta1, phi1, theta2, phi2, theta3, phi3);
    gMC->Gspos("TSPI", il, "TPSG", x, y, z0, idrotm[idr], "ONLY");
    
  }
  
  for (iu = 1; iu <= 12; ++iu) {
    phi1 = opu * .5 + (iu - 1) * 2. * opu;
    theta1 = 90.;
    phi2 = opu * .5 + 90. + (iu - 1) * 2. * opu;
    if (phi2 > 360.) {
      phi2 += -360.;
    }
    theta2 = 90.;
    phi3   = 0.;
    theta3 = 0.;
    
    alpha = phi1 * kDegrad;
    x     = x2 * TMath::Cos(alpha);
    y     = x2 * TMath::Sin(alpha);
    
    idr = iu + 42;
    
    AliMatrix(idrotm[idr], theta1, phi1, theta2, phi2, theta3, phi3);
    gMC->Gspos("TSPO", iu, "TPSG", x, y, z0, idrotm[idr], "ONLY");
  }
  // -------------------------------------------------------- 
  //       endcap cover (C, 0.86% X0) 
  // -------------------------------------------------------- 
  dm[0] = 76.+0.09776;
  dm[1] = 257.;
  dm[2] = 0.16168*0.5;
  
  gMC->Gsvolu("TCOV", "TUBE", idtmed[407], dm, 3);
  
  z0 = z_side - dm[2];
  
  gMC->Gspos("TCOV", 1, "TPSG", 0, 0, z0, 0, "ONLY");
  // -------------------------------------------------------- 
  //         put the readout chambers into the TPC 
  // -------------------------------------------------------- 
  theta1 = 90.;
  phi1   = 0.;
  theta2 = 90.;
  phi2   = 270.;
  theta3 = 180.;
  phi3   = 0.;
  
  AliMatrix(idrotm[55], theta1, phi1, theta2, phi2, theta3, phi3);
  
  z0 = z_side + 250.;
  
  gMC->Gspos("TPSG", 1, "TPC ", 0, 0, z0, 0, "ONLY");
  gMC->Gspos("TPSG", 2, "TPC ", 0, 0, -z0, idrotm[55], "ONLY");
  // --------------------------------------------------------- 
  //     outer gas insulation (CO2) 
  // --------------------------------------------------------- 
  dm[0] = 257.+0.09776;
  dm[1] = 278.-0.25004;
  dm[2] = 275.-23.;
  
  gMC->Gsvolu("TPOI", "TUBE", idtmed[406], dm, 3);
  
  gMC->Gspos("TPHV", 1, "TGAS", 0, 0, 0, 0, "ONLY");
  gMC->Gspos("TGAS", 1, "TPC ", 0, 0, 0, 0, "ONLY");
  gMC->Gspos("TPOI", 1, "TPC ", 0, 0, 0, 0, "ONLY");
  
  gMC->Gspos("TPC ", 1, "ALIC", 0, 0, 0, 0, "ONLY");
  // ====================================================== 
  //      all volumes below are positioned in ALIC 
  // ====================================================== 
  // ------------------------------------------------------ 
  //        the last parts of the smaller wheel (TSWS) 
  // ------------------------------------------------------ 
  dm[0] = 74.;
  dm[1] = 76.;
  dm[2] = 1.;
  
  z0 = 253.;
  
  gMC->Gsposp("TSWS", 10, "TPC ", 0, 0, z0, 0, "ONLY", dm, 3);
  gMC->Gsposp("TSWS", 11, "TPC ", 0, 0, -z0, 0, "ONLY", dm, 3);
  
  dm[0] = 70.;
  
  z0 += 21.;
  
  gMC->Gsposp("TSWS", 12, "TPC ", 0, 0, z0, 0, "ONLY", dm, 3);
  gMC->Gsposp("TSWS", 13, "TPC ", 0, 0, -z0, 0, "ONLY", dm, 3);
  // ---------------------------------------------------- 
  //             Inner vessel (PCON) 
  //   This volume is to be positioned directly in ALIC 
  // ---------------------------------------------------- 
  dm[0] = 0.;
  dm[1] = 360.;
  dm[2] = 4.;
  
  dm[3] = -250.;
  dm[4] = 75.;
  dm[5] = 76.;
  
  dm[6] = -64.5;
  dm[7] = 50.;
  dm[8] = 76.;
  
  dm[9] = 64.5;
  dm[10] = 50.;
  dm[11] = 76.;
  
  dm[12] = 250.;
  dm[13] = 75.;
  dm[14] = 76.;
  
  gMC->Gsvolu("TPIV", "PCON", idtmed[407], dm, 15);
  // -------------------------------------------------------- 
  //     fill the inner vessel with CO2, (HV kDegrader) 
  //     cone parts have different thickness 
  //     than the central barrel, according to the TP 
  // -------------------------------------------------------- 
  tana = 75./185.5;
  
  dm[0] = 0.;
  dm[1] = 360.;
  dm[2] = 6.;
  
  dm[3] = -(250.-0.2162);
  dm[4] = (185.5-0.2126)*tana+0.2126;
  dm[5] = 76-0.001;
  
  dm[6] = -64.5;
  dm[7] = 50.+0.2162;
  dm[8] = 76-0.001;
  
  dm[9]  = -64.5;
  dm[10] = 50+0.05076;
  dm[11] = 76-0.001;
  
  dm[12] = 64.5;
  dm[13] = 50+0.05076;
  dm[14] = 76-0.001;
  
  dm[15] = 64.5;
  dm[16] = 50.+0.2162;
  dm[17] = 76-0.001;
  
  dm[18] = (250.-0.2162);
  dm[19] = (185.5-0.2126)*tana+0.2126;
  dm[20] = 76-0.001;
  
  gMC->Gsvolu("TPVD", "PCON", idtmed[406], dm, 21);
  
  gMC->Gspos("TPVD", 1, "TPIV", 0, 0, 0, 0, "ONLY");
  
  gMC->Gspos("TPIV", 1, "ALIC", 0, 0, 0, 0, "ONLY");
  // --------------------------------------------------- 
  //               volumes ordering 
  // --------------------------------------------------- 
  gMC->Gsord("TPSG", 6);
}
 
 
//_____________________________________________________________________________
void AliTPCv3::DrawDetector()
{
  //
  // Draw a shaded view of the Time Projection Chamber version 1
  //


  // Set everything unseen
  gMC->Gsatt("*", "seen", -1);
  // 
  // Set ALIC mother transparent
  gMC->Gsatt("ALIC","SEEN",0);
  //
  // Set the volumes visible
  gMC->Gsatt("TPC","SEEN",0);
  gMC->Gsatt("TGAS","SEEN",0);
  gMC->Gsatt("TPSG","SEEN",0);
  gMC->Gsatt("TPHV","SEEN",1);
  gMC->Gsatt("TRCS","SEEN",1);
  gMC->Gsatt("TRCL","SEEN",1);
  gMC->Gsatt("TSWS","SEEN",1);
  gMC->Gsatt("TPW1","SEEN",1);
  gMC->Gsatt("TPW2","SEEN",1);
  gMC->Gsatt("TPW3","SEEN",1);
  gMC->Gsatt("TPW4","SEEN",1);
  gMC->Gsatt("TSPI","SEEN",1);
  gMC->Gsatt("TSP1","SEEN",0);
  gMC->Gsatt("TSPO","SEEN",1);
  gMC->Gsatt("TSP2","SEEN",0);
  gMC->Gsatt("TSWH","SEEN",1);
  gMC->Gsatt("TSW1","SEEN",1);
  gMC->Gsatt("TPOI","SEEN",1);
  gMC->Gsatt("TPIV","SEEN",1);
  gMC->Gsatt("TPVD","SEEN",1);
  //
  gMC->Gdopt("hide", "on");
  gMC->Gdopt("shad", "on");
  gMC->Gsatt("*", "fill", 7);
  gMC->SetClipBox(".");
  gMC->SetClipBox("*", 0, 1000, -1000, 1000, -1000, 1000);
  gMC->DefaultRange();
  gMC->Gdraw("alic", 40, 30, 0, 12, 9.5, .025, .025);
  gMC->Gdhead(1111, "Time Projection Chamber");
  gMC->Gdman(18, 4, "MAN");
  gMC->Gdopt("hide","off");
}

//_____________________________________________________________________________
void AliTPCv3::CreateMaterials()
{
  //
  // Define materials for version 2 of the Time Projection Chamber
  //


  //
  // Increase maximum number of steps
  gMC->SetMaxNStep(30000);
  //
  AliTPC::CreateMaterials();
}

//_____________________________________________________________________________
void AliTPCv3::Init()
{
  //
  // Initialises version 3 of the TPC after that it has been built
  //
  Int_t *idtmed = fIdtmed->GetArray()-399;

  AliTPC::Init();

  fIdSens1=gMC->VolId("TGAS"); // drift gas as a sensitive volume

  gMC->SetMaxNStep(30000); // max. number of steps increased

  gMC->Gstpar(idtmed[403],"LOSS",5);

  printf("*** TPC version 3 initialized ***\n");
  printf("Maximum number of steps = %d\n",gMC->GetMaxNStep());

  //
  
}

//_____________________________________________________________________________
void AliTPCv3::StepManager()
{
  //
  // Called for every step in the Time Projection Chamber
  //

  //
  // parameters used for the energy loss calculations
  //
  const Float_t prim = 14.35; // number of primary collisions per 1 cm
  const Float_t poti = 20.77e-9; // first ionization potential for Ne/CO2
  const Float_t w_ion = 35.97e-9; // energy for the ion-electron pair creation 
 
 
  const Float_t big = 1.e10;

  Int_t id,copy;
  TLorentzVector pos;
  Float_t hits[4];
  Int_t vol[2];  
  TClonesArray &lhits = *fHits;
  
  vol[1]=0;
  vol[0]=0;

  //

  gMC->SetMaxStep(big);
  
  if(!gMC->IsTrackAlive()) return; // particle has disappeared
  
  Float_t charge = gMC->TrackCharge();
  
  if(TMath::Abs(charge)<=0.) return; // take only charged particles
  
  
  id=gMC->CurrentVolID(copy);
  
  // Check the sensitive volume
  
  if (id != fIdSens1) return;
  
  //
  //  charged particle is in the sensitive volume
  //
  
  if(gMC->TrackStep() > 0) {

    
    Int_t nel = (Int_t)(((gMC->Edep())-poti)/w_ion) + 1;
    nel=TMath::Min(nel,300); // 300 electrons corresponds to 10 keV
    
    gMC->TrackPosition(pos);
    hits[0]=pos[0];
    hits[1]=pos[1];
    hits[2]=pos[2];

    //
    // check the selected side of the TPC
    //
 
    if(fSide && fSide*hits[2]<=0.) return;

    hits[3]=(Float_t)nel;
    
    // Add this hit
   
    new(lhits[fNhits++]) AliTPChit(fIshunt,gAlice->CurrentTrack(),vol,hits);
    
  } 
  
  // Stemax calculation for the next step
  
  Float_t pp;
  TLorentzVector mom;
  gMC->TrackMomentum(mom);
  Float_t ptot=mom.Rho();
  Float_t beta_gamma = ptot/gMC->TrackMass();
  
  if(gMC->TrackPid() <= 3 && ptot > 0.002)
    { 
      pp = prim*1.58; // electrons above 20 MeV/c are on the plateau!
    }
  else
    {
      pp=prim*BetheBloch(beta_gamma);    
      if(TMath::Abs(charge) > 1.) pp *= (charge*charge);
    }
  
  Float_t random[1];
  gMC->Rndm(random,1); // good, old GRNDM from Geant3
  
  Double_t rnd = (Double_t)random[0];
  
  gMC->SetMaxStep(-TMath::Log(rnd)/pp);
  
}

//_____________________________________________________________________________
Float_t AliTPCv3::BetheBloch(Float_t bg)
{
  //
  // Bethe-Bloch energy loss formula
  //
  const Double_t p1=0.76176e-1;
  const Double_t p2=10.632;
  const Double_t p3=0.13279e-4;
  const Double_t p4=1.8631;
  const Double_t p5=1.9479;

  Double_t dbg = (Double_t) bg;

  Double_t beta = dbg/TMath::Sqrt(1.+dbg*dbg);

  Double_t aa = TMath::Power(beta,p4);
  Double_t bb = TMath::Power(1./dbg,p5);

  bb=TMath::Log(p3+bb);
  
  return ((Float_t)((p2-aa-bb)*p1/aa));
}
