///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Time Projection Chamber version 2 -- detailed TPC and slow simulation    //
//                                                                           //
//Begin_Html
/*
<img src="gif/AliTPCv2Class.gif">
*/
//End_Html
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TMath.h>
#include <TGeometry.h>
#include "AliTPCv2.h"
#include "AliRun.h"
#include <iostream.h>
#include <fstream.h>

#include "AliMC.h"
#include "AliConst.h"
#include <stdlib.h>

#include "AliTPCParam.h"
#include "AliTPCD.h"

ClassImp(AliTPCv2)
 
//_____________________________________________________________________________
AliTPCv2::AliTPCv2(const char *name, const char *title) :
  AliTPC(name, title) 
{
  //
  // Standard constructor for Time Projection Chamber version 2
  //
  fIdSens1=0;
  fIdSens2=0;
  SetBufferSize(128000);
}
 
//_____________________________________________________________________________
void AliTPCv2::CreateGeometry()
{
  //
  // Create the geometry of Time Projection Chamber version 2
  //
  //Begin_Html
  /*
    <img src="gif/AliTPCv2.gif">
  */
  //End_Html
  //Begin_Html
  /*
    <img src="gif/AliTPCv2Tree.gif">
  */
  //End_Html

  AliMC* pMC = AliMC::GetMC();

  Int_t *idtmed = gAlice->Idtmed();

  AliTPCParam * fTPCParam = &(fDigParam->GetParam());

  Float_t tana;
  Int_t isll;
  Float_t rlsl, wlsl, rssl, rlsu, wssl, wlsu, rssu, wssu;
  Int_t i;
  Float_t alpha, x, y, z, sec_thick;
  
  Float_t r1, r2, x1, z0, z1, x2, theta1, theta2, theta3, dm[21];
  Int_t il, iu;
  Float_t z_side, zz;
  Int_t idrotm[100];
  
  Float_t x0l, x0u;
  Int_t idr;
  Float_t thl, opl;
  Int_t ils, iss;
  Float_t thu, opu;
  Int_t ifl1 = 0, ifl2 = 0;
  Float_t phi1, phi2, phi3;
  
  // --------------------------------------------------- 
  //        sector specification check 
  // --------------------------------------------------- 
  if (fSecAL >= 0) {
    ifl1 = 0;
    
    for (i = 0; i < 6; ++i) {
      if (fSecLows[i] > 0 && fSecLows[i] <25) {
	ifl1 = 1;
	printf("*** SECTOR %d selected\n",fSecLows[i]);
      }
    }

  } else {
    printf("*** ALL LOWER SECTORS SELECTED ***\n");
    ifl1 = 1;
  }

  if (fSecAU >= 0) {
    ifl2 = 0;
    
    for (i = 0; i < 12; ++i) {
      if (fSecUps[i] > 24 && fSecUps[i] < 73) {
	ifl2 = 1;
	printf("*** SECTOR %d selected\n",fSecUps[i]);
      }
    }
    
  } else {
    printf("*** ALL UPPER SECTORS SELECTED ***\n");
    ifl1 = 1;
  }
  
  if (ifl1 == 0 && ifl2 == 0) {
    printf("*** ERROR: AT LEAST ONE SECTOR MUST BE SPECIFIED ***\n");
    printf("!!! PROGRAM STOPPED !!!\n");
    exit(1);
  }
  
  if ((fSecAL < 0 || fSecAU < 0) && fSens >= 0) {
    printf("** ERROR: STRIPS CANNOT BE SPECIFIED FOR ALL SECTORS **\n");
    printf("!!! PROGRAM STOPPED !!!\n");
    exit(1);
  }
  // ---------------------------------------------------- 
  //          FIELD CAGE WITH ENDCAPS - CARBON FIBER 
  //          THIS IS ALSO A TPC MOTHER VOLUME 
  // ---------------------------------------------------- 
  dm[0] = 76.;
  dm[1] = 278.;
  dm[2] = 275.;
  
  pMC->Gsvolu("TPC ", "TUBE", idtmed[407], dm, 3);
  // ------------------------------------------------------- 
  //     drift gas Ne/CO2 (90/10 volume) - nonsensitive 
  //     field cage thickness = 0.52% X0 
  // ---------------------------------------------------- 
  dm[0] = 76.+0.09776;
  dm[1] = 257.;
  dm[2] = 250.;
  
  pMC->Gsvolu("TGAS", "TUBE", idtmed[402], dm, 3);
  // ------------------------------------------------------ 
  //     "side" gas volume (the same as drift gas), 
  //      here the readout chambers are positioned 
  // ------------------------------------------------------ 
  dm[2] = 0.5*(275.-250.);
  z_side = dm[2];

  pMC->Gsvolu("TPSG", "TUBE", idtmed[401], dm, 3);
  // ------------------------------------------------------ 
  //      HV midplane - 20 microns of mylar 
  // ----------------------------------------------------- 
  dm[2] = .001;
  
  pMC->Gsvolu("TPHV", "TUBE", idtmed[405], dm, 3);
  
  // ==================================================== 
  //   lower and upper readout chambers 
  // ==================================================== 
  //   sectros opening angles in degrees 
  // --------------------------------------------------- 
  opl = 30.;
  opu = 15.;
  thl = TMath::Tan(opl * .5 * kDegrad);
  thu = TMath::Tan(opu * .5 * kDegrad);
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
  //     S-sectors (lower sectors) 
  // --------------------------------------------------- 
  dm[0] = wssl * .5;
  dm[1] = wssu * .5;
  dm[2] = sec_thick * .5;
  dm[3] = (rssu - rssl) * .5;
  
  x0l = rssl + dm[3];
  
  pMC->Gsvolu("TRCS", "TRD1", idtmed[399], dm, 4);
  // ----------------------------------------------------- 
  //     S-sectors --> "gas sectors" - sensitive 
  // ----------------------------------------------------- 
  dm[2] = (250.-0.001)/2.;
  pMC->Gsvolu("TSGA", "TRD1", idtmed[403], dm, 4);
  // ------------------------------------------------------------- 
  //  Only for the debugging purpose and resolution calculation 
  //  Sensitive strips at the pad-row center 
  // ------------------------------------------------------------- 
  if (fSens >= 0) {
    pMC->Gsvolu("TSST", "TRD1", idtmed[403], dm, 0);
    dm[3] = .005;
  
    z0    = rssl + (rssu - rssl) * .5;
    
    for (iss = 0; iss < fTPCParam->GetNRowLow(); ++iss) {
      r1    = fTPCParam->GetPadRowRadiiLow(iss);
      r2    = r1 + dm[3] * 2.;
      dm[0] = r1 * thl - 2.63;
      dm[1] = r2 * thl - 2.63;
      
      zz    = -z0 + r1 +dm[3];

      pMC->Gsposp("TSST", iss+1, "TSGA", 0, 0, zz, 0, "ONLY", dm, 4);
    }
    pMC->Gsord("TSGA", 3);
  }
  // --------------------------------------------------- 
  //     L-sectors (upper sectors) 
  // --------------------------------------------------- 
  dm[0] = wlsl * .5;
  dm[1] = wlsu * .5;
  dm[2] = sec_thick * .5;
  dm[3] = (rlsu - rlsl) * .5;
  
  x0u = rlsl + dm[3];
  
  pMC->Gsvolu("TRCL", "TRD1", idtmed[399], dm, 4);
  // ----------------------------------------------------- 
  //     L-sectors - "gas sectors" - sensitive! 
  // ----------------------------------------------------- 
  dm[2] = (250.-0.001)/2.;
  pMC->Gsvolu("TLGA", "TRD1", idtmed[403], dm, 4);
  // ------------------------------------------------------------- 
  //  Only for the debugging purpose and resolution calculation 
  //  Sensitive strips at the pad-row center 
  // ------------------------------------------------------------- 
  if (fSens >= 0) {
    pMC->Gsvolu("TLST", "TRD1", idtmed[403], dm, 0);
    
    dm[3] = .005;

    z0   = rlsl+ (rlsu - rlsl) * .5;
    
    for (ils = 0; ils <fTPCParam->GetNRowUp(); ++ils) {
      r1    = fTPCParam->GetPadRowRadiiUp(ils);
      r2    = r1 + dm[3] * 2.;
      dm[0] = r1 * thu - 2.63;
      dm[1] = r2 * thu - 2.63;

      zz  = -z0 + r1 +dm[3]; 

      pMC->Gsposp("TLST", ils+1, "TLGA", 0, 0, zz, 0, "ONLY", dm, 4);
    }
    pMC->Gsord("TLGA", 3);
  }
  // ****************************************************** 
  // ------------------------------------------------ 
  //      positioning of lower sectors (1-12)*2 
  //          rotation matrices 1-12 
  
  //      the isec_al flag allows to select all (<0) or 
  //      only a few (up to 6) sectors 
  // ------------------------------------------------ 
  z  = (250.+0.001)/2.;
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
      AliMatrix(idrotm[idr], theta1, phi1, theta2, phi2, theta3, phi3);
      
      alpha = (il - 1) * opl * kDegrad;
      x     = x0l * TMath::Cos(alpha);
      y     = x0l * TMath::Sin(alpha);
      
      if (fSecAL < 0) {
	// ----------------------------------------------------------- 
	//       position ALL lower sectors 
	// ----------------------------------------------------------- 
	pMC->Gspos("TSGA", il, "TGAS", x, y, z, idrotm[idr], "ONLY");
	pMC->Gspos("TSGA",il+12 , "TGAS", x, y, -z, idrotm[idr], "ONLY");
      } else {
	// ----------------------------------------------------------- 
	//       position selected lower sectors 
	// ----------------------------------------------------------- 
	for (isll = 1; isll <= 6; ++isll) {
	  if (fSecLows[isll - 1] == il) {
	    pMC->Gspos("TSGA", il, "TGAS", x, y, z, idrotm[idr], "ONLY");
	  } else if (fSecLows[isll - 1] == il + 12) {
	    pMC->Gspos("TSGA",il+12 , "TGAS", x, y, -z, idrotm[idr],"ONLY");
	  }
	}
      }
      
      pMC->Gspos("TRCS", il, "TPSG", x, y, z1, idrotm[idr], "ONLY");
      
    }
    // ---------------------------------------------------- 
    //      positioning of upper sectors (1-24)*2 
    //          rotation matrices 13-36 
    //      the isec_au flag allows to select all (<0) or 
    //      only a few (up to 12) sectors 
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

      if (fSecAU < 0) {
	// ------------------------------------------------------------- 
	//           position ALL upper sectors 
	// ------------------------------------------------------------- 
	pMC->Gspos("TLGA", iu, "TGAS", x, y, z, idrotm[idr], "ONLY");
	pMC->Gspos("TLGA",iu+24 , "TGAS", x, y, -z, idrotm[idr], "ONLY");
      } else {
	// ------------------------------------------------------------- 
	//         position selected upper sectors 
	// ------------------------------------------------------------- 
	for (isll = 1; isll <= 12; ++isll) {
	  if (fSecUps[isll - 1] == iu + 24) {
	    pMC->Gspos("TLGA", iu, "TGAS", x, y, z, idrotm[idr], "ONLY");
	  } else if (fSecUps[isll - 1] == iu + 48) {
	    pMC->Gspos("TLGA",iu+24 , "TGAS", x, y, -z, idrotm[idr],"ONLY");
	  }
	}
      }
      
      pMC->Gspos("TRCL", iu, "TPSG", x, y, z1, idrotm[idr], "ONLY");
    }
    // -------------------------------------------------------- 
    //             Spoke wheel structures 
    // -------------------------------------------------------- 
    pMC->Gsvolu("TSWS", "TUBE", idtmed[399], dm, 0);
    
    z0 = -z_side + 2.;
    
    dm[0] = 82.;
    dm[1] = 86.;
    dm[2] = 1.;
    
    pMC->Gsposp("TSWS", 1, "TPSG", 0, 0, z0, 0, "ONLY", dm, 3);
    
    dm[0] = 253.;
    dm[1] = 257.;
    
    pMC->Gsposp("TSWS", 2, "TPSG", 0, 0, z0, 0, "ONLY", dm, 3);
    
    dm[0] = 140.9;
    dm[1] = 141.9;
    
    pMC->Gsposp("TSWS", 3, "TPSG", 0, 0, z0, 0, "ONLY", dm, 3);
    
    // ------------------------------------------------------- 
    //    this volumes are to avoid overlaping 
    // ------------------------------------------------------- 
    z0 = 253.;
    
    dm[0] = 76.;
    dm[1] = 76.+0.09776;

    pMC->Gsposp("TSWS", 4, "TPC ", 0, 0, z0, 0, "ONLY", dm, 3);
    pMC->Gsposp("TSWS", 5, "TPC ", 0, 0, -z0, 0, "ONLY", dm, 3);

    z0 += 21.;

    pMC->Gsposp("TSWS", 6, "TPC ", 0, 0, z0, 0, "ONLY", dm, 3);
    pMC->Gsposp("TSWS", 7, "TPC ", 0, 0, -z0, 0, "ONLY", dm, 3);

    dm[0] = 257.;
    dm[1] = 257.+0.09776;
    dm[2] = 11.5;

    z0 = 263.5;

    pMC->Gsposp("TSWS", 8, "TPC ", 0, 0, z0, 0, "ONLY", dm, 3);
    pMC->Gsposp("TSWS", 9, "TPC ", 0, 0, -z0, 0, "ONLY", dm, 3);
    // ========================================================== 
    //                  wheels 
    // ========================================================== 
    // ---------------------------------------------------------- 
    //       Large wheel -> positioned in the TPC 
    // ---------------------------------------------------------- 
    dm[0] = 257.+0.09776;
    dm[1] = 278.;
    dm[2] = 11.5;
    pMC->Gsvolu("TPW1", "TUBE", idtmed[399], dm, 3);
    
    dm[0] = 259.;
    dm[1] = 278.;
    dm[2] = 9.5;
    
    pMC->Gsvolu("TPW2", "TUBE", idtmed[498], dm, 3);
    
    pMC->Gspos("TPW2", 1, "TPW1", 0, 0, 0, 0, "ONLY");

    pMC->Gspos("TPW1", 1, "TPC ", 0, 0, z0, 0, "ONLY");
    pMC->Gspos("TPW1", 2, "TPC ", 0, 0, -z0, 0, "ONLY");
    // ----------------------------------------------------------- 
    //     Small wheel -> positioned in the TPSG 
    // ----------------------------------------------------------- 
    dm[0] = 76.+0.09776;
    dm[1] = 82.;
    dm[2] = 11.5;
    
    pMC->Gsvolu("TPW3", "TUBE", idtmed[399], dm, 3);
    
    dm[0] = 76.+0.09776;
    dm[1] = 80.;
    dm[2] = 9.5;
    
    pMC->Gsvolu("TPW4", "TUBE", idtmed[401], dm, 3);
    
    pMC->Gspos("TPW4", 1, "TPW3", 0, 0, 0, 0, "ONLY");
    
    z0 = 1.;
    
    pMC->Gspos("TPW3", 1, "TPSG", 0, 0, z0, 0, "ONLY");
    // --------------------------------------------------------- 
    //       spokes, inner and outer, also the inner ring 
    // --------------------------------------------------------- 
    dm[0] = 0.5*(135.9-82.1);
    dm[1] = 3.;
    dm[2] = 2.;
    
    x1 = dm[0] + 82.;
    
    pMC->Gsvolu("TSPI", "BOX ", idtmed[399], dm, 3);
    
    dm[1] = 2.;
    dm[2] = 1.;
    
    pMC->Gsvolu("TSP1", "BOX ", idtmed[498], dm, 3);
    
    pMC->Gspos("TSP1", 1, "TSPI", 0, 0, 0, 0, "ONLY");
    
    dm[0] = 0.5*(256.9-142.1);
    dm[1] = 3.;
    dm[2] = 2.;
    
    x2 = dm[0] + 142.;
    
    pMC->Gsvolu("TSPO", "BOX ", idtmed[399], dm, 3);
    
    dm[1] = 2.;
    dm[2] = 1.;
    
    pMC->Gsvolu("TSP2", "BOX ", idtmed[498], dm, 3);

    pMC->Gspos("TSP2", 1, "TSPO", 0, 0, 0, 0, "ONLY");
    // -------------------------------------------------------- 
    dm[0] = 136.;
    dm[1] = 142.;
    dm[2] = 2.;

    pMC->Gsvolu("TSWH", "TUBE", idtmed[399], dm, 3);

    dm[0] = 137.;
    dm[1] = 141.;
    dm[2] = 1.;

    pMC->Gsvolu("TSW1", "TUBE", idtmed[498], dm, 3);

    pMC->Gspos("TSW1", 1, "TSWH", 0, 0, 0, 0, "ONLY");

    z0 = z_side - .16168 - 2.;
    // -------------------------------------------------------- 
    pMC->Gspos("TSWH", 1, "TPSG", 0, 0, z0, 0, "ONLY");
    // ------------------------------------------------------- 
    //     posiioning of the inner spokes 
    // ------------------------------------------------------- 
    for (il = 1; il <= 6; ++il) {
      phi1   = opl * .5 + (il - 1) * 2. * opl;
      theta1 = 90.;
      phi2   = opl * .5 + 90. + (il - 1) * 2. * opl;
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
      pMC->Gspos("TSPI", il, "TPSG", x, y, z0, idrotm[idr], "ONLY");
      
    }
    
    for (iu = 1; iu <= 12; ++iu) {
      phi1   = opu * .5 + (iu - 1) * 2. * opu;
      theta1 = 90.;
      phi2   = opu * .5 + 90. + (iu - 1) * 2. * opu;
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
      pMC->Gspos("TSPO", iu, "TPSG", x, y, z0, idrotm[idr], "ONLY");
    }
    // -------------------------------------------------------- 
    //       endcap cover (C, 0.86% X0) 
    // -------------------------------------------------------- 
    dm[0] = 76.+0.09776;
    dm[1] = 257.;
    dm[2] = 0.16168*0.5;
    
    pMC->Gsvolu("TCOV", "TUBE", idtmed[407], dm, 3);
    
    z0 = z_side - dm[2];
    
    pMC->Gspos("TCOV", 1, "TPSG", 0, 0, z0, 0, "ONLY");
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
    
    pMC->Gspos("TPSG", 1, "TPC ", 0, 0, z0, 0, "ONLY");
    pMC->Gspos("TPSG", 2, "TPC ", 0, 0, -z0, idrotm[55], "ONLY");
    // --------------------------------------------------------- 
    //     outer gas insulation (CO2) 
    // --------------------------------------------------------- 
    dm[0] = 257.+0.09776;
    dm[1] = 278.-0.25004;
    dm[2] = 275.-23.;
    
    pMC->Gsvolu("TPOI", "TUBE", idtmed[406], dm, 3);
    
    pMC->Gspos("TPHV", 1, "TGAS", 0, 0, 0, 0, "ONLY");
    pMC->Gspos("TGAS", 1, "TPC ", 0, 0, 0, 0, "ONLY");
    pMC->Gspos("TPOI", 1, "TPC ", 0, 0, 0, 0, "ONLY");
    
    pMC->Gspos("TPC ", 1, "ALIC", 0, 0, 0, 0, "ONLY");
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
    
    pMC->Gsposp("TSWS", 10, "TPC ", 0, 0, z0, 0, "ONLY", dm, 3);
    pMC->Gsposp("TSWS", 11, "TPC ", 0, 0, -z0, 0, "ONLY", dm, 3);

    dm[0] = 70.;
    
    z0 += 21.;
    
    pMC->Gsposp("TSWS", 12, "TPC ", 0, 0, z0, 0, "ONLY", dm, 3);
    pMC->Gsposp("TSWS", 13, "TPC ", 0, 0, -z0, 0, "ONLY", dm, 3);
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

    dm[9]  = 64.5;
    dm[10] = 50.;
    dm[11] = 76.;

    dm[12] = 250.;
    dm[13] = 75.;
    dm[14] = 76.;

    pMC->Gsvolu("TPIV", "PCON", idtmed[407], dm, 15);
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
    
    dm[9] = -64.5;
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
    
    pMC->Gsvolu("TPVD", "PCON", idtmed[406], dm, 21);
    
    pMC->Gspos("TPVD", 1, "TPIV", 0, 0, 0, 0, "ONLY");
    
    pMC->Gspos("TPIV", 1, "ALIC", 0, 0, 0, 0, "ONLY");
    // --------------------------------------------------- 
    //               volumes ordering 
    // --------------------------------------------------- 
    pMC->Gsord("TGAS", 6);
    pMC->Gsord("TPSG", 6);
}
 
//_____________________________________________________________________________
void AliTPCv2::DrawDetector()
{
  //
  // Draw a shaded view of the Time Projection Chamber version 1
  //

  AliMC* pMC = AliMC::GetMC();

  // Set everything unseen
  pMC->Gsatt("*", "seen", -1);
  // 
  // Set ALIC mother transparent
  pMC->Gsatt("ALIC","SEEN",0);
  //
  // Set the volumes visible
  pMC->Gsatt("TPC","SEEN",0);
  pMC->Gsatt("TGAS","SEEN",0);
  pMC->Gsatt("TPSG","SEEN",0);
  pMC->Gsatt("TPHV","SEEN",1);
  pMC->Gsatt("TRCS","SEEN",1);
  pMC->Gsatt("TRCL","SEEN",1);
  pMC->Gsatt("TSWS","SEEN",1);
  pMC->Gsatt("TPW1","SEEN",1);
  pMC->Gsatt("TPW3","SEEN",1);
  pMC->Gsatt("TSPI","SEEN",1);
  pMC->Gsatt("TSPO","SEEN",1);
  pMC->Gsatt("TSWH","SEEN",1);
  pMC->Gsatt("TPOI","SEEN",1);
  pMC->Gsatt("TPIV","SEEN",1);
  pMC->Gsatt("TPVD","SEEN",1);
  //
  pMC->Gdopt("hide", "on");
  pMC->Gdopt("shad", "on");
  pMC->Gsatt("*", "fill", 7);
  pMC->SetClipBox(".");
  pMC->SetClipBox("*", 0, 1000, -1000, 1000, -1000, 1000);
  pMC->DefaultRange();
  pMC->Gdraw("alic", 40, 30, 0, 12, 9.5, .025, .025);
  pMC->Gdhead(1111, "Time Projection Chamber");
  pMC->Gdman(18, 4, "MAN");
  pMC->Gdopt("hide","off");
}

//_____________________________________________________________________________
void AliTPCv2::CreateMaterials()
{
  //
  // Define materials for version 2 of the Time Projection Chamber
  //

  AliMC* pMC = AliMC::GetMC();

  //
  // Increase maximum number of steps
  pMC->SetMaxNStep(30000);
  //
  AliTPC::CreateMaterials();
}

//_____________________________________________________________________________
void AliTPCv2::Init()
{
  //
  // Initialises version 2 of the TPC after that it has been built
  //
  AliMC* pMC=AliMC::GetMC();
  Int_t *idtmed = gAlice->Idtmed();
  AliTPC::Init();
  fIdSens1=pMC->VolId("TLGA"); // L-sector
  fIdSens2=pMC->VolId("TSGA"); // S-sector 
  fIdSens3=pMC->VolId("TSST"); // strip - S-sector (not always used)
  fIdSens4=pMC->VolId("TLST"); // strip - S-sector (not always used)

  pMC->SetMaxNStep(30000); // max. number of steps increased

  pMC->Gstpar(idtmed[403],"LOSS",5);

  printf("*** TPC version 2 initialized ***\n");
  printf("Maximum number of steps = %d\n",pMC->GetMaxNStep());

  //
  
}

//_____________________________________________________________________________
void AliTPCv2::StepManager()
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
  Float_t hits[4];
  Int_t vol[2];  
  TClonesArray &lhits = *fHits;
  AliMC* pMC=AliMC::GetMC();
  
  vol[1]=0;

  //

  pMC->SetMaxStep(big);
  
  if(!pMC->TrackAlive()) return; // particle has disappeared
  
  Float_t charge = pMC->TrackCharge();
  
  if(TMath::Abs(charge)<=0.) return; // take only charged particles
  
  
  id=pMC->CurrentVol(0, copy);
  
  // Check the sensitive volume
  
  if(id == fIdSens1) 
    {
      vol[0] = copy + 24; // L-sector number
    }
  else if(id == fIdSens2) 
    {
      vol[0] = copy; // S-sector number 
    }
  else if(id == fIdSens3 && pMC->TrackEntering())
    {
      vol[1] = copy;  // row number  
      id = pMC->CurrentVolOff(1,0,copy);
      vol[0] = copy; // sector number (S-sector)
      
      pMC->TrackPosition(hits);
      hits[3]=0.; // this hit has no energy loss
      new(lhits[fNhits++]) AliTPChit(fIshunt,gAlice->CurrentTrack(),vol,hits);
    }
  else if(id == fIdSens4 && pMC->TrackEntering())
    {
      vol[1] = copy; // row number 
      id = pMC->CurrentVolOff(1,0,copy);
      vol[0] = copy+24; // sector number (L-sector)
      
      pMC->TrackPosition(hits);
      hits[3]=0.; // this hit has no energy loss
      new(lhits[fNhits++]) AliTPChit(fIshunt,gAlice->CurrentTrack(),vol,hits);
    }
  else return;
  
  //
  //  charged particle is in the sensitive volume
  //
  
  if(pMC->TrackStep() > 0) {
    
    Int_t nel = (Int_t)(((pMC->Edep())-poti)/w_ion) + 1;
    nel=TMath::Min(nel,300); // 300 electrons corresponds to 10 keV
    
    pMC->TrackPosition(hits);
    hits[3]=(Float_t)nel;
    
    // Add this hit
    
    new(lhits[fNhits++]) AliTPChit(fIshunt,gAlice->CurrentTrack(),vol,hits);
    
  } 
  
  // Stemax calculation for the next step
  
  Float_t pp;
  Float_t vect[4];
  pMC->TrackMomentum(vect);
  Float_t ptot = vect[3];
  Float_t beta_gamma = ptot/(pMC->TrackMass());
  
  if(pMC->TrackPid() <= 3 && ptot > 0.002)
    { 
      pp = prim*1.58; // electrons above 20 MeV/c are on the plateau!
    }
  else
    {
      pp=prim*BetheBloch(beta_gamma);    
      if(TMath::Abs(charge) > 1.) pp *= (charge*charge);
    }
  
  Float_t random[1];
  pMC->Rndm(random,1); // good, old GRNDM from Geant3
  
  Double_t rnd = (Double_t)random[0];
  
  pMC->SetMaxStep(-TMath::Log(rnd)/pp);
  
}

//_____________________________________________________________________________
Float_t AliTPCv2::BetheBloch(Float_t bg)
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
