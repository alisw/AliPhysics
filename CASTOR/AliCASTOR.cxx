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
Revision 1.7  2000/01/19 17:16:41  fca
Introducing a list of lists of hits -- more hits allowed for detector now

Revision 1.6  1999/09/29 09:24:07  fca
Introduction of the Copyright and cvs Log

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  CASTOR                                                                   //
//  This class contains the description of the CASTOR detector               //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliCASTORClass.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>The responsible person for this module is
<a href="mailto:aris.angelis@cern.ch">Aris Angelis</a>.
</font>
<pre>
*/
//End_Html
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include "AliCASTOR.h"
#include <TNode.h>
#include <TPGON.h>
#include "TGeometry.h"
#include "AliMagF.h"
#include "AliRun.h"
#include "AliMC.h"
#include "AliConst.h"

ClassImp(AliCASTOR)
 
//_____________________________________________________________________________
AliCASTOR::AliCASTOR()
{
  //
  // Default constructor for CASTOR
  //
  fIshunt   = 0;
}
 
//_____________________________________________________________________________
AliCASTOR::AliCASTOR(const char *name, const char *title)
       : AliDetector(name,title)
{
  //
  // Standard constructor for CASTOR
  //

  //
  // Create a tree of castor hits
  fHits   = new TClonesArray("AliCASTORhit",  405);
  gAlice->AddHitList(fHits);
  
  fIshunt     =  0;
   
  SetMarkerColor(7);
  SetMarkerStyle(2);
  SetMarkerSize(0.4);
}
 
//_____________________________________________________________________________
void AliCASTOR::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{
  //
  // Add a CASTOR hit
  //
  TClonesArray &lhits = *fHits;
  new(lhits[fNhits++]) AliCASTORhit(fIshunt,track,vol,hits);
}

//_____________________________________________________________________________
void AliCASTOR::BuildGeometry()
{
  //
  // Build CASTOR ROOT TNode geometry for event display
  TNode *Node, *Top;
  TPGON *pgon;
  const int kColorCASTOR  = 4;
  //
  Top=gAlice->GetGeometry()->GetNode("alice");
  
  // CASTOR
  pgon = new TPGON("S_CASTOR","S_CASTOR","void",22.5,360,8,2);
  pgon->DefineSection(0,-69.05885,2.598121,12.86874);
  pgon->DefineSection(1,69.05885,2.787778,13.88912);
  new TRotMatrix("rotcas","rotcas",90,180,90,90,180,0);

  Top->cd();
  Node = new TNode("CASTOR","CASTOR","S_CASTOR",0,0,-1809.59,"rotcas");
  Node->SetLineColor(kColorCASTOR);
  fNodes->Add(Node);
}

//_____________________________________________________________________________
Int_t AliCASTOR::DistancetoPrimitive(Int_t , Int_t )
{
   return 9999;
}
 
 
ClassImp(AliCASTORv1)
 
//_____________________________________________________________________________
AliCASTORv1::AliCASTORv1() : AliCASTOR()
{
  //
  // Default constructor for CASTOR version 1
  //
  fOdFiber = 0;
  fOdCladding = 0;
  fOdAbsorber = 0;
  fOctants = 0;
  fLayersEM = 0;
  fLayersHad = 0;
  fPhiOct = 0;
  fRadCore = 0;
  fRadFactor = 0;
}
 
//_____________________________________________________________________________
AliCASTORv1::AliCASTORv1(const char *name, const char *title)
       : AliCASTOR(name,title)
{
  //
  // Standard constructor for CASTOR version 1
  //
  fOdFiber = 0;
  fOdCladding = 0;
  fOdAbsorber = 0;
  fOctants = 0;
  fLayersEM = 0;
  fLayersHad = 0;
  fPhiOct = 0;
  fRadCore = 0;
  fRadFactor = 0;
}
 
//_____________________________________________________________________________
void AliCASTORv1::CreateGeometry()
{
  //
  // Creation of the geometry of the CASTOR detector
  //
  //Begin_Html
  /*
    <img src="picts/AliCASTORTree.gif">
  */
  //End_Html
  //Begin_Html
  /*
    <img src="picts/AliCASTOR.gif">
  */
  //End_Html
  //
  //   28 March 1997   23 February 1998              Aris L. S. Angelis   * 
  // >--------------------------------------------------------------------<* 
  
  
  Float_t dhad[11], dcal[3], beta, doct[11], alfa1, alfa2, fact1, fact2,fact3;
  Float_t dclha[3], dcoha[3], dclem[3], dbxha[3], dcoem[3], dcalt[5], dcalv[5], dbxem[3];
  Float_t rzhig;
  Float_t s1, s2, s3, rxyin, rzlow, rxyut, facemd, facein, dlayha, dlayem, doctem, doctha, faceut, zendha, phicov;
  Float_t doctnt;
  Float_t zemhad;
  Int_t idrotm[100];
  Float_t thecen, xp, xxmdhi, zp, yp, rinbeg;
  Float_t rutbeg, xxinhi, rinend, rutend, xxmdlo;
  Float_t dztotl, xxinlo, xxuthi;
  Float_t xxutlo, dem[11], ang;
  Int_t nfx;
  Float_t rxy;
  // Angle (deg) of inclination of quartz fibres w.r.t. to beam (Cerenkov angle).
  const Float_t kBetaD = 45;
  //Rapidity range covered by the calorimeter.
  const Float_t kEtaLow  = 5.6;
  const Float_t kEtaHigh = 7.2;
  // Z position (cm) of beginning of calorimeter EM section (the tip.
  const Float_t kZbegem = 1740;
  // Number of azimuthal calorimeter sectors: octants.
  fOctants = 8;
  // Number of e-m and hadronic layers (each layer comprises a slice
  // of absorber material followed by a slice of active quartz fibres).
  //     DATA NLAYEM,NLAYHA /9,69/  ! 0.64 + 9.73 lambda_i
  fLayersEM  = 8;
  fLayersHad = 72;  // 0.57 + 10.15 lambda_i
  // Number of planes of quartz fibres within each active slice for
  // e-m and hadronic sections.
  const Int_t kFibersEM  = 2;
  const Int_t kFibersHad = 4;
  // Thickness (cm) of absorber material for e-m and hadronic layers.
  const Float_t kAbsorberEM  = 0.5;
  const Float_t kAbsorberHad = 1;
  // Diameter (cm) of fibre core and of fibre with cladding.
  const Float_t kDiamCore     = 0.043;
  const Float_t kDiamCladding = 0.045;

  Int_t i;
  static Int_t debugFlag = 0;
  
  Int_t *idtmed = fIdtmed->GetArray()-1499;

  
  // >--------------------------------------------------------------------<*
  // **> Note: ALICE frame XYZ, proper ref. frame of a trapezoid X'Y'Z'. 
  // --- Common which contains debug flags for the various detectors --- 
  // --- Also control flags (JPAWF,JOUTF) for each detector added --- 
  
  // **> Common containing some of the Castor FCAL geometry data. 
  
  //**> Angle (deg) of inclination of quartz fibres w.r.t. to beam
  //**> (Cerenkovangle).
  // **> Rapidity range covered by the calorimeter. 
  // **> Z position (cm) of beginning of calorimeter EM section (the tip. 
  // **> Number of planes of quartz fibres within each active slice for 
  // **> e-m and hadronic sections. 
  // **> Thickness (cm) of absorber material for e-m and hadronic layers. 
  // **> Diameter (cm) of fibre core and of fibre with cladding. 
  // **> E-M and hadronic sections of an octant and complete octant module 
  // **> (general trapezoids). 
  // **> Imaginary box to hold the complete calorimeter. 
  // **> Imaginary rectangular boxes containing the trapezoids of the 
  // **> EM and Hadronic sections of an Octant. 
  // **> Cylindrical volumes for clad fibres and fibre cores in the 
  // **> EM and Had sections. 
  //**> Narrow stainless steel conical beam tube traversing the calorimeter.
  // **> Print calorimeter parameters. 
  // **> Number of azimuthal calorimeter sectors: octants. 
  //      DATA NOCTS / 16 / 
  // **> Number of e-m and hadronic layers (each layer comprises a slice 
  // **> of absorber material followed by a slice of active quartz fibres). 
  //      DATA NLAYEM,NLAYHA /9,69/  ! 0.64 + 9.73 lambda_i 
  // 0.57 + 10.15 lambda_i 
  if (debugFlag > 0) {
    printf("----------------------------------\n");
    printf(" EtaLo = %f, EtaHigh = %f, ZbegEM =%f\n",kEtaLow, kEtaHigh,kZbegem);
    printf(" Nocts =%d, NlayEM=%d, NlayHad = %d\n",fOctants,fLayersEM,fLayersHad);
    printf("----------------------------------\n");
  }
  // **> Radius of sensitive fibre core. 
  fRadCore = kDiamCore/2;
  // **> Radius normalised to radius of 0.5 mm used in the calculation of 
  // **> the Cherenkov tables. 
  fRadFactor = fRadCore / .05;
  // **> Total number of sensitive QF plane layers. 
  //nqemly = fLayersEM*kFibersEM;
  //nqhaly = fLayersHad*kFibersHad;
  beta   = kBetaD*kDegrad; // **> Conversions to radians. 
  // **> Thickness of e-m and hadronic layers: 
  // **> Thickness = Thickness_of_Absorber + Thickness_of_N_Fibre_Planes 
  // **> For N pair: Thickness_of_N_Fibre_Planes = N/2 * [2+TMath::Sqrt(3)]*R_fibre
  // **> taking into account staggering of fibres in adjacent planes. 
  //**> For simplicity staggering not yet introduced, use TMath::Sqrt(4) temporarily.
  dlayem = kAbsorberEM +(0.5*kFibersEM )*(2+TMath::Sqrt(4.))*kDiamCladding/2;
  dlayha = kAbsorberHad+(0.5*kFibersHad)*(2+TMath::Sqrt(4.))*kDiamCladding/2;
  if (debugFlag > 0) {
    printf(" Layer Thickness. EM = %f, Had = %f\n",dlayem,dlayha);
  }
  // **> Thickness of complete octant, along the line perpendicular 
  // **> to the layers. 
  // **> Thickness = NlayerEM*DlayerEM + NlayerHad*DlayerHad (DeltaZ'). 
  doctem = fLayersEM*dlayem;
  doctha = fLayersHad*dlayha;
  doctnt = doctem + doctha;
  if (debugFlag > 0) {
    printf(" Octant Thickness. EM = %f, Had = %f, Total = %f\n",doctem,doctha,doctnt);
  }
  // **> Construct one octant module: general trapezoid, rotated such 
  // **> that the fibre planes are perpenicular to the Z axis of the 
  // **> proper reference frame (X'Y'Z' frame). 
  // **> Calculation of the length of the faces at +/- DeltaZ'/2 of an 
  // **> octant, projected onto the Y'Z' plane (see notes dated 4/4/97). 
  alfa1 = TMath::ATan(exp(-kEtaLow)) * 2.;
  alfa2 = TMath::ATan(exp(-kEtaHigh)) * 2.;
  fact1 = (TMath::Tan(alfa1) - TMath::Tan(alfa2)) * TMath::Cos(alfa1) / TMath::Sin(beta - alfa1);
  if (debugFlag > 0) {
    printf(" Beta =%f,Fact1 =%f\n",kBetaD, fact1);
    printf(" EtaLow=%f, EtaHigh=%f, Alfa1=%f, Alfa2=%f\n",kEtaLow,kEtaHigh,alfa1*kRaddeg,alfa2*kRaddeg);
  }
  // **> Face at entrance to E-M section (-DeltaZ'/2). 
  facein = fact1 * kZbegem;
  // **> Face at interface from E-M to Hadronic section. 
  facemd = (doctem / TMath::Sin(beta) + kZbegem) * fact1;
  // **> Face at exit of Hadronic section (+DeltaZ'/2). 
  faceut = (doctnt / TMath::Sin(beta) + kZbegem) * fact1;
  if (debugFlag > 0) {
    printf(" Octant Face Length. Front: %f, Back: %f, EM-Had: %f\n",facein,faceut,facemd);
  }
  // **> Angular coverage of octant (360./8) projected onto plane 
  // **> tilted at angle Beta (see notes dated 28/3/97). 
  //**> PhiTilted = 2*atan[TMath::Tan(phi/2)TMath::Cos(beta)] = 32.65 deg for beta=45,phi=22.5.
  fPhiOct = k2PI / fOctants;
  phicov = TMath::ATan(TMath::Tan(fPhiOct / 2.) * TMath::Cos(beta)) * 2.;
  if (debugFlag > 0) {
    printf(" FPhiOct =%f, PhiCov =%f\n",fPhiOct * kRaddeg,phicov * kRaddeg);
  }
  // **> Dimensions along X' of front and back faces of calorimeter 
  // **> (see notes dated 8/4/97). 
  fact2  = TMath::Tan(alfa2) / TMath::Sin(beta);
  fact3  = TMath::Cos(alfa2) / TMath::Sin(beta - alfa2);
  zendha = doctnt * fact3 + kZbegem;
  zemhad = doctem * fact3 + kZbegem;
  if (debugFlag > 0) {
    printf(" ZbegEM =%f, ZendHA =%f, ZEMHad =%f\n",kZbegem,zendha, zemhad);
    printf(" Fact2 =%f, Fact3 =%f\n",fact2,fact3);
  }
  // **> DeltaX' at -DeltaY'/2, -DeltaZ'/2. 
  xxinlo = fact2 * 2*kZbegem * TMath::Tan(phicov / 2.);
  // **> DeltaX' at +DeltaY'/2, -DeltaZ'/2. 
  xxinhi = (fact2 + fact1) * 2*kZbegem * TMath::Tan(phicov / 2.);
  // **> DeltaX' at -DeltaY'/2, +DeltaZ'/2. 
  xxutlo = zendha * 2. * fact2 * TMath::Tan(phicov / 2.);
  // **> DeltaX' at +DeltaY'/2, +DeltaZ'/2. 
  xxuthi = zendha * 2. * (fact2 + fact1) * TMath::Tan(phicov / 2.);
  // **> DeltaX' at -DeltaY'/2, at EM/Had interface. 
  xxmdlo = zemhad * 2. * fact2 * TMath::Tan(phicov / 2.);
  // **> DeltaX' at +DeltaY'/2, at EM/Had interface. 
  xxmdhi = zemhad * 2. * (fact2 + fact1) * TMath::Tan(phicov / 2.);
  if (debugFlag > 0) {
    printf(" XXinLo=%f, XXinHi=%f, XXutLo=%f, XXutHi=%f, XXmdLo=%f, XXmdHi=%f\n",
           xxinlo,xxinhi,xxutlo,xxuthi,xxmdlo,xxmdhi);
  }
  //**> Calculate the polar angle in the X'Y'Z' frame of the line joining the
  //**> centres of the front and back faces of the octant (see notes dated 9/4/97).
  s1  = (1. - fact2 * TMath::Cos(beta)) * kZbegem;
  s2  = (fact2 + fact1 / 2.) * kZbegem;
  s3  = TMath::Sqrt(s1 * s1 + s2 * s2 - s1 * s2 * TMath::Cos(kPI - beta));
  ang = TMath::ASin(sin(kPI - beta) * s2 / s3);
  thecen = kPI/2 - beta + ang;
  if (debugFlag > 0) {
    printf(" S1=%f, S2=%f, S3=%f, Ang=%f, TheCen=%f\n",s1,s2,s3,ang*kRaddeg,thecen*kRaddeg);
  }
  // **> Construct the octant volume. 
  doct[0] = 180*0.125;
  doct[1] = 360.;
  doct[2] = 8.;
  doct[3] = 2.;
  doct[4] = -(zendha - kZbegem + faceut * TMath::Cos(beta)) / 2.;
  doct[5] = TMath::Tan(alfa2) * kZbegem;
  doct[6] = TMath::Tan(alfa1) * kZbegem;
  doct[7] = (zendha - kZbegem + faceut * TMath::Cos(beta)) / 2.;
  doct[8] = zendha * TMath::Tan(alfa2);
  doct[9] = (faceut + zendha * fact2) * TMath::Sin(beta);
  
  if (debugFlag > 0) {
    printf("\n Doct(1-10) = ");
    for (i = 1; i <= 10; ++i) {
      printf("%f, ",doct[i - 1]);
    }
    printf("   \n");
  }
  gMC->Gsvolu("OCTA", "PGON", idtmed[fOdAbsorber - 1], doct, 10);
  gMC->Gsdvn("OCT ", "OCTA", 8, 2);
  // absorber material. 
  // **> Construct the E-M section volume. 
  dem[0]  = doctem / 2.;      // DeltaZ'/2 
  dem[1]  = thecen *kRaddeg;  // Theta[(Centre(-DeltaZ')--Centre(+DeltaZ' 
  dem[2]  = 90.;              // Phi[(Centre(-DeltaZ')--Centre(+DeltaZ')] 
  dem[3]  = facein / 2.;      // DeltaY'/2 at -DeltaZ'/2. 
  dem[4]  = xxinlo / 2.;      // DeltaX'/2 at -DeltaY'/2 at -DeltaZ'/2. 
  dem[5]  = xxinhi / 2.;      // DeltaX'/2 at +DeltaY'/2 at -DeltaZ'/2. 
  dem[6]  = 0.;               // Angle w.r.t. Y axis of line joining cent 
                                // at +/- DeltaY at -DeltaZ. // Angle w.r.t. Y axis of line joining cent 
  dem[7]  = facemd / 2.;      // DeltaY'/2 at +DeltaZ'. 
  dem[8]  = xxmdlo / 2.;      // DeltaX'/2 at -DeltaY'/2 at +DeltaZ'/2. 
  dem[9]  = xxmdhi / 2.;      // DeltaX'/2 at +DeltaY'/2 at +DeltaZ'/2. 
  dem[10] = 0.;               // Angle w.r.t. Y axis of line joining cent
                                // at +/- DeltaY at +DeltaZ. 
  
  if (debugFlag > 0) {
    printf("\n De-m(1-11) =");
    for (i = 1; i <= 11; ++i) {
      printf("%f, ",dem[i - 1]);
    }
    printf("   \n");
  }
  gMC->Gsvolu("EM  ", "TRAP", idtmed[fOdAbsorber - 1], dem, 11);
  // absorber material. 
  // **> Construct the Hadronic section volume. 
  // Fill with s 
  dhad[0]  = doctha / 2.;      // DeltaZ'/2 
  dhad[1]  = thecen *kRaddeg;  // Theta[(Centre(-DeltaZ')--Centre(+DeltaZ' 
  dhad[2]  = 90.;              // Phi[(Centre(-DeltaZ')--Centre(+DeltaZ')] 
  dhad[3]  = facemd / 2.;      // DeltaY'/2 at -DeltaZ'/2. 
  dhad[4]  = xxmdlo / 2.;      // DeltaX'/2 at -DeltaY'/2 at -DeltaZ'/2. 
  dhad[5]  = xxmdhi / 2.;      // DeltaX'/2 at +DeltaY'/2 at -DeltaZ'/2. 
  dhad[6]  = 0.;               // Angle w.r.t. Y axis of line joining cent
  // at +/- DeltaY at -DeltaZ. 
  dhad[7]  = faceut / 2.;      // DeltaY'/2 at +DeltaZ'. 
  dhad[8]  = xxutlo / 2.;      // DeltaX'/2 at -DeltaY'/2 at +DeltaZ'/2. 
  dhad[9]  = xxuthi / 2.;      // DeltaX'/2 at +DeltaY'/2 at +DeltaZ'/2. 
  dhad[10] = 0.;               // Angle w.r.t. Y axis of line joining cent
  // at +/- DeltaY at +DeltaZ. 
  
  if (debugFlag > 0) {
    printf("\n Dhad(1-11) = ");
    for (i = 1; i <= 11; ++i) {
      printf("%f, ",dhad[i - 1]);
    }
    printf("   \n");
  }
  gMC->Gsvolu("HAD ", "TRAP", idtmed[fOdAbsorber - 1], dhad, 11); // absorber material. 
  // **> Rotation matrix to rotate fibres verticaly to fit into holes. 
  // Fill with 
  AliMatrix(idrotm[0], 90., 0., 180., 0., 90., 90.);
  // **> Internal structure of the EM section starts here.  <--- 
  // **> Construct one sampling module 
  gMC->Gsdvn("SLEM", "EM  ", fLayersEM, 3);
  gMC->Gsatt("SLEM", "SEEN", 0);
  // **> Construct the (imaginary) rectangular box embedding the fibres 
  // **> Fill with air, make it invisible on the drawings. 
  dbxem[0] = xxmdhi / 2.;
  dbxem[2] = kFibersEM*kDiamCladding/2;
  dbxem[1] = facemd / 2. + dbxem[2] * TMath::Tan(thecen);
  if (debugFlag > 0) {
    printf(" DbxEM(1-3) =");
    for (i = 1; i <= 3; ++i) {
      printf("%f, ",dbxem[i - 1]);
    }
    printf("   \n");
  }
  gMC->Gsvolu("BXEM", "BOX ", idtmed[1501], dbxem, 3);
  gMC->Gsatt("BXEM", "SEEN", 0);
  // **> Divide along Z to obtain one layer 
  gMC->Gsdvn("RWEM", "BXEM", 2, 3);
  gMC->Gsatt("RWEM", "SEEN", 0);
  // **> Divide along X' to accomodate the maximum number of individual 
  //**> fibres packed along X', make the divisions invisible on the drawings.
  nfx = Int_t(xxmdhi / .045);
  if (debugFlag > 0) {
    printf(" NfxEM = %d\n",nfx);
  }
  gMC->Gsdvn("FXEM", "RWEM", nfx, 1);
  gMC->Gsatt("FXEM", "SEEN", 0);
  // **> Construct the fiber cladding 
  dclem[0] = 0.;
  dclem[1] = kDiamCladding/2;
  dclem[2] = dbxem[1];
  if (debugFlag > 0) {
    printf(" DclEM(1-3) = \n");
    for (i = 1; i <= 3; ++i) {
      printf("%f, ",dclem[i - 1]);
    }
    printf("   \n");
  }
  gMC->Gsvolu("CLEM", "TUBE", idtmed[fOdCladding - 1], dclem,3);
  gMC->Gsatt("CLEM", "SEEN", 0);
  //**> Construct the cylindrical volume for a fibre core in the EM section.
  //**> Fill with selected fibre material, make it invisible on the drawings.
  dcoem[0] = 0.;
  dcoem[1] = kDiamCore/2;
  dcoem[2] = dbxem[1];
  if (debugFlag > 0) {
    printf(" DcoEM(1-3) = ");
    for (i = 1; i <= 3; ++i) {
      printf("%f, ",dcoem[i - 1]);
    }
    printf("   \n");
  }
  gMC->Gsvolu("COEM", "TUBE", idtmed[fOdFiber - 1], dcoem,3);
  gMC->Gsatt("COEM", "SEEN", 0);
  // **> Position the volumes 
  // **> Put the air section inside one sampling module 
  // **> Use MANY to obtain clipping of protruding edges. 
  xp = 0.;
  zp = dlayem / 2. - 0.5*kFibersEM*kDiamCladding;
  yp = zp * TMath::Tan(thecen);
  gMC->Gspos("BXEM", 1, "SLEM", xp, yp, zp, 0, "MANY");
  // **> Place the core fibre in the clad 
  xp = 0.;
  yp = 0.;
  zp = 0.;
  gMC->Gspos("COEM", 1, "CLEM", xp, yp, zp, 0, "MANY");
  // **> Put the fiber in its air box 
  gMC->Gspos("CLEM", 1, "FXEM", xp, yp, zp, idrotm[0], "MANY");
  // **> Internal structure of the Hadronic section starts here.  <--- 
  gMC->Gsdvn("SLHA", "HAD ", fLayersHad, 3);
  gMC->Gsatt("SLHA", "SEEN", 0);
  // **> Construct the air section where the fibers are 
  dhad[0] = 0.5*kFibersEM*kDiamCladding;
  gMC->Gsvolu("AIHA", "TRAP", idtmed[1501], dhad, 11);
  // **> Divide along z in the appropriate number of layers 
  gMC->Gsdvn("SAHA", "AIHA", 4, 3);
  //**> Construct the (imaginary) rectangular box embedding one lauer of fibres
  // **> Fill with air, make it invisible on the drawings. 
  dbxha[0] = xxuthi / 2.;
  dbxha[2] = 0.5*kFibersHad*kDiamCladding;
  dbxha[1] = faceut / 2. + dbxha[2] * TMath::Tan(thecen);
  if (debugFlag > 0) {
    printf(" DbxHa(1-3) = ");
    for (i = 1; i <= 3; ++i) {
      printf("%f, ",dbxem[i - 1]);
    }
    printf("   \n");
  }
  gMC->Gsvolu("BXHA", "BOX ", idtmed[1501], dbxha, 3);
  gMC->Gsatt("BXHA", "SEEN", 0);
  // **> Divide along Z to obtain one layer 
  gMC->Gsdvn("RWHA", "BXHA", 4, 3);
  gMC->Gsatt("RWHA", "SEEN", 0);
  // **> Divide along X' to accomodate the maximum number of individual 
  //**> fibres packed along X', make the divisions invisible on the drawings.
  nfx = Int_t(xxuthi / .045);
  if (debugFlag > 0) {
    printf(" NfxHad = %d\n",nfx);
  }
  gMC->Gsdvn("FXHA", "RWHA", nfx, 1);
  gMC->Gsatt("FXHA", "SEEN", 0);
  // **> Construct one fiber cladding 
  dclha[0] = 0.;
  dclha[1] = 0.5*kDiamCladding;
  dclha[2] = dbxha[1];
  if (debugFlag > 0) {
    printf(" DclHa(1-3) = ");
    for (i = 1; i <= 3; ++i) {
      printf("%f, ",dclha[i - 1]);
    }
    printf("   \n");
  }
  gMC->Gsvolu("CLHA", "TUBE", idtmed[fOdCladding - 1], dclha,3);
  gMC->Gsatt("CLHA", "SEEN", 0);
  //**> Construct the cylindrical volume for a fibre core in the Had section.
  //**> Fill with selected fibre material, make it invisible on the drawings.
  dcoha[0] = 0.;
  dcoha[1] = 0.5*kDiamCore;
  dcoha[2] = dbxha[1];
  if (debugFlag > 0) {
    printf(" DcoHa(1-3) = ");
    for (i = 1; i <= 3; ++i) {
      printf("%f, ",dcoha[i - 1]);
    }
    printf("   \n");
  }
  gMC->Gsvolu("COHA", "TUBE", idtmed[fOdFiber - 1], dcoha,3);
  gMC->Gsatt("COHA", "SEEN", 0);
  // **> Position the volumes 
  // **> Put the air section inside one sampling module 
  // **> Use MANY to obtain clipping of protruding edges. 
  xp = 0.;
  zp = dlayha / 2. - 0.5*kFibersHad*kDiamCladding;
  yp = zp * TMath::Tan(thecen);
  gMC->Gspos("BXHA", 1, "SLHA", xp, yp, zp, 0, "MANY");
  // **> Place the core fibre in the clad 
  xp = 0.;
  yp = 0.;
  zp = 0.;
  gMC->Gspos("COHA", 1, "CLHA", xp, yp, zp, 0, "MANY");
  // **> Place the fibre in its air box 
  gMC->Gspos("CLHA", 1, "FXHA", xp, yp, zp, idrotm[0], "MANY");
  // **> Rotation matrices for consecutive calorimeter octants 
  // **> filling the imaginary box. 
  AliMatrix(idrotm[1], 90., -90., 45., 0., 45., 180.);
  // **> Place the EM and Hadronic sections inside the Octant. 
  rzlow = (doct[5] + doct[6]) * .5;
  rzhig = (doct[8] + doct[9]) * .5;
  zp = doct[7] - (faceut * TMath::Cos(beta) + doctha * fact3) * .5;
  yp = 0.;
  xp = rzlow + (rzhig - rzlow) * .5 * (zp - doct[4]) / doct[7];
  gMC->Gspos("HAD ", 1, "OCT ", xp, yp, zp, idrotm[1], "ONLY");
  yp = 0.;
  zp = doct[7] - faceut * TMath::Cos(beta) * .5 - doctha * fact3 - doctem * fact3 * .5;
  xp = rzlow + (rzhig - rzlow) * .5 * (zp - doct[4]) / doct[7];
  gMC->Gspos("EM  ", 1, "OCT ", xp, yp, zp, idrotm[1], "ONLY");
  // **> An imaginary box to hold the complete calorimeter. 
  dcal[0] = (faceut + zendha * fact2) * TMath::Sin(beta);
  dcal[1] = dcal[0];
  dcal[2] = (zendha - kZbegem + faceut * TMath::Cos(beta)) / 2.;
  if (debugFlag > 0) {
    printf(" Dcal(1-3) = ");
    for (i = 1; i <= 3; ++i) {
      printf("%f, ",dcal[i - 1]);
    }
    printf("   \n");
  }
  gMC->Gsvolu("CAL ", "BOX ", idtmed[1501], dcal, 3);
  // Fill with air 
  rinbeg = TMath::Tan(alfa2) * kZbegem;
  rutbeg = TMath::Tan(alfa1) * kZbegem;
  dztotl = dcal[2] * 2.;
  rinend = (dztotl + kZbegem) * TMath::Tan(alfa2);
  rutend = (dztotl + kZbegem) * TMath::Tan(alfa1);
  if (debugFlag > 0) {
    printf(" RinBeg=%f, RoutBeg=%f\n",rinbeg,rutbeg);
    printf(" RinEnd=%f, RoutEnd=%f\n",rinend,rutend);
    printf(" DeltaZtotal = %f\n",dztotl);
  }
  // **> Build the calorimeter inside the imaginary box. 
  rxyin = (fact2 + fact1 / 2.) * kZbegem; // Radius to centre of octant in X'Y' 
  // plane at calorimeter entrance. 
  rxyut = zendha * (fact2 + fact1 / 2.);  // Radius to centre of octant in X'Y'
  // plane at calorimeter exit. 
  rxy   = (rxyin + rxyut) / 2.;           // Radius to geometrical centre of octant in 
  rxy  *= TMath::Sin(beta);               // projected to the XY plane. 
  if (debugFlag > 0) {
    printf(" \n");
  }
  gMC->Gspos("OCTA", 1, "CAL ", 0., 0., 0., 0, "ONLY");
  //**> Construct the narrow stainless steel conical beam tube traversing the
  // **> calorimeter and its vacuum filling:  WallThickness = 0.1 cm, 
  // **> Router = touching the inner side of the calorimeter, 
  // **> DeltaZ = all through the calorimeter box. 
  dcalt[0] = dcal[2];
  dcalt[2] = TMath::Tan(alfa2) * kZbegem;
  dcalt[1] = dcalt[2] - .1 / TMath::Cos(alfa2);
  dcalt[4] = (dcalt[0] * 2. + kZbegem) * TMath::Tan(alfa2);
  dcalt[3] = dcalt[4] - .1 / TMath::Cos(alfa2);
  dcalv[0] = dcalt[0];
  dcalv[2] = dcalt[1];
  dcalv[1] = 0.;
  dcalv[4] = dcalt[3];
  dcalv[3] = 0.;
  gMC->Gsvolu("CALT", "CONE", idtmed[1506], dcalt, 5);
  // Fe (steel a 
  gMC->Gsvolu("CALV", "CONE", idtmed[1500], dcalv, 5);
  // Vacuum. 
  gMC->Gsatt("CALV", "SEEN", 0);
  // **> Position at centre of calorimeter box. 
  zp = 0.;
  gMC->Gspos("CALT", 1, "CAL ", 0., 0., zp, 0, "ONLY");
  gMC->Gspos("CALV", 1, "CAL ", 0., 0., zp, 0, "ONLY");
  if (debugFlag > 0) {
    printf(" Dcalt,Zp,-/+ = ");
    for (i = 1; i <= 5; ++i) {
      printf("%f, ",dcalt[i - 1]);
    }
    printf("%f, %f, %f\n",zp, zp - dcalt[0], zp + dcalt[0]);
    printf(" Dcalt,Zp,-/+ = ");
    for (i = 1; i <= 5; ++i) {
      printf("%f, ",dcalt[i - 1]);
    }
    printf("%f, %f, %f\n",zp, zp - dcalt[0], zp + dcalt[0]);
  }
  // **> Rotate the imaginary box carrying the calorimeter and place it 
  // **> in the ALICE volume on the -Z side. 
  xp = 0.;
  yp = 0.;
  zp = dcal[2] + kZbegem;
  AliMatrix(idrotm[2], 90., 180., 90., 90., 180., 0.);
  // -X theta and phi w.r.t. to box XYZ. 
  //  Y theta and phi w.r.t. to box XYZ. 
  // -Z theta and phi w.r.t. to box XYZ. 
  gMC->Gspos("CAL ", 1, "ALIC", xp, yp, -zp, idrotm[2], "ONLY");
  if (debugFlag > 0) {
    printf(" Dcal,Zp,-/+ = ");
    for (i = 1; i <= 3; ++i) {
      printf("%f, ",dcal[i - 1]);
    }
    printf("%f, %f, %f\n",zp, zp - dcal[2], zp + dcal[2]);
  }
}

//_____________________________________________________________________________
void AliCASTORv1::DrawModule()
{
  //
  // Draw a shaded view of CASTOR version 1
  //

  
  gMC->Gsatt("*", "seen", -1);
  gMC->Gsatt("alic", "seen", 0);
  //
  // Set visibility of elements
  gMC->Gsatt("OCTA","seen",0);
  gMC->Gsatt("EM  ","seen",0);
  gMC->Gsatt("HAD ","seen",0);
  gMC->Gsatt("CAL ","seen",0);
  gMC->Gsatt("CALT","seen",1);
  gMC->Gsatt("OCT ","seen",0);
  gMC->Gsatt("SLEM","seen",1);
  gMC->Gsatt("SLHA","seen",1);
  gMC->Gsatt("SAHA","seen",1);
  //
  gMC->Gdopt("hide", "on");
  gMC->Gdopt("shad", "on");
  gMC->Gsatt("*", "fill", 7);
  gMC->SetClipBox(".");
  gMC->SetClipBox("*", 0, 20, -20, 20, -1900, -1700);
  gMC->DefaultRange();
  gMC->Gdraw("alic", 40, 30, 0, -191.5, -78, .19, .19);
  gMC->Gdhead(1111, "CASTOR Version 1");
  gMC->Gdman(15,-2, "MAN");
  gMC->Gdopt("hide", "off");
}

//_____________________________________________________________________________
void AliCASTORv1::CreateMaterials()
{
  //
  // Create materials for CASTOR version 1
  //
  //   30 March 1997   27 November 1997              Aris L. S. Angelis   * 
  // >--------------------------------------------------------------------<* 
  Int_t   ISXFLD = gAlice->Field()->Integ();
  Float_t SXMGMX = gAlice->Field()->Max();
  
  Int_t *idtmed = fIdtmed->GetArray()-1499;
  
  Float_t cute, ubuf[1], cutg, epsil, awmix[3], dwmix, stmin;
  Int_t isvol;
  Float_t wwmix[3], zwmix[3], aq[2], dq, zq[2], wq[2];
  Float_t tmaxfd, stemax, deemax;
  Int_t kod;
  
  
  // **> Quartz and Wmixture. 
  // **> UBUF is the value of r0, used for calculation of the radii of 
  // **> the nuclei and the Woods-Saxon potential. 
  ubuf[0] = .68;
  AliMaterial(1, "Vacuum$", 1e-16, 1e-16, 1e-16, 1e16, 1e16, ubuf, 1);
  ubuf[0] = .68;
  AliMaterial(2, "Air   $", 14.61, 7.3, .001205, 30420., 67500., ubuf, 1);
  //**> Quartz (SiO2) and fluorinated (?) quartz for cladding (insensitive).
  dq    = 2.64;
  aq[0] = 28.086;
  aq[1] = 15.9994;
  zq[0] = 14.;
  zq[1] = 8.;
  wq[0] = 1.;
  wq[1] = 2.;
  AliMixture(3, "Quartz$", aq, zq, dq, -2, wq);
  // After a call with ratios by number (negative number of elements), 
  // the ratio array is changed to the ratio by weight, so all successive 
  // calls with the same array must specify the number of elements as 
  // positive 
  AliMixture(4, "FQuartz$", aq, zq, dq, 2, wq);
  // **> W mixture (90% W + 7.5% Ni + 2.5% Cu). 
  awmix[0] = 183.85;
  zwmix[0] = 74.;
  wwmix[0] = .9;
  awmix[1] = 58.69;
  zwmix[1] = 28.;
  wwmix[1] = .075;
  awmix[2] = 63.55;
  zwmix[2] = 29.;
  wwmix[2] = .025;
  dwmix    = 17.2;
  // **> (Pure W and W mixture are given the same material number 
  // **> so that they can be used interchangeably). 
  ubuf[0] = 1.1;
  AliMixture(5, "W Mix $", awmix, zwmix, dwmix, 3, wwmix);
  // **> Lead. 
  ubuf[0] = 1.12;
  AliMaterial(6, "Pb208 $", 207.19, 82., 11.35, .56, 18.5, ubuf, 1);
  // **> Iron. 
  ubuf[0] = .99;
  AliMaterial(7, "Fe56  $", 55.85, 26., 7.87, 1.76, 16.7, ubuf, 1);
  // **> Copper. 
  ubuf[0] = 1.01;
  AliMaterial(8, "Cu63  $", 63.54, 29., 8.96, 1.43, 15., ubuf, 1);
  // **> Debug Printout. 
  //      CALL GPRINT('MATE',0) 
  // **> (Negative values for automatic calculation in case of AUTO=0). 
  isvol  = 0;    // Sensitive volume flag. 
  tmaxfd = .1;   // Max allowed angular deviation in 1 step due to field 
  stemax = -.5;  // Maximum permitted step size (cm). 
  deemax = -.2;  // Maximum permitted fractional energy loss. 
  epsil  = .01;  // Boundary crossing precision (cm). 
  stmin  = -.1;  // Minimum permitted step size inside absorber (cm). 
  AliMedium(1, "Vacuum$", 1, isvol, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(2, "Air   $", 2, isvol, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  
  // **> Options for Cherenkov fibres and cladding. 
  isvol = 1;    // Declare fibre core as sensitive. 
  AliMedium(3, "Quartz$", 3, isvol, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  isvol = 0;    // Declare fibre cladding as not sensitive. 
  AliMedium(4, "FQuartz$", 4, isvol, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  
  // **> Options for absorber material (not sensitive). 
  isvol  = 0;   // Sensitive volume flag. 
  stemax = .5;  // Maximum permitted step size (cm). 
  deemax = .5;  // Maximum permitted fractional energy loss. 
  stmin  = .1;  // Minimum permitted step size inside absorber (cm). 
  AliMedium(5, "W Mix $",  5, isvol, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(6, "Pb208 $",  6, isvol, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(7, "Fe56  $ ", 7, isvol, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  AliMedium(8, "Cu63  $ ", 8, isvol, ISXFLD, SXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
  
  // **> Select material for the Cherenkov fibres. 
  fOdFiber    = 1503;
  //      CALL GPTMED(IDTMED(KODFBR)) 
  // **> Select material for the fibre cladding. 
  // Quartz. 
  fOdCladding = 1504;
  //      CALL GPTMED(IDTMED(KODCLD)) 
  // **> Select absorber material. 
  // FQuartz. 
  fOdAbsorber = 1505;  // W184/Mix 
  //      KODABS=1506   ! Pb208. 
  //      KODABS=1507   ! Fe56. 
  //      KODABS=1508   ! Cu63. 
  //      CALL GPTMED(IDTMED(KODABS)) 
  // **> Set by default all interactions and decays explicitly ON 
  // **> and redefine the kinetic energy cutoffs: 
  //      CUTE=0.0031       ! Allow beta >= 0.99 only. 
  cute = 7e-4;  // Allow beta >= 0.67 only. 
  cutg = cute * 1.33;
  
  // **> Inside the absorber material, 
  for (kod = 1505; kod <= 1508; ++kod) {
    Int_t absorber = idtmed[kod - 1];
    gMC->Gstpar(absorber, "CUTELE", cute);  // Allow beta >= 0.xx 
    gMC->Gstpar(absorber, "CUTGAM", cutg);  // = 1.33 cutele. 
    gMC->Gstpar(absorber, "CUTNEU", .01);   // Default. 
    gMC->Gstpar(absorber, "CUTHAD", .01);   // Default. 
    gMC->Gstpar(absorber, "CUTMUO", .01);   // Default. 
    gMC->Gstpar(absorber, "BCUTE", cutg);   // = cutgam. 
    gMC->Gstpar(absorber, "BCUTM", cutg);   // = cutgam. 
    gMC->Gstpar(absorber, "DCUTE", cute);   // = cutele. 
    gMC->Gstpar(absorber, "DCUTM", cute);   // = cutele. 
    gMC->Gstpar(absorber, "PPCUTM", cutg);  // = 1.33 cutele. 
    gMC->Gstpar(absorber, "DCAY", 1.);
    gMC->Gstpar(absorber, "MULS", 1.);
    gMC->Gstpar(absorber, "PFIS", 1.);
    gMC->Gstpar(absorber, "MUNU", 1.);
    gMC->Gstpar(absorber, "LOSS", 1.);
    gMC->Gstpar(absorber, "PHOT", 1.);
    gMC->Gstpar(absorber, "COMP", 1.);
    gMC->Gstpar(absorber, "PAIR", 1.);
    gMC->Gstpar(absorber, "BREM", 1.);
    gMC->Gstpar(absorber, "RAYL", 1.);
    gMC->Gstpar(absorber, "DRAY", 1.);
    gMC->Gstpar(absorber, "ANNI", 1.);
    gMC->Gstpar(absorber, "HADR", 1.);
    gMC->Gstpar(absorber, "LABS", 1.);
  }
  // **> Inside the cladding, 
  Int_t cladding = idtmed[fOdCladding - 1];
  gMC->Gstpar(cladding, "CUTELE", cute);  // Allow beta >= 0.xx 
  gMC->Gstpar(cladding, "CUTGAM", cutg);  // = 1.33 cutele. 
  gMC->Gstpar(cladding, "CUTNEU", .01);   // Default. 
  gMC->Gstpar(cladding, "CUTHAD", .01);   // Default. 
  gMC->Gstpar(cladding, "CUTMUO", .01);   // Default. 
  gMC->Gstpar(cladding, "BCUTE", cutg);   // = cutgam. 
  gMC->Gstpar(cladding, "BCUTM", cutg);   // = cutgam. 
  gMC->Gstpar(cladding, "DCUTE", cute);   // = cutele. 
  gMC->Gstpar(cladding, "DCUTM", cute);   // = cutele. 
  gMC->Gstpar(cladding, "PPCUTM", cutg);  // = 1.33 cutele. 
  gMC->Gstpar(cladding, "DCAY", 1.);
  gMC->Gstpar(cladding, "MULS", 1.);
  gMC->Gstpar(cladding, "PFIS", 1.);
  gMC->Gstpar(cladding, "MUNU", 1.);
  gMC->Gstpar(cladding, "LOSS", 1.);
  gMC->Gstpar(cladding, "PHOT", 1.);
  gMC->Gstpar(cladding, "COMP", 1.);
  gMC->Gstpar(cladding, "PAIR", 1.);
  gMC->Gstpar(cladding, "BREM", 1.);
  gMC->Gstpar(cladding, "RAYL", 1.);
  gMC->Gstpar(cladding, "DRAY", 1.);
  gMC->Gstpar(cladding, "ANNI", 1.);
  gMC->Gstpar(cladding, "HADR", 1.);
  gMC->Gstpar(cladding, "LABS", 1.);
  
  // **> and Inside the Cherenkov fibres, 
  Int_t fiber = idtmed[fOdFiber - 1];
  gMC->Gstpar(fiber, "CUTELE", cute);  // Allow beta >= 0.xx 
  gMC->Gstpar(fiber, "CUTGAM", cutg);  // = 1.33 cutele. 
  gMC->Gstpar(fiber, "CUTNEU", .01);   // Default. 
  gMC->Gstpar(fiber, "CUTHAD", .01);   // Default. 
  gMC->Gstpar(fiber, "CUTMUO", .01);   // Default. 
  gMC->Gstpar(fiber, "BCUTE", cutg);   // = cutgam. 
  gMC->Gstpar(fiber, "BCUTM", cutg);   // = cutgam. 
  gMC->Gstpar(fiber, "DCUTE", cute);   // = cutele. 
  gMC->Gstpar(fiber, "DCUTM", cute);   // = cutele. 
  gMC->Gstpar(fiber, "PPCUTM", cutg);  // = 1.33 cutele. 
  gMC->Gstpar(fiber, "DCAY", 1.);
  gMC->Gstpar(fiber, "MULS", 1.);
  gMC->Gstpar(fiber, "PFIS", 1.);
  gMC->Gstpar(fiber, "MUNU", 1.);
  gMC->Gstpar(fiber, "LOSS", 1.);
  gMC->Gstpar(fiber, "PHOT", 1.);
  gMC->Gstpar(fiber, "COMP", 1.);
  gMC->Gstpar(fiber, "PAIR", 1.);
  gMC->Gstpar(fiber, "BREM", 1.);
  gMC->Gstpar(fiber, "RAYL", 1.);
  gMC->Gstpar(fiber, "DRAY", 1.);
  gMC->Gstpar(fiber, "ANNI", 1.);
  gMC->Gstpar(fiber, "HADR", 1.);
  gMC->Gstpar(fiber, "LABS", 1.);
}

//_____________________________________________________________________________
void AliCASTORv1::StepManager()
{
  //
  // Called at every step in CASTOR
  //
}

//_____________________________________________________________________________
void AliCASTORv1::Init()
{
  //
  // Initialise CASTOR detector after it has been built
  //
  Int_t i;
  //
  printf("\n");
  for(i=0;i<35;i++) printf("*");
  printf(" CASTOR_INIT ");
  for(i=0;i<35;i++) printf("*");
  printf("\n");
  //
  // Here the ABSO initialisation code (if any!)
  for(i=0;i<80;i++) printf("*");
  printf("\n");
}

ClassImp(AliCASTORhit)

//_____________________________________________________________________________
AliCASTORhit::AliCASTORhit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits):
AliHit(shunt, track)
{
  //
  // Store a CASTOR hit
  //
  fVolume  = vol[0];
  fX=hits[0];
  fY=hits[1];
  fZ=hits[2];
}
 
 
