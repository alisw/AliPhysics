///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Inner Traking System version 1                                           //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliITSv1Class.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>The responsible person for this module is
<a href="mailto:roberto.barbera@ct.infn.it">Roberto Barbera</a>.
</font>
<pre>
*/
//End_Html
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
 
#include <TMath.h>
#include <TRandom.h>
#include <TVector.h>
#include "AliITSv1.h"
#include "AliRun.h"

#include "AliMC.h"
#include "AliConst.h"

ClassImp(AliITSv1)
 
//_____________________________________________________________________________
AliITSv1::AliITSv1() : AliITS() 
{
  //
  // Default constructor for the ITS
  //
}
 
//_____________________________________________________________________________
AliITSv1::AliITSv1(const char *name, const char *title)
  : AliITS(name, title)
{ 
  //
  // Standard constructor for the ITS
  //
}
 
//_____________________________________________________________________________
void AliITSv1::CreateGeometry()
{
  //
  // Create geometry for version 1 of the ITS
  //
  //
  // Create Geometry for ITS version 0
  //
  //Begin_Html
  /*
    <img src="picts/AliITSv1Tree.gif">
  */
  //End_Html
  //Begin_Html
  /*
    <img src="picts/AliITSv1.gif">
  */
  //End_Html
  

  AliMC* pMC = AliMC::GetMC();
  
  Float_t drcer[6] = { 0.,0.,.08,.08,0.,0. };           //CERAMICS THICKNESS
  Float_t drepx[6] = { 0.,0.,0.,0.,.5357,.5357 };       //EPOXY THICKNESS
  Float_t drpla[6] = { 0.,0.,0.,0.,.1786,.1786 };       //PLASTIC THICKNESS
  Float_t dzb[6]   = { 0.,0.,15.,15.,4.,4. };           //LENGTH OF BOXES
  Float_t dphi[6]  = { 72.,72.,72.,72.,50.6,45. };      //COVERED PHI-RANGE FOR LAYERS 1-6
  Float_t rl[6]    = { 3.9,7.6,14.,24.,40.,45. };       //SILICON LAYERS INNER RADIUS
  Float_t drl[6]   = { .755,.755,.809,.809,.7,.7 };     //THICKNESS OF LAYERS (in % radiation length)
  Float_t dzl[6]   = { 12.67,16.91,20.85,29.15,45.11,50.975 };//HALF LENGTH OF LAYERS
  Float_t drpcb[6] = { 0.,0.,.06,.06,0.,0. };           //PCB THICKNESS
  Float_t drcu[6]  = { 0.,0.,.0504,.0504,.0357,.0357 }; //COPPER THICKNESS
  Float_t drsi[6]  = { 0.,0.,.006,.006,.3571,.3571 };   //SILICON THICKNESS

  Float_t drca = 0, dzfc;
  Int_t i, nsec;
  Float_t rend, drca_tpc, dzco, zend, dits[3], rlim, drsu, zmax;
  Float_t zpos, dzco1, dzco2;
  Float_t drcac[6], acone, dphii;
  Float_t pcits[15], xltpc;
  Float_t rzcone, rstep, r0, z0, acable, fp, dz, zi, ri;
  Int_t idrotm[399];
  Float_t dgh[15];
  
  Int_t *idtmed = fIdtmed->GetArray()-199;
  
  //     CONVERT INTO CM (RL(SI)=9.36 CM) 
  for (i = 0; i < 6; ++i) {
    drl[i] = drl[i] / 100. * 9.36;
  }
  
  //     SUPPORT ENDPLANE THICKNESS 
  drsu = 2.*0.06+1./20;  // 1./20. is 1 cm of honeycomb (1/20 carbon density);
  
  //     CONE BELOW TPC 
  
  drca_tpc = 1.2/4.;
  
  //     CABLE THICKNESS (CONICAL CABLES CONNECTING THE LAYERS) 

  
  //     ITS CONE ANGLE 
  
  acone  = 45.;
  acone *= kDegrad;
  
  //     CONE RADIUS AT 1ST LAYER 
  
  rzcone = 30.;
  
  //     FIELD CAGE HALF LENGTH 
  
  dzfc  = 64.5;
  rlim  = 48.;
  zmax  = 80.;
  xltpc = 275.;
  
  
  //     PARAMETERS FOR SMALL (1/2) ITS 

  for (i = 0; i < 6; ++i) {
    dzl[i] /= 2.;
    dzb[i] /= 2.;
  }
  drca     /= 2.;
  acone    /= 2.;
  drca_tpc /= 2.;
  rzcone   /= 2.;
  dzfc     /= 2.;
  zmax     /= 2.;
  xltpc    /= 2.;
  acable    = 15.;
  
  
  
  //     EQUAL DISTRIBUTION INTO THE 6 LAYERS 
  rstep = drca_tpc / 6.;
  for (i = 0; i < 6; ++i) {
    drcac[i] = (i+1) * rstep;
  }

  //     NUMBER OF PHI SECTORS 
  
  nsec = 5;
  
  //     PACK IN PHI AS MUCH AS POSSIBLE 
  //     NOW PACK USING THICKNESS 
  
  for (i = 0; i < 6; ++i) {
    
//     PACKING FACTOR 
    fp = rl[5] / rl[i];
    
    //      PHI-PACKING NOT SUFFICIENT ? 
    
    if (dphi[i]/45 < fp) {
      drcac[i] = drcac[i] * fp * 45/dphi[i];
    }
  }
  
  
  // --- Define ghost volume containing the six layers and fill it with air 
  
  dgh[0] = 3.5;
  dgh[1] = 50.;
  dgh[2] = zmax;
  pMC->Gsvolu("ITSV", "TUBE", idtmed[275], dgh, 3);
  
  // --- Place the ghost volume in its mother volume (ALIC) and make it 
  //     invisible 
  
  pMC->Gspos("ITSV", 1, "ALIC", 0., 0., 0., 0, "ONLY");
  pMC->Gsatt("ITSV", "SEEN", 0);
  
  //     ITS LAYERS (SILICON) 
  
  dits[0] = rl[0];
  dits[1] = rl[0] + drl[0];
  dits[2] = dzl[0];
  pMC->Gsvolu("ITS1", "TUBE", idtmed[199], dits, 3);
  pMC->Gspos("ITS1", 1, "ITSV", 0., 0., 0., 0, "ONLY");
  
  dits[0] = rl[1];
  dits[1] = rl[1] + drl[1];
  dits[2] = dzl[1];
  pMC->Gsvolu("ITS2", "TUBE", idtmed[199], dits, 3);
  pMC->Gspos("ITS2", 1, "ITSV", 0., 0., 0., 0, "ONLY");
  
  dits[0] = rl[2];
  dits[1] = rl[2] + drl[2];
  dits[2] = dzl[2];
  pMC->Gsvolu("ITS3", "TUBE", idtmed[224], dits, 3);
  pMC->Gspos("ITS3", 1, "ITSV", 0., 0., 0., 0, "ONLY");
  
  dits[0] = rl[3];
  dits[1] = rl[3] + drl[3];
  dits[2] = dzl[3];
  pMC->Gsvolu("ITS4", "TUBE", idtmed[224], dits, 3);
  pMC->Gspos("ITS4", 1, "ITSV", 0., 0., 0., 0, "ONLY");
  
  dits[0] = rl[4];
  dits[1] = rl[4] + drl[4];
  dits[2] = dzl[4];
  pMC->Gsvolu("ITS5", "TUBE", idtmed[249], dits, 3);
  pMC->Gspos("ITS5", 1, "ITSV", 0., 0., 0., 0, "ONLY");
  
  dits[0] = rl[5];
  dits[1] = rl[5] + drl[5];
  dits[2] = dzl[5];
  pMC->Gsvolu("ITS6", "TUBE", idtmed[249], dits, 3);
  pMC->Gspos("ITS6", 1, "ITSV", 0., 0., 0., 0, "ONLY");
  
  //    ELECTRONICS BOXES 
  
  //     PCB (layer #3 and #4) 
  
  pMC->Gsvolu("IPCB", "TUBE", idtmed[233], dits, 0);
  for (i = 2; i < 4; ++i) {
    dits[0] = rl[i];
    dits[1] = dits[0] + drpcb[i];
    dits[2] = dzb[i] / 2.;
    zpos = dzl[i] + dits[2];
    pMC->Gsposp("IPCB", i-1, "ITSV", 0., 0., zpos, 0, "ONLY", dits, 3);
    pMC->Gsposp("IPCB", i+1, "ITSV", 0., 0.,-zpos, 0, "ONLY", dits, 3);
  }
  
  //     COPPER (layer #3 and #4) 
  
  pMC->Gsvolu("ICO2", "TUBE", idtmed[234], dits, 0);
  for (i = 2; i < 4; ++i) {
    dits[0] = rl[i] + drpcb[i];
    dits[1] = dits[0] + drcu[i];
    dits[2] = dzb[i] / 2.;
    zpos = dzl[i] + dits[2];
    pMC->Gsposp("ICO2", i-1, "ITSV", 0., 0., zpos, 0, "ONLY", dits, 3);
    pMC->Gsposp("ICO2", i+1, "ITSV", 0., 0.,-zpos, 0, "ONLY", dits, 3);
  }
  
  //     CERAMICS (layer #3 and #4) 
  
  pMC->Gsvolu("ICER", "TUBE", idtmed[235], dits, 0);
  for (i = 2; i < 4; ++i) {
    dits[0] = rl[i] + drpcb[i] + drcu[i];
    dits[1] = dits[0] + drcer[i];
    dits[2] = dzb[i] / 2.;
    zpos = dzl[i] + dits[2];
    pMC->Gsposp("ICER", i-1, "ITSV", 0., 0., zpos, 0, "ONLY", dits, 3);
    pMC->Gsposp("ICER", i+1, "ITSV", 0., 0.,-zpos, 0, "ONLY", dits, 3);
  }
  
  //     SILICON (layer #3 and #4) 
  
  pMC->Gsvolu("ISI2", "TUBE", idtmed[226], dits, 0);
  for (i = 2; i < 4; ++i) {
    dits[0] = rl[i] + drpcb[i] + drcu[i] + drcer[i];
    dits[1] = dits[0] + drsi[i];
    dits[2] = dzb[i] / 2.;
    zpos = dzl[i] + dits[2];
    pMC->Gsposp("ISI2", i-1, "ITSV", 0., 0., zpos, 0, "ONLY", dits, 3);
    pMC->Gsposp("ISI2", i+1, "ITSV", 0., 0.,-zpos, 0, "ONLY", dits, 3);
  }
  
  //     PLASTIC (G10FR4) (layer #5 and #6) 
  
  pMC->Gsvolu("IPLA", "TUBE", idtmed[262], dits, 0);
  for (i = 4; i < 6; ++i) {
    dits[0] = rl[i];
    dits[1] = dits[0] + drpla[i];
    dits[2] = dzb[i] / 2.;
    zpos = dzl[i] + dits[2];
    pMC->Gsposp("IPLA", i-1, "ITSV", 0., 0., zpos, 0, "ONLY", dits, 3);
    pMC->Gsposp("IPLA", i+1, "ITSV", 0., 0.,-zpos, 0, "ONLY", dits, 3);
  }
  
  //     COPPER (layer #5 and #6) 
  
  pMC->Gsvolu("ICO3", "TUBE", idtmed[259], dits, 0);
  for (i = 4; i < 6; ++i) {
    dits[0] = rl[i] + drpla[i];
    dits[1] = dits[0] + drcu[i];
    dits[2] = dzb[i] / 2.;
    zpos = dzl[i] + dits[2];
    pMC->Gsposp("ICO3", i-1, "ITSV", 0., 0., zpos, 0, "ONLY", dits, 3);
    pMC->Gsposp("ICO3", i+1, "ITSV", 0., 0.,-zpos, 0, "ONLY", dits, 3);
  }
  
  //     EPOXY (layer #5 and #6) 
  
  pMC->Gsvolu("IEPX", "TUBE", idtmed[262], dits, 0);
  for (i = 4; i < 6; ++i) {
    dits[0] = rl[i] + drpla[i] + drcu[i];
    dits[1] = dits[0] + drepx[i];
    dits[2] = dzb[i] / 2.;
    zpos = dzl[i] + dits[2];
    pMC->Gsposp("IEPX", i-1, "ITSV", 0., 0., zpos, 0, "ONLY", dits, 3);
    pMC->Gsposp("IEPX", i+1, "ITSV", 0., 0.,-zpos, 0, "ONLY", dits, 3);
  }
  
  //     SILICON (layer #5 and #6) 
  
  pMC->Gsvolu("ISI3", "TUBE", idtmed[251], dits, 0);
  for (i = 4; i < 6; ++i) {
    dits[0] = rl[i] + drpla[i] + drcu[i] + drepx[i];
    dits[1] = dits[0] + drsi[i];
    dits[2] = dzb[i] / 2.;
    zpos = dzl[i] + dits[2];
    pMC->Gsposp("ISI3", i-1, "ITSV", 0., 0., zpos, 0, "ONLY", dits, 3);
    pMC->Gsposp("ISI3", i+1, "ITSV", 0., 0.,-zpos, 0, "ONLY", dits, 3);
  }
  
  //    SUPPORT 
  
  pMC->Gsvolu("ISUP", "TUBE", idtmed[274], dits, 0);
  for (i = 0; i < 6; ++i) {
    dits[0] = rl[i];
    if (i < 5) dits[1] = rl[i];
    else       dits[1] = rlim;
    dits[2] = drsu / 2.;
    zpos = dzl[i] + dzb[i] + dits[2];
    pMC->Gsposp("ISUP", i+1, "ITSV", 0., 0., zpos, 0, "ONLY", dits, 3);
    pMC->Gsposp("ISUP", i+7, "ITSV", 0., 0.,-zpos, 0, "ONLY", dits, 3);
  }
  
  // CABLES (HORIZONTAL) 
  
  pMC->Gsvolu("ICHO", "TUBE", idtmed[278], dits, 0);
  for (i = 0; i < 6; ++i) {
    dits[0] = rl[i];
    dits[1] = dits[0] + drca;
    dits[2] = (rzcone + TMath::Tan(acone) * (rl[i] - rl[0]) - (dzl[i]+ dzb[i] + drsu)) / 2.;
    zpos = dzl[i - 1] + dzb[i] + drsu + dits[2];
    pMC->Gsposp("ICHO", i+1, "ITSV", 0., 0., zpos, 0, "ONLY", dits, 3);
    pMC->Gsposp("ICHO", i+7, "ITSV", 0., 0.,-zpos, 0, "ONLY", dits, 3);
  }
  //    DEFINE A CONICAL GHOST VOLUME FOR THE PHI SEGMENTATION 
  pcits[0] = 0.;
  pcits[1] = 360.;
  pcits[2] = 2.;
  pcits[3] = rzcone;
  pcits[4] = 3.5;
  pcits[5] = rl[0];
  pcits[6] = pcits[3] + TMath::Tan(acone) * (rlim - rl[0]);
  pcits[7] = rlim - rl[0] + 3.5;
  pcits[8] = rlim;
  pMC->Gsvolu("ICMO", "PCON", idtmed[275], pcits, 9);
  AliMatrix(idrotm[200], 90., 0., 90., 90., 180., 0.);
  pMC->Gspos("ICMO", 1, "ITSV", 0., 0., 0., 0, "ONLY");
  pMC->Gspos("ICMO", 2, "ITSV", 0., 0., 0., idrotm[200], "ONLY");
  
  //     DIVIDE INTO NSEC PHI-SECTIONS 
  
  pMC->Gsdvn("ICMD", "ICMO", nsec, 2);
  pMC->Gsatt("ICMO", "SEEN", 0);
  pMC->Gsatt("ICMD", "SEEN", 0);
  
  //     CONICAL CABLES 
  
  pcits[2] = 2.;
  pMC->Gsvolu("ICCO", "PCON", idtmed[278], pcits, 0);
  for (i = 1; i < 6; ++i) {
    pcits[0] = -dphi[i] / 2.;
    pcits[1] = dphi[i];
    if (i < 5) {
      dzco = TMath::Tan(acone) * (rl[i+1] - rl[i]);
    } else {
      dzco1 = zmax - (rzcone + TMath::Tan(acone) * (rl[5] - rl[0])) -2.;
      dzco2 = (rlim - rl[5]) * TMath::Tan(acone);
      if (rl[5] + dzco1 / TMath::Tan(acone) < rlim) {
	dzco = dzco1;
      } else {
	dzco = dzco2;
      }
    }
    pcits[3] = rzcone + TMath::Tan(acone) * (rl[i] - rl[0]);
    pcits[4] = rl[i] - drcac[i] / TMath::Sin(acone);
    pcits[5] = rl[i];
    pcits[6] = pcits[3] + dzco;
    pcits[7] = rl[i] + dzco / TMath::Tan(acone) - drcac[i] / TMath::Sin(acone);
    pcits[8] = rl[i] + dzco / TMath::Tan(acone);
    
    pMC->Gsposp("ICCO", i, "ICMD", 0., 0., 0., 0, "ONLY", pcits, 9);
    
  }
  zend = pcits[6];
  rend = pcits[8];
  
  //  CONICAL CABLES BELOW TPC 
  
  //    DEFINE A CONICAL GHOST VOLUME FOR THE PHI SEGMENTATION 
  pcits[0] = 0.;
  pcits[1] = 360.;
  pcits[2] = 2.;
  pcits[3] = zend;
  pcits[5] = rend;
  pcits[4] = pcits[5] - drca_tpc;
  pcits[6] = xltpc;
  pcits[8] = pcits[4] + (pcits[6] - pcits[3]) * TMath::Tan(acable * kDegrad);
  pcits[7] = pcits[8] - drca_tpc;
  AliMatrix(idrotm[200], 90., 0., 90., 90., 180., 0.);
  pMC->Gsvolu("ICCM", "PCON", idtmed[275], pcits, 9);
  pMC->Gspos("ICCM", 1, "ALIC", 0., 0., 0., 0, "ONLY");
  pMC->Gspos("ICCM", 2, "ALIC", 0., 0., 0., idrotm[200], "ONLY");
  pMC->Gsdvn("ITMD", "ICCM", nsec, 2);
  pMC->Gsatt("ITMD", "SEEN", 0);
  pMC->Gsatt("ICCM", "SEEN", 0);
  
  //     NOW PLACE SEGMENTS WITH DECREASING PHI SEGMENTS INTO THE 
  //     GHOST-VOLUME 
  
  pcits[2] = 2.;
  pMC->Gsvolu("ITTT", "PCON", idtmed[278], pcits, 0);
  r0 = rend;
  z0 = zend;
  dz = (xltpc - zend) / 9.;
  for (i = 0; i < 9; ++i) {
    zi = z0 + i*dz + dz / 2.;
    ri = r0 + (zi - z0) * TMath::Tan(acable * kDegrad);
    dphii = dphi[5] * r0 / ri;
    pcits[0] = -dphii / 2.;
    pcits[1] = dphii;
    pcits[3] = zi - dz / 2.;
    pcits[5] = r0 + (pcits[3] - z0) * TMath::Tan(acable * kDegrad);
    pcits[4] = pcits[5] - drca_tpc;
    pcits[6] = zi + dz / 2.;
    pcits[8] = r0 + (pcits[6] - z0) * TMath::Tan(acable * kDegrad);
    pcits[7] = pcits[8] - drca_tpc;
    
    pMC->Gsposp("ITTT", i+1, "ITMD", 0., 0., 0., 0, "ONLY", pcits, 9);
  }
  
  // --- Outputs the geometry tree in the EUCLID/CAD format 
  
  if (fEuclidOut) {
    pMC->WriteEuclid("ITSgeometry", "ITSV", 1, 5);
  }
}

//_____________________________________________________________________________
void AliITSv1::CreateMaterials()
{
  //
  // Create the materials for ITS
  //
  AliITS::CreateMaterials();
}

//_____________________________________________________________________________
void AliITSv1::Init()
{
  //
  // Initialise the ITS after it has been built
  //
  AliITS::Init();
}  
 
//_____________________________________________________________________________
void AliITSv1::DrawModule()
{ 
  //
  // Draw a shaded view of the FMD version 1
  //

  AliMC* pMC = AliMC::GetMC();
  
  // Set everything unseen
  pMC->Gsatt("*", "seen", -1);
  // 
  // Set ALIC mother visible
  pMC->Gsatt("ALIC","SEEN",0);
  //
  // Set the volumes visible
  pMC->Gsatt("ITSV","SEEN",0);
  pMC->Gsatt("ITS1","SEEN",1);
  pMC->Gsatt("ITS2","SEEN",1);
  pMC->Gsatt("ITS3","SEEN",1);
  pMC->Gsatt("ITS4","SEEN",1);
  pMC->Gsatt("ITS5","SEEN",1);
  pMC->Gsatt("ITS6","SEEN",1);

  pMC->Gsatt("IPCB","SEEN",1);
  pMC->Gsatt("ICO2","SEEN",1);
  pMC->Gsatt("ICER","SEEN",0);
  pMC->Gsatt("ISI2","SEEN",0);
  pMC->Gsatt("IPLA","SEEN",0);
  pMC->Gsatt("ICO3","SEEN",0);
  pMC->Gsatt("IEPX","SEEN",0);
  pMC->Gsatt("ISI3","SEEN",1);
  pMC->Gsatt("ISUP","SEEN",0);
  pMC->Gsatt("ICHO","SEEN",0);
  pMC->Gsatt("ICMO","SEEN",0);
  pMC->Gsatt("ICMD","SEEN",0);
  pMC->Gsatt("ICCO","SEEN",1);
  pMC->Gsatt("ICCM","SEEN",0);
  pMC->Gsatt("ITMD","SEEN",0);
  pMC->Gsatt("ITTT","SEEN",1);

  //
  pMC->Gdopt("hide", "on");
  pMC->Gdopt("shad", "on");
  pMC->Gsatt("*", "fill", 7);
  pMC->SetClipBox(".");
  pMC->SetClipBox("*", 0, 300, -300, 300, -300, 300);
  pMC->DefaultRange();
  pMC->Gdraw("alic", 40, 30, 0, 11, 10, .07, .07);
  pMC->Gdhead(1111, "Inner Tracking System Version 1");
  pMC->Gdman(17, 6, "MAN");
}

//_____________________________________________________________________________
void AliITSv1::StepManager()
{ 
  //
  // Called at every step in the ITS
  //
  Int_t         copy, id;
  Float_t       hits[7];
  Int_t         vol[3];
  Float_t       position[3];
  Float_t       momentum[4];
  TClonesArray &lhits = *fHits;
  AliMC* pMC = AliMC::GetMC();
  //
  if(pMC->TrackCharge() && pMC->Edep()) {
    //
    // Only entering charged tracks
    if((id=pMC->CurrentVol(0,copy))==fIdSens1) {  
      vol[0]=1;
      id=pMC->CurrentVolOff(1,0,copy);      
      vol[1]=copy;
      id=pMC->CurrentVolOff(2,0,copy);
      vol[2]=copy;                       
    } else if(id==fIdSens2) {
      vol[0]=2;
      id=pMC->CurrentVolOff(1,0,copy);       
      vol[1]=copy;
      id=pMC->CurrentVolOff(2,0,copy);
      vol[2]=copy;                    
    } else if(id==fIdSens3) {
      vol[0]=3;
      vol[1]=copy;
      id=pMC->CurrentVolOff(1,0,copy);
      vol[2]=copy;             
    } else if(id==fIdSens4) {
      vol[0]=4;
      vol[1]=copy;
      id=pMC->CurrentVolOff(1,0,copy);
      vol[2]=copy;                  
    } else if(id==fIdSens5) {
      vol[0]=5;
      vol[1]=copy;
      id=pMC->CurrentVolOff(1,0,copy);
      vol[2]=copy;               
    } else if(id==fIdSens6) {
      vol[0]=6;
      vol[1]=copy;
      id=pMC->CurrentVolOff(1,0,copy);
      vol[2]=copy;                      
    } else return;
    pMC->TrackPosition(position);
    pMC->TrackMomentum(momentum);
    hits[0]=position[0];
    hits[1]=position[1];
    hits[2]=position[2];          
    hits[3]=momentum[0]*momentum[3];
    hits[4]=momentum[1]*momentum[3];
    hits[5]=momentum[2]*momentum[3];        
    hits[6]=pMC->Edep();
    new(lhits[fNhits++]) AliITShit(fIshunt,gAlice->CurrentTrack(),vol,hits);
  }      
}
