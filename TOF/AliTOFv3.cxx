///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Time Of Flight                                                           //
//  This class contains the functions for version 3 of the Time Of Flight    //
//  detector.                                                                //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliTOFv3Class.gif">
*/
//End_Html
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTOFv3.h"
#include <TNode.h>
#include <TTUBE.h>
#include "AliRun.h"
#include "AliMC.h"
 
ClassImp(AliTOFv3)
 
//_____________________________________________________________________________
AliTOFv3::AliTOFv3()
{
  //
  // Default constructor
  //
}
 
//_____________________________________________________________________________
AliTOFv3::AliTOFv3(const char *name, const char *title)
       : AliTOF(name,title)
{
  //
  // Standard constructor
  //
}
 
//_____________________________________________________________________________
void AliTOFv3::CreateGeometry()
{
  //
  // Create geometry for Time Of Flight version 2
  //
  //Begin_Html
  /*
    <img src="picts/AliTOFv3.gif">
  */
  //End_Html
  //

  //
  // Create common geometry between version 2 and 3
  //
  AliTOF::CreateGeometry();
}
 
//_____________________________________________________________________________
void AliTOFv3::TOFpc(Float_t xm, Float_t ym, Float_t zm0,
		     Float_t zm1, Float_t zm2)
{
  //
  // Definition of the Time Of Fligh Parallel Plate Chambers
  //

  Int_t inum;
  Float_t xcor, zcor, ytop;
  Int_t inum1;
  Float_t xcor1, xcor2, ycoor;
  Float_t stepx, stepz, dx, dy, dz, xp, yp, zp, shiftx, shiftz, ywidth;
  Float_t shiftx1, shiftx2, xad, zad;
  Int_t ink;
  Float_t par[10];
  Int_t inz;
  Float_t xzd;
  Int_t nxp, npx, npz;
  Float_t xsz, ysz, zsz;
  Int_t nzp0, nzp1, nzp2;
  
  Int_t *idtmed = fIdtmed->GetArray()-499;
  
  // X size of PPC plate 
  xsz = 54.;
  // Y size of PPC plate 
  ysz = .2;
  // Z size of PPC plate 
  zsz = 48.;
  // First return additional shift along X 
  xad = 1.5;
  // Second return additional shift along X 
  xzd = .5;
  // Return additional shift along Z 
  zad = .25;
  // Width of DME box 
  ywidth = 4.;
  // X size of PPC chamber 
  xp = 5.7;
  // Y size of PPC chamber 
  yp = .32;
  // Z size of PPC chamber 
  zp = 5.7;
  // Frame width along X,Y and Z axis of PPC chambers 
  dx = .2;
  dy = .1;
  dz = .2;
  // No sensitive volumes with DME 
  par[0] = xm / 2.;
  par[1] = ywidth / 2.;
  par[2] = zm0 / 2.;
  ycoor = ym / 3. - ywidth / 2.;
  gMC->Gsvolu("FBT1", "BOX ", idtmed[505], par, 3);
  gMC->Gspos("FBT1", 0, "FTO1", 0., 0., 0., 0, "ONLY");
  par[2] = zm1 / 2.;
  gMC->Gsvolu("FBT2", "BOX ", idtmed[505], par, 3);
  gMC->Gspos("FBT2", 1, "FTO2", 0., 0., 0., 0, "ONLY");
  par[2] = zm2 / 2.;
  gMC->Gsvolu("FBT3", "BOX ", idtmed[505], par, 3);
  gMC->Gspos("FBT3", 2, "FTO3", 0., 0., 0., 0, "ONLY");
  // Electronic plate 
  par[1] = ysz / 2.;
  par[2] = zm0 / 2.;
  ycoor = ywidth / 2. - ysz / 2.;
  gMC->Gsvolu("FPE1", "BOX ", idtmed[504], par, 3);
  gMC->Gspos("FPE1", 0, "FBT1", 0., ycoor, 0., 0, "ONLY");
  gMC->Gspos("FPE1", 1, "FBT1", 0., -ycoor, 0., 0, "ONLY");
  par[2] = zm1 / 2.;
  gMC->Gsvolu("FPE2", "BOX ", idtmed[504], par, 3);
  gMC->Gspos("FPE2", 0, "FBT2", 0., ycoor, 0., 0, "ONLY");
  gMC->Gspos("FPE2", 1, "FBT2", 0., -ycoor, 0., 0, "ONLY");
  par[2] = zm2 / 2.;
  gMC->Gsvolu("FPE3", "BOX ", idtmed[504], par, 3);
  gMC->Gspos("FPE3", 0, "FBT3", 0., ycoor, 0., 0, "ONLY");
  gMC->Gspos("FPE3", 1, "FBT3", 0., -ycoor, 0., 0, "ONLY");
  // Electronic insensitive volumes 
  par[1] = yp / 2.;
  par[2] = zm0 / 2.;
  ytop = ywidth / 2. - (ysz * 2 + yp) / 2.;
  gMC->Gsvolu("FST1", "BOX ", idtmed[505], par, 3);
  gMC->Gsvolu("FLT1", "BOX ", idtmed[505], par, 3);
  gMC->Gspos("FST1", 0, "FBT1", 0., ytop, 0., 0, "ONLY");
  gMC->Gspos("FLT1", 0, "FBT1", 0., -ytop, 0., 0, "ONLY");
  par[2] = zm1 / 2.;
  gMC->Gsvolu("FST2", "BOX ", idtmed[505], par, 3);
  gMC->Gsvolu("FLT2", "BOX ", idtmed[505], par, 3);
  gMC->Gspos("FST2", 0, "FBT2", 0., ytop, 0., 0, "ONLY");
  gMC->Gspos("FLT2", 0, "FBT2", 0., -ytop, 0., 0, "ONLY");
  par[2] = zm2 / 2.;
  gMC->Gsvolu("FST3", "BOX ", idtmed[505], par, 3);
  gMC->Gsvolu("FLT3", "BOX ", idtmed[505], par, 3);
  gMC->Gspos("FST3", 0, "FBT3", 0., ytop, 0., 0, "ONLY");
  gMC->Gspos("FLT3", 0, "FBT3", 0., -ytop, 0., 0, "ONLY");
  // PPC-plate number along X axis 
  nxp = Int_t (xm / xsz);
  // PPC-plate number along Z axis 
  nzp0 = Int_t (zm0 / zsz);
  nzp1 = Int_t (zm1 / zsz);
  nzp2 = Int_t (zm2 / zsz);
  // Position of big PPC-plate 
  par[0] = xm * .5 / nxp;
  par[2] = zm0 * .5 / nzp0;
  gMC->Gsvolu("FSK1", "BOX ", idtmed[505], par, 3);
  gMC->Gsvolu("FLK1", "BOX ", idtmed[505], par, 3);
  inum = 0;
  for (ink = 1; ink <= nxp; ++ink) {
    xcor = xm * .5 * ((ink * 2 - 1) / (Float_t) nxp - 1.);
    for (inz = 1; inz <= nzp0; ++inz) {
      zcor = zm0 * .5 * ((inz * 2 - 1) / (Float_t) nzp0 - 1.);
      ++inum;
      gMC->Gspos("FSK1", inum, "FST1", xcor, 0., zcor, 0, "ONLY");
      gMC->Gspos("FLK1", inum, "FLT1", xcor, 0., zcor, 0, "ONLY");
    }
    for (inz = 1; inz <= nzp1; ++inz) {
      zcor = zm1 * .5 * ((inz * 2 - 1) / (Float_t) nzp1 - 1.);
      ++inum;
      gMC->Gspos("FSK1", inum, "FST2", xcor, 0., zcor, 0, "ONLY");
      gMC->Gspos("FLK1", inum, "FLT2", xcor, 0., zcor, 0, "ONLY");
    }
    for (inz = 1; inz <= nzp2; ++inz) {
      zcor = zm2 * .5 * ((inz * 2 - 1) / (Float_t) nzp2 - 1.);
      ++inum;
      gMC->Gspos("FSK1", inum, "FST3", xcor, 0., zcor, 0, "ONLY");
      gMC->Gspos("FLK1", inum, "FLT3", xcor, 0., zcor, 0, "ONLY");
    }
  }
  par[0] = xsz / 2.;
  par[1] = yp / 2.;
  par[2] = zsz / 2.;
  gMC->Gsvolu("FSL1", "BOX ", idtmed[505], par, 3);
  gMC->Gsvolu("FLL1", "BOX ", idtmed[505], par, 3);
  shiftx = (xp / 2. + xad / 2.) / 2.;
  shiftz = (zm0 / nzp0 - zsz) / 2.;
  gMC->Gspos("FSL1", 0, "FSK1", -shiftx, 0., -shiftz, 0, "ONLY");
  gMC->Gspos("FLL1", 0, "FLK1", shiftx, 0., shiftz, 0, "ONLY");
  // PPC position on PPC-plate 
  npx = 4;
  npz = 8;
  par[0] = xp / 2.;
  par[1] = yp / 2.;
  par[2] = zp / 2.;
  stepx = (xad + xzd + xp * 2) / 2.;
  stepz = (zp + zad) / 2.;
  shiftz = npz * (zad + zp) / 2.;
  shiftx = npx * (xp * 2 + xad + xzd) / 2.;
  shiftx1 = (xp * 2 + xzd + xad) / 2. - xp / 2.;
  shiftx2 = (xp * 2 + xzd + xad) / 2. - xp / 2. - xzd;
  gMC->Gsvolu("FPG1", "BOX ", idtmed[507], par, 3);
  for (ink = 1; ink <= npx; ++ink) {
    xcor1 = -shiftx + stepx * (ink * 2 - 1) - shiftx1;
    xcor2 = -shiftx + stepx * (ink * 2 - 1) + shiftx2;
    for (inz = 1; inz <= npz; ++inz) {
      zcor = -shiftz + stepz * (inz * 2 - 1);
      ++inum;
      inum1 = npx * npz + inum;
      gMC->Gspos("FPG1", inum, "FSL1", xcor1, 0., zcor, 0, "ONLY");
      gMC->Gspos("FPG1", inum1, "FSL1", xcor2, 0., zcor, 0, "ONLY");
      gMC->Gspos("FPG1", inum, "FLL1", xcor1, 0., zcor, 0, "ONLY");
      gMC->Gspos("FPG1", inum1, "FLL1", xcor2, 0., zcor, 0, "ONLY");
    }
  }
  par[0] = xp / 2. - dx;
  par[1] = yp / 2. - dy;
  par[2] = zp / 2. - dz;
  gMC->Gsvolu("FPG2", "BOX ", idtmed[509], par, 3);
  gMC->Gspos("FPG2", 0, "FPG1", 0., 0., 0., 0, "ONLY");
}

//_____________________________________________________________________________
void AliTOFv3::DrawModule()
{
  //
  // Draw a shaded view of the Time Of Flight version 3
  //

  // Set everything unseen
  gMC->Gsatt("*", "seen", -1);
  // 
  // Set ALIC mother transparent
  gMC->Gsatt("ALIC","SEEN",0);
  //
  // Set the volumes visible
  gMC->Gsatt("ALIC","SEEN",0);
  gMC->Gsatt("FBAR","SEEN",0);
  gMC->Gsatt("FTO1","SEEN",0);
  gMC->Gsatt("FTO2","SEEN",0);
  gMC->Gsatt("FTO3","SEEN",0);
  gMC->Gsatt("FBT1","SEEN",0);
  gMC->Gsatt("FBT2","SEEN",0);
  gMC->Gsatt("FBT3","SEEN",0);
  gMC->Gsatt("FST1","SEEN",0);
  gMC->Gsatt("FLT1","SEEN",0);
  gMC->Gsatt("FST2","SEEN",0);
  gMC->Gsatt("FLT2","SEEN",0);
  gMC->Gsatt("FST3","SEEN",0);
  gMC->Gsatt("FLT3","SEEN",0);
  gMC->Gsatt("FSK1","SEEN",0);
  gMC->Gsatt("FLK1","SEEN",0);
  gMC->Gsatt("FSL1","SEEN",1);
  gMC->Gsatt("FLL1","SEEN",1);
  //
  gMC->Gdopt("hide", "on");
  gMC->Gdopt("shad", "on");
  gMC->Gsatt("*", "fill", 7);
  gMC->SetClipBox(".");
  gMC->SetClipBox("*", 0, 1000, -1000, 1000, -1000, 1000);
  gMC->DefaultRange();
  gMC->Gdraw("alic", 40, 30, 0, 12, 9.5, .02, .02);
  gMC->Gdhead(1111, "Time Of Flight");
  gMC->Gdman(18, 4, "MAN");
  gMC->Gdopt("hide","off");
}

//_____________________________________________________________________________
void AliTOFv3::CreateMaterials()
{
  //
  // Define materials for the Time Of Flight
  //
   AliTOF::CreateMaterials();
}
 
//_____________________________________________________________________________
void AliTOFv3::Init()
{
  //
  // Initialise the detector after the geometry has been defined
  //

  AliTOF::Init();
  fIdFTO2=gMC->VolId("FTO2");
  fIdFTO3=gMC->VolId("FTO3");
  fIdFLT1=gMC->VolId("FLT1");
  fIdFLT2=gMC->VolId("FLT2");
  fIdFLT3=gMC->VolId("FLT3");
}
 
//_____________________________________________________________________________
void AliTOFv3::StepManager()
{
  //
  // Procedure called at each step in the Time Of Flight
  //
  Float_t hits[8];
  Int_t vol[3];
  Int_t copy, id;
  Int_t *idtmed = fIdtmed->GetArray()-499;
  if(gMC->GetMedium()==idtmed[510-1] && 
     gMC->TrackEntering() && gMC->TrackCharge()
     && gMC->CurrentVol(0,copy)==fIdSens) {
    TClonesArray &lhits = *fHits;
    //
    // Record only charged tracks at entrance
    gMC->CurrentVolOff(1,0,copy);
    vol[2]=copy;
    gMC->CurrentVolOff(3,0,copy);
    vol[1]=copy;
    id=gMC->CurrentVolOff(6,0,copy);
    vol[0]=copy;
    if(id==fIdFTO3) {
      vol[0]+=22;
      id=gMC->CurrentVolOff(4,0,copy);
      if(id==fIdFLT3) vol[1]+=6;
    } else if (id==fIdFTO2) {
      vol[0]+=20;
      id=gMC->CurrentVolOff(4,0,copy);
      if(id==fIdFLT2) vol[1]+=8;
    } else {
      id=gMC->CurrentVolOff(4,0,copy);
      if(id==fIdFLT1) vol[1]+=14;
    }
    gMC->TrackPosition(hits);
    gMC->TrackMomentum(&hits[3]);
    hits[7]=gMC->TrackTime();
    new(lhits[fNhits++]) AliTOFhit(fIshunt,gAlice->CurrentTrack(),vol,hits);
  }
}
