///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Time Of Flight: design of P.Fonte                                        //
//  This class contains the functions for version 2 of the Time Of Flight    //
//  detector.                                                                //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliTOFv2Class.gif">
*/
//End_Html
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTOFv2.h"
#include "AliRun.h"
#include "AliConst.h"
 
ClassImp(AliTOFv2)
 
//_____________________________________________________________________________
AliTOFv2::AliTOFv2()
{
  //
  // Default constructor
  //
}
 
//_____________________________________________________________________________
AliTOFv2::AliTOFv2(const char *name, const char *title)
       : AliTOF(name,title)
{
  //
  // Standard constructor
  //
}
 
//_____________________________________________________________________________
void AliTOFv2::CreateGeometry()
{
  //
  // Create geometry for Time Of Flight version 0
  //
  //Begin_Html
  /*
    <img src="picts/AliTOFv2.gif">
  */
  //End_Html
  //
  //
  // Create common geometry
  //
  AliTOF::CreateGeometry();
}
 
//_____________________________________________________________________________
void AliTOFv2::TOFpc(Float_t xm, Float_t ym, Float_t zm0,
		     Float_t zm1, Float_t zm2)
{
  //
  // Definition of the Time Of Fligh Resistive Plate Chambers
  // xm, ym, zm - sizes of TOF modules (large)
  
  Int_t inum;
  Float_t xcor, ycor, zcor, ycoor;
  Float_t zazor, dx, dy, dz, xp, yp, zp, ywidth;
  Int_t ink;
  Float_t par[10];
  Int_t inz, nxp, npx, npz;
  Float_t xsz, ysz, zsz;
  Int_t nzp0, nzp1, nzp2;
  
  Int_t *idtmed = fIdtmed->GetArray()-499;
  
  // X size of small RPC plate G10 
  xsz = 60.;
  // Y size (thickness) of large && small RPC plate G10 
  ysz = .26;
  // Z size of small RPC plate G10
  zsz = 50.;
  // Width of CO2 box with RPC
  ywidth = 4.;
  // Frame width along X,Y and Z axis of RPC chambers 
  dx = 0.;
  dy = .3; //this is 1mm(ceramic) + 1mm(Al) + 1mm(polyethelene) 
  dz = 0.;
  // gap in RPC chamber 
  zazor = .03;
  // Sizes of RPC chamber 
  xp = 3.06; //small pixel
//xp = 3.9; //large pixel 
  yp = zazor + dy * 2; //=0.83cm total thickness of RPC
  zp = 3.06; //small pixel
//zp = 4.1; //large pixel
  // Large not sensitive volumes with CO2 
  par[0] = xm / 2.;
  par[1] = ywidth / 2.;
  par[2] = zm0 / 2.;
  gMC->Gsvolu("FBT1", "BOX ", idtmed[506], par, 3); // CO2
  gMC->Gspos("FBT1", 0, "FTO1", 0., 0., 0., 0, "ONLY");
  par[2] = zm1 / 2.;
  gMC->Gsvolu("FBT2", "BOX ", idtmed[506], par, 3); // CO2
  gMC->Gspos("FBT2", 1, "FTO2", 0., 0., 0., 0, "ONLY");
  par[2] = zm2 / 2.;
  gMC->Gsvolu("FBT3", "BOX ", idtmed[506], par, 3); // CO2
  gMC->Gspos("FBT3", 2, "FTO3", 0., 0., 0., 0, "ONLY");
  // Large electronic plate (G10) after
  par[0] = xm / 2.;
  par[1] = ysz / 2.;
  par[2] = zm0 / 2.;
  ycoor = yp + par[1];
  gMC->Gsvolu("FPE1", "BOX ", idtmed[504], par, 3); // G10
  gMC->Gspos("FPE1", 0, "FBT1", 0., ycoor, 0., 0, "ONLY");
  par[2] = zm1 / 2.;
  gMC->Gsvolu("FPE2", "BOX ", idtmed[504], par, 3); // G10
  gMC->Gspos("FPE2", 0, "FBT2", 0., ycoor, 0., 0, "ONLY");
  par[2] = zm2 / 2.;
  gMC->Gsvolu("FPE3", "BOX ", idtmed[504], par, 3); // G10
  gMC->Gspos("FPE3", 0, "FBT3", 0., ycoor, 0., 0, "ONLY");
  // Electronics (5cm thick) after
  //first - Cu (0.02574cm thick - 1.8% X0)
  par[0] = xm / 2.;
  par[1] = 0.02574 / 2.;
  par[2] = zm0 / 2.;
  ycoor = yp + ysz + 5/2 - par[1];
  gMC->Gsvolu("FEC1", "BOX ", idtmed[501], par, 3); // Cu
  gMC->Gspos("FEC1", 0, "FBT1", 0., ycoor, 0., 0, "ONLY");
  par[2] = zm1 / 2.;
  gMC->Gsvolu("FEC2", "BOX ", idtmed[501], par, 3); // Cu
  gMC->Gspos("FEC2", 0, "FBT2", 0., ycoor, 0., 0, "ONLY");
  par[2] = zm2 / 2.;
  gMC->Gsvolu("FEC3", "BOX ", idtmed[501], par, 3); // Cu
  gMC->Gspos("FEC3", 0, "FBT3", 0., ycoor, 0., 0, "ONLY");
  //second - G10 (0.2328cm thick - 1.2% X0)
  par[0] = xm / 2.;
  par[1] = 0.2328 / 2.;
  par[2] = zm0 / 2.;
  ycoor = yp + ysz + 5/2 + par[1];
  gMC->Gsvolu("FEG1", "BOX ", idtmed[504], par, 3); // G10
  gMC->Gspos("FEG1", 0, "FBT1", 0., ycoor, 0., 0, "ONLY");
  par[2] = zm1 / 2.;
  gMC->Gsvolu("FEG2", "BOX ", idtmed[504], par, 3); // G10
  gMC->Gspos("FEG2", 0, "FBT2", 0., ycoor, 0., 0, "ONLY");
  par[2] = zm2 / 2.;
  gMC->Gsvolu("FEG3", "BOX ", idtmed[504], par, 3); // G10
  gMC->Gspos("FEG3", 0, "FBT3", 0., ycoor, 0., 0, "ONLY");
  // Al support (5mm thick) after
  par[0] = xm / 2.;
  par[1] = 0.5 / 2.;
  par[2] = zm0 / 2.;
  ycoor = yp + ysz + par[1];
  gMC->Gsvolu("FSP1", "BOX ", idtmed[508], par, 3); // Al
  gMC->Gspos("FSP1", 0, "FBT1", 0., ycoor, 0., 0, "ONLY");
  par[2] = zm1 / 2.;
  gMC->Gsvolu("FSP2", "BOX ", idtmed[508], par, 3); // Al
  gMC->Gspos("FSP2", 0, "FBT2", 0., ycoor, 0., 0, "ONLY");
  par[2] = zm2 / 2.;
  gMC->Gsvolu("FSP3", "BOX ", idtmed[508], par, 3); // Al
  gMC->Gspos("FSP3", 0, "FBT3", 0., ycoor, 0., 0, "ONLY");
  // Mylar layer in front 0.5mm thick at 5mm from detector
  par[0] = xm / 2.;
  par[1] = 0.05 / 2;
  par[2] = zm0 / 2.;
  ycoor = -yp - 0.5 - par[1];
  gMC->Gsvolu("FMY1", "BOX ", idtmed[511], par, 3); // G10
  gMC->Gspos("FMY1", 0, "FBT1", 0., ycoor, 0., 0, "ONLY");
  par[2] = zm1 / 2.;
  gMC->Gsvolu("FMY2", "BOX ", idtmed[511], par, 3); // G10
  gMC->Gspos("FMY2", 0, "FBT2", 0., ycoor, 0., 0, "ONLY");
  par[2] = zm2 / 2.;
  gMC->Gsvolu("FMY3", "BOX ", idtmed[511], par, 3); // G10
  gMC->Gspos("FMY3", 0, "FBT3", 0., ycoor, 0., 0, "ONLY");
  //  insensitive volumes - large box for RPCs 
  par[1] = yp; // two times thicker than RPC
  par[2] = zm0 / 2.;
  gMC->Gsvolu("FLT1", "BOX ", idtmed[512], par, 3); //Freon not senc.
  gMC->Gspos("FLT1", 0, "FBT1", 0., 0., 0., 0, "ONLY");
  par[2] = zm1 / 2.;
  gMC->Gsvolu("FLT2", "BOX ", idtmed[512], par, 3); //Freon not senc.
  gMC->Gspos("FLT2", 0, "FBT2", 0., 0., 0., 0, "ONLY");
  par[2] = zm2 / 2.;
  gMC->Gsvolu("FLT3", "BOX ", idtmed[512], par, 3); //Freon not senc.
  gMC->Gspos("FLT3", 0, "FBT3", 0., 0., 0., 0, "ONLY");
  // RPC box (small plate) number along X axis 
  nxp = Int_t (xm / xsz);
  // RPC box (small plate) number along Z axis 
  nzp0 = Int_t (zm0 / zsz);
  nzp1 = Int_t (zm1 / zsz);
  nzp2 = Int_t (zm2 / zsz);
  // (small) box (plate) for RPC size with insencitive Freon
  par[0] = xm * .5 / nxp;
  par[1] = yp; // two times thicker than RPC
  par[2] = zm0 * .5 / nzp0;
  gMC->Gsvolu("FLK0", "BOX ", idtmed[512], par, 3); //Freon not sencitive
  // Position of (small) RPC boxes 
  inum = 0;
  for (ink = 1; ink <= nxp; ++ink) {
    xcor = xm * .5 * ((ink * 2 - 1) / (Float_t) nxp -  1.);
    for (inz = 1; inz <= nzp0; ++inz) {
      zcor = zm0 * .5 * ((inz * 2 - 1) / (Float_t) nzp0 - 1.);
      ++inum;
      gMC->Gspos("FLK0", inum, "FLT1", xcor, 0., zcor, 0, "ONLY");
    }
    for (inz = 1; inz <= nzp1; ++inz) {
      zcor = zm1 * .5 * ((inz * 2 - 1) / (Float_t) nzp1 - 1.);
      ++inum;
      gMC->Gspos("FLK0", inum, "FLT2", xcor, 0., zcor, 0, "ONLY");
    }
    for (inz = 1; inz <= nzp2; ++inz) {
      zcor = zm2 * .5 * ((inz * 2 - 1) / (Float_t) nzp2 - 1.);
      ++inum;
      gMC->Gspos("FLK0", inum, "FLT3", xcor, 0., zcor, 0, "ONLY");
    }
  }
  // Polyethilene boxes for RPC cell
  npx = 19; //number of small pixels along X
  npz = 16; //number of small pixels along Z
  //  npx = 15; //large pixel
  //  npz = 12; //large pixel
  par[0] = xsz * .5 / npx;
  par[1] = yp/2; 
  par[2] = zsz * .5 / npz;
  gMC->Gsvolu("FPP0", "BOX ", idtmed[503], par, 3); // Polyethilene
  inum = 0;
  for (ink = 1; ink <= npx; ++ink) {
    xcor = xsz * .5 * ((ink * 2 - 1) / (Float_t) npx - 1.);
    if (ink%2 != 0) ycor=yp/2; else ycor=-yp/2;
    for (inz = 1; inz <= npz; ++inz) {
      zcor = zsz * .5 * ((inz * 2 - 1) / (Float_t) npz - 1.);
      ++inum;
      gMC->Gspos("FPP0", inum, "FLK0", xcor, ycor, zcor, 0, "ONLY");
      ycor=-ycor;
    }
  }
  //Al RPC geometry 
  par[0] = xp / 2.;
  par[1] = yp / 2. - 0.1; //minus 1mm of poliethelene
  par[2] = zp / 2.;
  gMC->Gsvolu("FPA0", "BOX ", idtmed[508], par, 3);// Al
  gMC->Gspos("FPA0", inum, "FPP0", 0., 0., 0., 0, "ONLY");
  //Ceramic RPC geometry 
  par[0] = xp / 2.;
  par[1] = par[1] - 0.1; //minus 1mm of Al
  par[2] = zp / 2.;
  gMC->Gsvolu("FPC0", "BOX ", idtmed[507], par, 3);// Ceramic
  gMC->Gspos("FPC0", inum, "FPA0", 0., 0., 0., 0, "ONLY");
  // Freon gas sencitive volume
  par[0] = xp / 2. - dx;
  par[1] = yp / 2. - dy;
  par[2] = zp / 2. - dz;
  gMC->Gsvolu("FPG0", "BOX ", idtmed[513], par, 3);// Freon 
  gMC->Gspos("FPG0", 0, "FPC0", 0., 0., 0., 0, "ONLY");
}

//_____________________________________________________________________________
void AliTOFv2::DrawModule()
{
  //
  // Draw a shaded view of the Time Of Flight version 2
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
  gMC->Gsatt("FLT1","SEEN",0);
  gMC->Gsatt("FLT2","SEEN",0);
  gMC->Gsatt("FLT3","SEEN",0);
  gMC->Gsatt("FLK0","SEEN",1);
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
void AliTOFv2::CreateMaterials()
{
  //
  // Define materials for the Time Of Flight
  //
  AliTOF::CreateMaterials();
}
 
//_____________________________________________________________________________
void AliTOFv2::Init()
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
void AliTOFv2::StepManager()
{
  //
  // Procedure called at each step in the Time Of Flight
  //
  TLorentzVector mom, pos;
  Float_t hits[8];
  Int_t vol[3];
  Int_t copy, id, i;
  Int_t *idtmed = fIdtmed->GetArray()-499;
  if(gMC->GetMedium()==idtmed[514-1] && 
     gMC->IsTrackEntering() && gMC->TrackCharge()
     && gMC->CurrentVolID(copy)==fIdSens) {
    TClonesArray &lhits = *fHits;
    //
    // Record only charged tracks at entrance
    gMC->CurrentVolOffID(1,copy);
    vol[2]=copy;
    gMC->CurrentVolOffID(3,copy);
    vol[1]=copy;
    id=gMC->CurrentVolOffID(6,copy);
    vol[0]=copy;
    if(id==fIdFTO3) {
      vol[0]+=22;
      id=gMC->CurrentVolOffID(4,copy);
      if(id==fIdFLT3) vol[1]+=4;
    } else if (id==fIdFTO2) {
      vol[0]+=20;
      id=gMC->CurrentVolOffID(4,copy);
      if(id==fIdFLT2) vol[1]+=8;
    } else {
      id=gMC->CurrentVolOffID(4,copy);
      if(id==fIdFLT1) vol[1]+=14;
    }
    gMC->TrackPosition(pos);
    gMC->TrackMomentum(mom);
    //
    Double_t ptot=mom.Rho();
    Double_t norm=1/ptot;
    for(i=0;i<3;++i) {
      hits[i]=pos[i];
      hits[i+3]=mom[i]*norm;
    }
    hits[6]=ptot;
    hits[7]=pos[3];
    new(lhits[fNhits++]) AliTOFhit(fIshunt,gAlice->CurrentTrack(),vol,hits);
  }
}
