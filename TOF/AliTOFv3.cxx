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
*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Time Of Flight: as for version 1 but full coverage                       //
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
#include "AliRun.h"
#include "AliConst.h"
 
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
  // Create geometry for Time Of Flight version 0
  //
  //Begin_Html
  /*
    <img src="picts/AliTOFv3.gif">
  */
  //End_Html
  //
  //
  // Create common geometry
  //
  AliTOF::CreateGeometry();
}
 
//_____________________________________________________________________________
void AliTOFv3::TOFpc(Float_t xm, Float_t ym, Float_t zm0,
		     Float_t zm1, Float_t zm2)
{
  //
  // Definition of the Time Of Fligh Resistive Plate Chambers
  // xm, ym, zm - sizes of TOF modules (large)
  
  Float_t  ycoor;
  Float_t zazor, xp, yp, zp;
  Float_t par[10];
  
  Int_t *idtmed = fIdtmed->GetArray()-499;
  
  // gap in RPC chamber 
  zazor = .03;
  // Sizes of RPC chamber 
  xp = 3.0; //small pixel
//xp = 3.9; //large pixel 
  yp = 12.3*0.05; // 5% X0 of glass 
  zp = 3.0; //small pixel
//zp = 4.1; //large pixel
  // Large not sensitive volumes with CO2 
  par[0] = xm/2;
  par[1] = ym/2;
  par[2] = zm0/2;
  gMC->Gsvolu("FBT1", "BOX ", idtmed[506], par, 3); // CO2
  gMC->Gspos("FBT1", 0, "FTO1", 0., 0., 0., 0, "ONLY");
  gMC->Gsdvn("FDT1", "FBT1", 2, 3); // 2 large modules along Z
  par[2] = zm1 / 2;
  gMC->Gsvolu("FBT2", "BOX ", idtmed[506], par, 3); // CO2
  gMC->Gspos("FBT2", 1, "FTO2", 0., 0., 0., 0, "ONLY");
  gMC->Gsdvn("FDT2", "FBT2", 2, 3); // 2 (PHOS) modules along Z
  par[2] = zm2 / 2;
  gMC->Gsvolu("FBT3", "BOX ", idtmed[506], par, 3); // CO2
  gMC->Gspos("FBT3", 2, "FTO3", 0., 0., 0., 0, "ONLY");
  gMC->Gsdvn("FDT3", "FBT3", 1, 3); // 1 (RICH) module along Z
  //
  //  subtraction of dead boundaries in X=2 cm and Z=7/2 cm 
  par[0] = par[0]-2.;
  Int_t nz0, nz1, nz2, nx; //- numbers of pixels
  nx = Int_t (par[0]*2/xp);
  cout <<"************************* TOF geometry **************************"<<endl;
  cout<< "nx = "<< nx << " x size = "<< par[0]*2/nx << endl;
  par[1] = -1;
  par[2] = (zm0 / 2.)/2.; //this is half size of module after division by 2
  par[2]=par[2]-7/2.;
  nz0 = Int_t (par[2]*2/zp);
cout<< "nz0 = "<< nz0 << " z0 size = "<< par[2]*2/nz0 << endl;
  gMC->Gsvolu("FLT1", "BOX ", idtmed[506], par, 3); // CO2
  gMC->Gspos("FLT1", 0, "FDT1", 0., 0., 0., 0, "ONLY");
  par[2] = (zm1 / 2.)/2.; //this is half size of module after division by 2
  par[2]=par[2]-7/2.;
  nz1 = Int_t (par[2]*2/zp);
cout<< "nz1 = "<< nz1 << " z1 size = "<< par[2]*2/nz1 << endl;
  gMC->Gsvolu("FLT2", "BOX ", idtmed[506], par, 3); // CO2
  gMC->Gspos("FLT2", 0, "FDT2", 0., 0., 0., 0, "ONLY");
  par[2] = (zm2 / 2.); //this is half size of module after division by 1
  par[2]=par[2]-7/2.;
  nz2 = Int_t (par[2]*2/zp);
cout<< "nz2 = "<< nz2 << " z2 size = "<< par[2]*2/nz2 << endl;
  gMC->Gsvolu("FLT3", "BOX ", idtmed[506], par, 3); // CO2
  gMC->Gspos("FLT3", 0, "FDT3", 0., 0., 0., 0, "ONLY");
  //
////////// Layers before detector ////////////////////
// Mylar layer in front 0.5mm thick at the beginning
  par[0] = -1;
  par[1] = 0.05 / 2;
  par[2] = -1;
  ycoor = -ym/2 + par[1];
  gMC->Gsvolu("FMY1", "BOX ", idtmed[511], par, 3); // Mylar
  gMC->Gspos("FMY1", 0, "FLT1", 0., ycoor, 0., 0, "ONLY");
  gMC->Gsvolu("FMY2", "BOX ", idtmed[511], par, 3); // Mylar
  gMC->Gspos("FMY2", 0, "FLT2", 0., ycoor, 0., 0, "ONLY");
  gMC->Gsvolu("FMY3", "BOX ", idtmed[511], par, 3); // Mylar
  gMC->Gspos("FMY3", 0, "FLT3", 0., ycoor, 0., 0, "ONLY");
// Honeycomb layer (1cm of special!!! polyethilene)
  ycoor = ycoor + par[1];
  par[0] = -1;
  par[1] = 1. / 2;
  par[2] = -1;
  ycoor = ycoor + par[1];
  gMC->Gsvolu("FPL1", "BOX ", idtmed[503], par, 3); // Hony
  gMC->Gspos("FPL1", 0, "FLT1", 0., ycoor, 0., 0, "ONLY");
  gMC->Gsvolu("FPL2", "BOX ", idtmed[503], par, 3); // Hony
  gMC->Gspos("FPL2", 0, "FLT2", 0., ycoor, 0., 0, "ONLY");
  gMC->Gsvolu("FPL3", "BOX ", idtmed[503], par, 3); // Hony
  gMC->Gspos("FPL3", 0, "FLT3", 0., ycoor, 0., 0, "ONLY");
  //
///////////////// Detector itself //////////////////////
  par[0] = -1;
  par[1] = yp/2; // 5 %X0 thick of glass  
  par[2] = -1;
  ycoor = -ym/2 + 2;
  gMC->Gsvolu("FLD1", "BOX ", idtmed[514], par, 3); // Glass
  gMC->Gspos("FLD1", 0, "FLT1", 0., ycoor, 0., 0, "ONLY");
  gMC->Gsvolu("FLD2", "BOX ", idtmed[514], par, 3); // Glass
  gMC->Gspos("FLD2", 0, "FLT2", 0., ycoor, 0., 0, "ONLY");
  gMC->Gsvolu("FLD3", "BOX ", idtmed[514], par, 3); // Glass
  gMC->Gspos("FLD3", 0, "FLT3", 0., ycoor, 0., 0, "ONLY");
  //
  gMC->Gsdvn("FLZ1", "FLD1", nz0, 3); //pixel size xp=zp=3
  gMC->Gsdvn("FLZ2", "FLD2", nz1, 3); 
  gMC->Gsdvn("FLZ3", "FLD3", nz2, 3); 
  gMC->Gsdvn("FLX1", "FLZ1", nx, 1); 
  gMC->Gsdvn("FLX2", "FLZ2", nx, 1); 
  gMC->Gsdvn("FLX3", "FLZ3", nx, 1); 
  // RPC pixel itself 
  par[0] = -1;//xp/2;
  par[1] = -1;//yp/2; // 5 %X0 thick of glass  
  par[2] = -1;//zp/2;
  gMC->Gsvolu("FPA0", "BOX ", idtmed[514], par, 3);// Glass
  gMC->Gspos("FPA0", 1, "FLX1", 0., 0., 0., 0, "ONLY");
  gMC->Gspos("FPA0", 2, "FLX2", 0., 0., 0., 0, "ONLY");
  gMC->Gspos("FPA0", 3, "FLX3", 0., 0., 0., 0, "ONLY");
  // Freon gas sencitive volume
  par[0] = -1;
  par[1] = zazor/2;
  par[2] = -1;
  gMC->Gsvolu("FPG0", "BOX ", idtmed[513], par, 3);// Freon 
  gMC->Gspos("FPG0", 0, "FPA0", 0., 0., 0., 0, "ONLY");
  //
////////// Layers after detector ////////////////////
  // Honeycomb layer after (3cm)
  par[0] = -1;
  par[1] = 1.2 / 2.;
  par[2] = -1;
  ycoor = -ym/2 + 6. - par[1];
  gMC->Gsvolu("FPE1", "BOX ", idtmed[503], par, 3); // Hony
  gMC->Gspos("FPE1", 0, "FLT1", 0., ycoor, 0., 0, "ONLY");
  gMC->Gsvolu("FPE2", "BOX ", idtmed[503], par, 3); // Hony
  gMC->Gspos("FPE2", 0, "FLT2", 0., ycoor, 0., 0, "ONLY");
  gMC->Gsvolu("FPE3", "BOX ", idtmed[503], par, 3); // Hony
  gMC->Gspos("FPE3", 0, "FLT3", 0., ycoor, 0., 0, "ONLY");
  // Electronics (Cu) after
  par[0] = -1;
  par[1] = 1.43*0.05 / 2.; // 5% of X0
  par[2] = -1;
  ycoor = -ym/2 + 6.+par[1];
  gMC->Gsvolu("FEC1", "BOX ", idtmed[501], par, 3); // Cu
  gMC->Gspos("FEC1", 0, "FLT1", 0., ycoor, 0., 0, "ONLY");
  gMC->Gsvolu("FEC2", "BOX ", idtmed[501], par, 3); // Cu
  gMC->Gspos("FEC2", 0, "FLT2", 0., ycoor, 0., 0, "ONLY");
  gMC->Gsvolu("FEC3", "BOX ", idtmed[501], par, 3); // Cu
  gMC->Gspos("FEC3", 0, "FLT3", 0., ycoor, 0., 0, "ONLY");
 // Cooling water after
  ycoor = ycoor+par[1];
  par[0] = -1;
  par[1] = 36.1*0.02 / 2.; // 2% of X0
  par[2] = -1;
  ycoor = ycoor+par[1];
  gMC->Gsvolu("FWA1", "BOX ", idtmed[515], par, 3); // Water
  gMC->Gspos("FWA1", 0, "FLT1", 0., ycoor, 0., 0, "ONLY");
  gMC->Gsvolu("FWA2", "BOX ", idtmed[515], par, 3); // Water
  gMC->Gspos("FWA2", 0, "FLT2", 0., ycoor, 0., 0, "ONLY");
  gMC->Gsvolu("FWA3", "BOX ", idtmed[515], par, 3); // Water
  gMC->Gspos("FWA3", 0, "FLT3", 0., ycoor, 0., 0, "ONLY");
  //back plate honycomb (2cm)
  par[0] = -1;
  par[1] = 2 / 2.;
  par[2] = -1;
  ycoor = ym/2 - par[1];
  gMC->Gsvolu("FEG1", "BOX ", idtmed[503], par, 3); // Hony
  gMC->Gspos("FEG1", 0, "FLT1", 0., ycoor, 0., 0, "ONLY");
  gMC->Gsvolu("FEG2", "BOX ", idtmed[503], par, 3); // Hony
  gMC->Gspos("FEG2", 0, "FLT2", 0., ycoor, 0., 0, "ONLY");
  gMC->Gsvolu("FEG3", "BOX ", idtmed[503], par, 3); // Hony
  gMC->Gspos("FEG3", 0, "FLT3", 0., ycoor, 0., 0, "ONLY");
}

//_____________________________________________________________________________
void AliTOFv3::DrawModule()
{
  //
  // Draw a shaded view of the Time Of Flight version 0
  //
  // Set everything unseen
  gMC->Gsatt("*", "seen", -1);
  // 
  // Set ALIC mother transparent
  gMC->Gsatt("ALIC","SEEN",0);
  //
  // Set the volumes visible
  gMC->Gsatt("ALIC","SEEN",0);
  gMC->Gsatt("FBAR","SEEN",1);
  gMC->Gsatt("FTO1","SEEN",1);
  gMC->Gsatt("FTO2","SEEN",1);
  gMC->Gsatt("FTO3","SEEN",1);
  gMC->Gsatt("FBT1","SEEN",1);
  gMC->Gsatt("FBT2","SEEN",1);
  gMC->Gsatt("FBT3","SEEN",1);
  gMC->Gsatt("FDT1","SEEN",1);
  gMC->Gsatt("FDT2","SEEN",1);
  gMC->Gsatt("FDT3","SEEN",1);
  gMC->Gsatt("FLT1","SEEN",1);
  gMC->Gsatt("FLT2","SEEN",1);
  gMC->Gsatt("FLT3","SEEN",1);
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
    id=gMC->CurrentVolOffID(8,copy);
    vol[0]=copy;
    if(id==fIdFTO3) {
      vol[0]+=22;
      id=gMC->CurrentVolOffID(5,copy);
      if(id==fIdFLT3) vol[1]+=6;
    } else if (id==fIdFTO2) {
      vol[0]+=20;
      id=gMC->CurrentVolOffID(5,copy);
      if(id==fIdFLT2) vol[1]+=8;
    } else {
      id=gMC->CurrentVolOffID(5,copy);
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
