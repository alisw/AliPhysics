///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Time Of Flight                                                           //
//  This class contains the functions for version 0 of the Time Of Flight    //
//  detector.                                                                //
//                                                                           //
//Begin_Html
/*
<img src="gif/AliTOFv0Class.gif">
*/
//End_Html
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTOFv0.h"
#include "AliRun.h"
#include "AliMC.h"
#include "AliConst.h"
  
ClassImp(AliTOFv0)
 
//_____________________________________________________________________________
AliTOFv0::AliTOFv0() : AliTOF()
{
  //
  // Default constructor
  //
}
 
//_____________________________________________________________________________
AliTOFv0::AliTOFv0(const char *name, const char *title)
       : AliTOF(name,title)
{
  //
  // Standard constructor for version 0 of the Time Of Flight
  //
}
 
//_____________________________________________________________________________
void AliTOFv0::CreateGeometry()
{
  //
  // Definition of Geometry
  // Authors :   Maxim Martemianov, Boris Zagreev (ITEP)   18/09/98 
  //Begin_Html
  /*
    <img src="gif/AliTOFv0.gif">
  */
  //End_Html
  //

  AliMC* pMC = AliMC::GetMC();
  
  Float_t fil_rich;
  Int_t lmax;
  Float_t phos_phi, zcor2, zcor3, ztof0, ztof1, ztof2;
  Float_t zl, phos_r;
  Int_t idrotm[101];
  Float_t phos_x;
  Float_t rp1, rp2;
  Float_t par[10], fil_min, fil_max, ysz, fil0;
  //
  Int_t *idtmed = gAlice->Idtmed();
  //
  // barrel size along Z axis 
  rp1 = 360.;
  rp2 = 372.;
  zl  = 720.;
  //
  // TOF width along radius of barrel 
  //xtof  = rp2 - rp1;
  ztof0 = 350.;
  ztof1 = 200.;
  ztof2 = 150.;
  //
  // Plate width 
  ysz = .6;
  // PHOS and RICH angles 
  phos_x = 214.6;
  phos_r = 467.;
  //phos_z = 260.;
  //rich_z = 472.5;
  fil_rich = 30.;
  lmax = 19;
  zcor2 = ztof0 - ztof1 / 2.;
  zcor3 = ztof0 - ztof2 / 2.;
  phos_phi = TMath::ATan(phos_x / (phos_r * 2.));
  fil_min = (kPI - phos_phi * 4.) * kRaddeg - 180. / lmax;
  fil_max = (phos_phi * 4. + kPI) * kRaddeg + 180. / lmax;
  // barrel radius in ALIC 
  par[0] = rp1;
  par[1] = rp2;
  par[2] = zl / 2.;
  pMC->Gsvolu("FBAR", "TUBE", idtmed[500], par, 3);
  pMC->Gspos("FBAR", 1, "ALIC", 0., 0., 0., 0, "ONLY");
  pMC->Gsatt("FBAR", "SEEN", 0);
  // First Block
  par[0] = (rp1+rp2-ysz)/2.;
  par[1] = (rp1+rp2+ysz)/2.;
  par[2] = ztof0;
  par[3] = 90. - fil_min;
  par[4] = 90. - fil_rich;
  fil0 = 180. - (par[3] + par[4]);
  pMC->Gsvolu("FBT1", "TUBS", idtmed[507], par, 5);
  AliMatrix(idrotm[1], 90., fil0, 90., fil0 + 90., 0., 0.);
  pMC->Gspos("FBT1", 0, "FBAR", 0., 0., 0., 0, "ONLY");
  pMC->Gspos("FBT1", 1, "FBAR", 0., 0., 0., idrotm[1], "ONLY");
  // --- Second block 
  par[2] = ztof1 / 2.;
  par[3] = 90. - fil_max;
  par[4] = 90. - fil_min;
  pMC->Gsvolu("FBT2", "TUBS", idtmed[507], par, 5);
  pMC->Gspos("FBT2", 0, "FBAR", 0., 0., zcor2, 0, "ONLY");
  pMC->Gspos("FBT2", 1, "FBAR", 0., 0.,-zcor2, 0, "ONLY");
  // --- Third block 
  par[2] = ztof2 / 2.;
  par[3] = 90. - fil_rich;
  par[4] = fil_rich + 90.;
  pMC->Gsvolu("FBT3", "TUBS", idtmed[507], par, 5);
  pMC->Gspos("FBT3", 0, "FBAR", 0., 0., zcor3, 0, "ONLY");
  pMC->Gspos("FBT3", 1, "FBAR", 0., 0., -zcor3, 0, "ONLY");
}
 
//_____________________________________________________________________________
void AliTOFv0::DrawModule()
{
  //
  // Draw a shaded view of the common part of the TOF geometry
  // for versions 2 and 3
  //

  AliMC* pMC = AliMC::GetMC();
  
  // Set everything unseen
  pMC->Gsatt("*", "seen", -1);
  // 
  // Set ALIC mother transparent
  pMC->Gsatt("ALIC","SEEN",0);
  //
  // Set the volumes visible
  pMC->Gsatt("FBAR","SEEN",0);
  pMC->Gsatt("FBT1","SEEN",1);
  pMC->Gsatt("FBT2","SEEN",1);
  pMC->Gsatt("FBT3","SEEN",1);
  //
  pMC->Gdopt("hide", "on");
  pMC->Gdopt("shad", "on");
  pMC->Gsatt("*", "fill", 7);
  //
  pMC->SetClipBox(".");
  pMC->SetClipBox("*", 0, 1000, -1000, 1000, -1000, 1000);
  pMC->DefaultRange();
  pMC->Gdraw("alic", 40, 30, 0, 12, 9.5, .02, .02);
  pMC->Gdhead(1111, "Time Of Flight");
  pMC->Gdman(18, 4, "MAN");
  pMC->Gdopt("hide","off");
}

//_____________________________________________________________________________
void AliTOFv0::CreateMaterials()
{
  //
  // Define materials for the Time Of Flight
  //
  AliTOF::CreateMaterials();
}
 
//_____________________________________________________________________________
void AliTOFv0::Init()
{
  //
  // Initialise detector after that it has been built
  //

  AliMC* pMC = AliMC::GetMC();
  
  AliTOF::Init();
  fIdFBT2=pMC->VolId("FBT2");
  fIdFBT3=pMC->VolId("FBT3");
}
 
//_____________________________________________________________________________
void AliTOFv0::StepManager()
{
  //
  // Procedure called at each step in the Time Of Flight
  Float_t hits[8];
  Int_t vol[3];
  Int_t copy, id;
  //
  // Get the pointer to the MonteCarlo
  AliMC *pMC= AliMC::GetMC();
  Int_t *idtmed = gAlice->Idtmed();
  if(pMC->GetMedium()==idtmed[510-1] && 
     pMC->TrackEntering() && pMC->TrackCharge()
     && (id=pMC->CurrentVol(0,copy))==fIdSens) {
    TClonesArray &lhits = *fHits;
    //
    // Record only charged tracks at entrance
    vol[2]=copy;
    vol[1]=pMC->CurrentVolOff(1,0,copy);
    if(id==fIdFBT2) copy+=2; else 
      if(id==fIdFBT2) copy+=4;
    vol[0]=1;
    pMC->TrackPosition(hits);
    pMC->TrackMomentum(&hits[3]);
    hits[7]=pMC->TrackTime();
    new(lhits[fNhits++]) AliTOFhit(fIshunt,gAlice->CurrentTrack(),vol,hits);
  }
}

 
