///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Time Of Flight                                                           //
//  This class contains the functions for version 1 of the Time Of Flight    //
//  detector.                                                                //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliTOFv1Class.gif">
*/
//End_Html
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTOFv1.h"
#include "AliRun.h"
#include "AliConst.h"
 
ClassImp(AliTOFv1)
 
//_____________________________________________________________________________
AliTOFv1::AliTOFv1()
{
  //
  // Default constructor
  //
}
 
//_____________________________________________________________________________
AliTOFv1::AliTOFv1(const char *name, const char *title)
       : AliTOF(name,title)
{
  //
  // Standard constructor for version 1 of Time Of Flight
  //
}
 
//_____________________________________________________________________________
void AliTOFv1::CreateGeometry()
{
  //
  // Create geometry for version 1 of Time Of Flight
  // Authors :   Maxim Martemianov, Boris Zagreev (ITEP)   18/09/98 
  //Begin_Html
  /*
    <img src="picts/AliTOFv1.gif">
  */
  //End_Html
  //

  Float_t fil_rich;
  Int_t lmax;
  Float_t phos_phi, zcor2, zcor3, ztof0, ztof1, ztof2;
  Float_t zl, phos_r, zazor;
  Int_t idrotm[101];
  Float_t phos_x;
  Float_t rp1, rp2;
  Float_t par[10], fil_min, fil_max, ysz, fil0;
  //
  Int_t *idtmed = fIdtmed->GetArray()-499;
  //
  // barrel size along Z axis 
  //
  // Temporary fix TOF people should really check this!!
  //  rp1 = 360.;
  // rp2 = 372.;
  rp1 = 370;
  rp2 = rp1 + 12;
  zl  = 720.;
  //
  // TOF width along radius of barrel 
  //xtof  = rp2 - rp1;
  ztof0 = 350.;
  ztof1 = 200.;
  ztof2 = 150.;
  //
  // Plate width 
  ysz = .3;
  //
  // DME barrel width
  zazor = 0.03;
  //
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
  //
  // barrel radius in ALIC 
  par[0] = rp1;
  par[1] = rp2;
  par[2] = zl / 2.;
  gMC->Gsvolu("FBAR", "TUBE", idtmed[500], par, 3);
  gMC->Gspos("FBAR", 1, "ALIC", 0., 0., 0., 0, "ONLY");
  gMC->Gsatt("FBAR", "SEEN", 0);
  //
  // --- First Block
  par[0] = (rp1+rp2-zazor)/2.;
  par[1] = (rp1+rp2+zazor)/2.;
  par[2] = ztof0;
  par[3] = 90. - fil_min;
  par[4] = 90. - fil_rich;
  fil0 = 180. - (par[3] + par[4]);
  //
  // --- Sensitive volume 
  gMC->Gsvolu("FBT1", "TUBS", idtmed[509], par, 5);
  AliMatrix(idrotm[1], 90., fil0, 90., fil0 + 90., 0., 0.);
  gMC->Gspos("FBT1", 0, "FBAR", 0., 0., 0., 0, "ONLY");
  gMC->Gspos("FBT1", 1, "FBAR", 0., 0., 0., idrotm[1], "ONLY");
  //
  // ALUMINA
  par[0] = (rp1+rp2+zazor)/2.;
  par[1] = (rp1+rp2+zazor)/2.+ysz;
  gMC->Gsvolu("FPE1","TUBS",idtmed[507], par, 0);
  gMC->Gsposp("FPE1",2,"FBAR", 0., 0., 0., 0,         "ONLY", par, 5);
  gMC->Gsposp("FPE1",3,"FBAR", 0., 0., 0., idrotm[1], "ONLY", par, 5);
  //
  par[0] = (rp1+rp2-zazor)/2.-ysz;
  par[1] = (rp1+rp2-zazor)/2.;
  gMC->Gsposp("FPE1",4,"FBAR", 0., 0., 0., 0,         "ONLY", par, 5);
  gMC->Gsposp("FPE1",5,"FBAR", 0., 0., 0., idrotm[1], "ONLY", par, 5);
  // --- Second block 
  par[0] = (rp1+rp2-zazor)/2.;
  par[1] = (rp1+rp2+zazor)/2.;
  par[2] = ztof1 / 2.;
  par[3] = 90. - fil_max;
  par[4] = 90. - fil_min;
  gMC->Gsvolu("FBT2", "TUBS", idtmed[509], par, 5);
  gMC->Gspos("FBT2", 0, "FBAR", 0., 0., zcor2, 0, "ONLY");
  gMC->Gspos("FBT2", 1, "FBAR", 0., 0., -zcor2, 0, "ONLY");
  //
  par[0]=(rp1+rp2+zazor)/2.;
  par[1]=(rp1+rp2+zazor)/2.+ysz;
  gMC->Gsvolu("FPE2","TUBS",idtmed[507], par, 0);
  gMC->Gsposp("FPE2",2,"FBAR",0.,0.,zcor2,0,"ONLY",par,5);
  gMC->Gsposp("FPE2",3,"FBAR",0.,0.,-zcor2,0,"ONLY",par,5);
  //
  par[0]=(rp1+rp2-zazor)/2.-ysz;
  par[1]=(rp1+rp2-zazor)/2.;
  gMC->Gsposp("FPE2",4,"FBAR",0.,0.,zcor2,0,"ONLY",par,5);
  gMC->Gsposp("FPE2",5,"FBAR",0.,0.,-zcor2,0,"ONLY",par,5);
  //
  // --- Third block 
  par[0] = (rp1+rp2-zazor)/2.;
  par[1] = (rp1+rp2-zazor)/2.;
  par[2] = ztof2 / 2.;
  par[3] = 90. - fil_rich;
  par[4] = fil_rich + 90.;
  gMC->Gsvolu("FBT3", "TUBS", idtmed[509], par, 5);
  gMC->Gspos("FBT3", 0, "FBAR", 0., 0., zcor3, 0, "ONLY");
  gMC->Gspos("FBT3", 1, "FBAR", 0., 0., -zcor3, 0, "ONLY");
  //
  par[0]=(rp1+rp2+zazor)/2.;
  par[1]=(rp1+rp2+zazor)/2.+ysz;
  gMC->Gsvolu("FPE3","TUBS",idtmed[507], par, 0);
  gMC->Gsposp("FPE3",2,"FBAR",0.,0.,zcor3,0,"ONLY",par,5);
  gMC->Gsposp("FPE3",3,"FBAR",0.,0.,-zcor3,0,"ONLY",par,5);
  //
  par[0]=(rp1+rp2-zazor)/2.-ysz;
  par[1]=(rp1+rp2-zazor)/2.;
  gMC->Gsposp("FPE3",4,"FBAR",0.,0.,zcor3,0,"ONLY",par,5);
  gMC->Gsposp("FPE3",5,"FBAR",0.,0.,-zcor3,0,"ONLY",par,5);
  //
}
 
//_____________________________________________________________________________
void AliTOFv1::DrawModule()
{
  //
  // Draw a shaded view of the Time Of Flight version 1
  //

  // Set everything unseen
  gMC->Gsatt("*", "seen", -1);
  // 
  // Set ALIC mother transparent
  gMC->Gsatt("ALIC","SEEN",0);
  //
  // Set the volumes visible
  gMC->Gsatt("FBAR","SEEN",0);
  gMC->Gsatt("FPE1","SEEN",1);
  gMC->Gsatt("FPE2","SEEN",1);
  gMC->Gsatt("FPE3","SEEN",1);
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

//___________________________________________
void AliTOFv1::CreateMaterials()
{
  //
  // Define materials for Time Of Flight
  //
  AliTOF::CreateMaterials();
}
 
//______________________________________________________________________________
void AliTOFv1::Init()
{
  //
  // Initialise Time Of Flight after it has been built

  AliTOF::Init();
  fIdFBT2=gMC->VolId("FBT2");
  fIdFBT3=gMC->VolId("FBT3");
}
 
//______________________________________________________________________________
void AliTOFv1::StepManager()
{
  TLorentzVector mom, pos;
  Float_t hits[8];
  Int_t vol[3];
  Int_t copy, id, i;
  Int_t *idtmed = fIdtmed->GetArray()-499;
  if(gMC->GetMedium()==idtmed[510-1] && 
     gMC->IsTrackEntering() && gMC->TrackCharge()
     && (id=gMC->CurrentVolID(copy))==fIdSens) {
    TClonesArray &lhits = *fHits;
//
// Record only charged tracks at entrance
    vol[2]=copy;
    vol[1]=gMC->CurrentVolOffID(1,copy);
    if(id==fIdFBT2) copy+=2; else 
      if(id==fIdFBT2) copy+=4;
    vol[0]=1;
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
 
