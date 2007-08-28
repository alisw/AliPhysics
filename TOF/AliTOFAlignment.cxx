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
***************************************************************************/

/*
$Log$
Revision 1.17  2007/06/06 16:26:46  arcelli
remove fall-back call to local CDB storage

Revision 1.16  2007/05/15 16:25:44  cvetan
Moving the alignment-related static methods from AliAlignObj to the new geometry steering class AliGeomManager (macro from Raffaele)

Revision 1.15  2007/05/03 09:25:10  decaro
Coding convention: RN13 violation -> suppression

Revision 1.14  2007/04/18 14:49:54  arcelli
Some code cleanup, added more debug info

Revision 1.13  2007/04/17 16:38:36  arcelli
Include Methods to derive TOF AlignObjs from Survey Data

Revision 1.12  2007/02/28 18:09:23  arcelli
Add protection against failed retrieval of the CDB cal object

Revision 1.11  2006/09/19 14:31:26  cvetan
Bugfixes and clean-up of alignment object classes. Introduction of so called symbolic names used to identify the alignable volumes (Raffaele and Cvetan)

Revision 1.10  2006/08/22 13:26:05  arcelli
removal of effective c++ warnings (C.Zampolli)

Revision 1.9  2006/08/10 14:46:54  decaro
TOF raw data format: updated version

Revision 1.8  2006/05/04 19:41:42  hristov
Possibility for partial TOF geometry (S.Arcelli)

Revision 1.7  2006/04/27 13:13:29  hristov
Moving the destructor to the implementation file

Revision 1.6  2006/04/20 22:30:49  hristov
Coding conventions (Annalisa)

Revision 1.5  2006/04/16 22:29:05  hristov
Coding conventions (Annalisa)

Revision 1.4  2006/04/05 08:35:38  hristov
Coding conventions (S.Arcelli, C.Zampolli)

Revision 1.3  2006/03/31 13:49:07  arcelli
Removing some junk printout

Revision 1.2  2006/03/31 11:26:30  arcelli
 changing CDB Ids according to standard convention

Revision 1.1  2006/03/28 14:54:48  arcelli
class for TOF alignment

author: Silvia Arcelli, arcelli@bo.infn.it
*/  

/////////////////////////////////////////////////////////
//                                                     //
//            Class for alignment procedure            //
//                                                     //
//                                                     //
//                                                     //
/////////////////////////////////////////////////////////

#include <Rtypes.h>

#include "TMath.h"
#include "TFile.h"
#include "TRandom.h"

#include "AliLog.h"
#include "AliAlignObj.h"
#include "AliAlignObjParams.h"
#include "AliAlignObjMatrix.h"
#include "AliCDBManager.h"
#include "AliCDBMetaData.h"
#include "AliCDBId.h"
#include "AliCDBEntry.h"
#include "AliTOFAlignment.h"


ClassImp(AliTOFAlignment)
const Double_t AliTOFAlignment::fgkXsizeTOF  = 124.5; // x size of the TOF ext. volume, cm
const Double_t AliTOFAlignment::fgkYsizeTOF  = 29.0;  // y size of the TOF ext. volume, cm
const Double_t AliTOFAlignment::fgkZsizeTOF  = 913.8; // z size of the TOF ext. volume, cm
const Double_t AliTOFAlignment::fgkRorigTOF  = 384.5; // Mean Radius of the TOF ext. volume, cm
const Double_t AliTOFAlignment::fgkXFM = 38.0; //x pos of FM in the LRS, cm 
const Double_t AliTOFAlignment::fgkYFM = 11.2; //y pos of FM in the LRS, cm
const Double_t AliTOFAlignment::fgkZFM = 457.3;//z pos of FM in the LRS, cm
const Double_t AliTOFAlignment::fgkZsizeTOFSens=741.2; //z size of the TOF sensitive volume, cm

//_____________________________________________________________________________
AliTOFAlignment::AliTOFAlignment():
  TTask("AliTOFAlignment",""),
  fNTOFAlignObj(0),
  fTOFmgr(0x0),
  fTOFAlignObjArray(0x0)
 { 
   //AliTOFalignment main Ctor
   for(Int_t ism=0;ism<18;ism++){
     for(Int_t iFM=0;iFM<4;iFM++){
       for(Int_t iFMc=0;iFMc<3;iFMc++){
	 fTOFSurveyFM[ism][iFM][iFMc]=-1.;
       }
     }
   }
}
//_____________________________________________________________________________
AliTOFAlignment::AliTOFAlignment(const AliTOFAlignment &t):
  TTask("AliTOFAlignment",""),
  fNTOFAlignObj(0),
  fTOFmgr(0x0),
  fTOFAlignObjArray(0x0)
{ 
  //AliTOFAlignment copy Ctor

  fNTOFAlignObj=t.fNTOFAlignObj;
  fTOFAlignObjArray=t.fTOFAlignObjArray;
  //AliTOFalignment main Ctor
  for(Int_t iSM=0;iSM<18;iSM++){
     for(Int_t iFM=0;iFM<4;iFM++){
       for(Int_t iFMc=0;iFMc<3;iFMc++){
	 fTOFSurveyFM[iSM][iFM][iFMc]=-1.;
       }
     }
  }  
}
//_____________________________________________________________________________
AliTOFAlignment& AliTOFAlignment::operator=(const AliTOFAlignment &t){ 
  //AliTOFAlignment assignment operator

  this->fNTOFAlignObj=t.fNTOFAlignObj;
  this->fTOFmgr=t.fTOFmgr;
  this->fTOFAlignObjArray=t.fTOFAlignObjArray;
  return *this;

}
//_____________________________________________________________________________
AliTOFAlignment::~AliTOFAlignment() {
  delete fTOFAlignObjArray;
  delete fTOFmgr;
}

//_____________________________________________________________________________
void AliTOFAlignment::Smear( Float_t *tr, Float_t *rot)
{
  //Introduce Random Offset/Tilts
  fTOFAlignObjArray = new TObjArray(kMaxAlignObj);
  Float_t dx, dy, dz;  // shifts
  Float_t dpsi, dtheta, dphi; // angular displacements
  TRandom *rnd   = new TRandom(1567);
 
  Int_t nSMTOF = 18;
  AliGeomManager::ELayerID iLayer = AliGeomManager::kInvalidLayer;
  UShort_t iIndex=0; //dummy volume index
  //  AliGeomManager::ELayerID iLayer = AliGeomManager::kTOF;
  //  Int_t iIndex=1; //dummy volume index
  UShort_t dvoluid = AliGeomManager::LayerToVolUID(iLayer,iIndex); //dummy volume identity 
  Int_t i;
  for (i = 0; i<nSMTOF ; i++) {
    Char_t  path[100];
    sprintf(path,"/ALIC_1/B077_1/BSEGMO%i_1/BTOF%i_1",i,i);

    dx = (rnd->Gaus(0.,1.))*tr[0];
    dy = (rnd->Gaus(0.,1.))*tr[1];
    dz = (rnd->Gaus(0.,1.))*tr[2];
    dpsi   = rot[0];
    dtheta = rot[1];
    dphi   = rot[2];
    AliAlignObjParams *o =new AliAlignObjParams(path, dvoluid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
    fTOFAlignObjArray->Add(o);
  }

  fNTOFAlignObj=fTOFAlignObjArray->GetEntries();
  AliInfo(Form("Number of Alignable Volumes: %d",fNTOFAlignObj));
  delete rnd;
}

//_____________________________________________________________________________
void AliTOFAlignment::Align( Float_t *tr, Float_t *rot)
{
  //Introduce Offset/Tilts

  fTOFAlignObjArray = new TObjArray(kMaxAlignObj);
  Float_t dx, dy, dz;  // shifts
  Float_t dpsi, dtheta, dphi; // angular displacements


  Int_t nSMTOF = 18;
  AliGeomManager::ELayerID iLayer = AliGeomManager::kInvalidLayer;
  UShort_t iIndex=0; //dummy volume index
  UShort_t dvoluid = AliGeomManager::LayerToVolUID(iLayer,iIndex); //dummy volume identity 
  Int_t i;
  for (i = 0; i<nSMTOF ; i++) {

    Char_t  path[100];
    sprintf(path,"/ALIC_1/B077_1/BSEGMO%i_1/BTOF%i_1",i,i);
    dx = tr[0];
    dy = tr[1];
    dz = tr[2];
    dpsi   = rot[0];
    dtheta = rot[1];
    dphi   = rot[2];
    
    AliAlignObjParams *o =new AliAlignObjParams(path, dvoluid, dx, dy, dz, dpsi, dtheta, dphi, kTRUE);
    fTOFAlignObjArray->Add(o);
  }
  fNTOFAlignObj=fTOFAlignObjArray->GetEntries();
  AliInfo(Form("Number of Alignable Volumes: %d",fNTOFAlignObj));
}
//_____________________________________________________________________________
void AliTOFAlignment::WriteParOnCDB(Char_t *sel, Int_t minrun, Int_t maxrun)
{
  //Write Align Par on CDB
  AliCDBManager *man = AliCDBManager::Instance();
  Char_t *sel1 = "AlignPar" ;
  Char_t  out[100];
  sprintf(out,"%s/%s",sel,sel1); 
  AliCDBId idTOFAlign(out,minrun,maxrun);
  AliCDBMetaData *mdTOFAlign = new AliCDBMetaData();
  mdTOFAlign->SetResponsible("TOF");
  AliInfo(Form("Number of Alignable Volumes: %d",fNTOFAlignObj));
  man->Put(fTOFAlignObjArray,idTOFAlign,mdTOFAlign);
}
//_____________________________________________________________________________
void AliTOFAlignment::ReadParFromCDB(Char_t *sel, Int_t nrun)
{
  //Read Align Par from CDB
  AliCDBManager *man = AliCDBManager::Instance();
  Char_t *sel1 = "AlignPar" ;
  Char_t  out[100];
  sprintf(out,"%s/%s",sel,sel1); 
  AliCDBEntry *entry = man->Get(out,nrun);
  if (!entry) { 
    AliError(Form("Failed to get entry: %s",out));
    return; 
  }
  fTOFAlignObjArray=(TObjArray*)entry->GetObject();
  fNTOFAlignObj=fTOFAlignObjArray->GetEntries();
  AliInfo(Form("Number of Alignable Volumes from CDB: %d",fNTOFAlignObj));

}
//_____________________________________________________________________________
void AliTOFAlignment::WriteSimParOnCDB(Char_t *sel, Int_t minrun, Int_t maxrun)
{
  //Write Sim Align Par on CDB
  AliCDBManager *man = AliCDBManager::Instance();
  Char_t *sel1 = "AlignSimPar" ;
  Char_t  out[100];
  sprintf(out,"%s/%s",sel,sel1); 
  AliCDBId idTOFAlign(out,minrun,maxrun);
  AliCDBMetaData *mdTOFAlign = new AliCDBMetaData();
  mdTOFAlign->SetResponsible("TOF");
  AliInfo(Form("Number of Alignable Volumes: %d",fNTOFAlignObj));
  man->Put(fTOFAlignObjArray,idTOFAlign,mdTOFAlign);
}
//_____________________________________________________________________________
void AliTOFAlignment::ReadSimParFromCDB(Char_t *sel, Int_t nrun){
  //Read Sim Align Par from CDB
  AliCDBManager *man = AliCDBManager::Instance();
  Char_t *sel1 = "AlignSimPar" ;
  Char_t  out[100];
  sprintf(out,"%s/%s",sel,sel1); 
  AliCDBEntry *entry = man->Get(out,nrun);
  fTOFAlignObjArray=(TObjArray*)entry->GetObject();
  fNTOFAlignObj=fTOFAlignObjArray->GetEntries();
  AliInfo(Form("Number of Alignable Volumes from CDB: %d",fNTOFAlignObj));

}
//_____________________________________________________________________________
void AliTOFAlignment::WriteOnCDBforDC()
{
  //Write Align Par on CDB for DC06
  AliCDBManager *man = AliCDBManager::Instance();
  AliCDBId idTOFAlign("TOF/Align/Data",0,0);
  AliCDBMetaData *mdTOFAlign = new AliCDBMetaData();
  mdTOFAlign->SetComment("Alignment objects for ideal geometry, i.e. applying them to TGeo has to leave geometry unchanged");
  mdTOFAlign->SetResponsible("TOF");
  AliInfo(Form("Number of Alignable Volumes: %d",fNTOFAlignObj));
  man->Put(fTOFAlignObjArray,idTOFAlign,mdTOFAlign);
}
//_____________________________________________________________________________
void AliTOFAlignment::ReadFromCDBforDC()
{
  //Read Sim Align Par from CDB for DC06
  AliCDBManager *man = AliCDBManager::Instance();
  AliCDBEntry *entry = man->Get("TOF/Align/Data",0);
  fTOFAlignObjArray=(TObjArray*)entry->GetObject();
  fNTOFAlignObj=fTOFAlignObjArray->GetEntries();
  AliInfo(Form("Number of Alignable Volumes from CDB: %d",fNTOFAlignObj));

}
//_____________________________________________________________________________
void AliTOFAlignment::BuildGeomForSurvey()
{

  //Generates the ideal TOF structure with four Fiducial Marks in each 
  //supermodule (two on each z side) in their expected position. 
  //Highly inspired to Raffaele's example... 

  fTOFmgr = new TGeoManager("Geom","survey to alignment for TOF");
  TGeoMedium *medium = 0;
  TGeoVolume *top = fTOFmgr->MakeBox("TOP",medium,1000,1000,1000);
  fTOFmgr->SetTopVolume(top);
  // make shape components:  
  // This is the big box containing the TOF master sensitive volume+services  
  TGeoBBox *sbox0  = new TGeoBBox(fgkXsizeTOF*0.5,fgkYsizeTOF*0.5,fgkZsizeTOF*0.5);
  TGeoVolume* box0[18];
  // This is the big box containing the TOF master sensitive volume  
  TGeoBBox *sbox1  = new TGeoBBox(fgkXsizeTOF*0.5,fgkYsizeTOF*0.5,fgkZsizeTOFSens*0.5);
  TGeoVolume* box1 = new TGeoVolume("B1",sbox1);
  box1->SetLineColor(3);//green

  // Now four fiducial marks on SM, expressed in local coordinates
  // They are positioned at x=+/- 38 cm, y=11.2, z=+/- 456.94 cm

  TGeoBBox *fmbox  = new TGeoBBox(1,1,1);
  TGeoVolume* fm = new TGeoVolume("FM",fmbox);
  fm->SetLineColor(2);//color

  TGeoTranslation* mAtr = new TGeoTranslation("mAtr",-fgkXFM, fgkYFM ,fgkZFM);
  TGeoTranslation* mBtr = new TGeoTranslation("mBtr", fgkXFM, fgkYFM, fgkZFM);
  TGeoTranslation* mCtr = new TGeoTranslation("mCtr", fgkXFM, fgkYFM,-fgkZFM);
  TGeoTranslation* mDtr = new TGeoTranslation("mDtr",-fgkXFM, fgkYFM,-fgkZFM);

  // position all this stuff in the global ALICE frame

  char name[16];
  Double_t smX = 0.;
  Double_t smY = 0.;
  Double_t smZ = 0.;
  Float_t  smR = fgkRorigTOF;

  for (Int_t iSM = 0; iSM < 18; iSM++) {
    Int_t mod = iSM + 13;
    if (mod > 17) mod -= 18;
    sprintf(name, "BTOF%d",mod);
    box0[iSM] = new TGeoVolume(name,sbox0);
    Float_t phi  = iSM * 20.;
    Float_t phirot = 180 + phi;    
    smX =  TMath::Sin(phi*TMath::Pi()/180.)*smR;
    smY = -TMath::Cos(phi*TMath::Pi()/180.)*smR;
    smZ = 0.;
    TGeoRotation* smRot = new TGeoRotation("smRot",phirot,0,0.);    
    TGeoCombiTrans trans = *(new TGeoCombiTrans(smX,smY,smZ, smRot));
    TGeoMatrix* id = new TGeoHMatrix();
    TGeoHMatrix  transMat = *id * trans;
    TGeoHMatrix  *smTrans = new TGeoHMatrix(transMat);
    box0[iSM]->SetVisDaughters();
    box0[iSM]->SetLineColor(1); //black
    top->AddNode(box0[iSM],1,smTrans); //place the extended SM volume
    box0[iSM]->AddNode(box1,1); //place the inner SM volume
    box0[iSM]->AddNode(fm,1,mAtr);
    box0[iSM]->AddNode(fm,2,mBtr);
    box0[iSM]->AddNode(fm,3,mCtr);
    box0[iSM]->AddNode(fm,4,mDtr);
  }  

  fTOFmgr->CloseGeometry();
  fTOFmgr->GetTopVolume()->Draw();
  fTOFmgr->SetVisOption(0);
  fTOFmgr->SetVisLevel(6);

  // Now Store the "Ideal" Matrices for later use....

  for (Int_t iSM = 0; iSM < 18; iSM++) {

    sprintf(name, "TOP_1/BTOF%d_1", iSM);
    printf("\n\n*****************  TOF SuperModule:  %s ****************** \n",name);
    TGeoPhysicalNode* pn3 = fTOFmgr->MakePhysicalNode(name);
    fTOFMatrixId[iSM] = pn3->GetMatrix(); //save "ideal" global matrix 
    printf("\n\n***************  The Ideal Matrix in GRS *****************\n");
    fTOFMatrixId[iSM]->Print();

  }
}
//_____________________________________________________________________________
void AliTOFAlignment::InsertMisAlignment( Float_t *mis)
{
  // Now Apply the Displacements and store the misaligned FM positions...

  Double_t lA[3]={-fgkXFM,fgkYFM, fgkZFM};
  Double_t lB[3]={ fgkXFM,fgkYFM, fgkZFM};
  Double_t lC[3]={ fgkXFM,fgkYFM,-fgkZFM};
  Double_t lD[3]={-fgkXFM,fgkYFM,-fgkZFM};

  for(Int_t iSM=0;iSM<18;iSM++){
  // ************* get ideal global matrix *******************
    char name[16];
    sprintf(name, "TOP_1/BTOF%d_1", iSM);
    fTOFmgr->cd(name);
    printf("\n\n******Misaligning TOF SuperModule ************** %s \n",name);

  // ************* get ideal local matrix *******************
    TGeoHMatrix g3 = *fTOFmgr->GetCurrentMatrix(); 
    TGeoNode* n3 = fTOFmgr->GetCurrentNode();
    TGeoMatrix* l3 = n3->GetMatrix(); 

    Double_t gA[3], gB[3], gC[3], gD[3]; // ideal FM point coord., global RS
    g3.LocalToMaster(lA,gA);
    g3.LocalToMaster(lB,gB);
    g3.LocalToMaster(lC,gC);
    g3.LocalToMaster(lD,gD);


    // We apply a delta transformation to the surveyed vol to represent
    // its real position, given below by ng3 nl3, which differs from its
    // ideal position saved above in g3 and l3


    Double_t dx     = mis[0]; // shift along x
    Double_t dy     = mis[1]; // shift along y
    Double_t dz     = mis[2]; // shift along z
    Double_t dphi   = mis[3]; // rot around z
    Double_t dtheta = mis[4]; // rot around x'
    Double_t dpsi   = mis[5]; // rot around z'

    TGeoRotation* rrot = new TGeoRotation("rot",dphi,dtheta,dpsi);
    TGeoCombiTrans localdelta = *(new TGeoCombiTrans(dx,dy,dz, rrot));
  // new local matrix, representing real position
    TGeoHMatrix nlocal = *l3 * localdelta;
    TGeoHMatrix* nl3 = new TGeoHMatrix(nlocal);
    TGeoPhysicalNode* pn3 = fTOFmgr->MakePhysicalNode(name);

    pn3->Align(nl3);    //Align....    
    
    TGeoHMatrix* ng3 = pn3->GetMatrix(); //"real" global matrix, what survey sees 
    printf("\n\n*************  The Misaligned Matrix in GRS **************\n");
    ng3->Print();
    Double_t ngA[3], ngB[3], ngC[3], ngD[3];// real FM point coord., global RS
    ng3->LocalToMaster(lA,ngA);
    ng3->LocalToMaster(lB,ngB);
    ng3->LocalToMaster(lC,ngC);
    ng3->LocalToMaster(lD,ngD);    

    for(Int_t iFM=0;iFM<3;iFM++){
      fTOFSurveyFM[iSM][0][iFM]=ngA[iFM];
      fTOFSurveyFM[iSM][1][iFM]=ngB[iFM];
      fTOFSurveyFM[iSM][2][iFM]=ngC[iFM];
      fTOFSurveyFM[iSM][3][iFM]=ngD[iFM];
    }
  }
}

//_____________________________________________________________________________
void AliTOFAlignment::AlignFromSurvey()
{
  //From Survey data, derive the needed transformations to get the 
  //Alignment Objects. 
  //Again, highly "inspired" to Raffaele's example... 

  fTOFAlignObjArray = new TObjArray(kMaxAlignObj);
  Int_t index=0; //let all SM modules have index=0
  AliGeomManager::ELayerID layer = AliGeomManager::kInvalidLayer;
  UShort_t dvoluid = AliGeomManager::LayerToVolUID(layer,index); //dummy vol id 
  
  for(Int_t iSM=0;iSM<18;iSM++){

    printf("\n\n******Survey analysis for TOF SuperModule ************** %i \n",iSM);

    Double_t ngA[3], ngB[3], ngC[3], ngD[3];// real FM point coord., global RS
 
   // Get the 'realistic' input from the Survey Matrix
    for(Int_t iFM=0;iFM<3;iFM++){
      ngA[iFM]=   fTOFSurveyFM[iSM][0][iFM];
      ngB[iFM]=   fTOFSurveyFM[iSM][1][iFM];
      ngC[iFM]=   fTOFSurveyFM[iSM][2][iFM];
      ngD[iFM]=   fTOFSurveyFM[iSM][3][iFM];
    }

    // From the new fiducial marks coordinates derive back the
    // new global position of the surveyed volume
    //*** What follows is the actual survey-to-alignment procedure
    
    Double_t ab[3], bc[3], n[3];
    Double_t plane[4], s=1.;
    
    // first vector on the plane of the fiducial marks
    for(Int_t i=0;i<3;i++){
      ab[i] = (ngB[i] - ngA[i]);
    }
    
    // second vector on the plane of the fiducial marks
    for(Int_t i=0;i<3;i++){
      bc[i] = (ngC[i] - ngB[i]);
    }
    
    // vector normal to the plane of the fiducial marks obtained
    // as cross product of the two vectors on the plane d0^d1
    n[0] = (ab[1] * bc[2] - ab[2] * bc[1]);
    n[1] = (ab[2] * bc[0] - ab[0] * bc[2]);
    n[2] = (ab[0] * bc[1] - ab[1] * bc[0]);
    
    Double_t sizen = TMath::Sqrt( n[0]*n[0] + n[1]*n[1] + n[2]*n[2] );
    if(sizen>1.e-8){
      s = Double_t(1.)/sizen ; //normalization factor
    }else{
      AliInfo("Problem in normalizing the vector");
    }
    
    // plane expressed in the hessian normal form, see:
    // http://mathworld.wolfram.com/HessianNormalForm.html
    // the first three are the coordinates of the orthonormal vector
    // the fourth coordinate is equal to the distance from the origin
  
    for(Int_t i=0;i<3;i++){
      plane[i] = n[i] * s;
    }
    plane[3] = ( plane[0] * ngA[0] + plane[1] * ngA[1] + plane[2] * ngA[2] );
    
    // The center of the square with fiducial marks as corners
    // as the middle point of one diagonal - md
    // Used below to get the center - orig - of the surveyed box

    Double_t orig[3], md[3];
    for(Int_t i=0;i<3;i++){
      md[i] = (ngA[i] + ngC[i]) * 0.5;
    }
    
    // The center of the box, gives the global translation
    for(Int_t i=0;i<3;i++){
      orig[i] = md[i] - plane[i]*fgkYFM;
    }
    
    // get local directions needed to write the global rotation matrix
    // for the surveyed volume by normalising vectors ab and bc
    Double_t sx = TMath::Sqrt(ab[0]*ab[0] + ab[1]*ab[1] + ab[2]*ab[2]);
    if(sx>1.e-8){
      for(Int_t i=0;i<3;i++){
	ab[i] /= sx;
      }
    }
    Double_t sy = TMath::Sqrt(bc[0]*bc[0] + bc[1]*bc[1] + bc[2]*bc[2]);
    if(sy>1.e-8){
      for(Int_t i=0;i<3;i++){
	bc[i] /= sy;
      }
    }
    Double_t rot[9] = {ab[0],plane[0],bc[0],ab[1],plane[1],-bc[1],ab[2],plane[2],-bc[2]}; // the rotation matrix

    // the Aligned matrix for the current TOF SMS in the Global RS, as derived from Survey:
    TGeoHMatrix ng;
    ng.SetTranslation(orig);
    ng.SetRotation(rot);
    printf("\n\n**** The Misaligned Matrix in GRS, as from Survey data ***\n");
    ng.Print();

    // Calculate the delta transformation wrt Ideal geometry
    // (Should be gdelta.rot ==I and gdelta.tr=0 if no misalignment is applied.)
    printf("\n\n**** The ideal matrix ***\n");
    fTOFMatrixId[iSM]->Print(); 
    TGeoHMatrix gdelta =fTOFMatrixId[iSM]->Inverse();
    printf("\n\n**** The inverse of the ideal matrix ***\n");
    gdelta.Print(); 
    gdelta.MultiplyLeft(&ng);
    printf("\n\n**** The Delta Matrix in GRS, as from Survey data ***\n");
    gdelta.Print(); 
    
    // Now Write the Alignment Objects....
    TString symname(Form("TOF/sm%02d",iSM));
    AliAlignObjMatrix* o = new AliAlignObjMatrix(symname.Data(),dvoluid,gdelta,kTRUE);
    fTOFAlignObjArray->Add(o);
  }
  // saving TOF AligObjs from survey on a file, for the moment.. 
  fNTOFAlignObj=fTOFAlignObjArray->GetEntries();
  AliInfo(Form("Number of Alignable Volumes: %d",fNTOFAlignObj));
  TFile f("TOFAlignFromSurvey.root","RECREATE");
  f.cd();
  f.WriteObject(fTOFAlignObjArray,"TOFAlignObjs","kSingleKey");
  f.Close();
}
