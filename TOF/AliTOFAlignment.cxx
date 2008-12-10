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
Revision 1.19  2007/10/02 09:46:08  arcelli
add methods to retrieve real survey data, and make some analysis (by B. Guerzoni)

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

#include "TGeoMatrix.h"
#include "TMath.h"
#include "TFile.h"
#include "TRandom.h"
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoBBox.h"
#include "TGeoTrd1.h"
#include "TGeoPhysicalNode.h"
#include "TGeoNode.h"
#include "TObjString.h"

#include "AliLog.h"
//#include "AliAlignObj.h"
#include "AliAlignObjParams.h"
#include "AliAlignObjMatrix.h"
#include "AliCDBManager.h"
#include "AliCDBMetaData.h"
#include "AliCDBId.h"
#include "AliCDBEntry.h"
#include "AliTOFAlignment.h"
#include "AliSurveyObj.h"
#include "AliSurveyPoint.h"

ClassImp(AliTOFAlignment)

const Double_t AliTOFAlignment::fgkRorigTOF  = 384.5; // Mean Radius of the TOF ext. volume, cm
const Double_t AliTOFAlignment::fgkX1BTOF = 124.5;    //x1 size of BTOF
const Double_t AliTOFAlignment::fgkX2BTOF = 134.7262; //x2 size of BTOF
const Double_t AliTOFAlignment::fgkYBTOF = 747.2;     //y size of BTOF
const Double_t AliTOFAlignment::fgkZBTOF = 29.0;      //z size of BTOF
const Double_t AliTOFAlignment::fgkXFM = 38.0;     //x pos of FM in BTOF, cm 
const Double_t AliTOFAlignment::fgkYFM = 457.3;    //y pos of FM in BTOF, cm
const Double_t AliTOFAlignment::fgkZFM = 11.2;     //z pos of FM in BTOF, cm

//_____________________________________________________________________________
AliTOFAlignment::AliTOFAlignment():
  TTask("AliTOFAlignment",""),
  fNTOFAlignObj(0),
  fTOFmgr(0x0),
  fTOFAlignObjArray(0x0)
 { 
   //AliTOFalignment main Ctor
   for(Int_t i=0; i<18;i++)
     for(Int_t j=0; j<5; j++)
       fNFMforSM[i][j]=0;
   for(Int_t i=0; i<72; i++)
    for (Int_t j=0; j<6; j++)
      fCombFMData[i][j]=0;
}
//_____________________________________________________________________________
AliTOFAlignment::AliTOFAlignment(const AliTOFAlignment &t):
  TTask(t),
  fNTOFAlignObj(t.fNTOFAlignObj),
  fTOFmgr(0x0),
  fTOFAlignObjArray(t.fTOFAlignObjArray)
{ 
  //AliTOFAlignment copy Ctor

  //AliTOFalignment main Ctor
  for(Int_t i=0; i<18;i++)
     for(Int_t j=0; j<5; j++)
       fNFMforSM[i][j]=t.fNFMforSM[i][j];
  for(Int_t i=0; i<72; i++)
    for (Int_t j=0; j<6; j++)
      fCombFMData[i][j]=t.fCombFMData[i][j]; 
}
//_____________________________________________________________________________
AliTOFAlignment& AliTOFAlignment::operator=(const AliTOFAlignment &t){ 
  //AliTOFAlignment assignment operator

  if (&t == this)
    return *this;

  TTask::operator=(t);
  fNTOFAlignObj=t.fNTOFAlignObj;
  fTOFmgr=t.fTOFmgr;
  fTOFAlignObjArray=t.fTOFAlignObjArray;
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
void AliTOFAlignment::WriteParOnCDB(const Char_t *sel, Int_t minrun, Int_t maxrun)
{
  //Write Align Par on CDB
  AliCDBManager *man = AliCDBManager::Instance();
  const Char_t *sel1 = "AlignPar" ;
  Char_t  out[100];
  sprintf(out,"%s/%s",sel,sel1); 
  AliCDBId idTOFAlign(out,minrun,maxrun);
  AliCDBMetaData *mdTOFAlign = new AliCDBMetaData();
  mdTOFAlign->SetResponsible("TOF");
  AliInfo(Form("Number of Alignable Volumes: %d",fNTOFAlignObj));
  man->Put(fTOFAlignObjArray,idTOFAlign,mdTOFAlign);
}
//_____________________________________________________________________________
void AliTOFAlignment::ReadParFromCDB(const Char_t *sel, Int_t nrun)
{
  //Read Align Par from CDB
  AliCDBManager *man = AliCDBManager::Instance();
  const Char_t *sel1 = "AlignPar" ;
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
void AliTOFAlignment::WriteSimParOnCDB(const Char_t *sel, Int_t minrun, Int_t maxrun)
{
  //Write Sim Align Par on CDB
  AliCDBManager *man = AliCDBManager::Instance();
  const Char_t *sel1 = "AlignSimPar" ;
  Char_t  out[100];
  sprintf(out,"%s/%s",sel,sel1); 
  AliCDBId idTOFAlign(out,minrun,maxrun);
  AliCDBMetaData *mdTOFAlign = new AliCDBMetaData();
  mdTOFAlign->SetResponsible("TOF");
  AliInfo(Form("Number of Alignable Volumes: %d",fNTOFAlignObj));
  man->Put(fTOFAlignObjArray,idTOFAlign,mdTOFAlign);
}
//_____________________________________________________________________________
void AliTOFAlignment::ReadSimParFromCDB(const Char_t *sel, Int_t nrun){
  //Read Sim Align Par from CDB
  AliCDBManager *man = AliCDBManager::Instance();
  const Char_t *sel1 = "AlignSimPar" ;
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
  //Make BTOF

  fTOFmgr = new TGeoManager("Geom","survey to alignment for TOF");
  TGeoMedium *medium = 0;
  TGeoVolume *top = fTOFmgr->MakeBox("TOP",medium,1000,1000,1000);
  fTOFmgr->SetTopVolume(top);
  // make shape components:  
  // This is the BTOF containing the FTOA  
  TGeoTrd1 *strd1  = new TGeoTrd1(fgkX1BTOF*0.5,fgkX2BTOF*0.5, fgkYBTOF*0.5,fgkZBTOF*0.5);
  TGeoVolume* trd1[18];

  // Now four fiducial marks on SM, expressed in local coordinates
  // They are positioned at x=+/- 38 cm, y=+/- 457.3 cm, z=11.2 cm
  
  TGeoBBox *fmbox  = new TGeoBBox(1,1,1);
  TGeoVolume* fm = new TGeoVolume("FM",fmbox);
  fm->SetLineColor(2);
  

  TGeoTranslation* mAtr = new TGeoTranslation("mAtr",-fgkXFM, -fgkYFM ,fgkZFM);
  TGeoTranslation* mBtr = new TGeoTranslation("mBtr",fgkXFM, -fgkYFM ,fgkZFM );
  TGeoTranslation* mCtr = new TGeoTranslation("mCtr",fgkXFM, fgkYFM ,fgkZFM );
  TGeoTranslation* mDtr = new TGeoTranslation("mDtr",-fgkXFM, fgkYFM ,fgkZFM );

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
    trd1[iSM] = new TGeoVolume(name,strd1);
    Float_t phi  = iSM * 20.;
    Float_t phi2 = 270 + phi;
    if (phi2 >= 360.) phi2 -= 360.;
    smX =  TMath::Sin(phi*TMath::Pi()/180.)*smR;
    smY = -TMath::Cos(phi*TMath::Pi()/180.)*smR;
    smZ = 0.;  
    TGeoRotation* bTOFRot = new TGeoRotation("bTOFRot",phi,90,0.);
    TGeoCombiTrans trans = *(new TGeoCombiTrans(smX,smY,smZ, bTOFRot));
    TGeoMatrix* id = new TGeoHMatrix();
    TGeoHMatrix  transMat = *id * trans;
    TGeoHMatrix  *smTrans = new TGeoHMatrix(transMat);
    
    trd1[iSM]->AddNode(fm,1,mAtr);        //place FM in BTOF
    trd1[iSM]->AddNode(fm,2,mBtr);
    trd1[iSM]->AddNode(fm,3,mCtr);
    trd1[iSM]->AddNode(fm,4,mDtr);
    top->AddNode(trd1[iSM],1,smTrans);    //place BTOF_iSM in ALICE
    trd1[iSM]->SetVisDaughters();
    trd1[iSM]->SetLineColor(iSM);         //black
    
  }  

  fTOFmgr->CloseGeometry();
  fTOFmgr->GetTopVolume()->Draw();
  fTOFmgr->SetVisOption(0);
  fTOFmgr->SetVisLevel(6);

  // Now Store the "Ideal"  Global Matrices (local to global) for later use
  
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
void AliTOFAlignment::InsertMisAlignment(Float_t *mis)
{
  // Now Apply the Displacements and store the misaligned FM positions...
  //
  //

  Double_t lA[3]={-fgkXFM, -fgkYFM ,fgkZFM};
  Double_t lB[3]={fgkXFM, -fgkYFM ,fgkZFM};
  Double_t lC[3]={fgkXFM, fgkYFM ,fgkZFM};
  Double_t lD[3]={-fgkXFM, fgkYFM ,fgkZFM};

  for(Int_t iSM=0;iSM<18;iSM++){
     char name[16];
     sprintf(name, "TOP_1/BTOF%d_1", iSM);
     fTOFmgr->cd(name);
     printf("\n\n******Misaligning TOF SuperModule ************** %s \n",name);

    // ************* get ideal global matrix *******************
    TGeoHMatrix g3 = *fTOFmgr->GetCurrentMatrix(); 
    AliInfo(Form("This is the ideal global trasformation of SM %i",iSM));
    g3.Print(); // g3 is the local(BTOF) to global (ALICE) matrix and is the same of fTOFMatrixId
    TGeoNode* n3 = fTOFmgr->GetCurrentNode(); 
    TGeoMatrix* l3 = n3->GetMatrix(); 
    
    Double_t gA[3], gB[3], gC[3], gD[3]; // ideal global FM point coord.
    g3.LocalToMaster(lA,gA);
    g3.LocalToMaster(lB,gB);
    g3.LocalToMaster(lC,gC);
    g3.LocalToMaster(lD,gD);
   
    //  We apply a delta transformation to the surveyed vol to represent
    //  its real position, given below by ng3 nl3, which differs from its
    //  ideal position saved above in g3 and l3

    //we have to express the displacements as regards the old local RS (non misaligned BTOF)
    Double_t dx     = mis[0]; // shift along x 
    Double_t dy     = mis[1]; // shift along y 
    Double_t dz     = mis[2]; // shift along z 
    Double_t dphi   = mis[3]; // rot around z 
    Double_t dtheta = mis[4]; // rot around x' 
    Double_t dpsi   = mis[5]; // rot around z''

    TGeoRotation* rrot = new TGeoRotation("rot",dphi,dtheta,dpsi);
    TGeoCombiTrans localdelta = *(new TGeoCombiTrans(dx,dy,dz, rrot));
    AliInfo(Form("This is the local delta trasformation for SM %i \n",iSM));
    localdelta.Print();
    TGeoHMatrix nlocal = *l3 * localdelta;
    TGeoHMatrix* nl3 = new TGeoHMatrix(nlocal); // new matrix, representing real position (from new local mis RS to the global one)
   
    TGeoPhysicalNode* pn3 = fTOFmgr->MakePhysicalNode(name);

    pn3->Align(nl3);   
    
    TGeoHMatrix* ng3 = pn3->GetMatrix(); //"real" global matrix, what survey sees 
    printf("\n\n*************  The Misaligned Matrix in GRS **************\n");
    ng3->Print();
    Double_t ngA[3], ngB[3], ngC[3], ngD[3];// real FM point coord., global RS 
    ng3->LocalToMaster(lA,ngA);
    ng3->LocalToMaster(lB,ngB);
    ng3->LocalToMaster(lC,ngC);
    ng3->LocalToMaster(lD,ngD);    

    for(Int_t coord=0;coord<3;coord++){
      fCombFMData[iSM*4][2*coord]=ngA[coord];
      fCombFMData[iSM*4][2*coord+1]=1;
      fCombFMData[iSM*4+1][2*coord]=ngB[coord];
      fCombFMData[iSM*4+1][2*coord+1]=1;
      fCombFMData[iSM*4+2][2*coord]=ngC[coord];
      fCombFMData[iSM*4+2][2*coord+1]=1;
      fCombFMData[iSM*4+3][2*coord]=ngD[coord];
      fCombFMData[iSM*4+3][2*coord+1]=1;
      }
    }

}

//____________________________________________________________________________
void AliTOFAlignment::WriteCombData(const Char_t *nomefile, Int_t option)
{
  // 1 for simulated data; 0 for data from survey file
  // write combined data on a file
  //

  FILE *data;
  /* Open file in text mode: */
  if( (data = fopen( nomefile, "w+t" )) != NULL ){
    if (option==1){
      fprintf( data, "simulated data\n" );} else {
	fprintf( data, "survey data\n" );}
    if (option==1){
      fprintf( data, "data from InsertMisAlignmentBTOF method\n");}
    else {fprintf( data, "real survey data from text file (coordinate in global RS)\n");}
    fprintf( data, "Point Name,XPH,YPH,ZPH,PrecisionX(mm),PrecisionY(mm),PrecisionZ(mm)\n");
    fprintf( data, "> Data:\n");
    for(Int_t i=0;i<72;i++){
      if (fCombFMData[i][0]!=0){
	fprintf( data, "SM%02iFM%i %f %f %f M Y %f %f %f\n", (i-i%4)/4, i%4, fCombFMData[i][0],fCombFMData[i][2],fCombFMData[i][4],fCombFMData[i][1]*10,fCombFMData[i][3]*10,fCombFMData[i][5]*10); 
      }
    }
    fclose( data );
   }
  else{
    printf(  "Problem opening the file\n" );
  }
  
  return;  
}

//____________________________________________________________________________
void AliTOFAlignment::WriteSimSurveyData(const Char_t *nomefile)
{
  // write sim data in standard format
  //
  //

  FILE *data;
  /* Open file in text mode: */
  if( (data = fopen( nomefile, "w+t" )) != NULL )
   {
      fprintf( data, "> Title:\n" );
      fprintf( data, "simulated data\n" );
      fprintf( data, "> Date:\n" );
      fprintf( data, "24.09.2007\n" );
      fprintf( data, "> Subdetector:\n" );
      fprintf( data, "TOF\n" );
      fprintf( data, "> Report URL:\n" );
      fprintf( data, "https://edms.cern.ch/document/835615\n" );
      fprintf( data, "> Version:\n" );
      fprintf( data, "1\n");
      fprintf( data, "> General Observations:\n"); 
      fprintf( data, "data from InsertMisAlignmentBTOF method\n");
      fprintf( data, "> Coordinate System:\n");
      fprintf( data, "\\ALICEPH\n");
      fprintf( data, "> Units:\n");
      fprintf( data, "cm\n");
      fprintf( data, "> Nr Columns:\n");
      fprintf( data, "9\n");
      fprintf( data, "> Column Names:\n");
      fprintf( data, "Point Name,XPH,YPH,ZPH,Point Type,Target Used,PrecisionX(mm),PrecisionY(mm),PrecisionZ(mm)\n");
      fprintf( data, "> Data:\n");
      for(Int_t i=0;i<72;i++)
        if (fCombFMData[i][0]!=0)
 	  fprintf( data, "SM%02iFM%i %f %f %f M Y %f %f %f\n", (i-i%4)/4, i%4, fCombFMData[i][0],fCombFMData[i][2],fCombFMData[i][4],fCombFMData[i][1],fCombFMData[i][3],fCombFMData[i][5]); 
      
       fclose( data );
   }
   else
     printf(  "Problem opening the file\n" );
}

//____________________________________________________________________________
void AliTOFAlignment::MakeDefData(const Int_t nf,TString namefiles[])
{
  //this method combines survey data from different files (namefiles[]) 
  //
  // 
 
  Float_t data[72][6][100];
  for (Int_t i=0;i<72;i++)
    for (Int_t j=0; j<6; j++)
      for(Int_t k=0; k<100; k++)
        data[i][j][k]=0;
  Int_t nfm=0;
  Int_t nsm=0;
  Long64_t totdata[72]={0};

  for (Int_t ii=0;ii<nf; ii++)
    {
      AliSurveyObj *so = new AliSurveyObj();
      const Char_t *nome=namefiles[ii];
      so->FillFromLocalFile(nome);
      TObjArray *points = so->GetData();
      Int_t nSurveyPoint=points->GetEntries();
      for(Int_t jj=0;jj<nSurveyPoint;jj++){
        const char* pointName= ((AliSurveyPoint *) points->At(jj))->GetPointName().Data();
        nfm=atoi(&pointName[6]);
        nsm=atoi(&pointName[2]);
        data[nsm*4+nfm][0][totdata[nsm*4+nfm]]=((AliSurveyPoint *) points->At(jj))->GetX();
        data[nsm*4+nfm][2][totdata[nsm*4+nfm]]=((AliSurveyPoint *) points->At(jj))->GetY();
        data[nsm*4+nfm][4][totdata[nsm*4+nfm]]=((AliSurveyPoint *) points->At(jj))->GetZ();
        data[nsm*4+nfm][1][totdata[nsm*4+nfm]]=((AliSurveyPoint *) points->At(jj))->GetPrecisionX();
        data[nsm*4+nfm][3][totdata[nsm*4+nfm]]=((AliSurveyPoint *) points->At(jj))->GetPrecisionY();
        data[nsm*4+nfm][5][totdata[nsm*4+nfm]]=((AliSurveyPoint *) points->At(jj))->GetPrecisionZ();
        totdata[nsm*4+nfm]=totdata[nsm*4+nfm]+1;
      } 
      delete so;
    }

  
  for(Int_t i=0; i<72 ;i++){
    Float_t numx=0, numy=0,numz=0, comodox=0, comodoy=0, comodoz=0,denx=0, deny=0, denz=0;
    if(totdata[i]!=0){    
      for(Int_t j=0; j<totdata[i]; j++){
        comodox=1/(data[i][1][j]/10*data[i][1][j]/10);//precision in mm, position in cm
        numx=numx+data[i][0][j]*comodox;
        denx=denx+comodox;
        comodoy=1/(data[i][3][j]/10*data[i][3][j]/10);
        numy=numy+data[i][2][j]*comodoy;
        deny=deny+comodoy;
        comodoz=1/(data[i][5][j]/10*data[i][5][j]/10);
        numz=numz+data[i][4][j]*comodoz;
        denz=denz+comodoz;
        }
      fCombFMData[i][1]=TMath::Sqrt(1/denx); //error for x position
      fCombFMData[i][3]=TMath::Sqrt(1/deny); //error for y position
      fCombFMData[i][5]=TMath::Sqrt(1/denz); //error for z position
      fCombFMData[i][0]=numx/denx;           //combined survey data for x position of FM
      fCombFMData[i][2]=numy/deny;           //combined survey data for y position of FM
      fCombFMData[i][4]=numz/denz;           //combined survey data for z position of FM
      } else continue;
    }

  for(Int_t i=0;i<72;i++)
    if (fCombFMData[i][0]!=0){
      fNFMforSM[(i-i%4)/4][i%4]=1;
      fNFMforSM[(i-i%4)/4][4]=fNFMforSM[(i-i%4)/4][4]+1;
    }
}

//_____________________________________________________________________________
void AliTOFAlignment::ReadSurveyDataAndAlign(){
  //
  // read the survey data and, if we know the positions of at least 3 FM 
  //for a SM, call the right Alignement procedure  

  fTOFAlignObjArray = new TObjArray(kMaxAlignObj);

  Float_t deltaFM0=0, deltaFM1=0, deltaFM2=0, deltaFM3=0;

  for(Int_t i=0; i<18; i++){
    switch(fNFMforSM[i][4]){
    case 0:
      printf("we don't know the position of any FM of SM %i\n",i);
      break;
    case 1:
      printf("we know the position of only one FM for SM %i\n",i);
     
      break;
    case 2:
      printf("we know the position of only 2 FM for SM %i\n",i);
      
      break;
    case 3:
      if (fNFMforSM[i][0]==1 && fNFMforSM[i][1]==1 && fNFMforSM[i][2]==1){
        printf("we know the position of FM A B C for SM %i\n",i);
        AliTOFAlignment::AlignFromSurveyABC(i);};

        
      if (fNFMforSM[i][0]==1 && fNFMforSM[i][1]==1 && fNFMforSM[i][3]==1){
        printf("we know the position of FM A B D for SM %i\n",i);
        AliTOFAlignment::AlignFromSurveyABD(i);};

        
      if (fNFMforSM[i][0]==1 && fNFMforSM[i][2]==1 && fNFMforSM[i][3]==1){
        printf("we know the position of FM A C D for SM %i\n",i);
        AliTOFAlignment::AlignFromSurveyACD(i);};

        
      if (fNFMforSM[i][1]==1 && fNFMforSM[i][2]==1 && fNFMforSM[i][3]==1){
        printf("we know the position of FM B C D for SM %i\n",i);
        AliTOFAlignment::AlignFromSurveyBCD(i);};

        
      break;
    case 4:
      printf("we know the position of all the 4 FM for SM %i\n",i);
      //check the precision of the measurement

      deltaFM0=fCombFMData[i*4][1]/TMath::Abs(fCombFMData[i*4][0])+fCombFMData[i*4][3]/TMath::Abs(fCombFMData[i*4][2])+fCombFMData[i*4][5]/TMath::Abs(fCombFMData[i*4][4]);
      deltaFM1=fCombFMData[i*4+1][1]/TMath::Abs(fCombFMData[i*4+1][0])+fCombFMData[i*4+1][3]/TMath::Abs(fCombFMData[i*4+1][2])+fCombFMData[i*4+1][5]/TMath::Abs(fCombFMData[i*4+1][4]);
      deltaFM2=fCombFMData[i*4+2][1]/TMath::Abs(fCombFMData[i*4+2][0])+fCombFMData[i*4+2][3]/TMath::Abs(fCombFMData[i*4+2][2])+fCombFMData[i*4+2][5]/TMath::Abs(fCombFMData[i*4+2][4]);
      deltaFM3=fCombFMData[i*4+3][1]/TMath::Abs(fCombFMData[i*4+3][0])+fCombFMData[i*4+3][3]/TMath::Abs(fCombFMData[i*4+3][2])+fCombFMData[i*4+3][5]/TMath::Abs(fCombFMData[i*4+3][4]);

      //to AlignFromSurvey we use the 3 FM whose positions are known with greatest precision
      if(deltaFM0>=deltaFM1 && deltaFM0>=deltaFM2 && deltaFM0>=deltaFM3){
	printf("to Align we use FM B,C,D");
	AliTOFAlignment::AlignFromSurveyBCD(i);} else
	  if(deltaFM1>=deltaFM0 && deltaFM1>=deltaFM2 && deltaFM1>=deltaFM3){
           printf("to Align we use FM A,C,D");
           AliTOFAlignment::AlignFromSurveyACD(i);} else
	     if(deltaFM2>=deltaFM0 && deltaFM2>=deltaFM1 && deltaFM2>=deltaFM3){
               printf("to Align we use FM A,B,D");
	       AliTOFAlignment::AlignFromSurveyABD(i);} else{
                 printf("to Align we use FM A,B,C");
	         AliTOFAlignment::AlignFromSurveyABC(i);}
     
      break;
    }
  
  }

    // saving TOF AligObjs from survey on a file, for the moment.. 
  fNTOFAlignObj=fTOFAlignObjArray->GetEntries();
  AliInfo(Form("Number of Alignable Volumes: %d",fNTOFAlignObj));
  TFile f("TOFAlignFromSurvey.root","RECREATE");
  f.cd();
  f.WriteObject(fTOFAlignObjArray,"TOFAlignObjs","kSingleKey");
  f.Close();
  

}

//_____________________________________________________________________________
void AliTOFAlignment::AlignFromSurveyABC(Int_t iSM)
{

  //From Survey data, derive the needed transformations to get the 
  //Alignment Objects. 
  //Again, highly "inspired" to Raffaele's example... 
  //we use FM A,B,C
    
    Double_t ngA[3], ngB[3], ngC[3]; // real FM point coord., global RS
    // Get the 'realistic' input from the Survey Matrix
      for(Int_t coord=0;coord<3;coord++){
      ngA[coord]=   fCombFMData[iSM*4][coord*2];
      ngB[coord]=   fCombFMData[iSM*4+1][coord*2];
      ngC[coord]=   fCombFMData[iSM*4+2][coord*2];
      }

    printf("\n\n******Survey analysis for TOF SuperModule ************** %i \n",iSM);

    // From the real fiducial marks coordinates derive back the
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
      orig[i] = md[i] - plane[i]*fgkZFM;
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
    Double_t rot[9] = {ab[0],bc[0],plane[0],ab[1],bc[1],plane[1],ab[2],bc[2],plane[2]}; // the rotation matrix
    // the Aligned matrix for the current TOF SM in the Global RS, as derived from Survey
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
    gdelta.Print(); //this is the global delta trasformation
    
    // Now Write the Alignment Objects....
    Int_t index=0; //let all SM modules have index=0
    AliGeomManager::ELayerID layer = AliGeomManager::kInvalidLayer;
    UShort_t dvoluid = AliGeomManager::LayerToVolUID(layer,index); //dummy vol id 
    TString symname(Form("TOF/sm%02d",iSM));
    AliAlignObjMatrix* o = new AliAlignObjMatrix(symname.Data(),dvoluid,gdelta,kTRUE);
    fTOFAlignObjArray->Add(o);

  }


//_____________________________________________________________________________
void AliTOFAlignment::AlignFromSurveyABD(Int_t iSM)
{
  
  //From Survey data, derive the needed transformations to get the 
  //Alignment Objects. 
  //Again, highly "inspired" to Raffaele's example... 
  //we use FM A,B,D
    
  Double_t ngA[3], ngB[3], ngD[3];// real FM point coord., global RS
    
   // Get the 'realistic' input from the Survey Matrix
      for(Int_t coord=0;coord<3;coord++){
      ngA[coord]=   fCombFMData[iSM*4][coord*2];
      ngB[coord]=   fCombFMData[iSM*4+1][coord*2];
      ngD[coord]=   fCombFMData[iSM*4+3][coord*2];
      }

    printf("\n\n******Survey analysis for TOF SuperModule ************** %i \n",iSM);

    // From the new fiducial marks coordinates derive back the
    // new global position of the surveyed volume
    //*** What follows is the actual survey-to-alignment procedure
    
    Double_t ab[3], ad[3], n[3];
    Double_t plane[4], s=1.;
    
    // first vector on the plane of the fiducial marks
    for(Int_t i=0;i<3;i++){
      ab[i] = (ngB[i] - ngA[i]);
    }
    
    // second vector on the plane of the fiducial marks
    for(Int_t i=0;i<3;i++){
      ad[i] = (ngD[i] - ngA[i]);
    }
    
    // vector normal to the plane of the fiducial marks obtained
    // as cross product of the two vectors on the plane d0^d1
    n[0] = (ab[1] * ad[2] - ab[2] * ad[1]);
    n[1] = (ab[2] * ad[0] - ab[0] * ad[2]);
    n[2] = (ab[0] * ad[1] - ab[1] * ad[0]);
    
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
      md[i] = (ngB[i] + ngD[i]) * 0.5;
    }
    
    // The center of the box, gives the global translation
    for(Int_t i=0;i<3;i++){
      orig[i] = md[i] - plane[i]*fgkZFM;
    }
    
    // get local directions needed to write the global rotation matrix
    // for the surveyed volume by normalising vectors ab and bc
    Double_t sx = TMath::Sqrt(ab[0]*ab[0] + ab[1]*ab[1] + ab[2]*ab[2]);
    if(sx>1.e-8){
      for(Int_t i=0;i<3;i++){
	ab[i] /= sx;
      }
    }
    Double_t sy = TMath::Sqrt(ad[0]*ad[0] + ad[1]*ad[1] + ad[2]*ad[2]);
    if(sy>1.e-8){
      for(Int_t i=0;i<3;i++){
	ad[i] /= sy;
      }
    }
    Double_t rot[9] = {ab[0],ad[0],plane[0],ab[1],ad[1],plane[1],ab[2],ad[2],plane[2]};
    // the Aligned matrix for the current TOF SM in the Global RS, as derived from Survey:
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
    gdelta.Print();  //global delta trasformation
    
    // Now Write the Alignment Objects....
    Int_t index=0; //let all SM modules have index=0
    AliGeomManager::ELayerID layer = AliGeomManager::kInvalidLayer;
    UShort_t dvoluid = AliGeomManager::LayerToVolUID(layer,index); //dummy vol id 
    TString symname(Form("TOF/sm%02d",iSM));
    AliAlignObjMatrix* o = new AliAlignObjMatrix(symname.Data(),dvoluid,gdelta,kTRUE);
    fTOFAlignObjArray->Add(o);

  }
//_____________________________________________________________________________
void AliTOFAlignment::AlignFromSurveyACD(Int_t iSM)
{
  //From Survey data, derive the needed transformations to get the 
  //Alignment Objects. 
  //Again, highly "inspired" to Raffaele's example... 
  //we use FM A,C,D
  
    
    Double_t ngA[3], ngC[3], ngD[3];// real FM point coord., global RS
    
   // Get the 'realistic' input from the Survey Matrix
      for(Int_t coord=0;coord<3;coord++){
      ngA[coord]=   fCombFMData[iSM*4][coord*2];
      ngC[coord]=   fCombFMData[iSM*4+2][coord*2];
      ngD[coord]=   fCombFMData[iSM*4+3][coord*2];
      }

    printf("\n\n******Survey analysis for TOF SuperModule ************** %i \n",iSM);

    // From the new fiducial marks coordinates derive back the
    // new global position of the surveyed volume
    //*** What follows is the actual survey-to-alignment procedure
    
    Double_t cd[3], ad[3], n[3];
    Double_t plane[4], s=1.;
    
    // first vector on the plane of the fiducial marks
    for(Int_t i=0;i<3;i++){
      cd[i] = (ngC[i] - ngD[i]);
    }
    
    // second vector on the plane of the fiducial marks
    for(Int_t i=0;i<3;i++){
      ad[i] = (ngD[i] - ngA[i]);
    }
    
    // vector normal to the plane of the fiducial marks obtained
    // as cross product of the two vectors on the plane d0^d1
    n[0] = (ad[1] * cd[2] - ad[2] * cd[1]);
    n[1] = (ad[2] * cd[0] - ad[0] * cd[2]);
    n[2] = (ad[0] * cd[1] - ad[1] * cd[0]);
    
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
      orig[i] = md[i] + plane[i]*fgkZFM;
    }
    
    // get local directions needed to write the global rotation matrix
    // for the surveyed volume by normalising vectors ab and bc
    Double_t sx = TMath::Sqrt(ad[0]*ad[0] + ad[1]*ad[1] + ad[2]*ad[2]);
    if(sx>1.e-8){
      for(Int_t i=0;i<3;i++){
	ad[i] /= sx;
      }
    }
    Double_t sy = TMath::Sqrt(cd[0]*cd[0] + cd[1]*cd[1] + cd[2]*cd[2]);
    if(sy>1.e-8){
      for(Int_t i=0;i<3;i++){
	cd[i] /= sy;
      }
    }
    Double_t rot[9] = {cd[0],ad[0],-plane[0],cd[1],ad[1],-plane[1],cd[2],ad[2],-plane[2]};
    // the Aligned matrix for the current TOF SM in the Global RS, as derived from Survey:
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
    gdelta.Print(); //global delta trasformation
    
    // Now Write the Alignment Objects....
    Int_t index=0; //let all SM modules have index=0
    AliGeomManager::ELayerID layer = AliGeomManager::kInvalidLayer;
    UShort_t dvoluid = AliGeomManager::LayerToVolUID(layer,index); //dummy vol id 
    TString symname(Form("TOF/sm%02d",iSM));
    AliAlignObjMatrix* o = new AliAlignObjMatrix(symname.Data(),dvoluid,gdelta,kTRUE);
    fTOFAlignObjArray->Add(o);
  }

//___________________________________________________________________________
void AliTOFAlignment::AlignFromSurveyBCD(Int_t iSM)
{
  //From Survey data, derive the needed transformations to get the 
  //Alignment Objects. 
  //Again, highly "inspired" to Raffaele's example... 
  //we use FM B,C,D

    Double_t ngB[3], ngC[3], ngD[3];// real FM point coord., global RS
    
    
   // Get the 'realistic' input from the Survey Matrix
      for(Int_t coord=0;coord<3;coord++){
      ngB[coord]=   fCombFMData[iSM*4+1][coord*2];
      ngC[coord]=   fCombFMData[iSM*4+2][coord*2];
      ngD[coord]=   fCombFMData[iSM*4+3][coord*2];
      }

    printf("\n\n******Survey analysis for TOF SuperModule ************** %i \n",iSM);

    // From the new fiducial marks coordinates derive back the
    // new global position of the surveyed volume
    //*** What follows is the actual survey-to-alignment procedure
    
    Double_t cd[3], bc[3], n[3];
    Double_t plane[4], s=1.;
    
    // first vector on the plane of the fiducial marks
    for(Int_t i=0;i<3;i++){
      cd[i] = (ngC[i] - ngD[i]);
    }
    
    // second vector on the plane of the fiducial marks
    for(Int_t i=0;i<3;i++){
      bc[i] = (ngC[i] - ngB[i]);
    }
    
    // vector normal to the plane of the fiducial marks obtained
    // as cross product of the two vectors on the plane d0^d1
    n[0] = (bc[1] * cd[2] - bc[2] * cd[1]);
    n[1] = (bc[2] * cd[0] - bc[0] * cd[2]);
    n[2] = (bc[0] * cd[1] - bc[1] * cd[0]);
    
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
    plane[3] = ( plane[0] * ngB[0] + plane[1] * ngB[1] + plane[2] * ngB[2] );
    
    // The center of the square with fiducial marks as corners
    // as the middle point of one diagonal - md
    // Used below to get the center - orig - of the surveyed box

    Double_t orig[3], md[3];
    for(Int_t i=0;i<3;i++){
      md[i] = (ngB[i] + ngD[i]) * 0.5;
    }
    
    // The center of the box, gives the global translation
    for(Int_t i=0;i<3;i++){
      orig[i] = md[i] + plane[i]*fgkZFM;
    }
    
    // get local directions needed to write the global rotation matrix
    // for the surveyed volume by normalising vectors ab and bc
    Double_t sx = TMath::Sqrt(cd[0]*cd[0] + cd[1]*cd[1] + cd[2]*cd[2]);
    if(sx>1.e-8){
      for(Int_t i=0;i<3;i++){
	cd[i] /= sx;
      }
    }
    Double_t sy = TMath::Sqrt(bc[0]*bc[0] + bc[1]*bc[1] + bc[2]*bc[2]);
    if(sy>1.e-8){
      for(Int_t i=0;i<3;i++){
	bc[i] /= sy;
      }
    }
    Double_t rot[9] = {cd[0],bc[0],-plane[0],cd[1],bc[1],-plane[1],cd[2],bc[2],-plane[2]};
    // the Aligned matrix for the current TOF SM in the Global RS, as derived from Survey:
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
    gdelta.Print(); //global delta trasformation
    
    // Now Write the Alignment Objects....
    Int_t index=0; //let all SM modules have index=0
    AliGeomManager::ELayerID layer = AliGeomManager::kInvalidLayer;
    UShort_t dvoluid = AliGeomManager::LayerToVolUID(layer,index); //dummy vol id 
    TString symname(Form("TOF/sm%02d",iSM));
    AliAlignObjMatrix* o = new AliAlignObjMatrix(symname.Data(),dvoluid,gdelta,kTRUE);
    fTOFAlignObjArray->Add(o);
  }


