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
$Log: AliFITAlign.cxx,v $
 Revision   2008/01/30
Removing code violations 

 Version 1.1  2006/10
Preliminary test version (T.Malkiewicz)
*/

#include "AliFITAlign.h"
#include "TROOT.h"
#include "Riostream.h"
#include "TFile.h"
#include "TMath.h"
#include "TSystem.h"
#include "TString.h"
#include "AliSurveyObj.h"
#include "AliAlignObjParams.h"
#include "AliCDBStorage.h"
#include <TClonesArray.h>
#include <TFile.h>
#include "AliLog.h"
#include "AliCDBManager.h"
#include "AliSurveyPoint.h" 

// Class creating the FIT aligmnent objects 
// from the surveys done by surveyers at Point2.
// Survey results are fetched from 
// Survey Depot, based on survey results 
// position of FIT alignment objects is computed.

using std::endl;
using std::cout;
ClassImp(AliFITAlign)

AliFITAlign::AliFITAlign() :
  TObject(),
  fFileGlob(0x0),
  fFITAAlignObj(0x0),
  fFITCAlignObj(0x0),
  fDebug(0),
  fXPosC(0.),
  fYPosC(0.),
  fXPosA(0.),
  fYPosA(0.),
  fRepLoc(0),
  fRepGlob(0),
  fSide(0x0),
  fUser(0x0)
{
  //
  //  default constructor
  //
}   
//________________________________________________________________________
AliFITAlign::AliFITAlign(Int_t reportloc, Int_t side, Int_t reportglob) :
  TObject(),
  fFileGlob(0x0),
  fFITAAlignObj(0x0),
  fFITCAlignObj(0x0),
  fDebug(0),
  fXPosC(0.),
  fYPosC(0.),
  fXPosA(0.),
  fYPosA(0.),
  fRepLoc(0),
  fRepGlob(0),
  fSide(0x0),
  fUser(0x0)
{
  //
  // constructor - defines data files
  //
  fRepLoc = reportloc;
  fRepGlob = reportglob;
  fSide = side;
  // Char_t path[50];
  TString path = Form("%s",gSystem->Getenv("ALICE_ROOT")) ;
  // fFileGlob = new Char_t[80];
  //  fUser = new Char_t[10];
  fFileGlob = Form("%s/FIT/Survey_%d_V0.txt",path.Data(),reportglob);
  fUser = Form("%s/FIT/Survey_%d_V0.txt",path.Data(),reportglob);
  // sprintf(path,gSystem->Getenv("ALICE_ROOT")); 
  //
  // sprintf(fFileLoc,"%s/FIT/Survey_%d_FIT.txt",path,reportloc);
  // sprintf(fFileGlob,"%s/FIT/Survey_%d_V0.txt",path,reportglob);
  //
  // sprintf(fUser,gSystem->Getenv("alien_API_USER"));
}
//_________________________________________________________________________
AliFITAlign::AliFITAlign(const AliFITAlign &align) :
  TObject(),
  fFileGlob(0x0),
  fFITAAlignObj(0x0),
  fFITCAlignObj(0x0),
  fDebug(0),
  fXPosC(0.),
  fYPosC(0.),
  fXPosA(0.),
  fYPosA(0.),
  fRepLoc(0),
  fRepGlob(0),
  fSide(0x0),
  fUser(0x0)  
{
  //
  //  copy constructor - dummy
  //
  ((AliFITAlign &) align).Copy(*this);

}
//__________________________________________________________________________
AliFITAlign & AliFITAlign::operator =(const AliFITAlign & align)
{
  //
  // assignment operator - dummy
  //
  if (this != &align) ((AliFITAlign &) align).Copy(*this);

   return (*this);
}

//__________________________________________________________________________
AliFITAlign::~AliFITAlign()
{
  //
  // destructor
  //
  if(fFITAAlignObj) delete fFITAAlignObj;
  if(fFITCAlignObj) delete fFITCAlignObj;
  if(fFileGlob) delete[] fFileGlob;
  if(fUser) delete[] fUser;
}
//__________________________________________________________________________
Bool_t AliFITAlign::LoadSurveyData()
{
  //
  // Create a new survey object and fill it.
 
 AliSurveyObj * s1 = new AliSurveyObj();
 const int numberPoints = 2;
 TString pointNames[numberPoints]={"Flange_0","C67_6_Beamcircle"}; 
 
 if(fRepLoc == 0) 
 { 
   //
   // Filling from DCDB (via GRID)
   //
   s1->SetGridUser(fUser);
   if(fSide == 0) 
   {
     s1->Fill("FIT", fRepGlob, fUser);
   }
   else if(fSide == 1)
   {
     s1->Fill("FIT", fRepGlob, fUser);
   }
   else
   {
     cout<<"Enter the side properly: '0'- A side, '1'- C side'" <<endl;
     return 0;
   }
 }
 else
 //
 // Filling from local file
 //
 {
   s1->FillFromLocalFile(fFileGlob);
 }
 //
 Float_t surveyedPoints [numberPoints][2];
 AliSurveyPoint *currPoint;

 //
 for(Int_t i=0;i<numberPoints;i++)
 {
   currPoint=0;
   currPoint = (AliSurveyPoint *) s1->GetData()->FindObject(pointNames[i]);
   //
   if(currPoint)
   {
     surveyedPoints[i][0]=currPoint->GetX();
     surveyedPoints[i][1]=currPoint->GetY();
   //  surveyedPoints[i]=currPoint->GetZ();
     if(fDebug)
     Printf("INFO: Point %s coordinates read.\n", pointNames[i].Data() ) ;
   }
   else 
   {
     if(fDebug)
     {
       Printf("ERROR: Essential point missing: %s\n", pointNames[i].Data() ) ;
       return 1;
     }
   }  
 }
 if(fSide == 0)
 {
   fXPosA = surveyedPoints[0][0];
   fYPosA = surveyedPoints[0][1];
 }
 else if(fSide == 1)
 {
   fXPosC = surveyedPoints[1][0]; 
   fYPosC = surveyedPoints[1][1];
 }
 //
 delete s1;
 //
 return 0;
}
//_________________________________________________________________

Double_t AliFITAlign::ComputePosition()
{
 //  Float_t fZPos, shift;
 //  fZPos = surveyedPoints[3] - shift;
  return 0;
}
//_______________________________________________________________________
void AliFITAlign::CreateAlignObj()
{
  //
  //  TClonesArray *array = new TClonesArray("AliAlignObjParams",2);
  // TClonesArray &alobj = *array;
  
  Double_t dx=0., dy=0., dz=0., dpsi=0., dtheta=0., dphi=0.;
  dx=fXPosA;
  dy=fYPosA;
  fFITAAlignObj = new AliAlignObjParams("ALIC_1/0STL_1",0,dx,dy,dz,dpsi,dtheta,dphi,kTRUE);
  
  dx=fXPosC;
  dy=fYPosC;
  // dz=surveyedPoints[2];
  dz=0.;
  fFITCAlignObj = new AliAlignObjParams("ALIC_1/0STR_1",0,dx,dy,dz,dpsi,dtheta,dphi,kTRUE);
  
}

//______________________________________________________________________
void AliFITAlign::Run()
{
  //
  // runs the full chain
  //
  // SetDebug(0);
  Bool_t flag = LoadSurveyData();
    if(flag) 
  {
    cout<<"Missing points"<<endl;
    return;
  }
  // ComputePosition();
  CreateAlignObj();
  StoreAlignObj();
}
//_________________________________________________________________________

void AliFITAlign::StoreAlignObj()
{
 //
 // Storing FIT alignment objects 
 //
 AliCDBManager* cdb = AliCDBManager::Instance();
 if(!cdb->IsDefaultStorageSet()) cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
 //
 TClonesArray *array = new TClonesArray("AliAlignObjParams",2);
//
 Double_t shifts[3];
 Double_t rots[3];
 //
 fFITAAlignObj->GetPars(shifts,rots);
 new((*array)[0]) AliAlignObjParams(fFITAAlignObj->GetSymName(),0,shifts[0],
                   shifts[1],shifts[2],rots[0],rots[1],rots[2],kTRUE);
 fFITCAlignObj->GetPars(shifts,rots);
 new((*array)[1]) AliAlignObjParams(fFITCAlignObj->GetSymName(),0,shifts[0],
                   shifts[1],shifts[2],rots[0],rots[1],rots[2],kTRUE);


//
// storing either in the OCDB or local file
//
  if( TString(gSystem->Getenv("TOCDB")) != TString("kTRUE") ){
    // save on file
    const char* filename = "FITSurveyMisalignment.root";
    //  Char_t fullname[80];
    //  sprintf(fullname,"%s/FIT/Align/Data/%s",gSystem->Getenv("ALICE_ROOT"),filename);
    TString fullname = Form("%s/FIT/Align/Data/%s",gSystem->Getenv("ALICE_ROOT"), filename);
    TFile *f = new TFile(fullname.Data(),"RECREATE");
    if(!f){
      AliError("cannot open file for output\n");
      return;
    }
    AliInfo(Form("Saving alignment objects to the file %s", filename));
    f->cd();
    f->WriteObject(array,"FITAlignObjs","kSingleKey");
    f->Close();
  }else{
    // save in CDB storage
    AliCDBStorage* storage;
    //
   TString Storage = gSystem->Getenv("STORAGE");
    if(!Storage.BeginsWith("local://") && !Storage.BeginsWith("alien://")) {
      AliError(Form(
      "STORAGE variable set to %s is not valid. Exiting\n",Storage.Data()));
      return;
    }
    storage = cdb->GetStorage(Storage.Data());
    if(!storage){
      AliError(Form("Unable to open storage %s\n",Storage.Data()));
      return;
    }
    //
    AliCDBMetaData* md = new AliCDBMetaData();
    md->SetResponsible("Tomasz Malkiewicz");
    md->SetComment("Position of FIT-A and FIT-C from survey");
    AliCDBId id("FIT/Align/Data",0,AliCDBRunRange::Infinity());
    storage->Put(array,id,md);
  }
}

