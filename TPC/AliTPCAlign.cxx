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
Revision 1.3  2007/10/22 14:20:25  kowal2
Hardwired "ideal" coordinates replaced with those from the DCDB.
Corrected bug in SetDebug method

Revision 1.2  2007/10/03 09:27:50  marian
Extra semicolon removed (Marian)

Revision 1.1  2007/10/01 14:12:45  kowal2
Class creating the aligmnent object fro the surveyor measurements.

*/ 

//
//  Creates the TPC align object
//

#include "AliTPCAlign.h"
#include "AliSurveyPoint.h"
#include "AliSurveyPoint.h"
//
#include "TROOT.h"
#include "Riostream.h"
#include "TFile.h"
#include "TMath.h"
#include "TSystem.h"
#include "AliSurveyObj.h"
#include "AliSurveyPoint.h"
#include "AliAlignObjParams.h"
#include "AliCDBStorage.h"
#include <TClonesArray.h>
#include <TFile.h>
#include "AliLog.h"
#include "AliCDBManager.h"

ClassImp(AliTPCAlign)

AliTPCAlign::AliTPCAlign() :
  TObject(),
  fFileLoc(0x0),
  fFileGlob(0x0),
  fTPCAlignObj(0x0),
  fX(),
  fA(),
  fY(),
  fDebug(0)
{
  //
  //  default constructor
  //
}   
//________________________________________________________________________
AliTPCAlign::AliTPCAlign(Int_t reportloc, Int_t reportglob) :
  TObject(),
  fFileLoc(0x0),
  fFileGlob(0x0),
  fTPCAlignObj(0x0),
  fX(6,1),
  fA(24,6),
  fY(24,1),
  fDebug(0)
{
  //
  // constructor - defines data files
  //
  fFileLoc = new Char_t[80];
  fFileGlob = new Char_t[80];
  Char_t path[50];
  //sprintf(path,gSystem->Getenv("ALICE_ROOT")); 
  snprintf(path,50,gSystem->Getenv("ALICE_ROOT"));
  //
  //sprintf(fFileLoc,"%s/TPC/Survey_%d_TPC.txt",path,reportloc);
  snprintf(fFileLoc,80,"%s/TPC/Survey_%d_TPC.txt",path,reportloc);
  //sprintf(fFileGlob,"%s/TPC/Survey_%d_TPC.txt",path,reportglob);
  snprintf(fFileGlob,80,"%s/TPC/Survey_%d_TPC.txt",path,reportglob);
  //

}
//_________________________________________________________________________
AliTPCAlign::AliTPCAlign(const AliTPCAlign &align) :
  TObject(),
  fFileLoc(0x0),
  fFileGlob(0x0),
  fTPCAlignObj(0x0),
  fX(),
  fA(),
  fY(),
  fDebug(0)
{
  //
  //  copy constructor - dummy
  //
  fDebug = align.fDebug;
}
//__________________________________________________________________________
AliTPCAlign & AliTPCAlign::operator =(const AliTPCAlign & align)
{
  //
  // assignment operator - dummy
  //
  if(this!=&align){
    fDebug=align.fDebug;
  }
  return (*this);
}

//__________________________________________________________________________
AliTPCAlign::~AliTPCAlign(){
  //
  // destructor
  //
  if(fTPCAlignObj) delete fTPCAlignObj;
}
//__________________________________________________________________________
Bool_t AliTPCAlign::LoadSurveyData(){
  //
  // for a time being it loads from the local file the surveyed point
  // and has the ideal points hardwired. I am waiting until Ricardo
  // completes his job
  // 

 AliSurveyObj * s1 = new AliSurveyObj();
 s1->FillFromLocalFile(fFileGlob);
 //
 Int_t numberPoints = 8;
 //
 TString pointNames[8] = {"T1Final_R04241","T1Final_R05241","T1Final_R06241",
                          "T1Final_R07241","T1Final_R08241","T1Final_R10241",
                          "T1Final_R11241","T1Final_R12241"};
 //
 Float_t surveyedPoints[8][3];
 AliSurveyPoint *currPoint;
 //
 for(Int_t i=0;i<numberPoints;i++){
   currPoint=0;
   currPoint = (AliSurveyPoint *) s1->GetData()->FindObject(pointNames[i]);
   //
   if(currPoint){
     surveyedPoints[i][0]=currPoint->GetX();
     surveyedPoints[i][1]=currPoint->GetY();
     surveyedPoints[i][2]=currPoint->GetZ();
     if(fDebug)
     Printf(Form("INFO: Point \"%s\" coordinates read.", pointNames[i].Data()));
   }
   else {
     if(fDebug){
    Printf(Form("ERROR: Essential point missing: \"%s\"", pointNames[i].Data()));
    return 1;
     }
   }  
 }
 //
 //  Ideal points
 //
 Float_t idealPoints[8][3];
 //
 AliSurveyObj * s2 = new AliSurveyObj();
 s2->FillFromLocalFile(fFileLoc);
 //
 TString pointNamesL[8] ={"R04","R05","R06","R07","R08","R10","R11","R12"};
 //
  AliSurveyPoint *currPointL;
  //
 for(Int_t i=0;i<numberPoints;i++){
   currPointL=0;
   currPointL = (AliSurveyPoint *) s2->GetData()->FindObject(pointNamesL[i]);
   if(currPointL){
     idealPoints[i][0]=currPointL->GetY();
     idealPoints[i][1]=currPointL->GetZ();
     idealPoints[i][2]=currPointL->GetX();
          if(fDebug)
     Printf(Form("INFO: Point \"%s\" coordinates read.", pointNamesL[i].Data()));
   }
   else{
     if(fDebug){
    Printf(Form("ERROR: Essential point missing: \"%s\"", pointNamesL[i].Data()));
    return 1; }
   }
 }
 //
 // Create and fill matrices a & y
 //
 for(Int_t i = 0;i<numberPoints;i++){
   for(Int_t j=0;j<3;j++){
     fY(i*3+j,0)=surveyedPoints[i][j]-idealPoints[i][j];

   }
 }
 fA.Zero(); 
 //
 //
 // setting matrix a
 //
 for(Int_t i=0;i<numberPoints;i++){
   fA(3*i,0)=    -idealPoints[i][1];
   fA(3*i,1)=    idealPoints[i][2];
   fA(3*i,3)=1.; 
   fA(3*i+1,0)=    idealPoints[i][0];
   fA(3*i+1,2)=    -idealPoints[i][2];
   fA(3*i+1,4)=    1.;
   fA(3*i+2,1)=   -idealPoints[i][0];
   fA(3*i+2,2)=   idealPoints[i][1];
   fA(3*i+2,5)=1.; 
 }
 //
 delete s1;
 //
 return 0;
}
//_________________________________________________________________
Double_t AliTPCAlign::ComputeTransform(){
  //
  // Here simple matrix operations for the linear least square
  // The trigonometric function sin is approximated by the angle
  // and the cos by "1", because angles are very small (Y-convention!)
  // This secures the linearity of the problem
  // This method returns a sum of squares of residuals
  //
  TMatrixD tt1(TMatrixD::kInverted,(TMatrixD(fA,TMatrixD::kTransposeMult,fA)));
  fX=TMatrixD(tt1,TMatrixD::kMult,TMatrixD(fA,TMatrixD::kTransposeMult,fY));
  //
 TMatrixD Unit(24,24);
 Unit.UnitMatrix();
 TMatrixD xxt1(TMatrixD::kInverted,TMatrixD(fA,TMatrixD::kTransposeMult,fA));
 TMatrixD t(fA,TMatrixD::kMult,TMatrixD(xxt1,TMatrixD::kMultTranspose,fA));
 TMatrixD t2(Unit,TMatrixD::kMinus,t);
 TMatrixD chi2(fY,TMatrixD::kTransposeMult,TMatrixD(t2,TMatrixD::kMult,fY));
 //

 return chi2(0,0);
 	    
}
//_______________________________________________________________________
void AliTPCAlign::CreateAlignObj(){
  //
  // This method creates AliAlignObj and fills it with Euler angles
  // and shifts. The units are degrees and cm.
  // 
 fTPCAlignObj = new AliAlignObjParams();
 fTPCAlignObj->SetSymName("ALIC_1/TPC_M_1");
 fTPCAlignObj->SetVolUID(0);
 Double_t raddeg = TMath::RadToDeg();
 Double_t phi,theta,psi;
 //
 phi=fX(0,0)*raddeg;
 theta=fX(1,0)*raddeg;
 psi=fX(2,0)*raddeg;
 //
 Double_t dx,dy,dz;
 dx=fX(3,0)*100.;
 dy=fX(4,0)*100.;
 dz=fX(5,0)*100.;
 fTPCAlignObj->SetRotation(psi,theta,phi);
 fTPCAlignObj->SetTranslation(dx,dy,dz);
 fTPCAlignObj->Print("");
}
//______________________________________________________________________
void AliTPCAlign::Run(){
  //
  // runs the full chain
  //
  SetDebug(0);
  Bool_t flag = LoadSurveyData();
  if(flag) {
    cout<<"Missing points"<<endl;
    return;
  }
  Double_t chi2 = ComputeTransform();
  if(chi2>0.01) return;
  CreateAlignObj();
  StoreAlignObj();
  //
}
//_________________________________________________________________________
void AliTPCAlign::StoreAlignObj(){
  //
AliCDBManager* cdb = AliCDBManager::Instance();
 if(!cdb->IsDefaultStorageSet()) cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
 //
TClonesArray *array = new TClonesArray("AliAlignObjParams",1);
//
 Double_t shifts[3];
 Double_t rots[3];
 //
 fTPCAlignObj->GetPars(shifts,rots);
 new((*array)[0]) AliAlignObjParams(fTPCAlignObj->GetSymName(),0,shifts[0],
                   shifts[1],shifts[2],rots[0],rots[1],rots[2],kTRUE);

//
// storing either in the OCDB or local file
//
  if( TString(gSystem->Getenv("TOCDB")) != TString("kTRUE") ){
    // save on file
    const char* filename = "TPCSurveyMisalignment.root";
    Char_t fullname[80];
    //sprintf(fullname,"%s/TPC/%s",gSystem->Getenv("ALICE_ROOT"),filename);
    snprintf(fullname,80,"%s/TPC/%s",gSystem->Getenv("ALICE_ROOT"),filename);
    TFile *f = new TFile(fullname,"RECREATE");
    if(!f){
      AliError("cannot open file for output\n");
      return;
    }
    AliInfo(Form("Saving alignment objects to the file %s", filename));
    f->cd();
    f->WriteObject(array,"TPCAlignObjs","kSingleKey");
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
    md->SetResponsible("Marek Kowalski");
    md->SetComment("Full misalignment of entire TPC from surveyors");
    AliCDBId id("TPC/Align/Data",0,AliCDBRunRange::Infinity());
    storage->Put(array,id,md);
  }

}
