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

#include "AliACORDEAlign.h"
#include "TROOT.h"
#include "Riostream.h"
#include "TFile.h"
#include "TMath.h"
#include "TSystem.h"
#include "AliSurveyObj.h"
#include "AliAlignObjParams.h"
#include "AliCDBStorage.h"
#include <TClonesArray.h>
#include <TFile.h>
#include "AliLog.h"
#include "AliCDBManager.h"
#include "AliSurveyPoint.h" 
#include "AliACORDEConstants.h" 

// Class creating the ACORDE aligmnent objects 
// from the surveys done by surveyers at Point2.
// Survey results are fetched from 
// Survey Depot, based on survey results 
// position of ACORDE alignment objects is computed.


ClassImp(AliACORDEAlign)


//________________________________________________________________________
AliACORDEAlign::AliACORDEAlign() :
  TObject(),
  fFileGlob(0x0),
  fRepLoc(0),
  fRepGlob(0),
  fUser(0x0),  
  fX(),
  fAlignACORDEObjArray(0x0),
  fDebug(0)
{
  //
  //  default constructor
  //
}

AliACORDEAlign::AliACORDEAlign(/*Int_t reportloc,*/Int_t reportglob):
  TObject(),
  fFileGlob(0x0),
  fRepLoc(0),
  fRepGlob(0),
  fUser(0x0),
  fX(120,4),
  fAlignACORDEObjArray(0x0),
  fDebug(0)
{
  
//
  // constructor
  //fRepLoc = new reportloc[80];
  //fRepGlob = new reportglob[80];
  Char_t path[50];
   fFileGlob = new Char_t[80];
   fUser = new Char_t[10];
  snprintf(path,50,"%s",gSystem->Getenv("ALICE_ROOT")); 
  // 
  snprintf(fFileGlob,80,"%s/ACORDE/Survey_%d_ACORDE.txt",path,reportglob);
  //
 snprintf(fUser,10,"%s",gSystem->Getenv("alien_API_USER"));

}




AliACORDEAlign::AliACORDEAlign(const AliACORDEAlign &align):
  TObject(),
  fFileGlob(0x0),
  fRepLoc(0),
  fRepGlob(0),
  fUser(0x0),
  fX(),
  fAlignACORDEObjArray(0x0),
  fDebug(0)
{
  //
  //  default copy constructor
fDebug = align.fDebug;
}

//__________________________________________________________________________
AliACORDEAlign & AliACORDEAlign::operator =(const AliACORDEAlign &align)
 

{
  //
  // assignment operator - dummy

 //
fDebug = align.fDebug;
  return (*this);
}

//__________________________________________________________________________
AliACORDEAlign::~AliACORDEAlign(){
  //
  // destructor
  //
if(fAlignACORDEObjArray) delete fAlignACORDEObjArray;
  if(fFileGlob) delete[] fFileGlob;
  if(fUser) delete[] fUser;
}


void AliACORDEAlign::LoadSurveyData()
{

//
// Create a new survey object and fill it.
 
AliSurveyObj * s1 = new AliSurveyObj(); 

if(fRepLoc != 0) 
 { 
 // Filling from DCDB (via GRID)
 s1->SetGridUser(fUser);
 s1->Fill("ACORDE",1014872,1,fUser); 
 }
 else
 {
   s1->FillFromLocalFile(fFileGlob);
 }


 //s1->GetEntries();
 //s1->GetUnits();
 //TObjArray* arr = s1->GetData();
 //cout<< "number of entries " << arr->GetEntries() <<endl;
 //arr->UncheckedAt(0)->ClassName();
 //AliSurveyPoint *sp0 = (AliSurveyPoint*) arr->UncheckedAt(0);   
 //cout << "point name " << sp0->GetPointName() << endl  ;
  

//
TString ML= "M" ;
//TString PL= "P";
TString underscore =  "_";
TString  endInner =  "_I";
TString  endOuter =  "_O";
TString  endCenter = "_P";
//
TString surveyname;
TString surveynameInner;
TString surveynameOuter;
TString surveynameCenter;

//TString surveynameAngles;
// 
TString pointNamesInner[60];
TString pointNamesOuter[60];
TString pointNamesCenter[60];

//
Int_t  nid=0;
//
 //for regular modules 
 for (Int_t ncolum=0; ncolum<6; ncolum++)
   {
     for (Int_t nrow=0; nrow<10; nrow++)
       {	
         
         surveyname=ML;
         surveyname+=ncolum;
         surveyname+=underscore;
         surveyname+=nrow;
  
         surveynameInner=surveyname; 
         surveynameInner+=endInner;

         surveynameOuter=surveyname;
         surveynameOuter+=endOuter; 
                          
         surveynameCenter=surveyname;
         surveynameCenter+=endCenter;
  
        pointNamesInner[nid] =  surveynameInner;
        pointNamesOuter[nid] = surveynameOuter; 
        pointNamesCenter[nid] = surveynameCenter;
  	 ++nid; 
       }
   }


//Read  two points 
AliSurveyPoint  *InnerPoint;
AliSurveyPoint *OuterPoint
; 
AliSurveyPoint  *CenterPoint;


 for(Int_t i=0;i<60;i++)
 {

   InnerPoint=0;
   OuterPoint=0; 
   CenterPoint=0;
 
   InnerPoint = (AliSurveyPoint *) s1->GetData()->FindObject(pointNamesInner[i]);
   OuterPoint = (AliSurveyPoint *) s1->GetData()->FindObject(pointNamesOuter[i]);
   CenterPoint = (AliSurveyPoint *) s1->GetData()->FindObject(pointNamesCenter[i]);



  if(InnerPoint && OuterPoint)
   {     
     //Use center if it is available
     if(CenterPoint)
       { 
         fX(i+60,0) =  CenterPoint->GetX()*100;  
         fX(i+60,1) =  CenterPoint->GetY()*100; 
         fX(i+60,2) =  CenterPoint->GetZ()*100;
        }
      else
        {
      //calculate center point 
         fX(i+60,0) = 100*(InnerPoint->GetX() + OuterPoint->GetX())/2.0;  
         fX(i+60,1) = 100*(InnerPoint->GetY() + OuterPoint->GetY())/2.0; 
         fX(i+60,2) = 100*(InnerPoint->GetZ() + OuterPoint->GetZ())/2.0;
        } 
   
      fX(i,0) =  OuterPoint->GetX()*100;  
      fX(i,1) =  OuterPoint->GetY()*100;
      fX(i,2) =  OuterPoint->GetZ()*100;
   }
   else 
   {
       if(InnerPoint && CenterPoint) 
         {

          fX(i+60,0) =  CenterPoint->GetX()*100;  
          fX(i+60,1) =  CenterPoint->GetY()*100; 
          fX(i+60,2) =  CenterPoint->GetZ()*100;

          fX(i,0) =  InnerPoint->GetX()*100;  
          fX(i,1) =  InnerPoint->GetY()*100;
          fX(i,2) =  InnerPoint->GetZ()*100;   
         } 
        else
        {
          if(OuterPoint && CenterPoint)
            { 
        
             fX(i+60,0) =  CenterPoint->GetX()*100;  
             fX(i+60,1) =  CenterPoint->GetY()*100; 
             fX(i+60,2) =  CenterPoint->GetZ()*100;

             fX(i,0) =  OuterPoint->GetX()*100;  
             fX(i,1) =  OuterPoint->GetY()*100;
             fX(i,2) =  OuterPoint->GetZ()*100;   
            }
          else
            { 

             fX(i+60,0) = -99.0;  
             fX(i+60,1) = -99.0; 
             fX(i+60,2) = -99.0;

             fX(i,0) =  -99.0;  
             fX(i,1) =  -99.0;
             fX(i,2) =  -99.0;   
            
            }
        }  
   } 


 }//ends  for

 delete s1;

}

void  AliACORDEAlign::ComputePosition()
{


//Residuals for rotations

Double_t theta;
Double_t resphi[60]; 
Double_t resiphi;

for (Int_t imod=0; imod<60; imod++)
  {
   if(TMath::Abs(fX(imod+60,0)-fX(imod,0))>=0.000001)
   {
   theta = (fX(imod+60,1)-fX(imod,1))/(fX(imod+60,0)-fX(imod,0));
   resiphi = TMath::ATan(theta)*(180.0/TMath::Pi());
   // calculate the residuals  special  modules 
   if(imod==0 || imod==9 || imod==50 || imod==59 )
    {    
    resphi[imod] = 0.0-resiphi;
    continue; 
    }
   // for module with no measurements 
   if(imod == 42 )
    {
     resphi[imod]= 0.0;
    continue;
    } 
   //face A
   if(imod>0 && imod <20)
    {
    resphi[imod] = resiphi + 45.0;  
    }
   //face B
   if(imod>=20 && imod <40)
    {
    resphi[imod] = -resiphi;  
    }
   //face C
   if(imod>=40 && imod <60)
    {
    resphi[imod] = resiphi - 45.0;  
    }
 
   }

}


//Get the  residuals for translations 

AliCDBManager* cdb = AliCDBManager::Instance();
if(!cdb->IsDefaultStorageSet()) cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
cdb->SetRun(0);

//AliCDBStorage* storage;
AliGeomManager::LoadGeometry(); 


TString symname;
TString basename = "ACORDE/Array";


// Get the ideal directly from the geometry 
 Double_t *tr;
 TGeoHMatrix *matrix;  
  for (Int_t imod=0; imod<60; imod++)
  {
    symname = basename;
    symname += imod; 
    cout<< symname << endl;
    matrix = AliGeomManager::GetMatrix(symname);
    tr=matrix->GetTranslation();  

    if(imod == 42)
      {
    fX(imod+60,0) = 0.0;  
    fX(imod+60,1) = 0.0;  
    fX(imod+60,2) = 0.0; 
    fX(imod,0) = 0.0;  
    fX(imod,1) = 0.0;  
    fX(imod,2) = 0.0; 
    continue;
      }

    fX(imod+60,0)=fX(imod+60,0)- tr[0];
    fX(imod+60,1)=fX(imod+60,1)- tr[1]- 4.0;
    fX(imod+60,2)=fX(imod+60,2)- tr[2];  
  
    fX(imod,0) = resphi[imod];  
    fX(imod,1) = 0.0;  
    fX(imod,2) = 0.0; 

   }
  



}

//______________________________________________________________________
void AliACORDEAlign::Run(){
  //
  // runs the full chain
  //
  
  //if(!LoadSurveyFromAlienFile("ACORDE",999999,1))
  //{
   // cout<<"Missing points"<<endl;
    //return;
  //}
  //else 
  //{
    //LoadSurveyfromLocalFile("ACORDE",99999,1);
  //} 


  LoadSurveyData();
  ComputePosition();
  //CreateACORDEAlignObjs();
  StoreAlignObj();

}

//_________________________________________________________________________

void AliACORDEAlign::StoreAlignObj()
{
  //
  // Storing ACORDE alignment objects 
  //

  AliCDBManager* cdb = AliCDBManager::Instance();
  if(!cdb->IsDefaultStorageSet()) cdb->SetDefaultStorage("local://$ALICE_ROOT");

  TClonesArray *array = new TClonesArray("AliAlignObjParams",60);
  //
  // storing either in the OCDB or local file
  //

  TString symname;
  TString basename = "ACORDE/Array";
  Int_t iIndex=0; 

  AliGeomManager::ELayerID iLayer = AliGeomManager::kInvalidLayer;
  UShort_t volid = AliGeomManager::LayerToVolUID(iLayer,iIndex);

  Double_t dx=0., dy=0., dz=0., dpsi=0., dtheta=0., dphi=0.0;
  
   for (Int_t imod=0; imod<60; imod++)
     {
       
       dphi = fX(imod,0);
       dtheta = fX(imod,1);
       dpsi = fX(imod,2);
       dx = fX(imod+60,0);
       dy = fX(imod+60,1);
       dz = fX(imod+60,2);  
       symname = basename;
       symname +=  imod;         
       new((*array)[imod]) AliAlignObjParams(symname,volid,dx,dy,dz,dpsi,dtheta,dphi,kFALSE);     
     }


  if( TString(gSystem->Getenv("TOCDB"))!= TString("kTRUE") )
   {
   
    

 // save on file
    const char* filename = "ACORDESurveyMisalignment.root";
    Char_t fullname[80];

    

    snprintf(fullname,80,"%s",filename);
       
   
    TFile *f = new TFile(fullname,"RECREATE");

    
    if(!f)
      {
	AliError("cannot open file for output\n");
	return;
      }
    AliInfo(Form("Saving alignment objects to the file %s", filename));
    f->cd();
    f->WriteObject(array,"ACORDEAlignObjs","kSingleKey");
    f->Close();
  }
  else
    {
      // save in CDB storage
      AliCDBStorage* storage;
      //
       TString Storage = gSystem->Getenv("STORAGE");
       if(!Storage.BeginsWith("local://") && !Storage.BeginsWith("alien://"))
	 {
	   AliError(Form("STORAGE variable set to %s is not valid. Exiting\n",Storage.Data()));
	   return;
	 }
       storage = cdb->GetStorage(Storage.Data());
       if(!storage)
	 {
	   AliError(Form("Unable to open storage %s\n",Storage.Data()));
	   return;
	 }
       //
       AliCDBMetaData* md = new AliCDBMetaData();
       md->SetResponsible("Pedro Podesta");
       md->SetComment("Full misalignment of ACORDE from surveyors");
       AliCDBId id("ACORDE/Align/Data",0,AliCDBRunRange::Infinity());
       storage->Put(fAlignACORDEObjArray,id,md);
    }

}
