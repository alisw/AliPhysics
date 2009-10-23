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
//
// Collection class for histograms
// Stores either histograms or vectors of histograms
// 
// Author:
// Matus Kalisky  <matus.kalisky@cern.ch>
//

//#include <iostream>

#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TString.h>
#include <TBrowser.h>
#include <TIterator.h>

#include "AliLog.h"
#include "AliHFEcollection.h"

using namespace std;


ClassImp(AliHFEcollection)

//___________________________________________________________________
AliHFEcollection::AliHFEcollection():
  TNamed()
  , fListE(0x0)
{
  //
  // default constructor
  //

  fListE = new TList();
  if(!fListE){
    AliError("Initialization of the list failed");
  }
  else{
    // list is owner of the objects. Once list is deleted, the objects
    // it contains will be deleted too
    fListE->SetOwner(kTRUE);
  }
  //Printf("%s:%d,%p",(char*)__FILE__,__LINE__,fInstance);
  
}
//___________________________________________________________________
AliHFEcollection::AliHFEcollection(char* name, char* title):
  TNamed(name, title)
  , fListE(0x0)
{
 
  //
  // constructor
  //
 
  fListE = new TList();
  if(!fListE){
    AliError("Initialization of the list failed");
  }
  else{
    // list is owner of the objects. Once list is deleted, the objects
    // it contains will be deleted too
    fListE->SetOwner(kTRUE);
  }
}
//___________________________________________________________________
AliHFEcollection::AliHFEcollection(const AliHFEcollection &c) :
  TNamed(c)
  , fListE(0x0)
{

  //
  // copy operator
  // 
  
  c.Copy(*this);
}
//___________________________________________________________________
AliHFEcollection &AliHFEcollection::operator=(const AliHFEcollection &ref)
{
  //
  // Assignment operator
  //
  
  if(this != &ref){
    ref.Copy(*this);
  }
  return *this;
}
//___________________________________________________________________
void AliHFEcollection::Copy(TObject &ref) const {
  //
  // Performs the copying of the object
  //
  AliHFEcollection &target = dynamic_cast<AliHFEcollection &>(ref);

  target.fListE = fListE;          
}
//___________________________________________________________________
AliHFEcollection::~AliHFEcollection(){
  //
  // Destructor
  //
  AliInfo("DESTRUCTOR");
}
//___________________________________________________________________
Bool_t AliHFEcollection::CreateTH1F(const char* name, const char* title, Int_t nBin, Float_t nMin, Float_t nMax){
  //
  // Creates a TH1F histogram for the collection
  //
  if(!fListE){
    AliError("No TList pointer ! ");
    return kFALSE;
  }
  else{
    fListE->Add(new TH1F(name, title, nBin, nMin, nMax));
    return CheckObject(name);
  }
}
//___________________________________________________________________
Bool_t AliHFEcollection::CreateTH2F(const char* name, const char* title, Int_t nBinX, Float_t nMinX, Float_t nMaxX, Int_t nBinY, Float_t nMinY, Float_t nMaxY){
  //
  // Creates a TH2F histogram for the collection
  //
  if(!fListE){
    AliError("No TList pointer ! ");
    return kFALSE;
  }
  fListE->Add(new TH2F(name, title, nBinX, nMinX, nMaxX, nBinY, nMinY, nMaxY));
  return CheckObject(name); 
}
//___________________________________________________________________
Bool_t AliHFEcollection::CreateTH1Fvector1(Int_t X, const char* name, const char* title, Int_t nBin, Float_t nMin, Float_t nMax){
  //
  // create a 1 dimensional array of size [X]
  //
  if(!fListE){
    AliError("No TList pointer ! ");
    return kFALSE;
  }
  if(X <=0){
    AliError("can not create array with negative or zero size ");
    return kFALSE;
  }
  TString hname;
  for(Int_t i=0; i<X; ++i){
    hname = "";
    hname.Append(Form("%s_[%d]", name, i));
    //cout<<" -D: name: "<<name.str().c_str()<<endl;
    //cout<<" -D: nBin: "<<_nBin<<" ,Min: "<<_nMin<<" , Max: "<<_nMax<<endl;
    CreateTH1F(hname.Data(), title, nBin, nMin, nMax);
    if(!CheckObject(hname.Data())){
      AliError(Form("Not possible to create object: ", hname.Data()));
      return kFALSE;
    }    
  }
  return kTRUE;  
}
//___________________________________________________________________
Bool_t AliHFEcollection::CreateTH2Fvector1(Int_t X, const char* name, const char* title, Int_t nBinX, Float_t nMinX, Float_t nMaxX, Int_t nBinY, Float_t nMinY, Float_t nMaxY){
  //
  // create a 1 dimensinal array of TH2F histograms with size [X]
  //
  if(!fListE){
    AliError("No TList pointer !");
    return kFALSE;
  }
  if(X <=0){
    AliError("can not create array with negative or zero size ");
    return kFALSE;
  }
  TString hname;
  for(Int_t i=0; i<X; ++i){
    hname = "";
    hname.Append(Form("%s_[%d]", name, i));
    //cout<<" -D: name: "<<name<<endl;
    //cout<<" -D: nBin: "<<_nBin<<" ,Min: "<<_nMin<<" , Max: "<<_nMax<<endl;
    CreateTH2F(hname.Data(), title, nBinX, nMinX, nMaxX, nBinY, nMinY, nMaxY);
    if(!CheckObject(hname.Data())){
      AliError(Form("Not possible to create object: %s", hname.Data()));
      return kFALSE;
    }    
  }
  return kTRUE;  
}
//___________________________________________________________________
Bool_t AliHFEcollection::CreateTH1Fvector2(Int_t X, Int_t Y, const char* name, const char* title, Int_t nBin, Float_t nMin, Float_t nMax){
  //
  // create a 2 dimensional array of histograms of size [X, Y]
  //
  if(!fListE){
    AliError("No TList pointer ! ");
    return kFALSE;
  }
  if(X <=0 || Y <=0){
    AliError("can not create array with negative or zero size ");
    return kFALSE;
  }
  TString hname;
  for(Int_t i=0; i<X; ++i){
    for(Int_t j=0; j<Y; ++j){
      hname = "";
      hname.Append(Form("%s_[%d][%d]", name, i, j));
      //cout<<" -D: name: "<<name.str().c_str()<<endl;
      //cout<<" -D: nBin: "<<_nBin<<" ,Min: "<<_nMin<<" , Max: "<<_nMax<<endl;
      CreateTH1F(hname.Data(), title, nBin, nMin, nMax);
      if(!CheckObject(hname.Data())){
	      AliError(Form("Not possible to create object: %s", hname.Data()));
	      return kFALSE;
      }
    }
  }
  return kTRUE;
  
  
}
//___________________________________________________________________
TObject* AliHFEcollection::Get(const char* name, Int_t X){
  //
  // Get the histogram from the vector
  // Paramter: 
  //  name - vector name
  //  X - Number of the desired histogram
  //
  TString hname = name;
  hname.Append(Form("_[%d]", X));
  if(!CheckObject(hname.Data())){
    AliError("No such object found in the list");
    return 0x0;
  }
  else{
    return Get(hname.Data());
  }
}
//___________________________________________________________________
TObject* AliHFEcollection::Get(const char* name, Int_t X, Int_t Y){
  //
  // Get histogram from the 2D vector
  // Parameters:
  //  name - Name of the vector
  //  X,Y - Indices of the histogram
  //
  TString hname = name;
  hname.Append(Form("_[%d][%d]", X, Y));
  if(!CheckObject(hname.Data())){
    AliError("No such object found in the list");
    AliError(Form("name: %s", hname.Data()));
    return 0x0;
  }
  else{
    return Get(hname.Data());
  }
}
//___________________________________________________________________
Bool_t AliHFEcollection::CheckObject(const char* name){
  //
  // check wheter the creation of the histogram was succesfull
  //
  
  if(!fListE){
    AliError("No TList pointer ! ");
    return kFALSE;
  }
  
  if(!fListE->FindObject(name)){
    AliError("Creating or Finding the object failed");
    return kFALSE;
  }
  return kTRUE;
}
//___________________________________________________________________
TObject* AliHFEcollection::Get(const char* name){ 
  //
  // Get histogram with the required name
  // 
  return fListE->FindObject(name); 
}
//___________________________________________________________________
Long64_t AliHFEcollection::Merge(TCollection *list){
  //
  // Merge the collections
  //
  if(!fListE){
    AliError("AliHFEcollection::Merge : No TList pointer ! ");
    return 0;
  }

  return fListE->Merge(list);
  
}
//____________________________________________________________________
void AliHFEcollection::Browse(TBrowser *b)
{
  //
  // Browse the content of the directory.
  //

   if (b) {
      TObject *obj = 0;
      TIter nextin(fListE);

      //Add objects that are only in memory
      while ((obj = nextin())) {
         b->Add(obj, obj->GetName());
      }
   }
}
