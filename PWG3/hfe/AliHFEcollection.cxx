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

#include <iostream>

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
  AliInfo("DESTRUCTOR");
}
//___________________________________________________________________
Bool_t AliHFEcollection::CreateTH1F(const char* _name, const char* _title, Int_t _nBin, Float_t _nMin, Float_t _nMax){
  
  if(!fListE){
    AliError("No TList pointer ! ");
    return kFALSE;
  }
  else{
    fListE->Add(new TH1F(_name, _title, _nBin, _nMin, _nMax));
     return CheckObject(_name);
  }
}
//___________________________________________________________________
Bool_t AliHFEcollection::CreateTH2F(const char* _name, const char* _title, Int_t _nBinX, Float_t _nMinX, Float_t _nMaxX, Int_t _nBinY, Float_t _nMinY, Float_t _nMaxY){
  
  if(!fListE){
    AliError("No TList pointer ! ");
    return kFALSE;
  }
  fListE->Add(new TH2F(_name, _title, _nBinX, _nMinX, _nMaxX, _nBinY, _nMinY, _nMaxY));
  return CheckObject(_name); 
}
//___________________________________________________________________
Bool_t AliHFEcollection::CreateTH1Fvector1(Int_t _X, const char* _name, const char* _title, Int_t _nBin, Float_t _nMin, Float_t _nMax){

  // create a 1 dimensional array of size [_X]
  
  if(!fListE){
    AliError("No TList pointer ! ");
    return kFALSE;
  }
  if(_X <=0){
    AliError("can not create array with negative or zero size ");
    return kFALSE;
  }
  TString name;
  for(Int_t i=0; i<_X; ++i){
    name = "";
    name.Append(Form("%s_[%d]", _name, i));
    //cout<<" -D: name: "<<name.str().c_str()<<endl;
    //cout<<" -D: nBin: "<<_nBin<<" ,Min: "<<_nMin<<" , Max: "<<_nMax<<endl;
    CreateTH1F(name.Data(), _title, _nBin, _nMin, _nMax);
    if(!CheckObject(name.Data())){
      AliError(Form("Not possible to create object: ", name.Data()));
      return kFALSE;
    }    
  }
  return kTRUE;  
}
//___________________________________________________________________
Bool_t AliHFEcollection::CreateTH2Fvector1(Int_t _X, const char* _name, const char* _title, Int_t _nBinX, Float_t _nMinX, Float_t _nMaxX, Int_t _nBinY, Float_t _nMinY, Float_t _nMaxY){

  // create a 1 dimensinal array of TH2F histograms with size [_X]
  if(!fListE){
    AliError("No TList pointer !");
    return kFALSE;
  }
  if(_X <=0){
    AliError("can not create array with negative or zero size ");
    return kFALSE;
  }
  TString name;
  for(Int_t i=0; i<_X; ++i){
    name = "";
    name.Append(Form("%s_[%d]", _name, i));
    //cout<<" -D: name: "<<name<<endl;
    //cout<<" -D: nBin: "<<_nBin<<" ,Min: "<<_nMin<<" , Max: "<<_nMax<<endl;
    CreateTH2F(name.Data(), _title, _nBinX, _nMinX, _nMaxX, _nBinY, _nMinY, _nMaxY);
    if(!CheckObject(name.Data())){
      AliError(Form("Not possible to create object: %s", name.Data()));
      return kFALSE;
    }    
  }
  return kTRUE;  
}
//___________________________________________________________________
Bool_t AliHFEcollection::CreateTH1Fvector2(Int_t _X, Int_t _Y, const char* _name, const char* _title, Int_t _nBin, Float_t _nMin, Float_t _nMax){

  // create a 2 dimensional array of histograms of size [_X, _Y]
  if(!fListE){
    AliError("No TList pointer ! ");
    return kFALSE;
  }
  if(_X <=0 || _Y <=0){
    AliError("can not create array with negative or zero size ");
    return kFALSE;
  }
  TString name;
  for(Int_t i=0; i<_X; ++i){
    for(Int_t j=0; j<_Y; ++j){
      name = "";
      name.Append(Form("%s_[%d][%d]", _name, i, j));
      //cout<<" -D: name: "<<name.str().c_str()<<endl;
      //cout<<" -D: nBin: "<<_nBin<<" ,Min: "<<_nMin<<" , Max: "<<_nMax<<endl;
      CreateTH1F(name.Data(), _title, _nBin, _nMin, _nMax);
      if(!CheckObject(name.Data())){
	      AliError(Form("Not possible to create object: %s", name.Data()));
	      return kFALSE;
      }
    }
  }
  return kTRUE;
  
  
}
//___________________________________________________________________
TObject* AliHFEcollection::Get(const char* _name, Int_t _X){
  
  TString name = _name;
  name.Append(Form("_[%d]", _X));
  if(!CheckObject(name.Data())){
    AliError("No such object found in the list");
    return 0x0;
  }
  else{
    return Get(name.Data());
  }
}
//___________________________________________________________________
TObject* AliHFEcollection::Get(const char* _name, Int_t _X, Int_t _Y){
  
  TString name = _name;
  name.Append(Form("_[%d][%d]", _X, _Y));
  if(!CheckObject(name.Data())){
    AliError("No such object found in the list");
    AliError(Form("name: %s", name.Data()));
    return 0x0;
  }
  else{
    return Get(name.Data());
  }
}
//___________________________________________________________________
Bool_t AliHFEcollection::CheckObject(const char* _name){

  // check wheter the creation of the histogram was succesfull
  
  if(!fListE){
    AliError("No TList pointer ! ");
    return kFALSE;
  }
  
  if(!fListE->FindObject(_name)){
    AliError("Creating or Finding the object failed");
    return kFALSE;
  }
  return kTRUE;
}
//___________________________________________________________________
TObject* AliHFEcollection::Get(const char* _name){ 
  return fListE->FindObject(_name); 
}
//___________________________________________________________________
Long64_t AliHFEcollection::Merge(TCollection *list){

  if(!fListE){
    AliError("AliHFEcollection::Merge : No TList pointer ! ");
    return 0;
  }

  return fListE->Merge(list);
  
}
//____________________________________________________________________
void AliHFEcollection::Browse(TBrowser *b)
{
   // Browse the content of the directory.

   if (b) {
      TObject *obj = 0;
      TIter nextin(fListE);

      //Add objects that are only in memory
      while ((obj = nextin())) {
         b->Add(obj, obj->GetName());
      }
   }
}
