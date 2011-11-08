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

#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TProfile.h>
#include <TString.h>
#include <TBrowser.h>
#include <TMath.h>

#include "AliLog.h"

#include "AliHFEcollection.h"

using namespace std;


ClassImp(AliHFEcollection)

//___________________________________________________________________
AliHFEcollection::AliHFEcollection():
  TNamed()
  , fList(NULL)
{

  //
  // default constructor
  //
}
//___________________________________________________________________
AliHFEcollection::AliHFEcollection(const char* name, const char* title):
  TNamed(name, title)
  , fList(NULL)
{
 
  //
  // constructor
  //
 
  fList = new THashList();
  if(fList){
    fList->SetOwner();
    fList->SetName(Form("list_%s", name));
  }
}
//___________________________________________________________________
AliHFEcollection::AliHFEcollection(const AliHFEcollection &c) :
  TNamed(c)
  , fList(NULL)
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

  // Clone List Content
  target.fList = new THashList();          
  target.fList->SetOwner();
  for(Int_t ien = 0; ien < fList->GetEntries(); ien++)
    target.fList->Add(fList->At(ien)->Clone());
}
//___________________________________________________________________
AliHFEcollection::~AliHFEcollection(){

  //
  // Destructor
  //
  delete fList;
  AliDebug(1, "DESTRUCTOR");
}
//___________________________________________________________________
Bool_t AliHFEcollection::CreateTH1F(const char* name, const char* title, Int_t nBin, Float_t nMin, Float_t nMax, Int_t logAxis){

  //
  // Creates a TH1F histogram for the collection
  //

  if(!fList){
    AliError("No TList pointer ! ");
    return kFALSE;
  }
  else{
    fList->Add(new TH1F(name, title, nBin, nMin, nMax));
    if(logAxis >= 0){
      BinLogAxis(name, logAxis);
    }
    return CheckObject(name);
  }
}

//___________________________________________________________________
Bool_t AliHFEcollection::CreateTH1Farray(const char* name, const char* title, Int_t nBin, const Double_t* xbins){

  //
  // Creates a TH1F histogram for the collection 2nd version 
  //

  if(!fList){
    AliError("No TList pointer ! ");
    return kFALSE;
  }
  else{
    fList->Add(new TH1F(name, title, nBin, xbins));
    return CheckObject(name);
  }
}

//___________________________________________________________________
Bool_t AliHFEcollection::CreateTH2F(const char* name, const char* title, Int_t nBinX, Float_t nMinX, Float_t nMaxX, Int_t nBinY, Float_t nMinY, Float_t nMaxY, Int_t logAxis){

  //
  // Creates a TH2F histogram for the collection
  //

  if(!fList){
    AliError("No TList pointer ! ");
    return kFALSE;
  }
  fList->Add(new TH2F(name, title, nBinX, nMinX, nMaxX, nBinY, nMinY, nMaxY));
  if(logAxis >= 0){
    BinLogAxis(name, logAxis);
  }
  return CheckObject(name); 
}
//___________________________________________________________________
Bool_t AliHFEcollection::CreateTH1Fvector1(Int_t X, const char* name, const char* title, Int_t nBin, Float_t nMin, Float_t nMax, Int_t logAxis){

  //
  // create a 1 dimensional array of size [X]
  //

  if(!fList){
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
    CreateTH1F(hname.Data(), title, nBin, nMin, nMax, logAxis);
    if(!CheckObject(hname.Data())){
      AliError(Form("Not possible to create object: %s", hname.Data()));
      return kFALSE;
    }    
  }
  return kTRUE;  
}
//___________________________________________________________________
Bool_t AliHFEcollection::CreateTH2Fvector1(Int_t X, const char* name, const char* title, Int_t nBinX, Float_t nMinX, Float_t nMaxX, Int_t nBinY, Float_t nMinY, Float_t nMaxY, Int_t logAxis){

  //
  // create a 1 dimensinal array of TH2F histograms with size [X]
  //

  if(!fList){
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
    CreateTH2F(hname.Data(), title, nBinX, nMinX, nMaxX, nBinY, nMinY, nMaxY, logAxis);
    if(!CheckObject(hname.Data())){
      AliError(Form("Not possible to create object: %s", hname.Data()));
      return kFALSE;
    }    
  }
  return kTRUE;  
}
//___________________________________________________________________
Bool_t AliHFEcollection::CreateTH1Fvector2(Int_t X, Int_t Y, const char* name, const char* title, Int_t nBin, Float_t nMin, Float_t nMax, Int_t logAxis){

  //
  // create a 2 dimensional array of histograms of size [X, Y]
  //

  if(!fList){
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
      CreateTH1F(hname.Data(), title, nBin, nMin, nMax, logAxis);
      if(!CheckObject(hname.Data())){
	      AliError(Form("Not possible to create object: %s", hname.Data()));
	      return kFALSE;
      }
    }
  }
  return kTRUE;  
}
//___________________________________________________________________
Bool_t AliHFEcollection::CreateTH3F(const char* name, const char* title, Int_t nBinX, Float_t nMinX, Float_t nMaxX, Int_t nBinY, Float_t nMinY, Float_t nMaxY, Int_t nBinZ, Float_t nMinZ, Float_t nMaxZ, Int_t logAxis){

  //
  // Creates a TH2F histogram for the collection
  //

  if(!fList){
    AliError("No TList pointer ! ");
    return kFALSE;
  }
  fList->Add(new TH3F(name, title, nBinX, nMinX, nMaxX, nBinY, nMinY, nMaxY, nBinZ, nMinZ, nMaxZ));
  if(logAxis >= 0){
    BinLogAxis(name, logAxis);
  }
  return CheckObject(name); 
}
//___________________________________________________________________
Bool_t AliHFEcollection::CreateProfile(const char* name, const char* title, Int_t nbins, Double_t xmin, Double_t xmax){
  
  //
  // create a simple TProfile
  //

  if(!fList){
    AliError("No TList pointer ! ");
    return kFALSE;
  }
  fList->Add(new TProfile(name, title, nbins, xmin, xmax));
  return CheckObject(name);

}
//___________________________________________________________________
Bool_t AliHFEcollection::CreateTHnSparse(const char* name, const char* title, Int_t dim, const Int_t* nbins, const Double_t* xmin, const Double_t* xmax){

  //
  // create 'dim' dimensional THnSparse
  //

  if(!fList){
    AliError("No TList pointer ! ");
    return kFALSE;
  }
  fList->Add(new THnSparseF(name, title, dim, nbins, xmin, xmax));
  return CheckObject(name);

}
//___________________________________________________________________
TObject* AliHFEcollection::Get(const char* name){ 

  //
  // Get histogram with the required name
  // 
  

  if(!CheckObject(name)){
    AliWarning(Form("Not possible to return pointer to the object '%s'\n", name));
    return 0;
  }

  return fList->FindObject(name);
  
}
//___________________________________________________________________
Bool_t AliHFEcollection::Fill(const char* name, Double_t v){

  //
  // fill function for one TH1 histograms
  //

  if(!CheckObject(name)){
    AliError(Form("Not possible to fill the object '%s', the object does not exist\n", name));
    return kFALSE;
  }

  TH1 *htmp = dynamic_cast<TH1F*>(fList->FindObject(name));
  // chack the possible object types
  if(htmp){
    htmp->Fill(v);
    return kTRUE;
  }
  
  return kFALSE;

}
//___________________________________________________________________
Bool_t AliHFEcollection::Fill(const char* name, Int_t v){

  //
  // fill function for one TH1 histograms for integer numbers
  //

  return Fill(name, v*1.0);
}
//___________________________________________________________________
Bool_t AliHFEcollection::Fill(const char* name, Int_t X, Double_t v){

  //
  // fill function for one dimension arrays of TH1
  //

  const char* n = Form("%s_[%d]", name, X);
  TObject *o = Get(n);
  if(!o){
    return kFALSE;
  }
  Fill(o->GetName(), v);
  return kTRUE;
}
//___________________________________________________________________
Bool_t AliHFEcollection::Fill(const char* name, Int_t X, Int_t Y, Double_t v){

  //
  // Fill function fir 2 dimensional TH1 arrays
  //
  
  const char* n = Form("%s_[%d][%d]", name, X, Y);
  TObject *o = Get(n);
  if(!o){
    return kFALSE;
  }
  Fill(o->GetName(), v);
  return kTRUE;
}
//___________________________________________________________________
Bool_t AliHFEcollection::Fill(const char* name, Int_t X, Double_t v1, Double_t v2){

  //
  // fill function for one dimension array of TH2
  //

  const char* n = Form("%s_[%d]", name, X);
  TObject *o = Get(n);
  if(!o){
    return kFALSE;
  }
  Fill(o->GetName(), v1, v2);
  
  return kTRUE;
}
//___________________________________________________________________
Bool_t AliHFEcollection::Fill(const char* name, Double_t v1, Double_t v2){

  //
  // fill function for TH2 objects
  //

  if(!CheckObject(name)){
    AliError(Form("Not possible to fill the object '%s', the object does not exist\n", name));
    return kFALSE;
  }

  // chack the possible object types
  if(fList->FindObject(name)->InheritsFrom("TH2")){
    TH2 *h2 = dynamic_cast<TH2F*>(fList->FindObject(name));
    if(h2) h2->Fill(v1, v2);
    return kTRUE;
  }  
  if(fList->FindObject(name)->InheritsFrom("TProfile")){
    TProfile *pr = dynamic_cast<TProfile*>(fList->FindObject(name));
    if(pr) pr->Fill(v1, v2);
    return kTRUE;
  }  
  
  return kFALSE;
  
}
//___________________________________________________________________
Bool_t AliHFEcollection::Fill(const char* name, Double_t v1, Double_t v2, Double_t v3){

  //
  // fill function for TH3 objects
  //

  if(!CheckObject(name)){
    AliError(Form("Not possible to fill the object '%s', the object does not exist\n", name));
    return kFALSE;
  }

  // chack the possible object types
  TH3 *h3 = dynamic_cast<TH3F*>(fList->FindObject(name));
  if(h3){ 
    h3->Fill(v1, v2, v3);
    return kTRUE;
  }

  return kFALSE;
  
}
//___________________________________________________________________
Bool_t AliHFEcollection::Fill(const char* name, Double_t* entry, Double_t weight){
  //
  // Fill a THnSparse object
  //

  if(!CheckObject(name)){
    AliError(Form("Not possible to fill the object '%s', the object does not exist\n", name));
    return kFALSE;
  }
  
  THnSparseF *htmp = dynamic_cast<THnSparseF*>(fList->FindObject(name));
   if(htmp){
     htmp->Fill(entry, weight);
     return kTRUE;
   }
   return kFALSE;

}
//___________________________________________________________________
Bool_t AliHFEcollection::CheckObject(const char* name){

  //
  // check wheter the creation of the histogram was succesfull
  //
  
  if(!fList){
    AliError("No TList pointer ! ");
    return kFALSE;
  }
  
  if(!fList->FindObject(name)){
    AliWarning(Form("Creating or Finding the object '%s' failed\n", name));
    return kFALSE;
  }
  return kTRUE;
}
//___________________________________________________________________
Bool_t AliHFEcollection::Sumw2(const char* name){
  //
  // Set Sumw2 for the given object
  //
  if(!CheckObject(name)){
    return kFALSE;
  }

  TObject *o = Get(name);
  THnSparse *htmp = dynamic_cast<THnSparse*>(o);
  if(htmp){
    htmp->Sumw2();
  }
  return kTRUE;
}
//___________________________________________________________________
Bool_t AliHFEcollection::BinLogAxis(const char* name, Int_t dim){

  // 
  // converts the axis (defined by the dimension) of THx or THnSparse
  // object to Log scale. Number of bins and bin min and bin max are preserved
  //


  if(!CheckObject(name)){
    return kFALSE;
  }

  TObject *o = Get(name);
  TAxis *axis = NULL;
  TString type(o->IsA()->GetName()); 
  if(type.Contains("TH1")){ // 1D histogram
    TH1 *h1 = dynamic_cast<TH1F*>(o);
    if(h1) axis = h1->GetXaxis();
  } else if(type.Contains("TH2")){
    TH2 *h2 = dynamic_cast<TH2F*>(o);
    if(h2){
      if(0 == dim){
        axis = h2->GetXaxis();
      }
      else if(1 == dim){
        axis = h2->GetYaxis();
      }
      else{
         AliError("Only dim = 0 or 1 possible for TH2F");
      }
    }
  } else if(type.Contains("TH3")){
    TH3 *h3 = dynamic_cast<TH3F*>(o);
    if(h3){
      if(0 == dim){
        axis = h3->GetXaxis();
      }
      else if(1 == dim){
        axis = h3->GetYaxis();
      }
      else if(2 == dim){
        axis = h3->GetZaxis();
      }
      else{
         AliError("Only dim = 0, 1 or 2 possible for TH3F");
      }
    }
  } else if(type.Contains("THnSparse")){
    THnSparse *hs = dynamic_cast<THnSparse*>(o);
    if(hs) axis = hs->GetAxis(dim);
  }

  if(!axis){
    AliError(Form("Axis '%d' could not be identified in the object '%s'\n", dim, name));
    return kFALSE;
  }
  
  Int_t bins = axis->GetNbins();

  Double_t from = axis->GetXmin();
  if(from <= 0){
    AliError(Form(" Log binning not possible for object '%s'because the '%d' axis starts from '%f\n'", name, dim, from));
    return kFALSE;
  }
  Double_t to = axis->GetXmax();
  Double_t *newBins = new Double_t[bins+1];
  newBins[0] = from;
  Double_t factor = TMath::Power(to/from, 1./bins);
  for(Int_t i=1; i<=bins; ++i){
    newBins[i] = factor * newBins[i-1];
  }
  axis->Set(bins, newBins);
  delete[] newBins;

  return kTRUE;


}
//___________________________________________________________________
Long64_t AliHFEcollection::Merge(const TCollection *list){

  //
  // Merge the collections
  //
  if(!list)
    return 0;
  if(list->IsEmpty())
    return 1;
  
  TIter it(list);
  TObject *o = NULL;
  Int_t index = 0;
  TList templist;       // Create temporary list containing all the lists to merge
  while((o = it())){
    AliHFEcollection *coll = dynamic_cast<AliHFEcollection *>(o);
    if(!coll) continue; 
    templist.Add(coll->fList);
    index++;
  }
  fList->Merge(&templist);
  return index + 1;
}
//____________________________________________________________________
void AliHFEcollection::Browse(TBrowser *b)
{

  //
  // Browse the content of the directory.
  //

   if (b) {
      TObject *obj = 0;
      TIter nextin(fList);

      //Add objects that are only in memory
      while ((obj = nextin())) {
         b->Add(obj, obj->GetName());
      }
   }
}
//____________________________________________________________________
void AliHFEcollection::Print(Option_t *) const{
  //
  // Print content of the collection
  //
  TIter histIter(fList);
  TObject *o = NULL;
  Int_t nHistos = 0;
  printf("Collection %s\n", GetName());
  printf("Content of the collection:\n=========================================\n");
  while((o = histIter())){
    printf("Histo %s, Type %s\n", o->GetName(), o->IsA()->GetName());
    nHistos++;
  }
  printf("Number of histos in the collection: %d\n", nHistos);
  printf("\n");
}

