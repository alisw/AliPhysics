/*************************************************************************
* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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
***********************************************************

 Event plane framework with helper functions wrapped in a namespace
 Contact: Jaap Onderwaater, j.onderwaater@gsi.de, jacobus.onderwaater@cern.ch
 
   Based on:
   AliDielectronHelper
   Dielectron helper functions wrapped in a namespace
   Authors: 
   Jens Wiechula <Jens.Wiechula@cern.ch>
   Frederick Kramer <Frederick.Kramer@cern.ch> 
   Julian Book <Julian.Book@cern.ch>
***********************************************************
*/



#ifndef ALIEVENTPLANEHELPER_H
#include "AliEventPlaneHelper.h"
#endif

#include <TError.h>
#include <TMath.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <TVectorD.h>
#include <TF1.h>
#include <TRandom.h>
#include <TProfile.h>
#include <TChain.h>
#include <TKey.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TDirectoryFile.h>
#include <TList.h>
#include <THashList.h>
#include <TClonesArray.h>


#include <iostream>
#include <fstream>

ClassImp(AliEventPlaneHelper)
TDirectoryFile* AliEventPlaneHelper::fgHistCali=0x0;  // main directory for a standard tree analysis output (used for calibration, plotting etc.)
TFile* AliEventPlaneHelper::fgHistCaliFile=0x0;      // pointer to a TFile opened for reading

AliEventPlaneHelper::AliEventPlaneHelper() :
//   TCollection(),
  TNamed("AliEventPlaneHelper","Reduced Histogram Container")
{
  //
  // Default constructor
  //
}

//_____________________________________________________________________________
AliEventPlaneHelper::~AliEventPlaneHelper()
{
  //
  // Destructor
  //
  fgHistCali->Clear("C");
}




//_____________________________________________________________________________
TVectorD* AliEventPlaneHelper::MakeLogBinning(Int_t nbinsX, Double_t xmin, Double_t xmax)
{
  //
  // Make logarithmic binning
  // the user has to delete the array afterwards!!!
  //
  
  //check limits
  if (xmin<1e-20 || xmax<1e-20){
    //Error("AliEventPlaneHelper::MakeLogBinning","For Log binning xmin and xmax must be > 1e-20. Using linear binning instead!");
    return AliEventPlaneHelper::MakeLinBinning(nbinsX, xmin, xmax);
  }
  if (xmax<xmin){
    Double_t tmp=xmin;
    xmin=xmax;
    xmax=tmp;
  }
  TVectorD *binLim=new TVectorD(nbinsX+1);
  Double_t first=xmin;
  Double_t last=xmax;
  Double_t expMax=TMath::Log(last/first);
  for (Int_t i=0; i<nbinsX+1; ++i){
    (*binLim)[i]=first*TMath::Exp(expMax/nbinsX*(Double_t)i);
  }
  return binLim;
}

//_____________________________________________________________________________
TVectorD* AliEventPlaneHelper::MakeLinBinning(Int_t nbinsX, Double_t xmin, Double_t xmax)
{
  //
  // Make linear binning
  // the user has to delete the array afterwards!!!
  //
  if (xmax<xmin){
    Double_t tmp=xmin;
    xmin=xmax;
    xmax=tmp;
  }
  TVectorD *binLim=new TVectorD(nbinsX+1);
  Double_t first=xmin;
  Double_t last=xmax;
  Double_t binWidth=(last-first)/nbinsX;
  for (Int_t i=0; i<nbinsX+1; ++i){
    (*binLim)[i]=first+binWidth*(Double_t)i;
  }
  return binLim;
}

//_____________________________________________________________________________
TVectorD* AliEventPlaneHelper::MakeArbitraryBinning(const char* bins)
{
  //
  // Make arbitrary binning, bins separated by a ','
  //
  TString limits(bins);
  if (limits.IsNull()){
    //Error("AliEventPlaneHelper::MakeArbitraryBinning","Bin Limit string is empty, cannot add the variable");
    return 0x0;
  }
  
  TObjArray *arr=limits.Tokenize(",");
  Int_t nLimits=arr->GetEntries();
  if (nLimits<2){
    //Error("AliEventPlaneHelper::MakeArbitraryBinning","Need at leas 2 bin limits, cannot add the variable");
    delete arr;
    return 0x0;
  }
  
  TVectorD *binLimits=new TVectorD(nLimits);
  for (Int_t iLim=0; iLim<nLimits; ++iLim){
    (*binLimits)[iLim]=(static_cast<TObjString*>(arr->At(iLim)))->GetString().Atof();
  }
  
  delete arr;
  return binLimits;
}



//____________________________________________________________________________________
void AliEventPlaneHelper::InitFile(const Char_t* filename) {
  //
  // Open an existing ROOT file containing lists of histograms and initialize the global list pointer
  //
  //if(fgHistCaliFile) fgHistCaliFile->Close();
  TString histfilename="";
  if(fgHistCaliFile) histfilename = fgHistCaliFile->GetName();
  if(!histfilename.Contains(filename)) {
    fgHistCaliFile = TFile::Open(filename);    // open file only if not already open
    //fgHistCaliFile = new TFile(filename);    // open file only if not already open
  
    if(!fgHistCaliFile) {
      //cout << "GlobalMacros.C::GetHistogram() : File " << filename << " not opened!!" << endl;
      return;
    }
    if(fgHistCaliFile->IsZombie()) {
      //cout << "GlobalMacros.C::GetHistogram() : File " << filename << " not opened!!" << endl;
      return;
    }
    TList* list1 = fgHistCaliFile->GetListOfKeys();
    TKey* key1 = (TKey*)list1->At(0);
    fgHistCali = (TDirectoryFile*)key1->ReadObj();
  }
}




//____________________________________________________________________________________
void AliEventPlaneHelper::CloseFile() {
  //
  // Close the opened file
  //
  fgHistCali = 0x0;
  if(fgHistCaliFile && fgHistCaliFile->IsOpen()) fgHistCaliFile->Close();
}



//_________________________________________________________________
TChain* AliEventPlaneHelper::GetChain(const Char_t* filename, Int_t howMany, Int_t offset, Long64_t& entries,
                                  TChain* friendChain/*=0x0*/, const Char_t* friendChainFile/*=0x0*/) {
  //
  // read an ascii file containing a list of root files with reduced trees
  // and build a TChain
  //
  //cout << "Creating the data chain from " << filename << " ..." << flush; 
  TChain* chain = new TChain("DstTree");
  ifstream inBuf;
  inBuf.open(filename);
  Int_t index = 0;
  while(inBuf.good()) {
    Char_t str[512];
    inBuf.getline(str,512,'\n');
    
    if(index<offset) {++index; continue;}
    if(index>=offset+howMany) break;
    
    TString strstr = str;
    if(!strstr.Contains(".root")) continue;
    //cout << endl << "Adding file " << str << endl;
    chain->Add(str);
    if(friendChain) {
      TObjArray* arr = strstr.Tokenize("/");
      strstr.ReplaceAll(arr->At(arr->GetEntries()-1)->GetName(),friendChainFile);
      friendChain->Add(strstr.Data());
    }    
    ++index;
  }
  inBuf.close();
  entries = chain->GetEntries();
  Long64_t entriesFriend = (friendChain ? friendChain->GetEntries() : 0);
  //cout << "AliEventPlaneHelper::GetChain() Chain entries = " << entries << endl;
  if(friendChain)
    //cout << "AliEventPlaneHelper::GetChain() Friend chain entries = " << entriesFriend << endl;
  if(friendChain && (entries!=entriesFriend)) {
    //cout << "AliEventPlaneHelper::GetChain() The friend chain does not have the same number of entries as the main chain!!!" << endl;
    //cout << "                            Check it out and retry !!" << endl;
    return 0x0;
  }
  //cout << " done" << endl;
  return chain;
}



//_________________________________________________________________
TChain* AliEventPlaneHelper::GetChain(const Char_t* filename, Int_t howMany, Int_t offset, Long64_t& entries,
                                  const Char_t* friendpath, TChain* friendChain/*=0x0*/, const Char_t* friendChainFile/*=0x0*/) {
  //
  // read an ascii file containing a list of root files with reduced trees
  // and build a TChain
  //
  //cout << "Creating the data chain from " << filename << " ..." << flush;
  TChain* chain = new TChain("DstTree");
  ifstream inBuf;
  inBuf.open(filename);
  Int_t index = 0;
  while(inBuf.good()) {
    Char_t str[512];
    inBuf.getline(str,512,'\n');

    if(index<offset) {++index; continue;}
    if(index>=offset+howMany) break;

    TString strstr = str;
    if(!strstr.Contains(".root")) continue;
    //cout << endl << "Adding file " << str << endl;
    chain->Add(str);
    if(friendChain) {
      TObjArray* arr = strstr.Tokenize("/");
      friendChain->Add(Form("%s/%s/%s/%s", friendpath, arr->At(arr->GetEntries()-3)->GetName(), arr->At(arr->GetEntries()-2)->GetName(), friendChainFile));
    }
    ++index;
  }
  inBuf.close();
  entries = chain->GetEntries();
  Long64_t entriesFriend = (friendChain ? friendChain->GetEntries() : 0);
  //cout << "AliEventPlaneHelper::GetChain() Chain entries = " << entries << endl;
  if(friendChain)
    //cout << "AliEventPlaneHelper::GetChain() Friend chain entries = " << entriesFriend << endl;
  if(friendChain && (entries!=entriesFriend)) {
    //cout << "AliEventPlaneHelper::GetChain() The friend chain does not have the same number of entries as the main chain!!!" << endl;
    //cout << "                            Check it out and retry !!" << endl;
    return 0x0;
  }
  //cout << " done" << endl;
  return chain;
}





//__________________________________________________________________
Double_t * AliEventPlaneHelper::GetElements(TClonesArray* array, Int_t vector){

  //
  // Retrieve fElements from given TVectorT<double> (vector) in TClonesArray (array)
  //

  if(vector==-1) return 0x0;

  TVectorT<double>* vec  = (TVectorT<double>*) array->At(vector);

  return vec->GetMatrixArray();

}


//__________________________________________________________________
TArrayD * AliEventPlaneHelper::ArrayConversion(TClonesArray* array, Int_t nvectors){

  //
  // Convert a TClonesArray of TVector<double> to a TArrayD of TArrayDs 
  // 

  TArrayD * binLimits;
  binLimits = new TArrayD [nvectors];
  for(Int_t idim=0; idim<nvectors; idim++) {
    TVectorT<double> * tmp = (TVectorT<double>*) array->At(idim);
    binLimits[idim] = TArrayD(tmp->GetNrows(), tmp->GetMatrixArray());
  }


  return binLimits;
}


//__________________________________________________________________
TArrayD * AliEventPlaneHelper::AppendArray(TArrayD* array, TArrayD arr, Int_t length){

  //
  // Append a TArrayD to an array of TArrayDs 
  // 

    
  TArrayD * binLimits;
  binLimits = new TArrayD [length];
  
  for(Int_t iarr = 0; iarr < length-1; iarr++){
    binLimits[iarr] = array[iarr];
  }
  binLimits[length-1] = arr;

  return binLimits;
}



//____________________________________________________________________________________
Double_t * AliEventPlaneHelper::MakeBins(Int_t nbins, Double_t min, Double_t max) {

  Double_t * newbins;
  newbins = new Double_t [nbins+1];
  
  newbins[0] = min;
  
  Double_t xwidth = (max-min)/(Double_t (nbins));
  
  for(Int_t ibin=1; ibin<nbins+1; ibin++) 
    newbins[ibin] = newbins[ibin-1] + xwidth;
  
  
  return newbins;
}

//____________________________________________________________________________________
TObject* AliEventPlaneHelper::GetHistogram(const Char_t* listname, const Char_t* hname) {
  //
  // Retrieve a histogram from the list hlist
  //
  if(!fgHistCali) {
    //cout << "GlobalMacros.C::GetHistogram() : The main list was not initialized." << endl;
    //cout << "                   A ROOT file must pe initialized first!!" << endl;
    return 0x0;
  }
  if(fgHistCali->FindObject(listname)) return fgHistCali->FindObject(listname)->FindObject(hname);
  else return 0x0;
  //TKey* listKey = fgHistCali->FindKey(listname);
  //TDirectoryFile* hlist = (TDirectoryFile*)listKey->ReadObj();
  //TKey* key = hlist->FindKey(hname);
  //if(key==0) {/*cout<<"Histogram does not exist"<<endl;*/ return 0x0;}
  //return key->ReadObj();
}


//____________________________________________________________________________________
TList* AliEventPlaneHelper::GetHistogramList(const Char_t* listname) {
  //
  // Retrieve a histogram from the list hlist
  //
  if(!fgHistCali) {
    //cout << "GlobalMacros.C::GetHistogram() : The main list was not initialized." << endl;
    //cout << "                   A ROOT file must pe initialized first!!" << endl;
    return 0x0;
  }
  TKey* listKey = fgHistCali->FindKey(listname);
  TDirectoryFile* hdir = (TDirectoryFile*)listKey->ReadObj();
  return hdir->GetListOfKeys();
}


//_________________________________________________________________
THnF* AliEventPlaneHelper::AddHistogram( const Char_t* name, const Char_t* title, Int_t nDimensions,TAxis* binLimits){
  //
  // create a multi-dimensional histogram THnF with equal or variable bin widths
  //
  if(!binLimits) return 0x0;
  TString hname = name;

  TString titleStr(title);
  TObjArray* arr=titleStr.Tokenize(";");

  Double_t* xmin = new Double_t[nDimensions];
  Double_t* xmax = new Double_t[nDimensions];
  Int_t* nBins = new Int_t[nDimensions];
  for(Int_t idim=0;idim<nDimensions;++idim) {
    nBins[idim] = binLimits[idim].GetNbins();
    xmin[idim] = binLimits[idim].GetBinLowEdge(0);
    xmax[idim] = binLimits[idim].GetBinUpEdge(nBins[idim]);
  }

  //THnF* h=new THnF(hname.Data(),arr->At(0)->GetName(),nDimensions,nBins,xmin,xmax);
  THnF* h=new THnF(hname.Data(),title,nDimensions,nBins,xmin,xmax);
  for(Int_t idim=0;idim<nDimensions;++idim) {
    TAxis* axis=h->GetAxis(idim);
    axis->Set(nBins[idim], binLimits[idim].GetXbins()->GetArray());
  }

  h->Sumw2();

  delete [] xmin;
  delete [] xmax;
  delete [] nBins;
  //delete [] binLimits;

  return h;
}



//_________________________________________________________________
THnF* AliEventPlaneHelper::AddHistogram( const Char_t* name, const Char_t* title, Int_t nDimensions,TArrayD* binLimits){
  //
  // create a multi-dimensional histogram THnF with equal or variable bin widths
  //
  TString hname = name;

  TString titleStr(title);
  TObjArray* arr=titleStr.Tokenize(";");

  Double_t* xmin = new Double_t[nDimensions];
  Double_t* xmax = new Double_t[nDimensions];
  Int_t* nBins = new Int_t[nDimensions];
  for(Int_t idim=0;idim<nDimensions;++idim) {
    nBins[idim] = binLimits[idim].GetSize()-1;
    xmin[idim] = binLimits[idim][0];
    xmax[idim] = binLimits[idim][nBins[idim]];
  }

  //THnF* h=new THnF(hname.Data(),arr->At(0)->GetName(),nDimensions,nBins,xmin,xmax);
  THnF* h=new THnF(hname.Data(),title,nDimensions,nBins,xmin,xmax);
  for(Int_t idim=0;idim<nDimensions;++idim) {
    TAxis* axis=h->GetAxis(idim);
    axis->Set(nBins[idim], binLimits[idim].GetArray());
  }

  h->Sumw2();

  delete [] xmin;
  delete [] xmax;
  delete [] nBins;
  //delete [] binLimits;

  return h;
}

