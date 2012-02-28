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

//
// Generic Histogram container with support for groups and filling of groups by passing
// a vector of data
//
// Authors: 
//   Jens Wiechula <Jens.Wiechula@cern.ch> 
// 

#include <iostream>
#include <TH1.h>
#include <TH1F.h>
#include <TH2.h>
#include <TH3.h>
#include <TTree.h>
#include <TCollection.h>
#include <THashList.h>
#include <TString.h>
#include <TObjString.h>
#include <TObjArray.h>
#include <TFile.h>
#include <TError.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TROOT.h>
#include <TLegend.h>
#include <TKey.h>
#include <TAxis.h>
#include <TVirtualPS.h>
#include <TVectorD.h>

#include "AliDielectronHelper.h"
#include "AliDielectronHistosTaku.h"
#include "AliDielectronVarManager.h"

ClassImp(AliDielectronHistosTaku)


AliDielectronHistosTaku::AliDielectronHistosTaku() :
//   TCollection(),
  TNamed("AliDielectronHistosTaku","Dielectron Histogram Container"),
  fHistoList(),
  fTreeList(),
  fList(0x0),
  fTree(0x0),
  fReservedWords(new TString)
{
  //
  // Default constructor
  //
  fHistoList.SetOwner(kTRUE);
  fHistoList.SetName("Dielectron_Histos");
  fTreeList.SetOwner(kTRUE);
  fTreeList.SetName("Dielectron_Tree");
}

//_____________________________________________________________________________
AliDielectronHistosTaku::AliDielectronHistosTaku(const char* name, const char* title) :
//   TCollection(),
  TNamed(name, title),
  fHistoList(),
  fTreeList(),
  fList(0x0),
  fTree(0x0),
  fReservedWords(new TString)
{
  //
  // TNamed constructor
  //
  fHistoList.SetOwner(kTRUE);
  fHistoList.SetName(name);
  fTreeList.SetOwner(kTRUE);
  char tname[100];
  sprintf(tname,"%s_tree", name);
  fTreeList.SetName(tname);
}

//_____________________________________________________________________________
AliDielectronHistosTaku::~AliDielectronHistosTaku()
{
  //
  // Destructor
  //
  fHistoList.Clear();
  if (fList) fList->Clear();
  delete fReservedWords;
  fTreeList.Clear();
  if (fTree) fTree->Clear();
}

//_____________________________________________________________________________
void AliDielectronHistosTaku::UserHistogram(const char* histClass,const char *name, const char* title,
                                        Int_t nbinsX, Double_t xmin, Double_t xmax,
                                        UInt_t valTypeX, Bool_t logBinX)
{
  //
  // Default histogram creation 1D case
  //

  TVectorD *binLimX=0x0;
  
  if (logBinX) {
    binLimX=AliDielectronHelper::MakeLogBinning(nbinsX, xmin, xmax);
  } else {
    binLimX=AliDielectronHelper::MakeLinBinning(nbinsX, xmin, xmax);
  }

  UserHistogram(histClass,name,title,binLimX,valTypeX);
}

//_____________________________________________________________________________
void AliDielectronHistosTaku::UserHistogram(const char* histClass,const char *name, const char* title,
                                        Int_t nbinsX, Double_t xmin, Double_t xmax,
                                        Int_t nbinsY, Double_t ymin, Double_t ymax,
                                        UInt_t valTypeX, UInt_t valTypeY,
                                        Bool_t logBinX, Bool_t logBinY)
{
  //
  // Default histogram creation 2D case
  //
  if (!IsHistogramOk(histClass,name)) return;

  TVectorD *binLimX=0x0;
  TVectorD *binLimY=0x0;
  
  if (logBinX) {
    binLimX=AliDielectronHelper::MakeLogBinning(nbinsX, xmin, xmax);
  } else {
    binLimX=AliDielectronHelper::MakeLinBinning(nbinsX, xmin, xmax);
  }
  if (logBinY) {
    binLimY=AliDielectronHelper::MakeLogBinning(nbinsY, ymin, ymax);
  } else {
    binLimY=AliDielectronHelper::MakeLinBinning(nbinsY, ymin, ymax);
  }
  
  UserHistogram(histClass,name,title,binLimX,binLimY,valTypeX,valTypeY);
}


//_____________________________________________________________________________
void AliDielectronHistosTaku::UserHistogram(const char* histClass,const char *name, const char* title,
                                        Int_t nbinsX, Double_t xmin, Double_t xmax,
                                        Int_t nbinsY, Double_t ymin, Double_t ymax,
                                        Int_t nbinsZ, Double_t zmin, Double_t zmax,
                                        UInt_t valTypeX, UInt_t valTypeY, UInt_t valTypeZ,
                                        Bool_t logBinX, Bool_t logBinY, Bool_t logBinZ)
{
  //
  // Default histogram creation 3D case
  //
  if (!IsHistogramOk(histClass,name)) return;

  TVectorD *binLimX=0x0;
  TVectorD *binLimY=0x0;
  TVectorD *binLimZ=0x0;
  
  if (logBinX) {
    binLimX=AliDielectronHelper::MakeLogBinning(nbinsX, xmin, xmax);
  } else {
    binLimX=AliDielectronHelper::MakeLinBinning(nbinsX, xmin, xmax);
  }
  
  if (logBinY) {
    binLimY=AliDielectronHelper::MakeLogBinning(nbinsY, ymin, ymax);
  } else {
    binLimY=AliDielectronHelper::MakeLinBinning(nbinsY, ymin, ymax);
  }
  
  if (logBinZ) {
    binLimZ=AliDielectronHelper::MakeLogBinning(nbinsZ, zmin, zmax);
  } else {
    binLimZ=AliDielectronHelper::MakeLinBinning(nbinsZ, zmin, zmax);
  }

  UserHistogram(histClass,name,title,binLimX,binLimY,binLimZ,valTypeX,valTypeY,valTypeZ);
}

//_____________________________________________________________________________
void AliDielectronHistosTaku::UserHistogram(const char* histClass,const char *name, const char* title,
                                        const char* binning,
                                        UInt_t valTypeX)
{
  //
  // Histogram creation 1D case with arbitraty binning
  //

  TVectorD *binLimX=AliDielectronHelper::MakeArbitraryBinning(binning);
  UserHistogram(histClass,name,title,binLimX,valTypeX);
}

//_____________________________________________________________________________
void AliDielectronHistosTaku::UserHistogram(const char* histClass,const char *name, const char* title,
                                        const TVectorD * const binsX,
                                        UInt_t valTypeX/*=kNoAutoFill*/)
{
  //
  // Histogram creation 1D case with arbitraty binning X
  // the TVectorD is assumed to be surplus after the creation and will be deleted!!!
  //

  Bool_t isOk=kTRUE;
  isOk&=IsHistogramOk(histClass,name);
  isOk&=(binsX!=0x0);

  if (isOk){
    TH1* hist=new TH1F(name,title,binsX->GetNrows()-1,binsX->GetMatrixArray());
  
    Bool_t isReserved=fReservedWords->Contains(histClass);
    if (isReserved)
      UserHistogramReservedWords(histClass, hist, valTypeX);
    else
      UserHistogram(histClass, hist, valTypeX);
  }
  
  delete binsX;
}

//_____________________________________________________________________________
void AliDielectronHistosTaku::UserHistogram(const char* histClass,const char *name, const char* title,
                                        const TVectorD * const binsX, const TVectorD * const binsY,
                                        UInt_t valTypeX/*=kNoAutoFill*/, UInt_t valTypeY/*=0*/)
{
  //
  // Histogram creation 1D case with arbitraty binning X
  // the TVectorD is assumed to be surplus after the creation and will be deleted!!!
  //

  Bool_t isOk=kTRUE;
  isOk&=IsHistogramOk(histClass,name);
  isOk&=(binsX!=0x0);
  isOk&=(binsY!=0x0);

  if (isOk){
    TH1* hist=new TH2F(name,title,
                       binsX->GetNrows()-1,binsX->GetMatrixArray(),
                       binsY->GetNrows()-1,binsY->GetMatrixArray());
  
    Bool_t isReserved=fReservedWords->Contains(histClass);
    if (isReserved)
      UserHistogramReservedWords(histClass, hist, valTypeX+100*valTypeY);
    else
      UserHistogram(histClass, hist, valTypeX+100*valTypeY);
  }
  
  delete binsX;
  delete binsY;
  
}

//_____________________________________________________________________________
void AliDielectronHistosTaku::UserHistogram(const char* histClass,const char *name, const char* title,
                                        const TVectorD * const binsX, const TVectorD * const binsY, const TVectorD * const binsZ,
                                        UInt_t valTypeX/*=kNoAutoFill*/, UInt_t valTypeY/*=0*/, UInt_t valTypeZ/*=0*/)
{
  //
  // Histogram creation 1D case with arbitraty binning X
  // the TVectorD is assumed to be surplus after the creation and will be deleted!!!
  //

  Bool_t isOk=kTRUE;
  isOk&=IsHistogramOk(histClass,name);
  isOk&=(binsX!=0x0);
  isOk&=(binsY!=0x0);
  isOk&=(binsZ!=0x0);
  
  if (isOk){
    TH1* hist=new TH3F(name,title,
                       binsX->GetNrows()-1,binsX->GetMatrixArray(),
                       binsY->GetNrows()-1,binsY->GetMatrixArray(),
                       binsZ->GetNrows()-1,binsZ->GetMatrixArray());
  
    Bool_t isReserved=fReservedWords->Contains(histClass);
    if (isReserved)
      UserHistogramReservedWords(histClass, hist, valTypeX+100*valTypeY+10000*valTypeZ);
    else
      UserHistogram(histClass, hist, valTypeX+100*valTypeY+10000*valTypeZ);
  }
  
  delete binsX;
  delete binsY;
  delete binsZ;
}

//_____________________________________________________________________________
void AliDielectronHistosTaku::UserHistogram(const char* histClass, TH1* hist, UInt_t valTypes)
{
  //
  // Add any type of user histogram
  //

  //special case for the calss Pair. where histograms will be created for all pair classes
  Bool_t isReserved=fReservedWords->Contains(histClass);
  if (isReserved) {
    UserHistogramReservedWords(histClass, hist, valTypes);
    return;
  }
  
  if (!IsHistogramOk(histClass,hist->GetName())) return;
  
  THashList *classTable=(THashList*)fHistoList.FindObject(histClass);
  hist->SetDirectory(0);
  hist->SetUniqueID(valTypes);
  classTable->Add(hist);
}

//_____________________________________________________________________________
void AliDielectronHistosTaku::AddClass(const char* histClass)
{
  //
  // Add a class of histograms
  // Several classes can be added by separating them by a ';' e.g. 'class1;class2;class3'
  //
  TString hists(histClass);
  TObjArray *arr=hists.Tokenize(";");
  TIter next(arr);
  TObject *o=0;
  while ( (o=next()) ){
    if (fHistoList.FindObject(o->GetName())){
      Warning("AddClass","Cannot create class '%s' it already exists.",histClass);
      continue;
    }
    if (fReservedWords->Contains(o->GetName())){
      Error("AddClass","Pair is a reserved word, please use another name");
      continue;
    }
    THashList *table=new THashList;
    table->SetOwner(kTRUE);
    table->SetName(o->GetName());
    fHistoList.Add(table);
  }
  delete arr;
}

//_____________________________________________________________________________
void AliDielectronHistosTaku::Fill(const char* histClass, const char* name, Double_t xval)
{
  //
  // Fill function 1D case
  //
  THashList *classTable=(THashList*)fHistoList.FindObject(histClass);
  TH1* hist=0;
  if (!classTable || !(hist=(TH1*)classTable->FindObject(name)) ){
    Warning("Fill","Cannot fill histogram. Either class '%s' or histogram '%s' not existing.",histClass,name);
    return;
  }
  hist->Fill(xval);
}

//_____________________________________________________________________________
void AliDielectronHistosTaku::Fill(const char* histClass, const char* name, Double_t xval, Double_t yval)
{
  //
  // Fill function 2D case
  //
  THashList *classTable=(THashList*)fHistoList.FindObject(histClass);
  TH2* hist=0;
  if (!classTable || !(hist=(TH2*)classTable->FindObject(name)) ){
    Warning("UserHistogram","Cannot fill histogram. Either class '%s' or histogram '%s' not existing.",histClass,name);
    return;
  }
  hist->Fill(xval,yval);
}

//_____________________________________________________________________________
void AliDielectronHistosTaku::Fill(const char* histClass, const char* name, Double_t xval, Double_t yval, Double_t zval)
{
  //
  // Fill function 3D case
  //
  THashList *classTable=(THashList*)fHistoList.FindObject(histClass);
  TH3* hist=0;
  if (!classTable || !(hist=(TH3*)classTable->FindObject(name)) ){
    Warning("UserHistogram","Cannot fill histogram. Either class '%s' or histogram '%s' not existing.",histClass,name);
    return;
  }
  hist->Fill(xval,yval,zval);
}

//_____________________________________________________________________________
void AliDielectronHistosTaku::FillClass(const char* histClass, Int_t nValues, const Double_t *values)
{
  //
  // Fill class 'histClass' (by name)
  //
  
  THashList *classTable=(THashList*)fHistoList.FindObject(histClass);
  if (!classTable){
    Warning("FillClass","Cannot fill class '%s' its not defined.",histClass);
    return;
  }
  
  TIter nextHist(classTable);
  TH1 *hist=0;
  while ( (hist=(TH1*)nextHist()) ){
    UInt_t valueTypes=hist->GetUniqueID();
    if (valueTypes==(UInt_t)kNoAutoFill) continue;
    UInt_t value1=valueTypes%100;        //last two digits
    UInt_t value2=valueTypes/100%100;    //second last two digits
    UInt_t value3=valueTypes/10000%100;  //third last two digits
    if (value1>=(UInt_t)nValues||value2>=(UInt_t)nValues||value3>=(UInt_t)nValues) {
      Warning("FillClass","One of the values is out of range. Not filling histogram '%s/%s'.", histClass, hist->GetName());
      continue;
    }
    switch (hist->GetDimension()){
    case 1:
      hist->Fill(values[value1]);
      break;
    case 2:
      ((TH2*)hist)->Fill(values[value1],values[value2]);
      break;
    case 3:
      ((TH3*)hist)->Fill(values[value1],values[value2],values[value3]);
      break;
    }
  }
}

//_____________________________________________________________________________
// void AliDielectronHistosTaku::FillClass(const char* histClass, const TVectorD &vals)
// {
//   //
//   //
//   //
//   FillClass(histClass, vals.GetNrows(), vals.GetMatrixArray());
// }

//_____________________________________________________________________________
void AliDielectronHistosTaku::UserHistogramReservedWords(const char* histClass, const TH1 *hist, UInt_t valTypes)
{
  //
  // Creation of histogram for all pair types
  //
  TString title(hist->GetTitle());
  // Same Event Like Sign
  TIter nextClass(&fHistoList);
  THashList *l=0;
  while ( (l=static_cast<THashList*>(nextClass())) ){
    TString name(l->GetName());
    if (name.Contains(histClass)){
      TH1 *h=static_cast<TH1*>(hist->Clone());
      h->SetDirectory(0);
      h->SetTitle(Form("%s %s",title.Data(),l->GetName()));
      UserHistogram(l->GetName(),h,valTypes);
    }
  }
  delete hist;
}

//_____________________________________________________________________________
void AliDielectronHistosTaku::DumpToFile(const char* file)
{
  //
  // Dump the histogram list to a newly created root file
  //
  TFile f(file,"recreate");
  fHistoList.Write(fHistoList.GetName(),TObject::kSingleKey);
  f.Close();
}

//_____________________________________________________________________________
TH1* AliDielectronHistosTaku::GetHistogram(const char* histClass, const char* name) const
{
  //
  // return histogram 'name' in 'histClass'
  //
  THashList *classTable=(THashList*)fHistoList.FindObject(histClass);
  if (!classTable) return 0x0;
  return (TH1*)classTable->FindObject(name);
}

//_____________________________________________________________________________
TH1* AliDielectronHistosTaku::GetHistogram(const char* cutClass, const char* histClass, const char* name) const
{
  //
  // return histogram from list of list of histograms
  // this function is thought for retrieving histograms if a list of AliDielectronHistosTaku is set
  //
  
  if (!fList) return 0x0;
  THashList *h=dynamic_cast<THashList*>(fList->FindObject(cutClass));
  if (!h)return 0x0;
  THashList *classTable=dynamic_cast<THashList*>(h->FindObject(histClass));
  if (!classTable) return 0x0;
  return (TH1*)classTable->FindObject(name);
}

//_____________________________________________________________________________
void AliDielectronHistosTaku::Draw(const Option_t* option)
{
  //
  // Draw histograms
  //

  TString drawStr(option);
  TObjArray *arr=drawStr.Tokenize(";");
  arr->SetOwner();
  TIter nextOpt(arr);

  TString drawClasses;
  TObjString *ostr=0x0;

  TString currentOpt;
  TString testOpt;
  while ( (ostr=(TObjString*)nextOpt()) ){
    currentOpt=ostr->GetString();
    currentOpt.Remove(TString::kBoth,'\t');
    currentOpt.Remove(TString::kBoth,' ');

    testOpt="classes=";
    if ( currentOpt.Contains(testOpt.Data()) ){
      drawClasses=currentOpt(testOpt.Length(),currentOpt.Length());
    }
  }

  delete arr;
  drawStr.ToLower();
  //optionsfList
//   Bool_t same=drawOpt.Contains("same"); //FIXME not yet implemented

  TCanvas *c=0x0;
  if (gVirtualPS) {
    if (!gPad){
      Error("Draw","When writing to a file you have to create a canvas before opening the file!!!");
      return;
    }
    c=gPad->GetCanvas();
    c->cd();
//     c=new TCanvas;
  }
  
  TIter nextClass(&fHistoList);
  THashList *classTable=0;
//   Bool_t first=kTRUE;
  while ( (classTable=(THashList*)nextClass()) ){
    //test classes option
    if (!drawClasses.IsNull() && !drawClasses.Contains(classTable->GetName())) continue;
    //optimised division
    Int_t nPads = classTable->GetEntries();
    Int_t nCols = (Int_t)TMath::Ceil( TMath::Sqrt(nPads) );
    Int_t nRows = (Int_t)TMath::Ceil( (Double_t)nPads/(Double_t)nCols );

    //create canvas
    if (!gVirtualPS){
      TString canvasName;
      canvasName.Form("c%s_%s",GetName(),classTable->GetName());
      c=(TCanvas*)gROOT->FindObject(canvasName.Data());
      if (!c) c=new TCanvas(canvasName.Data(),Form("%s: %s",GetName(),classTable->GetName()));
      c->Clear();
    } else {
//       if (first){
//         first=kFALSE;
//         if (nPads>1) gVirtualPS->NewPage();
//       } else {
        if (nPads>1) c->Clear();
//       }
    }
    if (nCols>1||nRows>1) c->Divide(nCols,nRows);
    
    //loop over histograms and draw them
    TIter nextHist(classTable);
    Int_t iPad=0;
    TH1 *h=0;
    while ( (h=(TH1*)nextHist()) ){
      TString drawOpt;
      if ( (h->InheritsFrom(TH2::Class())) ) drawOpt="colz";
      if (nCols>1||nRows>1) c->cd(++iPad);
      if ( TMath::Abs(h->GetXaxis()->GetBinWidth(1)-h->GetXaxis()->GetBinWidth(2))>1e-10 ) gPad->SetLogx();
      if ( TMath::Abs(h->GetYaxis()->GetBinWidth(1)-h->GetYaxis()->GetBinWidth(2))>1e-10 ) gPad->SetLogy();
      if ( TMath::Abs(h->GetZaxis()->GetBinWidth(1)-h->GetZaxis()->GetBinWidth(2))>1e-10 ) gPad->SetLogz();
      TString histOpt=h->GetOption();
      histOpt.ToLower();
      if (histOpt.Contains("logx")) gPad->SetLogx();
      if (histOpt.Contains("logy")) gPad->SetLogy();
      if (histOpt.Contains("logz")) gPad->SetLogz();
      histOpt.ReplaceAll("logx","");
      histOpt.ReplaceAll("logy","");
      histOpt.ReplaceAll("logz","");
      h->Draw(drawOpt.Data());
    }
    if (gVirtualPS) {
      c->Update();
    }
    
  }
//   if (gVirtualPS) delete c;
}

//_____________________________________________________________________________
void AliDielectronHistosTaku::Print(const Option_t* option) const
{
  //
  // Print classes and histograms
  //
  TString optString(option);

  if (optString.IsNull()) PrintStructure();



}

//_____________________________________________________________________________
void AliDielectronHistosTaku::PrintStructure() const
{
  //
  // Print classes and histograms in the class to stdout
  //
  if (!fList){
    TIter nextClass(&fHistoList);
    THashList *classTable=0;
    while ( (classTable=(THashList*)nextClass()) ){
      TIter nextHist(classTable);
      TObject *o=0;
      printf("+ %s\n",classTable->GetName());
      while ( (o=nextHist()) )
        printf("| ->%s\n",o->GetName());
    }
  } else {
    TIter nextCutClass(fList);
    THashList *cutClass=0x0;
    while ( (cutClass=(THashList*)nextCutClass()) ) {
      printf("+ %s\n",cutClass->GetName());
      TIter nextClass(cutClass);
      THashList *classTable=0;
      while ( (classTable=(THashList*)nextClass()) ){
        TIter nextHist(classTable);
        TObject *o=0;
        printf("|  + %s\n",classTable->GetName());
        while ( (o=nextHist()) )
          printf("|  | ->%s\n",o->GetName());
      }
      
    }
  }
}

//_____________________________________________________________________________
void AliDielectronHistosTaku::SetHistogramList(THashList &list, Bool_t setOwner/*=kTRUE*/)
{
  //
  // set histogram classes and histograms to this instance. It will take onwnership!
  //
  ResetHistogramList();
  TString name(GetName());
  if (name == "AliDielectronHistosTaku") SetName(list.GetName());
  TIter next(&list);
  TObject *o;
  while ( (o=next()) ){
    fHistoList.Add(o);
  }
  if (setOwner){
    list.SetOwner(kFALSE);
    fHistoList.SetOwner(kTRUE);
  } else {
    fHistoList.SetOwner(kFALSE);
  }
}

//_____________________________________________________________________________
Bool_t AliDielectronHistosTaku::SetCutClass(const char* cutClass)
{
  //
  // Assign histogram list according to cutClass
  //

  if (!fList) return kFALSE;
  ResetHistogramList();
  THashList *h=dynamic_cast<THashList*>(fList->FindObject(cutClass));
  if (!h) {
    Warning("SetCutClass","cutClass '%s' not found", cutClass);
    return kFALSE;
  }
  SetHistogramList(*h,kFALSE);
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliDielectronHistosTaku::IsHistogramOk(const char* histClass, const char* name)
{
  //
  // check whether the histogram class exists and the histogram itself does not exist yet
  //
  Bool_t isReserved=fReservedWords->Contains(histClass);
  if (!fHistoList.FindObject(histClass)&&!isReserved){
    Warning("IsHistogramOk","Cannot create histogram. Class '%s' not defined. Please create it using AddClass before.",histClass);
    return kFALSE;
  }
  if (GetHistogram(histClass,name)){
    Warning("IsHistogramOk","Cannot create histogram '%s' in class '%s': It already exists!",name,histClass);
    return kFALSE;
  }
  return kTRUE;
}

// //_____________________________________________________________________________
// TIterator* AliDielectronHistosTaku::MakeIterator(Bool_t dir) const
// {
//   //
//   //
//   //
//   return new TListIter(&fHistoList, dir);
// }

//_____________________________________________________________________________
void AliDielectronHistosTaku::ReadFromFile(const char* file)
{
  //
  // Read histos from file
  //
  TFile f(file);
  TIter nextKey(f.GetListOfKeys());
  TKey *key=0;
  while ( (key=(TKey*)nextKey()) ){
    TObject *o=f.Get(key->GetName());
    THashList *list=dynamic_cast<THashList*>(o);
    if (!list) continue;
    SetHistogramList(*list);
    break;
  }
  f.Close();
}

//_____________________________________________________________________________
void AliDielectronHistosTaku::DrawSame(const char* histName, const Option_t *opt)
{
  //
  // Draw all histograms with the same name into one canvas
  // if option contains 'leg' a legend will be created with the class name as caption
  // if option contains 'can' a new canvas is created
  //

  TString optString(opt);
  optString.ToLower();
  Bool_t optLeg=optString.Contains("leg");
  Bool_t optCan=optString.Contains("can");

  TLegend *leg=0;
  TCanvas *c=0;
  if (optCan){
    c=(TCanvas*)gROOT->FindObject(Form("c%s",histName));
    if (!c) c=new TCanvas(Form("c%s",histName),Form("All '%s' histograms",histName));
    c->Clear();
    c->cd();
  }

  if (optLeg) leg=new TLegend(.8,.3,.99,.9);
  
  Int_t i=0;
  TIter next(&fHistoList);
  THashList *classTable=0;
  Double_t max=-1e10;
  TH1 *hFirst=0x0;
  while ( (classTable=(THashList*)next()) ){
    if ( TH1 *h=(TH1*)classTable->FindObject(histName) ){
      if (i==0) hFirst=h;
      h->SetLineColor(i+1);
      h->SetMarkerColor(i+1);
      h->Draw(i>0?"same":"");
      if (leg) leg->AddEntry(h,classTable->GetName(),"lp");
      ++i;
      max=TMath::Max(max,h->GetMaximum());
    }
  }
  if (leg){
    leg->SetFillColor(10);
    leg->SetY1(.9-i*.05);
    leg->Draw();
  }
  if (hFirst&&(hFirst->GetYaxis()->GetXmax()<max)){
    hFirst->SetMaximum(max);
  }
}

//_____________________________________________________________________________
void AliDielectronHistosTaku::SetReservedWords(const char* words)
{
  //
  // set reserved words
  //
  
  (*fReservedWords)=words;
}

//_____________________________________________________________________________
void AliDielectronHistosTaku::UserTree(const char* name, const char *title)
{
  fTree = new TTree(name, title);
}

//_____________________________________________________________________________
void AliDielectronHistosTaku::SetReserveVariableInTree(UInt_t valTypes) 
{
  switch(valTypes){
  case AliDielectronVarManager::kPx:       
    fTree->Branch("kPx",&fgData[AliDielectronVarManager::kPx],"kPx/D"); break;
  case AliDielectronVarManager::kPy:       
    fTree->Branch("kPy",&fgData[AliDielectronVarManager::kPy],"kPy/D"); break;
  case AliDielectronVarManager::kPz:       
    fTree->Branch("kPz",&fgData[AliDielectronVarManager::kPz],"kPz/D"); break;
  case AliDielectronVarManager::kPt:       
    fTree->Branch("kPt",&fgData[AliDielectronVarManager::kPt],"kPt/D"); break;
  case AliDielectronVarManager::kP:        
    fTree->Branch("kP",&fgData[AliDielectronVarManager::kP],"kP/D"); break;
  case AliDielectronVarManager::kXv:       
    fTree->Branch("kXv",&fgData[AliDielectronVarManager::kXv],"kXv/D"); break;
  case AliDielectronVarManager::kYv:       
    fTree->Branch("kYv",&fgData[AliDielectronVarManager::kYv],"kYv/D"); break;
  case AliDielectronVarManager::kZv:       
    fTree->Branch("kZv",&fgData[AliDielectronVarManager::kZv],"kZv/D"); break;
  case AliDielectronVarManager::kOneOverPt:        
    fTree->Branch("kOneOverPt",&fgData[AliDielectronVarManager::kOneOverPt],"kOneOverPt/D"); break;
  case AliDielectronVarManager::kPhi:      
    fTree->Branch("kPhi",&fgData[AliDielectronVarManager::kPhi],"kPhi/D"); break;
  case AliDielectronVarManager::kTheta:    
    fTree->Branch("kTheta",&fgData[AliDielectronVarManager::kTheta],"kTheta/D"); break;
  case AliDielectronVarManager::kEta:      
    fTree->Branch("kEta",&fgData[AliDielectronVarManager::kEta],"kEta/D"); break;
  case AliDielectronVarManager::kY:        
    fTree->Branch("kY",&fgData[AliDielectronVarManager::kY],"kY/D"); break;
  case AliDielectronVarManager::kE:        
    fTree->Branch("kE",&fgData[AliDielectronVarManager::kE],"kE/D"); break;
  case AliDielectronVarManager::kM:        
    fTree->Branch("kM",&fgData[AliDielectronVarManager::kM],"kM/D"); break;
  case AliDielectronVarManager::kCharge:   
    fTree->Branch("kCharge",&fgData[AliDielectronVarManager::kCharge],"kCharge/D"); break;
  case AliDielectronVarManager::kNclsITS:  
    fTree->Branch("kNclsITS",&fgData[AliDielectronVarManager::kNclsITS],"kNclsITS/D"); break;
  case AliDielectronVarManager::kNclsTPC:  
    fTree->Branch("kNclsTPC",&fgData[AliDielectronVarManager::kNclsTPC],"kNclsTPC/D"); break;
  case AliDielectronVarManager::kNclsTPCiter1:     
    fTree->Branch("kNclsTPCiter1",&fgData[AliDielectronVarManager::kNclsTPCiter1],"kNclsTPCiter1/D"); break;
  case AliDielectronVarManager::kNFclsTPC:         
    fTree->Branch("kNFclsTPC",&fgData[AliDielectronVarManager::kNFclsTPC],"kNFclsTPC/D"); break;
  case AliDielectronVarManager::kNFclsTPCr:        
    fTree->Branch("kNFclsTPCr",&fgData[AliDielectronVarManager::kNFclsTPCr],"kNFclsTPCr/D"); break;
  case AliDielectronVarManager::kNFclsTPCrFrac:    
    fTree->Branch("kNFclsTPCrFrac",&fgData[AliDielectronVarManager::kNFclsTPCrFrac],"kNFclsTPCrFrac/D"); break;
  case AliDielectronVarManager::kTPCsignalN:       
    fTree->Branch("kTPCsignalN",&fgData[AliDielectronVarManager::kTPCsignalN],"kTPCsignalN/D"); break;
  case AliDielectronVarManager::kTPCsignalNfrac:   
    fTree->Branch("kTPCsignalNfrac",&fgData[AliDielectronVarManager::kTPCsignalNfrac],"kTPCsignalNfrac/D"); break;
  case AliDielectronVarManager::kTPCchi2Cl:        
    fTree->Branch("kTPCchi2Cl",&fgData[AliDielectronVarManager::kTPCchi2Cl],"kTPCchi2Cl/D"); break;
  case AliDielectronVarManager::kTrackStatus:      
    fTree->Branch("kTrackStatus",&fgData[AliDielectronVarManager::kTrackStatus],"kTrackStatus/D"); break;
  case AliDielectronVarManager::kNclsTRD:  
    fTree->Branch("kNclsTRD",&fgData[AliDielectronVarManager::kNclsTRD],"kNclsTRD/D"); break;
  case AliDielectronVarManager::kTRDntracklets:    
    fTree->Branch("kTRDntracklets",&fgData[AliDielectronVarManager::kTRDntracklets],"kTRDntracklets/D"); break;
  case AliDielectronVarManager::kTRDpidQuality:    
    fTree->Branch("kTRDpidQuality",&fgData[AliDielectronVarManager::kTRDpidQuality],"kTRDpidQuality/D"); break;
  case AliDielectronVarManager::kTRDprobEle:       
    fTree->Branch("kTRDprobEle",&fgData[AliDielectronVarManager::kTRDprobEle],"kTRDprobEle/D"); break;
  case AliDielectronVarManager::kTRDprobPio:       
    fTree->Branch("kTRDprobPio",&fgData[AliDielectronVarManager::kTRDprobPio],"kTRDprobPio/D"); break;
  case AliDielectronVarManager::kImpactParXY:      
    fTree->Branch("kImpactParXY",&fgData[AliDielectronVarManager::kImpactParXY],"kImpactParXY/D"); break;
  case AliDielectronVarManager::kImpactParZ:       
    fTree->Branch("kImpactParZ",&fgData[AliDielectronVarManager::kImpactParZ],"kImpactParZ/D"); break;
  case AliDielectronVarManager::kTrackLength:      
    fTree->Branch("kTrackLength",&fgData[AliDielectronVarManager::kTrackLength],"kTrackLength/D"); break;
  case AliDielectronVarManager::kPdgCode:  
    fTree->Branch("kPdgCode",&fgData[AliDielectronVarManager::kPdgCode],"kPdgCode/D"); break;
  case AliDielectronVarManager::kPdgCodeMother:    
    fTree->Branch("kPdgCodeMother",&fgData[AliDielectronVarManager::kPdgCodeMother],"kPdgCodeMother/D"); break;
  case AliDielectronVarManager::kPdgCodeGrandMother:       
    fTree->Branch("kPdgCodeGrandMother",&fgData[AliDielectronVarManager::kPdgCodeGrandMother],"kPdgCodeGrandMother/D"); break;
  case AliDielectronVarManager::kNumberOfDaughters:        
    fTree->Branch("kNumberOfDaughters",&fgData[AliDielectronVarManager::kNumberOfDaughters],"kNumberOfDaughters/D"); break;
  case AliDielectronVarManager::kHaveSameMother:   
    fTree->Branch("kHaveSameMother",&fgData[AliDielectronVarManager::kHaveSameMother],"kHaveSameMother/D"); break;
  case AliDielectronVarManager::kIsJpsiPrimary:    
    fTree->Branch("kIsJpsiPrimary",&fgData[AliDielectronVarManager::kIsJpsiPrimary],"kIsJpsiPrimary/D"); break;
  case AliDielectronVarManager::kITSsignal:        
    fTree->Branch("kITSsignal",&fgData[AliDielectronVarManager::kITSsignal],"kITSsignal/D"); break;
  case AliDielectronVarManager::kITSsignalSSD1:    
    fTree->Branch("kITSsignalSSD1",&fgData[AliDielectronVarManager::kITSsignalSSD1],"kITSsignalSSD1/D"); break;
  case AliDielectronVarManager::kITSsignalSSD2:    
    fTree->Branch("kITSsignalSSD2",&fgData[AliDielectronVarManager::kITSsignalSSD2],"kITSsignalSSD2/D"); break;
  case AliDielectronVarManager::kITSsignalSDD1:    
    fTree->Branch("kITSsignalSDD1",&fgData[AliDielectronVarManager::kITSsignalSDD1],"kITSsignalSDD1/D"); break;
  case AliDielectronVarManager::kITSsignalSDD2:    
    fTree->Branch("kITSsignalSDD2",&fgData[AliDielectronVarManager::kITSsignalSDD2],"kITSsignalSDD2/D"); break;
  case AliDielectronVarManager::kITSclusterMap:    
    fTree->Branch("kITSclusterMap",&fgData[AliDielectronVarManager::kITSclusterMap],"kITSclusterMap/D"); break;
  case AliDielectronVarManager::kITSnSigmaEle:     
    fTree->Branch("kITSnSigmaEle",&fgData[AliDielectronVarManager::kITSnSigmaEle],"kITSnSigmaEle/D"); break;
  case AliDielectronVarManager::kITSnSigmaPio:     
    fTree->Branch("kITSnSigmaPio",&fgData[AliDielectronVarManager::kITSnSigmaPio],"kITSnSigmaPio/D"); break;
  case AliDielectronVarManager::kITSnSigmaMuo:     
    fTree->Branch("kITSnSigmaMuo",&fgData[AliDielectronVarManager::kITSnSigmaMuo],"kITSnSigmaMuo/D"); break;
  case AliDielectronVarManager::kITSnSigmaKao:     
    fTree->Branch("kITSnSigmaKao",&fgData[AliDielectronVarManager::kITSnSigmaKao],"kITSnSigmaKao/D"); break;
  case AliDielectronVarManager::kITSnSigmaPro:     
    fTree->Branch("kITSnSigmaPro",&fgData[AliDielectronVarManager::kITSnSigmaPro],"kITSnSigmaPro/D"); break;
  case AliDielectronVarManager::kPIn:      
    fTree->Branch("kPIn",&fgData[AliDielectronVarManager::kPIn],"kPIn/D"); break;
  case AliDielectronVarManager::kTPCsignal:        
    fTree->Branch("kTPCsignal",&fgData[AliDielectronVarManager::kTPCsignal],"kTPCsignal/D"); break;
  case AliDielectronVarManager::kTOFsignal:        
    fTree->Branch("kTOFsignal",&fgData[AliDielectronVarManager::kTOFsignal],"kTOFsignal/D"); break;
  case AliDielectronVarManager::kTOFbeta:  
    fTree->Branch("kTOFbeta",&fgData[AliDielectronVarManager::kTOFbeta],"kTOFbeta/D"); break;
  case AliDielectronVarManager::kTPCnSigmaEle:     
    fTree->Branch("kTPCnSigmaEle",&fgData[AliDielectronVarManager::kTPCnSigmaEle],"kTPCnSigmaEle/D"); break;
  case AliDielectronVarManager::kTPCnSigmaPio:     
    fTree->Branch("kTPCnSigmaPio",&fgData[AliDielectronVarManager::kTPCnSigmaPio],"kTPCnSigmaPio/D"); break;
  case AliDielectronVarManager::kTPCnSigmaMuo:     
    fTree->Branch("kTPCnSigmaMuo",&fgData[AliDielectronVarManager::kTPCnSigmaMuo],"kTPCnSigmaMuo/D"); break;
  case AliDielectronVarManager::kTPCnSigmaKao:     
    fTree->Branch("kTPCnSigmaKao",&fgData[AliDielectronVarManager::kTPCnSigmaKao],"kTPCnSigmaKao/D"); break;
  case AliDielectronVarManager::kTPCnSigmaPro:     
    fTree->Branch("kTPCnSigmaPro",&fgData[AliDielectronVarManager::kTPCnSigmaPro],"kTPCnSigmaPro/D"); break;
  case AliDielectronVarManager::kTOFnSigmaEle:     
    fTree->Branch("kTOFnSigmaEle",&fgData[AliDielectronVarManager::kTOFnSigmaEle],"kTOFnSigmaEle/D"); break;
  case AliDielectronVarManager::kTOFnSigmaPio:     
    fTree->Branch("kTOFnSigmaPio",&fgData[AliDielectronVarManager::kTOFnSigmaPio],"kTOFnSigmaPio/D"); break;
  case AliDielectronVarManager::kTOFnSigmaMuo:     
    fTree->Branch("kTOFnSigmaMuo",&fgData[AliDielectronVarManager::kTOFnSigmaMuo],"kTOFnSigmaMuo/D"); break;
  case AliDielectronVarManager::kTOFnSigmaKao:     
    fTree->Branch("kTOFnSigmaKao",&fgData[AliDielectronVarManager::kTOFnSigmaKao],"kTOFnSigmaKao/D"); break;
  case AliDielectronVarManager::kTOFnSigmaPro:     
    fTree->Branch("kTOFnSigmaPro",&fgData[AliDielectronVarManager::kTOFnSigmaPro],"kTOFnSigmaPro/D"); break;
  case AliDielectronVarManager::kKinkIndex0:       
    fTree->Branch("kKinkIndex0",&fgData[AliDielectronVarManager::kKinkIndex0],"kKinkIndex0/D"); break;
  case AliDielectronVarManager::kChi2NDF:  
    fTree->Branch("kChi2NDF",&fgData[AliDielectronVarManager::kChi2NDF],"kChi2NDF/D"); break;
  case AliDielectronVarManager::kDecayLength:      
    fTree->Branch("kDecayLength",&fgData[AliDielectronVarManager::kDecayLength],"kDecayLength/D"); break;
  case AliDielectronVarManager::kR:        
    fTree->Branch("kR",&fgData[AliDielectronVarManager::kR],"kR/D"); break;
  case AliDielectronVarManager::kOpeningAngle:     
    fTree->Branch("kOpeningAngle",&fgData[AliDielectronVarManager::kOpeningAngle],"kOpeningAngle/D"); break;
  case AliDielectronVarManager::kThetaHE:  
    fTree->Branch("kThetaHE",&fgData[AliDielectronVarManager::kThetaHE],"kThetaHE/D"); break;
  case AliDielectronVarManager::kPhiHE:    
    fTree->Branch("kPhiHE",&fgData[AliDielectronVarManager::kPhiHE],"kPhiHE/D"); break;
  case AliDielectronVarManager::kThetaCS:  
    fTree->Branch("kThetaCS",&fgData[AliDielectronVarManager::kThetaCS],"kThetaCS/D"); break;
  case AliDielectronVarManager::kPhiCS:    
    fTree->Branch("kPhiCS",&fgData[AliDielectronVarManager::kPhiCS],"kPhiCS/D"); break;
  case AliDielectronVarManager::kLegDist:  
    fTree->Branch("kLegDist",&fgData[AliDielectronVarManager::kLegDist],"kLegDist/D"); break;
  case AliDielectronVarManager::kLegDistXY:        
    fTree->Branch("kLegDistXY",&fgData[AliDielectronVarManager::kLegDistXY],"kLegDistXY/D"); break;
  case AliDielectronVarManager::kDeltaEta:         
    fTree->Branch("kDeltaEta",&fgData[AliDielectronVarManager::kDeltaEta],"kDeltaEta/D"); break;
  case AliDielectronVarManager::kDeltaPhi:         
    fTree->Branch("kDeltaPhi",&fgData[AliDielectronVarManager::kDeltaPhi],"kDeltaPhi/D"); break;
  case AliDielectronVarManager::kMerr:     
    fTree->Branch("kMerr",&fgData[AliDielectronVarManager::kMerr],"kMerr/D"); break;
  case AliDielectronVarManager::kDCA:      
    fTree->Branch("kDCA",&fgData[AliDielectronVarManager::kDCA],"kDCA/D"); break;
  case AliDielectronVarManager::kPairType:         
    fTree->Branch("kPairType",&fgData[AliDielectronVarManager::kPairType],"kPairType/D"); break;
  case AliDielectronVarManager::kPseudoProperTime:         
    fTree->Branch("kPseudoProperTime",&fgData[AliDielectronVarManager::kPseudoProperTime],"kPseudoProperTime/D"); break;
  case AliDielectronVarManager::kXvPrim:  
    fTree->Branch("kXvPrim",&fgData[AliDielectronVarManager::kXvPrim],"kXvPrim=kPairMax/D"); break;
  case AliDielectronVarManager::kYvPrim:   
    fTree->Branch("kYvPrim",&fgData[AliDielectronVarManager::kYvPrim],"kYvPrim/D"); break;
  case AliDielectronVarManager::kZvPrim:   
    fTree->Branch("kZvPrim",&fgData[AliDielectronVarManager::kZvPrim],"kZvPrim/D"); break;
  case AliDielectronVarManager::kXRes:     
    fTree->Branch("kXRes",&fgData[AliDielectronVarManager::kXRes],"kXRes/D"); break;
  case AliDielectronVarManager::kYRes:     
    fTree->Branch("kYRes",&fgData[AliDielectronVarManager::kYRes],"kYRes/D"); break;
  case AliDielectronVarManager::kZRes:     
    fTree->Branch("kZRes",&fgData[AliDielectronVarManager::kZRes],"kZRes/D"); break;
  case AliDielectronVarManager::kNTrk:     
    fTree->Branch("kNTrk",&fgData[AliDielectronVarManager::kNTrk],"kNTrk/D"); break;
  case AliDielectronVarManager::kTracks:   
    fTree->Branch("kTracks",&fgData[AliDielectronVarManager::kTracks],"kTracks/D"); break;
  case AliDielectronVarManager::kNacc:     
    fTree->Branch("kNacc",&fgData[AliDielectronVarManager::kNacc],"kNacc/D"); break;
  case AliDielectronVarManager::kNaccTrcklts:      
    fTree->Branch("kNaccTrcklts",&fgData[AliDielectronVarManager::kNaccTrcklts],"kNaccTrcklts/D"); break;
  case AliDielectronVarManager::kNch:      
    fTree->Branch("kNch",&fgData[AliDielectronVarManager::kNch],"kNch/D"); break;
  case AliDielectronVarManager::kCentrality:       
    fTree->Branch("kCentrality",&fgData[AliDielectronVarManager::kCentrality],"kCentrality/D"); break;
  case AliDielectronVarManager::kNevents:  
    fTree->Branch("kNevents",&fgData[AliDielectronVarManager::kNevents],"kNevents/D"); break;
  }
  fTreeList.Add(fTree);
}

//_____________________________________________________________________________
void AliDielectronHistosTaku::FillTree(Int_t nValues, const Double_t *values)
{

  for(int i=0;i<nValues;i++){
    fgData[i]=values[i];
  }
  
  //std::cout<<" --> "<<fgData[AliDielectronVarManager::kPIn]<<" "<<fgData[AliDielectronVarManager::kTPCsignal]<<" "<<fgData[AliDielectronVarManager::kTOFbeta]<<std::endl;
  //fTree->Fill();
}

//_____________________________________________________________________________
void AliDielectronHistosTaku::SetTreeList(THashList &list, Bool_t setOwner/*=kTRUE*/)
{
  //
  // set histogram classes and histograms to this instance. It will take onwnership!
  //
  ResetHistogramList();
  TString name(GetName());
  if (name == "AliDielectronHistosTaku") SetName(list.GetName());
  TIter next(&list);
  TObject *o;
  while ( (o=next()) ){
    fTreeList.Add(o);
  }
  if (setOwner){
    list.SetOwner(kFALSE);
    fTreeList.SetOwner(kTRUE);
  } else {
    fTreeList.SetOwner(kFALSE);
  }
}
