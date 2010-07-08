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

#include <TH1.h>
#include <TH1F.h>
#include <TH2.h>
#include <TH3.h>
#include <TCollection.h>
#include <THashList.h>
#include <TString.h>
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
// #include <TVectorD.h>

#include "AliDielectronHistos.h"


ClassImp(AliDielectronHistos)


AliDielectronHistos::AliDielectronHistos() :
//   TCollection(),
  TNamed("AliDielectronHistos","Dielectron Histogram Container"),
  fHistoList(),
  fReservedWords(new TString)
{
  //
  // Default constructor
  //
  fHistoList.SetOwner(kTRUE);
  fHistoList.SetName("Dielectron_Histos");
}

//_____________________________________________________________________________
AliDielectronHistos::AliDielectronHistos(const char* name, const char* title) :
//   TCollection(),
  TNamed(name, title),
  fHistoList(),
  fReservedWords(new TString)
{
  //
  // TNamed constructor
  //
  fHistoList.SetOwner(kTRUE);
  fHistoList.SetName(name);
}

//_____________________________________________________________________________
AliDielectronHistos::~AliDielectronHistos()
{
  //
  // Destructor
  //
  fHistoList.Delete();
  delete fReservedWords;
}

//_____________________________________________________________________________
void AliDielectronHistos::UserHistogram(const char* histClass,const char *name, const char* title,
                                        Int_t nbinsX, Double_t xmin, Double_t xmax,
                                        UInt_t valTypeX, Bool_t logBinX)
{
  //
  // Default histogram creation 1D case
  //
  if (!IsHistogramOk(histClass,name)) return;

  Double_t *binLimX=0x0;
  
  if (logBinX) {
    binLimX=MakeLogBinning(nbinsX, xmin, xmax);
  } else {
    binLimX=MakeLinBinning(nbinsX, xmin, xmax);
  }

  TH1* hist=new TH1F(name,title,nbinsX,binLimX);

  delete [] binLimX;
  
  Bool_t isReserved=fReservedWords->Contains(histClass);
  if (isReserved)
    UserHistogramReservedWords(histClass, hist, valTypeX);
  else
    UserHistogram(histClass, hist, valTypeX);
}

//_____________________________________________________________________________
void AliDielectronHistos::UserHistogram(const char* histClass,const char *name, const char* title,
                                        Int_t nbinsX, Double_t xmin, Double_t xmax,
                                        Int_t nbinsY, Double_t ymin, Double_t ymax,
                                        UInt_t valTypeX, UInt_t valTypeY,
                                        Bool_t logBinX, Bool_t logBinY)
{
  //
  // Default histogram creation 2D case
  //
  if (!IsHistogramOk(histClass,name)) return;

  Double_t *binLimX=0x0;
  Double_t *binLimY=0x0;
  
  if (logBinX) {
    binLimX=MakeLogBinning(nbinsX, xmin, xmax);
  } else {
    binLimX=MakeLinBinning(nbinsX, xmin, xmax);
  }
  if (logBinY) {
    binLimY=MakeLogBinning(nbinsY, ymin, ymax);
  } else {
    binLimY=MakeLinBinning(nbinsY, ymin, ymax);
  }
  
  TH1* hist=new TH2F(name,title,nbinsX,binLimX,nbinsY,binLimY);
  
  delete [] binLimX;
  delete [] binLimY;
  
  
  Bool_t isReserved=fReservedWords->Contains(histClass);
  if (isReserved)
    UserHistogramReservedWords(histClass, hist, valTypeX+100*valTypeY);
  else
    UserHistogram(histClass, hist, valTypeX+100*valTypeY);
}


//_____________________________________________________________________________
void AliDielectronHistos::UserHistogram(const char* histClass,const char *name, const char* title,
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

  Double_t *binLimX=0x0;
  Double_t *binLimY=0x0;
  Double_t *binLimZ=0x0;
  
  if (logBinX) {
    binLimX=MakeLogBinning(nbinsX, xmin, xmax);
  } else {
    binLimX=MakeLinBinning(nbinsX, xmin, xmax);
  }
  
  if (logBinY) {
    binLimY=MakeLogBinning(nbinsY, ymin, ymax);
  } else {
    binLimY=MakeLinBinning(nbinsY, ymin, ymax);
  }
  
  if (logBinZ) {
    binLimZ=MakeLogBinning(nbinsZ, zmin, zmax);
  } else {
    binLimZ=MakeLinBinning(nbinsZ, zmin, zmax);
  }
  
  TH1* hist=new TH3F(name,title,nbinsX,binLimX,nbinsY,binLimY,nbinsZ,binLimZ);
  
  delete [] binLimX;
  delete [] binLimY;
  delete [] binLimZ;
  
  Bool_t isReserved=fReservedWords->Contains(histClass);
  if (isReserved)
    UserHistogramReservedWords(histClass, hist, valTypeX+100*valTypeY+10000*valTypeZ);
  else
    UserHistogram(histClass, hist, valTypeX+100*valTypeY+10000*valTypeZ);
}

//_____________________________________________________________________________
void AliDielectronHistos::UserHistogram(const char* histClass, TH1* hist, UInt_t valTypes)
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
void AliDielectronHistos::AddClass(const char* histClass)
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
void AliDielectronHistos::Fill(const char* histClass, const char* name, Double_t xval)
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
void AliDielectronHistos::Fill(const char* histClass, const char* name, Double_t xval, Double_t yval)
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
void AliDielectronHistos::Fill(const char* histClass, const char* name, Double_t xval, Double_t yval, Double_t zval)
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
void AliDielectronHistos::FillClass(const char* histClass, Int_t nValues, const Double_t *values)
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
// void AliDielectronHistos::FillClass(const char* histClass, const TVectorD &vals)
// {
//   //
//   //
//   //
//   FillClass(histClass, vals.GetNrows(), vals.GetMatrixArray());
// }

//_____________________________________________________________________________
void AliDielectronHistos::UserHistogramReservedWords(const char* histClass, TH1 *hist, UInt_t valTypes)
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
void AliDielectronHistos::DumpToFile(const char* file)
{
  //
  // Dump the histogram list to a newly created root file
  //
  TFile f(file,"recreate");
  fHistoList.Write(fHistoList.GetName(),TObject::kSingleKey);
  f.Close();
}

//_____________________________________________________________________________
TH1* AliDielectronHistos::GetHistogram(const char* histClass, const char* name) const
{
  //
  // return histogram 'name' in 'histClass'
  //
  THashList *classTable=(THashList*)fHistoList.FindObject(histClass);
  if (!classTable) return 0x0;
  return (TH1*)classTable->FindObject(name);
}

//_____________________________________________________________________________
void AliDielectronHistos::Draw(const Option_t* option)
{
  //
  // Draw histograms
  //

  TString drawStr(option);
  drawStr.ToLower();
  //options
//   Bool_t same=drawOpt.Contains("same"); //FIXME not yet implemented

  TCanvas *c=0x0;
  if (gVirtualPS) {
    if (!gPad){
      Error("Draw","When writing to a file you have to create a canvas before opening the file!!!");
      return;
    }
    c=gPad->GetCanvas();
//     c=new TCanvas;
  }
  
  TIter nextClass(&fHistoList);
  THashList *classTable=0;
  Bool_t first=kTRUE;
  while ( (classTable=(THashList*)nextClass()) ){
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
      if (first){
        first=kFALSE;
      } else {
        c->Clear();
      }
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
    if (gVirtualPS) c->Update();
  }
//   if (gVirtualPS) delete c;
}

//_____________________________________________________________________________
void AliDielectronHistos::Print(const Option_t* option) const
{
  //
  // Print classes and histograms
  //
  TString optString(option);

  if (optString.IsNull()) PrintStructure();



}

//_____________________________________________________________________________
void AliDielectronHistos::PrintStructure() const
{
  //
  // Print classes and histograms in the class to stdout
  //
  TIter nextClass(&fHistoList);
  THashList *classTable=0;
  while ( (classTable=(THashList*)nextClass()) ){
    TIter nextHist(classTable);
    TObject *o=0;
    printf("+ %s\n",classTable->GetName());
    while ( (o=nextHist()) )
      printf("| ->%s\n",o->GetName());
  }
  
}

//_____________________________________________________________________________
void AliDielectronHistos::SetHistogramList(THashList &list)
{
  //
  // set histogram classes and histograms to this instance. It will take onwnership!
  //
  ResetHistogramList();
  SetName(list.GetName());
  TIter next(&list);
  TObject *o;
  while ( (o=next()) ){
    fHistoList.Add(o);
  }
  list.SetOwner(kFALSE);
  fHistoList.SetOwner(kTRUE);
}

//_____________________________________________________________________________
Bool_t AliDielectronHistos::IsHistogramOk(const char* histClass, const char* name)
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
// TIterator* AliDielectronHistos::MakeIterator(Bool_t dir) const
// {
//   //
//   //
//   //
//   return new TListIter(&fHistoList, dir);
// }

//_____________________________________________________________________________
void AliDielectronHistos::ReadFromFile(const char* file)
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
void AliDielectronHistos::DrawSame(const char* histName, const Option_t *opt)
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
void AliDielectronHistos::SetReservedWords(const char* words)
{
  //
  // set reserved words
  //
  
  (*fReservedWords)=words;
}

//_____________________________________________________________________________
Double_t* AliDielectronHistos::MakeLogBinning(Int_t nbinsX, Double_t xmin, Double_t xmax) const
{
  //
  // Make logarithmic binning
  // the user has to delete the array afterwards!!!
  //

  //check limits
  if (xmin<1e-20 || xmax<1e-20){
    Error("MakeLogBinning","For Log binning xmin and xmax must be > 1e-20. Using linear binning instead!");
    return MakeLinBinning(nbinsX, xmin, xmax);
  }
  Double_t *binLim=new Double_t[nbinsX+1];
  Double_t first=xmin;
  Double_t last=xmax;
  Double_t expMax=TMath::Log(last/first);
  for (Int_t i=0; i<nbinsX+1; ++i){
    binLim[i]=first*TMath::Exp(expMax/nbinsX*(Double_t)i);
  }
  return binLim;
}

//_____________________________________________________________________________
Double_t* AliDielectronHistos::MakeLinBinning(Int_t nbinsX, Double_t xmin, Double_t xmax) const
{
  //
  // Make logarithmic binning
  // the user has to delete the array afterwards!!!
  //
  Double_t *binLim=new Double_t[nbinsX+1];
  Double_t first=xmin;
  Double_t last=xmax;
  Double_t binWidth=(last-first)/nbinsX;
  for (Int_t i=0; i<nbinsX+1; ++i){
    binLim[i]=first+binWidth*(Double_t)i;
  }
  return binLim;
}

