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
//   Julian Book   <Julian.Book@cern.ch> 
// 

#include <TH1.h>
#include <TH1F.h>
#include <TH2.h>
#include <TH3.h>
#include <THnBase.h>
#include <THn.h>
#include <THnSparse.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TProfile3D.h>
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
#include "AliDielectronVarManager.h"
#include "AliDielectronHistos.h"

ClassImp(AliDielectronHistos)


AliDielectronHistos::AliDielectronHistos() :
//   TCollection(),
  TNamed("AliDielectronHistos","Dielectron Histogram Container"),
  fHistoList(),
  fList(0x0),
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
  fList(0x0),
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
  fHistoList.Clear();
  if (fList) fList->Clear();
  delete fReservedWords;
}

//_____________________________________________________________________________
void AliDielectronHistos::UserProfile(const char* histClass,const char *name, const char* title,
				      UInt_t valTypeP,
				      Int_t nbinsX, Double_t xmin, Double_t xmax,
				      UInt_t valTypeX, Bool_t logBinX, TString option)
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
  UserProfile(histClass,name,title,valTypeP,binLimX,valTypeX,option);
}

//_____________________________________________________________________________
void AliDielectronHistos::UserProfile(const char* histClass,const char *name, const char* title,
				      UInt_t valTypeP,
				      Int_t nbinsX, Double_t xmin, Double_t xmax,
				      Int_t nbinsY, Double_t ymin, Double_t ymax,
				      UInt_t valTypeX, UInt_t valTypeY,
				      Bool_t logBinX, Bool_t logBinY, TString option)
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
  UserProfile(histClass,name,title,valTypeP,binLimX,binLimY,valTypeX,valTypeY,option);
}


//_____________________________________________________________________________
void AliDielectronHistos::UserProfile(const char* histClass,const char *name, const char* title,
				      UInt_t valTypeP,
				      Int_t nbinsX, Double_t xmin, Double_t xmax,
				      Int_t nbinsY, Double_t ymin, Double_t ymax,
				      Int_t nbinsZ, Double_t zmin, Double_t zmax,
				      UInt_t valTypeX, UInt_t valTypeY, UInt_t valTypeZ,
				      Bool_t logBinX, Bool_t logBinY, Bool_t logBinZ, TString option)
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

  UserProfile(histClass,name,title,valTypeP,binLimX,binLimY,binLimZ,valTypeX,valTypeY,valTypeZ,option);
}

//_____________________________________________________________________________
void AliDielectronHistos::UserProfile(const char* histClass,const char *name, const char* title,
				      UInt_t valTypeP,
				      const char* binning,
				      UInt_t valTypeX, TString option)
{
  //
  // Histogram creation 1D case with arbitraty binning
  //

  TVectorD *binLimX=AliDielectronHelper::MakeArbitraryBinning(binning);
  UserProfile(histClass,name,title,valTypeP,binLimX,valTypeX,option);
}

//_____________________________________________________________________________
void AliDielectronHistos::UserProfile(const char* histClass,const char *name, const char* title,
				      UInt_t valTypeP,
				      const TVectorD * const binsX,
				      UInt_t valTypeX/*=kNoAutoFill*/, TString option)
{
  //
  // Histogram creation 1D case with arbitraty binning X
  // the TVectorD is assumed to be surplus after the creation and will be deleted!!!
  //

  Bool_t isOk=kTRUE;
  isOk&=IsHistogramOk(histClass,name);
  isOk&=(binsX!=0x0);
  TH1 *hist=0x0;

  if (isOk){
    if(valTypeP==999)
      hist=new TH1F(name,title,binsX->GetNrows()-1,binsX->GetMatrixArray());
    else {
      TString opt=""; Double_t pmin=0., pmax=0.;
      if(!option.IsNull()) {
	TObjArray *arr=option.Tokenize(";");
	arr->SetOwner();
	opt=((TObjString*)arr->At(0))->GetString();
	if(arr->GetEntriesFast()>1) pmin=(((TObjString*)arr->At(1))->GetString()).Atof();
	if(arr->GetEntriesFast()>2) pmax=(((TObjString*)arr->At(2))->GetString()).Atof();
	delete arr;
      }
      hist=new TProfile(name,title,binsX->GetNrows()-1,binsX->GetMatrixArray());
      ((TProfile*)hist)->BuildOptions(pmin,pmax,opt.Data());
      //      printf(" name %s PROFILE options: pmin %.1f pmax %.1f err %s \n",name,((TProfile*)hist)->GetYmin(),((TProfile*)hist)->GetYmax(),((TProfile*)hist)->GetErrorOption() );
    }

    // store variales in axes
    UInt_t valType[20] = {0};
    valType[0]=valTypeX;     valType[1]=valTypeP;
    StoreVariables(hist, valType);

    // adapt the name and title of the histogram in case they are empty
    AdaptNameTitle(hist, histClass);

    Bool_t isReserved=fReservedWords->Contains(histClass);
    if(valTypeX==kNoAutoFill) hist->SetUniqueID(valTypeX);
    if (isReserved)
      UserHistogramReservedWords(histClass, hist, 999);
    else
      UserHistogram(histClass, hist, 999);
  }
  
  delete binsX;
}

//_____________________________________________________________________________
void AliDielectronHistos::UserProfile(const char* histClass,const char *name, const char* title,
				      UInt_t valTypeP,
				      const TVectorD * const binsX, const TVectorD * const binsY,
				      UInt_t valTypeX/*=kNoAutoFill*/, UInt_t valTypeY/*=0*/, TString option)
{
  //
  // Histogram creation 2D case with arbitraty binning X
  // the TVectorD is assumed to be surplus after the creation and will be deleted!!!
  //

  Bool_t isOk=kTRUE;
  isOk&=IsHistogramOk(histClass,name);
  isOk&=(binsX!=0x0);
  isOk&=(binsY!=0x0);
  TH1 *hist=0x0;

  if (isOk){
    if(valTypeP==999) {
      hist=new TH2F(name,title,
                    binsX->GetNrows()-1,binsX->GetMatrixArray(),
                    binsY->GetNrows()-1,binsY->GetMatrixArray()); 
    }
    else  {
      TString opt=""; Double_t pmin=0., pmax=0.;
      if(!option.IsNull()) {
	TObjArray *arr=option.Tokenize(";");
	arr->SetOwner();
	opt=((TObjString*)arr->At(0))->GetString();
	if(arr->GetEntriesFast()>1) pmin=(((TObjString*)arr->At(1))->GetString()).Atof();
	if(arr->GetEntriesFast()>2) pmax=(((TObjString*)arr->At(2))->GetString()).Atof();
	delete arr;
      }
      hist=new TProfile2D(name,title,
                          binsX->GetNrows()-1,binsX->GetMatrixArray(),
                          binsY->GetNrows()-1,binsY->GetMatrixArray());
      ((TProfile2D*)hist)->BuildOptions(pmin,pmax,opt.Data());
      //      printf(" name %s PROFILE options: pmin %.1f pmax %.1f err %s \n",name,((TProfile*)hist)->GetYmin(),((TProfile*)hist)->GetYmax(),((TProfile*)hist)->GetErrorOption() );
    }

    // store variales in axes
    UInt_t valType[20] = {0};
    valType[0]=valTypeX;     valType[1]=valTypeY;     valType[2]=valTypeP;
    StoreVariables(hist, valType);

    // adapt the name and title of the histogram in case they are empty
    AdaptNameTitle(hist, histClass);

    Bool_t isReserved=fReservedWords->Contains(histClass);
    if(valTypeX==kNoAutoFill) hist->SetUniqueID(valTypeX);
    if (isReserved)
      UserHistogramReservedWords(histClass, hist, 999);
    else
      UserHistogram(histClass, hist, 999);
  }
  
  delete binsX;
  delete binsY;
  
}

//_____________________________________________________________________________
void AliDielectronHistos::UserProfile(const char* histClass,const char *name, const char* title,
				      UInt_t valTypeP,
				      const TVectorD * const binsX, const TVectorD * const binsY, const TVectorD * const binsZ,
				      UInt_t valTypeX/*=kNoAutoFill*/, UInt_t valTypeY/*=0*/, UInt_t valTypeZ/*=0*/, TString option)
{
  //
  // Histogram creation 3D case with arbitraty binning X
  // the TVectorD is assumed to be surplus after the creation and will be deleted!!!
  //

  Bool_t isOk=kTRUE;
  isOk&=IsHistogramOk(histClass,name);
  isOk&=(binsX!=0x0);
  isOk&=(binsY!=0x0);
  isOk&=(binsZ!=0x0);
  TH1 *hist=0x0;

  if (isOk) {
    if(valTypeP==999) {
      hist=new TH3F(name,title,
		    binsX->GetNrows()-1,binsX->GetMatrixArray(),
		    binsY->GetNrows()-1,binsY->GetMatrixArray(),
		    binsZ->GetNrows()-1,binsZ->GetMatrixArray());
    }
    else {
      TString opt=""; Double_t pmin=0., pmax=0.;
      if(!option.IsNull()) {
	TObjArray *arr=option.Tokenize(";");
	arr->SetOwner();
	opt=((TObjString*)arr->At(0))->GetString();
	if(arr->GetEntriesFast()>1) pmin=(((TObjString*)arr->At(1))->GetString()).Atof();
	if(arr->GetEntriesFast()>2) pmax=(((TObjString*)arr->At(2))->GetString()).Atof();
	delete arr;
      }
      hist=new TProfile3D(name,title,
			  binsX->GetNrows()-1,binsX->GetMatrixArray(),
			  binsY->GetNrows()-1,binsY->GetMatrixArray(),
			  binsZ->GetNrows()-1,binsZ->GetMatrixArray());
      ((TProfile3D*)hist)->BuildOptions(pmin,pmax,opt.Data());
      //      printf(" name %s PROFILE options: pmin %.1f pmax %.1f err %s \n",name,((TProfile*)hist)->GetYmin(),((TProfile*)hist)->GetYmax(),((TProfile*)hist)->GetErrorOption() );
    }

    // store variales in axes
    UInt_t valType[20] = {0};
    valType[0]=valTypeX;     valType[1]=valTypeY;     valType[2]=valTypeZ;     valType[3]=valTypeP;
    StoreVariables(hist, valType);

    // adapt the name and title of the histogram in case they are empty
    AdaptNameTitle(hist, histClass);

    Bool_t isReserved=fReservedWords->Contains(histClass);
    if(valTypeX==kNoAutoFill) hist->SetUniqueID(valTypeX);
    if (isReserved)
      UserHistogramReservedWords(histClass, hist, 999);
    else
      UserHistogram(histClass, hist, 999);
  }
  
  delete binsX;
  delete binsY;
  delete binsZ;
}

//_____________________________________________________________________________
void AliDielectronHistos::UserHistogram(const char* histClass, Int_t ndim, Int_t *bins, Double_t *mins, Double_t *maxs, UInt_t *vars)
{
  //
  // Histogram creation 4-n dimension only with linear binning
  //

  Bool_t isOk=kTRUE;
  isOk&=(ndim<21 && ndim>3);
  if(!isOk) { Warning("UserHistogram","Array sizes should be between 3 and 20. Not adding Histogram to '%s'.", histClass); return; }

  // set automatic histo name
  TString name;
  for(Int_t iv=0; iv < ndim; iv++)
    name+=Form("%s_",AliDielectronVarManager::GetValueName(vars[iv]));
  name.Resize(name.Length()-1);

  isOk&=IsHistogramOk(histClass,name);

  THnD *hist;
  if (isOk) {
    hist=new THnD(name.Data(),"", ndim, bins, mins, maxs);

    // store variales in axes
    StoreVariables(hist, vars);

    Bool_t isReserved=fReservedWords->Contains(histClass);
    if (isReserved)
      UserHistogramReservedWords(histClass, hist, 999);
    else
      UserHistogram(histClass, hist, 999);

  }
}

//_____________________________________________________________________________
void AliDielectronHistos::UserHistogram(const char* histClass, Int_t ndim, TObjArray *limits, UInt_t *vars)
{
  //
  // Histogram creation n>3 dimension only with non-linear binning
  //

  Bool_t isOk=kTRUE;
  isOk&=(ndim<21 && ndim>3);
  if(!isOk) { Warning("UserHistogram","Array sizes should be between 3 and 20. Not adding Histogram to '%s'.", histClass); return; }
  isOk&=(ndim==limits->GetEntriesFast());
  if(!isOk) return;

  // set automatic histo name
  TString name;
  for(Int_t iv=0; iv < ndim; iv++)
    name+=Form("%s_",AliDielectronVarManager::GetValueName(vars[iv]));
  name.Resize(name.Length()-1);

  isOk&=IsHistogramOk(histClass,name);

  THnD *hist;
  Int_t bins[ndim];
  if (isOk) {
    // get number of bins
    for(Int_t idim=0 ;idim<ndim; idim++) {
      TVectorD *vec = (TVectorD*) limits->At(idim);
      bins[idim]=vec->GetNrows()-1;
    }

    hist=new THnD(name.Data(),"", ndim, bins, 0x0, 0x0);

    // set binning
    for(Int_t idim=0 ;idim<ndim; idim++) {
      TVectorD *vec = (TVectorD*) limits->At(idim);
      hist->SetBinEdges(idim,vec->GetMatrixArray());
    }

    // store variales in axes
    StoreVariables(hist, vars);

    Bool_t isReserved=fReservedWords->Contains(histClass);
    if (isReserved)
      UserHistogramReservedWords(histClass, hist, 999);
    else
      UserHistogram(histClass, hist, 999);

  }
}

//_____________________________________________________________________________
void AliDielectronHistos::UserSparse(const char* histClass, Int_t ndim, Int_t *bins, Double_t *mins, Double_t *maxs, UInt_t *vars)
{
  //
  // THnSparse creation with linear binning
  //

  Bool_t isOk=kTRUE;

  // set automatic histo name
  TString name;
  for(Int_t iv=0; iv < ndim; iv++)
    name+=Form("%s_",AliDielectronVarManager::GetValueName(vars[iv]));
  name.Resize(name.Length()-1);

  isOk&=IsHistogramOk(histClass,name);

  THnSparseD *hist;
  if (isOk) {
    hist=new THnSparseD(name.Data(),"", ndim, bins, mins, maxs);

    // store variales in axes
    StoreVariables(hist, vars);

    Bool_t isReserved=fReservedWords->Contains(histClass);
    if (isReserved)
      UserHistogramReservedWords(histClass, hist, 999);
    else
      UserHistogram(histClass, hist, 999);

  }
}

//_____________________________________________________________________________
void AliDielectronHistos::UserSparse(const char* histClass, Int_t ndim, TObjArray *limits, UInt_t *vars)
{
  //
  // THnSparse creation with non-linear binning
  //

  Bool_t isOk=kTRUE;
  isOk&=(ndim==limits->GetEntriesFast());
  if(!isOk) return;

  // set automatic histo name
  TString name;
  for(Int_t iv=0; iv < ndim; iv++)
    name+=Form("%s_",AliDielectronVarManager::GetValueName(vars[iv]));
  name.Resize(name.Length()-1);

  isOk&=IsHistogramOk(histClass,name);

  THnSparseD *hist;
  Int_t bins[ndim];
  if (isOk) {
    // get number of bins
    for(Int_t idim=0 ;idim<ndim; idim++) {
      TVectorD *vec = (TVectorD*) limits->At(idim);
      bins[idim]=vec->GetNrows()-1;
    }

    hist=new THnSparseD(name.Data(),"", ndim, bins, 0x0, 0x0);

    // set binning
    for(Int_t idim=0 ;idim<ndim; idim++) {
      TVectorD *vec = (TVectorD*) limits->At(idim);
      hist->SetBinEdges(idim,vec->GetMatrixArray());
    }

    // store variales in axes
    StoreVariables(hist, vars);

    Bool_t isReserved=fReservedWords->Contains(histClass);
    if (isReserved)
      UserHistogramReservedWords(histClass, hist, 999);
    else
      UserHistogram(histClass, hist, 999);

  }
}

//_____________________________________________________________________________
void AliDielectronHistos::UserHistogram(const char* histClass, TObject* hist, UInt_t valTypes)
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
  //  hist->SetDirectory(0);

  // store variables axis
  UInt_t valType[20] = {0};
  // incase valTypes is given old way of extracting variables
  if(valTypes!=999) {
    valType[0]=valTypes%1000;          //last three digits
    valType[1]=valTypes/1000%1000;     //second last three digits
    valType[2]=valTypes/1000000%1000;  //third last three digits
    hist->SetUniqueID(valTypes);
  }
  else {
    // extract variables from axis
    FillVarArray(hist, valType);
  }
  StoreVariables(hist, valType);

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
    Warning("FillClass","Cannot fill class '%s' its not defined. nValues %d",histClass,nValues);
    return;
  }

  TIter nextHist(classTable);
  TObject *obj=0;
  while ( (obj=(TObject*)nextHist()) )  FillValues(obj, values);

  return;
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
void AliDielectronHistos::UserHistogramReservedWords(const char* histClass, const TObject *hist, UInt_t valTypes)
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
      TObject *h=hist->Clone();
      // Tobject has no function SetDirectory, didn't we need this???
      //      h->SetDirectory(0);
      ((TH1*)h)->SetTitle(Form("%s %s",title.Data(),l->GetName()));

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
TObject* AliDielectronHistos::GetHist(const char* histClass, const char* name) const
{
  //
  // return object 'name' in 'histClass'
  //
  THashList *classTable=(THashList*)fHistoList.FindObject(histClass);
  if (!classTable) return 0x0;
  return classTable->FindObject(name);
}

//_____________________________________________________________________________
TH1* AliDielectronHistos::GetHistogram(const char* histClass, const char* name) const
{
  //
  // return histogram 'name' in 'histClass'
  //
  return ((TH1*) GetHist(histClass, name));
}

//_____________________________________________________________________________
TObject* AliDielectronHistos::GetHist(const char* cutClass, const char* histClass, const char* name) const
{
  //
  // return object from list of list of histograms
  // this function is thought for retrieving histograms if a list of AliDielectronHistos is set
  //
  
  if (!fList) return 0x0;
  THashList *h=dynamic_cast<THashList*>(fList->FindObject(cutClass));
  if (!h)return 0x0;
  THashList *classTable=dynamic_cast<THashList*>(h->FindObject(histClass));
  if (!classTable) return 0x0;
  return classTable->FindObject(name);
}

//_____________________________________________________________________________
TH1* AliDielectronHistos::GetHistogram(const char* cutClass, const char* histClass, const char* name) const
{
  //
  // return histogram from list of list of histograms
  // this function is thought for retrieving histograms if a list of AliDielectronHistos is set
  //
  return ((TH1*) GetHist(cutClass, histClass, name));
}

//_____________________________________________________________________________
void AliDielectronHistos::Draw(const Option_t* option)
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
void AliDielectronHistos::SetHistogramList(THashList &list, Bool_t setOwner/*=kTRUE*/)
{
  //
  // set histogram classes and histograms to this instance. It will take onwnership!
  //
  ResetHistogramList();
  TString name(GetName());
  if (name == "AliDielectronHistos") SetName(list.GetName());
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
Bool_t AliDielectronHistos::SetCutClass(const char* cutClass)
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
  if (GetHist(histClass,name)){
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
void AliDielectronHistos::StoreVariables(TObject *obj, UInt_t valType[20])
{
  //
  //
  //
  if (!obj) return;
  if      (obj->InheritsFrom(TH1::Class()))         StoreVariables(static_cast<TH1*>(obj), valType);
  else if (obj->InheritsFrom(THnBase::Class()))         StoreVariables(static_cast<THnBase*>(obj), valType);
  //  else if (obj->InheritsFrom(THnSparse::Class()))   StoreVariables(static_cast<THnSparse*>(obj), valType);

  return;

}


//_____________________________________________________________________________
void AliDielectronHistos::StoreVariables(TH1 *obj, UInt_t valType[20])
{
  //
  // store variables in the axis (special for TProfile3D)
  //

  Int_t dim   = obj->GetDimension();

  // dimension correction for profiles
  if(obj->IsA() == TProfile::Class() || obj->IsA() == TProfile2D::Class() || obj->IsA() == TProfile3D::Class()) {
    dim++;
  }

  switch( dim ) {
  case 4:
    obj->SetUniqueID(valType[3]); // Tprofile3D variable
  case 3:
    obj->GetZaxis()->SetUniqueID(valType[2]);
  case 2:
    obj->GetYaxis()->SetUniqueID(valType[1]);
  case 1:
    obj->GetXaxis()->SetUniqueID(valType[0]);
  }

  return;
}

//_____________________________________________________________________________
void AliDielectronHistos::StoreVariables(THnBase *obj, UInt_t valType[20])
{
  //
  // store variables in the axis
  //

  Int_t dim = obj->GetNdimensions();

  for(Int_t it=0; it<dim; it++) {
    obj->GetAxis(it)->SetUniqueID(valType[it]);
    obj->GetAxis(it)->SetTitle(Form("%s %s", AliDielectronVarManager::GetValueLabel(valType[it]), AliDielectronVarManager::GetValueUnit(valType[it])));
  }
  obj->Sumw2();
  return;
}

//_____________________________________________________________________________
void AliDielectronHistos::FillValues(TObject *obj, const Double_t *values)
{
  //
  //
  //
  if (!obj) return;
  if      (obj->InheritsFrom(TH1::Class()))   FillValues(static_cast<TH1*>(obj), values);
  else if (obj->InheritsFrom(THn::Class()))   FillValues(static_cast<THn*>(obj), values);
  else if (obj->InheritsFrom(THnSparse::Class()))   FillValues(static_cast<THnSparse*>(obj), values);

  return;

}

//_____________________________________________________________________________
void AliDielectronHistos::FillValues(TH1 *obj, const Double_t *values)
{
  //
  // fill values for TH1 inherted classes
  //

  Int_t dim   = obj->GetDimension();
  Bool_t bprf = kFALSE;
  //  UInt_t nValues = (UInt_t) AliDielectronVarManager::kNMaxValues;

  UInt_t valueTypes=obj->GetUniqueID();
  if (valueTypes==(UInt_t)kNoAutoFill) return;

  if(obj->IsA() == TProfile::Class() || obj->IsA() == TProfile2D::Class() || obj->IsA() == TProfile3D::Class())
    bprf=kTRUE;
  
  UInt_t value1=obj->GetXaxis()->GetUniqueID();
  UInt_t value2=obj->GetYaxis()->GetUniqueID();
  UInt_t value3=obj->GetZaxis()->GetUniqueID();
  UInt_t value4=obj->GetUniqueID();            // get profile var stored in the unique ID

  //  if (value1>=(UInt_t)nValues||value2>=(UInt_t)nValues||value3>=(UInt_t)nValues||(value4>=(UInt_t)nValues && value4!=999)) {
  //  Warning("FillClass","One of the values is out of range. Not filling Histogram '%s/%s'.", histClass, obj->GetName());
  //  return;
  // }

  switch ( dim ) {
  case 1:
    if(!bprf) obj->Fill(values[value1]);                    // histograms
    else ((TProfile*)obj)->Fill(values[value1],values[value2]);   // profiles
    break;
  case 2:
    if(!bprf) ((TH2*)obj)->Fill(values[value1],values[value2]);
    else ((TProfile2D*)obj)->Fill(values[value1],values[value2],values[value3]);
    break;
  case 3:
    if(!bprf) ((TH3*)obj)->Fill(values[value1],values[value2],values[value3]);
    else ((TProfile3D*)obj)->Fill(values[value1],values[value2],values[value3],values[value4]);
    break;
  }
  
  return;
}

//_____________________________________________________________________________
void AliDielectronHistos::FillValues(THnBase *obj, const Double_t *values)
{
  //
  // fill values for THn inherted classes
  //

  const Int_t dim   = obj->GetNdimensions();

  UInt_t valueTypes=obj->GetUniqueID();
  if (valueTypes==(UInt_t)kNoAutoFill) return;

  Double_t fill[dim];
  for(Int_t it=0; it<dim; it++)   fill[it] = values[obj->GetAxis(it)->GetUniqueID()];
  obj->Fill(fill);

  return;
}

//_____________________________________________________________________________
void AliDielectronHistos::FillVarArray(TObject *obj, UInt_t *valType)
{
  //
  // extract variables stored in the axis (special for TProfile3D)
  //


  if (!obj) return;
  //  printf(" fillvararray %s \n",obj->GetName());

  if (obj->InheritsFrom(TH1::Class())) {
    valType[0]=((TH1*)obj)->GetXaxis()->GetUniqueID();
    valType[1]=((TH1*)obj)->GetYaxis()->GetUniqueID();
    valType[2]=((TH1*)obj)->GetZaxis()->GetUniqueID();
    valType[3]=((TH1*)obj)->GetUniqueID();  // tprofile var stored in unique ID
  }
  else if (obj->InheritsFrom(THnBase::Class())) {
    for(Int_t it=0; it<((THn*)obj)->GetNdimensions(); it++)
      valType[it]=((THn*)obj)->GetAxis(it)->GetUniqueID();
  }
  return;
}

//_____________________________________________________________________________
void AliDielectronHistos::AdaptNameTitle(TH1 *hist, const char* histClass) {

  //
  // adapt name and title of the histogram
  //

  Int_t dim            = hist->GetDimension();
  TString currentName  = hist->GetName();
  TString currentTitle = hist->GetTitle();


  Bool_t bname  = (currentName.IsNull());
  Bool_t btitle = (currentTitle.IsNull());
  Bool_t bprf   = kFALSE;
  if(hist->IsA() == TProfile::Class() || hist->IsA() == TProfile2D::Class() || hist->IsA() == TProfile3D::Class())
    bprf=kTRUE;

  // tprofile options
  Double_t pmin=0., pmax=0.;
  TString option = "", calcrange="";
  Bool_t bStdOpt=kTRUE;
  if(bprf) {
    switch( dim ) {
    case 3:
      option = ((TProfile3D*)hist)->GetErrorOption();
      pmin   = ((TProfile3D*)hist)->GetTmin();
      pmax   = ((TProfile3D*)hist)->GetTmax();
      break;
    case 2:
      option = ((TProfile2D*)hist)->GetErrorOption();
      pmin   = ((TProfile2D*)hist)->GetZmin();
      pmax   = ((TProfile2D*)hist)->GetZmax();
      break;
    case 1:
      option = ((TProfile*)hist)->GetErrorOption();
      pmin   = ((TProfile*)hist)->GetYmin();
      pmax   = ((TProfile*)hist)->GetYmax();
      break;
    }
    if(option.Contains("s",TString::kIgnoreCase)) bStdOpt=kFALSE;
    if(pmin!=pmax) calcrange=Form("#cbar_{%+.*f}^{%+.*f}",GetPrecision(pmin),pmin,GetPrecision(pmax),pmax);
  }

  UInt_t varx = hist->GetXaxis()->GetUniqueID();
  UInt_t vary = hist->GetYaxis()->GetUniqueID();
  UInt_t varz = hist->GetZaxis()->GetUniqueID();
  UInt_t varp = hist->GetUniqueID();

  // store titles in the axis
  if(btitle) {
    switch( dim ) {
    case 3:
      hist->GetXaxis()->SetTitle(Form("%s %s",
				      AliDielectronVarManager::GetValueLabel(varx), 
				      AliDielectronVarManager::GetValueUnit(varx))
				 );
      hist->GetYaxis()->SetTitle(Form("%s %s",
				      AliDielectronVarManager::GetValueLabel(vary), 
				      AliDielectronVarManager::GetValueUnit(vary))
				 );
      hist->GetZaxis()->SetTitle(Form("%s %s",
				      AliDielectronVarManager::GetValueLabel(varz), 
				      AliDielectronVarManager::GetValueUnit(varz))
				 );
      if(bprf)
	hist->SetTitle(Form("%s  %s%s%s%s %s",
			    hist->GetTitle(),
			    (bStdOpt ? "#LT" : "RMS("),
			    AliDielectronVarManager::GetValueLabel(varp), 
			    (bStdOpt ? "#GT" : ")"),
			    calcrange.Data(),
			    AliDielectronVarManager::GetValueUnit(varp))
		       );
      break;
    case 2:
      hist->GetXaxis()->SetTitle(Form("%s %s",
				      AliDielectronVarManager::GetValueLabel(varx), 
				      AliDielectronVarManager::GetValueUnit(varx))
				 );
      hist->GetYaxis()->SetTitle(Form("%s %s",
				      AliDielectronVarManager::GetValueLabel(vary), 
				      AliDielectronVarManager::GetValueUnit(vary))
				 );
      hist->GetZaxis()->SetTitle(Form("#%ss",histClass));
      if(bprf)
	hist->GetZaxis()->SetTitle(Form("%s%s%s%s %s",
					(bStdOpt ? "#LT" : "RMS("),
					AliDielectronVarManager::GetValueLabel(varz), 
					(bStdOpt ? "#GT" : ")"),
					calcrange.Data(),
					AliDielectronVarManager::GetValueUnit(varz))
				   );
      break;
    case 1:
      hist->GetXaxis()->SetTitle(Form("%s %s",
				      AliDielectronVarManager::GetValueLabel(varx), 
				      AliDielectronVarManager::GetValueUnit(varx))
				 );
      hist->GetYaxis()->SetTitle(Form("#%ss",histClass));
      if(bprf)
	hist->GetYaxis()->SetTitle(Form("%s%s%s%s %s",
					(bStdOpt ? "#LT" : "RMS("),
					AliDielectronVarManager::GetValueLabel(vary), 
					(bStdOpt ? "#GT" : ")"),
					calcrange.Data(),
					AliDielectronVarManager::GetValueUnit(vary))
				   );
      break;
    }

    // create an unique name
    if(bname)
      switch(dim) {
      case 3:
	currentName+=Form("%s_",AliDielectronVarManager::GetValueName(varx));
	currentName+=Form("%s_",AliDielectronVarManager::GetValueName(vary));
	currentName+=Form("%s",AliDielectronVarManager::GetValueName(varz));
	if(bprf) currentName+=Form("-%s%s",AliDielectronVarManager::GetValueName(varp),(bStdOpt ? "avg" : "rms"));
	break;
      case 2:
	currentName+=Form("%s_",AliDielectronVarManager::GetValueName(varx));
	currentName+=Form("%s",AliDielectronVarManager::GetValueName(vary));
	if(bprf) currentName+=Form("-%s%s",AliDielectronVarManager::GetValueName(varz),(bStdOpt ? "avg" : "rms"));
	break;
      case 1:
	currentName+=Form("%s",AliDielectronVarManager::GetValueName(varx));
	if(bprf) currentName+=Form("-%s%s",AliDielectronVarManager::GetValueName(vary),(bStdOpt ? "avg" : "rms"));
	break;
      }
    // to differentiate btw. leg and pair histos
    if(!strcmp(histClass,"Pair")) currentName.Prepend("p");
    hist->SetName(currentName.Data());
  }

}

Int_t AliDielectronHistos::GetPrecision(Double_t value) {

  //
  // computes the precision of a double
  // usefull for axis ranges etc
  //

  Bool_t bfnd     = kFALSE;
  Int_t precision = 0;

  while(!bfnd) {
    //    printf(" value %f precision %d bfnd %d \n",TMath::Abs(value*TMath::Power(10,precision)), precision, bfnd);
    bfnd = (TMath::Abs(value*TMath::Power(10,precision))  -  TMath::Floor(TMath::Abs(value*TMath::Power(10,precision))) != 0.0
	    ? kFALSE
	    : kTRUE);
    if(!bfnd) precision++;
  }

  //  printf("precision for %f found to be %d \n", value, precision);
  return precision;

}



