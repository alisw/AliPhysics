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

///////////////////////////////////////////////////////////////////////////
//       Dielectron Correction framework draw helper                     //
//                                                                       //
/*









*/
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include <TSeqCollection.h>
#include <TObjArray.h>
#include <TKey.h>
#include <TList.h>
#include <TObject.h>
#include <TFile.h>
#include <TString.h>
#include <TObjString.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TPad.h>
#include <TCanvas.h>

#include <AliLog.h>

#include "AliDielectronCFdraw.h"

ClassImp(AliDielectronCFdraw)

AliDielectronCFdraw::AliDielectronCFdraw() :
  TNamed(),
  fCfContainer(0x0)
{
  //
  // Ctor
  //
}

//________________________________________________________________
AliDielectronCFdraw::AliDielectronCFdraw(const char*name, const char* title) :
  TNamed(name,title),
  fCfContainer(0x0)
{
  //
  // Named Ctor
  //
  
}

//________________________________________________________________
void AliDielectronCFdraw::SetCFContainers(const TSeqCollection *arr)
{
  //
  // Merge CF Container out of several containers
  //

  TIter next(arr);
  TObject *o=0x0;

  Int_t nstep=0;
  while ( (o=next()) ){
    AliCFContainer *cf=dynamic_cast<AliCFContainer*>(o);
    if (!o) continue;
    nstep+=cf->GetNStep();
  }
  Int_t nbins[1]={1};
  fCfContainer=new AliCFContainer(GetName(), GetTitle(), nstep, 1, nbins);

  //delete unneeded steps
  for (Int_t istep=0; istep<nstep; ++istep) delete fCfContainer->GetGrid(istep);

  //add step to the new container
  Int_t istep=0;
  for (Int_t icf=0; icf<arr->GetEntries(); ++icf){
    AliCFContainer *cf=dynamic_cast<AliCFContainer*>(arr->At(icf));
    if (!cf) continue;
    for (Int_t istepCurr=0; istepCurr<cf->GetNStep(); ++istepCurr){
      fCfContainer->SetGrid(istep, cf->GetGrid(istepCurr));
      fCfContainer->SetStepTitle(istep,Form("%s: %s",cf->GetName(),cf->GetStepTitle(istepCurr)));
      ++istep;
    }
  }
}

//________________________________________________________________
void AliDielectronCFdraw::SetCFContainers(const char* filename)
{
  //
  // get CF containers from file
  //

  TFile f(filename);
  TList *l=f.GetListOfKeys();
  Int_t entries=l->GetEntries();
  if (entries==0) return;
  
  TKey *k=(TKey*)l->At(0);
  TSeqCollection *arr=dynamic_cast<TSeqCollection*>(k->ReadObj());
  if (!arr) return;
  
  SetCFContainers(arr);
}

//________________________________________________________________
void AliDielectronCFdraw::SetRangeUser(Int_t ivar, Double_t min, Double_t max, const char* slices)
{
  //
  // Set range of cut steps defined in slices
  // Steps may be separated by one the the characteres ,;:
  //
  TObjArray *arr=TString(slices).Tokenize(",:;");

  if (arr->GetEntriesFast()==0){
    // all slices in case of 0 entries
    for (Int_t istep=0; istep<fCfContainer->GetNStep(); ++istep){
      fCfContainer->SetRangeUser(ivar,min,max,istep);
    }
  } else {
    TIter next(arr);
    TObjString *ostr=0x0;
    while ( (ostr=static_cast<TObjString*>(next())) ) {
      Int_t istep=ostr->GetString().Atoi();
      fCfContainer->SetRangeUser(ivar,min,max,istep);
    }
  }
  delete arr;
}

//________________________________________________________________
void AliDielectronCFdraw::UnsetRangeUser(Int_t ivar, const char* slices)
{
  //
  // Unset range of cut steps defined in slices
  // Steps may be separated by one the the characteres ,;:
  //
  TObjArray *arr=TString(slices).Tokenize(",:;");
  
  if (arr->GetEntriesFast()==0){
    // all slices in case of 0 entries
    for (Int_t istep=0; istep<fCfContainer->GetNStep(); ++istep){
      fCfContainer->SetRangeUser(ivar,fCfContainer->GetAxis(ivar,istep)->GetXmin(),
                                 fCfContainer->GetAxis(ivar,istep)->GetXmax(),istep);
    }
  } else {
    TIter next(arr);
    TObjString *ostr=0x0;
    while ( (ostr=static_cast<TObjString*>(next())) ) {
      Int_t istep=ostr->GetString().Atoi();
      fCfContainer->SetRangeUser(ivar,fCfContainer->GetAxis(ivar,istep)->GetXmin(),
                                 fCfContainer->GetAxis(ivar,istep)->GetXmax(),istep);
    }
  }
  delete arr;
}

//________________________________________________________________
void AliDielectronCFdraw::Draw(const Option_t* varnames, const char* slices, const char* opt)
{
  //
  // Draw 'variables' of 'slices'
  // for multidimensional draw variables may be separated by a ':'
  // slice numbers may be separated by any of ,:;
  //
  // variables may be called by either their name or number
  //

  TObjArray *arrVars=TString(varnames).Tokenize(":");
  Int_t entries=arrVars->GetEntriesFast();
  if (entries<1||entries>3){
    AliError("Wrong number of variables, supported are 1 - 3 dimensions");
    delete arrVars;
    return;
  }
  
  TIter next(arrVars);
  TObjString *ostr=0x0;
  Int_t ivar[3]={-1,-1,-1};
  for (Int_t i=0; i<entries; ++i){
    ostr=static_cast<TObjString*>(next());
    if (ostr->GetString().IsDigit()){
      ivar[i]=ostr->GetString().Atoi();
    } else {
      ivar[i]=fCfContainer->GetVar(ostr->GetName());
    }
  }

  switch (entries){
  case 1:
    Draw(ivar[0],slices,opt);
    break;
  case 2:
    Draw(ivar[0],ivar[1],slices,opt);
    break;
  case 3:
    Draw(ivar[0],ivar[1],ivar[2],slices,opt);
    break;
  }
  delete arrVars;
}

//________________________________________________________________
void AliDielectronCFdraw::Draw(Int_t var, const char* opt, const char* slices)
{
  //
  // Draw variable var for all slices
  // slices may be divided by and of ,;:
  //
  // if opt contains 'same' all histograms are drawn in the same pad
  // otherwise the pad will be divided in sub pads and the histograms
  // are drawn in each sub pad
  //

  const Int_t ndim=1;
  Int_t vars[ndim]={var};
  TObjArray *arr=CollectHistos(ndim,vars,slices);
  Draw(arr,opt);
  delete arr; 
}

//________________________________________________________________
void AliDielectronCFdraw::Draw(Int_t var0, Int_t var1, const char* opt, const char* slices)
{
  //
  // Draw 2D case
  //
  const Int_t ndim=2;
  Int_t vars[ndim]={var0,var1};
  TObjArray *arr=CollectHistos(ndim,vars,slices);
  Draw(arr,opt);
  delete arr;
}

//________________________________________________________________
void AliDielectronCFdraw::Draw(Int_t var0, Int_t var1, Int_t var2, const char* opt, const char* slices)
{
  //
  // Draw 3D case
  //
  const Int_t ndim=3;
  Int_t vars[ndim]={var0,var1,var2};
  TObjArray *arr=CollectHistos(ndim,vars,slices);
  Draw(arr,opt);
  delete arr;
}

//________________________________________________________________
void AliDielectronCFdraw::Draw(const TObjArray *arr, const char* opt)
{
  //
  // draw all objects in arr
  //
  TString optStr(opt);
  optStr.ToLower();
  Bool_t drawSame=optStr.Contains("same");
  optStr.ReplaceAll("same","");
  if (!gPad) new TCanvas;
  TCanvas *c=gPad->GetCanvas();
  c->Clear();
  if (!drawSame){
    //optimised division
    Int_t nPads = arr->GetEntries();
    Int_t nCols = (Int_t)TMath::Ceil( TMath::Sqrt(nPads) );
    Int_t nRows = (Int_t)TMath::Ceil( (Double_t)nPads/(Double_t)nCols );
    c->Divide(nCols,nRows);
    for (Int_t i=0; i<nPads; ++i){
      c->cd(i+1);
      arr->UncheckedAt(i)->Draw(optStr.Data());
    }
  }
}

//________________________________________________________________
TObjArray* AliDielectronCFdraw::CollectHistos(Int_t dim, Int_t *vars, const char* slices)
{
  //
  // Collect histos with 'dim'ension of the 'slices' separated by one of "':;'"
  // in a TObjArray and return it
  //
  TObjArray *arr=TString(slices).Tokenize(",:;");
  TObjArray *arrHists=0x0;
  if (arr->GetEntriesFast()==0){
    // all slices in case of 0 entries
    arrHists=new TObjArray(fCfContainer->GetNStep());
    for (Int_t istep=0; istep<fCfContainer->GetNStep(); ++istep){
      arrHists->Add(Project(dim,vars,istep));
    }
  } else {
    arrHists=new TObjArray(arr->GetEntriesFast());
    TIter next(arr);
    TObjString *ostr=0x0;
    while ( (ostr=static_cast<TObjString*>(next())) ) {
      Int_t istep=ostr->GetString().Atoi();
      arrHists->Add(Project(dim,vars,istep));
    }
  }
  delete arr;

  return arrHists;
}

//________________________________________________________________
TH1* AliDielectronCFdraw::Project(Int_t ndim, Int_t *vars, Int_t slice)
{
  //
  // Do an nim projection
  //
  switch (ndim){
  case 1:
    return fCfContainer->Project(vars[0],slice);
    break;
  case 2:
    return fCfContainer->Project(vars[0],vars[1],slice);
    break;
  case 3:
    return fCfContainer->Project(vars[0],vars[1],vars[2],slice);
    break;
  }
  return 0x0;
}

