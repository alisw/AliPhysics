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

/* $Id$ */

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
#include <TClass.h>
#include <TObject.h>
#include <TVirtualPS.h>
#include <TFile.h>
#include <TString.h>
#include <TObjString.h>
#include <TVectorD.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <AliCFEffGrid.h>

#include <AliLog.h>

#include "AliDielectronCFdraw.h"

ClassImp(AliDielectronCFdraw)

AliDielectronCFdraw::AliDielectronCFdraw() :
  TNamed(),
  fCfContainer(0x0),
  fEffGrid(0x0),
  fVdata(0)
{
  //
  // Ctor
  //
}

//________________________________________________________________
AliDielectronCFdraw::AliDielectronCFdraw(const char*name, const char* title) :
  TNamed(name,title),
  fCfContainer(0x0),
  fEffGrid(0x0),
  fVdata(0)
{
  //
  // Named Ctor
  //
  
}

//________________________________________________________________
AliDielectronCFdraw::AliDielectronCFdraw(AliCFContainer *cont) :
  TNamed(cont->GetName(), cont->GetTitle()),
  fCfContainer(cont),
  fEffGrid(new AliCFEffGrid("eff","eff",*cont)),
  fVdata(0)
{
  //
  // directly set the CF container
  //

}

//________________________________________________________________
AliDielectronCFdraw::AliDielectronCFdraw(const char* filename) :
  TNamed(),
  fCfContainer(0x0),
  fEffGrid(0x0),
  fVdata(0)
{
  //
  // get CF container(s) from file 'filename'
  //
  SetCFContainers(filename);
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
    if (!cf) continue;
    nstep+=cf->GetNStep();
  }
  if (nstep==0) return;
  Int_t nbins[1]={1};
  fCfContainer=new AliCFContainer(GetName(), GetTitle(), nstep, 1, nbins);

  //delete unneeded steps
//   for (Int_t istep=0; istep<nstep; ++istep) delete fCfContainer->GetGrid(istep);

  //add step to the new container
  Int_t istep=0;
  for (Int_t icf=0; icf<arr->GetEntries(); ++icf){
    AliCFContainer *cf=dynamic_cast<AliCFContainer*>(arr->At(icf));
    if (!cf) continue;
    for (Int_t istepCurr=0; istepCurr<cf->GetNStep(); ++istepCurr){
      fCfContainer->SetGrid(istep, cf->GetGrid(istepCurr));
      fCfContainer->SetStepTitle(istep,Form("%s, Pair: %s",cf->GetTitle(),cf->GetStepTitle(istepCurr)));
      ++istep;
    }
  }
  if (fEffGrid) delete fEffGrid;
  fEffGrid=new AliCFEffGrid("eff","eff",*fCfContainer);
}

//________________________________________________________________
void AliDielectronCFdraw::SetCFContainers(const char* filename)
{
  //
  // get CF containers from file
  //

  TFile f(filename);
  TList *l=f.GetListOfKeys();
  TIter nextKey(l);
  TKey *k=0x0;
  while ( (k=static_cast<TKey*>(nextKey())) ){
    TObject *o=k->ReadObj();
    if (o->IsA()->InheritsFrom(TSeqCollection::Class())){
      TSeqCollection *arr=static_cast<TSeqCollection*>(o);
      SetCFContainers(arr);
    } else if (o->IsA()==AliCFContainer::Class()){
      fCfContainer=static_cast<AliCFContainer*>(o);
      if (fEffGrid) delete fEffGrid;
      fEffGrid=new AliCFEffGrid("eff","eff",*fCfContainer);
    }
  }
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
      fCfContainer->GetGrid(istep)->SetRangeUser(ivar,min,max);
    }
  } else {
    TIter next(arr);
    TObjString *ostr=0x0;
    while ( (ostr=static_cast<TObjString*>(next())) ) {
      Int_t istep=ostr->GetString().Atoi();
      fCfContainer->GetGrid(istep)->SetRangeUser(ivar,min,max);
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
      fCfContainer->GetAxis(ivar,istep)->SetRange(0,0);
    }
  } else {
    TIter next(arr);
    TObjString *ostr=0x0;
    while ( (ostr=static_cast<TObjString*>(next())) ) {
      Int_t istep=ostr->GetString().Atoi();
      fCfContainer->GetAxis(ivar,istep)->SetRange(0,0);
    }
  }
  delete arr;
}

//________________________________________________________________
void AliDielectronCFdraw::Draw(const Option_t* varnames, const char* opt, const char* slices)
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
  for (Int_t i=entries-1; i>=0; --i){
    ostr=static_cast<TObjString*>(next());
    if (ostr->GetString().IsDigit()){
      ivar[i]=ostr->GetString().Atoi();
    } else {
      ivar[i]=fCfContainer->GetVar(ostr->GetName());
    }
  }

  Draw(ivar[0],ivar[1],ivar[2],opt,slices);
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

  Int_t vars[3]={var,-1,-1};
  TObjArray *arr=CollectHistosProj(vars,slices);
  Draw(arr,opt);
  delete arr; 
}

//________________________________________________________________
void AliDielectronCFdraw::Draw(Int_t var0, Int_t var1, const char* opt, const char* slices)
{
  //
  // Draw 2D case
  //
  Int_t vars[3]={var0,var1,-1};
  TObjArray *arr=CollectHistosProj(vars,slices);
  Draw(arr,opt);
  delete arr;
}

//________________________________________________________________
void AliDielectronCFdraw::Draw(Int_t var0, Int_t var1, Int_t var2, const char* opt, const char* slices)
{
  //
  // Draw 3D case
  //
  Int_t vars[3]={var0,var1,var2};
  TObjArray *arr=CollectHistosProj(vars,slices);
  Draw(arr,opt);
  delete arr;
}

//________________________________________________________________
TObjArray* AliDielectronCFdraw::CollectHistosProj(const Int_t vars[3], const char* slices)
{
  //
  // Collect histos with up to 3 dimension of the 'slices' separated by one of "':;'"
  // in a TObjArray and return it
  // if a dimension is not used it must be set to -1
  //
  TObjArray *arr=TString(slices).Tokenize(",:;");
  TObjArray *arrHists=0x0;
  if (arr->GetEntriesFast()==0){
    // all slices in case of 0 entries
    arrHists=new TObjArray(fCfContainer->GetNStep());
    for (Int_t istep=0; istep<fCfContainer->GetNStep(); ++istep){
      TH1 *hproj=Project(vars,istep);
      if (!hproj) continue;
      hproj->SetName(Form("proj_%02d",istep));
      hproj->SetTitle(fCfContainer->GetStepTitle(istep));
      arrHists->Add(hproj);
    }
  } else {
    arrHists=new TObjArray(arr->GetEntriesFast());
    TIter next(arr);
    TObjString *ostr=0x0;
    while ( (ostr=static_cast<TObjString*>(next())) ) {
      Int_t istep=ostr->GetString().Atoi();
      TH1 *hproj=Project(vars,istep);
      if (!hproj) continue;
      hproj->SetName(Form("proj_%02d",istep));
      hproj->SetTitle(fCfContainer->GetStepTitle(istep));
      arrHists->Add(hproj);
    }
  }
  delete arr;

  return arrHists;
}

//________________________________________________________________
TH1* AliDielectronCFdraw::Project(const Int_t *vars, Int_t slice)
{
  //
  // Do an ndim projection
  //
  return fCfContainer->Project(slice,vars[0],vars[1],vars[2]);
}

//________________________________________________________________
TH1* AliDielectronCFdraw::Project(const Option_t* var, Int_t slice)
{
  //
  // translate variable names and do projection
  //
  TObjArray *arrVars=TString(var).Tokenize(":");
  Int_t entries=arrVars->GetEntriesFast();
  if (entries<1||entries>3){
    AliError("Wrong number of variables, supported are 1 - 3 dimensions");
    delete arrVars;
    return 0x0;
  }
  
  TObjString *ostr=0x0;
  Int_t ivar[3]={-1,-1,-1};
  for (Int_t i=entries-1; i>=0; --i){
    ostr=static_cast<TObjString*>(arrVars->At(i));
    if (ostr->GetString().IsDigit()){
      ivar[i]=ostr->GetString().Atoi();
    } else {
      ivar[i]=fCfContainer->GetVar(ostr->GetName());
    }
  }
  if (ivar[0]==-1) return 0x0;
  delete arrVars;
  return fCfContainer->Project(slice,ivar[0],ivar[1],ivar[2]);
}

//________________________________________________________________
void AliDielectronCFdraw::DrawEfficiency(const char* varnames, const char* numerators, Int_t denominator, const char* opt)
{
  //
  // plot efficiencies for variables. Variable may be up to 3 dim, separated by a ':'
  // you may have several numerators, sparated by one of ',:;'
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

  Int_t type=0;
  TString optStr(opt);
  if (optStr.Contains("2")) type=1;
  
  DrawEfficiency(ivar[2],ivar[1],ivar[0],numerators, denominator,opt,type);
  delete arrVars;
}

//________________________________________________________________
void AliDielectronCFdraw::DrawEfficiency(Int_t var, const char* numerators, Int_t denominator, const char* opt, Int_t type)
{
  //
  // Draw Efficiencies for all numerators
  // numerators may be divided by and of ,;:
  //
  // if opt contains 'same' all histograms are drawn in the same pad
  // otherwise the pad will be divided in sub pads and the histograms
  // are drawn in each sub pad
  //
  
  Int_t vars[3]={var,-1,-1};
  TObjArray *arr=CollectHistosEff(vars,numerators,denominator,type);
  TString drawOpt=opt;
  drawOpt+="eff";
  Draw(arr,drawOpt);
  delete arr;
}

//________________________________________________________________
void AliDielectronCFdraw::DrawEfficiency(Int_t var0, Int_t var1, const char* numerators, Int_t denominator, const char* opt, Int_t type)
{
  //
  // Draw 2D case
  //
  Int_t vars[3]={var0,var1,-1};
  TObjArray *arr=CollectHistosEff(vars,numerators,denominator,type);
  TString drawOpt=opt;
  drawOpt+="eff";
  Draw(arr,drawOpt);
  delete arr;
}

//________________________________________________________________
void AliDielectronCFdraw::DrawEfficiency(Int_t var0, Int_t var1, Int_t var2, const char* numerators, Int_t denominator, const char* opt, Int_t type)
{
  //
  // Draw 3D case
  //
  Int_t vars[3]={var0,var1,var2};
  TObjArray *arr=CollectHistosEff(vars,numerators,denominator,type);
  TString drawOpt=opt;
  drawOpt+="eff";
  Draw(arr,drawOpt);
  delete arr;
}

//________________________________________________________________
TObjArray* AliDielectronCFdraw::CollectHistosEff(const  Int_t vars[3], const char* numerators, Int_t denominator, Int_t type)
{
  //
  // Collect histos with 'dim'ension of the 'slices' separated by one of "':;'"
  // in a TObjArray and return it
  //
  TObjArray *arr=TString(numerators).Tokenize(",:;");
  TObjArray *arrHists=0x0;

  if (type==0){
    if (arr->GetEntriesFast()==0){
    // all slices in case of 0 entries
      arrHists=new TObjArray(fCfContainer->GetNStep());
      fVdata.ResizeTo(arrHists->GetSize());
      for (Int_t istep=0; istep<fCfContainer->GetNStep(); ++istep){
        fEffGrid->CalculateEfficiency(istep,denominator);
        TH1 *hproj=ProjectEff(vars);
        if (!hproj) continue;
        Float_t eff=fEffGrid->GetAverage();
        fVdata(istep)=eff;
        hproj->SetName(Form("eff_%02d/%02d",istep,denominator));
        hproj->SetTitle(Form("%s (%.3f)",fCfContainer->GetStepTitle(istep),eff));
        arrHists->Add(hproj);
      }
    } else {
      arrHists=new TObjArray(arr->GetEntriesFast());
      TIter next(arr);
      TObjString *ostr=0x0;
      fVdata.ResizeTo(arrHists->GetSize());
      Int_t count=0;
      while ( (ostr=static_cast<TObjString*>(next())) ) {
        Int_t istep=ostr->GetString().Atoi();
        fEffGrid->CalculateEfficiency(istep,denominator);
        TH1 *hproj=ProjectEff(vars);
        if (!hproj) continue;
        Float_t eff=fEffGrid->GetAverage();
        fVdata(count++)=eff;
        hproj->SetName(Form("eff_%02d/%02d",istep,denominator));
        hproj->SetTitle(Form("%s (%.3f)",fCfContainer->GetStepTitle(istep),eff));
        arrHists->Add(hproj);
      }
    }
  }

  //second approach
  if (type==1){
    TH1 *hDen=Project(vars,denominator);
    Double_t entriesDen=hDen->GetEffectiveEntries();
    if (arr->GetEntriesFast()==0){
    // all slices in case of 0 entries
      arrHists=new TObjArray(fCfContainer->GetNStep());
      fVdata.ResizeTo(arrHists->GetSize());
      for (Int_t istep=0; istep<fCfContainer->GetNStep(); ++istep){
        TH1 *hproj=Project(vars,istep);
        if (!hproj) continue;
        Float_t eff=0;
        if (entriesDen>0) eff=hproj->GetEffectiveEntries()/entriesDen;
        fVdata(istep)=eff;
        hproj->Divide(hproj,hDen,1,1,"B");
        hproj->SetName(Form("eff_%02d/%02d",istep,denominator));
        hproj->SetTitle(Form("%s (%.3f)",fCfContainer->GetStepTitle(istep),eff));
        arrHists->Add(hproj);
      }
    } else {
      arrHists=new TObjArray(arr->GetEntriesFast());
      fVdata.ResizeTo(arrHists->GetSize());
      TIter next(arr);
      TObjString *ostr=0x0;
      Int_t count=0;
      while ( (ostr=static_cast<TObjString*>(next())) ) {
        Int_t istep=ostr->GetString().Atoi();
        TH1 *hproj=Project(vars,istep);
        if (!hproj) continue;
        Float_t eff=0;
        if (entriesDen>0) eff=hproj->GetEffectiveEntries()/entriesDen;
        fVdata(count++)=eff;
        hproj->Divide(hproj,hDen,1,1,"B");
        hproj->SetName(Form("eff_%02d/%02d",istep,denominator));
        hproj->SetTitle(Form("%s (%.3f)",fCfContainer->GetStepTitle(istep),eff));
        arrHists->Add(hproj);
      }
    }
    delete hDen;
  }
  

  delete arr;
  return arrHists;
}

//________________________________________________________________
TH1* AliDielectronCFdraw::ProjectEff(const Int_t vars[3])
{
  //
  // Do an nim projection
  //
  return fEffGrid->Project(vars[0],vars[1],vars[2]);
}

//________________________________________________________________
void AliDielectronCFdraw::Draw(const TObjArray *arr, const char* opt)
{
  //
  // draw all objects in arr
  //
  TString optStr(opt);
  optStr.ToLower();
  Bool_t drawSame     = optStr.Contains("same");
  Bool_t drawSamePlus = optStr.Contains("same+");
  Bool_t drawEff      = optStr.Contains("eff");
  Bool_t optLeg       = optStr.Contains("leg");
  Bool_t optScaleMax  = optStr.Contains("max");
  
  if (!drawSamePlus) optStr.ReplaceAll("same","");
  
  optStr.ReplaceAll("+","");
  optStr.ReplaceAll("eff","");
  optStr.ReplaceAll("leg","");
  optStr.ReplaceAll("max","");
  
  if (!gPad) new TCanvas;
  
  Int_t nPads = arr->GetEntriesFast();
  if (nPads==0) return;
  
//   if (nPads==1){
//     arr->UncheckedAt(0)->Draw(optStr.Data());
//     return;
//   }
  
  TCanvas *c=gPad->GetCanvas();
  if (!gVirtualPS&&!drawSamePlus&&!drawSame&&nPads>1) c->Clear();
  
  if (!drawSame&&nPads>1){
    //optimised division
    Int_t nCols = (Int_t)TMath::Ceil( TMath::Sqrt(nPads) );
    Int_t nRows = (Int_t)TMath::Ceil( (Double_t)nPads/(Double_t)nCols );
    c->Divide(nCols,nRows);
    for (Int_t i=0; i<nPads; ++i){
      c->cd(i+1);
      arr->UncheckedAt(i)->Draw(optStr.Data());
    }
  } else {
    TLegend *leg=0;
    if (drawSamePlus){
      //find already existing legend to attach entries to it
      TIter nextPrimitive(gPad->GetListOfPrimitives());
      TObject *o=0x0;
      while ((o=nextPrimitive())) if (o->IsA()==TLegend::Class()) leg=(TLegend*)o;
    }
    
    if (optLeg&&!leg) leg=new TLegend(.2,.3,.99,.9);
    Int_t addColor=0;
    if (leg) addColor=leg->GetNRows();
    Int_t offset=20;
    if (nPads<7) offset=24;
    for (Int_t i=0; i<nPads; ++i){
      if (i==1&&!drawSamePlus) optStr+="same";
      TH1 *hist=static_cast<TH1*>(arr->UncheckedAt(i));
      hist->SetLineColor(i+1+addColor);
      hist->SetLineWidth(2);
      hist->SetMarkerColor(i+1+addColor);
      hist->SetMarkerStyle(offset+i+addColor);
      hist->Draw(optStr.Data());
      
      if (drawEff&&i==0&&!drawSamePlus) {
        hist->GetYaxis()->SetRangeUser(0,2);
        hist->SetYTitle("Rec. Signal / Gen. Signal");
      }
      
      if (leg) leg->AddEntry(hist,hist->GetTitle(),"lp");
    }
    if (leg){
      leg->SetFillColor(10);
      leg->SetY1(.9-leg->GetNRows()*.05);
      leg->SetY1NDC(.9-leg->GetNRows()*.05);
      leg->SetMargin(.1);
      if (!drawSamePlus) leg->Draw();
    }
    if (optScaleMax){
      TIter nextPrimitive(gPad->GetListOfPrimitives());
      TObject *o=0x0;
      TH1 *hfirst=0x0;
      Float_t max=0;
      while ((o=nextPrimitive())) if (o->InheritsFrom(TH1::Class())){
        TH1 *h=(TH1*)o;
        if (!hfirst) hfirst=h;
        if (h->GetMaximum()>max) max=h->GetMaximum();
      }
      max*=1.1;
      hfirst->SetMaximum(max);
    }
  }
  
}

//________________________________________________________________
Double_t AliDielectronCFdraw::GetAverageEfficiency(Int_t numerator, Int_t denominator, Double_t &effErr)
{
  //
  // Extract the mean efficiency of the numerator and denominator
  //

  //use variable 0 as default, since for the average it doesn't matter
  TH1 *hDen=fCfContainer->Project(denominator,0);
  Double_t entriesDen=hDen->GetEffectiveEntries();
  TH1 *hproj=fCfContainer->Project(numerator,0);
  if (!hproj) return -1.;
  Double_t entriesNum=hproj->GetEffectiveEntries();
  if (entriesDen<1||entriesNum<1) return -1;
  
  Double_t eff=-1.;
  if (entriesDen>0) eff=entriesNum/entriesDen;
  effErr=TMath::Sqrt(1./entriesNum+1./entriesDen)*eff;
  delete hDen;
  delete hproj;
  return eff;
}
