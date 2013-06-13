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
//       Dielectron Histogram framework helper                           //
//                                                                       //
/*

A helper class to extract objects(histograms and/or profiles) from
a AliDielctronHF array of objects.


How to use it:

  AliDielectronHFhelper *hf = new AliDielectronHFhelper("path/to/the/output/file.root", "ConfigName");
  // print the structure
  hf->Print();

  //apply some cuts and print them
  hf->SetRangeUser("cut1name",cutmin1,cutmax1);
  hf->SetRangeUser(AliDielectronVarManager::kPt,ptmin,ptmax);
  hf->PrintCuts();

  // collect 1-,2- or 3-dim histograms or profiles with error option (default:"")
  TObjArray *arrHists = hf->CollectHistos(AliDielectronVarManager::kM);
  TObjArray *arrProfs = hf->CollectProfiles("",AliDielectronVarManager::kM,AliDielectronVarManager::kPt);

  // then you are left with an array of histograms for all pair types or MC signals

*/
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include <TObjArray.h>
#include <TKey.h>
#include <TList.h>
#include <TClass.h>
#include <TObject.h>
#include <TFile.h>
#include <TString.h>
#include <TObjString.h>
#include <TVectorD.h>
#include <TMath.h>
#include <TH1.h>

#include <AliLog.h>

#include "AliDielectron.h"
#include "AliDielectronHFhelper.h"
#include "AliDielectronHF.h"

ClassImp(AliDielectronHFhelper)

//const char* AliDielectronHFhelper::fCutVars[AliDielectronHFhelper::kMaxCuts] = {"","","","","","","","","",""};

//________________________________________________________________
AliDielectronHFhelper::AliDielectronHFhelper(const char* filename, const char* container) :
  TNamed(),
  fMainArr(0x0),
  fCutVars(0x0),
  fCutLowLimits(0),
  fCutUpLimits(0)
{
  //
  // get HF container(s) from file 'filename'
  //
  SetHFArray(filename, container);

}

//________________________________________________________________
AliDielectronHFhelper::~AliDielectronHFhelper()
{
  //
  // dtor
  //
  if(fMainArr) delete fMainArr;
  if(fCutVars)     delete fCutVars;

}

//________________________________________________________________
void AliDielectronHFhelper::SetHFArray(const char* filename, const char* container)
{
  //
  // get HF container from file
  //

  TFile *f = TFile::Open(filename);

  TList *l=f->GetListOfKeys();
  TIter nextKey(l);
  TKey *k=0x0;
  while ( (k=static_cast<TKey*>(nextKey())) ){

    TObject *o=k->ReadObj();
    if (o->IsA()==TList::Class()){

      TList *tlist=(TList*)o;

      TIter next(tlist);
      TObject *obj=0x0;
      while ((obj = next())) {
	TString objname(obj->GetName());

	if( objname.Contains(Form("%s_HF",container)) && obj->IsA()==TObjArray::Class()) {
	  fMainArr = new TObjArray( *(dynamic_cast<TObjArray*>(obj)) );
	  //fMainArr->Print();
	  return;
	}
      }
    }
  }

}
//________________________________________________________________
void AliDielectronHFhelper::SetRangeUser(const char *varname, Double_t min, Double_t max, Bool_t leg)
{
  //
  // Set range from variable name
  //

  Int_t size=fCutLowLimits.GetNrows();

  // check if cut is already set
  for(Int_t icut=0; icut<size; icut++) {
    TString cutName = fCutVars->At(icut)->GetName();
    if(!cutName.CompareTo(Form("%s%s",(leg?"Leg":""),varname))) {
      UnsetRangeUser(varname,leg);
      SetRangeUser(varname, min, max, leg);
      return;
    }
  }

  if(size>=kMaxCuts) return;

  // arrays
  if(!fCutVars) {
    fCutVars = new TObjArray();
    fCutVars->SetOwner();
  }
  fCutLowLimits.ResizeTo(size+1);
  fCutUpLimits.ResizeTo(size+1);

  // fill
  TObjString *str = new TObjString(Form("%s%s",(leg?"Leg":""),varname));
  fCutVars->Add(str);
  fCutLowLimits(size) = min;
  fCutUpLimits(size)  = max;
  AliWarning(Form(" %s [%.2f,%.2f]",fCutVars->At(size)->GetName(),fCutLowLimits(size),fCutUpLimits(size)));
}

//________________________________________________________________
void AliDielectronHFhelper::SetRangeUser(AliDielectronVarManager::ValueTypes type, Double_t min, Double_t max, Bool_t leg)
{
  //
  // Set range from AliDielectronVarManager
  //
  SetRangeUser(AliDielectronVarManager::GetValueName(type), min, max, leg);
}

//________________________________________________________________
void AliDielectronHFhelper::UnsetRangeUser(const char *varname, Bool_t leg)
{
  //
  // unset range from variable name
  //

  Int_t size=fCutLowLimits.GetNrows();
  TVectorD newlow;
  TVectorD newup;

  // find cut and build new vectors w/o it
  Int_t ientries = 0;
  for(Int_t icut=0; icut<size; icut++) {

    TString cutName = fCutVars->At(icut)->GetName();
    if(!cutName.CompareTo(Form("%s%s",(leg?"Leg":""),varname))) {
      fCutVars->AddAt(0x0,icut);
      continue;
    }

    // fill new vectors
    newlow.ResizeTo(ientries+1);
    newup.ResizeTo(ientries+1);
    newlow(ientries) = fCutLowLimits(icut);
    newup(ientries)  = fCutUpLimits(icut);

    ientries++;
  }

  // adapt new arrays/vectors
  fCutVars->Compress();

  fCutLowLimits.ResizeTo(ientries);
  fCutUpLimits.ResizeTo(ientries);
  for(Int_t icut=0; icut<ientries; icut++) {
    fCutLowLimits(icut) = newlow(icut);
    fCutUpLimits(icut)  = newup(icut);
  }
  // PrintCuts();

}

//________________________________________________________________
void AliDielectronHFhelper::UnsetRangeUser(AliDielectronVarManager::ValueTypes type, Bool_t leg)
{
  //
  // Unset range from AliDielectronVarManager
  //
  UnsetRangeUser(AliDielectronVarManager::GetValueName(type), leg);
}

//________________________________________________________________
TObjArray* AliDielectronHFhelper::CollectProfiles(TString option,
						  AliDielectronVarManager::ValueTypes varx,
						  AliDielectronVarManager::ValueTypes vary,
						  AliDielectronVarManager::ValueTypes varz,
						  AliDielectronVarManager::ValueTypes vart)
{
  //
  // collect 1-3 dimensional TProfiles for all kind of pair types or sources
  //

  // reconstruct the histogram/profile name resp. key
  Int_t dim    = 0;
  if(varx < AliDielectronVarManager::kNMaxValues) dim++;
  if(vary < AliDielectronVarManager::kNMaxValues) dim++;
  if(varz < AliDielectronVarManager::kNMaxValues) dim++;
  Bool_t bPairClass=0;
  if( varx < AliDielectronVarManager::kPairMax ||
      vary < AliDielectronVarManager::kPairMax ||
      varz < AliDielectronVarManager::kPairMax ||
      vart < AliDielectronVarManager::kPairMax  ) bPairClass=kTRUE;
  Bool_t bprf = !(option.Contains("hist",TString::kIgnoreCase));
  Bool_t bStdOpt=kTRUE;
  if(option.Contains("S",TString::kIgnoreCase) || option.Contains("rms",TString::kIgnoreCase)) bStdOpt=kFALSE;

  TString key = "";
  if(bprf) dim--;
  switch(dim) {
  case 3:
    key+=Form("%s_",AliDielectronVarManager::GetValueName(varx));
    key+=Form("%s_",AliDielectronVarManager::GetValueName(vary));
    key+=Form("%s",AliDielectronVarManager::GetValueName(varz));
    if(bprf) key+=Form("-%s%s",AliDielectronVarManager::GetValueName(vart),(bStdOpt ? "avg" : "rms"));
    break;
  case 2:
    key+=Form("%s_",AliDielectronVarManager::GetValueName(varx));
    key+=Form("%s",AliDielectronVarManager::GetValueName(vary));
    if(bprf) key+=Form("-%s%s",AliDielectronVarManager::GetValueName(varz),(bStdOpt ? "avg" : "rms"));
    break;
  case 1:
    key+=Form("%s",AliDielectronVarManager::GetValueName(varx));
    if(bprf) key+=Form("-%s%s",AliDielectronVarManager::GetValueName(vary),(bStdOpt ? "avg" : "rms"));
    break;
  }
  // to differentiate btw. leg and pair histos
  if(bPairClass) key.Prepend("p");
  // prepend HF
  key.Prepend("HF_");

  //TODO:  printf("--------> KEY: %s \n",key.Data());
  // get the requested object from the arrays
  TObjArray *collection = new TObjArray(fMainArr->GetEntriesFast());

  TObjArray *cloneArr = (TObjArray*) fMainArr->Clone("tmpArr");
  if(!cloneArr) return 0x0;
  cloneArr->SetOwner(kTRUE);

  // loop over all pair types or sources
  for(Int_t i=0; i<cloneArr->GetEntriesFast(); i++) {
    if(!cloneArr->At(i))                             continue;
    if(!((TObjArray*)cloneArr->At(i))->GetEntries()) continue;
    TString stepName = cloneArr->At(i)->GetName();

    // find histogram of interest from array of objects
    TObjArray *arr = (TObjArray*) GetObject(cloneArr->At(i)->GetName(),cloneArr);
    if(arr) {
      collection->AddAt(arr->FindObject(key.Data()), i);
      
      // modify the key name
      stepName.ReplaceAll("(","");
      stepName.ReplaceAll(")","");
      stepName.ReplaceAll(": ","");
      stepName.ReplaceAll(" ","_");
      stepName.ReplaceAll("Signal","");
      ((TH1*)collection->At(i))->SetName(Form("%s_%s",key.Data(),stepName.Data()));
    }

  }

  // clean up the clone
  if(cloneArr) {
    delete cloneArr;
    cloneArr=0;
  }

  return collection;
}

//________________________________________________________________
TObject* AliDielectronHFhelper::GetObject(const char *step, TObjArray *histArr)
{
  //
  // main function to recieve a pair type or MC signal array
  //

  AliDebug(1,Form(" Step %s selected",step));
  // TODO: check memory

  TObjArray *stepArr = 0x0; // this is the requested step
  TObject *hist      = 0x0;
  if(!histArr) {
    stepArr = (TObjArray*) fMainArr->FindObject(step)->Clone("tmpArr");
  }
  else {
    stepArr = (TObjArray*) histArr->FindObject(step);
  }

  if(stepArr) hist   = FindObjects(stepArr);
  return hist;

}

//________________________________________________________________
TObject* AliDielectronHFhelper::FindObjects(TObjArray *stepArr)
{
  //
  // apply cuts and exclude objects from the array
  //

  // debug
  // TString title    = stepArr->At(0)->GetTitle();
  // TObjArray* vars  = title.Tokenize(":");
  // AliDebug(1,Form(" number of cuts/vars: %d/%d",fCutLowLimits.GetNrows(),vars->GetEntriesFast()));

  // check for missing cuts
  CheckCuts(stepArr);

  // loop over all cuts
  for(Int_t icut=0; icut<fCutLowLimits.GetNrows(); icut++) {

    Bool_t bFndBin      = kFALSE; // bin with exact limits found
    const char *cutvar  = fCutVars->At(icut)->GetName(); // cut variable
    Double_t min        = fCutLowLimits(icut);           // lower limit
    Double_t max        = fCutUpLimits(icut);            // upper limit
    AliDebug(5,Form(" Cut %d: %s [%.2f,%.2f]",icut,cutvar,min,max));

    // loop over the full grid of given step
    for(Int_t i=0; i<stepArr->GetEntriesFast(); i++) {

      // continue if already empty
      if(!stepArr->At(i)) continue;

      // collect bins from the name
      TString title    = stepArr->At(i)->GetName();
      if(title.IsNull()) continue;
      AliDebug(5,Form(" %03d object name: %s",i,title.Data()));

      // loop over all variables
      TObjArray *vars  = title.Tokenize(":");
      for(Int_t ivar=0; ivar<vars->GetEntriesFast(); ivar++) {
	TString binvar = vars->At(ivar)->GetName();
	AliDebug(10,Form(" --> %d check bin var %s",ivar,binvar.Data()));

	// check cuts set by the user, and compare to the bin limits
	if(binvar.Contains(cutvar)) {
	  // bin limits
	  TObjArray *limits = binvar.Tokenize("#");
	  Double_t binmin = atof(limits->At(1)->GetName()); // lower bin limit
	  Double_t binmax = atof(limits->At(2)->GetName()); // upper bin limit
	  AliDebug(10,Form(" bin %s var %s [%.2f,%.2f]",binvar.Data(),limits->At(0)->GetName(),binmin,binmax));
	  if(limits) delete limits;

	  // cut and remove objects from the array
	  if(binmin < min || binmax < min || binmin > max || binmax > max ) {
	    AliDebug(10,Form(" removed, out of range! lower bin,cut: %.2f,%.2f or upper bin,cut: %.2f>%.2f",binmin,min,binmax,max));
	    stepArr->AddAt(0x0,i);
	  }
	  if(bFndBin && !(binmin == min && binmax == max)) {
	    stepArr->AddAt(0x0,i);
	    AliDebug(10,Form(" removed, within range! lower bin,cut: %.2f,%.2f or upper bin,cut: %.2f,%.2f",binmin,min,binmax,max));
	  }

	  // did we found a bin with exact the cut ranges
	  // this can happen only once per variable
	  if(binmin==min && binmax==max) bFndBin=kTRUE;

	}

      }
      // clean up
      if(vars)   delete vars;

    }

  }

  // compress the array by removing all empty entries
  stepArr->Compress();
  AliDebug(1,Form(" Compression: %d objects left",stepArr->GetEntriesFast()));

  // merge left objects
  TObject* hist = Merge(stepArr);
  //  if(hist) AliDebug(1,Form(" Merging: %e  entries",hist->GetEntries()));
  return hist;
}

//________________________________________________________________
TObject* AliDielectronHFhelper::Merge(TObjArray *arr)
{
  //
  // merge left objects to a single one
  //

  if(arr->GetEntriesFast()<1) { AliError(" No more objects left!"); return 0x0; }

  TObject *final=arr->At(0)->Clone();
  if(!final) return 0x0;

  TList listH;
  TString listHargs;
  listHargs.Form("((TCollection*)0x%lx)", (ULong_t)&listH);
  Int_t error = 0;

  //  final->Reset("CE");
  //  final->SetTitle(""); //TODO: change in future
  for(Int_t i=1; i<arr->GetEntriesFast(); i++) {
    listH.Add(arr->At(i));
    //   printf("%d: ent %.0f \n",i,((TH1*)((TObjArray*)arr->At(i))->At(0))->GetEntries());
    //    final->Add((TH1*)arr->At(i));
  }
  //  arr->Clear();

  final->Execute("Merge", listHargs.Data(), &error);
  return final;
}

//________________________________________________________________
void AliDielectronHFhelper::CheckCuts(TObjArray *arr)
{
  //
  // Compare binning and cut variables. Add necessary cuts (full range, no exclusion)
  //

  // build array with bin variable, minimum and maximum bin values
  TString titleFIRST       = arr->First()->GetName();
  TString titleLAST        = arr->Last()->GetName();
  TObjArray* binvarsF  = titleFIRST.Tokenize(":#");
  TObjArray* binvarsL  = titleLAST.Tokenize(":#");
  Double_t binmin[kMaxCuts]= {0.0};
  Double_t binmax[kMaxCuts]= {0.0};
  for(Int_t ivar=0; ivar<binvarsF->GetEntriesFast(); ivar++) {

    TString elementF=binvarsF->At(ivar)->GetName();
    TString elementL=binvarsL->At(ivar)->GetName();
    AliDebug(1,Form(" binvar %d: %s,%s",ivar,elementF.Data(),elementL.Data()));

    switch(ivar%3) {
    case 0: continue; break;
    case 1: binmin[(int)ivar/3]=atof(elementF.Data()); break;
    case 2: binmax[(int)ivar/3]=atof(elementL.Data()); break;
    }

    binvarsF->AddAt(0x0,ivar);
  }
  binvarsF->Compress();

  // loop over all vars and cuts, check for missing stuff
  for(Int_t ivar=0; ivar<binvarsF->GetEntriesFast(); ivar++) {

    TString binvar=binvarsF->At(ivar)->GetName();
    Bool_t selected=kFALSE;

    AliDebug(1,Form(" check cuts %d %s [%.2f,%.2f]",ivar,binvar.Data(),binmin[ivar],binmax[ivar]));
    // loop over all cuts and check for missing stuff
    for(Int_t icut=0; icut<fCutLowLimits.GetNrows(); icut++) {
      if(binvar.Contains(fCutVars->At(icut)->GetName())) { selected=kTRUE; break; }
      //      else break;
    }

    // add missing cut with max limits
    if(!selected) {
      AliWarning(Form(" Bin variable %s not covered. Add cut!",binvar.Data()));
      Bool_t leg = binvar.BeginsWith("Leg");
      if(leg) binvar.Remove(0,3);
      SetRangeUser(binvar.Data(),binmin[ivar],binmax[ivar], leg);
    }

  }

  // clean up
  if(binvarsF) delete binvarsF;
  if(binvarsL) delete binvarsL;
}

//________________________________________________________________
void AliDielectronHFhelper::Print(const Option_t* /*option*/) const
{

  //
  // Print out object contents
  //
  AliInfo(Form(" Container:               %s",fMainArr->GetName()));

  // pairtypes, steps and sources
  Int_t stepLast=0;
  AliInfo(Form(" Number of filled steps:  %d",fMainArr->GetEntries()));
  for(Int_t istep=0; istep<fMainArr->GetEntriesFast(); istep++) {
    if(!fMainArr->At(istep))                             continue;
    if(!((TObjArray*)fMainArr->At(istep))->GetEntries()) continue;
    AliInfo(Form(" step %d: %s",istep,fMainArr->At(istep)->GetName()));
    stepLast=istep;
  }

  AliInfo(Form(" Number of objects:    %d",
	       ((TObjArray*) ((TObjArray*)fMainArr->At(stepLast)) ->First())->GetEntriesFast()));

  TString title       = ((TObjArray*)fMainArr->At(stepLast))->First()->GetName();
  TObjArray* binvars  = title.Tokenize(":");
  AliInfo(Form(" Number of variables:     %d",binvars->GetEntriesFast()));
  delete binvars;

  TObjArray* binvars2  = title.Tokenize(":#");
  for(Int_t ivar=0; ivar<binvars2->GetEntriesFast(); ivar++) {
    if(ivar%3) continue;
    AliInfo(Form(" variable %.0f: %s",((Double_t)ivar)/3+1,binvars2->At(ivar)->GetName()));
  }
  delete binvars2;

}

//________________________________________________________________
void AliDielectronHFhelper::PrintCuts()
{

  //
  // Print cuts
  //

  // loop over all cuts
  AliInfo(" Selected cuts:");
  for(Int_t icut=0; icut<fCutLowLimits.GetNrows(); icut++)
    AliInfo(Form(" %d: %s [%.2f,%.2f]",icut,fCutVars->At(icut)->GetName(),fCutLowLimits(icut),fCutUpLimits(icut)));

}

