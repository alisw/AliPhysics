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
//       Dielectron Histogram framework helper                     //
//                                                                       //
/*









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

//ClassImp(AliDielectronHFhelper)

//const char* AliDielectronHFhelper::fCutVars[AliDielectronHFhelper::kMaxCuts] = {"","","","","","","","","",""};

//________________________________________________________________
AliDielectronHFhelper::AliDielectronHFhelper(const char* filename) :
  TNamed(),
  fArrPairType(0x0),
  fCutVars(0x0),
  fCutLowLimits(0),
  fCutUpLimits(0)
{
  //
  // get HF container(s) from file 'filename'
  //
  SetHFArray(filename);
  
}

//________________________________________________________________
AliDielectronHFhelper::~AliDielectronHFhelper()
{
  //
  // dtor
  //
  if(fArrPairType) delete fArrPairType;
  if(fCutVars)     delete fCutVars;
  
}

//________________________________________________________________
void AliDielectronHFhelper::SetHFArray(const char* filename)
{
  //
  // get HF containers from file
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

	if( objname.Contains("_HF") && obj->IsA()==TObjArray::Class()) {
	  fArrPairType = new TObjArray( *(dynamic_cast<TObjArray*>(obj)) );
	  //fArrPairType->Print();
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
  //  Int_t size=sizeof(fCutVars)/sizeof(const char*);
  Int_t size=fCutLowLimits.GetNrows();

  if(size>=kMaxCuts) return;

  // arrays
  if(!fCutVars) {
    fCutVars = new TObjArray();
    fCutVars->SetOwner();
  }
  fCutLowLimits.ResizeTo(size+1);
  fCutUpLimits.ResizeTo(size+1);

  // fill
  //  fCutVars[size]=Form("%s%s",(leg?"Leg":""),varname);
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
TObjArray* AliDielectronHFhelper::CollectHistos() 
{
  //
  // collect histograms for all kind of pair types or sources
  //
  
  TObjArray *collection = new TObjArray(AliDielectron::kEv1PMRot+1);

  // loop over max. available pair types
  for(Int_t i=0; i<AliDielectron::kEv1PMRot+1; i++) {

    collection->AddAt(GetHistogram(AliDielectron::PairClassName(i)),i);
    
  }

  return collection;
}

//________________________________________________________________
TH1F* AliDielectronHFhelper::GetHistogram(const char *step)
{
  //
  // main function to recive a single histogram
  //
  
  TObjArray *histos = (TObjArray*) fArrPairType->FindObject(step);
  AliInfo(Form(" Step %s selected",step));
  if(!histos) return 0x0;

  TH1F *hist = FindHistograms(histos);
  return hist;
}

//________________________________________________________________
TH1F* AliDielectronHFhelper::FindHistograms(TObjArray *histos)
{
  //
  // exclude histograms
  //
  
  // debug
  // TString title    = histos->At(0)->GetTitle();
  // TObjArray* vars  = title.Tokenize(":");
  // AliDebug(1,Form(" number of cuts/vars: %d/%d",fCutLowLimits.GetNrows(),vars->GetEntriesFast()));  
  
  // check for missing cuts
  CheckCuts(histos);
  
  // loop over all cuts
  for(Int_t icut=0; icut<fCutLowLimits.GetNrows(); icut++) {

    Bool_t bFndBin = kFALSE; // exact bin found
    const char *cutvar   = fCutVars->At(icut)->GetName();
    Double_t min   = fCutLowLimits(icut);
    Double_t max   = fCutUpLimits(icut);
    AliDebug(1,Form(" Cut %d: %s [%.2f,%.2f]",icut,cutvar,min,max));

    // loop over all histograms
    for(Int_t i=0; i<histos->GetEntriesFast(); i++) {

      // continue if already empty
      if(!histos->At(i)) continue;

      // collect binning from histo title
      TString title    = histos->At(i)->GetTitle();
      if(title.IsNull()) continue;
      AliDebug(1,Form(" histo title: %s",title.Data()));

      TObjArray *vars  = title.Tokenize(":");
      for(Int_t ivar=0; ivar<vars->GetEntriesFast(); ivar++) {
	TString binvar = vars->At(ivar)->GetName();
	AliDebug(1,Form(" Check ivar %d binvar %s",ivar,binvar.Data()));
	
	// check for cuts and ranges by the user
	if(binvar.Contains(cutvar)) {
	  TObjArray *limits = binvar.Tokenize("#");
	  
	  Double_t binmin = atof(limits->At(1)->GetName());
	  Double_t binmax = atof(limits->At(2)->GetName());
	  AliDebug(1,Form(" bin %s var %s [%.2f,%.2f]",binvar.Data(),limits->At(0)->GetName(),binmin,binmax));
	  
	  // remove histogram from array
	  if(binmin < min || binmax < min || binmin > max || binmax > max ) {
	    AliDebug(1,Form(" removed, out of range min %.2f,%.2f  max %.2f,%.2f",binmin,min,binmax,max));
	    histos->AddAt(0x0,i);
	  }
	  if(bFndBin && !(binmin == min && binmax == max)) {
	      histos->AddAt(0x0,i);
	      AliDebug(1,Form(" removed, within range min %.2f,%.2f  max %.2f,%.2f",binmin,min,binmax,max));
	    }
	  // clean up
	  if(limits) delete limits;

	  // do we have found an exact bin
	  if(binmin==min && binmax==max) bFndBin=kTRUE;

	}

      }
      // clean up
      if(vars)   delete vars;
      
    }
    
  }
  
  // compress the array by removing all empty histos
  histos->Compress();  
  AliDebug(1,Form(" Compression: %d histograms left",histos->GetEntriesFast()));

  // merge histograms
  TH1F* hist = MergeHistos(histos);
  return hist;
}

//________________________________________________________________
TH1F* AliDielectronHFhelper::MergeHistos(TObjArray *arr)
{
  //
  // merge histos to one single histogram
  //
  
  if(arr->GetEntriesFast()<1) { AliError("[E] No more histosgrams left!"); return 0x0; }

  TH1F *final=(TH1F*) arr->At(0)->Clone();  
  if(!final) return 0x0;
  
  final->Reset("CE");
  final->SetTitle(""); //TODO: change in future
  for(Int_t i=0; i<arr->GetEntriesFast(); i++) {
    final->Add((TH1F*)arr->At(i));
  }
  arr->Clear();

  return final;
}

//________________________________________________________________
void AliDielectronHFhelper::CheckCuts(TObjArray *arr)
{
  //
  // Compare histo binning and cut variables. Add necessary cuts (largest limit)
  //


  // build array with bin variable, minimum and maximum bin values
  TString titleFIRST       = arr->First()->GetTitle();
  TString titleLAST        = arr->Last()->GetTitle();
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
      AliWarning("Not all bin variables covered. Add additional cut!");
      SetRangeUser(binvar.Data(),binmin[ivar],binmax[ivar]);
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

  // pairtypes, steps and sources
  AliInfo(Form(" Number of steps: \t %d",fArrPairType->GetEntriesFast()));
  for(Int_t i=0; i<fArrPairType->GetEntriesFast(); i++) {
    AliInfo(Form(" Step %d:",0));
  }

}

//________________________________________________________________
void AliDielectronHFhelper::PrintCuts()
{

  //
  // Print out object contents
  //

  // loop over all cuts
  AliInfo(" Selected cuts:");
  for(Int_t icut=0; icut<fCutLowLimits.GetNrows(); icut++)
    AliInfo(Form(" %d: %s [%.2f,%.2f]",icut,fCutVars->At(icut)->GetName(),fCutLowLimits(icut),fCutUpLimits(icut)));
  
}

