
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
//                Add general description                                 //
//                                                                       //
//                                                                       //
/*
Detailed description


*/
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include <TMath.h>
#include <TVectorT.h>
#include <TH1.h>
#include <TROOT.h>
#include <TCanvas.h>

#include <AliCFContainer.h>
#include <AliCFGridSparse.h>

#include "AliDielectronSignalBase.h"
#include "AliDielectronSignalFunc.h"

#include "AliDielectronSpectrum.h"

ClassImp(AliDielectronSpectrum)

AliDielectronSpectrum::AliDielectronSpectrum() :
  TNamed(),
  fCFSignal(0x0),
  fCFCorrection(0x0),
  fCFSpectrum(0x0),
  fCFCorrMatrix(0x0),
  fStepSignal(kTRUE),
  fStepSignificance(kFALSE),
  fStepSOB(kFALSE),
  fStepMass(kFALSE),
  fStepMassWidth(kFALSE),
  fSignalStep(-1),
  fCorrNom(-1),
  fCorrDenom(-1),
  fSignalMethods(0),
  fVariables("Pt"),
  fOwnerSpectrum(kTRUE),
  fVisualDebug(kFALSE),
  fNvars(0),
  fVars(0x0),
  fNbins(0x0),
  fCurrentBins(0x0),
  fCurrentPositions(0x0)
{
  //
  // Default Constructor
  //
  fSignalMethods.SetOwner(kFALSE);
}

//______________________________________________
AliDielectronSpectrum::AliDielectronSpectrum(const char* name, const char* title) :
  TNamed(name, title),
  fCFSignal(0x0),
  fCFCorrection(0x0),
  fCFSpectrum(0x0),
  fCFCorrMatrix(0x0),
  fStepSignal(kTRUE),
  fStepSignificance(kFALSE),
  fStepSOB(kFALSE),
  fStepMass(kFALSE),
  fStepMassWidth(kFALSE),
  fSignalStep(-1),
  fCorrNom(-1),
  fCorrDenom(-1),
  fSignalMethods(0),
  fVariables("Pt"),
  fOwnerSpectrum(kTRUE),
  fVisualDebug(kFALSE),
  fNvars(0),
  fVars(0x0),
  fNbins(0x0),
  fCurrentBins(0x0),
  fCurrentPositions(0x0)
{
  //
  // Named Constructor
  //
  fSignalMethods.SetOwner(kFALSE);
}

//______________________________________________
AliDielectronSpectrum::~AliDielectronSpectrum()
{
  //
  // Default Destructor
  //

  if (fNbins) delete [] fNbins;
  if (fVars) delete [] fVars;
  if (fCFCorrMatrix) delete fCFCorrMatrix;
  if ( fOwnerSpectrum && fCFSpectrum ) delete fCFSpectrum;
}

//______________________________________________
void AliDielectronSpectrum::Process()
{
  //
  // Extract signal and perform correction in the specified bins
  //

  //sanity checks
  if (!fCFSignal){
    AliError("Cannot perform signal extraction, no signal container set");
    return;
  }
  
  if (fSignalMethods.GetEntries()==0){
    AliWarning("No Signal extraction method specified, using a default one");
    AliDielectronSignalFunc *func=new AliDielectronSignalFunc("gaus+exp","Gauss for Signal and Exponential background");
    func->SetDefaults(1);
    fSignalMethods.Add(func);
    fSignalMethods.SetOwner();
  }

  //setup configured variables
  if (!SetupVariables()) return;
  
  //create container for the spectrum
  CreateCFSpectrum();

  if (!fCFSpectrum){
    AliError("Could not create the Spectrum container");
    return;
  }

  //get efficiency map if correction container is available
  if (fCFCorrection&&!fCFCorrMatrix){
    CreateCorrectionMatrix();
  }

  //loop over all configured bins and extract the signal
  fCurrentBins=new Int_t[fNvars];
  fCurrentPositions=new Double_t[fNvars];
  ExtractSignalInBins();
  delete [] fCurrentBins;
  delete [] fCurrentPositions;
  if (fSignalMethods.IsOwner()) {
    fSignalMethods.Delete();
    fSignalMethods.SetOwner(kFALSE);
  }
  
}

//______________________________________________
Bool_t AliDielectronSpectrum::SetupVariables()
{
  //
  // Setup the variables arrays
  //
  
  TObjArray *arr=fVariables.Tokenize(":");
  fNvars=arr->GetEntries();
  fVars=new Int_t[fNvars];
  fNbins=new Int_t[fNvars];
  
  for (Int_t iVar=0; iVar<fNvars; ++iVar){
    fVars[iVar]=fCFSignal->GetVar(arr->UncheckedAt(iVar)->GetName());
    if (fVars[iVar]==-1){
      AliError(Form("Variable '%s' not found in Signal container!",arr->UncheckedAt(iVar)->GetName()));
      delete [] fVars;
      fVars=0x0;
      delete [] fNbins;
      fNbins=0x0;
      delete arr;
      return kFALSE;
    }
    
    fNbins[iVar]=fCFSignal->GetNBins(fVars[iVar]);
  }
  delete arr;
  return kTRUE;
}

//______________________________________________
void AliDielectronSpectrum::CreateCFSpectrum()
{
  //
  // Create CF container for the spectrum
  //

  Int_t nAddStep=0;
  if (fStepSignal)       nAddStep+=2;
  if (fStepSignificance) ++nAddStep;
  if (fStepSOB)          ++nAddStep;
  if (fStepMass)         ++nAddStep;
  if (fStepMassWidth)    ++nAddStep;
  
  Int_t nStep=nAddStep*(fSignalMethods.GetEntries());
  if (fSignalMethods.GetEntries()>1) nStep+=nAddStep;

  fCFSpectrum = new AliCFContainer(GetName(), GetTitle(), nStep, fNvars, fNbins);

  // initialize the variables and their bin limits
  for (Int_t iVar=0; iVar<fNvars; iVar++) {
    fCFSpectrum->SetBinLimits(iVar, fCFSignal->GetBinLimits(fVars[iVar]));
    fCFSpectrum->SetVarTitle(iVar, fCFSignal->GetVarTitle(fVars[iVar]));
  }

  // setup step titles
  Int_t steps=0;
  for (Int_t iMethod=0; iMethod<fSignalMethods.GetEntries(); ++iMethod){
    TString name(fSignalMethods.UncheckedAt(iMethod)->GetName());
    if (fStepSignal){
      fCFSpectrum->SetStepTitle(steps++,(name+" (Siganl)").Data());
      fCFSpectrum->SetStepTitle(steps++,(name+" (Corrected Siganl)").Data());
    }
    if (fStepSignificance){
      fCFSpectrum->SetStepTitle(steps++,(name+" (Significance)").Data());
    }
    if (fStepSOB){
      fCFSpectrum->SetStepTitle(steps++,(name+" (S/B)").Data());
    }
    if (fStepMass){
      fCFSpectrum->SetStepTitle(steps++,(name+" (Mass)").Data());
    }
    if (fStepMassWidth){
      fCFSpectrum->SetStepTitle(steps++,(name+" (Mass width)").Data());
    }
  }

  if (fSignalMethods.GetEntries()>1){
    fCFSpectrum->SetStepTitle(steps++,"Mean of methods");
    fCFSpectrum->SetStepTitle(steps++,"Mean of methods (Corrected)");
  }

  if (nStep!=steps){
    AliError("Something went wrong in the step creation");
    delete fCFSpectrum;
    fCFSpectrum=0x0;
  }
}

//______________________________________________
void AliDielectronSpectrum::CreateCorrectionMatrix()
{
  //
  // Get the correction matrix for the corresponding variables
  //

  if (!fCFCorrection) return;
  
  TObjArray *arr=fVariables.Tokenize(":");
  Int_t nvars=arr->GetEntries();
  Int_t *vars=new Int_t[nvars];
  
  for (Int_t iVar=0; iVar<fNvars; ++iVar){
    vars[iVar]=fCFCorrection->GetVar(arr->UncheckedAt(iVar)->GetName());
    if (vars[iVar]==-1){
      AliError(Form("Variable '%s' not found in Correction container!",arr->UncheckedAt(iVar)->GetName()));
      delete [] vars;
      delete arr;
      return;
    }
  }
  delete arr;
  
  fCFCorrMatrix  =fCFCorrection->GetGrid(fCorrNom)->Project(nvars,vars,0,0);
  AliCFGridSparse *hnDeNom=fCFCorrection->GetGrid(fCorrDenom)->Project(nvars,vars,0,0);
  fCFCorrMatrix->Divide(hnDeNom);
  delete hnDeNom;
}

//______________________________________________
void AliDielectronSpectrum::ExtractSignalInBins(Int_t variable)
{
  //
  // recursively loop over bins and extract signal
  //

  Int_t varPairType=fCFSignal->GetVar("PairType");
  Int_t varMass=fCFSignal->GetVar("M");
  
  for (Int_t ibin=0; ibin<fNbins[variable]; ++ibin){
    Int_t bin=ibin+1;
    fCFSignal->GetGrid(fSignalStep)->GetGrid()->GetAxis(fVars[variable])->SetRange(bin,bin);
    fCurrentBins[variable]=bin;
    fCurrentPositions[variable]=fCFSignal->GetGrid(fSignalStep)->GetBinCenter(fVars[variable],bin);
    
    if (variable != fNvars-1) ExtractSignalInBins(variable+1);
    
    TObjArray arrPairTypes(10);
    arrPairTypes.SetOwner();
  
    for (Int_t itype=0; itype<3; ++itype){
//     Int_t itype=1;
      fCFSignal->SetRangeUser(varPairType,itype,itype,fSignalStep);
      TH1 *h=fCFSignal->Project(varMass,fSignalStep);
      h->SetDirectory(0);
      arrPairTypes.AddAt(h,itype);
    }
    AliInfo(Form("Processing bin: %d (%.2f)",ibin, fCurrentPositions[variable]));
    //loop over all signal extraction methods and retrieve signals
    for (Int_t iMethod=0; iMethod<fSignalMethods.GetEntries(); ++iMethod){
      AliDielectronSignalBase *signalMethod=(AliDielectronSignalBase*)fSignalMethods.At(iMethod);
      signalMethod->Process(&arrPairTypes);
      if (fVisualDebug){
        TCanvas *c=(TCanvas*)gROOT->GetListOfCanvases()->FindObject("SpectrumVisualDebug");
        if (!c) c=new TCanvas("SpectrumVisualDebug","SpectrumVisualDebug");
        c->Clear();
        TString binsProc;
        for (Int_t ivar=0; ivar<fNvars; ++ivar) binsProc+=Form("_%02d",fCurrentBins[ivar]);
        signalMethod->Draw("stat");
        binsProc.Append(".png");
        binsProc.Prepend("SpectrumVisualDebug");
        c->Update();
        c->Print(binsProc.Data());
      }
      const TVectorD &val=signalMethod->GetValues();
      const TVectorD &err=signalMethod->GetErrors();
      
      //Fill sparse containers
      Int_t step=0;
      if (fStepSignal){
        //signal
        fCFSpectrum->GetGrid(step)->SetElement(fCurrentBins,val(0));
        fCFSpectrum->GetGrid(step)->SetElementError(fCurrentBins,err(0));
        ++step;

        //corrected signal
        if (fCFCorrMatrix){
          Float_t corrFactor = fCFCorrMatrix->GetElement(fCurrentPositions);
          Float_t corrError  = fCFCorrMatrix->GetElementError(fCurrentPositions);

          Float_t value=val(0)*corrFactor;
          fCFSpectrum->GetGrid(step)->SetElement(fCurrentBins,value);
          Float_t error=TMath::Sqrt( (err(0)/val(0))*(err(0)/val(0)) +
                                     (corrError/corrFactor)*(corrError/corrFactor)
                                   )*value;
          fCFSpectrum->GetGrid(step)->SetElementError(fCurrentBins,error);
//           printf("corrFactor: %f+-%f\n",value,error);
        }
        ++step;
      }

      if (fStepSignificance) {
        fCFSpectrum->GetGrid(step)->SetElement(fCurrentBins,val(2));
        fCFSpectrum->GetGrid(step)->SetElementError(fCurrentBins,err(2));
        ++step;
      }
      
      if (fStepSOB) {
        fCFSpectrum->GetGrid(step)->SetElement(fCurrentBins,val(3));
        fCFSpectrum->GetGrid(step)->SetElementError(fCurrentBins,err(3));
        ++step;
      }

      if (fStepMass) {
        fCFSpectrum->GetGrid(step)->SetElement(fCurrentBins,val(4));
        fCFSpectrum->GetGrid(step)->SetElementError(fCurrentBins,err(4));
        ++step;
      }
      
      if (fStepMassWidth) {
        fCFSpectrum->GetGrid(step)->SetElement(fCurrentBins,val(5));
        fCFSpectrum->GetGrid(step)->SetElementError(fCurrentBins,err(5));
        ++step;
      }
    }// end method loop
    
    arrPairTypes.Delete();
    
  }// end bin loop
}

