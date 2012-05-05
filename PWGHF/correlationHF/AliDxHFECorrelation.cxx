// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE Project            * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  Sedat Altinpinar <Sedat.Altinpinar@cern.ch>           *
//*                  Hege Erdal       <hege.erdal@gmail.com>               *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/// @file   AliDxHFECorrelation.cxx
/// @author Sedat Altinpinar, Hege Erdal, Matthias Richter
/// @date   2012-04-25
/// @brief  Worker class for D0-HF electron correlation
///

#include "AliDxHFECorrelation.h"
#include "AliVParticle.h"
#include "AliLog.h"
#include "TObjArray.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"
#include <iostream>
#include <cerrno>
#include <memory>

using namespace std;

ClassImp(AliDxHFECorrelation)

AliDxHFECorrelation::AliDxHFECorrelation(const char* name)
  : TNamed(name?name:"AliDxHFECorrelation", "")
  , fHistograms(NULL)
{
  // default constructor
  // 
  //

}

AliDxHFECorrelation::~AliDxHFECorrelation()
{
  // destructor
  //
  //
  if (fHistograms) {
    delete fHistograms;
    fHistograms=NULL;
  }
}

int AliDxHFECorrelation::Init()
{
  /// init class and create histograms
  if (fHistograms) delete fHistograms;
  fHistograms=new TObjArray;
  if (!fHistograms) return -ENOMEM;
  fHistograms->SetOwner(kTRUE);
  AliInfo(Form("initializing %s ", GetName()));

  // TODO: This is just a mockup to illustrate the functionality
  // the specific histograms need to be defined

  // avoid the objects to be added to the global directory 
  bool statusAddDirectory=TH1::AddDirectoryStatus();
  TH1::AddDirectory(false);
  const double Pii=TMath::Pi();
  TObjArray* a=fHistograms;
  a->AddAt(new TH1D("hD0pT"       , "D0 pT"              ,100,0.0,10.0)   , khD0pT        );
  a->AddAt(new TH1D("hD0Phi"      , "D0 phi"             ,180,0.0,2*Pii)  , khD0Phi       );
  a->AddAt(new TH1D("hD0Eta"      , "D0 eta"             ,100,-10.0,10.0) , khD0Eta       );
  a->AddAt(new TH1D("hElectronpT" , "Electron pT"        ,100,0.0,10.0)   , khElectronpT  );
  a->AddAt(new TH1D("hElectronPhi", "Electron phi"       ,180,0.0,2*Pii)  , khElectronPhi );
  a->AddAt(new TH1D("hElectronEta", "Electron eta"       ,100,-10.0,10.0) , khElectronEta );
  a->AddAt(new TH1D("hDeltaPhi"   , "#Delta#Phi D0 - HFE",180,0.0,2*Pii)  , khDeltaPhi    );

  TH1::AddDirectory(statusAddDirectory);
  return 0;
}

int AliDxHFECorrelation::Fill(const TObjArray* candidatesD0, const TObjArray* candidatesElectron)
{
  /// fill histograms from array of AliVParticle objects
  if (!candidatesD0 || !candidatesElectron) return -EINVAL;
  if (!fHistograms) {
    Init();
  }
  if (!fHistograms) {
    AliError("Initialisation failed, can not fill histograms");
  }

  const double Pii=TMath::Pi();

  TIter itrigger(candidatesD0);
  TObject* otrigger=NULL;
  int ctrigger=-1;
  while ((otrigger=itrigger())!=NULL) {
    // loop over trigger D0 particle
    ctrigger++;
    AliVParticle* ptrigger=reinterpret_cast<AliVParticle*>(otrigger);
    if (!ptrigger) continue;
    ((TH1D*)fHistograms->At(khD0pT)) ->Fill(ptrigger->Pt());
    ((TH1D*)fHistograms->At(khD0Phi))->Fill(ptrigger->Phi());
    ((TH1D*)fHistograms->At(khD0Eta))->Fill(ptrigger->Eta());
    // TODO: add further correlation specific cuts here, e.g acceptance
    // which are no primarily part of the particle selection

    TIter iElectron(candidatesElectron);
    TObject* oElectron=NULL;
    int cElectron=-1;
    while ((oElectron=iElectron())!=NULL) {
      // loop over electrons
      cElectron++;
      AliVParticle* pElectron=reinterpret_cast<AliVParticle*>(oElectron);
      if (!pElectron) continue;
      ((TH1D*)fHistograms->At(khElectronpT)) ->Fill(pElectron->Pt());
      ((TH1D*)fHistograms->At(khElectronPhi))->Fill(pElectron->Phi());
      ((TH1D*)fHistograms->At(khElectronEta))->Fill(pElectron->Eta());

      // phi difference to trigger particle
      Double_t DeltaPhi = ptrigger->Phi() - pElectron->Phi();
      if (DeltaPhi<-0.5*Pii) DeltaPhi += 2*Pii;
      if (DeltaPhi>1.5*Pii)  DeltaPhi -= 2*Pii;
      ((TH1D*)fHistograms->At(khDeltaPhi))->Fill(DeltaPhi);

    } // loop over electrons
  } // loop over D0 trigger particle

  return 0;
}

void AliDxHFECorrelation::Clear(Option_t * /*option*/)
{
  /// overloaded from TObject: cleanup

  // nothing to be done so far
  return TObject::Clear();
}

void AliDxHFECorrelation::Print(Option_t */*option*/) const
{
  /// overloaded from TObject: print info
  cout << "====================================================================" << endl;
  TObject::Print();
  if (fHistograms) {
    fHistograms->Print();
  }
}

void AliDxHFECorrelation::Draw(Option_t */*option*/)
{
  /// overloaded from TObject: draw histograms
  if (fHistograms) {
    TString name;
    int canvasno=1;
    int padno=1;
    const char* drawoption="";
    name.Form("%s_%d", GetName(), canvasno++);
    TCanvas* c=new TCanvas(name);
    c->SetWindowSize(1600,800);
    c->SetTitle(Form("%s: particle properties", GetName()));
    c->Divide(3,2);
    padno=1;
    c->cd(padno++); fHistograms->At(khD0pT)       ->Draw(drawoption);
    c->cd(padno++); fHistograms->At(khD0Phi)      ->Draw(drawoption);
    c->cd(padno++); fHistograms->At(khD0Eta)      ->Draw(drawoption);
    c->cd(padno++); fHistograms->At(khDeltaPhi)   ->Draw(drawoption);
    c->cd(padno++); fHistograms->At(khElectronpT) ->Draw(drawoption);
    c->cd(padno++); fHistograms->At(khElectronPhi)->Draw(drawoption);
    c->cd(padno++); fHistograms->At(khElectronEta)->Draw(drawoption);
    c->Print(".png");
  }
}

TObject* AliDxHFECorrelation::FindObject(const char *name) const
{
  /// overloaded from TObject: find object by name
  if (fHistograms) {
    return fHistograms->FindObject(name);
  }
  return NULL;
}

TObject* AliDxHFECorrelation::FindObject(const TObject *obj) const
{
  /// overloaded from TObject: find object by pointer
  if (fHistograms) {
    return fHistograms->FindObject(obj);
  }
  return NULL;
}

void     AliDxHFECorrelation::SaveAs(const char *filename,Option_t */*option*/) const
{
  /// overloaded from TObject: save to file
  std::auto_ptr<TFile> output(TFile::Open(filename, "RECREATE"));
  if (!output.get() || output->IsZombie()) {
    AliError(Form("can not open file %s from writing", filename));
    return;
  }
  output->cd();
  if (fHistograms) fHistograms->Write();
  output->Close();
}

AliDxHFECorrelation& AliDxHFECorrelation::operator+=(const AliDxHFECorrelation& other)
{
  /// add histograms from another instance
  if (!fHistograms || !other.fHistograms) return *this;
  
  for (int i=0; i<kNofHistograms; i++) {
    if (fHistograms->At(i)==NULL || other.fHistograms->At(i)==NULL) continue;
    TH1* target=reinterpret_cast<TH1*>(fHistograms->At(i));
    TH1* source=reinterpret_cast<TH1*>(other.fHistograms->At(i));
    if (!target || !source) continue;
    TString name(fHistograms->At(i)->GetName());
    if (name.CompareTo(target->GetName())!=0) {
      AliWarning(Form("skipping incompatible objects at position %d: %s vs %s", i, source->GetName(), target->GetName()));
      continue;
    }
    if (source->IsA()!=target->IsA()) {
      AliWarning(Form("skipping incompatible classes at position %d: %s vs %s", i, source->ClassName(), target->ClassName()));
      continue;
    }
    target->Add(source);
  }
  return *this;
}
