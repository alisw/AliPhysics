/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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
/************************************************************************
 *                                                                      *
 * Class for TRD PID                                                    *
 * Implements the abstract base class AliHFEpidbase        *
 * Make PID does the PID decision                                       *
 * Class further contains TRD specific cuts and QA histograms           *
 *                                                                      *
 * Authors:                                                             *
 *   Markus Fasel <M.Fasel@gsi.de>                                      *
 *                                                                      *
 ************************************************************************/
#include <TAxis.h>
#include <TFile.h>
#include <TH1F.h>
#include <TIterator.h>
#include <TKey.h>
#include <TMap.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TString.h>
#include <TROOT.h>

#include "AliESDtrack.h"
#include "AliPID.h"

#include "AliHFEpidTRD.h"

ClassImp(AliHFEpidTRD)

//___________________________________________________________________
AliHFEpidTRD::AliHFEpidTRD(const char* name) :
    AliHFEpidBase(name)
  , fThresholdFile("TRD.PIDthresholds.root")
  , fPIDMethod(kNN)
  , fTRDthresholds(0x0)
  , fTRDelectronEfficiencies(0x0)
{
  //
  // default  constructor
  // 
  fTRDthresholds = new TMap();
}

//___________________________________________________________________
AliHFEpidTRD::AliHFEpidTRD(const AliHFEpidTRD &ref):
    AliHFEpidBase("")
  , fThresholdFile("")
  , fPIDMethod(kLQ)
  , fTRDthresholds(0x0)
  , fTRDelectronEfficiencies(0x0)
{
  //
  // Copy constructor
  //
  ref.Copy(*this);
}

//___________________________________________________________________
AliHFEpidTRD &AliHFEpidTRD::operator=(const AliHFEpidTRD &ref){
  //
  // Assignment operator
  //
  if(this != &ref){
    ref.Copy(*this);
  }
  return *this;
}

//___________________________________________________________________
void AliHFEpidTRD::Copy(TObject &ref) const {
  //
  // Performs the copying of the object
  //
  AliHFEpidTRD &target = dynamic_cast<AliHFEpidTRD &>(ref);

  target.fThresholdFile = fThresholdFile;
  target.fPIDMethod = fPIDMethod;
  if(fTRDthresholds){
    target.fTRDthresholds = dynamic_cast<TMap *>(fTRDthresholds->Clone());
  }
  if(fTRDelectronEfficiencies){
    target.fTRDelectronEfficiencies = new TAxis(*fTRDelectronEfficiencies);
  }
  AliHFEpidBase::Copy(ref);
}

//___________________________________________________________________
AliHFEpidTRD::~AliHFEpidTRD(){
  //
  // Destructor
  //
  if(fTRDthresholds){
    fTRDthresholds->Delete();
    delete fTRDthresholds;
  }
  if(fTRDelectronEfficiencies) delete fTRDelectronEfficiencies;
}

//______________________________________________________
Bool_t AliHFEpidTRD::InitializePID(){
  //
  // InitializePID: Load TRD thresholds and create the electron efficiency axis
  // to navigate 
  //
  LoadTRDthresholds();
  // Fill the electron efficiencies axis for the TRD alone PID
  Int_t nEffs = fTRDthresholds->GetEntries() + 1;
  TIterator* thres =fTRDthresholds->MakeIterator();
  TObjString *key = 0x0;
  Double_t *tmp = new Double_t[nEffs], *titer = tmp;
  while((key = dynamic_cast<TObjString *>((*thres)()))){
    (*titer++) = static_cast<Double_t>(key->String().Atoi())/100.;
  }
  delete thres;
  *titer = 1.;
  // Sort the electron efficiencies and put them into the TAxis for later navigation
  Int_t *ind = new Int_t[nEffs], *iiter = ind;
  TMath::Sort(nEffs, tmp, ind, kFALSE);
  Double_t *eleffs = new Double_t[nEffs], *eiter = eleffs;
  while(eiter < &eleffs[nEffs]) *(eiter++) = tmp[*(iiter++)];
  // print the content
  Int_t cnt = 0;
  if(GetDebugLevel() > 1){
    printf("Printing electron efficiency bins:\n");
    eiter = eleffs;
    thres = fTRDthresholds->MakeIterator();
    TObject *object = 0x0;
    while(eiter < &eleffs[nEffs - 1]){
      printf("eleffs[%d] = %f", cnt++, *(eiter++));
      key = dynamic_cast<TObjString *>(thres->Next());
      object = fTRDthresholds->GetValue(key->String().Data());
      printf(", Content: %p\n", (void *)object);
    }
  }
  delete[] tmp; delete[] ind;
  fTRDelectronEfficiencies = new TAxis(nEffs - 1, eleffs);
  if(GetDebugLevel() > 1){
    printf("Printing axis content:\n");
    for(Int_t ibin = fTRDelectronEfficiencies->GetFirst(); ibin <= fTRDelectronEfficiencies->GetLast(); ibin++)
      printf("%d.) minimum: %f, maximum? %f\n", ibin, fTRDelectronEfficiencies->GetBinLowEdge(ibin), fTRDelectronEfficiencies->GetBinUpEdge(ibin));
  }
  delete[] eleffs;
  return kTRUE;
}

//______________________________________________________
Int_t AliHFEpidTRD::IsSelected(AliVParticle *track){
  //
  // Does PID for TRD alone:
  // PID thresholds based on 90% Electron Efficiency level approximated by a linear 
  // step function
  //
  AliESDtrack *esdTrack = 0x0;
  if(!(esdTrack = dynamic_cast<AliESDtrack *>(track))) return kFALSE;
  Double_t p = esdTrack->GetOuterParam() ? esdTrack->GetOuterParam()->P() : esdTrack->P();
  if(p < 0.6) return 0;

  // Get the Histograms
  TH1 *threshist = GetTRDthresholds(0.91);
  Int_t bin = 0;
  if(p > threshist->GetXaxis()->GetXmax()) 
    bin = threshist->GetXaxis()->GetLast();
  else if(p < threshist->GetXaxis()->GetXmin()) 
    bin = threshist->GetXaxis()->GetFirst();
  else
    bin = threshist->GetXaxis()->FindBin(p);

  Double_t pidProbs[AliPID::kSPECIES];
  esdTrack->GetTRDpid(pidProbs);
  if(pidProbs[AliPID::kElectron] > threshist->GetBinContent(bin)) return 11;
  return 0;
}

//___________________________________________________________________
void AliHFEpidTRD::LoadTRDthresholds(){
  //
  // Load TRD threshold histograms from File
  //
  TFile *mythresholds = TFile::Open(fThresholdFile);
  TKey *object = 0x0;
  TString electron_eff;
  TObjArray *histos = 0x0;
  Float_t eff;
  TH1F *refhist = 0x0;
  TIterator *keyIterator = mythresholds->GetListOfKeys()->MakeIterator();
  TString histnames[2] = {"fHistThreshLQ", "fHistThreshNN"};
  gROOT->cd();
  while((object = dynamic_cast<TKey *>((*keyIterator)()))){
    // Get the electron efficiency bin this histogram was taken with
    electron_eff = object->GetName();
    electron_eff = electron_eff.Remove(0,3);
    eff = static_cast<Float_t>(electron_eff.Atoi())/100.;

    // Get the threshold according to the selected 
    histos = dynamic_cast<TObjArray *>(object->ReadObj());
    refhist = dynamic_cast<TH1F *>(histos->FindObject(histnames[fPIDMethod].Data()));
    SetTRDthresholds(refhist, eff);
    histos->Delete();
    delete histos;
  }
  delete keyIterator;
  mythresholds->Close();
  delete mythresholds;
}

//___________________________________________________________________
void AliHFEpidTRD::SetTRDthresholds(TH1F *thresholds, Float_t electronEff){
  //
  // Set the threshold histogram for the TRD pid
  //
  fTRDthresholds->Add(new TObjString(Form("%d", TMath::Nint(electronEff * 100.))), new TH1F(*thresholds));
}

//___________________________________________________________________
TH1F *AliHFEpidTRD::GetTRDthresholds(Float_t electronEff){ 
  Int_t bin = 0;
  if(electronEff < fTRDelectronEfficiencies->GetXmin()) bin = fTRDelectronEfficiencies->GetFirst();
  else if(electronEff > fTRDelectronEfficiencies->GetXmax()) bin = fTRDelectronEfficiencies->GetLast();
  else bin = fTRDelectronEfficiencies->FindBin(electronEff);
 TObjString keyname = Form("%d", TMath::Nint(fTRDelectronEfficiencies->GetBinLowEdge(bin)* 100.));
/*  printf("Key: %s\n", keyname.String().Data());*/
  TH1F *thresholds = dynamic_cast<TH1F *>((dynamic_cast<TPair *>(fTRDthresholds->FindObject(&keyname)))->Value());
/*  printf("thresholds: %p\n", thresholds);*/
  return thresholds;
}
