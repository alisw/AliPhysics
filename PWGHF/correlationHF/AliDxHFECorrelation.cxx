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
//#include "AliAnalysisCuts.h"         // required dependency libANALYSISalice.so
//#include "AliFlowTrackSimple.h"      // required dependency libPWGflowBase.so
//#include "AliFlowCandidateTrack.h"   // required dependency libPWGflowTasks.so
//#include "AliCFContainer.h"          // required dependency libCORRFW.so
#include "AliAODRecoDecayHF2Prong.h"   // libPWGHFvertexingHF
#include "AliRDHFCutsD0toKpi.h"
#include "TObjArray.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THnSparse.h"
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include <iostream>
#include <cerrno>
#include <memory>

using namespace std;

ClassImp(AliDxHFECorrelation)

AliDxHFECorrelation::AliDxHFECorrelation(const char* name)
  : TNamed(name?name:"AliDxHFECorrelation", "")
  , fHistograms(NULL)  
  , fControlObjects(NULL)
  , fCorrProperties(NULL)
  , fhEventControlCorr(NULL)
  , fCuts(NULL)
  , fUseMC(kFALSE)
{
  // default constructor
  // 
  //

}

const char* AliDxHFECorrelation::fgkEventControlBinNames[]={
  "nEventsAll",
  "nEventsSelected",
  "nEventsD0",
  "nEventsD0e"
};

// TODO: maybe delete PtD0 in the future, note: there are more places!
const char* AliDxHFECorrelation::fgkCorrControlBinNames[]={
  "D0InvMass",
  "PtD0",
  "PhiD0",
  "PtBinD0",
  "dPhi",
  "Pt electron"
};

AliDxHFECorrelation::~AliDxHFECorrelation()
{
  // destructor
  //
  //
  if (fHistograms) delete fHistograms;
  fHistograms=NULL;

  // NOTE: fControlObjects owns the object, and they are deleted in the
  // destructor of TList
  if (fControlObjects) delete fControlObjects;
  fControlObjects=NULL;
  fCorrProperties=NULL;
  fhEventControlCorr=NULL;

  // NOTE: the external object is deleted elsewhere
  fCuts=NULL;
}

int AliDxHFECorrelation::Init()
{
  AliInfo("Setting up THnSparse");
  // Are using THnSparse instead of histograms to store information
  // At the moment use old setup for THnSparse, change to new later

  /*
    TString name;
    const int thnSize = 7;
    const double pi=TMath::Pi();
    // 			         0      1     2      3      4      5    6
    // 			      mass      Pt   Phi    ePt    ePhi  eEta  DeltaPhi
    int    thnBins[thnSize] = {   200,   1000,  100,   100,   100,  100,  180  };
    double thnMin [thnSize] = {  1.5648,   0,    0,     0,    0.0, -1.0,  0.0  };
    double thnMax [thnSize] = {  2.1648, 100,  (2*pi), 100, (2*pi), 1.0, (2*pi)};

    name.Form("%s info", GetName());
    std::auto_ptr<THnSparseF> CorrProperties(new THnSparseF(name, name, thnSize, thnBins, thnMin, thnMax));

    if (CorrProperties.get()==NULL) {
    return -ENOMEM;
    }
    int axis=0;
    CorrProperties->GetAxis(axis++)->SetTitle("D0 Inv Mass");
    CorrProperties->GetAxis(axis++)->SetTitle("D0 Pt");
    CorrProperties->GetAxis(axis++)->SetTitle("D0 Phi"); 
    CorrProperties->GetAxis(axis++)->SetTitle("electron Pt"); 
    CorrProperties->GetAxis(axis++)->SetTitle("electron Phi"); 
    CorrProperties->GetAxis(axis++)->SetTitle("electron Eta"); 
    CorrProperties->GetAxis(axis++)->SetTitle("#Delta#Phi D0 -HFE"); */


  // TODO: think about removing PtD0, but remember to change this in all places!
  static const int sizeEventdphi = 6;  
  Double_t minPhi= -TMath::Pi()/2;
  Double_t maxPhi= 3*TMath::Pi()/2;
  // 			                    0         1      2      3       4       5
  // 			                  D0invmass PtD0   PhiD0  PtbinD0  dphi   Pte
  int    binsEventdphi[sizeEventdphi] = {   200,    1000,    100,     21,  100,    1000};
  double minEventdphi [sizeEventdphi] = { 1.5648,    0,      0,      0,    minPhi,   0 };
  double maxEventdphi [sizeEventdphi] = {  2.1648,  100, 2*(TMath::Pi()), 20,  maxPhi, 100 };

  TString name;
  name.Form("%s info", GetName());
  std::auto_ptr<THnSparseF> CorrProperties(new THnSparseF(name, name, sizeEventdphi, binsEventdphi,minEventdphi , maxEventdphi));
  if (CorrProperties.get()==NULL) {
    return -ENOMEM;
  }
  int iLabel=0;

  for (iLabel=0; iLabel<sizeEventdphi; iLabel++)
    CorrProperties->GetAxis(iLabel)->SetTitle(fgkCorrControlBinNames[iLabel]);

  //----------------------------------------------
  // Histogram for storing event information

  std::auto_ptr<TH1D> hEventControl(new TH1D("hEventControlCorr", "hEventControlCorr", 10, 0, 10));
  if (!hEventControl.get()) {
    return -ENOMEM;
  }

  for (iLabel=0; iLabel<kNEventControlLabels; iLabel++)
    hEventControl->GetXaxis()->SetBinLabel(iLabel+1, fgkEventControlBinNames[iLabel]);

  fCorrProperties=CorrProperties.release();
  AddControlObject(fCorrProperties);
  fhEventControlCorr=hEventControl.release();
  AddControlObject(fhEventControlCorr);

  return 0;
}

int AliDxHFECorrelation::AddControlObject(TObject* pObj)
{
  AliInfo("Adding object");
  /// add control object to list, the base class becomes owner of the object
  if (!pObj) return -EINVAL;
  if (!fControlObjects) {
    fControlObjects=new TList;
    if (!fControlObjects) return -ENOMEM;
    fControlObjects->SetOwner();
  }
  if (fControlObjects->FindObject(pObj->GetName())) {
    AliError(Form("ignoring duplicate object '%s' of type %s", pObj->GetName(), pObj->ClassName()));
    return -EEXIST;
  }
  fControlObjects->Add(pObj);
  return 0;
}

int AliDxHFECorrelation::HistogramEventProperties(int bin)
{
  /// histogram event properties
  if (!fhEventControlCorr) return 0;
  fhEventControlCorr->Fill(bin);

  return 0;
}

int AliDxHFECorrelation::Fill(const TObjArray* triggerCandidates, const TObjArray* associatedTracks)
{
  /// fill ThnSparse from array of AliVParticle objects
  if (!triggerCandidates || !associatedTracks) return -EINVAL;
  if (!fControlObjects) {
    Init();
  }
  if (!fControlObjects) {
    AliError("Initialisation failed, can not fill THnSparse");
  }

  const double Pii=TMath::Pi();

  TIter itrigger(triggerCandidates);
  TObject* otrigger=NULL;
  int ctrigger=-1;

  // For the moment this is very specific to D0-electron correlation. Should be 
  // changed to be more specific. 
  while ((otrigger=itrigger())!=NULL) {
    // loop over trigger D0 particle
    ctrigger++;
    AliAODRecoDecayHF2Prong *d0 = dynamic_cast<AliAODRecoDecayHF2Prong*>(otrigger);
    if (!d0) continue;

    TIter iElectron(associatedTracks);
    TObject* oElectron=NULL;
    int cElectron=-1;
    while ((oElectron=iElectron())!=NULL) {
      // loop over electrons

      cElectron++;
      AliVParticle* pElectron=dynamic_cast<AliVParticle*>(oElectron);
      if (!pElectron) continue;

      //Calculating dPhi using TLorentzVectors
      Double_t mPDG=TDatabasePDG::Instance()->GetParticle(421)->Mass();
      TLorentzVector D0vector(0.,0.,0.,0.);
      TLorentzVector evector(0.,0.,0.,0.);
      D0vector.SetXYZM(d0->Px(),d0->Py(),d0->Pz(),mPDG);  
      evector.SetXYZM(pElectron->Px(),pElectron->Py(),pElectron->Pz(),0.000511);
      Double_t DeltaPhi=D0vector.DeltaPhi(evector);
      if(DeltaPhi<-TMath::PiOver2()) DeltaPhi=DeltaPhi+(2*Pii);

      /*Double_t CorrelationArray[]={d0->InvMassD0(),d0->Pt(),d0->Phi(),
	pElectron->Pt(),pElectron->Phi(), pElectron->Eta(),
	DeltaPhi};*/

      // TODO: think about a method to retrieve the pt bin from the
      // selection object
      Int_t ptbin=fCuts->PtBin(d0->Pt());
      Double_t CorrelationArray[]={d0->InvMassD0(),d0->Pt(),d0->Phi(),ptbin,
				   DeltaPhi, pElectron->Pt() };
      fCorrProperties->Fill(CorrelationArray);

    } // loop over associated tracks
  } // loop over trigger particle

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
}

TObject* AliDxHFECorrelation::FindObject(const char *name) const
{
  /// overloaded from TObject: find object by name
  if (fControlObjects) {
    return fControlObjects->FindObject(name);
  }
  return NULL;
}

TObject* AliDxHFECorrelation::FindObject(const TObject *obj) const
{
  /// overloaded from TObject: find object by pointer
  if (fControlObjects) {
    return fControlObjects->FindObject(obj);
  }
  return NULL;
}

void AliDxHFECorrelation::SaveAs(const char *filename, Option_t */*option*/) const
{
  /// overloaded from TObject: save to file
  std::auto_ptr<TFile> output(TFile::Open(filename, "RECREATE"));
  if (!output.get() || output->IsZombie()) {
    AliError(Form("can not open file %s from writing", filename));
    return;
  }
  output->cd();
  if (fControlObjects) fControlObjects->Write();
  output->Close();
}

AliDxHFECorrelation& AliDxHFECorrelation::operator+=(const AliDxHFECorrelation& other)
{
  /// add histograms from another instance
  // TODO - need to change this to ThnSparse?
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
