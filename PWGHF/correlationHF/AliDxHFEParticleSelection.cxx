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

/// @file   AliDxHFEParticleSelection.cxx
/// @author Sedat Altinpinar, Hege Erdal, Matthias Richter
/// @date   2012-03-19
/// @brief  Base class for particle selection
///

#include "AliDxHFEParticleSelection.h"
#include "AliLog.h"
#include "AliVEvent.h"
#include "AliVParticle.h"
#include "TObjArray.h"
#include "TList.h"
#include "TMath.h"
#include "TH1D.h"
#include "THnSparse.h"
#include "AliReducedParticle.h"
#include "TFile.h"
#include <iostream>
#include <cerrno>
#include <memory>

using namespace std;

/// ROOT macro for the implementation of ROOT specific class methods
ClassImp(AliDxHFEParticleSelection)

AliDxHFEParticleSelection::AliDxHFEParticleSelection(const char* name, const char* opt)
  : TNamed(name?name:"AliDxHFEParticleSelection", name?name:"AliDxHFEParticleSelection")
  , fOption(opt)
  , fSelectedTracks(NULL)
  , fControlObjects(NULL)
  , fhEventControl(NULL)
  , fhTrackControl(NULL)
  , fUseMC(0)
  , fVerbosity(0)
  , fDimThn(-1)
  , fParticleProperties(NULL)
{
  // constructor
  // 
  // 
  // 
  // 
}

const char* AliDxHFEParticleSelection::fgkEventControlBinNames[]={
  "nEventsAll",
  "nEventsSelected",
  "nEventsWParticle"
};

const char* AliDxHFEParticleSelection::fgkTrackControlBinNames[]={
  "nTrackAll",
  "nTrackSelected",
};

AliDxHFEParticleSelection::~AliDxHFEParticleSelection()
{
  // destructor
  if (fParticleProperties) delete[] fParticleProperties;
  fParticleProperties=NULL;
  if (fSelectedTracks) delete fSelectedTracks;
  fSelectedTracks=NULL;
  if (fControlObjects) delete fControlObjects;
  fControlObjects=NULL;
  fhEventControl=NULL;
  fhTrackControl=NULL;
}

int AliDxHFEParticleSelection::Init()
{
  //
  // Init part sel. Calls InitControlObjects()
  //

  InitControlObjects();
  return 0;
}


int AliDxHFEParticleSelection::InitControlObjects()
{
  //
  // init control objects
  // TODO: Change to private now that have Init()?
  //
  if (fVerbosity>0) {
    AliInfo("Setting up control objects");
  }

  /// init the control objects, can be overloaded by childs which should
  /// call AliDxHFEParticleSelection::InitControlObjects() explicitly
  std::auto_ptr<TH1D> hEventControl(new TH1D("hEventControl", "hEventControl", 10, 0, 10));
  std::auto_ptr<TH1D> hTrackControl(new TH1D("hTrackControl", "hTrackControl", 10, 0, 10));

  fhEventControl=hEventControl.release();
  for (int iLabel=0; iLabel<kNEventPropertyLabels; iLabel++)
    fhEventControl->GetXaxis()->SetBinLabel(iLabel+1, fgkEventControlBinNames[iLabel]);
  AddControlObject(fhEventControl);
  fhTrackControl=hTrackControl.release();
  for (int iLabel=0; iLabel<kNTrackPropertyLabels; iLabel++)
    fhTrackControl->GetXaxis()->SetBinLabel(iLabel+1, fgkTrackControlBinNames[iLabel]);
  AddControlObject(fhTrackControl);

  return 0;
}

THnSparse* AliDxHFEParticleSelection::CreateControlTHnSparse(const char* name,
							     int thnSize,
							     int* thnBins,
							     double* thnMin,
							     double* thnMax,
							     const char** binLabels) const
{
  //
  // Creates THnSparse. 
  //

  AliInfo("Setting up THnSparse");

  std::auto_ptr<THnSparseF> th(new THnSparseF(name, name, thnSize, thnBins, thnMin, thnMax));
  if (th.get()==NULL) {
    return NULL;
  }
  for (int iLabel=0; iLabel<thnSize; iLabel++) {
    th->GetAxis(iLabel)->SetTitle(binLabels[iLabel]);    
   
  }
 return th.release();

}

THnSparse* AliDxHFEParticleSelection::DefineTHnSparse() 
{
  //
  // Defines the THnSparse. For now, only calls CreatControlTHnSparse

  //TODO: Should make it more general. Or maybe one can use this here, and skip in PartSelEl?

  // here is the only place to change the dimension
  const int thnSize = 3;
  InitTHnSparseArray(thnSize);

  const double Pi=TMath::Pi();
  TString name;
  //     		            0    1       2
  // 	 	                    Pt   Phi    Eta
  int         thnBins [thnSize] = { 1000,  200, 500};
  double      thnMin  [thnSize] = {    0,    0, -1.};
  double      thnMax  [thnSize] = {  100, 2*Pi,  1.};
  const char* thnNames[thnSize] = { "Pt","Phi","Eta"};

  name.Form("%s info", GetName());

  return CreateControlTHnSparse(name,thnSize,thnBins,thnMin,thnMax,thnNames);
}

int AliDxHFEParticleSelection::AddControlObject(TObject* pObj)
{
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
  if (GetVerbosity()>0) {
    AliInfo(Form("Adding object '%s' of type %s",pObj->GetName(),pObj->ClassName()));
  }
  fControlObjects->Add(pObj);
  return 0;
}

int AliDxHFEParticleSelection::HistogramEventProperties(int bin)
{
  /// histogram event properties
  if (!fControlObjects) return 0;

  // TODO: use enums for the bins of the control histogram
  // for now: 0=all, 1=events with D0s, 2=events with correlated D0s
  fhEventControl->Fill(bin);
  return 0;
}

int AliDxHFEParticleSelection::HistogramParticleProperties(AliVParticle* p, int selected)
{
  /// histogram particle properties

  if (!p) return -EINVAL;
  if (!fControlObjects) return 0;

  // TODO: use enums for the bins of the control histogram
  fhTrackControl->Fill(kTrackAll);
  if (selected) fhTrackControl->Fill(kTrackSel);
  return 0;
}

int AliDxHFEParticleSelection::FillParticleProperties(AliVParticle* p, Double_t* data, int dimension) const
{
  // fill the data array from the particle data
  if (!data) return -EINVAL;
  int i=0;
  if (dimension!=GetDimTHnSparse()) {
    // TODO: think about filling only the available data and throwing a warning
    return -ENOSPC;
  }
  data[i++]=p->Pt();
  data[i++]=p->Phi();
  data[i++]=p->Eta();
  return i;
}

TObjArray* AliDxHFEParticleSelection::Select(const AliVEvent* pEvent)
{
  /// create selection from 'Tracks' member of the event,
  /// array contains only pointers but does not own the objects
  /// object array needs to be deleted by caller
  if (!pEvent) return NULL;
  TObjArray* selectedTracks=new TObjArray;
  if (!selectedTracks) return NULL;
  selectedTracks->SetOwner(kFALSE); // creating new track objects below
  int nofTracks=pEvent->GetNumberOfTracks();
  for (int itrack=0; itrack<nofTracks; itrack++) {
    AliVParticle* track=pEvent->GetTrack(itrack);
    int selectionCode=IsSelected(track,pEvent);
    HistogramParticleProperties(track, selectionCode);
    if (selectionCode==0) continue;
    selectedTracks->Add(CreateParticle(track));
  }
  return selectedTracks;
}

TObjArray* AliDxHFEParticleSelection::Select(TObjArray* pParticles, const AliVEvent* pEvent)
{
  /// create selection from the array of particles,
  /// array contains only pointers but does not own the objects
  /// object array needs to be deleted by caller
  if (!pParticles) return NULL;
  TObjArray* selectedTracks=new TObjArray;
  if (!selectedTracks) return NULL;
  selectedTracks->SetOwner(); // creating new track objects below
  TIter next(pParticles);
  TObject* pObj=NULL;
  while ((pObj=next())) {
    AliVParticle* track=dynamic_cast<AliVParticle*>(pObj);
    if (!track) continue;
    int selectionCode=IsSelected(track, pEvent);
    HistogramParticleProperties(track, selectionCode);
    if (selectionCode ==0) continue;
    selectedTracks->Add(CreateParticle(track));
  }
  return selectedTracks;
}

int AliDxHFEParticleSelection::CheckAndAdd(AliVParticle* /*p*/)
{
  /// check and add track to internal array
  /// TODO: check if needed
  return -ENOSYS;
}

int AliDxHFEParticleSelection::IsSelected(AliVParticle* /*p*/, const AliVEvent* /*e*/)
{
  /// check particle if it passes the selection criteria
  /// childs can overload, by default all tracks are selected
  return 1;
}

void AliDxHFEParticleSelection::AliDxHFEParticleSelection::Clear(Option_t * /*option*/)
{
  /// inherited from TObject: cleanup
}

void AliDxHFEParticleSelection::Print(Option_t */*option*/) const
{
  /// inherited from TObject: print info
  cout << "====================================================================" << endl;
  TNamed::Print();
  if (fControlObjects) fControlObjects->Print();
}
 
void AliDxHFEParticleSelection::SaveAs(const char* filename, Option_t */*option*/) const
{
  /// inherited from TObject: save selection criteria
  TString fileoption;
  // TODO: options recreate
  fileoption="RECREATE";
  //else fileoption="UPDATE";

  std::auto_ptr<TFile> output(TFile::Open(filename,fileoption));
  if (!output.get() || output->IsZombie()) {
    AliError(Form("can not open file %s from writing", filename));
    return;
  }
  output->cd();
  if (fControlObjects) fControlObjects->Write();
  output->Close();
}

void AliDxHFEParticleSelection::Draw(Option_t* /*option*/)
{
  /// inherited from TObject: draw content

  // TODO: implement drawing code
  // - create canvas objects
  // - plot internal objects
  // - optionally save canvases to file
  //
  // It might be appropriate to have another Draw function taking a
  // TList as argument and implementing the actual drawing. If this
  // function is 'static', it can be used stand-alone also from macros
}

TObject* AliDxHFEParticleSelection::FindObject(const char* name) const
{
  /// inherited from TObject: find object by name

  if (fControlObjects) {
    return fControlObjects->FindObject(name);
  }
  return NULL;
}

TObject* AliDxHFEParticleSelection::FindObject(const TObject* obj) const
{
  /// inherited from TObject: find object by pointer
  if (fControlObjects) {
    return fControlObjects->FindObject(obj);
  }
  return NULL;
}

AliVParticle *AliDxHFEParticleSelection::CreateParticle(AliVParticle* track)
{
  // Creating object with reduced particle properties
  AliReducedParticle *part = new AliReducedParticle(track->Eta(), track->Phi(), track->Pt(), track->Charge(), 0);

  return part;

}
