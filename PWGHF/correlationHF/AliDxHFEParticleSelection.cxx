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
{
  // constructor
  // 
  // 
  // 
  // 
}

AliDxHFEParticleSelection::~AliDxHFEParticleSelection()
{
  // destructor
  if (fSelectedTracks) delete fSelectedTracks;
  fSelectedTracks=NULL;
  if (fControlObjects) delete fControlObjects;
  fControlObjects=NULL;
}

int AliDxHFEParticleSelection::InitControlObjects()
{
  /// init the control objects, can be overloaded by childs which should
  /// call AliDxHFEParticleSelection::InitControlObjects() explicitly
  std::auto_ptr<TH1D> hEventControl(new TH1D("hEventControl", "hEventControl", 10, 0, 10));
  std::auto_ptr<TH1D> hTrackControl(new TH1D("hTrackControl", "hTrackControl", 10, 0, 10));

  fhEventControl=hEventControl.release();
  AddControlObject(fhEventControl);
  fhTrackControl=hTrackControl.release();
  AddControlObject(fhTrackControl);

  return 0;
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
  fControlObjects->Add(pObj);
  return 0;
}

int AliDxHFEParticleSelection::HistogramParticleProperties(AliVParticle* p, bool selected)
{
  /// histogram particle properties
  if (!p) return -EINVAL;
  if (!fControlObjects) return 0;

  // TODO: use enums for the bins of the control histogram
  fhTrackControl->Fill(0);
  if (selected) fhTrackControl->Fill(1);
  return 0;
}

TObjArray* AliDxHFEParticleSelection::Select(const AliVEvent* pEvent)
{
  /// create selection, array contains only pointers but does not own the objects
  /// object array needs to be deleted by caller
  if (!pEvent) return NULL;
  TObjArray* selectedTracks=new TObjArray;
  if (!selectedTracks) return NULL;
  int nofTracks=pEvent->GetNumberOfTracks();
  for (int itrack=0; itrack<nofTracks; itrack++) {
    AliVParticle* track=pEvent->GetTrack(itrack);
    bool selected=IsSelected(track);
    HistogramParticleProperties(track, selected);
    if (!selected) continue;
    selectedTracks->Add(track);
  }
  return selectedTracks;
}

TObjArray* AliDxHFEParticleSelection::Select(TObjArray* pTracks)
{
  /// create selection, array contains only pointers but does not own the objects
  /// object array needs to be deleted by caller
  if (!pTracks) return NULL;
  TObjArray* selectedTracks=new TObjArray;
  if (!selectedTracks) return NULL;
  TIter itrack(pTracks);
  TObject* pObj=NULL;
  while ((pObj=itrack())!=NULL) {
    AliVParticle* track=dynamic_cast<AliVParticle*>(pObj);
    if (!track) continue;
    bool selected=IsSelected(track);
    HistogramParticleProperties(track, selected);
    if (!selected) continue;
    selectedTracks->Add(track);
  }
  return selectedTracks;
}

int AliDxHFEParticleSelection::CheckAndAdd(AliVParticle* /*p*/)
{
  /// check and add track to internal array
  /// TODO: check if needed
  return -ENOSYS;
}

bool AliDxHFEParticleSelection::IsSelected(AliVParticle* /*p*/)
{
  /// check particle if it passes the selection criteria
  /// childs can overload, by default all tracks are selected
  return true;
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
 
void AliDxHFEParticleSelection::SaveAs(const char* filename,Option_t */*option*/) const
{
  /// inherited from TObject: save selection criteria
  std::auto_ptr<TFile> output(TFile::Open(filename, "RECREATE"));
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
