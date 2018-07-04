/*
 * AliFemtoPairCleanerHists.cxx
 *
 *  Created on: 22 Dec 2017
 *      Author: bernhardhohlweger
 */

#include "AliFemtoDreamPairCleanerHists.h"
ClassImp(AliFemtoDreamPairCleanerHists)
AliFemtoDreamPairCleanerHists::AliFemtoDreamPairCleanerHists()
:fTrackDecays(nullptr)
,fDecayDecays(nullptr)
,fPairInvMass(nullptr)
,fOutput(0)
{
}

AliFemtoDreamPairCleanerHists::AliFemtoDreamPairCleanerHists(
    int nTrackDecays,int nDecayDecays,int nInvMassPairs)
{
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("PairCleaner");

  fTrackDecays=new TH1F*[nTrackDecays];
  for (int i=0;i<nTrackDecays;++i) {
    TString histName=Form("DaugthersSharedTracks_%d",i);
    fTrackDecays[i]=new TH1F(histName.Data(),histName.Data(),20,0,20);
    fOutput->Add(fTrackDecays[i]);
  }
  fDecayDecays=new TH1F*[nDecayDecays];
  for (int i=0;i<nTrackDecays;++i) {
    TString histName=Form("DaugthersSharedDaughters_%d",i);
    fDecayDecays[i]=new TH1F(histName.Data(),histName.Data(),20,0,20);
    fOutput->Add(fDecayDecays[i]);
  }
  fPairInvMass=new TH1F*[nInvMassPairs];
  for (int i=0;i<nInvMassPairs;++i) {
    TString histName=Form("InvMassPair_%d",i);
    //this is tuned to look for the H Dibaryon, if neccessary setters need to be
    //introduced.
    fPairInvMass[i]=new TH1F(histName.Data(),histName.Data(),1000,2.25,2.5);
    fPairInvMass[i]->Sumw2();
    fOutput->Add(fPairInvMass[i]);
  }
}

AliFemtoDreamPairCleanerHists::~AliFemtoDreamPairCleanerHists() {
  // TODO Auto-generated destructor stub
}

