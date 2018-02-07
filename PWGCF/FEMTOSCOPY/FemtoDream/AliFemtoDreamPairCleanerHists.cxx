/*
 * AliFemtoPairCleanerHists.cxx
 *
 *  Created on: 22 Dec 2017
 *      Author: bernhardhohlweger
 */

#include "AliFemtoDreamPairCleanerHists.h"
ClassImp(AliFemtoDreamPairCleanerHists)
AliFemtoDreamPairCleanerHists::AliFemtoDreamPairCleanerHists()
:fTrackDecays(0)
,fDecayDecays(0)
,fOutput(0)
{
}

AliFemtoDreamPairCleanerHists::AliFemtoDreamPairCleanerHists(
    int nTrackDecays,int nDecayDecays)
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
}

AliFemtoDreamPairCleanerHists::~AliFemtoDreamPairCleanerHists() {
  // TODO Auto-generated destructor stub
}

