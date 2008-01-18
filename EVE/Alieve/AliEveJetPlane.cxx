// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveJetPlane.h"
#include <TString.h>
#include <TBuffer3D.h>
#include <TBuffer3DTypes.h>
#include <TVirtualPad.h>
#include <TVirtualViewer3D.h>


//______________________________________________________________________
// AliEveJetPlane
//


Bool_t AliEveJetPlane::fgOneMomentumXYZ      = kFALSE;
Bool_t AliEveJetPlane::fgOneMomentumPhiTheta = kFALSE;
Bool_t AliEveJetPlane::fgOneEta     = kFALSE;
Bool_t AliEveJetPlane::fgOneE       = kFALSE;
Bool_t AliEveJetPlane::fgOneChgMass = kFALSE;


ClassImp(AliEveJetPlane)

AliEveJetPlane::AliEveJetPlane(Int_t iev) :
  TEveElementList(Form("AliEveJetPlane %i",iev), Form("%i",iev)),

  fMinEta (-1.5 ),
  fMaxEta ( 1.5 ),
  fMinPhi (-TMath::Pi() ),
  fMaxPhi ( TMath::Pi() ),

  fNEtaDiv(30),
  fNPhiDiv(30),

  fEtaScale(350/1.5),
  fPhiScale(350/TMath::Pi()),
  fEnergyScale(100.0),

  fEnergyColorScale (0.),

  fGridColor(5),

  fRnrJets (kTRUE),
  fRnrTracks (kTRUE),

  fOneSelection (kTRUE),
  fTwoSelection (kFALSE),

  fSelectionFlag (1)
{
  SetMainColorPtr(&fGridColor);
}

/**************************************************************************/

void AliEveJetPlane::AddJet(AliAODJet jet)
{
  fJets.push_back(jet);
}

/**************************************************************************/

void AliEveJetPlane::AddTrack(AliAODTrack track)
{
  fTracks.push_back(track);
}


/**************************************************************************/

void AliEveJetPlane::ComputeBBox()
{
  BBoxInit();
  BBoxCheckPoint(-350, -350, -20);
  BBoxCheckPoint( 350,  350,  20);
}

void AliEveJetPlane::Paint(Option_t* /*option*/)
{
  TBuffer3D buff(TBuffer3DTypes::kGeneric);

  // Section kCore
  buff.fID           = this;
  buff.fColor        = fGridColor;
  buff.fTransparency = 0;
  fHMTrans.SetBuffer3D(buff);
  buff.SetSectionsValid(TBuffer3D::kCore);

  Int_t reqSections = gPad->GetViewer3D()->AddObject(buff);
  if (reqSections == TBuffer3D::kNone) {
    // printf("AliEveJetPlane::Paint viewer was happy with Core buff3d.\n");
    return;
  }
}
