// $Header$

#include "JetPlane.h"
#include <TString.h>
#include <TBuffer3D.h>
#include <TBuffer3DTypes.h>
#include <TVirtualPad.h>
#include <TVirtualViewer3D.h>

using namespace Alieve;

//______________________________________________________________________
// JetPlane
//


Bool_t JetPlane::fgOneMomentumXYZ      = kFALSE;
Bool_t JetPlane::fgOneMomentumPhiTheta = kFALSE;
Bool_t JetPlane::fgOneEta     = kFALSE;
Bool_t JetPlane::fgOneE       = kFALSE;
Bool_t JetPlane::fgOneChgMass = kFALSE;


ClassImp(JetPlane)

JetPlane::JetPlane(Int_t iev) :
  TEveElementList(Form("JetPlane %i",iev), Form("%i",iev)),

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

void JetPlane::AddJet(AliAODJet jet)
{
  fJets.push_back(jet);
}

/**************************************************************************/

void JetPlane::AddTrack(AliAODTrack track)
{
  fTracks.push_back(track);
}


/**************************************************************************/

void JetPlane::ComputeBBox()
{
  BBoxInit();
  BBoxCheckPoint(-350, -350, -20);
  BBoxCheckPoint( 350,  350,  20);
}

void JetPlane::Paint(Option_t* /*option*/)
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
    // printf("JetPlane::Paint viewer was happy with Core buff3d.\n");
    return;
  }
}
