// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveJetPlane.h"

#include <TEveTrans.h>
#include <TEveArrow.h>

#include <TBuffer3D.h>
#include <TBuffer3DTypes.h>
#include <TVirtualPad.h>
#include <TVirtualViewer3D.h>

//______________________________________________________________________________
//
// Show jets and tracks in eta-phi plane.
//
// 

ClassImp(AliEveJetPlane)

Bool_t AliEveJetPlane::fgOneMomentumXYZ      = kFALSE;
Bool_t AliEveJetPlane::fgOneMomentumPhiTheta = kFALSE;
Bool_t AliEveJetPlane::fgOneEta              = kFALSE;
Bool_t AliEveJetPlane::fgOneE                = kFALSE;
Bool_t AliEveJetPlane::fgOneChgMass          = kFALSE;


AliEveJetPlane::AliEveJetPlane(Int_t iev) :
  TEveElementList(Form("AliEveJetPlane %i",iev), Form("%i",iev)),

  fMinEta (-1.5 ),
  fMaxEta ( 1.5 ),
  fMinPhi ( 0.0 ),
  fMaxPhi ( 2.0 * TMath::Pi() ),

  fNEtaDiv(30),
  fNPhiDiv(30),

  fEtaScale(350/1.5),
  fPhiScale(350/(TMath::Pi())),
  fEnergyScale(100.0),

  fEnergyColorScale (0.),

  fGridColor(5),

  fJets(),
  fTracks(),

  fRnrJets (kTRUE),
  fRnrTracks (kTRUE),

  fOneSelection (kTRUE),
  fTwoSelection (kFALSE),

  fJet1(0), fJet2(0), fTrack1(0), fTrack2(0),

  fSelectionFlag (1)
{
  SetMainColorPtr(&fGridColor);
  InitMainTrans();
}

/******************************************************************************/

void AliEveJetPlane::AddJet(AliAODJet* jet) 
{
  // Add a jet for display.
  
  fJets.push_back(*jet);
}

/******************************************************************************/

void AliEveJetPlane::AddTrack(AliAODTrack* track)
{
  // Add a track for display.

  fTracks.push_back(*track);
}

void AliEveJetPlane::CreateArrows()
{
  // Create arrows according to current state.

  DestroyElements();

  Double_t eta, phi, e, x, y;

  if (fRnrJets)
  {
    UInt_t jetid = 0;
    std::vector<AliAODJet>::iterator j = fJets.begin();
    while (j != fJets.end())
    {
      eta = j->Eta();
      phi = j->Phi();
      e   = j->E();
	
      x = eta*(fEtaScale);
      y = phi*(fPhiScale) - 350;
	
      //       printf("Jet 4-momentum: %f, %f, %f, %f \n", v.Px(),v.Py(),v.Pz(),v.Pt() );
      //       printf("Eta-Phi values: %f, %f\n", v.Eta(), v.Phi());	
	
      TEveArrow *a = new TEveArrow(0, 0 , 100, x, y, 0);
      a->SetSourceObject(&*j);
      a->SetElementName (Form("Jet %d", jetid));
      a->SetElementTitle(Form("Jet 4-momentum: %f, %f, %f, %f \n, Eta-Phi values: %f, %f, %f \n",  j->Px(), j->Py(), j->Pz(), j->Pt(), eta, phi, e ));
	
      a->SetPickable(kTRUE);
      a->SetMainColor(kRed);
      a->SetTubeR(0.081);
      a->SetConeR(0.211);
      a->SetConeL(0.481);
      AddElement(a);

      ++j; ++jetid;
    }
  }

  if (fRnrTracks)
  {
    UInt_t trackid = 0;
    std::vector<AliAODTrack>::iterator k = fTracks.begin();  
    while (k != fTracks.end())
    {
      eta = k->Eta();
      phi = k->Phi();
      e   = k->E();
	
      x = eta*(fEtaScale);
      y = phi*(fPhiScale) - 350;
	
      TEveArrow *a = new TEveArrow(0, 0 , 50, x, y, 0);
      a->SetSourceObject(&*k);
      a->SetElementName (Form("Track %d", trackid));
      a->SetElementTitle("Tooltip");
      a->SetPickable(kTRUE);
      a->SetMainColor(kOrange);
      a->SetTubeR(0.061);
      a->SetConeR(0.211);
      a->SetConeL(0.481);
      AddElement(a);

      ++k; ++trackid;
    }
  }
}

/******************************************************************************/

void AliEveJetPlane::ComputeBBox()
{
  // Calculate bounding-box.

  BBoxInit();
  BBoxCheckPoint(-350, -350, -20);
  BBoxCheckPoint( 350, 350,  20);
}


void AliEveJetPlane::Paint(Option_t* /*option*/)
{
  // Paint the object.

  TBuffer3D buff(TBuffer3DTypes::kGeneric);

  // Section kCore
  buff.fID           = this;
  buff.fColor        = GetMainColor();
  buff.fTransparency = GetMainTransparency();
  if (HasMainTrans()) RefMainTrans().SetBuffer3D(buff);
  buff.SetSectionsValid(TBuffer3D::kCore);

  Int_t reqSections = gPad->GetViewer3D()->AddObject(buff);
  if (reqSections == TBuffer3D::kNone) {
    // printf("AliEveJetPlane::Paint viewer was happy with Core buff3d.\n");
    return;
  }
}
