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
#include <TEveSelection.h>
#include <TEveManager.h>

#include <TBuffer3D.h>
#include <TBuffer3DTypes.h>
#include <TVirtualPad.h>
#include <TVirtualViewer3D.h>

#include <TColor.h>
#include <TStyle.h>
#include <TROOT.h>

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
  fEnergyScale(50.0),

  fArrowJetScale (2.0),
	fArrowTrackScale (0.7),

  fGridColor(5),

  fJets(),
  fTracks(),

  fRnrJets (kTRUE),
  fRnrTracks (kTRUE),

  fOneSelection (kTRUE),
  fTwoSelection (kFALSE),
  fSelConnected (kFALSE),

  fJet1(0), fJet2(0), fTrack1(0), fTrack2(0),

  fSelectionFlag (1)
{
  // Constructor.

  SetMainColorPtr(&fGridColor);
  InitMainTrans();
}

AliEveJetPlane::~AliEveJetPlane()
{
  // Destructor.

  if (fSelConnected)
  {
    gEve->GetSelection()->Disconnect("SelectionAdded(TEveElement*)", this);
  }
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

  TEveManager::TRedrawDisabler noRedraw(gEve);

  DestroyElements();

  // Finding the maximum energy
  Double_t eJetMax = 0., eTrackMax = 0., eMax;
  {
    std::vector<AliAODTrack>::iterator k = fTracks.begin();
    std::vector<AliAODJet>::iterator   j = fJets.begin();

    while (j != fJets.end())
    {
      if (j->E() > eJetMax) eJetMax = j->E();
      ++j;
    }

    while (k != fTracks.end())
    {
      if (k->E() > eTrackMax) eTrackMax = k->E();
      ++k;
    }

    eMax = eJetMax > eTrackMax ? eJetMax : eTrackMax;
  }

  // Colors
  Int_t    nCol = gStyle->GetNumberOfColors();

  Double_t eta, phi, e, x, y, h;

  if (fRnrJets)
  {
    UInt_t jetid = 0;
    std::vector<AliAODJet>::iterator j = fJets.begin();
    while (j != fJets.end())
    {
      eta = j->Eta();
      phi = j->Phi();
      e   = j->E();
      h   = TMath::Log(e + 1.) * fEnergyScale;

      x = eta*(fEtaScale);
      y = phi*(fPhiScale) - 350;

      Int_t colBin = TMath::Min((Int_t) ((nCol-2)*TMath::Log(e + 1.)/(TMath::Log(eMax + 1.))),nCol-2);
      Int_t colIdx = gStyle->GetColorPalette(colBin);

      TEveArrow *a = new TEveArrow(0, 0 , h, x, y, 0);
      a->SetSourceObject(&*j);
      a->SetElementName (Form("Jet %d", jetid));
      a->SetElementTitle(Form("Jet 4-momentum: %f, %f, %f, %f \n Pt-Eta-Phi values: %f, %f, %f",
                              j->Px(), j->Py(), j->Pz(), e, j->Pt(), eta, phi ));
      a->SetPickable(kTRUE);
      a->SetMainColor(colIdx);
      a->SetTubeR(0.016*fArrowJetScale);
      a->SetConeR(0.049*fArrowJetScale);
      a->SetConeL(0.170*fArrowJetScale);
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
      h   = TMath::Log(e + 1.) * fEnergyScale;

      if (e < 0)
      {
        Warning("CreateArrows()",
                "Track %d has negative energy - NOT DISPLAYED.", trackid);
	++k; ++trackid;
	continue;
      }

      x = eta*(fEtaScale);
      y = phi*(fPhiScale) - 350;

      Int_t colBin = TMath::Min((Int_t) ((nCol-2)*TMath::Log(e + 1.)/(TMath::Log(eMax + 1.))),nCol-2);
      Int_t colIdx = gStyle->GetColorPalette(colBin);

      TEveArrow *a = new TEveArrow(0, 0 , h, x, y, 0);
      a->SetSourceObject(&*k);
      a->SetElementName (Form("Track %d", trackid));
      a->SetElementTitle(Form("Track 4-momentum: %f, %f, %f, %f \n Pt-Eta-Phi values: %f, %f, %f",
                              k->Px(), k->Py(), k->Pz(), e, k->Pt(), eta, phi ));
      a->SetPickable(kTRUE);
      a->SetMainColor(colIdx);
      a->SetTubeR(0.015*fArrowTrackScale);
      a->SetConeR(0.040*fArrowTrackScale);
      a->SetConeL(0.130*fArrowTrackScale);
      AddElement(a);

      ++k; ++trackid;
    }
  }

  if ( ! fSelConnected)
  {
    gEve->GetSelection()->Connect("SelectionAdded(TEveElement*)",
       "AliEveJetPlane", this, "SelectionAdded(TEveElement*)");
    fSelConnected = kTRUE;
  }
}

/******************************************************************************/

void AliEveJetPlane::SelectionAdded(TEveElement* el)
{
  // Slot called when EVE selection gets a new element.


  if (fOneSelection)
  {
    if (HasChild(el))
    {
      printf("\n\nNOW SELECTED: %s\n", el->GetElementName());

      TObject *src = el->GetSourceObject();

      AliAODTrack *k = dynamic_cast<AliAODTrack*>(src);
      AliAODJet   *j = dynamic_cast<AliAODJet*>  (src);

      //		printf ("Track %p --- Jet %p\n", (void*)t, (void*)j);
      if (k) printf("Track 4-momentum: %f, %f, %f, %f \nPt-Eta-Phi values: %f, %f, %f \n",
                    k->Px(), k->Py(), k->Pz(), k->E(), k->Pt(), k->Eta(), k->Phi() );

      if (j) printf("Jet 4-momentum: %f, %f, %f, %f \nPt-Eta-Phi values: %f, %f, %f \n",
                    j->Px(), j->Py(), j->Pz(), j->E(), j->Pt(), j->Eta(), j->Phi() );
    }
  }

  TEveSelection *sel = gEve->GetSelection();
  if (fTwoSelection && sel->NumChildren() == 2)
  {
    Int_t          numvps   = 0;
    AliVParticle  *vpart[2] = { 0 };
    TEveElement   *elmnt[2] = { 0 };

    for (List_i i = sel->BeginChildren(); i != sel->EndChildren(); ++i)
    {
      if (HasChild(*i))
      {
        TEveElement  *chld = *i;
        TObject      *src  = chld->GetSourceObject();
	AliVParticle *p    = dynamic_cast<AliVParticle*>(src);

        if (p != 0)
        {
          vpart[numvps] = p;
          elmnt[numvps] = chld;
          ++numvps;
          if (numvps >= 2)
            break;
        }
      }
    }

    Double_t eta0, eta1, phi0, phi1, d;

    if (numvps == 2)
    {
      eta0 = vpart[0]->Eta();
      eta1 = vpart[1]->Eta();
      phi0 = vpart[0]->Phi();
      phi1 = vpart[1]->Phi();

      d = TMath::Sqrt(TMath::Power(eta1-eta0,2) + TMath::Power(phi1-phi0,2));
      printf("\n\nNOW SELECTED: '%s' and '%s'\n",
             elmnt[0]->GetElementName(), elmnt[1]->GetElementName());
      printf("Eta-Phi Distance: %f\n", d);
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
