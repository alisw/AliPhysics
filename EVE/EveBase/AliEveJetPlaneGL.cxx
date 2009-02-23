// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#include "AliEveJetPlaneGL.h"
#include "AliEveJetPlane.h"

#include <TGLRnrCtx.h>
#include <TGLSelectRecord.h>
#include <TGLIncludes.h>

#include "TGLUtil.h"
#include "TGLAxis.h"

#include <TColor.h>
#include <TStyle.h>
#include <TROOT.h>

//==============================================================================
//==============================================================================
// AliEveJetPlaneGL
//==============================================================================

//______________________________________________________________________________
//
// GL renderer for AliEveJetPlane.

ClassImp(AliEveJetPlaneGL)

AliEveJetPlaneGL::AliEveJetPlaneGL() : TGLObject(), fM(0)
{
  // Constructor.

  fDLCache = kFALSE; // Disable display list -- axis pain.
}

/******************************************************************************/

Bool_t AliEveJetPlaneGL::SetModel(TObject* obj, const Option_t* /*opt*/)
{
  // Set model object.

  if (SetModelCheckClass(obj, AliEveJetPlane::Class())) {
    fM = dynamic_cast<AliEveJetPlane*>(obj);
    return kTRUE;
  }
  return kFALSE;
}

void AliEveJetPlaneGL::SetBBox()
{
  // Set bounding box.

  SetAxisAlignedBBox(((AliEveJetPlane*)fExternalObj)->AssertBBox());
}

/******************************************************************************/

void AliEveJetPlaneGL::DirectDraw(TGLRnrCtx& rnrCtx) const
{
  // Render the object.

  Float_t minEta = (fM->fMinEta)*(fM->fEtaScale);
  Float_t maxEta = (fM->fMaxEta)*(fM->fEtaScale);
  Float_t minPhi = (fM->fMinPhi)*(fM->fPhiScale) - 350;
  Float_t maxPhi = (fM->fMaxPhi)*(fM->fPhiScale) - 350;
  Float_t phiCoord, etaCoord, dPhi, dEta;
  Double_t eta, phi, e, x, y;

  // Show frame for Eta-Phi coordinates

  glBegin(GL_LINE_LOOP);
  glVertex3f(minEta, minPhi, 0);
  glVertex3f(maxEta, minPhi, 0);
  glVertex3f(maxEta, maxPhi, 0);
  glVertex3f(minEta, maxPhi, 0);
  glEnd();

  if (rnrCtx.Selection() == kFALSE && rnrCtx.Highlight() == kFALSE)
  {

    // Show grid in Eta-Phi coordinates

    dPhi = (maxPhi-minPhi)/(fM->fNPhiDiv - 1);
    dEta = (maxEta-minEta)/(fM->fNEtaDiv - 1);

    for (Int_t count = 1; count < fM->fNPhiDiv-1; ++count)
    {
      phiCoord = minPhi + count*dPhi;
      glBegin(GL_LINES);
      glVertex3f( minEta, phiCoord, 0);
      glVertex3f( maxEta, phiCoord, 0);
      glEnd();
    }

    for (Int_t count = 1; count < fM->fNEtaDiv-1; ++count)
    {
      etaCoord = minEta + count*dEta;
      glBegin(GL_LINES);
      glVertex3f(etaCoord, minPhi, 0);
      glVertex3f(etaCoord, maxPhi, 0);
      glEnd();
    }

    // Show axis tick marks and labels

    {
      TGLCapabilitySwitch lightsOff(GL_LIGHTING, false);

      TGLAxis ap;
      ap.SetLineColor(fM->fGridColor);
      ap.SetTextColor(fM->fGridColor);
      TGLVector3 start, end;

      start.Set(minEta, minPhi, 0);
      end.Set(maxEta, minPhi, 0);
      ap.PaintGLAxis(start.CArr(), end.CArr(), fM->fMinEta, fM->fMaxEta, 205);

      start.Set(maxEta, minPhi, 0);
      end.Set(maxEta, maxPhi, 0);
      ap.PaintGLAxis(start.CArr(), end.CArr(), fM->fMinPhi, fM->fMaxPhi, 205);
    }

  }

//   // Finding the maximum energy
// 
//   std::vector<AliAODTrack>::iterator k = fM->fTracks.begin();
//   std::vector<AliAODJet>::iterator   j = fM->fJets.begin();
// 
//   Double_t eJetMax = 0., eTrackMax = 0., eMax;
// 
//   while (j != fM->fJets.end())
//   {
//     if (j->E() > eJetMax) eJetMax = j->E();
//     ++j;
//   }
// 
//   while (k != fM->fTracks.end())
//   {
//     if (k->E() > eTrackMax) eTrackMax = k->E();
//     ++k;
//   }
// 
//   eMax = eJetMax > eTrackMax ? eJetMax : eTrackMax;
// 
//   // Rendering Jets and Tracks in the Eta-Phi plane
// 
//   Int_t     nCol = gStyle->GetNumberOfColors();
//   Float_t col[4] = { 0., 0., 0., 0.75};
// 
//   glBlendFunc(GL_SRC_ALPHA,GL_ONE);
//   TGLUtil::SetDrawQuality(6);
// 
// 
//   glPushName(0);
//   if (fM->fRnrJets)
//   {
//     glEnable(GL_BLEND);		   // Turn Blending On
//     glDisable(GL_DEPTH_TEST); // Turn Depth Testing Off
//     UInt_t jetid = 0;
// 
// 
//     j = fM->fJets.begin();
// 
//     while (j != fM->fJets.end())
//     {
//       eta = j->Eta();
//       phi = j->Phi();
//       e   = j->E();
// 
//       x = eta*(fM->fEtaScale);
//       y = phi*(fM->fPhiScale) - 350;
// 
//       Int_t colBin = TMath::Min((Int_t) ((nCol-2)*e*
// 					 TMath::Power(10.,fM->fEnergyColorScale)/(eMax)),nCol-2);
//       Int_t colIdx = gStyle->GetColorPalette(colBin);
//       TColor* c    = gROOT->GetColor(colIdx);
// 
//       if(c)
//       {
// 	col[0] = c->GetRed();
// 	col[1] = c->GetGreen();
// 	col[2] = c->GetBlue();
//       }
// 
//       glLoadName(jetid);
//       TGLUtil::DrawLine(TGLVertex3(x,y,0.),
// 			TGLVector3(0.,0.,TMath::Log(e + 1.)*fM->fEnergyScale),
// 			TGLUtil::kLineHeadArrow, 25.0, col);
//       ++j; ++jetid;
//     }
//   }

//   col[3] = 1.0;
// 
//   glPushName(1);
//   if(fM->fRnrTracks)
//   {
//     glDisable(GL_BLEND);		   // Turn Blending Off
//     glEnable(GL_DEPTH_TEST);	 // Turn Depth Testing On
//     UInt_t trackid = 0;
// 
//     k = fM->fTracks.begin();
// 
//     while (k != fM->fTracks.end())
//     {
//       eta = k->Eta();
//       phi = k->Phi();
//       e   = k->E();
// 
//       if (e < 0.)
//       {
// 	//			printf(" WARNING: Particle with negative energy has been found.\n");
// 	//			printf(" PARTICLE NOT DISPLAYED. TrackID: %i\n", trackid);
// 	++k; ++trackid;
// 	continue;
//       }
// 
//       x = eta*(fM->fEtaScale);
//       y = phi*(fM->fPhiScale) - 350;
// 
//       Int_t colBin = TMath::Min((Int_t) ((nCol-2)*e*
// 					 TMath::Power(10.,fM->fEnergyColorScale)/(eMax)),nCol-2);
//       Int_t colIdx = gStyle->GetColorPalette(colBin);
//       TColor* c    = gROOT->GetColor(colIdx);
// 
//       if(c)
//       {
// 	col[0] = c->GetRed();
// 	col[1] = c->GetGreen();
// 	col[2] = c->GetBlue();
//       }
// 
//       glLoadName(trackid);
//       TGLUtil::DrawLine( TGLVertex3(x,y,0.),
// 			 TGLVector3(0.,0.,TMath::Log(e + 1.)*fM->fEnergyScale),
// 			 TGLUtil::kLineHeadArrow, 5.0, col);
// 
//       ++k; ++trackid;
//     }
//   }

//  glPopName();
//  TGLUtil::ResetDrawQuality();
}

// /******************************************************************************/
// 
void AliEveJetPlaneGL::ProcessSelection(TGLRnrCtx & /*rnrCtx*/, TGLSelectRecord & rec)
{
//   // Process selection and print jet information.
// 
//   //   printf("beep %u\n", rec.GetN());
//   //   rec.Print();
//   static Int_t jet1State;
//   static Int_t jet2State;
//   static Int_t track1State;
//   static Int_t track2State;
// 
//   if (fM->fOneSelection)
//   {
//     printf("\n");
// 
//     if (rec.GetN() == 2)
//     {
//       AliAODJet v = fM->fJets[rec.GetItem(1)];
//       printf("Jet 4-momentum: %f, %f, %f, %f \n", v.Px(),v.Py(),v.Pz(),v.Pt() );
//       printf("Eta-Phi values: %f, %f\n", v.Eta(), v.Phi());
//     }
// 
//     if (rec.GetN() == 3)
//     {
//       AliAODTrack v = fM->fTracks[rec.GetItem(2)];
//       printf("TEveTrack 4-momentum: %f, %f, %f, %f \n", v.Px(),v.Py(),v.Pz(),v.Pt() );
//       printf("Eta-Phi values: %f, %f\n", v.Eta(), v.Phi());
//     }
//   }
// 
//   if (fM->fTwoSelection)
//   {
// 
//     if ( fM->fSelectionFlag == 1)
//     {
//       if (rec.GetN() == 2)
//       {
// 	fM->SetJet1(&(fM->fJets[rec.GetItem(1)]));
//         jet1State = 1;
// 	track1State = 0;
//       }
// 
//       if (rec.GetN() == 3)
//       {
// 	fM->SetTrack1(&(fM->fTracks[rec.GetItem(1)]));
// 	jet1State = 0;
// 	track1State = 1;
//       }
// 
//       fM->SetSelectionFlag(2);
// 
//       return;
//     }
// 
//     if ( fM->fSelectionFlag == 2)
//     {
//       printf("\n");
//       if (rec.GetN() == 2)
//       {
// 	fM->SetJet2(&(fM->fJets[rec.GetItem(1)]));
//         jet2State = 1;
// 	track2State = 0;
//       }
// 
//       if (rec.GetN() == 3)
//       {
// 	fM->SetTrack2(&(fM->fTracks[rec.GetItem(1)]));
// 	jet2State = 0;
// 	track2State = 1;
//       }
// 
//       printf("Jet: %i, TEveTrack: %i \n", jet1State, track1State);
//       printf("Jet: %i, TEveTrack: %i \n\n", jet2State, track2State);
// 
//       if(jet1State && jet2State)
//       {
// 	Double_t eta1, eta2, phi1, phi2, d;
// 
// 	eta1 = (fM->GetJet1()).Eta();
// 	eta2 = (fM->GetJet2()).Eta();
// 	phi1 = (fM->GetJet1()).Phi();
//         phi2 = (fM->GetJet2()).Phi();
// 
// 	d = TMath::Sqrt(TMath::Power(eta2-eta1,2) + TMath::Power(phi2-phi1,2));
// 
// 	printf("Eta-Phi: %f, %f\n", eta1, phi1);
// 	printf("Eta-Phi: %f, %f\n", eta2, phi2);
// 	printf("Eta-Phi Distance: %f\n", d);
//       }
// 
//       fM->SetSelectionFlag(1);
//     }
// 
//   }
// 
}




