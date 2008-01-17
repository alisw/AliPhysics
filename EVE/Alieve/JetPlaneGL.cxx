// $Header$

#include "JetPlaneGL.h"
#include <Alieve/JetPlane.h>

#include <TGLRnrCtx.h>
#include <TGLSelectRecord.h>
#include <TGLIncludes.h>

#include "TGLUtil.h"
#include "TGLAxis.h"

#include <TColor.h>
#include <TStyle.h>
#include <TROOT.h>
using namespace Alieve;

//______________________________________________________________________
// JetPlaneGL
//

ClassImp(JetPlaneGL)

JetPlaneGL::JetPlaneGL() : TGLObject(), fM(0)
{
  fDLCache = kFALSE; // Disable display list -- axis pain.
}

JetPlaneGL::~JetPlaneGL()
{}

/**************************************************************************/

Bool_t JetPlaneGL::SetModel(TObject* obj, const Option_t* /*opt*/)
{
  if(SetModelCheckClass(obj, Alieve::JetPlane::Class())) {
    fM = dynamic_cast<JetPlane*>(obj);
    return kTRUE;
  }
  return kFALSE;
}

void JetPlaneGL::SetBBox()
{
  // !! This ok if master sub-classed from TAttBBox
  SetAxisAlignedBBox(((JetPlane*)fExternalObj)->AssertBBox());
}

/**************************************************************************/

void JetPlaneGL::DirectDraw(TGLRnrCtx & /*rnrCtx*/) const
{

  Float_t minEta = (fM->fMinEta)*(fM->fEtaScale);
  Float_t maxEta = (fM->fMaxEta)*(fM->fEtaScale);
  Float_t minPhi = (fM->fMinPhi)*(fM->fPhiScale);
  Float_t maxPhi = (fM->fMaxPhi)*(fM->fPhiScale);
  Float_t phiCoord, etaCoord, dPhi, dEta;
  Double_t eta, phi, E, x, y;

  // Show frame for Eta-Phi coordinates

  glBegin(GL_LINE_LOOP);
  glVertex3f( minEta, minPhi, 0);
  glVertex3f( maxEta, minPhi, 0);
  glVertex3f( maxEta, maxPhi, 0);
  glVertex3f( minEta, maxPhi, 0);
  glEnd();

  // Show grid in Eta-Phi coordinates

  dPhi = (maxPhi-minPhi)/(fM->fNPhiDiv - 1);
  dEta = (maxEta-minEta)/(fM->fNEtaDiv - 1);

  for (int count=1; count < fM->fNPhiDiv-1; count++)
  {
    phiCoord = minPhi + count*dPhi;
    glBegin(GL_LINES);
    glVertex3f( minEta, phiCoord, 0);
    glVertex3f( maxEta, phiCoord, 0);
    glEnd();
  }

  for (int count=1; count < fM->fNEtaDiv-1; count++)
  {
    etaCoord = minEta + count*dEta;
    glBegin(GL_LINES);
    glVertex3f( etaCoord, minPhi, 0);
    glVertex3f( etaCoord, maxPhi, 0);
    glEnd();
  }

  // Show axis tick marks and labels

  {
    TGLCapabilitySwitch lights_off(GL_LIGHTING, false);

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

  // Finding the maximum energy

  std::vector<AliAODTrack>::iterator k = fM->fTracks.begin();
  std::vector<AliAODJet>::iterator   j = fM->fJets.begin();

  Double_t eJetMax = 0., eTrackMax = 0., eMax;

  while (j != fM->fJets.end())
  {
    if (j->E() > eJetMax) eJetMax = j->E();
    ++j;
  }

  while (k != fM->fTracks.end())
  {
    if (k->E() > eTrackMax) eTrackMax = k->E();
    ++k;
  }

  eMax = eJetMax > eTrackMax ? eJetMax : eTrackMax;

  // Rendering Jets and Tracks in the Eta-Phi plane

  Int_t     nCol = gStyle->GetNumberOfColors();
  Float_t col[4] = { 0., 0., 0., 0.75};

  glBlendFunc(GL_SRC_ALPHA,GL_ONE);
  TGLUtil::SetDrawQuality(6);


  glPushName(0);
  if (fM->fRnrJets)
  {
    glEnable(GL_BLEND);		   // Turn Blending On
    glDisable(GL_DEPTH_TEST); // Turn Depth Testing Off
    UInt_t jetid = 0;


    j = fM->fJets.begin();

    while (j != fM->fJets.end())
    {
      eta = j->Eta();
      phi = j->Phi();
      E   = j->E();

      x = eta*(fM->fEtaScale);
      y = phi*(fM->fPhiScale);

      Int_t colBin = TMath::Min((Int_t) ((nCol-2)*E*
					 TMath::Power(10.,fM->fEnergyColorScale)/(eMax)),nCol-2);
      Int_t colIdx = gStyle->GetColorPalette(colBin);
      TColor* c    = gROOT->GetColor(colIdx);

      if(c)
      {
	col[0] = c->GetRed();
	col[1] = c->GetGreen();
	col[2] = c->GetBlue();
      }

      glLoadName(jetid);
      TGLUtil::DrawLine( TGLVertex3(x,y,0.),
			 TGLVector3(0.,0.,TMath::Log(E + 1.)*fM->fEnergyScale),
			 TGLUtil::kLineHeadArrow, 25.0, col);

      ++j; ++jetid;
    }
  }

  col[3] = 1.0;

  glPushName(1);
  if(fM->fRnrTracks)
  {
    glDisable(GL_BLEND);		   // Turn Blending Off
    glEnable(GL_DEPTH_TEST);	 // Turn Depth Testing On
    UInt_t trackid = 0;


    k = fM->fTracks.begin();

    while (k != fM->fTracks.end())
    {
      eta = k->Eta();
      phi = k->Phi();
      E   = k->E();

      if (E < 0.)
      {
	//			printf(" WARNING: Particle with negative energy has been found.\n");
	//			printf(" PARTICLE NOT DISPLAYED. TrackID: %i\n", trackid);
	++k; ++trackid;
	continue;
      }

      x = eta*(fM->fEtaScale);
      y = phi*(fM->fPhiScale);

      Int_t colBin = TMath::Min((Int_t) ((nCol-2)*E*
					 TMath::Power(10.,fM->fEnergyColorScale)/(eMax)),nCol-2);
      Int_t colIdx = gStyle->GetColorPalette(colBin);
      TColor* c    = gROOT->GetColor(colIdx);

      if(c)
      {
	col[0] = c->GetRed();
	col[1] = c->GetGreen();
	col[2] = c->GetBlue();
      }

      glLoadName(trackid);
      TGLUtil::DrawLine( TGLVertex3(x,y,0.),
			 TGLVector3(0.,0.,TMath::Log(E + 1.)*fM->fEnergyScale),
			 TGLUtil::kLineHeadArrow, 5.0, col);

      ++k; ++trackid;
    }
  }

  glPopName();
  TGLUtil::ResetDrawQuality();



}

/**************************************************************************/

void JetPlaneGL::ProcessSelection(TGLRnrCtx & /*rnrCtx*/, TGLSelectRecord & rec)
{
  //   printf("beep %u\n", rec.GetN());
  //   rec.Print();
  static Int_t jet1State;
  static Int_t jet2State;
  static Int_t track1State;
  static Int_t track2State;

  if (fM->fOneSelection)
  {
    printf("\n");

    if (rec.GetN() == 2)
    {
      AliAODJet v = fM->fJets[rec.GetItem(1)];
      printf("Jet 4-momentum: %f, %f, %f, %f \n", v.Px(),v.Py(),v.Pz(),v.Pt() );
      printf("Eta-Phi values: %f, %f\n", v.Eta(), v.Phi());
    }

    if (rec.GetN() == 3)
    {
      AliAODTrack v = fM->fTracks[rec.GetItem(2)];
      printf("TEveTrack 4-momentum: %f, %f, %f, %f \n", v.Px(),v.Py(),v.Pz(),v.Pt() );
      printf("Eta-Phi values: %f, %f\n", v.Eta(), v.Phi());
    }
  }

  if (fM->fTwoSelection)
  {

    if ( fM->fSelectionFlag == 1)
    {
      if (rec.GetN() == 2)
      {
	fM->SetJet1(&(fM->fJets[rec.GetItem(1)]));
        jet1State = 1;
	track1State = 0;
      }

      if (rec.GetN() == 3)
      {
	fM->SetTrack1(&(fM->fTracks[rec.GetItem(1)]));
	jet1State = 0;
	track1State = 1;
      }

      fM->SetSelectionFlag(2);

      return;
    }

    if ( fM->fSelectionFlag == 2)
    {
      printf("\n");
      if (rec.GetN() == 2)
      {
	fM->SetJet2(&(fM->fJets[rec.GetItem(1)]));
        jet2State = 1;
	track2State = 0;
      }

      if (rec.GetN() == 3)
      {
	fM->SetTrack2(&(fM->fTracks[rec.GetItem(1)]));
	jet2State = 0;
	track2State = 1;
      }

      printf("Jet: %i, TEveTrack: %i \n", jet1State, track1State);
      printf("Jet: %i, TEveTrack: %i \n\n", jet2State, track2State);

      if(jet1State && jet2State)
      {
	Double_t Eta1, Eta2, Phi1, Phi2, d;

	Eta1 = (fM->GetJet1()).Eta();
	Eta2 = (fM->GetJet2()).Eta();
	Phi1 = (fM->GetJet1()).Phi();
        Phi2 = (fM->GetJet2()).Phi();

	d = TMath::Sqrt(TMath::Power(Eta2-Eta1,2) + TMath::Power(Phi2-Phi1,2));

	printf("Eta-Phi: %f, %f\n", Eta1, Phi1 );
	printf("Eta-Phi: %f, %f\n", Eta2, Phi2 );
	printf("Eta-Phi Distance: %f\n", d);
      }

      fM->SetSelectionFlag(1);
    }

  }

}




