/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
/*
$Log$
Revision 1.1.2.6  2000/06/09 21:57:40  morsch
Most coding rule violations corrected.

Revision 1.1.2.5  2000/05/05 11:33:56  morsch
Log inside comments.

Revision 1.1.2.4  2000/05/05 10:12:09  morsch
Log messages included
*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  This class contains the points for the ALICE event display               //
//                                                                           //
//Begin_Html
/*
<img src="gif/AliMUONPointsClass.gif">
*/
//End_Html
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliMUONPoints.h"
#include "AliMUONDisplay.h"
#include "AliRun.h"
#include "AliMUON.h"
#include "AliMUONChamber.h"
#include "AliMUONSegmentation.h"
#include "AliMUONResponse.h"
#include "AliMUONHit.h"
#include "AliMUONPadHit.h"
#include "AliMUONDigit.h"
#include "AliMUONRawCluster.h"


#include <TPolyMarker3D.h>
#include <TMatrix.h>
#include <TPad.h>
#include <TVirtualPad.h>
#include <TPolyLine3D.h>
#include <TPaveText.h>
#include <TView.h>
#include <TMath.h>

const Int_t kMaxNipx=400, kMaxNipy=800;
 
ClassImp(AliMUONPoints)

//_____________________________________________________________________________
AliMUONPoints::AliMUONPoints()
{
  //
  // Default constructor
  //
  fHitIndex = 0;
  fTrackIndex = 0;
  fDigitIndex = 0;
  fMarker[0] = fMarker[1] = fMarker[2]=0;
  fMatrix = 0;
}

//_____________________________________________________________________________
AliMUONPoints::AliMUONPoints(Int_t npoints)
  :AliPoints(npoints)
{
  //
  // Standard constructor
  //
  fHitIndex = 0;
  fTrackIndex = 0;
  fDigitIndex = 0;
  fMarker[0] = fMarker[1] = fMarker[2]=0;
  fMatrix = 0;
}

AliMUONPoints::AliMUONPoints(const AliMUONPoints& points)
{
// Copy constructor
}
	 
//_____________________________________________________________________________
AliMUONPoints::~AliMUONPoints()
{
  //
  // Default destructor
  //
  fHitIndex = 0;
  fTrackIndex = 0;
  fDigitIndex = 0;
  for (Int_t i=0;i<3;i++){
      if (fMarker[i]) delete fMarker[i];
  }
  fMatrix = 0;
}

//_____________________________________________________________________________
void AliMUONPoints::DumpHit()
{
  //
  //   Dump hit corresponding to this point
  //
  AliMUONHit *hit = GetHit();
  if (hit) hit->Dump();
}

//_____________________________________________________________________________
void AliMUONPoints::DumpDigit()
{
  //
  //   Dump digit corresponding to this point
  //
  AliMUONDigit *digit = GetDigit();
  if (digit) digit->Dump();
}

//_____________________________________________________________________________
void AliMUONPoints::InspectHit()
{
  //
  //   Inspect hit corresponding to this point
  //

  if (fHitIndex < 0 ) return;
  TVirtualPad *padsav = gPad;
  AliMUONHit *hit = GetHit();
  if (hit) hit->Inspect();
  TVirtualPad *padinspect = (TVirtualPad*)(gROOT->GetListOfCanvases())->FindObject("inspect");
   padinspect->cd();
   Float_t xmin = gPad->GetX1();
   Float_t xmax = gPad->GetX2();
   Float_t ymin = gPad->GetY1();
   Float_t ymax = gPad->GetY2();
   Float_t dy   = ymax-ymin;

      TPaveText *pad = new TPaveText(xmin, ymin+0.1*dy, xmax, ymin+0.15*dy);
      pad->SetBit(kCanDelete);
      pad->SetFillColor(42);
      pad->Draw();
      char ptitle[100];
      sprintf(ptitle," %s , fTrack: %d  fTrackIndex: %d ",GetName(),fIndex,fTrackIndex);
      pad->AddText(ptitle);
      padinspect->cd();
      padinspect->Update();
  if (padsav) padsav->cd();

}

//_____________________________________________________________________________
void AliMUONPoints::InspectDigit()
{
  //
  //   Inspect digit corresponding to this point
  //
  if (fDigitIndex < 0) return;
  TVirtualPad *padsav = gPad;
  AliMUONDigit *digit = GetDigit();
  if (digit) digit->Inspect();
  TVirtualPad *padinspect = (TVirtualPad*)(gROOT->GetListOfCanvases())->FindObject("inspect");
   padinspect->cd();
   Float_t xmin = gPad->GetX1();
   Float_t xmax = gPad->GetX2();
   Float_t ymin = gPad->GetY1();
   Float_t ymax = gPad->GetY2();
   Float_t dy   = ymax-ymin;

      TPaveText *pad = new TPaveText(xmin, ymin+0.1*dy, xmax, ymin+0.25*dy);
      pad->SetBit(kCanDelete);
      pad->SetFillColor(42);
      pad->Draw();
      char ptitle[11][100];
      //      sprintf(ptitle[11],"Tracks making this digit");
      //      pad->AddText(ptitle[11]);
  for (int i=0;i<10;i++) {
      if (digit->fTracks[i] == 0) continue;  
      sprintf(ptitle[i],"fTrackIndex: %d  Charge: %d",digit->fTracks[i],digit->fTcharges[i]);
      pad->AddText(ptitle[i]);
  }
      padinspect->cd();
      padinspect->Update();
  if (padsav) padsav->cd();
      
}

//_____________________________________________________________________________
Int_t AliMUONPoints::GetTrackIndex()
{
  //
  //   Dump digit corresponding to this point
  //

  this->Inspect();
  /*
  if (fDigitIndex != 0) {
    Int_t ncol=this->fMatrix->GetNcols();
    for (int i=0;i<ncol;i++) {
        printf(" track charge %f %f \n",(*(this->fMatrix))(0,i),(*(this->fMatrix))(1,i));
    }
  }
  */
  return fTrackIndex;
}

//_____________________________________________________________________________
AliMUONHit *AliMUONPoints::GetHit() const
{
  //
  //   Returns pointer to hit index in AliRun::fParticles
  //
  AliMUON *pMUON  = (AliMUON*)gAlice->GetModule("MUON");
  gAlice->TreeH()->GetEvent(fTrackIndex);
  TClonesArray *muonHits  = pMUON->Hits();
  Int_t nhits = muonHits->GetEntriesFast();
  if (fHitIndex < 0 || fHitIndex >= nhits) return 0;
  return (AliMUONHit*)muonHits->UncheckedAt(fHitIndex);
}

//_____________________________________________________________________________
AliMUONDigit *AliMUONPoints::GetDigit() const
{
  //
  //   Returns pointer to digit index in AliRun::fParticles
  //

  AliMUONDisplay *display=(AliMUONDisplay*)gAlice->Display();
  Int_t chamber=display->GetChamber();
  Int_t cathode=display->GetCathode();
   
  AliMUON *pMUON  = (AliMUON*)gAlice->GetModule("MUON");
  TClonesArray *muonDigits  = pMUON->DigitsAddress(chamber-1);
  Int_t nent=(Int_t)gAlice->TreeD()->GetEntries();
  gAlice->TreeD()->GetEvent(nent-2+cathode-1);
  //gAlice->TreeD()->GetEvent(cathode);
  Int_t ndigits = muonDigits->GetEntriesFast();
  if (fDigitIndex < 0 || fDigitIndex >= ndigits) return 0;
  return (AliMUONDigit*)muonDigits->UncheckedAt(fDigitIndex);
}
//_____________________________________________________________________________

AliMUONPoints& AliMUONPoints::operator= (const AliMUONPoints& rhs)
{
    return *this;
}
