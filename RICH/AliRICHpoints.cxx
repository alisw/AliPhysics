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
*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  This class contains the points for the ALICE event display               //
//                                                                           //
//Begin_Html
/*
<img src="gif/AliRICHpointsClass.gif">
*/
//End_Html
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
#include "AliRICHdisplay.h"
#include "AliRICHpoints.h"
#include "AliRun.h"
#include "TPad.h"
#include "TView.h"
#include "TMath.h"

const Int_t MAX_Nipx=400, MAX_Nipy=800;
 
ClassImp(AliRICHpoints)

//_____________________________________________________________________________
AliRICHpoints::AliRICHpoints()
{
  //
  // Default constructor
  //
  fHitIndex = 0;
  fTrackIndex = 0;
  fDigitIndex = 0;
}

//_____________________________________________________________________________
AliRICHpoints::AliRICHpoints(Int_t npoints)
  :AliPoints(npoints)
{
  //
  // Standard constructor
  //
  fHitIndex = 0;
  fTrackIndex = 0;
  fDigitIndex = 0;
}
	 
//_____________________________________________________________________________
AliRICHpoints::~AliRICHpoints()
{
  //
  // Default destructor
  //
  fHitIndex = 0;
  fTrackIndex = 0;
  fDigitIndex = 0;
}

//_____________________________________________________________________________
void AliRICHpoints::DumpHit()
{
  //
  //   Dump hit corresponding to this point
  //
  AliRICHhit *hit = GetHit();
  if (hit) hit->Dump();
}

//_____________________________________________________________________________
void AliRICHpoints::DumpDigit()
{
  //
  //   Dump digit corresponding to this point
  //
  AliRICHdigit *digit = GetDigit();
  if (digit) digit->Dump();
}

//_____________________________________________________________________________
void AliRICHpoints::InspectHit()
{
  //
  //   Inspect hit corresponding to this point
  //
  AliRICHhit *hit = GetHit();
  if (hit) hit->Inspect();
}

//_____________________________________________________________________________
void AliRICHpoints::InspectDigit()
{
  //
  //   Inspect digit corresponding to this point
  //
  AliRICHdigit *digit = GetDigit();
  if (digit) digit->Inspect();
}

//_____________________________________________________________________________
Int_t AliRICHpoints::GetTrackIndex()
{
  //
  //   Dump digit corresponding to this point
  //
  printf("GetTrackIndex - fTrackIndex %d \n",fTrackIndex);
  this->Inspect();
  return fTrackIndex;
}

//_____________________________________________________________________________
AliRICHhit *AliRICHpoints::GetHit() const
{
  //
  //   Returns pointer to hit index in AliRun::fParticles
  //
  AliRICH *RICH  = (AliRICH*)gAlice->GetDetector("RICH");
  gAlice->TreeH()->GetEvent(fTrackIndex);
  TClonesArray *RICHhits  = RICH->Hits();
  Int_t nhits = RICHhits->GetEntriesFast();
  if (fHitIndex < 0 || fHitIndex >= nhits) return 0;
  return (AliRICHhit*)RICHhits->UncheckedAt(fHitIndex);
}

//_____________________________________________________________________________
AliRICHdigit *AliRICHpoints::GetDigit() const
{
  //
  //   Returns pointer to digit index in AliRun::fParticles
  //

  AliRICHdisplay *display=(AliRICHdisplay*)gAlice->Display();
  Int_t chamber=display->GetChamber();
  Int_t cathode=display->GetCathode();
   
  AliRICH *RICH  = (AliRICH*)gAlice->GetDetector("RICH");
  TClonesArray *RICHdigits  = RICH->DigitsAddress(chamber-1);
  gAlice->TreeD()->GetEvent(cathode);
  Int_t ndigits = RICHdigits->GetEntriesFast();
  if (fDigitIndex < 0 || fDigitIndex >= ndigits) return 0;
  return (AliRICHdigit*)RICHdigits->UncheckedAt(fDigitIndex);
}
//_____________________________________________________________________________
struct Bin {
   const AliRICHdigit *dig;
   int idx;
   Bin() {dig=0; idx=-1;}
};

struct PreCluster : public AliRICHreccluster {
   const AliRICHdigit* summit;
   int idx;
   int cut;
   int npeaks;
   PreCluster() : AliRICHreccluster() {cut=npeaks=0;}
};
//_____________________________________________________________________________

static void FindCluster(AliRICHchamber *iChamber, AliRICHsegmentation *segmentation, int i, int j, Bin bins[MAX_Nipx][MAX_Nipy], PreCluster &c) 

{

  //
  // Find clusters
  //

  printf("I'm in FindCluster \n"); 

  Bin& b=bins[i][j];
  Int_t q=b.dig->fSignal;

  printf("FindCluster - i j q %d %d %d\n",i,j,q);
  
  if (q<0) { 
    q=-q;
    c.cut=1;
  } 
  if (b.idx >= 0 && b.idx != c.idx) {
    c.idx=b.idx;
    c.npeaks++;
  }
  
  if (q > TMath::Abs(c.summit->fSignal)) c.summit=b.dig;

  Int_t npx  = segmentation->Npx();
  Int_t npy  = segmentation->Npy();
  Float_t x,y;
  segmentation->GetPadCxy(i-npx, j-npy, x,y);
  printf("FindCluster - x  y %f %f \n",x,y);


  c.fX += q*x;
  c.fY += q*y;
  c.fQ += q;
  
  b.dig = 0;  b.idx = c.idx;
  
  if (bins[i-1][j].dig) FindCluster(iChamber,segmentation,i-1,j,bins,c);
  if (bins[i][j-1].dig) FindCluster(iChamber,segmentation,i,j-1,bins,c);
  if (bins[i+1][j].dig) FindCluster(iChamber,segmentation,i+1,j,bins,c);
  if (bins[i][j+1].dig) FindCluster(iChamber,segmentation,i,j+1,bins,c);

}

//_____________________________________________________________________________

void AliRICHpoints::GetCenterOfGravity()
{
  //
  // simple RICH cluster finder from digits -- finds neighbours and 
  // calculates center of gravity for the cluster
  //
  const Int_t MAX_nipx=400, MAX_nipy=800;
  printf("\n Hallo world");
  AliRICHdisplay *display=(AliRICHdisplay*)gAlice->Display();
  Int_t chamber=display->GetChamber();
  Int_t cathode=display->GetCathode();
   
  AliRICH *RICH  = (AliRICH*)gAlice->GetDetector("RICH");
  AliRICHchamber *iChamber;
  AliRICHsegmentation *segmentation;
  iChamber =&(RICH->Chamber(chamber-1));
  segmentation=iChamber->GetSegmentationModel(cathode);
  Int_t npx  = segmentation->Npx();
  Int_t npy  = segmentation->Npy();
  Float_t zpos=iChamber->ZPosition();
  
  TClonesArray *RICHdigits  = RICH->DigitsAddress(chamber-1);
  gAlice->TreeD()->GetEvent(cathode);
  Int_t ndigits = RICHdigits->GetEntriesFast();
  if (fDigitIndex < 0 || fDigitIndex >= ndigits) return;

  AliRICHdigit  *dig;
  dig=(AliRICHdigit*)RICHdigits->UncheckedAt(fDigitIndex);
  Int_t ipx=dig->fPadX;
  Int_t ipy=dig->fPadY;
  Bin bins[MAX_nipx][MAX_nipy]; 
  bins[ipx+npx][ipy+npy].dig=dig;
    
  int ndig;
  int ncls=0;
  for (ndig=0; ndig<ndigits; ndig++) {
      dig = (AliRICHdigit*)RICHdigits->UncheckedAt(ndig);
      int i=dig->fPadX, j=dig->fPadY;
      bins[i+npx][j+npy].dig=dig;
  }

  PreCluster c; c.summit=bins[ipx+npx][ipy+npy].dig; c.idx=ncls;
  FindCluster(iChamber,segmentation,ipx+npx, ipy+npy, bins, c);
  if (c.npeaks>1) {
      printf("GetCenterOfGravity -- more than one peak");
  }
  c.fX /= c.fQ;
  c.fY /= c.fQ;
  printf("GetCenterOfGravity - c.fX c.fY c.fQ %f %f %d \n",c.fX,c.fY,c.fQ);
  
  c.fTracks[0]=c.summit->fTracks[0];
  c.fTracks[1]=c.summit->fTracks[1];
  c.fTracks[2]=c.summit->fTracks[2];
  ncls++;
  AliRICHpoints *points = 0;
  points = new AliRICHpoints(1);
  points->SetMarkerColor(kYellow);
  points->SetMarkerStyle(5);
  points->SetMarkerSize(1.);
  points->SetPoint(0,c.fX,c.fY,zpos);
  points->Draw();
  
  printf("GetCenterOfGravity -- ncls %d \n",ncls);

}



