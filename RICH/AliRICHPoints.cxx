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
<img src="gif/AliRICHPointsClass.gif">
*/
//End_Html
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
#include "AliRICHDisplay.h"
#include "AliRICHPoints.h"
#include "AliRun.h"
#include "TPad.h"
#include "TView.h"
#include "TMath.h"

const Int_t MAX_Nipx=400, MAX_Nipy=800;
 
ClassImp(AliRICHPoints)

//_____________________________________________________________________________
AliRICHPoints::AliRICHPoints()
{
  //
  // Default constructor
  //
  fHitIndex = 0;
  fTrackIndex = 0;
  fDigitIndex = 0;
  fMarker[0] = fMarker[1] = fMarker[2]=0;
}

//_____________________________________________________________________________
AliRICHPoints::AliRICHPoints(Int_t npoints)
  :AliPoints(npoints)
{
  //
  // Standard constructor
  //
  fHitIndex = 0;
  fTrackIndex = 0;
  fDigitIndex = 0;
  fMarker[0] = fMarker[1] = fMarker[2]=0;
}
	 
//_____________________________________________________________________________
AliRICHPoints::~AliRICHPoints()
{
  //
  // Default destructor
  //
  fHitIndex = 0;
  fTrackIndex = 0;
  fDigitIndex = 0;
}

//_____________________________________________________________________________
void AliRICHPoints::DumpHit()
{
  //
  //   Dump hit corresponding to this point
  //
  AliRICHHit *hit = GetHit();
  if (hit) hit->Dump();
}

//_____________________________________________________________________________
void AliRICHPoints::DumpDigit()
{
  //
  //   Dump digit corresponding to this point
  //
  AliRICHDigit *digit = GetDigit();
  if (digit) digit->Dump();
}

//_____________________________________________________________________________
void AliRICHPoints::InspectHit()
{
  //
  //   Inspect hit corresponding to this point
  //
  AliRICHHit *hit = GetHit();
  if (hit) hit->Inspect();
}

//_____________________________________________________________________________
void AliRICHPoints::InspectDigit()
{
  //
  //   Inspect digit corresponding to this point
  //
  AliRICHDigit *digit = GetDigit();
  if (digit) digit->Inspect();
}

//_____________________________________________________________________________
Int_t AliRICHPoints::GetTrackIndex()
{
  //
  //   Dump digit corresponding to this point
  //
  printf("GetTrackIndex - fTrackIndex %d \n",fTrackIndex);
  this->Inspect();
  return fTrackIndex;
}

//_____________________________________________________________________________
AliRICHHit *AliRICHPoints::GetHit() const
{
  //
  //   Returns pointer to hit index in AliRun::fParticles
  //
  AliRICH *RICH  = (AliRICH*)gAlice->GetDetector("RICH");
  gAlice->TreeH()->GetEvent(fTrackIndex);
  TClonesArray *RICHhits  = RICH->Hits();
  Int_t nhits = RICHhits->GetEntriesFast();
  if (fHitIndex < 0 || fHitIndex >= nhits) return 0;
  return (AliRICHHit*)RICHhits->UncheckedAt(fHitIndex);
}

//_____________________________________________________________________________
AliRICHDigit *AliRICHPoints::GetDigit() const
{
  //
  //   Returns pointer to digit index in AliRun::fParticles
  //

  AliRICHDisplay *display=(AliRICHDisplay*)gAlice->Display();
  Int_t chamber=display->GetChamber();
  Int_t cathode=display->GetCathode();
   
  AliRICH *RICH  = (AliRICH*)gAlice->GetDetector("RICH");
  TClonesArray *RICHdigits  = RICH->DigitsAddress(chamber-1);
  gAlice->TreeD()->GetEvent(cathode);
  Int_t ndigits = RICHdigits->GetEntriesFast();
  if (fDigitIndex < 0 || fDigitIndex >= ndigits) return 0;
  return (AliRICHDigit*)RICHdigits->UncheckedAt(fDigitIndex);
}
//_____________________________________________________________________________
struct Bin {
   const AliRICHDigit *dig;
   int idx;
   Bin() {dig=0; idx=-1;}
};

struct PreCluster : public AliRICHRawCluster {
   const AliRICHDigit* summit;
   int idx;
   int cut;
   int npeaks;
   PreCluster() : AliRICHRawCluster() {cut=npeaks=0;}
};
//_____________________________________________________________________________

static void FindCluster(AliRICHChamber *iChamber, AliRICHSegmentation *segmentation, int i, int j, Bin bins[MAX_Nipx][MAX_Nipy], PreCluster &c) 

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

void AliRICHPoints::GetCenterOfGravity()
{
  //
  // simple RICH cluster finder from digits -- finds neighbours and 
  // calculates center of gravity for the cluster
  //
  const Int_t MAX_nipx=400, MAX_nipy=800;
  printf("\n Hallo world");
  AliRICHDisplay *display=(AliRICHDisplay*)gAlice->Display();
  Int_t chamber=display->GetChamber();
  Int_t cathode=display->GetCathode();
   
  AliRICH *RICH  = (AliRICH*)gAlice->GetDetector("RICH");
  AliRICHChamber *iChamber;
  AliRICHSegmentation *segmentation;
  iChamber =&(RICH->Chamber(chamber-1));
  segmentation=iChamber->GetSegmentationModel(cathode);
  Int_t npx  = segmentation->Npx();
  Int_t npy  = segmentation->Npy();
  Float_t zpos=iChamber->ZPosition();
  
  TClonesArray *RICHdigits  = RICH->DigitsAddress(chamber-1);
  gAlice->TreeD()->GetEvent(cathode);
  Int_t ndigits = RICHdigits->GetEntriesFast();
  if (fDigitIndex < 0 || fDigitIndex >= ndigits) return;

  AliRICHDigit  *dig;
  dig=(AliRICHDigit*)RICHdigits->UncheckedAt(fDigitIndex);
  Int_t ipx=dig->fPadX;
  Int_t ipy=dig->fPadY;
  Bin bins[MAX_nipx][MAX_nipy]; 
  bins[ipx+npx][ipy+npy].dig=dig;
    
  int ndig;
  int ncls=0;
  for (ndig=0; ndig<ndigits; ndig++) {
      dig = (AliRICHDigit*)RICHdigits->UncheckedAt(ndig);
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
  AliRICHPoints *points = 0;
  points = new AliRICHPoints(1);
  points->SetMarkerColor(kYellow);
  points->SetMarkerStyle(5);
  points->SetMarkerSize(1.);
  points->SetPoint(0,c.fX,c.fY,zpos);
  points->Draw();
  
  printf("GetCenterOfGravity -- ncls %d \n",ncls);

}

//_____________________________________________________________________________
const Text_t *AliRICHPoints::GetName() const
{
  //
  // Return name of the Geant3 particle corresponding to this point
  //
  TParticle *particle = GetParticle();
  if (!particle) return "Particle";
  return particle->GetName();
}

//_____________________________________________________________________________
Text_t *AliRICHPoints::GetObjectInfo(Int_t, Int_t)
{
  //
  //   Redefines TObject::GetObjectInfo.
  //   Displays the info (particle,etc
  //   corresponding to cursor position px,py
  //
  static char info[64];
  sprintf(info,"%s %d",GetName(),fIndex);
  return info;
}



