///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  This class contains the points for the ALICE event display               //
//                                                                           //
//Begin_Html
/*
<img src="gif/AliMUONpointsClass.gif">
*/
//End_Html
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
#include "AliMUONdisplay.h"
#include "AliMUONpoints.h"
#include "AliRun.h"
#include "TPad.h"
#include "TView.h"
#include "TMath.h"

const Int_t MAX_Nipx=400, MAX_Nipy=800;
 
ClassImp(AliMUONpoints)

//_____________________________________________________________________________
AliMUONpoints::AliMUONpoints()
{
  //
  // Default constructor
  //
  fHitIndex = 0;
  fTrackIndex = 0;
  fDigitIndex = 0;
}

//_____________________________________________________________________________
AliMUONpoints::AliMUONpoints(Int_t npoints)
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
AliMUONpoints::~AliMUONpoints()
{
  //
  // Default destructor
  //
  fHitIndex = 0;
  fTrackIndex = 0;
  fDigitIndex = 0;
}

//_____________________________________________________________________________
void AliMUONpoints::DumpHit()
{
  //
  //   Dump hit corresponding to this point
  //
  AliMUONhit *hit = GetHit();
  if (hit) hit->Dump();
}

//_____________________________________________________________________________
void AliMUONpoints::DumpDigit()
{
  //
  //   Dump digit corresponding to this point
  //
  AliMUONdigit *digit = GetDigit();
  if (digit) digit->Dump();
}

//_____________________________________________________________________________
void AliMUONpoints::InspectHit()
{
  //
  //   Inspect hit corresponding to this point
  //
  AliMUONhit *hit = GetHit();
  if (hit) hit->Inspect();
}

//_____________________________________________________________________________
void AliMUONpoints::InspectDigit()
{
  //
  //   Inspect digit corresponding to this point
  //
  AliMUONdigit *digit = GetDigit();
  if (digit) digit->Inspect();
}

//_____________________________________________________________________________
Int_t AliMUONpoints::GetTrackIndex()
{
  //
  //   Dump digit corresponding to this point
  //
  printf("GetTrackIndex - fTrackIndex %d \n",fTrackIndex);
  this->Inspect();
  return fTrackIndex;
}

//_____________________________________________________________________________
AliMUONhit *AliMUONpoints::GetHit() const
{
  //
  //   Returns pointer to hit index in AliRun::fParticles
  //
  AliMUON *MUON  = (AliMUON*)gAlice->GetModule("MUON");
  gAlice->TreeH()->GetEvent(fTrackIndex);
  TClonesArray *MUONhits  = MUON->Hits();
  Int_t nhits = MUONhits->GetEntriesFast();
  if (fHitIndex < 0 || fHitIndex >= nhits) return 0;
  return (AliMUONhit*)MUONhits->UncheckedAt(fHitIndex);
}

//_____________________________________________________________________________
AliMUONdigit *AliMUONpoints::GetDigit() const
{
  //
  //   Returns pointer to digit index in AliRun::fParticles
  //

  AliMUONdisplay *display=(AliMUONdisplay*)gAlice->Display();
  Int_t chamber=display->GetChamber();
  Int_t cathode=display->GetCathode();
   
  AliMUON *MUON  = (AliMUON*)gAlice->GetModule("MUON");
  TClonesArray *MUONdigits  = MUON->DigitsAddress(chamber-1);
  gAlice->TreeD()->GetEvent(cathode);
  Int_t ndigits = MUONdigits->GetEntriesFast();
  if (fDigitIndex < 0 || fDigitIndex >= ndigits) return 0;
  return (AliMUONdigit*)MUONdigits->UncheckedAt(fDigitIndex);
}
//_____________________________________________________________________________
struct Bin {
   const AliMUONdigit *dig;
   int idx;
   Bin() {dig=0; idx=-1;}
};

struct PreCluster : public AliMUONreccluster {
   const AliMUONdigit* summit;
   int idx;
   int cut;
   int npeaks;
   PreCluster() : AliMUONreccluster() {cut=npeaks=0;}
};
//_____________________________________________________________________________

static void FindCluster(AliMUONchamber *iChamber, AliMUONsegmentation *segmentation, int i, int j, Bin bins[MAX_Nipx][MAX_Nipy], PreCluster &c) 

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

void AliMUONpoints::GetCenterOfGravity()
{
  //
  // simple MUON cluster finder from digits -- finds neighbours and 
  // calculates center of gravity for the cluster
  //
  const Int_t MAX_Nipx=400, MAX_Nipy=800;
  printf("\n Hallo world");
  AliMUONdisplay *display=(AliMUONdisplay*)gAlice->Display();
  Int_t chamber=display->GetChamber();
  Int_t cathode=display->GetCathode();
   
  AliMUON *MUON  = (AliMUON*)gAlice->GetModule("MUON");
  AliMUONchamber *iChamber;
  AliMUONsegmentation *segmentation;
  iChamber =&(MUON->Chamber(chamber-1));
  segmentation=iChamber->GetSegmentationModel(cathode);
  Int_t npx  = segmentation->Npx();
  Int_t npy  = segmentation->Npy();
  Float_t zpos=iChamber->ZPosition();
  
  TClonesArray *MUONdigits  = MUON->DigitsAddress(chamber-1);
  gAlice->TreeD()->GetEvent(cathode);
  Int_t ndigits = MUONdigits->GetEntriesFast();
  if (fDigitIndex < 0 || fDigitIndex >= ndigits) return;

  AliMUONdigit  *dig;
  dig=(AliMUONdigit*)MUONdigits->UncheckedAt(fDigitIndex);
  Int_t ipx=dig->fPadX;
  Int_t ipy=dig->fPadY;
  Bin bins[MAX_Nipx][MAX_Nipy]; 
  bins[ipx+npx][ipy+npy].dig=dig;
    
  int ndig;
  int ncls=0;
  for (ndig=0; ndig<ndigits; ndig++) {
      dig = (AliMUONdigit*)MUONdigits->UncheckedAt(ndig);
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
  AliMUONpoints *points = 0;
  points = new AliMUONpoints(1);
  points->SetMarkerColor(kYellow);
  points->SetMarkerStyle(5);
  points->SetMarkerSize(1.);
  points->SetPoint(0,c.fX,c.fY,zpos);
  points->Draw();
  
  printf("GetCenterOfGravity -- ncls %d \n",ncls);

}



