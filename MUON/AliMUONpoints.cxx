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
<img src="gif/AliMUONpointsClass.gif">
*/
//End_Html
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliMUONpoints.h"
#include "AliMUONdisplay.h"
#include "AliRun.h"
#include "TPad.h"
#include "TVirtualPad.h"
#include "TPolyLine3D.h"
#include "TPaveText.h"
#include "TView.h"
#include "TMath.h"

//const Int_t MAX_Nipx=1026, MAX_Nipy=1026;
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
  fMarker[0] = fMarker[1] = fMarker[2]=0;
  fMatrix = 0;
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
  fMarker[0] = fMarker[1] = fMarker[2]=0;
  fMatrix = 0;
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
  for (Int_t i=0;i<3;i++){
      if (fMarker[i]) delete fMarker[i];
  }
  fMatrix = 0;
}

//_____________________________________________________________________________
//void AliMUONpoints::ExecuteEvent(Int_t event, Int_t px, Int_t py)
//{
  //
  //*-*-*-*-*-*-*-*-*-*Execute action corresponding to one event*-*-*-*-*-*-*-*
  //*-*                =========================================
  //*-*
  //*-*  This member function must be implemented to realize the action
  //*-*  corresponding to the mouse click on the object in the window
  //*-*
  //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//  gPad->SetCursor(kCross);
  
//}
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

  if (fHitIndex < 0 ) return;
  TVirtualPad *padsav = gPad;
  AliMUONhit *hit = GetHit();
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
void AliMUONpoints::InspectDigit()
{
  //
  //   Inspect digit corresponding to this point
  //
  if (fDigitIndex < 0) return;
  TVirtualPad *padsav = gPad;
  AliMUONdigit *digit = GetDigit();
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
Int_t AliMUONpoints::GetTrackIndex()
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
  Int_t nent=(Int_t)gAlice->TreeD()->GetEntries();
  gAlice->TreeD()->GetEvent(nent-2+cathode-1);
  //gAlice->TreeD()->GetEvent(cathode);
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
   int npoly;
   float xpoly[100];
   float ypoly[100];
   float zpoly[100];
   PreCluster() : AliMUONreccluster() {cut=npeaks=npoly=0; 
                                       for (int k=0;k<100;k++) {
                                         xpoly[k]=ypoly[k]=zpoly[k]=0;
				       }
   }
                              
};
//_____________________________________________________________________________

static void FindCluster(AliMUONchamber *iChamber, AliMUONsegmentation *segmentation, int i, int j, Bin bins[][MAX_Nipy], PreCluster &c) 

{

  //
  // Find clusters
  //


  Bin& b=bins[i][j];
  Int_t q=b.dig->fSignal;

  if (q<0) { 
    q=-q;
    c.cut=1;
  } 
  //  if (b.idx >= 0 && b.idx != c.idx) {
  if (b.idx >= 0 && b.idx > c.idx) {
    c.idx=b.idx;
    c.npeaks++;
  }
  
  if (q > TMath::Abs(c.summit->fSignal)) c.summit=b.dig;

  Float_t zpos=iChamber->ZPosition();

  Int_t npx  = segmentation->Npx();
  Int_t npy  = segmentation->Npy();


  // get pad coordinates and prepare the up and down steps   
  Int_t jup  =(j-npy > 0) ? j+1 : j-1;
  Int_t jdown=(j-npy > 0) ? j-1 : j+1;

  Float_t x, y;
  segmentation->GetPadCxy(i-npx, j-npy,x, y);
  Int_t isec0=segmentation->Sector(i-npx,j-npy);

  Float_t dpy  = segmentation->Dpy(isec0);
  Float_t dpx  = segmentation->Dpx(isec0)/16;
  Int_t ixx, iyy;
  Float_t absx=TMath::Abs(x);
  // iup
  (y >0) ? segmentation->GetPadIxy(absx+dpx,y+dpy,ixx,iyy) : segmentation->GetPadIxy(absx+dpx,y-dpy,ixx,iyy);

  //  Int_t jtest=TMath::Abs(iyy)-npy-1;

  Int_t iup=(x >0) ? ixx+npx : -ixx+npx;
  // idown
  (y >0) ? segmentation->GetPadIxy(absx+dpx,y-dpy,ixx,iyy) : segmentation->GetPadIxy(absx+dpx,y+dpy,ixx,iyy);

  Int_t idown=(x >0) ? ixx+npx : -ixx+npx;
  if (bins[idown][jdown].dig == 0) {
     (y >0) ? segmentation->GetPadIxy(absx-dpx,y-dpy,ixx,iyy) : segmentation->GetPadIxy(absx-dpx,y+dpy,ixx,iyy);
     idown=(x >0) ? ixx+npx : -ixx+npx;
  }
 

  // calculate center of gravity
  c.npoly++;
  if (c.npoly > 100 ) {
    printf("FindCluster - npoly >100,  npoly %d \n",c.npoly);
    c.npoly=100;
  }
  c.xpoly[c.npoly-1]=x;
  c.ypoly[c.npoly-1]=y;
  c.zpoly[c.npoly-1]=zpos;

  c.fX += q*x;
  c.fY += q*y;
  c.fQ += q;
  
  b.dig = 0;  b.idx = c.idx;

  // left and right  
  if (bins[i-1][j].dig) FindCluster(iChamber,segmentation,i-1,j,bins,c);
  if (bins[i+1][j].dig) FindCluster(iChamber,segmentation,i+1,j,bins,c);
  // up and down
  if (bins[iup][jup].dig) FindCluster(iChamber,segmentation,iup,jup,bins,c);
  if (bins[idown][jdown].dig) FindCluster(iChamber,segmentation,idown,jdown,bins,c);

}

//_____________________________________________________________________________

void AliMUONpoints::GetCenterOfGravity()
{
  //
  // simple MUON cluster finder from digits -- finds neighbours and 
  // calculates center of gravity for the cluster
  //

  //  const Int_t MAX_Nipx=1026, MAX_Nipy=1026;

  Bin bins[MAX_Nipx][MAX_Nipy]; 

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
  Int_t nent=(Int_t)gAlice->TreeD()->GetEntries();
  gAlice->TreeD()->GetEvent(nent-2+cathode-1);
  //gAlice->TreeD()->GetEvent(cathode);
  Int_t ndigits = MUONdigits->GetEntriesFast();
  if (fDigitIndex < 0 || fDigitIndex >= ndigits) return;

  AliMUONdigit  *dig;
  dig=(AliMUONdigit*)MUONdigits->UncheckedAt(fDigitIndex);
  Int_t ipx=dig->fPadX;
  Int_t ipy=dig->fPadY;
  bins[ipx+npx][ipy+npy].dig=dig;
    
  int ndig;
  int ncls=0;
  for (ndig=0;ndig<ndigits;ndig++) {
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
   points->SetParticle(-1);
   points->Draw();

   TPolyLine3D *pline=0;
   /*
   pline=new TPolyLine3D(c.npoly);
   Int_t np=c.npoly;
   TVector *xp=new TVector(c.npoly);
   TVector *yp=new TVector(c.npoly);
   TVector *zp=new TVector(c.npoly);
   for (int i=0;i<np;i++) {
     (*xp)(i)=c.xpoly[i];
     (*yp)(i)=c.ypoly[i];
     (*zp)(i)=c.zpoly[i];
     pline->SetPoint(i,(*xp)(i),(*yp)(i),(*zp)(i));
     //printf("np, i, xp, yp, zp %d %d %f %f %f \n",np,i,(*xp)(i),(*yp)(i),(*zp)(i));
   }
   */
   pline=new TPolyLine3D(c.npoly,c.xpoly,c.ypoly,c.zpoly);
   pline->SetLineColor(kWhite);
   pline->Draw();
   /*
   delete xp;
   delete yp;
   delete zp;
   */  
   for (int k=0;k<c.npoly;k++) {
     c.xpoly[k]=c.ypoly[k]=c.zpoly[k]=0;
   }
   c.npoly=0;


}

