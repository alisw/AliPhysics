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


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  This class contains the points for the ALICE event display               //
//                                                                           //
//Begin_Html
/*
<img src="gif/AliITSpointsClass.gif">
*/
//End_Html
//                                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TPad.h>
#include <TVirtualPad.h>
#include <TPolyLine3D.h>
#include <TPaveText.h>
#include <TView.h>
#include <TMath.h>
#include <TH1.h>
#include <TF1.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TSpectrum.h>

#include "AliRun.h"
#include "AliITSpoints.h"
#include "AliITSdisplay.h"
#include "AliITSMap.h"
#include "AliITSRecPointSDDnew.h"


//const Int_t MAXNipx=1026, MAXNipy=1026;
const Int_t MAXNipx=752, MAXNipy=256;

static AliITSMapA2 *sMap = 0; 
static Int_t  sModule=0; 
 
ClassImp(AliITSpoints)

//_____________________________________________________________________________
AliITSpoints::AliITSpoints()
{
  //
  // Default constructor
  //
  fModuleIndex = -1;
  fHitIndex = -1;
  fTrackIndex = -1;
  fDigitIndex = -1;
  fMarker[0] = fMarker[1] = fMarker[2]=0;
  fMatrix = 0;

  //fConnect=kFALSE;
}

//_____________________________________________________________________________
AliITSpoints::AliITSpoints(Int_t npoints)
  :AliPoints(npoints)
{
  //
  // Standard constructor
  //
  fModuleIndex = -1;
  fHitIndex = -1;
  fTrackIndex = -1;
  fDigitIndex = -1;
  fMarker[0] = fMarker[1] = fMarker[2]=0;
  fMatrix = 0;

  //fConnect=kFALSE;

}
	 
//_____________________________________________________________________________
AliITSpoints::~AliITSpoints()
{
  //
  // Default destructor
  //
  fHitIndex = 0;
  fTrackIndex = 0;
  fDigitIndex = 0;
  fModuleIndex = 0;
  for (Int_t i=0;i<3;i++){
      if (fMarker[i]) delete fMarker[i];
  }
  fMatrix = 0;
}

//_____________________________________________________________________________
void AliITSpoints::ExecuteEvent(Int_t event, Int_t px, Int_t py)
{
  //
  //*-*-*-*-*-*-*-*-*-*Execute action corresponding to one event*-*-*-*-*-*-*-*
  //*-*                =========================================
  //*-*
  //*-*  This member function must be implemented to realize the action
  //*-*  corresponding to the mouse click on the object in the window
  //*-*
  //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

  if (!gPad) return;

  gPad->SetCursor(kCross);

  TObject *select = gPad->GetSelected();
  if(!select) {printf("no select \n"); return;}
  if (!select->InheritsFrom("AliITSpoints")) {gPad->SetUniqueID(0); return;}

   //erase old position and draw a line at current position
   switch (event) {


   case kButton1Down:
        return;
   case kButton1Motion:
       AnodeProjection(px,py);
       return;
   case kButton1Up:
       return;


   case kButton1Double:
       TimeProjection(px,py);
       return;

       /* why button2 does not work ?
        it does not know kButton2Down, kButton2Motion, kButton2Up
   case kButton2Down:
       printf("Button2Down \n");
       return;
   case kButton2Motion:
       printf("Button2Motion \n");
       TimeProjection(px,py);
       return;
   case kButton2Up:
       return;
   */
   }
   
}

//_____________________________________________________________________________
void AliITSpoints::AnodeProjection(Int_t px, Int_t py) 
{
  //if(!fDrawHist) return;

    gPad->SetCursor(kMove);
    gPad->GetCanvas()->FeedbackMode(kTRUE);

   Float_t xmin = gPad->GetX1();
   Float_t xmax = gPad->GetX2();
   Float_t ymin = gPad->GetY1();
   Float_t ymax = gPad->GetY2();
   //printf("xmin,xmax,ymin,ymax %f %f %f %f\n",xmin,xmax,ymin,ymax);
    Float_t x,y;
   // x,y in the -1,1 scale
    x = gPad->AbsPixeltoX(px);
    y = gPad->AbsPixeltoY(py);

    Float_t scale[3],center[3];
    Int_t irep;
    gPad->GetView()->FindScope(scale,center,irep);

    AliITSdisplay* display=((AliITSdisplay*)(gAlice->Display()));
     Int_t anode, timebin;
    display->GetPadIxy(anode,timebin,x*scale[0],y*scale[0]);
    
    Int_t amin, amax, tmin, tmax;
    display->GetPadIxy(amin,tmin,xmin*scale[0],ymin*scale[0]);
    display->GetPadIxy(amax,tmax,xmax*scale[0],ymax*scale[0]);

    if (xmin < 0 && xmax > 0) {
       if (0-xmin > xmax) amax-=374;
       else {amin+=374; tmin=tmax;}
       tmax=256;
    } 
    Int_t tminold=tmin;
    Int_t tmaxold=tmax;
    tmin=TMath::Min(tminold,tmaxold);
    tmax=TMath::Max(tminold,tmaxold);

    Int_t nbiny = tmax-tmin;

    //printf("amin amax tmin tmax anode, timebin %d %d %d %d %d %d\n",amin,amax,tmin,tmax,anode,timebin);

    //create or set the new canvas c2
   TVirtualPad *padsav = gPad;
   //printf("AnodeProj: gPad %p\n",gPad);
   TCanvas *c2 = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("c2");
   if(c2) delete c2->GetPrimitive("Anode Projection");
   if(c2) delete c2->GetPrimitive("Time Projection");
   else     c2=new TCanvas("c2");
   c2->cd();
   c2->SetFillColor(0);

   char title[30];
   sprintf(title,"anode=%d",anode);

   TH1F *h1=new TH1F("Anode Projection","",nbiny,(float)tmin,(float)tmax);

   h1->SetTitle(title);
   h1->SetTitleSize(2.);
   h1->SetTitleOffset(2.);
   h1->SetLabelSize(0.05);
   h1->SetYTitle(title);
   //h1->SetStats(0);

   AliITSMapA2 *sMap=display->GetMap();

   for (int j=tmin;j<tmax;j++) {
     //printf("an tbin signal %d %d %f\n",anode,j,(float)sMap->GetSignal(anode-1,j-1));
      h1->Fill((float)j,(float)sMap->GetSignal(anode-1,j-1));
   }
   h1->Smooth(1);
   h1->Fit("gaus","ql");
   h1->Draw();

   c2->Update();
   padsav->cd();
   //printf("AnodeProj: gPad padsav %p %p\n",gPad,padsav);

}

//_____________________________________________________________________________
void AliITSpoints::TimeProjection(Int_t px, Int_t py)
{
  //if(!fDrawHist) return;
  //printf("TimeProjection \n");

    gPad->SetCursor(kRightSide);
    gPad->GetCanvas()->FeedbackMode(kTRUE);

   Float_t xmin = gPad->GetX1();
   Float_t xmax = gPad->GetX2();
   Float_t ymin = gPad->GetY1();
   Float_t ymax = gPad->GetY2();
   //printf("xmin,xmax,ymin,ymax %f %f %f %f\n",xmin,xmax,ymin,ymax);
    Float_t x,y;
   // x,y in the -1,1 scale
    x = gPad->AbsPixeltoX(px);
    y = gPad->AbsPixeltoY(py);

    Float_t scale[3],center[3];
    Int_t irep;
    gPad->GetView()->FindScope(scale,center,irep);

    AliITSdisplay* display=((AliITSdisplay*)(gAlice->Display()));
     Int_t anode, timebin;
    display->GetPadIxy(anode,timebin,x*scale[0],y*scale[2]);
    
    Int_t amin, amax, tmin, tmax;
    display->GetPadIxy(amin,tmin,xmin*scale[0],ymin*scale[2]);
    display->GetPadIxy(amax,tmax,xmax*scale[0],ymax*scale[2]);

    if (xmin < 0 && xmax > 0) {
       if (0-xmin > xmax) amax-=374;
       else {amin+=374;tmin=tmax;}
       tmax=256;
    } 

    Int_t aminold=amin;
    Int_t amaxold=amax;
    amin=TMath::Min(aminold,amaxold);
    amax=TMath::Max(aminold,amaxold);

    Int_t nbinx = amax-amin;

    //printf("amin amax tmin tmax anode, timebin %d %d %d %d %d %d\n",amin,amax,tmin,tmax,anode,timebin);

    //create or set the new canvas c2
   TVirtualPad *padsav = gPad;
   TCanvas *c2 = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("c2");
   if(c2) delete c2->GetPrimitive("Anode Projection");
   if(c2) delete c2->GetPrimitive("Time Projection");
   else     c2=new TCanvas("c2");
   c2->cd();
   c2->SetFillColor(0);

   char title[30];
   //printf("title: anode, timebin %d %d\n",anode,timebin);
   sprintf(title,"time sample=%d",timebin);

   TH1F *h2=new TH1F("Time Projection","",nbinx,(float)amin,(float)amax);

   h2->SetTitle(title);
   h2->SetTitleSize(2.);
   h2->SetTitleOffset(2.);
   h2->SetLabelSize(0.05);
   h2->SetYTitle(title);
   //h2->SetStats(0);

   AliITSMapA2 *sMap=display->GetMap();

   for (int i=amin;i<amax;i++) {
     //printf("an tbin signal %d %d %f\n",i,timebin,(float)sMap->GetSignal(i-1,timebin-1));
      h2->Fill((float)i,(float)sMap->GetSignal(i-1,timebin-1));
   }
   h2->Smooth(1);
   h2->Fit("gaus","ql");
   h2->Draw();

   c2->Update();
   padsav->cd();

}
//_____________________________________________________________________________
void AliITSpoints::DisplayModule()
{
  /*
 
   if (!((AliITSDisplay*)(gAlice->Display()))->GetDrawTracksOpt()) return;
   // ????
   ((AliITSDisplay*)(gAlice->Display()))->SetModuleNumber(fModuleIndex);
   Int_t event=((AliITSDisplay*)(gAlice->Display()))->GetEvent();
   ((AliITSDisplay*)(gAlice->Display()))->SetEvent(event);
   ((AliITSDisplay*)(gAlice->Display()))->SetRange(4,4);
   //((AliITSDisplay*)(gAlice->Display()))->DrawModule(fModuleIndex);
  */
}

//_____________________________________________________________________________
const Text_t *AliITSpoints::GetName() const
{
  //
  // Return name of the Geant3 particle corresponding to this point
  //
  TParticle *particle = GetParticle();
  if (particle) return particle->GetName();
  else  return IsA()->GetName();
  //if (!particle) return "Particle";
}

//_____________________________________________________________________________
Text_t *AliITSpoints::GetObjectInfo(Int_t px, Int_t py)
{
  //
  //   Redefines TObject::GetObjectInfo.
  //   Displays the info (particle,etc
  //   corresponding to cursor position px,py
  //
   if (!gPad) return (char*)"";
  static char info[64];
  char an[6], tbin[9], track[6];
  an="Anode"; tbin="Time bin"; track="Track"; 
  if(strcmp(GetName(),"AliITSpoints ")) {
        AliITSdigitSDD *dig=GetDigit();
        if (!dig) {sprintf(info,"%s %s %d",GetName(),track,fIndex); return info;} 
        int anode=dig->fCellX;
        int time=dig->fCellY;
        sprintf(info,"%s %d %s %d %s %d",an,anode+1,tbin,time+1,track,fIndex);
  } else {
        if(strcmp(GetName(),"TView ")) {
            Float_t x = gPad->AbsPixeltoX(px);
            Float_t y = gPad->AbsPixeltoY(py);
            sprintf(info,"%s x=%.3g, y=%.3g",GetName(),gPad->PadtoX(x),gPad->PadtoY(y));
	} else {
            sprintf(info,"%s %s %d",GetName(),track,fIndex); 
	}

  }

  return info;
}

//_____________________________________________________________________________
void AliITSpoints::DumpHit()
{
  //
  //   Dump hit corresponding to this point
  //
  AliITShit *hit = GetHit();
  if (hit) hit->Dump();

}

//_____________________________________________________________________________
void AliITSpoints::DumpDigit()
{
  //
  //   Dump digit corresponding to this point
  //
  AliITSdigitSDD *digit = GetDigit();
  if (digit) digit->Dump();
}

//_____________________________________________________________________________
void AliITSpoints::InspectHit()
{
  //
  //   Inspect hit corresponding to this point
  //

  if (fHitIndex < 0 ) return;
  TVirtualPad *padsav = gPad;
  AliITShit *hit = GetHit();
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
void AliITSpoints::InspectDigit()
{
  //
  //   Inspect digit corresponding to this point
  //

  if (fDigitIndex < 0) return;
  TVirtualPad *padsav = gPad;
  AliITSdigitSDD *digit = GetDigit();
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
  for (int i=0;i<3;i++) {
      if (digit->fTracks[i] == 0) continue;  
      sprintf(ptitle[i],"fTrackIndex: %d  Charge: %f",digit->fTracks[i],digit->fTcharges[i]);
      pad->AddText(ptitle[i]);
  }
      padinspect->cd();
      padinspect->Update();
  if (padsav) padsav->cd();
    
}

//_____________________________________________________________________________
Int_t AliITSpoints::GetTrackIndex()
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
AliITShit *AliITSpoints::GetHit() const
{
  //
  //   Returns pointer to hit index in AliRun::fParticles
  //
  AliITS *ITS  = (AliITS*)gAlice->GetModule("ITS");
  gAlice->TreeH()->GetEvent(fTrackIndex);
  TClonesArray *ITShits  = ITS->Hits();
  Int_t nhits = ITShits->GetEntriesFast();
  if (fHitIndex < 0 || fHitIndex >= nhits) return 0;
  return (AliITShit*)ITShits->UncheckedAt(fHitIndex);
}

//_____________________________________________________________________________
AliITSdigitSDD *AliITSpoints::GetDigit() const
{
  //
  //   Returns pointer to digit index in AliRun::fParticles
  //

  AliITSdisplay *display=(AliITSdisplay*)gAlice->Display();
  Int_t module=display->GetModule();
  Int_t layer=display->GetLayer();
  Int_t id=(Int_t)((layer-1)/2);
   
  AliITS *ITS  = (AliITS*)gAlice->GetModule("ITS");
  TClonesArray *ITSdigits  = ITS->DigitsAddress(id);
  gAlice->TreeD()->GetEvent(module);
  Int_t ndigits = ITSdigits->GetEntriesFast();
  if (fDigitIndex < 0 || fDigitIndex >= ndigits) return 0;
  return (AliITSdigitSDD*)ITSdigits->UncheckedAt(fDigitIndex);

  /*
  // have smth like ??
  AliITSgeom *geom  = ITS->GetITSgeom;
  Int_t lastSPD=geom->GetLastSPD();
  Int_t lastSDD=geom->GetLastSDD();

  Int_t nent=(Int_t)gAlice->TreeD()->GetEntries();
  if (nent < lastSDD) {
       // take the appropriate action ...using lastSPD, lastSDD ...
       // or better introduce an Option_t in the display : All, SPD, SDD, SSD
       // or combination
       printf("Digitisation wasn't done for all det types !\n");
  }
  */
  // have smth like 
  // if (module < lastSPD) return (AliITSdigitSPD*)ITSdigits->UncheckedAt(fDigitIndex)
  // else  (module < lastSDD) return (AliITSdigitSDD*)ITSdigits->UncheckedAt(fDigitIndex)  .... ??

  //return 0;
}

//_____________________________________________________________________________


struct Bin {
   AliITSdigitSDD *dig;
   int idx;
   int digidx;
   Bin() {dig=0; idx=-1; digidx=-1;}
};


struct PreCluster : public AliITSRecPointSDDnew {
   AliITSdigitSDD* summit;
   int idx;   
   int npeaks;
   int npoly;
   float xpoly[300];
   float ypoly[300];
   float zpoly[300];
   PreCluster() : AliITSRecPointSDDnew() {npeaks=npoly=0; 
                                       for (int k=0;k<300;k++) {
                                         xpoly[k]=ypoly[k]=zpoly[k]=0;
				       }
   }
                              
};


//_____________________________________________________________________________

//static void FindCluster(int i, int j, Bin bins[][MAXNipy], PreCluster &c) 
static void FindCluster(int i, int j, Bin *bins, PreCluster &c, int thresh) 

{

  //
  // Find clusters
  //

  //Bin& b=bins[i][j];
  Bin& b=bins[i*MAXNipy+j];
  Int_t q=b.dig->fSignal - thresh;

  printf("i j q %d %d %d\n",i,j,q);

  //  if (b.idx >= 0 && b.idx != c.idx) {
  if (b.idx >= 0 && b.idx > c.idx) {
    c.idx=b.idx;
    c.npeaks++;
  }
  
  if (q > TMath::Abs(c.summit->fSignal)) c.summit=b.dig;

  // get pad coordinates and prepare the up and down steps   

  Float_t xl[3], dx[3];
  AliITSdisplay *display=(AliITSdisplay*)gAlice->Display();
  display->GetPadCxy(i,j,xl,dx);

  // calculate center of gravity
  c.npoly++;
  c.fMulDigits++;
  c.fDigitsList[c.npoly-1]=b.digidx;

  //printf("npoly c.fDigitsList[c.npoly-1] %d %d \n",c.npoly,c.fDigitsList[c.npoly-1]);

  if (c.npoly > 300 ) {
    printf("FindCluster - npoly >300,  npoly %d \n",c.npoly);
    c.npoly=300;
  }
  c.xpoly[c.npoly-1]=xl[0];
  c.ypoly[c.npoly-1]=xl[1];
  c.zpoly[c.npoly-1]=xl[2];

  c.fX += q*xl[0];
  c.fZ += q*xl[2];
  c.fY =  xl[1];
  c.fQ += q;
  
  b.dig = 0;  b.idx = c.idx;

  /*
  // left and right  
  if (i && bins[i-1][j].dig) FindCluster(i-1,j,bins,c);
  if (i < MAXNipx && bins[i+1][j].dig) FindCluster(i+1,j,bins,c);
  // up and down
  if (j+1 < MAXNipy && bins[i][j+1].dig) FindCluster(i,j+1,bins,c);
  if (j && bins[i][j-1].dig) FindCluster(i,j-1,bins,c);
  */
  // left and right  
  if (i && bins[(i-1)*MAXNipy+j].dig) FindCluster(i-1,j,bins,c,thresh);
  if (i < MAXNipx && bins[(i+1)*MAXNipy+j].dig) FindCluster(i+1,j,bins,c,thresh);
  // up and down
  if (j+1 < MAXNipy && bins[i*MAXNipy+(j+1)].dig) FindCluster(i,j+1,bins,c,thresh);
  if (j && bins[i*MAXNipy+(j-1)].dig) FindCluster(i,j-1,bins,c,thresh);

}

//_____________________________________________________________________________
void AliITSpoints::GetCenterOfGravity()
{
  //
  // simple ITS cluster finder from digits -- finds neighbours and 
  // calculates center of gravity for the cluster
  //

  const int THRESHOLD=20;

  Bin *sBin=new Bin[MAXNipx*MAXNipy];
  //else memset(sBin,0,sizeof(Bin)*MAXNipx*MAXNipy);
  //printf("sBin %p\n",sBin);


  //Bin bins[MAXNipx][MAXNipy]; 

  AliITSdisplay *display=(AliITSdisplay*)gAlice->Display();
  Int_t module=display->GetModule(); // or module=fModuleIndex;
  Int_t layer=display->GetLayer();
  Int_t id=(Int_t)((layer-1)/2);
   
  AliITS *ITS  = (AliITS*)gAlice->GetModule("ITS");

  TClonesArray *ITSdigits  = ITS->DigitsAddress(id);
  Int_t nent=(Int_t)gAlice->TreeD()->GetEntries();
  //gAlice->TreeD()->GetEvent(nent-nofmodules+module-1);
  gAlice->TreeD()->GetEvent(module);
  Int_t ndigits = ITSdigits->GetEntriesFast();
  if (fDigitIndex < 0 || fDigitIndex >= ndigits) return;

  AliITSdigitSDD  *dig;
  dig=(AliITSdigitSDD*)ITSdigits->UncheckedAt(fDigitIndex);
  Int_t ipx=dig->fCellX;
  Int_t ipy=dig->fCellY;
  //bins[ipx][ipy].dig=dig;
  sBin[ipx*MAXNipy+ipy].dig=dig;
  //printf("ipx, ipy, Sbin.dig %d %d %p\n",ipx,ipy,sBin[ipx*MAXNipy+ipy].dig);
    
  int ndig;
  int ncls=0;
  for (ndig=0;ndig<ndigits;ndig++) {
    dig = (AliITSdigitSDD*)ITSdigits->UncheckedAt(ndig);
    int i=dig->fCellX; int j=dig->fCellY;
    //printf("ndig i j %d %d %d \n",ndig,i,j);
    //bins[i][j].dig=dig;
    sBin[i*MAXNipy+j].dig=dig;
    sBin[i*MAXNipy+j].digidx=ndig;
  }
  //PreCluster c; c.summit=bins[ipx][ipy].dig; c.idx=ncls;
  // FindCluster(ipx, ipy, bins, c);
   PreCluster c; c.summit=sBin[ipx*MAXNipy+ipy].dig; c.idx=ncls;
   FindCluster(ipx, ipy, sBin, c, 0);
    
   if (c.npeaks>1) {
      printf("GetCenterOfGravity -- more than one peak");
   }
   c.fX /= c.fQ;
   c.fZ /= c.fQ;
   ncls++;

   AliITSpoints *points = 0;
   points = new AliITSpoints(1);
   points->SetMarkerColor(kYellow);
   points->SetMarkerStyle(5);
   points->SetMarkerSize(1.);
   points->SetPoint(0,c.fX,c.fY,c.fZ);
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
   /*
   delete xp;
   delete yp;
   delete zp;
   */ 

   //if (fConnect) {
      pline=new TPolyLine3D(c.npoly,c.xpoly,c.ypoly,c.zpoly);
      pline->SetLineColor(kWhite);
      pline->Draw();
      //}



   for (int k=0;k<c.npoly;k++) {
     c.xpoly[k]=c.ypoly[k]=c.zpoly[k]=0;
   }
   c.npoly=0;

}

//_____________________________________________________________________________
void AliITSpoints::FindLocalMaxima()
{
  //

  Bin *lBin=new Bin[MAXNipx*MAXNipy];

  AliITSdisplay *display=(AliITSdisplay*)gAlice->Display();
  Int_t module=display->GetModule(); // or module=fModuleIndex;
  Int_t layer=display->GetLayer();
  Int_t id=(Int_t)((layer-1)/2);
   
  AliITS *ITS  = (AliITS*)gAlice->GetModule("ITS");

  TClonesArray *ITSdigits  = ITS->DigitsAddress(id);
  Int_t nent=(Int_t)gAlice->TreeD()->GetEntries();
  //gAlice->TreeD()->GetEvent(nent-nofmodules+module-1);
  gAlice->TreeD()->GetEvent(module);
  Int_t ndigits = ITSdigits->GetEntriesFast();
  if (fDigitIndex < 0 || fDigitIndex >= ndigits) return;

  AliITSdigitSDD  *dig;
  dig=(AliITSdigitSDD*)ITSdigits->UncheckedAt(fDigitIndex);
  Int_t ipx=dig->fCellX;
  Int_t ipy=dig->fCellY;
  lBin[ipx*MAXNipy+ipy].dig=dig;
    
  int ndig;
  int ncls=0;

  const int kMaxNeighbours=4;
  
  for (ndig=0;ndig<ndigits;ndig++) {
    dig = (AliITSdigitSDD*)ITSdigits->UncheckedAt(ndig);
    int i=dig->fCellX; int j=dig->fCellY;
    //printf("ndig i j %d %d %d \n",ndig,i,j);
    lBin[i*MAXNipy+j].dig=dig;
    lBin[i*MAXNipy+j].digidx=ndig;
  }
   PreCluster c; c.summit=lBin[ipx*MAXNipy+ipy].dig; c.idx=ncls;
   FindCluster(ipx, ipy, lBin, c, 0);
    
   printf("finished FindCluster !!!\n");

   Int_t fMul=c.npoly;
   printf("fMul %d \n",fMul);
   AliITSdigitSDD* digt;

  int fIx[fMul];
  int fIy[fMul];
  int fSortIx[fMul];
  int fSortIy[fMul];
  int fQ[fMul];
  int fQini[fMul];
  int fIndLocal[fMul];
  int fIndSaddle[fMul];
  AliITSdigitSDD*    fDig[fMul];        // current list of digits 

  for (int i=0; i<fMul; i++)
  {
        int idx=c.fDigitsList[i];
        digt = (AliITSdigitSDD*)ITSdigits->UncheckedAt(idx);
        fDig[i]=digt;
	fIx[i]= digt->fCellX;
	fIy[i]= digt->fCellY;
	fSortIx[i]= digt->fCellX;
	fSortIy[i]= digt->fCellY;
	fQini[i] = digt->fSignal;
  }

    Int_t idx_sort[fMul];

    TMath::Sort(fMul,fSortIy,idx_sort,0);
    Int_t tmin=fIy[idx_sort[0]];
    Int_t tmax=fIy[idx_sort[fMul-1]];

    TMath::Sort(fMul,fSortIx,idx_sort,0);
    Int_t amin=fIx[idx_sort[0]];
    Int_t amax=fIx[idx_sort[fMul-1]];
    //printf("tmin tmax amin amax %d %d %d %d\n",tmin,tmax,amin,amax);


   AliITSMapA2 *sMap=display->GetMap();

   Int_t nt=tmax-tmin;
   Int_t nbinx=amax-amin;

   float **source = new float* [nbinx+1];
   for(int i=0;i<nbinx+1;i++) source[i] = new float[nt+1];

   Double_t Qsmooth[nbinx+1][nt+1];
   Double_t temp[nt+1];
   TH1F *h1=new TH1F("Anode Projection","",nbinx+1,(float)amin,(float)amax+1);

   for (int i=amin;i<=amax;i++) {
     for (int j=tmin;j<=tmax;j++) {
         source[i-amin][j-tmin]=0;
         Qsmooth[i-amin][j-tmin]=0.;
         source[i-amin][j-tmin]=(float)sMap->GetSignal(i,j);
         Qsmooth[i-amin][j-tmin]=sMap->GetSignal(i,j);
	 //printf("i j Qsmooth signal map %d %d %f %f\n",i,j,Qsmooth[i-amin][j-tmin],sMap->GetSignal(i,j));
         temp[j-tmin]=sMap->GetSignal(i,j);
     }
     h1->SmoothArray(nt+1,temp,1);
     for (int j=tmin;j<=tmax;j++) {
         Qsmooth[i-amin][j-tmin]=temp[j-tmin];
	 //printf("i j Qsmooth temp %d %d %f %f\n",i,j,Qsmooth[i-amin][j-tmin],temp[j-tmin]);
     }

   }


//
//  Find local maxima//
    int fNLocal=0;
    int fNSaddle=0;
    Bool_t IsLocal[fMul];
    Bool_t IsSaddle[fMul];
    Int_t AssocPeak[fMul];
    Int_t nn;
    Int_t Qneighb[2*kMaxNeighbours];
    Int_t X[2*kMaxNeighbours], Y[2*kMaxNeighbours];
    for (int i=0; i<fMul; i++) {
        fQ[i]=(Int_t)Qsmooth[fIx[i]-amin][fIy[i]-tmin];
	Neighbours(fIx[i], fIy[i], &nn, X, Y);
	IsLocal[i]=kTRUE;
	IsSaddle[i]=kTRUE;
	for (int j=0; j<nn; j++) {
            Qneighb[j]=0;
            if (X[j]<amin || X[j]>amax || Y[j]<tmin || Y[j]>tmax) continue;
            Qneighb[j]=(Int_t)Qsmooth[X[j]-amin][Y[j]-tmin];
	}
        int nlzero=0;
	for (int j=0; j<nn; j++) {
	  //if (!fQ[i] || !Qneighb[j] || !fQini[i]) {nlzero++;continue;}
	    if (!fQ[i] || Qneighb[j]<4 || !fQini[i]) {nlzero++;continue;}
	    //if (Qneighb[j] > fQ[i] || nlzero>=2) { //???
	    if (Qneighb[j] > fQ[i] ) {
	      //printf("i j %d %d  local kFALSE Qneighb fQ %d %d\n",i,j,Qneighb[j],fQ[i]);
		IsLocal[i]=kFALSE;
		break;
//
// handle special case of neighbouring pads with equal signal
	     } else if (Qneighb[j] == fQ[i]) {
		if (fNLocal >0) {
		    for (Int_t k=0; k<fNLocal; k++) {
			if (X[j]==fIx[fIndLocal[k]] && (Y[j]==fIy[fIndLocal[k]] || Y[j]==fIy[fIndLocal[k]+1] || Y[j]==fIy[fIndLocal[k]-1]))
			{
			    IsLocal[i]=kFALSE;
			} 
                        if(j-2 >= tmin && j+2 <= tmax && (fQ[i] == Qsmooth[X[j]-amin][Y[j-2]-tmin]) || (fQ[i] == Qsmooth[X[j]-amin][Y[j+2]-tmin]))
			{
			    IsLocal[i]=kFALSE;
			} 
		    } // loop over local maxima
		} // are there are already local maxima
	    } // IsLocal
	} // loop over neighb

	// find saddle points
        int nzero=0; int nlower=0;
	for (int j=0; j<nn; j++) {
            if (X[j]<amin || X[j]>amax || Y[j]<tmin || Y[j]>tmax) continue;
            Qneighb[j]=(int)source[X[j]-amin][Y[j]-tmin];
	    if (!Qneighb[j] || !fQ[i] || !fQini[i]) {nzero++;continue;}
	    //if (!fQ[i] || Qneighb[j]<4 || !fQini[i]) {nzero++;continue;}
	    if (Qneighb[j] < fQ[i] ) {
	      //printf("fIx fIy X[j] Y[j] %d %d %d %d  saddle kFALSE Qneighb fQ %d %d\n",fIx[i],fIy[i],X[j],Y[j],Qneighb[j],fQ[i]);
              nlower++;
              
	      IsSaddle[i]=kFALSE;
	      break;
//
// handle special case of neighbouring pads with equal signal
	      // this is dangerous for saddle points ! better out !!
	    } else if (Qneighb[j] == fQ[i]) {
		if (fNSaddle >0) {
		    for (Int_t k=0; k<fNSaddle; k++) {
			if ((X[j]==fIx[fIndSaddle[k]] && (Y[j]==fIy[fIndSaddle[k]] || Y[j]==fIy[fIndSaddle[k]-1] || Y[j]==fIy[fIndSaddle[k]+1]))) 

			  //|| (X[j]==fIx[fIndSaddle[k]-1] && (Y[j]==fIy[fIndSaddle[k]] || Y[j]==fIy[fIndSaddle[k]-1] ||Y[j]==fIy[fIndSaddle[k]+1])) || (X[j]==fIx[fIndSaddle[k]+1] && (Y[j]==fIy[fIndSaddle[k]] || Y[j]==fIy[fIndSaddle[k]-1] ||Y[j]==fIy[fIndSaddle[k]+1])))
			{
			    IsSaddle[i]=kFALSE;
			} 

		    } // loop over saddle points
		} // are there are already saddle points
	    } // IsSaddle
	} // loop over next neighbours

	// Maxima should not be on the edge - therefore cond on nlzero ?
	if (IsLocal[i] && fQ[i]>0) {
	  //printf(" i fNLocal nlzero %d %d %d\n",i,fNLocal,nlzero);
	    fIndLocal[fNLocal]=i;
	    fNLocal++;
	} else fQ[i]=0;
	// Maxima should not be on the edge
	//if (IsSaddle[i] && fQ[i]>0 && !nzero ) {
	if (IsSaddle[i] && fQini[i]>0) {
	  //printf(" i fNSaddle nzero %d %d %d\n",i,fNSaddle,nzero);
          if (fIy[i]-1 < tmin || fIy[i]+1 > tmax) continue;
	  if (source[fIx[i]-amin][fIy[i]-1-tmin] && source[fIx[i]-amin][fIy[i]+1-tmin]) {  
	    fIndSaddle[fNSaddle]=i;
	    fNSaddle++;
	  }
	}
    } // loop over all digits
    printf("fNLocal %d  fNSaddle %d fMul %d \n",fNLocal,fNSaddle,fMul);


    // take the highest local maxima and spread the charge, i.e. define the
    // response 

    Int_t idx=TMath::LocMax(fMul,fQ);
    int thresh=0;
    int qsaddle[fNSaddle];
    if (fNSaddle) {
       for (int i=0;i<fNSaddle;i++) {
               qsaddle[i]=fQini[fIndSaddle[i]];
	       printf("isaddle qsaddle %d %d \n",i,qsaddle[i]);
       }
       Int_t idmin=TMath::LocMin(fNSaddle,qsaddle);
       thresh=qsaddle[idmin];
    }
    printf("idx fQ[idx] thresh %d %d %d\n",idx,fQ[idx],thresh);

    // this part should be recursive for fNLocal>2 and should apply only to
    // clusters with fNLocal==2

    Float_t xc[fNLocal], yc[fNLocal], zc[fNLocal];

    if (fNLocal > 1 && fNSaddle) {
      for (int i=0; i<fMul; i++) {
	  printf("before i %d\n",i);
          if (source[fIx[i]-amin][fIy[i]-tmin]  > thresh+3) lBin[fIx[i]*MAXNipy+fIy[i]].dig=fDig[i];
          else lBin[fIx[i]*MAXNipy+fIy[i]].dig=0;
	  printf("fDig %p \n",fDig[i]);

	  printf("i j source %d %d %f lBin %p \n",fIx[i],fIy[i],source[fIx[i]-amin][fIy[i]-tmin],lBin[fIx[i]*MAXNipy+fIy[i]].dig);
      }

      for (int i=0; i<fNLocal; i++) {
        int ipxl=fIx[fIndLocal[i]]; int ipyl=fIy[fIndLocal[i]];
	printf("ipxl ipyl %d %d\n",ipxl,ipyl);
	//printf("q %p %d \n",fDig[fIndLocal[i]],fDig[fIndLocal[i]]->fSignal);
	printf("%p  \n",lBin[ipxl*MAXNipy+ipyl].dig);

	PreCluster cnew; cnew.summit=lBin[ipxl*MAXNipy+ipyl].dig; cnew.idx=ncls;
	printf("q %p %d \n",fDig[fIndLocal[i]],lBin[ipx*MAXNipy+ipy].dig->fSignal);
	//FindCluster(ipx, ipy, lBin, cnew, 2*thresh+1); // it's further away
	FindCluster(ipxl, ipyl, lBin, cnew, 0);
	printf("finish FindCluster i %d\n",i);
	cnew.fX /= cnew.fQ;
	cnew.fZ /= cnew.fQ;
	ncls++;
	printf("i cnew.fX cnew.fY cnew.fZ %d %f %f %f\n",i,cnew.fX, cnew.fY, cnew.fZ);

	xc[i]=cnew.fX;
	yc[i]=cnew.fY;
	zc[i]=cnew.fZ;

	
        Int_t fMulnew=cnew.npoly;
	printf("fMulnew %d \n",fMulnew);

	for (int j=0; j<fMulnew; j++)
	  {
	    int idx=cnew.fDigitsList[j];
	    //printf("fMul idx %d  %d \n",i,idx);
	    digt = (AliITSdigitSDD*)ITSdigits->UncheckedAt(idx);
	    fDig[j]=digt;
	    //printf("digt %p \n",digt);
	    fIx[j]= digt->fCellX;
	    fIy[j]= digt->fCellY;
	    fSortIx[j]= digt->fCellX;
	    fSortIy[j]= digt->fCellY;
	    fQini[j] = digt->fSignal;
	    // printf("i fIx fIy fQini %d %d %d %d\n",j,fIx[j],fIy[j],fQini[j]);
	  }

	Int_t idx_sort[fMulnew];

	TMath::Sort(fMulnew,fSortIy,idx_sort,0);
	Int_t tmin=fIy[idx_sort[0]];
	Int_t tmax=fIy[idx_sort[fMulnew-1]];
    
        if (ipy-tmin != tmax-ipy) {
             printf("non-symetric cluster in time direction\n");
	     // do something - take Local maxima coord
	} else {
            // take cog coord
	}
        printf("ipy-tmin tmax-ipy %d %d\n",ipy-tmin, tmax-ipy);

	TMath::Sort(fMulnew,fSortIx,idx_sort,0);
	Int_t amin=fIx[idx_sort[0]];
	Int_t amax=fIx[idx_sort[fMulnew-1]];

        if (ipxl-amin != amax-ipxl) {
             printf("non-symetric cluster in anode direction\n");
	     // do something - take Local maxima coord
	} else {
            // take cog coord
	}
        printf("ipxl-amin amax-ipxl %d %d\n",ipxl-amin, amax-ipxl);

	printf("inside loop over local maxima: i ipxl ipyl %d %d %d\n",i,ipxl,ipyl);
      } // loop over local maxima

    } // if fNLocal


    if (fNLocal == 1) {
        int ipx1=fIx[fIndLocal[0]]; int ipy1=fIy[fIndLocal[0]];
        if (ipy1-tmin != tmax-ipy1) {
             printf("non-symetric cluster in time direction\n");
	     // do something - take Local maxima coord or fit
	} else {
            // take cog coord or fit
	}
        printf("ipy1-tmin tmax-ipy1 %d %d\n",ipy1-tmin, tmax-ipy1);

        if (ipx1-amin != amax-ipx1) {
             printf("non-symetric cluster in anode direction\n");
	     // do something - take Local maxima coord or fit
	} else {
            // take cog coord or fit
	}
        printf("ipx1-amin amax-ipx1 %d %d\n",ipx1-amin, amax-ipx1);
    }


    printf("Here! fNLocal fNSaddle %d %d\n",fNLocal,fNSaddle);
   
   AliITSpoints *points = 0;
   for (int i=0; i<fNLocal; i++) {
	points = new AliITSpoints();
	points->SetMarkerColor(kYellow);
	points->SetMarkerStyle(5);
	points->SetMarkerSize(1.);
	points->SetPoint(0,xc[i],yc[i],zc[i]);
	points->SetParticle(-1);
	points->Draw();
   }


   //points = new AliITSpoints(fNLocal);
   for (int i=0; i<fNLocal; i++) {
       points = new AliITSpoints();
       int idx=fIndLocal[i];
       points->SetMarkerColor(kGreen);
       points->SetMarkerStyle(5);
       points->SetMarkerSize(1.5);
       points->SetPoint(0,c.xpoly[idx],c.ypoly[idx],c.zpoly[idx]);
       points->SetParticle(-1);
       points->Draw();
   }
   //points = new AliITSpoints(fNSaddle);
   for (int i=0; i<fNSaddle; i++) {
       points = new AliITSpoints();
       int idx=fIndSaddle[i];
       points->SetMarkerColor(kWhite);
       points->SetMarkerStyle(5);
       points->SetMarkerSize(1.5);
       points->SetPoint(0,c.xpoly[idx],c.ypoly[idx],c.zpoly[idx]);
       points->SetParticle(-1);
       points->Draw();
   }


   for (int k=0;k<c.npoly;k++) {
     c.xpoly[k]=c.ypoly[k]=c.zpoly[k]=0;
   }
   c.npoly=0;

   delete [] lBin;
}

//_____________________________________________________________________________
static Double_t sinoid(Double_t *x, Double_t *par)
{
    Double_t arg = -2*TMath::Pi()*x[0];
    Double_t fitval= par[0]*TMath::Sin(arg)+
	par[1]*TMath::Sin(2*arg)+
	par[2]*TMath::Sin(3*arg)+
	par[3]*TMath::Sin(4*arg)+
	par[4]*TMath::Sin(5*arg);
    return fitval;
 }



//_____________________________________________________________________________
void AliITSpoints::
Neighbours(Int_t iX, Int_t iY, Int_t* Nlist, Int_t Xlist[8], Int_t Ylist[8])
{

    *Nlist=8;
    Xlist[0]=Xlist[1]=iX;
    if(iX) Xlist[2]=iX-1;
    else Xlist[2]=iX;
    if (iX < MAXNipx) Xlist[3]=iX+1;
    else Xlist[3]=iX;
    if(iY) Ylist[0]=iY-1;
    else Ylist[0]=iY;
    if (iY < MAXNipy) Ylist[1]=iY+1;
    else Ylist[1]=iY;
    Ylist[2]=Ylist[3]=iY;

    // diagonal elements

    Xlist[4]=Xlist[5]=Xlist[2];
    Xlist[6]=Xlist[7]=Xlist[3];

    Ylist[4]=Ylist[7]=Ylist[1];
    Ylist[5]=Ylist[6]=Ylist[0];
    

}

