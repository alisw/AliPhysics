//Author:        Anders Strand Vestbo
//Last Modified: 12.01.2001

#include <TCanvas.h>
#include <TView.h>
#include <TPolyMarker3D.h>
#include <TPolyLine3D.h>
#include <TNode.h>
#include <TGeometry.h>
#include <TShape.h>
#include <TFile.h>

#include "AliL3Transform.h"
#include "AliL3HoughTrack.h"
#include "AliL3TrackArray.h"
#include "AliL3Logging.h"
#include "AliL3HoughDisplay.h"
#include "AliL3Defs.h"

ClassImp(AliL3HoughDisplay)


AliL3HoughDisplay::AliL3HoughDisplay()
{
  //Ctor. Specify which slices you want to look at.

  fTracks = 0;

  Init();
}

AliL3HoughDisplay::~AliL3HoughDisplay()
{
  if(fTransform)
    delete fTransform;
}

void AliL3HoughDisplay::Init()
{
  TFile *file = TFile::Open("/prog/alice/data/GEO/alice.geom");
  if(!file->IsOpen())
    LOG(AliL3Log::kError,"AliL3HoughDisplay::AliL3HoughDisplay","File Open")
      <<"Geometry file alice.geom does not exist"<<ENDLOG;
  fGeom = (TGeometry*)file->Get("AliceGeom");
  file->Close();
  fTransform = new AliL3Transform();
}


void AliL3HoughDisplay::GenerateHits(AliL3HoughTrack *track,Float_t *x,Float_t *y,Float_t *z,Int_t &n)
{
  n=0;
  Float_t xyz[3];
  Int_t slice = track->GetSlice();
  for(Int_t i=track->GetFirstRow(); i<track->GetLastRow(); i++)
    {
      if(track->GetCrossingPoint(i,xyz))
	{
	  fTransform->Local2Global(xyz,slice);
	  x[n] = xyz[0];
	  y[n] = xyz[1];
	  z[n] = xyz[2];
	  n++;
	}
      else
	break;
    }
}

void AliL3HoughDisplay::DisplayTracks()
{
  //Display the found tracks.
  
  if(!fTracks)
    {
      printf("AliL3HoughDisplay::DisplayTracks() : No tracks\n");
      return;
    }
  
  TCanvas *c1 = new TCanvas("c1","",700,700);
  c1->cd();
  
  TView *v = new TView(1);
  v->SetRange(-430,-560,-430,430,560,1710);
  
  c1->Clear();
  c1->SetFillColor(1);
  c1->SetTheta(90.);
  c1->SetPhi(0.);
    
  Int_t ntracks = fTracks->GetNTracks();
  TPolyLine3D *line = new TPolyLine3D[ntracks];
  
  Int_t n;
  Float_t x[176],y[176],z[176];
  for(Int_t j=0; j<ntracks; j++)
    {
      AliL3HoughTrack *track = (AliL3HoughTrack*)fTracks->GetCheckedTrack(j); 
      if(!track) continue;        
      GenerateHits(track,x,y,z,n);
      TPolyMarker3D *pm = new TPolyMarker3D(n);
      for(Int_t h=0; h<n; h++)
	pm->SetPoint(h,x[h],y[h],z[h]);
      
      pm->SetMarkerColor(2);
      pm->Draw();
      TPolyLine3D *current_line = &(line[j]);
      current_line = new TPolyLine3D(n,x,y,z,"");
      current_line->SetLineColor(1);
      current_line->Draw("same");
      
    }
  
  fGeom->Draw("same");
  c1->x3d();
}
