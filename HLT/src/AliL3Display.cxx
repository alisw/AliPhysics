//Author:        Anders Strand Vestbo
//Last Modified: 12.01.2001

#include <TCanvas.h>
#include <TView.h>
#include <TPolyMarker3D.h>
#include <TPolyLine3D.h>
#include <TH2.h>
#include <TTree.h>
#include <TNode.h>
#include <TGeometry.h>
#include <TShape.h>
#include <TParticle.h>

#include "AliL3Evaluate.h"
#include "AliRun.h"
#include "AliSimDigits.h"
#include "AliTPCParam.h"
#include "AliL3Transform.h"
#include "AliL3Track.h"
#include "AliL3TrackArray.h"
#include "AliL3SpacePointData.h"
#include "AliL3FileHandler.h"
#include "AliL3Logging.h"
#include "AliL3Display.h"

// AliL3Display
// Simple display class for the Level3 tracker.

ClassImp(AliL3Display)

AliL3Display::AliL3Display()
{
  fGeom = NULL;
  fTracks = NULL;
}

AliL3Display::AliL3Display(Int_t *slice)
{
  //Ctor. Specify which slices you want to look at.

  TFile *file = new TFile("/prog/alice/data/GEO/alice.geom");
  if(!file) printf("NO FILE\n");
  if(file->IsOpen())
    LOG(AliL3Log::kError,"AliL3Display::AliL3Display","File Open")
      <<"Geometry file alice.geom does not exist"<<ENDLOG;
  
  fGeom = (TGeometry*)file->Get("AliceGeom");
  fMinSlice = slice[0];
  fMaxSlice = slice[1];

  fTransform = new AliL3Transform();
  file->Close();
  delete file;
}

AliL3Display::~AliL3Display()
{

  if(fTracks)
    delete fTracks;
  if(fTransform)
    delete fTransform;
}

void AliL3Display::Setup(Char_t *trackfile,Char_t *path)
{
  //Read in the hit and track information from produced files.
  
  Char_t fname[256];
  AliL3FileHandler *clusterfile[36][5];
  for(Int_t s=fMinSlice; s<=fMaxSlice; s++)
    {
      for(Int_t p=0; p<5; p++)
	{
	  clusterfile[s][p] = new AliL3FileHandler();
	  sprintf(fname,"%spoints_%d_%d.raw",path,s,p);
	  if(!clusterfile[s][p]->SetBinaryInput(fname))
	    {
	      LOG(AliL3Log::kError,"AliL3Evaluation::Setup","File Open")
		<<"Inputfile "<<fname<<" does not exist"<<ENDLOG; 
	      delete clusterfile[s][p];
              clusterfile[s][p] = 0; 
	      continue;
	    }
	  fClusters[s][p] = (AliL3SpacePointData*)clusterfile[s][p]->Allocate();
	  clusterfile[s][p]->Binary2Memory(fNcl[s][p],fClusters[s][p]);
	  clusterfile[s][p]->CloseBinaryInput();
	}
    }
  
  
  AliL3FileHandler *tfile = new AliL3FileHandler();
  if(!tfile->SetBinaryInput(trackfile))
    {
      LOG(AliL3Log::kError,"AliL3Evaluation::Setup","File Open")
	<<"Inputfile "<<trackfile<<" does not exist"<<ENDLOG; 
      return;
    }
  fTracks = new AliL3TrackArray();
  tfile->Binary2TrackArray(fTracks);
  tfile->CloseBinaryInput();
  delete tfile;

}

void AliL3Display::DisplayTracks(Int_t min_hits)
{
  //Display the found tracks.

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
  Float_t xcl[176];
  Float_t ycl[176];
  Float_t zcl[176];
  
  for(Int_t j=0; j<ntracks; j++)
    {
      AliL3Track *gtrack = fTracks->GetCheckedTrack(j); 
      if(!gtrack) continue;        
      Int_t nHits = gtrack->GetNHits();
      UInt_t *hitnum = gtrack->GetHitNumbers();
      if(nHits < min_hits) continue;
      TPolyMarker3D *pm = new TPolyMarker3D(nHits);
      
      for(Int_t h=0; h<nHits; h++)
	{
	  UInt_t id=hitnum[h];
	  Int_t slice = (id>>25) & 0x7f;
	  Int_t patch = (id>>22) & 0x7;
	  UInt_t pos = id&0x3fffff;	      
	  AliL3SpacePointData *points = fClusters[slice][patch];
	  if(slice < fMinSlice || slice > fMaxSlice)
	    continue;

	  if(!points) {
	    LOG(AliL3Log::kError,"AliL3Display::DisplayTracks","Clusterarray")
	      <<"No points at slice "<<slice<<" patch "<<patch<<" pos "<<pos<<ENDLOG;
	    continue;
	  }
	  if(pos>=fNcl[slice][patch]) {printf("Error \n"); continue;}
	  Float_t xyz_tmp[3];
	  xyz_tmp[0] = points[pos].fX;
	  xyz_tmp[1] = points[pos].fY;
	  xyz_tmp[2] = points[pos].fZ;
	  	  
	  xcl[h] = xyz_tmp[0];
	  ycl[h] = xyz_tmp[1];
	  zcl[h] = xyz_tmp[2];
	  
	  pm->SetPoint(h,xcl[h],ycl[h],zcl[h]);
	}
      pm->SetMarkerColor(2);
      pm->Draw();
      TPolyLine3D *current_line = &(line[j]);
      current_line = new TPolyLine3D(nHits,xcl,ycl,zcl,"");
      
      current_line->SetLineColor(1);
      current_line->Draw("same");
            
    }
  
  /*Take this if you want black&white display for printing.
    Char_t fname[256];
    Int_t i;
    Int_t color = 1;
    for(i=0; i<10; i++)
    {
    sprintf(fname,"LS0%d",i);
    geom->GetNode(fname)->SetLineColor(color);
    sprintf(fname,"US0%d",i);
    geom->GetNode(fname)->SetLineColor(color);
    }
    for(i=10; i<18; i++)
    {
    sprintf(fname,"LS%d",i);
    geom->GetNode(fname)->SetLineColor(color);
    sprintf(fname,"US%d",i);
    geom->GetNode(fname)->SetLineColor(color);
    }
  */
  fGeom->Draw("same");
  
  c1->x3d();
  
}

void AliL3Display::DisplayClusters()
{
  //Display all clusters.
  
  TCanvas *c1 = new TCanvas("c1","",700,700);
  c1->cd();
  TView *v = new TView(1);
  v->SetRange(-430,-560,-430,430,560,1710);
  c1->Clear();
  c1->SetFillColor(1);
  c1->SetTheta(90.);
  c1->SetPhi(0.);
  
  for(Int_t s=fMinSlice; s<=fMaxSlice; s++)
    {
      for(Int_t p=0;p<5;p++)
	{
	  AliL3SpacePointData *points = fClusters[s][p];
	  if(!points) continue;
	  Int_t npoints = fNcl[s][p];
	  TPolyMarker3D *pm = new TPolyMarker3D(npoints);
	  
	  Float_t xyz[3];
	  for(Int_t i=0; i<npoints; i++)
	    {
	      
	      xyz[0] = points[i].fX;
	      xyz[1] = points[i].fY;
	      xyz[2] = points[i].fZ;
	      
	      pm->SetPoint(i,xyz[0],xyz[1],xyz[2]); 
	      
	    }
	  pm->SetMarkerColor(2);
	  pm->Draw("");
	}
    }
  fGeom->Draw("same");
  
  c1->x3d(); 
}


void AliL3Display::DisplayAll(Int_t min_hits)
{
  //Display tracks & all hits.

  
  TCanvas *c1 = new TCanvas("c1","",700,700);
  c1->cd();
  TView *v = new TView(1);
  v->SetRange(-430,-560,-430,430,560,1710);
  c1->Clear();
  c1->SetFillColor(1);
  c1->SetTheta(90.);
  c1->SetPhi(0.);
  
  for(Int_t s=fMinSlice; s<=fMaxSlice; s++)
    {
      for(Int_t p=0;p<1;p++)
	{
	  AliL3SpacePointData *points = fClusters[s][p];
	  if(!points) continue;
	  Int_t npoints = fNcl[s][p];
	  TPolyMarker3D *pm = new TPolyMarker3D(npoints);
	  
	  Float_t xyz[3];
	  for(Int_t i=0; i<npoints; i++){
	    xyz[0] = points[i].fX;
	    xyz[1] = points[i].fY;
	    xyz[2] = points[i].fZ;
	    
	    pm->SetPoint(i,xyz[0],xyz[1],xyz[2]); 
	    
	  }
	  pm->SetMarkerColor(2);
	  pm->Draw("");
	}
    }
  
  Int_t ntracks = fTracks->GetNTracks();
  TPolyLine3D *line = new TPolyLine3D[ntracks];
  Float_t xcl[176];
  Float_t ycl[176];
  Float_t zcl[176];
  
  
  for(Int_t j=0; j<ntracks; j++)
    {
      AliL3Track *gtrack = fTracks->GetCheckedTrack(j); 
      if(!gtrack) continue;        
      Int_t nHits = gtrack->GetNHits();
      UInt_t *hitnum = gtrack->GetHitNumbers();
      if(nHits < min_hits) continue;
      TPolyMarker3D *pm = new TPolyMarker3D(nHits);
      for(Int_t h=0; h<nHits; h++)
	{
	  UInt_t id=hitnum[h];
	  Int_t slice = (id>>25) & 0x7f;
	  Int_t patch = (id>>22) & 0x7;
	  UInt_t pos = id&0x3fffff;	      
	  
	  AliL3SpacePointData *points = fClusters[slice][patch];

	  if(!points) {
	    LOG(AliL3Log::kError,"AliL3Display::DisplayTracks","Clusterarray")
	      <<"No points at slice "<<slice<<" patch "<<patch<<" pos "<<pos<<ENDLOG;
	    continue;
	  }
	  if(pos>=fNcl[slice][patch]) {printf("Error \n"); continue;}
	  xcl[h] = points[pos].fX;
	  ycl[h] = points[pos].fY;
	  zcl[h] = points[pos].fZ;
	  pm->SetPoint(h,xcl[h],ycl[h],zcl[h]);
	}
      pm->SetMarkerColor(3);
      pm->Draw();
      TPolyLine3D *current_line = &(line[j]);
      current_line = new TPolyLine3D(nHits,xcl,ycl,zcl,"");
      
      current_line->SetLineColor(4);
      current_line->Draw("same");
    }

  fGeom->Draw("same");
  
  c1->x3d();
  
}

void AliL3Display::DisplayClusterRow(Int_t slice,Int_t padrow,Char_t *digitsFile,Char_t *type)
{
  //Display the found clusters on this row together with the raw data.
  
  
  TFile *file = new TFile(digitsFile);
  AliTPCParam *param = (AliTPCParam*)file->Get("75x40_100x60");
  TTree *TD=(TTree*)file->Get("TreeD_75x40_100x60");
  AliSimDigits da, *digits=&da;
  TD->GetBranch("Segment")->SetAddress(&digits); //Return pointer to branch segment.
  AliL3Transform *transform = new AliL3Transform();
  
  Int_t sector,row;
  transform->Slice2Sector(slice,padrow,sector,row);
  Int_t npads = param->GetNPads(sector,row);
  Int_t ntimes = param->GetMaxTBin();
  TH2F *histdig = new TH2F("histdig","",npads,0,npads-1,ntimes,0,ntimes-1);
  TH2F *histfast = new TH2F("histfast","",npads,0,npads-1,ntimes,0,ntimes-1);
  TH2F *histpart = new TH2F("histpart","",npads,0,npads-1,ntimes,0,ntimes-1);

  
  Int_t sectors_by_rows=(Int_t)TD->GetEntries();
  Int_t i;
  for (i=0; i<sectors_by_rows; i++) {
    if (!TD->GetEvent(i)) continue;
    Int_t sec,ro;
    param->AdjustSectorRow(digits->GetID(),sec,ro);
    
    if(sec != sector) continue;
    if(ro < row) continue;
    if(ro != row) break;
    printf("sector %d row %d\n",sec,ro);
    digits->First();
    while (digits->Next()) {
      Int_t it=digits->CurrentRow(), ip=digits->CurrentColumn();
      Short_t dig = digits->GetDigit(it,ip);
      if(dig<=param->GetZeroSup()) continue;
      if(it < param->GetMaxTBin()-1 && it > 0)
	if(digits->GetDigit(it+1,ip) <= param->GetZeroSup()
	   && digits->GetDigit(it-1,ip) <= param->GetZeroSup())
	  continue;
      
      histdig->Fill(ip,it,dig);
    }
  }
  
  /*file->cd();
  AliRun *gAlice = (AliRun*)file->Get("gAlice");
  gAlice->GetEvent(0);
  TClonesArray *fParticles=gAlice->Particles(); 
  TParticle *part = (TParticle*)fParticles->UncheckedAt(0);
  AliL3Evaluate *eval = new AliL3Evaluate();
  Float_t xyz_cross[3];
  */
  /*
  for(Int_t p=0;p<5;p++)
    {
      AliL3SpacePointData *points = fClusters[slice][p];
      if(!points) continue;

      Int_t npoints = fNcl[slice][p];     
      Float_t xyz[3];
      for(Int_t i=0; i<npoints; i++)
	{
	  if(points[i].fPadRow != padrow) continue;
	  xyz[0] = points[i].fX;
	  xyz[1] = points[i].fY;
	  xyz[2] = points[i].fZ;
	  transform->Local2Raw(xyz,sector,row);
	  histfast->Fill(xyz[1],xyz[2],1);
	  
	  
	  }
	  
	  }
  */
  TCanvas *c1 = new TCanvas("c1","",900,900);
  c1->cd();
  histdig->Draw();
  return;
  histfast->SetMarkerColor(2);
  histfast->SetMarkerStyle(4);
  histpart->SetMarkerColor(2);
  histpart->SetMarkerStyle(3);

  histdig->GetXaxis()->SetTitle("Pad #");
  histdig->GetYaxis()->SetTitle("Timebin #");
  histdig->Draw(type);
  histfast->Draw("psame");
  //histpart->Draw("psame");
}
