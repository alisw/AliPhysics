// @(#) $Id$

/** \class AliHLTTPCDisplay
<pre>
//_____________________________________________________________
// AliHLTTPCDisplay
//
// Simple display class for the HLT tracks.
</pre>
*/
// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group 

#include "AliHLTTPCStandardIncludes.h"
#include <TView.h>
#include <TPolyMarker3D.h>
#include <TPolyLine3D.h>
#include <TH2.h>
#include <TTree.h>
#include <TNode.h>
#include <TGeometry.h>
#include <TShape.h>
#include <TParticle.h>
#include <TFile.h>
#ifdef use_aliroot
#include <TClonesArray.h>
#include <AliRun.h>
#include <AliSimDigits.h>
#include <AliTPCParam.h>
#endif

#include "AliHLTTPCLogging.h"
#include "AliHLTTPCDisplay.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCTrack.h"
#include "AliHLTTPCTrackArray.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCMemHandler.h"

#if __GNUC__ == 3
using namespace std;
#endif


ClassImp(AliHLTTPCDisplay)

AliHLTTPCDisplay::AliHLTTPCDisplay()
{
  //constructor
  fGeom = NULL;
  fTracks = NULL;
  fc1 = new TCanvas("c1","",700,700);
  memset(fClusters,0,36*6*sizeof(AliHLTTPCSpacePointData*));
  memset(fNcl, 0, 36*6*sizeof(UInt_t));
}

AliHLTTPCDisplay::AliHLTTPCDisplay(Int_t *slice,Char_t *gfile)
{
  //ctor. Specify which slices you want to look at.
  LoadGeometrie(gfile);
  if (slice) {
    SetSlices(slice[0], slice[1]);
  }

  fc1 = new TCanvas("c1","",700,700);
  memset(fClusters,0,36*6*sizeof(AliHLTTPCSpacePointData*));
  memset(fNcl, 0, 36*6*sizeof(UInt_t));
}

AliHLTTPCDisplay::~AliHLTTPCDisplay()
{
  //destructor
  if(fTracks)
    delete fTracks;
  if (fc1)
    delete fc1;
}

Bool_t AliHLTTPCDisplay::SetSlices(Int_t minslice, Int_t maxslice) {
  fMinSlice = minslice;
  fMaxSlice = maxslice;
  return kTRUE;
}

Bool_t AliHLTTPCDisplay::LoadGeometrie(Char_t *gfile)
{
  if (gfile) {
  TFile *file = TFile::Open(gfile);
  if(!file)
    {
      LOG(AliHLTTPCLog::kError,"AliHLTTPCDisplay::AliHLTTPCDisplay","File Open")
	<<"Geometry file " << gfile << " does not exist!"<<ENDLOG;
      return kFALSE;
    }
  
  fGeom = (TGeometry*)file->Get("AliceGeom");

  file->Close();
  delete file;
  }
  return kTRUE;
}

void AliHLTTPCDisplay::SetupClusterDataForPatch(Int_t slice, Int_t patch, UInt_t nofClusters, AliHLTTPCSpacePointData* data)
{
  if (data && slice>=0 && slice<36 && patch>=0 && patch<AliHLTTPCTransform::GetNPatches()) {
    if (fClusters[slice][patch]!=NULL) {
      delete(fClusters[slice][patch]);
      fClusters[slice][patch]=NULL;
    }
    Int_t arraysize=nofClusters*sizeof(AliHLTTPCSpacePointData);
    fClusters[slice][patch] = (AliHLTTPCSpacePointData*)new Byte_t[arraysize];
    if (fClusters[slice][patch]) {
      memcpy(fClusters[slice][patch], data, arraysize);
      fNcl[slice][patch]=nofClusters;
    } else {
      fNcl[slice][patch]=nofClusters;
      LOG(AliHLTTPCLog::kError,"AliHLTTPCDisplay::SetupClusterDataForPatch","memory allocation")
	<<"mmemory allocation failed "<<ENDLOG; 
    }
  } else {
    LOG(AliHLTTPCLog::kError,"AliHLTTPCDisplay::SetupClusterDataForPatch","argument check")
      <<"invalid argument "<<ENDLOG; 
  } 
}

void AliHLTTPCDisplay::Setup(Char_t *trackfile,Char_t *path,Int_t event,Bool_t sp)
{
  //Read in the hit and track information from produced files.
  
  Char_t fname[256];
  AliHLTTPCMemHandler *clusterfile[36][6];
  memset(fClusters,0,36*6*sizeof(AliHLTTPCSpacePointData*));
  for(Int_t s=fMinSlice; s<=fMaxSlice; s++)
    {
      for(Int_t p=0; p<AliHLTTPCTransform::GetNPatches(); p++)
	{
	  Int_t patch;
	  if(sp==kTRUE)
	    patch=-1;
	  else
	    patch=p;
	  clusterfile[s][p] = new AliHLTTPCMemHandler();
	  if(event<0)
	    sprintf(fname,"%s/points_%d_%d.raw",path,s,patch);
	  else
	    sprintf(fname,"%s/points_%d_%d_%d.raw",path,event,s,patch);
	  if(!clusterfile[s][p]->SetBinaryInput(fname))
	    {
	      LOG(AliHLTTPCLog::kError,"AliHLTTPCEvaluation::Setup","File Open")
		<<"Inputfile "<<fname<<" does not exist"<<ENDLOG; 
	      delete clusterfile[s][p];
              clusterfile[s][p] = 0; 
	      continue;
	    }
	  fClusters[s][p] = (AliHLTTPCSpacePointData*)clusterfile[s][p]->Allocate();
	  clusterfile[s][p]->Binary2Memory(fNcl[s][p],fClusters[s][p]);
	  clusterfile[s][p]->CloseBinaryInput();
	  if(sp==kTRUE)
	    break;
	}
    }
  
  if(!trackfile) return;
  AliHLTTPCMemHandler *tfile = new AliHLTTPCMemHandler();
  if(!tfile->SetBinaryInput(trackfile))
    {
      LOG(AliHLTTPCLog::kError,"AliHLTTPCEvaluation::Setup","File Open")
	<<"Inputfile "<<trackfile<<" does not exist"<<ENDLOG; 
      return;
    }
  fTracks = new AliHLTTPCTrackArray();
  tfile->Binary2TrackArray(fTracks);
  tfile->CloseBinaryInput();
  delete tfile;

}

void AliHLTTPCDisplay::DisplayTracks(Int_t minhits,Bool_t x3don,Float_t thr)
{
  //Display the found tracks.

  if (!fc1) return;
  fc1->cd();
  
  TView *v = new TView(1);
  v->SetRange(-430,-560,-430,430,560,1710);
  fc1->Clear();
  fc1->SetFillColor(1);
  fc1->SetTheta(45.);
  fc1->SetPhi(0.);
    
  Int_t ntracks = fTracks->GetNTracks();
  TPolyLine3D *line = new TPolyLine3D[ntracks];
  Float_t xcl[176];
  Float_t ycl[176];
  Float_t zcl[176];
  
  for(Int_t j=0; j<ntracks; j++)
    {
      AliHLTTPCTrack *gtrack = fTracks->GetCheckedTrack(j); 
      if(!gtrack) continue;
      if((thr>=0)&&(gtrack->GetPt()<thr)) continue;        
      Int_t nHits = gtrack->GetNHits();
      UInt_t *hitnum = gtrack->GetHitNumbers();
      if(nHits < minhits) continue;
      TPolyMarker3D *pm = new TPolyMarker3D(nHits);
      Int_t hitcount=0;
      for(Int_t h=0; h<nHits; h++)
	{

	  UInt_t id=hitnum[h];
	  Int_t slice = (id>>25) & 0x7f;
	  Int_t patch = (id>>22) & 0x7;
	  UInt_t pos = id&0x3fffff;	      
	  //cout << h << " id " << pos << endl;
	  AliHLTTPCSpacePointData *points = fClusters[slice][patch];
	  if(slice < fMinSlice || slice > fMaxSlice)
	    continue;

	  if(!points) {
	    LOG(AliHLTTPCLog::kError,"AliHLTTPCDisplay::DisplayTracks","Clusterarray")
	      <<"No points at slice "<<slice<<" patch "<<patch<<" pos "<<pos<<ENDLOG;
	    continue;
	  }
	  if(pos>=fNcl[slice][patch]){
	    LOG(AliHLTTPCLog::kError,"AliHLTTPCDisplay::DisplayTracks","Clusterarray")
	      <<"Pos is too large: pos "<<pos <<" ncl "<<fNcl[slice][patch]<<ENDLOG;
	    continue;
	  }

	  Float_t xyztmp[3];
	  xyztmp[0] = points[pos].fX;
	  xyztmp[1] = points[pos].fY;
	  xyztmp[2] = points[pos].fZ;
	  	  
	  xcl[h] = xyztmp[0];
	  ycl[h] = xyztmp[1];
	  zcl[h] = xyztmp[2];
	  
	  pm->SetPoint(h,xcl[h],ycl[h],zcl[h]);
	  hitcount++;
	}
      if(hitcount==0) continue;
      pm->SetMarkerColor(2);
      pm->Draw();
      TPolyLine3D *currentline = &(line[j]);
      currentline = new TPolyLine3D(nHits,xcl,ycl,zcl,"");
      
      currentline->SetLineColor(4);
      currentline->Draw("same");
            
    }
  
  //Take this if you want black&white display for printing.
  Char_t fname[256];
  Int_t i;
  Int_t color = 1;
  fc1->SetFillColor(10);
  for(i=0; i<10; i++)
    {
      sprintf(fname,"LS0%d",i);
      fGeom->GetNode(fname)->SetLineColor(color);
      sprintf(fname,"US0%d",i);
      fGeom->GetNode(fname)->SetLineColor(color);
    }
  for(i=10; i<18; i++)
    {
      sprintf(fname,"LS%d",i);
      fGeom->GetNode(fname)->SetLineColor(color);
      sprintf(fname,"US%d",i);
      fGeom->GetNode(fname)->SetLineColor(color);
    }
  
  fGeom->Draw("same");
  
  if(x3don) fc1->x3d();
  
}

void AliHLTTPCDisplay::DisplayClusters(Bool_t x3don, Float_t* etaRange)
{
  //Display all clusters.
  
  if (!fc1) return;
  fc1->cd();

  TView *v = new TView(1);
  v->SetRange(-430,-560,-430,430,560,1710);
  fc1->Clear();
  fc1->SetFillColor(1);
  fc1->SetTheta(90.);
  fc1->SetPhi(0.);
  
  Int_t processed = 0, discarded = 0;
  for(Int_t s=fMinSlice; s<=fMaxSlice; s++)
    {
      for(Int_t p=0;p<6;p++)
	{
	  AliHLTTPCSpacePointData *points = fClusters[s][p];
	  if(!points) continue;
	  Int_t npoints = fNcl[s][p];
	  TPolyMarker3D *pm = new TPolyMarker3D(npoints);
	  
	  Float_t xyz[3];
	  for(Int_t i=0; i<npoints; i++)
	    {
	      xyz[0] = points[i].fX;
	      xyz[1] = points[i].fY;
	      xyz[2] = points[i].fZ;
	      if ( etaRange )
		  {
		  Double_t pointEta = AliHLTTPCTransform::GetEta( xyz );
		  if ( pointEta<etaRange[0] || pointEta>etaRange[1] )
		      {
		      discarded++;
		      continue;
		      }
		  }
	      processed++;
	      //AliHLTTPCTransform::Local2Global(xyz,s);
	      pm->SetPoint(i,xyz[0],xyz[1],xyz[2]); 
 	    }
	  pm->SetMarkerColor(2);
	  pm->Draw("");
	}
    }
  printf( "Processed: %d - Discarded: %d\n", processed, discarded );
  fGeom->Draw("same");
  fc1->Draw();
  
  if(x3don) fc1->x3d(); 
  fc1->Modified();
  fc1->Update();
}


void AliHLTTPCDisplay::DisplayAll(Int_t minhits,Bool_t x3don, Float_t* etaRange)
{
  //Display tracks & all hits.

  if (!fc1) return;
  fc1->cd();
  TView *v = new TView(1);
  v->SetRange(-430,-560,-430,430,560,1710);
  fc1->Clear();
  fc1->SetFillColor(1);
  fc1->SetTheta(90.);
  fc1->SetPhi(0.);
  
  Int_t processed = 0, discarded = 0;
  for(Int_t s=fMinSlice; s<=fMaxSlice; s++)
    {
      for(Int_t p=0;p<6;p++)
	{
	  AliHLTTPCSpacePointData *points = fClusters[s][p];
	  if(!points) continue;
	  Int_t npoints = fNcl[s][p];
	  TPolyMarker3D *pm = new TPolyMarker3D(npoints);
	  
	  Float_t xyz[3];
	  for(Int_t i=0; i<npoints; i++){
	    xyz[0] = points[i].fX;
	    xyz[1] = points[i].fY;
	    xyz[2] = points[i].fZ;
	    if ( etaRange )
		{
		Double_t pointEta = AliHLTTPCTransform::GetEta( xyz );
		if ( pointEta<etaRange[0] || pointEta>etaRange[1] )
		      {
		      discarded++;
		      continue;
		      }
		}
	    processed++;
	    
	    pm->SetPoint(i,xyz[0],xyz[1],xyz[2]); 
	    
	  }
	  pm->SetMarkerColor(2);
	  pm->Draw("");
	}
    }
  printf( "Processed: %d - Discarded: %d\n", processed, discarded );
  
  Int_t ntracks = fTracks->GetNTracks();
  TPolyLine3D *line = new TPolyLine3D[ntracks];
  Float_t xcl[176];
  Float_t ycl[176];
  Float_t zcl[176];
  
  for(Int_t j=0; j<ntracks; j++)
    {
      AliHLTTPCTrack *gtrack = fTracks->GetCheckedTrack(j); 
      if(!gtrack) continue;        
      Int_t nHits = gtrack->GetNHits();
      UInt_t *hitnum = gtrack->GetHitNumbers();
      if(nHits < minhits) continue;
      TPolyMarker3D *pm = new TPolyMarker3D(nHits);
      Int_t hitcount=0;
      for(Int_t h=0; h<nHits; h++)
	{
	  UInt_t id=hitnum[h];
	  Int_t slice = (id>>25) & 0x7f;
	  Int_t patch = (id>>22) & 0x7;
	  UInt_t pos = id&0x3fffff;	      
	  if(slice < fMinSlice || slice > fMaxSlice)
	    continue;
	  
	  AliHLTTPCSpacePointData *points = fClusters[slice][patch];
	  if(!points) {
	    LOG(AliHLTTPCLog::kError,"AliHLTTPCDisplay::DisplayAll","Clusterarray")
	      <<"No points at slice "<<slice<<" patch "<<patch<<" pos "<<pos<<ENDLOG;
	    continue;
	  }
	  if(pos>=fNcl[slice][patch]) {
	    LOG(AliHLTTPCLog::kError,"AliHLTTPCDisplay::DisplayAll","Clusterarray")
	      <<"Pos is too large: pos "<<pos <<" ncl "<<fNcl[slice][patch]<<ENDLOG;
	    continue;
	  }
	  xcl[h] = points[pos].fX;
	  ycl[h] = points[pos].fY;
	  zcl[h] = points[pos].fZ;
	  pm->SetPoint(h,xcl[h],ycl[h],zcl[h]);
	  hitcount++;
	}
      if(hitcount==0) continue;
      pm->SetMarkerColor(3);
      pm->Draw();
      TPolyLine3D *currentline = &(line[j]);
      currentline = new TPolyLine3D(nHits,xcl,ycl,zcl,"");
      currentline->SetLineColor(4);
      currentline->SetLineWidth(2);
      currentline->Draw("same");
    }
  
  Char_t fname[256];
  Int_t i;
  Int_t color = 1;
  fc1->SetFillColor(10);
  for(i=0; i<10; i++)
    {
      sprintf(fname,"LS0%d",i);
      fGeom->GetNode(fname)->SetLineColor(color);
      sprintf(fname,"US0%d",i);
      fGeom->GetNode(fname)->SetLineColor(color);
    }
  for(i=10; i<18; i++)
    {
      sprintf(fname,"LS%d",i);
      fGeom->GetNode(fname)->SetLineColor(color);
      sprintf(fname,"US%d",i);
      fGeom->GetNode(fname)->SetLineColor(color);
    }
    
  fGeom->Draw("same");
  
  if(x3don) fc1->x3d();
}

void AliHLTTPCDisplay::DisplayClusterRow(Int_t slice,Int_t padrow,Char_t *digitsFile,Char_t *type)
{
  //Display the found clusters on this row together with the raw data.
  
  if (!fc1) return;
#ifdef use_aliroot
  TFile *file = new TFile(digitsFile);
  AliTPCParam *param = (AliTPCParam*)file->Get(AliHLTTPCTransform::GetParamName());

  Char_t dname[100];
  sprintf(dname,"TreeD_%s_0",AliHLTTPCTransform::GetParamName());
  TTree *td=(TTree*)file->Get(dname);
  AliSimDigits da, *digits=&da;
  td->GetBranch("Segment")->SetAddress(&digits); //Return pointer to branch segment.
  
  Int_t sector,row;
  AliHLTTPCTransform::Slice2Sector(slice,padrow,sector,row);
  Int_t npads = param->GetNPads(sector,row);
  Int_t ntimes = param->GetMaxTBin();
  TH2F *histdig = new TH2F("histdig","",npads,0,npads-1,ntimes,0,ntimes-1);
  TH2F *histfast = new TH2F("histfast","",npads,0,npads-1,ntimes,0,ntimes-1);
  TH2F *histpart = new TH2F("histpart","",npads,0,npads-1,ntimes,0,ntimes-1);

  
  Int_t sectorsbyrows=(Int_t)td->GetEntries();
  Int_t i;
  for (i=0; i<sectorsbyrows; i++) {
    if (!td->GetEvent(i)) continue;
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
      /*
      if(it < param->GetMaxTBin()-1 && it > 0)
	if(digits->GetDigit(it+1,ip) <= param->GetZeroSup()
	   && digits->GetDigit(it-1,ip) <= param->GetZeroSup())
	  continue;
      */
      histdig->Fill(ip,it,dig);
    }
  }
  
  /*file->cd();
  AliRun *gAlice = (AliRun*)file->Get("gAlice");
  gAlice->GetEvent(0);
  TClonesArray *fParticles=gAlice->Particles(); 
  TParticle *part = (TParticle*)fParticles->UncheckedAt(0);
  AliHLTTPCEvaluate *eval = new AliHLTTPCEvaluate();
  Float_t xyzcross[3];
  */
  
  for(Int_t p=0;p<6;p++)
    {
      AliHLTTPCSpacePointData *points = fClusters[slice][p];
      if(!points) continue;
      
      Int_t npoints = fNcl[slice][p];     
      Float_t xyz[3];
      for(Int_t i=0; i<npoints; i++)
	{
	  if(points[i].fPadRow != padrow) continue;
	  xyz[0] = points[i].fX;
	  xyz[1] = points[i].fY;
	  xyz[2] = points[i].fZ;
	  AliHLTTPCTransform::Global2Raw(xyz,sector,row);
	  //AliHLTTPCTransform::Local2Raw(xyz,sector,row);
	  histfast->Fill(xyz[1],xyz[2],1);
	  
	  
	}
      
    }
  
  fc1->cd();
  histdig->Draw();
  histfast->SetMarkerColor(2);
  histfast->SetMarkerStyle(4);
  histpart->SetMarkerColor(2);
  histpart->SetMarkerStyle(3);

  histdig->GetXaxis()->SetTitle("Pad #");
  histdig->GetYaxis()->SetTitle("Timebin #");
  histdig->Draw(type);
  histfast->Draw("psame");
  //histpart->Draw("psame");

#endif
  return;
}
