//$Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ASV 

#include "AliL3StandardIncludes.h"
#include <TFile.h>
#include <TH1.h>
#include <TParticle.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TNtupleD.h>

#include <AliRun.h>
#include <AliSimDigits.h>
#include <AliTPC.h>
#include <AliTPCcluster.h>
#include <AliTPCClustersArray.h>
#include <AliTPCClustersRow.h>
#include <AliTPCParam.h>
#include <AliComplexCluster.h>

#if GCCVERSION == 3
#include <fstream>
#include <iosfwd>
#else
#include <fstream.h>
#endif

#include "AliL3Logging.h"
#include "AliL3Transform.h"
#include "AliL3SpacePointData.h"
#include "AliL3Track.h"
#include "AliL3FileHandler.h"
#include "AliL3TrackArray.h"
#include "AliL3Evaluate.h"

#if GCCVERSION == 3
using namespace std;
#endif

//_____________________________________________________________
// AliL3Evaluate
//
// Evaluation class for tracking. Plots efficiencies etc..
//


ClassImp(AliL3Evaluate)

AliL3Evaluate::AliL3Evaluate()
{
  fTracks = 0;
  fNFastPoints = 0;
  fMcIndex = 0;
  fMcId = 0;
  fMinSlice=0;
  fMaxSlice=0;
  fGoodTracks=0;
}

AliL3Evaluate::AliL3Evaluate(Char_t *datapath,Int_t min_clusters,Int_t minhits,Double_t minpt,Int_t *slice)
{
  
  if(slice)
    {
      fMinSlice=slice[0];
      fMaxSlice=slice[1];
    }
  else
    {
      fMinSlice=0;
      fMaxSlice=35;
    }
  fMaxFalseClusters = 0.1;
  fGoodFound = 0;
  fGoodGen = 0;
  fMinPointsOnTrack = min_clusters;
  fMinHitsFromParticle = minhits;
  fMinGoodPt = minpt;
  Char_t fname[1024];
  AliL3FileHandler *clusterfile[36][6];
  for(Int_t s=fMinSlice; s<=fMaxSlice; s++)
    {
      for(Int_t p=0; p<AliL3Transform::GetNPatches(); p++)
	{
	  fClusters[s][p] = 0;
	  clusterfile[s][p] = new AliL3FileHandler();
	  sprintf(fname,"%s/points_%d_%d.raw",datapath,s,p);
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

  sprintf(fname,"%s/tracks.raw",datapath);
  AliL3FileHandler *tfile = new AliL3FileHandler();
  if(!tfile->SetBinaryInput(fname)){
    LOG(AliL3Log::kError,"AliL3Evaluation::Setup","File Open")
      <<"Inputfile "<<fname<<" does not exist"<<ENDLOG; 
    return;
  }
  fTracks = new AliL3TrackArray();
  tfile->Binary2TrackArray(fTracks);
  tfile->CloseBinaryInput();
  delete tfile;
}


AliL3Evaluate::~AliL3Evaluate()
{
  if(fGoodTracks) delete fGoodTracks;
  if(fDigitsTree) fDigitsTree->Delete();
  if(fTracks) delete fTracks;
  if(fPtRes) delete fPtRes;
  if(fNGoodTracksPt) delete fNGoodTracksPt;
  if(fNFoundTracksPt) delete fNFoundTracksPt;
  if(fNFakeTracksPt) delete fNFakeTracksPt;
  if(fTrackEffPt) delete fTrackEffPt;
  if(fFakeTrackEffPt) delete fFakeTrackEffPt;
  if(fNGoodTracksEta) delete fNGoodTracksEta;
  if(fNFoundTracksEta) delete fNFoundTracksEta;
  if(fNFakeTracksEta) delete fNFakeTracksEta;
  if(fTrackEffEta) delete fTrackEffEta;
  if(fFakeTrackEffEta) delete fFakeTrackEffEta;
  if(fMcIndex) delete [] fMcIndex;
  if(fMcId)    delete [] fMcId;
  if(fNtuppel) delete fNtuppel;
}

  
void AliL3Evaluate::AssignIDs()
{
  //Assign MC id to the tracks.
#ifndef do_mc
  cerr<<"AliL3Evaluate::AssignIDs() : You need to compile with the do_mc flag!"<<endl;
  return;
#endif
  
  fTracks->QSort();
  LOG(AliL3Log::kDebug,"AliL3Evaluate::AssignIDs","Track Loop")
    <<"Assigning MC id to the found tracks...."<<ENDLOG;
  for(Int_t i=0; i<fTracks->GetNTracks(); i++)
    {
      AliL3Track *track = (AliL3Track*)fTracks->GetCheckedTrack(i);
      if(!track) continue; 
      if(track->GetNumberOfPoints() < fMinPointsOnTrack) break;
      
      fGoodFound++;
      Int_t tID = GetMCTrackLabel(track);
      track->SetMCid(tID);
    }
  cout<<"Found "<<fGoodFound<<" good tracks "<<endl;
}


struct S {Int_t lab; Int_t max;};
Int_t AliL3Evaluate::GetMCTrackLabel(AliL3Track *track){ 
  //Returns the MCtrackID of the belonging clusters.
  //If MCLabel < 0, means that track is fake.
  //Fake track means that more than 10 percent of clusters are assigned incorrectly.
  
#ifdef do_mc
  Int_t num_of_clusters = track->GetNumberOfPoints();
  S *s=new S[num_of_clusters];
  Int_t i;
  for (i=0; i<num_of_clusters; i++) s[i].lab=s[i].max=0;
  UInt_t *hitnum = track->GetHitNumbers();  
  UInt_t id;
    
  Int_t lab=123456789;
  for (i=0; i<num_of_clusters; i++) 
    {
      //Tricks to get the clusters belonging to this track:
      id = hitnum[i];
      Int_t slice = (id>>25) & 0x7f;
      Int_t patch = (id>>22) & 0x7;
      UInt_t pos = id&0x3fffff;	      
      
      AliL3SpacePointData *points = fClusters[slice][patch];
      if(!points) continue; 
      if(pos>=fNcl[slice][patch]) 
	{
	  LOG(AliL3Log::kError,"AliL3Evaluate::GetMCTrackLabel","Clusterarray")
	    <<AliL3Log::kDec<<"ERROR"<<ENDLOG;
	  continue;
	}
      
      //Get the label of the cluster:
      
      lab=abs(points[pos].fTrackID[0]);
      if(lab < 0)
	cout<<"Track had negative id : "<<lab<<" padrow "<<(Int_t)points[pos].fPadRow<<" nhits "<<num_of_clusters<<" pt "<<track->GetPt()<<endl;
      
      Int_t j;
      for (j=0; j<num_of_clusters; j++)
        if (s[j].lab==lab || s[j].max==0) break;
      s[j].lab=lab;
      s[j].max++;
    }
  
  Int_t max=0;
  for (i=0; i<num_of_clusters; i++) 
    if (s[i].max>max) {max=s[i].max; lab=s[i].lab;}
  
  delete[] s;
  
  for (i=0; i<num_of_clusters; i++) 
    {
      id = hitnum[i];
      Int_t slice = (id>>25) & 0x7f;
      Int_t patch = (id>>22) & 0x7;
      UInt_t pos = id&0x3fffff;	      
      
      AliL3SpacePointData *points = fClusters[slice][patch];
      if(!points) continue; 
      if(pos>=fNcl[slice][patch]) 
	{
	  LOG(AliL3Log::kError,"AliL3Evaluate::GetMCTrackLabel","Clusterarray")
	    <<AliL3Log::kDec<<"ERROR"<<ENDLOG;
	  continue;
	}
      
      if (abs(points[pos].fTrackID[1]) == lab || 
	  abs(points[pos].fTrackID[2]) == lab ) max++;
    }
  
  if (1.-Float_t(max)/num_of_clusters > fMaxFalseClusters) 
    {
      return -lab;
    }
  
  return lab;
#else //If we are running with mc_ids or not
  return 0;
#endif

}

void AliL3Evaluate::GetFastClusterIDs(Char_t *path)
{
  //Get the MC id of space points in case of using the fast simulator. 
  char fname[256];
  sprintf(fname,"%s/point_mc.dat",path);
  FILE *infile = fopen(fname,"r");
  if(!infile) return;
  Int_t hitid,hitmc,i;
  
  for(i=0; ; i++)
    if(fscanf(infile,"%d %d",&hitid,&hitmc)==EOF) break;
  rewind(infile);
  fNFastPoints = i;
  fMcId = new Int_t[fNFastPoints];
  fMcIndex = new UInt_t[fNFastPoints];
  
  for(i=0; i<fNFastPoints; i++)
    {
      if(fscanf(infile,"%d %d",&hitid,&hitmc)==EOF) break;
      fMcId[i] = hitmc;
      fMcIndex[i] = hitid;
    }
  fclose(infile);
}

void AliL3Evaluate::CreateHistos(Int_t nbin,Float_t xlow,Float_t xup)
{
  //Create the histograms 
  
  LOG(AliL3Log::kInformational,"AliL3Evaluate::CreateHistos","Allocating")
    <<"Creating histograms..."<<ENDLOG;
  
  fNtuppel = new TNtuple("fNtuppel","Pt resolution","pt_gen:pt_found:nHits");
  fNtuppel->SetDirectory(0);
  fPtRes = new TH1F("fPtRes","Relative Pt resolution",30,-10.,10.); 
  fNGoodTracksPt = new TH1F("fNGoodTracksPt","Good tracks vs pt",nbin,xlow,xup);    
  fNFoundTracksPt = new TH1F("fNFoundTracksPt","Found tracks vs pt",nbin,xlow,xup);
  fNFakeTracksPt = new TH1F("fNFakeTracksPt","Fake tracks vs pt",nbin,xlow,xup);
  fTrackEffPt = new TH1F("fTrackEffPt","Tracking efficiency vs pt",nbin,xlow,xup);
  fFakeTrackEffPt = new TH1F("fFakeTrackEffPt","Efficiency for fake tracks vs pt",nbin,xlow,xup);
  
  fNGoodTracksEta = new TH1F("fNGoodTracksEta","Good tracks vs eta",20,-50,50);
  fNFoundTracksEta = new TH1F("fNFoundTracksEta","Found tracks vs eta",20,-50,50);
  fNFakeTracksEta = new TH1F("fNFakeTracksEta","Fake tracks vs eta",20,-50,50);
  fTrackEffEta = new TH1F("fTrackEffEta","Tracking efficienct vs eta",20,-50,50);
  fFakeTrackEffEta = new TH1F("fFakeTrackEffEta","Efficiency for fake tracks vs eta",20,-50,50);

}

void AliL3Evaluate::GetGoodParticles(Char_t *path,Bool_t sector)
{
  //Read the good particles from file. This file should already have been
  //generated by macro AliTPCComparison.C.
  
  Char_t filename[1024];
  if(!sector)
    sprintf(filename,"%s/good_tracks_tpc",path);
  else
    sprintf(filename,"%s/good_tracks_tpc_sector",path);//Sectorwise comparison.
  ifstream in(filename);
  if(!in)
    {
      cerr<<"AliL3Evaluate::GetGoodParticles : Problems opening file :"<<filename<<endl;
      return;
    }
  Int_t MaxTracks=20000;
  fGoodTracks = new GoodTrack[MaxTracks];
  
  while (in>>fGoodTracks[fGoodGen].label>>fGoodTracks[fGoodGen].code>>
	 fGoodTracks[fGoodGen].px>>fGoodTracks[fGoodGen].py>>fGoodTracks[fGoodGen].pz>>
	 fGoodTracks[fGoodGen].x>>fGoodTracks[fGoodGen].y >>fGoodTracks[fGoodGen].z>>fGoodTracks[fGoodGen].nhits>>fGoodTracks[fGoodGen].sector) 

    {
      fGoodGen++;
      if (fGoodGen==MaxTracks) 
	{
	  cerr<<"AliL3Evaluate::GetGoodParticles : Too many good tracks !\n";
	  break;
	}
    }  
}


void AliL3Evaluate::FillEffHistos()
{  
  if(!fGoodTracks)
    {
      cerr<<"AliL3Evaluate::FillEffHistos : No good tracks"<<endl;
      return;
    }
  cout<<"Comparing "<<fGoodGen<<" good tracks ..."<<endl;
  for(Int_t i=0; i<fGoodGen; i++)
    {
      //cout<<"Checking particle "<<i<<endl;
      if(fGoodTracks[i].nhits < fMinHitsFromParticle) continue;
      Float_t ptg = TMath::Sqrt(fGoodTracks[i].px*fGoodTracks[i].px + fGoodTracks[i].py*fGoodTracks[i].py);
      if(ptg < fMinGoodPt) continue;
      Float_t pzg=fGoodTracks[i].pz;
      Float_t dipangle=TMath::ATan2(pzg,ptg)*180./TMath::Pi();

      fNGoodTracksPt->Fill(ptg);
      fNGoodTracksEta->Fill(dipangle);
      Int_t found = 0;
      
      for(Int_t k=0; k<fTracks->GetNTracks(); k++)
	{
	  AliL3Track *track = fTracks->GetCheckedTrack(k);
	  if(!track) continue;
	  Int_t nHits = track->GetNumberOfPoints();
	  if(nHits < fMinPointsOnTrack) break;
	  Int_t tracklabel;
	  tracklabel = track->GetMCid();
	  
	  if(TMath::Abs(tracklabel) != fGoodTracks[i].label) continue;
	  found=1;
	  if(tracklabel == fGoodTracks[i].label) {fNFoundTracksPt->Fill(ptg); fNFoundTracksEta->Fill(dipangle);}
	  else {fNFakeTracksPt->Fill(ptg); fNFakeTracksEta->Fill(dipangle);}
	  Float_t pt=track->GetPt();
	  fPtRes->Fill((pt-ptg)/ptg*100.);
	  fNtuppel->Fill(ptg,pt,nHits);
	  break;
	  
	}
      if(!found)
	cout<<"Track "<<fGoodTracks[i].label<<" was not found"<<endl;
    }
}

void AliL3Evaluate::FillEffHistosNAIVE()
{  
  //Fill the efficiency histograms.
  
  cout<<endl<<"Note: Doing NAIVE evaluation "<<endl;
  for(Int_t i=0; i<fGoodGen; i++)
    {
      Double_t ptg=TMath::Sqrt(fGoodTracks[i].px*fGoodTracks[i].px + fGoodTracks[i].py*fGoodTracks[i].py);
      Double_t pzg=fGoodTracks[i].pz;
      Float_t dipangle=TMath::ATan2(pzg,ptg)*180./TMath::Pi();
      //printf("filling particle with pt %f and dipangle %f\n",ptg,dipangle);
      fNGoodTracksPt->Fill(ptg);
      fNGoodTracksEta->Fill(dipangle);
      
    }
  
  for(Int_t k=0; k<fTracks->GetNTracks(); k++)
    {
      AliL3Track *track = fTracks->GetCheckedTrack(k);
      if(!track) continue;
      Int_t nHits = track->GetNumberOfPoints();
      if(nHits < fMinPointsOnTrack) break;
      if(track->GetPt()<fMinGoodPt) continue;
      if(fabs(track->GetPseudoRapidity())>0.9) continue;

      fNFoundTracksPt->Fill(track->GetPt()); fNFoundTracksEta->Fill(track->GetPseudoRapidity());
      //Float_t pt=track->GetPt();
      //fPtRes->Fill((pt-ptg)/ptg*100.);
      //fNtuppel->Fill(ptg,pt,nHits);
            
    }
}

void AliL3Evaluate::CalcEffHistos(){  
  
  Stat_t ngood=fNGoodTracksPt->GetEntries();
  Stat_t nfound=fNFoundTracksPt->GetEntries();
  Stat_t nfake=fNFakeTracksPt->GetEntries();
  
  LOG(AliL3Log::kInformational,"AliL3Evaluate::FillEffHistos","Efficiency")
    <<AliL3Log::kDec<<"There was "<<ngood<<" generated good tracks"<<ENDLOG;
  LOG(AliL3Log::kInformational,"AliL3Evaluate::FillEffHistos","Efficiency")
    <<AliL3Log::kDec<<"Found "<<nfound<<" tracks"<<ENDLOG;
  LOG(AliL3Log::kInformational,"AliL3Evaluate::FillEffHistos","Efficiency")
    <<AliL3Log::kDec<<"Integral efficiency is about "<<nfound/ngood*100<<ENDLOG;
  LOG(AliL3Log::kInformational,"AliL3Evaluate::FillEffHistos","Efficiency")
    <<AliL3Log::kDec<<"Fake tracks relative is about "<<nfake/ngood*100<<ENDLOG;
  
  LOG(AliL3Log::kInformational,"AliL3Evaluate::FillEffHistos","Efficiency")
    <<AliL3Log::kDec<<"Naive efficiency "<<(Double_t)fGoodFound/(Double_t)fGoodGen<<ENDLOG;

  
  fNFoundTracksPt->Sumw2(); fNGoodTracksPt->Sumw2();
  fTrackEffPt->Divide(fNFoundTracksPt,fNGoodTracksPt,1,1.,"b");
  fFakeTrackEffPt->Divide(fNFakeTracksPt,fNGoodTracksPt,1,1.,"b");
  fTrackEffPt->SetMaximum(1.4);
  fTrackEffPt->SetXTitle("P_{T} [GeV]");
  fTrackEffPt->SetLineWidth(2);
  fFakeTrackEffPt->SetFillStyle(3013);
  fTrackEffPt->SetLineColor(4);
  fFakeTrackEffPt->SetFillColor(2);
  
  fNFoundTracksEta->Sumw2(); fNGoodTracksEta->Sumw2();
  fTrackEffEta->Divide(fNFoundTracksEta,fNGoodTracksEta,1,1.,"b");
  fFakeTrackEffEta->Divide(fNFakeTracksEta,fNGoodTracksEta,1,1.,"b");
  fTrackEffEta->SetMaximum(1.4);
  fTrackEffEta->SetXTitle("#lambda [degrees]");
  fTrackEffEta->SetLineWidth(2);
  fFakeTrackEffEta->SetFillStyle(3013);
  fTrackEffEta->SetLineColor(4);
  fFakeTrackEffEta->SetFillColor(2);
       
}

void AliL3Evaluate::Write2File(Char_t *outputfile)
{
  //Write histograms to file:
  
  TFile *of = TFile::Open(outputfile,"RECREATE");
  if(!of->IsOpen())
    {
      LOG(AliL3Log::kError,"AliL3Evaluate::Write2File","File Open")
	<<"Problems opening rootfile"<<ENDLOG;
      return;
    }
  
  of->cd();
  fNtuppel->Write();
  fPtRes->Write();
  fNGoodTracksPt->Write();
  fNFoundTracksPt->Write();
  fNFakeTracksPt->Write();
  fTrackEffPt->Write();
  fFakeTrackEffPt->Write();
  fNGoodTracksEta->Write();
  fNFoundTracksEta->Write();
  fNFakeTracksEta->Write();
  fTrackEffEta->Write();
  fFakeTrackEffEta->Write();
  
  of->Close();

}

TNtupleD *AliL3Evaluate::CalculateResiduals(Char_t *datapath)
{

  TNtupleD *ntuppel=new TNtupleD("ntuppel","Residuals","residual_trans:residual_long:zHit:pt:dipangle:beta:padrow:nHits");
  ntuppel->SetDirectory(0);
  
  for(int f=fMinSlice; f<=fMaxSlice; f++)
    {
      AliL3FileHandler *tfile = new AliL3FileHandler();
      char fname[256];
      sprintf(fname,"%s/tracks_tr_%d_0.raw",datapath,f);
      if(!tfile->SetBinaryInput(fname)){
	LOG(AliL3Log::kError,"AliL3Evaluation::Setup","File Open")
	  <<"Inputfile "<<fname<<" does not exist"<<ENDLOG; 
	return 0;
      }
      fTracks = new AliL3TrackArray();
      tfile->Binary2TrackArray(fTracks);
      tfile->CloseBinaryInput();
      delete tfile;
      printf("Looking in slice %d\n",f);
      for(Int_t i=0; i<fTracks->GetNTracks(); i++)
	{
	  
	  AliL3Track *track = (AliL3Track*)fTracks->GetCheckedTrack(i);
	  if(!track) continue;
	  if(track->GetNHits() < fMinPointsOnTrack) continue;
	  
	  track->CalculateHelix();
	  UInt_t *hitnum = track->GetHitNumbers();
	  UInt_t id;
	  
	  Float_t xyz[3];
	  Int_t padrow;
	  for(Int_t j=0; j<track->GetNumberOfPoints()-1; j++)
	    {
	      id = hitnum[j];
	      Int_t slice = (id>>25) & 0x7f;
	      Int_t patch = (id>>22) & 0x7;
	      UInt_t pos = id&0x3fffff;	      
	      
	      //if(slice!=1) continue;
	      
	      AliL3SpacePointData *points = fClusters[slice][patch];
	      
	      if(!points) 
		{
		  LOG(AliL3Log::kError,"AliL3Evaluate::CalculateResiduals","Clusterarray")
		    <<"No points at slice "<<slice<<" patch "<<patch<<" pos "<<pos<<ENDLOG;
		  continue;
		}
	      if(pos>=fNcl[slice][patch]) 
		{
		  LOG(AliL3Log::kError,"AliL3Evaluate::CalculateResiduals","Clusterarray")
		    <<AliL3Log::kDec<<"ERROR"<<ENDLOG;
		  continue;
		}
	      
	      xyz[0] = points[pos].fX;
	      xyz[1] = points[pos].fY;
	      xyz[2] = points[pos].fZ;
	      padrow = points[pos].fPadRow;
	      AliL3Transform::Global2Local(xyz,slice);
	      
	      Float_t xyz_cross[3];
	      track->GetCrossingPoint(padrow,xyz_cross);
 	      Double_t beta = track->GetCrossingAngle(padrow);
	      
	      Double_t yres = xyz_cross[1] - xyz[1];
	      Double_t zres = xyz_cross[2] - xyz[2];
	      
	      Double_t dipangle = atan(track->GetTgl());
	      ntuppel->Fill(yres,zres,xyz_cross[2],track->GetPt(),dipangle,beta,padrow,track->GetNumberOfPoints());
	      
	    }
	}
      if(fTracks)
	delete fTracks;
    }
  return ntuppel;
}

TNtuple *AliL3Evaluate::EvaluatePoints(Char_t *path)
{
  //Compare points to the exact crossing points of track and padrows.
  //The input file to this function, contains the exact clusters calculated
  //in AliTPC::Hits2ExactClusters.
  
  cout<<"Evaluating points"<<endl;
  TNtuple *ntuppel = new TNtuple("ntuppel","residuals","slice:padrow:resy:resz:zHit:pt");
  ntuppel->SetDirectory(0);
  
  Char_t filename[1024];
  sprintf(filename,"%s/alirunfile.root",path);
  TFile *exfile = TFile::Open(filename);
  if(!exfile)
    {
      cerr<<"Error opening rootfile "<<filename<<endl;
      return 0;
    }
  gAlice = (AliRun*)exfile->Get("gAlice");
  if (!gAlice) 
    {
      LOG(AliL3Log::kError,"AliL3Evaluate::InitMC","gAlice")
	<<"AliRun object non existing on file"<<ENDLOG;
      return false;
    }
  
  gAlice->GetEvent(0);
  AliTPCParam *param = (AliTPCParam*)exfile->Get(AliL3Transform::GetParamName());
  
  //Get the exact clusters from file:
  AliTPCClustersArray *arr = new AliTPCClustersArray;
  arr->Setup(param);
  arr->SetClusterType("AliComplexCluster");
  char treeName[500];
  sprintf(treeName,"TreeCExact_%s",param->GetTitle());
  Bool_t clusterok = arr->ConnectTree(treeName);//Segment Tree (for offline clusters)
  if(!clusterok) {printf("AliL3Evaluate::EvaluatePoints : Error in clusterloading\n"); return 0;}
  
  cout<<"Entering loop with "<<(Int_t)arr->GetTree()->GetEntries()<<endl;
  for(Int_t i=0; i<arr->GetTree()->GetEntries(); i++)
    {
      //Get the exact clusters for this row:
      Int_t cursec,currow;
      AliSegmentID *s = arr->LoadEntry(i);
      param->AdjustSectorRow(s->GetID(),cursec,currow);
      
      AliTPCClustersRow *ro = (AliTPCClustersRow *)arr->GetRow(cursec,currow);
      TClonesArray *clusters = ro->GetArray();
      int num_of_offline=clusters->GetEntriesFast();
      
      //Get the found clusters:
      Int_t slice,padrow;
      AliL3Transform::Sector2Slice(slice,padrow,cursec,currow);
      if(slice<fMinSlice || slice>fMaxSlice) continue;
      AliL3SpacePointData *points = fClusters[slice][0];
      if(!points)
	{
	  cerr<<"AliL3Evaluate::EvalutePoints : Error getting clusters "<<endl;
	  return 0;
	}
      printf("Checking slice %d padrow %d with %d clusters\n",slice,padrow,num_of_offline);
      cout<<"There are "<<fNcl[slice][0]<<" clusters here"<<endl;
      for(UInt_t c=0; c<fNcl[slice][0]; c++)
	{
	  if((Int_t)points[c].fPadRow!=padrow) continue;
	  Float_t xyz_cl[3] = {points[c].fX,points[c].fY,points[c].fZ};
	  Float_t xyz_ex[3];
	  AliL3Transform::Global2Local(xyz_cl,cursec);

	  for(Int_t m=0; m<num_of_offline; m++)
	    {
	      AliComplexCluster *cluster = (AliComplexCluster *)clusters->UncheckedAt(m);
	      Int_t mcId = cluster->fTracks[0];

	      if(mcId <0) continue;
	      TParticle *part = gAlice->Particle(mcId);
	      if(points[c].fTrackID[0]!=mcId &&
		 points[c].fTrackID[1]!=mcId &&
		 points[c].fTrackID[2]!=mcId)
		continue;

	      AliL3Transform::Raw2Local(xyz_ex,cursec,currow,cluster->fY,cluster->fX);
	      
	      //In function AliTPC::Hits2ExactClusters the time offset is not included,
	      //so we have to substract it again here.
	      xyz_ex[2]-=AliL3Transform::GetZOffset();
	      
	      Float_t resy = xyz_cl[1] - xyz_ex[1];//cluster->GetY()
	      Float_t resz = xyz_cl[2] - xyz_ex[2];//cluster->GetZ()

	      ntuppel->Fill(slice,padrow,resy,resz,xyz_ex[2],part->Pt());
	    }
	}      
      arr->ClearRow(cursec,currow);
    }

  return ntuppel;
}

void AliL3Evaluate::GetCFeff(Char_t *outfile)
{
  
  TNtuple *ntuppel = new TNtuple("ntuppel","Cluster finder efficiency","row:ncrossings:nclusters");
  
  AliTPC *TPC = (AliTPC*)gAlice->GetDetector("TPC");
  
  TPC->SetParam(fParam);
  
  Int_t ver = TPC->IsVersion();
  LOG(AliL3Log::kInformational,"AliL3Evaluate::GetCFeff","TPC version")
    <<"TPC version "<<ver<<" found on file"<<ENDLOG;
  
  Int_t zero=TPC->GetParam()->GetZeroSup();
  
  Int_t np = gAlice->GetNtrack();
  
  Int_t crossed,recs;
  Int_t *count = new Int_t[np]; //np number of particles.
  Int_t i;
  Float_t xyz[3];
  for (i=0; i<np; i++) count[i]=0;
  for(Int_t sl=fMinSlice; sl<=fMaxSlice; sl++)
    {
      for (i=0; i<=175; i++) 
	{
	  crossed=0;
	  recs=0;
	  Int_t index = fRowid[sl][i];
	  if (!fDigitsTree->GetEvent(index)) continue;
	  Int_t sec,row;
	  fParam->AdjustSectorRow(fDigits->GetID(),sec,row);
	  fDigits->First();
	  do {
	    Int_t it=fDigits->CurrentRow(), ip=fDigits->CurrentColumn();
	    Short_t dig = fDigits->GetDigit(it,ip);
	    
	    if(dig<=fParam->GetZeroSup()) continue;
	    if(it < fParam->GetMaxTBin()-1 && it > 0)
	      if(fDigits->GetDigit(it+1,ip) <= fParam->GetZeroSup()
		 && fDigits->GetDigit(it-1,ip) <= fParam->GetZeroSup())
		continue;
	    
	    AliL3Transform::Raw2Local(xyz,sec,row,ip,it);
	    if(fParam->GetPadRowRadii(sec,row)<230./250.*fabs(xyz[2]))
	      continue;
	    
	    
	    Int_t idx0=fDigits->GetTrackID(it,ip,0); 
	    Int_t idx1=fDigits->GetTrackID(it,ip,1);
	    Int_t idx2=fDigits->GetTrackID(it,ip,2);
	    
	    if (idx0>=0 && dig>=zero) count[idx0]+=1;
	    if (idx1>=0 && dig>=zero) count[idx1]+=1;
	    if (idx2>=0 && dig>=zero) count[idx2]+=1;
	  } while (fDigits->Next());
	  
	  for (Int_t j=0; j<np; j++) 
	    {
	      if (count[j]>1) //at least two digits at this padrow 
		{
		  crossed++;
		  count[j]=0;
		}
	    }
	  AliL3SpacePointData *points = fClusters[sl][0];
	  for(UInt_t k=0; k<fNcl[sl][0]; k++)
	    {
	      if(points[k].fPadRow!=i) continue;
	      recs++;
	    }
	  ntuppel->Fill(i,crossed,recs);
	}
      
    }
  delete[] count;
  
  TFile *file = TFile::Open(outfile,"RECREATE");
  ntuppel->Write();
  file->Close();
  
}

