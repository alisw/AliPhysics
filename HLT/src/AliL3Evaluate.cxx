//Author:        Anders Strand Vestbo
//Last Modified: 5.01.2001

#include <stdio.h>
#include <TFile.h>
#include <TH1.h>
#include <TParticle.h>
#include <TTree.h>
#include <TClonesArray.h>

#include "AliRun.h"
#include "AliSimDigits.h"
#include "AliTPC.h"
#include "AliTPCClustersArray.h"
#include "AliTPCClustersRow.h"
#include "AliTPCcluster.h"
#include "AliTPCParam.h"

#include "AliL3Transform.h"
#include "AliL3SpacePointData.h"
#include "AliL3Track.h"
#include "AliL3FileHandler.h"
#include "AliL3TrackArray.h"
#include "AliL3Evaluate.h"
#include "AliL3Logging.h"

//AliL3Evaluate
//Class for tracking evaluation.

ClassImp(AliL3Evaluate)

AliL3Evaluate::AliL3Evaluate()
{
  fMCFile = NULL;
  fTracks = NULL;
  fMCclusterfile = NULL;
}

AliL3Evaluate::AliL3Evaluate(Char_t *mcfile,Int_t *slice)
{
  //Normal constructor. Input are the rootfile containing the 
  //original MC information of the simulated event. 

  fMCFile = new TFile(mcfile,"READ");
  if(!fMCFile->IsOpen())
    {
      LOG(AliL3Log::kError,"AliL3Evaluation::AliL3Evaluation","File Open")
	<<"Inputfile "<<mcfile<<" does not exist"<<ENDLOG;
      return;
    }
  
  fParam = (AliTPCParam*)fMCFile->Get("75x40_100x60");
  fTransform = new AliL3Transform();
  
  fMinSlice = slice[0];
  fMaxSlice = slice[1];

}

AliL3Evaluate::~AliL3Evaluate()
{
  if(fDigitsTree) 
    fDigitsTree->Delete();
  

  /*if(fMCFile) 
    {
      fMCFile->Close();
      delete fMCFile;
    }
  */
  if(fTransform)
    delete fTransform;
  if(fTracks)
    delete fTracks;
  if(fPtRes)
    delete fPtRes;
  if(fNGoodTracksPt)
    delete fNGoodTracksPt;
  if(fNFoundTracksPt)
    delete fNFoundTracksPt;
  if(fNFakeTracksPt)
    delete fNFakeTracksPt;
  if(fTrackEffPt)
    delete fTrackEffPt;
  if(fFakeTrackEffPt)
    delete fFakeTrackEffPt;
  if(fNGoodTracksEta)
    delete fNGoodTracksEta;
  if(fNFoundTracksEta)
    delete fNFoundTracksEta;
  if(fNFakeTracksEta)
    delete fNFakeTracksEta;
  if(fTrackEffEta)
    delete fTrackEffEta;
  if(fFakeTrackEffEta)
    delete fFakeTrackEffEta;
}

void AliL3Evaluate::Setup(Char_t *trackfile,Char_t *clustfile)
{
  //Setup for using the slow simulator.
  //Read in the hit and track information from produced files.
  
  fIsSlow = true;

  Char_t fname[256];
  AliL3FileHandler *clusterfile[36][5];
  for(Int_t s=fMinSlice; s<=fMaxSlice; s++)
    {
      for(Int_t p=0; p<5; p++)
	{
	  clusterfile[s][p] = new AliL3FileHandler();
	  sprintf(fname,"points_%d_%d.raw",s,p);
	  if(!clusterfile[s][p]->SetBinaryInput(fname))
	    {
	      LOG(AliL3Log::kError,"AliL3Evaluation::Setup","File Open")
		<<"Inputfile "<<fname<<" does not exist"<<ENDLOG; 
	      return;
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
  
  if(!SetDigitsTree())
    LOG(AliL3Log::kError,"AliL3Evaluation::Setup","Digits Tree")
      <<"Error setting up digits tree"<<ENDLOG;
  if(!SetMCParticleArray())
    LOG(AliL3Log::kError,"AliL3Evaluation::Setup","Particle array")
      <<"Error setting up particle array"<<ENDLOG;
  
}

void AliL3Evaluate::SetupFast(Char_t *trackfile,Char_t *clustfile,Char_t *mcClusterfile)
{
  //Setup for using the fast simulator.

  fIsSlow = false;
  
  fMCclusterfile = new TFile(mcClusterfile);
  if(!fMCclusterfile->IsOpen())
    LOG(AliL3Log::kError,"AliL3Evaluation::SetupFast","File Open")
      <<"Inputfile "<<mcClusterfile<<" does not exist"<<ENDLOG; 

  Char_t fname[256];
  AliL3FileHandler *clusterfile[36][5];
  for(Int_t s=fMinSlice; s<=fMaxSlice; s++)
    {
      for(Int_t p=0; p<5; p++)
	{
	  clusterfile[s][p] = new AliL3FileHandler();
	  sprintf(fname,"points_%d_%d.raw",s,p);
	  if(!clusterfile[s][p]->SetBinaryInput(fname))
	    {
	      LOG(AliL3Log::kError,"AliL3Evaluation::Setup","File Open")
		<<"Inputfile "<<fname<<" does not exist"<<ENDLOG; 
	      return;
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
  
  if(!SetMCParticleArray())
    LOG(AliL3Log::kError,"AliL3Evaluation::SetupFast","Particle array")
      <<"Error setting up particle array"<<ENDLOG;
}

void AliL3Evaluate::CreateHistos(Int_t nbin,Int_t xlow,Int_t xup)
{
  //Create the histograms 
  
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

Bool_t AliL3Evaluate::SetDigitsTree()
{

  fMCFile->cd();
  fDigitsTree = (TTree*)fMCFile->Get("TreeD_75x40_100x60");
  if(!fDigitsTree) return false;
  fDigitsTree->GetBranch("Segment")->SetAddress(&fDigits);
  for(Int_t i=0; i<fDigitsTree->GetEntries(); i++)
    {
      if(!fDigitsTree->GetEvent(i)) continue;
      Int_t se,ro,slice,slicerow;
      fParam->AdjustSectorRow(fDigits->GetID(),se,ro);
      fTransform->Sector2Slice(slice,slicerow,se,ro);
      fRowid[slice][slicerow] = i;
    }
  
  return true;
}

Bool_t AliL3Evaluate::SetMCParticleArray()
{
  fMCFile->cd();
  AliRun *gAlice = (AliRun*)fMCFile->Get("gAlice");
  if (!gAlice) 
    {
      LOG(AliL3Log::kError,"AliL3Evaluate::SetParticleArray","gAlice")
	<<"AliRun object non existing on file"<<ENDLOG;
      return false;
    }
  
  gAlice->GetEvent(0);
  fParticles=gAlice->Particles(); 
  return true;
}


TObjArray *AliL3Evaluate::DefineGoodTracks(Int_t slice,Int_t *padrow,Int_t good_number,Int_t *particle_id)
{
  //Loop over MC particles, and mark the good ones
  //(which the tracker should find...)

  Int_t np=fParticles->GetEntriesFast();
  AliTPC *TPC = (AliTPC*)gAlice->GetDetector("TPC");
  TPC->SetParam(fParam);
  Int_t ver = TPC->IsVersion();
  LOG(AliL3Log::kInformational,"AliL3Evaluate::DefineGoodTracks","TPC version")
    <<"TPC version "<<ver<<" found on file"<<ENDLOG;
  
  //Int_t nrow_up=TPC->GetParam()->GetNRowUp();
  //Int_t nrows=TPC->GetParam()->GetNRowLow()+nrow_up;
  Int_t zero=TPC->GetParam()->GetZeroSup();
  
  //Int_t number_of_rows = padrow[1] - padrow[0] + 1;
  
  Int_t *good = new Int_t[np];
  for(Int_t ii=0; ii<np; ii++)
    good[ii] = 0;

  if(ver==1)
    {
      if(fIsSlow)
	LOG(AliL3Log::kError,"AliL3Evaluate::DefineGoodTracks","TPC version")
	  <<"TPC version "<<ver<<" does not match."<<ENDLOG;
      fMCclusterfile->cd();
      AliTPCClustersArray carray;
      carray.Setup(fParam);
      carray.SetClusterType("AliTPCcluster");
      Bool_t clusterok = carray.ConnectTree("Segment Tree");
      if(!clusterok) 
	LOG(AliL3Log::kError,"AliL3Evaluate::DefineGoodTracks","Cluster Array")
	  <<"Error loading clusters from rootfile"<<ENDLOG;
      
      for(Int_t i=padrow[0]; i<=padrow[1]; i++)
	{
	  Int_t sec,row,sl,lr;
	  AliSegmentID *s = carray.LoadEntry(i);
	  fParam->AdjustSectorRow(s->GetID(),sec,row);
	  fTransform->Sector2Slice(sl,lr,sec,row);
	  
	  if(sl != slice) {carray.ClearRow(sec,row); continue;}
	  if(lr != i) {carray.ClearRow(sec,row); continue;}
	  AliTPCClustersRow *cRow = carray.GetRow(sec,row);
	  for(Int_t j=0; j<cRow->GetArray()->GetEntriesFast(); j++)
	    {
	      AliTPCcluster *cluster=(AliTPCcluster*)(*cRow)[j];
	      Int_t lab=cluster->GetLabel(0);
	      if(lab<0) continue;
	      lab=TMath::Abs(lab);
	      good[lab]++;
	    }
	  if(carray.GetRow(sec,row)) 
	    carray.ClearRow(sec,row);
	}
    }
  else if(ver==2)
    {
      if(!fIsSlow)
	LOG(AliL3Log::kError,"AliL3Evaluate::DefineGoodTracks","TPC version")
	  <<"TPC version "<<ver<<" does not match."<<ENDLOG;
      Int_t *count = new Int_t[np]; //np number of particles.
      Int_t i;
      for (i=0; i<np; i++) count[i]=0;
      for (i=padrow[0]; i<=padrow[1]; i++) {
	Int_t index = fRowid[slice][i];
	if (!fDigitsTree->GetEvent(index)) continue;
	Int_t sec,row;
	fParam->AdjustSectorRow(fDigits->GetID(),sec,row);
	fDigits->First();
	while (fDigits->Next()) {
	  Int_t it=fDigits->CurrentRow(), ip=fDigits->CurrentColumn();
	  Short_t dig = fDigits->GetDigit(it,ip);
	  Int_t idx0=fDigits->GetTrackID(it,ip,0); 
	  Int_t idx1=fDigits->GetTrackID(it,ip,1);
	  Int_t idx2=fDigits->GetTrackID(it,ip,2);
	  if (idx0>=0 && dig>=zero) count[idx0]+=1;
	  if (idx1>=0 && dig>=zero) count[idx1]+=1;
	  if (idx2>=0 && dig>=zero) count[idx2]+=1;
	}
	for (Int_t j=0; j<np; j++) {
	  if (count[j]>1) {//at least two digits at this padrow 
	    good[j]++;
	  }
	  count[j]=0;
	}
      }
      delete[] count;
    }
  else 
    {
      LOG(AliL3Log::kError,"AliL3Evaluation::FillEffHistos","TPC version")
	<<"No valid TPC version found"<<ENDLOG;
      return 0;
    }
  
  Float_t torad=TMath::Pi()/180;
  
  Float_t phi_min = slice*20 - 10;
  Float_t phi_max = slice*20 + 10;
  TObjArray *good_part = new TObjArray();
  
  for(Int_t i=0; i<fParticles->GetEntriesFast(); i++)
   {
     TParticle *p = (TParticle*)fParticles->UncheckedAt(i);
     if(p->GetFirstMother()>0) continue; //secondary particle
     if(good[i] < good_number) {continue;}
     
     Double_t ptg=p->Pt(),pxg=p->Px(),pyg=p->Py(),pzg=p->Pz();
     Double_t phi_part = TMath::ATan2(pyg,pxg);
     if (phi_part < 0) phi_part += 2*TMath::Pi();
     

     if(phi_part < phi_min*torad || phi_part > phi_max*torad) {continue;}
     if(ptg<0.100) continue;
     if(fabs(pzg/ptg)>0.999) {continue;}
     Int_t entries = good_part->GetEntriesFast();
     good_part->AddLast(p);
     particle_id[entries] = i;
   }
  delete [] good;
  return good_part;
}

void AliL3Evaluate::EvaluatePatch(Int_t slice,Int_t patch,Int_t min_points,Int_t good_number)
{
  //Make efficiency plots for tracking on patch level (before any merging).
  
  Int_t row[5][2] = {{ 0, 45},{46,77},{78,109},{110,141},{142,173}};
  Int_t *particle_id = new Int_t[fParticles->GetEntriesFast()];
  TObjArray *good_particles = DefineGoodTracks(slice,row[patch],good_number,particle_id);
  SetMinPoints(min_points);
  AssignIDs();
  FillEffHistos(good_particles,particle_id);
  delete good_particles;
  delete [] particle_id;
}

void AliL3Evaluate::EvaluateSlice(Int_t slice,Int_t min_points,Int_t good_number)
{
  //Make efficiency plots for tracking on a slice (after merging).
  //min_points = minimum points on track to be considered for evaluation
  //good_number = minimum hits (padrows) produced by simulated track for consideration.

  Int_t row[174] = {0,173};
  Int_t *particle_id = new Int_t[fParticles->GetEntriesFast()];
  TObjArray *good_particles = DefineGoodTracks(slice,row,good_number,particle_id);
  SetMinPoints(min_points);
  AssignIDs();
  FillEffHistos(good_particles,particle_id);
  delete good_particles;
  delete [] particle_id;
}

void AliL3Evaluate::EvaluateGlobal()
{
  //Make efficiency plots for tracking on several slices.
  /*
  Int_t row[174] = {0,173};
  TObjArray *good_particles = DefineGoodTracks(slice,row,good_number);
  FillEffHistos(good_particles);
  delete good_particles;
  */
}

void AliL3Evaluate::FillEffHistos(TObjArray *good_particles,Int_t *particle_id)
{  
  //Fill the efficiency histograms.

  
  CreateHistos();

  for(Int_t i=0; i<good_particles->GetEntriesFast(); i++)
    {
      TParticle *p = (TParticle*)good_particles->UncheckedAt(i);
      Double_t ptg=p->Pt(),pzg=p->Pz();
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
	  
	  if(TMath::Abs(tracklabel) != particle_id[i]) continue;
	  found=1;
	  if(tracklabel == particle_id[i]) {fNFoundTracksPt->Fill(ptg); fNFoundTracksEta->Fill(dipangle);}
	  else {fNFakeTracksPt->Fill(ptg); fNFakeTracksEta->Fill(dipangle);}
	  Float_t pt=track->GetPt();
	  fPtRes->Fill((pt-ptg)/ptg*100.);
	  break;
	  
	}
    }
  
  
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

void AliL3Evaluate::AssignIDs()
{
  //Assign MC id to the tracks.
  
  UInt_t *index=0,tmp_ind=0;
  Int_t *pID=0,npoints=0;
  
  if(!fIsSlow)
    {
      pID = GetFastIDs(tmp_ind,npoints);
      index = (UInt_t*)tmp_ind;
    }
  
  fTracks->QSort();
  LOG(AliL3Log::kDebug,"AliL3Evaluate::AssignIDs","Track Loop")
    <<"Assigning MC id to the found tracks...."<<ENDLOG;
  for(Int_t i=0; i<fTracks->GetNTracks(); i++)
    {
      AliL3Track *track = (AliL3Track*)fTracks->GetCheckedTrack(i);
      if(!track) 
	{
	  LOG(AliL3Log::kWarning,"AliL3Evaluate::AssignIDs","Track Loop")
	    <<AliL3Log::kDec<<"No track in track array, index "<<i<<ENDLOG;
	  continue;
	}
      if(track->GetNumberOfPoints() < fMinPointsOnTrack) break;
      Int_t tID;
      if(!fIsSlow)
	tID = GetMCTrackLabel(track,index,pID,npoints);
      else
	tID = GetMCTrackLabel(track);
      track->SetMCid(tID);
      printf("track %i id %d\n",i,tID);
    }
  
  if(!fIsSlow)
    {
      delete [] pID;
      delete [] index;
    }
 
}


struct S {Int_t lab; Int_t max;};
Int_t AliL3Evaluate::GetMCTrackLabel(AliL3Track *track,UInt_t *index,Int_t *pID,Int_t npoints) 
{
  //Returns the MCtrackID of the belonging clusters.
  //If MCLabel < 0, means that track is fake.
  //Fake track means that more than 10 percent of clusters are assigned incorrectly.
  
  Int_t num_of_clusters = track->GetNumberOfPoints();
  S *s=new S[num_of_clusters];
  Int_t i;
  for (i=0; i<num_of_clusters; i++) s[i].lab=s[i].max=0;
  UInt_t *hitnum = track->GetHitNumbers();  
  UInt_t id;

  Int_t **trackID;
    
  if(fIsSlow)
    trackID = GetClusterIDs(track);
  else
    trackID = GetClusterIDs(track,index,pID,npoints);
  

  Int_t lab=123456789;
  for (i=0; i<num_of_clusters; i++) 
    {
      //Tricks to get the clusters belonging to this track:
      id = hitnum[i];
      Int_t slice = (id>>25) & 0x7f;
      Int_t patch = (id>>22) & 0x7;
      UInt_t pos = id&0x3fffff;	      
      
      AliL3SpacePointData *points = fClusters[slice][patch];
      
      if(!points) 
	continue;
      
      if(pos>=fNcl[slice][patch]) 
	{
	  LOG(AliL3Log::kError,"AliL3Evaluate::GetMCTrackLabel","Clusterarray")
	    <<AliL3Log::kDec<<"ERROR"<<ENDLOG;
	  continue;
	}
      
      //Get the label of the cluster:
      //printf("label %d %d %d\n",trackID[i][0],trackID[i][1],trackID[i][2]);
      lab=abs(trackID[i][0]);
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
      
      if(!points) 
	continue;
      
      if(pos>=fNcl[slice][patch]) 
	{
	  LOG(AliL3Log::kError,"AliL3Evaluate::GetMCTrackLabel","Clusterarray")
	    <<AliL3Log::kDec<<"ERROR"<<ENDLOG;
	  continue;
	}
      
      if (abs(trackID[i][1]) == lab || 
	  abs(trackID[i][2]) == lab ) max++;
    }
  
  //check if more than 10% of the clusters are incorrectly assigned (fake track):
  if (1.-Float_t(max)/num_of_clusters > 0.10) 
    {
      return -lab;
    }
  
  delete [] trackID;
  return lab;
}


Int_t **AliL3Evaluate::GetClusterIDs(AliL3Track *track,UInt_t *index,Int_t *pID,Int_t npoints)
{
  //Return the MC information of all clusters belonging to track.
  
  Int_t num_of_clusters = track->GetNumberOfPoints();
  Int_t **trackID = new Int_t*[num_of_clusters];
  
  UInt_t *hitnum = track->GetHitNumbers();  
  UInt_t id;

  
  Float_t xyz[3];
  Int_t sector,padrow;
  for(Int_t i=0; i<num_of_clusters; i++)
    {
      id = hitnum[i];
      Int_t slice = (id>>25) & 0x7f;
      Int_t patch = (id>>22) & 0x7;
      UInt_t pos = id&0x3fffff;	      
      
      AliL3SpacePointData *points = fClusters[slice][patch];
      
      if(!points) 
	{
	  LOG(AliL3Log::kError,"AliL3Evaluate::GetClusterIDs","Clusterarray")
	    <<"No points at slice "<<slice<<" patch "<<patch<<" pos "<<pos<<ENDLOG;
	  continue;
	}
      if(pos>=fNcl[slice][patch]) 
	{
	  LOG(AliL3Log::kError,"AliL3Evaluate::GetClusterIDs","Clusterarray")
	    <<AliL3Log::kDec<<"ERROR"<<ENDLOG;
	  continue;
	}
      
      xyz[0] = points[pos].fX;
      xyz[1] = points[pos].fY;
      xyz[2] = points[pos].fZ;
      //sector = points[pos].fSector;
      padrow = points[pos].fPadRow;
      Int_t se,ro;
      fTransform->Slice2Sector(slice,padrow,se,ro);
      fTransform->Global2Raw(xyz,se,ro);
      
      if(fIsSlow)
	{
	  Int_t p = fRowid[slice][padrow];
	  
	  if(!fDigitsTree->GetEvent(p)) 
	    LOG(AliL3Log::kError,"AliL3Evaluate::GetClusterIDs","Digits Tree")
	      <<"Error reading digits tree"<<ENDLOG;
	  
	  trackID[i] = new Int_t[3];
	  trackID[i][0] = fDigits->GetTrackID((Int_t)xyz[2],(Int_t)xyz[1],0);
	  trackID[i][1] = fDigits->GetTrackID((Int_t)xyz[2],(Int_t)xyz[1],1);
	  trackID[i][2] = fDigits->GetTrackID((Int_t)xyz[2],(Int_t)xyz[1],2);
	}
      else
	{
	  Int_t tmp_pid=0;
	  for(Int_t ii=0; ii<npoints; ii++)
	    {
	      tmp_pid = pID[ii];
	      if(index[ii] == id) break;
	    }
	  trackID[i] = new Int_t[3];
	  trackID[i][0] = tmp_pid;
	  trackID[i][1] = -1;
	  trackID[i][2] = -1;
	}
    }
  return trackID;
}

Int_t *AliL3Evaluate::GetFastIDs(UInt_t &tmp_ind,Int_t &npoints)
{
  //Get the MC id of space points in case of using the fast simulator. 

  FILE *infile = fopen("point_mc.dat","r");
  if(!infile) return 0;
  Int_t hitid,hitmc,i;
  
  for(i=0; ; i++)
    if(fscanf(infile,"%d %d",&hitid,&hitmc)==EOF) break;
  npoints = i;
  rewind(infile);
  Int_t *pID = new Int_t[npoints];
  UInt_t *ind = new UInt_t[npoints];
  tmp_ind = (UInt_t)ind;
  
  for(i=0; i<npoints; i++)
    {
      if(fscanf(infile,"%d %d",&hitid,&hitmc)==EOF) break;
      pID[i] = hitmc;
      ind[i] = hitid;
    }
  fclose(infile);

  return pID;

}

void AliL3Evaluate::Write2File(Char_t *outputfile)
{
  //Write histograms to file:
  
  TFile *of = new TFile(outputfile,"RECREATE");
  if(!of->IsOpen())
    {
      LOG(AliL3Log::kError,"AliL3Evaluate::Write2File","File Open")
	<<"Problems opening rootfile"<<ENDLOG;
      return;
    }
  
  of->cd();
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
  delete of;

}

