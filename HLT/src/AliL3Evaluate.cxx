//Author:        Anders Strand Vestbo
//Last Modified: 5.01.2001

#include <stdio.h>
#include <math.h>
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
#include "AliComplexCluster.h"

#include "AliL3Defs.h"
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
  fDigitsFile = NULL;
  fTracks = NULL;
  fMCclusterfile = NULL;
  fNFastPoints = 0;
  fMcIndex = 0;
  fMcId = 0;
  fMinSlice=0;
  fMaxSlice=0;
  fTransform = new AliL3Transform();
}

AliL3Evaluate::AliL3Evaluate(Char_t *mcfile,Char_t *digitsfile,Int_t *slice)
{
  //Normal constructor. Input are the rootfile containing the 
  //original MC information of the simulated event. 

  fDigitsFile = new TFile(digitsfile,"READ");
  fEventFile = new TFile(mcfile,"READ");
  if(!fDigitsFile->IsOpen())
    {
      LOG(AliL3Log::kError,"AliL3Evaluation::AliL3Evaluation","File Open")
	<<"Inputfile "<<mcfile<<" does not exist"<<ENDLOG;
      return;
    }
  
  fParam = (AliTPCParam*)fEventFile->Get("75x40_100x60");
  fTransform = new AliL3Transform();
  
  fMinSlice = slice[0];
  fMaxSlice = slice[1];
  fMinGoodPt = 0.1;
  fNoOverlap = kFALSE;
}

AliL3Evaluate::AliL3Evaluate(Int_t *slice)
{
  //ctor to use if you do not need any rootfile.


  fMinSlice = slice[0];
  fMaxSlice = slice[1];
  fTransform = new AliL3Transform();

}

AliL3Evaluate::~AliL3Evaluate()
{
  if(fDigitsTree) fDigitsTree->Delete();
  if(fDigitsFile) {
    fDigitsFile->Close();
    delete fDigitsFile;
  }
  if(fTransform) delete fTransform;
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

void AliL3Evaluate::Setup(Char_t *trackfile,Char_t *path)
{
  //Read in the hit and track information from produced files.
  
  fGoodFound = 0;
  fGoodGen = 0;
  Char_t fname[256];
  AliL3FileHandler *clusterfile[36][5];
  for(Int_t s=fMinSlice; s<=fMaxSlice; s++)
    {
      for(Int_t p=0; p<5; p++)
	{
	  fClusters[s][p] = 0;
	  clusterfile[s][p] = new AliL3FileHandler();
	  sprintf(fname,"%s/points_%d_%d.raw",path,s,p);
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
  if(!tfile->SetBinaryInput(trackfile)){
    LOG(AliL3Log::kError,"AliL3Evaluation::Setup","File Open")
    <<"Inputfile "<<trackfile<<" does not exist"<<ENDLOG; 
    return;
  }
  fTracks = new AliL3TrackArray();
  tfile->Binary2TrackArray(fTracks);
  tfile->CloseBinaryInput();
  delete tfile;
}

void AliL3Evaluate::SetupSlow(Char_t *trackfile,Char_t *path)
{
  //Setup for using the slow simulator.
  
  fIsSlow = true;
  Setup(trackfile,path);
  
  if(!InitMC())
    LOG(AliL3Log::kError,"AliL3Evaluation::SetupSlow","Digits Tree")
    <<"Error setting up digits tree"<<ENDLOG;
  
}

void AliL3Evaluate::SetupFast(Char_t *trackfile,Char_t *mcClusterfile,Char_t *path)
{
  //Setup for using the fast simulator.

  fIsSlow = false;
  GetFastClusterIDs(path);
  
  fMCclusterfile = new TFile(mcClusterfile);
  if(!fMCclusterfile->IsOpen())
    LOG(AliL3Log::kError,"AliL3Evaluation::SetupFast","File Open")
      <<"Inputfile "<<mcClusterfile<<" does not exist"<<ENDLOG; 

  Setup(trackfile,path);
  InitMC();
}

Bool_t AliL3Evaluate::InitMC()
{
  if(fIsSlow)
    {
      fDigitsFile->cd();
      fDigitsTree = (TTree*)fDigitsFile->Get("TreeD_75x40_100x60");
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
    }
  
  fEventFile->cd();
  AliRun *gAlice = (AliRun*)fEventFile->Get("gAlice");
  if (!gAlice) 
    {
      LOG(AliL3Log::kError,"AliL3Evaluate::SetParticleArray","gAlice")
	<<"AliRun object non existing on file"<<ENDLOG;
      return false;
    }
  
  gAlice->GetEvent(0);
  
  return true;
  
}



void AliL3Evaluate::DefineGoodTracks(Int_t *slice,Int_t *padrow,Int_t good_number)
{
  //Loop over MC particles, and mark the good ones
  //(which the tracker should find...)
  
  AliTPC *TPC = (AliTPC*)gAlice->GetDetector("TPC");
  
  TPC->SetParam(fParam);
  
  Int_t ver = TPC->IsVersion();
  LOG(AliL3Log::kInformational,"AliL3Evaluate::DefineGoodTracks","TPC version")
    <<"TPC version "<<ver<<" found on file"<<ENDLOG;
  
  Int_t zero=TPC->GetParam()->GetZeroSup();
  
  Int_t np = gAlice->GetNtrack();
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
         
      for(Int_t i=0; i<carray.GetTree()->GetEntries(); i++)  
	{
	  Int_t sec,row,sl,lr;
	  AliSegmentID *s = carray.LoadEntry(i);
	  fParam->AdjustSectorRow(s->GetID(),sec,row);
	  fTransform->Sector2Slice(sl,lr,sec,row);
	  
	  if(sl != slice[0]) {carray.ClearRow(sec,row); continue;}
	  if(lr < padrow[0]) {carray.ClearRow(sec,row); continue;}
	  if(lr > padrow[1]) {carray.ClearRow(sec,row); continue;}
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
      for(Int_t sl=slice[0]; sl<=slice[1]; sl++)
	{
	  for (i=padrow[0]; i<=padrow[1]; i++) {
	    Int_t index = fRowid[sl][i];
	    if (!fDigitsTree->GetEvent(index)) continue;
	    Int_t sec,row;
	    fParam->AdjustSectorRow(fDigits->GetID(),sec,row);
	    fDigits->First();
	    do {
	      Int_t it=fDigits->CurrentRow(), ip=fDigits->CurrentColumn();
	      Short_t dig = fDigits->GetDigit(it,ip);
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
		  good[j]++;
		
		count[j]=0;
	      }
	  }
	}
      delete[] count;
    }
  
  else 
    {
      LOG(AliL3Log::kError,"AliL3Evaluation::FillEffHistos","TPC version")
	<<"No valid TPC version found"<<ENDLOG;
      return;
    }
  
  
  //The following code has been taken from offlinemacro->AliTPCComparison.C
  
  TTree *TH=gAlice->TreeH();
  Int_t npart=(Int_t)TH->GetEntries();
  Int_t nt=0;
  Int_t max = 15000;
  while (npart--) {
    AliTPChit *hit0=0;
    
      TPC->ResetHits();
      TH->GetEvent(npart);
      AliTPChit * hit = (AliTPChit*) TPC->FirstHit(-1);
      while (hit){
        if (hit->fQ==0.) break;
        hit =  (AliTPChit*) TPC->NextHit();
      }
      if (hit) {
        hit0 = new AliTPChit(*hit); //Make copy of hit
        hit = hit0;
      }
      else continue;
      AliTPChit *hit1=(AliTPChit*)TPC->NextHit();       
      if (hit1==0) continue;
      if (hit1->fQ != 0.) continue;
      Int_t i=hit->Track();
      TParticle *p = (TParticle*)gAlice->Particle(i);
      
      printf("Checking particle %d with code %d\n",i,p->GetPdgCode());
      if (p->GetFirstMother()>=0) continue;  //secondary particle
      if (good[i] < good_number) continue;
      if (p->Pt()<fMinGoodPt) continue;
      if (TMath::Abs(p->Pz()/p->Pt())>0.999) continue;
      printf("Checking particle %d, nHits %d\n",i,good[i]);
      fGoodGen++;
      fGoodTracks[nt].label=i;
      fGoodTracks[nt].code=p->GetPdgCode();
      //**** px py pz - in global coordinate system, x y z - in local !
      fGoodTracks[nt].px=hit->X(); fGoodTracks[nt].py=hit->Y(); fGoodTracks[nt].pz=hit->Z();
      
      nt++;     
      if (hit0) delete hit0;
      if (nt==max) {cerr<<"Too many good tracks !n"; break;}
   }
  
  delete [] good;
    
}

void AliL3Evaluate::EvaluatePatch(Int_t slice,Int_t patch,Int_t min_points,Int_t good_number)
{
  //Make efficiency plots for tracking on patch level (before any merging).
  
  Int_t row[5][2] = {{ 0, 45},{46,77},{78,109},{110,141},{142,175}};
  Int_t sl[2] ={slice,slice};
  DefineGoodTracks(sl,row[patch],good_number);
  SetMinPoints(min_points);
  AssignIDs();
  CreateHistos();
  FillEffHistos();
  CalcEffHistos();
}

void AliL3Evaluate::EvaluateSlice(Int_t slice,Int_t min_points,Int_t good_number)
{
  //Make efficiency plots for tracking on a slice (after merging).
  //min_points = minimum points on track to be considered for evaluation
  //good_number = minimum hits (padrows) produced by simulated track for consideration.

  Int_t row[2] = {0,175};
  Int_t sl[2] ={slice,slice};
  DefineGoodTracks(sl,row,good_number);

  SetMinPoints(min_points);
  
  AssignIDs();
  CreateHistos();
  FillEffHistos();
  CalcEffHistos();
}

void AliL3Evaluate::EvaluateGlobal(Int_t min_points,Int_t good_number)
{
  //Make efficiency plots for tracking on several slices.
  
  Int_t row[2] = {0,175};
  SetMinPoints(min_points);
  Int_t slice[2] = {fMinSlice,fMaxSlice};
  DefineGoodTracks(slice,row,good_number);
  AssignIDs();
  CreateHistos(20,0,5);
  FillEffHistos();
  CalcEffHistos();
}

void AliL3Evaluate::AssignIDs()
{
  //Assign MC id to the tracks.
  
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
      printf("track %i id %d nHits %d\n",i,tID,track->GetNumberOfPoints());
    }
}


struct S {Int_t lab; Int_t max;};
Int_t AliL3Evaluate::GetMCTrackLabel(AliL3Track *track){ 
  //Returns the MCtrackID of the belonging clusters.
  //If MCLabel < 0, means that track is fake.
  //Fake track means that more than 10 percent of clusters are assigned incorrectly.

  Int_t num_of_clusters = track->GetNumberOfPoints();
  S *s=new S[num_of_clusters];
  Int_t i;
  for (i=0; i<num_of_clusters; i++) s[i].lab=s[i].max=0;
  UInt_t *hitnum = track->GetHitNumbers();  
  UInt_t id;
    
  Int_t **trackID = GetClusterIDs(track);

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
      if(!points) continue; 
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
  
  for(Int_t j=0; j<track->GetNumberOfPoints(); j++)
    delete [] trackID[j];
  delete [] trackID;
  return lab;
}


Int_t **AliL3Evaluate::GetClusterIDs(AliL3Track *track)
{
  //Return the MC information of all clusters belonging to track.

  Int_t num_of_clusters = track->GetNumberOfPoints();
  Int_t **trackID = new Int_t*[num_of_clusters];
  
  UInt_t *hitnum = track->GetHitNumbers();  
  UInt_t id;
  
  Float_t xyz[3];
  Int_t padrow;
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
	  trackID[i][0] = fDigits->GetTrackID((Int_t)rint(xyz[2]),(Int_t)rint(xyz[1]),0);
	  trackID[i][1] = fDigits->GetTrackID((Int_t)rint(xyz[2]),(Int_t)rint(xyz[1]),1);
	  trackID[i][2] = fDigits->GetTrackID((Int_t)rint(xyz[2]),(Int_t)rint(xyz[1]),2);
	  //if(trackID[i][0]==6 || trackID[i][0]==32)
	  //printf("trackID %d, padrow %d pad %d time %d\n",fDigits->GetTrackID((Int_t)rint(xyz[2]),(Int_t)rint(xyz[1]),0),padrow,(int)rint(xyz[1]),(int)rint(xyz[2]));
	  /*if(trackID[i][0]<0)
	    {
	      
	    printf("trackID %d, padrow %d pad %d time %d\n",trackID[i][0],padrow,(int)rint(xyz[1]),(int)rint(xyz[2]));
	    printf("on the side %d %d %d %d\n",fDigits->GetTrackID(((int)rint(xyz[2])-1),((int)rint(xyz[1])),0),fDigits->GetTrackID(((int)rint(xyz[2])+1),((int)rint(xyz[1])),0),fDigits->GetTrackID(((int)rint(xyz[2])),((int)rint(xyz[1])-1),0),fDigits->GetTrackID(((int)rint(xyz[2])),((int)rint(xyz[1])+1),0));
	    }*/
	}
      else
	{
	  Int_t tmp_pid=0;
	  for(Int_t ii=0; ii<fNFastPoints; ii++)
	    {
	      tmp_pid = fMcId[ii];
	      if(fMcIndex[ii] == id) break;
	    }
	  trackID[i] = new Int_t[3];
	  trackID[i][0] = tmp_pid;
	  trackID[i][1] = -1;
	  trackID[i][2] = -1;
	}
    }

  return trackID;
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

void AliL3Evaluate::CreateHistos(Int_t nbin,Int_t xlow,Int_t xup)
{
  //Create the histograms 
  
  fNtuppel = new TNtuple("fNtuppel","Pt resolution","pt_gen:pt_found:nHits");

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

void AliL3Evaluate::FillEffHistos()
{  
  //Fill the efficiency histograms.

  for(Int_t i=0; i<fGoodGen; i++)
    {
      Double_t ptg=fGoodTracks[i].pt,pzg=fGoodTracks[i].pz;
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
  
  TFile *of = new TFile(outputfile,"RECREATE");
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
  delete of;

}

TNtuple *AliL3Evaluate::CalculateResiduals()
{

  TNtuple *ntuppel=new TNtuple("ntuppel","Residuals","residual_trans:residual_long:zHit:pt:dipangle:beta:padrow:nHits");

  for(int f=fMinSlice; f<=fMaxSlice; f++)
    {
      AliL3FileHandler *tfile = new AliL3FileHandler();
      char fname[256];
      sprintf(fname,"tracks_tr_%d_0.raw",f);
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
	      //fTransform->Global2Local(xyz,slice);
	      
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

TNtuple *AliL3Evaluate::EvaluatePoints(Char_t *rootfile)
{
  //Compare points to the exact crossing points of track and padrows.
  //The input file to this function, contains the exact clusters calculated
  //in AliTPC::Hits2ExactClusters.
    
 
  TNtuple *ntuppel = new TNtuple("ntuppel","residuals","slice:padrow:resy:resz:zHit:pt");
  
  TFile *exfile = new TFile(rootfile);
  
  AliTPCParam *param = (AliTPCParam*)exfile->Get("75x40_100x60");
  
  //Get the exact clusters from file:
  AliTPCClustersArray *arr = new AliTPCClustersArray;
  arr->Setup(param);
  arr->SetClusterType("AliComplexCluster");
  char treeName[500];
  sprintf(treeName,"TreeCExact_%s",param->GetTitle());
  Bool_t clusterok = arr->ConnectTree(treeName);//Segment Tree (for offline clusters)
  if(!clusterok) {printf("AliL3Evaluate::EvaluatePoints : Error in clusterloading\n"); return 0;}
  
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
      fTransform->Sector2Slice(slice,padrow,cursec,currow);
      if(slice<fMinSlice || slice>fMaxSlice) continue;
      AliL3SpacePointData *points = fClusters[slice][0];
      
      Int_t index = fRowid[slice][padrow];
      if(!fDigitsTree->GetEvent(index))
	printf("AliL3Evaluate::EvaluatePoints : ERROR IN DIGITSTREE\n");
      printf("Checking slice %d padrow %d with %d clusters\n",slice,padrow,num_of_offline);
      
      for(UInt_t c=0; c<fNcl[slice][0]; c++)
	{
	  if(points[c].fPadRow!=padrow) continue;
	  for(Int_t m=0; m<num_of_offline; m++)
	    {
	      AliComplexCluster *cluster = (AliComplexCluster *)clusters->UncheckedAt(m);
	      Int_t mcId = cluster->fTracks[0];
	      //Int_t mcId = cluster->GetLabel(0);
	      if(mcId <0) continue;
	      TParticle *part = gAlice->Particle(mcId);
	      
	      Float_t xyz_cl[3] = {points[c].fX,points[c].fY,points[c].fZ};
	      Float_t xyz_ex[3];
	      fTransform->Global2Raw(xyz_cl,cursec,currow);
	      if(fDigits->GetTrackID((Int_t)rint(xyz_cl[2]),(Int_t)rint(xyz_cl[1]),0)!=mcId &&
		 fDigits->GetTrackID((Int_t)rint(xyz_cl[2]),(Int_t)rint(xyz_cl[1]),1)!=mcId &&
		 fDigits->GetTrackID((Int_t)rint(xyz_cl[2]),(Int_t)rint(xyz_cl[1]),2)!=mcId)
		continue;
	      fTransform->Raw2Local(xyz_ex,cursec,currow,cluster->fY,cluster->fX);
	      fTransform->Raw2Local(xyz_cl,cursec,currow,xyz_cl[1],xyz_cl[2]);
	      Float_t resy = xyz_cl[1] - xyz_ex[1];//cluster->GetY()
	      Float_t resz = xyz_cl[2] - xyz_ex[2];//cluster->GetZ()
	      
	      ntuppel->Fill(slice,padrow,resy,resz,xyz_ex[2],part->Pt());
	    }
	}      
      arr->ClearRow(cursec,currow);
    }

  return ntuppel;
}
/*
TNtuple *AliL3Evaluate::EvaluateGEANT()
{
  TNtuple *ntuppel = new TNtuple("ntuppel","residuals","resy:ptgen:padrow:zHit:slice");
  
  Int_t good_number=175;
  Int_t row[2] = {0,175};
      
  Float_t xyz_cross[3];
  Float_t xyz_cl[3];
  
  for(Int_t slice=fMinSlice; slice<=fMaxSlice; slice++)
    {
      Int_t *particle_id = new Int_t[fParticles->GetEntriesFast()];
      TObjArray *good_particles = DefineGoodTracks(slice,row,good_number,particle_id);
      
      if(fMCclusterfile)
	{
	  fMCclusterfile->Close();
	  delete fMCclusterfile;
	  fMCclusterfile = new TFile("/nfs/david/subatom/alice/data/V3.04/fast/clusters/hg_8k_v0_s1-3_e0_cl.root");
	  fMCclusterfile->cd();
	}
      AliTPCClustersArray arr;
      arr.Setup(fParam);
      arr.SetClusterType("AliTPCcluster");
      Bool_t clusterok = arr.ConnectTree("Segment Tree");
      if(!clusterok) 
	LOG(AliL3Log::kError,"AliL3Evaluate::DefineGoodTracks","Cluster Array")
	  <<"Error loading clusters from rootfile"<<ENDLOG;
      
      for(Int_t ci=0; ci<arr.GetTree()->GetEntries(); ci++)  
	{
	  Int_t sec,row,sl,lr;
	  AliSegmentID *s = arr.LoadEntry(ci);
	  fParam->AdjustSectorRow(s->GetID(),sec,row);
	  fTransform->Sector2Slice(sl,lr,sec,row);
	  if(sl != slice) {arr.ClearRow(sec,row); continue;}
	  //if(lr!=pr) {carray.ClearRow(sec,row); continue;}
	  AliTPCClustersRow *cRow = arr.GetRow(sec,row);
	  printf("evaluating slice %d row %d\n",slice,lr);
	  for(Int_t j=0; j<cRow->GetArray()->GetEntriesFast(); j++)
	    {
	      
	      AliTPCcluster *cluster=(AliTPCcluster*)(*cRow)[j];
	      xyz_cl[1] = cluster->GetY();
	      xyz_cl[2] = cluster->GetZ();
	      Int_t lab=cluster->GetLabel(0);
	      if(lab<0) continue;
	      
	      for(Int_t pa=0; pa<good_particles->GetEntriesFast(); pa++)
		{
		  if(particle_id[pa]!=lab) continue;
		  TParticle *part = (TParticle*)good_particles->UncheckedAt(pa);
		  GetParticleCrossingPoint(part,slice,lr,xyz_cross);
		  Double_t yres = xyz_cl[1] - xyz_cross[1];
		  ntuppel->Fill(yres,part->Pt(),lr,xyz_cl[2],slice);
		}
	    }
	  if(arr.GetRow(sec,row)) 
	    arr.ClearRow(sec,row);
	}
      delete [] particle_id;
      delete good_particles;
    }
  return ntuppel;
}
*/
Bool_t AliL3Evaluate::GetParticleCrossingPoint(TParticle *part,Int_t slice,Int_t padrow,Float_t *xyz)
{
  //Calcluate the crossing point between a generated particle and given padrow.
  

  Double_t kappa = 0.2*0.0029980/part->Pt();
  
  Double_t radius = 1/fabs(kappa);
  if(part->GetPdgCode() > 0) kappa = -kappa;
  
  Int_t charg;
  if(kappa>0)
    charg=-1;
  else
    charg=1;
  
  Float_t angl[1] = {part->Phi()};
  
  fTransform->Global2LocalAngle(angl,slice);
  
  Double_t charge = -1.*kappa;
  Double_t trackPhi0 = angl[0] + charge*0.5*Pi/fabs(charge);
  
  Double_t x0=0;
  Double_t y0=0;
  Double_t xc = x0 - radius * cos(trackPhi0);
  Double_t yc = y0 - radius * sin(trackPhi0);
  
  //printf("radius %f xc %f yc %f\n",radius,xc,yc);

  Double_t xHit = fTransform->Row2X(padrow);
  xyz[0] = xHit;
  Double_t aa = (xHit - xc)*(xHit - xc);
  Double_t r2 = radius*radius;
  if(aa > r2)
    return false;

  Double_t aa2 = sqrt(r2 - aa);
  Double_t y1 = yc + aa2;
  Double_t y2 = yc - aa2;
  xyz[1] = y1;
  if(fabs(y2) < fabs(y1)) xyz[1] = y2;

  Double_t tgl = part->Pz()/part->Pt();

  Double_t yHit = xyz[1];
  Double_t angle1 = atan2((yHit - yc),(xHit - xc));
  if(angle1 < 0) angle1 += 2.*Pi;
  Double_t angle2 = atan2((0 - yc),(0 - xc));
  if(angle2 < 0) angle2 += 2.*Pi;
  Double_t diff_angle = angle1 - angle2;
  diff_angle = fmod(diff_angle,2*Pi);
  if((charg*diff_angle) > 0) diff_angle = diff_angle - charg*2.*Pi;
  Double_t s_tot = fabs(diff_angle)*radius;
  Double_t zHit = 0 + s_tot*tgl;
  xyz[2] = zHit;

  return true;
}

