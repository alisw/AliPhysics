//Author:        Anders Strand Vestbo
//Last Modified: 28.6.01

#include <TClonesArray.h>
#include <TH2.h>
#include <TH1.h>
#include <TFile.h>
#include <AliRun.h>
#include <TParticle.h>
#include <TTree.h>

#include "AliTPCParam.h"
#include "AliSimDigits.h"
#include "AliL3TrackArray.h"
#include "AliL3Transform.h"
#include "AliL3HoughTransformer.h"
#include "AliL3Defs.h"
#include "AliL3HoughTrack.h"
#include "AliL3HoughEval.h"

ClassImp(AliL3HoughEval)


AliL3HoughEval::AliL3HoughEval()
{
  //Default constructor
  fTransform = new AliL3Transform();
  fHoughTransformer = 0;
  fNumOfRowsToMiss = 1;
  fNumOfPadsToLook = 1;
}


AliL3HoughEval::AliL3HoughEval(AliL3HoughTransformer *transformer)
{
  //Constructor
  fHoughTransformer = transformer;
  fTransform = new AliL3Transform();
  fNumOfRowsToMiss = 1;
  fNumOfPadsToLook = 1;
}


AliL3HoughEval::~AliL3HoughEval()
{
  //Destructor
  if(fTransform)
    delete fTransform;
  if(fMcTrackTable)
    delete [] fMcTrackTable;
}


Bool_t AliL3HoughEval::LookInsideRoad(AliL3HoughTrack *track,Int_t eta_index,Bool_t remove)
{
  //Look at rawdata along the road specified by the track candidates.
  
  if(!fHoughTransformer)
    {
      printf("\nAliL3HoughEval: No transformer object\n");
      return kFALSE;
    }
  
  Int_t patch = fHoughTransformer->fPatch;
  Int_t slice = fHoughTransformer->fSlice;

  
  Int_t sector,row;
  Int_t lut_index;
  Char_t **track_lut = fHoughTransformer->fTrackTable;
  Int_t nrow=0,npixs=0;
  Float_t xyz[3];
    
  for(Int_t padrow = NRows[patch][0]; padrow <= NRows[patch][1]; padrow++)
    {
      Int_t prow = padrow - NRows[patch][0];
      if(!track->GetCrossingPoint(padrow,xyz))  
	{
	  printf("AliL3HoughEval::LookInsideRoad : Track does not cross line!!\n");
	  continue;
	}
      
      fTransform->Slice2Sector(slice,padrow,sector,row);
      fTransform->Local2Raw(xyz,sector,row);
      npixs=0;
      //Look at both sides of the crossing point:
      for(Int_t p=(Int_t)xyz[1]-fNumOfPadsToLook; p<=(Int_t)xyz[1]+fNumOfPadsToLook; p++)
	{
	  if(p<0 || p>fTransform->GetNPads(padrow)) continue;
	  lut_index = (prow<<8) + p;
	  if(track_lut[eta_index][lut_index]>0) //There was a signal here
	    npixs++;
	}
      
      if(npixs > 0)
	{
	  nrow++;
	}	  
    }
  if(nrow >= NRows[patch][1]-NRows[patch][0]-fNumOfRowsToMiss)//this was a good track
    {
      track->SetEtaIndex(eta_index);
      if(remove)
	RemoveTrackFromImage(track,eta_index);
      return kTRUE;
    }
  else
    return kFALSE;
}

void AliL3HoughEval::LookInsideRawRoad(AliL3TrackArray *tracks,Int_t eta_index,Bool_t remove)
{
  //Evalaute the track candidates by looking along the trajectory.
  //If remove is on, the pixels along the track will be removed. 


  if(!fHoughTransformer)
    {
      printf("AliL3HoughEval: No transformer object\n");
      return;
    }
  
  Int_t patch = fHoughTransformer->fPatch;
  Int_t slice = fHoughTransformer->fSlice;

  
  Int_t sector,row;
  Int_t lut_index;
  Char_t **track_lut = fHoughTransformer->fTrackTable;
  Int_t nrow,npixs;
  Float_t xyz[3];
  
  for(Int_t i=0; i<tracks->GetNTracks(); i++)
    {
      AliL3HoughTrack *track = (AliL3HoughTrack*)tracks->GetCheckedTrack(i);
      if(!track) {printf("No track\n"); return;}
      
      nrow = 0;
      
      for(Int_t padrow = NRows[patch][0]; padrow <= NRows[patch][1]; padrow++)
	{
	  Int_t prow = padrow - NRows[patch][0];
	  if(!track->GetCrossingPoint(padrow,xyz))  
	    {
	      printf("AliL3HoughEval::LookInsideRoad : Track does not cross line!!\n");
	      continue;
	    }
	  	  
	  fTransform->Slice2Sector(slice,padrow,sector,row);
	  fTransform->Local2Raw(xyz,sector,row);
	  npixs=0;
	  //Look at both sides of the crossing point:
	  for(Int_t p=(Int_t)xyz[1]-fNumOfPadsToLook; p<=(Int_t)xyz[1]+fNumOfPadsToLook; p++)
	    {
	      if(p<0 || p>fTransform->GetNPads(padrow)) continue;
	      lut_index = (prow<<8) + p;
	      if(track_lut[eta_index][lut_index]>0) //There was a signal here
		npixs++;
	    }
	  
	  if(npixs > 0)
	    {
	      nrow++;
	    }	  
	}
      if(nrow < NRows[patch][1]-NRows[patch][0]-fNumOfRowsToMiss)
	tracks->Remove(i); //this was not a good enough track
      else if(remove)
	RemoveTrackFromImage(track,eta_index); //this was a good track, so remove it from the image.
    }
  
  tracks->Compress();
}

void AliL3HoughEval::RemoveTrackFromImage(AliL3HoughTrack *track,Int_t eta_index)
{
  //Remove the pixels along the track in the image. 
  
  Int_t patch = fHoughTransformer->fPatch;
  Int_t slice = fHoughTransformer->fSlice;

  Int_t maxnpads = fTransform->GetNPads(NRows[patch][1]);
  
  Int_t lut_index;
  Char_t **track_lut = fHoughTransformer->fTrackTable;
  Int_t sector,row;
  Float_t xyz[3];
  for(Int_t padrow = NRows[patch][0]; padrow<=NRows[patch][1]; padrow++)
    {
      Int_t prow = padrow - NRows[patch][0];
      if(!track->GetCrossingPoint(padrow,xyz))  
	{
	  printf("AliL3HoughEval::LookInsideRoad : Track does not cross line!!\n");
	  continue;
	}
      
      fTransform->Slice2Sector(slice,padrow,sector,row);
      fTransform->Local2Raw(xyz,sector,row);
      
      for(Int_t p=(Int_t)xyz[1]-fNumOfPadsToLook; p<=(Int_t)xyz[1]+fNumOfPadsToLook; p++)
	{
	  if(p<0 || p>fTransform->GetNPads(padrow)) continue;
	  lut_index = (prow<<8) + p;
	  track_lut[eta_index][lut_index] = -1; //remove it!!
	}
    }
  
  
}

void AliL3HoughEval::DisplaySlice(TH2F *hist)
{

  AliL3Digits *pixel;
  
  for(Int_t padrow = 0; padrow < 174; padrow++)
    {
      
      for(pixel=(AliL3Digits*)fHoughTransformer->fRowContainer[padrow].first; pixel!=0; pixel=(AliL3Digits*)pixel->nextRowPixel)
	{
	  
	  Int_t sector,row;
	  Float_t xyz_pix[3];
	  fTransform->Slice2Sector(fHoughTransformer->fSlice,padrow,sector,row);
	  fTransform->Raw2Local(xyz_pix,sector,row,pixel->fPad,pixel->fTime); //y alone pad
	  
	  if(hist)
	    {
	      hist->Fill(xyz_pix[0],xyz_pix[1],pixel->fCharge);
	    }    
	  
	}
    }
  
}

void AliL3HoughEval::DefineGoodParticles(Char_t *rootfile,Double_t pet)
{
  //define the particles that produce good enough signals to be recognized in the transform

  Int_t num_eta_segments = fHoughTransformer->fNumEtaSegments;
  Double_t eta_max = fHoughTransformer->fEtaMax;
  fMcTrackTable = new Int_t[num_eta_segments];
  for(Int_t i=0; i<num_eta_segments; i++)
    fMcTrackTable[i]=0;
  Double_t etaslice = eta_max/num_eta_segments;

  Int_t patch = fHoughTransformer->fPatch;
  Int_t slice = fHoughTransformer->fSlice;
  
  TFile *file = new TFile(rootfile);
  file->cd();

  AliRun *gAlice = (AliRun*)file->Get("gAlice");
  gAlice->GetEvent(0);
  
  TClonesArray *particles=gAlice->Particles();
  Int_t np = particles->GetEntriesFast();
  
  Int_t row_to_miss = fNumOfRowsToMiss;
  AliTPCParam *param = (AliTPCParam*)file->Get("75x40_100x60");
  Int_t zero = param->GetZeroSup();
  TTree *TD=(TTree*)gDirectory->Get("TreeD_75x40_100x60");
  AliSimDigits da, *digits=&da;
  TD->GetBranch("Segment")->SetAddress(&digits); //Return pointer to branch segment.
  Int_t *good=new Int_t[np];
  Int_t *count = new Int_t[np]; //np number of particles.
  Int_t good_number = NRows[patch][1]-NRows[patch][0]-row_to_miss;
  Int_t i;
  for (i=0; i<np; i++) count[i]=0;
  Int_t sectors_by_rows=(Int_t)TD->GetEntries();
  for (i=0; i<sectors_by_rows; i++) 
    {
      if (!TD->GetEvent(i)) continue;
      Int_t sec,row;
      param->AdjustSectorRow(digits->GetID(),sec,row);
      Int_t ss,sr;
      fTransform->Sector2Slice(ss,sr,sec,row);
      if(ss!=slice) continue;
      if(sr < NRows[patch][0]) continue;
      if(sr > NRows[patch][1]) break;
      digits->First();
      while (digits->Next()) {
	Int_t it=digits->CurrentRow(), ip=digits->CurrentColumn();
	Short_t dig = digits->GetDigit(it,ip);
	Int_t idx0=digits->GetTrackID(it,ip,0); 
	Int_t idx1=digits->GetTrackID(it,ip,1);
	Int_t idx2=digits->GetTrackID(it,ip,2);
	if (idx0>=0 && dig>=zero) count[idx0]+=1;
	if (idx1>=0 && dig>=zero) count[idx1]+=1;
	if (idx2>=0 && dig>=zero) count[idx2]+=1;
      }
      for (Int_t j=0; j<np; j++) 
	{
	  if (count[j]>1) {//at least two digits at this padrow
	    good[j]++;
	  }
	  count[j]=0;
	}
    }
  delete[] count;
  
  Int_t good_one=0;
  //TObjArray *part = new TObjArray(0,0);
  for(i=0; i<np; i++)
   {
     TParticle *p = (TParticle*)particles->UncheckedAt(i);
     if(p->GetFirstMother()>0) continue; //secondary particle
     if(good[i] < good_number) continue; //too few padrows
     Double_t ptg=p->Pt();
     if(ptg<pet) continue;
     Int_t eta_index = (Int_t)(p->Eta()/etaslice);
     if(eta_index < 0 || eta_index > num_eta_segments)
       continue;
     //fMcTrackTable[eta_index]++;
     //part->AddLast(p);
     good_one++;
   }
  
  printf("nparticles %d\n",good_one);
  file->Close();
  delete [] good;
  delete file;
  //return part;
}


void AliL3HoughEval::CompareMC(Char_t *rootfile,AliL3TrackArray *merged_tracks,Float_t *eta)
{
  
  Int_t slice = fHoughTransformer->fSlice;
  
  TFile *file = new TFile(rootfile);
  file->cd();
  
  AliRun *gAlice = (AliRun*)file->Get("gAlice");
  gAlice->GetEvent(0);
  
  TClonesArray *particles=gAlice->Particles();  
  Int_t n=particles->GetEntriesFast();
  Float_t torad=TMath::Pi()/180;
  Float_t phi_min = slice*20 - 10;
  Float_t phi_max = slice*20 + 10;
  
  for (Int_t j=0; j<n; j++) {
    TParticle *p=(TParticle*)particles->UncheckedAt(j);
    if (p->GetFirstMother()>=0) continue;  //secondary particle
    if(p->Eta() < eta[0] || p->Eta() > eta[1]) continue;
    Double_t ptg=p->Pt(),pxg=p->Px(),pyg=p->Py();//,pzg=p->Pz();
    Double_t phi_part = TMath::ATan2(pyg,pxg);
    if (phi_part < 0) phi_part += 2*TMath::Pi();
    
    if(phi_part < phi_min*torad || phi_part > phi_max*torad) {continue;}
    if(ptg<0.100) continue;
    
    AliL3HoughTrack *sel_track=0;
    
    Double_t min_dist = 10000;
    
    for(Int_t t=0; t<merged_tracks->GetNTracks(); t++)
      {
	AliL3HoughTrack *track = (AliL3HoughTrack*)merged_tracks->GetCheckedTrack(t);
	if(!track) {printf("AliL3HoughEval: NO TRACK\n"); continue;}
	Float_t phi0[2] = {track->GetPhi0(),0};
	fTransform->Local2GlobalAngle(phi0,slice);
	Double_t dPt = ptg - track->GetPt();
	Double_t dPhi0 = phi_part - phi0[0];
	Double_t distance = pow(dPt,2) + pow(dPhi0,2);
	if(distance < min_dist)
	  {
	    min_dist = distance;
	    sel_track = track;
	  }      
	
      }
    if(sel_track)
      sel_track->SetBestMCid(j,min_dist);
    
    Double_t dpt = fabs(ptg-sel_track->GetPt());
    Double_t dphi0 = fabs(phi_part-sel_track->GetPhi0());
    //printf("Found match, min_dist %f dPt %f dPhi0 %f\n",min_dist,dpt,dphi0);
  }
  file->Close();
  delete file;
  
}
