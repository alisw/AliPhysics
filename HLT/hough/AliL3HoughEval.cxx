#include <TClonesArray.h>
#include <TH2.h>
#include <TH1.h>
#include <TFile.h>
#include <AliRun.h>
#include <TParticle.h>

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
  fHoughTransformer = NULL;
}


AliL3HoughEval::AliL3HoughEval(AliL3HoughTransformer *transformer)
{
  //Constructor
  fHoughTransformer = transformer;
  fTransform = new AliL3Transform();
}


AliL3HoughEval::~AliL3HoughEval()
{
  //Destructor
  if(fTransform)
    delete fTransform;
 
}


void AliL3HoughEval::LookInsideRoad(AliL3TrackArray *tracks,TH2F *hist,TH2F *fake)
{
  //Look at rawdata along the road specified by the track candidates.

  if(!fHoughTransformer)
    {
      printf("AliL3HoughEval: No transformer object\n");
      return;
    }
  
  AliL3Digits *pixel;
  Int_t patch = fHoughTransformer->fPatch;
  Int_t slice = fHoughTransformer->fSlice;
  Int_t num_of_pads_to_look = 1;
  Int_t rows_to_miss = 1;

  for(Int_t i=0; i<tracks->GetNTracks(); i++)
    {
      AliL3HoughTrack *track = (AliL3HoughTrack*)tracks->GetCheckedTrack(i);
      if(!track) {printf("No track\n"); return;}
      
      Int_t nrow = 0;
      
      for(Int_t padrow = NRows[patch][0]; padrow <= NRows[patch][1]; padrow++)
	{
	  
	  Float_t xyz[3];
	  //Get the crossing point of the track and current padrow:
	  //if(!track->GetCrossingPoint(slice,padrow,xyz))
	  if(!track->GetCrossingPoint(padrow,xyz))  
	    {
	      printf("AliL3HoughEval::LookInsideRoad : Track does not cross line!!\n");
	      continue;
	    }
	  
	  if(fake)
	    fake->Fill(xyz[0],xyz[1],1);
	  Int_t npixs = 0;
	  
	  //Get the pixels along the track candidate
	  for(pixel=(AliL3Digits*)fHoughTransformer->fRowContainer[padrow].first; pixel!=0; pixel=(AliL3Digits*)pixel->nextRowPixel)
	    {
	      
	      Int_t sector,row;
	      Float_t xyz_pix[3];
	      fTransform->Slice2Sector(slice,padrow,sector,row);
	      fTransform->Raw2Local(xyz_pix,sector,row,pixel->fPad,pixel->fTime); //y alone pad
	      
	      
	      //check if we are inside road
	      if(fabs(xyz_pix[1] - xyz[1]) > num_of_pads_to_look*fTransform->GetPadPitchWidthLow()) continue; 
	      npixs++;
	      
	      if(hist)
		{
		  //fTransform->Local2Global(xyz_pix,slice);
		  hist->Fill(xyz_pix[0],xyz_pix[1],pixel->fCharge);
		}    
	      
	    }
	  if(npixs > 0)
	    {
	      nrow++;
	    }	  
	}
      
      if(nrow < NRows[patch][1]-NRows[patch][0]-rows_to_miss)
	tracks->Remove(i);
	      
    }
  
  tracks->Compress();

  
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

TClonesArray *AliL3HoughEval::GetParticles(Char_t *rootfile)
{
  
  TFile *file = new TFile(rootfile);
  file->cd();

  AliRun *gAlice = (AliRun*)file->Get("gAlice");
  gAlice->GetEvent(0);
  
  TClonesArray *particles=gAlice->Particles();
  return particles;

  file->Close();
  delete file;
}


void AliL3HoughEval::CompareMC(Char_t *rootfile,AliL3TrackArray *merged_tracks,Float_t *eta)
{
  
  Int_t slice = fSlice;
  
  TClonesArray *particles = GetParticles(rootfile);
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
    
}
