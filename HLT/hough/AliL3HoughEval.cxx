// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliL3StandardIncludes.h"

#ifdef use_root
#include <TH1.h>
#include <TFile.h>
#endif

#include "AliL3Logging.h"
#include "AliL3HoughEval.h"
#include "AliL3MemHandler.h"
#include "AliL3TrackArray.h"
#include "AliL3HoughBaseTransformer.h"
#include "AliL3DigitData.h"
#include "AliL3HoughTrack.h"
#include "AliL3Transform.h"
#include "AliL3Histogram.h"
#include "AliL3Histogram1D.h"

#if __GNUC__ == 3
using namespace std;
#endif

/** /class AliL3HoughEval
//<pre>
//_____________________________________________________________
// AliL3HoughEval
//
// Evaluation class for tracklets produced by the Hough transform.
//
</pre>
*/

ClassImp(AliL3HoughEval)

AliL3HoughEval::AliL3HoughEval()
{
  
  fRemoveFoundTracks = kFALSE;
  fNumOfPadsToLook = 1;
  fNumOfRowsToMiss = 1;
  fEtaHistos=0;
  fRowPointers = 0;
}


AliL3HoughEval::~AliL3HoughEval()
{
  fHoughTransformer = 0;
  if(fRowPointers)
    {
      for(Int_t i=0; i<fNrows; i++)
	fRowPointers[i] = 0;
      delete [] fRowPointers;
    }
}

void AliL3HoughEval::InitTransformer(AliL3HoughBaseTransformer *transformer)
{
  fHoughTransformer = transformer;
  fSlice = fHoughTransformer->GetSlice();
  fPatch = fHoughTransformer->GetPatch();
  fNrows = AliL3Transform::GetLastRow(fPatch) - AliL3Transform::GetFirstRow(fPatch) + 1;
  fNEtaSegments = fHoughTransformer->GetNEtaSegments();
  fEtaMin = fHoughTransformer->GetEtaMin();
  fEtaMax = fHoughTransformer->GetEtaMax();
  fZVertex = fHoughTransformer->GetZVertex();
  GenerateLUT();
}

void AliL3HoughEval::GenerateLUT()
{
  //Generate a Look-up table, to limit the access to raw data
  
  if(!fRowPointers)
    fRowPointers = new AliL3DigitRowData*[fNrows];

  AliL3DigitRowData *tempPt = (AliL3DigitRowData*)fHoughTransformer->GetDataPointer();
  if(!tempPt)
    printf("\nAliL3HoughEval::GenerateLUT : Zero data pointer\n");
  
  for(Int_t i=AliL3Transform::GetFirstRow(fPatch); i<=AliL3Transform::GetLastRow(fPatch); i++)
    {
      Int_t prow = i - AliL3Transform::GetFirstRow(fPatch);
      fRowPointers[prow] = tempPt;
      AliL3MemHandler::UpdateRowPointer(tempPt);
    }
  
}

Bool_t AliL3HoughEval::LookInsideRoad(AliL3HoughTrack *track,Int_t &nrows_crossed,Int_t *rowrange,Bool_t remove)
{
  //Look at rawdata along the road specified by the track candidates.
  //If track is good, return true, if not return false.
  
  Int_t sector,row;
  
  Int_t nrow=0,npixs=0;//,rows_crossed=0;
  Float_t xyz[3];
  
  Int_t total_charge=0;//total charge along the road
  
  //for(Int_t padrow = AliL3Transform::GetFirstRow(fPatch); padrow <= AliL3Transform::GetLastRow(fPatch); padrow++)
  for(Int_t padrow = rowrange[0]; padrow<=rowrange[1]; padrow++)
    {
      Int_t prow = padrow - AliL3Transform::GetFirstRow(fPatch);
      if(track->IsHelix())
	{
	  if(!track->GetCrossingPoint(padrow,xyz))  
	    {
	      continue;
	    }
	}
      else
	{
	  track->GetLineCrossingPoint(padrow,xyz);
	  xyz[0] += AliL3Transform::Row2X(track->GetFirstRow());
	  Float_t R = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]);
	  xyz[2] = R*track->GetTgl();
	}
      
      AliL3Transform::Slice2Sector(fSlice,padrow,sector,row);
      AliL3Transform::Local2Raw(xyz,sector,row);

      npixs=0;
      
      //Get the timebins for this pad
      AliL3DigitRowData *tempPt = fRowPointers[prow];
      if(!tempPt) 
	{
	  printf("AliL3HoughEval::LookInsideRoad : Zero data pointer\n");
	  continue;
	}
      
      //Look at both sides of the pad:
      for(Int_t p=(Int_t)rint(xyz[1])-fNumOfPadsToLook; p<=(Int_t)rint(xyz[1])+fNumOfPadsToLook; p++)
	{
	  AliL3DigitData *digPt = tempPt->fDigitData;
	  for(UInt_t j=0; j<tempPt->fNDigit; j++)
	    {
	      Int_t pad = digPt[j].fPad;
	      Int_t charge = digPt[j].fCharge;
	      if(charge <= fHoughTransformer->GetLowerThreshold()) continue;
	      if(pad < p) continue;
	      if(pad > p) break;
	      UShort_t time = digPt[j].fTime;
	      Double_t eta = AliL3Transform::GetEta(fSlice,padrow,pad,time);
	      Int_t pixel_index = fHoughTransformer->GetEtaIndex(eta);
	      if(pixel_index != track->GetEtaIndex()) continue;
	      total_charge += digPt[j].fCharge;
	      if(remove)
		digPt[j].fCharge = 0; //Erease the track from image
	      npixs++;
	    }
	}
            
      if(npixs > 1)//At least 2 digits on this padrow
	{
	  nrow++;
	}	  
    }
  if(remove)
    return kTRUE;
  
  nrows_crossed += nrow; //Update the number of rows crossed.
  
  if(nrow >= rowrange[1]-rowrange[0]+1 - fNumOfRowsToMiss)//this was a good track
    {
      if(fRemoveFoundTracks)
	{
	  Int_t dummy=0;
	  LookInsideRoad(track,dummy,rowrange,kTRUE);
	}
      return kTRUE;
    }
  else
    return kFALSE;
}

void AliL3HoughEval::FindEta(AliL3TrackArray *tracks)
{
  
  Int_t sector,row;
  Float_t xyz[3];
  
  Int_t ntracks = tracks->GetNTracks();
  fEtaHistos = new AliL3Histogram1D*[ntracks];
  
  Char_t hname[100];
  for(Int_t i=0; i<ntracks; i++)
    {
      sprintf(hname,"etahist_%d",i);
      fEtaHistos[i] = new AliL3Histogram1D(hname,hname,100,0,1);
    }
  Double_t etaslice = (fEtaMax - fEtaMin)/fNEtaSegments;
  
  for(Int_t ntr=0; ntr<ntracks; ntr++)
    {
      AliL3HoughTrack *track = (AliL3HoughTrack*)tracks->GetCheckedTrack(ntr);
      if(!track) continue;
      for(Int_t padrow = AliL3Transform::GetFirstRow(fPatch); padrow <= AliL3Transform::GetLastRow(fPatch); padrow++)
	{
	  Int_t prow = padrow - AliL3Transform::GetFirstRow(fPatch);
	  
	  if(!track->GetCrossingPoint(padrow,xyz))  
	    {
	      printf("AliL3HoughEval::LookInsideRoad : Track does not cross line!!\n");
	      continue;
	    }
	  
	  AliL3Transform::Slice2Sector(fSlice,padrow,sector,row);
	  AliL3Transform::Local2Raw(xyz,sector,row);
	  
	  //Get the timebins for this pad
	  AliL3DigitRowData *tempPt = fRowPointers[prow];
	  if(!tempPt) 
	    {
	      printf("AliL3HoughEval::LookInsideRoad : Zero data pointer\n");
	      continue;
	    }
	  
	  //Look at both sides of the pad:
	  for(Int_t p=(Int_t)rint(xyz[1])-fNumOfPadsToLook; p<=(Int_t)rint(xyz[1])+fNumOfPadsToLook; p++)
	    {
	      AliL3DigitData *digPt = tempPt->fDigitData;
	      for(UInt_t j=0; j<tempPt->fNDigit; j++)
		{
		  UChar_t pad = digPt[j].fPad;
		  Int_t charge = digPt[j].fCharge;
		  if(charge <= fHoughTransformer->GetLowerThreshold()) continue;
		  if(pad < p) continue;
		  if(pad > p) break;
		  UShort_t time = digPt[j].fTime;
		  Double_t eta = AliL3Transform::GetEta(fSlice,padrow,pad,time);
		  Int_t pixel_index = (Int_t)(eta/etaslice);
		  if(pixel_index > track->GetEtaIndex()+1) continue;
		  if(pixel_index < track->GetEtaIndex()-1) break;
		  fEtaHistos[ntr]->Fill(eta,digPt[j].fCharge);
		}
	    }
	}
    }
  
  for(Int_t i=0; i<ntracks; i++)
    {
      AliL3Histogram1D *hist = fEtaHistos[i];
      Int_t max_bin = hist->GetMaximumBin();
      Double_t max_value = hist->GetBinContent(max_bin);
      AliL3HoughTrack *track = (AliL3HoughTrack*)tracks->GetCheckedTrack(i);
      if(!track) continue;
      if(hist->GetBinContent(max_bin-1)<max_value && hist->GetBinContent(max_bin+1)<max_value)
	{
	  track->SetWeight((Int_t)max_value,kTRUE); 
	  track->SetEta(hist->GetBinCenter(max_bin));
	  track->SetNHits(track->GetWeight());
	}
      else
	{
	  track->SetWeight(0);
	  tracks->Remove(i); //remove this track, because it was not a peak
	}    
    }
  tracks->Compress();
  
  //for(Int_t i=0; i<ntracks; i++)
  //delete fEtaHistos[i];
  //delete [] fEtaHistos;
}

void AliL3HoughEval::DisplayEtaSlice(Int_t eta_index,AliL3Histogram *hist)
{
  //Display the current raw data inside the (slice,patch)

  if(!hist)
    {
      printf("AliL3HoughEval::DisplayEtaSlice : No input histogram!\n");
      return;
    }
  
  for(Int_t padrow = AliL3Transform::GetFirstRow(fPatch); padrow <= AliL3Transform::GetLastRow(fPatch); padrow++)
    {
      Int_t prow = padrow - AliL3Transform::GetFirstRow(fPatch);
                  
      AliL3DigitRowData *tempPt = fRowPointers[prow];
      if(!tempPt) 
	{
	  printf("AliL3HoughEval::DisplayEtaSlice : Zero data pointer\n");
	  continue;
	}
      
      AliL3DigitData *digPt = tempPt->fDigitData;
      if((Int_t)tempPt->fRow != padrow)
	{
	  printf("\nAliL3HoughEval::DisplayEtaSlice : Mismatching padrows!!!\n");
	  return;
	}
      for(UInt_t j=0; j<tempPt->fNDigit; j++)
	{
	  UChar_t pad = digPt[j].fPad;
	  UChar_t charge = digPt[j].fCharge;
	  UShort_t time = digPt[j].fTime;
	  if((Int_t)charge <= fHoughTransformer->GetLowerThreshold() || (Int_t)charge >= fHoughTransformer->GetUpperThreshold()) continue;
	  Float_t xyz[3];
	  Int_t sector,row;
	  AliL3Transform::Slice2Sector(fSlice,padrow,sector,row);
	  AliL3Transform::Raw2Local(xyz,sector,row,pad,time);
	  xyz[2] -= fZVertex;
	  Double_t eta = AliL3Transform::GetEta(xyz);
	  Int_t pixel_index = fHoughTransformer->GetEtaIndex(eta);//(Int_t)(eta/etaslice);
	  if(pixel_index != eta_index) continue;
	  hist->Fill(xyz[0],xyz[1],charge);
	}
    }
  
}

#ifdef use_root
void AliL3HoughEval::CompareMC(AliL3TrackArray */*tracks*/,Char_t */*trackfile*/,Int_t /*threshold*/)
{
  /*  
  struct GoodTrack goodtracks[15000];
  Int_t nt=0;
  ifstream in(trackfile);
  if(in)
    {
      printf("Reading good tracks from file %s\n",trackfile);
      while (in>>goodtracks[nt].label>>goodtracks[nt].code>>
	     goodtracks[nt].px>>goodtracks[nt].py>>goodtracks[nt].pz>>
	     goodtracks[nt].pt>>goodtracks[nt].eta>>goodtracks[nt].nhits) 
	{
	  nt++;
	  if (nt==15000) 
	    {
	      cerr<<"Too many good tracks"<<endl;
	      break;
	    }
	}
      if (!in.eof())
	{
	  LOG(AliL3Log::kError,"AliL3HoughEval::CompareMC","Input file")
	    <<"Error in file reading"<<ENDLOG;
	  return;
	}
    }
  else
    {
      LOG(AliL3Log::kError,"AliL3HoughEval::CompareMC","Input")
	<<"No input trackfile "<<trackfile<<ENDLOG;
    }
  
  Int_t *particles = new Int_t[fNEtaSegments];
  Int_t *ftracks = new Int_t[fNEtaSegments];
  for(Int_t i=0; i<fNEtaSegments; i++)
    {
      particles[i]=0;
      ftracks[i]=0;
    }
  
  TH1F *ptgood = new TH1F("ptgood","ptgood",5,0,2);
  TH1F *ptfound = new TH1F("ptfound","ptgood",5,0,2);
  TH1F *pteff = new TH1F("pteff","pteff",5,0,2);
  TH1F *etafound = new TH1F("etafound","etafound",5,0,1);
  TH1F *etagood = new TH1F("etagood","etagood",5,0,1);
  TH1F *etaeff = new TH1F("etaeff","etaeff",5,0,1);
  
  Double_t etaslice = (fEtaMax - fEtaMin)/fNEtaSegments;
  for(Int_t i=0; i<tracks->GetNTracks(); i++)
    {
      AliL3HoughTrack *tr = (AliL3HoughTrack*)tracks->GetCheckedTrack(i);
      if(!tr) continue;
      if(tr->GetWeight()<threshold) continue;
      Int_t trackindex = tr->GetEtaIndex();
      if(trackindex <0 || trackindex >= fNEtaSegments) continue;
      ftracks[trackindex]++;
      ptfound->Fill(tr->GetPt());
      etafound->Fill(tr->GetEta());
    }
  for(Int_t i=0; i<nt; i++)
    {
      if(goodtracks[i].nhits < 174) continue;
      if(goodtracks[i].pt < 0.2) continue;
      Int_t particleindex = (Int_t)(goodtracks[i].eta/etaslice);
      if(particleindex < 0 || particleindex >= fNEtaSegments) continue;
      particles[particleindex]++;
      ptgood->Fill(goodtracks[i].pt);
      etagood->Fill(goodtracks[i].eta);
    }
  
  Double_t found=0;
  Double_t good =0;
  for(Int_t i=0; i<fNEtaSegments; i++)
    {
      //printf("Slice %d : Found tracks %d, good tracks %d\n",i,ftracks[i],particles[i]);
      found += ftracks[i];
      good += particles[i];
    }
  printf("And the total efficiency was: %f\n",found/good);

  ptgood->Sumw2(); ptfound->Sumw2();
  etagood->Sumw2(); etafound->Sumw2();
  pteff->Divide(ptfound,ptgood,1,1,"b");
  etaeff->Divide(etafound,etagood,1,1,"b");
  TFile *file = TFile::Open("eff.root","RECREATE");
  ptgood->Write();
  ptfound->Write();
  pteff->Write();
  etafound->Write();
  etagood->Write();
  etaeff->Write();
  file->Close();
  
  delete [] particles;
  delete [] ftracks;
  */  
}

#endif
