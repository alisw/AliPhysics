// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliHLTStandardIncludes.h"

#ifdef use_root
#include <TH1.h>
#include <TFile.h>
#endif

#include "AliHLTLogging.h"
#include "AliHLTHoughEval.h"
#include "AliHLTMemHandler.h"
#include "AliHLTTrackArray.h"
#include "AliHLTHoughBaseTransformer.h"
#include "AliHLTDigitData.h"
#include "AliHLTHoughTrack.h"
#include "AliHLTTransform.h"
#include "AliHLTHistogram.h"
#include "AliHLTHistogram1D.h"

#if __GNUC__ == 3
using namespace std;
#endif

/** /class AliHLTHoughEval
//<pre>
//_____________________________________________________________
// AliHLTHoughEval
//
// Evaluation class for tracklets produced by the Hough transform.
//
</pre>
*/

ClassImp(AliHLTHoughEval)

AliHLTHoughEval::AliHLTHoughEval()
{
  //default ctor  
  fRemoveFoundTracks = kFALSE;
  fNumOfPadsToLook = 1;
  fNumOfRowsToMiss = 1;
  fEtaHistos=0;
  fRowPointers = 0;
}


AliHLTHoughEval::~AliHLTHoughEval()
{
  //dtor
  fHoughTransformer = 0;
  if(fRowPointers)
    {
      for(Int_t i=0; i<fNrows; i++)
	fRowPointers[i] = 0;
      delete [] fRowPointers;
    }
}

void AliHLTHoughEval::InitTransformer(AliHLTHoughBaseTransformer *transformer)
{
  //Init hough transformer
  fHoughTransformer = transformer;
  fSlice = fHoughTransformer->GetSlice();
  fPatch = fHoughTransformer->GetPatch();
  fNrows = AliHLTTransform::GetLastRow(fPatch) - AliHLTTransform::GetFirstRow(fPatch) + 1;
  fNEtaSegments = fHoughTransformer->GetNEtaSegments();
  fEtaMin = fHoughTransformer->GetEtaMin();
  fEtaMax = fHoughTransformer->GetEtaMax();
  fZVertex = fHoughTransformer->GetZVertex();
  GenerateLUT();
}

void AliHLTHoughEval::GenerateLUT()
{
  //Generate a Look-up table, to limit the access to raw data
  
  if(!fRowPointers)
    fRowPointers = new AliHLTDigitRowData*[fNrows];

  AliHLTDigitRowData *tempPt = (AliHLTDigitRowData*)fHoughTransformer->GetDataPointer();
  if(!tempPt)
    printf("\nAliHLTHoughEval::GenerateLUT : Zero data pointer\n");
  
  for(Int_t i=AliHLTTransform::GetFirstRow(fPatch); i<=AliHLTTransform::GetLastRow(fPatch); i++)
    {
      Int_t prow = i - AliHLTTransform::GetFirstRow(fPatch);
      fRowPointers[prow] = tempPt;
      AliHLTMemHandler::UpdateRowPointer(tempPt);
    }
  
}

Bool_t AliHLTHoughEval::LookInsideRoad(AliHLTHoughTrack *track,Int_t &nrowscrossed,Int_t *rowrange,Bool_t remove)
{
  //Look at rawdata along the road specified by the track candidates.
  //If track is good, return true, if not return false.
  
  Int_t sector,row;
  
  Int_t nrow=0,npixs=0;//,rows_crossed=0;
  Float_t xyz[3];
  
  Int_t totalcharge=0;//total charge along the road
  
  //for(Int_t padrow = AliHLTTransform::GetFirstRow(fPatch); padrow <= AliHLTTransform::GetLastRow(fPatch); padrow++)
  for(Int_t padrow = rowrange[0]; padrow<=rowrange[1]; padrow++)
    {
      Int_t prow = padrow - AliHLTTransform::GetFirstRow(fPatch);
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
	  xyz[0] += AliHLTTransform::Row2X(track->GetFirstRow());
	  Float_t r = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]);
	  xyz[2] = r*track->GetTgl();
	}
      
      AliHLTTransform::Slice2Sector(fSlice,padrow,sector,row);
      AliHLTTransform::Local2Raw(xyz,sector,row);

      npixs=0;
      
      //Get the timebins for this pad
      AliHLTDigitRowData *tempPt = fRowPointers[prow];
      if(!tempPt) 
	{
	  printf("AliHLTHoughEval::LookInsideRoad : Zero data pointer\n");
	  continue;
	}
      
      //Look at both sides of the pad:
      for(Int_t p=(Int_t)rint(xyz[1])-fNumOfPadsToLook; p<=(Int_t)rint(xyz[1])+fNumOfPadsToLook; p++)
	{
	  AliHLTDigitData *digPt = tempPt->fDigitData;
	  for(UInt_t j=0; j<tempPt->fNDigit; j++)
	    {
	      Int_t pad = digPt[j].fPad;
	      Int_t charge = digPt[j].fCharge;
	      if(charge <= fHoughTransformer->GetLowerThreshold()) continue;
	      if(pad < p) continue;
	      if(pad > p) break;
	      UShort_t time = digPt[j].fTime;
	      Double_t eta = AliHLTTransform::GetEta(fSlice,padrow,pad,time);
	      Int_t pixelindex = fHoughTransformer->GetEtaIndex(eta);
	      if(pixelindex != track->GetEtaIndex()) continue;
	      totalcharge += digPt[j].fCharge;
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
  
  nrowscrossed += nrow; //Update the number of rows crossed.
  
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

void AliHLTHoughEval::FindEta(AliHLTTrackArray *tracks)
{
  //Find the corresponding eta slice hough space  
  Int_t sector,row;
  Float_t xyz[3];
  
  Int_t ntracks = tracks->GetNTracks();
  fEtaHistos = new AliHLTHistogram1D*[ntracks];
  
  Char_t hname[100];
  for(Int_t i=0; i<ntracks; i++)
    {
      sprintf(hname,"etahist_%d",i);
      fEtaHistos[i] = new AliHLTHistogram1D(hname,hname,100,0,1);
    }
  Double_t etaslice = (fEtaMax - fEtaMin)/fNEtaSegments;
  
  for(Int_t ntr=0; ntr<ntracks; ntr++)
    {
      AliHLTHoughTrack *track = (AliHLTHoughTrack*)tracks->GetCheckedTrack(ntr);
      if(!track) continue;
      for(Int_t padrow = AliHLTTransform::GetFirstRow(fPatch); padrow <= AliHLTTransform::GetLastRow(fPatch); padrow++)
	{
	  Int_t prow = padrow - AliHLTTransform::GetFirstRow(fPatch);
	  
	  if(!track->GetCrossingPoint(padrow,xyz))  
	    {
	      printf("AliHLTHoughEval::LookInsideRoad : Track does not cross line!!\n");
	      continue;
	    }
	  
	  AliHLTTransform::Slice2Sector(fSlice,padrow,sector,row);
	  AliHLTTransform::Local2Raw(xyz,sector,row);
	  
	  //Get the timebins for this pad
	  AliHLTDigitRowData *tempPt = fRowPointers[prow];
	  if(!tempPt) 
	    {
	      printf("AliHLTHoughEval::LookInsideRoad : Zero data pointer\n");
	      continue;
	    }
	  
	  //Look at both sides of the pad:
	  for(Int_t p=(Int_t)rint(xyz[1])-fNumOfPadsToLook; p<=(Int_t)rint(xyz[1])+fNumOfPadsToLook; p++)
	    {
	      AliHLTDigitData *digPt = tempPt->fDigitData;
	      for(UInt_t j=0; j<tempPt->fNDigit; j++)
		{
		  UChar_t pad = digPt[j].fPad;
		  Int_t charge = digPt[j].fCharge;
		  if(charge <= fHoughTransformer->GetLowerThreshold()) continue;
		  if(pad < p) continue;
		  if(pad > p) break;
		  UShort_t time = digPt[j].fTime;
		  Double_t eta = AliHLTTransform::GetEta(fSlice,padrow,pad,time);
		  Int_t pixelindex = (Int_t)(eta/etaslice);
		  if(pixelindex > track->GetEtaIndex()+1) continue;
		  if(pixelindex < track->GetEtaIndex()-1) break;
		  fEtaHistos[ntr]->Fill(eta,digPt[j].fCharge);
		}
	    }
	}
    }
  
  for(Int_t i=0; i<ntracks; i++)
    {
      AliHLTHistogram1D *hist = fEtaHistos[i];
      Int_t maxbin = hist->GetMaximumBin();
      Double_t maxvalue = hist->GetBinContent(maxbin);
      AliHLTHoughTrack *track = (AliHLTHoughTrack*)tracks->GetCheckedTrack(i);
      if(!track) continue;
      if(hist->GetBinContent(maxbin-1)<maxvalue && hist->GetBinContent(maxbin+1)<maxvalue)
	{
	  track->SetWeight((Int_t)maxvalue,kTRUE); 
	  track->SetEta(hist->GetBinCenter(maxbin));
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

void AliHLTHoughEval::DisplayEtaSlice(Int_t etaindex,AliHLTHistogram *hist)
{
  //Display the current raw data inside the (slice,patch)

  if(!hist)
    {
      printf("AliHLTHoughEval::DisplayEtaSlice : No input histogram!\n");
      return;
    }
  
  for(Int_t padrow = AliHLTTransform::GetFirstRow(fPatch); padrow <= AliHLTTransform::GetLastRow(fPatch); padrow++)
    {
      Int_t prow = padrow - AliHLTTransform::GetFirstRow(fPatch);
                  
      AliHLTDigitRowData *tempPt = fRowPointers[prow];
      if(!tempPt) 
	{
	  printf("AliHLTHoughEval::DisplayEtaSlice : Zero data pointer\n");
	  continue;
	}
      
      AliHLTDigitData *digPt = tempPt->fDigitData;
      if((Int_t)tempPt->fRow != padrow)
	{
	  printf("\nAliHLTHoughEval::DisplayEtaSlice : Mismatching padrows!!!\n");
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
	  AliHLTTransform::Slice2Sector(fSlice,padrow,sector,row);
	  AliHLTTransform::Raw2Local(xyz,sector,row,pad,time);
	  xyz[2] -= fZVertex;
	  Double_t eta = AliHLTTransform::GetEta(xyz);
	  Int_t pixelindex = fHoughTransformer->GetEtaIndex(eta);//(Int_t)(eta/etaslice);
	  if(pixelindex != etaindex) continue;
	  hist->Fill(xyz[0],xyz[1],charge);
	}
    }
  
}

#ifdef use_root
void AliHLTHoughEval::CompareMC(AliHLTTrackArray */*tracks*/,Char_t */*trackfile*/,Int_t /*threshold*/)
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
	  LOG(AliHLTLog::kError,"AliHLTHoughEval::CompareMC","Input file")
	    <<"Error in file reading"<<ENDLOG;
	  return;
	}
    }
  else
    {
      LOG(AliHLTLog::kError,"AliHLTHoughEval::CompareMC","Input")
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
      AliHLTHoughTrack *tr = (AliHLTHoughTrack*)tracks->GetCheckedTrack(i);
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
