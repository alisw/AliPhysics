#include <string.h>

#include <TH2.h>

#include "AliL3TrackArray.h"
#include "AliL3HoughTrack.h"
#include "AliL3HoughMaxFinder.h"

ClassImp(AliL3HoughMaxFinder)

  
AliL3HoughMaxFinder::AliL3HoughMaxFinder()
{
  //Default constructor
  fThreshold = 0;
  //fTracks = 0;
  fHistoType=0;
}


AliL3HoughMaxFinder::AliL3HoughMaxFinder(Char_t *histotype)
{
  //Constructor

  //fTracks = new AliL3TrackArray("AliL3HoughTrack");
  if(strcmp(histotype,"KappaPhi")==0) fHistoType='c';
  if(strcmp(histotype,"DPsi")==0) fHistoType='l';
  
}


AliL3HoughMaxFinder::~AliL3HoughMaxFinder()
{
  //Destructor

}


AliL3TrackArray *AliL3HoughMaxFinder::FindMaxima(TH2F *hist,Int_t *rowrange,Int_t ref_row)
{
  //Locate all the maxima in input histogram.
  //Maxima is defined as bins with more entries than the
  //immediately neighbouring bins. 
  
  Int_t xmin = hist->GetXaxis()->GetFirst();
  Int_t xmax = hist->GetXaxis()->GetLast();
  Int_t ymin = hist->GetYaxis()->GetFirst();
  Int_t ymax = hist->GetYaxis()->GetLast();
  Int_t bin[9],track_counter=0;
  Stat_t value[9];
  
  AliL3TrackArray *tracks = new AliL3TrackArray("AliL3HoughTrack");
  AliL3HoughTrack *track;

  for(Int_t xbin=xmin+1; xbin<xmax-1; xbin++)
    {
      for(Int_t ybin=ymin+1; ybin<ymax-1; ybin++)
	{
	  bin[0] = hist->GetBin(xbin-1,ybin-1);
	  bin[1] = hist->GetBin(xbin,ybin-1);
	  bin[2] = hist->GetBin(xbin+1,ybin-1);
	  bin[3] = hist->GetBin(xbin-1,ybin);
	  bin[4] = hist->GetBin(xbin,ybin);
	  bin[5] = hist->GetBin(xbin+1,ybin);
	  bin[6] = hist->GetBin(xbin-1,ybin+1);
	  bin[7] = hist->GetBin(xbin,ybin+1);
	  bin[8] = hist->GetBin(xbin+1,ybin+1);
	  value[0] = hist->GetBinContent(bin[0]);
	  value[1] = hist->GetBinContent(bin[1]);
	  value[2] = hist->GetBinContent(bin[2]);
	  value[3] = hist->GetBinContent(bin[3]);
	  value[4] = hist->GetBinContent(bin[4]);
	  value[5] = hist->GetBinContent(bin[5]);
	  value[6] = hist->GetBinContent(bin[6]);
	  value[7] = hist->GetBinContent(bin[7]);
	  value[8] = hist->GetBinContent(bin[8]);
	  
	  if(value[4] <= fThreshold) continue;//central bin below threshold
	  
	  if(value[4]>value[0] && value[4]>value[1] && value[4]>value[2]
	     && value[4]>value[3] && value[4]>value[5] && value[4]>value[6]
	     && value[4]>value[7] && value[4]>value[8])
	    {
	      //Found a local maxima
	      Float_t max_x = hist->GetXaxis()->GetBinCenter(xbin);
	      Float_t max_y = hist->GetYaxis()->GetBinCenter(ybin);
	      
	      track = (AliL3HoughTrack*)tracks->NextTrack();
	      
	      
	      switch(fHistoType)
		{
		case 'c':
		  track->SetTrackParameters(max_x,max_y,(Int_t)value[4]);
		  break;
		case 'l':
		  track->SetLineParameters(max_x,max_y,(Int_t)value[4],rowrange,ref_row);
		  break;
		default:
		  printf("AliL3HoughMaxFinder: Error in tracktype\n");
		}
	      
	      track_counter++;
	      

	    }
	  else
	    continue; //not a maxima
	}
    }
  tracks->QSort();
  return tracks;
}
