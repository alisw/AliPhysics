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

AliL3TrackArray *AliL3HoughMaxFinder::FindBigMaxima(TH2F *hist)
{
  Int_t xmin = hist->GetXaxis()->GetFirst();
  Int_t xmax = hist->GetXaxis()->GetLast();
  Int_t ymin = hist->GetYaxis()->GetFirst();
  Int_t ymax = hist->GetYaxis()->GetLast();
  Int_t bin[25],bin_index;
  Stat_t value[25];
  
  AliL3TrackArray *tracks = new AliL3TrackArray("AliL3HoughTrack");
  AliL3HoughTrack *track;
  
  for(Int_t xbin=xmin+2; xbin<xmax-3; xbin++)
    {
      for(Int_t ybin=ymin+2; ybin<ymax-3; ybin++)
	{
	  bin_index=0;
	  for(Int_t xb=xbin-2; xb<xbin+3; xb++)
	    {
	      for(Int_t yb=ybin-2; yb<ybin+3; yb++)
		{
		  bin[bin_index]=hist->GetBin(xb,yb);
		  value[bin_index]=hist->GetBinContent(bin[bin_index]);
		  bin_index++;
		}
	      
	    }
	  if(value[12]==0) continue;
	  Int_t b=0;
	  while(1)
	    {
	      if(value[b] > value[12] || b==bin_index) break;
	      b++;
	      printf("b %d\n",b);
	    }
	  if(b == bin_index)
	    {
	      //Found maxima
	      Double_t max_x = hist->GetXaxis()->GetBinCenter(xbin);
	      Double_t max_y = hist->GetYaxis()->GetBinCenter(ybin);
	      track = (AliL3HoughTrack*)tracks->NextTrack();
	      track->SetTrackParameters(max_x,max_y,value[12]);
	    }
	}
    }
  
  tracks->QSort();
  return tracks;
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


AliL3TrackArray *AliL3HoughMaxFinder::LookForPeaks(TH2F *hist,Int_t nbins)
{
  
  AliL3TrackArray *tracks = new AliL3TrackArray("AliL3HoughTrack");
  AliL3HoughTrack *track;
  Int_t xmin = hist->GetXaxis()->GetFirst();
  Int_t xmax = hist->GetXaxis()->GetLast();
  Int_t ymin = hist->GetYaxis()->GetFirst();
  Int_t ymax = hist->GetYaxis()->GetLast();
  
  Int_t weight_loc;
  Int_t candidate=0;
  for(Int_t xbin=xmin+nbins; xbin <= xmax - nbins; xbin++)
    {
      for(Int_t ybin=ymin+nbins; ybin <= ymax - nbins; ybin++)
	{
	  weight_loc=0;
	  for(Int_t xbin_loc = xbin-nbins; xbin_loc <= xbin+nbins; xbin_loc++)
	    {
	      for(Int_t ybin_loc = ybin-nbins; ybin_loc <= ybin+nbins; ybin_loc++)
		{
		  Int_t bin_loc = hist->GetBin(xbin_loc,ybin_loc);
		  weight_loc += (Int_t)hist->GetBinContent(bin_loc);
		}
	    }
	  
	  if(weight_loc > 0)
	    {
	      track = (AliL3HoughTrack*)tracks->NextTrack();
	      track->SetTrackParameters(hist->GetXaxis()->GetBinCenter(xbin),hist->GetYaxis()->GetBinCenter(ybin),weight_loc);
	    }
	  
	}
    }
  tracks->QSort();
  
  AliL3HoughTrack *track1,*track2;

  for(Int_t i=1; i<tracks->GetNTracks(); i++)
    {
      track1 = (AliL3HoughTrack*)tracks->GetCheckedTrack(i);
      if(!track1) continue;
      Int_t xbin1 = hist->GetXaxis()->FindBin(track1->GetKappa());
      Int_t ybin1 = hist->GetYaxis()->FindBin(track1->GetPhi0());
      for(Int_t j=0; j < i; j++)
	{
	  track2 = (AliL3HoughTrack*)tracks->GetCheckedTrack(j);
	  if(!track2) continue;
	  Int_t xbin2 = hist->GetXaxis()->FindBin(track2->GetKappa());
	  Int_t ybin2 = hist->GetYaxis()->FindBin(track2->GetPhi0());
	  if(abs(xbin1-xbin2) < 10 && abs(ybin1-ybin2) < 10)
	    {
	      tracks->Remove(i);
	      break;
	    }
	}

    }
  tracks->Compress();
  return tracks;
}

AliL3TrackArray *AliL3HoughMaxFinder::LookInWindows(TH2F *hist,Int_t nbins,Int_t t1,Double_t t2,Int_t t3)
{

  AliL3TrackArray *tracks = new AliL3TrackArray("AliL3HoughTrack");
  AliL3HoughTrack *track = (AliL3HoughTrack*)tracks->NextTrack();

  Int_t xmin = hist->GetXaxis()->GetFirst();
  Int_t xmax = hist->GetXaxis()->GetLast();
  Int_t ymin = hist->GetYaxis()->GetFirst();
  Int_t ymax = hist->GetYaxis()->GetLast();
  
  for(Int_t xbin=xmin+nbins; xbin < xmax - nbins; xbin += nbins)
    {
      for(Int_t ybin=ymin+nbins; ybin<ymax-nbins; ybin += nbins)
	{
	  //Int_t bin = hist->GetBin(xbin,ybin);
	  //if((Int_t)hist->GetBinContent(bin)==0) continue;

	  Int_t xrange[2] = {xbin-nbins,xbin+nbins};
	  Int_t yrange[2] = {ybin-nbins,ybin+nbins};
	  if(LocatePeak(hist,track,xrange,yrange,t1,t2,t3))
	    track = (AliL3HoughTrack*)tracks->NextTrack();
	  else
	    continue;
	}
    }
  tracks->QSort();
  return tracks;
  
}

Bool_t AliL3HoughMaxFinder::LocatePeak(TH2F *hist,AliL3HoughTrack *track,Int_t *xrange,Int_t *yrange,Int_t t1,Double_t t2,Int_t t3)
{

  /*  Int_t xmin = hist->GetXaxis()->FindBin(track->GetKappa()) - nbins;
  Int_t xmax = hist->GetXaxis()->FindBin(track->GetKappa()) + nbins;
  Int_t ymin = hist->GetYaxis()->FindBin(track->GetPhi0()) - nbins;
  Int_t ymax = hist->GetYaxis()->FindBin(track->GetPhi0()) + nbins;
  Int_t nbinsx = nbins*2 + 1;

  if(xmin < 0 || xmax > hist->GetXaxis()->GetNbins() || ymin < 0 || ymax > hist->GetYaxis()->GetNbins())
    return false;
  */
  Int_t xmin = xrange[0];
  Int_t xmax = xrange[1];
  Int_t ymin = yrange[0];
  Int_t ymax = yrange[1];
  
  if(xmin < 0)
    xmin=0;
  if(xmax > hist->GetXaxis()->GetNbins())
    xmax = hist->GetXaxis()->GetNbins();
  if(ymin < 0)
    ymin=0;
  if(ymax > hist->GetYaxis()->GetNbins())
    ymax = hist->GetYaxis()->GetNbins();

  printf("Defined xrange %d %d yrange %d %d\n",xmin,xmax,ymin,ymax);


  Int_t totbins = hist->GetXaxis()->GetNbins();
  Int_t *m = new Int_t[totbins];
  Int_t *m_low = new Int_t[totbins];
  Int_t *m_up = new Int_t[totbins];
  
  
  for(Int_t i=0; i<totbins; i++)
    {
      m[i]=0;
      m_low[i]=0;
      m_up[i]=0;
    }

  Int_t max_x=0,sum=0,max_xbin=0,bin;

  for(Int_t xbin=xmin; xbin<=xmax; xbin++)
    {
      for(Int_t ybin=ymin; ybin < ymax - t1; ybin++)
	{
	  sum = 0;
	  for(Int_t y=ybin; y < ybin+t1; y++)
	    {
	      //Inside window
	      bin = hist->GetBin(xbin,y);
	      sum += (Int_t)hist->GetBinContent(bin);
	      
	    }
	  if(sum > m[xbin]) //Max value locally in this xbin
	    {
	      m[xbin]=sum;
	      m_low[xbin]=ybin;
	      m_up[xbin]=ybin + t1 - 1;
	    }
	  
	}
      
      if(m[xbin] > max_x) //Max value globally in x-direction
	{
	  max_xbin = xbin;
	  max_x = m[xbin];//sum;
	}
    }
  
  if(max_x == 0) return false;
  
  //printf("max_xbin %d max_x %d m_low %d m_up %d\n",max_xbin,max_x,m_low[max_xbin],m_up[max_xbin]);
  //printf("ylow %f yup %f\n",hist->GetYaxis()->GetBinCenter(m_low[max_xbin]),hist->GetYaxis()->GetBinCenter(m_up[max_xbin]));

  //Determine a width in the x-direction
  Int_t x_low=0,x_up=0;
  for(Int_t xbin=max_xbin-1; xbin >= xmin; xbin--)
    {
      if(m[xbin]*t2 < max_x)
	{
	  x_low = xbin+1;
	  break;
	}
    }
  for(Int_t xbin = max_xbin+1; xbin <=xmax; xbin++)
    {
      if(m[xbin]*t2 < max_x)
	{
	  x_up = xbin-1;
	  break;
	}
    }
  
  Double_t top=0,butt=0,value,x_peak;
  
  // printf("xlow %f xup %f\n",hist->GetXaxis()->GetBinCenter(x_low),hist->GetXaxis()->GetBinCenter(x_up));
  //printf("Spread in x %d\n",x_up-x_low +1);

  //Now, calculate the center of mass in x-direction
  for(Int_t xbin=x_low; xbin <= x_up; xbin++)
    {
      value = hist->GetXaxis()->GetBinCenter(xbin);
      top += value*m[xbin];
      butt += m[xbin];
    }
  x_peak = top/butt;
  
  //printf("FOund x_peak %f\n",x_peak);

  //Find the peak in y direction:
  Int_t x_l = hist->GetXaxis()->FindBin(x_peak);
  if(hist->GetXaxis()->GetBinCenter(x_l) > x_peak)
    x_l--;

  Int_t x_u = x_l + 1;
  
  if(hist->GetXaxis()->GetBinCenter(x_l) > x_peak || hist->GetXaxis()->GetBinCenter(x_u) <= x_peak)
    printf("\nAliL3HoughMaxFinder::FindPeak : Wrong xrange %f %f %f\n\n",hist->GetXaxis()->GetBinCenter(x_l),x_peak,hist->GetXaxis()->GetBinCenter(x_u));
  
  //printf("\nxlow %f xup %f\n",hist->GetXaxis()->GetBinCenter(x_l),hist->GetXaxis()->GetBinCenter(x_u));

  value=top=butt=0;
  
  //printf("ylow %f yup %f\n",hist->GetYaxis()->GetBinCenter(m_low[x_l]),hist->GetYaxis()->GetBinCenter(m_up[x_l]));
  //printf("ylow %f yup %f\n",hist->GetYaxis()->GetBinCenter(m_low[x_u]),hist->GetYaxis()->GetBinCenter(m_up[x_u]));
  
  for(Int_t ybin=m_low[x_l]; ybin <= m_up[x_l]; ybin++)
    {
      value = hist->GetYaxis()->GetBinCenter(ybin);
      bin = hist->GetBin(x_l,ybin);
      top += value*hist->GetBinContent(bin);
      butt += hist->GetBinContent(bin);
    }
  Double_t y_peak_low = top/butt;
  
  //printf("y_peak_low %f\n",y_peak_low);

  value=top=butt=0;
  for(Int_t ybin=m_low[x_u]; ybin <= m_up[x_u]; ybin++)
    {
      value = hist->GetYaxis()->GetBinCenter(ybin);
      bin = hist->GetBin(x_u,ybin);
      top += value*hist->GetBinContent(bin);
      butt += hist->GetBinContent(bin);
    }
  Double_t y_peak_up = top/butt;
  
  //printf("y_peak_up %f\n",y_peak_up);

  Double_t x_value_up = hist->GetXaxis()->GetBinCenter(x_u);
  Double_t x_value_low = hist->GetXaxis()->GetBinCenter(x_l);

  Double_t y_peak = (y_peak_low*(x_value_up - x_peak) + y_peak_up*(x_peak - x_value_low))/(x_value_up - x_value_low);

  //Find the weight:
  bin = hist->FindBin(x_peak,y_peak);
  Int_t weight = (Int_t)hist->GetBinContent(bin);

  if(weight==0) return false;
  
  track->SetTrackParameters(x_peak,y_peak,weight);
  
  delete [] m;
  delete [] m_up;
  delete [] m_low;
  return true;
}


AliL3TrackArray *AliL3HoughMaxFinder::FindPeak(TH2F *hist,Int_t t1,Double_t t2,Int_t t3)
{
  //Attempt of a more sophisticated peak finder.
  
  AliL3TrackArray *tracks = new AliL3TrackArray("AliL3HoughTrack",1);

  Int_t xmin = hist->GetXaxis()->GetFirst();
  Int_t xmax = hist->GetXaxis()->GetLast();
  Int_t ymin = hist->GetYaxis()->GetFirst();
  Int_t ymax = hist->GetYaxis()->GetLast();
  Int_t nbinsx = hist->GetXaxis()->GetNbins()+1;
  
  Int_t *m = new Int_t[nbinsx];
  Int_t *m_low = new Int_t[nbinsx];
  Int_t *m_up = new Int_t[nbinsx];
  
  
 recompute:  //this is a goto.
  
  for(Int_t i=0; i<nbinsx; i++)
    {
      m[i]=0;
      m_low[i]=0;
      m_up[i]=0;
    }

  Int_t max_x=0,sum=0,max_xbin=0,bin;

  for(Int_t xbin=xmin; xbin<=xmax; xbin++)
    {
      for(Int_t ybin=ymin; ybin < ymax - t1; ybin++)
	{
	  sum = 0;
	  for(Int_t y=ybin; y < ybin+t1; y++)
	    {
	      //Inside window
	      bin = hist->GetBin(xbin,y);
	      sum += (Int_t)hist->GetBinContent(bin);
	      
	    }
	  if(sum > m[xbin]) //Max value locally in this xbin
	    {
	      m[xbin]=sum;
	      m_low[xbin]=ybin;
	      m_up[xbin]=ybin + t1 - 1;
	    }
	  
	}
      
      if(m[xbin] > max_x) //Max value globally in x-direction
	{
	  max_xbin = xbin;
	  max_x = m[xbin];//sum;
	}
    }
  //printf("max_xbin %d max_x %d m_low %d m_up %d\n",max_xbin,max_x,m_low[max_xbin],m_up[max_xbin]);
  //printf("ylow %f yup %f\n",hist->GetYaxis()->GetBinCenter(m_low[max_xbin]),hist->GetYaxis()->GetBinCenter(m_up[max_xbin]));

  //Determine a width in the x-direction
  Int_t x_low=0,x_up=0;
  for(Int_t xbin=max_xbin-1; xbin >= xmin; xbin--)
    {
      if(m[xbin]*t2 < max_x)
	{
	  x_low = xbin+1;
	  break;
	}
    }
  for(Int_t xbin = max_xbin+1; xbin <=xmax; xbin++)
    {
      if(m[xbin]*t2 < max_x)
	{
	  x_up = xbin-1;
	  break;
	}
    }
  
  Double_t top=0,butt=0,value,x_peak;
  if(x_up - x_low + 1 > t3)
    {
      t1 -= 1;
      printf("\nxrange out if limit x_up %d x_low %d\n\n",x_low,x_up);
      if(t1 > 1)
	goto recompute;
      else
	{
	  x_peak = hist->GetXaxis()->GetBinCenter(max_xbin);
	  goto moveon;
	}
    }
  
  //printf("xlow %f xup %f\n",hist->GetXaxis()->GetBinCenter(x_low),hist->GetXaxis()->GetBinCenter(x_up));
  //printf("Spread in x %d\n",x_up-x_low +1);

  //Now, calculate the center of mass in x-direction
  for(Int_t xbin=x_low; xbin <= x_up; xbin++)
    {
      value = hist->GetXaxis()->GetBinCenter(xbin);
      top += value*m[xbin];
      butt += m[xbin];
    }
  x_peak = top/butt;
  
 moveon:
  
  //Find the peak in y direction:
  Int_t x_l = hist->GetXaxis()->FindBin(x_peak);
  if(hist->GetXaxis()->GetBinCenter(x_l) > x_peak)
    x_l--;

  Int_t x_u = x_l + 1;
  
  if(hist->GetXaxis()->GetBinCenter(x_l) > x_peak || hist->GetXaxis()->GetBinCenter(x_u) <= x_peak)
    printf("\nAliL3HoughMaxFinder::FindPeak : Wrong xrange %f %f %f\n\n",hist->GetXaxis()->GetBinCenter(x_l),x_peak,hist->GetXaxis()->GetBinCenter(x_u));
    
    //printf("\nxlow %f xup %f\n",hist->GetXaxis()->GetBinCenter(x_l),hist->GetXaxis()->GetBinCenter(x_u));

  value=top=butt=0;
  
  //printf("ylow %f yup %f\n",hist->GetYaxis()->GetBinCenter(m_low[x_l]),hist->GetYaxis()->GetBinCenter(m_up[x_l]));
  //printf("ylow %f yup %f\n",hist->GetYaxis()->GetBinCenter(m_low[x_u]),hist->GetYaxis()->GetBinCenter(m_up[x_u]));
  
  for(Int_t ybin=m_low[x_l]; ybin <= m_up[x_l]; ybin++)
    {
      value = hist->GetYaxis()->GetBinCenter(ybin);
      bin = hist->GetBin(x_l,ybin);
      top += value*hist->GetBinContent(bin);
      butt += hist->GetBinContent(bin);
    }
  Double_t y_peak_low = top/butt;
  
  //printf("y_peak_low %f\n",y_peak_low);

  value=top=butt=0;
  for(Int_t ybin=m_low[x_u]; ybin <= m_up[x_u]; ybin++)
    {
      value = hist->GetYaxis()->GetBinCenter(ybin);
      bin = hist->GetBin(x_u,ybin);
      top += value*hist->GetBinContent(bin);
      butt += hist->GetBinContent(bin);
    }
  Double_t y_peak_up = top/butt;
  
  //printf("y_peak_up %f\n",y_peak_up);

  Double_t x_value_up = hist->GetXaxis()->GetBinCenter(x_u);
  Double_t x_value_low = hist->GetXaxis()->GetBinCenter(x_l);

  Double_t y_peak = (y_peak_low*(x_value_up - x_peak) + y_peak_up*(x_peak - x_value_low))/(x_value_up - x_value_low);


  //Find the weight:
  bin = hist->FindBin(x_peak,y_peak);
  Int_t weight = (Int_t)hist->GetBinContent(bin);

  AliL3HoughTrack *track = (AliL3HoughTrack*)tracks->NextTrack();
  track->SetTrackParameters(x_peak,y_peak,weight);
  
  
  delete [] m;
  delete [] m_low;
  delete [] m_up;
  
  return tracks;

    
}

