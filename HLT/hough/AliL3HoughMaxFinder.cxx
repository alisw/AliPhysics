//Author:        Anders Strand Vestbo
//Last Modified: 28.6.01

#include <string.h>

#include "AliL3Histogram.h"
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


AliL3HoughMaxFinder::AliL3HoughMaxFinder(Char_t *histotype,AliL3Histogram *hist)
{
  //Constructor

  //fTracks = new AliL3TrackArray("AliL3HoughTrack");
  if(strcmp(histotype,"KappaPhi")==0) fHistoType='c';
  if(strcmp(histotype,"DPsi")==0) fHistoType='l';
  
  if(hist)
    fCurrentHisto = hist;
}


AliL3HoughMaxFinder::~AliL3HoughMaxFinder()
{
  //Destructor

}

Int_t *AliL3HoughMaxFinder::FindAbsMaxima()
{
  if(!fCurrentHisto)
    {
      printf("AliL3HoughMaxFinder::FindAbsMaxima : No histogram!\n");
      return 0;
    }
  AliL3Histogram *hist = fCurrentHisto;

  Int_t xmin = hist->GetFirstXbin();
  Int_t xmax = hist->GetLastXbin();
  Int_t ymin = hist->GetFirstYbin();
  Int_t ymax = hist->GetLastYbin();  
  Int_t bin;
  Int_t *max_bin = new Int_t[2];
  Stat_t value,max_value=0;

  for(Int_t xbin=xmin; xbin<=xmax; xbin++)
    {
      for(Int_t ybin=ymin; ybin<=ymax; ybin++)
	{
	  bin = hist->GetBin(xbin,ybin);
	  value = hist->GetBinContent(bin);
	  if(value>max_value)
	    {
	      max_value = value;
	      max_bin[0] = xbin;
	      max_bin[1] = ybin;
	    }
	}
    }
  return max_bin;
}

AliL3TrackArray *AliL3HoughMaxFinder::FindBigMaxima(AliL3Histogram *hist)
{
  
  Int_t xmin = hist->GetFirstXbin();
  Int_t xmax = hist->GetLastXbin();
  Int_t ymin = hist->GetFirstYbin();
  Int_t ymax = hist->GetLastYbin();
  Int_t bin[25],bin_index;
  Stat_t value[25];
  
  AliL3TrackArray *tracks = new AliL3TrackArray("AliL3HoughTrack");
  AliL3HoughTrack *track;
  
  for(Int_t xbin=xmin+2; xbin<xmax-3; xbin++)
    {
      for(Int_t ybin=ymin+2; ybin<ymax-3; ybin++)
	{
	  bin_index=0;
	  for(Int_t xb=xbin-2; xb<=xbin+2; xb++)
	    {
	      for(Int_t yb=ybin-2; yb<=ybin+2; yb++)
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
	      Double_t max_x = hist->GetBinCenterX(xbin);
	      Double_t max_y = hist->GetBinCenterY(ybin);
	      track = (AliL3HoughTrack*)tracks->NextTrack();
	      track->SetTrackParameters(max_x,max_y,value[12]);
	    }
	}
    }
  
  tracks->QSort();
  return tracks;
}


AliL3TrackArray *AliL3HoughMaxFinder::FindMaxima(AliL3Histogram *hist,Int_t *rowrange,Int_t ref_row)
{
  //Locate all the maxima in input histogram.
  //Maxima is defined as bins with more entries than the
  //immediately neighbouring bins. 
  
  Int_t xmin = hist->GetFirstXbin();
  Int_t xmax = hist->GetLastXbin();
  Int_t ymin = hist->GetFirstYbin();
  Int_t ymax = hist->GetLastYbin();
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
	      Float_t max_x = hist->GetBinCenterX(xbin);
	      Float_t max_y = hist->GetBinCenterY(ybin);
	      
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


AliL3TrackArray *AliL3HoughMaxFinder::LookForPeaks(AliL3Histogram *hist,Int_t nbins)
{
  
  AliL3TrackArray *tracks = new AliL3TrackArray("AliL3HoughTrack");
  AliL3HoughTrack *track;
  Int_t xmin = hist->GetFirstXbin();
  Int_t xmax = hist->GetLastXbin();
  Int_t ymin = hist->GetFirstYbin();
  Int_t ymax = hist->GetLastYbin();
  
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
	      track->SetTrackParameters(hist->GetBinCenterX(xbin),hist->GetBinCenterY(ybin),weight_loc);
	    }
	  
	}
    }
  tracks->QSort();
  
  AliL3HoughTrack *track1,*track2;

  for(Int_t i=1; i<tracks->GetNTracks(); i++)
    {
      track1 = (AliL3HoughTrack*)tracks->GetCheckedTrack(i);
      if(!track1) continue;
      Int_t xbin1 = hist->FindXbin(track1->GetKappa());
      Int_t ybin1 = hist->FindXbin(track1->GetPhi0());
      for(Int_t j=0; j < i; j++)
	{
	  track2 = (AliL3HoughTrack*)tracks->GetCheckedTrack(j);
	  if(!track2) continue;
	  Int_t xbin2 = hist->FindXbin(track2->GetKappa());
	  Int_t ybin2 = hist->FindYbin(track2->GetPhi0());
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

AliL3TrackArray *AliL3HoughMaxFinder::LookInWindows(AliL3Histogram *hist,Int_t nbins,Int_t t1,Double_t t2,Int_t t3)
{

  AliL3TrackArray *tracks = new AliL3TrackArray("AliL3HoughTrack");
  AliL3HoughTrack *track = (AliL3HoughTrack*)tracks->NextTrack();

  Int_t xmin = hist->GetFirstXbin();
  Int_t xmax = hist->GetLastXbin();
  Int_t ymin = hist->GetFirstYbin();
  Int_t ymax = hist->GetLastYbin();
  
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
  tracks->RemoveLast();
  return tracks;
  
}

AliL3HoughTrack *AliL3HoughMaxFinder::CalculatePeakInWindow(Int_t *maxbin,Int_t t0,Int_t t1,Double_t t2,Int_t t3)
{
  //Try to expand the area around the maxbin +- t0

  if(!fCurrentHisto)
    {
      printf("AliL3HoughMaxFinder::LocatePeak : No histogram\n");
      return 0;
    }
  AliL3Histogram *hist = fCurrentHisto;
  Int_t xmin = hist->GetFirstXbin();
  Int_t xmax = hist->GetLastXbin();
  Int_t ymin = hist->GetFirstYbin();
  Int_t ymax = hist->GetLastYbin();

  Int_t xlow = maxbin[0]-t0;
  if(xlow < xmin)
    xlow = xmin;
  Int_t xup = maxbin[0]+t0;
  if(xup > xmax)
    xup = xmax;
  Int_t ylow = maxbin[1]-t0;
  if(ylow < ymin)
    ylow = ymin;
  Int_t yup = maxbin[1]+t0;
  if(yup > ymax)
    yup = ymax;
  
  Int_t nbinsx = hist->GetNbinsX()+1;
  
  Int_t *m = new Int_t[nbinsx];
  Int_t *m_low = new Int_t[nbinsx];
  Int_t *m_up = new Int_t[nbinsx];
  
 recompute:

  Int_t max_x=0,sum=0,max_xbin=0,bin;
  
  for(Int_t xbin=xlow; xbin<=xup; xbin++)
    {
      for(Int_t ybin=ylow; ybin <= yup; ybin++)
	{
	  sum = 0;
	  for(Int_t y=ybin; y <= ybin+t1; y++)
	    {
	      if(y>yup) break; //reached the upper limit in y.
	      //Inside window
	      bin = hist->GetBin(xbin,y);
	      sum += (Int_t)hist->GetBinContent(bin);
	      
	    }
	  if(sum > m[xbin]) //Max value locally in this xbin
	    {
	      m[xbin]=sum;
	      m_low[xbin]=ybin;
	      m_up[xbin]=ybin + t1;
	    }
	  
	}
      
      if(m[xbin] > max_x) //Max value globally in x-direction
	{
	  max_xbin = xbin;
	  max_x = m[xbin];//sum;
	}
    }
  //printf("max_xbin %d max_x %d m_low %d m_up %d\n",max_xbin,max_x,m_low[max_xbin],m_up[max_xbin]);
  //printf("ylow %f yup %f\n",hist->GetBinCenterY(m_low[max_xbin]),hist->GetBinCenterY(m_up[max_xbin]));

  //Determine a width in the x-direction
  Int_t x_low=0,x_up=0;
  for(Int_t xbin=max_xbin-1; xbin >= xmin; xbin--)
    {
      if(m[xbin] < max_x*t2)
	{
	  x_low = xbin+1;
	  break;
	}
    }
  for(Int_t xbin = max_xbin+1; xbin <=xmax; xbin++)
    {
      if(m[xbin] < max_x*t2)
	{
	  x_up = xbin-1;
	  break;
	}
    }
  printf("x_low %d x_up %d\n",x_low,x_up);

  Double_t top=0,butt=0,value,x_peak;
  if(x_up - x_low + 1 > t3)
    {
      t1 -= 1;
      printf("\nxrange out if limit x_up %d x_low %d\n\n",x_low,x_up);
      if(t1 > 1)
	goto recompute;
      else
	{
	  x_peak = hist->GetBinCenterX(max_xbin);
	  goto moveon;
	}
    }
  
  //printf("xlow %f xup %f\n",hist->GetBinCenterX(x_low),hist->GetBinCenterX(x_up));
  //printf("Spread in x %d\n",x_up-x_low +1);

  //Now, calculate the center of mass in x-direction
  for(Int_t xbin=x_low; xbin <= x_up; xbin++)
    {
      value = hist->GetBinCenterX(xbin);
      top += value*m[xbin];
      butt += m[xbin];
    }
  x_peak = top/butt;
  
 moveon:
  
  //Find the peak in y direction:
  Int_t x_l = hist->FindXbin(x_peak);
  if(hist->GetBinCenterX(x_l) > x_peak)
    x_l--;

  Int_t x_u = x_l + 1;
  
  if(hist->GetBinCenterX(x_l) > x_peak || hist->GetBinCenterX(x_u) <= x_peak)
    printf("\nAliL3HoughMaxFinder::FindPeak : Wrong xrange %f %f %f\n\n",hist->GetBinCenterX(x_l),x_peak,hist->GetBinCenterX(x_u));
    
    //printf("\nxlow %f xup %f\n",hist->GetBinCenterX(x_l),hist->GetBinCenterX(x_u));

  value=top=butt=0;
  
  //printf("ylow %f yup %f\n",hist->GetBinCenterY(m_low[x_l]),hist->GetBinCenterY(m_up[x_l]));
  //printf("ylow %f yup %f\n",hist->GetBinCenterY(m_low[x_u]),hist->GetBinCenterY(m_up[x_u]));
  
  for(Int_t ybin=m_low[x_l]; ybin <= m_up[x_l]; ybin++)
    {
      value = hist->GetBinCenterY(ybin);
      bin = hist->GetBin(x_l,ybin);
      top += value*hist->GetBinContent(bin);
      butt += hist->GetBinContent(bin);
    }
  Double_t y_peak_low = top/butt;
  
  //printf("y_peak_low %f\n",y_peak_low);

  value=top=butt=0;
  for(Int_t ybin=m_low[x_u]; ybin <= m_up[x_u]; ybin++)
    {
      value = hist->GetBinCenterY(ybin);
      bin = hist->GetBin(x_u,ybin);
      top += value*hist->GetBinContent(bin);
      butt += hist->GetBinContent(bin);
    }
  Double_t y_peak_up = top/butt;
  
  //printf("y_peak_up %f\n",y_peak_up);

  Double_t x_value_up = hist->GetBinCenterX(x_u);
  Double_t x_value_low = hist->GetBinCenterX(x_l);

  Double_t y_peak = (y_peak_low*(x_value_up - x_peak) + y_peak_up*(x_peak - x_value_low))/(x_value_up - x_value_low);


  //Find the weight:
  bin = hist->FindBin(x_peak,y_peak);
  Int_t weight = (Int_t)hist->GetBinContent(bin);

  AliL3HoughTrack *track = new AliL3HoughTrack();
  track->SetTrackParameters(x_peak,y_peak,weight);
  
  //Reset area around peak
  for(Int_t xbin=x_low; xbin<=x_up; xbin++)
    {
      for(Int_t ybin=m_low[xbin]; ybin<=m_up[xbin]; ybin++)
	{
	  bin = hist->GetBin(xbin,ybin);
	  hist->SetBinContent(bin,0);
	}
    }
  
  delete [] m;
  delete [] m_low;
  delete [] m_up;
  
  return track;

  
}


AliL3HoughTrack *AliL3HoughMaxFinder::FindPeak(Int_t t1,Double_t t2,Int_t t3)
{
  //Attempt of a more sophisticated peak finder.
  
  if(!fCurrentHisto)
    {
      printf("AliL3HoughMaxFinder::FindPeak : No histogram!!\n");
      return 0;
    }
  AliL3Histogram *hist = fCurrentHisto;

  Int_t xmin = hist->GetFirstXbin();
  Int_t xmax = hist->GetLastXbin();
  Int_t ymin = hist->GetFirstYbin();
  Int_t ymax = hist->GetLastYbin();
  Int_t nbinsx = hist->GetNbinsX()+1;
  
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
      for(Int_t ybin=ymin; ybin <= ymax - t1; ybin++)
	{
	  sum = 0;
	  for(Int_t y=ybin; y <= ybin+t1; y++)
	    {
	      if(y>ymax) break;
	      //Inside window
	      bin = hist->GetBin(xbin,y);
	      sum += (Int_t)hist->GetBinContent(bin);
	      
	    }
	  if(sum > m[xbin]) //Max value locally in this xbin
	    {
	      m[xbin]=sum;
	      m_low[xbin]=ybin;
	      m_up[xbin]=ybin + t1;
	    }
	  
	}
      
      if(m[xbin] > max_x) //Max value globally in x-direction
	{
	  max_xbin = xbin;
	  max_x = m[xbin];//sum;
	}
    }
  //printf("max_xbin %d max_x %d m_low %d m_up %d\n",max_xbin,max_x,m_low[max_xbin],m_up[max_xbin]);
  //printf("ylow %f yup %f\n",hist->GetBinCenterY(m_low[max_xbin]),hist->GetBinCenterY(m_up[max_xbin]));

  //Determine a width in the x-direction
  Int_t x_low=0,x_up=0;
  
  for(Int_t xbin=max_xbin-1; xbin >= xmin; xbin--)
    {
      if(m[xbin] < max_x*t2)
	{
	  x_low = xbin+1;
	  break;
	}
    }
  for(Int_t xbin = max_xbin+1; xbin <=xmax; xbin++)
    {
      if(m[xbin] < max_x*t2)
	{
	  x_up = xbin-1;
	  break;
	}
    }
  
  Double_t top=0,butt=0,value,x_peak;
  if(x_up - x_low + 1 > t3)
    {
      t1 -= 1;
      printf("\nxrange out if limit x_up %d x_low %d t1 %d\n\n",x_low,x_up,t1);
      if(t1 > 1)
	goto recompute;
      else
	{
	  x_peak = hist->GetBinCenterX(max_xbin);
	  goto moveon;
	}
    }
  
  //printf("xlow %f xup %f\n",hist->GetBinCenterX(x_low),hist->GetBinCenterX(x_up));
  //printf("Spread in x %d\n",x_up-x_low +1);

  //Now, calculate the center of mass in x-direction
  for(Int_t xbin=x_low; xbin <= x_up; xbin++)
    {
      value = hist->GetBinCenterX(xbin);
      top += value*m[xbin];
      butt += m[xbin];
    }
  x_peak = top/butt;
  
 moveon:
  
  //Find the peak in y direction:
  Int_t x_l = hist->FindXbin(x_peak);
  if(hist->GetBinCenterX(x_l) > x_peak)
    x_l--;

  Int_t x_u = x_l + 1;
  
  if(hist->GetBinCenterX(x_l) > x_peak || hist->GetBinCenterX(x_u) <= x_peak)
    printf("\nAliL3HoughMaxFinder::FindPeak : Wrong xrange %f %f %f\n\n",hist->GetBinCenterX(x_l),x_peak,hist->GetBinCenterX(x_u));
    
    //printf("\nxlow %f xup %f\n",hist->GetBinCenterX(x_l),hist->GetBinCenterX(x_u));

  value=top=butt=0;
  
  //printf("ylow %f yup %f\n",hist->GetBinCenterY(m_low[x_l]),hist->GetBinCenterY(m_up[x_l]));
  //printf("ylow %f yup %f\n",hist->GetBinCenterY(m_low[x_u]),hist->GetBinCenterY(m_up[x_u]));
  
  for(Int_t ybin=m_low[x_l]; ybin <= m_up[x_l]; ybin++)
    {
      value = hist->GetBinCenterY(ybin);
      bin = hist->GetBin(x_l,ybin);
      top += value*hist->GetBinContent(bin);
      butt += hist->GetBinContent(bin);
    }
  Double_t y_peak_low = top/butt;
  
  //printf("y_peak_low %f\n",y_peak_low);

  value=top=butt=0;
  for(Int_t ybin=m_low[x_u]; ybin <= m_up[x_u]; ybin++)
    {
      value = hist->GetBinCenterY(ybin);
      bin = hist->GetBin(x_u,ybin);
      top += value*hist->GetBinContent(bin);
      butt += hist->GetBinContent(bin);
    }
  Double_t y_peak_up = top/butt;
  
  //printf("y_peak_up %f\n",y_peak_up);

  Double_t x_value_up = hist->GetBinCenterX(x_u);
  Double_t x_value_low = hist->GetBinCenterX(x_l);

  Double_t y_peak = (y_peak_low*(x_value_up - x_peak) + y_peak_up*(x_peak - x_value_low))/(x_value_up - x_value_low);


  //Find the weight:
  bin = hist->FindBin(x_peak,y_peak);
  Int_t weight = (Int_t)hist->GetBinContent(bin);

  AliL3HoughTrack *track = new AliL3HoughTrack();
  track->SetTrackParameters(x_peak,y_peak,weight);
  
  
  delete [] m;
  delete [] m_low;
  delete [] m_up;
  
  return track;

    
}

