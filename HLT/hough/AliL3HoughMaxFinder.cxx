//$Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ASV 

#include "AliL3StandardIncludes.h"

#ifndef no_root
#include <TNtuple.h>
#include <TFile.h>
#endif

#include "AliL3Logging.h"
#include "AliL3HoughMaxFinder.h"
#include "AliL3Histogram.h"
#include "AliL3TrackArray.h"
#include "AliL3HoughTrack.h"

#if GCCVERSION == 3
using namespace std;
#endif

//_____________________________________________________________
// AliL3HoughMaxFinder
//
// Maximum finder

ClassImp(AliL3HoughMaxFinder)

  
AliL3HoughMaxFinder::AliL3HoughMaxFinder()
{
  //Default constructor
  fThreshold = 0;
  fHistoType=0;
  fXPeaks=0;
  fYPeaks=0;
  fNPeaks=0;
  fNMax=0;
#ifndef no_root
  fNtuppel = 0;
#endif
}


AliL3HoughMaxFinder::AliL3HoughMaxFinder(Char_t *histotype,Int_t nmax,AliL3Histogram *hist)
{
  //Constructor

  //fTracks = new AliL3TrackArray("AliL3HoughTrack");
  if(strcmp(histotype,"KappaPhi")==0) fHistoType='c';
  if(strcmp(histotype,"DPsi")==0) fHistoType='l';
  
  if(hist)
    fCurrentHisto = hist;
  
  fNMax=nmax;
  fXPeaks = new Float_t[fNMax];
  fYPeaks = new Float_t[fNMax];
  fWeight = new Int_t[fNMax];
#ifndef no_root
  fNtuppel = 0;
#endif
  fThreshold=0;
}


AliL3HoughMaxFinder::~AliL3HoughMaxFinder()
{
  //Destructor
  if(fXPeaks)
    delete [] fXPeaks;
  if(fYPeaks)
    delete [] fYPeaks;
  if(fWeight)
    delete [] fWeight;
#ifndef no_root
  if(fNtuppel)
    delete fNtuppel;
#endif
}

void AliL3HoughMaxFinder::Reset()
{
  for(Int_t i=0; i<fNMax; i++)
    {
      fXPeaks[i]=0;
      fYPeaks[i]=0;
      fWeight[i]=0;
    }
  fNPeaks=0;
}

void AliL3HoughMaxFinder::CreateNtuppel()
{
#ifndef no_root
  //content#; neighbouring bins of the peak.
  fNtuppel = new TNtuple("ntuppel","Peak charateristics","kappa:phi0:weigth:content3:content5:content1:content7");
#endif  
}

void AliL3HoughMaxFinder::WriteNtuppel(Char_t *filename)
{
#ifndef no_root
  TFile *file = TFile::Open(filename,"RECREATE");
  if(!file)
    {
      cerr<<"AliL3HoughMaxFinder::WriteNtuppel : Error opening file "<<filename<<endl;
      return;
    }
  fNtuppel->Write();
  file->Close();
#endif
}

void AliL3HoughMaxFinder::FindAbsMaxima()
{
  
  if(!fCurrentHisto)
    {
      cerr<<"AliL3HoughMaxFinder::FindAbsMaxima : No histogram"<<endl;
      return;
    }
  AliL3Histogram *hist = fCurrentHisto;
  
  Int_t xmin = hist->GetFirstXbin();
  Int_t xmax = hist->GetLastXbin();
  Int_t ymin = hist->GetFirstYbin();
  Int_t ymax = hist->GetLastYbin();  
  Int_t bin;
  Double_t value,max_value=0;
  
  Int_t max_xbin=0,max_ybin=0;
  for(Int_t xbin=xmin; xbin<=xmax; xbin++)
    {
      for(Int_t ybin=ymin; ybin<=ymax; ybin++)
	{
	  bin = hist->GetBin(xbin,ybin);
	  value = hist->GetBinContent(bin);
	  if(value>max_value)
	    {
	      max_value = value;
	      max_xbin = xbin;
	      max_ybin = ybin;
	    }
	}
    }
  
  if(fNPeaks > fNMax)
    {
      cerr<<"AliL3HoughMaxFinder::FindAbsMaxima : Array out of range : "<<fNPeaks<<endl;
      return;
    }
  
      
  Double_t max_x = hist->GetBinCenterX(max_xbin);
  Double_t max_y = hist->GetBinCenterY(max_ybin);
  fXPeaks[fNPeaks] = max_x;
  fYPeaks[fNPeaks] = max_y;
  fWeight[fNPeaks] = (Int_t)max_value;
  fNPeaks++;
#ifndef no_root
  if(fNtuppel)
    {
      Int_t bin3 = hist->GetBin(max_xbin-1,max_ybin);
      Int_t bin5 = hist->GetBin(max_xbin+1,max_ybin);
      Int_t bin1 = hist->GetBin(max_xbin,max_ybin-1);
      Int_t bin7 = hist->GetBin(max_xbin,max_ybin+1);
      
      fNtuppel->Fill(max_x,max_y,max_value,hist->GetBinContent(bin3),hist->GetBinContent(bin5),hist->GetBinContent(bin1),hist->GetBinContent(bin7));
    }
#endif  
}

void AliL3HoughMaxFinder::FindBigMaxima()
{
  
  AliL3Histogram *hist = fCurrentHisto;
  Int_t xmin = hist->GetFirstXbin();
  Int_t xmax = hist->GetLastXbin();
  Int_t ymin = hist->GetFirstYbin();
  Int_t ymax = hist->GetLastYbin();
  Int_t bin[25],bin_index;
  Double_t value[25];
  
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
	      //printf("b %d\n",b);
	    }
	  if(b == bin_index)
	    {
	      //Found maxima
	      if(fNPeaks > fNMax)
		{
		  cerr<<"AliL3HoughMaxFinder::FindBigMaxima : Array out of range "<<fNPeaks<<endl;
		  return;
		}
	      
	      Double_t max_x = hist->GetBinCenterX(xbin);
	      Double_t max_y = hist->GetBinCenterY(ybin);
	      fXPeaks[fNPeaks] = max_x;
	      fYPeaks[fNPeaks] = max_y;
	      fNPeaks++;
	    }
	}
    }
}

void AliL3HoughMaxFinder::FindMaxima(Double_t grad_x,Double_t grad_y)
{
  //Locate all the maxima in input histogram.
  //Maxima is defined as bins with more entries than the
  //immediately neighbouring bins. 
  
  Int_t xmin = fCurrentHisto->GetFirstXbin();
  Int_t xmax = fCurrentHisto->GetLastXbin();
  Int_t ymin = fCurrentHisto->GetFirstYbin();
  Int_t ymax = fCurrentHisto->GetLastYbin();
  Int_t bin[9];
  Double_t value[9];
  
  for(Int_t xbin=xmin+1; xbin<xmax-1; xbin++)
    {
      for(Int_t ybin=ymin+1; ybin<ymax-1; ybin++)
	{
	  bin[0] = fCurrentHisto->GetBin(xbin-1,ybin-1);
	  bin[1] = fCurrentHisto->GetBin(xbin,ybin-1);
	  bin[2] = fCurrentHisto->GetBin(xbin+1,ybin-1);
	  bin[3] = fCurrentHisto->GetBin(xbin-1,ybin);
	  bin[4] = fCurrentHisto->GetBin(xbin,ybin);
	  bin[5] = fCurrentHisto->GetBin(xbin+1,ybin);
	  bin[6] = fCurrentHisto->GetBin(xbin-1,ybin+1);
	  bin[7] = fCurrentHisto->GetBin(xbin,ybin+1);
	  bin[8] = fCurrentHisto->GetBin(xbin+1,ybin+1);
	  value[0] = fCurrentHisto->GetBinContent(bin[0]);
	  value[1] = fCurrentHisto->GetBinContent(bin[1]);
	  value[2] = fCurrentHisto->GetBinContent(bin[2]);
	  value[3] = fCurrentHisto->GetBinContent(bin[3]);
	  value[4] = fCurrentHisto->GetBinContent(bin[4]);
	  value[5] = fCurrentHisto->GetBinContent(bin[5]);
	  value[6] = fCurrentHisto->GetBinContent(bin[6]);
	  value[7] = fCurrentHisto->GetBinContent(bin[7]);
	  value[8] = fCurrentHisto->GetBinContent(bin[8]);
	  
	  
	  
	  if(value[4]>value[0] && value[4]>value[1] && value[4]>value[2]
	     && value[4]>value[3] && value[4]>value[5] && value[4]>value[6]
	     && value[4]>value[7] && value[4]>value[8])
	    {
	      //Found a local maxima
	      Float_t max_x = fCurrentHisto->GetBinCenterX(xbin);
	      Float_t max_y = fCurrentHisto->GetBinCenterY(ybin);
	      
	      if((Int_t)value[4] <= fThreshold) continue;//central bin below threshold
	      if(fNPeaks >= fNMax)
		{
		  cout<<"AliL3HoughMaxFinder::FindMaxima : Array out of range "<<fNPeaks<<endl;
		  return;
		}

	      /*
	      //Check the gradient:
	      if(value[3]/value[4] > 1./grad_x || value[5]/value[4] < 1./grad_x ||
		 value[1]/value[4] < 1./grad_y || value[7]/value[4] < 1./grad_y)
		continue;
	      */


	      fXPeaks[fNPeaks] = max_x;
	      fYPeaks[fNPeaks] = max_y;
	      fWeight[fNPeaks] = (Int_t)value[4];
	      fNPeaks++;
	      
	      /*
	      //Check if the peak is overlapping with a previous:
	      Bool_t bigger = kFALSE;
	      for(Int_t p=0; p<entries; p++)
	      {
	        if(fabs(max_x - xpeaks[p]) < kappa_overlap && fabs(max_y - ypeaks[p]) < phi_overlap)
	      {
	      bigger = kTRUE;
		      if(value[4] > weight[p]) //this peak is bigger.
			{
			  xpeaks[p] = max_x;
			  ypeaks[p] = max_y;
			  weight[p] = (Int_t)value[4];
			}
		      else
			continue; //previous peak is bigger.
		    }
		}
	      if(!bigger) //there were no overlapping peaks.
		{
		xpeaks[entries] = max_x;
		  ypeaks[entries] = max_y;
		  weight[entries] = (Int_t)value[4];
		  entries++;
		}
	      */
	    }
	  else
	    continue; //not a maxima
	}
    }
  
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

AliL3HoughTrack *AliL3HoughMaxFinder::FindPeakLine(Double_t rho,Double_t theta)
{
  //Peak finder based on a second line transformation on kappa-phi space, 
  //to use as a baseline.

  if(!fCurrentHisto)
    {
      printf("AliL3HoughTransformer::FindPeakLine : No input histogram\n");
      return 0;
    }
  
  //get the line parameters:
  Double_t a = -1./tan(theta);
  Double_t b = rho/sin(theta);
  
  printf("rho %f theta %f\n",rho,theta);
  //now, start looking along the line.
  Int_t xmin = fCurrentHisto->GetFirstXbin();
  Int_t xmax = fCurrentHisto->GetLastXbin();
    
  Int_t max_weight=0;
  Double_t max_bin[2];
  for(Int_t xbin=xmin; xbin<=xmax; xbin++)
    {
      Double_t x = fCurrentHisto->GetBinCenterX(xbin);
      Double_t y = a*x + b;
      Int_t bin = fCurrentHisto->FindBin(x,y);
      //printf("x %f y %f weight %d\n",x,y,fCurrentHisto->GetBinContent(bin));
      if(fCurrentHisto->GetBinContent(bin) > max_weight)
	{
	  max_weight = (Int_t)fCurrentHisto->GetBinContent(bin);
	  max_bin[0] = x;
	  max_bin[1] = y;
	}
    }
  
  AliL3HoughTrack *track = new AliL3HoughTrack();
  track->SetTrackParameters(max_bin[0],max_bin[1],max_weight);
  return track;
}

void AliL3HoughMaxFinder::FindPeak1(Int_t y_window,Int_t x_bin_sides)
{
  //Testing mutliple peakfinding.
  //The algorithm searches the histogram for prepreaks by looking in windows
  //for each bin on the xaxis. The size of these windows is controlled by y_window.
  //Then the prepreaks are sorted according to their weight (sum inside window),
  //and the peak positions are calculated by taking the weighted mean in both
  //x and y direction. The size of the peak in x-direction is controlled by x_bin_sides.

  if(!fCurrentHisto)
    {
      printf("AliL3HoughMaxFinder::FindPeak1 : No input histogram\n");
      return;
    }  
  //Int_t y_window=2;
  //Int_t x_bin_sides=1;
  
  //Float_t max_kappa = 0.001;
  //Float_t max_phi0 = 0.05;
  
  Int_t max_sum=0;
  
  Int_t xmin = fCurrentHisto->GetFirstXbin();
  Int_t xmax = fCurrentHisto->GetLastXbin();
  Int_t ymin = fCurrentHisto->GetFirstYbin();
  Int_t ymax = fCurrentHisto->GetLastYbin();
  Int_t nbinsx = fCurrentHisto->GetNbinsX()+1;
  
  AxisWindow **windowPt = new AxisWindow*[nbinsx];
  AxisWindow **anotherPt = new AxisWindow*[nbinsx];
  
  for(Int_t i=0; i<nbinsx; i++)
    {
      windowPt[i] = new AxisWindow;
      bzero((void*)windowPt[i],sizeof(AxisWindow));
      anotherPt[i] = windowPt[i];
    }
  
  for(Int_t xbin=xmin; xbin<=xmax; xbin++)
    {
      max_sum = 0;
      for(Int_t ybin=ymin; ybin<=ymax-y_window; ybin++)
	{
	  Int_t sum_in_window=0;
	  for(Int_t b=ybin; b<ybin+y_window; b++)
	    {
	      //inside window
	      Int_t bin = fCurrentHisto->GetBin(xbin,b);
	      sum_in_window += (Int_t)fCurrentHisto->GetBinContent(bin);
	    }
	  
	  if(sum_in_window > max_sum)
	    {
	      max_sum = sum_in_window;
	      windowPt[xbin]->ymin = ybin;
	      windowPt[xbin]->ymax = ybin + y_window;
	      windowPt[xbin]->weight = sum_in_window;
	      windowPt[xbin]->xbin = xbin;
	    }
	}
    }

  //Sort the windows according to the weight
  SortPeaks(windowPt,0,nbinsx);
  
  Float_t top,butt;
  for(Int_t i=0; i<nbinsx; i++)
    {
      top=butt=0;
      Int_t xbin = windowPt[i]->xbin;
      
      if(xbin<xmin || xbin > xmax-1) continue;
      
      //Check if this is really a local maxima
      if(anotherPt[xbin-1]->weight > anotherPt[xbin]->weight ||
	 anotherPt[xbin+1]->weight > anotherPt[xbin]->weight)
	continue;

      for(Int_t j=windowPt[i]->ymin; j<windowPt[i]->ymax; j++)
	{
	  //Calculate the mean in y direction:
	  Int_t bin = fCurrentHisto->GetBin(windowPt[i]->xbin,j);
	  top += (fCurrentHisto->GetBinCenterY(j))*(fCurrentHisto->GetBinContent(bin));
	  butt += fCurrentHisto->GetBinContent(bin);
	}
      
      if(butt < fThreshold)
	continue;
      
      fXPeaks[fNPeaks] = fCurrentHisto->GetBinCenterX(windowPt[i]->xbin);
      fYPeaks[fNPeaks] = top/butt;
      fWeight[fNPeaks] = (Int_t)butt;
      //cout<<"mean in y "<<ypeaks[n]<<" on x "<<windowPt[i]->xbin<<" content "<<butt<<endl;
      fNPeaks++;
      if(fNPeaks==fNMax) break;
    }

  
  //Improve the peaks by including the region around in x.
  Float_t ytop,ybutt;
  Int_t prev;
  Int_t w;
  for(Int_t i=0; i<fNPeaks; i++)
    {
      Int_t xbin = fCurrentHisto->FindXbin(fXPeaks[i]);
      if(xbin - x_bin_sides < xmin || xbin + x_bin_sides > xmax) continue;
      top=butt=0;
      ytop=0,ybutt=0;	  
      w=0;
      prev = xbin - x_bin_sides+1;
      for(Int_t j=xbin-x_bin_sides; j<=xbin+x_bin_sides; j++)
	{
	  /*
	  //Check if the windows are overlapping:
	  if(anotherPt[j]->ymin > anotherPt[prev]->ymax) {prev=j; continue;}
	  if(anotherPt[j]->ymax < anotherPt[prev]->ymin) {prev=j; continue;}
	  prev = j;
	  */
	  
	  top += fCurrentHisto->GetBinCenterX(j)*anotherPt[j]->weight;
	  butt += anotherPt[j]->weight;
	  
	  for(Int_t k=anotherPt[j]->ymin; k<anotherPt[j]->ymax; k++)
	    {
	      Int_t bin = fCurrentHisto->GetBin(j,k);
	      ytop += (fCurrentHisto->GetBinCenterY(k))*(fCurrentHisto->GetBinContent(bin));
	      ybutt += fCurrentHisto->GetBinContent(bin);
	      w+=(Int_t)fCurrentHisto->GetBinContent(bin);
	    }
	}
      
      fXPeaks[i] = top/butt;
      fYPeaks[i] = ytop/ybutt;
      fWeight[i] = w;
      //cout<<"Setting weight "<<w<<" kappa "<<xpeaks[i]<<" phi0 "<<ypeaks[i]<<endl;
      
      /*
      //Check if this peak is overlapping with a previous:
      for(Int_t p=0; p<i-1; p++)
	{
	  if(fabs(fXPeaks[p] - fXPeaks[i]) < max_kappa ||
	     fabs(fYPeaks[p] - fYPeaks[i]) < max_phi0)
	    {
	      fWeight[i]=0;
	      break;
	    }
	}
      */
    }
  
  for(Int_t i=0; i<nbinsx; i++)
    delete windowPt[i];
  delete [] windowPt;
  delete [] anotherPt;
  
}

void AliL3HoughMaxFinder::SortPeaks(struct AxisWindow **a,Int_t first,Int_t last)
{
  //General sorting routine
  //Sort according to PeakCompare()

  static struct AxisWindow *tmp;
  static int i;           // "static" to save stack space
  int j;
  
  while (last - first > 1) {
    i = first;
    j = last;
    for (;;) {
      while (++i < last && PeakCompare(a[i], a[first]) < 0)
	;
      while (--j > first && PeakCompare(a[j], a[first]) > 0)
	;
      if (i >= j)
	break;
      
      tmp  = a[i];
      a[i] = a[j];
      a[j] = tmp;
    }
    if (j == first) {
      ++first;
      continue;
    }
    tmp = a[first];
    a[first] = a[j];
    a[j] = tmp;
    if (j - first < last - (j + 1)) {
      SortPeaks(a, first, j);
      first = j + 1;   // QSort(j + 1, last);
    } else {
      SortPeaks(a, j + 1, last);
      last = j;        // QSort(first, j);
    }
  }
  
}

Int_t AliL3HoughMaxFinder::PeakCompare(struct AxisWindow *a,struct AxisWindow *b)
{
  if(a->weight < b->weight) return 1;
  if(a->weight > b->weight) return -1;
  return 0;

}

void AliL3HoughMaxFinder::FindPeak(Int_t t1,Double_t t2,Int_t t3)
{
  //Attempt of a more sophisticated peak finder.
  //Finds the best peak in the histogram, and returns the corresponding
  //track object.

  if(!fCurrentHisto)
    {
      printf("AliL3HoughMaxFinder::FindPeak : No histogram!!\n");
      return;
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

  Int_t max_x=0,sum=0,max_xbin=0,bin=0;

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
  //bin = hist->FindBin(x_peak,y_peak);
  //Int_t weight = (Int_t)hist->GetBinContent(bin);

  //AliL3HoughTrack *track = new AliL3HoughTrack();
  //track->SetTrackParameters(x_peak,y_peak,weight);
  fXPeaks[fNPeaks]=x_peak;
  fYPeaks[fNPeaks]=y_peak;
  fWeight[fNPeaks]=(Int_t)hist->GetBinContent(bin);
  fNPeaks++;
  
  delete [] m;
  delete [] m_low;
  delete [] m_up;
  
  //return track;

    
}

