// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliL3StandardIncludes.h"

#include "AliL3Logging.h"
#include "AliL3MemHandler.h"
#include "AliL3Transform.h"
#include "AliL3DigitData.h"
#include "AliL3HistogramAdaptive.h"
#include "AliL3HoughTransformerRow.h"

#if __GNUC__ == 3
using namespace std;
#endif

//_____________________________________________________________
// AliL3HoughTransformerRow
//
// TPC rows Hough transformation class
//

ClassImp(AliL3HoughTransformerRow)

UChar_t **AliL3HoughTransformerRow::fRowCount = 0;
UChar_t **AliL3HoughTransformerRow::fGapCount = 0;
UChar_t **AliL3HoughTransformerRow::fCurrentRowCount = 0;
#ifdef do_mc
TrackIndex **AliL3HoughTransformerRow::fTrackID = 0;
#endif


AliL3HoughTransformerRow::AliL3HoughTransformerRow()
{
  //Default constructor
  fParamSpace = 0;
  fDoMC = kFALSE;;

  fLUT2sinphi0up=0;
  fLUT2sinphi0low=0;
  fLUT2cosphi0up=0;
  fLUT2cosphi0low=0;

  fLUTforwardZ=0;
  fLUTforwardZ2=0;
  fLUTbackwardZ=0;
  fLUTbackwardZ2=0;
}

AliL3HoughTransformerRow::AliL3HoughTransformerRow(Int_t slice,Int_t patch,Int_t n_eta_segments,Bool_t /*DoMC*/,Float_t zvertex) : AliL3HoughBaseTransformer(slice,patch,n_eta_segments,zvertex)
{
  //Normal constructor
  fParamSpace = 0;
  fDoMC = kFALSE;
  fDoMC=kFALSE;
#ifdef do_mc
  fDoMC = kTRUE;
#endif

  fLUT2sinphi0up=0;
  fLUT2sinphi0low=0;
  fLUT2cosphi0up=0;
  fLUT2cosphi0low=0;

  fLUTforwardZ=0;
  fLUTforwardZ2=0;
  fLUTbackwardZ=0;
  fLUTbackwardZ2=0;
}

AliL3HoughTransformerRow::~AliL3HoughTransformerRow()
{
  DeleteHistograms();
#ifdef do_mc
  if(fTrackID)
    {
      for(Int_t i=0; i<GetNEtaSegments(); i++)
	{
	  if(!fTrackID[i]) continue;
	  delete fTrackID[i];
	}
      delete [] fTrackID;
      fTrackID = 0;
    }
#endif

  if(fRowCount)
    {
      for(Int_t i=0; i<GetNEtaSegments(); i++)
	{
	  if(!fRowCount[i]) continue;
	  delete [] fRowCount[i];
	}
      delete [] fRowCount;
      fRowCount = 0;
    }
  if(fGapCount)
    {
      for(Int_t i=0; i<GetNEtaSegments(); i++)
	{
	  if(!fGapCount[i]) continue;
	  delete [] fGapCount[i];
	}
      delete [] fGapCount;
      fGapCount = 0;
    }
  if(fCurrentRowCount)
    {
      for(Int_t i=0; i<GetNEtaSegments(); i++)
	{
	  if(fCurrentRowCount[i])
	    delete [] fCurrentRowCount[i];
	}
      delete [] fCurrentRowCount;
      fCurrentRowCount = 0;
    }
}

void AliL3HoughTransformerRow::DeleteHistograms()
{
  delete[] fLUT2sinphi0up;
  delete[] fLUT2cosphi0up;
  delete[] fLUT2sinphi0low;
  delete[] fLUT2cosphi0low;

  fLUT2sinphi0up=0;
  fLUT2cosphi0up=0;
  fLUT2sinphi0low=0;
  fLUT2cosphi0low=0;

  delete[] fLUTforwardZ;
  delete[] fLUTforwardZ2;
  delete[] fLUTbackwardZ;
  delete[] fLUTbackwardZ2;

  fLUTforwardZ=0;
  fLUTforwardZ2=0;
  fLUTbackwardZ=0;
  fLUTbackwardZ2=0;

  if(!fParamSpace)
    return;
  for(Int_t i=0; i<GetNEtaSegments(); i++)
    {
      if(!fParamSpace[i]) continue;
      delete fParamSpace[i];
    }
  delete [] fParamSpace;
}

void AliL3HoughTransformerRow::CreateHistograms(Int_t nxbin,Float_t pt_min,
						Int_t nybin,Float_t phimin,Float_t phimax)
{
  //Create the histograms (parameter space).
  //These are 2D histograms, span by kappa (curvature of track) 
  //and phi0 (emission angle with x-axis).
  //The arguments give the range and binning; 
  //nxbin = #bins in kappa
  //nybin = #bins in phi0
  //pt_min = mimium Pt of track (corresponding to maximum kappa)
  //phi_min = mimimum phi0 
  //phi_max = maximum phi0 
    
  Double_t x = AliL3Transform::GetBFact()*AliL3Transform::GetBField()/pt_min;
  //Double_t torad = AliL3Transform::Pi()/180;
  CreateHistograms(nxbin,-1.*x,x,nybin,phimin/**torad*/,phimax/**torad*/);
}

void AliL3HoughTransformerRow::CreateHistograms(Int_t nxbin,Float_t xmin,Float_t xmax,
						Int_t nybin,Float_t ymin,Float_t ymax)
{
  
  fParamSpace = new AliL3Histogram*[GetNEtaSegments()];
  
  Char_t histname[256];
  for(Int_t i=0; i<GetNEtaSegments(); i++)
    {
      sprintf(histname,"paramspace_%d",i);
      fParamSpace[i] = new AliL3Histogram(histname,"",nxbin,xmin,xmax,nybin,ymin,ymax);
    }
  
#ifdef do_mc
  if(fDoMC)
    {
      AliL3Histogram *hist = fParamSpace[0];
      Int_t ncells = (hist->GetNbinsX()+2)*(hist->GetNbinsY()+2);
      if(!fTrackID)
	{
	  cout<<"Transformer: Allocating "<<GetNEtaSegments()*ncells*sizeof(TrackIndex)<<" bytes to fTrackID"<<endl;
	  fTrackID = new TrackIndex*[GetNEtaSegments()];
	  for(Int_t i=0; i<GetNEtaSegments(); i++)
	    fTrackID[i] = new TrackIndex[ncells];
	}
    }
#endif
  AliL3Histogram *hist = fParamSpace[0];
  Int_t ncells = (hist->GetNbinsX()+2)*(hist->GetNbinsY()+2);
  if(!fRowCount)
    {
      cout<<"Transformer: Allocating "<<GetNEtaSegments()*ncells*sizeof(UChar_t)<<" bytes to fRowCount"<<endl;
      fRowCount = new UChar_t*[GetNEtaSegments()];
      for(Int_t i=0; i<GetNEtaSegments(); i++)
	fRowCount[i] = new UChar_t[ncells];
    }
  if(!fGapCount)
    {
      cout<<"Transformer: Allocating "<<GetNEtaSegments()*ncells*sizeof(UChar_t)<<" bytes to fGapCount"<<endl;
      fGapCount = new UChar_t*[GetNEtaSegments()];
      for(Int_t i=0; i<GetNEtaSegments(); i++)
	fGapCount[i] = new UChar_t[ncells];
    }
  if(!fCurrentRowCount)
    {
      cout<<"Transformer: Allocating "<<GetNEtaSegments()*ncells*sizeof(UChar_t)<<" bytes to fCurrentRowCount"<<endl;
      fCurrentRowCount = new UChar_t*[GetNEtaSegments()];
      for(Int_t i=0; i<GetNEtaSegments(); i++)
	fCurrentRowCount[i] = new UChar_t[ncells];
    }

  //create lookup table for sin and cos
  Int_t nbins = hist->GetNbinsY()+2;
  fLUT2sinphi0up=new Float_t[nbins];
  fLUT2cosphi0up=new Float_t[nbins];
  fLUT2sinphi0low=new Float_t[nbins];
  fLUT2cosphi0low=new Float_t[nbins];
  Float_t hist_bin=hist->GetBinWidthY()/2.0;
  for(Int_t i=hist->GetFirstYbin(); i<=hist->GetLastYbin(); i++){
    Float_t phi0=hist->GetBinCenterY(i);
    fLUT2sinphi0low[i]=2.*sin(phi0+hist_bin);
    fLUT2cosphi0low[i]=2.*cos(phi0+hist_bin);
    fLUT2sinphi0up[i]=2.*sin(phi0-hist_bin);
    fLUT2cosphi0up[i]=2.*cos(phi0-hist_bin);
  }  

  //create lookup table for z of the digits
  Int_t ntimebins = AliL3Transform::GetNTimeBins();
  fLUTforwardZ = new Float_t[ntimebins];
  fLUTforwardZ2 = new Float_t[ntimebins];
  fLUTbackwardZ = new Float_t[ntimebins];
  fLUTbackwardZ2 = new Float_t[ntimebins];
  for(Int_t i=0; i<ntimebins; i++){
    Float_t z;
    z=AliL3Transform::GetZFast(0,i,GetZVertex());
    fLUTforwardZ[i]=z;
    fLUTforwardZ2[i]=z*z;
    z=AliL3Transform::GetZFast(18,i,GetZVertex());
    fLUTbackwardZ[i]=z;
    fLUTbackwardZ2[i]=z*z;
  }

}

void AliL3HoughTransformerRow::Reset()
{
  //Reset all the histograms

  if(!fParamSpace)
    {
      LOG(AliL3Log::kWarning,"AliL3HoughTransformer::Reset","Histograms")
	<<"No histograms to reset"<<ENDLOG;
      return;
    }
  
  for(Int_t i=0; i<GetNEtaSegments(); i++)
    fParamSpace[i]->Reset();
  
#ifdef do_mc
  if(fDoMC)
    {
      AliL3Histogram *hist = fParamSpace[0];
      Int_t ncells = (hist->GetNbinsX()+2)*(hist->GetNbinsY()+2);
      for(Int_t i=0; i<GetNEtaSegments(); i++)
	memset(fTrackID[i],0,ncells*sizeof(TrackIndex));
    }
#endif
  AliL3Histogram *hist = fParamSpace[0];
  Int_t ncells = (hist->GetNbinsX()+2)*(hist->GetNbinsY()+2);
  for(Int_t i=0; i<GetNEtaSegments(); i++)
    {
      memset(fRowCount[i],0,ncells*sizeof(UChar_t));
      memset(fGapCount[i],1,ncells*sizeof(UChar_t));
      memset(fCurrentRowCount[i],160,ncells*sizeof(UChar_t));
    }
}

Int_t AliL3HoughTransformerRow::GetEtaIndex(Double_t eta)
{
  //Return the histogram index of the corresponding eta. 

  Double_t etaslice = (GetEtaMax() - GetEtaMin())/GetNEtaSegments();
  Double_t index = (eta-GetEtaMin())/etaslice;
  return (Int_t)index;
}

inline AliL3Histogram *AliL3HoughTransformerRow::GetHistogram(Int_t eta_index)
{
  if(!fParamSpace || eta_index >= GetNEtaSegments() || eta_index < 0)
    return 0;
  if(!fParamSpace[eta_index])
    return 0;
  return fParamSpace[eta_index];
}

Double_t AliL3HoughTransformerRow::GetEta(Int_t eta_index,Int_t /*slice*/)
{
  Double_t eta_slice = (GetEtaMax()-GetEtaMin())/GetNEtaSegments();
  Double_t eta=0;
  eta=(Double_t)((eta_index+0.5)*eta_slice);
  return eta;
}

struct EtaRow {
  UChar_t start_pad;
  UChar_t end_pad;
  Bool_t found;
  Float_t start_y;
#ifdef do_mc
  Int_t mc_labels[MaxTrack];
#endif
};

void AliL3HoughTransformerRow::TransformCircle()
{
  Int_t n_eta_segments = GetNEtaSegments();
  Double_t eta_min = GetEtaMin();
  Double_t eta_slice = (GetEtaMax() - eta_min)/n_eta_segments;

  Int_t lower_threshold = GetLowerThreshold();

  //Assumes that all the etaslice histos are the same!
  AliL3Histogram *h = fParamSpace[0];
  Int_t first_bin = h->GetFirstYbin();
  Int_t last_bin = h->GetLastYbin();
  Float_t x_min = h->GetXmin();
  Float_t x_max = h->GetXmax();
  Float_t x_bin = (x_max-x_min)/h->GetNbinsX();
  Int_t first_binx = h->GetFirstXbin()-1;
  Int_t last_binx = h->GetLastXbin()+1;
  Int_t nbinx = h->GetNbinsX()+2;

  UChar_t last_pad;
  Int_t last_eta_index;
  EtaRow *eta_clust = new EtaRow[n_eta_segments];

  AliL3DigitRowData *tempPt = GetDataPointer();
  if(!tempPt)
    {
      LOG(AliL3Log::kError,"AliL3HoughTransformer::TransformCircle","Data")
	<<"No input data "<<ENDLOG;
      return;
    }

  Int_t ipatch = GetPatch();
  Double_t pad_pitch = AliL3Transform::GetPadPitchWidth(ipatch);
  Int_t islice = GetSlice();
  Float_t *fLUTz;
  Float_t *fLUTz2;
  if(islice < 18) {
    fLUTz = fLUTforwardZ;
    fLUTz2 = fLUTforwardZ2;
  }
  else {
    fLUTz = fLUTbackwardZ;
    fLUTz2 = fLUTbackwardZ2;
  }

  Int_t npads_in_patch = AliL3Transform::GetNPads(AliL3Transform::GetLastRow(ipatch))+2;
  Bool_t **IsPrevRow = new Bool_t*[n_eta_segments];
  Int_t nrows_in_patch = AliL3Transform::GetLastRow(ipatch)-AliL3Transform::GetFirstRow(ipatch)+2;
  for(Int_t i=0; i<n_eta_segments; i++) {
    IsPrevRow[i] = new Bool_t[npads_in_patch*nrows_in_patch];
    memset(IsPrevRow[i],1,npads_in_patch*nrows_in_patch*sizeof(Bool_t));
  }

  //Loop over the padrows:
  for(UChar_t i=AliL3Transform::GetFirstRow(ipatch); i<=AliL3Transform::GetLastRow(ipatch); i++)
    {
      Int_t npads = AliL3Transform::GetNPads((Int_t)i)-1;

      last_pad = 255;
      //Flush eta clusters array
      memset(eta_clust,0,n_eta_segments*sizeof(EtaRow));  

      Float_t x = AliL3Transform::Row2X((Int_t)i);
      Float_t x2 = x*x;
      Float_t y=0,r2=0;

      Int_t lrow = (i-AliL3Transform::GetFirstRow(ipatch))*npads_in_patch;

      //Get the data on this padrow:
      AliL3DigitData *digPt = tempPt->fDigitData;
      if((Int_t)i != (Int_t)tempPt->fRow)
	{
	  cerr<<"AliL3HoughTransform::TransformCircle : Mismatching padrow numbering "<<(Int_t)i<<" "<<(Int_t)tempPt->fRow<<endl;
	  continue;
	}
      //      cout<<" Starting row "<<i<<endl;
      //Loop over the data on this padrow:
      for(UInt_t j=0; j<tempPt->fNDigit; j++)
	{
	  UShort_t charge = digPt[j].fCharge;
	  if((Int_t)charge <= lower_threshold)
	    continue;
	  UChar_t pad = digPt[j].fPad;
	  UShort_t time = digPt[j].fTime;
	  if(pad != last_pad)
	    {
	      y = (pad-0.5*npads)*pad_pitch;
	      Float_t y2 = y*y;
	      r2 = x2 + y2;
	      last_eta_index = -1;
	    }

	  //Transform data to local cartesian coordinates:
	  Float_t z = fLUTz[(Int_t)time];
	  Float_t z2 = fLUTz2[(Int_t)time];
	  //Calculate the eta:
	  Double_t r = sqrt(r2+z2);
	  Double_t eta = 0.5 * log((r+z)/(r-z));
	  //Get the corresponding index, which determines which histogram to fill:
	  Int_t eta_index = (Int_t)((eta-eta_min)/eta_slice);

#ifndef do_mc
	  if(eta_index == last_eta_index) continue;
#endif
	  //	  cout<<" Digit at patch "<<ipatch<<" row "<<i<<" pad "<<(Int_t)pad<<" time "<<time<<" etaslice "<<eta_index<<endl;
	  
	  if(eta_index < 0 || eta_index >= n_eta_segments)
	    continue;

	  if(!eta_clust[eta_index].found)
	    {
	      eta_clust[eta_index].start_pad = pad;
	      eta_clust[eta_index].end_pad = pad;
	      eta_clust[eta_index].start_y = y - pad_pitch/2.0;
	      eta_clust[eta_index].found = 1;
#ifdef do_mc
	      if(fDoMC)
		{
		  for(Int_t t=0; t<3; t++)
		    {
		      Int_t label = digPt[j].fTrackID[t];
		      if(label < 0) break;
		      UInt_t c;
		      for(c=0; c<(MaxTrack-1); c++)
			if(eta_clust[eta_index].mc_labels[c] == label || eta_clust[eta_index].mc_labels[c] == 0)
			  break;
		      //		      if(c == MaxTrack-1) cerr<<"AliL3HoughTransformer::TransformCircle : Cluster array reached maximum!! "<<c<<endl;
		      eta_clust[eta_index].mc_labels[c] = label;
		    }
		}
#endif
	      continue;
	    }
	  else
	    {
	      if(pad <= (eta_clust[eta_index].end_pad + 1))
		{
		  eta_clust[eta_index].end_pad = pad;
#ifdef do_mc
		  if(fDoMC)
		    {
		      for(Int_t t=0; t<3; t++)
			{
			  Int_t label = digPt[j].fTrackID[t];
			  if(label < 0) break;
			  UInt_t c;
			  for(c=0; c<(MaxTrack-1); c++)
			    if(eta_clust[eta_index].mc_labels[c] == label || eta_clust[eta_index].mc_labels[c] == 0)
			      break;
			  //			  if(c == MaxTrack-1) cerr<<"AliL3HoughTransformer::TransformCircle : Cluster array reached maximum!! "<<c<<endl;
			  eta_clust[eta_index].mc_labels[c] = label;
			}
		    }
#endif
		}
	      else
		{

		  Bool_t fill_cluster = kFALSE;
		  Bool_t *IsPrevRow2 = &IsPrevRow[eta_index][lrow+eta_clust[eta_index].start_pad];
		  for(Int_t ipad = 0; ipad <= (eta_clust[eta_index].end_pad - eta_clust[eta_index].start_pad + 2); ipad++)
		    {
		      if(*IsPrevRow2)
			{
			  fill_cluster = kTRUE;
			  break;
			}
		      IsPrevRow2++;
		    }

		  //		  if(eta_index == 51)
		  //		    cout<<" Cluster to fill:"<<" row "<<(Int_t)i<<" pads from "<<(Int_t)eta_clust[eta_index].start_pad<<"("<<eta_clust[eta_index].start_y<<") to "<<(Int_t)eta_clust[eta_index].end_pad<<"("<<(eta_clust[eta_index].end_pad-0.5*(npads-1))*pad_pitch<<") fill "<<fill_cluster<<endl;

		  if(fill_cluster) {
		  //		  cout<<" Cluster found at etaslice "<<eta_index<<" from pad "<<(Int_t)eta_clust[eta_index].start_pad<<" to pad "<<(Int_t)eta_clust[eta_index].end_pad<<endl;
		  //		  Bool_t fill_next_row = kFALSE;

		  UChar_t *nrows = fRowCount[eta_index];
		  UChar_t *ngaps = fGapCount[eta_index];
		  UChar_t *currentrow = fCurrentRowCount[eta_index];

		  //Do the transformation:
		  Float_t start_y = eta_clust[eta_index].start_y;
		  if(eta_clust[eta_index].start_pad == 0)
		    start_y -= 3.0*pad_pitch;
		  Float_t R1 = x2 + start_y*start_y;
		  Float_t x_over_R1 = x/R1;
		  Float_t start_y_over_R1 = start_y/R1;
		  Float_t end_y = (eta_clust[eta_index].end_pad-0.5*(npads-1))*pad_pitch;
		  if(eta_clust[eta_index].end_pad == npads)
		    end_y += 3.0*pad_pitch;
		  Float_t R2 = x2 + end_y*end_y; 
		  Float_t x_over_R2 = x/R2;
		  Float_t end_y_over_R2 = end_y/R2;

		  //Fill the histogram along the phirange
		  for(Int_t b=first_bin; b<=last_bin; b++)
		    {
		      Float_t kappa1 = start_y_over_R1*fLUT2cosphi0low[b]-x_over_R1*fLUT2sinphi0low[b];
		      if(kappa1>x_max) continue;
		      Float_t kappa2 = end_y_over_R2*fLUT2cosphi0up[b]-x_over_R2*fLUT2sinphi0up[b];
		      if(kappa2<x_min) break;

		      Int_t binx1 = 1 + (Int_t)((kappa1-x_min)/x_bin);
		      if(binx1<first_binx) binx1 = first_binx;
		      Int_t binx2 = 1 + (Int_t)((kappa2-x_min)/x_bin);
		      if(binx2>last_binx) binx2 = last_binx;
		      //		      if(eta_index == 51)
		      //			cout<<"     phi bin "<<b<<" kappa1 "<<kappa1<<" kappa2 "<<kappa2<<" from bin "<<binx1<<" to "<<binx2<<endl;
		      Int_t temp_bin = b*nbinx + binx1;
		      UChar_t *nrows2 = nrows + temp_bin;
		      UChar_t *ngaps2 = ngaps + temp_bin;
		      UChar_t *currentrow2 = currentrow + temp_bin;
		      for(Int_t bin=binx1;bin<=binx2;bin++)
			{
			  //Assumes threshold about 32
			  if(*ngaps2 < 3) {
			    if(i != *currentrow2)
			      {
				(*nrows2)++;
				if(i > (*currentrow2 + 1))
				  (*ngaps2)++;
				*currentrow2=i;
			      }
#ifdef do_mc
			    if(fDoMC)
			      {
				for(UInt_t t=0;t<(MaxTrack-1); t++)
				  {
				    Int_t label = eta_clust[eta_index].mc_labels[t];
				    if(label == 0) break;
				    UInt_t c;
				    Int_t temp_bin2 = b*nbinx + bin;
				    for(c=0; c<(MaxTrack-1); c++)
				      if(fTrackID[eta_index][temp_bin2].fLabel[c] == label || fTrackID[eta_index][temp_bin2].fNHits[c] == 0)
					break;
				    //				    if(c == MaxTrack-1) cerr<<"AliL3HoughTransformer::TransformCircle : Array reached maximum!! "<<c<<endl;
				    fTrackID[eta_index][temp_bin2].fLabel[c] = label;
				    if(fTrackID[eta_index][temp_bin2].fCurrentRow[c] != i) {
				      fTrackID[eta_index][temp_bin2].fNHits[c]++;
				      fTrackID[eta_index][temp_bin2].fCurrentRow[c] = i;
				    }
				  }
			      }
#endif
			  }
			  nrows2++;
			  ngaps2++;
			  currentrow2++;
			}
		    }
		  //End of the transformation

		  /*		  if(!fill_next_row) {
		    Bool_t *IsCurrentRow = &IsPrevRow[eta_index][lrow+npads_in_patch];
		    memset(&IsCurrentRow[eta_clust[eta_index].start_pad + 1],0,eta_clust[eta_index].end_pad - eta_clust[eta_index].start_pad + 1);
		    } */
		  }
		  else {
		    Bool_t *IsCurrentRow = &IsPrevRow[eta_index][lrow+npads_in_patch+eta_clust[eta_index].start_pad+1];
		    memset(IsCurrentRow,0,eta_clust[eta_index].end_pad - eta_clust[eta_index].start_pad + 1);
		  }

		  Bool_t *IsCurrentRow = &IsPrevRow[eta_index][lrow+npads_in_patch+eta_clust[eta_index].end_pad+2];
		  memset(IsCurrentRow,0,pad - eta_clust[eta_index].end_pad - 1);

		  eta_clust[eta_index].start_pad = pad;
		  eta_clust[eta_index].end_pad = pad;
		  eta_clust[eta_index].start_y = y - pad_pitch/2.0;

#ifdef do_mc
		  if(fDoMC)
		    {
		      memset(eta_clust[eta_index].mc_labels,0,MaxTrack);
		      for(Int_t t=0; t<3; t++)
			{
			  Int_t label = digPt[j].fTrackID[t];
			  if(label < 0) break;
			  UInt_t c;
			  for(c=0; c<(MaxTrack-1); c++)
			    if(eta_clust[eta_index].mc_labels[c] == label || eta_clust[eta_index].mc_labels[c] == 0)
			      break;
			  //			  if(c == MaxTrack-1) cerr<<"AliL3HoughTransformer::TransformCircle : Cluster array reached maximum!! "<<c<<endl;
			  eta_clust[eta_index].mc_labels[c] = label;
			}
		    }
#endif
		}
	    }
	  last_pad = pad;
	  last_eta_index = eta_index;
	}
      //Write remaining clusters
      for(Int_t eta_index = 0;eta_index < n_eta_segments;eta_index++)
	{
	  //Check for empty row
	  if((eta_clust[eta_index].start_pad == 0) && (eta_clust[eta_index].end_pad == 0)) continue;

	  UChar_t *nrows = fRowCount[eta_index];
	  UChar_t *ngaps = fGapCount[eta_index];
	  UChar_t *currentrow = fCurrentRowCount[eta_index];

	  //Do the transformation:
	  Float_t start_y = eta_clust[eta_index].start_y;
	  if(eta_clust[eta_index].start_pad == 0)
	    start_y -= 3.0*pad_pitch;
	  Float_t R1 = x2 + start_y*start_y; 
	  Float_t x_over_R1 = x/R1;
	  Float_t start_y_over_R1 = start_y/R1;
	  Float_t end_y = (eta_clust[eta_index].end_pad-0.5*(npads-1))*pad_pitch;
	  if(eta_clust[eta_index].end_pad == npads)
	    end_y += 3.0*pad_pitch;
	  Float_t R2 = x2 + end_y*end_y; 
	  Float_t x_over_R2 = x/R2;
	  Float_t end_y_over_R2 = end_y/R2;

	  //Fill the histogram along the phirange
	  for(Int_t b=first_bin; b<=last_bin; b++)
	    {
	      Float_t kappa1 = start_y_over_R1*fLUT2cosphi0low[b]-x_over_R1*fLUT2sinphi0low[b];
	      if(kappa1>x_max) continue;
	      Float_t kappa2 = end_y_over_R2*fLUT2cosphi0up[b]-x_over_R2*fLUT2sinphi0up[b];
	      if(kappa2<x_min) break;

	      Int_t binx1 = 1 + (Int_t)((kappa1-x_min)/x_bin);
	      if(binx1<first_binx) binx1 = first_binx;
	      Int_t binx2 = 1 + (Int_t)((kappa2-x_min)/x_bin);
	      if(binx2>last_binx) binx2 = last_binx;

	      Int_t temp_bin = b*nbinx + binx1;
	      UChar_t *nrows2 = nrows + temp_bin;
	      UChar_t *ngaps2 = ngaps + temp_bin;
	      UChar_t *currentrow2 = currentrow + temp_bin;
	      for(Int_t bin=binx1;bin<=binx2;bin++)
		{
		  //Assumes threshold about 32
		  if(*ngaps2 < 3) {
		    if(i != *currentrow2)
		      {
			(*nrows2)++;
			if(i > (*currentrow2 + 1))
			  (*ngaps2)++;
			*currentrow2=i;
		      }
#ifdef do_mc
		    if(fDoMC)
		      {
			for(UInt_t t=0;t<(MaxTrack-1); t++)
			  {
			    Int_t label = eta_clust[eta_index].mc_labels[t];
			    if(label == 0) break;
			    UInt_t c;
			    Int_t temp_bin2 = b*nbinx + bin;
			    for(c=0; c<(MaxTrack-1); c++)
			      if(fTrackID[eta_index][temp_bin2].fLabel[c] == label || fTrackID[eta_index][temp_bin2].fNHits[c] == 0)
				break;
			    //			    if(c == MaxTrack-1) cout<<"AliL3HoughTransformer::TransformCircle : Array reached maximum!! "<<c<<endl;
			    fTrackID[eta_index][temp_bin2].fLabel[c] = label;
			    if(fTrackID[eta_index][temp_bin2].fCurrentRow[c] != i) {
			      fTrackID[eta_index][temp_bin2].fNHits[c]++;
			      fTrackID[eta_index][temp_bin2].fCurrentRow[c] = i;
			    }
			  }
		      }
#endif
		  }
		  nrows2++;
		  ngaps2++;
		  currentrow2++;
		}
	    }
	  //End of the transformation

	  if(eta_clust[eta_index].end_pad < npads)
	    {
	      Bool_t *IsCurrentRow = &IsPrevRow[eta_index][lrow+npads_in_patch+eta_clust[eta_index].end_pad+2];
	      memset(IsCurrentRow,0,npads - eta_clust[eta_index].end_pad);
	    }

	}

      //Move the data pointer to the next padrow:
      AliL3MemHandler::UpdateRowPointer(tempPt);
    }


  for(Int_t i=0; i<n_eta_segments; i++)
    {
      delete [] IsPrevRow[i];
    }
  delete [] IsPrevRow;

  delete [] eta_clust;
}

Int_t AliL3HoughTransformerRow::GetTrackID(Int_t eta_index,Double_t kappa,Double_t psi)
{
  if(!fDoMC)
    {
      cerr<<"AliL3HoughTransformer::GetTrackID : Flag switched off"<<endl;
      return -1;
    }
  
#ifdef do_mc
  if(eta_index < 0 || eta_index > GetNEtaSegments())
    {
      cerr<<"AliL3HoughTransformer::GetTrackID : Wrong etaindex "<<eta_index<<endl;
      return -1;
    }
  AliL3Histogram *hist = fParamSpace[eta_index];
  Int_t bin = hist->FindBin(kappa,psi);
  Int_t label=-1;
  Int_t max=0;
  for(UInt_t i=0; i<(MaxTrack-1); i++)
    {
      Int_t nhits=fTrackID[eta_index][bin].fNHits[i];
      if(nhits == 0) break;
      if(nhits > max)
	{
	  max = nhits;
	  label = fTrackID[eta_index][bin].fLabel[i];
	}
    }
  Int_t label2=-1;
  Int_t max2=0;
  for(UInt_t i=0; i<(MaxTrack-1); i++)
    {
      Int_t nhits=fTrackID[eta_index][bin].fNHits[i];
      if(nhits == 0) break;
      if(nhits > max2)
	{
	  if(fTrackID[eta_index][bin].fLabel[i]!=label) {
	    max2 = nhits;
	    label2 = fTrackID[eta_index][bin].fLabel[i];
	  }
	}
    }
  //  cout<<" TrackID label "<<label<<" max "<<max<<" label2 "<<label2<<" max2 "<<max2<<" "<<(Float_t)max2/(Float_t)max<<" "<<fTrackID[eta_index][bin].fLabel[MaxTrack-1]<<" "<<(Int_t)fTrackID[eta_index][bin].fNHits[MaxTrack-1]<<endl;
  return label;
#endif
  cout<<"AliL3HoughTransformer::GetTrackID : Compile with do_mc flag!"<<endl;
  return -1;
}

UChar_t *AliL3HoughTransformerRow::GetRowCount(Int_t eta_index)
{
  return fRowCount[eta_index];
}

UChar_t *AliL3HoughTransformerRow::GetGapCount(Int_t eta_index)
{
  return fGapCount[eta_index];
}
