// @(#) $Id$


// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliL3StandardIncludes.h"

#include "AliL3Logging.h"
#include "AliL3MemHandler.h"
#include "AliL3Transform.h"
#include "AliL3DigitData.h"
#include "AliL3HistogramAdaptive.h"
#include "AliL3HoughTrack.h"
#include "AliL3HoughTransformerRow.h"

#if GCCVERSION == 3
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
UChar_t *AliL3HoughTransformerRow::fTrackNRows = 0;
UChar_t *AliL3HoughTransformerRow::fTrackFirstRow = 0;
UChar_t *AliL3HoughTransformerRow::fTrackLastRow = 0;

Float_t AliL3HoughTransformerRow::fBeta1 = AliL3Transform::Row2X(79)/pow(AliL3Transform::Row2X(79),2);
Float_t AliL3HoughTransformerRow::fBeta2 = (AliL3Transform::Row2X(158)+6.0)/pow((AliL3Transform::Row2X(158)+6.0),2);

AliL3HoughTransformerRow::AliL3HoughTransformerRow()
{
  //Default constructor
  fParamSpace = 0;
  fDoMC = kFALSE;;

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
  if(fTrackNRows)
    {
      delete [] fTrackNRows;
      fTrackNRows = 0;
    }
  if(fTrackFirstRow)
    {
      delete [] fTrackFirstRow;
      fTrackFirstRow = 0;
    }
  if(fTrackLastRow)
    {
      delete [] fTrackLastRow;
      fTrackLastRow = 0;
    }
}

void AliL3HoughTransformerRow::DeleteHistograms()
{
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
	  LOG(AliL3Log::kInformational,"AliL3HoughTransformerRow::CreateHistograms()","")
	    <<"Transformer: Allocating "<<GetNEtaSegments()*ncells*sizeof(TrackIndex)<<" bytes to fTrackID"<<ENDLOG;
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
      LOG(AliL3Log::kInformational,"AliL3HoughTransformerRow::CreateHistograms()","")
	<<"Transformer: Allocating "<<GetNEtaSegments()*ncells*sizeof(UChar_t)<<" bytes to fRowCount"<<ENDLOG;
      fRowCount = new UChar_t*[GetNEtaSegments()];
      for(Int_t i=0; i<GetNEtaSegments(); i++)
	fRowCount[i] = new UChar_t[ncells];
    }
  if(!fGapCount)
    {
      LOG(AliL3Log::kInformational,"AliL3HoughTransformerRow::CreateHistograms()","")
	<<"Transformer: Allocating "<<GetNEtaSegments()*ncells*sizeof(UChar_t)<<" bytes to fGapCount"<<ENDLOG;
      fGapCount = new UChar_t*[GetNEtaSegments()];
      for(Int_t i=0; i<GetNEtaSegments(); i++)
	fGapCount[i] = new UChar_t[ncells];
    }
  if(!fCurrentRowCount)
    {
      LOG(AliL3Log::kInformational,"AliL3HoughTransformerRow::CreateHistograms()","")
	<<"Transformer: Allocating "<<GetNEtaSegments()*ncells*sizeof(UChar_t)<<" bytes to fCurrentRowCount"<<ENDLOG;
      fCurrentRowCount = new UChar_t*[GetNEtaSegments()];
      for(Int_t i=0; i<GetNEtaSegments(); i++)
	fCurrentRowCount[i] = new UChar_t[ncells];
    }

  if(!fTrackNRows)
    {
      LOG(AliL3Log::kInformational,"AliL3HoughTransformerRow::CreateHistograms()","")
	<<"Transformer: Allocating "<<ncells*sizeof(UChar_t)<<" bytes to fTrackNRows"<<ENDLOG;
      fTrackNRows = new UChar_t[ncells];
      LOG(AliL3Log::kInformational,"AliL3HoughTransformerRow::CreateHistograms()","")
	<<"Transformer: Allocating "<<ncells*sizeof(UChar_t)<<" bytes to fTrackFirstRow"<<ENDLOG;
      fTrackFirstRow = new UChar_t[ncells];
      LOG(AliL3Log::kInformational,"AliL3HoughTransformerRow::CreateHistograms()","")
	<<"Transformer: Allocating "<<ncells*sizeof(UChar_t)<<" bytes to fTrackLastRow"<<ENDLOG;
      fTrackLastRow = new UChar_t[ncells];

      AliL3HoughTrack track;
      Int_t xmin = hist->GetFirstXbin();
      Int_t xmax = hist->GetLastXbin();
      Int_t ymin = hist->GetFirstYbin();
      Int_t ymax = hist->GetLastYbin();
      Int_t nxbins = hist->GetNbinsX()+2;
      for(Int_t ybin=ymin; ybin<=ymax; ybin++)
	{
	  for(Int_t xbin=xmin; xbin<=xmax; xbin++)
	    {
	      //cvetan: we get strange warning on gcc-2.95
	      //warning: large integer implicitly truncated to unsigned type
	      fTrackNRows[xbin + ybin*nxbins] = 99999;
	      for(Int_t deltay = -1; deltay <= 1; deltay += 2) {
		for(Int_t deltax = -1; deltax <= 1; deltax += 2) {
		
		  Float_t xtrack = hist->GetPreciseBinCenterX((Float_t)xbin+0.5*(Float_t)deltax);
		  Float_t ytrack = hist->GetPreciseBinCenterY((Float_t)ybin+0.5*(Float_t)deltay);

		  Float_t psi = atan((xtrack-ytrack)/(fBeta1-fBeta2));
		  Float_t kappa = 2.0*(xtrack*cos(psi)-fBeta1*sin(psi));
		  track.SetTrackParameters(kappa,psi,1);
		  Bool_t first_row = kFALSE;
		  UInt_t maxfirstrow = 0;
		  UInt_t maxlastrow = 0;
		  UInt_t curfirstrow = 0;
		  UInt_t curlastrow = 0;
		  for(Int_t j=AliL3Transform::GetFirstRow(0); j<=AliL3Transform::GetLastRow(5); j++)
		    {
		      Float_t hit[3];
		      if(!track.GetCrossingPoint(j,hit)) continue;
		      AliL3Transform::LocHLT2Raw(hit,0,j);
		      if(hit[1]>=0 && hit[1]<AliL3Transform::GetNPads(j))
			{
			  if(!first_row) {
			    curfirstrow = j;
			    first_row = kTRUE;
			  }
			  curlastrow = j;
			}
		      else {
			if(first_row) {
			  first_row = kFALSE;
			  if((curlastrow-curfirstrow) >= (maxlastrow-maxfirstrow)) {
			    maxfirstrow = curfirstrow;
			    maxlastrow = curlastrow;
			  }
			}
		      }
		    }
		  if((curlastrow-curfirstrow) >= (maxlastrow-maxfirstrow)) {
		    maxfirstrow = curfirstrow;
		    maxlastrow = curlastrow;
		  }
		  if((maxlastrow-maxfirstrow) < fTrackNRows[xbin + ybin*nxbins]) {
		    fTrackNRows[xbin + ybin*nxbins] = maxlastrow-maxfirstrow;
		    fTrackFirstRow[xbin + ybin*nxbins] = maxfirstrow;
		    fTrackLastRow[xbin + ybin*nxbins] = maxlastrow;
		  }
		}
	      }
	      //	      cout<<" fTrackNRows "<<xbin<<" "<<ybin<<" "<<(Int_t)fTrackNRows[xbin + ybin*nxbins]<<" "<<endl;
	    }
	}
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
      //      memset(fCurrentRowCount[i],160,ncells*sizeof(UChar_t));
      memcpy(fCurrentRowCount[i],fTrackFirstRow,ncells*sizeof(UChar_t));
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
  Float_t beta1 = fBeta1;
  Float_t beta2 = fBeta2;
  Float_t beta1_minus_beta2 = fBeta1 - fBeta2;

  Int_t n_eta_segments = GetNEtaSegments();
  Double_t eta_min = GetEtaMin();
  Double_t eta_slice = (GetEtaMax() - eta_min)/n_eta_segments;

  Int_t lower_threshold = GetLowerThreshold();

  //Assumes that all the etaslice histos are the same!
  AliL3Histogram *h = fParamSpace[0];
  Float_t y_min = h->GetYmin();
  //Float_t y_max = h->GetYmax();
  Float_t hist_bin = h->GetBinWidthY();
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

      //Get the data on this padrow:
      AliL3DigitData *digPt = tempPt->fDigitData;
      if((Int_t)i != (Int_t)tempPt->fRow)
	{
	  LOG(AliL3Log::kError,"AliL3HoughTransformerRow::TransformCircle","Data")
	    <<"AliL3HoughTransform::TransformCircle : Mismatching padrow numbering "<<(Int_t)i<<" "<<(Int_t)tempPt->fRow<<ENDLOG;
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
		  Bool_t fill_cluster = kTRUE;
		  if(fill_cluster) {

		  UChar_t *nrows = fRowCount[eta_index];
		  UChar_t *ngaps = fGapCount[eta_index];
		  UChar_t *currentrow = fCurrentRowCount[eta_index];
		  UChar_t *lastrow = fTrackLastRow;

		  //Do the transformation:
		  Float_t start_y = eta_clust[eta_index].start_y;
		  if(eta_clust[eta_index].start_pad == 0)
		    start_y -= 0.0*pad_pitch;
		  Float_t R1 = x2 + start_y*start_y;
		  Float_t x_over_R1 = x/R1;
		  Float_t start_y_over_R1 = start_y/R1;
		  Float_t end_y = (eta_clust[eta_index].end_pad-0.5*(npads-1))*pad_pitch;
		  if(eta_clust[eta_index].end_pad == npads)
		    end_y += 0.0*pad_pitch;
		  Float_t R2 = x2 + end_y*end_y; 
		  Float_t x_over_R2 = x/R2;
		  Float_t end_y_over_R2 = end_y/R2;
		  Float_t A1 = beta1_minus_beta2/(x_over_R1-beta2);
		  Float_t B1 = (x_over_R1-beta1)/(x_over_R1-beta2);
		  Float_t A2 = beta1_minus_beta2/(x_over_R2-beta2);
		  Float_t B2 = (x_over_R2-beta1)/(x_over_R2-beta2);

		  Float_t kappa1 = (A1*start_y_over_R1+B1*y_min-x_min)/x_bin;
		  Float_t delta_kappa1 = B1*hist_bin/x_bin;
		  if(B1<0)
		    kappa1 += delta_kappa1;
		  Float_t kappa2 = (A2*end_y_over_R2+B2*y_min-x_min)/x_bin;
		  Float_t delta_kappa2 = B2*hist_bin/x_bin;
		  if(B2>=0)
		    kappa2 += delta_kappa2;

		  //Fill the histogram along the phirange
		  for(Int_t b=first_bin; b<=last_bin; b++, kappa1 += delta_kappa1, kappa2 += delta_kappa2)
		    {
		      Int_t binx1 = 1 + (Int_t)kappa1;
		      if(binx1>last_binx) continue;
		      if(binx1<first_binx) binx1 = first_binx;
		      Int_t binx2 = 1 + (Int_t)kappa2;
		      if(binx2<first_binx) continue;
		      if(binx2>last_binx) binx2 = last_binx;
#ifdef do_mc
		      if(binx2<binx1) {
			LOG(AliL3Log::kWarning,"AliL3HoughTransformerRow::TransformCircle()","")
			  <<"Wrong filling "<<binx1<<" "<<binx2<<" "<<i<<" "<<x<<" "<<start_y<<" "<<end_y<<ENDLOG;
		      }
#endif
		      Int_t temp_bin = b*nbinx;
		      UChar_t *nrows2 = nrows + temp_bin;
		      UChar_t *ngaps2 = ngaps + temp_bin;
		      UChar_t *currentrow2 = currentrow + temp_bin;
		      UChar_t *lastrow2 = lastrow + temp_bin;
		      for(Int_t bin=binx1;bin<=binx2;bin++)
			{
			  if(ngaps2[bin] < MAX_N_GAPS) {
			    if(i > (currentrow2[bin] + MAX_GAP_SIZE) && i < lastrow2[bin]) {
			      ngaps2[bin] = MAX_N_GAPS;
			      continue;
			    }
			    if(i > currentrow2[bin] && i < lastrow2[bin])
			      {
				nrows2[bin]++;
				if(i > (currentrow2[bin] + 1))
				  ngaps2[bin]++;
				currentrow2[bin]=i;
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
			}
		    }
		  //End of the transformation

		  }

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
	  UChar_t *lastrow = fTrackLastRow;

	  //Do the transformation:
	  Float_t start_y = eta_clust[eta_index].start_y;
	  if(eta_clust[eta_index].start_pad == 0)
	    start_y -= 0.0*pad_pitch;
	  Float_t R1 = x2 + start_y*start_y; 
	  Float_t x_over_R1 = x/R1;
	  Float_t start_y_over_R1 = start_y/R1;
	  Float_t end_y = (eta_clust[eta_index].end_pad-0.5*(npads-1))*pad_pitch;
	  if(eta_clust[eta_index].end_pad == npads)
	    end_y += 0.0*pad_pitch;
	  Float_t R2 = x2 + end_y*end_y; 
	  Float_t x_over_R2 = x/R2;
	  Float_t end_y_over_R2 = end_y/R2;
	  Float_t A1 = beta1_minus_beta2/(x_over_R1-beta2);
	  Float_t B1 = (x_over_R1-beta1)/(x_over_R1-beta2);
	  Float_t A2 = beta1_minus_beta2/(x_over_R2-beta2);
	  Float_t B2 = (x_over_R2-beta1)/(x_over_R2-beta2);

	  Float_t kappa1 = (A1*start_y_over_R1+B1*y_min-x_min)/x_bin;
	  Float_t delta_kappa1 = B1*hist_bin/x_bin;
	  if(B1<0)
	    kappa1 += delta_kappa1;
	  Float_t kappa2 = (A2*end_y_over_R2+B2*y_min-x_min)/x_bin;
	  Float_t delta_kappa2 = B2*hist_bin/x_bin;
	  if(B2>=0)
	    kappa2 += delta_kappa2;

	  //Fill the histogram along the phirange
	  for(Int_t b=first_bin; b<=last_bin; b++, kappa1 += delta_kappa1, kappa2 += delta_kappa2)
	    {
	      Int_t binx1 = 1 + (Int_t)kappa1;
	      if(binx1>last_binx) continue;
	      if(binx1<first_binx) binx1 = first_binx;
	      Int_t binx2 = 1 + (Int_t)kappa2;
	      if(binx2<first_binx) continue;
	      if(binx2>last_binx) binx2 = last_binx;
#ifdef do_mc
	      if(binx2<binx1) {
		LOG(AliL3Log::kWarning,"AliL3HoughTransformerRow::TransformCircle()","")
		  <<"Wrong filling "<<binx1<<" "<<binx2<<" "<<i<<" "<<x<<" "<<start_y<<" "<<end_y<<ENDLOG;
	      }
#endif
	      Int_t temp_bin = b*nbinx;
	      UChar_t *nrows2 = nrows + temp_bin;
	      UChar_t *ngaps2 = ngaps + temp_bin;
	      UChar_t *currentrow2 = currentrow + temp_bin;
	      UChar_t *lastrow2 = lastrow + temp_bin;
	      for(Int_t bin=binx1;bin<=binx2;bin++)
		{
		  if(ngaps2[bin] < MAX_N_GAPS) {
		    if(i > (currentrow2[bin] + MAX_GAP_SIZE) && i < lastrow2[bin]) {
		      ngaps2[bin] = MAX_N_GAPS;
		      continue;
		    }
		    if(i > currentrow2[bin] && i < lastrow2[bin])
		      {
			nrows2[bin]++;
			if(i > (currentrow2[bin] + 1))
			  ngaps2[bin]++;
			currentrow2[bin]=i;
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
		}
	    }
	  //End of the transformation

	}

      //Move the data pointer to the next padrow:
      AliL3MemHandler::UpdateRowPointer(tempPt);
    }

  delete [] eta_clust;
}

Int_t AliL3HoughTransformerRow::GetTrackID(Int_t eta_index,Double_t kappa,Double_t psi)
{
  if(!fDoMC)
    {
      LOG(AliL3Log::kWarning,"AliL3HoughTransformerRow::GetTrackID","Data")
	<<"Flag switched off"<<ENDLOG;
      return -1;
    }
  
#ifdef do_mc
  if(eta_index < 0 || eta_index > GetNEtaSegments())
    {
      LOG(AliL3Log::kWarning,"AliL3HoughTransformerRow::GetTrackID","Data")
	<<"Wrong etaindex "<<eta_index<<ENDLOG;
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
  LOG(AliL3Log::kWarning,"AliL3HoughTransformerRow::GetTrackID()","")
    <<"Compile with do_mc flag!"<<ENDLOG;
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

UChar_t *AliL3HoughTransformerRow::GetCurrentRowCount(Int_t eta_index)
{
  return fCurrentRowCount[eta_index];
}

UChar_t *AliL3HoughTransformerRow::GetTrackNRows()
{
  return fTrackNRows;
}


UChar_t *AliL3HoughTransformerRow::GetTrackFirstRow()
{
  return fTrackFirstRow;
}

UChar_t *AliL3HoughTransformerRow::GetTrackLastRow()
{
  return fTrackLastRow;
}
