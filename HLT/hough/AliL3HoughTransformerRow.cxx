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

UChar_t **AliL3HoughTransformerRow::fgRowCount = 0;
UChar_t **AliL3HoughTransformerRow::fgGapCount = 0;
UChar_t **AliL3HoughTransformerRow::fgCurrentRowCount = 0;
#ifdef do_mc
TrackIndex **AliL3HoughTransformerRow::fgTrackID = 0;
#endif
UChar_t *AliL3HoughTransformerRow::fgTrackNRows = 0;
UChar_t *AliL3HoughTransformerRow::fgTrackFirstRow = 0;
UChar_t *AliL3HoughTransformerRow::fgTrackLastRow = 0;

Float_t AliL3HoughTransformerRow::fgBeta1 = AliL3Transform::Row2X(79)/pow(AliL3Transform::Row2X(79),2);
Float_t AliL3HoughTransformerRow::fgBeta2 = (AliL3Transform::Row2X(158)+6.0)/pow((AliL3Transform::Row2X(158)+6.0),2);

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

AliL3HoughTransformerRow::AliL3HoughTransformerRow(Int_t slice,Int_t patch,Int_t netasegments,Bool_t /*DoMC*/,Float_t zvertex) : AliL3HoughBaseTransformer(slice,patch,netasegments,zvertex)
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
  //Destructor
  DeleteHistograms();
#ifdef do_mc
  if(fgTrackID)
    {
      for(Int_t i=0; i<GetNEtaSegments(); i++)
	{
	  if(!fgTrackID[i]) continue;
	  delete fgTrackID[i];
	}
      delete [] fgTrackID;
      fgTrackID = 0;
    }
#endif

  if(fgRowCount)
    {
      for(Int_t i=0; i<GetNEtaSegments(); i++)
	{
	  if(!fgRowCount[i]) continue;
	  delete [] fgRowCount[i];
	}
      delete [] fgRowCount;
      fgRowCount = 0;
    }
  if(fgGapCount)
    {
      for(Int_t i=0; i<GetNEtaSegments(); i++)
	{
	  if(!fgGapCount[i]) continue;
	  delete [] fgGapCount[i];
	}
      delete [] fgGapCount;
      fgGapCount = 0;
    }
  if(fgCurrentRowCount)
    {
      for(Int_t i=0; i<GetNEtaSegments(); i++)
	{
	  if(fgCurrentRowCount[i])
	    delete [] fgCurrentRowCount[i];
	}
      delete [] fgCurrentRowCount;
      fgCurrentRowCount = 0;
    }
  if(fgTrackNRows)
    {
      delete [] fgTrackNRows;
      fgTrackNRows = 0;
    }
  if(fgTrackFirstRow)
    {
      delete [] fgTrackFirstRow;
      fgTrackFirstRow = 0;
    }
  if(fgTrackLastRow)
    {
      delete [] fgTrackLastRow;
      fgTrackLastRow = 0;
    }
}

void AliL3HoughTransformerRow::DeleteHistograms()
{
  // Clean up
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

void AliL3HoughTransformerRow::CreateHistograms(Int_t nxbin,Float_t ptmin,
						Int_t nybin,Float_t phimin,Float_t phimax)
{
  //Create the histograms (parameter space).
  //These are 2D histograms, span by kappa (curvature of track) 
  //and phi0 (emission angle with x-axis).
  //The arguments give the range and binning; 
  //nxbin = #bins in kappa
  //nybin = #bins in phi0
  //ptmin = mimium Pt of track (corresponding to maximum kappa)
  //phi_min = mimimum phi0 
  //phi_max = maximum phi0 
    
  Double_t x = AliL3Transform::GetBFact()*AliL3Transform::GetBField()/ptmin;
  //Double_t torad = AliL3Transform::Pi()/180;
  CreateHistograms(nxbin,-1.*x,x,nybin,phimin/**torad*/,phimax/**torad*/);
}

void AliL3HoughTransformerRow::CreateHistograms(Int_t nxbin,Float_t xmin,Float_t xmax,
						Int_t nybin,Float_t ymin,Float_t ymax)
{
  //Create the histograms (parameter space)  
  //nxbin = #bins in X
  //nybin = #bins in Y
  //xmin xmax ymin ymax = histogram limits in X and Y
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
      Int_t ncellsx = (hist->GetNbinsX()+3)/2;
      Int_t ncellsy = (hist->GetNbinsY()+3)/2;
      Int_t ncells = ncellsx*ncellsy;
      if(!fgTrackID)
	{
	  LOG(AliL3Log::kInformational,"AliL3HoughTransformerRow::CreateHistograms()","")
	    <<"Transformer: Allocating "<<GetNEtaSegments()*ncells*sizeof(TrackIndex)<<" bytes to fgTrackID"<<ENDLOG;
	  fgTrackID = new TrackIndex*[GetNEtaSegments()];
	  for(Int_t i=0; i<GetNEtaSegments(); i++)
	    fgTrackID[i] = new TrackIndex[ncells];
	}
    }
#endif
  AliL3Histogram *hist = fParamSpace[0];
  Int_t ncells = (hist->GetNbinsX()+2)*(hist->GetNbinsY()+2);
  if(!fgRowCount)
    {
      LOG(AliL3Log::kInformational,"AliL3HoughTransformerRow::CreateHistograms()","")
	<<"Transformer: Allocating "<<GetNEtaSegments()*ncells*sizeof(UChar_t)<<" bytes to fgRowCount"<<ENDLOG;
      fgRowCount = new UChar_t*[GetNEtaSegments()];
      for(Int_t i=0; i<GetNEtaSegments(); i++)
	fgRowCount[i] = new UChar_t[ncells];
    }
  if(!fgGapCount)
    {
      LOG(AliL3Log::kInformational,"AliL3HoughTransformerRow::CreateHistograms()","")
	<<"Transformer: Allocating "<<GetNEtaSegments()*ncells*sizeof(UChar_t)<<" bytes to fgGapCount"<<ENDLOG;
      fgGapCount = new UChar_t*[GetNEtaSegments()];
      for(Int_t i=0; i<GetNEtaSegments(); i++)
	fgGapCount[i] = new UChar_t[ncells];
    }
  if(!fgCurrentRowCount)
    {
      LOG(AliL3Log::kInformational,"AliL3HoughTransformerRow::CreateHistograms()","")
	<<"Transformer: Allocating "<<GetNEtaSegments()*ncells*sizeof(UChar_t)<<" bytes to fgCurrentRowCount"<<ENDLOG;
      fgCurrentRowCount = new UChar_t*[GetNEtaSegments()];
      for(Int_t i=0; i<GetNEtaSegments(); i++)
	fgCurrentRowCount[i] = new UChar_t[ncells];
    }

  if(!fgTrackNRows)
    {
      LOG(AliL3Log::kInformational,"AliL3HoughTransformerRow::CreateHistograms()","")
	<<"Transformer: Allocating "<<ncells*sizeof(UChar_t)<<" bytes to fgTrackNRows"<<ENDLOG;
      fgTrackNRows = new UChar_t[ncells];
      LOG(AliL3Log::kInformational,"AliL3HoughTransformerRow::CreateHistograms()","")
	<<"Transformer: Allocating "<<ncells*sizeof(UChar_t)<<" bytes to fgTrackFirstRow"<<ENDLOG;
      fgTrackFirstRow = new UChar_t[ncells];
      LOG(AliL3Log::kInformational,"AliL3HoughTransformerRow::CreateHistograms()","")
	<<"Transformer: Allocating "<<ncells*sizeof(UChar_t)<<" bytes to fgTrackLastRow"<<ENDLOG;
      fgTrackLastRow = new UChar_t[ncells];

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
	      fgTrackNRows[xbin + ybin*nxbins] = 99999;
	      for(Int_t deltay = -1; deltay <= 1; deltay += 2) {
		for(Int_t deltax = -1; deltax <= 1; deltax += 2) {
		
		  Float_t xtrack = hist->GetPreciseBinCenterX((Float_t)xbin+0.5*(Float_t)deltax);
		  Float_t ytrack = hist->GetPreciseBinCenterY((Float_t)ybin+0.5*(Float_t)deltay);

		  Float_t psi = atan((xtrack-ytrack)/(fgBeta1-fgBeta2));
		  Float_t kappa = 2.0*(xtrack*cos(psi)-fgBeta1*sin(psi));
		  track.SetTrackParameters(kappa,psi,1);
		  Bool_t firstrow = kFALSE;
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
			  if(!firstrow) {
			    curfirstrow = j;
			    firstrow = kTRUE;
			  }
			  curlastrow = j;
			}
		      else {
			if(firstrow) {
			  firstrow = kFALSE;
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
		  if((maxlastrow-maxfirstrow) < fgTrackNRows[xbin + ybin*nxbins]) {
		    fgTrackNRows[xbin + ybin*nxbins] = maxlastrow-maxfirstrow;
		    fgTrackFirstRow[xbin + ybin*nxbins] = maxfirstrow;
		    fgTrackLastRow[xbin + ybin*nxbins] = maxlastrow;
		  }
		}
	      }
	      //	      cout<<" fgTrackNRows "<<xbin<<" "<<ybin<<" "<<(Int_t)fgTrackNRows[xbin + ybin*nxbins]<<" "<<endl;
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
  //Reset all the histograms. Should be done when processing new slice

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
      Int_t ncellsx = (hist->GetNbinsX()+3)/2;
      Int_t ncellsy = (hist->GetNbinsY()+3)/2;
      Int_t ncells = ncellsx*ncellsy;
      for(Int_t i=0; i<GetNEtaSegments(); i++)
	memset(fgTrackID[i],0,ncells*sizeof(TrackIndex));
    }
#endif
  AliL3Histogram *hist = fParamSpace[0];
  Int_t ncells = (hist->GetNbinsX()+2)*(hist->GetNbinsY()+2);
  for(Int_t i=0; i<GetNEtaSegments(); i++)
    {
      memset(fgRowCount[i],0,ncells*sizeof(UChar_t));
      memset(fgGapCount[i],1,ncells*sizeof(UChar_t));
      //      memset(fgCurrentRowCount[i],160,ncells*sizeof(UChar_t));
      memcpy(fgCurrentRowCount[i],fgTrackFirstRow,ncells*sizeof(UChar_t));
    }
}

Int_t AliL3HoughTransformerRow::GetEtaIndex(Double_t eta)
{
  //Return the histogram index of the corresponding eta. 

  Double_t etaslice = (GetEtaMax() - GetEtaMin())/GetNEtaSegments();
  Double_t index = (eta-GetEtaMin())/etaslice;
  return (Int_t)index;
}

inline AliL3Histogram *AliL3HoughTransformerRow::GetHistogram(Int_t etaindex)
{
  // Return a pointer to the histogram which contains etaindex eta slice
  if(!fParamSpace || etaindex >= GetNEtaSegments() || etaindex < 0)
    return 0;
  if(!fParamSpace[etaindex])
    return 0;
  return fParamSpace[etaindex];
}

Double_t AliL3HoughTransformerRow::GetEta(Int_t etaindex,Int_t /*slice*/)
{
  // Return eta calculated in the middle of the eta slice
  Double_t etaslice = (GetEtaMax()-GetEtaMin())/GetNEtaSegments();
  Double_t eta=0;
  eta=(Double_t)((etaindex+0.5)*etaslice);
  return eta;
}

struct AliL3EtaRow {
  UChar_t fStartPad; //First pad in the cluster
  UChar_t fEndPad; //Last pad in the cluster
  Bool_t fIsFound; //Is the cluster already found
  Float_t fStartY; //Y position of the first pad in the cluster
#ifdef do_mc
  Int_t fMcLabels[MaxTrack]; //Array to store mc labels inside cluster
#endif
};

void AliL3HoughTransformerRow::TransformCircle()
{
  //Do the Hough Transform
  Float_t beta1 = fgBeta1;
  Float_t beta2 = fgBeta2;
  Float_t beta1minusbeta2 = fgBeta1 - fgBeta2;

  Int_t netasegments = GetNEtaSegments();
  Double_t etamin = GetEtaMin();
  Double_t etaslice = (GetEtaMax() - etamin)/netasegments;

  Int_t lowerthreshold = GetLowerThreshold();

  //Assumes that all the etaslice histos are the same!
  AliL3Histogram *h = fParamSpace[0];
  Float_t ymin = h->GetYmin();
  //Float_t y_max = h->GetYmax();
  Float_t histbin = h->GetBinWidthY();
  Int_t firstbin = h->GetFirstYbin();
  Int_t lastbin = h->GetLastYbin();
  Float_t xmin = h->GetXmin();
  Float_t xmax = h->GetXmax();
  Float_t xbin = (xmax-xmin)/h->GetNbinsX();
  Int_t firstbinx = h->GetFirstXbin()-1;
  Int_t lastbinx = h->GetLastXbin()+1;
  Int_t nbinx = h->GetNbinsX()+2;

  UChar_t lastpad;
  Int_t lastetaindex;
  AliL3EtaRow *etaclust = new AliL3EtaRow[netasegments];

  AliL3DigitRowData *tempPt = GetDataPointer();
  if(!tempPt)
    {
      LOG(AliL3Log::kError,"AliL3HoughTransformer::TransformCircle","Data")
	<<"No input data "<<ENDLOG;
      return;
    }

  Int_t ipatch = GetPatch();
  Double_t padpitch = AliL3Transform::GetPadPitchWidth(ipatch);
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

      lastpad = 255;
      //Flush eta clusters array
      memset(etaclust,0,netasegments*sizeof(AliL3EtaRow));  

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
	  if((Int_t)charge <= lowerthreshold)
	    continue;
	  UChar_t pad = digPt[j].fPad;
	  UShort_t time = digPt[j].fTime;

	  if(pad != lastpad)
	    {
	      y = (pad-0.5*npads)*padpitch;
	      Float_t y2 = y*y;
	      r2 = x2 + y2;
	      lastetaindex = -1;
	    }

	  //Transform data to local cartesian coordinates:
	  Float_t z = fLUTz[(Int_t)time];
	  Float_t z2 = fLUTz2[(Int_t)time];
	  //Calculate the eta:
	  Double_t r = sqrt(r2+z2);
	  Double_t eta = 0.5 * log((r+z)/(r-z));
	  //Get the corresponding index, which determines which histogram to fill:
	  Int_t etaindex = (Int_t)((eta-etamin)/etaslice);

#ifndef do_mc
	  if(etaindex == lastetaindex) continue;
#endif
	  //	  cout<<" Digit at patch "<<ipatch<<" row "<<i<<" pad "<<(Int_t)pad<<" time "<<time<<" etaslice "<<etaindex<<endl;
	  
	  if(etaindex < 0 || etaindex >= netasegments)
	    continue;

	  if(!etaclust[etaindex].fIsFound)
	    {
	      etaclust[etaindex].fStartPad = pad;
	      etaclust[etaindex].fEndPad = pad;
	      etaclust[etaindex].fStartY = y - padpitch/2.0;
	      etaclust[etaindex].fIsFound = 1;
#ifdef do_mc
	      if(fDoMC)
		{
		  for(Int_t t=0; t<3; t++)
		    {
		      Int_t label = digPt[j].fTrackID[t];
		      if(label < 0) break;
		      UInt_t c;
		      for(c=0; c<(MaxTrack-1); c++)
			if(etaclust[etaindex].fMcLabels[c] == label || etaclust[etaindex].fMcLabels[c] == 0)
			  break;
		      //		      if(c == MaxTrack-1) cerr<<"AliL3HoughTransformer::TransformCircle : Cluster array reached maximum!! "<<c<<endl;
		      etaclust[etaindex].fMcLabels[c] = label;
		    }
		}
#endif
	      continue;
	    }
	  else
	    {
	      if(pad <= (etaclust[etaindex].fEndPad + 1))
		{
		  etaclust[etaindex].fEndPad = pad;
#ifdef do_mc
		  if(fDoMC)
		    {
		      for(Int_t t=0; t<3; t++)
			{
			  Int_t label = digPt[j].fTrackID[t];
			  if(label < 0) break;
			  UInt_t c;
			  for(c=0; c<(MaxTrack-1); c++)
			    if(etaclust[etaindex].fMcLabels[c] == label || etaclust[etaindex].fMcLabels[c] == 0)
			      break;
			  //			  if(c == MaxTrack-1) cerr<<"AliL3HoughTransformer::TransformCircle : Cluster array reached maximum!! "<<c<<endl;
			  etaclust[etaindex].fMcLabels[c] = label;
			}
		    }
#endif
		}
	      else
		{
		  Bool_t fillcluster = kTRUE;
		  if(fillcluster) {

		  UChar_t *nrows = fgRowCount[etaindex];
		  UChar_t *ngaps = fgGapCount[etaindex];
		  UChar_t *currentrow = fgCurrentRowCount[etaindex];
		  UChar_t *lastrow = fgTrackLastRow;

		  //Do the transformation:
		  Float_t starty = etaclust[etaindex].fStartY;
		  if(etaclust[etaindex].fStartPad == 0)
		    starty -= 0.0*padpitch;
		  Float_t r1 = x2 + starty*starty;
		  Float_t xoverr1 = x/r1;
		  Float_t startyoverr1 = starty/r1;
		  Float_t endy = (etaclust[etaindex].fEndPad-0.5*(npads-1))*padpitch;
		  if(etaclust[etaindex].fEndPad == npads)
		    endy += 0.0*padpitch;
		  Float_t r2 = x2 + endy*endy; 
		  Float_t xoverr2 = x/r2;
		  Float_t endyoverr2 = endy/r2;
		  Float_t a1 = beta1minusbeta2/(xoverr1-beta2);
		  Float_t b1 = (xoverr1-beta1)/(xoverr1-beta2);
		  Float_t a2 = beta1minusbeta2/(xoverr2-beta2);
		  Float_t b2 = (xoverr2-beta1)/(xoverr2-beta2);

		  Float_t kappa1 = (a1*startyoverr1+b1*ymin-xmin)/xbin;
		  Float_t deltaalpha1 = b1*histbin/xbin;
		  if(b1<0)
		    kappa1 += deltaalpha1;
		  Float_t kappa2 = (a2*endyoverr2+b2*ymin-xmin)/xbin;
		  Float_t deltaalpha2 = b2*histbin/xbin;
		  if(b2>=0)
		    kappa2 += deltaalpha2;

		  //Fill the histogram along the phirange
		  for(Int_t b=firstbin; b<=lastbin; b++, kappa1 += deltaalpha1, kappa2 += deltaalpha2)
		    {
		      Int_t binx1 = 1 + (Int_t)kappa1;
		      if(binx1>lastbinx) continue;
		      if(binx1<firstbinx) binx1 = firstbinx;
		      Int_t binx2 = 1 + (Int_t)kappa2;
		      if(binx2<firstbinx) continue;
		      if(binx2>lastbinx) binx2 = lastbinx;
#ifdef do_mc
		      if(binx2<binx1) {
			LOG(AliL3Log::kWarning,"AliL3HoughTransformerRow::TransformCircle()","")
			  <<"Wrong filling "<<binx1<<" "<<binx2<<" "<<i<<" "<<x<<" "<<starty<<" "<<endy<<ENDLOG;
		      }
#endif
		      Int_t tempbin = b*nbinx;
		      UChar_t *nrows2 = nrows + tempbin;
		      UChar_t *ngaps2 = ngaps + tempbin;
		      UChar_t *currentrow2 = currentrow + tempbin;
		      UChar_t *lastrow2 = lastrow + tempbin;
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
				    Int_t label = etaclust[etaindex].fMcLabels[t];
				    if(label == 0) break;
				    UInt_t c;
				    Int_t tempbin2 = ((Int_t)(b/2))*((Int_t)((nbinx+1)/2)) + (Int_t)(bin/2);
				    for(c=0; c<(MaxTrack-1); c++)
				      if(fgTrackID[etaindex][tempbin2].fLabel[c] == label || fgTrackID[etaindex][tempbin2].fNHits[c] == 0)
					break;
				    //				    if(c == MaxTrack-1) cerr<<"AliL3HoughTransformer::TransformCircle : Array reached maximum!! "<<c<<endl;
				    fgTrackID[etaindex][tempbin2].fLabel[c] = label;
				    if(fgTrackID[etaindex][tempbin2].fCurrentRow[c] != i) {
				      fgTrackID[etaindex][tempbin2].fNHits[c]++;
				      fgTrackID[etaindex][tempbin2].fCurrentRow[c] = i;
				    }
				  }
			      }
#endif
			  }
			}
		    }
		  //End of the transformation

		  }

		  etaclust[etaindex].fStartPad = pad;
		  etaclust[etaindex].fEndPad = pad;
		  etaclust[etaindex].fStartY = y - padpitch/2.0;

#ifdef do_mc
		  if(fDoMC)
		    {
		      memset(etaclust[etaindex].fMcLabels,0,MaxTrack);
		      for(Int_t t=0; t<3; t++)
			{
			  Int_t label = digPt[j].fTrackID[t];
			  if(label < 0) break;
			  UInt_t c;
			  for(c=0; c<(MaxTrack-1); c++)
			    if(etaclust[etaindex].fMcLabels[c] == label || etaclust[etaindex].fMcLabels[c] == 0)
			      break;
			  //			  if(c == MaxTrack-1) cerr<<"AliL3HoughTransformer::TransformCircle : Cluster array reached maximum!! "<<c<<endl;
			  etaclust[etaindex].fMcLabels[c] = label;
			}
		    }
#endif
		}
	    }
	  lastpad = pad;
	  lastetaindex = etaindex;
	}
      //Write remaining clusters
      for(Int_t etaindex = 0;etaindex < netasegments;etaindex++)
	{
	  //Check for empty row
	  if((etaclust[etaindex].fStartPad == 0) && (etaclust[etaindex].fEndPad == 0)) continue;

	  UChar_t *nrows = fgRowCount[etaindex];
	  UChar_t *ngaps = fgGapCount[etaindex];
	  UChar_t *currentrow = fgCurrentRowCount[etaindex];
	  UChar_t *lastrow = fgTrackLastRow;

	  //Do the transformation:
	  Float_t starty = etaclust[etaindex].fStartY;
	  if(etaclust[etaindex].fStartPad == 0)
	    starty -= 0.0*padpitch;
	  Float_t r1 = x2 + starty*starty; 
	  Float_t xoverr1 = x/r1;
	  Float_t startyoverr1 = starty/r1;
	  Float_t endy = (etaclust[etaindex].fEndPad-0.5*(npads-1))*padpitch;
	  if(etaclust[etaindex].fEndPad == npads)
	    endy += 0.0*padpitch;
	  Float_t r2 = x2 + endy*endy; 
	  Float_t xoverr2 = x/r2;
	  Float_t endyoverr2 = endy/r2;
	  Float_t a1 = beta1minusbeta2/(xoverr1-beta2);
	  Float_t b1 = (xoverr1-beta1)/(xoverr1-beta2);
	  Float_t a2 = beta1minusbeta2/(xoverr2-beta2);
	  Float_t b2 = (xoverr2-beta1)/(xoverr2-beta2);

	  Float_t kappa1 = (a1*startyoverr1+b1*ymin-xmin)/xbin;
	  Float_t deltaalpha1 = b1*histbin/xbin;
	  if(b1<0)
	    kappa1 += deltaalpha1;
	  Float_t kappa2 = (a2*endyoverr2+b2*ymin-xmin)/xbin;
	  Float_t deltaalpha2 = b2*histbin/xbin;
	  if(b2>=0)
	    kappa2 += deltaalpha2;

	  //Fill the histogram along the phirange
	  for(Int_t b=firstbin; b<=lastbin; b++, kappa1 += deltaalpha1, kappa2 += deltaalpha2)
	    {
	      Int_t binx1 = 1 + (Int_t)kappa1;
	      if(binx1>lastbinx) continue;
	      if(binx1<firstbinx) binx1 = firstbinx;
	      Int_t binx2 = 1 + (Int_t)kappa2;
	      if(binx2<firstbinx) continue;
	      if(binx2>lastbinx) binx2 = lastbinx;
#ifdef do_mc
	      if(binx2<binx1) {
		LOG(AliL3Log::kWarning,"AliL3HoughTransformerRow::TransformCircle()","")
		  <<"Wrong filling "<<binx1<<" "<<binx2<<" "<<i<<" "<<x<<" "<<starty<<" "<<endy<<ENDLOG;
	      }
#endif
	      Int_t tempbin = b*nbinx;
	      UChar_t *nrows2 = nrows + tempbin;
	      UChar_t *ngaps2 = ngaps + tempbin;
	      UChar_t *currentrow2 = currentrow + tempbin;
	      UChar_t *lastrow2 = lastrow + tempbin;
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
			    Int_t label = etaclust[etaindex].fMcLabels[t];
			    if(label == 0) break;
			    UInt_t c;
			    Int_t tempbin2 = ((Int_t)(b/2))*((Int_t)((nbinx+1)/2)) + (Int_t)(bin/2);
			    for(c=0; c<(MaxTrack-1); c++)
			      if(fgTrackID[etaindex][tempbin2].fLabel[c] == label || fgTrackID[etaindex][tempbin2].fNHits[c] == 0)
				break;
			    //			    if(c == MaxTrack-1) cout<<"AliL3HoughTransformer::TransformCircle : Array reached maximum!! "<<c<<endl;
			    fgTrackID[etaindex][tempbin2].fLabel[c] = label;
			    if(fgTrackID[etaindex][tempbin2].fCurrentRow[c] != i) {
			      fgTrackID[etaindex][tempbin2].fNHits[c]++;
			      fgTrackID[etaindex][tempbin2].fCurrentRow[c] = i;
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

  delete [] etaclust;
}

Int_t AliL3HoughTransformerRow::GetTrackID(Int_t etaindex,Double_t kappa,Double_t psi)
{
  // Returns the MC label for a given peak found in the Hough space
  if(!fDoMC)
    {
      LOG(AliL3Log::kWarning,"AliL3HoughTransformerRow::GetTrackID","Data")
	<<"Flag switched off"<<ENDLOG;
      return -1;
    }
  
#ifdef do_mc
  if(etaindex < 0 || etaindex > GetNEtaSegments())
    {
      LOG(AliL3Log::kWarning,"AliL3HoughTransformerRow::GetTrackID","Data")
	<<"Wrong etaindex "<<etaindex<<ENDLOG;
      return -1;
    }
  AliL3Histogram *hist = fParamSpace[etaindex];
  Int_t bin = hist->FindLabelBin(kappa,psi);
  if(bin==-1) {
    LOG(AliL3Log::kWarning,"AliL3HoughTransformerRow::GetTrackID()","")
      <<"Track candidate outside Hough space boundaries: "<<kappa<<" "<<psi<<ENDLOG;
    return -1;
  }
  Int_t label=-1;
  Int_t max=0;
  for(UInt_t i=0; i<(MaxTrack-1); i++)
    {
      Int_t nhits=fgTrackID[etaindex][bin].fNHits[i];
      if(nhits == 0) break;
      if(nhits > max)
	{
	  max = nhits;
	  label = fgTrackID[etaindex][bin].fLabel[i];
	}
    }
  Int_t label2=-1;
  Int_t max2=0;
  for(UInt_t i=0; i<(MaxTrack-1); i++)
    {
      Int_t nhits=fgTrackID[etaindex][bin].fNHits[i];
      if(nhits == 0) break;
      if(nhits > max2)
	{
	  if(fgTrackID[etaindex][bin].fLabel[i]!=label) {
	    max2 = nhits;
	    label2 = fgTrackID[etaindex][bin].fLabel[i];
	  }
	}
    }
  if(max2 !=0 ) {
    LOG(AliL3Log::kDebug,"AliL3HoughTransformerRow::GetTrackID()","")
      <<" TrackID"<<" label "<<label<<" max "<<max<<" label2 "<<label2<<" max2 "<<max2<<" "<<(Float_t)max2/(Float_t)max<<" "<<fgTrackID[etaindex][bin].fLabel[MaxTrack-1]<<" "<<(Int_t)fgTrackID[etaindex][bin].fNHits[MaxTrack-1]<<ENDLOG;
  }
  return label;
#endif
  LOG(AliL3Log::kWarning,"AliL3HoughTransformerRow::GetTrackID()","")
    <<"Compile with do_mc flag!"<<ENDLOG;
  return -1;
}

UChar_t *AliL3HoughTransformerRow::GetRowCount(Int_t etaindex)
{
  return fgRowCount[etaindex];
}

UChar_t *AliL3HoughTransformerRow::GetGapCount(Int_t etaindex)
{
  return fgGapCount[etaindex];
}

UChar_t *AliL3HoughTransformerRow::GetCurrentRowCount(Int_t etaindex)
{
  return fgCurrentRowCount[etaindex];
}

UChar_t *AliL3HoughTransformerRow::GetTrackNRows()
{
  return fgTrackNRows;
}


UChar_t *AliL3HoughTransformerRow::GetTrackFirstRow()
{
  return fgTrackFirstRow;
}

UChar_t *AliL3HoughTransformerRow::GetTrackLastRow()
{
  return fgTrackLastRow;
}
