// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliL3StandardIncludes.h"

#include "AliL3Logging.h"
#include "AliL3MemHandler.h"
#include "AliL3Transform.h"
#include "AliL3DigitData.h"
#include "AliL3HistogramAdaptive.h"
#include "AliL3HoughTransformerGap.h"

#if GCCVERSION == 3
using namespace std;
#endif

//_____________________________________________________________
// AliL3HoughTransformerGap
//
// Hough transformation class
//

ClassImp(AliL3HoughTransformerGap)

#ifdef COUNTINGROWS
UChar_t **AliL3HoughTransformerGap::fRowCount = 0;
UChar_t **AliL3HoughTransformerGap::fGapCount = 0;
UChar_t **AliL3HoughTransformerGap::fCurrentRowCount = 0;
#endif

AliL3HoughTransformerGap::AliL3HoughTransformerGap()
{
  //Default constructor
  fParamSpace = 0;
  fDoMC = kFALSE;;

#ifdef do_mc
  fTrackID = 0;
#endif
  fZVertex = 0.;
#ifdef SUBSLICES
  fEtaSpace = 0;
  fEtaSpaceSize = 0;
#endif
#ifdef COUNTINGROWS
  fLUT2sinphi0up=0;
  fLUT2sinphi0low=0;
  fLUT2cosphi0up=0;
  fLUT2cosphi0low=0;
#endif
}

AliL3HoughTransformerGap::AliL3HoughTransformerGap(Int_t slice,Int_t patch,Int_t n_eta_segments,Bool_t DoEtaOverlap,Bool_t DoMC,Float_t zvertex) : AliL3HoughBaseTransformer(slice,patch,n_eta_segments)
{
  //Normal constructor
  fParamSpace = 0;
  fDoMC = kFALSE;
  fDoMC=kFALSE;

#ifdef do_mc
  fTrackID = 0;
#endif
  fZVertex = zvertex;
#ifdef SUBSLICES
  fEtaSpace = 0;
  fEtaSpaceSize = 800/n_eta_segments;
#endif
#ifdef COUNTINGROWS
  fLUT2sinphi0up=0;
  fLUT2sinphi0low=0;
  fLUT2cosphi0up=0;
  fLUT2cosphi0low=0;
#endif
}

AliL3HoughTransformerGap::~AliL3HoughTransformerGap()
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
    }
#endif
#ifdef SUBSLICES
  if(fEtaSpace)
    {
      for(Int_t i=0; i<GetNEtaSegments(); i++)
	{
	  if(!fEtaSpace[i]) continue;
	  delete fEtaSpace[i];
	}
      delete [] fEtaSpace;
    }
#endif
#ifdef COUNTINGROWS
  if(fRowCount)
    {
      for(Int_t i=0; i<GetNEtaSegments(); i++)
	{
	  if(!fRowCount[i]) continue;
	  delete fRowCount[i];
	}
      delete [] fRowCount;
    }
  if(fGapCount)
    {
      for(Int_t i=0; i<GetNEtaSegments(); i++)
	{
	  if(!fGapCount[i]) continue;
	  delete fGapCount[i];
	}
      delete [] fGapCount;
    }
  if(fCurrentRowCount)
    {
      for(Int_t i=0; i<GetNEtaSegments(); i++)
	{
	  if(fCurrentRowCount[i])
	    delete fCurrentRowCount[i];
	}
      delete [] fCurrentRowCount;
    }
#endif
}

void AliL3HoughTransformerGap::DeleteHistograms()
{
#ifdef COUNTINGROWS
  delete[] fLUT2sinphi0up;
  delete[] fLUT2cosphi0up;
  delete[] fLUT2sinphi0low;
  delete[] fLUT2cosphi0low;

  fLUT2sinphi0up=0;
  fLUT2cosphi0up=0;
  fLUT2sinphi0low=0;
  fLUT2cosphi0low=0;
#endif
  if(!fParamSpace)
    return;
  for(Int_t i=0; i<GetNEtaSegments(); i++)
    {
      if(!fParamSpace[i]) continue;
      delete fParamSpace[i];
    }
  delete [] fParamSpace;
}

void AliL3HoughTransformerGap::CreateHistograms(Int_t nxbin,Float_t pt_min,
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

void AliL3HoughTransformerGap::CreateHistograms(Int_t nxbin,Float_t xmin,Float_t xmax,
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
      cout<<"Transformer: Allocating "<<GetNEtaSegments()*ncells*sizeof(TrackIndex)<<" bytes to fTrackID"<<endl;
      fTrackID = new TrackIndex*[GetNEtaSegments()];
      for(Int_t i=0; i<GetNEtaSegments(); i++)
	fTrackID[i] = new TrackIndex[ncells];
    }
#endif
#ifdef SUBSLICES
  AliL3Histogram *hist = fParamSpace[0];
  Int_t ncells = (hist->GetNbinsX()+2)*(hist->GetNbinsY()+2);
  cout<<"Transformer: Allocating "<<GetNEtaSegments()*ncells*fEtaSpaceSize*sizeof(Int_t)<<" bytes to fEtaSpace"<<endl;
  fEtaSpace = new Int_t*[GetNEtaSegments()];
  for(Int_t i=0; i<GetNEtaSegments(); i++)
    fEtaSpace[i] = new Int_t[ncells*fEtaSpaceSize];
#endif
#ifdef COUNTINGROWS
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
#endif

}

void AliL3HoughTransformerGap::Reset()
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
#ifdef SUBSLICES
  AliL3Histogram *hist = fParamSpace[0];
  Int_t ncells = (hist->GetNbinsX()+2)*(hist->GetNbinsY()+2);
  for(Int_t i=0; i<GetNEtaSegments(); i++)
    memset(fEtaSpace[i],0,ncells*fEtaSpaceSize*sizeof(Int_t));
#endif
#ifdef COUNTINGROWS
  AliL3Histogram *hist = fParamSpace[0];
  Int_t ncells = (hist->GetNbinsX()+2)*(hist->GetNbinsY()+2);
  for(Int_t i=0; i<GetNEtaSegments(); i++)
    {
      memset(fRowCount[i],0,ncells*sizeof(UChar_t));
      memset(fGapCount[i],1,ncells*sizeof(UChar_t));
      memset(fCurrentRowCount[i],160,ncells*sizeof(UChar_t));
    }
#endif
}

Int_t AliL3HoughTransformerGap::GetEtaIndex(Double_t eta)
{
  //Return the histogram index of the corresponding eta. 

  Double_t etaslice = (GetEtaMax() - GetEtaMin())/GetNEtaSegments();
  Double_t index = (eta-GetEtaMin())/etaslice;
  return (Int_t)index;
}

inline AliL3Histogram *AliL3HoughTransformerGap::GetHistogram(Int_t eta_index)
{
  if(!fParamSpace || eta_index >= GetNEtaSegments() || eta_index < 0)
    return 0;
  if(!fParamSpace[eta_index])
    return 0;
  return fParamSpace[eta_index];
}

Double_t AliL3HoughTransformerGap::GetEta(Int_t eta_index,Int_t slice)
{
  Double_t eta_slice = (GetEtaMax()-GetEtaMin())/GetNEtaSegments();
  Double_t eta=0;
  eta=(Double_t)((eta_index+0.5)*eta_slice);
  return eta;
}

#ifdef SUBSLICES
Int_t AliL3HoughTransformerGap::GetEtaSubIndex(Double_t eta,Int_t eta_index)
{
  //Return the eta space sub index of the corresponding eta. 

  Double_t sub_eta = GetEta(eta_index,GetSlice());
  Int_t eta_subindex = (Int_t)((eta-sub_eta+0.5/(Double_t)GetNEtaSegments())*800.);
  return eta_subindex;
}

Double_t AliL3HoughTransformerGap::GetMaxSubEta(Float_t kappa,Float_t phi0,Int_t eta_index,Int_t slice)
{
  AliL3Histogram *hist = fParamSpace[eta_index];
  Int_t bin = hist->FindBin(kappa,phi0);
  Int_t maxbin=fEtaSpaceSize/2;
  Int_t maxvalue=0;
  for(Int_t i=0;i<fEtaSpaceSize;i++) {
    if(fEtaSpace[eta_index][bin*fEtaSpaceSize+i]>maxvalue) {
      maxbin=i;
      maxvalue=fEtaSpace[eta_index][bin*fEtaSpaceSize+i];
    }
  }
  Double_t eta = GetEta(eta_index,slice)+((Double_t)maxbin-0.5*(Double_t)fEtaSpaceSize+0.5)/800.;
  cout<<"DEBUG maxetasubindex "<<eta<<" "<<eta_index<<" "<<maxbin<<endl;
  return eta;
}
#endif

#ifdef notused
void AliL3HoughTransformerGap::TransformCircle()
{

  AliL3DigitRowData *tempPt = GetDataPointer();
  if(!tempPt)
    {
      LOG(AliL3Log::kError,"AliL3HoughTransformer::TransformCircle","Data")
	<<"No input data "<<ENDLOG;
      return;
    }

  Int_t ipatch = GetPatch();
  //Loop over the padrows:
  for(Int_t i=AliL3Transform::GetFirstRow(GetPatch()); i<=AliL3Transform::GetLastRow(GetPatch()); i++)
    {
      //Get the data on this padrow:
      AliL3DigitData *digPt = tempPt->fDigitData;
      if(i != (Int_t)tempPt->fRow)
	{
	  cerr<<"AliL3HoughTransform::TransformCircle : Mismatching padrow numbering "<<i<<" "<<(Int_t)tempPt->fRow<<endl;
	  continue;
	}
      
      //Loop over the data on this padrow:
      for(UInt_t j=0; j<tempPt->fNDigit; j++)
	{
	  UShort_t charge = digPt[j].fCharge;
	  UChar_t pad = digPt[j].fPad;
	  UShort_t time = digPt[j].fTime;
	  if((Int_t)charge <= GetLowerThreshold())
	    continue;
	  
	  if((Int_t)charge > GetUpperThreshold())
	    charge = GetUpperThreshold();
	  
	  Int_t sector,row;
	  Float_t xyz[3];

	  //Transform data to local cartesian coordinates:
	  AliL3Transform::Slice2Sector(GetSlice(),i,sector,row);
	  AliL3Transform::Raw2Local(xyz,sector,row,(Int_t)pad,(Int_t)time);
	  	  
	  //Calculate the eta:
	  xyz[2] -= fZVertex;
	  Double_t eta = AliL3Transform::GetEta(xyz);
	  
	  //Get the corresponding index, which determines which histogram to fill:
	  Int_t eta_index = GetEtaIndex(eta);
	  	  
	  if(eta_index < 0 || eta_index >= GetNEtaSegments())
	    continue;
	  
	  //Get the correct histogrampointer:
	  AliL3Histogram *hist = fParamSpace[eta_index];
	  if(!hist)
	    {
	      cerr<<"AliL3HoughTransformer::TransformCircle : Error getting histogram in index "<<eta_index<<endl;
	      continue;
	    }
	  
	  //Do the transformation:
#ifndef COUNTINGROWS
	  Float_t R = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]); 
	  Float_t phi = AliL3Transform::GetPhi(xyz);
#else
	  xyz[1] -= AliL3Transform::GetPadPitchWidth(ipatch)/2.0;
	  Float_t R1 = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]); 
	  Float_t phi1 = AliL3Transform::GetPhi(xyz);
	  xyz[1] += AliL3Transform::GetPadPitchWidth(ipatch);
	  Float_t R2 = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]); 
	  Float_t phi2 = AliL3Transform::GetPhi(xyz);
#endif
#ifdef SUBSLICES	  
	  Int_t eta_subindex = GetEtaSubIndex(eta,eta_index);
#endif	  
	  //Fill the histogram along the phirange
	  for(Int_t b=hist->GetFirstYbin(); b<=hist->GetLastYbin(); b++)
	    {
	      Float_t phi0 = hist->GetBinCenterY(b);
#ifndef COUNTINGROWS
	      Float_t kappa = 2*sin(phi - phi0)/R;
#else
	      Float_t kappa1 = 2*sin(phi1 - phi0 - hist->GetBinWidthY()/2.)/R1;
	      Float_t kappa2 = 2*sin(phi2 - phi0 + hist->GetBinWidthY()/2.)/R2;
#endif
	      //hist->Fill(kappa,phi0,(int)rint(log((Float_t)charge)));
	      //	      hist->Fill(kappa,phi0,charge);
#ifndef COUNTINGROWS
	      hist->Fill(kappa,phi0,30);
	      //hist->Fill(kappa,phi0,1);
#else
	      Int_t binx1,binx2;
	      binx1 = hist->FindXbin(kappa1);
	      binx2 = hist->FindXbin(kappa2);

	      //	      cout<<"DEBUG counting rows "<<ipatch<<" "<<i<<" "<<" bins "<<binx1<<" "<<binx2<<" "<<" kappas "<<kappa1<<" "<<kappa2<<endl;

	      Int_t bin1 = hist->GetBin(binx1,b);
	      Int_t bin2 = hist->GetBin(binx2,b);
	      if((bin1==0) || (bin2==0)) continue;
	      //	      cout<<"DEBUG counting rows "<<ipatch<<" "<<i<<" "<<" bins "<<binx1<<" "<<binx2<<" "<<bin1<<" "<<bin2<<" kappas "<<kappa1<<" "<<kappa2<<endl;
	      for(Int_t bin=bin1;bin<=bin2;bin++) {
		fRowCount[eta_index][bin]=fRowCount[eta_index][bin]|(1<<(i-AliL3Transform::GetFirstRow(ipatch)));
	      }
#endif
#ifndef COUNTINGROWS
#ifdef SUBSLICES
	      //Fill Eta space histogram
	      Int_t bin = hist->FindBin(kappa,phi0);
	      //	      cout<<" DEBUG etasubindex "<<eta<<" "<<eta_index<<" "<<eta_subindex<<endl;
	      fEtaSpace[eta_index][bin*fEtaSpaceSize+eta_subindex]++;
#endif
#ifdef do_mc
	      if(fDoMC)
		{
		  Int_t bin = hist->FindBin(kappa,phi0);
		  for(Int_t t=0; t<3; t++)
		    {
		      Int_t label = digPt[j].fTrackID[t];
		      if(label < 0) break;
		      UInt_t c;
		      for(c=0; c<MaxTrack; c++)
			if(fTrackID[eta_index][bin].fLabel[c] == label || fTrackID[eta_index][bin].fNHits[c] == 0)
			  break;
		      if(c == MaxTrack-1) cerr<<"AliL3HoughTransformer::TransformCircle : Array reached maximum!! "<<c<<endl;
		      fTrackID[eta_index][bin].fLabel[c] = label;
		      fTrackID[eta_index][bin].fNHits[c]++;
		    }
		}
#endif
#endif
	    }
	}
      
      //Move the data pointer to the next padrow:
      AliL3MemHandler::UpdateRowPointer(tempPt);
    }
}
#endif

#ifdef COUNTINGROWS

struct EtaRow {
  UChar_t start_pad;
  UChar_t end_pad;
  Bool_t found;
  Float_t start_y;
};

void AliL3HoughTransformerGap::TransformCircle()
{
  //Assumes that all the etaslice histos are the same!
  AliL3Histogram *h = fParamSpace[0];
  //Float_t hist_bin = h->GetBinWidthY()/2.;
  Int_t first_bin = h->GetFirstYbin();
  Int_t last_bin = h->GetLastYbin();
  Float_t x_min = h->GetXmin();
  Float_t x_max = h->GetXmax();
  Float_t x_bin = (x_max-x_min)/h->GetNbinsX();
  Int_t first_binx = h->GetFirstXbin()-1;
  Int_t last_binx = h->GetLastXbin()+1;
  Int_t nbinx = h->GetNbinsX()+2;

  UChar_t last_pad;
  Int_t last_eta_index=-1;;
  EtaRow *eta_clust = new EtaRow[GetNEtaSegments()];

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
  //Loop over the padrows:
  for(Int_t i=AliL3Transform::GetFirstRow(ipatch); i<=AliL3Transform::GetLastRow(ipatch); i++)
    {
      Int_t npads = AliL3Transform::GetNPads(i)-1;

      last_pad = 255;
      //Flush eta clusters array
      memset(eta_clust,0,GetNEtaSegments()*sizeof(EtaRow));  

      Float_t x = AliL3Transform::Row2X(i);
      Float_t x2 = x*x;
      Float_t y=0.,r2=0.;

      UChar_t lrow = i-AliL3Transform::GetFirstRow(ipatch);

      //Get the data on this padrow:
      AliL3DigitData *digPt = tempPt->fDigitData;
      if(i != (Int_t)tempPt->fRow)
	{
	  cerr<<"AliL3HoughTransform::TransformCircle : Mismatching padrow numbering "<<i<<" "<<(Int_t)tempPt->fRow<<endl;
	  continue;
	}
      //      cout<<" Starting row "<<i<<endl;
      //Loop over the data on this padrow:
      for(UInt_t j=0; j<tempPt->fNDigit; j++)
	{
	  UShort_t charge = digPt[j].fCharge;
	  if((Int_t)charge <= GetLowerThreshold())
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
	  Float_t z = AliL3Transform::GetZFast(islice,(Int_t)time,fZVertex);
	  Float_t z2 = z*z;
	  //Calculate the eta:
	  Double_t r = sqrt(r2+z2);
	  Double_t eta = 0.5 * log((r+z)/(r-z));
	  //Get the corresponding index, which determines which histogram to fill:
	  Int_t eta_index = GetEtaIndex(eta);

	  if(eta_index == last_eta_index) continue;

	  //	  cout<<" Digit at patch "<<ipatch<<" row "<<i<<" pad "<<(Int_t)pad<<" time "<<time<<" etaslice "<<eta_index<<endl;
	  
	  if(eta_index < 0 || eta_index >= GetNEtaSegments())
	    continue;

	  if(!eta_clust[eta_index].found)
	    {
	      eta_clust[eta_index].start_pad = pad;
	      eta_clust[eta_index].end_pad = pad;
	      eta_clust[eta_index].start_y = y - pad_pitch/2.0;
	      eta_clust[eta_index].found = 1;
	      continue;
	    }
	  else
	    {
	      if(pad <= (eta_clust[eta_index].end_pad + 1))
		eta_clust[eta_index].end_pad = pad;
	      else
		{

		  //		  cout<<" Cluster found at etaslice "<<eta_index<<" from pad "<<(Int_t)eta_clust[eta_index].start_pad<<" to pad "<<(Int_t)eta_clust[eta_index].end_pad<<endl;
		  //Get the correct histogrampointer:
		  AliL3Histogram *hist = fParamSpace[eta_index];
		  if(!hist)
		    {
		      cerr<<"AliL3HoughTransformer::TransformCircle : Error getting histogram in index "<<eta_index<<endl;
		      continue;
		    }
		  UChar_t *nrows = fRowCount[eta_index];
		  UChar_t *ngaps = fGapCount[eta_index];
		  UChar_t *currentrow = fCurrentRowCount[eta_index];

		  //Do the transformation:
		  Float_t start_y = eta_clust[eta_index].start_y;
		  Float_t R1 = x2 + start_y*start_y; 
		  Float_t end_y = (eta_clust[eta_index].end_pad-0.5*npads)*pad_pitch + pad_pitch/2.0;
		  Float_t R2 = x2 + end_y*end_y; 

		  //Fill the histogram along the phirange
		  for(Int_t b=first_bin; b<=last_bin; b++)
		    {
#ifdef WOLUT
		      Float_t phi0 = hist->GetBinCenterY(b);
		      //		      Float_t kappa1 = 2*sin(phi1 - phi0 - hist_bin)/R1;
		      //		      Float_t kappa2 = 2*sin(phi2 - phi0 + hist_bin)/R2;
		      Float_t kappa1 = 2*(start_y*cos(phi0 + hist_bin)-x*sin(phi0 + hist_bin))/R1;
		      if(kappa1>x_max) continue;
		      Float_t kappa2 = 2*(end_y*cos(phi0 - hist_bin)-x*sin(phi0 - hist_bin))/R2;
		      if(kappa2<x_min) break;
#else
		      Float_t kappa1 = (start_y*fLUT2cosphi0low[b]-x*fLUT2sinphi0low[b])/R1;
		      if(kappa1>x_max) continue;
		      Float_t kappa2 = (end_y*fLUT2cosphi0up[b]-x*fLUT2sinphi0up[b])/R2;
		      if(kappa2<x_min) break;
#endif
		      Int_t binx1,binx2;
		      if(kappa1<x_min)
			binx1 = first_binx;
		      else
			binx1 = 1 + (Int_t)((kappa1-x_min)/x_bin);
		      if(kappa2>x_max)
			binx2 = last_binx;
		      else
			binx2 = 1 + (Int_t)((kappa2-x_min)/x_bin);
		      
		      Int_t temp_bin = b*nbinx + binx1;
		      UChar_t *nrows2 = nrows + temp_bin;
		      UChar_t *ngaps2 = ngaps + temp_bin;
		      UChar_t *currentrow2 = currentrow + temp_bin;
		      for(Int_t bin=binx1;bin<=binx2;bin++)
			{
			  //Assumes threshold about 32
			  if(*ngaps2 < 5) {
			    if(lrow != *currentrow2)
			      {
				(*nrows2)++;
				if(lrow > (*currentrow2 + 1))
				  (*ngaps2)++;
				*currentrow2=lrow;
			      }
			  }
			  nrows2++;
			  ngaps2++;
			  currentrow2++;
			}
		    }
		  //End of the transformation

		  eta_clust[eta_index].start_pad = pad;
		  eta_clust[eta_index].end_pad = pad;
		  eta_clust[eta_index].start_y = y - pad_pitch/2.0;
		}
	    }
	  last_pad = pad;
	  last_eta_index = eta_index;
	}
      //Write remaining clusters
      for(Int_t eta_index = 0;eta_index < GetNEtaSegments();eta_index++)
	{
	  //Check for empty row
	  if((eta_clust[eta_index].start_pad == 0) && (eta_clust[eta_index].end_pad == 0)) continue;

	  //Get the correct histogrampointer:
	  AliL3Histogram *hist = fParamSpace[eta_index];
	  if(!hist)
	    {
	      cerr<<"AliL3HoughTransformer::TransformCircle : Error getting histogram in index "<<eta_index<<endl;
	      continue;
	    }
	  UChar_t *nrows = fRowCount[eta_index];
	  UChar_t *ngaps = fGapCount[eta_index];
	  UChar_t *currentrow = fCurrentRowCount[eta_index];

	  //Do the transformation:
	  Float_t start_y = eta_clust[eta_index].start_y;
	  Float_t R1 = x2 + start_y*start_y; 
	  Float_t end_y = (eta_clust[eta_index].end_pad-0.5*npads)*pad_pitch + pad_pitch/2.0;
	  Float_t R2 = x2 + end_y*end_y; 

	  //Fill the histogram along the phirange
	  for(Int_t b=first_bin; b<=last_bin; b++)
	    {
#ifdef WOLUT
	      Float_t phi0 = hist->GetBinCenterY(b);
	      //		      Float_t kappa1 = 2*sin(phi1 - phi0 - hist_bin)/R1;
	      //		      Float_t kappa2 = 2*sin(phi2 - phi0 + hist_bin)/R2;
	      Float_t kappa1 = 2*(start_y*cos(phi0 + hist_bin)-x*sin(phi0 + hist_bin))/R1;
	      if(kappa1>x_max) continue;
	      Float_t kappa2 = 2*(end_y*cos(phi0 - hist_bin)-x*sin(phi0 - hist_bin))/R2;
	      if(kappa2<x_min) break;
#else
	      Float_t kappa1 = (start_y*fLUT2cosphi0low[b]-x*fLUT2sinphi0low[b])/R1;
	      if(kappa1>x_max) continue;
	      Float_t kappa2 = (end_y*fLUT2cosphi0up[b]-x*fLUT2sinphi0up[b])/R2;
	      if(kappa2<x_min) break;
#endif
	      Int_t binx1,binx2;
	      if(kappa1<x_min)
		binx1 = first_binx;
	      else
		binx1 = 1 + (Int_t)((kappa1-x_min)/x_bin);
	      if(kappa2>x_max)
		binx2 = last_binx;
	      else
		binx2 = 1 + (Int_t)((kappa2-x_min)/x_bin);

	      Int_t temp_bin = b*nbinx + binx1;
	      UChar_t *nrows2 = nrows + temp_bin;
	      UChar_t *ngaps2 = ngaps + temp_bin;
	      UChar_t *currentrow2 = currentrow + temp_bin;
	      for(Int_t bin=binx1;bin<=binx2;bin++)
		{
		  //Assumes threshold about 32
		  if(*ngaps2 < 5) {
		    if(lrow != *currentrow2)
		      {
			(*nrows2)++;
			if(lrow > (*currentrow2 + 1))
			  (*ngaps2)++;
			*currentrow2=lrow;
		      }
		  }
		  nrows2++;
		  ngaps2++;
		  currentrow2++;
		}
	    }
	  //End of the transformation
	}

      //Move the data pointer to the next padrow:
      AliL3MemHandler::UpdateRowPointer(tempPt);
    }
  delete [] eta_clust;
}
#endif

#ifdef COUNTINGROWS
/*
void AliL3HoughTransformer::CalculateNRows()
{
  Int_t ipatch = GetPatch();

  for(Int_t i=0; i<GetNEtaSegments(); i++) {
    AliL3Histogram *hist = fParamSpace[i];

    Int_t xmin = hist->GetFirstXbin();
    Int_t xmax = hist->GetLastXbin();
    Int_t ymin = hist->GetFirstYbin();
    Int_t ymax = hist->GetLastYbin();
    Int_t sum,ngaps;
    Bool_t isPrevRow;

    for(Int_t ybin=ymin; ybin<=ymax; ybin++)
      {
      for(Int_t xbin=xmin; xbin<=xmax; xbin++)
	{
	  sum = 0;
	  ngaps = 1;
	  isPrevRow = kFALSE;
	  Int_t bin = hist->GetBin(xbin,ybin);
	  for(Int_t irow=0;irow<=(AliL3Transform::GetLastRow(ipatch)-AliL3Transform::GetFirstRow(ipatch));irow++) {
	    if(((fRowCount[i][bin]>>irow)&1)==1) {
	      sum++;
	      isPrevRow = kTRUE;
	    }
	    else {
	      if(isPrevRow)
		ngaps++;
	      isPrevRow = kFALSE;
	    }
	  }
	  hist->SetBinContent(bin,sum/ngaps);
	}
      }
    hist->SetNEntries(1);
  }
}
*/

UChar_t *AliL3HoughTransformerGap::GetRowCount(Int_t eta_index)
{
  return fRowCount[eta_index];
}

UChar_t *AliL3HoughTransformerGap::GetGapCount(Int_t eta_index)
{
  return fGapCount[eta_index];
}

#endif
