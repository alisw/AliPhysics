//$Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ASV 

#include "AliL3StandardIncludes.h"

#include "AliL3Logging.h"
#include "AliL3HoughTransformer.h"
#include "AliL3MemHandler.h"
#include "AliL3Transform.h"
#include "AliL3DigitData.h"
#include "AliL3HistogramAdaptive.h"

#if GCCVERSION == 3
using namespace std;
#endif

//_____________________________________________________________
// AliL3HoughTransformer
//
// Hough transformation class
//

ClassImp(AliL3HoughTransformer)

AliL3HoughTransformer::AliL3HoughTransformer()
{
  //Default constructor
  fParamSpace = 0;
  fDoMC = kFALSE;
#ifdef do_mc
  fTrackID = 0;
#endif
}

AliL3HoughTransformer::AliL3HoughTransformer(Int_t slice,Int_t patch,Int_t n_eta_segments) : AliL3HoughBaseTransformer(slice,patch,n_eta_segments)
{
  //Normal constructor
  fParamSpace = 0;
  fDoMC = kFALSE;
#ifdef do_mc
  fTrackID = 0;
#endif
}

AliL3HoughTransformer::~AliL3HoughTransformer()
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
}

//void AliL3HoughTransformer::Init(Int_t slice=0,Int_t patch=0,Int_t n_eta_segments=100){}

void AliL3HoughTransformer::DeleteHistograms()
{
  if(!fParamSpace)
    return;
  for(Int_t i=0; i<GetNEtaSegments(); i++)
    {
      if(!fParamSpace[i]) continue;
      delete fParamSpace[i];
    }
  delete [] fParamSpace;
}

void AliL3HoughTransformer::CreateHistograms(Int_t nxbin,Double_t pt_min,
					     Int_t nybin,Double_t phimin,Double_t phimax)
{
  //Create the histograms (parameter space).
  //These are 2D histograms, span by kappa (curvature of track) and phi0 (emission angle with x-axis).
  //The arguments give the range and binning; 
  //nxbin = #bins in kappa
  //nybin = #bins in phi0
  //pt_min = mimium Pt of track (corresponding to maximum kappa)
  //phi_min = mimimum phi0 (degrees)
  //phi_max = maximum phi0 (degrees)
    
  Double_t x = AliL3Transform::GetBFact()*AliL3Transform::GetBField()/pt_min;
  Double_t torad = AliL3Transform::Pi()/180;
  CreateHistograms(nxbin,-1.*x,x,nybin,phimin*torad,phimax*torad);
}

void AliL3HoughTransformer::CreateHistograms(Int_t nxbin,Double_t xmin,Double_t xmax,
					     Int_t nybin,Double_t ymin,Double_t ymax)
{
  
  fParamSpace = new AliL3Histogram*[GetNEtaSegments()];
  
  Char_t histname[256];
  for(Int_t i=0; i<GetNEtaSegments(); i++)
    {
      sprintf(histname,"paramspace_%d",i);
      //fParamSpace[i] = new AliL3HistogramAdaptive(histname,0.1,1,0.05,64,ymin,ymax);
      fParamSpace[i] = new AliL3Histogram(histname,"",nxbin,xmin,xmax,nybin,ymin,ymax);
    }
  
#ifdef do_mc
  if(fDoMC)
    {
      AliL3Histogram *hist = fParamSpace[0];
      Int_t ncells = (hist->GetNbinsX()+2)*(hist->GetNbinsY()+2);
      cout<<"Allocating "<<GetNEtaSegments()*ncells*sizeof(TrackIndex)<<" bytes to fTrackID"<<endl;
      fTrackID = new TrackIndex*[GetNEtaSegments()];
      for(Int_t i=0; i<GetNEtaSegments(); i++)
	fTrackID[i] = new TrackIndex[ncells];
    }
#endif
}

void AliL3HoughTransformer::Reset()
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
}

Int_t AliL3HoughTransformer::GetEtaIndex(Double_t eta)
{
  //Return the histogram index of the corresponding eta. 

  Double_t etaslice = (GetEtaMax() - GetEtaMin())/GetNEtaSegments();
  Double_t index = (eta-GetEtaMin())/etaslice;
  return (Int_t)index;
}

inline AliL3Histogram *AliL3HoughTransformer::GetHistogram(Int_t eta_index)
{
  if(!fParamSpace || eta_index >= GetNEtaSegments() || eta_index < 0)
    return 0;
  if(!fParamSpace[eta_index])
    return 0;
  return fParamSpace[eta_index];
}

Double_t AliL3HoughTransformer::GetEta(Int_t eta_index,Int_t slice)
{
  Double_t eta_slice = (GetEtaMax()-GetEtaMin())/GetNEtaSegments();
  Double_t eta=(Double_t)((eta_index+0.5)*eta_slice);
  if(slice>17) eta*=-1;
  return eta;
}

void AliL3HoughTransformer::TransformCircle()
{
  //Transform the input data with a circle HT.
  //The function loops over all the data, and transforms each pixel with the equations:
  // 
  //kappa = 2/R*sin(phi - phi0)
  //
  //where R = sqrt(x*x +y*y), and phi = arctan(y/x)
  //
  //Each pixel then transforms into a curve in the (kappa,phi0)-space. In order to find
  //which histogram in which the pixel should be transformed, the eta-value is calcluated
  //and the proper histogram index is found by GetEtaIndex(eta).


  AliL3DigitRowData *tempPt = GetDataPointer();
  if(!tempPt)
    {
      LOG(AliL3Log::kError,"AliL3HoughTransformer::TransformCircle","Data")
	<<"No input data "<<ENDLOG;
      return;
    }
  
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
	  if((Int_t)charge <= GetLowerThreshold() || (Int_t)charge > GetUpperThreshold())
	    continue;
	  Int_t sector,row;
	  Float_t xyz[3];
	  
	  //Transform data to local cartesian coordinates:
	  AliL3Transform::Slice2Sector(GetSlice(),i,sector,row);
	  AliL3Transform::Raw2Local(xyz,sector,row,(Int_t)pad,(Int_t)time);
	  
	  //Calculate the eta:
	  Double_t eta = AliL3Transform::GetEta(xyz);
	  
	  //Get the corresponding index, which determines which histogram to fill:
	  Int_t eta_index = GetEtaIndex(eta);
	  if(eta_index < 0 || eta_index >= GetNEtaSegments())
	    continue;
	  
	  //Get the correct histogrampointer:
	  AliL3Histogram *hist = fParamSpace[eta_index];
	  if(!hist)
	    {
	      printf("AliL3HoughTransformer::TransformCircle : Error getting histogram in index %d\n",eta_index);
	      continue;
	    }

	  //Do the transformation:
	  Float_t R = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]); 
	  Float_t phi = AliL3Transform::GetPhi(xyz);
	  
	  //Fill the histogram along the phirange
	  for(Int_t b=hist->GetFirstYbin(); b<=hist->GetLastYbin(); b++)
	    {
	      Float_t phi0 = hist->GetBinCenterY(b);
	      Float_t kappa = 2*sin(phi - phi0)/R;
	      hist->Fill(kappa,phi0,charge);
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
	    }
	}
      
      //Move the data pointer to the next padrow:
      AliL3MemHandler::UpdateRowPointer(tempPt);
    }
}

void AliL3HoughTransformer::TransformCircleC(Int_t row_range)
{
  //Circle transform, using combinations of every 2 points lying
  //on different padrows and within the same etaslice.
  
  AliL3DigitRowData *tempPt = GetDataPointer();
  if(!tempPt)
    LOG(AliL3Log::kError,"AliL3HoughTransformer::TransformCircleC","Data")
      <<"No input data "<<ENDLOG;
 
  Int_t counter=0;
  for(Int_t i=AliL3Transform::GetFirstRow(GetPatch()); i<=AliL3Transform::GetLastRow(GetPatch()); i++)
    {
      counter += tempPt->fNDigit;
      AliL3MemHandler::UpdateRowPointer(tempPt);
    }
  
  struct Digit {
    Int_t row;
    Double_t r;
    Double_t phi;
    Int_t eta_index;
    Int_t charge;
    Int_t trackID[3];
  };
  
  Digit *digits = new Digit[counter];
  cout<<"Allocating "<<counter*sizeof(Digit)<<" bytes to digitsarray"<<endl;
  
  Int_t total_digits=counter;
  Int_t sector,row,tot_charge,pad,time,charge;
  Double_t r1,r2,phi1,phi2,eta,kappa,phi_0;
  Float_t xyz[3];
  
  counter=0;
  tempPt = GetDataPointer();
  
  for(Int_t i=AliL3Transform::GetFirstRow(GetPatch()); i<=AliL3Transform::GetLastRow(GetPatch()); i++)
    {
      AliL3DigitData *digPt = tempPt->fDigitData;
      for(UInt_t di=0; di<tempPt->fNDigit; di++)
	{
	  charge = digPt[di].fCharge;
	  pad = digPt[di].fPad;
	  time = digPt[di].fTime;
	  AliL3Transform::Slice2Sector(GetSlice(),i,sector,row);
	  AliL3Transform::Raw2Local(xyz,sector,row,(Int_t)pad,(Int_t)time);
	  eta = AliL3Transform::GetEta(xyz);
	  digits[counter].row = i;
	  digits[counter].r = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]);
	  digits[counter].phi = atan2(xyz[1],xyz[0]);
	  digits[counter].eta_index = GetEtaIndex(eta);
	  digits[counter].charge = charge;
#ifdef do_mc
	  if(fDoMC)
	    {
	      digits[counter].trackID[0] = digPt[di].fTrackID[0];
	      digits[counter].trackID[1] = digPt[di].fTrackID[1];
	      digits[counter].trackID[2] = digPt[di].fTrackID[2];
	    }
#endif
	  counter++;
	}
      AliL3MemHandler::UpdateRowPointer(tempPt);
    }
  
  for(Int_t i=0; i<total_digits; i++)
    {
      if(digits[i].eta_index < 0 || digits[i].eta_index >= GetNEtaSegments()) continue;
      Int_t ind = digits[i].eta_index;
      
      for(Int_t j=i+1; j<total_digits; j++)
	{
	  if(digits[i].row == digits[j].row) continue;
	  if(digits[i].eta_index != digits[j].eta_index) continue;
	  if(digits[i].row + row_range < digits[j].row) break;
	  
	  //Get the correct histogrampointer:
	  AliL3Histogram *hist = fParamSpace[ind];
	  if(!hist)
	    {
	      printf("AliL3HoughTransformer::TransformCircleC() : No histogram at index %d\n",ind);
	      continue;
	    }
	  
	  r1 = digits[i].r;
	  phi1 = digits[i].phi;
	  r2 = digits[j].r;
	  phi2 = digits[j].phi;
	  phi_0 = atan( (r2*sin(phi1)-r1*sin(phi2))/(r2*cos(phi1)-r1*cos(phi2)) );
	  kappa = 2*sin(phi2-phi_0)/r2;
	  tot_charge = digits[i].charge + digits[j].charge;
	  hist->Fill(kappa,phi_0,tot_charge);
#ifdef do_mc
	  if(fDoMC)
	    {
	      Int_t bin = hist->FindBin(kappa,phi_0);
	      for(Int_t l=0; l<3; l++)
		{
		  for(Int_t m=0; m<3; m++)
		    {
		      if(digits[i].trackID[l] == digits[j].trackID[m])
			{
			  Int_t label = digits[i].trackID[l];
			  if(label < 0) continue;
			  UInt_t c;
			  for(c=0; c<MaxTrack; c++)
			    if(fTrackID[ind][bin].fLabel[c] == label || fTrackID[ind][bin].fNHits[c] == 0)
			      break;
			  if(c == MaxTrack-1) cerr<<"AliL3HoughTransformer::TransformCircleC : Array reached maximum!! "<<c<<endl;
			  fTrackID[ind][bin].fLabel[c] = label;
			  fTrackID[ind][bin].fNHits[c]++;
			}
		    }
		}
	    }
#endif
	}
    }
  delete [] digits;
}

void AliL3HoughTransformer::TransformLine()
{
  //Do a line transform on the data.

  
  AliL3DigitRowData *tempPt = GetDataPointer();
  if(!tempPt)
    {
      LOG(AliL3Log::kError,"AliL3HoughTransformer::TransformLine","Data")
	<<"No input data "<<ENDLOG;
      return;
    }
    
  for(Int_t i=AliL3Transform::GetFirstRow(GetPatch()); i<=AliL3Transform::GetLastRow(GetPatch()); i++)
    {
      AliL3DigitData *digPt = tempPt->fDigitData;
      if(i != (Int_t)tempPt->fRow)
	{
	  printf("AliL3HoughTransform::TransformLine : Mismatching padrow numbering\n");
	  continue;
	}
      for(UInt_t j=0; j<tempPt->fNDigit; j++)
	{
	  UShort_t charge = digPt[j].fCharge;
	  UChar_t pad = digPt[j].fPad;
	  UShort_t time = digPt[j].fTime;
	  if(charge < GetLowerThreshold())
	    continue;
	  Int_t sector,row;
	  Float_t xyz[3];
	  AliL3Transform::Slice2Sector(GetSlice(),i,sector,row);
	  AliL3Transform::Raw2Local(xyz,sector,row,(Int_t)pad,(Int_t)time);
	  Float_t eta = AliL3Transform::GetEta(xyz);
	  Int_t eta_index = GetEtaIndex(eta);//(Int_t)(eta/etaslice);
	  if(eta_index < 0 || eta_index >= GetNEtaSegments())
	    continue;
	  
	  //Get the correct histogram:
	  AliL3Histogram *hist = fParamSpace[eta_index];
	  if(!hist)
	    {
	      printf("AliL3HoughTransformer::TransformLine : Error getting histogram in index %d\n",eta_index);
	      continue;
	    }
	  for(Int_t xbin=hist->GetFirstXbin(); xbin<hist->GetLastXbin(); xbin++)
	    {
	      Double_t theta = hist->GetBinCenterX(xbin);
	      Double_t rho = xyz[0]*cos(theta) + xyz[1]*sin(theta);
	      hist->Fill(theta,rho,charge);
	    }
	}
      AliL3MemHandler::UpdateRowPointer(tempPt);
    }
  
}

Int_t AliL3HoughTransformer::GetTrackID(Int_t eta_index,Double_t kappa,Double_t psi)
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
  for(UInt_t i=0; i<MaxTrack; i++)
    {
      Int_t nhits=fTrackID[eta_index][bin].fNHits[i];
      if(nhits == 0) break;
      if(nhits > max)
	{
	  max = nhits;
	  label = fTrackID[eta_index][bin].fLabel[i];
	}
    }
  return label;
#endif
  cout<<"AliL3HoughTransformer::GetTrackID : Compile with do_mc flag!"<<endl;
  return -1;
}

