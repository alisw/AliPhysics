// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliL3StandardIncludes.h"

#include "AliL3Logging.h"
#include "AliL3HoughTransformer.h"
#include "AliL3MemHandler.h"
#include "AliL3Transform.h"
#include "AliL3DigitData.h"
#include "AliL3HistogramAdaptive.h"

#if __GNUC__ == 3
using namespace std;
#endif

/** \class AliL3HoughTransformer
<pre>
//_____________________________________________________________
// AliL3HoughTransformer
//
// Hough transformation class
//
</pre>
*/

ClassImp(AliL3HoughTransformer)

AliL3HoughTransformer::AliL3HoughTransformer()
{
  //Default constructor
  fParamSpace = 0;
  fDoMC = kFALSE;;
  fEtaOverlap=kFALSE;
#ifdef do_mc
  fTrackID = 0;
#endif
}

AliL3HoughTransformer::AliL3HoughTransformer(Int_t slice,Int_t patch,Int_t n_eta_segments,Bool_t DoEtaOverlap,Bool_t /*DoMC*/) : AliL3HoughBaseTransformer(slice,patch,n_eta_segments)
{
  //Normal constructor
  fParamSpace = 0;
  fDoMC = kFALSE;
  fEtaOverlap = DoEtaOverlap;
  fDoMC=kFALSE;
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

void AliL3HoughTransformer::CreateHistograms(Float_t ptmin,Float_t ptmax,Float_t ptres,
					     Int_t nybin,Float_t psi)
{
  //Create histograms.
  //_Only_ to be used in case of the adaptive histograms!
  //phimax is given in radians!!
  
  if(ptmin > ptmax)
    {
      cerr<<"AliL3HoughTransformer::CreateHistograms: Error in ptrange "<<ptmin<<" "<<ptmax<<endl;
      return;
    }
  if(psi < 0)
    {
      cerr<<"AliL3HoughTransformer::CreateHistograms: Wrong psi-angle "<<psi<<endl;
      return;
    }
  
  fParamSpace = new AliL3Histogram*[GetNEtaSegments()];
  Char_t histname[256];
  Int_t i;
  for(i=0; i<GetNEtaSegments(); i++)
    {
      sprintf(histname,"paramspace_%d",i);
      fParamSpace[i] = new AliL3HistogramAdaptive(histname,ptmin,ptmax,ptres,nybin,-psi,psi);
    }
}

void AliL3HoughTransformer::CreateHistograms(Int_t nxbin,Float_t pt_min,
					     Int_t nybin,Float_t phimin,Float_t phimax)
{
  //Create the histograms (parameter space).
  //These are 2D histograms, span by kappa (curvature of track) and phi0 (emission angle with x-axis).
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

void AliL3HoughTransformer::CreateHistograms(Int_t nxbin,Float_t xmin,Float_t xmax,
					     Int_t nybin,Float_t ymin,Float_t ymax)
{
  
  fParamSpace = new AliL3Histogram*[GetNEtaSegments()];
  
  Char_t histname[256];
  for(Int_t i=0; i<GetNEtaSegments(); i++)
    {
      sprintf(histname,"paramspace_%d",i);
      //fParamSpace[i] = new AliL3HistogramAdaptive(histname,0.5,1.5,0.05,nybin,ymin,ymax);
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

void AliL3HoughTransformer::GetEtaIndexes(Double_t eta,Int_t *indexes)
{
  //Return histogram indexes in case of overlapping etaslices.
  
  Double_t etaslice = (GetEtaMax() - GetEtaMin())/GetNEtaSegments();
  Int_t index = (Int_t)((eta-GetEtaMin())/etaslice);
  if(index%2 == 0)
    {
      indexes[0] = index;
      indexes[1] = index - 1;
    }
  else
    {
      indexes[0] = index - 1;
      indexes[1] = index;
    }
}

inline AliL3Histogram *AliL3HoughTransformer::GetHistogram(Int_t eta_index)
{
  if(!fParamSpace || eta_index >= GetNEtaSegments() || eta_index < 0)
    return 0;
  if(!fParamSpace[eta_index])
    return 0;
  return fParamSpace[eta_index];
}

Double_t AliL3HoughTransformer::GetEta(Int_t eta_index,Int_t /*slice*/)
{
  Double_t eta_slice = (GetEtaMax()-GetEtaMin())/GetNEtaSegments();
  Double_t eta=0;
  if(fEtaOverlap)
    {
      Int_t index = eta_index + 1;
      eta=(Double_t)((index)*eta_slice);
    }
  else
    eta=(Double_t)((eta_index+0.5)*eta_slice);
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
  //which histogram in which the pixel should be transformed, the eta-value is calculated
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
	  Float_t R = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]); 
	  Float_t phi = AliL3Transform::GetPhi(xyz);
	  
	  
	  //Fill the histogram along the phirange
	  for(Int_t b=hist->GetFirstYbin(); b<=hist->GetLastYbin(); b++)
	    {
	      Float_t phi0 = hist->GetBinCenterY(b);
	      Float_t kappa = 2*sin(phi - phi0)/R;
	      //hist->Fill(kappa,phi0,(int)rint(log((Float_t)charge)));
	      hist->Fill(kappa,phi0,charge);
	      //hist->Fill(kappa,phi0,1);
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

struct Digit {
  Int_t row;
  Double_t r;
  Double_t phi;
  Int_t charge;
  Digit *next;
};

struct EtaContainer {
  Digit *first;
  Digit *last;
};

void AliL3HoughTransformer::TransformCircleC(Int_t *row_range,Int_t every)
{
  //Circle transform, using combinations of every 2 points lying
  //on different padrows and within the same etaslice.
  
  AliL3DigitRowData *tempPt = GetDataPointer();
  if(!tempPt)
    LOG(AliL3Log::kError,"AliL3HoughTransformer::TransformCircleC","Data")
      <<"No input data "<<ENDLOG;
  
  Int_t minrow = AliL3Transform::GetFirstRow(GetPatch());
  Int_t maxrow = AliL3Transform::GetLastRow(GetPatch());
  if(row_range)
    {
      minrow = row_range[0];
      maxrow = row_range[1];
      if(minrow < AliL3Transform::GetFirstRow(GetPatch()) || minrow >= AliL3Transform::GetLastRow(GetPatch()))
	minrow = AliL3Transform::GetFirstRow(GetPatch());
      if(maxrow < AliL3Transform::GetFirstRow(GetPatch()) || maxrow >= AliL3Transform::GetLastRow(GetPatch()))
	maxrow = AliL3Transform::GetLastRow(GetPatch());
      if(minrow > maxrow || maxrow==minrow)
	{
	  cerr<<"AliL3HoughTransformer::TransformCircleC : Bad row range "<<minrow<<" "<<maxrow<<endl;
	  return;
	}
    }
  else
    {
      minrow = AliL3Transform::GetFirstRow(GetPatch());
      maxrow = AliL3Transform::GetLastRow(GetPatch());
    }
      
  Int_t counter=0;
  for(Int_t i=AliL3Transform::GetFirstRow(GetPatch()); i<=AliL3Transform::GetLastRow(GetPatch()); i++)
    {
      counter += tempPt->fNDigit;
      AliL3MemHandler::UpdateRowPointer(tempPt);
    }
  
  Int_t bound = (GetNEtaSegments()+1)*(AliL3Transform::GetNRows(GetPatch())+1);
  EtaContainer *etaPt = new EtaContainer[bound];
  memset(etaPt,0,bound*sizeof(EtaContainer));  
  
  Digit *digits = new Digit[counter];
  cout<<"Allocating "<<counter*sizeof(Digit)<<" bytes to digitsarray"<<endl;
  memset(digits,0,counter*sizeof(Digit));

  Int_t sector,row,tot_charge,pad,time,charge;
  Double_t r1,r2,phi1,phi2,eta,kappa,phi_0;
  Float_t xyz[3];
  
  counter=0;
  tempPt = GetDataPointer();
  
  cout<<"Calculating digits in patch "<<GetPatch()<<endl;
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
	  digits[counter].charge = charge;

	  if(!fEtaOverlap)
	    {
	      Int_t eta_index = GetEtaIndex(eta);
	      
	      Int_t index = (GetNEtaSegments()+1)*(i-AliL3Transform::GetFirstRow(GetPatch())) + eta_index;
	      
	      if(index > 0 && index < bound) 
		{
		  if(etaPt[index].first == 0)
		    etaPt[index].first = &digits[counter];
		  else
		    (etaPt[index].last)->next = &digits[counter];
		  etaPt[index].last = &digits[counter];
		}
	    }
	  else
	    {
	      Int_t eta_index[2];
	      GetEtaIndexes(eta,eta_index);
	      Int_t index[2];
	      index[0] = (GetNEtaSegments()+1)*(i-AliL3Transform::GetFirstRow(GetPatch())) + eta_index[0];
	      index[1] = (GetNEtaSegments()+1)*(i-AliL3Transform::GetFirstRow(GetPatch())) + eta_index[1];
	      if(index[0] == index[1])
		{
		  cerr<<"Same etaindexes "<<index[0]<<" "<<index[1]<<endl;
		  exit(5);
		}
	      
	      Int_t ind = index[0];
	      if(ind > 0 && ind < bound)
		{
		  if(etaPt[ind].first == 0)
		    etaPt[ind].first = &digits[counter];
		  else
		    (etaPt[ind].last)->next = &digits[counter];
		  etaPt[ind].last = &digits[counter];
		}
	      
	      ind = index[1];
	      if(ind > 0 && ind < bound)
		{
		  if(etaPt[ind].first == 0)
		    etaPt[ind].first = &digits[counter];
		  else
		    (etaPt[ind].last)->next = &digits[counter];
		  etaPt[ind].last = &digits[counter];
		}
	    }

	  counter++;
	}
      AliL3MemHandler::UpdateRowPointer(tempPt);
    }
  
  cout<<"Doing the combinatorics"<<endl;
  
  Digit *dPt1,*dPt2;
  
  for(Int_t e=0; e<GetNEtaSegments(); e++)
    {
      for(Int_t i=minrow; i<=maxrow; i+=every)
	{
	  Int_t index1 = (GetNEtaSegments()+1)*(i-AliL3Transform::GetFirstRow(GetPatch())) + e;
	  
	  for(dPt1 = (Digit*)etaPt[index1].first; dPt1 != 0; dPt1 = (Digit*)dPt1->next)
	    {
	      for(Int_t j=i+every; j<=maxrow; j+=every)
		{
		  Int_t index2 = (GetNEtaSegments()+1)*(j-AliL3Transform::GetFirstRow(GetPatch())) + e;
		  
		  for(dPt2 = (Digit*)etaPt[index2].first; dPt2 != 0; dPt2 = (Digit*)dPt2->next)
		    {
		      if(dPt1->row == dPt2->row)
			{
			  cerr<<"same row; indexes "<<index1<<" "<<index2<<endl;
			  exit(5);
			}
		      
		      //Get the correct histogrampointer:
		      AliL3Histogram *hist = fParamSpace[e];
		      if(!hist)
			{
			  printf("AliL3HoughTransformer::TransformCircleC() : No histogram at index %d\n",i);
			  continue;
			}
		      
		      //Do the transform:
		      r1 = dPt1->r;
		      phi1 = dPt1->phi;
		      r2 = dPt2->r;
		      phi2 = dPt2->phi;
		      phi_0 = atan( (r2*sin(phi1)-r1*sin(phi2))/(r2*cos(phi1)-r1*cos(phi2)) );
		      kappa = 2*sin(phi2-phi_0)/r2;
		      tot_charge = dPt1->charge + dPt2->charge;
		      hist->Fill(kappa,phi_0,tot_charge);
		      
		    }
		}
	    }
	}
    }

  cout<<"done"<<endl;
  delete [] etaPt;
  delete [] digits;

}

void AliL3HoughTransformer::TransformLine(Int_t *row_range,Float_t *phirange)
{
  //Do a line transform on the data.


  AliL3DigitRowData *tempPt = GetDataPointer();
  if(!tempPt)
    {
      LOG(AliL3Log::kError,"AliL3HoughTransformer::TransformLine","Data")
	<<"No input data "<<ENDLOG;
      return;
    }
  
  Int_t minrow = AliL3Transform::GetFirstRow(GetPatch());
  Int_t maxrow = AliL3Transform::GetLastRow(GetPatch());
  if(row_range)
    {
      minrow = row_range[0];
      maxrow = row_range[1];
      if(minrow < AliL3Transform::GetFirstRow(GetPatch()) || minrow >= AliL3Transform::GetLastRow(GetPatch()))
	minrow = AliL3Transform::GetFirstRow(GetPatch());
      if(maxrow < AliL3Transform::GetFirstRow(GetPatch()) || maxrow >= AliL3Transform::GetLastRow(GetPatch()))
	maxrow = AliL3Transform::GetLastRow(GetPatch());
      if(minrow > maxrow || maxrow==minrow)
	{
	  cerr<<"AliL3HoughTransformer::TransformCircleC : Bad row range "<<minrow<<" "<<maxrow<<endl;
	  return;
	}
    }

  for(Int_t i=minrow; i<=maxrow; i++)
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
	  
	  if(phirange)
	    {
	      Float_t phi = AliL3Transform::GetPhi(xyz);
	      if(phi < phirange[0] || phi > phirange[1])
		continue;
	    }
	  Float_t eta = AliL3Transform::GetEta(xyz);
	  Int_t eta_index = GetEtaIndex(eta);//(Int_t)(eta/etaslice);
	  if(eta_index < 0 || eta_index >= GetNEtaSegments())
	    continue;
	  
	  xyz[0] = xyz[0] - AliL3Transform::Row2X(minrow);
	  
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

struct LDigit {
  Int_t row;
  Int_t charge;
  Float_t y;
  LDigit *next;
};
struct LEtaContainer {
  LDigit *first;
  LDigit *last;
};
void AliL3HoughTransformer::TransformLineC(Int_t *rowrange,Float_t *phirange)
{
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
  
  Int_t bound = (GetNEtaSegments()+1)*(AliL3Transform::GetNRows(GetPatch())+1);
  LEtaContainer *etaPt = new LEtaContainer[bound];
  memset(etaPt,0,bound*sizeof(LEtaContainer));  
  
  LDigit *digits = new LDigit[counter];
  cout<<"Allocating "<<counter*sizeof(LDigit)<<" bytes to digitsarray"<<endl;
  memset(digits,0,counter*sizeof(LDigit));

  Int_t sector,row;
  Float_t xyz[3];
  
  counter=0;
  tempPt = GetDataPointer();
  
  cout<<"Calculating digits in patch "<<GetPatch()<<endl;
  for(Int_t i=AliL3Transform::GetFirstRow(GetPatch()); i<=AliL3Transform::GetLastRow(GetPatch()); i++)
    {
      AliL3DigitData *digPt = tempPt->fDigitData;
      for(UInt_t di=0; di<tempPt->fNDigit; di++)
	{
	  Int_t charge = digPt[di].fCharge;
	  Int_t pad = digPt[di].fPad;
	  Int_t time = digPt[di].fTime;
	  AliL3Transform::Slice2Sector(GetSlice(),i,sector,row);
	  AliL3Transform::Raw2Local(xyz,sector,row,(Int_t)pad,(Int_t)time);
	  Double_t eta = AliL3Transform::GetEta(xyz);
	  
	  Float_t phi = atan2(xyz[1],xyz[0]);
	  if(phi < phirange[0] || phi > phirange[1]) continue;
	  
	  digits[counter].row = i;
	  digits[counter].y = xyz[1];
	  digits[counter].charge = charge;
	  
	  Int_t eta_index = GetEtaIndex(eta);
	  Int_t index = (GetNEtaSegments()+1)*(i-AliL3Transform::GetFirstRow(GetPatch())) + eta_index;
	  
	  if(index > 0 && index < bound) 
	    {
	      if(etaPt[index].first == 0)
		etaPt[index].first = &digits[counter];
	      else
		(etaPt[index].last)->next = &digits[counter];
	      etaPt[index].last = &digits[counter];
	    }
	  counter++;
	}
      AliL3MemHandler::UpdateRowPointer(tempPt);
    }
  
  cout<<"Doing the combinatorics"<<endl;
  
  LDigit *dPt1,*dPt2;
  
  for(Int_t e=0; e<GetNEtaSegments(); e++)
    {
      for(Int_t i=rowrange[0]; i<=rowrange[1]; i++)
	{
	  Int_t index1 = (GetNEtaSegments()+1)*(i-AliL3Transform::GetFirstRow(GetPatch())) + e;
	  
	  for(dPt1 = (LDigit*)etaPt[index1].first; dPt1 != 0; dPt1 = (LDigit*)dPt1->next)
	    {
	      for(Int_t j=i+1; j<=rowrange[1]; j++)
		{
		  Int_t index2 = (GetNEtaSegments()+1)*(j-AliL3Transform::GetFirstRow(GetPatch())) + e;
		  
		  for(dPt2 = (LDigit*)etaPt[index2].first; dPt2 != 0; dPt2 = (LDigit*)dPt2->next)
		    {
		      if(dPt1->row == dPt2->row)
			{
			  cerr<<"same row; indexes "<<index1<<" "<<index2<<endl;
			  exit(5);
			}
		      
		      //Get the correct histogrampointer:
		      AliL3Histogram *hist = fParamSpace[e];
		      if(!hist)
			{
			  printf("AliL3HoughTransformer::TransformCircleC() : No histogram at index %d\n",i);
			  continue;
			}
		      
		      //Do the transform:
		      float x1 = AliL3Transform::Row2X(dPt1->row) - AliL3Transform::Row2X(rowrange[0]);
		      float x2 = AliL3Transform::Row2X(dPt2->row) - AliL3Transform::Row2X(rowrange[0]);
		      float y1 = dPt1->y;
		      float y2 = dPt2->y;
		      float theta = atan2(x2-x1,y1-y2);
		      float rho = x1*cos(theta)+y1*sin(theta);
		      hist->Fill(theta,rho,1);//dPt1->charge+dPt2->charge);
		    }
		}
	    }
	}
    }

  cout<<"done"<<endl;
  delete [] etaPt;
  delete [] digits;
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
  //nhits = max;
  return label;
#endif
  cout<<"AliL3HoughTransformer::GetTrackID : Compile with do_mc flag!"<<endl;
  return -1;
}

