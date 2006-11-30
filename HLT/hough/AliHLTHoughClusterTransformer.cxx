// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliHLTStandardIncludes.h"

#include "AliHLTLogging.h"
#include "AliHLTHoughClusterTransformer.h"
#include "AliHLTMemHandler.h"
#include "AliHLTSpacePointData.h"
#include "AliHLTTransform.h"
#include "AliHLTDigitData.h"
#include "AliHLTHistogram.h"
#include "AliHLTClustFinderNew.h"

#if __GNUC__ >= 3
using namespace std;
#endif

//_____________________________________________________________
// AliHLTHoughClusterTransformer
//
// Hough transformation class. 
//
// This class performs the Hough transform on _clusters_, and not raw data.
// Therefore, it first finds the clusters using the HLT cluster finder, and
// then uses these as input for the transform.

ClassImp(AliHLTHoughClusterTransformer)

AliHLTHoughClusterTransformer::AliHLTHoughClusterTransformer()
{
  //Default constructor
  fParamSpace = 0;
  fNClusters = 0;
  fClusters = 0;
  fMemHandler = 0;
#ifdef do_mc
  fTrackID = 0;
#endif
}

AliHLTHoughClusterTransformer::AliHLTHoughClusterTransformer(Int_t slice,Int_t patch,Int_t netasegments) : AliHLTHoughBaseTransformer(slice,patch,netasegments)
{
  //Normal constructor
  fParamSpace = 0;
  fNClusters=0;
  fClusters=0;
  fMemHandler=0;
#ifdef do_mc
  fTrackID=0;
#endif
}

AliHLTHoughClusterTransformer::~AliHLTHoughClusterTransformer()
{
  //dtor
  DeleteHistograms();
  if(fMemHandler)
    delete fMemHandler;
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

//void AliHLTHoughClusterTransformer::Init(Int_t slice=0,Int_t patch=0,Int_t netasegments=100){}

void AliHLTHoughClusterTransformer::DeleteHistograms()
{
  //Delete hough space histograms
  if(fParamSpace)
    {
      for(Int_t i=0; i<GetNEtaSegments(); i++)
	{
	  if(!fParamSpace[i]) continue;
	  delete fParamSpace[i];
	}
      delete [] fParamSpace;
    }
}

void AliHLTHoughClusterTransformer::CreateHistograms(Int_t nxbin,Float_t ptmin,
					     Int_t nybin,Float_t phimin,Float_t phimax)
{
  //Create the histograms (parameter space).
  //These are 2D histograms, span by kappa (curvature of track) and phi0 (emission angle with x-axis).
  //The arguments give the range and binning; 
  //nxbin = #bins in kappa
  //nybin = #bins in phi0
  //ptmin = mimium Pt of track (corresponding to maximum kappa)
  //phimin = mimimum phi0 (degrees)
  //phimax = maximum phi0 (degrees)
    
  Double_t x = AliHLTTransform::GetBFact()*AliHLTTransform::GetBField()/ptmin;
  Double_t torad = AliHLTTransform::Pi()/180;
  CreateHistograms(nxbin,-1.*x,x,nybin,phimin*torad,phimax*torad);
}

void AliHLTHoughClusterTransformer::CreateHistograms(Int_t nxbin,Float_t xmin,Float_t xmax,
						    Int_t nybin,Float_t ymin,Float_t ymax)
{
  //Create histograms which contain hough space
  Char_t histname[256];
  
  fParamSpace = new AliHLTHistogram*[GetNEtaSegments()];
  for(Int_t i=0; i<GetNEtaSegments(); i++)
    {
      sprintf(histname,"paramspace_%d",i);
      fParamSpace[i] = new AliHLTHistogram(histname,"",nxbin,xmin,xmax,nybin,ymin,ymax);
    }
#ifdef do_mc
  Int_t ncells = (nxbin+2)*(nybin+2);
  cout<<"Allocating "<<GetNEtaSegments()*ncells*sizeof(AliHLTTrackIndex)<<" bytes to fTrackID"<<endl;
  fTrackID = new AliHLTTrackIndex*[GetNEtaSegments()];
  for(Int_t i=0; i<GetNEtaSegments(); i++)
    fTrackID[i] = new AliHLTTrackIndex[ncells];
#endif
  
}

void AliHLTHoughClusterTransformer::Reset()
{
  //Reset all the histograms

  if(!fParamSpace)
    {
      LOG(AliHLTLog::kWarning,"AliHLTHoughClusterTransformer::Reset","Histograms")
	<<"No histograms to reset"<<ENDLOG;
      return;
    }
  
  for(Int_t i=0; i<GetNEtaSegments(); i++)
    fParamSpace[i]->Reset();

  if(fMemHandler)
    delete fMemHandler;
  fNClusters=0;
  fClusters=0;
#ifdef do_mc
  AliHLTHistogram *hist = fParamSpace[0];
  Int_t ncells = (hist->GetNbinsX()+2)*(hist->GetNbinsY()+2);
  for(Int_t i=0; i<GetNEtaSegments(); i++)
    memset(fTrackID[i],0,ncells*sizeof(AliHLTTrackIndex));
#endif
}

Int_t AliHLTHoughClusterTransformer::GetEtaIndex(Double_t eta) const
{
  //Return the histogram index of the corresponding eta. 
  
  Double_t etaslice = (GetEtaMax() - GetEtaMin())/GetNEtaSegments();
  Double_t index = (eta-GetEtaMin())/etaslice;
  return (Int_t)index;
}

inline AliHLTHistogram *AliHLTHoughClusterTransformer::GetHistogram(Int_t etaindex)
{
  //Returns the histogram which correspond to a given eta slice
  if(!fParamSpace || etaindex >= GetNEtaSegments() || etaindex < 0)
    return 0;
  if(!fParamSpace[etaindex])
    return 0;
  return fParamSpace[etaindex];
}

Double_t AliHLTHoughClusterTransformer::GetEta(Int_t etaindex,Int_t slice) const
{
  //Returns eta associated with given eta slice
  Double_t etaslice = (GetEtaMax()-GetEtaMin())/GetNEtaSegments();
  Double_t eta=(Double_t)((etaindex+0.5)*etaslice);
  if(slice>17) eta*=-1;
  return eta;
}

void AliHLTHoughClusterTransformer::FindClusters()
{
  //Find the clusters
  if(!GetDataPointer())
    {
      cerr<<"AliHLTHoughClusterTransformer::FindClusters : Zero data pointer"<<endl;
      return;
    }

  const Int_t kMaxpoints=100000;
  const Int_t kPointsize = kMaxpoints * sizeof(AliHLTSpacePointData);

  fMemHandler = new AliHLTMemHandler();
  fClusters = (AliHLTSpacePointData*)fMemHandler->Allocate(kPointsize);
  AliHLTClustFinderNew *cf = new AliHLTClustFinderNew();
  cf->InitSlice(0,GetPatch(),AliHLTTransform::GetFirstRow(GetPatch()),AliHLTTransform::GetLastRow(GetPatch()),kMaxpoints);
  cf->SetDeconv(kFALSE);
  cf->SetOutputArray(fClusters);
  cf->Read(1,GetDataPointer());
  cf->ProcessDigits();
  fNClusters = cf->GetNumberOfClusters();
  delete cf;
}

void AliHLTHoughClusterTransformer::TransformCircle()
{
  //Transform the input data with a circle HT.
  //The function loops over all the data, and transforms each cluster with the equations:
  // 
  //kappa = 2/r*sin(phi - phi0)
  //
  //where r = sqrt(x*x +y*y), and phi = arctan(y/x)
  //
  //Each cluster then transforms into a curve in the (kappa,phi0)-space. In order to find
  //which histogram in which the pixel should be transformed, the eta-value is calcluated
  //and the proper histogram index is found by GetEtaIndex(eta).

  FindClusters();
  if(!fClusters)
    {
      LOG(AliHLTLog::kError,"AliHLTHoughClusterTransformer::TransformCircle","Data")
	<<"No input data "<<ENDLOG;
      return;
    }
  
  Float_t xyz[3];
  for(Int_t i=0; i<fNClusters; i++)
    {
      xyz[0] = fClusters[i].fX;
      xyz[1] = fClusters[i].fY;
      xyz[2] = fClusters[i].fZ;
      Double_t eta = AliHLTTransform::GetEta(xyz);
      Int_t etaindex = GetEtaIndex(eta);
      if(etaindex < 0 || etaindex >= GetNEtaSegments())
	continue;
      
      AliHLTHistogram *hist = fParamSpace[etaindex];
      
      //Do the transformation:
      Double_t r = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]); 
      Double_t phi = AliHLTTransform::GetPhi(xyz);
      
      //Fill the histogram along the phirange
      for(Int_t b=hist->GetFirstYbin(); b<=hist->GetLastYbin(); b++)
	{
	  Double_t phi0 = hist->GetBinCenterY(b);
	  Double_t kappa = 2*sin(phi - phi0)/r;
	  hist->Fill(kappa,phi0,1);
#ifdef do_mc
	  Int_t bin = hist->FindBin(kappa,phi0);
	  for(Int_t t=0; t<3; t++)
	    {
	      Int_t label = fClusters[i].fTrackID[t];
	      if(label < 0) break;
	      UInt_t c;
	      for(c=0; c<MaxTrack; c++)
		if(fTrackID[etaindex][bin].fLabel[c] == label || fTrackID[etaindex][bin].fNHits[c] == 0)
		  break;
	      if(c == MaxTrack-1) cerr<<"AliHLTHoughClusterTransformer::TransformCircle : Array reached maximum!! "<<c<<endl;
	      fTrackID[etaindex][bin].fLabel[c] = label;
	      fTrackID[etaindex][bin].fNHits[c]++;
	    }
#endif
	}
    }
}

void AliHLTHoughClusterTransformer::TransformCircleC(Int_t */*row_range*/,Int_t /*every*/)
{
  //Circle transform, using combinations of every 2 points lying
  //on different padrows and within the same etaslice.
  
  FindClusters();
  if(!fClusters)
    LOG(AliHLTLog::kError,"AliHLTHoughClusterTransformer::TransformCircleC","Data")
      <<"No input data "<<ENDLOG;
  
  Float_t xyz[3];
  Double_t eta,r1,phi1,r2,phi2,kappa,psi;
  Int_t index1,index2;
  for(Int_t i=0; i<fNClusters; i++)
    {
      xyz[0] = fClusters[i].fX;
      xyz[1] = fClusters[i].fY;
      xyz[2] = fClusters[i].fZ;
      r1 = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]);
      phi1 = atan2(xyz[1],xyz[0]);
      eta = AliHLTTransform::GetEta(xyz);
      index1 = GetEtaIndex(eta);
      if(index1 < 0 || index1 >= GetNEtaSegments()) continue;
      
      for(Int_t j=i+1; j<fNClusters; j++)
	{
	  if(fClusters[j].fPadRow == fClusters[i].fPadRow) continue;
	  xyz[0] = fClusters[j].fX;
	  xyz[1] = fClusters[j].fY;
	  xyz[2] = fClusters[j].fZ;
	  eta = AliHLTTransform::GetEta(xyz);
	  index2 = GetEtaIndex(eta);
	  if(index1 != index2) continue;
	  
	  AliHLTHistogram *hist = fParamSpace[index1];
	  
	  r2 = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]);
	  phi2 = atan2(xyz[1],xyz[0]);
	  
	  psi = atan( (r2*sin(phi1) - r1*sin(phi2)) / (r2*cos(phi1) - r1*cos(phi2)) );
	  kappa = 2*sin(phi1-psi) / r1;
	  hist->Fill(kappa,psi,1);
#ifdef do_mc
	  Int_t bin = hist->FindBin(kappa,psi);
	  for(Int_t l=0; l<3; l++)
	    {
	      for(Int_t m=0; m<3; m++)
		{
		  if(fClusters[i].fTrackID[l] == fClusters[j].fTrackID[m])
		    {
		      Int_t label = fClusters[i].fTrackID[l];
		      if(label < 0) continue;
		      UInt_t c;
		      for(c=0; c<MaxTrack; c++)
			if(fTrackID[index1][bin].fLabel[c] == label || fTrackID[index1][bin].fNHits[c] == 0)
			  break;
		      if(c == MaxTrack-1) cerr<<"AliHLTHoughClusterTransformer::TransformCircleC : Array reached maximum!! "<<c<<endl;
		      fTrackID[index1][bin].fLabel[c] = label;
		      fTrackID[index1][bin].fNHits[c]++;
		    }
		}
	    }
#endif
	}
    }
}


#ifdef do_mc
Int_t AliHLTHoughClusterTransformer::GetTrackID(Int_t etaindex,Double_t kappa,Double_t psi) const
{
  //Returns the mc label for a given bin in the hough space
  if(etaindex < 0 || etaindex > GetNEtaSegments())
    {
      cerr<<"AliHLTHoughClusterTransformer::GetTrackID : Wrong etaindex "<<etaindex<<endl;
      return -1;
    }
  AliHLTHistogram *hist = fParamSpace[etaindex];
  Int_t bin = hist->FindBin(kappa,psi);
  Int_t label=-1;
  Int_t max=0;
  for(UInt_t i=0; i<MaxTrack; i++)
    {
      Int_t nhits=fTrackID[etaindex][bin].fNHits[i];
      if(nhits == 0) break;
      if(nhits > max)
	{
	  max = nhits;
	  label = fTrackID[etaindex][bin].fLabel[i];
	}
    }
  return label;
#else
  Int_t AliHLTHoughClusterTransformer::GetTrackID(Int_t /*etaindex*/,Double_t /*kappa*/,Double_t /*psi*/) const
{
  // Does nothing if do_mc undefinde
  cout<<"AliHLTHoughClusterTransformer::GetTrackID : Compile with do_mc flag!"<<endl;
  return -1;
#endif
}

