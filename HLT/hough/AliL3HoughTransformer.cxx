//Author:        Anders Strand Vestbo
//Last Modified: 14.09.01

#include "AliL3Logging.h"
#include "AliL3HoughTransformer.h"
#include "AliL3Defs.h"
#include "AliL3Transform.h"
#include "AliL3DigitData.h"
#include "AliL3Histogram.h"

ClassImp(AliL3HoughTransformer)

AliL3HoughTransformer::AliL3HoughTransformer()
{
  //Default constructor
  
}

AliL3HoughTransformer::AliL3HoughTransformer(Int_t slice,Int_t patch,Int_t n_eta_segments)
{
  fSlice = slice;
  fPatch = patch;
  fNEtaSegments = n_eta_segments;
  fEtaMin = 0;
  fEtaMax = fSlice < 18 ? 0.9 : -0.9;
  fTransform = new AliL3Transform();
  fThreshold = 0;
  
}

AliL3HoughTransformer::~AliL3HoughTransformer()
{
  if(fTransform)
    delete fTransform;
  DeleteHistograms();
}

void AliL3HoughTransformer::DeleteHistograms()
{
  if(!fParamSpace)
    return;
  for(Int_t i=0; i<fNEtaSegments; i++)
    {
      if(!fParamSpace[i]) continue;
      delete fParamSpace[i];
    }
  delete [] fParamSpace;
}

void AliL3HoughTransformer::CreateHistograms(Int_t nxbin,Double_t pt_min,
					     Int_t nybin,Double_t phimin,Double_t phimax)
{
  //Set the minimum absolute pt value, and phi0 angles given in degrees.

  Double_t torad = 3.1415/180;
  Double_t bfact = 0.0029980;
  Double_t bfield = 0.2;
  Double_t x = bfact*bfield/pt_min;
  CreateHistograms(nxbin,-1.*x,x,nybin,phimin*torad,phimax*torad);
}

void AliL3HoughTransformer::CreateHistograms(Int_t nxbin,Double_t xmin,Double_t xmax,
					     Int_t nybin,Double_t ymin,Double_t ymax)
{
  
  fParamSpace = new AliL3Histogram*[fNEtaSegments];
  
  Char_t histname[256];
  for(Int_t i=0; i<fNEtaSegments; i++)
    {
      sprintf(histname,"paramspace_%d",i);
      fParamSpace[i] = new AliL3Histogram(histname,"",nxbin,xmin,xmax,nybin,ymin,ymax);
    }
}

void AliL3HoughTransformer::Reset()
{
  //Reset all the histograms, and memory

  if(!fParamSpace)
    {
      LOG(AliL3Log::kWarning,"AliL3HoughTransformer::Reset","Histograms")
	<<"No histograms to reset"<<ENDLOG;
      return;
    }
  
  for(Int_t i=0; i<fNEtaSegments; i++)
    fParamSpace[i]->Reset();
  fNDigitRowData = 0;
  fDigitRowData = 0;
}

void AliL3HoughTransformer::SetInputData(UInt_t ndigits,AliL3DigitRowData *ptr)
{
  fNDigitRowData = ndigits;
  fDigitRowData = ptr;
}

void AliL3HoughTransformer::UpdateDataPointer(AliL3DigitRowData *&tempPt)
{
  //Update the data pointer to the next padrow
  
  Byte_t *tmp = (Byte_t*)tempPt;
  Int_t size = sizeof(AliL3DigitRowData) + tempPt->fNDigit*sizeof(AliL3DigitData);
  tmp += size;
  tempPt = (AliL3DigitRowData*)tmp;
}

void AliL3HoughTransformer::TransformCircle()
{
  //Transform the input data with a circle HT.
  

  //Set pointer to the data
  AliL3DigitRowData *tempPt = (AliL3DigitRowData*)fDigitRowData;
  if(!tempPt || fNDigitRowData==0)
    {
      printf("\nAliL3HoughTransformer::TransformCircle : No input data!!!\n\n");
      return;
    }
  
  Double_t etaslice = (fEtaMax - fEtaMin)/fNEtaSegments;
  for(Int_t i=NRows[fPatch][0]; i<=NRows[fPatch][1]; i++)
    {
      AliL3DigitData *digPt = tempPt->fDigitData;
      if(i != (Int_t)tempPt->fRow)
	{
	  printf("AliL3HoughTransform::TransformCircle : Mismatching padrow numbering\n");
	  continue;
	}
      for(UInt_t j=0; j<tempPt->fNDigit; j++)
	{
	  UShort_t charge = digPt[j].fCharge;
	  UChar_t pad = digPt[j].fPad;
	  UShort_t time = digPt[j].fTime;
	  if(charge < fThreshold)
	    continue;
	  Int_t sector,row;
	  Float_t xyz[3];
	  fTransform->Slice2Sector(fSlice,i,sector,row);
	  fTransform->Raw2Local(xyz,sector,row,(Int_t)pad,(Int_t)time);
	  Double_t eta = fTransform->GetEta(xyz);
	  Int_t eta_index = (Int_t)(eta/etaslice);
	  if(eta_index < 0 || eta_index >= fNEtaSegments)
	    continue;
	  
	  //Get the correct histogrampointer:
	  AliL3Histogram *hist = fParamSpace[eta_index];
	  if(!hist)
	    {
	      printf("AliL3HoughTransformer::TransformCircle : Error getting histogram in index %d\n",eta_index);
	      continue;
	    }

	  //Start transformation
	  Float_t R = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]); // + xyz[2]*xyz[2]);
	  Float_t phi = fTransform->GetPhi(xyz);
	  
	  //Fill the histogram along the phirange
	  for(Int_t b=hist->GetFirstYbin(); b<=hist->GetLastYbin(); b++)
	    {
	      Float_t phi0 = hist->GetBinCenterY(b);
	      Float_t kappa = 2*sin(phi - phi0)/R;
	      hist->Fill(kappa,phi0,charge);
	    }
	}
      UpdateDataPointer(tempPt);
    }
}

void AliL3HoughTransformer::TransformLine()
{
  //Do a line transform on the data.

  
  AliL3DigitRowData *tempPt = (AliL3DigitRowData*)fDigitRowData;
  if(!tempPt || fNDigitRowData==0)
    {
      printf("\nAliL3HoughTransformer::TransformLine : No input data!!!\n\n");
      return;
    }
  
  Double_t etaslice = (fEtaMax - fEtaMin)/fNEtaSegments;
  for(Int_t i=NRows[fPatch][0]; i<=NRows[fPatch][1]; i++)
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
	  if(charge < fThreshold)
	    continue;
	  Int_t sector,row;
	  Float_t xyz[3];
	  fTransform->Slice2Sector(fSlice,i,sector,row);
	  fTransform->Raw2Local(xyz,sector,row,(Int_t)pad,(Int_t)time);
	  Float_t eta = fTransform->GetEta(xyz);
	  Int_t eta_index = (Int_t)(eta/etaslice);
	  if(eta_index < 0 || eta_index >= fNEtaSegments)
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
      UpdateDataPointer(tempPt);
    }
  
}
