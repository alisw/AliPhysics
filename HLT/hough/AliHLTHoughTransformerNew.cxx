// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliHLTStandardIncludes.h"

#include "AliHLTLogging.h"
#include "AliHLTHoughTransformerNew.h"
#include "AliHLTMemHandler.h"
#include "AliHLTTransform.h"
#include "AliHLTDigitData.h"
#include "AliHLTHistogramAdaptive.h"
#include "AliHLTTrackArray.h"
#include "AliHLTHoughTrack.h"
#include "AliHLTSpacePointData.h"
#include "AliHLTVertex.h"
#include "AliHLTFitter.h"

#if __GNUC__ >= 3
using namespace std;
#endif

//_____________________________________________________________
// AliHLTHoughTransformerNew
//
// Hough transformation class
//

ClassImp(AliHLTHoughTransformerNew)

AliHLTHoughTransformerNew::AliHLTHoughTransformerNew()
{
  //default ctor
  fParamSpace3D=0;
}

AliHLTHoughTransformerNew::AliHLTHoughTransformerNew(Int_t slice,Int_t patch,Int_t netasegments) : AliHLTHoughTransformer(slice,patch,netasegments)
{
  //normal ctor
  fParamSpace3D=0;
}
AliHLTHoughTransformerNew::~AliHLTHoughTransformerNew()
{
  //dtor
  if(fParamSpace3D)
    delete fParamSpace3D;
}

void AliHLTHoughTransformerNew::CreateHistograms(Int_t nxbins,Float_t xlow,Float_t xup,
						Int_t nybins,Float_t ylow,Float_t yup,
						Int_t nzbins,Float_t zlow,Float_t zup)
{
  //Create the histogram which contain the hough space
  Char_t name[1024];
  sprintf(name,"paramspace_%p",(void*)this);
  fParamSpace3D = new TH3F(name,"",nxbins,xlow,xup,nybins,ylow,yup,nzbins,zlow,zup);
}

void AliHLTHoughTransformerNew::Reset()
{
  //Reset Hough space
  fParamSpace3D->Reset();
}



void AliHLTHoughTransformerNew::TransformLine(Int_t *rowrange,Float_t *phirange)
{
  //Hough Transform
  AliHLTDigitRowData *tempPt = GetDataPointer();
  if(!tempPt)
    {
      LOG(AliHLTLog::kError,"AliHLTHoughTransformerNew::TransformLine","Data")
	<<"No input data "<<ENDLOG;
      return;
    }
  
  TH3 *hist = fParamSpace3D;
  for(Int_t i=AliHLTTransform::GetFirstRow(GetPatch()); i<=AliHLTTransform::GetLastRow(GetPatch()); i++)
    {
      AliHLTDigitData *digPt = tempPt->fDigitData;
      if(i < rowrange[0])
	{
	  AliHLTMemHandler::UpdateRowPointer(tempPt);
	  continue;
	}
      else if(i > rowrange[1])
	break;
      if(i != (Int_t)tempPt->fRow)
	{
	  printf("AliHLTHoughTransform::TransformLine : Mismatching padrow numbering\n");
	  continue;
	}
      for(UInt_t j=0; j<tempPt->fNDigit; j++)
	{
	  UShort_t charge = digPt[j].fCharge;
	  UChar_t pad = digPt[j].fPad;
	  UShort_t time = digPt[j].fTime;
	  if(charge < GetLowerThreshold() || charge > GetUpperThreshold()) continue;
	  
	  Int_t sector,row;
	  Float_t xyz[3];
	  AliHLTTransform::Slice2Sector(GetSlice(),i,sector,row);
	  AliHLTTransform::Raw2Local(xyz,sector,row,(Int_t)pad,(Int_t)time);
	  
	  Float_t phi = AliHLTTransform::GetPhi(xyz);
	  if(phi < phirange[0] || phi > phirange[1])
	    continue;

	  xyz[0] = xyz[0] - AliHLTTransform::Row2X(rowrange[0]);
	  Float_t x = xyz[0] + AliHLTTransform::Row2X(rowrange[0]);
	  Float_t r = sqrt(x*x + xyz[1]*xyz[1]);
	  Float_t delta = atan(xyz[2]/r);
	  for(Int_t xbin=hist->GetXaxis()->GetFirst(); xbin<=hist->GetXaxis()->GetLast(); xbin++)
	    {
	      Float_t theta = hist->GetXaxis()->GetBinCenter(xbin);
	      Float_t rho = xyz[0]*cos(theta) + xyz[1]*sin(theta);
	      hist->Fill(theta,rho,delta,charge);
	    }
	}
      AliHLTMemHandler::UpdateRowPointer(tempPt);
    }
  
}

struct AliHLTDigit {
  Int_t fRow;//padrow index
  Int_t fCharge;//digits charge
  Float_t fY;//Y coordinate of the digit
  Float_t fZ;//Z coordinate of the digit
  AliHLTDigit *fNext;//pointer to the next digit
};

struct AliHLTEtaContainer {
  void *fFirst;//pointer to the first digit in a sequence
  void *fLast;//pointer to the last digit in a sequence
};

void AliHLTHoughTransformerNew::TransformLineC(Int_t *rowrange,Float_t *phirange)
{
  //Hough Transform
  AliHLTDigitRowData *tempPt = GetDataPointer();
  if(!tempPt)
    LOG(AliHLTLog::kError,"AliHLTHoughTransformer::TransformCircleC","Data")
      <<"No input data "<<ENDLOG;
  
  Int_t counter=0;
  for(Int_t i=AliHLTTransform::GetFirstRow(GetPatch()); i<=AliHLTTransform::GetLastRow(GetPatch()); i++)
    {
      counter += tempPt->fNDigit;
      AliHLTMemHandler::UpdateRowPointer(tempPt);
    }
  
  Int_t bound = (GetNEtaSegments()+1)*(AliHLTTransform::GetNRows(GetPatch())+1);
  AliHLTEtaContainer *etaPt = new AliHLTEtaContainer[bound];
  memset(etaPt,0,bound*sizeof(AliHLTEtaContainer));  
  
  AliHLTDigit *digits = new AliHLTDigit[counter];
  cout<<"Allocating "<<counter*sizeof(AliHLTDigit)<<" bytes to digitsarray"<<endl;
  memset(digits,0,counter*sizeof(AliHLTDigit));

  Int_t sector,row;
  Float_t xyz[3];
  
  counter=0;
  tempPt = GetDataPointer();
  
  cout<<"Calculating digits in patch "<<GetPatch()<<endl;
  for(Int_t i=AliHLTTransform::GetFirstRow(GetPatch()); i<=AliHLTTransform::GetLastRow(GetPatch()); i++)
    {
      AliHLTDigitData *digPt = tempPt->fDigitData;
      for(UInt_t di=0; di<tempPt->fNDigit; di++)
	{
	  Int_t charge = digPt[di].fCharge;
	  Int_t pad = digPt[di].fPad;
	  Int_t time = digPt[di].fTime;
	  AliHLTTransform::Slice2Sector(GetSlice(),i,sector,row);
	  AliHLTTransform::Raw2Local(xyz,sector,row,(Int_t)pad,(Int_t)time);
	  Double_t eta = AliHLTTransform::GetEta(xyz);
	  
	  Float_t phi = atan2(xyz[1],xyz[0]);
	  if(phi < phirange[0] || phi > phirange[1]) continue;
	  
	  digits[counter].fRow = i;
	  digits[counter].fY = xyz[1];
	  digits[counter].fZ = xyz[2];
	  digits[counter].fCharge = charge;
	  
	  Int_t etaindex = GetEtaIndex(eta);
	  Int_t index = (GetNEtaSegments()+1)*(i-AliHLTTransform::GetFirstRow(GetPatch())) + etaindex;
	  
	  if(index > 0 && index < bound) 
	    {
	      if(etaPt[index].fFirst == 0)
		etaPt[index].fFirst = (void*)(&digits[counter]);
	      else
		((AliHLTDigit*)(etaPt[index].fLast))->fNext = &digits[counter];
	      etaPt[index].fLast = (void*)(&digits[counter]);
	    }
	  counter++;
	}
      AliHLTMemHandler::UpdateRowPointer(tempPt);
    }
  
  cout<<"Doing the combinatorics"<<endl;
  
  AliHLTDigit *dPt1,*dPt2;
  TH3 *hist = fParamSpace3D;
  for(Int_t e=0; e<GetNEtaSegments(); e++)
    {
      for(Int_t i=rowrange[0]; i<=rowrange[1]; i++)
	{
	  Int_t index1 = (GetNEtaSegments()+1)*(i-AliHLTTransform::GetFirstRow(GetPatch())) + e;
	  
	  for(dPt1 = (AliHLTDigit*)etaPt[index1].fFirst; dPt1 != 0; dPt1 = (AliHLTDigit*)dPt1->fNext)
	    {
	      for(Int_t j=i+1; j<=rowrange[1]; j++)
		{
		  Int_t index2 = (GetNEtaSegments()+1)*(j-AliHLTTransform::GetFirstRow(GetPatch())) + e;
		  
		  for(dPt2 = (AliHLTDigit*)etaPt[index2].fFirst; dPt2 != 0; dPt2 = (AliHLTDigit*)dPt2->fNext)
		    {
		      if(dPt1->fRow == dPt2->fRow)
			{
			  cerr<<"same row; indexes "<<index1<<" "<<index2<<endl;
			  exit(5);
			}
		      
		      //Do the transform:
		      Float_t x1 = AliHLTTransform::Row2X(dPt1->fRow) - AliHLTTransform::Row2X(rowrange[0]);
		      Float_t x2 = AliHLTTransform::Row2X(dPt2->fRow) - AliHLTTransform::Row2X(rowrange[0]);
		      Float_t y1 = dPt1->fY;
		      Float_t y2 = dPt2->fY;
		      Float_t theta = atan2(x2-x1,y1-y2);
		      Float_t rho = x1*cos(theta)+y1*sin(theta);
		      Float_t r1 = sqrt(pow(AliHLTTransform::Row2X(dPt1->fRow),2) + pow(y1,2));
		      Float_t delta1 = atan2(dPt1->fZ,r1);
		      Float_t r2 = sqrt(pow(AliHLTTransform::Row2X(dPt2->fRow),2) + pow(y2,2));
		      Float_t delta2 = atan2(dPt2->fZ,r2);
		      Float_t delta = (dPt1->fCharge*delta1+dPt2->fCharge*delta2)/(dPt1->fCharge+dPt2->fCharge);
		      hist->Fill(theta,rho,delta,dPt1->fCharge+dPt2->fCharge);
		    }
		}
	    }
	}
    }
  
  cout<<"done"<<endl;
  delete [] etaPt;
  delete [] digits;
  
}

