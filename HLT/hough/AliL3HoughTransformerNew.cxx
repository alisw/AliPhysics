// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliL3StandardIncludes.h"

#include "AliL3Logging.h"
#include "AliL3HoughTransformerNew.h"
#include "AliL3MemHandler.h"
#include "AliL3Transform.h"
#include "AliL3DigitData.h"
#include "AliL3HistogramAdaptive.h"
#include "AliL3TrackArray.h"
#include "AliL3HoughTrack.h"
#include "AliL3SpacePointData.h"
#include "AliL3Vertex.h"
#include "AliL3Fitter.h"

#if __GNUC__ == 3
using namespace std;
#endif

//_____________________________________________________________
// AliL3HoughTransformerNew
//
// Hough transformation class
//

ClassImp(AliL3HoughTransformerNew)

AliL3HoughTransformerNew::AliL3HoughTransformerNew()
{
  //default ctor
  fParamSpace3D=0;
}

AliL3HoughTransformerNew::AliL3HoughTransformerNew(Int_t slice,Int_t patch,Int_t netasegments) : AliL3HoughTransformer(slice,patch,netasegments)
{
  //normal ctor
  fParamSpace3D=0;
}
AliL3HoughTransformerNew::~AliL3HoughTransformerNew()
{
  //dtor
  if(fParamSpace3D)
    delete fParamSpace3D;
}

void AliL3HoughTransformerNew::CreateHistograms(Int_t nxbins,Float_t xlow,Float_t xup,
						Int_t nybins,Float_t ylow,Float_t yup,
						Int_t nzbins,Float_t zlow,Float_t zup)
{
  //Create the histogram which contain the hough space
  Char_t name[1024];
  sprintf(name,"paramspace_%d",(Int_t)this);
  fParamSpace3D = new TH3F(name,"",nxbins,xlow,xup,nybins,ylow,yup,nzbins,zlow,zup);
}

void AliL3HoughTransformerNew::Reset()
{
  //Reset Hough space
  fParamSpace3D->Reset();
}



void AliL3HoughTransformerNew::TransformLine(Int_t *rowrange,Float_t *phirange)
{
  //Hough Transform
  AliL3DigitRowData *tempPt = GetDataPointer();
  if(!tempPt)
    {
      LOG(AliL3Log::kError,"AliL3HoughTransformerNew::TransformLine","Data")
	<<"No input data "<<ENDLOG;
      return;
    }
  
  TH3 *hist = fParamSpace3D;
  for(Int_t i=AliL3Transform::GetFirstRow(GetPatch()); i<=AliL3Transform::GetLastRow(GetPatch()); i++)
    {
      AliL3DigitData *digPt = tempPt->fDigitData;
      if(i < rowrange[0])
	{
	  AliL3MemHandler::UpdateRowPointer(tempPt);
	  continue;
	}
      else if(i > rowrange[1])
	break;
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
	  if(charge < GetLowerThreshold() || charge > GetUpperThreshold()) continue;
	  
	  Int_t sector,row;
	  Float_t xyz[3];
	  AliL3Transform::Slice2Sector(GetSlice(),i,sector,row);
	  AliL3Transform::Raw2Local(xyz,sector,row,(Int_t)pad,(Int_t)time);
	  
	  Float_t phi = AliL3Transform::GetPhi(xyz);
	  if(phi < phirange[0] || phi > phirange[1])
	    continue;

	  xyz[0] = xyz[0] - AliL3Transform::Row2X(rowrange[0]);
	  Float_t x = xyz[0] + AliL3Transform::Row2X(rowrange[0]);
	  Float_t r = sqrt(x*x + xyz[1]*xyz[1]);
	  Float_t delta = atan(xyz[2]/r);
	  for(Int_t xbin=hist->GetXaxis()->GetFirst(); xbin<=hist->GetXaxis()->GetLast(); xbin++)
	    {
	      Float_t theta = hist->GetXaxis()->GetBinCenter(xbin);
	      Float_t rho = xyz[0]*cos(theta) + xyz[1]*sin(theta);
	      hist->Fill(theta,rho,delta,charge);
	    }
	}
      AliL3MemHandler::UpdateRowPointer(tempPt);
    }
  
}

struct AliL3Digit {
  Int_t fRow;//padrow index
  Int_t fCharge;//digits charge
  Float_t fY;//Y coordinate of the digit
  Float_t fZ;//Z coordinate of the digit
  AliL3Digit *fNext;//pointer to the next digit
};

struct AliL3EtaContainer {
  void *fFirst;//pointer to the first digit in a sequence
  void *fLast;//pointer to the last digit in a sequence
};

void AliL3HoughTransformerNew::TransformLineC(Int_t *rowrange,Float_t *phirange)
{
  //Hough Transform
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
  AliL3EtaContainer *etaPt = new AliL3EtaContainer[bound];
  memset(etaPt,0,bound*sizeof(AliL3EtaContainer));  
  
  AliL3Digit *digits = new AliL3Digit[counter];
  cout<<"Allocating "<<counter*sizeof(AliL3Digit)<<" bytes to digitsarray"<<endl;
  memset(digits,0,counter*sizeof(AliL3Digit));

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
	  
	  digits[counter].fRow = i;
	  digits[counter].fY = xyz[1];
	  digits[counter].fZ = xyz[2];
	  digits[counter].fCharge = charge;
	  
	  Int_t etaindex = GetEtaIndex(eta);
	  Int_t index = (GetNEtaSegments()+1)*(i-AliL3Transform::GetFirstRow(GetPatch())) + etaindex;
	  
	  if(index > 0 && index < bound) 
	    {
	      if(etaPt[index].fFirst == 0)
		etaPt[index].fFirst = (void*)(&digits[counter]);
	      else
		((AliL3Digit*)(etaPt[index].fLast))->fNext = &digits[counter];
	      etaPt[index].fLast = (void*)(&digits[counter]);
	    }
	  counter++;
	}
      AliL3MemHandler::UpdateRowPointer(tempPt);
    }
  
  cout<<"Doing the combinatorics"<<endl;
  
  AliL3Digit *dPt1,*dPt2;
  TH3 *hist = fParamSpace3D;
  for(Int_t e=0; e<GetNEtaSegments(); e++)
    {
      for(Int_t i=rowrange[0]; i<=rowrange[1]; i++)
	{
	  Int_t index1 = (GetNEtaSegments()+1)*(i-AliL3Transform::GetFirstRow(GetPatch())) + e;
	  
	  for(dPt1 = (AliL3Digit*)etaPt[index1].fFirst; dPt1 != 0; dPt1 = (AliL3Digit*)dPt1->fNext)
	    {
	      for(Int_t j=i+1; j<=rowrange[1]; j++)
		{
		  Int_t index2 = (GetNEtaSegments()+1)*(j-AliL3Transform::GetFirstRow(GetPatch())) + e;
		  
		  for(dPt2 = (AliL3Digit*)etaPt[index2].fFirst; dPt2 != 0; dPt2 = (AliL3Digit*)dPt2->fNext)
		    {
		      if(dPt1->fRow == dPt2->fRow)
			{
			  cerr<<"same row; indexes "<<index1<<" "<<index2<<endl;
			  exit(5);
			}
		      
		      //Do the transform:
		      Float_t x1 = AliL3Transform::Row2X(dPt1->fRow) - AliL3Transform::Row2X(rowrange[0]);
		      Float_t x2 = AliL3Transform::Row2X(dPt2->fRow) - AliL3Transform::Row2X(rowrange[0]);
		      Float_t y1 = dPt1->fY;
		      Float_t y2 = dPt2->fY;
		      Float_t theta = atan2(x2-x1,y1-y2);
		      Float_t rho = x1*cos(theta)+y1*sin(theta);
		      Float_t r1 = sqrt(pow(AliL3Transform::Row2X(dPt1->fRow),2) + pow(y1,2));
		      Float_t delta1 = atan2(dPt1->fZ,r1);
		      Float_t r2 = sqrt(pow(AliL3Transform::Row2X(dPt2->fRow),2) + pow(y2,2));
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

