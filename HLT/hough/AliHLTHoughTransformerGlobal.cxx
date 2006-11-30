// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliHLTStandardIncludes.h"

#include "AliHLTLogging.h"
#include "AliHLTHoughTransformerGlobal.h"
#include "AliHLTFileHandler.h"
#include "AliHLTTransform.h"
#include "AliHLTDigitData.h"
#include "AliHLTTrack.h"
#include "AliHLTHistogram.h"
#include "AliHLTTrackArray.h"
#include "AliHLTHoughMaxFinder.h"
#include "AliHLTHoughTrack.h"
#include "AliHLTSpacePointData.h"

#include <TH2.h>

#if __GNUC__ >= 3
using namespace std;
#endif

//_____________________________________________________________
// AliHLTHoughTransformerGlobal
//
// Hough transformation class
//

ClassImp(AliHLTHoughTransformerGlobal)

AliHLTHoughTransformerGlobal::AliHLTHoughTransformerGlobal()
{
  //default ctor
  fPadMin=0;
  fPadMax=0;
  fNActiveSlice=0;
  fMemHandler=0;
  fEvent=0;
  fPsi=0;
  fPtMin=0;
  fTracks=0;
  fPeakFinder=0;
  fSeedPadRow=-1;
}

AliHLTHoughTransformerGlobal::AliHLTHoughTransformerGlobal(Char_t *path,Int_t event)
{
  //normal ctor
  strcpy(fPath,path);
  fEvent=event;
  fMemHandler=0;
  fPadMin=0;
  fPadMax=0;
  fPsi=0;
  fPtMin=0;
  fNActiveSlice=0;
  fTracks = new AliHLTTrackArray("AliHLTHoughTrack");
  fPeakFinder = new AliHLTHoughMaxFinder("KappaPhi",1000);
}

AliHLTHoughTransformerGlobal::~AliHLTHoughTransformerGlobal()
{
  //dtor
  if(fPadMin)
    delete [] fPadMin;
  if(fPadMax)
    delete [] fPadMax;
  if(fTracks)
    delete fTracks;
  if(fPeakFinder)
    delete fPeakFinder;
  UnloadActiveSlices();
  
}

void AliHLTHoughTransformerGlobal::CreateHistograms(Float_t /*ptmin*/,Int_t nxbin,Int_t nybin)
{
  //create histograms for HT
  if(fPsi==0)
    {
      cerr<<"AliHLTHoughTransformerGlobal::CreateHistograms : Call DefineRegion first"<<endl;
      exit(5);
    }
  AliHLTHoughTransformer::CreateHistograms(nxbin,fPtMin,nybin,-fPsi,fPsi);
  //AliHLTHoughTransformer::CreateHistograms(fPtMin,ptmax,ptres,nybin,fPsi);
}

void AliHLTHoughTransformerGlobal::TransformCircleC()
{
  //Hough Transform
  AliHLTSeed *clusters = new AliHLTSeed[1000];
  Int_t nclusters = LoadClusterSeeds(clusters);
  
  Int_t i=0,sector,row,slice;
  UInt_t dummy=0;
  Float_t xyz[3],rpe[3],kappa,psi;
  
  for(Int_t sl=GetSlice()-fNActiveSlice; sl<=GetSlice()+fNActiveSlice; sl++)
    {
      if(sl < 0 || sl < 18 && GetSlice() >= 18)
	slice = sl + 18;
      else if(sl > 35 || sl > 17 && GetSlice() <= 17)
	slice = sl-18;
      else
	slice = sl;
      cout<<"Transforming in slice "<<slice<<endl;
      
      AliHLTDigitRowData *rowPt = (AliHLTDigitRowData*)fMemHandler[i++]->GetDataPointer(dummy);
      
      for(Int_t padrow=0; padrow<AliHLTTransform::GetNRows(); padrow++)
	{
	  if(padrow == fSeedPadRow) continue;
	  AliHLTDigitData *digits = (AliHLTDigitData*)rowPt->fDigitData;
	  
	  for(UInt_t j=0; j<rowPt->fNDigit; j++)
	    {
	      Int_t pad = digits[j].fPad;

	      if(i==1 && pad < fPadMin[padrow])
		continue;
	      if(i==fNActiveSlice*2+1 && pad>fPadMax[padrow])
		continue;

	      Int_t time = digits[j].fTime;
	      Int_t charge = digits[j].fCharge;
	      
	      if(charge > GetUpperThreshold() || charge < GetLowerThreshold())
		continue;
	      
	      AliHLTTransform::Slice2Sector(slice,padrow,sector,row);
	      AliHLTTransform::Raw2Local(xyz,sector,row,pad,time);
	      
	      Rotate(xyz,i-1-fNActiveSlice);
	      
	      AliHLTTransform::XYZtoRPhiEta(rpe,xyz);
	      rpe[0] = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]);
	      
	      Int_t etaindex = GetEtaIndex(rpe[2]);
	      
	      AliHLTHistogram *hist = GetHistogram(etaindex);
	      
	      for(Int_t l=0; l<nclusters; l++)
		{
		  if(clusters[l].fIndex != etaindex) continue;
		  psi = atan( (clusters[l].fRadius*sin(rpe[1])-rpe[0]*sin(clusters[l].fPhi))/
			      (clusters[l].fRadius*cos(rpe[1])-rpe[0]*cos(clusters[l].fPhi)) );
		  kappa = 2*sin(clusters[l].fPhi-psi)/clusters[l].fRadius;
		  if(fabs(kappa) < 1e-33)
  		    continue;
		  //hist->Fill(kappa,psi,(int)rint(log((Float_t)charge)));
		  hist->Fill(kappa,psi,charge);
		}
	    }
	  AliHLTMemHandler::UpdateRowPointer(rowPt);
	}
    }
  delete [] clusters;
}

void AliHLTHoughTransformerGlobal::TransformCircle()
{
  //Hough Transform  
  Int_t i=0,sector,row,slice;
  UInt_t dummy=0;
  Float_t xyz[3],rpe[3],kappa,psi;
  
  for(Int_t sl=GetSlice()-fNActiveSlice; sl<=GetSlice()+fNActiveSlice; sl++)
    {
      if(sl < 0 || sl < 18 && GetSlice() >= 18)
	slice = sl + 18;
      else if(sl > 35 || sl > 17 && GetSlice() <= 17)
	slice = sl-18;
      else
	slice = sl;
      cout<<"Transforming in slice "<<slice<<endl;
      
      AliHLTDigitRowData *rowPt = (AliHLTDigitRowData*)fMemHandler[i++]->GetDataPointer(dummy);
      
      for(Int_t padrow=0; padrow<AliHLTTransform::GetNRows(); padrow++)
	{

	  AliHLTDigitData *digits = (AliHLTDigitData*)rowPt->fDigitData;
	  
	  for(UInt_t j=0; j<rowPt->fNDigit; j++)
	    {
	      Int_t pad = digits[j].fPad;
	      
	      if(i==1 && pad < fPadMin[padrow])
		continue;
	      if(i==fNActiveSlice*2+1 && pad>fPadMax[padrow])
		continue;
	      
	      Int_t time = digits[j].fTime;
	      Int_t charge = digits[j].fCharge;
	      
	      if(charge > GetUpperThreshold() || charge < GetLowerThreshold())
		continue;
	      
	      AliHLTTransform::Slice2Sector(slice,padrow,sector,row);
	      AliHLTTransform::Raw2Local(xyz,sector,row,pad,time);
	      
	      Rotate(xyz,i-1-fNActiveSlice);
	      
	      AliHLTTransform::XYZtoRPhiEta(rpe,xyz);
	      rpe[0] = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]);
	      
	      Int_t etaindex = GetEtaIndex(rpe[2]);
	      
	      AliHLTHistogram *hist = GetHistogram(etaindex);
	      
	      for(Int_t b=hist->GetFirstYbin(); b<=hist->GetLastYbin(); b++)
		{
		  psi = hist->GetBinCenterY(b);
		  kappa = 2*sin(rpe[1] - psi)/rpe[0];
		  
		  if(fabs(kappa) < 1e-33)
		    continue;
		  
		  hist->Fill(kappa,psi,charge);
		  //hist->Fill(kappa,psi,(int)rint(log((Float_t)charge)));
		}
	    }
	  AliHLTMemHandler::UpdateRowPointer(rowPt);
	}
    }
}

void AliHLTHoughTransformerGlobal::VerifyTracks(AliHLTTrackArray *tracks,Int_t &index)
{
  //remove tracks which do not cross the seed padrow
  for(int i=index; i<tracks->GetNTracks(); i++)
    {
      AliHLTHoughTrack *tr = (AliHLTHoughTrack*)tracks->GetCheckedTrack(i);
      if(!tr) continue;

      Float_t angle=0;
      AliHLTTransform::Local2GlobalAngle(&angle,GetSlice());
      if(!tr->CalculateReferencePoint(angle,AliHLTTransform::Row2X(fSeedPadRow)))
	{
	  tracks->Remove(i);
	  continue;
	}
      Int_t sector,row;
      AliHLTTransform::Slice2Sector(GetSlice(),fSeedPadRow,sector,row);
      float xyz[3] = {tr->GetPointX(),tr->GetPointY(),tr->GetPointZ()};
      AliHLTTransform::Global2Raw(xyz,sector,row);
      if(xyz[1]<0 || xyz[1]>=AliHLTTransform::GetNPads(fSeedPadRow) ||
	 xyz[2]<0 || xyz[2]>=AliHLTTransform::GetNTimeBins())
	tracks->Remove(i);
    }
  tracks->Compress();
  index = tracks->GetNTracks();
}

void AliHLTHoughTransformerGlobal::FindPeaks(AliHLTHistogram *hist,Float_t eta)
{
  //Peak Finder  
  fPeakFinder->Reset();
  fPeakFinder->SetHistogram(hist);
  fPeakFinder->FindAbsMaxima();
  if(fPeakFinder->GetWeight(0) < 1000)
    return;
  
  AliHLTHoughTrack *track = (AliHLTHoughTrack*)fTracks->NextTrack();
  track->SetTrackParameters(fPeakFinder->GetXPeak(0),fPeakFinder->GetYPeak(0),fPeakFinder->GetWeight(0));
  track->SetEta(eta);
  track->SetRowRange(AliHLTTransform::GetFirstRow(0),AliHLTTransform::GetLastRow(5));
}

void AliHLTHoughTransformerGlobal::Rotate(Float_t *xyz,Int_t relslice)
{
  //Rotate coordinates from one slice to the slice relative to it;
  //-1 means lower
  //1 means upper

  Float_t angle = (Float_t)relslice*AliHLTTransform::Deg2Rad(20);
  Float_t x=xyz[0],y=xyz[1];
  xyz[0] = x*cos(angle) - y*sin(angle);
  xyz[1] = x*sin(angle) + y*cos(angle);
}


void AliHLTHoughTransformerGlobal::DefineRegion(Float_t minpt,Float_t /*linephi*/,Int_t seedpadrow)
{
  //Setup the region to be included in the transform
  //This function should be called only once, since it is the same for all slices.
  
  fSeedPadRow=seedpadrow;
  fPtMin = minpt;
  
  fPadMin = new Int_t[AliHLTTransform::GetNRows()];
  fPadMax = new Int_t[AliHLTTransform::GetNRows()];
  
  //Phirange of data which should be included in the transform
  //Here we assume a min pt to define it; fPtMin
  //Based on pt (kappa) and the two points; origo and middle point, we can get psi angle
  
  //Calculate the upper point, the lower point is symmetrical around the x-axis
  //The track in this calculation has a positive charge, and the kappa sign is the opposite:  
  
  Float_t xyz[3],phi,phi2;
  
  phi = CalculateBorder(xyz,-1); //upward bending track
  phi2 = CalculateBorder(xyz,1); //downward bending track
  if(phi2 > phi)
    phi = phi2;
  
  cout<<"Phiangle "<<phi*180/AliHLTTransform::Pi()<<" psi "<<fPsi*180/AliHLTTransform::Pi()<<" nslices "<<fNActiveSlice<<endl;
  
  Float_t rotangle = fNActiveSlice*20;
  
  //Calculate the LUT for min/max pad for every padrow, and check which slices we need.
  Int_t pad,sector,row;
  for(Int_t i=0; i<AliHLTTransform::GetNRows(); i++)
    {
      AliHLTTransform::Slice2Sector(0,i,sector,row);
      
      //Lower boundary:
      pad = AliHLTTransform::GetNPads(i)-1;
      while(pad >= 0)
	{
	  AliHLTTransform::Raw2Local(xyz,sector,row,pad,0);
	  if(AliHLTTransform::GetPhi(xyz) > -1.*phi + AliHLTTransform::Deg2Rad(rotangle))
	    fPadMin[i] = pad;
	  else
	    break;
	  pad--;
	}
      
      //Upper boundary
      pad = 0;
      while(pad < AliHLTTransform::GetNPads(i))
	{
	  AliHLTTransform::Raw2Local(xyz,sector,row,pad,0);
	  if(AliHLTTransform::GetPhi(xyz) < phi - AliHLTTransform::Deg2Rad(rotangle))
	    fPadMax[i] = pad;
	  else
	    break;
	  pad++;
	}
    }
  
  cout<<"Padmax "<<fPadMax[155]<<endl;

}

Float_t AliHLTHoughTransformerGlobal::CalculateBorder(Float_t *xyz,Int_t charge)
{
  //Define Hough space borders
  Double_t lineradius = sqrt(pow(AliHLTTransform::Row2X(fSeedPadRow),2) + pow(AliHLTTransform::GetMaxY(fSeedPadRow),2));

  Double_t kappa = -charge*AliHLTTransform::GetBField()*AliHLTTransform::GetBFact()/fPtMin;

  //Calculate the psi angle of the track = emission angle with x-axis in origo
  if( fabs(lineradius*kappa/2) > 1)
    {
      cerr<<"AliHLTTransformerGlobal::DefineRegion : Angle too big"<<endl;
      exit(5);
    }
  
  Int_t nslices=0;
  Float_t psi = AliHLTTransform::Deg2Rad(10) - asin(lineradius*kappa/2);
  if(charge > 0)
    fPsi=psi;
  
  cout<<"Calculated psi-angle "<<fPsi<<endl;
  
  AliHLTTrack track;
  track.SetFirstPoint(0,0,0);
  track.SetPt(fPtMin);
  track.SetPsi(psi);
  track.SetCharge(charge);
  track.CalculateHelix();
  
  Int_t crossingrow=0;
  if(charge < 0)
    crossingrow = AliHLTTransform::GetNRows()-1;
  
  Float_t rotangle;
  
 redefine:
  
  rotangle = nslices*20;
  
  if(nslices==0)//here we are in local slice, so we use the appropriate function. a mess of course...
    {
      while(!track.GetCrossingPoint(crossingrow,xyz))
	{
	  if(crossingrow==0)
	    {	  
	      cerr<<"AliHLTHoughTransformerGlobal::DefineRegion : Error calculating point1 on row "<<crossingrow<<endl;
	      exit(5);
	    }
	  crossingrow--;
	}
    }
  else
    {
      while(!track.CalculateReferencePoint(AliHLTTransform::Deg2Rad(rotangle),AliHLTTransform::Row2X(crossingrow)))
	{
	  if(crossingrow==0)
	    {
	      cerr<<"AliHLTHoughTransformerGlobal::DefineRegion : Error calculating point2 on row "<<crossingrow<<endl;
	      exit(5);
	    }
	  crossingrow--;
	}
      
      xyz[0] = track.GetPointX();
      xyz[1] = track.GetPointY();
      xyz[2] = track.GetPointZ();
    }
  
  Float_t phi = atan2(xyz[1],xyz[0]);
  
  //Rotate coordinates back to local coordinates:
  Rotate(xyz,-1*nslices);
  
  Int_t sector,row;
  AliHLTTransform::Slice2Sector(0,crossingrow,sector,row);
  AliHLTTransform::Local2Raw(xyz,sector,row);
  if(xyz[1] < 0)
    {
      if(crossingrow>0)
	{
	  cerr<<"AliHLTHoughTransformerGlobal::DefineRegion : Wrong pad, probably a deadzone... "<<xyz[1]<<endl;
	  exit(5);
	}
      else
	return 0;//here you only want the crossing point with the outer padrow
    }
  if(xyz[1] >= AliHLTTransform::GetNPads(crossingrow)) //Here the range is even one more slice away
    {
      cerr<<"One more slice with pad "<<xyz[1]<<endl;
      nslices++;
      if(charge < 0)
	crossingrow = AliHLTTransform::GetNRows()-1;
      goto redefine;
    }
  
  if(nslices > fNActiveSlice)
    fNActiveSlice = nslices;
  
  cout<<"Calculated phi "<<phi<<" and pad "<<xyz[1]<<" on row "<<crossingrow<<endl;
  return phi;
}

void AliHLTHoughTransformerGlobal::LoadActiveSlices(Bool_t binary)
{
  //Load active slices
  if(fMemHandler)
    UnloadActiveSlices();
  fMemHandler = new AliHLTMemHandler*[(fNActiveSlice*2 + 1)];
  Char_t filename[1024];
  UInt_t dummy;
  Int_t i=0,slice;
  for(Int_t sl=GetSlice()-fNActiveSlice; sl<=GetSlice()+fNActiveSlice; sl++)
    {
      if(sl < 0 || sl < 18 && GetSlice() >= 18)
	slice = sl + 18;
      else if(sl > 35 || sl > 17 && GetSlice() <= 17)
	slice = sl-18;
      else
	slice = sl;
      cout<<"Loading slice "<<slice<<endl;
      fMemHandler[i] = new AliHLTFileHandler();
      fMemHandler[i]->Init(slice,-1);
      if(binary)
	{
	  sprintf(filename,"%s/binaries/digits_%d_%d_%d.raw",fPath,fEvent,slice,-1);
	  fMemHandler[i]->SetBinaryInput(filename);
	  fMemHandler[i]->CompBinary2Memory(dummy);
	  fMemHandler[i++]->CloseBinaryInput();
	}
      else
	{
	  sprintf(filename,"%s/digitfile.root",fPath);
	  fMemHandler[i]->SetAliInput(filename);
	  fMemHandler[i]->AliAltroDigits2Memory(dummy,fEvent);
	  fMemHandler[i++]->CloseAliInput();
	}
    }
  
}

void AliHLTHoughTransformerGlobal::UnloadActiveSlices()
{
  //Unload active slices
  if(!fMemHandler)
    return;
  for(Int_t i=0; i<=fNActiveSlice*2; i++)
    {
      if(!fMemHandler[i]) continue;
      delete fMemHandler[i];
    }
  delete [] fMemHandler;
  fMemHandler=0;
}

Int_t AliHLTHoughTransformerGlobal::LoadClusterSeeds(AliHLTSeed *seeds)
{
  //Load cluster seeds
  Char_t filename[1024];
  UInt_t npoints;
  sprintf(filename,"%s/hough/points_%d_%d_%d.raw",fPath,fEvent,GetSlice(),-1);
  //sprintf(filename,"../tracker/points_%d_%d_%d.raw",fEvent,GetSlice(),-1);
  AliHLTMemHandler mem;
  if(!mem.SetBinaryInput(filename))
    {
      cerr<<"AliHLTHoughTransformerGlobal::LoadClusterSeeds : Cannot open file "<<filename<<endl;
      exit(5);
    }
  AliHLTSpacePointData *points = (AliHLTSpacePointData*)mem.Allocate();
  mem.Binary2Memory(npoints,points);
  mem.CloseBinaryInput();
  
  Int_t counter=0;
  for(UInt_t i=0; i<npoints; i++)
    {
      Int_t lrow = (Int_t)points[i].fPadRow;
      if(lrow < fSeedPadRow) continue;
      if(lrow > fSeedPadRow) break;
      if(fabs(points[i].fY) > AliHLTTransform::GetMaxY(lrow))
	{
	  cerr<<"AliHLTHoughTransformerGlobal::LoadClusterSeeds : Seeds seems not to be in local coordinates: "
	      <<points[i].fY<<endl;
	  exit(5);
	}
      Float_t xyz[3] = {AliHLTTransform::Row2X(lrow),points[i].fY,points[i].fZ};
      Double_t eta = AliHLTTransform::GetEta(xyz);
      seeds[counter].fPhi = atan2(xyz[1],xyz[0]);
      seeds[counter].fRadius = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]);
      seeds[counter].fIndex = GetEtaIndex(eta);
      seeds[counter++].fEta = eta;
    }
  cout<<"Loaded "<<counter<<" cluster seeds on slice "<<GetSlice()<<" padrow "<<fSeedPadRow<<endl;
  mem.Free();
  return counter;
}

void AliHLTHoughTransformerGlobal::DisplayActiveRegion(TH2F *hist,Int_t etaindex)
{
  //Fill the active region in a histogram
  
  Int_t i=0,sector,row,slice;
  UInt_t dummy=0;
  Float_t xyz[3];
  for(Int_t sl=GetSlice()-fNActiveSlice; sl<=GetSlice()+fNActiveSlice; sl++)
    {
      if(sl < 0 || sl < 18 && GetSlice() >= 18)
	slice = sl + 18;
      else if(sl > 35 || sl > 17 && GetSlice() <= 17)
	slice = sl-18;
      else
	slice = sl;
      cout<<"Displaying slice "<<slice<<endl;
      AliHLTDigitRowData *rowPt = (AliHLTDigitRowData*)fMemHandler[i++]->GetDataPointer(dummy);
      for(Int_t padrow=0; padrow<AliHLTTransform::GetNRows(); padrow++)
	{
	  AliHLTDigitData *digits = (AliHLTDigitData*)rowPt->fDigitData;
	  for(UInt_t j=0; j<rowPt->fNDigit; j++)
	    {
	      Int_t pad = digits[j].fPad;
	      if(i==1 && pad < fPadMin[padrow])
 		continue;
	      if(i==fNActiveSlice*2+1 && pad>fPadMax[padrow])
		continue;
	      Int_t time = digits[j].fTime;
	      AliHLTTransform::Slice2Sector(slice,padrow,sector,row);
	      AliHLTTransform::Raw2Local(xyz,sector,row,pad,time);
	     
	      Rotate(xyz,i-1-fNActiveSlice);
	      
	      Double_t eta = AliHLTTransform::GetEta(xyz);
	      if(GetEtaIndex(eta) == etaindex)
		hist->Fill(xyz[0],xyz[1]);
	    }
	  AliHLTMemHandler::UpdateRowPointer(rowPt);
	}
    }
}
