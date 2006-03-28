// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliL3StandardIncludes.h"

#include "AliL3Logging.h"
#include "AliL3HoughTransformerGlobal.h"
#include "AliL3FileHandler.h"
#include "AliL3Transform.h"
#include "AliL3DigitData.h"
#include "AliL3Track.h"
#include "AliL3Histogram.h"
#include "AliL3TrackArray.h"
#include "AliL3HoughMaxFinder.h"
#include "AliL3HoughTrack.h"
#include "AliL3SpacePointData.h"

#include <TH2.h>

#if __GNUC__ >= 3
using namespace std;
#endif

//_____________________________________________________________
// AliL3HoughTransformerGlobal
//
// Hough transformation class
//

ClassImp(AliL3HoughTransformerGlobal)

AliL3HoughTransformerGlobal::AliL3HoughTransformerGlobal()
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

AliL3HoughTransformerGlobal::AliL3HoughTransformerGlobal(Char_t *path,Int_t event)
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
  fTracks = new AliL3TrackArray("AliL3HoughTrack");
  fPeakFinder = new AliL3HoughMaxFinder("KappaPhi",1000);
}

AliL3HoughTransformerGlobal::~AliL3HoughTransformerGlobal()
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

void AliL3HoughTransformerGlobal::CreateHistograms(Float_t /*ptmin*/,Int_t nxbin,Int_t nybin)
{
  //create histograms for HT
  if(fPsi==0)
    {
      cerr<<"AliL3HoughTransformerGlobal::CreateHistograms : Call DefineRegion first"<<endl;
      exit(5);
    }
  AliL3HoughTransformer::CreateHistograms(nxbin,fPtMin,nybin,-fPsi,fPsi);
  //AliL3HoughTransformer::CreateHistograms(fPtMin,ptmax,ptres,nybin,fPsi);
}

void AliL3HoughTransformerGlobal::TransformCircleC()
{
  //Hough Transform
  AliL3Seed *clusters = new AliL3Seed[1000];
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
      
      AliL3DigitRowData *rowPt = (AliL3DigitRowData*)fMemHandler[i++]->GetDataPointer(dummy);
      
      for(Int_t padrow=0; padrow<AliL3Transform::GetNRows(); padrow++)
	{
	  if(padrow == fSeedPadRow) continue;
	  AliL3DigitData *digits = (AliL3DigitData*)rowPt->fDigitData;
	  
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
	      
	      AliL3Transform::Slice2Sector(slice,padrow,sector,row);
	      AliL3Transform::Raw2Local(xyz,sector,row,pad,time);
	      
	      Rotate(xyz,i-1-fNActiveSlice);
	      
	      AliL3Transform::XYZtoRPhiEta(rpe,xyz);
	      rpe[0] = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]);
	      
	      Int_t etaindex = GetEtaIndex(rpe[2]);
	      
	      AliL3Histogram *hist = GetHistogram(etaindex);
	      
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
	  AliL3MemHandler::UpdateRowPointer(rowPt);
	}
    }
  delete [] clusters;
}

void AliL3HoughTransformerGlobal::TransformCircle()
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
      
      AliL3DigitRowData *rowPt = (AliL3DigitRowData*)fMemHandler[i++]->GetDataPointer(dummy);
      
      for(Int_t padrow=0; padrow<AliL3Transform::GetNRows(); padrow++)
	{

	  AliL3DigitData *digits = (AliL3DigitData*)rowPt->fDigitData;
	  
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
	      
	      AliL3Transform::Slice2Sector(slice,padrow,sector,row);
	      AliL3Transform::Raw2Local(xyz,sector,row,pad,time);
	      
	      Rotate(xyz,i-1-fNActiveSlice);
	      
	      AliL3Transform::XYZtoRPhiEta(rpe,xyz);
	      rpe[0] = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]);
	      
	      Int_t etaindex = GetEtaIndex(rpe[2]);
	      
	      AliL3Histogram *hist = GetHistogram(etaindex);
	      
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
	  AliL3MemHandler::UpdateRowPointer(rowPt);
	}
    }
}

void AliL3HoughTransformerGlobal::VerifyTracks(AliL3TrackArray *tracks,Int_t &index)
{
  //remove tracks which do not cross the seed padrow
  for(int i=index; i<tracks->GetNTracks(); i++)
    {
      AliL3HoughTrack *tr = (AliL3HoughTrack*)tracks->GetCheckedTrack(i);
      if(!tr) continue;

      Float_t angle=0;
      AliL3Transform::Local2GlobalAngle(&angle,GetSlice());
      if(!tr->CalculateReferencePoint(angle,AliL3Transform::Row2X(fSeedPadRow)))
	{
	  tracks->Remove(i);
	  continue;
	}
      Int_t sector,row;
      AliL3Transform::Slice2Sector(GetSlice(),fSeedPadRow,sector,row);
      float xyz[3] = {tr->GetPointX(),tr->GetPointY(),tr->GetPointZ()};
      AliL3Transform::Global2Raw(xyz,sector,row);
      if(xyz[1]<0 || xyz[1]>=AliL3Transform::GetNPads(fSeedPadRow) ||
	 xyz[2]<0 || xyz[2]>=AliL3Transform::GetNTimeBins())
	tracks->Remove(i);
    }
  tracks->Compress();
  index = tracks->GetNTracks();
}

void AliL3HoughTransformerGlobal::FindPeaks(AliL3Histogram *hist,Float_t eta)
{
  //Peak Finder  
  fPeakFinder->Reset();
  fPeakFinder->SetHistogram(hist);
  fPeakFinder->FindAbsMaxima();
  if(fPeakFinder->GetWeight(0) < 1000)
    return;
  
  AliL3HoughTrack *track = (AliL3HoughTrack*)fTracks->NextTrack();
  track->SetTrackParameters(fPeakFinder->GetXPeak(0),fPeakFinder->GetYPeak(0),fPeakFinder->GetWeight(0));
  track->SetEta(eta);
  track->SetRowRange(AliL3Transform::GetFirstRow(0),AliL3Transform::GetLastRow(5));
}

void AliL3HoughTransformerGlobal::Rotate(Float_t *xyz,Int_t relslice)
{
  //Rotate coordinates from one slice to the slice relative to it;
  //-1 means lower
  //1 means upper

  Float_t angle = (Float_t)relslice*AliL3Transform::Deg2Rad(20);
  Float_t x=xyz[0],y=xyz[1];
  xyz[0] = x*cos(angle) - y*sin(angle);
  xyz[1] = x*sin(angle) + y*cos(angle);
}


void AliL3HoughTransformerGlobal::DefineRegion(Float_t minpt,Float_t /*linephi*/,Int_t seedpadrow)
{
  //Setup the region to be included in the transform
  //This function should be called only once, since it is the same for all slices.
  
  fSeedPadRow=seedpadrow;
  fPtMin = minpt;
  
  fPadMin = new Int_t[AliL3Transform::GetNRows()];
  fPadMax = new Int_t[AliL3Transform::GetNRows()];
  
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
  
  cout<<"Phiangle "<<phi*180/AliL3Transform::Pi()<<" psi "<<fPsi*180/AliL3Transform::Pi()<<" nslices "<<fNActiveSlice<<endl;
  
  Float_t rotangle = fNActiveSlice*20;
  
  //Calculate the LUT for min/max pad for every padrow, and check which slices we need.
  Int_t pad,sector,row;
  for(Int_t i=0; i<AliL3Transform::GetNRows(); i++)
    {
      AliL3Transform::Slice2Sector(0,i,sector,row);
      
      //Lower boundary:
      pad = AliL3Transform::GetNPads(i)-1;
      while(pad >= 0)
	{
	  AliL3Transform::Raw2Local(xyz,sector,row,pad,0);
	  if(AliL3Transform::GetPhi(xyz) > -1.*phi + AliL3Transform::Deg2Rad(rotangle))
	    fPadMin[i] = pad;
	  else
	    break;
	  pad--;
	}
      
      //Upper boundary
      pad = 0;
      while(pad < AliL3Transform::GetNPads(i))
	{
	  AliL3Transform::Raw2Local(xyz,sector,row,pad,0);
	  if(AliL3Transform::GetPhi(xyz) < phi - AliL3Transform::Deg2Rad(rotangle))
	    fPadMax[i] = pad;
	  else
	    break;
	  pad++;
	}
    }
  
  cout<<"Padmax "<<fPadMax[155]<<endl;

}

Float_t AliL3HoughTransformerGlobal::CalculateBorder(Float_t *xyz,Int_t charge)
{
  //Define Hough space borders
  Double_t lineradius = sqrt(pow(AliL3Transform::Row2X(fSeedPadRow),2) + pow(AliL3Transform::GetMaxY(fSeedPadRow),2));

  Double_t kappa = -charge*AliL3Transform::GetBField()*AliL3Transform::GetBFact()/fPtMin;

  //Calculate the psi angle of the track = emission angle with x-axis in origo
  if( fabs(lineradius*kappa/2) > 1)
    {
      cerr<<"AliL3TransformerGlobal::DefineRegion : Angle too big"<<endl;
      exit(5);
    }
  
  Int_t nslices=0;
  Float_t psi = AliL3Transform::Deg2Rad(10) - asin(lineradius*kappa/2);
  if(charge > 0)
    fPsi=psi;
  
  cout<<"Calculated psi-angle "<<fPsi<<endl;
  
  AliL3Track track;
  track.SetFirstPoint(0,0,0);
  track.SetPt(fPtMin);
  track.SetPsi(psi);
  track.SetCharge(charge);
  track.CalculateHelix();
  
  Int_t crossingrow=0;
  if(charge < 0)
    crossingrow = AliL3Transform::GetNRows()-1;
  
  Float_t rotangle;
  
 redefine:
  
  rotangle = nslices*20;
  
  if(nslices==0)//here we are in local slice, so we use the appropriate function. a mess of course...
    {
      while(!track.GetCrossingPoint(crossingrow,xyz))
	{
	  if(crossingrow==0)
	    {	  
	      cerr<<"AliL3HoughTransformerGlobal::DefineRegion : Error calculating point1 on row "<<crossingrow<<endl;
	      exit(5);
	    }
	  crossingrow--;
	}
    }
  else
    {
      while(!track.CalculateReferencePoint(AliL3Transform::Deg2Rad(rotangle),AliL3Transform::Row2X(crossingrow)))
	{
	  if(crossingrow==0)
	    {
	      cerr<<"AliL3HoughTransformerGlobal::DefineRegion : Error calculating point2 on row "<<crossingrow<<endl;
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
  AliL3Transform::Slice2Sector(0,crossingrow,sector,row);
  AliL3Transform::Local2Raw(xyz,sector,row);
  if(xyz[1] < 0)
    {
      if(crossingrow>0)
	{
	  cerr<<"AliL3HoughTransformerGlobal::DefineRegion : Wrong pad, probably a deadzone... "<<xyz[1]<<endl;
	  exit(5);
	}
      else
	return 0;//here you only want the crossing point with the outer padrow
    }
  if(xyz[1] >= AliL3Transform::GetNPads(crossingrow)) //Here the range is even one more slice away
    {
      cerr<<"One more slice with pad "<<xyz[1]<<endl;
      nslices++;
      if(charge < 0)
	crossingrow = AliL3Transform::GetNRows()-1;
      goto redefine;
    }
  
  if(nslices > fNActiveSlice)
    fNActiveSlice = nslices;
  
  cout<<"Calculated phi "<<phi<<" and pad "<<xyz[1]<<" on row "<<crossingrow<<endl;
  return phi;
}

void AliL3HoughTransformerGlobal::LoadActiveSlices(Bool_t binary)
{
  //Load active slices
  if(fMemHandler)
    UnloadActiveSlices();
  fMemHandler = new AliL3MemHandler*[(fNActiveSlice*2 + 1)];
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
      fMemHandler[i] = new AliL3FileHandler();
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

void AliL3HoughTransformerGlobal::UnloadActiveSlices()
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

Int_t AliL3HoughTransformerGlobal::LoadClusterSeeds(AliL3Seed *seeds)
{
  //Load cluster seeds
  Char_t filename[1024];
  UInt_t npoints;
  sprintf(filename,"%s/hough/points_%d_%d_%d.raw",fPath,fEvent,GetSlice(),-1);
  //sprintf(filename,"../tracker/points_%d_%d_%d.raw",fEvent,GetSlice(),-1);
  AliL3MemHandler mem;
  if(!mem.SetBinaryInput(filename))
    {
      cerr<<"AliL3HoughTransformerGlobal::LoadClusterSeeds : Cannot open file "<<filename<<endl;
      exit(5);
    }
  AliL3SpacePointData *points = (AliL3SpacePointData*)mem.Allocate();
  mem.Binary2Memory(npoints,points);
  mem.CloseBinaryInput();
  
  Int_t counter=0;
  for(UInt_t i=0; i<npoints; i++)
    {
      Int_t lrow = (Int_t)points[i].fPadRow;
      if(lrow < fSeedPadRow) continue;
      if(lrow > fSeedPadRow) break;
      if(fabs(points[i].fY) > AliL3Transform::GetMaxY(lrow))
	{
	  cerr<<"AliL3HoughTransformerGlobal::LoadClusterSeeds : Seeds seems not to be in local coordinates: "
	      <<points[i].fY<<endl;
	  exit(5);
	}
      Float_t xyz[3] = {AliL3Transform::Row2X(lrow),points[i].fY,points[i].fZ};
      Double_t eta = AliL3Transform::GetEta(xyz);
      seeds[counter].fPhi = atan2(xyz[1],xyz[0]);
      seeds[counter].fRadius = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1]);
      seeds[counter].fIndex = GetEtaIndex(eta);
      seeds[counter++].fEta = eta;
    }
  cout<<"Loaded "<<counter<<" cluster seeds on slice "<<GetSlice()<<" padrow "<<fSeedPadRow<<endl;
  mem.Free();
  return counter;
}

void AliL3HoughTransformerGlobal::DisplayActiveRegion(TH2F *hist,Int_t etaindex)
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
      AliL3DigitRowData *rowPt = (AliL3DigitRowData*)fMemHandler[i++]->GetDataPointer(dummy);
      for(Int_t padrow=0; padrow<AliL3Transform::GetNRows(); padrow++)
	{
	  AliL3DigitData *digits = (AliL3DigitData*)rowPt->fDigitData;
	  for(UInt_t j=0; j<rowPt->fNDigit; j++)
	    {
	      Int_t pad = digits[j].fPad;
	      if(i==1 && pad < fPadMin[padrow])
 		continue;
	      if(i==fNActiveSlice*2+1 && pad>fPadMax[padrow])
		continue;
	      Int_t time = digits[j].fTime;
	      AliL3Transform::Slice2Sector(slice,padrow,sector,row);
	      AliL3Transform::Raw2Local(xyz,sector,row,pad,time);
	     
	      Rotate(xyz,i-1-fNActiveSlice);
	      
	      Double_t eta = AliL3Transform::GetEta(xyz);
	      if(GetEtaIndex(eta) == etaindex)
		hist->Fill(xyz[0],xyz[1]);
	    }
	  AliL3MemHandler::UpdateRowPointer(rowPt);
	}
    }
}
