//$Id$

// Author: Anders Vestbo <mailto:vestbo$fi.uib.no>
//*-- Copyright &copy ASV

#include <stream.h>
#include <string.h>
#include <math.h>

#include "AliL3ModelTrack.h"
#include "AliL3Defs.h"
#include "AliL3Transform.h"

ClassImp(AliL3ModelTrack)

AliL3ModelTrack::AliL3ModelTrack()
{
  fNClusters = 0;
  fClusters = 0;
  fOverlap = -1;
  fPad=0;
  fTime=0;
  fClusterCharge=0;
  fTrackModel=0;
}


AliL3ModelTrack::~AliL3ModelTrack()
{
  if(fClusters)
    delete [] fClusters;
  if(fPad)
    delete [] fPad;
  if(fTime)
    delete [] fTime;
  if(fTrackModel)
    delete fTrackModel;
}

void AliL3ModelTrack::Init(Int_t slice,Int_t patch)
{
  fNClusters = 0;
  Int_t nrows = NumRows[patch];
  fClusters = new AliL3ClusterModel[nrows];
  memset((void*)fClusters,0,nrows*sizeof(AliL3ClusterModel));
  
  fPad = new Float_t[NRowsSlice];
  fTime = new Float_t[NRowsSlice];
  fTrackModel = new AliL3TrackModel;
  memset(fTrackModel,0,sizeof(AliL3TrackModel));
  
  fClusterCharge = 100;
  AliL3Transform transform;
  
  fXYResidualQ = 0.01/transform.GetPadPitchWidth(patch);
  fZResidualQ = 0.01/transform.GetPadPitchWidth(patch);
}


void AliL3ModelTrack::SetCluster(Float_t fpad,Float_t ftime,Float_t charge,Float_t sigmaY2,Float_t sigmaZ2)
{

  AliL3ClusterModel *cl = &fClusters[fNClusters];
  if(!charge)
    cl->fPresent = kFALSE;
  else
    {
      cl->fPresent = kTRUE;
      cl->fDTime = (ftime - GetTimeHit(fNClusters))/fXYResidualQ;
      cl->fDPad = (fpad - GetPadHit(fNClusters))/fZResidualQ;
      cl->fDCharge = charge;
      cl->fDSigmaY2 = sigmaY2;
      cl->fDSigmaZ2 = sigmaZ2;
    }
  //cout<<"Pad "<<fpad<<" dtime "<<ftime<<" charge "<<charge<<" sigmaY2 "<<sigmaY2<<" sigmaZ2 "<<sigmaZ2<<endl;
  fNClusters++;
}



void AliL3ModelTrack::FillModel()
{
  //fTrackModel = new AliL3TrackModel;
  fTrackModel->fKappa = GetKappa();
  fTrackModel->fFirstPointX = GetFirstPointX();
  fTrackModel->fFirstPointY = GetFirstPointY();
  fTrackModel->fFirstPointZ = GetFirstPointZ();
  fTrackModel->fTgl = GetTgl();
  fTrackModel->fPsi = GetPsi();
  fTrackModel->fLength = GetLength();
  fTrackModel->fClusterCharge = fClusterCharge;
  fTrackModel->fNClusters = fNClusters;

}

void AliL3ModelTrack::Print()
{
  //Print info

  cout<<"---------------------"<<endl;
  for(Int_t i=0; i<fNClusters; i++)
    {
      AliL3ClusterModel *cl = &fClusters[i];
      if(!cl->fPresent)
	cout<<i<<" Empty"<<endl;
      else
	{
	  cout<<i<<" Dpad "<<cl->fDPad<<" Dtime "<<cl->fDTime<<" Dcharge "<<cl->fDCharge;
	  cout<<" Padcrossing "<<GetPadHit(i)<<" Timecrossing "<<GetTimeHit(i)<<endl;
	}
    }
}

Bool_t AliL3ModelTrack::GetPadResidual(Int_t row,Float_t &res)
{
  res = fClusters[row].fDPad;
  return fClusters[row].fPresent;
}

Bool_t AliL3ModelTrack::GetTimeResidual(Int_t row,Float_t &res)
{
  res = fClusters[row].fDTime;
  return fClusters[row].fPresent;
}

//----------Code below taken from AliTPCTracker.cxx-----------------------
//Functions that give the expected cluster errors based on track parameters.
Double_t AliL3ModelTrack::GetParSigmaY2(Double_t r)//, Double_t tgl, Double_t pt)
{
  
  //
  // Parametrised error of the cluster reconstruction (pad direction)   
  //
  // Sigma rphi
  
  Double_t tgl = GetTgl();
  Double_t pt = GetPt();
  
  const Float_t kArphi=0.41818e-2;
  const Float_t kBrphi=0.17460e-4;
  const Float_t kCrphi=0.30993e-2;
  const Float_t kDrphi=0.41061e-3;
  
  pt=fabs(pt)*1000.;
  Double_t x=r/pt;
  tgl=fabs(tgl);
  Double_t s=kArphi - kBrphi*r*tgl + kCrphi*x*x + kDrphi*x;
  if (s<0.4e-3) s=0.4e-3;
  s*=1.3; //Iouri Belikov

  return s;
}

Double_t AliL3ModelTrack::GetParSigmaZ2(Double_t r)//, Double_t tgl) 
{
  //
  // Parametrised error of the cluster reconstruction (drift direction)
  //
  // Sigma z
  
  Double_t tgl = GetTgl();

  const Float_t kAz=0.39614e-2;
  const Float_t kBz=0.22443e-4;
  const Float_t kCz=0.51504e-1;
  

  tgl=fabs(tgl);
  Double_t s=kAz - kBz*r*tgl + kCz*tgl*tgl;
  if (s<0.4e-3) s=0.4e-3;
  s*=1.3; //Iouri Belikov

  return s;
}
