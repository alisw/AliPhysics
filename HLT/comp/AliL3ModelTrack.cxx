//$Id$

// Author: Anders Vestbo <mailto:vestbo$fi.uib.no>
//*-- Copyright &copy ASV

#include <stream.h>
#include <string.h>
#include <math.h>

#include "AliL3ModelTrack.h"
#include "AliL3Transform.h"

ClassImp(AliL3ModelTrack)

AliL3ModelTrack::AliL3ModelTrack()
{
  fNClusters = 0;
  fClusters = 0;
  fOverlap = 0;
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
  if(fOverlap)
    delete [] fOverlap;
  if(fTrackModel)
    delete fTrackModel;
}

void AliL3ModelTrack::Init(Int_t slice,Int_t patch)
{
  fNClusters = 0;
  fSlice=slice;
  fPatch=patch;
  Int_t nrows = NumRows[patch];
  fClusters = new AliL3ClusterModel[nrows];
  fPad = new Float_t[nrows];
  fTime = new Float_t[nrows];
  fTrackModel = new AliL3TrackModel;
  fOverlap = new Int_t[nrows];
  
  memset(fClusters,0,nrows*sizeof(AliL3ClusterModel));
  memset(fPad,0,nrows*sizeof(Float_t));
  memset(fTime,0,nrows*sizeof(Float_t));
  memset(fTrackModel,0,sizeof(AliL3TrackModel));
  for(Int_t i=0; i<nrows; i++)
    fOverlap[i]=-1;

  fClusterCharge = 260;
  AliL3Transform transform;
  
  fXYResidualQ = 0.1/transform.GetPadPitchWidth(patch);
  fZResidualQ = 0.1/transform.GetPadPitchWidth(patch);
  

  fXYWidthQ = 0.01;
  fZWidthQ = 0.01;
}


void AliL3ModelTrack::SetCluster(Int_t row,Float_t fpad,Float_t ftime,Float_t charge,Float_t sigmaY2,Float_t sigmaZ2)
{
  Int_t index = row - NRows[fPatch][0];
  if(index != fNClusters)
    cout<<"AliL3ModelTrack::SetCluster() : Mismatch ; index: "<<index<<" nclusters "<<fNClusters<<endl;
  
  if(index < 0 || index > NumRows[fPatch])
    {
      cerr<<"AliL3ModelTrack::SetCluster() : Wrong index: "<<index<<" row "<<row<<endl;
      return;
    }
  AliL3ClusterModel *cl = GetClusterModel(row);
  if(!charge)
    cl->fPresent = kFALSE;
  else
    {
      cl->fPresent = kTRUE;
      cl->fDTime = (ftime - GetTimeHit(row))/fXYResidualQ;
      cl->fDPad = (fpad - GetPadHit(row))/fZResidualQ;
      cl->fDCharge = charge - fClusterCharge;
      cl->fDSigmaY2 = sigmaY2/fXYWidthQ;
      cl->fDSigmaZ2 = sigmaZ2/fZWidthQ;
    }
  
  fNClusters++;
}


void AliL3ModelTrack::FillModel()
{
  //Fill the track structure
  
  if(!fTrackModel)
    {
      cerr<<"AliL3ModelTrack::FillModel() : No trackmodel "<<endl;
      return;
    }
  fTrackModel->fKappa = GetKappa();
  fTrackModel->fFirstPointX = GetFirstPointX();
  fTrackModel->fFirstPointY = GetFirstPointY();
  fTrackModel->fFirstPointZ = GetFirstPointZ();
  fTrackModel->fTgl = GetTgl();
  fTrackModel->fPsi = GetPsi();
  fTrackModel->fLength = (Short_t)GetLength();
  fTrackModel->fClusterCharge = fClusterCharge;
  fTrackModel->fNClusters = fNClusters;

}

void AliL3ModelTrack::FillTrack()
{
  //Fill the track parameters from the structure.
  
  if(!fTrackModel)
    {
      cerr<<"AliL3ModelTrack::FillTrack() : No data!!"<<endl;
      return;
    }
  SetKappa(fTrackModel->fKappa);
  SetCharge((-1*(Int_t)copysign(1.,GetKappa())));
  SetFirstPoint(fTrackModel->fFirstPointX,fTrackModel->fFirstPointY,fTrackModel->fFirstPointZ);
  SetTgl(fTrackModel->fTgl);
  SetPsi(fTrackModel->fPsi);
  SetLength(fTrackModel->fLength);
  fClusterCharge=fTrackModel->fClusterCharge;
  fNClusters = fTrackModel->fNClusters;
  SetPt((BFACT*BField)/fabs(GetKappa()));
  
  CalculateHelix();
  
  AliL3Transform transform;
  Float_t hit[3];
  Int_t sector,row;
  for(Int_t i=NRows[fPatch][0]; i<=NRows[fPatch][1]; i++)
    {
      AliL3ClusterModel *cl = GetClusterModel(i);
      if(!cl) continue;
      GetCrossingPoint(i,hit);
      transform.Slice2Sector(fSlice,i,sector,row);
      transform.Local2Raw(hit,sector,row);
      SetPadHit(i,hit[1]);
      SetTimeHit(i,hit[2]);
    }
}

void AliL3ModelTrack::SetPadHit(Int_t row,Float_t pad)
{
  Int_t index = row-NRows[fPatch][0];
  if(index < 0 || index > NumRows[fPatch])
    {
      cerr<<"AliL3ModelTrack::SetPadHit() : Wrong index: "<<index<<endl;
      return;
    }
  fPad[index]=pad;
  
}

void AliL3ModelTrack::SetTimeHit(Int_t row,Float_t time)
{
  Int_t index = row-NRows[fPatch][0];
  if(index < 0 || index > NumRows[fPatch])
    {
      cerr<<"AliL3ModelTrack::SetTimeHit() : Wrong index: "<<index<<endl;
      return;
    }
  fTime[index]=time;
}

void AliL3ModelTrack::SetOverlap(Int_t row,Int_t id)
{
  Int_t index = row-NRows[fPatch][0];
  if(index < 0 || index > NumRows[fPatch])
    {
      cerr<<"AliL3ModelTrack::SetOverlap() : Wrong index: "<<index<<endl;
      return;
    }
  fOverlap[index]=id;
}

Bool_t AliL3ModelTrack::GetPad(Int_t row,Float_t &pad)
{
  //(ftime - GetTimeHit(fNClusters))/fXYResidualQ;
  //(fpad - GetPadHit(fNClusters))/fZResidualQ;

  AliL3ClusterModel *cl = GetClusterModel(row);
  pad = cl->fDPad*fXYResidualQ + GetPadHit(row);

  return (Bool_t)cl->fPresent;
}

Bool_t AliL3ModelTrack::GetTime(Int_t row,Float_t &time)
{
  AliL3ClusterModel *cl = GetClusterModel(row);
  time = cl->fDTime*fZResidualQ + GetTimeHit(row);

  return (Bool_t)cl->fPresent;
}

Bool_t AliL3ModelTrack::GetClusterCharge(Int_t row,Int_t &charge)
{
  AliL3ClusterModel *cl = GetClusterModel(row);
  charge = (Int_t)cl->fDCharge + fClusterCharge;
  
  return (Bool_t)cl->fPresent;
}

Bool_t AliL3ModelTrack::GetXYWidth(Int_t row,Float_t &width)
{
  AliL3ClusterModel *cl = GetClusterModel(row);
  width = cl->fDSigmaY2*fXYWidthQ;
  
  return (Bool_t)cl->fPresent;
}

Bool_t AliL3ModelTrack::GetZWidth(Int_t row,Float_t &width)
{
  AliL3ClusterModel *cl = GetClusterModel(row);
  width = cl->fDSigmaZ2*fZWidthQ;
  
  return (Bool_t)cl->fPresent;
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

Float_t AliL3ModelTrack::GetPadHit(Int_t row)
{
  Int_t index = row-NRows[fPatch][0];
  if(index < 0 || index > NumRows[fPatch])
    {
      cerr<<"AliL3ModelTrack::GetPadHit() : Wrong index: "<<index<<" row "<<row<<endl;
      return 0;
    }
  return fPad[index];
}

Float_t AliL3ModelTrack::GetTimeHit(Int_t row)
{
  Int_t index = row-NRows[fPatch][0];
  if(index < 0 || index > NumRows[fPatch])
    {
      cerr<<"AliL3ModelTrack::GetTimeHit() : Wrong index: "<<index<<" row "<<row<<endl;
      return 0;
    }
  return fTime[index];
}

Int_t AliL3ModelTrack::GetOverlap(Int_t row)
{
  Int_t index = row-NRows[fPatch][0];
  if(index < 0 || index > NumRows[fPatch])
    {
      cerr<<"AliL3ModelTrack::GetOverlap() : Wrong index: "<<index<<endl;
      return 0;
    }
  return fOverlap[index];
}

AliL3ClusterModel *AliL3ModelTrack::GetClusterModel(Int_t row)
{
  if(!fClusters) return 0; 
  Int_t index = row-NRows[fPatch][0];
  if(index < 0 || index > NumRows[fPatch])
    {
      cerr<<"AliL3ModelTrack::GetClusterModel() : Wrong index: "<<index<<endl;
      return 0;
    }
  return &fClusters[index];
}

void AliL3ModelTrack::Print()
{
  //Print info

  cout<<"---------------------"<<endl;
  cout<<"First point "<<GetFirstPointX()<<" "<<GetFirstPointY()<<" "<<GetFirstPointZ()<<endl;
  cout<<"Last point "<<GetLastPointX()<<" "<<GetLastPointY()<<" "<<GetLastPointZ()<<endl;
  cout<<"Pt "<<GetPt()<<" kappa "<<GetKappa()<<" tgl "<<GetTgl()<<" psi "<<GetPsi()<<" charge "<<GetCharge()<<endl;
  cout<<"Center "<<GetCenterX()<<" "<<GetCenterY()<<endl<<endl;
  cout<<"NHits "<<GetNClusters()<<endl;
  cout<<"Clusters:"<<endl;

  for(Int_t i=NRows[fPatch][0]; i<=NRows[fPatch][1]; i++)
    {
      AliL3ClusterModel *cl = GetClusterModel(i);
      
      if(!cl->fPresent)
	cout<<i<<" Empty"<<" Padcrossing "<<GetPadHit(i)<<" Timecrossing "<<GetTimeHit(i)<<" ";
      else
	{
	  cout<<i<<" Dpad "<<cl->fDPad<<" Dtime "<<cl->fDTime<<" Dcharge "<<cl->fDCharge;
	  cout<<" Padcrossing "<<GetPadHit(i)<<" Timecrossing "<<GetTimeHit(i)<<" ";
	}
      cout<<"Overlapping index "<<GetOverlap(i)<<endl;
    }
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
