//$Id$

// Author: Anders Vestbo <mailto:vestbo$fi.uib.no>
//*-- Copyright &copy ASV

#include "AliL3StandardIncludes.h"

#include "AliL3Logging.h"
#include "AliL3ModelTrack.h"
#include "AliL3Transform.h"

#if GCCVERSION == 3
using namespace std;
#endif

//_____________________________________________________________
// AliL3ModelTrack
//
// 

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
  fLabel=0;
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
  Int_t nrows = AliL3Transform::GetNRows(fPatch);
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
#ifdef do_mc
  for(Int_t i=0; i<nrows; i++)
    fClusters[i].fTrackID[0]=fClusters[i].fTrackID[1]=fClusters[i].fTrackID[2]=-2;
#endif
  fClusterCharge = 100;
  
  // 100 micrometers:
  fXYResidualQ = 0.01/AliL3Transform::GetPadPitchWidth(patch);
  fZResidualQ = 0.01/AliL3Transform::GetPadPitchWidth(patch);
  
  fXYWidthQ = 0.005/AliL3Transform::GetPadPitchWidth(patch);
  fZWidthQ = 0.005/AliL3Transform::GetPadPitchWidth(patch);
}


void AliL3ModelTrack::SetCluster(Int_t row,Float_t fpad,Float_t ftime,Float_t charge,Float_t sigmaY2,Float_t sigmaZ2,Int_t npads)
{
  Int_t index = row - AliL3Transform::GetFirstRow(fPatch);
  if(index != fNClusters)
    cout<<"AliL3ModelTrack::SetCluster() : Mismatch ; index: "<<index<<" nclusters "<<fNClusters<<endl;
  
  if(index < 0 || index > AliL3Transform::GetNRows(fPatch))
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
      cl->fDSigmaY2 = (sigmaY2 - GetParSigmaY2(row))/fXYWidthQ;
      cl->fDSigmaZ2 = (sigmaZ2 - GetParSigmaZ2(row))/fZWidthQ;
      cl->fNPads = npads;
    }
  
  fNClusters++;
}

Int_t AliL3ModelTrack::CheckClustersQuality(UInt_t npads)
{

  //Check the quality of clusters,- remove clusters with less than
  //npads. 
  //Returns the number of good clusters left.

  Int_t count=0;

  for(Int_t i=AliL3Transform::GetFirstRow(fPatch); i<=AliL3Transform::GetLastRow(fPatch); i++)
    {
      AliL3ClusterModel *cl = GetClusterModel(i);
      if(cl->fNPads < npads)
	cl->fPresent = kFALSE;
      if(cl->fPresent)
	count++;
    }
  
  return count;
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
  SetPt((BFACT*AliL3Transform::GetBField())/fabs(GetKappa()));
  
  CalculateHelix();
  
  Float_t hit[3];
  Int_t sector,row;
  for(Int_t i=AliL3Transform::GetFirstRow(fPatch); i<=AliL3Transform::GetLastRow(fPatch); i++)
    {
      AliL3ClusterModel *cl = GetClusterModel(i);
      if(!cl) continue;
      GetCrossingPoint(i,hit);
      AliL3Transform::Slice2Sector(fSlice,i,sector,row);
      AliL3Transform::Local2Raw(hit,sector,row);
      SetPadHit(i,hit[1]);
      SetTimeHit(i,hit[2]);
    }
}


void AliL3ModelTrack::SetTrackID(Int_t row,Int_t *trackID)
{
#ifdef do_mc
  AliL3ClusterModel *cluster = GetClusterModel(row);
  cluster->fTrackID[0] = trackID[0];
  cluster->fTrackID[1] = trackID[1];
  cluster->fTrackID[2] = trackID[2];
  return;
#endif
  cerr<<"AliL3ModelTrack::SetTrackID : Compile with do_mc flag"<<endl;
}


void AliL3ModelTrack::SetPadHit(Int_t row,Float_t pad)
{
  Int_t index = row-AliL3Transform::GetFirstRow(fPatch);
  if(index < 0 || index > AliL3Transform::GetNRows(fPatch))
    {
      cerr<<"AliL3ModelTrack::SetPadHit() : Wrong index: "<<index<<endl;
      return;
    }
  fPad[index]=pad;
  
}

void AliL3ModelTrack::SetTimeHit(Int_t row,Float_t time)
{
  Int_t index = row-AliL3Transform::GetFirstRow(fPatch);
  if(index < 0 || index > AliL3Transform::GetNRows(fPatch))
    {
      cerr<<"AliL3ModelTrack::SetTimeHit() : Wrong index: "<<index<<endl;
      return;
    }
  fTime[index]=time;
}

void AliL3ModelTrack::SetOverlap(Int_t row,Int_t id)
{
  Int_t index = row-AliL3Transform::GetFirstRow(fPatch);
  if(index < 0 || index > AliL3Transform::GetNRows(fPatch))
    {
      cerr<<"AliL3ModelTrack::SetOverlap() : Wrong index: "<<index<<endl;
      return;
    }
  fOverlap[index]=id;
}


Int_t AliL3ModelTrack::GetTrackID(Int_t row,Int_t index)
{
  
#ifdef do_mc
  AliL3ClusterModel *cl = GetClusterModel(row);
  return cl->fTrackID[index];
#endif
  cerr<<"AliL3ModelTrack::GetTrackID : Compile with do_mc flag"<<endl;
}


Int_t AliL3ModelTrack::GetNPads(Int_t row)
{
  AliL3ClusterModel *cl = GetClusterModel(row);
  return cl->fNPads;
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
  width = cl->fDSigmaY2*fXYWidthQ + GetParSigmaY2(row);
  
  return (Bool_t)cl->fPresent;
}

Bool_t AliL3ModelTrack::GetZWidth(Int_t row,Float_t &width)
{
  AliL3ClusterModel *cl = GetClusterModel(row);
  width = cl->fDSigmaZ2*fZWidthQ + GetParSigmaZ2(row);
  
  return (Bool_t)cl->fPresent;
}

Bool_t AliL3ModelTrack::GetPadResidual(Int_t row,Float_t &res)
{
  AliL3ClusterModel *cl = GetClusterModel(row);
  res = cl->fDPad;
  return cl->fPresent;
}

Bool_t AliL3ModelTrack::GetTimeResidual(Int_t row,Float_t &res)
{
  AliL3ClusterModel *cl = GetClusterModel(row);
  res = cl->fDTime;
  return cl->fPresent;
}

Bool_t AliL3ModelTrack::GetXYWidthResidual(Int_t row,Float_t &res)
{
  AliL3ClusterModel *cl = GetClusterModel(row);
  res = cl->fDSigmaY2;
  return cl->fPresent;
}

Bool_t AliL3ModelTrack::GetZWidthResidual(Int_t row,Float_t &res)
{
  AliL3ClusterModel *cl = GetClusterModel(row);
  res = cl->fDSigmaZ2;
  return cl->fPresent;
}

Float_t AliL3ModelTrack::GetPadHit(Int_t row)
{
  Int_t index = row-AliL3Transform::GetFirstRow(fPatch);
  if(index < 0 || index > AliL3Transform::GetNRows(fPatch))
    {
      cerr<<"AliL3ModelTrack::GetPadHit() : Wrong index: "<<index<<" row "<<row<<endl;
      return 0;
    }
  return fPad[index];
}

Float_t AliL3ModelTrack::GetTimeHit(Int_t row)
{
  Int_t index = row-AliL3Transform::GetFirstRow(fPatch);
  if(index < 0 || index > AliL3Transform::GetNRows(fPatch))
    {
      cerr<<"AliL3ModelTrack::GetTimeHit() : Wrong index: "<<index<<" row "<<row<<endl;
      return 0;
    }
  return fTime[index];
}

Int_t AliL3ModelTrack::GetOverlap(Int_t row)
{
  Int_t index = row-AliL3Transform::GetFirstRow(fPatch);
  if(index < 0 || index > AliL3Transform::GetNRows(fPatch))
    {
      cerr<<"AliL3ModelTrack::GetOverlap() : Wrong index: "<<index<<endl;
      return 0;
    }
  return fOverlap[index];
}

AliL3ClusterModel *AliL3ModelTrack::GetClusterModel(Int_t row)
{
  if(!fClusters) return 0; 
  Int_t index = row-AliL3Transform::GetFirstRow(fPatch);
  if(index < 0 || index > AliL3Transform::GetNRows(fPatch))
    {
      cerr<<"AliL3ModelTrack::GetClusterModel() : Wrong index: "<<index<<endl;
      return 0;
    }
  return &fClusters[index];
}

void AliL3ModelTrack::Print()
{
  //Print info

  cout<<"----Slice "<<fSlice<<" Patch "<<fPatch<<"----"<<endl;
  cout<<"First point "<<GetFirstPointX()<<" "<<GetFirstPointY()<<" "<<GetFirstPointZ()<<endl;
  cout<<"Last point "<<GetLastPointX()<<" "<<GetLastPointY()<<" "<<GetLastPointZ()<<endl;
  cout<<"Pt "<<GetPt()<<" kappa "<<GetKappa()<<" tgl "<<GetTgl()<<" psi "<<GetPsi()<<" charge "<<GetCharge()<<endl;
  cout<<"Center "<<GetCenterX()<<" "<<GetCenterY()<<endl<<endl;
  cout<<"NHits "<<GetNClusters()<<endl;
  cout<<"Clusters:"<<endl;

  for(Int_t i=AliL3Transform::GetFirstRow(fPatch); i<=AliL3Transform::GetLastRow(fPatch); i++)
    {
      AliL3ClusterModel *cl = GetClusterModel(i);
      
      if(!cl->fPresent)
	cout<<i<<" Empty"<<" Padcrossing "<<GetPadHit(i)<<" Timecrossing "<<GetTimeHit(i)<<" ";
      else
	{
	  cout<<i<<" Dpad "<<cl->fDPad<<" Dtime "<<cl->fDTime<<" Dcharge "<<cl->fDCharge;
	  cout<<" DsigmaY2 "<<cl->fDSigmaY2<<" DsigmaZ2 "<<cl->fDSigmaZ2;
	  cout<<" Padcrossing "<<GetPadHit(i)<<" Timecrossing "<<GetTimeHit(i)<<" ";
	  cout<<"Number of pads "<<GetNPads(i)<<endl;
	}
      cout<<"Overlapping index "<<GetOverlap(i)<<endl;
    }
}

Double_t AliL3ModelTrack::GetParSigmaY2(Int_t row)
{
  //Calculate the expected cluster width, based on the track parameters and drift distance.

  Float_t pad,time;
  if(!GetTime(row,time) || !GetPad(row,pad))
    return -1;
  
  Float_t xyz[3];
  Int_t sector,padrow;
  AliL3Transform::Slice2Sector(fSlice,row,sector,padrow);
  AliL3Transform::Raw2Local(xyz,sector,padrow,pad,time);
  
  //Calculate the drift length:
  Double_t drift;
  if(xyz[2] > 0)
    drift = AliL3Transform::GetZLength() - xyz[2];
  else
    drift = AliL3Transform::GetZLength() + xyz[2];
  
  Double_t prf = AliL3Transform::GetPRFSigma(fPatch);
  Double_t diffT = AliL3Transform::GetDiffT();
  Double_t padlength = AliL3Transform::GetPadLength(fPatch);
  Double_t anode = AliL3Transform::GetAnodeWireSpacing();
  Double_t beta = GetCrossingAngle(row);
  
  Double_t sigmaY2 = prf*prf + diffT*diffT*drift + padlength*padlength*tan(beta)*tan(beta)/12 + anode*anode*pow( tan(beta)-0.15, 2)/12;
  
  //Convert back to raw coordinates.
  sigmaY2 = sigmaY2/pow(AliL3Transform::GetPadPitchWidth(fPatch),2);
  return sigmaY2;
}

Double_t AliL3ModelTrack::GetParSigmaZ2(Int_t row)
{
  //Calculate the expected cluster width, based on the track parameters and drift distance.
  
  Float_t pad,time;
  if(!GetTime(row,time) || !GetPad(row,pad))
    return -1;
  
  Float_t xyz[3];
  Int_t sector,padrow;
  AliL3Transform::Slice2Sector(fSlice,row,sector,padrow);
  AliL3Transform::Raw2Local(xyz,sector,padrow,pad,time);
  
  //Calculate the drift length:
  Double_t drift;
  if(xyz[2] > 0)
    drift = AliL3Transform::GetZLength() - xyz[2];
  else
    drift = AliL3Transform::GetZLength() + xyz[2];
  
  Double_t sigma0 = AliL3Transform::GetTimeSigma();
  Double_t diffL = AliL3Transform::GetDiffL();
  Double_t padlength = AliL3Transform::GetPadLength(fPatch);
  Double_t tanl = GetTgl();
  
  Double_t sigmaZ2 = sigma0*sigma0 + diffL*diffL*drift + padlength*padlength * tanl*tanl/12;
  
  //Convert back to raw coodinates:
  sigmaZ2 = sigmaZ2/pow(AliL3Transform::GetZWidth(),2);
  return sigmaZ2;
  
}

void AliL3ModelTrack::AssignTrackID(Float_t wrong)
{
  //Assign a track ID to the track, corresponding to the MC TParticle ID.
  //Can only be done if you compiled with do_mc flag, of course.
  //The function loops over the assigned clusters, and finds the label (ID)
  //of each clusters, and assigns the ID with the most hits to the track.
  //If there are more than wrong% clusters of a different ID, the track is
  //considered to be fake, and label will be assigned as negative.
  
#ifdef do_mc
  Int_t *lb = new Int_t[GetNClusters()];
  Int_t *mx = new Int_t[GetNClusters()];

  Int_t i,j;
  for(Int_t i=0; i<GetNClusters(); i++) 
    lb[i]=mx[i]=0;
  
  Int_t lab=123456789;
  
  for(i=0; i<GetNClusters(); i++) 
    {
      lab = abs(GetTrackID(i,0));
      for (j=0; j<GetNClusters(); j++) 
	if (lb[j]==lab || mx[j]==0) break;
      lb[j]=lab;
      (mx[j])++;
    }
  
  Int_t max=0;
  for (i=0; i<GetNClusters(); i++) 
    {
      if(mx[i] > max) 
	{
	  max=mx[i]; 
	  lab=lb[i];
	}
    }
  
  for (i=0; i<GetNClusters(); i++) 
    {
      if(abs(GetTrackID(i,1)) == lab ||
	 abs(GetTrackID(i,2)) == lab)
	max++;
    }

  if ((1.- Float_t(max)/GetNClusters()) > wrong) lab=-lab;

  SetLabel(lab);

  delete[] lb;
  delete[] mx;
  return;
#endif
  cerr<<"AliL3ModelTrack::AssignTrackID : Compile with do_mc flag"<<endl;
}
