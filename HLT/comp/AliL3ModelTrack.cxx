// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo$fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliL3StandardIncludes.h"

#include "AliL3Logging.h"
#include "AliL3ModelTrack.h"
#include "AliL3Transform.h"
#include "AliL3DataCompressor.h"

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
  fNoverlaps=0;
  fClusterCharge=0;
  fTrackModel=0;
  fCrossingAngle=0;
  fParSigmaY2=0;
  fParSigmaZ2=0;
}


AliL3ModelTrack::~AliL3ModelTrack()
{
  if(fClusters)
    delete [] fClusters;
  if(fPad)
    delete [] fPad;
  if(fTime)
    delete [] fTime;
  if(fCrossingAngle)
    delete [] fCrossingAngle;
  if(fParSigmaY2)
    delete [] fParSigmaY2;
  if(fParSigmaZ2)
    delete [] fParSigmaZ2;
  if(fTrackModel)
    delete fTrackModel;
  if(fNoverlaps)
    delete [] fNoverlaps;
  if(fOverlap)
    {
      for(Int_t i=0; i<AliL3Transform::GetNRows(fPatch); i++)
	delete [] fOverlap[i];
      delete [] fOverlap;
    }
}

void AliL3ModelTrack::Init(Int_t slice,Int_t patch)
{
  fNClusters=AliL3Transform::GetNRows(patch);//is not incremented in setcluster anymore
  //fSlice=slice;
  fPatch=patch;
  Int_t nrows = AliL3Transform::GetNRows(fPatch);
  fClusters = new AliL3ClusterModel[nrows];
  fPad = new Float_t[nrows];
  fTime = new Float_t[nrows];
  fCrossingAngle = new Float_t[nrows];
  fParSigmaY2 = new Float_t[nrows];
  fParSigmaZ2 = new Float_t[nrows];
  fTrackModel = new AliL3TrackModel;
  
  fOverlap = new Int_t*[nrows];
  fNoverlaps = new Int_t[nrows];
  fMaxOverlaps = 5;
  
  memset(fNoverlaps,0,nrows*sizeof(Int_t));
  memset(fClusters,0,nrows*sizeof(AliL3ClusterModel));
  memset(fPad,0,nrows*sizeof(Float_t));
  memset(fTime,0,nrows*sizeof(Float_t));
  memset(fCrossingAngle,0,nrows*sizeof(Float_t));
  memset(fParSigmaY2,0,nrows*sizeof(Float_t));
  memset(fParSigmaZ2,0,nrows*sizeof(Float_t));
  memset(fTrackModel,0,sizeof(AliL3TrackModel));
  for(Int_t i=0; i<nrows; i++)
    {
      fOverlap[i] = new Int_t[fMaxOverlaps];
      for(Int_t j=0; j<fMaxOverlaps; j++)
	fOverlap[i][j]=-1;
      fClusters[i].fSlice = -1;
    }

}

void AliL3ModelTrack::CalculateClusterWidths(Int_t row,Bool_t parametrize)
{
  //Cluster widths
  
  Float_t xyz[3];
  Int_t sr,lr;
  Int_t index = row - AliL3Transform::GetFirstRow(fPatch);
  if(index < 0 || index > AliL3Transform::GetNRows(fPatch))
    {
      cerr<<"AliL3ModelTrack::CalculcateClusterWidths : Wrond index "<<index<<" row "<<row<<endl;
      return;
    }
  Int_t patch = AliL3Transform::GetPatch(row);
  AliL3Transform::Slice2Sector(0,row,sr,lr);
  AliL3Transform::Raw2Local(xyz,sr,lr,GetPadHit(row),GetTimeHit(row));
  fParSigmaY2[index] = AliL3Transform::GetParSigmaY2(row,xyz[2],GetCrossingAngleLUT(row));
  fParSigmaZ2[index] = AliL3Transform::GetParSigmaZ2(row,xyz[2],GetTgl());
  
  if(parametrize)
    {
      fParSigmaY2[index] = (fParSigmaY2[index] + (1./12)*pow(AliL3Transform::GetPadPitchWidth(patch),2) );
      fParSigmaY2[index] *= 0.108;
      if(patch<2)
	fParSigmaY2[index] *= 2.07;
     
      fParSigmaZ2[index] = (fParSigmaZ2[index] + (1./12)*pow(AliL3Transform::GetZWidth(),2) );
      fParSigmaZ2[index] *= 0.169;
      if(patch<2)
	fParSigmaZ2[index] *= 1.77;
    }
  
  //convert to raw coordinates:
  fParSigmaY2[index] /= pow(AliL3Transform::GetPadPitchWidth(patch),2);
  fParSigmaZ2[index] /= pow(AliL3Transform::GetZWidth(),2);
}

void AliL3ModelTrack::SetCluster(Int_t row,Float_t fpad,Float_t ftime,Float_t charge,
				 Float_t sigmaY2,Float_t sigmaZ2,Int_t npads)
{
  AliL3ClusterModel *cl = GetClusterModel(row);
  
  //First bit: Cluster is present or not
  //Second bit: Cluster was set, meaning an fit attempt was done (if true)
  
  cl->fPresent |= 0x2; //set second bit to true, because a fit attempt has been made
  
  Int_t patch = AliL3Transform::GetPatch(row);
  if(!charge || npads == 1)
    {
      cl->fPresent &= ~0x1; //set first bit to false
    }
  else
    {
      cl->fPresent|=0x1;//set first bit to true
      cl->fDTime = (ftime - GetTimeHit(row))/(AliL3DataCompressor::GetZResidualStep(row)/AliL3Transform::GetZWidth());   
      cl->fDPad = (fpad - GetPadHit(row))/(AliL3DataCompressor::GetXYResidualStep(row)/AliL3Transform::GetPadPitchWidth(patch));
      cl->fDCharge = charge;// - AliL3DataCompressor::GetClusterCharge();
      //cl->fSlice = fSlice;
      if(sigmaY2==0 && sigmaZ2==0)
	{
	  cl->fDSigmaY2=0;//if width is zero, shape is not supposed to be written
	  cl->fDSigmaZ2=0;
	}
      else
	{
	  cl->fDSigmaY2 = (sigmaY2 - GetParSigmaY2(row))/(AliL3DataCompressor::GetXYWidthStep()/pow(AliL3Transform::GetPadPitchWidth(patch),2));
	  cl->fDSigmaZ2 = (sigmaZ2 - GetParSigmaZ2(row))/(AliL3DataCompressor::GetZWidthStep()/pow(AliL3Transform::GetZWidth(),2));
	}
      cl->fNPads = npads;
    }
  
  //fNClusters++;
}

void AliL3ModelTrack::Set(AliL3Track *tpt)
{
  AliL3ModelTrack *tr = (AliL3ModelTrack*)tpt;
  SetRowRange(tr->GetFirstRow(),tr->GetLastRow());
  SetPhi0(tr->GetPhi0());
  SetKappa(tr->GetKappa());
  SetFirstPoint(tr->GetFirstPointX(),tr->GetFirstPointY(),tr->GetFirstPointZ());
  SetLastPoint(tr->GetLastPointX(),tr->GetLastPointY(),tr->GetLastPointZ());
  SetPt(tr->GetPt());
  SetPsi(tr->GetPsi());
  SetTgl(tr->GetTgl());
  SetCharge(tr->GetCharge());
  
  if(fClusters)
    {
      cerr<<"AliL3ModelTrack::Set : Init has already been called for this object!"<<endl;
      return;
    }

  //Init(tr->fSlice,tr->fPatch);
  Init(0,tr->fPatch);
  memcpy(fClusters,tr->fClusters,AliL3Transform::GetNRows(fPatch)*sizeof(AliL3ClusterModel));
  memcpy(fPad,tr->fPad,AliL3Transform::GetNRows(fPatch)*sizeof(Float_t));
  memcpy(fTime,tr->fTime,AliL3Transform::GetNRows(fPatch)*sizeof(Float_t));
  memcpy(fParSigmaY2,tr->fParSigmaY2,AliL3Transform::GetNRows(fPatch)*sizeof(Float_t));
  memcpy(fParSigmaZ2,tr->fParSigmaZ2,AliL3Transform::GetNRows(fPatch)*sizeof(Float_t));
  memcpy(fCrossingAngle,tr->fCrossingAngle,AliL3Transform::GetNRows(fPatch)*sizeof(Float_t));
  memcpy(fTrackModel,tr->fTrackModel,sizeof(AliL3TrackModel));

}

Int_t AliL3ModelTrack::GetNPresentClusters()
{
  //Return the number of assigned clusters to the track.
  //Differs from fNClusters, which should be equal to the 
  //number of padrows in the present patch.
  
  Int_t count=0;

  for(Int_t i=AliL3Transform::GetFirstRow(fPatch); i<=AliL3Transform::GetLastRow(fPatch); i++)
    if(IsPresent(i))
      count++;

  return count;
}

void AliL3ModelTrack::FillModel()
{
  //Fill the track structure
  
  if(fNClusters != AliL3Transform::GetNRows(fPatch))
    {
      cout<<"AliL3ModelTrack::FillModel : fNClusters != nrows; beware, this could be caused by a bug!!!"<<endl;
      fNClusters = AliL3Transform::GetNRows(fPatch);
    }

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
  //fTrackModel->fLength = (Short_t)GetLength();
  //fTrackModel->fNClusters = fNClusters;
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
  //SetLength(fTrackModel->fLength);
  fNClusters = AliL3Transform::GetNRows(fPatch);//fTrackModel->fNClusters;
  SetPt((AliL3Transform::GetBFact()*AliL3Transform::GetBField())/fabs(GetKappa()));
    
  CalculateHelix();

  for(Int_t i=AliL3Transform::GetFirstRow(fPatch); i<=AliL3Transform::GetLastRow(fPatch); i++)
    {
      AliL3ClusterModel *cl = GetClusterModel(i);
      if(!cl) continue;

      if(cl->fSlice == -1)
	{
	  SetPadHit(i,-1);
	  SetTimeHit(i,-1);
	  continue;
	}
      if(cl->fSlice < 0 || cl->fSlice > 35)
	{
	  cerr<<"AliL3ModelTrack::FillTrack : Slice out of range "<<cl->fSlice<<" on row "<<i<<endl;
	  exit(5);
	}
      
      Float_t angle = 0;
      
      AliL3Transform::Local2GlobalAngle(&angle,cl->fSlice);
      if(!CalculateReferencePoint(angle,AliL3Transform::Row2X(i)))
	{
	  if(IsPresent(i))
	    {
	      cerr<<"AliL3ModelTrack::FillTrack : Track does not cross slice "<<cl->fSlice<<" row "<<i<<" Points "
		  <<GetPointX()<<" "<<GetPointY()<<" "<<GetPointZ()<<endl;
	      Print();
	      exit(5);
	    }
	  SetPadHit(i,-1);
	  SetTimeHit(i,-1);
	  continue;
	}
      Float_t hit[3] = {GetPointX(),GetPointY(),GetPointZ()};
      Int_t sector,row;
      AliL3Transform::Slice2Sector(cl->fSlice,i,sector,row);
      AliL3Transform::Global2Raw(hit,sector,row);

      SetPadHit(i,hit[1]);
      SetTimeHit(i,hit[2]);

      Float_t crossingangle = GetCrossingAngle(i,cl->fSlice);
      
      SetCrossingAngleLUT(i,crossingangle);
      CalculateClusterWidths(i,kTRUE);
      
    }
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

void AliL3ModelTrack::SetCrossingAngleLUT(Int_t row,Float_t angle)
{
  Int_t index = row-AliL3Transform::GetFirstRow(fPatch);
  if(index < 0 || index > AliL3Transform::GetNRows(fPatch))
    {
      cerr<<"AliL3ModelTrack::SetCrossingAngle() : Wrong index: "<<index<<endl;
      return;
    }
  fCrossingAngle[index]=angle;
}

void AliL3ModelTrack::SetOverlap(Int_t row,Int_t id)
{

  Int_t index = row-AliL3Transform::GetFirstRow(fPatch);
  if(index < 0 || index > AliL3Transform::GetNRows(fPatch))
    {
      cerr<<"AliL3ModelTrack::SetOverlap() : Wrong index: "<<index<<endl;
      return;
    }
  if(fNoverlaps[index] >= fMaxOverlaps) return;
  fOverlap[index][fNoverlaps[index]++] = id;
}

Bool_t AliL3ModelTrack::IsPresent(Int_t row)
{
  AliL3ClusterModel *cl = GetClusterModel(row);
  return (Bool_t)(cl->fPresent & 0x1);
}

Bool_t AliL3ModelTrack::IsSet(Int_t row)
{
  AliL3ClusterModel *cl = GetClusterModel(row);
  return (Bool_t)(cl->fPresent & 0x2);
}

Int_t AliL3ModelTrack::GetNPads(Int_t row)
{
  AliL3ClusterModel *cl = GetClusterModel(row);
  return cl->fNPads;
}

Bool_t AliL3ModelTrack::GetPad(Int_t row,Float_t &pad)
{
  //(ftime - GetTimeHit(fNClusters))/AliL3DataCompressor::GetXYResidualStep();
  //(fpad - GetPadHit(fNClusters))/AliL3DataCompressor::GetZResidualStep();

  AliL3ClusterModel *cl = GetClusterModel(row);
  Int_t patch = AliL3Transform::GetPatch(row);
  pad = cl->fDPad*(AliL3DataCompressor::GetXYResidualStep(row)/AliL3Transform::GetPadPitchWidth(patch)) + GetPadHit(row);
  
  return IsPresent(row);
}

Bool_t AliL3ModelTrack::GetTime(Int_t row,Float_t &time)
{
  AliL3ClusterModel *cl = GetClusterModel(row);
  time = cl->fDTime*(AliL3DataCompressor::GetZResidualStep(row)/AliL3Transform::GetZWidth()) + GetTimeHit(row);
  
  return IsPresent(row);
}

Bool_t AliL3ModelTrack::GetClusterCharge(Int_t row,Int_t &charge)
{
  AliL3ClusterModel *cl = GetClusterModel(row);
  charge = (Int_t)cl->fDCharge;// + AliL3DataCompressor::GetClusterCharge();
  
  return IsPresent(row);
}

Bool_t AliL3ModelTrack::GetXYWidth(Int_t row,Float_t &width)
{
  //cl->fDSigmaY2 = (sigmaY2 - GetParSigmaY2(row))/AliL3DataCompressor::GetXYWidthStep();
  AliL3ClusterModel *cl = GetClusterModel(row);
  Int_t patch = AliL3Transform::GetPatch(row);
  width = cl->fDSigmaY2*(AliL3DataCompressor::GetXYWidthStep()/pow(AliL3Transform::GetPadPitchWidth(patch),2)) + GetParSigmaY2(row);

  return IsPresent(row);
}

Bool_t AliL3ModelTrack::GetZWidth(Int_t row,Float_t &width)
{
  AliL3ClusterModel *cl = GetClusterModel(row);
  width = cl->fDSigmaZ2*(AliL3DataCompressor::GetZWidthStep()/pow(AliL3Transform::GetZWidth(),2)) + GetParSigmaZ2(row);

  return IsPresent(row);
}

Bool_t AliL3ModelTrack::GetPadResidual(Int_t row,Float_t &res)
{
  AliL3ClusterModel *cl = GetClusterModel(row);
  res = cl->fDPad;
  return IsPresent(row);
}

Bool_t AliL3ModelTrack::GetTimeResidual(Int_t row,Float_t &res)
{
  AliL3ClusterModel *cl = GetClusterModel(row);
  res = cl->fDTime;
  return IsPresent(row);
}

Bool_t AliL3ModelTrack::GetXYWidthResidual(Int_t row,Float_t &res)
{
  AliL3ClusterModel *cl = GetClusterModel(row);
  res = cl->fDSigmaY2;
  return IsPresent(row);
}

Bool_t AliL3ModelTrack::GetZWidthResidual(Int_t row,Float_t &res)
{
  AliL3ClusterModel *cl = GetClusterModel(row);
  res = cl->fDSigmaZ2;
  return IsPresent(row);
}

Int_t AliL3ModelTrack::GetSlice(Int_t row)
{
  AliL3ClusterModel *cl = GetClusterModel(row);
  return cl->fSlice;
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

Float_t AliL3ModelTrack::GetCrossingAngleLUT(Int_t row)
{
  Int_t index = row-AliL3Transform::GetFirstRow(fPatch);
  if(index < 0 || index > AliL3Transform::GetNRows(fPatch))
    {
      cerr<<"AliL3ModelTrack::GetCrossingAngleLUT() : Wrong index: "<<index<<" row "<<row<<endl;
      return 0;
    }
  return fCrossingAngle[index];
}

Float_t AliL3ModelTrack::GetParSigmaY2(Int_t row)
{
  Int_t index = row-AliL3Transform::GetFirstRow(fPatch);
  if(index < 0 || index > AliL3Transform::GetNRows(fPatch))
    {
      cerr<<"AliL3ModelTrack::GetParSigmaY2() : Wrong index: "<<index<<" row "<<row<<endl;
      return 0;
    }
  return fParSigmaY2[index];
}

Float_t AliL3ModelTrack::GetParSigmaZ2(Int_t row)
{
  Int_t index = row-AliL3Transform::GetFirstRow(fPatch);
  if(index < 0 || index > AliL3Transform::GetNRows(fPatch))
    {
      cerr<<"AliL3ModelTrack::GetParSigmaZ2() : Wrong index: "<<index<<" row "<<row<<endl;
      return 0;
    }
  return fParSigmaZ2[index];
}

Int_t AliL3ModelTrack::GetNOverlaps(Int_t row)
{
  Int_t index = row - AliL3Transform::GetFirstRow(fPatch);
  if(index < 0 || index > AliL3Transform::GetNRows(fPatch))
    {
      cerr<<"AliL3ModelTrack::GetOverlap() : Wrong index: "<<index<<endl;
      return 0;
    }
  return fNoverlaps[index];
}

Int_t *AliL3ModelTrack::GetOverlaps(Int_t row)
{
  Int_t index = row - AliL3Transform::GetFirstRow(fPatch);
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

void AliL3ModelTrack::Print(Bool_t everything)
{
  //Print info

  cout<<"First point "<<GetFirstPointX()<<" "<<GetFirstPointY()<<" "<<GetFirstPointZ()<<endl;
  cout<<"Last point "<<GetLastPointX()<<" "<<GetLastPointY()<<" "<<GetLastPointZ()<<endl;
  cout<<"Pt "<<GetPt()<<" kappa "<<GetKappa()<<" tgl "<<GetTgl()<<" psi "<<GetPsi()<<" charge "<<GetCharge()<<endl;
  cout<<"Center "<<GetCenterX()<<" "<<GetCenterY()<<endl<<endl;
  if(!everything)
    return;
  cout<<"NHits "<<GetNClusters()<<endl;

  cout<<"Clusters:"<<endl;
  Int_t origslice=-1,counter=0;
  Float_t fpad,ftime;
  for(Int_t i=AliL3Transform::GetFirstRow(fPatch); i<=AliL3Transform::GetLastRow(fPatch); i++)
    {
      AliL3ClusterModel *cl = GetClusterModel(i);
      
      if(!IsPresent(i))
	cout<<i<<" Empty"<<" Slice "<<cl->fSlice<<" Padcrossing "<<GetPadHit(i)<<" Timecrossing "<<GetTimeHit(i)<<" ";
      else
	{
	  GetPad(i,fpad);
	  GetTime(i,ftime);
	  if(counter==0)
	    origslice=cl->fSlice;
	  else if(cl->fSlice != origslice)
	    cout<<"Change in slice "<<cl->fSlice<<" "<<origslice<<endl;
	  cout<<i<<" Slice "<<cl->fSlice<<" Dpad "<<cl->fDPad<<" Dtime "<<cl->fDTime<<" Dcharge "<<cl->fDCharge;
	  cout<<" sigmaY2 "<<GetParSigmaY2(i)<<" sigmaZ2 "<<GetParSigmaZ2(i);
	  cout<<" Pad "<<fpad<<" padhit "<<GetPadHit(i)<<" Time "<<ftime<<" timehit "<<GetTimeHit(i)<<" ";
	  cout<<"Number of pads "<<GetNPads(i)<<" Overlaps "<<GetNOverlaps(i);
	  counter++;
	}
      cout<<endl;
    }
}

void AliL3ModelTrack::SetClusterLabel(Int_t row,Int_t *trackID)
{
#ifdef do_mc
  AliL3ClusterModel *cl = GetClusterModel(row);
  cl->fTrackID[0] = trackID[0];
  cl->fTrackID[1] = trackID[1];
  cl->fTrackID[2] = trackID[2];
#endif
  return;
}

void AliL3ModelTrack::GetClusterLabel(Int_t row,Int_t *trackID)
{
#ifdef do_mc
  AliL3ClusterModel *cl = GetClusterModel(row);
  trackID[0] = cl->fTrackID[0];
  trackID[1] = cl->fTrackID[1];
  trackID[2] = cl->fTrackID[2];
#endif
  return;
}

