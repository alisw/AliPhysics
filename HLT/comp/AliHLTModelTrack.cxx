// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo$fi.uib.no>
//*-- Copyright &copy ALICE HLT Group
//_____________________________________________________________
// AliHLTModelTrack
//
// 


#include "AliHLTStandardIncludes.h"

#include "AliHLTLogging.h"
#include "AliHLTTransform.h"
#include "AliHLTVertex.h"
#include "AliHLTDataCompressorHelper.h"

#include "AliHLTModelTrack.h"

#if __GNUC__ >= 3
using namespace std;
#endif

ClassImp(AliHLTModelTrack)

AliHLTModelTrack::AliHLTModelTrack() 
{
  // default constructor
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
  fArraysCreated=kFALSE;
}


AliHLTModelTrack::~AliHLTModelTrack()
{
  // destructor
  DeleteArrays();
}

void AliHLTModelTrack::DeleteArrays()
{
  // deletes all arrays
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
      for(Int_t i=0; i<AliHLTTransform::GetNRows(fPatch); i++)
	delete [] fOverlap[i];
      delete [] fOverlap;
    }
  fArraysCreated=kFALSE;
}

void AliHLTModelTrack::Init(Int_t /*slice*/,Int_t patch)
{
  // Initialization
  if(fArraysCreated)
    {               
      DeleteArrays();
    }
  fNClusters=AliHLTTransform::GetNRows(patch);
  fPatch=patch;
  Int_t nrows = AliHLTTransform::GetNRows(fPatch);
  fClusters = new AliHLTClusterModel[nrows];
  fPad = new Float_t[nrows];
  fTime = new Float_t[nrows];
  fCrossingAngle = new Float_t[nrows];
  fParSigmaY2 = new Float_t[nrows];
  fParSigmaZ2 = new Float_t[nrows];
  fTrackModel = new AliHLTTrackModel;
  
  fOverlap = new Int_t*[nrows];
  fNoverlaps = new Int_t[nrows];
  fMaxOverlaps = 5;
  
  memset(fNoverlaps,0,nrows*sizeof(Int_t));
  memset(fClusters,0,nrows*sizeof(AliHLTClusterModel));
  memset(fPad,0,nrows*sizeof(Float_t));
  memset(fTime,0,nrows*sizeof(Float_t));
  memset(fCrossingAngle,0,nrows*sizeof(Float_t));
  memset(fParSigmaY2,0,nrows*sizeof(Float_t));
  memset(fParSigmaZ2,0,nrows*sizeof(Float_t));
  memset(fTrackModel,0,sizeof(AliHLTTrackModel));
  for(Int_t i=0; i<nrows; i++)
    {
      fOverlap[i] = new Int_t[fMaxOverlaps];
      for(Int_t j=0; j<fMaxOverlaps; j++)
	fOverlap[i][j]=-1;
      fClusters[i].fSlice = -1;
    }
  fArraysCreated=kTRUE;
}



void AliHLTModelTrack::CalculateClusterWidths(Int_t row,Bool_t parametrize)
{
  //Cluster widths
  
  Float_t xyz[3];
  Int_t sr,lr;
  Int_t index = row - AliHLTTransform::GetFirstRow(fPatch);
  if(index < 0 || index > AliHLTTransform::GetNRows(fPatch))
    {
      cerr<<"AliHLTModelTrack::CalculcateClusterWidths : Wrond index "<<index<<" row "<<row<<endl;
      return;
    }
  Int_t patch = AliHLTTransform::GetPatch(row);
  AliHLTTransform::Slice2Sector(0,row,sr,lr);
  AliHLTTransform::Raw2Local(xyz,sr,lr,GetPadHit(row),GetTimeHit(row));
  fParSigmaY2[index] = AliHLTTransform::GetParSigmaY2(row,xyz[2],GetCrossingAngleLUT(row));
  fParSigmaZ2[index] = AliHLTTransform::GetParSigmaZ2(row,xyz[2],GetTgl());
  
  if(parametrize)
    {
      fParSigmaY2[index] = (fParSigmaY2[index] + (1./12)*pow(AliHLTTransform::GetPadPitchWidth(patch),2) );
      fParSigmaY2[index] *= 0.108;
      if(patch<2)
	fParSigmaY2[index] *= 2.07;
     
      fParSigmaZ2[index] = (fParSigmaZ2[index] + (1./12)*pow(AliHLTTransform::GetZWidth(),2) );
      fParSigmaZ2[index] *= 0.169;
      if(patch<2)
	fParSigmaZ2[index] *= 1.77;
    }
  
  //convert to raw coordinates:
  fParSigmaY2[index] /= pow(AliHLTTransform::GetPadPitchWidth(patch),2);
  fParSigmaZ2[index] /= pow(AliHLTTransform::GetZWidth(),2);
}

void AliHLTModelTrack::SetCluster(Int_t row,Float_t fpad,Float_t ftime,Float_t charge,
				 Float_t sigmaY2,Float_t sigmaZ2,Int_t npads)
{
  AliHLTClusterModel *cl = GetClusterModel(row);
  
  //First bit: Cluster is present or not
  //Second bit: Cluster was set, meaning an fit attempt was done (if true)
  
  cl->fPresent |= 0x2; //set second bit to true, because a fit attempt has been made
  
  Int_t patch = AliHLTTransform::GetPatch(row);
  if(!charge || npads == 1)
    {
      cl->fPresent &= ~0x1; //set first bit to false
    }
  else
    {
      cl->fPresent|=0x1;//set first bit to true
      cl->fDPad = (fpad - GetPadHit(row))/(AliHLTDataCompressorHelper::GetXYResidualStep(row)/AliHLTTransform::GetPadPitchWidth(patch));
      cl->fDTime = (ftime - GetTimeHit(row))/(AliHLTDataCompressorHelper::GetZResidualStep(row)/AliHLTTransform::GetZWidth());
      cl->fDCharge = charge;
      if(sigmaY2==0 && sigmaZ2==0)
	{
	  cl->fDSigmaY=0;//if width is zero, shape is not supposed to be written
	  cl->fDSigmaZ=0;
	}
      else
	{
	  //cl->fDSigmaY2 = (sigmaY2 - GetParSigmaY2(row))/(pow(AliHLTDataCompressorHelper::GetXYWidthStep(),2)/pow(AliHLTTransform::GetPadPitchWidth(patch),2));
	  //cl->fDSigmaZ2 = (sigmaZ2 - GetParSigmaZ2(row))/(pow(AliHLTDataCompressorHelper::GetZWidthStep(),2)/pow(AliHLTTransform::GetZWidth(),2));
	  cl->fDSigmaY = (sqrt(sigmaY2) - sqrt(GetParSigmaY2(row)))/(AliHLTDataCompressorHelper::GetXYWidthStep()/AliHLTTransform::GetPadPitchWidth(patch));
	  cl->fDSigmaZ = (sqrt(sigmaZ2) - sqrt(GetParSigmaZ2(row)))/(AliHLTDataCompressorHelper::GetZWidthStep()/AliHLTTransform::GetZWidth());
	}
      cl->fNPads = npads;
    }
}


void AliHLTModelTrack::Set(AliHLTTrack *tpt)
{
  // Sets track and does initialization
  AliHLTModelTrack *tr = (AliHLTModelTrack*)tpt;
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
      cerr<<"AliHLTModelTrack::Set : Init has already been called for this object!"<<endl;
      return;
    }

  //Init(tr->fSlice,tr->fPatch);
  Init(0,tr->fPatch);
  memcpy(fClusters,tr->fClusters,AliHLTTransform::GetNRows(fPatch)*sizeof(AliHLTClusterModel));
  memcpy(fPad,tr->fPad,AliHLTTransform::GetNRows(fPatch)*sizeof(Float_t));
  memcpy(fTime,tr->fTime,AliHLTTransform::GetNRows(fPatch)*sizeof(Float_t));
  memcpy(fParSigmaY2,tr->fParSigmaY2,AliHLTTransform::GetNRows(fPatch)*sizeof(Float_t));
  memcpy(fParSigmaZ2,tr->fParSigmaZ2,AliHLTTransform::GetNRows(fPatch)*sizeof(Float_t));
  memcpy(fCrossingAngle,tr->fCrossingAngle,AliHLTTransform::GetNRows(fPatch)*sizeof(Float_t));
  memcpy(fTrackModel,tr->fTrackModel,sizeof(AliHLTTrackModel));

}

Int_t AliHLTModelTrack::GetNPresentClusters()
{
  //Return the number of assigned clusters to the track.
  //Differs from fNClusters, which should be equal to the 
  //number of padrows in the present patch.
  
  Int_t count=0;

  for(Int_t i=AliHLTTransform::GetFirstRow(fPatch); i<=AliHLTTransform::GetLastRow(fPatch); i++)
    if(IsPresent(i))
      count++;

  return count;
}

void AliHLTModelTrack::FillModel()
{
  //Fill the track structure
  
  if(fNClusters != AliHLTTransform::GetNRows(fPatch))
    {
      cout<<"AliHLTModelTrack::FillModel : fNClusters != nrows; beware, this could be caused by a bug!!!"<<endl;
      fNClusters = AliHLTTransform::GetNRows(fPatch);
    }

  if(!fTrackModel)
    {
      cerr<<"AliHLTModelTrack::FillModel() : No trackmodel "<<endl;
      return;
    }
  Double_t impact[3];
  AliHLTVertex vertex;
  CalculateHelix();
  GetClosestPoint(&vertex,impact[0],impact[1],impact[2]);
  fTrackModel->fKappa = GetKappa();
  fTrackModel->fPhi = atan2(impact[1],impact[0]);
  fTrackModel->fD = sqrt(impact[0]*impact[0] + impact[1]*impact[1]);
  fTrackModel->fZ0 = impact[2];
  fTrackModel->fTgl = GetTgl();
  
  //We have to check on which of the vertex the track fit is lying
  //This we need to encode the azimuthal angle coordinate of the center of curvature.
  if(GetRadius() < sqrt(GetCenterX()*GetCenterX()+GetCenterY()*GetCenterY()))
    fTrackModel->fD *=-1;
  
}

void AliHLTModelTrack::FillTrack()
{
  //Fill the track parameters from the structure.
  
  if(!fTrackModel)
    {
      cerr<<"AliHLTModelTrack::FillTrack() : No data!!"<<endl;
      return;
    }
  SetKappa(fTrackModel->fKappa);
  Double_t impact[3],psi;
  Float_t trackPhi0 = fTrackModel->fPhi;
  if(fTrackModel->fD < 0)
    trackPhi0 += AliHLTTransform::Pi();
  Int_t charge = -1*(Int_t)copysign(1.,GetKappa());
  impact[0] = fabs(fTrackModel->fD)*cos(fTrackModel->fPhi);
  impact[1] = fabs(fTrackModel->fD)*sin(fTrackModel->fPhi);
  impact[2] = fTrackModel->fZ0;

  psi = trackPhi0 - charge*0.5*AliHLTTransform::Pi();
  if(psi < 0) 
    psi += 2*AliHLTTransform::Pi();

  SetCharge(charge);
  SetFirstPoint(impact[0],impact[1],impact[2]);
  SetPsi(psi);
  SetTgl(fTrackModel->fTgl);
  SetPt((AliHLTTransform::GetBFact()*AliHLTTransform::GetBField())/fabs(GetKappa()));
  fNClusters = AliHLTTransform::GetNRows(fPatch);
  CalculateHelix();
  
  for(Int_t i=AliHLTTransform::GetFirstRow(fPatch); i<=AliHLTTransform::GetLastRow(fPatch); i++)
    {
      AliHLTClusterModel *cl = GetClusterModel(i);
      if(!cl) continue;

      if(cl->fSlice == -1)
	{
	  SetPadHit(i,-1);
	  SetTimeHit(i,-1);
	  continue;
	}
      if(cl->fSlice < 0 || cl->fSlice > 35)
	{
	  cerr<<"AliHLTModelTrack::FillTrack : Slice out of range "<<cl->fSlice<<" on row "<<i<<endl;
	  exit(5);
	}
      
      Float_t angle = 0;
      
      AliHLTTransform::Local2GlobalAngle(&angle,cl->fSlice);
      if(!CalculateReferencePoint(angle,AliHLTTransform::Row2X(i)))
	{
	  if(IsPresent(i))
	    {
	      cerr<<"AliHLTModelTrack::FillTrack : Track does not cross slice "<<cl->fSlice<<" row "<<i<<" Points "
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
      AliHLTTransform::Slice2Sector(cl->fSlice,i,sector,row);
      AliHLTTransform::Global2Raw(hit,sector,row);

      SetPadHit(i,hit[1]);
      SetTimeHit(i,hit[2]);

      Float_t crossingangle = GetCrossingAngle(i,cl->fSlice);
      
      SetCrossingAngleLUT(i,crossingangle);
      CalculateClusterWidths(i,kTRUE);
      
    }
}

void AliHLTModelTrack::SetPadHit(Int_t row,Float_t pad)
{
  // sets pad hit
  Int_t index = row-AliHLTTransform::GetFirstRow(fPatch);
  if(index < 0 || index > AliHLTTransform::GetNRows(fPatch))
    {
      cerr<<"AliHLTModelTrack::SetPadHit() : Wrong index: "<<index<<endl;
      return;
    }
  fPad[index]=pad;
}

void AliHLTModelTrack::SetTimeHit(Int_t row,Float_t time)
{
  // sets time hit
  Int_t index = row-AliHLTTransform::GetFirstRow(fPatch);
  if(index < 0 || index > AliHLTTransform::GetNRows(fPatch))
    {
      cerr<<"AliHLTModelTrack::SetTimeHit() : Wrong index: "<<index<<endl;
      return;
    }
  fTime[index]=time;
}

void AliHLTModelTrack::SetCrossingAngleLUT(Int_t row,Float_t angle)
{
  // sets LUT for crossing angle
  Int_t index = row-AliHLTTransform::GetFirstRow(fPatch);
  if(index < 0 || index > AliHLTTransform::GetNRows(fPatch))
    {
      cerr<<"AliHLTModelTrack::SetCrossingAngle() : Wrong index: "<<index<<endl;
      return;
    }
  fCrossingAngle[index]=angle;
}

void AliHLTModelTrack::SetOverlap(Int_t row,Int_t id)
{
  // sets overlap
  Int_t index = row-AliHLTTransform::GetFirstRow(fPatch);
  if(index < 0 || index > AliHLTTransform::GetNRows(fPatch))
    {
      cerr<<"AliHLTModelTrack::SetOverlap() : Wrong index: "<<index<<endl;
      return;
    }
  if(fNoverlaps[index] >= fMaxOverlaps) return;
  fOverlap[index][fNoverlaps[index]++] = id;
}

Bool_t AliHLTModelTrack::IsPresent(Int_t row)
{
  AliHLTClusterModel *cl = GetClusterModel(row);
  return (Bool_t)(cl->fPresent & 0x1);
}

Bool_t AliHLTModelTrack::IsSet(Int_t row)
{
  // checks if row was set
  AliHLTClusterModel *cl = GetClusterModel(row);
  return (Bool_t)(cl->fPresent & 0x2);
}

Int_t AliHLTModelTrack::GetNPads(Int_t row)
{
  // gets number of pads
  AliHLTClusterModel *cl = GetClusterModel(row);
  return cl->fNPads;
}

Bool_t AliHLTModelTrack::GetPad(Int_t row,Float_t &pad)
{
  // gets pad
  //(fpad - GetPadHit(row))/(AliHLTDataCompressorHelper::GetXYResidualStep(row)/AliHLTTransform::GetPadPitchWidth(patch));
  AliHLTClusterModel *cl = GetClusterModel(row);
  Int_t patch = AliHLTTransform::GetPatch(row);
  pad = cl->fDPad*(AliHLTDataCompressorHelper::GetXYResidualStep(row)/AliHLTTransform::GetPadPitchWidth(patch)) + GetPadHit(row);
  return IsPresent(row);
}

Bool_t AliHLTModelTrack::GetTime(Int_t row,Float_t &time)
{
  // gets time
  AliHLTClusterModel *cl = GetClusterModel(row);
  time = cl->fDTime*(AliHLTDataCompressorHelper::GetZResidualStep(row)/AliHLTTransform::GetZWidth()) + GetTimeHit(row);
  return IsPresent(row);
}

Bool_t AliHLTModelTrack::GetClusterCharge(Int_t row,Int_t &charge)
{
  // gets cluster's charge
  AliHLTClusterModel *cl = GetClusterModel(row);
  charge = (Int_t)cl->fDCharge;// + AliHLTDataCompressorHelperHelper::GetClusterCharge();
  return IsPresent(row);
}

Bool_t AliHLTModelTrack::GetSigmaY2(Int_t row,Float_t &sigma2)
{
  // gets SigmaY2
  //cl->fDSigmaY = (sqrt(sigmaY2) - sqrt(GetParSigmaY2(row)))/(AliHLTDataCompressorHelper::GetXYWidthStep()/AliHLTTransform::GetPadPitchWidth(patch));
  AliHLTClusterModel *cl = GetClusterModel(row);
  Int_t patch = AliHLTTransform::GetPatch(row);
  Float_t sigma = cl->fDSigmaY*(AliHLTDataCompressorHelper::GetXYWidthStep()/AliHLTTransform::GetPadPitchWidth(patch)) + sqrt(GetParSigmaY2(row));
  sigma2 = sigma*sigma;
  return IsPresent(row);
}

Bool_t AliHLTModelTrack::GetSigmaZ2(Int_t row,Float_t &sigma2)
{
  // gets SigmaZ2
  //cl->fDSigmaZ = (sqrt(sigmaZ2) - sqrt(GetParSigmaZ2(row)))/(AliHLTDataCompressorHelper::GetZWidthStep()/AliHLTTransform::GetZWidth());
  AliHLTClusterModel *cl = GetClusterModel(row);
  Float_t sigma = cl->fDSigmaZ*(AliHLTDataCompressorHelper::GetZWidthStep()/AliHLTTransform::GetZWidth()) + sqrt(GetParSigmaZ2(row));
  sigma2 = sigma*sigma;
  return IsPresent(row);
}

Bool_t AliHLTModelTrack::GetPadResidual(Int_t row,Float_t &res)
{
  // gets pad residual
  AliHLTClusterModel *cl = GetClusterModel(row);
  Int_t patch = AliHLTTransform::GetPatch(row);
  res = cl->fDPad*(AliHLTDataCompressorHelper::GetXYResidualStep(row)/AliHLTTransform::GetPadPitchWidth(patch));
  return IsPresent(row);
}

Bool_t AliHLTModelTrack::GetTimeResidual(Int_t row,Float_t &res)
{
  // gets time residual
  AliHLTClusterModel *cl = GetClusterModel(row);
  res = cl->fDTime*(AliHLTDataCompressorHelper::GetZResidualStep(row)/AliHLTTransform::GetZWidth());
  return IsPresent(row);
}

Bool_t AliHLTModelTrack::GetSigmaYResidual(Int_t row,Float_t &res)
{
  // gets SigmaY residual (?)
  AliHLTClusterModel *cl = GetClusterModel(row);
  Int_t patch = AliHLTTransform::GetPatch(row);
  res = cl->fDSigmaY*(AliHLTDataCompressorHelper::GetXYWidthStep()/AliHLTTransform::GetPadPitchWidth(patch));
  return IsPresent(row);
}

Bool_t AliHLTModelTrack::GetSigmaZResidual(Int_t row,Float_t &res)
{
  // gets SigmaZ resigual (?)
  AliHLTClusterModel *cl = GetClusterModel(row);
  res = cl->fDSigmaZ*(AliHLTDataCompressorHelper::GetZWidthStep()/AliHLTTransform::GetZWidth());
  return IsPresent(row);
}

Int_t AliHLTModelTrack::GetSlice(Int_t row)
{
  // Gets slice
  AliHLTClusterModel *cl = GetClusterModel(row);
  return cl->fSlice;
}

Float_t AliHLTModelTrack::GetPadHit(Int_t row)
{
  // Gets pad hit
  Int_t index = row-AliHLTTransform::GetFirstRow(fPatch);
  if(index < 0 || index > AliHLTTransform::GetNRows(fPatch))
    {
      cerr<<"AliHLTModelTrack::GetPadHit() : Wrong index: "<<index<<" row "<<row<<endl;
      return 0;
    }
  return fPad[index];
}

Float_t AliHLTModelTrack::GetTimeHit(Int_t row)
{
  // Gets time hit
  Int_t index = row-AliHLTTransform::GetFirstRow(fPatch);
  if(index < 0 || index > AliHLTTransform::GetNRows(fPatch))
    {
      cerr<<"AliHLTModelTrack::GetTimeHit() : Wrong index: "<<index<<" row "<<row<<endl;
      return 0;
    }
  return fTime[index];
}

Float_t AliHLTModelTrack::GetCrossingAngleLUT(Int_t row)
{
  // gets LUT for crossing angle
  Int_t index = row-AliHLTTransform::GetFirstRow(fPatch);
  if(index < 0 || index > AliHLTTransform::GetNRows(fPatch))
    {
      cerr<<"AliHLTModelTrack::GetCrossingAngleLUT() : Wrong index: "<<index<<" row "<<row<<endl;
      return 0;
    }
  return fCrossingAngle[index];
}

Float_t AliHLTModelTrack::GetParSigmaY2(Int_t row)
{
  // gets par SigmaY2 (?)
  Int_t index = row-AliHLTTransform::GetFirstRow(fPatch);
  if(index < 0 || index > AliHLTTransform::GetNRows(fPatch))
    {
      cerr<<"AliHLTModelTrack::GetParSigmaY2() : Wrong index: "<<index<<" row "<<row<<endl;
      return 0;
    }
  return fParSigmaY2[index];
}

Float_t AliHLTModelTrack::GetParSigmaZ2(Int_t row)
{
  // gets par SigmaZ2 (?)
  Int_t index = row-AliHLTTransform::GetFirstRow(fPatch);
  if(index < 0 || index > AliHLTTransform::GetNRows(fPatch))
    {
      cerr<<"AliHLTModelTrack::GetParSigmaZ2() : Wrong index: "<<index<<" row "<<row<<endl;
      return 0;
    }
  return fParSigmaZ2[index];
}

Int_t AliHLTModelTrack::GetNOverlaps(Int_t row)
{
  // gets number of overlaps
  Int_t index = row - AliHLTTransform::GetFirstRow(fPatch);
  if(index < 0 || index > AliHLTTransform::GetNRows(fPatch))
    {
      cerr<<"AliHLTModelTrack::GetOverlap() : Wrong index: "<<index<<endl;
      return 0;
    }
  return fNoverlaps[index];
}

Int_t *AliHLTModelTrack::GetOverlaps(Int_t row)
{
  // gets overlaps
  Int_t index = row - AliHLTTransform::GetFirstRow(fPatch);
  if(index < 0 || index > AliHLTTransform::GetNRows(fPatch))
    {
      cerr<<"AliHLTModelTrack::GetOverlap() : Wrong index: "<<index<<endl;
      return 0;
    }
  return fOverlap[index];
}

AliHLTClusterModel *AliHLTModelTrack::GetClusterModel(Int_t row)
{
  // gets cluster model
  if(!fClusters) return 0; 
  Int_t index = row-AliHLTTransform::GetFirstRow(fPatch);
  if(index < 0 || index > AliHLTTransform::GetNRows(fPatch))
    {
      cerr<<"AliHLTModelTrack::GetClusterModel() : Wrong index: "<<index<<endl;
      return 0;
    }
  return &fClusters[index];
}

void AliHLTModelTrack::Print(Bool_t everything)
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
  Float_t fpad,ftime,sigmaY2,sigmaZ2;
  for(Int_t i=AliHLTTransform::GetFirstRow(fPatch); i<=AliHLTTransform::GetLastRow(fPatch); i++)
    {
      AliHLTClusterModel *cl = GetClusterModel(i);
      
      if(!IsPresent(i))
	{
	  cout<<i<<" Empty"<<" Slice "<<cl->fSlice<<" Padcrossing "<<GetPadHit(i)<<" Timecrossing "<<GetTimeHit(i)<<" ";
	  //AliHLTTransform::RawHLT2Global(xyz,cl->fSlice,i,GetPadHit(i),GetTimeHit(i));
	  //cout<<i<<" slice "<<cl->fSlice<<" x "<<xyz[0]<<" y "<<xyz[1]<<" z "<<xyz[2];
	}
      else
	{
	  GetPad(i,fpad);
	  GetTime(i,ftime);
	  GetSigmaY2(i,sigmaY2);
	  GetSigmaZ2(i,sigmaZ2);
	  if(counter==0)
	    origslice=cl->fSlice;
	  else if(cl->fSlice != origslice)
	    cout<<"Change in slice "<<cl->fSlice<<" "<<origslice<<endl;
	  cout<<i<<" Slice "<<cl->fSlice<<" Dpad "<<cl->fDPad<<" Dtime "<<cl->fDTime<<" Dcharge "<<cl->fDCharge;
	  cout<<" sigmaY2 "<<sigmaY2<<" sigmaZ2 "<<sigmaZ2;
	  cout<<" parsigmaY2 "<<GetParSigmaY2(i)<<" parsigmaZ2 "<<GetParSigmaZ2(i);
	  cout<<" Pad "<<fpad<<" padhit "<<GetPadHit(i)<<" Time "<<ftime<<" timehit "<<GetTimeHit(i)<<" ";
	  counter++;
	}
      cout<<endl;
    }
}

#ifdef do_mc
void AliHLTModelTrack::SetClusterLabel(Int_t row,Int_t *trackID)
{
  // sets cluster label
  AliHLTClusterModel *cl = GetClusterModel(row);
  cl->fTrackID[0] = trackID[0];
  cl->fTrackID[1] = trackID[1];
  cl->fTrackID[2] = trackID[2];
#else
  void AliHLTModelTrack::SetClusterLabel(Int_t /*row*/,Int_t */*trackID*/)
{
  // Does nothing if do_mc undefined
  return;
#endif
}

#ifdef do_mc
void AliHLTModelTrack::GetClusterLabel(Int_t row,Int_t *trackID)
{
  // gets cluster label
  AliHLTClusterModel *cl = GetClusterModel(row);
  trackID[0] = cl->fTrackID[0];
  trackID[1] = cl->fTrackID[1];
  trackID[2] = cl->fTrackID[2];
#else
  void AliHLTModelTrack::GetClusterLabel(Int_t /*row*/,Int_t */*trackID*/)
{
  // Does nothing if do_mc undefined
  return;
#endif
}

