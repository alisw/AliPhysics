// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo$fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliHLTTPCStandardIncludes.h"

#include "AliHLTTPCLogging.h"
#include "AliHLTTPCModelTrack.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCDataCompressor.h"

#if GCCVERSION == 3
using namespace std;
#endif

//_____________________________________________________________
// AliHLTTPCModelTrack
//
// 

ClassImp(AliHLTTPCModelTrack)

AliHLTTPCModelTrack::AliHLTTPCModelTrack() 
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


AliHLTTPCModelTrack::~AliHLTTPCModelTrack()
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
      for(Int_t i=0; i<AliHLTTPCTransform::GetNRows(fPatch); i++)
	delete [] fOverlap[i];
      delete [] fOverlap;
    }
}

void AliHLTTPCModelTrack::Init(Int_t slice,Int_t patch)
{
  fNClusters=AliHLTTPCTransform::GetNRows(patch);//is not incremented in setcluster anymore
  //fSlice=slice;
  fPatch=patch;
  Int_t nrows = AliHLTTPCTransform::GetNRows(fPatch);
  fClusters = new AliHLTTPCClusterModel[nrows];
  fPad = new Float_t[nrows];
  fTime = new Float_t[nrows];
  fCrossingAngle = new Float_t[nrows];
  fParSigmaY2 = new Float_t[nrows];
  fParSigmaZ2 = new Float_t[nrows];
  fTrackModel = new AliHLTTPCTrackModel;
  
  fOverlap = new Int_t*[nrows];
  fNoverlaps = new Int_t[nrows];
  fMaxOverlaps = 5;
  
  memset(fNoverlaps,0,nrows*sizeof(Int_t));
  memset(fClusters,0,nrows*sizeof(AliHLTTPCClusterModel));
  memset(fPad,0,nrows*sizeof(Float_t));
  memset(fTime,0,nrows*sizeof(Float_t));
  memset(fCrossingAngle,0,nrows*sizeof(Float_t));
  memset(fParSigmaY2,0,nrows*sizeof(Float_t));
  memset(fParSigmaZ2,0,nrows*sizeof(Float_t));
  memset(fTrackModel,0,sizeof(AliHLTTPCTrackModel));
  for(Int_t i=0; i<nrows; i++)
    {
      fOverlap[i] = new Int_t[fMaxOverlaps];
      for(Int_t j=0; j<fMaxOverlaps; j++)
	fOverlap[i][j]=-1;
      fClusters[i].fSlice = -1;
    }

}

void AliHLTTPCModelTrack::CalculateClusterWidths(Int_t row,Bool_t parametrize)
{
  //Cluster widths
  
  Float_t xyz[3];
  Int_t sr,lr;
  Int_t index = row - AliHLTTPCTransform::GetFirstRow(fPatch);
  if(index < 0 || index > AliHLTTPCTransform::GetNRows(fPatch))
    {
      std::cerr<<"AliHLTTPCModelTrack::CalculcateClusterWidths : Wrond index "<<index<<" row "<<row<<std::endl;
      return;
    }
  Int_t patch = AliHLTTPCTransform::GetPatch(row);
  AliHLTTPCTransform::Slice2Sector(0,row,sr,lr);
  AliHLTTPCTransform::Raw2Local(xyz,sr,lr,GetPadHit(row),GetTimeHit(row));
  fParSigmaY2[index] = AliHLTTPCTransform::GetParSigmaY2(row,xyz[2],GetCrossingAngleLUT(row));
  fParSigmaZ2[index] = AliHLTTPCTransform::GetParSigmaZ2(row,xyz[2],GetTgl());
  
  if(parametrize)
    {
      fParSigmaY2[index] = (fParSigmaY2[index] + (1./12)*pow(AliHLTTPCTransform::GetPadPitchWidth(patch),2) );
      fParSigmaY2[index] *= 0.108;
      if(patch<2)
	fParSigmaY2[index] *= 2.07;
     
      fParSigmaZ2[index] = (fParSigmaZ2[index] + (1./12)*pow(AliHLTTPCTransform::GetZWidth(),2) );
      fParSigmaZ2[index] *= 0.169;
      if(patch<2)
	fParSigmaZ2[index] *= 1.77;
    }
  
  //convert to raw coordinates:
  fParSigmaY2[index] /= pow(AliHLTTPCTransform::GetPadPitchWidth(patch),2);
  fParSigmaZ2[index] /= pow(AliHLTTPCTransform::GetZWidth(),2);
}

void AliHLTTPCModelTrack::SetCluster(Int_t row,Float_t fpad,Float_t ftime,Float_t charge,
				 Float_t sigmaY2,Float_t sigmaZ2,Int_t npads)
{
  AliHLTTPCClusterModel *cl = GetClusterModel(row);
  
  //First bit: Cluster is present or not
  //Second bit: Cluster was set, meaning an fit attempt was done (if true)
  
  cl->fPresent |= 0x2; //set second bit to true, because a fit attempt has been made
  
  Int_t patch = AliHLTTPCTransform::GetPatch(row);
  if(!charge || npads == 1)
    {
      cl->fPresent &= ~0x1; //set first bit to false
    }
  else
    {
      cl->fPresent|=0x1;//set first bit to true
      cl->fDTime = (ftime - GetTimeHit(row))/(AliHLTTPCDataCompressor::GetZResidualStep(row)/AliHLTTPCTransform::GetZWidth());   
      cl->fDPad = (fpad - GetPadHit(row))/(AliHLTTPCDataCompressor::GetXYResidualStep(row)/AliHLTTPCTransform::GetPadPitchWidth(patch));
      cl->fDCharge = charge;// - AliHLTTPCDataCompressor::GetClusterCharge();
      //cl->fSlice = fSlice;
      if(sigmaY2==0 && sigmaZ2==0)
	{
	  cl->fDSigmaY2=0;//if width is zero, shape is not supposed to be written
	  cl->fDSigmaZ2=0;
	}
      else
	{
	  cl->fDSigmaY2 = (sigmaY2 - GetParSigmaY2(row))/(AliHLTTPCDataCompressor::GetXYWidthStep()/pow(AliHLTTPCTransform::GetPadPitchWidth(patch),2));
	  cl->fDSigmaZ2 = (sigmaZ2 - GetParSigmaZ2(row))/(AliHLTTPCDataCompressor::GetZWidthStep()/pow(AliHLTTPCTransform::GetZWidth(),2));
	}
      cl->fNPads = npads;
    }
  
  //fNClusters++;
}

void AliHLTTPCModelTrack::Set(AliHLTTPCTrack *tpt)
{
  AliHLTTPCModelTrack *tr = (AliHLTTPCModelTrack*)tpt;
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
      std::cerr<<"AliHLTTPCModelTrack::Set : Init has already been called for this object!"<<std::endl;
      return;
    }

  //Init(tr->fSlice,tr->fPatch);
  Init(0,tr->fPatch);
  memcpy(fClusters,tr->fClusters,AliHLTTPCTransform::GetNRows(fPatch)*sizeof(AliHLTTPCClusterModel));
  memcpy(fPad,tr->fPad,AliHLTTPCTransform::GetNRows(fPatch)*sizeof(Float_t));
  memcpy(fTime,tr->fTime,AliHLTTPCTransform::GetNRows(fPatch)*sizeof(Float_t));
  memcpy(fParSigmaY2,tr->fParSigmaY2,AliHLTTPCTransform::GetNRows(fPatch)*sizeof(Float_t));
  memcpy(fParSigmaZ2,tr->fParSigmaZ2,AliHLTTPCTransform::GetNRows(fPatch)*sizeof(Float_t));
  memcpy(fCrossingAngle,tr->fCrossingAngle,AliHLTTPCTransform::GetNRows(fPatch)*sizeof(Float_t));
  memcpy(fTrackModel,tr->fTrackModel,sizeof(AliHLTTPCTrackModel));

}

Int_t AliHLTTPCModelTrack::GetNPresentClusters()
{
  //Return the number of assigned clusters to the track.
  //Differs from fNClusters, which should be equal to the 
  //number of padrows in the present patch.
  
  Int_t count=0;

  for(Int_t i=AliHLTTPCTransform::GetFirstRow(fPatch); i<=AliHLTTPCTransform::GetLastRow(fPatch); i++)
    if(IsPresent(i))
      count++;

  return count;
}

void AliHLTTPCModelTrack::FillModel()
{
  //Fill the track structure
  
  if(fNClusters != AliHLTTPCTransform::GetNRows(fPatch))
    {
#if 0
      std::cout<<"AliHLTTPCModelTrack::FillModel : fNClusters != nrows; beware, this could be caused by a bug!!!"<<std::endl;
#endif
      fNClusters = AliHLTTPCTransform::GetNRows(fPatch);
    }

  if(!fTrackModel)
    {
      std::cerr<<"AliHLTTPCModelTrack::FillModel() : No trackmodel "<<std::endl;
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

void AliHLTTPCModelTrack::FillTrack()
{
  //Fill the track parameters from the structure.
  
  if(!fTrackModel)
    {
      std::cerr<<"AliHLTTPCModelTrack::FillTrack() : No data!!"<<std::endl;
      return;
    }
  SetKappa(fTrackModel->fKappa);
  SetCharge((-1*(Int_t)copysign(1.,GetKappa())));
  SetFirstPoint(fTrackModel->fFirstPointX,fTrackModel->fFirstPointY,fTrackModel->fFirstPointZ);
  SetTgl(fTrackModel->fTgl);
  SetPsi(fTrackModel->fPsi);
  //SetLength(fTrackModel->fLength);
  fNClusters = AliHLTTPCTransform::GetNRows(fPatch);//fTrackModel->fNClusters;
  SetPt((AliHLTTPCTransform::GetBFact()*AliHLTTPCTransform::GetBField())/fabs(GetKappa()));
    
  CalculateHelix();

  for(Int_t i=AliHLTTPCTransform::GetFirstRow(fPatch); i<=AliHLTTPCTransform::GetLastRow(fPatch); i++)
    {
      AliHLTTPCClusterModel *cl = GetClusterModel(i);
      if(!cl) continue;

      if(cl->fSlice == -1)
	{
	  SetPadHit(i,-1);
	  SetTimeHit(i,-1);
	  continue;
	}
      if(cl->fSlice < 0 || cl->fSlice > 35)
	{
	  std::cerr<<"AliHLTTPCModelTrack::FillTrack : Slice out of range "<<cl->fSlice<<" on row "<<i<<std::endl;
	  exit(5);
	}
      
      Float_t angle = 0;
      
      AliHLTTPCTransform::Local2GlobalAngle(&angle,cl->fSlice);
      if(!CalculateReferencePoint(angle,AliHLTTPCTransform::Row2X(i)))
	{
	  if(IsPresent(i))
	    {
	      std::cerr<<"AliHLTTPCModelTrack::FillTrack : Track does not cross slice "<<cl->fSlice<<" row "<<i<<" Points "
		  <<GetPointX()<<" "<<GetPointY()<<" "<<GetPointZ()<<std::endl;
	      Print();
	      exit(5);
	    }
	  SetPadHit(i,-1);
	  SetTimeHit(i,-1);
	  continue;
	}
      Float_t hit[3] = {GetPointX(),GetPointY(),GetPointZ()};
      Int_t sector,row;
      AliHLTTPCTransform::Slice2Sector(cl->fSlice,i,sector,row);
      AliHLTTPCTransform::Global2Raw(hit,sector,row);

      SetPadHit(i,hit[1]);
      SetTimeHit(i,hit[2]);

      Float_t crossingangle = GetCrossingAngle(i,cl->fSlice);
      
      SetCrossingAngleLUT(i,crossingangle);
      CalculateClusterWidths(i,kTRUE);
      
    }
}

void AliHLTTPCModelTrack::SetPadHit(Int_t row,Float_t pad)
{
  Int_t index = row-AliHLTTPCTransform::GetFirstRow(fPatch);
  if(index < 0 || index > AliHLTTPCTransform::GetNRows(fPatch))
    {
      std::cerr<<"AliHLTTPCModelTrack::SetPadHit() : Wrong index: "<<index<<std::endl;
      return;
    }
  fPad[index]=pad;
}

void AliHLTTPCModelTrack::SetTimeHit(Int_t row,Float_t time)
{
  Int_t index = row-AliHLTTPCTransform::GetFirstRow(fPatch);
  if(index < 0 || index > AliHLTTPCTransform::GetNRows(fPatch))
    {
      std::cerr<<"AliHLTTPCModelTrack::SetTimeHit() : Wrong index: "<<index<<std::endl;
      return;
    }
  fTime[index]=time;
}

void AliHLTTPCModelTrack::SetCrossingAngleLUT(Int_t row,Float_t angle)
{
  Int_t index = row-AliHLTTPCTransform::GetFirstRow(fPatch);
  if(index < 0 || index > AliHLTTPCTransform::GetNRows(fPatch))
    {
      std::cerr<<"AliHLTTPCModelTrack::SetCrossingAngle() : Wrong index: "<<index<<std::endl;
      return;
    }
  fCrossingAngle[index]=angle;
}

void AliHLTTPCModelTrack::SetOverlap(Int_t row,Int_t id)
{

  Int_t index = row-AliHLTTPCTransform::GetFirstRow(fPatch);
  if(index < 0 || index > AliHLTTPCTransform::GetNRows(fPatch))
    {
      std::cerr<<"AliHLTTPCModelTrack::SetOverlap() : Wrong index: "<<index<<std::endl;
      return;
    }
  if(fNoverlaps[index] >= fMaxOverlaps) return;
  fOverlap[index][fNoverlaps[index]++] = id;
}

Bool_t AliHLTTPCModelTrack::IsPresent(Int_t row)
{
  AliHLTTPCClusterModel *cl = GetClusterModel(row);
  return (Bool_t)(cl->fPresent & 0x1);
}

Bool_t AliHLTTPCModelTrack::IsSet(Int_t row)
{
  AliHLTTPCClusterModel *cl = GetClusterModel(row);
  return (Bool_t)(cl->fPresent & 0x2);
}

Int_t AliHLTTPCModelTrack::GetNPads(Int_t row)
{
  AliHLTTPCClusterModel *cl = GetClusterModel(row);
  return cl->fNPads;
}

Bool_t AliHLTTPCModelTrack::GetPad(Int_t row,Float_t &pad)
{
  //(ftime - GetTimeHit(fNClusters))/AliHLTTPCDataCompressor::GetXYResidualStep();
  //(fpad - GetPadHit(fNClusters))/AliHLTTPCDataCompressor::GetZResidualStep();

  AliHLTTPCClusterModel *cl = GetClusterModel(row);
  Int_t patch = AliHLTTPCTransform::GetPatch(row);
  pad = cl->fDPad*(AliHLTTPCDataCompressor::GetXYResidualStep(row)/AliHLTTPCTransform::GetPadPitchWidth(patch)) + GetPadHit(row);
  
  return IsPresent(row);
}

Bool_t AliHLTTPCModelTrack::GetTime(Int_t row,Float_t &time)
{
  AliHLTTPCClusterModel *cl = GetClusterModel(row);
  time = cl->fDTime*(AliHLTTPCDataCompressor::GetZResidualStep(row)/AliHLTTPCTransform::GetZWidth()) + GetTimeHit(row);
  
  return IsPresent(row);
}

Bool_t AliHLTTPCModelTrack::GetClusterCharge(Int_t row,Int_t &charge)
{
  AliHLTTPCClusterModel *cl = GetClusterModel(row);
  charge = (Int_t)cl->fDCharge;// + AliHLTTPCDataCompressor::GetClusterCharge();
  
  return IsPresent(row);
}

Bool_t AliHLTTPCModelTrack::GetXYWidth(Int_t row,Float_t &width)
{
  //cl->fDSigmaY2 = (sigmaY2 - GetParSigmaY2(row))/AliHLTTPCDataCompressor::GetXYWidthStep();
  AliHLTTPCClusterModel *cl = GetClusterModel(row);
  Int_t patch = AliHLTTPCTransform::GetPatch(row);
  width = cl->fDSigmaY2*(AliHLTTPCDataCompressor::GetXYWidthStep()/pow(AliHLTTPCTransform::GetPadPitchWidth(patch),2)) + GetParSigmaY2(row);

  return IsPresent(row);
}

Bool_t AliHLTTPCModelTrack::GetZWidth(Int_t row,Float_t &width)
{
  AliHLTTPCClusterModel *cl = GetClusterModel(row);
  width = cl->fDSigmaZ2*(AliHLTTPCDataCompressor::GetZWidthStep()/pow(AliHLTTPCTransform::GetZWidth(),2)) + GetParSigmaZ2(row);

  return IsPresent(row);
}

Bool_t AliHLTTPCModelTrack::GetPadResidual(Int_t row,Float_t &res)
{
  AliHLTTPCClusterModel *cl = GetClusterModel(row);
  res = cl->fDPad;
  return IsPresent(row);
}

Bool_t AliHLTTPCModelTrack::GetTimeResidual(Int_t row,Float_t &res)
{
  AliHLTTPCClusterModel *cl = GetClusterModel(row);
  res = cl->fDTime;
  return IsPresent(row);
}

Bool_t AliHLTTPCModelTrack::GetXYWidthResidual(Int_t row,Float_t &res)
{
  AliHLTTPCClusterModel *cl = GetClusterModel(row);
  res = cl->fDSigmaY2;
  return IsPresent(row);
}

Bool_t AliHLTTPCModelTrack::GetZWidthResidual(Int_t row,Float_t &res)
{
  AliHLTTPCClusterModel *cl = GetClusterModel(row);
  res = cl->fDSigmaZ2;
  return IsPresent(row);
}

Int_t AliHLTTPCModelTrack::GetSlice(Int_t row)
{
  AliHLTTPCClusterModel *cl = GetClusterModel(row);
  return cl->fSlice;
}

Float_t AliHLTTPCModelTrack::GetPadHit(Int_t row)
{
  Int_t index = row-AliHLTTPCTransform::GetFirstRow(fPatch);
  if(index < 0 || index > AliHLTTPCTransform::GetNRows(fPatch))
    {
      std::cerr<<"AliHLTTPCModelTrack::GetPadHit() : Wrong index: "<<index<<" row "<<row<<std::endl;
      return 0;
    }
  return fPad[index];
}

Float_t AliHLTTPCModelTrack::GetTimeHit(Int_t row)
{
  Int_t index = row-AliHLTTPCTransform::GetFirstRow(fPatch);
  if(index < 0 || index > AliHLTTPCTransform::GetNRows(fPatch))
    {
      std::cerr<<"AliHLTTPCModelTrack::GetTimeHit() : Wrong index: "<<index<<" row "<<row<<std::endl;
      return 0;
    }
  return fTime[index];
}

Float_t AliHLTTPCModelTrack::GetCrossingAngleLUT(Int_t row)
{
  Int_t index = row-AliHLTTPCTransform::GetFirstRow(fPatch);
  if(index < 0 || index > AliHLTTPCTransform::GetNRows(fPatch))
    {
      std::cerr<<"AliHLTTPCModelTrack::GetCrossingAngleLUT() : Wrong index: "<<index<<" row "<<row<<std::endl;
      return 0;
    }
  return fCrossingAngle[index];
}

Float_t AliHLTTPCModelTrack::GetParSigmaY2(Int_t row)
{
  Int_t index = row-AliHLTTPCTransform::GetFirstRow(fPatch);
  if(index < 0 || index > AliHLTTPCTransform::GetNRows(fPatch))
    {
      std::cerr<<"AliHLTTPCModelTrack::GetParSigmaY2() : Wrong index: "<<index<<" row "<<row<<std::endl;
      return 0;
    }
  return fParSigmaY2[index];
}

Float_t AliHLTTPCModelTrack::GetParSigmaZ2(Int_t row)
{
  Int_t index = row-AliHLTTPCTransform::GetFirstRow(fPatch);
  if(index < 0 || index > AliHLTTPCTransform::GetNRows(fPatch))
    {
      std::cerr<<"AliHLTTPCModelTrack::GetParSigmaZ2() : Wrong index: "<<index<<" row "<<row<<std::endl;
      return 0;
    }
  return fParSigmaZ2[index];
}

Int_t AliHLTTPCModelTrack::GetNOverlaps(Int_t row)
{
  Int_t index = row - AliHLTTPCTransform::GetFirstRow(fPatch);
  if(index < 0 || index > AliHLTTPCTransform::GetNRows(fPatch))
    {
      std::cerr<<"AliHLTTPCModelTrack::GetOverlap() : Wrong index: "<<index<<std::endl;
      return 0;
    }
  return fNoverlaps[index];
}

Int_t *AliHLTTPCModelTrack::GetOverlaps(Int_t row)
{
  Int_t index = row - AliHLTTPCTransform::GetFirstRow(fPatch);
  if(index < 0 || index > AliHLTTPCTransform::GetNRows(fPatch))
    {
      std::cerr<<"AliHLTTPCModelTrack::GetOverlap() : Wrong index: "<<index<<std::endl;
      return 0;
    }
  return fOverlap[index];
}

AliHLTTPCClusterModel *AliHLTTPCModelTrack::GetClusterModel(Int_t row)
{
  if(!fClusters) return 0; 
  Int_t index = row-AliHLTTPCTransform::GetFirstRow(fPatch);
  if(index < 0 || index > AliHLTTPCTransform::GetNRows(fPatch))
    {
      std::cerr<<"AliHLTTPCModelTrack::GetClusterModel() : Wrong index: "<<index<<std::endl;
      return 0;
    }
  return &fClusters[index];
}

void AliHLTTPCModelTrack::Print(Bool_t everything)
{
  //Print info

#if 0
  std::cout<<"First point "<<GetFirstPointX()<<" "<<GetFirstPointY()<<" "<<GetFirstPointZ()<<std::endl;
  std::cout<<"Last point "<<GetLastPointX()<<" "<<GetLastPointY()<<" "<<GetLastPointZ()<<std::endl;
  std::cout<<"Pt "<<GetPt()<<" kappa "<<GetKappa()<<" tgl "<<GetTgl()<<" psi "<<GetPsi()<<" charge "<<GetCharge()<<std::endl;
  std::cout<<"Center "<<GetCenterX()<<" "<<GetCenterY()<<std::endl<<std::endl;
#endif
  if(!everything)
    return;
#if 0
  std::cout<<"NHits "<<GetNClusters()<<std::endl;
#endif

#if 0
  std::cout<<"Clusters:"<<std::endl;#
#endif
  Int_t origslice=-1,counter=0;
  Float_t fpad,ftime;
  for(Int_t i=AliHLTTPCTransform::GetFirstRow(fPatch); i<=AliHLTTPCTransform::GetLastRow(fPatch); i++)
    {
      AliHLTTPCClusterModel *cl = GetClusterModel(i);
      
      if(!IsPresent(i))
	  {
#if 0
	std::cout<<i<<" Empty"<<" Slice "<<cl->fSlice<<" Padcrossing "<<GetPadHit(i)<<" Timecrossing "<<GetTimeHit(i)<<" ";
#endif
	  }
      else
	{
	  GetPad(i,fpad);
	  GetTime(i,ftime);
	  if(counter==0)
	    origslice=cl->fSlice;
	  else if(cl->fSlice != origslice)
	      {
#if 0
	    std::cout<<"Change in slice "<<cl->fSlice<<" "<<origslice<<std::endl;
#endif
	      }
#if 0
	  std::cout<<i<<" Slice "<<cl->fSlice<<" Dpad "<<cl->fDPad<<" Dtime "<<cl->fDTime<<" Dcharge "<<cl->fDCharge;
	  std::cout<<" sigmaY2 "<<GetParSigmaY2(i)<<" sigmaZ2 "<<GetParSigmaZ2(i);
	  std::cout<<" Pad "<<fpad<<" padhit "<<GetPadHit(i)<<" Time "<<ftime<<" timehit "<<GetTimeHit(i)<<" ";
	  std::cout<<"Number of pads "<<GetNPads(i)<<" Overlaps "<<GetNOverlaps(i);
#endif
	  counter++;
	}
#if 0
      std::cout<<std::endl;
#endif
    }
}

void AliHLTTPCModelTrack::SetClusterLabel(Int_t row,Int_t *trackID)
{
#ifdef do_mc
  AliHLTTPCClusterModel *cl = GetClusterModel(row);
  cl->fTrackID[0] = trackID[0];
  cl->fTrackID[1] = trackID[1];
  cl->fTrackID[2] = trackID[2];
#endif
  return;
}

void AliHLTTPCModelTrack::GetClusterLabel(Int_t row,Int_t *trackID)
{
#ifdef do_mc
  AliHLTTPCClusterModel *cl = GetClusterModel(row);
  trackID[0] = cl->fTrackID[0];
  trackID[1] = cl->fTrackID[1];
  trackID[2] = cl->fTrackID[2];
#endif
  return;
}

