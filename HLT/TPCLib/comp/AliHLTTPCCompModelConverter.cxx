// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Timm Steinbeck <timm@kip.uni-heidelberg.de>           *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** @file   AliHLTTPCCompModelConverter.cxx
    @author Timm Steinbeck
    @author changed by J. Wagner
    @date   17-11-2007
    @brief  A copy processing component for the HLT. */

#if __GNUC__ >= 3
using namespace std;
#endif

#include "AliHLTTPCCompModelConverter.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCTrack.h"
#include "AliHLTTPCModelTrack.h"
#include "AliHLTTPCCompDataCompressorHelper.h"
#include <cerrno>

AliHLTTPCCompModelConverter::AliHLTTPCCompModelConverter():
  fInputTrackArray(),
  fOutputTrackArray("AliHLTTPCModelTrack"),
  fModelAnalysisInstance(NULL),
  fMinHits(0)
    {
      // see header file for class documentation
      for ( UInt_t slice=0; slice<36; slice++ )
	for ( UInt_t patch=0; patch<6; patch++ )
	  {
	    fClusterUsedSizes[slice][patch] = 0;
	    fClusterUsed[slice][patch] = NULL;
	  }
      Init();
      fMinHits = 5;
    }

AliHLTTPCCompModelConverter::AliHLTTPCCompModelConverter(AliHLTTPCCompModelAnalysis* modelanalysis):
  fInputTrackArray(),
  fOutputTrackArray("AliHLTTPCModelTrack"),
  fModelAnalysisInstance(modelanalysis),
  fMinHits(0)
    {
      // see header file for class documentation
      for ( UInt_t slice=0; slice<36; slice++ )
	for ( UInt_t patch=0; patch<6; patch++ )
	  {
	    fClusterUsedSizes[slice][patch] = 0;
	    fClusterUsed[slice][patch] = NULL;
	  }
      Init();
      fMinHits = 5;
    }

AliHLTTPCCompModelConverter::~AliHLTTPCCompModelConverter()
    {
      // see header file for class documentation
      for ( UInt_t slice=0; slice<36; slice++ )
	for ( UInt_t patch=0; patch<6; patch++ )
	  {
	    if ( fClusterUsed[slice][patch] )
	      {
		delete [] fClusterUsed[slice][patch];
		fClusterUsed[slice][patch] = NULL;
	      }
	  }
    }

int AliHLTTPCCompModelConverter::Init()
    {
      // see header file for class documentation
      fInputTrackArray.Reset();
      fOutputTrackArray.Reset();
      for ( UInt_t slice=0; slice<36; slice++ )
	for ( UInt_t patch=0; patch<6; patch++ )
	  fClusters[slice][patch] = NULL;
      
      return 0;
    }

int AliHLTTPCCompModelConverter::SetInputTracks( AliHLTTPCTrackletData* tracklets )
    {
      // see header file for class documentation
      HLTDebug( "Filling %u tracks", (unsigned)tracklets->fTrackletCnt );
      fInputTrackArray.FillTracks( tracklets->fTrackletCnt, tracklets->fTracklets );
      return 0;
    }

int AliHLTTPCCompModelConverter::SetInputClusters( AliHLTTPCClusterData* clusters, UInt_t slice, UInt_t patch )
    {
      // see header file for class documentation
      if ( slice>=36 || patch>=6 )
	return EINVAL;
      if ( fClusters[slice][patch] )
	return EBUSY;
      fClusters[slice][patch] = clusters;
      if ( fClusterUsedSizes[slice][patch]<clusters->fSpacePointCnt ||
	   fClusterUsedSizes[slice][patch]>clusters->fSpacePointCnt*8 )
	{
	  delete [] fClusterUsed[slice][patch];
	  fClusterUsed[slice][patch] = NULL;
	}
      if ( !fClusterUsed[slice][patch] )
	{
	  fClusterUsed[slice][patch] = new bool[clusters->fSpacePointCnt];
	  if ( !fClusterUsed[slice][patch] )
	    {
	      HLTDebug( "Out of memory trying to allocate usage data for  %u clusters", (unsigned)clusters->fSpacePointCnt );
	      return ENOMEM;
	    }
	}
      for ( unsigned long nn=0; nn<clusters->fSpacePointCnt; nn++ )
	fClusterUsed[slice][patch][nn]=false;
      HLTDebug( "Filling %u clusters", (unsigned)clusters->fSpacePointCnt );
      return 0;
    }

void AliHLTTPCCompModelConverter::Convert()
    {
      // see header file for class documentation
      fInputTrackArray.QSort();
      for(Int_t i=0; i<fInputTrackArray.GetNTracks(); i++)
	{
	  AliHLTTPCTrack *intrack = fInputTrackArray.GetCheckedTrack(i);
	  
	  // NO WARNING IF intrack = NULL!
	  if(!intrack) continue;
	  
	  if((unsigned)intrack->GetNHits()<fMinHits) 
	    {
	      HLTDebug("Track %d with %d clusters is below minimum of %d clusters",i,intrack->GetNHits(),fMinHits);
	      break;
	    };
	  
	  // LOSS OF TRACKS due to following statement possible!
	  if(intrack->GetPt()<0.1)
	    {
	      HLTDebug("Discarding track with low pt.");
	      if(fModelAnalysisInstance)
		{
		  if(fModelAnalysisInstance->GetfModelAnalysis()) // analysis of model
		    {
		      fModelAnalysisInstance->MarkTrashTrack(intrack);
		    }
		}
	      
	      continue;
	    }
	  
	  intrack->CalculateHelix();
	  
	  AliHLTTPCModelTrack *outtrack = (AliHLTTPCModelTrack*)fOutputTrackArray.NextTrack();
	  outtrack->SetNHits(intrack->GetNHits());
	  outtrack->SetRowRange(intrack->GetFirstRow(),intrack->GetLastRow());
	  outtrack->SetFirstPoint(intrack->GetFirstPointX(),intrack->GetFirstPointY(),intrack->GetFirstPointZ());
	  outtrack->SetLastPoint(intrack->GetLastPointX(),intrack->GetLastPointY(),intrack->GetLastPointZ());
	  outtrack->SetPt(intrack->GetPt());
	  outtrack->SetPsi(intrack->GetPsi());
	  outtrack->SetTgl(intrack->GetTgl());
	  outtrack->SetCharge(intrack->GetCharge());
	  outtrack->CalculateHelix();
	  Int_t nhits = intrack->GetNHits();
	  UInt_t *hitids = intrack->GetHitNumbers();
	  Int_t origslice = AliHLTTPCSpacePointData::GetSlice(hitids[nhits-1]);
	  outtrack->Init(origslice,-1);
	  
	  for(Int_t j=nhits-1; j>=0; j--)
	    {
	      UInt_t id=hitids[j];
	      Int_t slice = AliHLTTPCSpacePointData::GetSlice(id);
	      Int_t patch = AliHLTTPCSpacePointData::GetPatch(id);
	      UInt_t pos = AliHLTTPCSpacePointData::GetNumber(id);

	      //UInt_t size;
	      if ( !fClusters[slice][patch] )
		{
		  //HLTWarning( "No clusters for slice %d, patch %d", slice, patch );
		  continue;
		}
	    if ( !fClusterUsed[slice][patch] )
	      {
		HLTWarning( "No cluster used data for slice %d, patch %d", slice, patch );
		continue;
	      }
	    if ( fClusters[slice][patch]->fSpacePointCnt<=pos )
		{
		  HLTWarning( "Clusters position %d too large in slice %d, patch %d (%u max.)", pos,
			      slice, patch, fClusters[slice][patch]->fSpacePointCnt );
		  continue;
		}
	    
	    AliHLTTPCSpacePointData *points = fClusters[slice][patch]->fSpacePoints;
	    bool* clustersUsed = fClusterUsed[slice][patch];
	    Float_t xyz[3] = {points[pos].fX,points[pos].fY,points[pos].fZ};
	    Int_t padrow = points[pos].fPadRow;
	    
	    //Calculate the crossing point between track and padrow
	    Float_t angle = 0; //Perpendicular to padrow in local coordinates
	    AliHLTTPCTransform::Local2GlobalAngle(&angle,slice);
	    if(!intrack->CalculateReferencePoint(angle,AliHLTTPCTransform::Row2X(padrow)))
	      {
		HLTError( "AliHLTDataCompressor::FillData : Error in crossing point calc on slice %d, padrow %d", slice, padrow );
		break;
		//outtrack->Print(kFALSE);
		//exit(5);
	      }
	    
	    Float_t xyzCross[3] = {intrack->GetPointX(),intrack->GetPointY(),intrack->GetPointZ()};
	    
	    Int_t sector,row;
	    AliHLTTPCTransform::Slice2Sector(slice,padrow,sector,row);
	    AliHLTTPCTransform::Global2Raw(xyzCross,sector,row);
#if 1
	    AliHLTTPCTransform::Local2Raw(xyz,sector,row);
#else
	    AliHLTTPCTransform::Global2Raw(xyz,sector,row);
#endif
	    
	    outtrack->SetPadHit(padrow,xyzCross[1]);
	    outtrack->SetTimeHit(padrow,xyzCross[2]);

	    outtrack->SetCrossingAngleLUT(padrow,intrack->GetCrossingAngle(padrow,slice));
	    outtrack->CalculateClusterWidths(padrow,kTRUE); // calculates parSigmas (with parametrisation) in raw coordinates
	    //HLTInfo("angle %f", outtrack->GetCrossingAngleLUT(padrow));
	    //HLTInfo("parsigma %f",outtrack->GetParSigmaY2(padrow));
	    patch = AliHLTTPCTransform::GetPatch(padrow);
	    // sigmay in units of pads (quantisation!) 
	    Float_t sigmaY2 = points[pos].fSigmaY2 / pow(AliHLTTPCTransform::GetPadPitchWidth(patch),2);
	    //HLTInfo("sigmaY conv.: %f", points[pos].fSigmaY2);

	    //HLTInfo("parSigmaY2 = %f",AliHLTTPCTransform::GetParSigmaY2(padrow, xyzCross[2],angle));
	    //Float_t testsigma = 0.0;
	    //outtrack->GetSigmaY2(padrow, testsigma);
	    //HLTInfo("DSigmaY2 = %f",testsigma);

	    //HLTInfo("sigmaY2 float: %f",sigmaY2);
	    Float_t sigmaZ2 = points[pos].fSigmaZ2 / pow(AliHLTTPCTransform::GetZWidth(),2);
	    outtrack->SetCluster(padrow,xyz[1],xyz[2],points[pos].fCharge,sigmaY2,sigmaZ2,3);
	    //AliHLTTPCClusterModel* test1 = outtrack->GetClusterModel(padrow);
	    //HLTInfo("Dsigma %f",test1->fDSigmaY);
	    
	    //IMPORTANT: Set the slice in which cluster is, you need it in AliHLTTPCModelTrack::FillTrack!
	    outtrack->GetClusterModel(padrow)->fSlice=slice;
#ifdef MODELDEBUG
	    outtrack->GetClusterModel(padrow)->fID=points[pos].fID;
	    HLTDebug( "Track %d cluster for padrow %d ID: %u (0x%08X) - fSlice: %u", i, padrow, 
		      outtrack->GetClusterModel(padrow)->fID, outtrack->GetClusterModel(padrow)->fID,
		      (unsigned)outtrack->GetClusterModel(padrow)->fSlice );
#endif
	    //points[pos].fCharge = 0;//Mark this cluster as used.
	    clustersUsed[pos] = true;//Mark this cluster as used.
	    //fNusedClusters++;
	    } //end of clusters for each track

	//outtrack->SetNClusters(AliHLTTPCTransform::GetNRows(-1)); // Equivalent call in ExpandTrackData
	} // end of track-loop
    ExpandTrackData();

    // validation test for clusternumbers of tracks:
    //for(unsigned long jj = 0; jj < (unsigned long) fOutputTrackArray.GetNTracks(); jj++)
    //  {
    //	AliHLTTPCModelTrack *track = (AliHLTTPCModelTrack*)fOutputTrackArray.GetCheckedTrack(jj);
    //	Int_t nhits = track->GetNHits();
    //	HLTInfo("Number of clusters for track %lu is %d",jj, nhits);
    //  }

    //comp->WriteFile(fOutputTrackArray);

    }

void AliHLTTPCCompModelConverter::ExpandTrackData()
    {
      // see header file for class documentation
      //Loop over tracks and try to assign unused clusters.
      //Only clusters which are closer than the max. residual are taken.
      
      HLTDebug( "Expanding %lu tracks", (unsigned long)fOutputTrackArray.GetNTracks() );
      for(Int_t i=0; i<fOutputTrackArray.GetNTracks(); i++)
	{
	  AliHLTTPCModelTrack *track = (AliHLTTPCModelTrack*)fOutputTrackArray.GetCheckedTrack(i);
	  
	  if(!track) continue;
	  
	  // tracks that hit every row already cannot be expanded in the current model!
	  if(track->GetNHits() == AliHLTTPCTransform::GetNRows()) continue;
	  
	  Int_t nhits = track->GetNHits();
	  
	  // validation test
	  //HLTInfo("Before expansion: track %u with number of clusters %d", i, nhits);
	  
	  Int_t lastSlice=-1;
	  for(Int_t padrow=AliHLTTPCTransform::GetNRows()-1; padrow>=0; padrow--)
	    {
	      if(track->IsPresent(padrow))
		{
		  lastSlice = track->GetClusterModel(padrow)->fSlice;
		  continue;
		}
	      
	      if(lastSlice < 0) //the outer cluster is missing, so skip it - it will be written anyhow.
		continue;
	      
	      //Check the slice of the next padrow:
	      Int_t nextPadrow = padrow-1;
	      Int_t nextSlice = -1;
	      while(nextPadrow >=0)
		{
		  if(track->IsPresent(nextPadrow))
		    {
		      nextSlice = track->GetClusterModel(nextPadrow)->fSlice;
		      break;
		    }
		  nextPadrow--;
		}
	      if(nextSlice>=0)
		if(nextSlice != lastSlice)//The track crosses a slice boundary here
		  continue;
	      
	      //UInt_t size;
	      if ( !fClusters[lastSlice][0] )
		{
		  HLTWarning( "No clusters for slice %d, patch %d", lastSlice, 0 );
		  continue;
		}
	      if ( !fClusterUsed[lastSlice][0] )
		{
		  HLTWarning( "No cluster used data for slice %d, patch %d", lastSlice, 0 );
		  continue;
		}
	      AliHLTTPCSpacePointData *points = fClusters[lastSlice][0]->fSpacePoints;//->GetDataPointer(size);
	      bool* clustersUsed = fClusterUsed[lastSlice][0];
	      
	      Float_t globalangle = 0;
	      AliHLTTPCTransform::Local2GlobalAngle(&globalangle,lastSlice);
	      if(!track->CalculateReferencePoint(globalangle,AliHLTTPCTransform::Row2X(padrow)))
		continue;
	      Float_t xyzCross[3] = {track->GetPointX(),track->GetPointY(),track->GetPointZ()};
	      AliHLTTPCTransform::Global2LocHLT(xyzCross,lastSlice);
	      Float_t mindist = 123456789;
	      AliHLTTPCSpacePointData *closest=0;
	      UInt_t closestJ=0;
	      for(UInt_t j=0; j<fClusters[lastSlice][0]->fSpacePointCnt; j++)
		{
		  //if(points[j].fCharge == 0) continue;// || points[j].fPadRow != padrow) continue;
		  if (clustersUsed[j]) continue; // Cluster already used
		  if(points[j].fPadRow < padrow) continue;
		  if(points[j].fPadRow > padrow) break;
		  Float_t xyz[3] = {points[j].fX,points[j].fY,points[j].fZ};
#if 1
#else
		  AliHLTTPCTransform::Global2LocHLT(xyz,lastSlice);
#endif
		  
		  //Check for overflow:
		  Int_t temp = (Int_t)rint((xyzCross[1]-xyz[1])/AliHLTTPCCompDataCompressorHelper::GetXYResidualStep(padrow));
		  if( abs(temp) > 1<<(AliHLTTPCCompDataCompressorHelper::GetNPadBits()-1))
		    continue;
		  
		  temp = (Int_t)rint((xyzCross[2]-xyz[2])/AliHLTTPCCompDataCompressorHelper::GetZResidualStep(padrow));
		  if( abs(temp) > 1<<(AliHLTTPCCompDataCompressorHelper::GetNTimeBits()-1))
		    continue;
		  
		  Float_t dist = sqrt( pow(xyzCross[1]-xyz[1],2) + pow(xyzCross[2]-xyz[2],2) );
		  if(dist < mindist)
		    {
		      closest = &points[j];
		      closestJ = j;
		      mindist = dist;
		    }
		}
	      if(closest) //there was a cluster assigned
		{
		  Int_t sector,row;
		  Float_t xyz[3] = {closest->fX,closest->fY,closest->fZ};
		  AliHLTTPCTransform::Slice2Sector(lastSlice,padrow,sector,row);
		  AliHLTTPCTransform::Local2Raw(xyzCross,sector,row);
#if 1
		  AliHLTTPCTransform::Local2Raw(xyz,sector,row);
#else
		  AliHLTTPCTransform::Global2Raw(xyz,sector,row);
#endif
		  
		  track->SetPadHit(padrow,xyzCross[1]);
		  track->SetTimeHit(padrow,xyzCross[2]);
		  
		  Float_t angle = track->GetCrossingAngle(padrow,lastSlice);
		  track->SetCrossingAngleLUT(padrow,angle);
		  track->CalculateClusterWidths(padrow,kTRUE);
		  Int_t patch = AliHLTTPCTransform::GetPatch(padrow);
		  Float_t sigmaY2 = closest->fSigmaY2 / pow(AliHLTTPCTransform::GetPadPitchWidth(patch),2);
		  Float_t sigmaZ2 = closest->fSigmaZ2 / pow(AliHLTTPCTransform::GetZWidth(),2);
		  track->SetCluster(padrow,xyz[1],xyz[2],closest->fCharge,sigmaY2,sigmaZ2,3);
		  //AliHLTTPCClusterModel* test1 = track->GetClusterModel(padrow);
		  //HLTInfo("Dsigma %f",test1->fDSigmaY);
		  
		  nhits++;
		  
		  //IMPORTANT: Set the slice in which cluster is, you need it in AliHLTTPCModelTrack::FillTrack!
		  track->GetClusterModel(padrow)->fSlice=lastSlice;
#ifdef MODELDEBUG
		  track->GetClusterModel(padrow)->fID=closest->fID;
		  HLTDebug( "Track %d cluster for padrow %d ID: %u (0x%08X) - fSlice: %u", i, padrow, 
			    track->GetClusterModel(padrow)->fID, track->GetClusterModel(padrow)->fID,
			    track->GetClusterModel(padrow)->fSlice );
#endif
		  //closest->fCharge = 0;//Mark this cluster as used.
		  clustersUsed[closestJ] = true;//Mark this cluster as used.
		}
	    }
	  track->SetNClusters(AliHLTTPCTransform::GetNRows());
	  //cout<<"Track was assigned "<<nhits<<" clusters"<<endl;
	  
	  // validation test
	  //HLTInfo( "After expansion: track %d with clusters %u", i, nhits);
	}
  
    }

unsigned long AliHLTTPCCompModelConverter::GetOutputModelDataSize()
    {
      // see header file for class documentation
      unsigned long dataSize=0;
      Short_t ntracks = fOutputTrackArray.GetNTracks();
      
      dataSize += sizeof(AliHLTUInt32_t);
      
      for(Int_t i=0; i<ntracks; i++)
	{
	  AliHLTTPCModelTrack *track = (AliHLTTPCModelTrack*)fOutputTrackArray.GetCheckedTrack(i);
	  if ( !track )
	    continue;
	  
	  dataSize += sizeof(AliHLTTPCTrackModel)+track->GetNClusters()*sizeof(AliHLTTPCClusterModel);
	}
      return dataSize;
    }

int AliHLTTPCCompModelConverter::OutputModelData( AliHLTUInt8_t* data )
    {
      // see header file for class documentation 
      unsigned long dataOffset=0;
      Short_t ntracks = fOutputTrackArray.GetNTracks();
      
      AliHLTTPCClusterModel *clusters=0;
      AliHLTTPCTrackModel *model=0;
      
      *(AliHLTUInt32_t*)data = 0; // Write format version number
      dataOffset += sizeof(AliHLTUInt32_t);
      
      for(Int_t i=0; i<ntracks; i++)
	{
	  AliHLTTPCModelTrack *track = (AliHLTTPCModelTrack*)fOutputTrackArray.GetCheckedTrack(i);
	  if ( !track )
	    continue;
	  
	  track->FillModel();
	  model = track->GetModel();
	  
	  clusters = track->GetClusters();
	  
	  // validation test
	  //HLTInfo( "Track %d clusters: %u", i, (unsigned)track->GetNPresentClusters() );
	  
	  for ( Int_t jj=0; jj<track->GetNClusters(); jj++ )
	    {
	      //HLTDebug( "  Cluster %d fPresent: %u", jj, (unsigned)clusters[jj].fPresent );
	    }
	  
	  memcpy( data+dataOffset, model, sizeof(AliHLTTPCTrackModel) );
	  dataOffset += sizeof(AliHLTTPCTrackModel);
	  
	  memcpy( data+dataOffset, clusters, track->GetNClusters()*sizeof(AliHLTTPCClusterModel) );
	  dataOffset += track->GetNClusters()*sizeof(AliHLTTPCClusterModel);
	}
      return 0;
    }

void AliHLTTPCCompModelConverter::SelectRemainingClusters()
    {
      // see header file for class documentation
      //Select which remaining clusters to write in addition to the compressed data.
      //In particular one can here make sure that "important" clusters are not missed:
      //The offline track finder perform seed finding in the outer padrows;
      //the first seeding is using pair of points on outermost padrow and 
      //0.125*nrows more rows towards the vertex. The second seeding uses pair
      //of points on the outermost padrow-0.5*0.125*nrows and 0.125*nrows + 0.5*0.125*nrows
      //more rows towards the vertex. In order to evaluate the seeds, the track offline
      //track finder checks whether a certain amount of possible clusters (padrows) is 
      //attached to the track, and then the kalman filtering starts.
      //To ensure a minimal loss off efficiency, all clusters in this region should be
      //intact.....
      
      Int_t nrows = AliHLTTPCTransform::GetNRows();
      Int_t gap=(Int_t)(0.125*nrows), shift=(Int_t)(0.5*gap);
      
      for(Int_t slice=0; slice<36; slice++)
	{
	  for(Int_t patch=0; patch<6; patch++)
	    {
	      if ( !fClusters[slice][patch] )
		continue;
	      AliHLTTPCSpacePointData *points = fClusters[slice][patch]->fSpacePoints;
	      bool* clustersUsed = fClusterUsed[slice][patch];
	      for(UInt_t i=0; i<fClusters[slice][patch]->fSpacePointCnt; i++)
		{
		  //if(points[i].fCharge == 0) continue; //Already removed
		  if (clustersUsed[i]) continue; //Already removed
		  Int_t padrow = (Int_t)points[i].fPadRow;
		  
		  //Check the widths (errors) of the cluster, and remove big bastards:
		  Float_t padw = sqrt(points[i].fSigmaY2) / AliHLTTPCTransform::GetPadPitchWidth(AliHLTTPCTransform::GetPatch(padrow));
		  Float_t timew = sqrt(points[i].fSigmaZ2) / AliHLTTPCTransform::GetZWidth();
		  if(padw >= 2.55 || timew >= 2.55)//Because we use 1 byte to store
		    {
		      //points[i].fCharge = 0;
		      clustersUsed[i] = true;
		      continue;
		    }
		  
		  Float_t xyz[3] = {points[i].fX,points[i].fY,points[i].fZ};
		  Int_t sector,row;
		  AliHLTTPCTransform::Slice2Sector(slice,padrow,sector,row);
		  AliHLTTPCTransform::Global2Raw(xyz,sector,row);
		  
		  if(padrow >= nrows-1-gap-shift) continue;//save all the clusters in this region
		  
		  //if(padrow >= nrows-1-shift) continue;
		  
		  //Save the clusters at the borders:
		  //if(xyz[1] < 3 || xyz[1] >= AliHLTTPCTransform::GetNPads(padrow)-4)
		  // continue;
		  
		  //Save clusters on padrows used for offline seeding:
		  if(padrow == nrows - 1 || padrow == nrows - 1 - gap ||                 //First seeding
		     padrow == nrows - 1 - shift || padrow == nrows - 1 - gap - shift)   //Second seeding
		    continue;
		  
		  //Cluster did not meet any of the above criteria, so disregard it:
		  //points[i].fCharge = 0;
		  clustersUsed[i] = true;
		}
	    }
	}
      
    }

unsigned long AliHLTTPCCompModelConverter::GetRemainingClustersOutputDataSize()
    {
      // see header file for class documentation
      int iResult=0;
#if 0
    for ( UInt_t slice=0; slice<36; slice++ )
	for ( UInt_t patch=0; patch<6; patch++ )
	    {
	    bool* clustersUsed = fClusterUsed[slice][patch];
	    if ( !clustersUsed || !fClusters[slice][patch] )
		continue;
	    for ( UInt_t pos=0; pos<fClusters[slice][patch]->fSpacePointCnt; pos++ )
		{
		if ( !clustersUsed[pos] )
		    clusterCnt++;
		}
	    }
    return clusterCnt*sizeof(AliHLTTPCClusterModel);
#else
    const Int_t nrows = AliHLTTPCTransform::GetNRows();
    Int_t * npoints = new Int_t[nrows];
    unsigned long dataWritten = 0;

    dataWritten += sizeof(AliHLTUInt32_t);

    // FIXME: get rid of hardcoded numbers
    for(Int_t slice=0; slice<35 && iResult>=0; slice++)
	{
	for(Int_t patch=0; patch < 6 && iResult>=0; patch++)
	    {
	    if ( !fClusters[slice][patch] )
		{
		dataWritten++;
		continue;
		}
	    AliHLTTPCSpacePointData *points = fClusters[slice][patch]->fSpacePoints;
	    bool* clustersUsed = fClusterUsed[slice][patch];
	    if ( !clustersUsed )
		continue;
	    memset(npoints,0,nrows*sizeof(Int_t));
	    Int_t nonZeroRows=0;
	  
	    for(UInt_t j=0; j<fClusters[slice][patch]->fSpacePointCnt; j++)
		{
		//if(points[j].fCharge == 0) continue; //has been used
		if ( clustersUsed[j] ) continue; //has been used
		if ( !npoints[points[j].fPadRow] )
		    nonZeroRows++;
		npoints[points[j].fPadRow]++;
		}

	    dataWritten++;

	    Int_t size =0;
	    Byte_t *data = 0;
	    AliHLTTPCRemainingRow *tempPt=0;
	    
	    Int_t lastRow = -2;
	    Int_t localcounter=0;
	    
	    for(UInt_t j=0; j<fClusters[slice][patch]->fSpacePointCnt; j++)
		{
		//if(points[j].fCharge == 0) continue; //has been used
		if ( clustersUsed[j] ) continue; //has been used
		
		Int_t padrow = points[j].fPadRow;
		if(padrow != lastRow)
		    {
		    if(lastRow != -2)
			{
			if(!tempPt)
			    {
			    HLTError( "Zero row pointer " );
			    iResult=-EINVAL;
			    break;
			    }
			if(localcounter != tempPt->fNClusters)
			    {
			    HLTError( "Mismatching clustercounter %lu - %d ", 
				      (unsigned long)localcounter, (Int_t)tempPt->fNClusters );
			    iResult=EINVAL;
			    break;
			    }
			dataWritten += size;
			}
		    if(data)
			delete [] data;
		    size = sizeof(AliHLTTPCRemainingRow) + npoints[padrow]*sizeof(AliHLTTPCRemainingCluster);
		    data = new Byte_t[size];
		    tempPt = reinterpret_cast<AliHLTTPCRemainingRow*>(data);
		    
		    localcounter=0;
		    tempPt->fPadRow = padrow;
		    tempPt->fNClusters = npoints[padrow];
		    lastRow = padrow;
		    }
		if(localcounter >= npoints[padrow])
		    {
		    HLTError( "Cluster counter out of range: %lu - %lu",
			      (unsigned long)localcounter, (unsigned long)npoints[padrow] );
		    iResult=-EINVAL;
		    break;
		    }
	      
		localcounter++;
		}
	    
	    //Write the last row:
	    if ( tempPt )
		{
		dataWritten += size;
		}
	    if(data)
	      delete [] data;
	    }
	}
    delete [] npoints;
    // FIXME check the caller and propagate an error condition
    if (iResult<0) return 0;
    return dataWritten;
#endif
    }

int AliHLTTPCCompModelConverter::GetRemainingClusters( AliHLTUInt8_t* const pTgt, unsigned long& dataSize )
    { 
      // see header file for class documentation
      int iResult=0;

      // FIXME: almost identical code to  GetRemainingClustersOutputDataSize
      // try to combine
      const Int_t nrows = AliHLTTPCTransform::GetNRows();
      Int_t * npoints = new Int_t[nrows];
      unsigned long dataWritten = 0;
      AliHLTUInt8_t* writePtr = pTgt;
      
      *(AliHLTUInt32_t*)writePtr = 0; // Write format version
      dataWritten += sizeof(AliHLTUInt32_t);
      writePtr += sizeof(AliHLTUInt32_t);

      for(Int_t slice=0; slice<=35 && iResult>=0; slice++)
	{
	  for(Int_t patch=0; patch < 6 && iResult>=0; patch++)
	    {
	      if ( !fClusters[slice][patch] )
		{
		  *writePtr = (AliHLTUInt8_t)0;
		  writePtr++;
		  dataWritten++;
		  continue;
		}
	      AliHLTTPCSpacePointData *points = fClusters[slice][patch]->fSpacePoints;
	      bool* clustersUsed = fClusterUsed[slice][patch];
	      if ( !clustersUsed )
		continue;
	      memset(npoints,0,nrows*sizeof(Int_t));
	      Int_t nonZeroRows=0;
	      
	      for(UInt_t j=0; j<fClusters[slice][patch]->fSpacePointCnt; j++)
		{
		  //if(points[j].fCharge == 0) continue; //has been used
		  if ( clustersUsed[j] ) continue; //has been used
		  if ( !npoints[points[j].fPadRow] )
		    nonZeroRows++;
		  npoints[points[j].fPadRow]++;
		}
	      
	      *writePtr = (AliHLTUInt8_t)nonZeroRows;
	      writePtr++;
	      dataWritten++;
	      
	      Int_t size =0;
	      Byte_t *data = 0;
	      AliHLTTPCRemainingRow *tempPt=0;
	      
	      Int_t lastRow = -2;
	      Int_t localcounter=0;
	    
	      for(UInt_t j=0; j<fClusters[slice][patch]->fSpacePointCnt; j++)
		{
		  //if(points[j].fCharge == 0) continue; //has been used
		  if ( clustersUsed[j] ) continue; //has been used
		  
		  Int_t padrow = points[j].fPadRow;
		  if(padrow != lastRow)
		    {
		      if(lastRow != -2)
			{
			  if(!tempPt)
			    {
			      HLTError( "Zero row pointer " );
			      iResult=-EINVAL;
			      break;
			    }
			  if(localcounter != tempPt->fNClusters)
			    {
			      HLTError( "Mismatching clustercounter %lu - %d ", 
					(unsigned long)localcounter, (Int_t)tempPt->fNClusters );
			      iResult=-EINVAL;
			      break;
			    }
			  //cout<<"Writing row "<<(int)tempPt->fPadRow<<" with "<<(int)tempPt->fNClusters<<" clusters"<<endl;
			  //fwrite(tempPt,size,1,outfile);
			  if ( dataWritten+size > dataSize )
			    {
			      HLTWarning( "Cannot write remaining clusters to output. Data size too large (exceeding %lu bytes)", (unsigned long)dataSize );
			      iResult=-ENOBUFS;
			      break;
			    }
			  memcpy( writePtr, tempPt, size );
			  dataWritten += size;
			  writePtr += size;
			}
		      if(data)
			delete [] data;
		      size = sizeof(AliHLTTPCRemainingRow) + npoints[padrow]*sizeof(AliHLTTPCRemainingCluster);
		      data = new Byte_t[size];
		      tempPt = (AliHLTTPCRemainingRow*)data;
		      
		      localcounter=0;
		      tempPt->fPadRow = padrow;
		      tempPt->fNClusters = npoints[padrow];
		      lastRow = padrow;
		    }
		  if(localcounter >= npoints[padrow])
		    {
		      HLTError( "Cluster counter out of range: %lu - %lu",
				(unsigned long)localcounter, (unsigned long)npoints[padrow] );
		      iResult=EINVAL;
		      break;
		    }
		  
		  Float_t xyz[3] = {points[j].fX,points[j].fY,points[j].fZ};
		  Int_t sector,row;
		  AliHLTTPCTransform::Slice2Sector(slice,padrow,sector,row);
#if 1
		  AliHLTTPCTransform::Local2Raw(xyz,sector,row);
#else
		  AliHLTTPCTransform::Global2Raw(xyz,sector,row);
#endif
		  
		  Float_t padw = points[j].fSigmaY2 / pow(AliHLTTPCTransform::GetPadPitchWidth(AliHLTTPCTransform::GetPatch(padrow)),2);
		  Float_t timew = points[j].fSigmaZ2 / pow(AliHLTTPCTransform::GetZWidth(),2);
		  tempPt->fClusters[localcounter].fPad = xyz[1];
		  tempPt->fClusters[localcounter].fTime = xyz[2];
		  tempPt->fClusters[localcounter].fCharge = points[j].fCharge;
		  tempPt->fClusters[localcounter].fSigmaY2 = padw;
		  tempPt->fClusters[localcounter].fSigmaZ2 = timew;
#ifdef MODELDEBUG
		  tempPt->fClusters[localcounter].fID = points[j].fID;
#endif
		  localcounter++;
		  if(fModelAnalysisInstance)
		    {
		      if(fModelAnalysisInstance->GetfModelAnalysis())
			{
			  fModelAnalysisInstance->MarkTrashCluster(fClusters[slice][patch], slice, patch);
			}
		    }
		}
	      
	      
	      //Write the last row:
	      if ( dataWritten+size > dataSize )
		{
		HLTWarning( "Cannot write remaining clusters to output. Data size too large (exceeding %lu bytes)", (unsigned long)dataSize );
		iResult=-ENOBUFS;
		if(data)
		  delete [] data;
		break;
		}
	      if ( tempPt )
		{
		  memcpy( writePtr, tempPt, size );
		  dataWritten += size;
		  writePtr += size;
		}
	      if(data)
		delete [] data;
	    }
	}
      dataSize = dataWritten;

      delete [] npoints;
      return iResult;
    }
