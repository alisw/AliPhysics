/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

// $Id$

#include "AliMUONClusterFinderCOG.h"

#include "AliLog.h"
#include "AliMUONCluster.h"
#include "AliMUONDigit.h"
#include "AliMUONPad.h"
#include "AliMUONPreClusterFinder.h"
#include "AliMpArea.h"
#include "TClonesArray.h"
#include "TVector2.h"

/// \class AliMUONClusterFinderCOG
///
/// A very basic (and mostly useless, probably) cluster finder.
/// 
/// We use AliMUONPreClusterFinder to actually build the cluster,
/// and then we simply use center-of-gravity to get the coordinates
/// of the cluster.
/// Only point to note is that we compute separately both
/// cathodes when doing, in order to take the positions from the 
/// direction with the better resolution.
///
/// \author Laurent Aphecetche

/// \nocond CLASSIMP
ClassImp(AliMUONClusterFinderCOG)
/// \endcond

//_____________________________________________________________________________
AliMUONClusterFinderCOG::AliMUONClusterFinderCOG()
: AliMUONVClusterFinder(),
fPreClusterFinder(0x0)
{
  /// ctor
}

//_____________________________________________________________________________
AliMUONClusterFinderCOG::~AliMUONClusterFinderCOG()
{
  /// dtor
  delete fPreClusterFinder;
}

//_____________________________________________________________________________
Bool_t 
AliMUONClusterFinderCOG::Prepare(const AliMpVSegmentation* segmentations[2],
                                       TClonesArray* digits[2])
{
  /// Prepare for clustering
  
  // Find out the DetElemId
  Int_t detElemId(-1);
  
  for ( Int_t i = 0; i < 2; ++i )
  {
    AliMUONDigit* d = static_cast<AliMUONDigit*>(digits[i]->First());
    if (d)
    {
      detElemId = d->DetElemId();
      break;
    }
  }
  
  if ( detElemId < 0 )
  {
    AliWarning("Could not find DE. Probably no digits at all ?");
    return kFALSE;
  }
  
  delete fPreClusterFinder;
  fPreClusterFinder = new AliMUONPreClusterFinder;
  return fPreClusterFinder->Prepare(segmentations,digits);
}

//_____________________________________________________________________________
AliMUONCluster* 
AliMUONClusterFinderCOG::NextCluster()
{
  /// Get next cluster
  
  if ( !fPreClusterFinder ) return 0x0;
  AliMUONCluster* cluster = fPreClusterFinder->NextCluster();
  if ( cluster )
  {
    ComputePosition(*cluster);

    if ( cluster->Charge() < 7 )
    {
      // skip that one
      return NextCluster();
    }    
  }
  return cluster;
}

//_____________________________________________________________________________
void 
AliMUONClusterFinderCOG::ComputePosition(AliMUONCluster& cluster)
{
  /// Compute a first estimate of cluster position by a basic center-of-gravity
  
  Double_t xmin = 1E9;
  Double_t ymin = 1E9;
  Double_t xmax = -1E9;
  Double_t ymax = -1E9;
  
  Double_t x[] = { 0.0, 0.0 };
  Double_t y[] = { 0.0, 0.0 };
  
  Double_t xsize[] = { 0.0, 0.0 } ;
  Double_t ysize[] = { 0.0, 0.0 } ;
  
  for ( Int_t cathode = 0; cathode < 2; ++cathode )
  {
    for ( Int_t i = 0; i < cluster.Multiplicity(); ++i )
    {
      AliMUONPad* pad = cluster.Pad(i);
      TVector2 padPosition = pad->Position();
      AliMpArea area(pad->Position(),pad->Dimensions());
      xmin = std::min(area.LeftBorder(),xmin);
      xmax = std::max(area.RightBorder(),xmax);
      ymin = std::min(area.DownBorder(),ymin);
      ymax = std::max(area.UpBorder(),ymax);
      if ( cathode == pad->Cathode() )
      {
        x[cathode] += padPosition.X()*pad->Charge();
        y[cathode] += padPosition.Y()*pad->Charge();
        xsize[cathode] += pad->Dimensions().X();
        ysize[cathode] += pad->Dimensions().Y();
      }
    }
    if ( cluster.Charge(cathode) )
    {
      x[cathode] /= cluster.Charge(cathode);
      y[cathode] /= cluster.Charge(cathode);
    }
    if ( cluster.Multiplicity(cathode) )
    {
      xsize[cathode] /= cluster.Multiplicity(cathode);
      ysize[cathode] /= cluster.Multiplicity(cathode);
    }
  }
  
  Double_t xCOG = 0;
  Double_t yCOG = 0;

  // take the positions from the direction with the better resolution
  xCOG = ( xsize[0] < xsize[1] ) ? x[0] : x[1];
  yCOG = ( ysize[0] < ysize[1] ) ? y[0] : y[1];
  
  AliDebug(1,Form("Cluster mult %d (x,y)=(%e,%e) boundaries=(xmin,ymin,xmax,ymax)=(%e,%e,%e,%e)"
                  " (x0,y0,x1,y1)=(%e,%e,%e,%e) ",
                  cluster.Multiplicity(),xCOG,yCOG,xmin,ymin,xmax,ymax,
                  x[0],y[0],x[1],y[1]));
  
  cluster.SetPosition(TVector2(xCOG,yCOG),cluster.Pad(0)->Dimensions()); // FIXME: what to put as an error here ?
}



