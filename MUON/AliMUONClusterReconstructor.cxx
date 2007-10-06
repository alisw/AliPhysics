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

/* $Id$ */

//-----------------------------------------------------------------------------
/// \class AliMUONClusterReconstructor
///
/// This class is just a steering class to loop over detection elements of tracking chambers,
/// extract the relevant digits, and pass them to the actual clusterizer class,
/// which operate on a single detection element at a time.
///
/// \author C. Finck and L. Aphecetche, Subatech
///
//-----------------------------------------------------------------------------

#include "AliMUONClusterReconstructor.h"

#include "AliMUONCluster.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONVCluster.h"
#include "AliMUONVClusterFinder.h"
#include "AliMUONVClusterStore.h"
#include "AliMUONVDigit.h"
#include "AliMUONVDigitStore.h"
#include "AliMUONPad.h"

#include "AliMpDEIterator.h"
#include "AliMpDEManager.h"
#include "AliMpSegmentation.h"

#include "AliLog.h"

#include <Riostream.h>
#include <float.h>

/// \cond CLASSIMP
ClassImp(AliMUONClusterReconstructor) // Class implementation in ROOT context
/// \endcond
 
//__________________________________________________________________________
AliMUONClusterReconstructor::AliMUONClusterReconstructor(AliMUONVClusterFinder* clusterFinder,
                                                         const AliMUONGeometryTransformer* transformer)
: TObject(),
  fClusterFinder(clusterFinder),
  fTransformer(transformer),
  fClusterStore(0x0),
  fNCluster(0)
{
    /// Standard Constructor
    /// Note that we adopt clusterFinder

  if (!transformer && clusterFinder)
  {
    AliFatal("I require a geometry transformer, otherwise I cannot compute "
             "global coordinates of the clusters !");    
  }
  
//  fRecModel->SetGhostChi2Cut(10);
}

//__________________________________________________________________________
AliMUONClusterReconstructor::~AliMUONClusterReconstructor(void)
{
  /// Destructor
  delete fClusterFinder;
}

//______________________________________________________________________________
void
AliMUONClusterReconstructor::ClusterizeOneDE(Int_t detElemId,
                                             const AliMUONVDigitStore& digitStore)
{
  /// Clusterize one detection element, which digits are in digitStore
  
//  AliDebug(1,Form("detElemId=%d, %d digits",detElemId,digitStore.GetSize()));
  
  if ( digitStore.IsEmpty() ) return;
  
  const AliMpVSegmentation* seg[2] = 
  { AliMpSegmentation::Instance()->GetMpSegmentation(detElemId,AliMp::kCath0),
    AliMpSegmentation::Instance()->GetMpSegmentation(detElemId,AliMp::kCath1)
  };
    
  if (!fClusterFinder->Prepare(seg,digitStore)) AliWarning(Form("No hit pad for DE %d ?",detElemId));
  
  AliMUONCluster* cluster;
  AliMUONVCluster *rawCluster;
  Int_t nPad;
  
  // Converts cluster objects into ones suitable for output
  while ( ( cluster = fClusterFinder->NextCluster() ) ) {
    
    // add new cluster to the store with information to build its ID
    // increment the number of clusters into the store
    rawCluster = fClusterStore->Add(AliMpDEManager::GetChamberId(detElemId), detElemId, fNCluster++);
    
    // fill array of Id of digits attached to this cluster
    nPad = cluster->Multiplicity();
    if (nPad < 1) AliWarning("no pad attached to the cluster");
    for (Int_t iPad=0; iPad<nPad; iPad++) {
      AliMUONPad *pad = cluster->Pad(iPad);
      rawCluster->AddDigitId(pad->GetUniqueID());
    }
    
    // fill charge and other cluster informations
    rawCluster->SetCharge(cluster->Charge());
    
    Double_t xg, yg, zg;
    fTransformer->Local2Global(detElemId, 
				cluster->Position().X(), cluster->Position().Y(), 
				0, xg, yg, zg);
    rawCluster->SetXYZ(xg, yg, zg);
    
    AliDebug(1,Form("Adding RawCluster detElemId %4d mult %2d charge %e (xl,yl,zl)=(%e,%e,%e) (xg,yg,zg)=(%e,%e,%e)",
		detElemId,nPad,cluster->Charge(),
		cluster->Position().X(),cluster->Position().Y(),0.0,
		xg,yg,zg));
  }
  
}

//____________________________________________________________________
void AliMUONClusterReconstructor::Digits2Clusters(const AliMUONVDigitStore& digitStore,
                                                  AliMUONVClusterStore& clusterStore)
{
  /// Clusterize the digitStore to produce a clusterStore
  
  fClusterStore = &clusterStore;
  fClusterStore->Clear();
  fNCluster = 0;
  
  AliMpDEIterator deIt;
  
  deIt.First();
  
  while (!deIt.IsDone())
  {
    AliMUONVDigitStore* deDigits = digitStore.Create();
    
    Int_t currentDE = deIt.CurrentDEId();
    AliMp::StationType stationType = AliMpDEManager::GetStationType(currentDE);
    if (stationType!=AliMp::kStationTrigger) 
    {
      TIter next(digitStore.CreateIterator(currentDE,currentDE));
      AliMUONVDigit* digit;

      while ( ( digit = static_cast<AliMUONVDigit*>(next()) ) )
      {
        if ( ! digit->Charge() > 0 ) continue; // skip void digits.
    
        deDigits->Add(*digit,AliMUONVDigitStore::kIgnore);
      }      
      ClusterizeOneDE(currentDE,*deDigits);
    }
    delete deDigits;
    deIt.Next();
  }
}
