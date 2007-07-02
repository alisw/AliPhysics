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

/// \class AliMUONClusterReconstructor
///
/// This class is just a steering class to loop over detection elements of tracking chambers,
/// extract the relevant digits, and pass them to the actual clusterizer class,
/// which operate on a single detection element at a time.
///
/// \author C. Finck and L. Aphecetche, Subatech
///

#include "AliMUONClusterReconstructor.h"

#include "AliLog.h"
#include "AliMUONCluster.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONRawCluster.h"
#include "AliMUONVClusterFinder.h"
#include "AliMUONVClusterStore.h"
#include "AliMUONVDigit.h"
#include "AliMUONVDigitStore.h"
#include "AliMpDEIterator.h"
#include "AliMpDEManager.h"
#include "AliMpSegmentation.h"
#include <Riostream.h>

/// \cond CLASSIMP
ClassImp(AliMUONClusterReconstructor) // Class implementation in ROOT context
/// \endcond
 
//__________________________________________________________________________
AliMUONClusterReconstructor::AliMUONClusterReconstructor(AliMUONVClusterFinder* clusterFinder,
                                                         const AliMUONGeometryTransformer* transformer)
: TObject(),
  fClusterFinder(clusterFinder),
  fTransformer(transformer),
  fClusterStore(0x0)
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
    
  Bool_t ok = fClusterFinder->Prepare(seg,digitStore);
  if ( !ok )
  {
    AliWarning(Form("No hit pad for DE %d ?",detElemId));
  }
  
  AliMUONCluster* cluster;
    
  while ( ( cluster = fClusterFinder->NextCluster() ) )
  {
    // Converts cluster objects into ones suitable for output
    //
    AliMUONRawCluster rawCluster;
    
    rawCluster.SetDetElemId(detElemId);
    
    for ( Int_t cathode = 0; cathode < 2; ++cathode )
    {
      rawCluster.SetMultiplicity(cathode,cluster->Multiplicity(cathode));
      rawCluster.SetCharge(cathode,cluster->Charge()); // both cathode get the total cluster charge
      Double_t xg, yg, zg;
      
      fTransformer->Local2Global(detElemId, 
                                 cluster->Position().X(), cluster->Position().Y(), 
                                 0, xg, yg, zg);
      
      if ( cathode == 0 )
      {
        AliDebug(1,Form("Adding RawCluster detElemId %4d mult %2d charge %e (xl,yl,zl)=(%e,%e,%e) (xg,yg,zg)=(%e,%e,%e)",
                        detElemId,cluster->Multiplicity(),cluster->Charge(),
                        cluster->Position().X(),cluster->Position().Y(),0.0,
                        xg,yg,zg));
      }
      rawCluster.SetX(cathode,xg);
      rawCluster.SetY(cathode,yg);
      rawCluster.SetZ(cathode,zg);      
    }
    fClusterStore->Add(rawCluster);
  }
}

//____________________________________________________________________
void AliMUONClusterReconstructor::Digits2Clusters(const AliMUONVDigitStore& digitStore,
                                                  AliMUONVClusterStore& clusterStore)
{
  /// Clusterize the digitStore to produce a clusterStore
  
  fClusterStore = &clusterStore;
  fClusterStore->Clear();
  
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
