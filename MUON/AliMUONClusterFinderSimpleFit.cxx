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

#include "AliMUONClusterFinderSimpleFit.h"

#include "AliLog.h"
#include "AliMpDEManager.h"
#include "AliMUONCluster.h"
#include "AliMUONConstants.h"
#include "AliMUONVDigit.h"
#include "AliMUONMathieson.h"
#include "AliMUONPad.h"
#include "AliMpArea.h"
#include "TClonesArray.h"
#include "TObjArray.h"
#include "TVector2.h"
#include "TVirtualFitter.h"
#include "TF1.h"
#include "AliMUONVDigitStore.h"
#include <Riostream.h>

//-----------------------------------------------------------------------------
/// \class AliMUONClusterFinderSimpleFit
///
/// Basic cluster finder 
/// 
/// We simply use AliMUONPreClusterFinder to get basic cluster,
/// and then we try to fit the charge repartition using a Mathieson
/// distribution, varying the position.
///
/// FIXME: this one is still at the developping stage...
///
/// \author Laurent Aphecetche
//-----------------------------------------------------------------------------

/// \cond CLASSIMP
ClassImp(AliMUONClusterFinderSimpleFit)
/// \endcond

namespace
{
  //___________________________________________________________________________
  void 
  FitFunction(Int_t& /*notused*/, Double_t* /*notused*/, 
              Double_t& f, Double_t* par, 
              Int_t /*notused*/)
  {
    /// Chi2 Function to minimize: Mathieson charge distribution in 2 dimensions
    
    TObjArray* userObjects = static_cast<TObjArray*>(TVirtualFitter::GetFitter()->GetObjectFit());
    
    AliMUONCluster* cluster = static_cast<AliMUONCluster*>(userObjects->At(0));
    AliMUONMathieson* mathieson = static_cast<AliMUONMathieson*>(userObjects->At(1));
    
    f = 0.0;
    Float_t qTot = cluster->Charge();
//    Float_t chargeCorrel[] = { cluster->Charge(0)/qTot, cluster->Charge(1)/qTot };
//    Float_t qRatio[] = { 1.0/par[2], par[2] };
    
    for ( Int_t i = 0 ; i < cluster->Multiplicity(); ++i )
    {
      AliMUONPad* pad = cluster->Pad(i);
      // skip pads w/ saturation or other problem(s)
      if ( pad->Status() ) continue; 
      TVector2 lowerLeft = TVector2(par[0],par[1]) - pad->Position() - pad->Dimensions();
      TVector2 upperRight(lowerLeft + pad->Dimensions()*2.0);
      Float_t estimatedCharge = mathieson->IntXY(lowerLeft.X(),lowerLeft.Y(),
                                                 upperRight.X(),upperRight.Y());
//      estimatedCharge *= 2/(1+qRatio[pad->Cathode()]);
      Float_t actualCharge = pad->Charge()/qTot;
      
      Float_t delta = (estimatedCharge - actualCharge);
      
      f += delta*delta;    
    }  
  }
}

//_____________________________________________________________________________
AliMUONClusterFinderSimpleFit::AliMUONClusterFinderSimpleFit(AliMUONVClusterFinder* clusterFinder)
: AliMUONVClusterFinder(),
fClusterFinder(clusterFinder),
fMathieson(0x0)
{
  /// ctor
}

//_____________________________________________________________________________
AliMUONClusterFinderSimpleFit::~AliMUONClusterFinderSimpleFit()
{
  /// dtor
  delete fClusterFinder;
  delete fMathieson;
}

//_____________________________________________________________________________
Bool_t 
AliMUONClusterFinderSimpleFit::Prepare(Int_t detElemId,
                                       TClonesArray* pads[2],
                                       const AliMpArea& area)
{
  /// Prepare for clustering

  // FIXME: should we get the Mathieson from elsewhere ?
  
  // Find out the DetElemId
  AliMq::Station12Type stationType = AliMpDEManager::GetStation12Type(detElemId);
  
  Float_t kx3 = AliMUONConstants::SqrtKx3();
  Float_t ky3 = AliMUONConstants::SqrtKy3();
  Float_t pitch = AliMUONConstants::Pitch();
  
  if ( stationType == AliMq::kStation1 )
  {
    kx3 = AliMUONConstants::SqrtKx3St1();
    ky3 = AliMUONConstants::SqrtKy3St1();
    pitch = AliMUONConstants::PitchSt1();
  }
  
  delete fMathieson;
  fMathieson = new AliMUONMathieson;
  
  fMathieson->SetPitch(pitch);
  fMathieson->SetSqrtKx3AndDeriveKx2Kx4(kx3);
  fMathieson->SetSqrtKy3AndDeriveKy2Ky4(ky3);

  return fClusterFinder->Prepare(detElemId,pads,area);
}

//_____________________________________________________________________________
AliMUONCluster* 
AliMUONClusterFinderSimpleFit::NextCluster()
{
  /// Returns next cluster
  
  if ( !fClusterFinder ) return 0x0;
  AliMUONCluster* cluster = fClusterFinder->NextCluster();
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
AliMUONClusterFinderSimpleFit::ComputePosition(AliMUONCluster& cluster)
{
  /// Compute the position of the given cluster, by fitting a Mathieson
  /// charge distribution to it
  
  TVirtualFitter* fitter = TVirtualFitter::Fitter(0,2);
  fitter->SetFCN(FitFunction);

  if ( cluster.Multiplicity() < 3 ) return;
  
  // We try a Mathieson fit, starting
  // with the center-of-gravity estimate as a first guess
  // for the cluster center.
  
  Double_t xCOG = cluster.Position().X();
  Double_t yCOG = cluster.Position().Y();
  
  Float_t stepX = 0.01; // cm
  Float_t stepY = 0.01; // cm
  
  Double_t arg(-1); // disable printout
  
  fitter->ExecuteCommand("SET PRINT",&arg,1);
  
  fitter->SetParameter(0,"cluster X position",xCOG,stepX,0,0);
  fitter->SetParameter(1,"cluster Y position",yCOG,stepY,0,0);
  
  TObjArray userObjects;
  
  userObjects.Add(&cluster);
  userObjects.Add(fMathieson);
  
  fitter->SetObjectFit(&userObjects);
  
  Int_t val = fitter->ExecuteCommand("MIGRAD",0,0);
  AliDebug(1,Form("ExecuteCommand returned value=%d",val));
  if ( val ) 
  {
    // fit failed. Using COG results, with big errors
    AliWarning("Fit failed. Using COG results for cluster=");
    StdoutToAliWarning(cluster.Print());
    cluster.SetPosition(TVector2(xCOG,yCOG),TVector2(TMath::Abs(xCOG),TMath::Abs(yCOG)));
    cluster.SetChi2(1E3);
  }
  
  Double_t results[] = { fitter->GetParameter(0),
    fitter->GetParameter(1) };
  
  Double_t errors[] = { fitter->GetParError(0),
    fitter->GetParError(1) };
  
  cluster.SetPosition(TVector2(results[0],results[1]),
                      TVector2(errors[0],errors[1]));
  
  Double_t amin, edm, errdef;
  Int_t nvpar, nparx;
  
  fitter->GetStats(amin, edm, errdef, nvpar, nparx);

  Double_t chi2 = amin;
  
  AliDebug(1,Form("Cluster fitted to (x,y)=(%e,%e) (xerr,yerr)=(%e,%e) \n chi2=%e ndf=%d",
                  results[0],results[1],
                  errors[0],errors[1],chi2,fitter->GetNumberFreeParameters()));
  cluster.SetChi2(chi2);
}



