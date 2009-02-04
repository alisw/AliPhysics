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
/// \class AliMUONClusterFinderPeakFit
/// 
/// Clusterizer class based on simple peak finder
///
/// Pre-clustering is handled by AliMUONPreClusterFinder
/// From a precluster a pixel array is built, and its local maxima are used
/// to get pads and make the fit with up to 3 hit candidates or compute pad 
/// centers of gravity for larger number of peaks.
///
/// \author Laurent Aphecetche (for the "new" C++ structure) and 
/// Alexander Zinchenko, JINR Dubna, for the hardcore of it ;-)
//-----------------------------------------------------------------------------

#include "AliMUONClusterFinderPeakFit.h"
#include "AliMUONCluster.h"
#include "AliMUONConstants.h"
#include "AliMUONPad.h"
#include "AliMUONMathieson.h"

#include "AliMpDEManager.h"
#include "AliMpPad.h"
#include "AliMpVSegmentation.h"

#include "AliLog.h"
#include "AliRunLoader.h"
//#include "AliCodeTimer.h"

#include <Riostream.h>
#include <TH2.h>
#include <TVirtualFitter.h>
#include <TMath.h>
//#include <TCanvas.h>

/// \cond CLASSIMP
ClassImp(AliMUONClusterFinderPeakFit)
/// \endcond
 
const Double_t AliMUONClusterFinderPeakFit::fgkZeroSuppression = 6; // average zero suppression value
//const Double_t AliMUONClusterFinderMLEM::fgkDistancePrecision = 1e-6; // (cm) used to check overlaps and so on
const Double_t AliMUONClusterFinderPeakFit::fgkDistancePrecision = 1e-3; // (cm) used to check overlaps and so on
const TVector2 AliMUONClusterFinderPeakFit::fgkIncreaseSize(-AliMUONClusterFinderPeakFit::fgkDistancePrecision,-AliMUONClusterFinderPeakFit::fgkDistancePrecision);
const TVector2 AliMUONClusterFinderPeakFit::fgkDecreaseSize(AliMUONClusterFinderPeakFit::fgkDistancePrecision,AliMUONClusterFinderPeakFit::fgkDistancePrecision);

// Status flags for pads
const Int_t AliMUONClusterFinderPeakFit::fgkZero = 0x0; ///< pad "basic" state
const Int_t AliMUONClusterFinderPeakFit::fgkMustKeep = 0x1; ///< do not kill (for pixels)
const Int_t AliMUONClusterFinderPeakFit::fgkUseForFit = 0x10; ///< should be used for fit
const Int_t AliMUONClusterFinderPeakFit::fgkOver = 0x100; ///< processing is over
const Int_t AliMUONClusterFinderPeakFit::fgkModified = 0x1000; ///< modified pad charge 
const Int_t AliMUONClusterFinderPeakFit::fgkCoupled = 0x10000; ///< coupled pad  

namespace
{
  //_____________________________________________________________________________
  Double_t Param2Coef(Int_t icand, Double_t coef, Double_t *par, Int_t nHits)
  {
    /// Extract hit contribution scale factor from fit parameters
  
    //Int_t nHits = TMath::Nint(par[8]);
    if (nHits == 1) return 1.;
    if (nHits == 2) return icand==0 ? par[2] : TMath::Max(1.-par[2],0.);
    if (icand == 0) return par[2];
    if (icand == 1) return TMath::Max((1.-par[2])*par[5], 0.);
    return TMath::Max(1.-par[2]-coef,0.);
  }

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
    AliMUONClusterFinderPeakFit* finder = 
      static_cast<AliMUONClusterFinderPeakFit*>(userObjects->At(2));
    
    f = 0.0;
    Int_t nHits = finder->GetNMax(), npads = cluster->Multiplicity();
    Double_t qTot = cluster->Charge(), coef = 0;
    //if (cluster->Multiplicity(0) == 0 || cluster->Multiplicity(1) == 0) qTot *= 2.;
    
    for ( Int_t i = 0 ; i < npads; ++i )
    {
      AliMUONPad* pad = cluster->Pad(i);
      // skip pads w/ saturation or other problem(s)
      //if ( pad->Status() ) continue; 
      if ( pad->IsSaturated() ) continue; 
      Double_t charge = 0.;
      for (Int_t j = 0; j < nHits; ++j) {
	// Sum over hits
	Int_t indx = 3 * j;
	TVector2 lowerLeft = TVector2(par[indx],par[indx+1]) - pad->Position() - pad->Dimensions();
	TVector2 upperRight(lowerLeft + pad->Dimensions()*2.0);
	Double_t estimatedCharge = mathieson->IntXY(lowerLeft.X(),lowerLeft.Y(),
						    upperRight.X(),upperRight.Y());
	coef = Param2Coef(j, coef, par, nHits);
	charge += estimatedCharge * coef;
      }
      charge *= qTot;
      
      Double_t delta = charge - pad->Charge();
      delta *= delta;
      delta /= pad->Charge();    
      f += delta;
    }
    f /= (qTot/npads);
    //cout << qTot << " " << par[0] << " " << par[1] << " " << f << endl;
  }
}

//_____________________________________________________________________________
AliMUONClusterFinderPeakFit::AliMUONClusterFinderPeakFit(Bool_t plot, AliMUONVClusterFinder* clusterFinder)
  : AliMUONVClusterFinder(),
fPreClusterFinder(clusterFinder),
fPreCluster(0x0),
fClusterList(),
fMathieson(0x0),
fEventNumber(0),
fDetElemId(-1),
fClusterNumber(0),
fNMax(0),
fHistAnode(0x0),
fPixArray(new TObjArray(20)),
fDebug(0),
fPlot(plot),
fNClusters(0),
fNAddVirtualPads(0)
{
  /// Constructor
 
  fSegmentation[1] = fSegmentation[0] = 0x0; 

  if (fPlot) fDebug = 1;
}

//_____________________________________________________________________________
AliMUONClusterFinderPeakFit::~AliMUONClusterFinderPeakFit()
{
/// Destructor
  delete fPixArray; fPixArray = 0;
  delete fPreClusterFinder;
  AliInfo(Form("Total clusters %d AddVirtualPad needed %d",
               fNClusters,fNAddVirtualPads));
  delete fMathieson;
}

//_____________________________________________________________________________
Bool_t 
AliMUONClusterFinderPeakFit::Prepare(Int_t detElemId, TClonesArray* pads[2],
				     const AliMpArea& area, const AliMpVSegmentation* seg[2])
{
  /// Prepare for clustering
//  AliCodeTimerAuto("")
  
  for ( Int_t i = 0; i < 2; ++i )
  {
    fSegmentation[i] = seg[i];
  }
  
  // Find out the DetElemId
  fDetElemId = detElemId;
  
  // find out current event number, and reset the cluster number
  AliRunLoader *runLoader = AliRunLoader::Instance();
  fEventNumber = runLoader ? runLoader->GetEventNumber() : 0;
  fClusterNumber = -1;
  fClusterList.Delete();
  
  AliDebug(3,Form("EVT %d DE %d",fEventNumber,fDetElemId));
  
  AliMq::Station12Type stationType = AliMpDEManager::GetStation12Type(fDetElemId);
  
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

  if ( fPreClusterFinder->NeedSegmentation() )
  {
    return fPreClusterFinder->Prepare(detElemId,pads,area,seg);
  }
  else
  {
    return fPreClusterFinder->Prepare(detElemId,pads,area);
  }
}

//_____________________________________________________________________________
AliMUONCluster* 
AliMUONClusterFinderPeakFit::NextCluster()
{
  /// Return next cluster
//  AliCodeTimerAuto("")
  
  // if the list of clusters is not void, pick one from there
  TObject* o = fClusterList.At(++fClusterNumber);
  if ( o != 0x0 ) return static_cast<AliMUONCluster*>(o);
  
  //FIXME : at this point, must check whether we've used all the digits
  //from precluster : if not, let the preclustering know about those unused
  //digits, so it can reuse them
  
  // if the cluster list is exhausted, we need to go to the next
  // pre-cluster and treat it

  fPreCluster = fPreClusterFinder->NextCluster();
  
  if (!fPreCluster)
  {
    // we are done
    return 0x0;
  }
    
  fClusterList.Delete(); // reset the list of clusters for this pre-cluster
  fClusterNumber = -1; 
  
  WorkOnPreCluster();

  // WorkOnPreCluster may have used only part of the pads, so we check that
  // now, and let the unused pads be reused by the preclustering...
  
  Int_t mult = fPreCluster->Multiplicity();
  for ( Int_t i = 0; i < mult; ++i )
  {
    AliMUONPad* pad = fPreCluster->Pad(i);
    if ( !pad->IsUsed() )
    {
      fPreClusterFinder->UsePad(*pad);
    }
  }
  
  return NextCluster();
}

//_____________________________________________________________________________
Bool_t
AliMUONClusterFinderPeakFit::WorkOnPreCluster()
{
  /// Starting from a precluster, builds a pixel array, and then
  /// extract clusters from this array
  
  //  AliCodeTimerAuto("")	

  if (fDebug) {
    cout << " *** Event # " << fEventNumber
	 << " det. elem.: " << fDetElemId << endl;
    for (Int_t j = 0; j < fPreCluster->Multiplicity(); ++j) {
      AliMUONPad* pad = fPreCluster->Pad(j);
      printf(" bbb %3d %1d %8.4f %8.4f %8.4f %8.4f %6.1f %3d %3d %2d %1d %1d \n",
	     j, pad->Cathode(), pad->Coord(0), pad->Coord(1), pad->DX()*2, pad->DY()*2,
             pad->Charge(), pad->Ix(), pad->Iy(), pad->Status(), pad->IsReal(), pad->IsSaturated());
    }
  }

  AliMUONCluster* cluster = CheckPrecluster(*fPreCluster);
  if (!cluster) return kFALSE;

  BuildPixArray(*cluster);
  
  if ( fPixArray->GetLast() < 0 )
  {
    AliDebug(1,"No pixel for the above cluster");
    delete cluster;
    return kFALSE;
  }
  
  Int_t nMax = 1, localMax[100], maxPos[100] = {0};
  Double_t maxVal[100];
  
  nMax = FindLocalMaxima(fPixArray, localMax, maxVal); // find local maxima

  if (nMax > 1) TMath::Sort(nMax, maxVal, maxPos, kTRUE); // in descending order
  
  if (nMax <= 3) {
    FindClusterFit(*cluster, localMax, maxPos, nMax);
  } else {
    for (Int_t i = 0; i < nMax; ++i) 
      {
	FindClusterCOG(*cluster, localMax, maxPos[i]);
      }
  }

  delete cluster;
  if (fPlot == 0) {
    delete fHistAnode;
    fHistAnode = 0x0;
  }
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t 
AliMUONClusterFinderPeakFit::Overlap(const AliMUONPad& pad, const AliMUONPad& pixel)
{
  /// Check if the pad and the pixel overlaps

  // make a fake pad from the pixel
  AliMUONPad tmp(pad.DetElemId(),pad.Cathode(),pad.Ix(),pad.Iy(),
                 pixel.Coord(0),pixel.Coord(1),
                 pixel.Size(0),pixel.Size(1),0);
  
  return AliMUONPad::AreOverlapping(pad,tmp,fgkDecreaseSize);
}

//_____________________________________________________________________________
AliMUONCluster* 
AliMUONClusterFinderPeakFit::CheckPrecluster(const AliMUONCluster& origCluster)
{
  /// Check precluster in order to attempt to simplify it (mostly for
  /// two-cathode preclusters)
    
  //  AliCodeTimerAuto("")

  // Disregard small clusters (leftovers from splitting or noise)
  if ((origCluster.Multiplicity()==1 || origCluster.Multiplicity()==2) &&
      origCluster.Charge(0)+origCluster.Charge(1) < 10) 
  { 
    return 0x0;
  }

  AliMUONCluster* cluster = new AliMUONCluster(origCluster);

  AliDebug(2,"Start of CheckPreCluster=");
  //StdoutToAliDebug(2,cluster->Print("full"));

  AliMUONCluster* rv(0x0);
  
  if (cluster->Multiplicity(0) && cluster->Multiplicity(1))
  { 
    rv = CheckPreclusterTwoCathodes(cluster);
  }
  else
  {
    rv = cluster;
  }
  return rv;
}

//_____________________________________________________________________________
AliMUONCluster*
AliMUONClusterFinderPeakFit::CheckPreclusterTwoCathodes(AliMUONCluster* cluster)
{
  /// Check two-cathode cluster
  
  Int_t npad = cluster->Multiplicity();
  Int_t* flags = new Int_t[npad];
  for (Int_t j = 0; j < npad; ++j) flags[j] = 0;
  
  // Check pad overlaps
  for ( Int_t i = 0; i < npad; ++i) 
  {
    AliMUONPad* padi = cluster->Pad(i);
    if ( padi->Cathode() != 0 ) continue;
    for (Int_t j = i+1; j < npad; ++j) 
    {
      AliMUONPad* padj = cluster->Pad(j);
      if ( padj->Cathode() != 1 ) continue;
      if ( !AliMUONPad::AreOverlapping(*padi,*padj,fgkDecreaseSize) ) continue;
      flags[i] = flags[j] = 1; // mark overlapped pads
    } 
  } 
  
  // Check if all pads overlap
  Int_t nFlags=0;
  for (Int_t i = 0; i < npad; ++i) 
  {
    if (!flags[i]) ++nFlags;
  }
  
  if (nFlags > 0) 
  {
    // not all pads overlap.
    if (fDebug) cout << " nFlags: " << nFlags << endl;
    TObjArray toBeRemoved;
    for (Int_t i = 0; i < npad; ++i) 
    {
      AliMUONPad* pad = cluster->Pad(i);
      if (flags[i]) continue;
      Int_t cath = pad->Cathode();
      Int_t cath1 = TMath::Even(cath);
      // Check for edge effect (missing pads on the _other_ cathode)
      AliMpPad mpPad = fSegmentation[cath1]->PadByPosition(pad->Position(),kFALSE);
      if (!mpPad.IsValid()) continue;
      //if (nFlags == 1 && pad->Charge() < fgkZeroSuppression * 3) continue;
      if (nFlags == 1 && pad->Charge() < 20) continue;
      AliDebug(2,Form("Releasing the following pad : de,cath,ix,iy %d,%d,%d,%d charge %e",
                      fDetElemId,pad->Cathode(),pad->Ix(),pad->Iy(),pad->Charge()));
      toBeRemoved.AddLast(pad);
      fPreCluster->Pad(i)->Release();
    }
    Int_t nRemove = toBeRemoved.GetEntriesFast();
    for ( Int_t i = 0; i < nRemove; ++i )
    {
      cluster->RemovePad(static_cast<AliMUONPad*>(toBeRemoved.UncheckedAt(i)));
    }
  } 
  
  // Check correlations of cathode charges
  if ( !cluster->IsSaturated() && cluster->ChargeAsymmetry() > 1 )
  {
    // big difference
    Int_t cathode = cluster->MaxRawChargeCathode();
    Int_t imin(-1);
    Int_t imax(-1);
    Double_t cmax(0);
    Double_t cmin(1E9);
    
    // get min and max pad charges on the cathode opposite to the 
    // max pad (given by MaxRawChargeCathode())
    //
    Int_t mult = cluster->Multiplicity();
    for ( Int_t i = 0; i < mult; ++i )
    {
      AliMUONPad* pad = cluster->Pad(i);
      if ( pad->Cathode() != cathode || !pad->IsReal() )
      {
        // only consider pads in the opposite cathode, and
        // only consider real pads (i.e. exclude the virtual ones)
        continue;
      }
      if ( pad->Charge() < cmin )
      {
        cmin = pad->Charge();
        imin = i;
	if (imax < 0) {
	  imax = imin;
	  cmax = cmin;
	}
      }
      else if ( pad->Charge() > cmax )
      {
        cmax = pad->Charge();
        imax = i;
      }      
    }
    AliDebug(2,Form("Pad imin,imax %d,%d cmin,cmax %e,%e",
                    imin,imax,cmin,cmax));
    //
    // arrange pads according to their distance to the max, normalized
    // to the pad size
    Double_t* dist = new Double_t[mult];
    Double_t dxMin(1E9);
    Double_t dyMin(1E9);
    Double_t dmin(0);
    
    AliMUONPad* padmax = cluster->Pad(imax);
    
    for ( Int_t i = 0; i < mult; ++i )
    {
      dist[i] = 0.0;
      if ( i == imax) continue;
      AliMUONPad* pad = cluster->Pad(i);
      if ( pad->Cathode() != cathode || !pad->IsReal() ) continue;
      Double_t dx = (pad->X()-padmax->X())/padmax->DX()/2.0;
      Double_t dy = (pad->Y()-padmax->Y())/padmax->DY()/2.0;
      dist[i] = TMath::Sqrt(dx*dx+dy*dy);
      if ( i == imin )
      {
        dmin = dist[i] + 1E-3; // distance to the pad with minimum charge
        dxMin = dx;
        dyMin = dy;
      }      
    }
    
    TMath::Sort(mult,dist,flags,kFALSE); // in ascending order
    Double_t xmax(-1), distPrev(999);
    TObjArray toBeRemoved;
    
    for ( Int_t i = 0; i < mult; ++i )
    {
      Int_t indx = flags[i];
      AliMUONPad* pad = cluster->Pad(indx);
      if ( pad->Cathode() != cathode || !pad->IsReal() ) continue;
      if ( dist[indx] > dmin )
      {
        // farther than the minimum pad
        Double_t dx = (pad->X()-padmax->X())/padmax->DX()/2.0;
        Double_t dy = (pad->Y()-padmax->Y())/padmax->DY()/2.0;
        dx *= dxMin;
        dy *= dyMin;
        if (dx >= 0 && dy >= 0) continue;
        if (TMath::Abs(dx) > TMath::Abs(dy) && dx >= 0) continue;
        if (TMath::Abs(dy) > TMath::Abs(dx) && dy >= 0) continue;        
      }
      if (dist[indx] > distPrev + 1) break; // overstepping empty pads
      if ( pad->Charge() <= cmax || TMath::Abs(dist[indx]-xmax) < 1E-3 )
      {
        // release pad
        if (TMath::Abs(dist[indx]-xmax) < 1.e-3) 
        {
          cmax = TMath::Max(pad->Charge(),cmax);
        }
        else
        {
          cmax = pad->Charge();
        }
        xmax = dist[indx];
	distPrev = dist[indx];
        AliDebug(2,Form("Releasing the following pad : de,cath,ix,iy %d,%d,%d,%d charge %e",
                        fDetElemId,pad->Cathode(),pad->Ix(),pad->Iy(),
                        pad->Charge()));
  
        toBeRemoved.AddLast(pad);
        fPreCluster->Pad(indx)->Release();
      }
    }
    Int_t nRemove = toBeRemoved.GetEntriesFast();
    for ( Int_t i = 0; i < nRemove; ++i )
    {
      cluster->RemovePad(static_cast<AliMUONPad*>(toBeRemoved.UncheckedAt(i)));
    }    
    delete[] dist;
  } // if ( !cluster->IsSaturated() && 
  
  delete[] flags;
  
  AliDebug(2,"End of CheckPreClusterTwoCathodes=");
  //StdoutToAliDebug(2,cluster->Print("full"));

  return cluster;    
}

//_____________________________________________________________________________
void
AliMUONClusterFinderPeakFit::CheckOverlaps()
{
  /// For debug only : check if some pixels overlap...
  
  Int_t nPix = fPixArray->GetLast()+1;
  Int_t dummy(0);
  
  for ( Int_t i = 0; i < nPix; ++i )
  {
    AliMUONPad* pixelI = Pixel(i);
    AliMUONPad pi(dummy,dummy,dummy,dummy,
                  pixelI->Coord(0),pixelI->Coord(1),
                  pixelI->Size(0),pixelI->Size(1),0.0);
    
    for ( Int_t j = i+1; j < nPix; ++j )
    {
      AliMUONPad* pixelJ = Pixel(j);
      AliMUONPad pj(dummy,dummy,dummy,dummy,
                    pixelJ->Coord(0),pixelJ->Coord(1),
                    pixelJ->Size(0),pixelJ->Size(1),0.0);  
      AliMpArea area;
      
      if ( AliMUONPad::AreOverlapping(pi,pj,fgkDecreaseSize,area) )
      {
        AliInfo(Form("The following 2 pixels (%d and %d) overlap !",i,j));
	/*
        StdoutToAliInfo(pixelI->Print();
                        cout << " Surface = " << pixelI->Size(0)*pixelI->Size(1)*4 << endl;
                        pixelJ->Print();
                        cout << " Surface = " << pixelJ->Size(0)*pixelJ->Size(1)*4 << endl;
                        cout << " Area surface = " << area.Dimensions().X()*area.Dimensions().Y()*4 << endl;
                        cout << "-------" << endl;
                        );
	*/        
      }
    }    
  }
}

//_____________________________________________________________________________
void AliMUONClusterFinderPeakFit::BuildPixArray(AliMUONCluster& cluster)
{
  /// Build pixel array 
  
  Int_t npad = cluster.Multiplicity();
  if (npad<=0) 
  {
    AliWarning("Got no pad at all ?!");
  }
  
  fPixArray->Delete();
  BuildPixArrayOneCathode(cluster);
  
//  StdoutToAliDebug(2,cout << "End of BuildPixelArray:" << endl;
//                   fPixArray->Print(););
  //CheckOverlaps();//FIXME : this is for debug only. Remove it.
}

//_____________________________________________________________________________
void AliMUONClusterFinderPeakFit::BuildPixArrayOneCathode(AliMUONCluster& cluster)
{
  /// Build the pixel array

//  AliDebug(2,Form("cluster.Multiplicity=%d",cluster.Multiplicity()));

  TVector2 dim = cluster.MinPadDimensions (-1, kFALSE);
  Double_t width[2] = {dim.X(), dim.Y()}, xy0[2];
  Int_t found[2] = {0}, mult = cluster.Multiplicity();

  for ( Int_t i = 0; i < mult; ++i) {
    AliMUONPad* pad = cluster.Pad(i);
    for (Int_t j = 0; j < 2; ++j) {
      if (found[j] == 0 && TMath::Abs(pad->Size(j)-width[j]) < fgkDistancePrecision) { 
	xy0[j] = pad->Coord(j);
	found[j] = 1;
      }
    }
    if (found[0] && found[1]) break;
  }

  Double_t min[2], max[2];
  Int_t cath0 = 0, cath1 = 1;
  if (cluster.Multiplicity(0) == 0) cath0 = 1;
  else if (cluster.Multiplicity(1) == 0) cath1 = 0;

  TVector2 leftDown = cluster.Area(cath0).LeftDownCorner();
  TVector2 rightUp = cluster.Area(cath0).RightUpCorner();
  min[0] = leftDown.X();
  min[1] = leftDown.Y();
  max[0] = rightUp.X();
  max[1] = rightUp.Y();
  if (cath1 != cath0) {
    leftDown = cluster.Area(cath1).LeftDownCorner();
    rightUp = cluster.Area(cath1).RightUpCorner();
    min[0] = TMath::Max (min[0], leftDown.X());
    min[1] = TMath::Max (min[1], leftDown.Y());
    max[0] = TMath::Min (max[0], rightUp.X());
    max[1] = TMath::Min (max[1], rightUp.Y());
  }

  // Adjust limits
  if (cath0 != cath1) {
    TVector2 dim0 = cluster.MinPadDimensions (0, -1, kFALSE);
    TVector2 dim1 = cluster.MinPadDimensions (1, -1, kFALSE);
    if (TMath::Abs(dim0.Y()-dim1.Y()) < fgkDistancePrecision) {
      // The same size of pads on both cathodes - check position
      AliMUONPad* pad0 = cluster.Pad(0);
      for ( Int_t i = 1; i < mult; ++i) {
	AliMUONPad* pad = cluster.Pad(i);
	if (pad->Cathode() == pad0->Cathode()) continue;
	Double_t dist = TMath::Abs (pad0->Coord(1) - pad->Coord(1));
	Double_t dd = dist - Int_t(dist/width[1]/2.) * width[1] * 2.;
	if (TMath::Abs(dd/width[1]/2.-0.5) < fgkDistancePrecision) { 
	  // Half pad shift between cathodes
	  width[0] /= 2.; 
	  width[1] /= 2.;
	}
	break;
      }
    }
  }

  Int_t nbins[2];
  for (Int_t i = 0; i < 2; ++i) {
    Double_t dist = (min[i] - xy0[i]) / width[i] / 2;
    if (TMath::Abs(dist) < 1.e-6) dist = -1.e-6;
    min[i] = xy0[i] + (TMath::Nint(dist-TMath::Sign(1.e-6,dist)) 
		       + TMath::Sign(0.5,dist)) * width[i] * 2;
    nbins[i] = TMath::Nint ((max[i] - min[i]) / width[i] / 2);
    if (nbins[i] == 0) ++nbins[i];
    max[i] = min[i] + nbins[i] * width[i] * 2;
    //cout << dist << " " << min[i] << " " << max[i] << " " << nbins[i] << endl;
  }

  // Book histogram
  TH2D *hist1 = new TH2D ("Grid", "", nbins[0], min[0], max[0], nbins[1], min[1], max[1]);
  TH2D *hist2 = new TH2D ("Entries", "", nbins[0], min[0], max[0], nbins[1], min[1], max[1]);
  TAxis *xaxis = hist1->GetXaxis();
  TAxis *yaxis = hist1->GetYaxis();

  // Fill histogram
  for ( Int_t i = 0; i < mult; ++i) {
    AliMUONPad* pad = cluster.Pad(i);
    Int_t ix0 = xaxis->FindBin(pad->X());
    Int_t iy0 = yaxis->FindBin(pad->Y());
    PadOverHist(0, ix0, iy0, pad, hist1, hist2);
  }

  // Store pixels
  for (Int_t i = 1; i <= nbins[0]; ++i) {
    Double_t x = xaxis->GetBinCenter(i);
    for (Int_t j = 1; j <= nbins[1]; ++j) {
      if (hist2->GetCellContent(i,j) < 0.1) continue;
      if (cath0 != cath1) {
	// Two-sided cluster
	Double_t cont = hist2->GetCellContent(i,j);
	if (cont < 999.) continue;
	if (cont-Int_t(cont/1000.)*1000. < 0.5) continue;
      }
      //if (hist2->GetCellContent(i,j) < 1.1 && cluster.Multiplicity(0) && 
      //  cluster.Multiplicity(1)) continue;
      Double_t y = yaxis->GetBinCenter(j);
      Double_t charge = hist1->GetCellContent(i,j);
      AliMUONPad* pixPtr = new AliMUONPad(x, y, width[0], width[1], charge);
      fPixArray->Add(pixPtr);
    }  
  }
  /*
  if (fPixArray->GetEntriesFast() == 1) {
    // Split pixel into 2
    AliMUONPad* pixPtr = static_cast<AliMUONPad*> (fPixArray->UncheckedAt(0));
    pixPtr->SetSize(0,width[0]/2.);
    pixPtr->Shift(0,-width[0]/4.);
    pixPtr = new AliMUONPad(pixPtr->X()+width[0], pixPtr->Y(), width[0]/2., width[1], pixPtr->Charge());
    fPixArray->Add(pixPtr);
  }
  */
  //fPixArray->Print();
  delete hist1;
  delete hist2;
}

//_____________________________________________________________________________
void AliMUONClusterFinderPeakFit::PadOverHist(Int_t idir, Int_t ix0, Int_t iy0, AliMUONPad *pad,
					      TH2D *hist1, TH2D *hist2)
{
  /// "Span" pad over histogram in the direction idir

  TAxis *axis = idir == 0 ? hist1->GetXaxis() : hist1->GetYaxis();
  Int_t nbins = axis->GetNbins(), cath = pad->Cathode();
  Double_t bin = axis->GetBinWidth(1), amask = TMath::Power(1000.,cath*1.);

  Int_t nbinPad = (Int_t)(pad->Size(idir)/bin*2+fgkDistancePrecision) + 1; // number of bins covered by pad

  for (Int_t i = 0; i < nbinPad; ++i) {
    Int_t ixy = idir == 0 ? ix0 + i : iy0 + i;
    if (ixy > nbins) break;
    Double_t lowEdge = axis->GetBinLowEdge(ixy);
    if (lowEdge + fgkDistancePrecision > pad->Coord(idir) + pad->Size(idir)) break;
    if (idir == 0) PadOverHist(1, ixy, iy0, pad, hist1, hist2); // span in the other direction
    else {
      // Fill histogram
      Double_t cont = pad->Charge();
      if (hist2->GetCellContent(ix0, ixy) > 0.1) 
	cont = TMath::Min (hist1->GetCellContent(ix0, ixy), cont) 
	     + TMath::Min (TMath::Max(hist1->GetCellContent(ix0, ixy),cont)*0.1, 10.);
      hist1->SetCellContent(ix0, ixy, cont);
      hist2->SetCellContent(ix0, ixy, hist2->GetCellContent(ix0, ixy)+amask);
    }
  }

  for (Int_t i = -1; i > -nbinPad; --i) {
    Int_t ixy = idir == 0 ? ix0 + i : iy0 + i;
    if (ixy < 1) break;
    Double_t upEdge = axis->GetBinUpEdge(ixy);
    if (upEdge - fgkDistancePrecision < pad->Coord(idir) - pad->Size(idir)) break;
    if (idir == 0) PadOverHist(1, ixy, iy0, pad, hist1, hist2); // span in the other direction
    else {
      // Fill histogram
      Double_t cont = pad->Charge();
      if (hist2->GetCellContent(ix0, ixy) > 0.1) 
	cont = TMath::Min (hist1->GetCellContent(ix0, ixy), cont)
	     + TMath::Min (TMath::Max(hist1->GetCellContent(ix0, ixy),cont)*0.1,10.);
      hist1->SetCellContent(ix0, ixy, cont);
      hist2->SetCellContent(ix0, ixy, hist2->GetCellContent(ix0, ixy)+amask);
    }
  }
}

//_____________________________________________________________________________
Int_t AliMUONClusterFinderPeakFit::FindLocalMaxima(TObjArray *pixArray, Int_t *localMax, Double_t *maxVal)
{
/// Find local maxima in pixel space 

  AliDebug(1,Form("nPix=%d",pixArray->GetLast()+1));

  //TH2D *hist = NULL;
  //delete ((TH2D*) gROOT->FindObject("anode"));
  //if (pixArray == fPixArray) hist = (TH2D*) gROOT->FindObject("anode");
  //else { hist = (TH2D*) gROOT->FindObject("anode1"); cout << hist << endl; }
  //if (hist) hist->Delete();
  delete fHistAnode;
 
  Double_t xylim[4] = {999, 999, 999, 999};

  Int_t nPix = pixArray->GetEntriesFast();
  AliMUONPad *pixPtr = 0;
  for (Int_t ipix = 0; ipix < nPix; ++ipix) {
    pixPtr = (AliMUONPad*) pixArray->UncheckedAt(ipix);
    for (Int_t i = 0; i < 4; ++i) 
         xylim[i] = TMath::Min (xylim[i], (i%2 ? -1 : 1)*pixPtr->Coord(i/2));
  }
  for (Int_t i = 0; i < 4; ++i) xylim[i] -= pixPtr->Size(i/2); 

  Int_t nx = TMath::Nint ((-xylim[1]-xylim[0])/pixPtr->Size(0)/2);
  Int_t ny = TMath::Nint ((-xylim[3]-xylim[2])/pixPtr->Size(1)/2);
  if (pixArray == fPixArray) fHistAnode = new TH2D("anode","anode",nx,xylim[0],-xylim[1],ny,xylim[2],-xylim[3]);
  else fHistAnode = new TH2D("anode1","anode1",nx,xylim[0],-xylim[1],ny,xylim[2],-xylim[3]);
  for (Int_t ipix = 0; ipix < nPix; ++ipix) {
    pixPtr = (AliMUONPad*) pixArray->UncheckedAt(ipix);
    fHistAnode->Fill(pixPtr->Coord(0), pixPtr->Coord(1), pixPtr->Charge());
  }
//  if (fDraw && pixArray == fPixArray) fDraw->DrawHist("c2", hist);

  Int_t nMax = 0, indx, nxy = ny * nx;
  Int_t *isLocalMax = new Int_t[nxy];
  for (Int_t i = 0; i < nxy; ++i) isLocalMax[i] = 0; 

  for (Int_t i = 1; i <= ny; ++i) {
    indx = (i-1) * nx;
    for (Int_t j = 1; j <= nx; ++j) {
      if (fHistAnode->GetCellContent(j,i) < 0.5) continue;
      //if (isLocalMax[indx+j-1] < 0) continue;
      if (isLocalMax[indx+j-1] != 0) continue;
      FlagLocalMax(fHistAnode, i, j, isLocalMax);
    }
  }

  for (Int_t i = 1; i <= ny; ++i) {
    indx = (i-1) * nx;
    for (Int_t j = 1; j <= nx; ++j) {
      if (isLocalMax[indx+j-1] > 0) { 
	localMax[nMax] = indx + j - 1; 
	maxVal[nMax++] = fHistAnode->GetCellContent(j,i);
	if (nMax > 99) break;
      }
    }
    if (nMax > 99) {
      AliError(" Too many local maxima !!!");
      break;
    }
  }
  if (fDebug) cout << " Local max: " << nMax << endl;
  delete [] isLocalMax; 
  return nMax;
}

//_____________________________________________________________________________
void AliMUONClusterFinderPeakFit::FlagLocalMax(TH2D *hist, Int_t i, Int_t j, Int_t *isLocalMax)
{
/// Flag pixels (whether or not local maxima)

  Int_t nx = hist->GetNbinsX();
  Int_t ny = hist->GetNbinsY();
  Int_t cont = TMath::Nint (hist->GetCellContent(j,i));
  Int_t cont1 = 0, indx = (i-1)*nx+j-1, indx1 = 0, indx2 = 0;

  Int_t ie = i + 2, je = j + 2;
  for (Int_t i1 = i-1; i1 < ie; ++i1) {
    if (i1 < 1 || i1 > ny) continue;
    indx1 = (i1 - 1) * nx;
    for (Int_t j1 = j-1; j1 < je; ++j1) {
      if (j1 < 1 || j1 > nx) continue;
      if (i == i1 && j == j1) continue;
      indx2 = indx1 + j1 - 1;
      cont1 = TMath::Nint (hist->GetCellContent(j1,i1));
      if (cont < cont1) { isLocalMax[indx] = -1; return; }
      else if (cont > cont1) isLocalMax[indx2] = -1;
      else { // the same charge
	isLocalMax[indx] = 1; 
	if (isLocalMax[indx2] == 0) {
	  FlagLocalMax(hist, i1, j1, isLocalMax);
	  if (isLocalMax[indx2] < 0) { isLocalMax[indx] = -1; return; }
	  else isLocalMax[indx2] = -1;
	}
      } 
    }
  }
  isLocalMax[indx] = 1; // local maximum
}

//_____________________________________________________________________________
void AliMUONClusterFinderPeakFit::FindClusterFit(AliMUONCluster& cluster, Int_t *localMax, 
						 Int_t *maxPos, Int_t nMax)
{
/// Fit pad charge distribution with nMax hit hypothesis

  //if (cluster.Multiplicity(0) == 0 || cluster.Multiplicity(1) == 0) cout << cluster.Multiplicity(0) << " " << cluster.Multiplicity(1) << " " << cluster.Charge() << " " << cluster.Charge(0) << " " << cluster.Charge(1) << " " << endl;

  //TH2D *hist = (TH2D*) gROOT->FindObject("anode");
  Int_t nx = fHistAnode->GetNbinsX();
  //Int_t ny = hist->GetNbinsY();
  Double_t xmin = fHistAnode->GetXaxis()->GetXmin(); //- hist->GetXaxis()->GetBinWidth(1);
  Double_t xmax = fHistAnode->GetXaxis()->GetXmax(); //+ hist->GetXaxis()->GetBinWidth(1);
  Double_t ymin = fHistAnode->GetYaxis()->GetXmin(); //- hist->GetYaxis()->GetBinWidth(1);
  Double_t ymax = fHistAnode->GetYaxis()->GetXmax(); //+ hist->GetYaxis()->GetBinWidth(1);

  TVirtualFitter* fitter = TVirtualFitter::Fitter(0,nMax*3);
  fitter->Clear("");
  fitter->SetFCN(FitFunction);

  Float_t stepX = 0.01; // cm
  Float_t stepY = 0.01; // cm
  Float_t stepQ = 0.01; // 
  
  Double_t args[10] = {-1.}; // disable printout
  
  fitter->ExecuteCommand("SET PRINT",args,1);
  fitter->ExecuteCommand("SET NOW",args,0); // no warnings
  
  Int_t indx = 0;
  fNMax = nMax;
  for (Int_t i = 0; i < nMax; ++i) {
    Int_t ic = localMax[maxPos[i]] / nx + 1;
    Int_t jc = localMax[maxPos[i]] % nx + 1;
    Double_t yc = fHistAnode->GetYaxis()->GetBinCenter(ic);
    Double_t xc = fHistAnode->GetXaxis()->GetBinCenter(jc);
    indx = 3 * i;
    fitter->SetParameter(indx,"Hit X position",xc,stepX,xmin,xmax);
    fitter->SetParameter(indx+1,"Hit Y position",yc,stepY,ymin,ymax);
    fitter->SetParameter(indx+2,"Hit contribution",0.6,stepQ,0.,1.);
  }
  fitter->SetParameter(indx+2,"Hit contribution",0.,0.,0,0);
  //fitter->SetParameter(8,"Number of hits",nMax,0.,0,0);

  TObjArray userObjects;
  
  userObjects.Add(&cluster);
  userObjects.Add(fMathieson);
  userObjects.Add(this);
  
  fitter->SetObjectFit(&userObjects);
  
  args[0] = 500.;
  args[1] = 1.;
  /*Int_t stat =*/ fitter->ExecuteCommand("MIGRAD",args,2);
  //if (stat) { cout << " stat = " << stat << " " << fDetElemId << endl; /*exit(0);*/ }
  //Int_t nvpar, nparx;
  //Double_t amin, edm, errdef; 
  //fitter->GetStats(amin, edm, errdef, nvpar, nparx);
  //cout << amin << endl;

  Double_t qTot = cluster.Charge(), par[9] = {0.}, err[9] = {0.}, coef = 0.;
  //par[8] = nMax;
  for (Int_t j = 0; j < nMax; ++j) {
    indx = 3 * j;
    par[indx+2] = fitter->GetParameter(indx+2);
    coef = Param2Coef(j, coef, par, nMax);
    par[indx] = fitter->GetParameter(indx);
    par[indx+1] = fitter->GetParameter(indx+1);
    err[indx] = fitter->GetParError(indx);
    err[indx+1] = fitter->GetParError(indx+1);
      
    if ( coef*qTot >= 14 )
    {
      AliMUONCluster* cluster1 = new AliMUONCluster(cluster);
      
      cluster1->SetCharge(coef*qTot,coef*qTot);
      
      cluster1->SetPosition(TVector2(par[indx],par[indx+1]),TVector2(err[indx],err[indx+1]));
      cluster1->SetChi2(0.);
      
      // FIXME: we miss some information in this cluster, as compared to 
      // the original AddRawCluster code.
      
      AliDebug(2,Form("Adding RawCluster detElemId %4d mult %2d charge %5d (xl,yl)=(%9.6g,%9.6g)",
		      fDetElemId,cluster1->Multiplicity(),(Int_t)cluster1->Charge(),
		      cluster1->Position().X(),cluster1->Position().Y()));
        
      fClusterList.Add(cluster1);
    }
  }
}

//_____________________________________________________________________________
void AliMUONClusterFinderPeakFit::FindClusterCOG(AliMUONCluster& cluster, 
						 Int_t *localMax, Int_t iMax)
{
/// Find COG of pad charge distribution around local maximum \a iMax 

  //TH2D *hist = (TH2D*) gROOT->FindObject("anode");
  /* Just for check
  TCanvas* c = new TCanvas("Anode","Anode",800,600);
  c->cd();
  hist->Draw("lego1Fb"); // debug
  c->Update();
  Int_t tmp;
  cin >> tmp;
  */
  Int_t nx = fHistAnode->GetNbinsX();
  //Int_t ny = hist->GetNbinsY();
  Int_t ic = localMax[iMax] / nx + 1;
  Int_t jc = localMax[iMax] % nx + 1;

  // Get min pad dimensions for the precluster
  Int_t nSides = 2;
  if (cluster.Multiplicity(0) == 0 || cluster.Multiplicity(1) == 0) nSides = 1;
  TVector2 dim0 = cluster.MinPadDimensions(0, -1, kFALSE);
  TVector2 dim1 = cluster.MinPadDimensions(1, -1, kFALSE);
  //Double_t width[2][2] = {{dim0.X(), dim0.Y()},{dim1.X(),dim1.Y()}};
  Int_t nonb[2] = {1, 0}; // coordinate index vs cathode
  if (nSides == 1 || dim0.X() < dim1.X() - fgkDistancePrecision) {
    nonb[0] = 0;
    nonb[1] = 1;
  }
  Bool_t samex = kFALSE, samey = kFALSE;
  if (TMath::Abs(dim0.X()-dim1.X()) < fgkDistancePrecision) samex = kTRUE; // the same X pad size on both planes 
  if (TMath::Abs(dim0.Y()-dim1.Y()) < fgkDistancePrecision) samey = kTRUE; // the same Y pad size on both planes 

  // Drop all pixels from the array - pick up only the ones from the cluster
  //fPixArray->Delete();

  Double_t wx = fHistAnode->GetXaxis()->GetBinWidth(1)/2; 
  Double_t wy = fHistAnode->GetYaxis()->GetBinWidth(1)/2;  
  Double_t yc = fHistAnode->GetYaxis()->GetBinCenter(ic);
  Double_t xc = fHistAnode->GetXaxis()->GetBinCenter(jc);
  Double_t cont = fHistAnode->GetCellContent(jc,ic);
  AliMUONPad* pixPtr = new AliMUONPad (xc, yc, wx, wy, cont);
  if (fDebug) pixPtr->Print("full"); 

  Int_t npad = cluster.Multiplicity();
  
  // Pick up pads which overlap with the maximum pixel and find pads with the max signal
  Double_t qMax[2] = {0}; 
  AliMUONPad *matrix[2][3] = {{0x0,0x0,0x0},{0x0,0x0,0x0}};
  for (Int_t j = 0; j < npad; ++j) 
  {
    AliMUONPad* pad = cluster.Pad(j);
    if ( Overlap(*pad,*pixPtr) )
    {
      if (fDebug) { cout << j << " "; pad->Print("full"); }
      if (pad->Charge() > qMax[pad->Cathode()]) {
	qMax[pad->Cathode()] = pad->Charge();
	matrix[pad->Cathode()][1] = pad;
	if (nSides == 1) matrix[!pad->Cathode()][1] = pad;
      }
    }
  }
  //if (nSides == 2 && (matrix[0][1] == 0x0 || matrix[1][1] == 0x0)) return; // ???

  // Find neighbours of maxima to have 3 pads per direction (if possible)
  for (Int_t j = 0; j < npad; ++j) 
  {
    AliMUONPad* pad = cluster.Pad(j);
    Int_t cath = pad->Cathode();
    if (pad == matrix[cath][1]) continue;
    Int_t nLoops = 3 - nSides;

    for (Int_t k = 0; k < nLoops; ++k) {
      Int_t cath1 = cath;
      if (k) cath1 = !cath;

      // Check the coordinate corresponding to the cathode (bending or non-bending case)
      Double_t dist = pad->Coord(nonb[cath1]) - matrix[cath][1]->Coord(nonb[cath1]);
      Double_t dir = TMath::Sign (1., dist);
      dist = TMath::Abs(dist) - pad->Size(nonb[cath1]) - matrix[cath][1]->Size(nonb[cath1]);

      if (TMath::Abs(dist) < fgkDistancePrecision) {
	// Check the other coordinate
	dist = pad->Coord(!nonb[cath1]) - matrix[cath1][1]->Coord(!nonb[cath1]);
	if (TMath::Abs(dist) > 
	    TMath::Max(pad->Size(!nonb[cath1]), matrix[cath1][1]->Size(!nonb[cath1])) - fgkDistancePrecision) break;
	Int_t idir = TMath::Nint (dir);
	if (matrix[cath1][1+idir] == 0x0) matrix[cath1][1+idir] = pad;
	else if (pad->Charge() > matrix[cath1][1+idir]->Charge()) matrix[cath1][1+idir] = pad; // diff. segmentation
	//cout << pad->Coord(nonb[cath1]) << " " << pad->Coord(!nonb[cath1]) << " " << pad->Size(nonb[cath1]) << " " << pad->Size(!nonb[cath1]) << " " << pad->Charge() << endl ;
	break;
      }
    }
  }

  Double_t coord[2] = {0.}, qAver = 0.;
  for (Int_t i = 0; i < 2; ++i) {
    Double_t q = 0.;
    Double_t coordQ = 0.;
    Int_t cath = matrix[i][1]->Cathode();
    if (i && nSides == 1) cath = !cath;
    for (Int_t j = 0; j < 3; ++j) {
      if (matrix[i][j] == 0x0) continue;
      Double_t dq = matrix[i][j]->Charge();
      q += dq;
      coordQ += dq * matrix[i][j]->Coord(nonb[cath]);
      //coordQ += (matrix[i][j]->Charge() * matrix[i][j]->Coord(nonb[cath]));
    }
    coord[cath] = coordQ / q;
    qAver = TMath::Max (qAver, q);
  }

  //qAver = TMath::Sqrt(qAver);
  if ( qAver >= 14 )
  {
    
    AliMUONCluster* cluster1 = new AliMUONCluster(cluster);
      
    cluster1->SetCharge(qAver,qAver);
    if (nonb[0] == 1) 
      cluster1->SetPosition(TVector2(coord[1],coord[0]),TVector2(0.,0.));
    else 
      cluster1->SetPosition(TVector2(coord[0],coord[1]),TVector2(0.,0.));

    cluster1->SetChi2(0.);
      
    // FIXME: we miss some information in this cluster, as compared to 
    // the original AddRawCluster code.
      
    AliDebug(2,Form("Adding RawCluster detElemId %4d mult %2d charge %5d (xl,yl)=(%9.6g,%9.6g)",
		    fDetElemId,cluster1->Multiplicity(),(Int_t)cluster1->Charge(),
		    cluster1->Position().X(),cluster1->Position().Y()));
        
    fClusterList.Add(cluster1);
  }
}

//_____________________________________________________________________________
AliMUONClusterFinderPeakFit&  
AliMUONClusterFinderPeakFit::operator=(const AliMUONClusterFinderPeakFit& rhs)
{
/// Protected assignement operator

  if (this == &rhs) return *this;

  AliFatal("Not implemented.");
    
  return *this;  
}    

//_____________________________________________________________________________
void AliMUONClusterFinderPeakFit::PadsInXandY(AliMUONCluster& cluster,
                                           Int_t &nInX, Int_t &nInY) const
{
  /// Find number of pads in X and Y-directions (excluding virtual ones and
  /// overflows)

  //Int_t statusToTest = 1;
  Int_t statusToTest = fgkUseForFit;
  
  //if ( nInX < 0 ) statusToTest = 0;
  if ( nInX < 0 ) statusToTest = fgkZero;
       
  Bool_t mustMatch(kTRUE);

  AliMpIntPair cn = cluster.NofPads(statusToTest,mustMatch);
  
  nInX = cn.GetFirst();
  nInY = cn.GetSecond();
}

//_____________________________________________________________________________
void AliMUONClusterFinderPeakFit::RemovePixel(Int_t i)
{
  /// Remove pixel at index i
  AliMUONPad* pixPtr = Pixel(i);
  fPixArray->RemoveAt(i); 
  delete pixPtr;
}

//_____________________________________________________________________________
AliMUONPad* 
AliMUONClusterFinderPeakFit::Pixel(Int_t i) const
{
  /// Returns pixel at index i
  return static_cast<AliMUONPad*>(fPixArray->UncheckedAt(i));
}

//_____________________________________________________________________________
void 
AliMUONClusterFinderPeakFit::Print(Option_t* what) const
{
  /// printout
  TString swhat(what);
  swhat.ToLower();
  if ( swhat.Contains("precluster") )
  {
    if ( fPreCluster) fPreCluster->Print();
  }
}


