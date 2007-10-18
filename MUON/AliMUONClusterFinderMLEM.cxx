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
/// \class AliMUONClusterFinderMLEM
/// 
/// Clusterizer class based on the Expectation-Maximization algorithm
///
/// Pre-clustering is handled by AliMUONPreClusterFinder
/// From a precluster a pixel array is built, and from this pixel array
/// a list of clusters is output (using AliMUONClusterSplitterMLEM).
///
/// \author Laurent Aphecetche (for the "new" C++ structure) and 
/// Alexander Zinchenko, JINR Dubna, for the hardcore of it ;-)
//-----------------------------------------------------------------------------

#include "AliMUONClusterFinderMLEM.h"
#include "AliLog.h"
#include "AliMUONCluster.h"
#include "AliMUONClusterSplitterMLEM.h"
#include "AliMUONVDigit.h"
#include "AliMUONPad.h"
#include "AliMUONPreClusterFinder.h"
#include "AliMpPad.h"
#include "AliMpVPadIterator.h"
#include "AliMpVSegmentation.h"
#include "AliRunLoader.h"
#include "AliMUONVDigitStore.h"
#include <Riostream.h>
#include <TH2.h>
#include <TMinuit.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TROOT.h>
#include "AliCodeTimer.h"

/// \cond CLASSIMP
ClassImp(AliMUONClusterFinderMLEM)
/// \endcond
 
const Double_t AliMUONClusterFinderMLEM::fgkZeroSuppression = 6; // average zero suppression value
const Double_t AliMUONClusterFinderMLEM::fgkSaturation = 3000; // average saturation level
//const Double_t AliMUONClusterFinderMLEM::fgkDistancePrecision = 1e-6; // (cm) used to check overlaps and so on
const Double_t AliMUONClusterFinderMLEM::fgkDistancePrecision = 1e-3; // (cm) used to check overlaps and so on
const TVector2 AliMUONClusterFinderMLEM::fgkIncreaseSize(-AliMUONClusterFinderMLEM::fgkDistancePrecision,-AliMUONClusterFinderMLEM::fgkDistancePrecision);
const TVector2 AliMUONClusterFinderMLEM::fgkDecreaseSize(AliMUONClusterFinderMLEM::fgkDistancePrecision,AliMUONClusterFinderMLEM::fgkDistancePrecision);

 TMinuit* AliMUONClusterFinderMLEM::fgMinuit = 0x0;

//_____________________________________________________________________________
AliMUONClusterFinderMLEM::AliMUONClusterFinderMLEM(Bool_t plot, AliMUONVClusterFinder* clusterFinder)
  : AliMUONVClusterFinder(),
fPreClusterFinder(clusterFinder),
fPreCluster(0x0),
fClusterList(),
fEventNumber(0),
fDetElemId(-1),
fClusterNumber(0),
fZpad(0.0),
fReco(1),
fCathBeg(0),
fPixArray(new TObjArray(20)),
fDebug(0),
fPlot(plot),
fSplitter(0x0),
fNClusters(0),
fNAddVirtualPads(0)
{
  /// Constructor
 
  fSegmentation[1] = fSegmentation[0] = 0x0; 

  fPadBeg[0] = fPadBeg[1] = fCathBeg;

  if (!fgMinuit) fgMinuit = new TMinuit(8);

  if (fPlot) fDebug = 1;
}

//_____________________________________________________________________________
AliMUONClusterFinderMLEM::~AliMUONClusterFinderMLEM()
{
/// Destructor
  delete fgMinuit; fgMinuit = 0; delete fPixArray; fPixArray = 0;
//  delete fDraw;
  delete fPreClusterFinder;
  delete fSplitter;
  AliInfo(Form("Total clusters %d AddVirtualPad needed %d",
               fNClusters,fNAddVirtualPads));
}

//_____________________________________________________________________________
Bool_t 
AliMUONClusterFinderMLEM::Prepare(const AliMpVSegmentation* segmentations[2],
                                  const AliMUONVDigitStore& digitStore)
{
  /// Prepare for clustering
//  AliCodeTimerAuto("")
  
  for ( Int_t i = 0; i < 2; ++i )
  {
    fSegmentation[i] = segmentations[i];
  }
  
  // Find out the DetElemId
  fDetElemId = -1;
  
  TIter next(digitStore.CreateIterator());
  AliMUONVDigit* d = static_cast<AliMUONVDigit*>(next());
  
  if (d)
  {
    fDetElemId = d->DetElemId();
  }
  
  if ( fDetElemId < 0 )
  {
    AliWarning("Could not find DE. Probably no digits at all ?");
    return kFALSE;
  }
  
  delete fSplitter;
  fSplitter = new AliMUONClusterSplitterMLEM(fDetElemId,fPixArray);
  fSplitter->SetDebug(fDebug);
    
  // find out current event number, and reset the cluster number
  fEventNumber = AliRunLoader::GetRunLoader()->GetEventNumber();
  fClusterNumber = -1;
  fClusterList.Delete();
  
  AliDebug(3,Form("EVT %d DE %d",fEventNumber,fDetElemId));
  
  return fPreClusterFinder->Prepare(segmentations,digitStore);
}

//_____________________________________________________________________________
AliMUONCluster* 
AliMUONClusterFinderMLEM::NextCluster()
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
  fClusterNumber = -1; //AZ
  
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
AliMUONClusterFinderMLEM::WorkOnPreCluster()
{
  /// Starting from a precluster, builds a pixel array, and then
  /// extract clusters from this array
  
  //  AliCodeTimerAuto("")	

  if (fDebug) {
    cout << " *** Event # " << AliRunLoader::GetRunLoader()->GetEventNumber() 
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
  
  // Use MLEM for cluster finder
  Int_t nMax = 1, localMax[100], maxPos[100];
  Double_t maxVal[100];
  
  if (cluster->Multiplicity() > 50) 
  {
    nMax = FindLocalMaxima(fPixArray, localMax, maxVal);
  }
  //nMax = 1; // just for test
  
  if (nMax > 1) 
  {
    TMath::Sort(nMax, maxVal, maxPos, kTRUE); // in decreasing order
  }
  
  Int_t iSimple = 0, nInX = -1, nInY;
  
  PadsInXandY(*cluster,nInX, nInY);
  
  if (nMax == 1 && nInX < 4 && nInY < 4) 
  {
    iSimple = 1; //1; // simple cluster
  }
  
  for (Int_t i = 0; i < nMax; ++i) 
  {
    if (nMax > 1) 
    {
      FindCluster(*cluster,localMax, maxPos[i]);
    }

    MainLoop(*cluster,iSimple);

    if (i < nMax-1) 
    {
      Int_t mult = cluster->Multiplicity();
      for (Int_t j = 0; j < mult; ++j) 
      {
        AliMUONPad* pad = cluster->Pad(j);
        if ( pad->Status() == 0 ) continue; // pad charge was not modified
        pad->SetStatus(0);
        pad->RevertCharge(); // use backup charge value
      }
    }
  } // for (Int_t i=0; i<nMax;
  if (nMax > 1) ((TH2D*) gROOT->FindObject("anode"))->Delete();
  //TH2D *mlem = (TH2D*) gROOT->FindObject("mlem");
  //if (mlem) mlem->Delete();
  delete cluster;
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t 
AliMUONClusterFinderMLEM::Overlap(const AliMUONPad& pad, const AliMUONPad& pixel)
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
AliMUONClusterFinderMLEM::CheckPrecluster(const AliMUONCluster& origCluster)
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

  AliMUONCluster* cluster = static_cast<AliMUONCluster*>(origCluster.Clone());

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
AliMUONClusterFinderMLEM::CheckPreclusterTwoCathodes(AliMUONCluster* cluster)
{
  /// Check two-cathode cluster
  
  Int_t i1 = cluster->Multiplicity(0) ? 0 : 1;
  Int_t i2 = cluster->Multiplicity(1) ? 1 : 0;
  
  Int_t npad = cluster->Multiplicity();
  Int_t* flags = new Int_t[npad];
  for (Int_t j = 0; j < npad; ++j) flags[j] = 0;
  
  // Check pad overlaps
  for ( Int_t i = 0; i < npad; ++i) 
  {
    AliMUONPad* padi = cluster->Pad(i);
    if ( padi->Cathode() != i1 ) continue;
    for (Int_t j = i+1; j < npad; ++j) 
    {
      AliMUONPad* padj = cluster->Pad(j);
      if ( padj->Cathode() != i2 ) continue;
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
      if (nFlags == 1 && pad->Charge() < fgkZeroSuppression * 3) continue;
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
    Int_t imin(0);
    Int_t imax(0);
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
      }
      if ( pad->Charge() > cmax )
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
  }
  
  delete[] flags;
  
  AliDebug(2,"End of CheckPreClusterTwoCathodes=");
  //StdoutToAliDebug(2,cluster->Print("full"));

  return cluster;    
}

//_____________________________________________________________________________
void
AliMUONClusterFinderMLEM::CheckOverlaps()
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
void AliMUONClusterFinderMLEM::BuildPixArray(AliMUONCluster& cluster)
{
  /// Build pixel array for MLEM method
  
  Int_t npad = cluster.Multiplicity();
  if (npad<=0) 
  {
    AliWarning("Got no pad at all ?!");
  }
  
  fPixArray->Delete();
  
  if ( !cluster.Multiplicity(0) || !cluster.Multiplicity(1) )
  {
    BuildPixArrayOneCathode(cluster);
  }
  else
  {
    //BuildPixArrayTwoCathodes(cluster); 
    BuildPixArrayOneCathode(cluster);
  }
  
  //fPixArray->Sort(); // FIXME : not really needed, only to compare with ClusterFinderAZ
  
  Int_t nPix = fPixArray->GetLast()+1;
  
//  AliDebug(2,Form("nPix after BuildPixArray=%d",nPix));
  
  Double_t xPadMin(1E9);
  Double_t yPadMin(1E9);
  
  //for ( Int_t i = 0; i < cluster.Multiplicity(); ++i ) 
  for ( Int_t i = 0; i < npad; ++i ) 
  {
    AliMUONPad* pad = cluster.Pad(i);
    xPadMin = TMath::Min (xPadMin, pad->DX());
    yPadMin = TMath::Min (yPadMin, pad->DY());
  }
  
  Double_t wxmin(1E9);
  Double_t wymin(1E9);
  
  for ( Int_t i = 0; i < nPix; ++i ) 
  {
    AliMUONPad* pixPtr = Pixel(i);
    wxmin = TMath::Min(wxmin, pixPtr->Size(0));
    wymin = TMath::Min(wymin, pixPtr->Size(1));
  }

  //wxmin = TMath::Abs (wxmin - xPadMin/2) > 0.001 ? xPadMin : xPadMin / 2;
  //wymin = TMath::Abs (wymin - yPadMin/2) > 0.001 ? yPadMin : yPadMin / 2;
  wxmin = xPadMin;
  wymin = yPadMin;
  
  /*
  // Check if small pixel X-size
  AdjustPixel(cluster,wxmin, 0);
  // Check if small pixel Y-size
  AdjustPixel(cluster,wymin, 1);
  // Check if large pixel size
  AdjustPixel(wxmin, wymin);
  */
  
  // Remove discarded pixels
  for (Int_t i = 0; i < nPix; ++i) 
  {
    AliMUONPad* pixPtr = Pixel(i);
    if (pixPtr->Charge() < 1) 
    { 
      AliDebug(2,Form("Removing pixel %d with charge<1 : ",i));
      //StdoutToAliDebug(2,pixPtr->Print());
      RemovePixel(i);
    }
  }
  
  fPixArray->Compress();
  nPix = fPixArray->GetEntriesFast();
  
//  AliDebug(2,Form("nPix after AdjustPixel=%d",nPix));

  //if ( nPix > cluster.Multiplicity() ) 
  if ( nPix > npad ) 
  {
//    AliDebug(2,Form("Will trim number of pixels to number of pads"));
    
    // Too many pixels - sort and remove pixels with the lowest signal
    fPixArray->Sort();
    for ( Int_t i = cluster.Multiplicity(); i < nPix; ++i ) 
    {
      RemovePixel(i);
    }
    fPixArray->Compress();
    nPix = fPixArray->GetEntriesFast();
  } // if (nPix > npad)

//  StdoutToAliDebug(2,cout << "End of BuildPixelArray:" << endl;
//                   fPixArray->Print(););
  //CheckOverlaps();//FIXME : this is for debug only. Remove it.
}

//_____________________________________________________________________________
void AliMUONClusterFinderMLEM::BuildPixArrayOneCathode(AliMUONCluster& cluster)
{
  /// Build the pixel array

//  AliDebug(2,Form("cluster.Multiplicity=%d",cluster.Multiplicity()));

  // Find min and max cluster dimensions
  Double_t minx[2] = {9999,9999}, maxx[2] = {-9999,-9999};
  Double_t miny[2] = {9999,9999}, maxy[2] = {-9999,-9999};

  TVector2 dim = cluster.MinPadDimensions (-1, kFALSE);
  Double_t width[2] = {dim.X(), dim.Y()}, xy0[2];
  Int_t found[2] = {0};

  for ( Int_t i = 0; i < cluster.Multiplicity(); ++i) {
    AliMUONPad* pad = cluster.Pad(i);
    Int_t cath = pad->Cathode();
    minx[cath] = TMath::Min (minx[cath], pad->Coord(0)-pad->Size(0));
    maxx[cath] = TMath::Max (maxx[cath], pad->Coord(0)+pad->Size(0));
    miny[cath] = TMath::Min (miny[cath], pad->Coord(1)-pad->Size(1));
    maxy[cath] = TMath::Max (maxy[cath], pad->Coord(1)+pad->Size(1));
    for (Int_t j = 0; j < 2; ++j) {
      if (found[j] == 0 && TMath::Abs(pad->Size(j)-width[j]) < fgkDistancePrecision) { 
	xy0[j] = pad->Coord(j);
	found[j] = 1;
      }
    }
  }
  /*
  TVector2 leftDown = cluster.Area(0).LeftDownCorner();
  TVector2 rightUp = cluster.Area(0).RightUpCorner();
  cout << leftDown.X() << " " << leftDown.Y() << " " << rightUp.X() << " " << rightUp.Y() << endl;
  leftDown = cluster.Area(1).LeftDownCorner();
  rightUp = cluster.Area(1).RightUpCorner();
  cout << leftDown.X() << " " << leftDown.Y() << " " << rightUp.X() << " " << rightUp.Y() << endl;
  */

  //cout << minx[0] << " " << maxx[0] << " " << minx[1] << " " << maxx[1] << endl;
  //cout << miny[0] << " " << maxy[0] << " " << miny[1] << " " << maxy[1] << endl;
  //cout << width[0] << " " << width[1] << endl;
  Double_t min[2], max[2];
  Int_t cath0 = 0, cath1 = 1;
  if (cluster.Multiplicity(0) == 0) cath0 = 1;
  else if (cluster.Multiplicity(1) == 0) cath1 = 0;
  min[0] = TMath::Max (minx[cath0], minx[cath1]);
  min[1] = TMath::Max (miny[cath0], miny[cath1]);
  max[0] = TMath::Min (maxx[cath0], maxx[cath1]);
  max[1] = TMath::Min (maxy[cath0], maxy[cath1]);

  // Adjust limits
  //width[0] /= 2; width[1] /= 2; // just for check
  Int_t nbins[2];
  for (Int_t i = 0; i < 2; ++i) {
    Double_t dist = (min[i] - xy0[i]) / width[i] / 2;
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
  Int_t mult = cluster.Multiplicity();
  for ( Int_t i = 0; i < mult; ++i) {
    AliMUONPad* pad = cluster.Pad(i);
    Int_t ix0 = xaxis->FindBin(pad->X());
    Int_t iy0 = yaxis->FindBin(pad->Y());
    PadOverHist(0, ix0, iy0, pad);
  }

  // Store pixels
  for (Int_t i = 1; i <= nbins[0]; ++i) {
    Double_t x = xaxis->GetBinCenter(i);
    for (Int_t j = 1; j <= nbins[1]; ++j) {
      if (hist2->GetCellContent(i,j) < 0.1) continue;
      if (hist2->GetCellContent(i,j) < 1.1 && cluster.Multiplicity(0) && 
	  cluster.Multiplicity(1)) continue;
      Double_t y = yaxis->GetBinCenter(j);
      Double_t charge = hist1->GetCellContent(i,j);
      AliMUONPad* pixPtr = new AliMUONPad(x, y, width[0], width[1], charge);
      fPixArray->Add(pixPtr);
    }  
  }
  //fPixArray->Print();
  delete hist1;
  delete hist2;
}

//_____________________________________________________________________________
void AliMUONClusterFinderMLEM::PadOverHist(Int_t idir, Int_t ix0, Int_t iy0, AliMUONPad *pad)
{
  /// "Span" pad over histogram in the direction idir

  TH2D *hist1 = static_cast<TH2D*> (gROOT->FindObject("Grid"));
  TH2D *hist2 = static_cast<TH2D*> (gROOT->FindObject("Entries"));
  TAxis *axis = idir == 0 ? hist1->GetXaxis() : hist1->GetYaxis();
  Int_t nbins = axis->GetNbins();
  Double_t bin = axis->GetBinWidth(1);

  Int_t nbinPad = (Int_t)(pad->Size(idir)/bin*2+fgkDistancePrecision) + 1; // number of bins covered by pad

  for (Int_t i = 0; i < nbinPad; ++i) {
    Int_t ixy = idir == 0 ? ix0 + i : iy0 + i;
    if (ixy > nbins) break;
    Double_t lowEdge = axis->GetBinLowEdge(ixy);
    if (lowEdge + fgkDistancePrecision > pad->Coord(idir) + pad->Size(idir)) break;
    if (idir == 0) PadOverHist(1, ixy, iy0, pad); // span in the other direction
    else {
      // Fill histogram
      Double_t cont = pad->Charge();
      if (hist2->GetCellContent(ix0, ixy) > 0.1) 
	cont = TMath::Min (hist1->GetCellContent(ix0, ixy), cont);
      hist1->SetCellContent(ix0, ixy, cont);
      hist2->SetCellContent(ix0, ixy, hist2->GetCellContent(ix0, ixy)+1);
    }
  }

  for (Int_t i = -1; i > -nbinPad; --i) {
    Int_t ixy = idir == 0 ? ix0 + i : iy0 + i;
    if (ixy < 1) break;
    Double_t upEdge = axis->GetBinUpEdge(ixy);
    if (upEdge - fgkDistancePrecision < pad->Coord(idir) - pad->Size(idir)) break;
    if (idir == 0) PadOverHist(1, ixy, iy0, pad); // span in the other direction
    else {
      // Fill histogram
      Double_t cont = pad->Charge();
      if (hist2->GetCellContent(ix0, ixy) > 0.1) 
	cont = TMath::Min (hist1->GetCellContent(ix0, ixy), cont);
      hist1->SetCellContent(ix0, ixy, cont);
      hist2->SetCellContent(ix0, ixy, hist2->GetCellContent(ix0, ixy)+1);
    }
  }
}

/*
//_____________________________________________________________________________
void AliMUONClusterFinderMLEM::BuildPixArrayOneCathode(AliMUONCluster& cluster)
{
  /// From a single-cathode cluster, build the pixel array

//  AliDebug(2,Form("cluster.Multiplicity=%d",cluster.Multiplicity()));

  for ( Int_t j=0; j<cluster.Multiplicity(); ++j) 
  {
    AliMUONPad* pad = cluster.Pad(j);
    AliMUONPad* pixPtr = new AliMUONPad(pad->Position(),pad->Dimensions(),
                                            pad->Charge());    
    fPixArray->Add(pixPtr);
  }  
}
*/

//_____________________________________________________________________________
void AliMUONClusterFinderMLEM::BuildPixArrayTwoCathodes(AliMUONCluster& cluster)
{
  /// From a two-cathodes cluster, build the pixel array
  
//  AliDebug(2,Form("cluster.Multiplicity=%d",cluster.Multiplicity()));
           
  Int_t i1 = cluster.Pad(0)->Cathode();
  Int_t i2 = TMath::Even(i1);
  
  for ( Int_t i = 0; i < cluster.Multiplicity(); ++i) 
  {
    AliMUONPad* padi = cluster.Pad(i);
    if (padi->Cathode() != i1) continue;

    for ( Int_t j = 1; j < cluster.Multiplicity(); ++j) 
    {
      AliMUONPad* padj = cluster.Pad(j);
      if (padj->Cathode() != i2) continue;

      AliMpArea overlap;

      if ( AliMUONPad::AreOverlapping(*padi,*padj,fgkDecreaseSize,overlap) )
      {      
        AliMUONPad* pixPtr = new AliMUONPad(overlap.Position(),
                                                overlap.Dimensions(),
                                                TMath::Min(padi->Charge(),padj->Charge()));
        if ( ( padi->Charge() <= padj->Charge() && padi->IsSaturated() ) ||
             ( padj->Charge() < padi->Charge() && padj->IsSaturated() ) )
        {
          // if smallest charge of the 2 pads is already above saturation, then
          // the pixel is saturated...
          pixPtr->SetSaturated(kTRUE);
        }
        pixPtr->SetReal(kFALSE);
        fPixArray->Add(pixPtr);
      }
    } 
  } 
}  

//_____________________________________________________________________________
void AliMUONClusterFinderMLEM::AdjustPixel(AliMUONCluster& /*cluster*/, 
                                           Float_t width, Int_t ixy)
{
  /// Check if some pixels have smaller size than others (adjust if necessary)

  AliDebug(2,Form("width=%e ixy=%d",width,ixy));
  
  AliMUONPad *pixPtr, *pixPtr1 = 0;
  Int_t ixy1 = !ixy;
  Int_t nPix = fPixArray->GetEntriesFast(), iOK = 1;
  
  Double_t xy0 = 0, minmax[2] = {9999,-9999}, dist = 0;
  // First, find a "normal" pixel
  for (Int_t i = 0; i < nPix; ++i) {
    pixPtr = Pixel(i);
    if (pixPtr->Charge() < 1) continue; // discarded pixel
    minmax[0] = TMath::Min (minmax[0], pixPtr->Size(ixy));
    minmax[1] = TMath::Max (minmax[1], pixPtr->Size(ixy));
    if (pixPtr->Size(ixy) - width < -fgkDistancePrecision) iOK = 0;
    if (TMath::Abs(pixPtr->Size(ixy)-width) > fgkDistancePrecision) continue;
    xy0 = pixPtr->Coord(ixy);
  }
  if (TMath::Abs(minmax[0]-minmax[1]) < fgkDistancePrecision) iOK = 1; // the same size
  if (iOK == 1) return; // all pixels have the same size in the direction IXY
  
  //cout << " --- " << xy0 << endl; fPixArray->Print();
  for (Int_t i = 0; i < nPix; ++i) 
  {
    pixPtr = Pixel(i);
    if (pixPtr->Charge() < 1) continue; // discarded pixel
    if (pixPtr->Size(ixy) - width < -fgkDistancePrecision) 
    {
      // try to merge 
      if (fDebug) cout << i << " Small X or Y: " << ixy << " " << pixPtr->Size(ixy) << " " << width << " " << pixPtr->Coord(0) << " " << pixPtr->Coord(1) << endl;
      for (Int_t j = i + 1; j < nPix; ++j) 
      {
        pixPtr1 = Pixel(j);
        if (pixPtr1->Charge() < 1) continue; // discarded pixel
        if (TMath::Abs(pixPtr1->Size(ixy)-width) < fgkDistancePrecision) continue; // right size 
        if (TMath::Abs(pixPtr1->Coord(ixy1)-pixPtr->Coord(ixy1)) > fgkDistancePrecision) continue; // different rows/columns
	dist = TMath::Abs (pixPtr1->Coord(ixy) - pixPtr->Coord(ixy));
        if (TMath::Abs(dist-pixPtr1->Size(ixy)-pixPtr->Size(ixy)) < fgkDistancePrecision) // neighbours 
        {
          // merge
          //Double_t dist = (pixPtr->Coord(ixy) + pixPtr1->Coord(ixy)) / 2;
	  //dist = TMath::Nint((dist-xy0)/width/2) * width * 2;
	  dist = (pixPtr->Coord(ixy)-xy0) / width / 2;
	  cout << j << " " << dist << endl;
	  dist = TMath::Nint(dist) * width * 2;
          pixPtr->SetCoord(ixy, xy0+dist);
          pixPtr->SetSize(ixy, width);
          pixPtr->SetCharge(TMath::Min (pixPtr->Charge(),pixPtr1->Charge()));
          pixPtr1->SetCharge(0);
          pixPtr1 = 0x0;
          break;
        }
      } // for (Int_t j = i + 1;
      if (pixPtr1 || i == nPix-1) {
        // edge pixel - just increase its size
	if (fDebug) cout << " No pair ..." << endl; 
	cout << (pixPtr->Coord(ixy)-xy0)/width/2 << endl;
	dist = (pixPtr->Coord(ixy) - xy0) / width / 2;
	dist = TMath::Nint(dist) * width * 2;
	pixPtr->SetCoord(ixy, xy0+dist);
	pixPtr->SetSize(ixy, width);
      }
    } // if (pixPtr->Size(ixy)-width < -fgkDistancePrecision)
  } // for (Int_t i = 0; i < nPix;
  //cout << " *** " << endl; fPixArray->Print();
}

//_____________________________________________________________________________
void AliMUONClusterFinderMLEM::AdjustPixel(Double_t wxmin, Double_t wymin)
{
/// Check if some pixels have large size (adjust if necessary)

  AliDebug(2,Form("wxmin=%e wymin=%e",wxmin,wymin));
  
  Int_t n2[2], iOK = 1, nPix = fPixArray->GetEntriesFast();
  AliMUONPad *pixPtr, pix;
  Double_t xy0[2] = {9999, 9999}, wxy[2], dist[2] = {0};

  // Check if large pixel size
  for (Int_t i = 0; i < nPix; i++) {
    pixPtr = (AliMUONPad*) fPixArray->UncheckedAt(i);
    if (pixPtr->Charge() < 1) continue; // discarded pixel
    if (pixPtr->Size(0) - wxmin < 1.e-4) {
      if (xy0[0] > 9998) xy0[0] = pixPtr->Coord(0); // position of a "normal" pixel
      if (pixPtr->Size(1) - wymin < 1.e-4) { 
	if (xy0[1] > 9998) xy0[1] = pixPtr->Coord(1); // position of a "normal" pixel
	continue;
      } else iOK = 0; // large pixel
    } else {
      iOK = 0; // large pixel
      if (xy0[1] > 9998 && pixPtr->Size(1) - wymin < 1.e-4) xy0[1] = pixPtr->Coord(1); // "normal" pixel
    }      
    if (xy0[0] < 9998 && xy0[1] < 9998) break;
  }
  if (iOK) return;

  wxy[0] = wxmin;
  wxy[1] = wymin;
  Int_t update[2] = {0};
  //cout << " --- " << endl; fPixArray->Print();
  cout << xy0[0] << " " << xy0[1] << endl;
  for (Int_t i = 0; i < nPix; i++) {
    pixPtr = (AliMUONPad*) fPixArray->UncheckedAt(i);
    if (pixPtr->Charge() < 1) continue; // discarded pixel
    n2[0] = n2[1] = 1;
    update[0] = update[1] = 0;
    for (Int_t j = 0; j < 2; j++) {
      if (pixPtr->Size(j) - wxy[j] < 1.e-4) continue;
      dist[j] = pixPtr->Coord(j) - xy0[j]; // distance to "normal" pixel
      // Go back to position of the first updated pixel
      dist[j] += (pixPtr->Size(j) - wxy[j]) * TMath::Sign(1.,-dist[j]);
      n2[j] = TMath::Nint (pixPtr->Size(j) / wxy[j]);
      update[j] = 1;
    }
    if (update[0] == 0 && update[1] == 0) continue;
    if (fDebug) cout << " Different " << pixPtr->Size(0) << " " << wxy[0] << " "
		     << pixPtr->Size(1) << " " << wxy[1] <<endl;
    
    pix = *pixPtr;
    pix.SetSize(0, wxy[0]); pix.SetSize(1, wxy[1]);
    //pixPtr->Print();
    for (Int_t ii = 0; ii < n2[0]; ii++) {
      if (update[0]) pix.SetCoord(0, xy0[0] + dist[0] + TMath::Sign(2.,dist[0]) * ii * wxy[0]);
      for (Int_t jj = 0; jj < n2[1]; jj++) {
	if (update[1]) pix.SetCoord(1, xy0[1] + dist[1] + TMath::Sign(2.,dist[1]) * jj * wxy[1]);
	fPixArray->Add(new AliMUONPad(pix));
	//pix.Print();
      }
    }
    pixPtr->SetCharge(0);
  } // for (Int_t i = 0; i < nPix;
  cout << " *** " << endl; fPixArray->Print();
}

//_____________________________________________________________________________
void
AliMUONClusterFinderMLEM::Plot(const char* basename)
{
  /// Make a plot and save it as png
  
  return; //AZ
  if (!fPlot) return;
  
  TCanvas* c = new TCanvas("MLEM","MLEM",800,600);
  c->Draw();
  Draw();
  c->Modified();
  c->Update();
  c->Print(Form("%s.EVT%d.DE%d.CLU%d.png",basename,fEventNumber,
                fDetElemId,fClusterNumber));
}

//_____________________________________________________________________________
void
AliMUONClusterFinderMLEM::ComputeCoefficients(AliMUONCluster& cluster,
                                              Double_t* coef,
                                              Double_t* probi)
{
  /// Compute coefficients needed for MLEM algorithm
  
  Int_t nPix = fPixArray->GetLast()+1;
  
  //memset(probi,0,nPix*sizeof(Double_t));
  for (Int_t j = 0; j < nPix; ++j) probi[j] = 0.;

  Int_t mult = cluster.Multiplicity();
  for ( Int_t j = 0; j < mult; ++j ) 
  {
    AliMUONPad* pad = cluster.Pad(j);
    Int_t indx = j*nPix;
  
    for ( Int_t ipix = 0; ipix < nPix; ++ipix ) 
    {
      Int_t indx1 = indx + ipix;
      if (pad->Status() < 0) 
      {   
        coef[indx1] = 0; 
        continue; 
      }
      AliMUONPad* pixPtr = Pixel(ipix);
      // coef is the charge (given by Mathieson integral) on pad, assuming
      // the Mathieson is center at pixel.
      coef[indx1] = fSplitter->ChargeIntegration(pixPtr->Coord(0), pixPtr->Coord(1), *pad);  
//      AliDebug(2,Form("pad=(%d,%d,%e,%e,%e,%e) pix=(%e,%e,%e,%e) coef %e",
//                      pad->Ix(),pad->Iy(),
//                      pad->X(),pad->Y(),
//                      pad->DX(),pad->DY(),
//                      pixPtr->Coord(0),pixPtr->Coord(1), 
//                      pixPtr->Size(0),pixPtr->Size(1),
//                      coef[indx1]));
      
      probi[ipix] += coef[indx1];
    } 
  } 
}

//_____________________________________________________________________________
Bool_t AliMUONClusterFinderMLEM::MainLoop(AliMUONCluster& cluster, Int_t iSimple)
{
  /// Repeat MLEM algorithm until pixel size becomes sufficiently small
  
  //  AliCodeTimerAuto("")
  
  Int_t nPix = fPixArray->GetLast()+1;

  AliDebug(2,Form("nPix=%d iSimple=%d, precluster=",nPix,iSimple));
  //StdoutToAliDebug(2,cluster.Print("full"););

  if ( nPix < 0 )
  {
    AliDebug(1,"No pixels, why am I here then ?");
  }
  
  AddVirtualPad(cluster); // add virtual pads if necessary
  
  Int_t npadTot = cluster.Multiplicity();
  Int_t npadOK = 0;
  for (Int_t i = 0; i < npadTot; ++i) 
  {
    if (cluster.Pad(i)->Status() == 0) ++npadOK;
  }

  TH2D* mlem(0x0);
  Double_t* coef(0x0);
  Double_t* probi(0x0);
  Int_t lc(0); // loop counter (for debug)
  
  //Plot("mlem.start");
  
  while (1) 
  {
    ++lc;
    mlem = (TH2D*) gROOT->FindObject("mlem");
    delete mlem;
    
    AliDebug(2,Form("lc %d nPix %d(%d) npadTot %d npadOK %d",lc,nPix,fPixArray->GetLast()+1,npadTot,npadOK));
    AliDebug(2,Form("EVT%d PixArray=",fEventNumber));
    //StdoutToAliDebug(2,fPixArray->Print("","full"));
        
    coef = new Double_t [npadTot*nPix];
    probi = new Double_t [nPix];

    // Calculate coefficients and pixel visibilities
    ComputeCoefficients(cluster,coef,probi);

    for (Int_t ipix = 0; ipix < nPix; ++ipix) 
    {
      if (probi[ipix] < 0.01) 
      {
        AliMUONPad* pixel = Pixel(ipix);
        AliDebug(2,Form("Setting the following pixel to invisible as its probi<0.01:"));
        //StdoutToAliDebug(2,cout << Form(" -- ipix %3d --- "); pixel->Print(););
        pixel->SetCharge(0); // "invisible" pixel
      }
    }
    
    // MLEM algorithm
    Mlem(cluster,coef, probi, 15);

    Double_t xylim[4] = {999, 999, 999, 999};
    AliMUONPad* pixPtr(0x0);
    
    for ( Int_t ipix = 0; ipix < nPix; ++ipix ) 
    {
      pixPtr = Pixel(ipix);
      for ( Int_t i = 0; i < 4; ++i ) 
      {
        xylim[i] = TMath::Min (xylim[i], (i%2 ? -1 : 1)*pixPtr->Coord(i/2));
      }
    }
    for (Int_t i = 0; i < 4; i++) 
    {
      xylim[i] -= pixPtr->Size(i/2); 
    }

    
    Int_t nx = TMath::Nint ((-xylim[1]-xylim[0])/pixPtr->Size(0)/2);
    Int_t ny = TMath::Nint ((-xylim[3]-xylim[2])/pixPtr->Size(1)/2);

    //StdoutToAliDebug(2,cout << "pixel used for nx,ny computation : "; pixPtr->Print(););
    AliDebug(2,Form("lc %d pixPtr size = %e,%e nx,ny=%d,%d xylim=%e,%e,%e,%e",
                    lc,pixPtr->Size(0),pixPtr->Size(1),nx,ny,
                    xylim[0],-xylim[1],xylim[2],-xylim[3]
                    ));
    
    mlem = new TH2D("mlem","mlem",nx,xylim[0],-xylim[1],ny,xylim[2],-xylim[3]);

    for (Int_t ipix = 0; ipix < nPix; ++ipix) 
    {
      AliMUONPad* pixPtr = Pixel(ipix);
      mlem->Fill(pixPtr->Coord(0),pixPtr->Coord(1),pixPtr->Charge());
    }

    // Check if the total charge of pixels is too low
    Double_t qTot = 0;
    for ( Int_t i = 0; i < nPix; ++i) 
    {
      qTot += Pixel(i)->Charge();
    }
    
    if ( qTot < 1.e-4 || ( npadOK < 3 && qTot < 7 ) ) 
    {
      AliDebug(1,Form("Deleting the above cluster (charge %e too low, npadOK=%d)",qTot,npadOK));
      delete [] coef; 
      delete [] probi; 
      fPixArray->Delete(); 
      for ( Int_t i = 0; i < npadTot; ++i) 
      {
        AliMUONPad* pad = cluster.Pad(i);
        if ( pad->Status() == 0) pad->SetStatus(-1);
      }
      return kFALSE; 
    }

    if (iSimple) 
    {
      // Simple cluster - skip further passes thru EM-procedure
      Simple(cluster);
      delete [] coef; 
      delete [] probi; 
      fPixArray->Delete(); 
      return kTRUE;
    }

    // Calculate position of the center-of-gravity around the maximum pixel
    Double_t xyCOG[2];
    FindCOG(mlem, xyCOG);

    if (TMath::Min(pixPtr->Size(0),pixPtr->Size(1)) < 0.07 && 
        pixPtr->Size(0) > pixPtr->Size(1)) break;

    // Sort pixels according to the charge
    fPixArray->Sort();
    Double_t pixMin = 0.01*Pixel(0)->Charge();
    pixMin = TMath::Min(pixMin,50.);

    // Decrease pixel size and shift pixels to make them centered at 
    // the maximum one
    Int_t indx = (pixPtr->Size(0)>pixPtr->Size(1)) ? 0 : 1;
    Int_t ix(1);
    Double_t width = 0;
    Double_t shift[2] = { 0.0, 0.0 };
    for (Int_t i = 0; i < 4; ++i) xylim[i] = 999;
    Int_t nPix1 = nPix; 
    nPix = 0;
    for (Int_t ipix = 0; ipix < nPix1; ++ipix) 
    {
      AliMUONPad* pixPtr = Pixel(ipix);
      if ( nPix >= npadOK  // too many pixels already
           ||
           pixPtr->Charge() < pixMin // too low charge
           ) 
      { 
        RemovePixel(ipix);
        continue;
      }
      for (Int_t i = 0; i < 2; ++i) 
      {
        if (!i) 
        {
          pixPtr->SetCharge(10);
          pixPtr->SetSize(indx, pixPtr->Size(indx)/2);
          width = -pixPtr->Size(indx);
          pixPtr->Shift(indx, width);
          // Shift pixel position
          if (ix) 
          {
            ix = 0;
            for (Int_t j = 0; j < 2; ++j) 
            {
              shift[j] = pixPtr->Coord(j) - xyCOG[j];
              shift[j] -= ((Int_t)(shift[j]/pixPtr->Size(j)/2))*pixPtr->Size(j)*2;
            }
          } // if (ix)
          pixPtr->Shift(0, -shift[0]);
          pixPtr->Shift(1, -shift[1]);
        } 
        else 
        {
          pixPtr = new AliMUONPad(*pixPtr);
          pixPtr->Shift(indx, -2*width);
          fPixArray->Add(pixPtr);
        } 
        for (Int_t i = 0; i < 4; ++i) 
        {
          xylim[i] = TMath::Min (xylim[i], (i%2 ? -1 : 1)*pixPtr->Coord(i/2));
        }
      } // for (Int_t i=0; i<2;
      nPix += 2;
    } // for (Int_t ipix=0;
    
    fPixArray->Compress();
    nPix = fPixArray->GetEntriesFast();

    AliDebug(2,Form("After shift:"));
    //StdoutToAliDebug(2,fPixArray->Print("","full"););
    //Plot(Form("mlem.lc%d",lc+1));
    
    AliDebug(2,Form(" xyCOG=%9.6f %9.6f xylim=%9.6f,%9.6f,%9.6f,%9.6f",
                    xyCOG[0],xyCOG[1],
                    xylim[0],xylim[1],
                    xylim[2],xylim[3]));

    // Remove excessive pixels
    if (nPix > npadOK) 
    {
      for (Int_t ipix = npadOK; ipix < nPix; ++ipix) 
      { 
        RemovePixel(ipix);
      }
    } 
    else 
    {
      AliMUONPad* pixPtr = Pixel(0);
      // add pixels if the maximum is at the limit of pixel area
      // start from Y-direction
      Int_t j = 0;
      for (Int_t i = 3; i > -1; --i) 
      {
        if (nPix < npadOK && 
            TMath::Abs((i%2 ? -1 : 1)*xylim[i]-xyCOG[i/2]) < pixPtr->Size(i/2)) 
        {
          AliMUONPad* p = static_cast<AliMUONPad*>(pixPtr->Clone());
          p->SetCoord(i/2, xyCOG[i/2]+(i%2 ? 2:-2)*pixPtr->Size(i/2));
          j = TMath::Even (i/2);
          p->SetCoord(j, xyCOG[j]);
          AliDebug(2,Form("Adding pixel on the edge (i=%d) ",i));
          //StdoutToAliDebug(2,cout << " ---- "; 
	  //               p->Print("corners"););
          fPixArray->Add(p);
          ++nPix;
        }
      }
    } 
    fPixArray->Compress();
    nPix = fPixArray->GetEntriesFast();
    delete [] coef; 
    delete [] probi; 
  } // while (1)

  AliDebug(2,Form("At the end of while loop nPix=%d : ",fPixArray->GetLast()+1));
  //StdoutToAliDebug(2,fPixArray->Print("","full"););

  // remove pixels with low signal or low visibility
  // Cuts are empirical !!!
  Double_t thresh = TMath::Max (mlem->GetMaximum()/100.,1.);
  thresh = TMath::Min (thresh,50.);
  Double_t charge = 0;

  // Mark pixels which should be removed
  for (Int_t i = 0; i < nPix; ++i) 
  {
    AliMUONPad* pixPtr = Pixel(i);
    charge = pixPtr->Charge();
    if (charge < thresh) 
    {
      pixPtr->SetCharge(-charge);
    }
  }

  // Move charge of removed pixels to their nearest neighbour (to keep total charge the same)
  Int_t near = 0;
  for (Int_t i = 0; i < nPix; ++i) 
  {
    AliMUONPad* pixPtr = Pixel(i);
    charge = pixPtr->Charge();
    if (charge > 0) continue;
    near = FindNearest(pixPtr);
    pixPtr->SetCharge(0);
    probi[i] = 0; // make it "invisible"
    AliMUONPad* pnear = Pixel(near);
    pnear->SetCharge(pnear->Charge() + (-charge));
  }
  Mlem(cluster,coef,probi,2);
  
  AliDebug(2,Form("Before splitting nPix=%d EVT %d DE %d",fPixArray->GetLast()+1,fEventNumber,fDetElemId));
  //StdoutToAliDebug(2,fPixArray->Print("","full"););
  //Plot("mlem.beforesplit");
  
           // Update histogram
  for (Int_t i = 0; i < nPix; ++i) 
  {
    AliMUONPad* pixPtr = Pixel(i);
    Int_t ix = mlem->GetXaxis()->FindBin(pixPtr->Coord(0));
    Int_t iy = mlem->GetYaxis()->FindBin(pixPtr->Coord(1));
    mlem->SetBinContent(ix, iy, pixPtr->Charge());
  }

  // Try to split into clusters
  Bool_t ok = kTRUE;
  if (mlem->GetSum() < 1) 
  {
    ok = kFALSE;
  }
  else 
  {
    fSplitter->Split(cluster,mlem,coef,fClusterList);
  }
  
  delete [] coef; 
  delete [] probi; 
  fPixArray->Delete(); 
  
  return ok;
}

//_____________________________________________________________________________
void AliMUONClusterFinderMLEM::Mlem(AliMUONCluster& cluster, 
                                    Double_t* coef, Double_t* probi, 
                                    Int_t nIter)
{
  /// Use MLEM to find pixel charges
  
  Int_t nPix = fPixArray->GetEntriesFast();

  Int_t npad = cluster.Multiplicity();

  Double_t* probi1 = new Double_t[nPix];
  Double_t probMax = TMath::MaxElement(nPix,probi);
  
  for (Int_t iter = 0; iter < nIter; ++iter) 
  {
    // Do iterations
    for (Int_t ipix = 0; ipix < nPix; ++ipix) 
    {
      Pixel(ipix)->SetChargeBackup(0);
      // Correct each pixel
      if (probi[ipix] < 0.01) continue; // skip "invisible" pixel
      Double_t sum = 0;
      probi1[ipix] = probMax;
      for (Int_t j = 0; j < npad; ++j) 
      {
        AliMUONPad* pad = cluster.Pad(j);
        if (pad->Status() < 0) continue; 
        Double_t sum1 = 0;
        Int_t indx1 = j*nPix;
        Int_t indx = indx1 + ipix;
        // Calculate expectation
        for (Int_t i = 0; i < nPix; ++i) 
        {
          sum1 += Pixel(i)->Charge()*coef[indx1+i];
	  //cout << i << " " << Pixel(i)->Charge() << " " << coef[indx1+i] << endl;
        } 
        if ( pad->IsSaturated() && sum1 > pad->Charge() ) 
        { 
          // correct for pad charge overflows
          probi1[ipix] -= coef[indx]; 
          continue; 
	  //sum1 = pad->Charge();
        } 

        if (sum1 > 1.e-6) sum += pad->Charge()*coef[indx]/sum1;
        //if (coef[indx] > 1.e-6) sum += pad->Charge()*coef[indx]/sum1;
      } // for (Int_t j=0;
      AliMUONPad* pixPtr = Pixel(ipix);
      if (probi1[ipix] > 1.e-6) 
      {
        //AZ pixPtr->SetCharge(pixPtr->Charge()*sum/probi1[ipix]);
        pixPtr->SetChargeBackup(pixPtr->Charge()*sum/probi1[ipix]);
      }
      //cout << " xxx " << ipix << " " << pixPtr->Charge() << " " << pixPtr->ChargeBackup() << " " << sum << " " << probi1[ipix] << endl;
    } // for (Int_t ipix=0;
    Double_t qTot = 0;
    for (Int_t i = 0; i < nPix; ++i) {
      AliMUONPad* pixPtr = Pixel(i);
      pixPtr->RevertCharge();
      qTot += pixPtr->Charge();
    }
    if (qTot < 1.e-6) {
      // Can happen in clusters with large number of overflows - speeding up 
      delete [] probi1;
      return;
    }
  } // for (Int_t iter=0;
  delete [] probi1;
}

//_____________________________________________________________________________
void AliMUONClusterFinderMLEM::FindCOG(TH2D *mlem, Double_t *xyc)
{
  /// Calculate position of the center-of-gravity around the maximum pixel

  Int_t ixmax, iymax, ix, nsumx=0, nsumy=0, nsum=0;
  Int_t i1 = -9, j1 = -9;
  mlem->GetMaximumBin(ixmax,iymax,ix);
  Int_t nx = mlem->GetNbinsX();
  Int_t ny = mlem->GetNbinsY();
  Double_t thresh = mlem->GetMaximum()/10;
  Double_t x, y, cont, xq=0, yq=0, qq=0;
  
  Int_t ie = TMath::Min(ny,iymax+1), je = TMath::Min(nx,ixmax+1);
  for (Int_t i = TMath::Max(1,iymax-1); i <= ie; ++i) {
    y = mlem->GetYaxis()->GetBinCenter(i);
    for (Int_t j = TMath::Max(1,ixmax-1); j <= je; ++j) {
      cont = mlem->GetCellContent(j,i);
      if (cont < thresh) continue;
      if (i != i1) {i1 = i; nsumy++;}
      if (j != j1) {j1 = j; nsumx++;}
      x = mlem->GetXaxis()->GetBinCenter(j);
      xq += x*cont;
      yq += y*cont;
      qq += cont;
      nsum++;
    }
  }
  
  Double_t cmax = 0;
  Int_t i2 = 0, j2 = 0;
  x = y = 0;
  if (nsumy == 1) {
    // one bin in Y - add one more (with the largest signal)
    for (Int_t i = TMath::Max(1,iymax-1); i <= ie; ++i) {
      if (i == iymax) continue;
      for (Int_t j = TMath::Max(1,ixmax-1); j <= je; ++j) {
        cont = mlem->GetCellContent(j,i);
        if (cont > cmax) {
          cmax = cont;
          x = mlem->GetXaxis()->GetBinCenter(j);
          y = mlem->GetYaxis()->GetBinCenter(i);
          i2 = i;
          j2 = j;
        }
      }
    }
    xq += x*cmax;
    yq += y*cmax;
    qq += cmax;
    if (i2 != i1) nsumy++;
    if (j2 != j1) nsumx++;
    nsum++;
  } // if (nsumy == 1)
  
  if (nsumx == 1) {
    // one bin in X - add one more (with the largest signal)
    cmax = x = y = 0;
    for (Int_t j = TMath::Max(1,ixmax-1); j <= je; ++j) {
      if (j == ixmax) continue;
      for (Int_t i = TMath::Max(1,iymax-1); i <= ie; ++i) {
        cont = mlem->GetCellContent(j,i);
        if (cont > cmax) {
          cmax = cont;
          x = mlem->GetXaxis()->GetBinCenter(j);
          y = mlem->GetYaxis()->GetBinCenter(i);
          i2 = i;
          j2 = j;
        }
      }
    }
    xq += x*cmax;
    yq += y*cmax;
    qq += cmax;
    if (i2 != i1) nsumy++;
    if (j2 != j1) nsumx++;
    nsum++;
  } // if (nsumx == 1)
  
  xyc[0] = xq/qq; xyc[1] = yq/qq;
  AliDebug(2,Form("x,y COG = %e,%e",xyc[0],xyc[1]));
}

//_____________________________________________________________________________
Int_t AliMUONClusterFinderMLEM::FindNearest(AliMUONPad *pixPtr0)
{
/// Find the pixel nearest to the given one
/// (algorithm may be not very efficient)

  Int_t nPix = fPixArray->GetEntriesFast(), imin = 0;
  Double_t rmin = 99999, dx = 0, dy = 0, r = 0;
  Double_t xc = pixPtr0->Coord(0), yc = pixPtr0->Coord(1);
  AliMUONPad *pixPtr;

  for (Int_t i = 0; i < nPix; ++i) {
    pixPtr = (AliMUONPad*) fPixArray->UncheckedAt(i);
    if (pixPtr == pixPtr0) continue;
    if (pixPtr->Charge() < 0.5) continue;
    dx = (xc - pixPtr->Coord(0)) / pixPtr->Size(0);
    dy = (yc - pixPtr->Coord(1)) / pixPtr->Size(1);
    r = dx *dx + dy * dy;
    if (r < rmin) { rmin = r; imin = i; }
  }
  return imin;
}

//_____________________________________________________________________________
void
AliMUONClusterFinderMLEM::Paint(Option_t*)
{
  /// Paint cluster and pixels
  
  AliMpArea area(fPreCluster->Area());
  
  gPad->Range(area.LeftBorder(),area.DownBorder(),area.RightBorder(),area.UpBorder());

  gVirtualX->SetFillStyle(1001);
  gVirtualX->SetFillColor(3);    
  gVirtualX->SetLineColor(3);
  
  Double_t s(1.0);
  
  for ( Int_t i = 0; i <= fPixArray->GetLast(); ++i)
  {
    AliMUONPad* pixel = Pixel(i);

    gPad->PaintBox(pixel->Coord(0)-pixel->Size(0)*s,
                   pixel->Coord(1)-pixel->Size(1)*s,
                   pixel->Coord(0)+pixel->Size(0)*s,
                   pixel->Coord(1)+pixel->Size(1)*s);

//    for ( Int_t sign = -1; sign < 2; sign +=2 )
//    {
//      gPad->PaintLine(pixel->Coord(0) - pixel->Size(0),
//                      pixel->Coord(1) + sign*pixel->Size(1),
//                      pixel->Coord(0) + pixel->Size(0),
//                      pixel->Coord(1) - sign*pixel->Size(1)
//                    );
//    }
  }      


  gVirtualX->SetFillStyle(0);
  
  fPreCluster->Paint();

  gVirtualX->SetLineColor(1);
  gVirtualX->SetLineWidth(2);
  gVirtualX->SetFillStyle(0);
  gVirtualX->SetTextColor(1);
  gVirtualX->SetTextAlign(22);
  
  for ( Int_t i = 0; i <= fPixArray->GetLast(); ++i)
  {
    AliMUONPad* pixel = Pixel(i);
    gPad->PaintBox(pixel->Coord(0)-pixel->Size(0),
                   pixel->Coord(1)-pixel->Size(1),
                   pixel->Coord(0)+pixel->Size(0),
                   pixel->Coord(1)+pixel->Size(1));    
    gVirtualX->SetTextSize(pixel->Size(0)*60);

    gPad->PaintText(pixel->Coord(0),pixel->Coord(1),Form("%d",(Int_t)(pixel->Charge())));
  }  
}

//_____________________________________________________________________________
Int_t AliMUONClusterFinderMLEM::FindLocalMaxima(TObjArray *pixArray, Int_t *localMax, Double_t *maxVal)
{
/// Find local maxima in pixel space for large preclusters in order to
/// try to split them into smaller pieces (to speed up the MLEM procedure)
/// or to find additional fitting seeds if clusters were not completely resolved  

  AliDebug(1,Form("nPix=%d",pixArray->GetLast()+1));

  TH2D *hist = NULL;
  //if (pixArray == fPixArray) hist = (TH2D*) gROOT->FindObject("anode");
  //else { hist = (TH2D*) gROOT->FindObject("anode1"); cout << hist << endl; }
  //if (hist) hist->Delete();
 
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
  if (pixArray == fPixArray) hist = new TH2D("anode","anode",nx,xylim[0],-xylim[1],ny,xylim[2],-xylim[3]);
  else hist = new TH2D("anode1","anode1",nx,xylim[0],-xylim[1],ny,xylim[2],-xylim[3]);
  for (Int_t ipix = 0; ipix < nPix; ++ipix) {
    pixPtr = (AliMUONPad*) pixArray->UncheckedAt(ipix);
    hist->Fill(pixPtr->Coord(0), pixPtr->Coord(1), pixPtr->Charge());
  }
//  if (fDraw && pixArray == fPixArray) fDraw->DrawHist("c2", hist);

  Int_t nMax = 0, indx, nxy = ny * nx;
  Int_t *isLocalMax = new Int_t[nxy];
  for (Int_t i = 0; i < nxy; ++i) isLocalMax[i] = 0; 

  for (Int_t i = 1; i <= ny; ++i) {
    indx = (i-1) * nx;
    for (Int_t j = 1; j <= nx; ++j) {
      if (hist->GetCellContent(j,i) < 0.5) continue;
      //if (isLocalMax[indx+j-1] < 0) continue;
      if (isLocalMax[indx+j-1] != 0) continue;
      FlagLocalMax(hist, i, j, isLocalMax);
    }
  }

  for (Int_t i = 1; i <= ny; ++i) {
    indx = (i-1) * nx;
    for (Int_t j = 1; j <= nx; ++j) {
      if (isLocalMax[indx+j-1] > 0) { 
	localMax[nMax] = indx + j - 1; 
	maxVal[nMax++] = hist->GetCellContent(j,i);
	if (nMax > 99) AliFatal(" Too many local maxima !!!");
      }
    }
  }
  if (fDebug) cout << " Local max: " << nMax << endl;
  delete [] isLocalMax; 
  if (nMax == 1) hist->Delete();
  return nMax;
}

//_____________________________________________________________________________
void AliMUONClusterFinderMLEM::FlagLocalMax(TH2D *hist, Int_t i, Int_t j, Int_t *isLocalMax)
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
void AliMUONClusterFinderMLEM::FindCluster(AliMUONCluster& cluster, 
                                           Int_t *localMax, Int_t iMax)
{
/// Find pixel cluster around local maximum \a iMax and pick up pads
/// overlapping with it

  TH2D *hist = (TH2D*) gROOT->FindObject("anode");
  /* Just for check
  TCanvas* c = new TCanvas("Anode","Anode",800,600);
  c->cd();
  hist->Draw("lego1Fb"); // debug
  c->Update();
  Int_t tmp;
  cin >> tmp;
  */
  Int_t nx = hist->GetNbinsX();
  Int_t ny = hist->GetNbinsY();
  Int_t ic = localMax[iMax] / nx + 1;
  Int_t jc = localMax[iMax] % nx + 1;
  Int_t nxy = ny * nx;
  Bool_t *used = new Bool_t[nxy];
  for (Int_t i = 0; i < nxy; ++i) used[i] = kFALSE; 

  // Drop all pixels from the array - pick up only the ones from the cluster
  fPixArray->Delete();

  Double_t wx = hist->GetXaxis()->GetBinWidth(1)/2; 
  Double_t wy = hist->GetYaxis()->GetBinWidth(1)/2;  
  Double_t yc = hist->GetYaxis()->GetBinCenter(ic);
  Double_t xc = hist->GetXaxis()->GetBinCenter(jc);
  Double_t cont = hist->GetCellContent(jc,ic);
  fPixArray->Add(new AliMUONPad (xc, yc, wx, wy, cont));
  used[(ic-1)*nx+jc-1] = kTRUE;
  AddBinSimple(hist, ic, jc);
  //fSplitter->AddBin(hist, ic, jc, 1, used, (TObjArray*)0); // recursive call

  Int_t nPix = fPixArray->GetEntriesFast();
  Int_t npad = cluster.Multiplicity();
  
  for (Int_t i = 0; i < nPix; ++i) 
  {
    AliMUONPad* pixPtr = Pixel(i);
    pixPtr->SetSize(0,wx);
    pixPtr->SetSize(1,wy);
  }

  // Pick up pads which overlap with found pixels
  for (Int_t i = 0; i < npad; ++i) 
  {
    cluster.Pad(i)->SetStatus(-1);
  }
  
  for (Int_t i = 0; i < nPix; ++i) 
  {
    AliMUONPad* pixPtr = Pixel(i);
    for (Int_t j = 0; j < npad; ++j) 
    {
      AliMUONPad* pad = cluster.Pad(j);
      if (pad->Status() == 0) continue;
      if ( Overlap(*pad,*pixPtr) )
      {
        pad->SetStatus(0);
	if (fDebug) { cout << j << " "; pad->Print("full"); }
      }
    }
  }

  delete [] used;
}

//_____________________________________________________________________________
void 
AliMUONClusterFinderMLEM::AddBinSimple(TH2D *mlem, Int_t ic, Int_t jc)
{
  /// Add adjacent bins (+-1 in X and Y) to the cluster
  
  Int_t nx = mlem->GetNbinsX();
  Int_t ny = mlem->GetNbinsY();
  Double_t cont1, cont = mlem->GetCellContent(jc,ic);
  AliMUONPad *pixPtr = 0;
  
  Int_t ie = TMath::Min(ic+1,ny), je = TMath::Min(jc+1,nx);
  for (Int_t i = TMath::Max(ic-1,1); i <= ie; ++i) {
    for (Int_t j = TMath::Max(jc-1,1); j <= je; ++j) {
      cont1 = mlem->GetCellContent(j,i);
      if (cont1 > cont) continue;
      if (cont1 < 0.5) continue;
      pixPtr = new AliMUONPad (mlem->GetXaxis()->GetBinCenter(j), 
			       mlem->GetYaxis()->GetBinCenter(i), 0, 0, cont1);
      fPixArray->Add(pixPtr);
    }
  }
}

//_____________________________________________________________________________
AliMUONClusterFinderMLEM&  
AliMUONClusterFinderMLEM::operator=(const AliMUONClusterFinderMLEM& rhs)
{
/// Protected assignement operator

  if (this == &rhs) return *this;

  AliFatal("Not implemented.");
    
  return *this;  
}    

//_____________________________________________________________________________
void
AliMUONClusterFinderMLEM::Neighbours(Int_t cathode, Int_t ix, Int_t iy,
                                     Int_t& n, Int_t* xList, Int_t* yList)
{
  /// Get the list of neighbours of pad at (cathode,ix,iy)
  n = 0;
  
  const AliMpVSegmentation* seg = fSegmentation[cathode];
  
  AliMpPad pad = seg->PadByIndices(AliMpIntPair(ix,iy),kTRUE);
	
  // Define the region to look into : a region slightly bigger
  // than the pad itself (5% bigger), in order to catch first neighbours.
  
  AliMpArea area(pad.Position(),pad.Dimensions()*1.05); 
		
  AliMpVPadIterator* it = seg->CreateIterator(area);
  it->First();
  while ( !it->IsDone() && n < 10 )
	{
		AliMpPad p = it->CurrentItem();
		if ( p != pad ) // skip self
		{
			xList[n] = p.GetIndices().GetFirst();
			yList[n] = p.GetIndices().GetSecond();
			++n;
		}
		it->Next();
	}
	delete it;
}

//_____________________________________________________________________________
void AliMUONClusterFinderMLEM::AddVirtualPad(AliMUONCluster& cluster)
{
  /// Add virtual pad (with small charge) to improve fit for clusters
  /// with number of pads == 2 per direction
  
  // Find out non-bending and bending planes
  Int_t nonb[2] = {1, 0}; // non-bending and bending cathodes

  TVector2 dim0 = cluster.MinPadDimensions(0, 0, kTRUE);
  TVector2 dim1 = cluster.MinPadDimensions(1, 0, kTRUE);
  if (dim0.X() < dim1.X() - fgkDistancePrecision) {
    nonb[0] = 0;
    nonb[1] = 1;
  } 

  Bool_t same = kFALSE;
  if (TMath::Abs(dim0.Y()-dim1.Y()) < fgkDistancePrecision) same = kTRUE; // the same pad size on both planes 

  AliMpIntPair cn;
  Bool_t check[2] = {kFALSE, kFALSE};
  Int_t nxy[2];
  nxy[0] = nxy[1] = 0;
  for (Int_t inb = 0; inb < 2; ++inb) {
    cn = cluster.NofPads(nonb[inb], 0, kTRUE);
    if (inb == 0 && cn.GetFirst() == 2) check[inb] = kTRUE; // check non-bending plane
    else if (inb == 1 && cn.GetSecond() == 2) check[inb] = kTRUE; // check bending plane
    if (same) {
      nxy[0] = TMath::Max (nxy[0], cn.GetFirst());
      nxy[1] = TMath::Max (nxy[1], cn.GetSecond());
      if (inb == 0 && nxy[0] < 2) nonb[inb] = !nonb[inb];
      else if (inb == 1 && cn.GetSecond() < 2) nonb[inb] = !nonb[inb];
    }
  }
  if (same) {
    if (nxy[0] > 2) check[0] = kFALSE;
    if (nxy[1] > 2) check[1] = kFALSE;
  }
  if (!check[0] && !check[1]) return;

  for (Int_t inb = 0; inb < 2; ++inb) {
    if (!check[inb]) continue;

    // Find pads with maximum and next to maximum charges 
    Int_t maxPads[2] = {-1, -1};
    Double_t amax[2] = {0};
    Int_t mult = cluster.Multiplicity();
    for (Int_t j = 0; j < mult; ++j) {
      AliMUONPad *pad = cluster.Pad(j);
      if (pad->Cathode() != nonb[inb]) continue;
      for (Int_t j2 = 0; j2 < 2; ++j2) {
	if (pad->Charge() > amax[j2]) {
	  if (j2 == 0) { amax[1] = amax[0]; maxPads[1] = maxPads[0]; }
	  amax[j2] = pad->Charge();
	  maxPads[j2] = j;
	  break;
	}
      }
    }

    // Find min and max dimensions of the cluster
    Double_t limits[2] = {9999, -9999};
    for (Int_t j = 0; j < mult; ++j) {
      AliMUONPad *pad = cluster.Pad(j);
      if (pad->Cathode() != nonb[inb]) continue;
      if (pad->Coord(inb) < limits[0]) limits[0] = pad->Coord(inb);
      if (pad->Coord(inb) > limits[1]) limits[1] = pad->Coord(inb);
    }

    // Loop over max and next to max pads
    Bool_t add = kFALSE;
    Int_t idirAdd = 0;
    for (Int_t j = 0; j < 2; ++j) {
      if (j == 1) {
	if (maxPads[j] < 0) continue;
	if (!add) break; 
	if (amax[1] / amax[0] < 0.5) break;
      }
      // Check if pad at the cluster limit
      AliMUONPad *pad = cluster.Pad(maxPads[j]);
      Int_t idir = 0;
      if (TMath::Abs(pad->Coord(inb)-limits[0]) < fgkDistancePrecision) idir = -1;
      else if (TMath::Abs(pad->Coord(inb)-limits[1]) < fgkDistancePrecision) idir = 1;
      else {
	//cout << " *** Pad not at the cluster limit: " << j << endl;
	break;
      }
      if (j == 1 && idir == idirAdd) break; // this pad has been already added

      // Add pad (if it exists)
      TVector2 pos;
      if (inb == 0) pos.Set (pad->X() + idir * (pad->DX()+fgkDistancePrecision), pad->Y());
      else pos.Set (pad->X(), pad->Y() + idir * (pad->DY()+fgkDistancePrecision));
      //AliMpPad mppad = fSegmentation[nonb[inb]]->PadByPosition(pos,kTRUE);
      AliMpPad mppad = fSegmentation[nonb[inb]]->PadByPosition(pos,kFALSE);
      if (!mppad.IsValid()) continue; // non-existing pad
      cn = mppad.GetIndices();
      AliMUONPad muonPad(fDetElemId, nonb[inb], cn.GetFirst(), cn.GetSecond(), 
			 mppad.Position().X(), mppad.Position().Y(), 
			 mppad.Dimensions().X(), mppad.Dimensions().Y(), 0);
      if (inb == 0) muonPad.SetCharge(TMath::Min (amax[j]/100, 5.));
      else muonPad.SetCharge(TMath::Min (amax[j]/15, fgkZeroSuppression));
      if (muonPad.Charge() < 1.) muonPad.SetCharge(1.);
      muonPad.SetReal(kFALSE);
      if (fDebug) printf(" ***** Add virtual pad in %d direction ***** %f %f %f %3d %3d %f %f \n",
			 inb, muonPad.Charge(), muonPad.X(), muonPad.Y(), muonPad.Ix(), 
			 muonPad.Iy(), muonPad.DX(), muonPad.DY());
      cluster.AddPad(muonPad); // add pad to the cluster
      add = kTRUE;
      idirAdd = idir;
    }
  }
}

//_____________________________________________________________________________
void AliMUONClusterFinderMLEM::PadsInXandY(AliMUONCluster& cluster,
                                           Int_t &nInX, Int_t &nInY) const
{
  /// Find number of pads in X and Y-directions (excluding virtual ones and
  /// overflows)

  Int_t statusToTest = 1;
  
  if ( nInX < 0 ) statusToTest = 0;
       
  Bool_t mustMatch(kTRUE);

  AliMpIntPair cn = cluster.NofPads(statusToTest,mustMatch);
  
  nInX = cn.GetFirst();
  nInY = cn.GetSecond();
}

//_____________________________________________________________________________
void AliMUONClusterFinderMLEM::RemovePixel(Int_t i)
{
  /// Remove pixel at index i
  AliMUONPad* pixPtr = Pixel(i);
  fPixArray->RemoveAt(i); 
  delete pixPtr;
}

//_____________________________________________________________________________
void AliMUONClusterFinderMLEM::Simple(AliMUONCluster& cluster)
{
/// Process simple cluster (small number of pads) without EM-procedure

  Int_t nForFit = 1, clustFit[1] = {0}, nfit;
  Double_t parOk[3] = {0.}; 
  TObjArray *clusters[1]; 
  clusters[0] = fPixArray;

  AliDebug(1,Form("nPix=%d",fPixArray->GetLast()+1));

  Int_t mult = cluster.Multiplicity();
  for (Int_t i = 0; i < mult; ++i) 
  {
    AliMUONPad* pad = cluster.Pad(i);
    if ( pad->IsSaturated()) 
    {
      pad->SetStatus(-9);
    }
    else 
    {
      pad->SetStatus(1);
    }
  }
  nfit = fSplitter->Fit(cluster,1, nForFit, clustFit, clusters, parOk, fClusterList);
}

//_____________________________________________________________________________
AliMUONPad* 
AliMUONClusterFinderMLEM::Pixel(Int_t i) const
{
  /// Returns pixel at index i
  return static_cast<AliMUONPad*>(fPixArray->UncheckedAt(i));
}

//_____________________________________________________________________________
void 
AliMUONClusterFinderMLEM::Print(Option_t* what) const
{
  /// printout
  TString swhat(what);
  swhat.ToLower();
  if ( swhat.Contains("precluster") )
  {
    if ( fPreCluster) fPreCluster->Print();
  }
}


