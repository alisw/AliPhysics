//-*- Mode: C++ -*-
// **************************************************************************
// This file is property of and copyright by the ALICE ITSU Project         *
// ALICE Experiment at CERN, All rights reserved.                           *
//                                                                          *
// Primary Author: Ruben Shahoyan                                           *
//                                                                          *
// Adapted to ITSU: Maximiliano Puccio <maximiliano.puccio@cern.ch>         *
//                  for the ITS Upgrade project                             *
//                                                                          *
// Permission to use, copy, modify and distribute this software and its     *
// documentation strictly for non-commercial purposes is hereby granted     *
// without fee, provided that the above copyright notice appears in all     *
// copies and that both the copyright notice and this permission notice     *
// appear in the supporting documentation. The authors make no claims       *
// about the suitability of this software for any purpose. It is            *
// provided "as is" without express or implied warranty.                    *
//                                                                          *
//***************************************************************************

#include "AliITSUCATrackingStation.h"
#include <TMath.h>
#include "AliITSUGeomTGeo.h"
#include "AliVertex.h"
#include <TRandom.h>
#include <TStopwatch.h>
#include <TString.h>
#include "AliITSMFTAux.h"
#include "AliITSURecoSens.h"
#include <Riostream.h>

using AliITSMFTAux::BringTo02Pi;

//_________________________________________________________________
AliITSUCATrackingStation::AliITSUCATrackingStation()
:fClusters(0x0)
,fVIDOffset(0)
,fNClusters(0)
,fZMin(0)
,fZMax(0)
,fDZInv(-1)
,fDPhiInv(-1)
,fNZBins(20)
,fNPhiBins(20)
,fQueryZBmin(-1)
,fQueryZBmax(-1)
,fQueryPhiBmin(-1)
,fQueryPhiBmax(-1)
,fBins(0)
,fOccBins(0)
,fNOccBins(0)
,fNFoundClusters(0)
,fFoundClusterIterator(0)
,fFoundBinIterator(0)
,fIndex()
,fFoundBins(0)
,fSortedClInfo(0)
,fDetectors()
{
  // def. c-tor
}

//_________________________________________________________________
AliITSUCATrackingStation::AliITSUCATrackingStation(int nzbins,int nphibins,
                                                   AliITSURecoLayer *lr, AliITSUGeomTGeo *geo)
:fClusters(0x0)
,fVIDOffset()
,fNClusters(0)
,fZMin(0.f)
,fZMax(0.f)
,fDZInv(-1)
,fDPhiInv(-1)
,fNZBins(nzbins)
,fNPhiBins(nphibins)
,fQueryZBmin(-1)
,fQueryZBmax(-1)
,fQueryPhiBmin(-1)
,fQueryPhiBmax(-1)
,fBins(0)
,fOccBins(0)
,fNOccBins(0)
,fNFoundClusters(0)
,fFoundClusterIterator(0)
,fFoundBinIterator(0)
,fIndex()
,fFoundBins(0)
,fSortedClInfo(0)
,fDetectors()
{
  // c-tor
  Init(lr,geo);
}

//_________________________________________________________________
AliITSUCATrackingStation::~AliITSUCATrackingStation()
{
  // d-tor
  delete[] fBins;
  delete[] fOccBins;
//  delete fClusters;
}

//_________________________________________________________________
void AliITSUCATrackingStation::Init(AliITSURecoLayer *lr, AliITSUGeomTGeo *geo)
{
  if (fClusters) {
    printf("Already initialized\n");
    return;
  }
  
  fClusters = *lr->GetClustersAddress();
  fZMin = lr->GetZMin();
  fZMax = lr->GetZMax();
  if (fNZBins < 1)   fNZBins = 2;
  if (fNPhiBins < 1) fNPhiBins = 1;
  fDZInv   = fNZBins / (fZMax - fZMin);
  fDPhiInv = fNPhiBins / TMath::TwoPi();
  //
  fBins = new ClBinInfo_t[fNZBins * fNPhiBins];
  fOccBins = new int[fNZBins * fNPhiBins];
  fNClusters = fClusters->GetEntriesFast();
  if(fNClusters == 0) return;
  fSortedClInfo.reserve(fClusters->GetEntriesFast());
  fVIDOffset = ((AliITSUClusterPix*)fClusters->UncheckedAt(0))->GetVolumeId();
  //
  // prepare detectors info
  int detID = -1;
  fIndex.resize(lr->GetNSensors(),-1);
  fDetectors.reserve(lr->GetNSensors());
  for (int iCl = 0; iCl < fClusters->GetEntriesFast(); ++iCl)
  { //Fill this layer with detectors
    AliITSUClusterPix* c = (AliITSUClusterPix*)fClusters->UncheckedAt(iCl);
    if (detID == c->GetVolumeId())
    {
      continue;
    }
    detID = c->GetVolumeId();
    ITSDetInfo_t det;
    det.index = iCl;
    //
    TGeoHMatrix m;
    geo->GetOrigMatrix(detID,m);
    //
    fIndex[detID - fVIDOffset] = fDetectors.size();
    AliITSURecoSens* sens = lr->GetSensorFromID(detID);
    const TGeoHMatrix *tm = geo->GetMatrixT2L(detID);
    m.Multiply(tm);
    double txyz[3] = {0.,0.,0.}, xyz[3] = {0.,0.,0.};
    m.LocalToMaster(txyz,xyz);
    det.xTF = sens->GetXTF(); // TMath::Sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1]);

    det.phiTF = sens->GetPhiTF();//TMath::ATan2(xyz[1],xyz[0]);
    det.sinTF = TMath::Sin(det.phiTF);
    det.cosTF = TMath::Cos(det.phiTF);
    //
    // compute the real radius (with misalignment)
    TGeoHMatrix mmisal(*(geo->GetMatrix(detID)));
    mmisal.Multiply(tm);
    xyz[0] = 0.;
    xyz[1] = 0.;
    xyz[2] = 0.;
    mmisal.LocalToMaster(txyz,xyz);
    det.xTFmisal = TMath::Sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1]);
//    if (!TMath::AreEqualAbs(det.xTFmisal, sens->GetXTF(), 1e-9)) {
//      cout << "Error: " << lr->GetActiveID() << " " << detID << " " << det.xTFmisal - sens->GetXTF()<< endl;
//    }
    
    fDetectors.push_back(det);
  } // end loop on detectors
}

//_________________________________________________________________
void AliITSUCATrackingStation::SortClusters(const AliVertex* vtx)
{
  // sort clusters and build fast lookup table
  //
  ClearSortedInfo();
  fSortedClInfo.reserve(fNClusters);
  //
  ClsInfo_t cl;
  for (int icl = fNClusters;icl--;) {
    AliITSUClusterPix* cluster = (AliITSUClusterPix*)fClusters->UncheckedAt(icl);
    cluster->GetGlobalXYZ( (float*)&cl );
    //
    if (vtx) { // phi and r will be computed wrt vertex
      cl.x -= vtx->GetX();
      cl.y -= vtx->GetY();
    }
    //
    cl.r = TMath::Sqrt(cl.x*cl.x + cl.y*cl.y);
    cl.phi = TMath::ATan2(cl.y,cl.x);
    BringTo02Pi(cl.phi);
    cl.index = icl;
    cl.zphibin = GetBinIndex(GetZBin(cl.z),GetPhiBin(cl.phi));
    cl.detid = cluster->GetVolumeId() - fVIDOffset;
    //
    fSortedClInfo.push_back( cl );
    //
  }
  sort(fSortedClInfo.begin(), fSortedClInfo.end()); // sort in phi, z
  //
  // fill cells in phi,z
  int currBin = -1;
  for (int icl = 0;icl < fNClusters; ++icl)
  {
    ClsInfo_t &t = fSortedClInfo[icl];
    if (t.zphibin > currBin)
    { // register new occupied bin
      currBin = t.zphibin;
      fBins[currBin].first = icl;
      fBins[currBin].index = fNOccBins;
      fOccBins[fNOccBins++] = currBin;
    }
    fBins[currBin].ncl++;
  }
  //  Print("clb"); //RS
}

//_________________________________________________________________
void AliITSUCATrackingStation::Clear(Option_t *)
{
  // clear cluster info
  ClearSortedInfo();
  fIndex.clear();
  fNClusters = 0;
  if (fClusters) fClusters = 0x0;
  //
}

//_________________________________________________________________
void AliITSUCATrackingStation::ClearSortedInfo()
{
  // clear cluster info
  fSortedClInfo.clear();
  memset(fBins,0,fNZBins * fNPhiBins * sizeof(ClBinInfo_t));
  memset(fOccBins,0,fNZBins * fNPhiBins * sizeof(int));
  fNOccBins = 0;
}

//_________________________________________________________________
void AliITSUCATrackingStation::Print(Option_t *opt) const
{
  // dump cluster bins info
  TString opts = opt;
  opts.ToLower();
  printf("Stored %d clusters in %d occupied bins\n",fNClusters,fNOccBins);
  //
  if (opts.Contains("c")) {
    printf("\nCluster info\n");
    for (int i = 0; i < fNClusters;i++)
    {
      const ClsInfo_t &t = fSortedClInfo[i];
      printf("#%5d Bin(phi/z):%03d/%03d Z:%+8.3f Phi:%+6.3f R:%7.3f Ind:%d ",
             i,t.zphibin/fNZBins,t.zphibin%fNZBins,t.z,t.phi,t.r,t.index);
      if (opts.Contains("l")) { // mc labels
        AliITSUClusterPix* rp = (AliITSUClusterPix*)fClusters->UncheckedAt(t.index);
        for (int l = 0;l < 3; l++) if (rp->GetLabel(l) >= 0) printf("| %d ",rp->GetLabel(l));
      }
      printf("\n");
    }
  }
  //
  if (opts.Contains("b")) {
    printf("\nBins info (occupied only)\n");
    for (int i=0;i<fNOccBins;i++) {
      printf("%4d %5d(phi/z: %03d/%03d) -> %3d cl from %d\n",i,fOccBins[i],fOccBins[i]/fNZBins,fOccBins[i]%fNZBins,
             fBins[fOccBins[i]].ncl,fBins[fOccBins[i]].first);
    }
  }
  //
}

//_____________________________________________________________
int AliITSUCATrackingStation::SelectClusters(float zmin,float zmax,float phimin,float phimax)
{
  // prepare occupied bins in the requested region
  //printf("Select: Z %f %f | Phi: %f %f\n",zmin,zmax,phimin,phimax);
  if (!fNOccBins) return 0;
  if (zmax < fZMin || zmin > fZMax || zmin > zmax) return 0;
  fFoundBins.clear();
  
  fQueryZBmin = GetZBin(zmin);
  if (fQueryZBmin < 0) fQueryZBmin = 0;
  fQueryZBmax = GetZBin(zmax);
  if (fQueryZBmax >= fNZBins) fQueryZBmax = fNZBins - 1;
  BringTo02Pi(phimin);
  BringTo02Pi(phimax);
  fQueryPhiBmin = GetPhiBin(phimin);
  fQueryPhiBmax = GetPhiBin(phimax);
  int dbz = 0;
  fNFoundClusters = 0;
  int nbcheck = fQueryPhiBmax - fQueryPhiBmin + 1; //TODO:(MP) check if a circular buffer is feasible
  if (nbcheck > 0)
  { // no wrapping around 0-2pi, fast case
    for (int ip = fQueryPhiBmin;ip <= fQueryPhiBmax;ip++) {
      int binID = GetBinIndex(fQueryZBmin,ip);
      if ( !(dbz = (fQueryZBmax-fQueryZBmin)) ) { // just one Z bin in the query range
        ClBinInfo_t& binInfo = fBins[binID];
        if (!binInfo.ncl) continue;
        fNFoundClusters += binInfo.ncl;
        fFoundBins.push_back(binID);
        continue;
      }
      int binMax = binID + dbz;
      for ( ; binID <= binMax; binID++) {
        ClBinInfo_t& binInfo = fBins[binID];
        if (!binInfo.ncl) continue;
        fNFoundClusters += binInfo.ncl;
        fFoundBins.push_back(binID);
      }
    }
  }
  else
  {  // wrapping
    nbcheck += fNPhiBins;
    for (int ip0 = 0;ip0 <= nbcheck;ip0++) {
      int ip = fQueryPhiBmin + ip0;
      if (ip >= fNPhiBins) ip -= fNPhiBins;
      int binID = GetBinIndex(fQueryZBmin,ip);
      if ( !(dbz = (fQueryZBmax - fQueryZBmin)) ) { // just one Z bin in the query range
        ClBinInfo_t& binInfo = fBins[binID];
        if (!binInfo.ncl) continue;
        fNFoundClusters += binInfo.ncl;
        fFoundBins.push_back(binID);
        continue;
      }
      int binMax = binID + dbz;
      for (;binID <= binMax;binID++) {
        ClBinInfo_t& binInfo = fBins[binID];
        if (!binInfo.ncl) continue;
        fNFoundClusters += binInfo.ncl;
        fFoundBins.push_back(binID);
      }
    }
  }
  fFoundClusterIterator = fFoundBinIterator = 0;
  /*
   //printf("Selected -> %d cl in %d bins\n",fNFoundClusters,(int)fFoundBins.size());
   for (int i=0;i<(int)fFoundBins.size();i++) {
   int bn = fFoundBins[i];
   ClBinInfo_t& bin=fBins[bn];
   printf("#%d b:%d 1st: %3d Ncl:%d\n",i,bn,bin.first,bin.ncl);
   }
   printf("\n");
   */
  return fNFoundClusters;
}

//_____________________________________________________________
int AliITSUCATrackingStation::GetNextClusterInfoID()
{
  if (fFoundBinIterator < 0) return 0;
  int currBin = fFoundBins[fFoundBinIterator];
  if (fFoundClusterIterator < fBins[currBin].ncl) { // same bin
    return fBins[currBin].first + fFoundClusterIterator++;
  }
  if (++fFoundBinIterator < int(fFoundBins.size())) {  // need to change bin
    currBin = fFoundBins[fFoundBinIterator];
    fFoundClusterIterator = 1;
    return fBins[currBin].first;
  }
  fFoundBinIterator = -1;
  return -1;
}

//_____________________________________________________________
void AliITSUCATrackingStation::ResetFoundIterator()
{
  // prepare for a new loop over found clusters
  if (fNFoundClusters)  fFoundClusterIterator = fFoundBinIterator = 0;
}
