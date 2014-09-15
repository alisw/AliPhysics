//-*- Mode: C++ -*-
// **************************************************************************
// This file is property of and copyright by the ALICE ITSU Project         *
// ALICE Experiment at CERN, All rights reserved.                           *
//                                                                          *
// Primary Author: Maximiliano Puccio <maximiliano.puccio@cern.ch>          *
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

#include <algorithm>
#include <TBranch.h>
#include <TMath.h>
#include <TTree.h>
#include <Riostream.h>
#include "AliLog.h"
#include "AliESDEvent.h"
#include "AliITSUClusterPix.h"
#include "AliITSUCATracker.h"
#include "AliITSUReconstructor.h"
#include "AliITSURecoDet.h"
#include "AliESDtrack.h"
#include <cassert>

using TMath::Abs;
using TMath::Sort;
using TMath::Sqrt;
using std::sort;
using std::cout;
using std::endl;
using std::flush;

ClassImp(AliITSUCATracker)


// tolerance for layer on-surface check
const Double_t AliITSUCATracker::fgkToler =  1e-6;
const Double_t AliITSUCATracker::fgkChi2Cut =  600.f;
const int AliITSUCATracker::fgkNumberOfIterations =  1;
const float AliITSUCATracker::fgkR[7] = {2.33959,3.14076,3.91924,19.6213,24.5597,34.388,39.3329};
//
const float kmaxDCAxy[5] = {0.05f,0.04f,0.05f,0.2f,0.4f};
const float kmaxDCAz[5] = {0.2f,0.4f,0.5f,0.6f,3.f};
const float kmaxDN[4] = {0.002f,0.009f,0.002f,0.005f};
const float kmaxDP[4] = {0.008f,0.0025f,0.003f,0.0035f};
const float kmaxDZ[6] = {0.1f,0.1f,0.3f,0.3f,0.3f,0.3f};
//const float kSigma2 = 0.0005f * 0.0005f;

//__________________________________________________________________________________________________
static inline float invsqrt(float _x  )
{
  //
  // The function calculates fast inverse sqrt. Google for 0x5f3759df.
  // Credits to ALICE HLT Project
  //

  union { float f; int i; } x = { _x };
  const float xhalf = 0.5f * x.f;
  x.i = 0x5f3759df - ( x.i >> 1 );
  x.f = x.f * ( 1.5f - xhalf * x.f * x.f );
  return x.f;
}

//__________________________________________________________________________________________________
static inline float Curvature(float x1, float y1, float x2, float y2, float x3, float y3)
{
  //
  // Initial approximation of the track curvature
  //
//  return - 2.f * ((x2 - x1) * (y3 - y2) - (x3 - x2) * (y2 - y1))
//               * invsqrt(((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1)) *
//                         ((x2 - x3) * (x2 - x3) + (y2 - y3) * (y2 - y3)) *
//                         ((x3 - x1) * (x3 - x1) + (y3 - y1) * (y3 - y1)));
  
  //calculates the curvature of track
  float den = (x3 - x1) * (y2 - y1) - (x2 - x1) * (y3 - y1);
  if(den * den < 1e-32) return 0.f;
  float a = ((y3-y1)*(x2*x2+y2*y2-x1*x1-y1*y1)-(y2-y1)*(x3*x3+y3*y3-x1*x1-y1*y1))/den;
  if((y2 - y1) * (y2 - y1) < 1e-32) return 0.f;
  float b = -(x2 * x2 - x1 * x1 + y2 * y2 - y1 * y1 + a * (x2 - x1)) / (y2 - y1);
  float c = -x1 * x1 - y1 * y1 - a * x1 - b * y1;
  float xc = -a / 2.f;
  
  if((a * a + b * b - 4 * c) < 0) return 0.f;
  float rad = TMath::Sqrt(a * a + b * b - 4 * c) / 2.f;
  if(rad * rad < 1e-32) return 1e16;
  
  if((x1 > 0.f && y1 > 0.f && x1 < xc)) rad *= -1.f;
  if((x1 < 0.f && y1 > 0.f && x1 < xc)) rad *= -1.f;
  //  if((x1<0 && y1<0 && x1<xc)) rad*=-1;
  // if((x1>0 && y1<0 && x1<xc)) rad*=-1;
  
  return 1.f/rad;

}

//__________________________________________________________________________________________________
static inline float TanLambda(float x1, float y1, float x2, float y2, float z1, float z2)
{
  //
  // Initial approximation of the tangent of the track dip angle
  //
  return (z1 - z2) * invsqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
}

//__________________________________________________________________________________________________
//static inline float XCenterOfCurvature(float x1, float y1, float x2, float y2, float x3, float y3)
//{
//    //
//    // Initial approximation of the x-coordinate of the center of curvature
//    //
//
//  const float k1 = (y2 - y1) / (x2 - x1), k2 = (y3 - y2) / (x3 - x2);
//  return TMath::Abs(k2 - k1) > kAlmost0 ?
//    0.5f * (k1 * k2 * (y1 - y3) + k2 * (x1 + x2) - k1 * (x2 + x3)) / (k2 - k1) : 1e12f;
//}

//__________________________________________________________________________________________________
static inline bool CompareAngles(float alpha, float beta, float tolerance)
{
	const float delta = TMath::Abs(alpha - beta);
	return (delta < tolerance || TMath::Abs(delta - TMath::TwoPi()) < tolerance);
}

//__________________________________________________________________________________________________
AliITSUCATracker::AliITSUCATracker(AliITSUReconstructor* rec) :
  fReconstructor(rec),
  fITS(0),
  fMatLUT(0),
  fUseMatLUT(kFALSE),
  fCurrMass(0.14),
  fLayer(),
  fUsedClusters(),
  fChi2Cut(fgkChi2Cut),
  fPhiCut(1),
  fZCut(0.5f),
  fCandidates()
{
  // This default constructor needs to be provided
  if (rec) Init(rec);
  for (int i = 0; i < 4; ++i)
  {
    fCandidates[i] = new TClonesArray("AliITSUTrackCooked",100000);
  }

#ifdef _TUNING_
  for (int i = 0; i < 6; ++i)
  {
    fGDZ[i] = new TH1F(Form("DGZ%i",i),Form("DZ%i;#Deltaz;Entries",i),500,0.f,5.f);
    fGDXY[i] = new TH1F(Form("DGPhi%i",i),Form("#Delta#Phi%i;#DeltaPhi;Entries",i),500,0.f,TMath::TwoPi());
    fFDZ[i] = new TH1F(Form("DFZ%i",i),Form("DZ%i;#Deltaz;Entries",i),500,0.f,5.f);
    fFDXY[i] = new TH1F(Form("DFPhi%i",i),Form("#Delta#Phi%i;#Delta#Phi;Entries",i),500,0.f,TMath::TwoPi());
  }
  for (int i = 0; i < 5; ++i)
  {
    fGDCAZ[i] = new TH1F(Form("DCAGZ%i",i),Form("DCAZ%i;#Deltaz;Entries",i),500,0.f,5.f);
    fGDCAXY[i] = new TH1F(Form("DCAGXY%i",i),Form("DCAXY%i;#Deltar;Entries",i),500,0.f,5.f);
    fFDCAZ[i] = new TH1F(Form("DCAFZ%i",i),Form("DCAZ%i;#Deltaz;Entries",i),500,0.f,5.f);
    fFDCAXY[i] = new TH1F(Form("DCAFXY%i",i),Form("DCAXY%i;#Deltar;Entries",i),500,0.f,5.f);
  }
  for(int i = 0; i < 4; ++i)
  {
    fGoodCombChi2[i] = new TH1F(Form("goodcombchi2%i",i),Form("%i;#chi^{2};Entries",i),2000,0,0.1);
    fFakeCombChi2[i] = new TH1F(Form("fakecombchi2%i",i),Form("%i;#chi^{2};Entries",i),2000,0,0.1);
    fGoodCombN[i] = new TH1F(Form("goodcombn%i",i),Form("%i;#Deltan;Entries",i),300,0.f,0.03f);
    fFakeCombN[i] = new TH1F(Form("fakecombn%i",i),Form("%i;#Deltan;Entries",i),300,0.f,0.03f);
  }
  fGoodCombChi2[4] = new TH1F("goodcomb4",";#chi^{2};Entries",200,0,500);
  fFakeCombChi2[4] = new TH1F("fakecomb4",";#chi^{2};Entries",200,0,500);
  fTan = new TH1F("tan","tan",2500,0,0.5);
  fTanF = new TH1F("tanF","tanF",2500,0,0.5);
  fPhi = new TH1F("phi","phi",2500,0,TMath::Pi());
  fPhiF = new TH1F("phi","phiF",2500,0,TMath::Pi());
  fNEntries = new TH1F("nentries","nentries",2001,-0.5,2000.5);
#endif
}

//__________________________________________________________________________________________________
AliITSUCATracker::AliITSUCATracker(const AliITSUCATracker &t): AliTracker(t),
  fReconstructor(t.fReconstructor),
  fITS(t.fITS),
  fMatLUT(t.fMatLUT),
  fUseMatLUT(t.fUseMatLUT),
  fCurrMass(t.fCurrMass),
  fLayer(),
  fUsedClusters(),
  fChi2Cut(fgkChi2Cut),
  fPhiCut(),
  fZCut(0.5f),
  fCandidates()
{
  // The copy constructor is protected
}

//__________________________________________________________________________________________________
AliITSUCATracker::~AliITSUCATracker()
{
  // d-tor
   for (int i = 0; i < 4; ++i)
  {
    if (fCandidates[i])
      delete fCandidates[i];
  }
  delete fMatLUT;
}

//__________________________________________________________________________________________________
void AliITSUCATracker::Init(AliITSUReconstructor* rec)
{
  // init with external reconstructor
  //
  fITS = rec->GetITSInterface();
  //
  // create material lookup table
  const int kNTest = 1000;
  const double kStepsPerCM = 5;
  fMatLUT  = new AliITSUMatLUT(fITS->GetRMin(),fITS->GetRMax(), \
    Nint(kStepsPerCM * (fITS->GetRMax() - fITS->GetRMin())));
  double zmn = 1e6;
  for (int ilr = fITS->GetNLayers(); ilr--;)
  {
    AliITSURecoLayer* lr = fITS->GetLayer(ilr);
    if (zmn>Abs(lr->GetZMin())) zmn = Abs(lr->GetZMin());
    if (zmn>Abs(lr->GetZMax())) zmn = Abs(lr->GetZMax());
  }
  fMatLUT->FillData(kNTest,-zmn,zmn);
  //
}

//__________________________________________________________________________________________________
void AliITSUCATracker::CellsTreeTraversal(vector<AliITSUCARoad> &roads,
                                          const int &iD, const int &doubl)
{

  if (doubl < 0) return;

  roads.back().AddElement(doubl,iD);
  //Road tmp = roads.back();
  const int currentN = roads.back().N;
  for (size_t iN = 0; iN < fCells[doubl][iD].NumberOfNeighbours(); ++iN)
  {
    const int currD = doubl - 1;
    const int neigh = fCells[doubl][iD](iN);

    if (iN > 0)
    {
      roads.push_back(roads.back());
      roads.back().N = currentN;
    }

    CellsTreeTraversal(roads,neigh,currD);
  }

  fCells[doubl][iD].SetLevel(0u); // Level = -1
}

//__________________________________________________________________________________________________
Int_t AliITSUCATracker::Clusters2Tracks(AliESDEvent *event)
{
  // This is the main tracking function
  // The clusters must already be loaded

  int ntrk = 0, ngood = 0;
  for (int iteration = 0; iteration < fgkNumberOfIterations; ++iteration)
  {

    fCandidates[0]->Clear();
    fCandidates[1]->Clear();
    fCandidates[2]->Clear();
    fCandidates[3]->Clear();

    FindTracksCA(iteration);

    for (int iL = 3; iL >= 0; --iL)
    {
      const int nCand = fCandidates[iL]->GetEntries();
      int index[nCand];
      float chi2[nCand];

      for (int iC = 0; iC < nCand; ++iC)
      {
        AliITSUTrackCooked *tr = (AliITSUTrackCooked*)fCandidates[iL]->At(iC);
//        Int_t clInfo[2 * AliITSUAux::kMaxLayers];
//        for (unsigned int i = 0; i < 2 * AliITSUAux::kMaxLayers; ++i)
//        {
//          clInfo[i] = -1;
//        }
//        for (int k = 0; k < tr->GetNumberOfClusters(); ++k)
//        {
//          const int layer = (tr->GetClusterIndex(k) & 0xf0000000) >> 28;
//          const int idx = (tr->GetClusterIndex(k) & 0x0fffffff);
//          clInfo[layer << 1] = idx;
//        }
        chi2[iC] = tr->GetChi2();//(RefitTrack(tr,clInfo,0.,-1));
        index[iC] = iC;
        CookLabel(tr,0.f);
      }

      TMath::Sort(nCand,chi2,index,false);

      for (int iUC = 0; iUC < nCand; ++iUC)
      {
        const int iC = index[iUC];
        if (chi2[iC] < 0.f)
        {
          continue;
        }
//        if (chi2[iC] > fChi2Cut)
//        {
//          break;
//        }

        AliITSUTrackCooked *tr = (AliITSUTrackCooked*)fCandidates[iL]->At(iC);

        bool sharingCluster = false;
        for (int k = 0; k < tr->GetNumberOfClusters(); ++k)
        {
          const int layer = (tr->GetClusterIndex(k) & 0xf0000000) >> 28;
          const int idx = (tr->GetClusterIndex(k) & 0x0fffffff);
          if (fUsedClusters[layer][idx])
          {
            sharingCluster = true;
            break;
          }
        }

        if (sharingCluster)
          continue;

        for (int k = 0; k < tr->GetNumberOfClusters(); ++k)
        {
          const int layer = (tr->GetClusterIndex(k) & 0xf0000000) >> 28;
          const int idx = (tr->GetClusterIndex(k) & 0x0fffffff);
          fUsedClusters[layer][idx] = true;
        }

        AliESDtrack outTrack;
        CookLabel(tr,0.f);
        ntrk++;
        if(tr->GetLabel() >= 0)
        {
          ngood++;
#ifdef _TUNING_
          fGoodCombChi2[4]->Fill(chi2[iC]);
        }
        else
        {
          fFakeCombChi2[4]->Fill(chi2[iC]);
#endif
        }

        outTrack.UpdateTrackParams(tr,AliESDtrack::kITSin);
        outTrack.SetLabel(tr->GetLabel());
        event->AddTrack(&outTrack);
      }
    }
  }
  Info("Clusters2Tracks","Reconstructed tracks: %d",ntrk);
  if (ntrk)
    Info("Clusters2Tracks","Good tracks/reconstructed: %f",Float_t(ngood)/ntrk);
//
  return 0;
}

//__________________________________________________________________________________________________
void AliITSUCATracker::FindTracksCA(int)
{

#ifdef _TUNING_
  unsigned int numberOfGoodDoublets = 0, totalNumberOfDoublets = 0;
  unsigned int numberOfGoodCells = 0, totalNumberOfCells = 0;
  unsigned int cellsCombiningSuccesses = 0, cellsWrongCombinations = 0;
  unsigned int totalNumberOfRoads = 0;
#endif
  
  vector<int> dLUT[5];
  for (int iL = 0; iL < 6; ++iL) {
    if (iL < 5) {
      dLUT[iL].resize(fLayer[iL + 1].GetNClusters(),0);
    }
    for (int iC = 0; iC < fLayer[iL].GetNClusters(); ++iC) {
      ClsInfo_t* cls = fLayer[iL].GetClusterInfo(iC);
      const float tanL = (cls->z - GetZ()) / cls->r;
      const float extz = tanL * (fgkR[iL + 1] - cls->r) + cls->z;
      const int nClust = fLayer[iL + 1].SelectClusters(extz - 5 * fZCut, extz + 5 * fZCut,
                                                       cls->phi - fPhiCut, cls->phi + fPhiCut);
      bool first = true;
   
      for (int iC2 = 0; iC2 < nClust; ++iC2) {
        const int iD2 = fLayer[iL + 1].GetNextClusterInfoID();
        ClsInfo_t* cls2 = fLayer[iL + 1].GetClusterInfo(iD2);
        const float dz = tanL * (cls2->r - cls->r) + cls->z - cls2->z;
        if (TMath::Abs(dz) < kmaxDZ[iL] && CompareAngles(cls->phi, cls2->phi, fPhiCut)) {
          if (first && iL > 0) {
            dLUT[iL - 1][iC] = fDoublets[iL].size();
            first = false;
          }
          const float dTanL = (cls->z - cls2->z) / (cls->r - cls2->r);
          const float phi = TMath::ATan2(cls->y - cls2->y, cls->x - cls2->x);
          fDoublets[iL].push_back(Doublets(iC,iD2,dTanL,phi));
#ifdef _TUNING_
          if (fLayer[iL].GetClusterSorted(iC)->GetLabel(0) == fLayer[iL + 1].GetClusterSorted(iD2)->GetLabel(0) &&
              fLayer[iL].GetClusterSorted(iC)->GetLabel(0) > 0) {
            numberOfGoodDoublets++;
            fGDZ[iL]->Fill(dz);
            fGDXY[iL]->Fill(fabs(cls->phi - cls2->phi));
          } else {
            fFDZ[iL]->Fill(dz);
            fFDXY[iL]->Fill(fabs(cls->phi - cls2->phi));
          }
          totalNumberOfDoublets++;
#endif
        }
      }
      fLayer[iL + 1].ResetFoundIterator();
    }
  }

  vector<int> tLUT[4];
  for (int iD = 0; iD < 5; ++iD)
  {
    if (iD < 4) {
      tLUT[iD].resize(fDoublets[iD + 1].size(),0);
    }
    for (size_t iD0 = 0; iD0 < fDoublets[iD].size(); ++iD0)
    {
      const int idx = fDoublets[iD][iD0].y;
      bool first = true;
      for (size_t iD1 = dLUT[iD][idx]; idx == fDoublets[iD + 1][iD1].x;++iD1)
      {
        //cout << dLUT[iD][fDoublets[iD][iD0].y] << " " << dLUT[iD][fDoublets[iD][iD0].y + 1] << " " << fLayer[iD + 1].GetNClusters() << endl;
        if (TMath::Abs(fDoublets[iD][iD0].tanL - fDoublets[iD + 1][iD1].tanL) < 0.025 && // TODO: cuts as parameters
            TMath::Abs(fDoublets[iD][iD0].phi - fDoublets[iD + 1][iD1].phi) < 0.14) {    // TODO: cuts as parameters
          const float tan = 0.5f * (fDoublets[iD][iD0].tanL + fDoublets[iD + 1][iD1].tanL);
          const float extz = -tan * fLayer[iD][fDoublets[iD][iD0].x]->r + fLayer[iD][fDoublets[iD][iD0].x]->z;
          if (fabs(extz - GetZ()) < kmaxDCAz[iD]) {
#ifdef _TUNING_
            fGood = (fLayer[iD].GetClusterSorted(fDoublets[iD][iD0].x)->GetLabel(0) == fLayer[iD + 1].GetClusterSorted(fDoublets[iD][iD0].y)->GetLabel(0) &&
                     fLayer[iD].GetClusterSorted(fDoublets[iD][iD0].x)->GetLabel(0) == fLayer[iD + 2].GetClusterSorted(fDoublets[iD + 1][iD1].y)->GetLabel(0) &&
                     fLayer[iD].GetClusterSorted(fDoublets[iD][iD0].x)->GetLabel(0) > 0);
#endif
            float curv, n[3];
            if (CellParams(iD, fLayer[iD][fDoublets[iD][iD0].x], fLayer[iD + 1][fDoublets[iD][iD0].y],
                           fLayer[iD + 2][fDoublets[iD + 1][iD1].y], curv, n)) {
              if (first && iD > 0) {
                tLUT[iD - 1][iD0] = fCells[iD].size();
                first = false;
              }
              fCells[iD].push_back(AliITSUCACell(fDoublets[iD][iD0].x,fDoublets[iD][iD0].y,
                                                 fDoublets[iD + 1][iD1].y,iD0,iD1,curv,n));
#ifdef _TUNING_
              if (fGood) {
                fTan->Fill(TMath::Abs(fDoublets[iD][iD0].tanL - fDoublets[iD + 1][iD1].tanL));
                fPhi->Fill(TMath::Abs(fDoublets[iD][iD0].phi - fDoublets[iD + 1][iD1].phi));
                fGDCAZ[iD]->Fill(fabs(extz-GetZ()));
                numberOfGoodCells++;
              } else {
                fTanF->Fill(TMath::Abs(fDoublets[iD][iD0].tanL - fDoublets[iD + 1][iD1].tanL));
                fPhiF->Fill(TMath::Abs(fDoublets[iD][iD0].phi - fDoublets[iD + 1][iD1].phi));
                fFDCAZ[iD]->Fill(fabs(extz - GetZ()));
              }
              totalNumberOfCells++;
#endif
            }
          }
        }
      }
    }
  }
  //
  for (int iD = 0; iD < 4; ++iD) {
    for (size_t c0 = 0; c0 < fCells[iD].size(); ++c0) {
      const int idx = fCells[iD][c0].d1();
      for (size_t c1 = tLUT[iD][idx]; idx == fCells[iD + 1][c1].d0(); ++c1) {
#ifdef _TUNING_
        fGood = (fLayer[iD].GetClusterSorted(fCells[iD][c0].x())->GetLabel(0) == fLayer[iD + 1].GetClusterSorted(fCells[iD][c0].y())->GetLabel(0) &&
                 fLayer[iD + 1].GetClusterSorted(fCells[iD][c0].y())->GetLabel(0) == fLayer[iD + 2].GetClusterSorted(fCells[iD][c0].z())->GetLabel(0) &&
                 fLayer[iD + 2].GetClusterSorted(fCells[iD][c0].z())->GetLabel(0) == fLayer[iD + 3].GetClusterSorted(fCells[iD + 1][c1].z())->GetLabel(0) &&
                 fLayer[iD].GetClusterSorted(fCells[iD][c0].x())->GetLabel(0) > 0);
#endif
        float *n0 = fCells[iD][c0].GetN();
        float *n1 = fCells[iD + 1][c1].GetN();
        const float dn2 = ((n0[0] - n1[0]) * (n0[0] - n1[0]) + (n0[1] - n1[1]) * (n0[1] - n1[1]) +
                           (n0[2] - n1[2]) * (n0[2] - n1[2]));
        const float dp = fabs(fCells[iD][c0].GetCurvature() - fCells[iD + 1][c1].GetCurvature());
        if (dn2 < kmaxDN[iD] && dp < kmaxDP[iD]) {
#ifdef _TUNING_
          assert(fCells[iD + 1][c1].Combine(fCells[iD][c0], c0));
          if (fGood) {
            fGoodCombChi2[iD]->Fill(dp);
            fGoodCombN[iD]->Fill(dn2);
            cellsCombiningSuccesses++;
          } else {
            fFakeCombChi2[iD]->Fill(dp);
            fFakeCombN[iD]->Fill(dn2);
            cellsWrongCombinations++;
          }
#else
          fCells[iD + 1][c1].Combine(fCells[iD][c0], c0);
#endif
        }
      }
    }
  }
  //
  for (int level = 5; level > 1; --level) {
    vector<AliITSUCARoad> roads;
    
    for (int iCL = 4; iCL >= level - 1; --iCL) {
      for (size_t iCell = 0; iCell < fCells[iCL].size(); ++iCell) {
        if (fCells[iCL][iCell].GetLevel() != level)
        {
          continue;
        }
        roads.push_back(AliITSUCARoad(iCL,iCell));
        for(size_t iN = 0; iN < fCells[iCL][iCell].NumberOfNeighbours(); ++iN) {
          const int currD = iCL - 1;
          const int neigh = fCells[iCL][iCell](iN);
          if(iN > 0)
          {
            roads.push_back(AliITSUCARoad(iCL,iCell));
          }
          CellsTreeTraversal(roads,neigh,currD);
        }
        fCells[iCL][iCell].SetLevel(0u); // Level = -1
      }
    }
    
    for (size_t iR = 0; iR < roads.size(); ++iR)
    {
      if (roads[iR].N != level)
        continue;
#ifdef _TUNING_
      ++totalNumberOfRoads;
#endif
      AliITSUTrackCooked tr;
      int first = -1,last = -1;
      ClsInfo_t *cl0 = 0x0,*cl1 = 0x0,*cl2 = 0x0;
      for(int i = 0; i < 5; ++i)
      {
        if (roads[iR][i] < 0)
          continue;
        
        if (first < 0)
        {
          cl0 = fLayer[i][fCells[i][roads[iR][i]].x()];
          tr.SetClusterIndex(i,fLayer[i][fCells[i][roads[iR][i]].x()]->index);
          tr.SetClusterIndex(i + 1,fLayer[i + 1][fCells[i][roads[iR][i]].y()]->index);
          first = i;
        }
        tr.SetClusterIndex(i + 2,fLayer[i + 2][fCells[i][roads[iR][i]].z()]->index);
        last = i;
      }
      AliITSUClusterPix* c = fLayer[last + 2].GetClusterSorted(fCells[last][roads[iR][last]].z());
      cl2 = fLayer[last + 2][fCells[last][roads[iR][last]].z()];
      first = (last + first) / 2;
//      if (roads[iR][first] < 0) {
//        printf("roads[%lu][%i] : %i %i %i %i %i\n",iR,first,roads[iR][0],roads[iR][1],roads[iR][2],roads[iR][3],roads[iR][4]);
//      }
      cl1 = fLayer[first + 1][fCells[first][roads[iR][first]].y()];
      // Init track parameters
      double cv = Curvature(cl0->x,cl0->y,cl1->x,cl1->y,cl2->x,cl2->y);
      double tgl = TanLambda(cl0->x,cl0->y,cl2->x,cl2->y,cl0->z,cl2->z);
      
      AliITSUCATrackingStation::ITSDetInfo_t det = fLayer[last + 2].GetDetInfo(cl2->detid);
      double x = det.xTF + c->GetX(); // I'd like to avoit using AliITSUClusterPix...
      double alp = det.phiTF;
      double par[5] = {c->GetY(),c->GetZ(),0,tgl,cv};
      double cov[15] = {
        5*5,
        0,  5*5,
        0,  0  , 0.7*0.7,
        0,  0,   0,       0.7*0.7,
        0,  0,   0,       0,       10
      };
      tr.Set(x,alp,par,cov);
      AliITSUTrackCooked tt(tr);
      tt.ResetClusters();
      if (RefitAt(2.1, &tt, &tr))
        new((*fCandidates[level - 2])[fCandidates[level - 2]->GetEntriesFast()]) AliITSUTrackCooked(tt);
    }
  }
#ifdef _TUNING_
//  Info("FindTracksCA","Iteration %i",iteration);
  Info("FindTracksCA","Good doublets: %d",numberOfGoodDoublets);
  Info("FindTracksCA","Number of doublets: %d",totalNumberOfDoublets);
  Info("FindTracksCA","Good cells: %d",numberOfGoodCells);
  Info("FindTracksCA","Number of cells: %d",totalNumberOfCells);
//  Info("FindTracksCA","Cells combining inefficiencies: %d",cellsCombiningInefficiencies);
  Info("FindTracksCA","Cells combining successes: %d",cellsCombiningSuccesses);
  Info("FindTracksCA","Cells wrong combinations: %d",cellsWrongCombinations);
//  Info("FindTracksCA","Number of found roads: %d",numberOfRoads);
//  Info("FindTracksCA","Number of FULL found roads: %d",numberOfFullRoads);
  Info("FindTracksCA","Roads survived after first selection: %d",totalNumberOfRoads);
#endif
}

//__________________________________________________________________________________________________
Int_t AliITSUCATracker::PropagateBack(AliESDEvent * event)
{
  // Here, we implement the Kalman smoother ?
  // The clusters must be already loaded
//  Int_t n=event->GetNumberOfTracks();
//  Int_t ntrk=0;
//  Int_t ngood=0;
//  for (Int_t i = 0; i < n; i++)
//  {
//    AliESDtrack *esdTrack=event->GetTrack(i);
//
//    if ((esdTrack->GetStatus()&AliESDtrack::kITSin) == 0)
//      continue;
//
//    AliITSUTrackCooked track(*esdTrack);
//
//    track.ResetCovariance(10.);
//
//    int points[2 * AliITSUAux::kMaxLayers];
//    for (UInt_t k = 0; k < 2 * AliITSUAux::kMaxLayers; k++)
//      points[k] = -1;
//    Int_t nc = track.GetNumberOfClusters();
//    for (Int_t k = 0; k < nc; k++)
//    {
//      const int layer = (track.GetClusterIndex(k) & 0xf0000000) >> 28;
//      const int idx = (track.GetClusterIndex(k) & 0x0fffffff);
//      points[layer << 1]=idx;
//    }
//
//    if (RefitTrack(&track,points,40,1) >= 0) {
//
//      CookLabel(&track, 0.); //For comparison only
//      Int_t label = track.GetLabel();
//      if (label > 0)
//        ngood++;
//
//      esdTrack->UpdateTrackParams(&track,AliESDtrack::kITSout);
//      ntrk++;
//    }
//  }
//
//  Info("PropagateBack","Back propagated tracks: %d",ntrk);
//  if (ntrk)
//    Info("PropagateBack","Good tracks/back propagated: %f",Float_t(ngood)/ntrk);
//
//  return 0;
  
  Int_t n=event->GetNumberOfTracks();
  Int_t ntrk=0;
  Int_t ngood=0;
  for (Int_t i=0; i<n; i++) {
    AliESDtrack *esdTrack=event->GetTrack(i);
    
    if ((esdTrack->GetStatus()&AliESDtrack::kITSin)==0) continue;
    
    AliITSUTrackCooked track(*esdTrack);
    AliITSUTrackCooked temp(*esdTrack);
    
    temp.ResetCovariance(10.);
    temp.ResetClusters();
    
    if (RefitAt(40., &temp, &track)) {
      
      CookLabel(&temp, 0.); //For comparison only
      Int_t label = temp.GetLabel();
      if (label > 0) ngood++;
      
      esdTrack->UpdateTrackParams(&temp,AliESDtrack::kITSout);
      //UseClusters(fTrackToFollow);
      ntrk++;
    }
  }
  
  Info("PropagateBack","Back propagated tracks: %d",ntrk);
  if (ntrk)
    Info("PropagateBack","Good tracks/back propagated: %f",Float_t(ngood)/ntrk);
  
  return 0;
}

//__________________________________________________________________________________________________
Int_t AliITSUCATracker::RefitInward(AliESDEvent * event)
{
  // Some final refit, after the outliers get removed by the smoother ?
  // The clusters must be loaded
//
//  Int_t n = event->GetNumberOfTracks();
//  Int_t ntrk = 0;
//  Int_t ngood = 0;
//  for (Int_t i = 0; i < n; i++) {
//    AliESDtrack *esdTrack=event->GetTrack(i);
//
//    if ((esdTrack->GetStatus()&AliESDtrack::kITSout) == 0) continue;
//
//    AliITSUTrackCooked track(*esdTrack);
//
//    track.ResetCovariance(10.);
//
//    int points[2 * AliITSUAux::kMaxLayers];
//    for (UInt_t k = 0; k < 2 * AliITSUAux::kMaxLayers; k++)
//      points[k] = -1;
//    Int_t nc = track.GetNumberOfClusters();
//    for (Int_t k = 0; k < nc; k++)
//    {
//      const int layer = (track.GetClusterIndex(k) & 0xf0000000) >> 28;
//      const int idx = (track.GetClusterIndex(k) & 0x0fffffff);
//      points[layer << 1] = idx;
//    }
//
//    if (RefitTrack(&track,points,1.8,1) >= 0) { //2.1,1)>=0) {
//
//      //if (!track.PropagateTo(1.8, 2.27e-3, 35.28*1.848)) continue;
//      CookLabel(&track, 0.); //For comparison only
//      Int_t label=track.GetLabel();
//      if (label>0) ngood++;
//
//      //cout << esdTrack->GetStatus() << " ";
//      esdTrack->UpdateTrackParams(&track,AliESDtrack::kITSrefit);
//      //cout << esdTrack->GetStatus() << endl;
//      ntrk++;
//    }
//  }
//
//  Info("RefitInward","Refitted tracks: %d",ntrk);
//  if (ntrk)
//    Info("RefitInward","Good tracks/refitted: %f",Float_t(ngood) / ntrk);
//
//  return 0;
  Int_t n = event->GetNumberOfTracks();
  Int_t ntrk = 0;
  Int_t ngood = 0;
  for (Int_t i = 0; i < n; i++) {
    AliESDtrack *esdTrack = event->GetTrack(i);
    
    if ((esdTrack->GetStatus() & AliESDtrack::kITSout) == 0) continue;
    
    AliITSUTrackCooked track(*esdTrack);
    AliITSUTrackCooked temp(*esdTrack);
    
    temp.ResetCovariance(10.);
    temp.ResetClusters();
  
    if (!RefitAt(2.1, &temp, &track)) continue;
    //Cross the beam pipe
    if (!temp.PropagateTo(1.8, 2.27e-3, 35.28 * 1.848)) continue;
    
    CookLabel(&temp, 0.); //For comparison only
    Int_t label = temp.GetLabel();
    if (label > 0) ngood++;
    
    esdTrack->UpdateTrackParams(&temp,AliESDtrack::kITSrefit);
    //esdTrack->RelateToVertex(event->GetVertex(),GetBz(),33.);
    //UseClusters(fTrackToFollow);
    ntrk++;
  }
  
  Info("RefitInward","Refitted tracks: %d",ntrk);
  if (ntrk)
    Info("RefitInward","Good tracks/refitted: %f",Float_t(ngood)/ntrk);
  
  return 0;
}

//__________________________________________________________________________________________________
Int_t AliITSUCATracker::LoadClusters(TTree *cluTree)
{
  // This function reads the ITSU clusters from the tree,
  // sort them, distribute over the internal tracker arrays, etc

  fITS->LoadClusters(cluTree);
  fITS->ProcessClusters();
  //
  // I consider a single vertex event for the moment.
  double vertex[3] = {GetX(),GetY(),GetZ()};
  for (int iL = 0; iL < 7; ++iL) {
    fLayer[iL].Init(fITS->GetLayerActive(iL),fITS->GetGeom());
    AliVertex v(vertex,1,1);
    fLayer[iL].SortClusters(&v);
    fUsedClusters[iL].resize(fLayer[iL].GetNClusters(),false);
  }
  return 0;
}

//__________________________________________________________________________________________________
void AliITSUCATracker::UnloadClusters()
{
  // This function unloads ITSU clusters from the memory
  for (int i = 0;i < 7;++i) {
    fLayer[i].Clear();
    fUsedClusters[i].clear();
  }
  for (int i = 0; i < 6; ++i) {
    fDoublets[i].clear();
  }
  for (int i = 0; i < 5; ++i) {
    fCells[i].clear();
  }
  for (int i = 0; i < 4; ++i)
  {
    fCandidates[i]->Clear("C");
  }
  // for(int i=0;i<5;++i)
  //   fCells[i].clear();

}

//__________________________________________________________________________________________________
Double_t AliITSUCATracker::RefitTrack(AliExternalTrackParam* trc,
                                      Int_t clInfo[2*AliITSUAux::kMaxLayers],
                                      Double_t rDest, Int_t stopCond)
{
  // refit track till radius rDest.
  // if stopCond<0 : propagate till last cluster then stop
  // if stopCond==0: propagate till last cluster then try to go till limiting
  //                 rDest, don't mind if fail
  // if stopCond>0 : rDest must be reached
  //
  // The clList should provide the indices of clusters at corresponding layer
  // (as stored in the layer TClonesArray, with convention (allowing for up to 2
  // clusters per layer due to the overlaps): if there is a cluster on given
  // layer I, then it should be stored at clInfo[2*I-1] if there is an
  // additional cluster on this layer, it goes to clInfo[2*I],
  // -1 means no cluster
  //
  double rCurr = Sqrt(trc->GetX()*trc->GetX() + trc->GetY()*trc->GetY());
  int dir,lrStart,lrStop;
  //
  dir = rCurr<rDest ? 1 : -1;
  lrStart = fITS->FindFirstLayerID(rCurr,dir);
  lrStop  = fITS->FindLastLayerID(rDest,dir);//lrid before which we have to stop
  //
  if (lrStop<0 || lrStart<0)
  {
    AliFatal(Form("Failed to find start(%d) or last(%d) layers. "
             "Track from %.3f to %.3f",lrStart,lrStop,rCurr,rDest));
  }
  //
  int nCl = 0;
  for (int i=2*fITS->GetNLayersActive();i--;)  {
    if (clInfo[i]<0)
      continue;
    nCl++;
  }
  //
  AliExternalTrackParam tmpTr(*trc);
  double chi2 = 0;
  int iclLr[2],nclLr;
  int nclFit = 0;
  //

  int lrStop1 = lrStop+dir;
  for (int ilr=lrStart;ilr!=lrStop1;ilr+=dir) {
    AliITSURecoLayer* lr = fITS->GetLayer(ilr);
    if ( dir*(rCurr-lr->GetR(dir))>0) continue; // this layer is already passed
    int ilrA2,ilrA = lr->GetActiveID();
    // passive layer or active w/o hits will be traversed on the way to next
    // cluster
    if (!lr->IsActive() || clInfo[ilrA2=(ilrA<<1)]<0) continue;
    //
    // select the order in which possible 2 clusters (in case of the overlap)
    // will be traversed and fitted
    nclLr=0;
    if (dir>0) { // clusters are stored in increasing radius order
      iclLr[nclLr++]=clInfo[ilrA2++];
      if (clInfo[ilrA2]>=0) iclLr[nclLr++]=clInfo[ilrA2];
    }
    else {
      if ( clInfo[ilrA2+1]>=0 ) iclLr[nclLr++]=clInfo[ilrA2+1];
      iclLr[nclLr++]=clInfo[ilrA2];
    }
    //
    Bool_t transportedToLayer = kFALSE;
    for (int icl=0;icl<nclLr;icl++) {
      AliITSUClusterPix* clus =  fLayer[ilrA].GetClusterUnSorted(iclLr[icl]);
      AliITSURecoSens* sens = lr->GetSensorFromID(clus->GetVolumeId());
      if (!tmpTr.Rotate(sens->GetPhiTF())) return -1;
      //
      double xClus = sens->GetXTF()+clus->GetX();
      if (!transportedToLayer) {
        if (ilr!=lrStart && !TransportToLayerX(&tmpTr,lrStart,ilr,xClus)) {
          return -1; // go to the entrance to the layer
        }
        lrStart = ilr;
        transportedToLayer = kTRUE;
      }
      //
      if (!PropagateSeed(&tmpTr,xClus,fCurrMass)) return -1;
      //
      Double_t p[2]={clus->GetY(), clus->GetZ()};
      Double_t cov[3]= \
        {clus->GetSigmaY2(), clus->GetSigmaYZ(), clus->GetSigmaZ2()};
      double chi2cl = tmpTr.GetPredictedChi2(p,cov);
      chi2 += chi2cl;
      //
      if ( !tmpTr.Update(p,cov) ) return -1;
      if (++nclFit==nCl && stopCond<0) {
        *trc = tmpTr;
        return chi2; // it was requested to not propagate after last update
      }
    }
    //
  }
  // All clusters were succesfully fitted. Even if the track does not reach
  // rDest, this is enough to validate it. Still, try to go as close as possible
  // to rDest.
  //
  if (lrStart!=lrStop) {
    if (!TransportToLayer(&tmpTr,lrStart,lrStop))
      return (stopCond>0) ? -chi2 : chi2; // rDest was obligatory
    if (!GoToExitFromLayer(&tmpTr,fITS->GetLayer(lrStop),dir))
      return (stopCond>0) ? -chi2 : chi2; // rDest was obligatory
  }
  // go to the destination radius. Note that here we don't select direction
  // to avoid precision problems
  if (!tmpTr.GetXatLabR(rDest,rDest,GetBz(),0) || \
   !PropagateSeed(&tmpTr,rDest,fCurrMass, 100, kFALSE)) {
    return (stopCond>0) ? -chi2 : chi2; // rDest was obligatory
  }
  *trc = tmpTr;

  return chi2;
}

//__________________________________________________________________________________________________
Bool_t AliITSUCATracker::PropagateSeed(AliExternalTrackParam *seed, Double_t xToGo, Double_t mass,
                                       Double_t maxStep, Bool_t matCorr)
{
  // propagate seed to given x applying material correction if requested
  const Double_t kEpsilon = 1e-5;
  Double_t xpos     = seed->GetX();
  Int_t dir         = (xpos < xToGo) ? 1 : -1;
  Double_t xyz0[3],xyz1[3];
  //
  Bool_t updTime = dir > 0 && seed->IsStartedTimeIntegral();
  if (matCorr || updTime) seed->GetXYZ(xyz1);   //starting global position
  while ( (xToGo - xpos) * dir > kEpsilon){
    Double_t step = dir * TMath::Min(TMath::Abs(xToGo-xpos), maxStep);
    Double_t x    = xpos + step;
    Double_t bz = GetBz();   // getting the local Bz
    if (!seed->PropagateTo(x,bz))  return kFALSE;
    double ds = 0;
    if (matCorr || updTime) {
      xyz0[0] = xyz1[0]; // global pos at the beginning of step
      xyz0[1] = xyz1[1];
      xyz0[2] = xyz1[2];
      seed->GetXYZ(xyz1);    //  // global pos at the end of step
      //
      if (matCorr)
      {
      	Double_t xrho,xx0;
      	ds = GetMaterialBudget(xyz0,xyz1,xx0,xrho);
        if (dir>0) xrho = -xrho; // outward should be negative
        if (!seed->CorrectForMeanMaterial(xx0,xrho,mass)) return kFALSE;
      }
      else
      { // matCorr is not requested but time integral is
      	double d0 = xyz1[0] - xyz0[0];
      	double d1 = xyz1[1] - xyz0[1];
      	double d2 = xyz1[2] - xyz0[2];
      	ds = TMath::Sqrt(d0 * d0 + d1 * d1 + d2 * d2);
      }
    }
    if (updTime) seed->AddTimeStep(ds);
    //
    xpos = seed->GetX();
  }
  return kTRUE;
}

//__________________________________________________________________________________________________
Bool_t AliITSUCATracker::TransportToLayer(AliExternalTrackParam* seed, Int_t lFrom, Int_t lTo,
                                          Double_t rLim)
{
  // transport track from layerFrom to the entrance of layerTo or to rLim (if>0), wathever is closer
  //
  if (lTo == lFrom) AliFatal(Form("was called with lFrom=%d lTo=%d",lFrom,lTo));
  //
  int dir = lTo > lFrom ? 1 : -1;
  // lrFr can be 0 when extrapolation from TPC to ITS is requested
  AliITSURecoLayer* lrFr = fITS->GetLayer(lFrom);
  Bool_t checkFirst = kTRUE;
  Bool_t limReached = kFALSE;
  while(lFrom != lTo)
  {
    if (lrFr)
    {
      if (!GoToExitFromLayer(seed,lrFr,dir,checkFirst))
      {
        return kFALSE; // go till the end of current layer
      }
      checkFirst = kFALSE;
    }
    AliITSURecoLayer* lrTo =  fITS->GetLayer( (lFrom += dir) );
    if (!lrTo) AliFatal(Form("Layer %d does not exist",lFrom));
    //
    // go the entrance of the layer, assuming no materials in between
    double xToGo = lrTo->GetR(-dir);
    if (rLim > 0)
    {
      if (dir > 0)
      {
      	if (rLim < xToGo)
        {
          xToGo = rLim;
          limReached = kTRUE;
        }
      }
      else
      {
      	if (rLim > xToGo)
        {
          xToGo = rLim;
          limReached = kTRUE;
        }
      }
    }
    //    double xts = xToGo;
    if (!seed->GetXatLabR(xToGo,xToGo,GetBz(),dir))
    {
      //      printf("FailHere1: %f %f %d\n",xts,xToGo,dir);
      //      seed->Print("etp");
      return kFALSE;
    }
    if (!PropagateSeed(seed,xToGo,fCurrMass,100, kFALSE ))
    {
      //printf("FailHere2: %f %f %d\n",xts,xToGo,dir);
      //seed->Print("etp");
      return kFALSE;
    }
    lrFr = lrTo;
    if (limReached) break;
  }
  return kTRUE;
  //
}

//__________________________________________________________________________________________________
Bool_t AliITSUCATracker::TransportToLayerX(AliExternalTrackParam* seed, Int_t lFrom, Int_t lTo,
	                                         Double_t xStop)
{
  // transport track from layerFrom to the entrance of layerTo but do not pass
  // control parameter X
  if (lTo==lFrom) AliFatal(Form("was called with lFrom=%d lTo=%d",lFrom,lTo));
  //
  int dir = lTo > lFrom ? 1 : -1;
  // lrFr can be 0 when extrapolation from TPC to ITS is requested
  AliITSURecoLayer* lrFr = fITS->GetLayer(lFrom);
  Bool_t checkFirst = kTRUE;
  while(lFrom != lTo)
  {
    if (lrFr) {
      if (!GoToExitFromLayer(seed,lrFr,dir,checkFirst))
      {
        return kFALSE; // go till the end of current layer
      }
      checkFirst = kFALSE;
    }
    AliITSURecoLayer* lrTo =  fITS->GetLayer( (lFrom += dir) );
    if (!lrTo) AliFatal(Form("Layer %d does not exist",lFrom));
    //
    // go the entrance of the layer, assuming no materials in between
    double xToGo = lrTo->GetR(-dir); // R of the entrance to layer
    //
    //    double xts = xToGo;
    if (!seed->GetXatLabR(xToGo,xToGo,GetBz(),dir)) {
      //      printf("FailHere1: %f %f %d\n",xts,xToGo,dir);
      //      seed->Print("etp");
      return kFALSE;
    }
    if ( (dir > 0 && xToGo > xStop) || (dir < 0 && xToGo < xStop) ) xToGo = xStop;
    //
	#ifdef _ITSU_DEBUG_
    AliDebug(2,Form("go in dir=%d to R=%.4f(X:%.4f)",dir,lrTo->GetR(-dir), xToGo));
	#endif
    if (!PropagateSeed(seed,xToGo,fCurrMass,100, kFALSE )) {
      //printf("FailHere2: %f %f %d\n",xts,xToGo,dir);
      //seed->Print("etp");
      return kFALSE;
    }
    lrFr = lrTo;
  }
  return kTRUE;
  //
}

//__________________________________________________________________________________________________
Bool_t AliITSUCATracker::GoToExitFromLayer(AliExternalTrackParam* seed, AliITSURecoLayer* lr,
	                                         Int_t dir, Bool_t check)
{
  // go to the exit from lr in direction dir, applying material corrections in steps specific
  // for this layer.
  // If check is requested, do this only provided the track has not exited the layer already
  double xToGo = lr->GetR(dir);
  if (check) { // do we need to track till the surface of the current layer ?
    double curR2 = seed->GetX() * seed->GetX() + seed->GetY() * seed->GetY(); // current radius
    if (dir > 0) {
    	if (curR2 - xToGo*xToGo > -fgkToler) return kTRUE;
    } // on the surface or outside of the layer
    else if (dir < 0) {
    	if (xToGo * xToGo - curR2 > -fgkToler) return kTRUE;
    } // on the surface or outside of the layer
  }
  if (!seed->GetXatLabR(xToGo,xToGo,GetBz(),dir)) return kFALSE;
  // go via layer to its boundary, applying material correction.
  if (!PropagateSeed(seed,xToGo,fCurrMass, lr->GetMaxStep())) return kFALSE;
  //
  return kTRUE;
  //
}

//__________________________________________________________________________________________________
Bool_t AliITSUCATracker::GoToEntranceToLayer(AliExternalTrackParam* seed, AliITSURecoLayer* lr,
	                                           Int_t dir, Bool_t check)
{
  // go to the entrance of lr in direction dir, w/o applying material corrections.
  // If check is requested, do this only provided the track did not reach the layer already
  double xToGo = lr->GetR(-dir);
  if (check) { // do we need to track till the surface of the current layer ?
    double curR2 = seed->GetX() * seed->GetX() + seed->GetY() * seed->GetY(); // current radius
    if      (dir>0) { if (curR2 - xToGo * xToGo > -fgkToler) return kTRUE; } // already passed
    else if (dir<0) { if (xToGo * xToGo - curR2 > -fgkToler) return kTRUE; } // already passed
  }
  if (!seed->GetXatLabR(xToGo,xToGo,GetBz(),dir)) return kFALSE;
  // go via layer to its boundary, applying material correction.
  if (!PropagateSeed(seed,xToGo,fCurrMass, 100, kFALSE)) return kFALSE;
  return kTRUE;
  //
}

//__________________________________________________________________________________________________
Double_t AliITSUCATracker::GetMaterialBudget(const double* pnt0,const double* pnt1, double& x2x0,
	                                           double& rhol) const
{
  double par[7];
  if (fUseMatLUT && fMatLUT) {
    double d = fMatLUT->GetMatBudget(pnt0,pnt1,par);
    x2x0 = par[AliITSUMatLUT::kParX2X0];
    rhol = par[AliITSUMatLUT::kParRhoL];
    return d;
  }
  else {
    MeanMaterialBudget(pnt0,pnt1,par);
    x2x0 = par[1];
    rhol = par[0]*par[4];
    return par[4];
  }
}

//__________________________________________________________________________________________________
inline AliCluster* AliITSUCATracker::GetCluster(Int_t index) const
{
	const Int_t l=(index & 0xf0000000) >> 28;
	const Int_t c=(index & 0x0fffffff);
	return (AliCluster*)fLayer[l].GetClusterUnSorted(c) ;
}

//__________________________________________________________________________________________________
bool AliITSUCATracker::CellParams(int l, ClsInfo_t* c1, ClsInfo_t* c2, ClsInfo_t* c3,
                                  float &curv, float n[3])
{
  // Calculation of cell params and filtering using a DCA cut wrt beam line position.
  
  // Mapping the hits
  const float mHit0[3] = {c1->x, c1->y, c1->r * c1->r};
  const float mHit1[3] = {c2->x, c2->y, c2->r * c2->r};
  const float mHit2[3] = {c3->x, c3->y, c3->r * c3->r};
  // Computing the deltas
  const float mD10[3] = {mHit1[0] - mHit0[0],mHit1[1] - mHit0[1],mHit1[2] - mHit0[2]};
  const float mD20[3] = {mHit2[0] - mHit0[0],mHit2[1] - mHit0[1],mHit2[2] - mHit0[2]};
  // External product of the deltas -> n
  n[0] = (mD10[1] * mD20[2]) - (mD10[2] * mD20[1]);
  n[1] = (mD10[2] * mD20[0]) - (mD10[0] * mD20[2]);
  n[2] = (mD10[0] * mD20[1]) - (mD10[1] * mD20[0]);
  // Normalisation
  const float norm = TMath::Sqrt((n[0] * n[0]) + (n[1] * n[1]) + (n[2] * n[2]));
  if (norm < 1e-20f || fabs(n[2]) < 1e-20f)
    return false;
  n[0] /= norm;
  n[1] /= norm;
  n[2] /= norm;
  // Center of the circle
  const float c[2] = {-0.5f * n[0] / n[2], -0.5f * n[1] / n[2]};
  // Constant
  const float k = - n[0] * mHit1[0] - n[1] * mHit1[1] - n[2] * mHit1[2];
  // Radius of the circle
  curv = TMath::Sqrt((1.f - n[2] * n[2] - 4.f * k * n[2]) / (4.f * n[2] * n[2]));
  // Distance of closest approach to the beam line
  const float dca = fabs(curv - sqrt(c[0] * c[0] + c[1] * c[1]));
  // Cut on the DCA
  if (dca > kmaxDCAxy[l]) {
    return false;
  }
#ifdef _TUNING_
  if (fGood) {
    fGDCAXY[l]->Fill(dca);
  } else {
    fFDCAXY[l]->Fill(dca);
  }
#endif
  
  curv = 1.f / curv;
  return true;
}

//__________________________________________________________________________________________________
Bool_t AliITSUCATracker::RefitAt(Double_t xx, AliITSUTrackCooked *t, const AliITSUTrackCooked *c) {
  //--------------------------------------------------------------------
  // This function refits the track "t" at the position "x" using
  // the clusters from "c"
  //--------------------------------------------------------------------
  
  const int nLayers = 7;
  Int_t index[nLayers];
  Int_t k;
  for (k = 0; k < nLayers; k++) index[k] = -1;
  Int_t nc = c->GetNumberOfClusters();
  for (k = 0; k < nc; k++) {
    Int_t idx = c->GetClusterIndex(k), nl = (idx&0xf0000000)>>28;
    index[nl] = idx;
  }
  
  Int_t from, to, step;
  if (xx > t->GetX()) {
    from = 0;
    to = nLayers;
    step = +1;
  } else {
    from = nLayers - 1;
    to = -1;
    step = -1;
  }
  
  for (Int_t i = from; i != to; i += step) {
    Int_t idx = index[i];
    if (idx >= 0) {
      const AliCluster *cl = GetCluster(idx);
      Float_t xr,ar;
      cl->GetXAlphaRefPlane(xr, ar);
      if (!t->Propagate(Double_t(ar), Double_t(xr), GetBz())) {
        //Warning("RefitAt","propagation failed !\n");
        return kFALSE;
      }
      Double_t chi2 = t->GetPredictedChi2(cl);
//      if (chi2 < 100)
      t->Update(cl, chi2, idx);
    } else {
      Double_t r = fgkR[i];
      Double_t phi,z;
      if (!t->GetPhiZat(r,phi,z)) {
        //Warning("RefitAt","failed to estimate track !\n");
        return kFALSE;
      }
      if (!t->Propagate(phi, r, GetBz())) {
        //Warning("RefitAt","propagation failed !\n");
        return kFALSE;
      }
    }
    Double_t xx0 = (i > 2) ? 0.008 : 0.003;  // Rough layer thickness
    Double_t x0  = 9.36; // Radiation length of Si [cm]
    Double_t rho = 2.33; // Density of Si [g/cm^3]
    Double_t mass = t->GetMass();
    t->CorrectForMeanMaterial(xx0, - step * xx0 * x0 * rho, mass, kTRUE);
  }
  
  if (!t->PropagateTo(xx,0.,0.)) return kFALSE;
  return kTRUE;
}

