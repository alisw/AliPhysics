#include "AliITSUVertexerFast.h"
#include "AliITSUClusterPix.h"

#include <TBranch.h>
#include <TClonesArray.h>
#include <TString.h>
#include <TTree.h>

#include <TGraph.h>
#include <TFile.h>

#include <sys/time.h>
#include <algorithm>
// #include <math.h>

#ifdef _TUNE
#include <iostream>
using std::cout;
using std::endl;
#include <TCanvas.h>
#include <TH2F.h>
#endif


ClassImp(AliITSUVertexerFast)

AliITSUVertexerFast::AliITSUVertexerFast()
:TObject()
,fGranularityPhi(64) 
,fGranularityZ(256)
,fSizePhi(fGranularityPhi / 6.28318530718f)
,fSizeZ()
,fHalfLayerLength()
,fRLayer()          
,fLUT()
,fClusters()
// ,fLines() 
{

}

AliITSUVertexerFast::AliITSUVertexerFast(int granularityPhi, int granularityZ, float layerSize)
:TObject()
,fGranularityPhi(granularityPhi) 
,fGranularityZ(granularityZ)
,fSizePhi(granularityPhi / 6.28318530718f) 
,fSizeZ(granularityZ / layerSize)
,fHalfLayerLength(0.5f * layerSize)
,fRLayer()
,fClusters()           
,fLUT()
// ,fLines() 
{

}

void AliITSUVertexerFast::LoadClusters(TClonesArray* clusters[3]) {
  fClusters[0].clear();
  fClusters[1].clear();
  fClusters[2].clear(); 

  fLUT[0].clear();
  fLUT[1].clear();
  fLUT[2].clear();
  // fLines.clear();
  
  #ifdef _THCLK
  cout<<"Total threads declared: "<<omp_get_max_threads()<<endl;
  SetTimesVecSize();
  #endif

#pragma omp parallel for
 
  for (int iL = 0; iL < 3; ++iL) {
    #ifdef _THCLK
    StartTime( omp_get_thread_num() );
    #endif

    #ifdef _TUNE
    fL[iL].reserve(clusters[iL]->GetEntriesFast());
    #endif

    for (int iC = 0; iC < clusters[iL]->GetEntriesFast(); ++iC) {
      float xyz[3];
      AliITSUClusterPix* c = (AliITSUClusterPix*)clusters[iL]->At(iC);
      c->GetGlobalXYZ(xyz);
      fClusters[iL].push_back({ xyz[0], xyz[1], xyz[2] });

      #ifdef _TUNE
      fL[iL].push_back(c->GetLabel(0));
      #endif
    }

    std::sort(fClusters[iL].begin(), fClusters[iL].end());

    float size = 6.28318530718f / fGranularityPhi;
    vector<int> &tLUT = fLUT[iL];
    tLUT.reserve(fGranularityPhi+1);
    tLUT.push_back(0);
    for (size_t iC = 0; iC < fClusters[iL].size(); ++iC) {
        fClusters[iL][iC].fID = int(fClusters[iL][iC].fP * fSizePhi);
        while (fClusters[iL][iC].fP > size * tLUT.size()) {
            tLUT.push_back(iC);
        }
    }

    tLUT.resize(fGranularityPhi+1, fClusters[iL].size());
    #ifdef _THCLK
    EndTime( omp_get_thread_num() );
    #endif
  }

  #ifdef _THCLK
  for(int a=0; a<GetMaxNumberOfThreads(); a++) {
      cout<<" -> Thread "<<a+1<<" took "<<ElapsedTicks(a)<<" ticks | "<<ElapsedTime(a)<<" sec"<<endl;
  }
  #endif
}

void AliITSUVertexerFast::FindVertex(float* xyz) {
  /// Trackleting
  cout << "Trackleting" << endl;
  int nThread = 1;
  #ifdef _TUNE
  int nGood = 0;
  #endif
  AliITSUVertexCandidate vertexCandidate;
  int nTot = 0;
  #pragma omp parallel
  {
    nThread = omp_get_num_threads();
    int tid = omp_get_thread_num();
    int n = 0;
    AliITSUVertexCandidate candidate;

    #pragma omp for
    for (size_t iC0 = 0; iC0 < fClusters[0].size(); ++iC0) {
      // int phiI = int(fClusters[0][iC0].fP * fSizePhi);
      for (int iC1 = fLUT[1][fClusters[0][iC0].fID]; iC1 < fLUT[1][fClusters[0][iC0].fID + 1]; ++iC1) {
        if (fabs(fClusters[0][iC0].fP - fClusters[1][iC1].fP) < 0.002) { 
          Line l;
          l.x[0] = fClusters[0][iC0].fX;
          l.x[1] = fClusters[0][iC0].fY;        
          l.x[2] = fClusters[0][iC0].fZ;
          l.c[0] = fClusters[1][iC1].fX - fClusters[0][iC0].fX;
          l.c[1] = fClusters[1][iC1].fY - fClusters[0][iC0].fY;
          l.c[2] = fClusters[1][iC1].fZ - fClusters[0][iC0].fZ;
          const float z = IntersectCylinder(l,fRLayer[2]);
          const int zI = int((z + fHalfLayerLength) * fSizeZ);
          for (int iC2 = fLUT[2][fClusters[0][iC0].fID]; iC2 < fLUT[2][fClusters[0][iC0].fID + 1]; ++iC2) {
            if (fabs(fClusters[2][iC2].fP - fClusters[1][iC1].fP) < 0.002 && fabs(fClusters[2][iC2].fZ - z) < 0.03) {
              candidate.Add(l);
              #ifdef _TUNE
                if (fL[0][iC0] == fL[1][iC1] && fL[1][iC1] == fL[2][iC2]) n++;
              #endif
            }
          }  
        }
      }
    }
    double localVert[3];

    #pragma omp critical
    {
      // fLines.insert(fLines.end(),lines.begin(),lines.end())
      vertexCandidate.Add(candidate);      
      nGood += n;
    }
  }

  vertexCandidate.ComputeClusterCentroid();
  vertexCandidate.GetVertex(xyz);

  #ifdef _TUNE
  cout << vertexCandidate.GetSize() <<" "<< nGood <<" "<< nThread << endl;
  #endif
}

float AliITSUVertexerFast::IntersectCylinder(Line &l, float r) {
  const float a = l.c[0] * l.c[0] + l.c[1] * l.c[1];
  const float b = l.c[0] * l.x[0] + l.c[1] * l.x[1];
  const float c = l.x[0] * l.x[0] + l.x[1] * l.x[1] - r * r;
  const float t = (sqrt(b * b - a * c) - b) / a;
  return l.x[2] + l.c[2] * t;
}
