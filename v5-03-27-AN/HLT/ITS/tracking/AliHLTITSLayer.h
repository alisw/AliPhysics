// ************************************************************************
// This file is property of and copyright by the ALICE HLT Project        *
// ALICE Experiment at CERN, All rights reserved.                         *
// See cxx source for full Copyright notice                               *
//                                                                        *
//*************************************************************************

#ifndef ALIHLTITSLAYER_H
#define ALIHLTITSLAYER_H

#include "AliHLTITSDetector.h"
#include "AliITSRecoParam.h"

/**
 * @class AliHLTITSLayer
 *
 * The AliHLTTPCCAHit class is the internal representation
 * of the TPC clusters for the AliHLTTPCCATracker algorithm.
 *
 */
class AliHLTITSLayer 
{
 public:
  AliHLTITSLayer();
  AliHLTITSLayer(Double_t r, Double_t p, Double_t z, Int_t nl, Int_t nd);
  ~AliHLTITSLayer();
  Int_t InsertCluster(AliITSRecPoint *c);
  void  SortClusters();
  void ResetClusters(); 
  void SelectClusters(Double_t zmin,Double_t zmax,Double_t ymin,Double_t ymax);
  const AliITSRecPoint *GetNextCluster(Int_t &ci,Bool_t test=kFALSE);
  void ResetRoad();
  Double_t GetRoad() const {return fRoad;}
  Double_t GetR() const {return fR;}
  Int_t FindClusterIndex(Float_t z) const;
  AliITSRecPoint *GetCluster(Int_t i) const {return i<fN? fClusters[i]:0;} 
  AliHLTITSDetector &GetDetector(Int_t n) const { return fDetectors[n]; }
  Int_t FindDetectorIndex(Double_t phi, Double_t z) const;
  Double_t GetThickness(Double_t y, Double_t z, Double_t &x0) const;
  Int_t InRoad() const ;
  Int_t GetNumberOfClusters() const {return fN;}
  Int_t GetNladders() const {return fNladders;}
  Int_t GetNdetectors() const {return fNdetectors;}
  Int_t GetSkip() const {return fSkip;}
  void  SetSkip(Int_t skip){fSkip=skip;}
  void IncAccepted(){fAccepted++;}
  Int_t GetAccepted() const {return fAccepted;}    

 protected:

  AliHLTITSLayer(const AliHLTITSLayer& layer);
  AliHLTITSLayer & operator=(const AliHLTITSLayer& layer){
    this->~AliHLTITSLayer();new(this) AliHLTITSLayer(layer);
    return *this;}
  Double_t fR;                // mean radius of this layer
  Double_t fPhiOffset;        // offset of the first detector in Phi
  Int_t fNladders;            // number of ladders
  Double_t fZOffset;          // offset of the first detector in Z
  Int_t fNdetectors;          // detectors/ladder
  AliHLTITSDetector *fDetectors; // array of detectors
  Int_t fN;                   // number of clusters
  AliITSRecPoint *fClusters[AliITSRecoParam::kMaxClusterPerLayer]; // pointers to clusters
  Int_t        fClusterIndex[AliITSRecoParam::kMaxClusterPerLayer]; // pointers to clusters
  Float_t fY[AliITSRecoParam::kMaxClusterPerLayer];                // y position of the clusters      
  Float_t fZ[AliITSRecoParam::kMaxClusterPerLayer];                // z position of the clusters      
  Float_t fYB[2];                                       // ymin and ymax
  //
  AliITSRecPoint *fClusters5[6][AliITSRecoParam::kMaxClusterPerLayer5]; // pointers to clusters -     slice in y
  Int_t        fClusterIndex5[6][AliITSRecoParam::kMaxClusterPerLayer5]; // pointers to clusters -     slice in y    
  Float_t fY5[6][AliITSRecoParam::kMaxClusterPerLayer5];                // y position of the clusters  slice in y    
  Float_t fZ5[6][AliITSRecoParam::kMaxClusterPerLayer5];                // z position of the clusters  slice in y 
  Int_t fN5[6];                                       // number of cluster in slice
  Float_t fDy5;                                       //delta y
  Float_t fBy5[6][2];                                    //slice borders
  //
  AliITSRecPoint *fClusters10[11][AliITSRecoParam::kMaxClusterPerLayer10]; // pointers to clusters -     slice in y
  Int_t        fClusterIndex10[11][AliITSRecoParam::kMaxClusterPerLayer10]; // pointers to clusters -     slice in y    
  Float_t fY10[11][AliITSRecoParam::kMaxClusterPerLayer10];                // y position of the clusters  slice in y    
  Float_t fZ10[11][AliITSRecoParam::kMaxClusterPerLayer10];                // z position of the clusters  slice in y 
  Int_t fN10[11];                                       // number of cluster in slice
  Float_t fDy10;                                        // delta y
  Float_t fBy10[11][2];                                 // slice borders
  //
  AliITSRecPoint *fClusters20[21][AliITSRecoParam::kMaxClusterPerLayer20]; // pointers to clusters -     slice in y
  Int_t        fClusterIndex20[21][AliITSRecoParam::kMaxClusterPerLayer20]; // pointers to clusters -     slice in y    
  Float_t fY20[21][AliITSRecoParam::kMaxClusterPerLayer20];                // y position of the clusters  slice in y    
  Float_t fZ20[21][AliITSRecoParam::kMaxClusterPerLayer20];                // z position of the clusters  slice in y 
  Int_t fN20[21];                                       // number of cluster in slice
  Float_t fDy20;                                        //delta y 
  Float_t fBy20[21][2];                                 //slice borders
  //
  AliITSRecPoint** fClustersCs;                         //clusters table in current slice
  Int_t   *fClusterIndexCs;                             //cluster index in current slice 
  Float_t *fYcs;                                        //y position in current slice
  Float_t *fZcs;                                        //z position in current slice
  Int_t    fNcs;                                        //number of clusters in current slice    
  Int_t fCurrentSlice;                                  //current slice
  //

  Float_t fZmax;      //    edges
  Float_t fYmin;      //   of  the
  Float_t fYmax;      //   "window"
  Int_t fI;            // index of the current cluster within the "window"
  Int_t fImax;            // index of the last cluster within the "window"    
  Int_t fSkip;     // indicates possibility to skip cluster
  Int_t fAccepted;     // accept indicator 
  Double_t fRoad;      // road defined by the cluster density
};


#endif
