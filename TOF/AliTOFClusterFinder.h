#ifndef ALITOFCLUSTERFINDER_H
#define ALITOFCLUSTERFINDER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
// AliTOFClusterFinder Class
// Task: Transform digits/raw data to TOF Clusters, to fill TOF RecPoints
// and feed TOF tracking 

#include "AliRawReader.h"

class TClonesArray;
class TFile;
class TTree;

class AliLoader;
class AliRunLoader;

class AliTOFGeometry;

class AliTOFClusterFinder : public TObject
{

  enum {kTofMaxCluster=77777}; //maximal number of the TOF clusters

 public:

  AliTOFClusterFinder();
  AliTOFClusterFinder(AliRunLoader* runLoader);
  virtual ~AliTOFClusterFinder();

  void Digits2RecPoints(Int_t ievt);
  void Digits2RecPoints(AliRawReader *rawReader, TTree *clustersTree);
  void Digits2RecPoints(Int_t ievt, AliRawReader *rawReader);
  void Raw2Digits(Int_t ievt, AliRawReader *rawReader); // temporary solution
  void FillRecPoint();
  void ResetRecpoint();
  void Load();
  void LoadClusters();
  void UnLoad();
  void UnLoadClusters();

 protected:
  AliRunLoader *fRunLoader;      // Pointer to Run Loader
  AliLoader    *fTOFLoader;      // Pointer to specific detector loader

  TTree        *fTreeD;          // Digits tree
  TTree        *fTreeR;          // Reconstructed points

  AliTOFcluster *fTofClusters[kTofMaxCluster];  // pointers to the TOF clusters

  AliTOFGeometry  *fTOFGeometry; // Pointer to TOF geometry
  TClonesArray *fDigits;         // List of digits
  TClonesArray *fRecPoints;      // List of reconstructed points

  Int_t fNumberOfTofClusters;    // Number of TOF Clusters

 private:

  //Int_t InsertCluster(Int_t *aa, Double_t *bb);    // Fills TofClusters Array
  //Int_t InsertCluster(Int_t *aa, Double_t *bb, Int_t *cc, Int_t d); // Fills TofClusters Array
  Int_t InsertCluster(AliTOFcluster *tofCluster);    // Fills TofClusters Array
  Int_t FindClusterIndex(Double_t z) const; // Returns cluster index 
  void  CalibrateRecPoint(); // Apply calibration pars to Clusters


  ClassDef(AliTOFClusterFinder,1) // To run TOF clustering
};
#endif

