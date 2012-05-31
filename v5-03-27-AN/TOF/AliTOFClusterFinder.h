#ifndef ALITOFCLUSTERFINDER_H
#define ALITOFCLUSTERFINDER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
// AliTOFClusterFinder Class
// Task: Transform digits/raw data to TOF Clusters, to fill TOF RecPoints
// and feed TOF tracking 

#include "TNamed.h"

#include "AliTOFRawStream.h"

class TClonesArray;
class TFile;
class TTree;

class AliLoader;
class AliRunLoader;
class AliRawReader;

class AliTOFGeometry;
class AliTOFcluster;
class AliTOFcalib;

class AliTOFClusterFinder : public TNamed
{

  enum {kTofMaxCluster=77777}; //maximal number of the TOF clusters

 public:

  AliTOFClusterFinder(AliTOFcalib *calib);
  AliTOFClusterFinder(AliRunLoader* runLoader, AliTOFcalib *calib);
  AliTOFClusterFinder(const AliTOFClusterFinder &source); // copy constructor
  AliTOFClusterFinder& operator=(const AliTOFClusterFinder &source); // ass. op.
  virtual ~AliTOFClusterFinder();

  void Digits2RecPoints(TTree* digitsTree, TTree* clusterTree);
  void Digits2RecPoints(Int_t ievt);
  void Digits2RecPoints(AliRawReader *rawReader, TTree *clustersTree);
  void Digits2RecPoints(Int_t ievt, AliRawReader *rawReader);
  void Raw2Digits(Int_t ievt, AliRawReader *rawReader); // temporary solution
  void Raw2Digits(AliRawReader *rawReader, TTree* digitsTree); 
  void FillRecPoint();
  void ResetRecpoint();
  void Load();
  void LoadClusters();
  void UnLoad();
  void UnLoadClusters();
  void SetVerbose(Int_t Verbose){fVerbose=Verbose;} // To set the verbose level
  void SetDecoderVersion(Int_t version){fDecoderVersion=version;} // To set the decoder version
  Int_t GetDecoderVersion() const {return fDecoderVersion;} // To get the decoder version
  UShort_t  GetClusterVolIndex(const Int_t * const ind) const; //Volume Id getter
  void GetClusterPars(Int_t *ind, Double_t *pos, Double_t *cov) const; //cluster par getter

 protected:
  AliRunLoader *fRunLoader;      // Pointer to Run Loader
  AliLoader    *fTOFLoader;      // Pointer to specific detector loader

  TTree        *fTreeD;          // Digits tree
  TTree        *fTreeR;          // Reconstructed points

  AliTOFcluster *fTofClusters[kTofMaxCluster];  // pointers to the TOF clusters

  TClonesArray *fDigits;         // List of digits
  TClonesArray *fRecPoints;      // List of reconstructed points

  Int_t fNumberOfTofClusters;    // Number of TOF Clusters

 private:

  //Int_t InsertCluster(Int_t *aa, Double_t *bb);    // Fills TofClusters Array
  //Int_t InsertCluster(Int_t *aa, Double_t *bb, Int_t *cc, Int_t d); // Fills TofClusters Array
  Int_t InsertCluster(AliTOFcluster *tofCluster);    // Fills TofClusters Array
  Int_t FindClusterIndex(Double_t z) const; // Returns cluster index 
  void  CalibrateRecPoint(UInt_t timestamp = 0); // Apply calibration pars to Clusters

  Int_t fVerbose;           // Verbose level (0:no msg,
                            //  1:msg, 2:digits in txt files)
  Int_t fDecoderVersion;   //setting whether to use the new decoder version 
  AliTOFcalib *fTOFcalib;         // pointer to the TOF calibration info
  AliTOFRawStream fTOFRawStream; // AliTOFRawStream variable

  ClassDef(AliTOFClusterFinder,7) // To run TOF clustering
};
#endif

