#ifndef ALIPMDCLUSTERFINDER_H
#define ALIPMDCLUSTERFINDER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//-----------------------------------------------------//
//                                                     //
//  Date   : August 05 2003                            //
//  This reads the file PMD.digits.root(TreeD),        //
//  calls the Clustering algorithm and stores the      //
//  clustering output in PMD.RecPoints.root(TreeR)     // 
//                                                     //
//-----------------------------------------------------//

class TClonesArray;
class TFile;
class TTree;

class AliLoader;
class AliRunLoader;
class AliRawReader;

class AliPMDClusterFinder
{

 public:

  AliPMDClusterFinder(AliRunLoader* runLoader);
  virtual ~AliPMDClusterFinder();

  void Digits2RecPoints(Int_t ievt);
  void Digits2RecPoints(Int_t ievt, AliRawReader *rawReader);
  void SetCellEdepCut(Float_t ecut);
  void SetDebug(Int_t idebug);
  void AddRecPoint(Int_t idet, Int_t ismn, Float_t * clusdata);
  void ResetCellADC();
  void ResetRecpoint();
  void Load();
  void LoadClusters();
  void UnLoad();
  void UnLoadClusters();

 protected:
  AliRunLoader *fRunLoader; // Pointer to Run Loader
  AliLoader    *fPMDLoader; // Pointer to specific detector loader

  TTree        *fTreeD;     // Digits tree
  TTree        *fTreeR;     // Reconstructed points

  TClonesArray *fDigits;    // List of digits
  TClonesArray *fRecpoints; // List of reconstructed points

  Int_t   fNpoint;          // 
  Int_t   fDetNo;           // Detector Number (0:PRE, 1:CPV)
  Int_t   fDebug;           // Debugging switch (0:NO, 1:YES)
  Float_t fEcut;            // Energy/ADC cut per cell

  static const Int_t fgkRow = 48; // Total number of rows in one unitmodule
  static const Int_t fgkCol = 96; // Total number of cols in one unitmodule
  Double_t fCellADC[fgkRow][fgkCol]; // Array containing individual cell ADC

  ClassDef(AliPMDClusterFinder,6) // To run PMD clustering
};
#endif

