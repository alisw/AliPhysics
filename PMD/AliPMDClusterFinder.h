#ifndef PMDClusterFinder_H
#define PMDClusterFinder_H
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

#include <Riostream.h>
#include <stdlib.h>
#include <math.h>
#include <TMath.h>

class TClonesArray;
class TFile;
class TObjArray;
class TTree;
class TNtuple;

class AliLoader;
class AliRunLoader;
class AliRun;
class AliDetector;
class AliHeader;

class AliPMDdigit;
class AliPMDClustering;
class AliPMDcluster;
class AliPMDrecpoint1;

class AliPMDClusterFinder
{
 protected:
  AliRunLoader *fRunLoader;
  AliRun       *gAlice;
  AliDetector  *PMD;      /* Get pointers to Alice detectors 
			     and Hits containers */
  AliLoader    *pmdloader;

  TTree        *treeD;
  TTree        *treeR;

  TClonesArray *fDigits;
  TClonesArray *fRecpoints;

  Int_t fNpoint;
  Int_t fDetNo;
  Int_t fDebug;
  Float_t fEcut;

  static const Int_t fRow = 48;
  static const Int_t fCol = 96;
  Double_t fCellADC[fRow][fCol];

 public:

  AliPMDClusterFinder();
  virtual ~AliPMDClusterFinder();

  void OpengAliceFile(char * /* galice.root */, Option_t * /* option */);

  void Digits2RecPoints(Int_t /* ievt */);
  void SetCellEdepCut(Float_t /* ecut */);
  void SetDebug(Int_t /* idebug */);
  void AddRecPoint(Float_t * /* clusdata */);
  void ResetCellADC();
  void ResetRecpoint();
  void UnLoad(Option_t * /* option */);

  ClassDef(AliPMDClusterFinder,2)
};
#endif

