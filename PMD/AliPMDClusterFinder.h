#ifndef PMDClusterFinder_H
#define PMDClusterFinder_H
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
#include <vector>
#include <algorithm>

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
class AliPMDContainer;
class AliPMDcluster;
class AliPMDrecpoint;

using namespace std;

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

  static const Int_t fTotSM = 27;
  static const Int_t fNCell = 72;
  Float_t fCPV[fTotSM][fNCell][fNCell];
  Float_t fPMD[fTotSM][fNCell][fNCell];

 public:

  AliPMDClusterFinder();
  virtual ~AliPMDClusterFinder();

  void OpengAliceFile(char * /* galice.root */, Option_t * /* option */);

  void Digits2RecPoints(Int_t /* ievt */);
  void AddRecPoint(Float_t * /* clusdata */);
  void ResetCellADC();
  void ResetRecpoint();
  void UnLoad(Option_t * /* option */);

  ClassDef(AliPMDClusterFinder,1)
};
#endif

