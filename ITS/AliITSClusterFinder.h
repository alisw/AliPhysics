#ifndef ALIITSCLUSTERFINDER_H
#define ALIITSCLUSTERFINDER_H


////////////////////////////////////////////////
//  ITS Cluster Finder Class                  //
////////////////////////////////////////////////

#include <TObject.h>
#include <TClonesArray.h>

class AliITSMap;
class AliITSresponse;
class AliITSsegmentation;
class AliITSRawCluster;
class AliITS;
class AliITSdigit;
class AliITSRecPoint;

//---------------------------------------------------------------
class AliITSClusterFinder :public TObject
{
public:
  
  AliITSClusterFinder(AliITSsegmentation *seg, AliITSresponse *resp, TClonesArray *digits);
  AliITSClusterFinder();
  virtual ~AliITSClusterFinder(){
    // destructor
  }
  AliITSClusterFinder(const AliITSClusterFinder &source); // copy constructor
  AliITSClusterFinder& operator=(const AliITSClusterFinder &source); // assignment operator

  virtual void SetResponse(AliITSresponse *response) {
    // set response
    fResponse=response;
  }
  virtual void SetSegmentation(AliITSsegmentation *segmentation) {
    // set segmentation
    fSegmentation=segmentation;
  }
  virtual void SetDigits(TClonesArray *ITSdigits) {
    // set digits
    fDigits=ITSdigits;
    fNdigits = fDigits->GetEntriesFast();
  }
  virtual AliITSdigit* GetDigit(Int_t i){
      return (AliITSdigit*) fDigits->UncheckedAt(i);
  }
  virtual TClonesArray* Digits(){
      return fDigits;
  }
  virtual Int_t   NDigits() {
    // Get Number of Digits
    return fNdigits;
  }
  
  AliITSMap   *Map()  {
    // map
    return fMap;
  }
  //
  virtual void AddCluster(Int_t branch, AliITSRawCluster *c);
  virtual void AddCluster(Int_t branch, AliITSRawCluster *c, AliITSRecPoint &rp);
  
  virtual void FindRawClusters() {
    // Search for raw clusters
  }
  virtual void FindCluster(Int_t i, Int_t j, AliITSRawCluster *c) {
    // find cluster
  }
  
  virtual void Decluster(AliITSRawCluster *cluster) {
    // Decluster
  }
  virtual void SetNperMax(Int_t npermax=3) {
    // Set max. Number of cells per local cluster
    fNperMax = npermax;
  }
  virtual void SetDeclusterFlag(Int_t flag=1) {
    // Decluster ?
    fDeclusterFlag =flag;
  }
  virtual void SetClusterSize(Int_t clsize=3) {
    // Set max. cluster size ; bigger clusters will be rejected
    fClusterSize = clsize;
  }
  virtual void CalibrateCOG() {
    // Self Calibration of COG 
  }
  virtual void CorrectCOG(){
    // correct COG
  }
  
  virtual Bool_t Centered(AliITSRawCluster *cluster) {
    // cluster
    return kTRUE;
  }
  virtual void   SplitByLocalMaxima(AliITSRawCluster *cluster) {
    // split by local maxima
  }
  virtual void   FillCluster(AliITSRawCluster *cluster, Int_t) {
    // fiil cluster
  }
  virtual void   FillCluster(AliITSRawCluster *cluster) {
    // fill cluster
    FillCluster(cluster,1);
  }
  
  // set the fitting methods in the derived classes 

  // data members

  TClonesArray           *fDigits;      // digits
  Int_t                   fNdigits;     // num of digits
  
protected:
  AliITSresponse         *fResponse;      // response
  AliITSsegmentation     *fSegmentation;  //segmentation
  
  Int_t                   fNRawClusters;  // in case we split the cluster
                                          // and want to keep track of 
                                          // the cluster which was splitted
  AliITSMap              *fMap;           // map
  
  Int_t                   fNperMax;       // NperMax
  Int_t                   fDeclusterFlag; // DeclusterFlag
  Int_t                   fClusterSize;   // ClusterSize
  Int_t                   fNPeaks;        // NPeaks
  
  
  ClassDef(AliITSClusterFinder,1) //Class for clustering and reconstruction of space points
    };
#endif







