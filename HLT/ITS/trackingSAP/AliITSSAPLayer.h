#ifndef ALIITSSAPLAYER_H
#define ALIITSSAPLAYER_H

#include <algorithm>
#include <vector>
#include <TObject.h>
#include <TObjArray.h>
#include "AliITSRecPoint.h"
class AliVertex;
class AliITSSAPLayer 
{

 public:
  struct ClsInfo   // cluster info, optionally XY origin at vertex
  { 
    float x,y,z,phi,r;    // lab params
    int   zphibin; // bins is z,phi
    int   index;          // index in RecPoints array
    int   detid;          // detector index //RS ??? Do we need it?
    bool operator<(const ClsInfo &rhs) const {return zphibin<rhs.zphibin;}
    //
  };
  typedef struct ClsInfo ClsInfo_t;
  //
  struct ClBinInfo  // info on bin clusters start, number of clusters
  {
    unsigned short ncl;    // number of clusters
    unsigned short first;  // entry of 1st cluster in sorted vector of ClsInfo
    int            index;  // index in the vector containing cells with non-0 occupancy
  };
  typedef struct ClBinInfo ClBinInfo_t;
  //
  struct ITSDetInfo  // info on sensor
  {
    int index; // sensor vid
    float xTF,xTFmisal,phiTF,sinTF,cosTF; //tracking frame parameters of the detector
  };
  typedef struct ITSDetInfo ITSDetInfo_t;


  AliITSSAPLayer();
  AliITSSAPLayer(int id, float zspan,int nzbins,int nphibins, int buffer=100);
  virtual ~AliITSSAPLayer();
  //
  int     GetVIDOffset()                const {return fVIDOffset;}
  int     GetNClusters()                const {return fNClusters;}
  int     GetNZBins()                   const {return fNZBins;}
  int     GetNPhiBins()                 const {return fNPhiBins;}
  float   GetZMin()                     const {return fZMin;}
  float   GetZMax()                     const {return fZMax;}
  //
  void    SetNZBins(int v)                    {fNZBins = v;}
  void    SetNPhiBins(int v)                  {fNPhiBins = v;}
  void    SetZMin(float v)                    {fZMin = v;}
  void    SetZMax(float v)                    {fZMax = v;}
  //
  void Init(int buffer=100);
  //
  void AddCluster(AliITSRecPoint *cl)         {fClusters->AddAtAndExpand(cl,fNClusters++);}
  //
  void SortClusters(const AliVertex* vtx=0);
  int  GetPhiBin(float phi)             const {return phi*fDPhiInv;}
  int  GetZBin  (float z)               const {return (z-fZMin)*fDZInv;}
  int  GetBinIndex(int iz, int iphi)    const {return iphi*fNZBins + iz;}
  int  GetBinZ(int ipz)                 const {return ipz%fNZBins;}
  int  GetBinPhi(int ipz)               const {return ipz/fNZBins;}
  void GetBinZPhi(int ipz,int &iz,int &iphi) const {iz = GetBinZ(ipz); iphi=GetBinPhi(ipz);}
  //
  int  SelectClusters(float zmin,float zmax,float phimin,float phimax);
  int  GetNFoundBins()                  const {return fFoundBins.size();}
  int  GetFoundBin(int i)               const {return fFoundBins[i];}
  int  GetFoundBinClusters(int i, int &first)  const;
  void ResetFoundIterator();
  AliITSSAPLayer::ClsInfo_t* GetClusterInfo(int i) const {return (AliITSSAPLayer::ClsInfo_t*)&fSortedClInfo[i];}
  AliITSSAPLayer::ClsInfo_t* GetNextClusterInfo();
  int                     GetNextClusterInfoID();
  AliITSRecPoint*         GetNextCluster();
  AliITSRecPoint*         GetClusterSorted(int i)   const {return (AliITSRecPoint*)fClusters->UncheckedAt(fSortedClInfo[i].index);}
  AliITSRecPoint*         GetClusterUnSorted(int i) const {return (AliITSRecPoint*)fClusters->UncheckedAt(i);}
  //
  AliITSSAPLayer::ITSDetInfo_t& GetDetInfo(int id)     const {return (ITSDetInfo_t&)fDetectors[id];}
  Int_t                   GetNDetectors()           const {return fDetectors.size();}

  void         ClearSortedInfo();
  virtual void Clear(Option_t *opt="");
  virtual void Print(Option_t *opt="")  const;

 protected:
  TObjArray* fClusters;       // externally supplied clusters
  int   fLrID;                // layer id
  int   fVIDOffset;           // offset of VID for detectors of this layer
  int   fNClusters;           // N clusters
  //
  float fZMin;                // Zmin
  float fZMax;                // Zmax
  float fDZInv;               // inverse size of Z bin
  float fDPhiInv;             // inverse size of Phi bin
  int   fNZBins;             // N cells in Z
  int   fNPhiBins;           // N cells in Phi
  //
  int   fQueryZBmin;         // min bin in Z of the query
  int   fQueryZBmax;         // max bin in Z of the query
  int   fQueryPhiBmin;       // min bin in phi of the query
  int   fQueryPhiBmax;       // max bin in phi of the query
  ClBinInfo_t* fBins;           // 2D (z,phi) grid of clusters binned in z,phi
  int* fOccBins;              // id's of bins with non-0 occupancy
  int  fNOccBins;             // number of occupied bins
  int  fNFoundClusters;       // number of found clusters in the query zone
  int  fFoundClusterIterator; // at which cluster within the bin we are?
  int  fFoundBinIterator;     // at which foune bin we are?
  std::vector<int>     fFoundBins;    // occupied bins satisfying to query
  std::vector<ClsInfo_t> fSortedClInfo; // processed cluster info
  std::vector<ITSDetInfo_t> fDetectors; // detector params
  //
};

//_____________________________________________________ 
inline int AliITSSAPLayer::GetFoundBinClusters(int i, int &first)  const {
  // set the entry of the first cl.info in the fSortedClInfo 
  // and return n clusters in the bin
  ClBinInfo_t& bin=fBins[GetFoundBin(i)];
  first = bin.first;
  return bin.ncl;
}

//_____________________________________________________ 
inline AliITSRecPoint* AliITSSAPLayer::GetNextCluster() {
  // return next cluster
  ClsInfo_t* cli=GetNextClusterInfo(); 
  return cli ? (AliITSRecPoint*)fClusters->UncheckedAt(cli->index) : 0;
}

//_____________________________________________________________
inline AliITSSAPLayer::ClsInfo_t* AliITSSAPLayer::GetNextClusterInfo()
{
  // return cluster info for next matching cluster
  int id = GetNextClusterInfoID();
  return id<0 ? 0 : (AliITSSAPLayer::ClsInfo_t*)&fSortedClInfo[id];
}

#endif
