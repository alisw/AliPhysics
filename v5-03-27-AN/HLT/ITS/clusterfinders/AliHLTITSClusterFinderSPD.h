#ifndef ALIHLTITSCLUSTERFINDERSPD_H
#define ALIHLTITSCLUSTERFINDERSPD_H
//--------------------------------------------------------------
//                       ITS clusterer V2 for SPD
//
//   This can be a "wrapping" for the V1 cluster finding classes
//   if compiled with uncommented "#define V1" line 
//   in the AliITSclustererV2.cxx file.
//
//   Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//--------------------------------------------------------------
#include "AliITSDetTypeRec.h"
class TClonesArray;
class AliRawReader;
class AliITSRawStream;
class AliITSRawStreamSPD;
class AliITSDetTypeRec;
class AliITSRecoParam;

class AliHLTITSClusterFinderSPD : public TObject 
{
public:
  AliHLTITSClusterFinderSPD(AliITSDetTypeRec* dettyp);
 
  virtual ~AliHLTITSClusterFinderSPD(){
    delete[] fSignal2Bin;
    delete[] fBin2Signal;
  }
  
  void RawdataToClusters(AliRawReader* rawReader,std::vector<AliITSRecPoint> & clusters);  

 protected:
  
  void FindClustersSPD(AliITSRawStreamSPD* input,std::vector<AliITSRecPoint> & clusters);
  Int_t ClustersSPD( std::vector<AliITSRecPoint> & clusters,Int_t iModule );
  void FindCluster(Int_t k,Int_t &n,Int_t *idx);

  AliITSRecoParam *fRecoParam;
  AliITSDetTypeRec* fDetTypeRec; //ITS object for reconstruction
 // Data members needed to fill AliCluster objects
  Int_t fNdet[2200];           // detector index  
  Int_t fNlayer[2200];         // detector layer  
  Int_t fNModules;             // total number of modules    

  Int_t fLastSPD1;       //index of the last SPD1 detector
  Int_t fNySPD;          //number of pixels in Y
  Int_t fNzSPD;          //number of pixels in Z
  Int_t fNzBins;
  Int_t fNyBins;
  Int_t fMaxBin;

  Float_t fYpitchSPD;    //pixel size in Y
  Float_t fZ1pitchSPD,fZ2pitchSPD;    //pixel sizes in Z
  Float_t fHwSPD;        //half width of the SPD detector
  Float_t fHlSPD;        //half length of the SPD detector
  Float_t fYSPD[260];    //Y-coordinates of pixel centers
  Float_t fZSPD[170];    //Z-coordinates of pixel centers

  Int_t fNSignals;
  UShort_t *fSignal2Bin; // array of input SPD signals
  UShort_t *fBin2Signal; // 2D map of SPD signals;

 private:

  AliHLTITSClusterFinderSPD(const AliHLTITSClusterFinderSPD&);
  AliHLTITSClusterFinderSPD &operator=( const AliHLTITSClusterFinderSPD &);

  ClassDef(AliHLTITSClusterFinderSPD,0)  // ITS cluster finder V2 for SPD
};

#endif
