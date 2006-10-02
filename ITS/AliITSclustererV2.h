#ifndef ALIITSCLUSTERERV2_H
#define ALIITSCLUSTERERV2_H
//--------------------------------------------------------------
//                       ITS clusterer V2
//
//   This can be a "wrapping" for the V1 cluster finding classes
//   if compiled with uncommented "#define V1" line 
//   in the AliITSclustererV2.cxx file.
//
//   Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//--------------------------------------------------------------
#include <TObject.h>

class TFile;
class TTree;
class TClonesArray;

class AliITSgeom;
class AliITSclusterV2;
class AliRawReader;
class AliITSRawStream;

class AliITSclustererV2 : public TObject {
public:
  AliITSclustererV2();
  AliITSclustererV2(const AliITSgeom *geom);

  void SetEvent(Int_t event) { fEvent=event; }
  Int_t Digits2Clusters(TTree *in, TTree *out);
  void Digits2Clusters(AliRawReader* rawReader);
  void FindClustersSPD(const TClonesArray *dig, TClonesArray *cls);
  void FindClustersSPD(AliITSRawStream* input, TClonesArray** clusters);
  void FindClustersSDD(const TClonesArray *dig, TClonesArray *cls);
  void FindClustersSDD(AliITSRawStream* input, TClonesArray** clusters);
  void FindClustersSSD(const TClonesArray *dig, TClonesArray *cls);
  void FindClustersSSD(AliITSRawStream* input, TClonesArray** clusters);

  void RecPoints2Clusters(const TClonesArray *p, Int_t idx, TClonesArray *c);

private:
  class Ali1Dcluster {
  private:
    Float_t fY; //cluster position
    Float_t fQ; //cluster charge
    Int_t fNd;  //number of digits
    Int_t fLab[3]; //track label
  public:
    void SetY(Float_t y) {fY=y;}
    void SetQ(Float_t q) {fQ=q;}
    void SetNd(Int_t n)  {fNd=n;}
    void SetLabels(Int_t *lab) {fLab[0]=lab[0];fLab[1]=lab[1];fLab[2]=lab[2];}
    Float_t GetY() const {return fY;}
    Float_t GetQ() const {return fQ;}
    Int_t GetNd()const {return fNd;}
    Int_t GetLabel(Int_t lab) const { return fLab[lab]; }
  };
  class AliBin {
  public:
    AliBin():fIndex(0),fMask(0xFFFFFFFE),fQ(0) {}
    void SetIndex(UInt_t idx) {fIndex=idx;}
    void SetQ(UShort_t q)  {fQ=q;}
    void SetMask(UInt_t m) {fMask=m;}

    void Use() {fMask&=0xFFFFFFFE;}
    Bool_t IsNotUsed() const {return (fMask&1);}
    Bool_t IsUsed() const {return !(IsNotUsed());}

    UInt_t   GetIndex() const {return fIndex;}
    UShort_t GetQ()     const {return fQ;}
    UInt_t   GetMask()  const {return fMask;}
  private:
    UInt_t fIndex; //digit index
    UInt_t fMask; //peak mask
    UShort_t fQ;  //signal
  };
  static Bool_t IsMaximum(Int_t k, Int_t max, const AliBin *bins);
  static void FindPeaks(Int_t k,Int_t m,AliBin*b,Int_t*idx,UInt_t*msk,Int_t&n);
  static void MarkPeak(Int_t k, Int_t max, AliBin *bins, UInt_t m);
  static void MakeCluster(Int_t k,Int_t max,AliBin *bins,UInt_t m,
   AliITSclusterV2 &c);
  void FindClustersSDD(AliBin* bins[2], Int_t nMaxBin, Int_t nMaxZ,
		       const TClonesArray *dig, TClonesArray *cls);
  void FindClustersSSD(Ali1Dcluster* neg, Int_t nn, 
		       Ali1Dcluster* pos, Int_t np,
		       TClonesArray *clusters);


  static void FindCluster(Int_t k,Int_t maxz,AliBin *bins,Int_t &n,Int_t *idx);

protected:
  Int_t fNModules;             // total number of modules

private:
  Int_t fEvent;           //event number

  Int_t fI;                     //index of the current subdetector
  Float_t fYshift[2200];       //y-shifts of detector local coor. systems 
  Float_t fZshift[2200];       //z-shifts of detector local coor. systems 
  Int_t fNdet[2200];            //detector index  
  Int_t fNlayer[2200];          //detector layer

  //SPD related values:
  Int_t fLastSPD1;       //index of the last SPD1 detector
  Int_t fNySPD;          //number of pixels in Y
  Int_t fNzSPD;          //number of pixels in Z
  Float_t fYpitchSPD;    //pixel size in Y
  Float_t fZ1pitchSPD,fZ2pitchSPD;    //pixel sizes in Z
  Float_t fHwSPD;        //half width of the SPD detector
  Float_t fHlSPD;        //half length of the SPD detector
  Float_t fYSPD[260];    //Y-coordinates of pixel centers
  Float_t fZSPD[170];    //Z-coordinates of pixel centers

  //SDD related values:
  Int_t fNySDD;           //number of "pixels" in Y
  Int_t fNzSDD;           //number of "pixels" in Z
  Float_t fYpitchSDD;     //"pixel size" in Y (drift direction)
  Float_t fZpitchSDD;     //"pixel sizes" in Z
  Float_t fHwSDD;         //half width of the SDD detector
  Float_t fHlSDD;         //half length of the SDD detector
  Float_t fYoffSDD;       //some delay in the drift channel   

  //SSD related values:
  Int_t fLastSSD1;        //index of the last SSD1 detector   
  Float_t fYpitchSSD;     //strip pitch (cm)
  Float_t fHwSSD;         //half-width of an SSD detector (cm)
  Float_t fHlSSD;         //half-length of an SSD detector (cm)
  Float_t fTanP;          //tangent of the stereo angle on the P side
  Float_t fTanN;          //tangent of the stereo angle on the N side

  ClassDef(AliITSclustererV2,1)  // ITS cluster finder V2
};

#endif
