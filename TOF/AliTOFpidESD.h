#ifndef ALITOFPIDESD_H
#define ALITOFPIDESD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------
//                    TOF PID class
// A very naive design... Should be made better by the detector experts...
//   Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//-------------------------------------------------------

#include <TObject.h>

class AliESD;
class TFile;
class TTree;

const Int_t kMaxCluster=77777; //maximal number of the TOF clusters

class AliTOFpidESD : public TObject {
public:
  AliTOFpidESD(){fR=376.; fDy=2.5; fDz=3.5; fN=0; fEventN=0;}
  AliTOFpidESD(Double_t *param) throw (const Char_t *);
  ~AliTOFpidESD(){UnloadClusters();}

  Int_t MakePID(AliESD *event);
  Int_t LoadClusters(const TFile *f);
  Int_t LoadClusters(TTree *f);
  void  UnloadClusters();
  void SetEventNumber(Int_t n) {fEventN=n;}

  Int_t GetEventNumber() const {return fEventN;}

public:
  class AliTOFcluster {
  public:
    AliTOFcluster(Double_t *h, Int_t *l) {
      fR=h[0]; fPhi=h[1]; fZ=h[2]; fTDC=h[3]; fADC=h[4];
      fLab[0]=l[0]; fLab[1]=l[1]; fLab[2]=l[2];
    }
    void Use() {fADC=-fADC;}

    Double_t GetR() const {return fR;}
    Double_t GetPhi() const {return fPhi;}
    Double_t GetZ()   const {return fZ;}
    Double_t GetTDC() const {return fTDC;}
    Double_t GetADC() const {return TMath::Abs(fADC);}
    Int_t IsUsed() const {return (fADC<0) ? 1 : 0;}
    Int_t GetLabel(Int_t n) const {return fLab[n];}
  private:
    Int_t fLab[3];
    Double_t fR;
    Double_t fPhi;
    Double_t fZ;
    Double_t fTDC;
    Double_t fADC;
  };

private:
  Int_t InsertCluster(AliTOFcluster*);
  Int_t FindClusterIndex(Double_t z) const;

  Int_t fEventN;          //event number

  Double_t fR;            // mean readius of the TOF barrel
  Double_t fDy;           // size of the TOF cell in R*Phi
  Double_t fDz;           // size of the TOF cell in Z

  Double_t fSigma;        // intrinsic TOF resolution
  Double_t fRange;        // one particle type PID range (in sigmas)

  Int_t fN;                  // number of the TOF clusters
  AliTOFcluster *fClusters[kMaxCluster];  // pointers to the TOF clusters

  ClassDef(AliTOFpidESD,1)   // TOF PID class
};

#endif


