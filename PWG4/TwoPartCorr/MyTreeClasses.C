// $Id:$

#ifndef MyTreeClasses_cxx
#define MyTreeClasses_cxx

#include <Riostream.h>
#include <TTree.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>

class MyHeader
{
public:
  MyHeader() : fRun(0), fOrbit(0), fTime(0), fPeriod(0), fBx(0),
               fNChips1(0), fNChips2(0), fNTracks(0), fNSelTracks(0), fNTracklets(0), 
               fVz(0), fVc(-1), fIsPileupSPD(0), fNPileupSPD(0), fNPileup(0),
               fTrClassMask(0), fTrCluster(0), fEvNumberInFile(-1), fFileId(-1),
               fVzSPD(0), fVcSPD(-1) {;}
  virtual ~MyHeader() {;}
  ULong64_t     GetEventId() const {
                  return ((ULong64_t)fBx+
                          (ULong64_t)fOrbit*3564+
                          (ULong64_t)fPeriod*16777215*3564);
                }

  Int_t         fRun;
  UInt_t        fOrbit; 
  UInt_t        fTime;   
  UInt_t        fPeriod;
  UShort_t      fBx;
  Short_t       fNChips1;
  Short_t       fNChips2;
  Short_t       fNTracks;
  Short_t       fNSelTracks;
  Short_t       fNTracklets;
  Double_t      fVz; //[0,0,16]
  Double_t      fVc; //[0,0,16]
  Bool_t        fIsPileupSPD; 
  Char_t        fNPileupSPD;
  Char_t        fNPileup;
  ULong64_t     fTrClassMask;
  UChar_t       fTrCluster;
  Int_t         fEvNumberInFile;
  Short_t       fFileId;
  Double_t      fVzSPD; //[0,0,16]
  Double_t      fVcSPD; //[0,0,16]

  ClassDef(MyHeader,3) // My header class
};

class MyPart : public TObject
{
public:
  MyPart(ULong_t st=0, Char_t c=0, Double_t pt=0, Double_t eta=0, Double_t phi=0) :
    TObject(), fSt(st), fC(c), fPt(pt), fEta(eta), fPhi(phi) {;}
  Char_t      Charge() const { return fC;   }
  Double_t    Eta()    const { return fEta; }
  Double_t    Phi()    const { return fPhi; }
  Double_t    Pt()     const { return fPt;  }
  Double_t    Px()     const { return fPt*TMath::Cos(fPhi);  }
  Double_t    Py()     const { return fPt*TMath::Sin(fPhi);  }
  Double_t    Pz()     const { return fPt*TMath::SinH(fEta); }
  ULong_t     Status() const { return fSt;  }
protected:
  ULong_t     fSt;
  Char_t      fC;
  Double32_t  fPt;  //[0,0,14]
  Double32_t  fEta; //[0,0,10]
  Double32_t  fPhi; //[0,0,10]
  ClassDef(MyPart,1) // My particle class in cylindrical coordinates
};

class MyTracklet : public TObject
{
public:
  MyTracklet(Double_t dphi=0, Double_t dth=0, Double_t eta=0, Double_t phi=0) :
    TObject(), fDPhi(dphi), fDTh(dth), fEta(eta), fPhi(phi) {;}
  Double_t    DPhi()   const { return fDPhi; }
  Double_t    DTh()    const { return fDTh;  }
  Double_t    Eta()    const { return fEta;  }
  Double_t    Phi()    const { return fPhi;  }
protected:
  Double32_t  fDPhi; //[0,0,10]
  Double32_t  fDTh;  //[0,0,10]
  Double32_t  fEta;  //[0,0,10]
  Double32_t  fPhi;  //[0,0,10]
  ClassDef(MyTracklet,1) // My tracklet class
};
#endif
