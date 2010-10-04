// $Id$

#ifndef TreeClasses_h
#define TreeClasses_h

#include <Riostream.h>
#include <TTree.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>

class Noti: public TObject
{
public:
  Noti() : fc(0) {;}
  virtual ~Noti() {;}
  Bool_t Notify()         { fc=1; return 1; }
  Bool_t Notified() const { return fc;      }
  void   Reset()          { fc=0;           }
protected:
  Bool_t fc; //=1 when file changed
  ClassDef (Noti,0) // Use to be notified when file in chain changes
};

class MyHeader
{
public:
  MyHeader() : fRun(0), fOrbit(0), fTime(0), fPeriod(0), fBx(0), fL0(0), fL1(0), fL2(0),
               fNChips1(0), fNChips2(0), fNTracks(0), fNSelTracks(0), fNTracklets(0), 
               fVx(0), fVy(0), fVz(0), fVc(-1), fIsPileupSPD(0), fNPileupSPD(0), fNPileup(0),
               fTrClassMask(0), fTrCluster(0), fEvNumberInFile(-1), fFileId(-1),
               fVxSPD(0), fVySPD(0), fVzSPD(0), fVcSPD(-1) {;}
  virtual ~MyHeader() {;}
  ULong64_t    GetEventId() const {
                 return (((ULong64_t)fPeriod << 36) |
                         ((ULong64_t)fOrbit  << 12) |
                         (ULong64_t)fBx); 
               }

public:
  Int_t         fRun;
  UInt_t        fOrbit; 
  UInt_t        fTime;   
  UInt_t        fPeriod;
  UShort_t      fBx;
  UInt_t        fL0;
  UInt_t        fL1;
  UShort_t      fL2;
  Short_t       fNChips1;
  Short_t       fNChips2;
  Short_t       fNTracks;
  Short_t       fNSelTracks;
  Short_t       fNTracklets;
  Double_t      fVx;             //[0,0,16]
  Double_t      fVy;             //[0,0,16]
  Double_t      fVz;             //[0,0,16]
  Double_t      fVc;             //[0,0,16]
  Bool_t        fIsPileupSPD; 
  Char_t        fNPileupSPD;
  Char_t        fNPileup;
  ULong64_t     fTrClassMask;
  UChar_t       fTrCluster;
  Int_t         fEvNumberInFile;
  Short_t       fFileId;
  Double_t      fVxSPD;          //[0,0,16]
  Double_t      fVySPD;          //[0,0,16]
  Double_t      fVzSPD;          //[0,0,16]
  Double_t      fVcSPD;          //[0,0,16]
  Double_t      fVxTPC;          //[0,0,16]
  Double_t      fVyTPC;          //[0,0,16]
  Double_t      fVzTPC;          //[0,0,16]
  Double_t      fVcTPC;          //[0,0,16]

  ClassDef(MyHeader,4) // My header class
};

class MyPart : public TObject
{
public:
  enum { /*from AliESDtrack.h */
    kITSin=0x0001,kITSout=0x0002,kITSrefit=0x0004,kITSpid=0x0008,
    kTPCin=0x0010,kTPCout=0x0020,kTPCrefit=0x0040,kTPCpid=0x0080,
    kTRDin=0x0100,kTRDout=0x0200,kTRDrefit=0x0400,kTRDpid=0x0800,
    kTOFin=0x1000,kTOFout=0x2000,kTOFrefit=0x4000,kTOFpid=0x8000,
    kHMPIDout=0x10000,kHMPIDpid=0x20000,
    kEMCALmatch=0x40000,
    kPHOSmatch=0x200000,
    kTRDbackup =0x80000,
    kTRDStop=0x20000000,
    kESDpid=0x40000000,
    kTIME=0x80000000,
    kGlobalMerge=0x08000000,
    kITSpureSA=0x10000000,
    kMultInV0=0x2000000,    //BIT(25): assumed to be belong to V0 in multiplicity estimates
    kMultSec=0x4000000     //BIT(26): assumed to be secondary (due to the DCA) in multiplicity estimates
  }; 

  MyPart(ULong_t st=0, Char_t c=0, Double_t pt=0, Double_t eta=0, Double_t phi=0,
         Short_t ncltpc=0, Short_t ncltpc1=0, Short_t ncltpcs=0,Char_t nclits=0, Double_t chi2tpc=0,
         Double_t chi2tpc1=0, Double_t chi2its=0, Double_t d=0, Double_t z=0, Double_t dtpc=0, Double_t ztpc=0) :
    TObject(), fSt(st), fC(c), fPt(pt), fEta(eta), fPhi(phi), fNClTPC(ncltpc), fNClTPC1(ncltpc1), 
    fNClTPCShared(ncltpcs),fNClITS(nclits), fChi2TPC(chi2tpc), fChi2TPC1(chi2tpc1), fChi2ITS(chi2its), 
    fD(d), fZ(z), fDTPC(dtpc), fZTPC(ztpc) {;}


  Double_t    Px()         const { return fPt*TMath::Cos(fPhi);  }
  Double_t    Py()         const { return fPt*TMath::Sin(fPhi);  }
  Double_t    Pz()         const { return fPt*TMath::SinH(fEta); }
  Bool_t      IsITSRefit() const { return (fSt&(ULong64_t)kITSrefit);    }
  Bool_t      IsTPCIn()    const { return (fSt&(ULong64_t)kTPCin);       }
  Bool_t      IsTPCRefit() const { return (fSt&(ULong64_t)kITSrefit);    }

public:
  ULong64_t   fSt;
  Char_t      fC;
  Double32_t  fPt;           //[0,0,16]
  Double32_t  fEta;          //[0,0,10]
  Double32_t  fPhi;          //[0,0,10]
  Short_t     fNClTPC;
  Short_t     fNClTPC1;
  Short_t     fNClTPCShared;
  Char_t      fNClITS;
  Double32_t  fChi2TPC;      //[0,0,10]
  Double32_t  fChi2TPC1;     //[0,0,10]
  Double32_t  fChi2ITS;      //[0,0,10]
  Double32_t  fD;            //[0,0,16]
  Double32_t  fZ;            //[0,0,16]
  Double32_t  fDTPC;         //[0,0,16]
  Double32_t  fZTPC;         //[0,0,16]

  ClassDef(MyPart,3) // My particle class in cylindrical coordinates
};

class MyTracklet : public TObject
{
public:
  MyTracklet(Double_t dphi=0, Double_t dth=0, Double_t eta=0, Double_t phi=0) :
    TObject(), fDPhi(dphi), fDTh(dth), fEta(eta), fPhi(phi) {;}

public:
  Double32_t  fDPhi; //[0,0,10]
  Double32_t  fDTh;  //[0,0,10]
  Double32_t  fEta;  //[0,0,10]
  Double32_t  fPhi;  //[0,0,10]

  ClassDef(MyTracklet,1) // My tracklet class
};
#endif
