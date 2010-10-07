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
  ULong64_t     GetEventId() const {
                  return (((ULong64_t)fPeriod << 36) |
                          ((ULong64_t)fOrbit  << 12) |
                          (ULong64_t)fBx); }
  Bool_t        IsL0Fired(Int_t bit)      const { return fL0 & (1<<bit);          }
  Bool_t        IsTrClassFired(Int_t bit) const { return fTrClassMask & (1<<bit); }
public:
  Int_t         fRun;            //         run number
  UInt_t        fOrbit;          //         orbit number
  UInt_t        fTime;           //         timestamp from DAQ
  UInt_t        fPeriod;         //         period number
  UShort_t      fBx;             //         bunch crossing id
  UInt_t        fL0;             //         l0 trigger bits
  UInt_t        fL1;             //         l1 trigger bits
  UShort_t      fL2;             //         l2 trigger bits
  Short_t       fNChips1;        //         number of chips layer 1
  Short_t       fNChips2;        //         number of chips layer 2
  Short_t       fNTracks;        //         number of tracks before cuts
  Short_t       fNSelTracks;     //         number of stored tracks
  Short_t       fNTracklets;     //         number of tracklets
  Double_t      fVx;             //[0,0,16] global vertex x
  Double_t      fVy;             //[0,0,16] global vertex y
  Double_t      fVz;             //[0,0,16] global vertex z
  Double_t      fVc;             //[0,0,16] number of contributors to global vertex
  Bool_t        fIsPileupSPD;    //         is pileup according to standard definition
  Char_t        fNPileupSPD;     //         number of additional vertices in SPD
  Char_t        fNPileup;        //         number of additional vertices in TPC
  ULong64_t     fTrClassMask;    //         trigger class mask
  UChar_t       fTrCluster;      //         trigger cluster mask
  Int_t         fEvNumberInFile; //         event number in ESD file
  Short_t       fFileId;         //         file id to access path name
  Double_t      fVxSPD;          //[0,0,16] SPD vertex x
  Double_t      fVySPD;          //[0,0,16] SPD vertex y
  Double_t      fVzSPD;          //[0,0,16] SPD vertex z
  Double_t      fVcSPD;          //[0,0,16] number of contributors to SPD vertex
  Double_t      fVxTPC;          //[0,0,16] TPC vertex x
  Double_t      fVyTPC;          //[0,0,16] TPC vertex y
  Double_t      fVzTPC;          //[0,0,16] TPC vertex z
  Double_t      fVcTPC;          //[0,0,16] number of contributors to TPC vertex

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
    kMultInV0=0x2000000,
    kMultSec=0x4000000
  }; 

  MyPart(ULong_t st=0, Char_t c=0, Double_t pt=0, Double_t eta=0, Double_t phi=0,
         Short_t ncltpc=0, Short_t ncltpc1=0, Short_t ncltpcs=0,Char_t nclits=0, Double_t chi2tpc=0,
         Double_t chi2tpc1=0, Double_t chi2its=0, Double_t d=0, Double_t z=0, Double_t dtpc=0, 
         Double_t ztpc=0, Double_t sigma=0.) :
    TObject(), fSt(st), fC(c), fPt(pt), fEta(eta), fPhi(phi), fNClTPC(ncltpc), fNClTPC1(ncltpc1), 
    fNClTPCShared(ncltpcs),fNClITS(nclits), fChi2TPC(chi2tpc), fChi2TPC1(chi2tpc1), fChi2ITS(chi2its), 
    fD(d), fZ(z), fDTPC(dtpc), fZTPC(ztpc), fSigma(sigma) {;}


  Double_t    Px()         const { return fPt*TMath::Cos(fPhi);  }
  Double_t    Py()         const { return fPt*TMath::Sin(fPhi);  }
  Double_t    Pz()         const { return fPt*TMath::SinH(fEta); }
  Bool_t      IsITSRefit() const { return (fSt&(ULong64_t)kITSrefit);    }
  Bool_t      IsTPCIn()    const { return (fSt&(ULong64_t)kTPCin);       }
  Bool_t      IsTPCRefit() const { return (fSt&(ULong64_t)kITSrefit);    }

public:
  ULong64_t   fSt;           //         status flag
  Double32_t  fC;            //[-1,1,2] charge
  Double32_t  fPt;           //[0,0,16] pt  from constrained fit
  Double32_t  fEta;          //[0,0,10] eta from constrained fit
  Double32_t  fPhi;          //[0,0,10] phi from constrained fit
  Short_t     fNClTPC;       //         number of clusters in TPC
  Short_t     fNClTPC1;      //         number of clusters in TPC if standalone
  Short_t     fNClTPCShared; //         number of shared clusters in TPC
  Char_t      fNClITS;       //         number of clusters in ITS
  Double32_t  fChi2TPC;      //[0,0,10] chi2 of TPC
  Double32_t  fChi2TPC1;     //[0,0,10] chi2 of TPC if standalone
  Double32_t  fChi2ITS;      //[0,0,10] chi2 of ITC
  Double32_t  fD;            //[0,0,16] transverse   DCA
  Double32_t  fZ;            //[0,0,16] longitudinal DCA
  Double32_t  fDTPC;         //[0,0,16] transverse   DCA TPC
  Double32_t  fZTPC;         //[0,0,16] longitudinal DCA TPC
  Double32_t  fSigma;        //[0,0,16] sigma to vertex

  ClassDef(MyPart,5) // My particle class in cylindrical coordinates
};

class MyTracklet : public TObject
{
public:
  MyTracklet(Double_t dphi=0, Double_t dth=0, Double_t eta=0, Double_t phi=0) :
    TObject(), fDPhi(dphi), fDTh(dth), fEta(eta), fPhi(phi) {;}

public:
  Double32_t  fDPhi; //[0,0,10] delta phi
  Double32_t  fDTh;  //[0,0,10] delta theta
  Double32_t  fEta;  //[0,0,10] eta
  Double32_t  fPhi;  //[0,0,10] phi

  ClassDef(MyTracklet,1) // My tracklet class
};

ULong_t flagValue[32] = 
  {
    MyPart::kITSin,
    MyPart::kITSout,
    MyPart::kITSrefit,
    MyPart::kITSpid,
    MyPart::kTPCin,
    MyPart::kTPCout,
    MyPart::kTPCrefit,
    MyPart::kTPCpid,
    MyPart::kTRDin,
    MyPart::kTRDout,
    MyPart::kTRDrefit,
    MyPart::kTRDpid,
    MyPart::kTOFin,
    MyPart::kTOFout,
    MyPart::kTOFrefit,
    MyPart::kTOFpid,
    MyPart::kHMPIDout,
    MyPart::kHMPIDpid,
    MyPart::kEMCALmatch,
    MyPart::kTRDbackup,
    0x0,     // 20 missing
    MyPart::kPHOSmatch,
    0x0,     // 22-24 missing
    0x0,
    0x0,
    MyPart::kMultInV0,
    MyPart::kMultSec,
    MyPart::kGlobalMerge,
    MyPart::kITSpureSA,
    MyPart::kTRDStop,
    MyPart::kESDpid,
    MyPart::kTIME
  };

TString flagLabel[32] =
  {
    "kITSin",
    "kITSout",
    "kITSrefit",
    "kITSpid",
    "kTPCin",
    "kTPCout",
    "kTPCrefit",
    "kTPCpid",
    "kTRDin",
    "kTRDout",
    "kTRDrefit",
    "kTRDpid",
    "kTOFin",
    "kTOFout",
    "kTOFrefit",
    "kTOFpid",
    "kHMPIDout",
    "kHMPIDpid",
    "kEMCALmatch",
    "kTRDbackup",
    "--",
    "kPHOSmatch",
    "--",
    "--",
    "--",
    "kMultInV0",
    "kMultSec",
    "kGlobalMerge",
    "kITSpureSA",
    "kTRDStop",
    "kESDpid",
    "kTIME"
  };
#endif
