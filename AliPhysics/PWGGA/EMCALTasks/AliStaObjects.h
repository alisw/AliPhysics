#ifndef AliStaObjects_H
#define AliStaObjects_H

// $Id$

#include <TLorentzVector.h>

class AliStaHeader
{
 public:
  AliStaHeader() : fRun(0), fOrbit(0), fPeriod(0), fBx(0), fL0(0), fL1(0), fL2(0),
                   fTrClassMask(0), fTrCluster(0), fOffTriggers(0), fFiredTriggers(),
                   fTcls(0), fV0And(0), fIsHT(0), fIsPileup(0), fIsPileup2(0), fIsPileup4(0), fIsPileup8(0), 
                   fNSpdVertices(0), fNTpcVertices(0), fV0Cent(0), fV0(0), fCl1Cent(0), fCl1(0), fTrCent(0), 
                   fTr(0), fCqual(-1), fPsi(0), fPsiRes(0), fNSelTr(0), fNSelPrimTr(0), fNSelPrimTr1(0),
                   fNSelPrimTr2(0), fNCells(0), fNCells0(0), fNCells01(0), fNCells03(0), 
                   fNCells1(0), fNCells2(0), fNCells5(0), fNClus(0), fNClus1(0), fNClus2(0), fNClus5(0), 
                   fMaxCellE(0), fMaxClusE(0), fMaxTrE(0), fNcSM0(0), fNcSM1(0), fNcSM2(0), fNcSM3(0), 
                   fNcSM4(0), fNcSM5(0), fNcSM6(0),fNcSM7(0),fNcSM8(0),fNcSM9(0) {;}
  virtual ~AliStaHeader() {;}

  ULong64_t     GetEventId() const {
                  return (((ULong64_t)fPeriod << 36) |
                          ((ULong64_t)fOrbit  << 12) |
                          (ULong64_t)fBx); 
                }

 public:
  Int_t         fRun;            //         run number
  UInt_t        fOrbit;          //         orbit number
  UInt_t        fPeriod;         //         period number
  UShort_t      fBx;             //         bunch crossing id
  UInt_t        fL0;             //         l0 trigger bits
  UInt_t        fL1;             //         l1 trigger bits
  UShort_t      fL2;             //         l2 trigger bits
  ULong64_t     fTrClassMask;    //         trigger class mask
  UChar_t       fTrCluster;      //         trigger cluster mask
  UInt_t        fOffTriggers;    //         fired offline triggers for this event
  TString       fFiredTriggers;  //         string with fired triggers
  UInt_t        fTcls;           //         custom trigger definition
  Bool_t        fV0And;          //         if V0AND (from AliTriggerAnalysis)
  Bool_t        fIsHT;           //         if EMCAL L0 (from AliTriggerAnalysis)
  Bool_t        fIsPileup;       //         indicate pileup from IsPileupFromSPD with 0.8 minzdist
  Bool_t        fIsPileup2;      //         indicate pileup from IsPileupFromSPD with 0.4 minzdist
  Bool_t        fIsPileup4;      //         indicate pileup from IsPileupFromSPD with 0.2 minzdist
  Bool_t        fIsPileup8;      //         indicate pileup from IsPileupFromSPD with 0.1 minzdist
  UShort_t      fNSpdVertices;   //         number of pileup vertices (spd)
  UShort_t      fNTpcVertices;   //         number of pileup vertices (tpc)
  Double32_t    fV0Cent;         //[0,0,16] v0 cent
  Double32_t    fV0;             //[0,0,16] v0 result used for cent 
  Double32_t    fCl1Cent;        //[0,0,16] cl1 cent
  Double32_t    fCl1;            //[0,0,16] cl1 result used for cent 
  Double32_t    fTrCent;         //[0,0,16] tr cent
  Double32_t    fTr;             //[0,0,16] tr result used for cent 
  Int_t         fCqual;          //         centrality quality
  Double32_t    fPsi;            //[0,0,16] event-plane angle
  Double32_t    fPsiRes;         //[0,0,16] event-plane ange resolution
  UShort_t      fNSelTr;         //         # selected tracks         
  UShort_t      fNSelPrimTr;     //         # selected tracks (primary)
  UShort_t      fNSelPrimTr1;    //         # selected tracks (primary) pt > 1 GeV/c
  UShort_t      fNSelPrimTr2;    //         # selected tracks (primary) pt > 2 GeV/c
  UShort_t      fNCells;         //         # cells
  UShort_t      fNCells0;        //         # cells > 0.45 GeV
  UShort_t      fNCells01;       //         # cells > 0.1  GeV
  UShort_t      fNCells03;       //         # cells > 0.3  GeV
  UShort_t      fNCells1;        //         # cells > 1    GeV
  UShort_t      fNCells2;        //         # cells > 2    GeV
  UShort_t      fNCells5;        //         # cells > 5    GeV
  UShort_t      fNClus;          //         # clus
  UShort_t      fNClus1;         //         # clus > 1 GeV
  UShort_t      fNClus2;         //         # clus > 2 GeV
  UShort_t      fNClus5;         //         # clus > 5 GeV
  Double32_t    fMaxCellE;       //[0,0,16] maximum cell energy
  Double32_t    fMaxClusE;       //[0,0,16] maximum clus energy
  Double32_t    fMaxTrE;         //[0,0,16] maximum trigger energy
  UShort_t      fNcSM0;          //         # cells > 0.1  GeV in SM 0
  UShort_t      fNcSM1;          //         # cells > 0.1  GeV in SM 1
  UShort_t      fNcSM2;          //         # cells > 0.1  GeV in SM 2
  UShort_t      fNcSM3;          //         # cells > 0.1  GeV in SM 3
  UShort_t      fNcSM4;          //         # cells > 0.1  GeV in SM 4
  UShort_t      fNcSM5;          //         # cells > 0.1  GeV in SM 5
  UShort_t      fNcSM6;          //         # cells > 0.1  GeV in SM 6
  UShort_t      fNcSM7;          //         # cells > 0.1  GeV in SM 7
  UShort_t      fNcSM8;          //         # cells > 0.1  GeV in SM 8
  UShort_t      fNcSM9;          //         # cells > 0.1  GeV in SM 9

  ClassDef(AliStaHeader,6) // Header class
};

class AliStaVertex
{
 public:
  AliStaVertex(Double_t x=0, Double_t y=0, Double_t z=0) : fVx(x), fVy(y), fVz(z), fVc(-1), fDisp(0), fZres(0),
                                                           fChi2(0), fSt(0), fIs3D(0), fIsZ(0) {;}
  virtual ~AliStaVertex() {;}

 public:
  Double_t      fVx;          //[0,0,16] vertex x
  Double_t      fVy;          //[0,0,16] vertex y
  Double_t      fVz;          //[0,0,16] vertex z
  Double_t      fVc;          //[0,0,16] number of contributors to vertex
  Double_t      fDisp;        //[0,0,16] dispersion
  Double_t      fZres;        //[0,0,16] z-resolution
  Double_t      fChi2;        //[0,0,16] chi2 of fit
  Bool_t        fSt;          //         status bit
  Bool_t        fIs3D;        //         is vertex from 3D
  Bool_t        fIsZ;         //         is vertex from Z only

  ClassDef(AliStaVertex,1) // Vertex class
};

class AliStaCluster : public TObject
{
 public:
  AliStaCluster() : TObject(), 
                    fE(0), fR(0), fEta(0), fPhi(0), fN(0), fN1(0), fN3(0), fIdMax(-1), fSM(-1), fEmax(0), fE2max(0), fEcross(0),
                    fTmax(0), fDbc(-1), fDisp(-1), fM20(-1), fM02(-1), fEcc(-1), fSig(-1), fSigEtaEta(-1), fSigPhiPhi(-1),
                    fIsTrackM(0), fTrDz(0), fTrDr(-1), fTrEp(0), fTrDedx(0), fTrIso(0), fTrIso1(0), fTrIso2(0),  
                    fTrIsoD1(0), fTrIso1D1(0), fTrIso2D1(0), fTrIsoD3(0), fTrIso1D3(0), fTrIso2D3(0),
                    fTrIsoD4(0), fTrIso1D4(0), fTrIso2D4(0), fTrIsoStrip(0), fCeIso(0), fCeIso1(0), 
                    fCeIso3(0), fCeIso4(0), fCeIso3x3(0), fCeIso4x4(0), fCeIso5x5(0), fCeCore(0), fCeIso3x22(0), 
                    fIsShared(0), fTrigId(-1), fTrigE(0), fMcLabel(-1), fEmbE(0) {;}

  void          GetMom(TLorentzVector& p, Double_t *vertex=0);
  void          GetMom(TLorentzVector& p, AliStaVertex *vertex);

 public:
  Double32_t    fE;                //[0,0,16] energy
  Double32_t    fR;                //[0,0,16] radius (cylinder)
  Double32_t    fEta;              //[0,0,16] eta
  Double32_t    fPhi;              //[0,0,16] phi
  UChar_t       fN;                //         number of cells
  UChar_t       fN1;               //         number of cells > 100 MeV
  UChar_t       fN3;               //         number of cells > 300 MeV
  Short_t       fIdMax;            //         id maximum cell
  Char_t        fSM;               //         super module number (from maximum cell)
  Double32_t    fEmax;             //[0,0,16] energy of maximum cell
  Double32_t    fE2max;            //[0,0,16] energy of second maximum cell
  Double32_t    fEcross;           //[0,0,16] energy of the 4 adjacent cells around the seed
  Double32_t    fTmax;             //[0,0,16] time of maximum cell
  Double32_t    fDbc;              //[0,0,16] distance to nearest bad channel
  Double32_t    fDisp;             //[0,0,16] cluster dispersion, for shape analysis
  Double32_t    fM20;              //[0,0,16] 2-nd moment along the main eigen axis
  Double32_t    fM02;              //[0,0,16] 2-nd moment along the second eigen axis
  Double32_t    fEcc;              //[0,0,16] eccentricity
  Double32_t    fSig;              //[0,0,16] sigma
  Double32_t    fSigEtaEta;        //[0,0,16] sigma eta-eta
  Double32_t    fSigPhiPhi;        //[0,0,16] sigma phi-phi
  Bool_t        fIsTrackM;         //         if true then track values are set
  Double32_t    fTrDz;             //[0,0,16] dZ to nearest track
  Double32_t    fTrDr;             //[0,0,16] dR to nearest track (in x,y)
  Double32_t    fTrEp;             //[0,0,16] E/P to nearest track 
  Double32_t    fTrDedx;           //[0,0,16] dE/dx (TPC signal) to nearest track 
  Double32_t    fTrIso;            //[0,0,16] track isolation
  Double32_t    fTrIso1;           //[0,0,16] track isolation (pt>1GeV/c)
  Double32_t    fTrIso2;           //[0,0,16] track isolation (pt>2GeV/c)
  Double32_t    fTrIsoD1;          //[0,0,16] track isolation, iso dist 0.25
  Double32_t    fTrIso1D1;         //[0,0,16] track isolation (pt>1GeV/c), iso dist 0.1
  Double32_t    fTrIso2D1;         //[0,0,16] track isolation (pt>2GeV/c), iso dist 0.1
  Double32_t    fTrIsoD3;          //[0,0,16] track isolation, iso dist 0.3
  Double32_t    fTrIso1D3;         //[0,0,16] track isolation (pt>1GeV/c), iso dist 0.3
  Double32_t    fTrIso2D3;         //[0,0,16] track isolation (pt>2GeV/c), iso dist 0.3
  Double32_t    fTrIsoD4;          //[0,0,16] track isolation, iso dist 0.4
  Double32_t    fTrIso1D4;         //[0,0,16] track isolation (pt>1GeV/c), iso dist 0.4
  Double32_t    fTrIso2D4;         //[0,0,16] track isolation (pt>2GeV/c), iso dist 0.4
  Double32_t    fTrIsoStrip;       //[0,0,16] track isolation strip, dEtaXdPhi=0.015x0.3
  Double32_t    fCeIso;            //[0,0,16] cell isolation in R=0.20
  Double32_t    fCeIso1;           //[0,0,16] cell isolation in R=0.10
  Double32_t    fCeIso3;           //[0,0,16] cell isolation in R=0.30
  Double32_t    fCeIso4;           //[0,0,16] cell isolation in R=0.40
  Double32_t    fCeIso3x3;         //[0,0,16] cell isolation in 3x3 cells
  Double32_t    fCeIso4x4;         //[0,0,16] cell isolation in 4x4 cells
  Double32_t    fCeIso5x5;         //[0,0,16] cell isolation in 5x5 cells
  Double32_t    fCeCore;           //[0,0,16] cell content in R=0.05 
  Double32_t    fCeIso3x22;        //[0,0,16] cell isolation in rectangular strip of dEtaXdPhi=0.042x0.308
  Bool_t        fIsShared;         //         =true then extends across more than one super module
  Short_t       fTrigId;           //         index of matched trigger tower
  Double32_t    fTrigE;            //[0,0,16] energy (FEE) of matched trigger tower
  Short_t       fMcLabel;          //         index of closest MC particle
  Double32_t    fEmbE;             //[0,0,16] sum of energy of embedded (MC) cells in cluster

  ClassDef(AliStaCluster,10) // Cluster class
};

class AliStaTrigger : public TObject
{
 public:
  AliStaTrigger() : TObject(), fE(0), fEta(0), fPhi(0), fIdMax(-1) {}

 public:
  Double32_t    fE;                //[0,0,16] energy
  Double32_t    fEta;              //[0,0,16] eta
  Double32_t    fPhi;              //[0,0,16] phi
  Short_t       fIdMax;            //         id maximum cell

  ClassDef(AliStaTrigger,2) // Trigger class
};

class AliStaPart : public TObject
{
 public:
  AliStaPart() : TObject(), fPt(0), fEta(0), fPhi(0), fVR(0), fVEta(0), fVPhi(0), fPid(0), fMo(-1), fDet(-2), 
                 fLab(-1), fNs(0) { memset(fDs,-1,sizeof(Short_t)*99); }

  Int_t         OnEmcal() const { return (fDet==8);  }
  Int_t         IsSim()   const { return (fDet!=-2); }
    
 public:
  Double32_t    fPt;               //[0,0,16] pt
  Double32_t    fEta;              //[0,0,16] eta
  Double32_t    fPhi;              //[0,0,16] phi
  Double32_t    fVR;               //[0,0,16] prod r (cylinder)
  Double32_t    fVEta;             //[0,0,16] prod eta
  Double32_t    fVPhi;             //[0,0,16] prod phi
  Short_t       fPid;              //         pid
  Short_t       fMo;               //         index of mother
  Short_t       fDet;              //         detector in which particle left trace (8 for EMCAL, see AliTrackReference.h)
    // the following must be filled before first usage
  Short_t       fLab;              //!        label (index in array)
  Short_t       fNs;               //!        number of daughters
  Short_t       fDs[99];           //!        daughters

  ClassDef(AliStaPart,1) // Particle class
};
#endif
