#ifndef Cumulants_h
#define Cumulants_h

#include <vector>
#include <TNamed.h>
#include <TList.h>
#include <TObjArray.h>
#include <TParticle.h>
#include <TComplex.h>
#include <AliVParticle.h>

class CPart;
class VPart;
class TH1;

class Cumulants : public TNamed {
 public:
  Cumulants(const char *name="cumulants", Int_t mbins=350, Int_t minM=10);
  virtual ~Cumulants() {;}
  void                  AddQC4withEG(Double_t etagap);
  void                  EnableEG();
  void                  EnableQC();
  void                  EnableQC4with4NL(Int_t mn=50,  Double_t etamin=0.0);
  void                  EnableQC4with3NL(Int_t mn=100, Double_t etamin=0.0);
  void                  EnableQC4withEG();
  TList                *GetList() const { return fList; }
  Int_t                 GetM()    const { return fM; }
  void                  RunAll();
  void                  SetTracks(TObjArray &trks, Bool_t doKinCuts=1);
  void                  SetKine(Double_t etamin, Double_t etamax, Double_t ptmin, Double_t ptmax) { fEtaMin=etamin; fEtaMax=etamax; fPtMin=ptmin; fPtMax=ptmax; }
  void                  SetDebug(Bool_t b) { fDoDebug = b; }
  void                  SetPrint(Bool_t b) { fDoPrint = b; }
 protected:
  void                  RunEG();
  void                  RunQC();  
  void                  RunQC4with4NL();
  void                  RunQC4with3NL();
  void                  RunQC4withEG();  
  void                  RunQC4withEG(Double_t etagap, Int_t &nn, Int_t &np, Double_t &val2, Double_t &val4);
  Int_t                 fCumMBins;    //number of M bins
  Int_t                 fMinM;        //minimum number of M
  Double_t              fEtaMin;      //min eta cut
  Double_t              fEtaMax;      //max eta cut
  Double_t              fPtMin;       //min pt cut
  Double_t              fPtMax;       //max pt cut
  Bool_t                fDoEtaGap;    //=true do eta gap method
  Bool_t                fDoCharge;    //=true do charge
  Bool_t                fDoQC;        //=true do QC method
  Bool_t                fDoQC44;      //=true do QC with 4 nested loops
  Int_t                 fMaxNL4;      //maximum M for which 4NL will be executed
  Double_t              fEGminNL4;    //minimum eta gap between particles in 4NL
  Bool_t                fDoQC43;      //=true do QC with 3 nested loops
  Int_t                 fMaxNL3;      //maximum M for which 3NL will be executed
  Double_t              fEGminNL3;    //minimum eta gap between particles in 3NL
  Bool_t                fDoQC4withEG; //=true do QC with 2 eta gaps
  Bool_t                fDoDebug;     //=true do debug statements
  Bool_t                fDoPrint;     //=true do print statements
  std::vector<Double_t> fEGCuts;      //eta gap cuts for two particle correlations
  std::vector<Double_t> fEGQCCuts;    //eta gap cuts for QC4
  std::vector<Double_t> fEGC2;        //!eta gap c2
  std::vector<Double_t> fEGC3;        //!eta gap c3
  std::vector<Double_t> fEGS2;        //!eta gap s2
  std::vector<Double_t> fEGS3;        //!eta gap s3
  std::vector<Int_t>    fEGCounts;    //!eta counts per gap
  TComplex              fQC[7];       //!QC values
  Int_t                 fM;           //!number of particles
  std::vector<CPart>    fParts;        //!vector with particles
  TList                *fList;        //!list with histograms
  TH1                  *fHists[999];  //!histograms
  ClassDef(Cumulants,1) //Cumulant class (CL)
};

class CPart {
 public:
  CPart() : fPt(0), fEta(0), fPhi(0), fCharge(0) {;}
  CPart(Double_t pt, Double_t eta, Double_t phi, Short_t ch) : fPt(pt), fEta(eta), fPhi(TVector2::Phi_0_2pi(phi)), fCharge(ch) {;}
  CPart(const CPart &p) : fPt(p.fPt), fEta(p.fEta), fPhi(p.fPhi), fCharge(p.fCharge) {;}
  CPart(const AliVParticle &p) : fPt(p.Pt()), fEta(p.Eta()), fPhi(TVector2::Phi_0_2pi(p.Phi())), fCharge(p.Charge()) {;}
  CPart(const TParticle &p);
  Double_t Pt()         const {return fPt;}
  Double_t Eta()        const {return fEta;}
  Double_t Phi()        const {return fPhi;}
  Short_t  Charge()     const {return fCharge;}
 protected:
  Double_t fPt;     //pt
  Double_t fEta;    //eta
  Double_t fPhi;    //phi
  Short_t  fCharge; //charge
};

class VPart: public AliVParticle {
 public:
  VPart() : AliVParticle(), fMom(), fCharge(0) {;}
  VPart(const VPart &p) : AliVParticle(p), fMom(p.fMom), fCharge(p.fCharge) {;} 
  VPart& operator=(const VPart &p) { return *this;}
  virtual ~VPart() {;}
  virtual Double_t Px() const {return 0;}
  virtual Double_t Py() const {return 0;}
  virtual Double_t Pz() const {return 0;}
  virtual Double_t Pt() const {return 0;}
  virtual Double_t P()  const {return 0;}
  virtual Bool_t   PxPyPz(Double_t p[3]) const {return 0;};
  virtual Double_t Xv() const {return 0;};
  virtual Double_t Yv() const {return 0;};
  virtual Double_t Zv() const {return 0;};
  virtual Bool_t   XvYvZv(Double_t x[3]) const {return 0;};  
  virtual Double_t OneOverPt()  const {return 0;}
  virtual Double_t Phi()        const {return 0;}
  virtual Double_t Theta()      const {return 0;}
  virtual Double_t E()          const {return 0;}
  virtual Double_t M()          const {return 0;}
  virtual Double_t Eta()        const {return 0;}
  virtual Double_t Y()          const {return 0;}
  virtual Short_t Charge()      const {return 0;}
  virtual Int_t   GetLabel()    const {return 0;}
  virtual Int_t   PdgCode()     const {return 0;}       
  virtual const Double_t *PID() const {return 0;};
 protected:
  TVector3 fMom;    //momentum
  Short_t  fCharge; //charge
  ClassDef(VPart, 1) // Virtual particle (CL)
};

#endif
