#ifndef CLINFO_H
#define CLINFO_H

#include <TMath.h>
#include <TVector3.h>

class CLInfo {
 public:
 CLInfo() : fTrig(0), fEvAcc(0), fSpdPu1(0), fSpdPu2(0), fSpdPu(0), fRun(0), fVz(0), fVcon(0), 
            fNtracks(0), fMult(0), fSkippedV(0), fSkippedC(0) {;}
  virtual ~CLInfo() {;}
  UInt_t        fTrig;         // trigger bits
  Bool_t        fEvAcc;         // if true, event accepted by standard EventCuts
  Bool_t        fSpdPu1;       // if true, spd pileup with low settings
  Bool_t        fSpdPu2;       // if true, spd pileup with high settings
  Bool_t        fSpdPu;        // if true, spd pileup with mult settings
  Int_t         fRun;          // run number
  Double32_t    fVz;           // [-32,32,16] vertex z
  Short_t       fVcon;         // number of contributors to vertex
  Short_t       fNtracks;      // number of aod tracks
  Short_t       fMult;         // mult in -0.5<eta<0.5
  Int_t         fSkippedV;     // how many events were skipped before the current one due to vertex
  Int_t         fSkippedC;     // how many events were skipped before the current one due to cuts
  ClassDef(CLInfo,1) // ClInfo (header) class
};

class CLPart : public TObject {
 public:
  CLPart() : TObject(), fEta(0), fPhi(0), fPt(0), fId(0) {;}
  virtual ~CLPart() {;}
  TVector3      mom() const { TVector3 r; r.SetPtEtaPhi(fPt,fEta,fPhi); return r; }
  Double32_t    fEta;        // [-1,1,10] eta
  Double32_t    fPhi;        // [0,6.2,10] phi
  Double32_t    fPt;         // [0,250,12] pT
  Int_t         fId;         // identifier
  ClassDef(CLPart,1) // CLPart (particle) class
};

class CLTrack : public CLPart {
 public:
  CLTrack() : CLPart(), fCh(0), fIsGlobal(0), fIsHybGlo(0), fIsTcon(0), fIsPrim(0), fIsTpcOk(0), fIsTofOk(0), fTsig(0), fNsd(0), fM2tof(0), fTfrac(0), fDcar(0), fDcaz(0) {;}
  virtual ~CLTrack() {;}
  Short_t       fCh;         // charge
  Bool_t        fIsGlobal;   // is global track
  Bool_t        fIsHybGlo;   // is hybrid-global track
  Bool_t        fIsTcon;     // is TPC constrained track
  Bool_t        fIsPrim;     // is primary candidate
  Bool_t        fIsTpcOk;    // is TPC pid ok
  Bool_t        fIsTofOk;    // is TOF pid ok
  Double32_t    fTsig;       //[0,0,16] dEdx signal
  Double32_t    fNsd;        //[0,99,16] number of sigma from deuteron dEdx
  Double32_t    fM2tof;      //[0,99,16] m2 tof
  Double32_t    fTfrac;      //[0,1,8] fraction of TPC clusters
  Double32_t    fDcar;       //[0,0,12] dca in transverse direction 
  Double32_t    fDcaz;       //[0,0,12] dca in beam direction 
  ClassDef(CLTrack,1) // CLTrack (track) class
};

class CLMCpart : public CLPart {
 public:
  CLMCpart() : CLPart(), fIsPrim(0) {;}
  virtual ~CLMCpart() {;}
  Bool_t        fIsPrim;      // primary particle flag
  ClassDef(CLMCpart,1) // CLMCTrack (mc track) class
};

class CLClus : public CLPart {
 public:
  CLClus() : CLPart(), fM02(0), fNmax(0), fM12(0), fE1(0), fE2(0) {;}
  virtual ~CLClus() {;}
  Double32_t         fM02;    //[0,0,16] long axis value
  Short_t            fNmax;   // number of maxima
  Double32_t         fM12;    //[0,10,12] inv mass of leading clusters
  Double32_t         fE1;     //[0,200,10] leading energy 1
  Double32_t         fE2;     //[0,200,10] leading energy 2
  ClassDef(CLClus,1) // CLClus (cluster) class
};
#endif

#ifndef ALICLTASK_H
#define ALICLTASK_H

class TClonesArray;
class TString;
class TH1;
class TH1F;
class TH2F;
class TH3F;
class TNtuple;
class TNtupleD;
class TTree;
class AliAODTrack;
class AliPIDResponse;
class AliCaloPID;
class AliCalorimeterUtils;
class AliEventCuts;

#include <AliAnalysisTaskSE.h>

class AliCLTask : public AliAnalysisTaskSE {
 public:
  AliCLTask(const char *name=0);
  virtual ~AliCLTask();
  void                        SetDoClust(Bool_t b)     { fDoClust = b; }      
  void                        SetDoSkip(Bool_t b)      { fDoSkip = b; }      
  void                        SetMCMode(Bool_t mc)     { fMcMode = mc; }
  void                        SetPtCut(Double_t ptcut) { fPtCut = ptcut; }
  void                        SetM2Cut(Double_t m2cut) { fM2Cut = m2cut; }      
 protected:
  void                        UserCreateOutputObjects();
  void                        UserExec(Option_t */*option*/="");
  void                        Terminate(Option_t */*option*/="") {;}
  Double_t                    GetM2tof(AliAODTrack *track) const;
  Bool_t                      fMcMode;    //  if true then only accept events with at least one deuteron
  Bool_t                      fDoClust;   //  if true the do clusters
  Bool_t                      fDoSkip;    //  if true then skip non-selected events
  Double_t                    fPtCut;     //  only store tracks/clusters above ptcut
  Double_t                    fM2Cut;     //  only store deuteron candidates above m2cut
  Int_t                       fCounter;   //! count skipped events due to vertex
  Int_t                       fCounterC;  //! count skipped events due to cuts
  AliPIDResponse             *fPidRes;    //! PID response
  AliCaloPID                 *fPidCalo;   //! PID calo 
  AliCalorimeterUtils        *fCaloUtils; //! Calo utils
  TTree                      *fMyTree;    //! tree to store
  CLInfo                     *fMyInfo;    //! info class
  TClonesArray               *fMyTracks;  //! tracks
  TClonesArray               *fMyClus;    //! clus
  TClonesArray               *fMyParts;   //! mc particles
  AliEventCuts               *fEventCuts; //! event cuts
 private:
  AliCLTask(const AliCLTask&);            // not implemented
  AliCLTask &operator=(const AliCLTask&); // not implemented
  ClassDef(AliCLTask, 1) // Constantin's Task
};
#endif
