// $Id$

#ifndef ALITKCHARGEDJETFINDER_H
#define ALITKCHARGEDJETFINDER_H

//if defined use Torstens version
//#define DOCHARGED

//if defined produce some debugging histos
#ifdef DOCHARGED
//#define DOHISTOS
#endif

#include <Riostream.h>
#include <list>
#include <map>

#include "AliTkEtaPhiVector.h"
#ifdef DOCHARGED
#include "TkChargedJet.h"
#else
#include "AliTkConeJetEvent.h"
#endif

#ifdef ALICEINTERFACE
#include <AliJetEventParticles.h>
#endif

class TTree;
class TFile;
class TClonesArray;
class TParticle;
class TObjArray;

//------------------------------------------------------------------------
// Faster Helper Class (charged jets)
//-------------------------------------------------------------------------
class jet {

 public:
  jet();

  Float_t Eta() const {return fEta;}
  Float_t Phi() const {return fPhi;}
  AliTkEtaPhiVector getCentroid()    const;

  Double_t getPt()                const;
  Int_t  getNParticles()          const; 
  TParticle *getParticle(Int_t i) const;
  TClonesArray *getParticles()    const;
  Double_t getDiff(TParticle *particle)   const;
  Double_t getDiffSq(TParticle *particle) const;

  void addParticle(TParticle *particle);
  friend ostream& operator<<(ostream& s, jet& j);

 private:
  list<TParticle *> fParticles;
  Double_t fPt;
  Double_t fEta;
  Double_t fPhi;
  Int_t fNParticles;
};

//------------------------------------------------------------------------
// Finder Class (charged jets)
//-------------------------------------------------------------------------

class AliTkChargedJetFinder : public TObject {

 public:
  AliTkChargedJetFinder();
  virtual ~AliTkChargedJetFinder();

  void defaultSettings();
  void setSettings(Int_t /*phibins*/,Int_t /*etabins*/) {;}
  void init();
  void initEvent(TClonesArray *newParticles,Int_t type = 1);
#ifdef ALICEINTERFACE
  void initEvent(const AliJetEventParticles *p,TString desc);
#endif
  void run();
  void finishEvent();
  void finish();

  void setEtMinJet(Float_t et)  { fMinJetPt=et; }  //minimum jet energy required
  void setEtCut(Float_t et)     { setPtSeed(et); } //min et for seedpoints     
  void setPtCut(Float_t pt)     { fPtCut=pt; }     //set pt cut 
  void setOutput(Bool_t out)    { fOutput=out; }
  void setRadius(Float_t r=0.7) { setFinderR(r); }

  void setFinderR(Float_t r) { fR=r;fRSq = r*r; }
  Float_t getFinderR() const { return fR; }
  void setEtaMin(Float_t eta) { fEtaMin = eta; }
  Float_t getEtaMin() const { return fEtaMin; }
  void setEtaMax(Float_t eta) { fEtaMax = eta; }
  Float_t getEtaMax() const { return fEtaMax; }
  void setPtSeed(Float_t s) { fPtSeed=s;}
  Float_t getPtSeed() const { return fPtSeed; }
  void setMinJetPt(Float_t s) { fMinJetPt=s;}
  Float_t getMinJetPt() const { return fMinJetPt; }
  
  void setEvOutFilename(const Char_t *filename);
  const Char_t *getEvOutFilename() { return fEvout_name; }
#ifdef DOHISTOS
  void setHistFilename(const Char_t *filename);
  const Char_t *getHistFilename() { return fOutput_name; }
#endif

  Bool_t  getOutput()   const {return fOutput;}
  Int_t   getNTowers()  const {return 0;}
  Int_t   getNEtaBins() const {return 0;}
  Float_t getEtaWidth() const {return 0;}
  Int_t   getPhiBins()  const {return 0;}
  Float_t getPhiMin()   const {return 0;} 
  Float_t getPhiMax()   const {return 2 * TMath::Pi();}
  Float_t getPhiWidth() const {return 0;}
  Float_t getEtCut()    const {return fPtSeed;};  
  Float_t getPtCut()    const {return fPtCut;};  
  Float_t getEtMinJet() const {return fMinJetPt;}
  Float_t getRadius()   const {return fR;}

 private:
  Int_t fOutput;
  list<TParticle *> fParticles; //!
  list<jet> fJets; //!

  // parameter for jetfinder
  Float_t fR;
  Float_t fRSq;
  Float_t fEtaMin;
  Float_t fEtaMax;
  Float_t fPtCut;
  Float_t fPtSeed;
  Float_t fMinJetPt;
  
#ifdef DOCHARGED
  TkChargedJet *fMyTJet; //!
#else
  AliTkConeJetEvent *fEvoutevent; //!
#endif

#ifdef DOHISTOS
  // histograms...
  Char_t *fOutput_name; //!
  TObjArray *fHistos; //!
  TFile *fHistoutfile; //!
#endif

  Char_t *fEvout_name; //!
  TFile *fEvoutfile; //!
  TTree *fEvouttree; //!

#ifdef ALICEINTERFACE
  TClonesArray *fAliParticles; //!
#endif

  bool isTParticleAccepted(TParticle *particle);
  void addParticle(TParticle *particle);
  list<jet>::iterator findHighestJet();
  void checkJets();

  ClassDef(AliTkChargedJetFinder,2)
};
//-------------------------------------------------------------------------

#endif
