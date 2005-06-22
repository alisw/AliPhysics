// $Id$

#ifndef ALIJFTOWERH
#define ALIJFTOWERH

#include <TMath.h>

class AliJFTower 
{
 public:
  AliJFTower(){};
  virtual ~AliJFTower(){};

  /*
  AliJFTower(const tower& t);
  AliJFTower(Float_t phimin, Float_t phimax, Float_t etamin, Float_t etamax);

  AliJFTower& operator+=(const Float_t E);
  */

 protected:
  Float_t fEtaMin;
  Float_t fEtaMax;
  Float_t fEtaCenter;
  Float_t fPhiMin;
  Float_t fPhiMax;
  Float_t fPhiCenter;
  Float_t fEt;

  ClassDef(AliJFTower,1)
};

#endif /*JFTOWERH*/

#if 0
// includes
#include "TROOT.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TFile.h"

#include <vector>
#include <list>
#include <map>

// classes
class tower {
 public:
  tower();
  tower(const tower& t);
  tower(Float_t phimin, Float_t phimax, Float_t etamin, Float_t etamax);
  Float_t eta_min;
  Float_t eta_max;
  Float_t eta_center;
  Float_t phi_min;
  Float_t phi_max;
  Float_t phi_center;
  Float_t Et;
  tower& operator+=(const Float_t E);
};

ostream& operator<<(ostream& s,const tower& t);

class seedpoint {
 public:
  Float_t eta;
  Float_t phi;
};

ostream& operator<<(ostream& s,const seedpoint& t);

class protojet {
 public:
  protojet() { eta = -999; phi = -999; Et = -999; }
  protojet(const protojet& p) { eta=p.eta; phi=p.phi; Et=p.Et; }
  protojet(const protojet * p) { eta=p->eta; phi=p->phi; Et=p->Et; }
  Float_t eta;
  Float_t phi;
  Float_t Et;
  bool operator<(const protojet &p1) const {return (Et < p1.Et);}
  bool operator==(const protojet &p1);
};

ostream& operator<<(ostream& s,const protojet& p);

class TkConeJetFinderV2 : public TObject {
 public:
  // some test stuff - I delete it later...
  void test();
  void testAddTowers();
  void testAddProtojet();
  void testAssocProtojetTower();

  // constructor
  TkConeJetFinderV2() : TObject() {}
  
  // run control
  void defaultSettings();
  void init();
  void initEvent(TClonesArray *particles,Int_t type = 1);
  void run();
  void finishEvent();
  void finish();

  // real physics functions...
  void createTowers(Int_t nPhiTower,Float_t phiMin,Float_t phiMax,
		    Int_t nEtaTower,Float_t etaMin,Float_t etaMax);
  void fillTowersFromTParticles(TClonesArray *particles);
  void createSeedPoints(Float_t EtCut);
  void findProtojets(Float_t radius);
  void findJets();

  // analysis functions...
  void createHistos();
  void createEventHistos();
  void clearEventHistos();
  void writeEventHistos();
  void writeHistos();

 protected:
 private:
  // containers to save info...
  vector<tower> *towers;
  Int_t mEtaBins;
  Float_t mEtaMin;
  Float_t mEtaMax;
  Float_t mEtaWidth;
  Int_t mPhiBins;
  Float_t mPhiMin;
  Float_t mPhiMax;
  Float_t mPhiWidth;

  //vector<tower>::iterator 
  Int_t findTower(Float_t phi, Float_t eta);
  Bool_t isTParticleAccepted(TParticle *particle);
  
  list<seedpoint> seedPoints;
  Float_t mEtCut;

  list<protojet> protojets;
  
  multimap<protojet,tower> assocTowers;
  list<protojet> jets;

  TObjArray mHistos;
  TObjArray mEventHistos;
  TFile *histFile;
  void fillTowerHist();

  ClassDef(TkConeJetFinderV2,1)
};
#endif


