// $Id$

// Version 2 of my Cone Jet Finder
// improved storage (see also version 1)

#ifndef ALITKCONJETFINDERV2_H
#define ALITKCONJETFINDERV2_H

#include <Riostream.h>
#include <TParticle.h>
#include <list>

// helper classes
class tower {

 public:
  tower();
  tower(const tower& t);
  tower(Float_t phimin, Float_t phimax, Float_t etamin, Float_t etamax);
  ~tower();

  tower& operator+=(const Float_t E);
  tower& operator+=(const TParticle *part);

  Float_t getEtaMin() const {return fEta_min;}
  Float_t getEtaMax() const {return fEta_max;}
  Float_t getEta()    const {return fEta_center;}
  Float_t getPhiMin() const {return fPhi_min;}
  Float_t getPhiMax() const {return fPhi_max;}
  Float_t getPhi() const {return fPhi_center;}
  Float_t getEt()  const {return fEt;}
  Float_t getEtaWidth() const {return (fEta_max-fEta_min);}
  Float_t getPhiWidth() const {return (fPhi_max-fPhi_min);}

  void setEtaMin(Float_t f) {fEta_min=f;}
  void setEtaMax(Float_t f) {fEta_max=f;}
  void setEta(Float_t f)    {fEta_center=f;}
  void setPhiMin(Float_t f) {fPhi_min=f;}
  void setPhiMax(Float_t f) {fPhi_max=f;}
  void setPhi(Float_t f)    {fPhi_center=f;}

  void addParticle(const TParticle* part) {fParticles.push_back(part);}
  void clearParticles() {fParticles.erase(fParticles.begin(),fParticles.end());}
  void clear() {fEt = 0; clearParticles();}

  list<const TParticle*> *getParticles();

 private:
  Float_t fEta_min;
  Float_t fEta_max;
  Float_t fEta_center;
  Float_t fPhi_min;
  Float_t fPhi_max;
  Float_t fPhi_center;
  Float_t fEt;

  list<const TParticle*> fParticles;
};

ostream& operator<<(ostream& s,const tower& t);

//-----------------------------------------------------

#include "AliTkEtaPhiVector.h"

class protojet {

 public:
  protojet();
  protojet(const protojet& p);
  protojet(const protojet *p);
  ~protojet(){eraseTowers();}

  Float_t Eta() const {return fCentroid.Eta();}
  Float_t Phi() const {return fCentroid.Phi();}
  Float_t getEt() const {
    if(fUpdate) cerr << "Eta update needed" << endl;
    return fEt;
  }  
  Float_t getEt() {
    if(fUpdate) update();
    return fEt;
  }

  Int_t getNTowers() const { return fTowers.size(); }

  void setCentroidPosition(AliTkEtaPhiVector center) {fCentroid = center;}

  AliTkEtaPhiVector getCentroidPosition()   const {
    return fCentroid; 
  }
  AliTkEtaPhiVector getEtWeightedPosition() const {
    if(fUpdate) cerr << "Weighted update needed" << endl;
    return fEtWeightedCentroid;
  }
  AliTkEtaPhiVector getEtWeightedPosition() {
    if(fUpdate) update();
    return fEtWeightedCentroid;
  }

  void addTower(tower *tower) {fTowers.push_back(tower);fUpdate=kTRUE;}
  void eraseTowers() {fTowers.erase(fTowers.begin(),fTowers.end()); }
  void eraseTower(tower *pTower) {fTowers.remove(pTower);}
  void clear(){eraseTowers();}

  list<tower*> getTowerList() const;
  list<tower*> getSharedTowerList(const protojet *other) const;

  bool hasTower(tower *pTower) const;
  bool shareTowers(const protojet *other) const;
  Float_t diffToCenter(tower *pTower) const;
  bool operator<(const protojet &p1) const {return (getEt() < p1.getEt());}
  bool operator==(const protojet &p1);
  void update();

 private:
  AliTkEtaPhiVector fCentroid;
  AliTkEtaPhiVector fEtWeightedCentroid;
  Float_t fEt;
  Bool_t fUpdate;
  list<tower*> fTowers;
};

ostream& operator<<(ostream& s,const protojet& p);

//-----------------------------------------------------
// Real ConeJetFinderClass
//-----------------------------------------------------

#include <TObject.h>
#include <TClonesArray.h>
#include <TTree.h>
#include <vector>
#include "AliTkConeJet.h"
#include "AliTkConeJetEvent.h"
#ifdef ALICEINTERFACE
#include <AliJetEventParticles.h>
#endif

//if defined produce some debugging histos
//#define DOHISTOS

class AliTkConeJetFinderV2 : public TObject {
 public:
  // constructor
  AliTkConeJetFinderV2();
  virtual ~AliTkConeJetFinderV2();
  
  // run control
  void defaultSettings();
  void setSettings(Int_t phibins,Int_t etabins);
  void setEtMinJet(Float_t et){fEtMinJet=et;} //minimum jet energy required
  void setEtCut(Float_t et){fEtCut=et;}       //min et for seedpoints     
  void setPtCut(Float_t pt){fPtCut=pt;}       //set pt cut 
  void setOutput(Bool_t out){fOutput=out;}
  void setRadius(Float_t r=0.7) {fRadius=r;}

  void init();
  void initEvent(const TClonesArray *particles,Int_t type);
  void initEvent(const TClonesArray *particles,Int_t type,TString desc); 
  void run();
  void finishEvent();
  void finish();

  // real physics functions...
  void createTowers();
  void fillTowersFromTParticles(const TClonesArray *particles);
#ifdef ALICEINTERFACE
  void fillTowersFromAliParticles(const TClonesArray *particles);
  void initEvent(const AliJetEventParticles *p,TString desc);
#endif

  void createSeedPoints();
  void findProtojets();
  void findJets();

#ifdef DOHISTOS
  // analysis functions...
  void createHistos();
  void createEventHistos();
  void clearEventHistos();
  void writeEventHistos();
  void writeHistos();
#endif

  // evout function
  void setEvOutFilename(const Char_t *filename);
  Char_t *getEvOutFilename() {return fEvout_name;}

  void addProtojet(protojet *pj);
  void dumpProtojets(Float_t etmin = 5.0);

  void splitMergeJets(protojet *jet1, protojet *jet2);
  void splitJets(protojet *jet1,protojet *jet2);
  void mergeJets(protojet *jet1,protojet *jet2);
  void addJet(protojet *pj) {fJets.push_back(pj);}
  void dumpJets();

  Bool_t isJetEnergy(Float_t min,Float_t max);
  Float_t maxJetEnergy(); //of protojets (before calling finishEvent)
  AliTkConeJetEvent* getEvent() const {return fEvoutevent;}

  //getters
  Bool_t getOutput()    const {return fOutput;}
  Int_t getNTowers()    const {return fNTowers;}
  Int_t getNEtaBins()   const {return fEtaBins;}
  Float_t getEtaMin()   const {return fEtaMin;}
  Float_t getEtaMax()   const {return fEtaMax;}
  Float_t getEtaWidth() const {return fEtaWidth;}
  Int_t getPhiBins()    const {return fPhiBins;}
  Float_t getPhiMin()   const {return fPhiMin;}
  Float_t getPhiMax()   const {return fPhiMax;}
  Float_t getPhiWidth() const {return fPhiWidth;}
  Float_t getEtCut()    const {return fEtCut;};  
  Float_t getPtCut()    const {return fPtCut;};  
  Float_t getEtMinJet() const {return fEtMinJet;}
  Float_t getRadius()   const {return fRadius;}

 protected:
  void initEvent_(const TClonesArray *particles,Int_t type);
  void fillTowerHist();
  Int_t findTower(Float_t phi, Float_t eta);
  Bool_t isTParticleAccepted(TParticle *particle);

  Bool_t fOutput; //outpug some information if true
  Int_t fNTowers; //# of towers
  Int_t fEtaBins; //eta-phi grid
  Float_t fEtaMin;
  Float_t fEtaMax;
  Float_t fEtaWidth;
  Int_t fPhiBins;
  Float_t fPhiMin;
  Float_t fPhiMax;
  Float_t fPhiWidth;
  Float_t fEtCut;    // cut for seedpoints
  Float_t fEtMinJet; // cut for jets
  Float_t fPtCut; // cut for particles
  Float_t fRadius;

  //towers, seeds, protojets and jets
  vector<tower> *fTowers; //!
  list<AliTkEtaPhiVector> fSeedPointsNew; //!
  list<protojet*> fProtojets; //!
  list<protojet*> fJets; //!

  //store results
  TFile *fEvoutfile; //!
  Char_t *fEvout_name; //!
  AliTkConeJetEvent *fEvoutevent; //!
  TTree *fEvouttree; //!

#ifdef ALICEINTERFACE
  TClonesArray *fAliParticles; //!
#endif

#ifdef DOHISTOS
  //histograms
  TObjArray fHistos;
  TObjArray fEventHistos;
  TFile *fHistFile; //!
#endif

  ClassDef(AliTkConeJetFinderV2,1)
};

inline Bool_t AliTkConeJetFinderV2::AliTkConeJetFinderV2::isTParticleAccepted(TParticle * /* particle */) 
{
  // check if particle is accepted
  // makes sense to write this into a own class, but now I'm lazy

  // check if particle is stable -> a detectable particle
#ifndef ALICEINTERFACE 
  if (particle->GetStatusCode()%100 != 1) {
    return kFALSE;
  }
#endif

  // default case: accept
  return kTRUE;
}

#endif
