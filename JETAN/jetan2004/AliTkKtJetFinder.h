// $Id$

#ifndef ALITKKTJETFINDER_H
#define ALITKKTJETFINDER_H

#include <TObject.h>
#include <Riostream.h>
#include <list>

struct TPreCluster;
struct TPreCluster {
  Int_t id;
  Float_t E;
  Float_t px;
  Float_t py;
  Float_t pz;
  TPreCluster *next;
};

struct Tdlist;
struct Tdlist {
  Float_t d;
  Int_t id1;
  Int_t id2;
  Tdlist *prev;
  Tdlist *next;
  bool operator==(const Tdlist &p1) const {return (d==p1.d);}
  bool operator<(const Tdlist &p1) const {return (d > p1.d);}
};

struct TJetList;
struct TJetList {
  Float_t E;
  Float_t px;
  Float_t py;
  Float_t pz;
  TJetList *next;
};

class AliTkKtJetFinder : public TObject {
 public:
  AliTkKtJetFinder() : TObject() {status = 0;}
  ~AliTkKtJetFinder() {}
  
  // run functions - these ones should be called by the user...
  void init();
  void makeTParticles(TClonesArray *particles);
  void clear();
  void finish();

  // main components of the jet finder - normally called by make()
  // public to allow easier testing/timing analysis in macro
  void initEvent();
  void preclusterTParticles(TClonesArray *particles);
  void findJets();

  // precluster options
  void setPhiMin(Float_t fPhiMin) {phiMin = fPhiMin;}
  void setPhiMax(Float_t fPhiMax) {phiMax = fPhiMax;}
  void setNPhiBins(Int_t nPhiBins) {phiBins = nPhiBins;}
  void setPhiBins(Int_t nPhiBins, Float_t fPhiMin, Float_t fPhiMax);

  void setThetaMin(Float_t fThetaMin) {thetaMin = fThetaMin;}
  void setThetaMax(Float_t fThetaMax) {thetaMax = fThetaMax;}
  void setNThetaBins(Int_t nThetaBins) {thetaBins = nThetaBins;}
  void setThetaBins(Int_t nThetaBins, Float_t fThetaMin, Float_t fThetaMax);
  
  // jet finder options
  void setD(Float_t fD) {finder_D = fD;}
  void setDCut(Float_t fDCut) {finder_DCut = fDCut;}

  // some "standard" options
  void setDefaultOptions();

  // debug options
  void setDebugLevel(Int_t nDebugLevel);
  void setDebugFilename(Char_t *sFilename) {debugFilename = sFilename;}

 private:
  Int_t status;

  // precluster parameters
  Float_t phiMin;
  Float_t phiMax;
  Int_t phiBins;

  Float_t thetaMin;
  Float_t thetaMax;
  Int_t thetaBins;

  // jet finder parameters
  Float_t finder_D;
  Float_t finder_DCut;

  // debug parameters
  Int_t debugLevel;
  Char_t *debugFilename;

  // data variables
  // list of preclusters with pointer to last precluster
  Int_t preClusterUID;
  TPreCluster *firstPreCluster;
  TPreCluster *lastPreCluster;

  // array of pointers to preclusters
  // makes access faster
  TPreCluster **preClusterArray;

  // functions related to preclusters...
  void addPreCluster(TPreCluster *precluster);
  void deletePreCluster(Int_t UID);
  void dumpPreClusters();
  void dumpPreClusterArray();

  // double linked sorted list of (relativ) transverse momenta 
  list<Tdlist *> myDList;
  //priority_queue<Tdlist *> myNewDList;
  //priority_queue<Tdlist *> myHeap;
  // priority queues are not known to CINT... :-(

  // functions related to (relativ) transverse momenta list...
  void addD(Tdlist *newD);
  void buildNewDList();
  Float_t calcD(TPreCluster *p1, TPreCluster *p2 = NULL);

  // list of final jets
  TJetList *firstJet;

  // private member functions
  Bool_t isTParticleAccepted(TParticle *particle);
  
  // debug/output functions
  void DebugOutput(Char_t *output);
  void TimingOutput(Char_t *output);
  void DetailedOutput(Char_t *output);
  
  ClassDef(AliTkKtJetFinder,1)
};
#endif

