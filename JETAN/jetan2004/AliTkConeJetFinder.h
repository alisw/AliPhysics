// $Id$

/*************************************************************************
 * AliTkConeJetFinder.h                                                     *
 * Thorsten Kollegger <kollegge@ikf.physik.uni-frankfurt.de>             *
 * Jet finder based on a cone algorithm with seeds and addition of       *
 * midpoints, for a description of the algorithm see hep-ex/0005012      *
 ************************************************************************/

/*************************************************************************
 * Some remarks:                                                         *
 * ---                                                                   *
 * This is a prototype finder - it's clearly not optimized for speed     *
 * e.g. one could win much by using TClonesArrays instead of TObjArrays  *
 * ---                                                                   *
 * for MC studies, it keeps track of the particle/proto-jet/jet origin   *
 ************************************************************************/

// NEXT STEPS: 09/03/02 12:30am
// Implement list of found stable protojets
// Implement function isProtoJetAlreadyFound()
// -> Finishes findProtoJets
// Implement merge/split process
// Implement seed + midpoint addition

// includes
#include <Riostream.h>
#include "TROOT.h"
#include "TObject.h"
#include "TH2.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "TParticle.h"


//forward declarations
class AliTkTower;
class AliTkProtoJet;

class AliTkConeJetFinder : public TObject {
 public:
  AliTkConeJetFinder();
  ~AliTkConeJetFinder();

  void setDefaultSettings();
  void Init();

  void InitEvent();

  void FillEtHistFromTParticles(TClonesArray *particles);
  Bool_t isTParticleAccepted(TParticle * particle);
  void AddEtHist(Float_t eta, Float_t phi, Double_t Et);
  
  
  Int_t run();

  /*
  Int_t CreateSeedList();
  Int_t CreateMidPoints();
  */
  Int_t FindProtoJets();
  Int_t CalculateConeCentroid(Float_t *eta, Float_t *phi, Float_t *Et);
  
  Int_t FindJets();
  /*
  Bool_t isProtoJetConeShared();
  Int_t FindNeighborCone(Int_t nProtoJet);
  Int_t SplitProtoJetCone(Int_t nProtoJet1, Int_t nProtoJet2);
  Int_t MergeProtoJetCone(Int_t nProtoJet1, Int_t nProtoJet2);
  */

  // eta binning
  void setEtaNBins(Int_t nbins);
  Int_t getEtaNBins();
  void setEtaRange(Float_t min, Float_t max);
  Float_t getEtaRangeMin();
  Float_t getEtaRangeMax();
  void setEtaGrid(Int_t nbins, Float_t min, Float_t max);

  // phi binnig
  void setPhiNBins(Int_t nbins);
  Int_t getPhiNBins();
  void setPhiRange(Float_t min, Float_t max);
  Float_t getPhiRangeMin();
  Float_t getPhiRangeMax();
  void setPhiGrid(Int_t nbins, Float_t min, Float_t max);

  // jet radius
  void setJetConeRadius(Float_t r);
  Float_t getJetConeRadius();
  /*
  // seed tower threshold
  void setTowerSeedEt(Float_t minEt);
  Float_t getTowerSeedEt();
  */

  // protojetlist
  TObjArray *getProtoJetList();

 protected:

 private:
  //----------------------------------------------------------------------
  // member variables
  TH2D *hEtHist;
  TH2D *hTowerEt;
  Int_t nEtaBins;
  Float_t nEtaMin;
  Float_t nEtaMax;
  Int_t nPhiBins;
  Float_t nPhiMin;
  Float_t nPhiMax;

  Float_t jetRadius;
  TH2D *hEtConeHist;
  TH2D *hStableProtoJetHist;

  // number of "towers" - calcualted in Init()
  Int_t nTower;
  // array of towers
  AliTkTower *towers;
  // seed list - array of ints with tower id's
  Int_t *SeedTower;
  // will not work this way - this requires more towers than tower due to the additon of midpoint

  // protojet list
  TObjArray *protojets;

  // private functions
  Float_t calcPhiDiff(Float_t phi1, Float_t phi2);
  Bool_t isProtoJetStable(Float_t etaDiff, Float_t phiDiff);
  Bool_t isProtoJetInTower(Int_t etaBin,Int_t phiBin, 
			   Float_t eta, Float_t phi);
  Bool_t isProtoJetInTower(Int_t tower, Float_t eta, Float_t phi);

  Int_t findTower(Float_t eta, Float_t phi);

  ClassDef(AliTkConeJetFinder,1)
};

class AliTkTower : public TObject {
 public:
  AliTkTower() : TObject() { etaMin = 0; etaMax = 0; phiMin = 0; phiMax = 0;}
  ~AliTkTower() { }
  Int_t uid;
  Float_t etaMin;
  Float_t etaMax;
  Float_t phiMin;
  Float_t phiMax;
  Float_t Et;
 private:

  ClassDef(AliTkTower,1)
};

class AliTkProtoJet : public TObject {
 public:
  AliTkProtoJet() : TObject() { }
  ~AliTkProtoJet() { }
  Float_t eta;
  Float_t phi;
  Float_t Et;

  Bool_t IsEqual(AliTkProtoJet *other);
 private:


  ClassDef(AliTkProtoJet,1)
};

struct SProtoJet;
struct SProtoJet {
  Float_t eta;
  Float_t phi;
  Float_t Et;
  bool operator==(const SProtoJet &s1);
};

