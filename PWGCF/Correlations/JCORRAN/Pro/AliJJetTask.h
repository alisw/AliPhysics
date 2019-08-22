#ifndef ALIJJETTASK_H
#define ALIJJETTASK_H
class TH1;
class TH1D;
class TH2;
class TH3;
class TGraphErrors;
class TProfile;
class AliJetContainer;
class AliParticleContainer;
class AliClusterContainer;


#include "AliAnalysisTaskEmcalJet.h"
#include <TClonesArray.h>
#include <TString.h>


class AliJJetTask : public AliAnalysisTaskEmcalJet {
 public:

  enum { kJUndefined = -1 , kJRecoTrack, kJMCParticle };
  AliJJetTask();
  AliJJetTask(const char *name, const int nJetFinder); // FIXME: Is camelCase good for args?
  virtual ~AliJJetTask();
  
  void UserCreateOutputObjects();

  void SetTrackArrayName( char *c ) { fTrackArrayName = c; }
  vector<TClonesArray> *GetAliJJetCollection() {return &fJJets;}
  TObjArray*            GetAliJJetList(int i) {return &fJJets[i]; }
  TClonesArray*         GetJTracks()  {return &fJTracks;}
  TClonesArray*         GetMCJTracks()  {return &fJMCTracks;}
  int GetNumberOfJetCollections() {return fNJetFinder;}

  int GetTaskEntry()        {return fTaskEntry;}
  void                        Terminate(Option_t *option);

  void SetTrackOrMCParticle( UInt_t i, int v ){ fTrackOrMCParticle[i] = v; }
  int  GetTrackOrMCParticle( UInt_t i ){ return fTrackOrMCParticle.at( i ); }
  void SetDebug(int n) {debug = n; }
  int  GetDebug(){ return debug; }
  void SetMC(int mc) {fIsMC = mc;}
  void SetnR(int nR) {fnR = nR;}
  void SetACside(int flag) {fACside = flag;}
  void SetIncludeFullJets(int full) {fDoFullJets = full;}
  int GetIncludeFullJets() {return fDoFullJets;}
  int GetnR() {return fnR;}
  void SetConeSize(UInt_t i, double radius) {fConeSizes[i] = radius;}
  double  GetConeSize( UInt_t i ){ return fConeSizes.at( i ); }
  void Setnkt(int nkt) {fnkt = nkt;}
  int Getnkt() {return fnkt;}
  vector<TString> &GetJetFinderString() { return fJetFinderString;}

 protected:
  void                        ExecOnce();
  Bool_t                      FillHistograms()   ;
  Bool_t                      Run()              ;
  void                        CheckClusTrackMatching();

  vector<AliJetContainer*>            fJetsCont;          //!Jets

 private:
  AliJJetTask(const AliJJetTask& ap);            // not implemented
  AliJJetTask &operator=(const AliJJetTask&); // not implemented

  //==== FOR AliJ<Object>
  // FIXME: We assume that we have only one of each. be carefull. This must be fixed later
  // FIXME: What for CaloCluster? Not needed?
  TClonesArray              fJTracks;     //! tracks array
  TClonesArray              fJMCTracks;   //! mc tracks array
  TClonesArray              fJClusters;   //! Clusters array

  vector<TClonesArray>      fJJets;       //! Jets array

  UInt_t                    fTaskEntry;
  vector<int>               fTrackOrMCParticle;
  vector<double>            fConeSizes;

  vector<TString>           fJetFinderString;

  TString    fTrackArrayName; // track constituents array name

  Int_t fNJetFinder;
  Int_t debug;
  Int_t fIsMC;
  Int_t fDoFullJets;
  Int_t fnR;
  Int_t fACside;
  Int_t fnkt;

  ClassDef(AliJJetTask, 5) 
    // 5 : add fTrackOrMCParticle
    // 5 : del fTracksCont
    // 5 : del fCaloClusterCont
    // 5 : Add fJClusters
};


#endif
