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

  AliJJetTask();
  AliJJetTask(const char *name, const int NJetFinder);
  //AliJJetTask(const char * name, int njetfinder );
  virtual ~AliJJetTask();

  
  void                        UserCreateOutputObjects();
  void SetTrackArrayName( char *c ) { fTrackArrayName = c; }
  vector<TClonesArray>* GetAliJJetCollection() {return &fJJets;}
  TObjArray*           GetAliJJetList(int i) {return &fJJets[i]; }
  TClonesArray*        GetJTracks()  {return &fJTracks;}
  int GetNumberOfJetCollections() {return fNJetFinder;}
  //void SetNumberOfJetCollections(const int n) {
  //    fJetsCont.resize(n); 
  //    fJJets.resize(n);
  //    fTracksCont.resize(n); 
  //    fCaloClustersCont.resize(n); 
  //    fNJetFinder=n;
  //}

  int GetTaskEntry()        {return fTaskEntry;}
  void                        Terminate(Option_t *option);

  void SetDebug(int n) {debug = n; }
  void SetMC(int mc) {fIsMC = mc;} 
  int  GetDebug(){ return debug; }
  vector<TString> &GetJetFinderString() { return fJetFinderString;}


 protected:
  void                        ExecOnce();
  Bool_t                      FillHistograms()   ;
  Bool_t                      Run()              ;
  void                        CheckClusTrackMatching();






  vector<AliJetContainer*>            fJetsCont;                   //!Jets
  //AliJetContainer            **fJetsConts;              //!Jets
  //Int_t 		     nJetsConts;
  vector<AliParticleContainer*>       fTracksCont;                 //!Tracks
  vector<AliClusterContainer*>        fCaloClustersCont;           //!Clusters  

 private:
  AliJJetTask(const AliJJetTask& ap);            // not implemented
  AliJJetTask &operator=(const AliJJetTask&); // not implemented

  //TVector *EtaGapThresholds;
  //TVector *RGapThresholds;
  //TVector *KlongBorders;
  //TVector *XeBorders;

  TClonesArray              fJTracks;  //! tracks array
  vector<TClonesArray>      fJJets;  //! tracks array
  UInt_t                    fTaskEntry;

  vector<TString> fJetFinderString;



  TString    fTrackArrayName; // track constituents array name


  Int_t fNJetFinder;
  Int_t debug;
  Int_t fIsMC;

  ClassDef(AliJJetTask, 4) // jet sample analysis task



};


#endif
