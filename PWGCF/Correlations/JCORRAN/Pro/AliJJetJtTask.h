#ifndef ALIJJETJTTASK_H
#define ALIJJETJTTASK_H

/// \class AliJJetJtTask
/// \brief Wrapper class for AliJJetJtAnalysis
///
/// Provides an AliAnalysisTask wrapper for AliJJetJtAnalysis class.  
/// When the task is created an AliJJetJtAnalysis object is created and correct parameters for the analysis are set.
/// 
/// 
/// \author Tomas Snellman, tsnellma@cern.ch
/// \author Beomkyu Kim
/// \author Dongjo Kim
/// \date Aug 3, 2017

#include "AliAnalysisTaskSE.h"
#include "AliJJetJtAnalysis.h"
#include "AliAnalysisUtils.h"
#include "AliJJetTask.h"
#include "AliJCard.h"
#include "AliJJet.h"
#include "AliJRunTable.h"
#include "AliParticleContainer.h"
#include "AliAODMCParticle.h"

class AliMultSelection;
class AliJJetJtAnalysis;
class AliJEfficiency;
//==============================================================

using namespace std;

class AliJJetJtTask : public AliAnalysisTaskSE {

 public:
  AliJJetJtTask();
  AliJJetJtTask(const char *name,  TString inputformat);
  AliJJetJtTask(const AliJJetJtTask& ap);   
  AliJJetJtTask& operator = (const AliJJetJtTask& ap);
  virtual ~AliJJetJtTask();

  // methods to fill from AliAnalysisTaskSE
  virtual void UserCreateOutputObjects(); 
  virtual void Init();   
  virtual void LocalInit() { Init(); }
  virtual void UserExec(Option_t *option);
          bool IsGoodEvent(AliVEvent *event);
  virtual void Terminate(Option_t *);
  virtual void FinishTaskOutput();
  virtual Bool_t UserNotify() { std::cout<<"DEBUG UserNotify"<<std::endl; return kTRUE;}

  void SetJetTaskName(TString name){ fJetTaskName=name; }
  void SetMCJetTaskName(TString name){ fMCJetTaskName=name; }
  void SetSelector(TString name) {fSelector =name;}
  void SetJJetJtAnalysis(AliJJetJtAnalysis * jco){ fJJetJtAnalysis=jco; }
  void SetCard( AliJCard * card ){ fCard=card; }
  void SetMC(int mc) {fDoMC = mc;};
  void SetLog(bool doLog) {fDoLog = doLog;}
  void SetNrandom( int Nrand) { NRandom = Nrand;}
  void SetMoveJet( int move) { moveJet = move;}
  void SetLeadingJets(int leading){fLeadingJets = leading;}
  void SetMaxDeltaRCorr(double maxDeltaR){fmaxDeltaRCorr = maxDeltaR;}
  void SetCentCut(double cent){fCentCut = cent;}
  void FindDaughters(AliJJet * jet, AliAODMCParticle * track, AliMCParticleContainer * mcTracksCont);
  void SetSide(int side) {fSide = side;} // 0 = both sides, -1 A side, +1 C side

 private:

    // TODO new Task - AliJJetTask?
    AliJJetTask           * fJetTask; ///< Pointer to jet finder task
    AliJJetTask           * fMCJetTask; ///< Pointer to MC jet finder task (obsolete)
    TString                 fJetTaskName; ///< Name of the jet finder task
    TString                 fMCJetTaskName; ///< Name of the MC jet finder task
    TString                 fSelector; ///< Name of multiplicity selector
    AliJJetJtAnalysis     * fJJetJtAnalysis; ///< Pointer to the jT analysis class
    TClonesArray           *fJMCTracks; ///< List of MC tracks
    TDirectory            * fOutput; 
    AliJCard              * fCard; ///< Pointer to the configuration card
    Bool_t fFirstEvent; ///< True if this is the first event analyzed
    int cBin; ///< Comment needed
    int zBin; ///< Comment needed
    int NRandom; ///< Number of times a random background is generated for each track
    int moveJet; ///< Comment needed
    int fDoMC; ///< Whether or not MC analysis is performed
    int fSide; ///< Possible to use only one side for analysis, 0: both sides, -1: A side, +1, C side
    int fLeadingJets; ///< Do only leading jets if >0
    double  fmaxDeltaRCorr; ///<
    double fCentCut; ///<
    double zVert; ///< Vertex position
    bool fDoLog; ///< Whether or not logarithmic histograms should be filled/created
    AliAnalysisUtils *fAnaUtils;
    AliJRunTable *fRunTable;
    AliJTH1D  fEventHist; /// Histogram to store event numbers
    AliMultSelection*               sel=nullptr;//!

    ClassDef(AliJJetJtTask, 1); 
};
#endif // ALIJJETJTTASK_H
