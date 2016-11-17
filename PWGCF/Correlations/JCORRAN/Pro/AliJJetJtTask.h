#ifndef ALIJJETJTTASK_H
#define ALIJJETJTTASK_H

/// \class AliJJetJtTask
/// \brief Wrapper class for AliJJetJtAnalysis
///
/// This is a longer description of the class. This longer description is
/// formatted using the Markdown syntax (see below) and can span on multiple
/// lines.
/// 
/// \author Tomas Snellman, tsnellma@cern.ch
/// \author Beomkyu Kim
/// \author Dongjo Kim
/// \date Nov 11, 2016

#include "AliAnalysisTaskSE.h"
#include "AliJJetJtAnalysis.h"
#include "AliAnalysisUtils.h"
#include "AliJJetTask.h"
#include "AliJCard.h"
#include "AliJJet.h"
#include "AliJRunTable.h"
#include "AliParticleContainer.h"
#include "AliAODMCParticle.h"

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
  void SetJJetJtAnalysis(AliJJetJtAnalysis * jco){ fJJetJtAnalysis=jco; }
  void SetCard( AliJCard * card ){ fCard=card; }
  void SetMC(int mc) {fDoMC = mc;};
  void SetNrandom( int Nrand) { NRandom = Nrand;}
  void SetMoveJet( int move) { moveJet = move;}
  void FindDaughters(AliJJet * jet, AliAODMCParticle * track, AliMCParticleContainer * mcTracksCont);

 private:

    // TODO new Task - AliJJetTask?
    AliJJetTask           * fJetTask; ///< Pointer to jet finder task
    AliJJetTask           * fMCJetTask; ///< Pointer to MC jet finder task (obsolete)
    TString                 fJetTaskName; ///< Name of the jet finder task
    TString                 fMCJetTaskName; ///< Name of the MC jet finder task
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
    double zVert; ///< Vertex position
    AliAnalysisUtils *fAnaUtils;
    AliJRunTable *fRunTable;
    TH1D * fEventHist;

    ClassDef(AliJJetJtTask, 1); 
};
#endif // ALIJJETJTTASK_H
