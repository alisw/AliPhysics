/// \class AliAnalysisTaskCutStudies

#ifndef AliAnalysisTaskCutStudies_H
#define AliAnalysisTaskCutStudies_H

#include "AliAnalysisTaskMKBase.h"
#include "Hist.h"

class AliESDtrackCuts;
class AliVEvent;
class AliESDEvent;
class AliAODEvent;
class AliMCEvent;
class AliStack;
class AliHeader;
class AliGenEventHeader;
class AliESDtrack;
class AliMCParticle;

class AliAnalysisTaskCutStudies : public AliAnalysisTaskMKBase
{
public:
  AliAnalysisTaskCutStudies();
  AliAnalysisTaskCutStudies(const char *name);
  virtual ~AliAnalysisTaskCutStudies();

  virtual void AddOutput();                     //called at the beginning
  virtual Bool_t IsEventSelected();             //called for each event
  virtual void AnaEvent();                      //called once for every selected event
  virtual void AnaTrack(Int_t flag = 0);        //called once for every track in DATA+MC event
  virtual void AnaTrackMC(Int_t flag = 0);      //called once for every track in DATA event
  virtual void AnaParticleMC(Int_t flag = 0);   //called once for every track in MC event

  static AliAnalysisTaskCutStudies* AddTaskCutStudies(const char* name = "TaskCutStudies", const char* outfile = 0);
  
private:
  typedef AnalysisHelpers::Hist<THnSparseI> Hist;
  
  AliAnalysisTaskCutStudies(const AliAnalysisTaskCutStudies&); // not implemented
  AliAnalysisTaskCutStudies& operator=(const AliAnalysisTaskCutStudies&); // not implemented
        
  Hist myHist;
  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskCutStudies, 1);
  /// \endcond
};

#endif
