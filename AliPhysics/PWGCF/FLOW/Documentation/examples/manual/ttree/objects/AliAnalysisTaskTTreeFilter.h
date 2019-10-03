#ifndef ALIANALYSISTTREEFILTER_H
#define ALIANALYSISTTREEFILTER_H

// AliRoot includes
#include "AliAnalysisTaskSE.h"

// forward declarations
class AliVEvent;
class AliVTrack;
class TTree;
class TClonesArray;
class AliFlowTTreeEvent;

class AliAnalysisTaskTTreeFilter : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskTTreeFilter();
  AliAnalysisTaskTTreeFilter(const char *name);

  virtual ~AliAnalysisTaskTTreeFilter();

  virtual void  UserCreateOutputObjects();
  virtual void  UserExec(Option_t *option);
  virtual void  Terminate(Option_t *); 

 private:
  Bool_t        ParseEvent(AliVEvent* event);
  void          ParseTracks(AliVEvent* event);
  void          PushToTTree();
  Bool_t        PassesCuts(AliVEvent* event);
  Bool_t        PassesCuts(AliVTrack* track);

  // Output objects
  TTree*                fTree;              //! output data
  AliFlowTTreeEvent*    fEvent;             //! custom event
  TClonesArray*         fTrackArray;        //! custom tracks

  ClassDef(AliAnalysisTaskTTreeFilter, 1);
};

#endif
