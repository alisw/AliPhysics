#ifndef ALIFEMTOK0ANALYSIS_H
#define ALIFEMTOK0ANALYSIS_H

//
// Class AliFemtoK0Analysis
//
// AliFemtoK0Analysis
// author:
//        Matthew Steinpreis (matthew.steinpreis@cern.ch)
//

class TH1F;
class TH1D;
class TH2D;
class TH3D;
class TProfile;
class TRandom3;

class AliESDEvent;
class AliAODEvent;
class AliESDtrackCuts;
class AliESDpid;

#include "AliAnalysisTask.h"
#include "AliAnalysisTaskSE.h"
#include "AliFemtoK0EventCollection.h"
#include "AliAODpidUtil.h"
#include "AliESDpid.h"

class AliFemtoK0Analysis : public AliAnalysisTaskSE {
 public:
  AliFemtoK0Analysis();
  AliFemtoK0Analysis(const char *name);
  virtual ~AliFemtoK0Analysis();
  AliFemtoK0Analysis(const AliFemtoK0Analysis&);
  AliFemtoK0Analysis& operator=(const AliFemtoK0Analysis&);

 private:
  
  virtual void   UserCreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);  

  void MyInit();
  void GetGlobalPositionAtGlobalRadiiThroughTPC(const AliESDtrack *track, const Float_t bfield, Float_t globalPositionsAtRadii[9][3]);
  
  enum 
  {
    kCentBins    = 16,
    kZVertexBins = 10,
    kEventsToMix =  5,
    kMultLimit   = 300,              //maximum number of v0s, array size
 
    ncthetabins = 36,
    nphibins    = 72
  };

  int fEventCount;

  AliFemtoK0EventCollection ***fEC; //!
  AliFemtoK0Event *fEvt; //!

  TRandom3* fRandomNumber; //!
  
  const char     *fName;
  AliAODEvent    *fAOD; //!    // AOD object
  TList          *fOutputList; //! Compact Output list
  AliAODpidUtil  *fPidAOD; //!
  AliESDpid      *fPidESD; //!
  AliESDtrack    *fPosDaughter1;//!
  AliESDtrack    *fPosDaughter2;//!
  AliESDtrack    *fNegDaughter1;//!
  AliESDtrack    *fNegDaughter2;//!
  
  ClassDef(AliFemtoK0Analysis, 1); 
};

#endif
