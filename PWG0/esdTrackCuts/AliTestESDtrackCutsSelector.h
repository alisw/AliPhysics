/* $Id$ */

#ifndef ALITESTESDTRACKCUTSSELECTOR_H
#define ALITESTESDTRACKCUTSSELECTOR_H

#include "AliSelectorRL.h"

class AliESDtrackCuts;

class AliTestESDtrackCutsSelector : public AliSelectorRL {
  public:
    AliTestESDtrackCutsSelector();
    virtual ~AliTestESDtrackCutsSelector();

    virtual void    Begin(TTree* tree);
    virtual void    SlaveBegin(TTree *tree);
    virtual void    Init(TTree *tree);
    virtual Bool_t  Process(Long64_t entry);
    virtual void    SlaveTerminate();
    virtual void    Terminate();

 protected:
    void ReadUserObjects(TTree* tree);

    AliESDtrackCuts*  fEsdTrackCutsAll;  // esd track cuts for all tracks   

    AliESDtrackCuts*  fEsdTrackCutsPri;  // cuts for tracks from primary particles
    AliESDtrackCuts*  fEsdTrackCutsSec;  // cuts for tracks from secondary particles
    AliESDtrackCuts*  fEsdTrackCutsPlusZ;  // cuts for tracks that go to z > 0
    AliESDtrackCuts*  fEsdTrackCutsMinusZ;  // cuts for tracks that go to z < 0
    AliESDtrackCuts*  fEsdTrackCutsPos;  // cuts for tracks from positive particles
    AliESDtrackCuts*  fEsdTrackCutsNeg;  // cuts for tracks from negative particles

 private:
    AliTestESDtrackCutsSelector(const AliTestESDtrackCutsSelector&);
    AliTestESDtrackCutsSelector& operator=(const AliTestESDtrackCutsSelector&);

  ClassDef(AliTestESDtrackCutsSelector, 0);
};

#endif
