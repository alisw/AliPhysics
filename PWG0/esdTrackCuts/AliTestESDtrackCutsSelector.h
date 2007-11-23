/* $Id$ */

#ifndef ALITESTESDTRACKCUTSSELECTOR_H
#define ALITESTESDTRACKCUTSSELECTOR_H

#include "AliSelectorRL.h"

class AliESDtrackCuts;
class TH1F;
class TH3F;

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
    AliESDtrackCuts*  fEsdTrackCutsNoVtx;  // all cuts except vtx

    AliESDtrackCuts*  fEsdTrackCutsPri;  // cuts for tracks from primary particles
    AliESDtrackCuts*  fEsdTrackCutsSec;  // cuts for tracks from secondary particles
    AliESDtrackCuts*  fEsdTrackCutsPlusZ;  // cuts for tracks that go to z > 0
    AliESDtrackCuts*  fEsdTrackCutsMinusZ;  // cuts for tracks that go to z < 0
    AliESDtrackCuts*  fEsdTrackCutsPos;  // cuts for tracks from positive particles
    AliESDtrackCuts*  fEsdTrackCutsNeg;  // cuts for tracks from negative particles

    TH1F*             fPIDAfterCutNoVtx;      // true PID of tracks that passed all cuts except vtx
    TH1F*             fPIDAfterCutAll;        // true PID of tracks that passed all cuts incl. vtx

    TH3F*             fVertex;                // originating vertex of specific particles


 private:
    AliTestESDtrackCutsSelector(const AliTestESDtrackCutsSelector&);
    AliTestESDtrackCutsSelector& operator=(const AliTestESDtrackCutsSelector&);

  ClassDef(AliTestESDtrackCutsSelector, 0);
};

#endif
