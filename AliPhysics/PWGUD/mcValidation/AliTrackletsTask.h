/* $Id:$ */

// -----------------------------------------------
// Task to extract distributions 
// for traclets paramters
// for a quick comparison 
// between MC and data
// -----------------------------------------------

#ifndef ALITRACKLETSTASK_H
#define ALITRACKLETSTASK_H

#include "AliAnalysisTask.h"

class TH1D;
class TH2D;
class TH3D;
class TH1I;
class AliESDEvent;

class AliTrackletsTask : public AliAnalysisTask {
  public:
    AliTrackletsTask();
    virtual ~AliTrackletsTask();

    virtual void   ConnectInputData(Option_t *);
    virtual void   CreateOutputObjects();
    virtual void   Exec(Option_t*);
    virtual void   Terminate(Option_t*);

 protected:
    AliESDEvent *fESD;      //! ESD object
    TList* fOutput;         //! list send on output slot 0
    TH1I* fNtracks;         //! nunmber of tracks
    TH1D* fPhi;             //! phi distribution
    TH2D* fEtaPhi;          //! phi vs eta distribution
    TH1D* fDeltaPhi;        //! deltaPhi distribution
    TH1D* fDeltaTheta;      //! deltaTheta distribution
    TH1D* fVtxX;            //! x of the SPD vertex distribution
    TH1D* fVtxY;            //! y of the SPD vertex distribution
    TH1D* fVtxZ;            //! z of the SPD vertex distribution
    TH3D* fVtx;             //! SPD vertex distribution
    TH3D* fVtxContributors; //! SPD vertex distribution with N contributors > 0

 private:
    AliTrackletsTask(const AliTrackletsTask&);
    AliTrackletsTask& operator=(const AliTrackletsTask&);

  ClassDef(AliTrackletsTask, 1);
};

#endif
