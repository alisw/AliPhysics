#ifndef ALITRDRECOTASK_H
#define ALITRDRECOTASK_H


// Author: Alexandru Bercuci, 10/09/2008

#ifndef ALIANALYSISTASK_H
#include "AliAnalysisTask.h"
#endif

class TObjArray;
class TTreeSRedirector;
class AliTRDrecoTask : public AliAnalysisTask 
{
public:
  enum AliTRDrecoSteeringBits{
    kMCdata       = BIT(20)
    ,kFriends     = BIT(21)
    ,kPostProcess = BIT(22)
  };
  AliTRDrecoTask(const char *name, const char *title);
  virtual ~AliTRDrecoTask();
  
  
  void           ConnectInputData(Option_t *);
  TObjArray*     Container() const {return fContainer;}
  virtual void   CreateOutputObjects() = 0;
  virtual void   Exec(Option_t *) = 0;

  Int_t          GetDebugLevel() const { return fDebugLevel;}
  Int_t          GetNRefFigures() const { return fNRefFigures; } 
  virtual void   GetRefFigure(Int_t ifig, Int_t &first, Int_t &last, Option_t *opt="lp");

  Bool_t         HasFriends() const {return TestBit(kFriends);};
  Bool_t         HasMCdata() const {return TestBit(kMCdata);};
  Bool_t         HasPostProcess() const {return TestBit(kPostProcess);};

  virtual Bool_t Load(const Char_t *filename);
  virtual Bool_t PostProcess();
  virtual void   SetDebugLevel(Int_t level);
  virtual void   SetFriends(Bool_t fr = kTRUE) {SetBit(kFriends, fr);}
  virtual void   SetMCdata(Bool_t mc = kTRUE) {SetBit(kMCdata, mc);}
  virtual void   SetPostProcess(Bool_t pp = kTRUE) {SetBit(kPostProcess, pp);}
  virtual void   Terminate(Option_t *) = 0;

private:
  AliTRDrecoTask(const AliTRDrecoTask&);
  AliTRDrecoTask& operator=(const AliTRDrecoTask&);

protected:
  UChar_t   fNRefFigures;  //! no of reference figures reported by task
  UChar_t   fDebugLevel;   //! Debug level 
  TObjArray *fContainer;   //! container to store results
  TObjArray *fTracks;      //! Array of tracks
  TTreeSRedirector *fDebugStream;  //! Debug stream 

  ClassDef(AliTRDrecoTask, 0) // base TRD reconstruction task
};

#endif

