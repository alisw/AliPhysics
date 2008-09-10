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
    kHasMCdata = BIT(14)
  };
  AliTRDrecoTask(const char *name, const char *title);
  virtual ~AliTRDrecoTask();

  void           ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects() = 0;
  virtual void   Exec(Option_t *) = 0;

  Int_t          GetDebugLevel() const { return fDebugLevel;}
  Int_t          GetNRefFigures() const { return fNRefFigures; } 
  virtual void   GetRefFigure(Int_t ifig, Int_t &first, Int_t &last);
  Bool_t         HasMCdata() const {return TestBit(kHasMCdata);};
  virtual Bool_t Load(Char_t *filename);
  virtual Bool_t PostProcess();
  virtual void   SetDebugLevel(Int_t level);
  void           SetMCdata(Bool_t mcdata) {SetBit(kHasMCdata, mcdata);}
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

