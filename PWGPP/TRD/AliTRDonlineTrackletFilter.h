#ifndef ALITRDONLINETRACKLETFILTER
#define ALITRDONLINETRACKLETFILTER

#include "AliAnalysisTask.h"

class TList;

class AliInputEventHandler;
class AliVEvent;
class AliAODEvent;
class AliMCEvent;

class AliTRDgeometry;
class AliTRDpadPlane;
class AliTRDtrackletMCM;
class AliTRDtrackletWord;

class AliTRDonlineTrackletFilter : public AliAnalysisTask
{
 public:
  AliTRDonlineTrackletFilter(const char *name);
  ~AliTRDonlineTrackletFilter();

  void ConnectInputData(const Option_t *option);
  void CreateOutputObjects();
  void Exec(const Option_t *option);
  void LocalInit();
  void Terminate(const Option_t *option);

  Bool_t Notify(); 
  Bool_t LoadEvent();


 protected:
  AliESDEvent *fESD;                    //!

  AliInputEventHandler *fInputHandler;  //!
  AliVEvent            *fInputEvent;    //!
  AliAODEvent          *fOutputAOD;     //!
  AliMCEvent           *fMCEvent;       //!

  TClonesArray         *fTrackletsRaw;  //!
  TClonesArray         *fTrackletsSim;  //!

  // ----- output objects -----
  TTree                *fTrackletTree;  //!

  // ----- internal use -----
  AliTRDgeometry       *fGeo; //! TRD geometry

  Int_t fNevent;

  TString fPath; //!
  TFile *fTrackletFile; //!
  Int_t fNEventsPerFile; //!
  Int_t fEvent;  //!
  Int_t fFileNumber; //!
  TTree *fTrackletTreeSim;  //!
  TTree *fTrackletTreeRaw; //!

 private:
  AliTRDonlineTrackletFilter(const AliTRDonlineTrackletFilter&); // not implemented
  AliTRDonlineTrackletFilter& operator=(const AliTRDonlineTrackletFilter&); // not implemented

  ClassDef(AliTRDonlineTrackletFilter, 0);
};

#endif
