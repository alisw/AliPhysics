#ifndef AliPHOSTracker_h
#define AliPHOSTracker_h

//-------------------------------------------------------------------------
//                          PHOS tracker.
// Matches ESD tracks with the PHOS and makes the PID.  
// Currently, has only one function implemented : PropagateBack(AliESD*).
//-------------------------------------------------------------------------

#include <AliTracker.h>
#include <AliLog.h>

class AliCluster;
class AliESD;
class TTree;
class AliRunLoader;

class AliPHOSTracker : public AliTracker
{
public:
  AliPHOSTracker():AliTracker()   {fRunLoader=0;}
  AliPHOSTracker(AliRunLoader *loader):AliTracker() {fRunLoader=loader;}
  virtual ~AliPHOSTracker()       {AliDebug(1,"Start.");}

  Int_t Clusters2Tracks(AliESD *) {AliDebug(1,"Start.");return 0;}
  Int_t RefitInward(AliESD *)     {AliDebug(1,"Start.");return 0;}
  void UnloadClusters()           {AliDebug(1,"Start.");}
  AliCluster *GetCluster(Int_t ) const {AliDebug(1,"Start.");return 0;}
  Int_t PropagateBack(AliESD *);
  Int_t LoadClusters(TTree *) {AliDebug(1,"Start.");return 0;}

  static void                SetDebug()   { fgDebug = kTRUE ; }
  static void                ResetDebug() { fgDebug = kFALSE ; }
  static Bool_t              Debug() { return fgDebug ; }

private:
  static Bool_t fgDebug ;    //! Verbosity controller
  AliRunLoader *fRunLoader;  //! Pointer to the run loader
  ClassDef(AliPHOSTracker,0)
};

#endif
