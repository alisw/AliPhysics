#ifndef ALIEMCALTRACKPROPAGATORTASK_H
#define ALIEMCALTRACKPROPAGATORTASK_H

// $Id$

class TClonesArray;
class AliEMCALRecoUtils;
class AliESDEvent;
class AliESDtrack;

#include "AliAnalysisTaskSE.h"

class AliEmcalTrackPropagatorTask : public AliAnalysisTaskSE {
 public:
  AliEmcalTrackPropagatorTask();
  AliEmcalTrackPropagatorTask(const char *name);
  virtual ~AliEmcalTrackPropagatorTask();

  void UserCreateOutputObjects();
  void UserExec(Option_t *option);
   
  void SetDist(Double_t d)                 { fDist       = d;    }
  void SetMinPt(Double_t pt)               { fMinPtCut = pt;     }
  void SetRecoUtils(AliEMCALRecoUtils *ru) { fRecoUtils  = ru;   }
  void SetTracksName(const char *name)     { fTracksName = name; }

 protected:
  AliEMCALRecoUtils *fRecoUtils;         // esd reco utils
  TString            fTracksName;        // name of tracks 
  Double_t           fDist;              // distance to surface (430cm default)
  Double_t           fMinPtCut;          // minimum track pt cut (500 MeV/c default)
  AliESDEvent       *fEsdEv;             //!esd event
  TClonesArray      *fTracks;            //!track array

 private:
  AliEmcalTrackPropagatorTask(const AliEmcalTrackPropagatorTask&);            // not implemented
  AliEmcalTrackPropagatorTask &operator=(const AliEmcalTrackPropagatorTask&); // not implemented

  ClassDef(AliEmcalTrackPropagatorTask, 1); // Class to propagate and store track parameters at EMCAL surface
};

#endif
