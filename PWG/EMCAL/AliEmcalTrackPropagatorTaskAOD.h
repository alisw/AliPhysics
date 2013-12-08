#ifndef ALIEMCALTRACKPROPAGATORTASKAOD_H
#define ALIEMCALTRACKPROPAGATORTASKAOD_H

// $Id$

class TClonesArray;
class AliEMCALRecoUtils;
class AliAODEvent;
class AliAODtrack;

#include "AliAnalysisTaskSE.h"

class AliEmcalTrackPropagatorTaskAOD : public AliAnalysisTaskSE {
 public:
  AliEmcalTrackPropagatorTaskAOD();
  AliEmcalTrackPropagatorTaskAOD(const char *name);
  virtual ~AliEmcalTrackPropagatorTaskAOD();

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
  AliAODEvent       *fAodEv;             //!aod event
  TClonesArray      *fTracks;            //!track array

 private:
  AliEmcalTrackPropagatorTaskAOD(const AliEmcalTrackPropagatorTaskAOD&);            // not implemented
  AliEmcalTrackPropagatorTaskAOD &operator=(const AliEmcalTrackPropagatorTaskAOD&); // not implemented

  ClassDef(AliEmcalTrackPropagatorTaskAOD, 1); // Class to propagate and store track parameters at EMCAL surface
};

#endif
