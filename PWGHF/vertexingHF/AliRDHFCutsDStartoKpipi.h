#ifndef ALIRDHFCUTSDSTARTOKPIPI_H
#define ALIRDHFCUTSDSTARTOKPIPI_H
/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//***********************************************************
/// \class Class AliRDHFCutsDStartoKpipi
/// \brief class for cuts on AOD reconstructed DStar->Kpipi
/// \author Author: A.Grelli, alessandro.grelli@uu.nl
/// PID method implemented by   Y.Wang, yifei@physi.uni-heidelberg.de
//***********************************************************

#include "AliRDHFCuts.h"

class AliAODEvent;
class AliAODRecoCascadeHF;
class AliAODRecoDecayHF;

class AliRDHFCutsDStartoKpipi : public AliRDHFCuts 
{
 public:

  AliRDHFCutsDStartoKpipi(const char* name="CutsDStartoKpipi");
  
  virtual ~AliRDHFCutsDStartoKpipi();

  AliRDHFCutsDStartoKpipi(const AliRDHFCutsDStartoKpipi& source);
  AliRDHFCutsDStartoKpipi& operator=(const AliRDHFCutsDStartoKpipi& source); 
 
  using AliRDHFCuts::GetCutVarsForOpt;
  virtual void GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgdaughters);

  using AliRDHFCuts::IsSelected;
  virtual Int_t IsSelected(TObject* obj, Int_t selectionLevel, AliAODEvent* aod);
  virtual Int_t IsSelected(TObject* obj, Int_t selectionLevel) {return IsSelected(obj,selectionLevel,0);}

  Int_t IsD0FromDStarSelected(Double_t pt, TObject* obj,Int_t selectionLevel, AliAODEvent* aod) const;
  Int_t IsD0FromDStarSelected(Double_t pt, TObject* obj,Int_t selectionLevel) {return IsD0FromDStarSelected(pt,obj,selectionLevel,0);};

  virtual Int_t PreSelect(TObjArray aodtracks);
  virtual Int_t IsSelectedPID(AliAODRecoDecayHF *rd);
  Int_t IsSelectedPID(Double_t pt, TObjArray aodTracks);
  virtual Int_t SelectPID(AliAODTrack *track, Int_t type);
  virtual Bool_t IsInFiducialAcceptance(Double_t pt,Double_t y) const;
  Float_t GetMassCut(Int_t iPtBin=0) const { return (GetCuts() ? fCutsRD[GetGlobalIndex(9,iPtBin)] : 1.e6);} // for the Dstar
  Float_t GetDCACut(Int_t iPtBin=0) const { return (GetCuts() ? fCutsRD[GetGlobalIndex(1,iPtBin)] : 1.e6);} // for the D0

  // standard cuts
  virtual void SetStandardCutsPP2010();
  virtual void SetStandardCutsPbPb2010();
  virtual void SetStandardCutsPbPb2011();

  // standard cuts 
  void SetStandardCutsPbPb2011DStar(TH1F *hfl);
  void SetStandardCutsPP2010DStarMult(Bool_t rec = kFALSE);

  
  void SetMaxPtPid(Float_t maxPt){fMaxPtPid = maxPt;}

  void SetOffHighPtPIDinTPC(Float_t TPCrem =999.){fTPCflag = TPCrem;}

  void AddTrackCutsSoftPi(const AliESDtrackCuts *cuts) 
     {fTrackCutsSoftPi=new AliESDtrackCuts(*cuts); return;}
  virtual AliESDtrackCuts *GetTrackCutsSoftPi() const {return fTrackCutsSoftPi;}

  Double_t GetCircRadius() { return fCircRadius; }
  void SetCircRadius(Double_t radius) { fCircRadius = radius; }

  void SetUseTPCtrackCutsOnD0Daughters(Bool_t flag=kTRUE) { fUseTPCtrackCutsOnD0Daughters = flag; }
  void SetUseTPCtrackCutsOnSoftPion(Bool_t flag=kFALSE) { fUseTPCtrackCutsOnSoftPion = flag; }
  Bool_t GetUseTPCtrackCutsOnD0Daughters() const { return fUseTPCtrackCutsOnD0Daughters; }
  Bool_t GetUseTPCtrackCutsOnSoftPion() const { return fUseTPCtrackCutsOnSoftPion; }

 protected:

  AliESDtrackCuts *fTrackCutsSoftPi; /// cuts for soft pion (AOD converted to ESD on the flight!)
  Float_t fMaxPtPid; /// maximum Dstar Pt for using PID
  Float_t fTPCflag;   ///
  Double_t fCircRadius; /// Radius for circular PID nsigma cut
  Bool_t fUseTPCtrackCutsOnD0Daughters; ///flag to apply TPC track quality cuts on D0 daughter from Dstar decay (used for different strategies for soft pion and D0daughters from Dstar decay) 
  Bool_t fUseTPCtrackCutsOnSoftPion; ///flag to apply TPC track quality cuts on soft pion from Dstar decay (used for different strategies for soft pion and D0daughters from Dstar decay)

  /// \cond CLASSIMP    
  ClassDef(AliRDHFCutsDStartoKpipi,8);  /// class for cuts on AOD reconstructed D0->Kpipi
  /// \endcond
};

#endif
