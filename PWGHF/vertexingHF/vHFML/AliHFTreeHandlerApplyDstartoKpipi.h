#ifndef ALIHFTREEHANDLERAPPLYDSTARTOKPIPI_H
#define ALIHFTREEHANDLERAPPLYDSTARTOKPIPI_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliHFTreeHandlerApplyDstartoKpipi
// \brief helper class to handle a tree for Dstar cut optimisation and MVA analyses
// \authors:
// L. Vermunt, luuk.vermunt@cern.ch
/////////////////////////////////////////////////////////////

#include "AliHFTreeHandlerApply.h"

class AliHFTreeHandlerApplyDstartoKpipi : public AliHFTreeHandlerApply
{
public:

  AliHFTreeHandlerApplyDstartoKpipi();
  AliHFTreeHandlerApplyDstartoKpipi(int PIDopt);
  virtual ~AliHFTreeHandlerApplyDstartoKpipi();

  AliHFTreeHandlerApplyDstartoKpipi(const AliHFTreeHandlerApplyDstartoKpipi &source) = delete;
  AliHFTreeHandlerApplyDstartoKpipi& operator=(const AliHFTreeHandlerApplyDstartoKpipi &source) = delete;

  virtual TTree* BuildTree(TString name="tree", TString title="tree");
  virtual bool SetVariables(int runnumber, int eventID, int eventID_Ext, Long64_t eventID_Long, float ptgen, float mlprob, AliAODRecoDecayHF* cand, float bfield, int masshypo=0, AliPIDResponse *pidrespo=nullptr, AliAODPidHF *pidhf=nullptr);
  void SetDstarMCMatchedLabel(int label) { fMCLabel = label; }

private:

  int fIDProng[knMaxProngs];         /// prong track ID
  int fChargeProng[knMaxProngs];     /// prong track charge
  float fImpParProng[knMaxProngs];   /// prong impact parameter
  int fCharge;                       /// Dstar charge
  float fCosThetaStar;               /// candidate cos theta star
  float fImpParProd;                 /// D0 product of impact parameter
  float fNormd0MeasMinusExp;         /// candidate topomatic variable
  float fAngleD0dkpPisoft;           /// angle between D0 decay plane and soft pion
  float fDeltaInvMassD0;             /// delta mass of D0

  /// \cond CLASSIMP
  ClassDef(AliHFTreeHandlerApplyDstartoKpipi,1); ///
  /// \endcond
};
#endif
