#ifndef ALIHFTREEHANDLERAPPLYDSTOKKPI_H
#define ALIHFTREEHANDLERAPPLYDSTOKKPI_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliHFTreeHandlerApplyDstoKKpi
// \brief helper class to handle a tree for Ds cut optimisation and MVA analyses
// \authors:
// L. Vermunt, luuk.vermunt@cern.ch
/////////////////////////////////////////////////////////////

#include "AliHFTreeHandlerApply.h"

class AliHFTreeHandlerApplyDstoKKpi : public AliHFTreeHandlerApply
{
public:
  
  enum massKKopt {kMassKK,kDeltaMassKKPhi};
  
  static const int kDplustoKKpi = BIT(11);
  
  AliHFTreeHandlerApplyDstoKKpi();
  AliHFTreeHandlerApplyDstoKKpi(int PIDopt);
  virtual ~AliHFTreeHandlerApplyDstoKKpi();
  
  AliHFTreeHandlerApplyDstoKKpi(const AliHFTreeHandlerApplyDstoKKpi &source) = delete;
  AliHFTreeHandlerApplyDstoKKpi& operator=(const AliHFTreeHandlerApplyDstoKKpi &source) = delete;
  
  virtual TTree* BuildTree(TString name="tree", TString title="tree");
  virtual bool SetVariables(int runnumber, unsigned int eventID, float ptgen, float mlprob, AliAODRecoDecayHF* cand, float bfield, int masshypo=0, AliPIDResponse *pidrespo=nullptr);
  
  void SetMassKKOption(int opt) {fMassKKOpt=opt;}
  void SetIsDplustoKKpi(bool isDplus) {
    if(isDplus) fCandType |= kDplustoKKpi;
    else fCandType &= ~kDplustoKKpi;
  }
  
  static bool IsDplustoKKpi(int candtype) {
    if(candtype>>11&1) return true;
    return false;
  }
  
private:
  
  float fImpParProng[knMaxProngs];       /// prong impact parameter
  float fSigmaVertex;                    /// candidate sigma vertex
  float fMassKK;                         /// candidate massKK
  float fCosPiDs;                        /// candidate cos3piDs
  float fCosPiKPhi;                      /// candidate cospiKphi
  float fNormd0MeasMinusExp;             /// candidate topomatic variable
  int fMassKKOpt;                        /// option for massKK variable (mass or delta mass wrt phi)
  
  /// \cond CLASSIMP
  ClassDef(AliHFTreeHandlerApplyDstoKKpi,1); ///
  /// \endcond
};
#endif
