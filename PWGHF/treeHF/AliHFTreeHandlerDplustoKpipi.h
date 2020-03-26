#ifndef ALIHFTREEHANDLERDPLUSTOKPIPI_H
#define ALIHFTREEHANDLERDPLUSTOKPIPI_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliHFTreeHandlerDplustoKpipi
// \brief helper class to handle a tree for D+ cut optimisation and MVA analyses
// \authors:
// F. Catalano, fabio.catalano@cern.ch
// A. Festanti, andrea.festanti@cern.ch
// F. Grosa, fabrizio.grosa@cern.ch
// G. Innocenti, gian.michele.innocenti@cern.ch
// F. Prino, prino@to.infn.it
// L. Vermunt, luuk.vermunt@cern.ch
// L. van Doremalen, lennart.van.doremalen@cern.ch
// J. Norman, jaime.norman@cern.ch
// G. Luparello, grazia.luparello@cern.ch
/////////////////////////////////////////////////////////////

#include "AliHFTreeHandler.h"

class AliHFTreeHandlerDplustoKpipi : public AliHFTreeHandler
{
  public:
    AliHFTreeHandlerDplustoKpipi();
    AliHFTreeHandlerDplustoKpipi(int PIDopt);

    virtual ~AliHFTreeHandlerDplustoKpipi();

    virtual TTree* BuildTree(TString name="tree", TString title="tree");
    virtual bool SetVariables(int runnumber, int eventID, int eventID_Ext, Long64_t eventID_Long, float ptgen, AliAODRecoDecayHF* cand, float bfield, int masshypo=0, AliPIDResponse *pidrespo=nullptr);

  private:

    float fImpParProng[knMaxProngs]; ///prong impact parameter
    float fSigmaVertex; /// candidate sigma vertex
    float fNormd0MeasMinusExp; ///candidate topomatic variable

    /// \cond CLASSIMP
    ClassDef(AliHFTreeHandlerDplustoKpipi,2); /// 
    /// \endcond
};
#endif
