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
/////////////////////////////////////////////////////////////

#include "AliHFTreeHandler.h"

using std::vector;

class AliHFTreeHandlerDplustoKpipi : public AliHFTreeHandler
{
  public:
    AliHFTreeHandlerDplustoKpipi();
    AliHFTreeHandlerDplustoKpipi(int PIDopt);

    virtual ~AliHFTreeHandlerDplustoKpipi();

    virtual TTree* BuildTree(TString name="tree", TString title="tree");
    virtual bool SetVariables(AliAODRecoDecayHF* cand, float bfield, int masshypo=0, AliAODPidHF* pidHF=0x0);
    virtual void FillTree();

  private:

    vector<float> fImpParProng[knMaxProngs]; ///vectors of prong impact parameter
    vector<float> fSigmaVertex; /// vector of candidate sigma vertex
    vector<float> fNormd0MeasMinusExp; ///vector of candidate topomatic variable

    /// \cond CLASSIMP
    ClassDef(AliHFTreeHandlerDplustoKpipi,1); /// 
    /// \endcond
};

#endif
