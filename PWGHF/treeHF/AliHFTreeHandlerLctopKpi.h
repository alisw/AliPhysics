#ifndef ALIHFTREEHANDLERLCTOPKPI_H
#define ALIHFTREEHANDLERLCTOPKPI_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliHFTreeHandlerLctopKpi
// \brief helper class to handle a tree for Lc->pKpi cut optimisation and MVA analyses
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

class AliHFTreeHandlerLctopKpi : public AliHFTreeHandler
{
  public:
    AliHFTreeHandlerLctopKpi();
    AliHFTreeHandlerLctopKpi(int PIDopt);

    virtual ~AliHFTreeHandlerLctopKpi();

    virtual TTree* BuildTree(TString name="tree", TString title="tree");
    virtual bool SetVariables(int runnumber, unsigned int eventID, float ptgen, AliAODRecoDecayHF* cand, float bfield, int masshypo=0, AliPIDResponse *pidrespo=0x0);

  private:

    float fImpParProng[knMaxProngs]; ///prong impact parameter
    float fSigmaVertex; ///candidate sigma vertex
    float fDist12toPrim; ///candidate distance between track 1-2 vertex to primary vertex
    float fDist23toPrim; ///candidate distance between track 2-3 vertex to primary vertex
    float fNormd0MeasMinusExp; ///candidate topomatic variable

    /// \cond CLASSIMP
    ClassDef(AliHFTreeHandlerLctopKpi,3); /// 
    /// \endcond
};
#endif
