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
#include <TRandom3.h>

using std::vector;

class AliHFTreeHandlerLctopKpi : public AliHFTreeHandler
{
  public:
    AliHFTreeHandlerLctopKpi();
    AliHFTreeHandlerLctopKpi(int PIDopt);

    virtual ~AliHFTreeHandlerLctopKpi();

    virtual TTree* BuildTree(TString name="tree", TString title="tree");
    virtual bool SetVariables(int runnumber, unsigned int eventID, AliAODRecoDecayHF* cand, float bfield, int masshypo=0, AliPIDResponse *pidrespo=0x0);
    virtual void FillTree();

  private:

    vector<float> fImpParProng[knMaxProngs]; ///vectors of prong impact parameter
    vector<float> fSigmaVertex; /// vector of candidate sigma vertex
    vector<float> fDist12toPrim; /// vector of candidate distance between track 1-2 vertex to primary vertex
    vector<float> fDist23toPrim; /// vector of candidate distance between track 2-3 vertex to primary vertex
    vector<float> fNormd0MeasMinusExp; ///vector of candidate topomatic variable
    TRandom3 *fRandom;

    /// \cond CLASSIMP
    ClassDef(AliHFTreeHandlerLctopKpi,2); /// 
    /// \endcond
};

#endif
