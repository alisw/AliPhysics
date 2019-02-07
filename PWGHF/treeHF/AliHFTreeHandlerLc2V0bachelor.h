#ifndef ALIHFTREEHANDLERLC2V0BACHELOR_H
#define ALIHFTREEHANDLERLC2V0BACHELOR_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliHFTreeHandlerLc2V0bachelor
// \brief helper class to handle a tree for Lc cut optimisation and MVA analyses
// \authors:
// C. Zampolli
/////////////////////////////////////////////////////////////

#include "AliHFTreeHandler.h"

using std::vector;

class AliHFTreeHandlerLc2V0bachelor : public AliHFTreeHandler
{
  public:
    AliHFTreeHandlerLc2V0bachelor();
    AliHFTreeHandlerLc2V0bachelor(int PIDopt);

    virtual ~AliHFTreeHandlerLc2V0bachelor();

    virtual TTree* BuildTree(TString name = "tree", TString title = "tree");
    virtual bool SetVariables(AliAODRecoDecayHF* cand, float bfield, int masshypo = 0, AliPIDResponse* pidrespo = 0x0);
    virtual void FillTree();

  private:

    vector<float> fImpParProng[knMaxProngs]; ///vectors of prong impact parameter
    vector<float> fInvMassK0s; /// invariant mass of K0s
    vector<float> fcTauK0s; /// vector of cTau of the K0s
    vector<float> fV0PointingAngle; ///vector of K0s pointing angle
    vector<float> fCosThetaStar; /// cos theta star
    vector<float> fsignd0; // signed d0

    /// \cond CLASSIMP
    ClassDef(AliHFTreeHandlerLc2V0bachelor, 1); /// 
    /// \endcond
};

#endif
