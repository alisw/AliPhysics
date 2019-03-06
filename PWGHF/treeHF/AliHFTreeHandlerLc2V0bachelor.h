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

    //Standard kSelected of AliHFTreeHandler is Lc->pK0s, but keep possibility to enable also Lc->Lpi (and charge conjugate together)
    enum isLctoLpi {
      kLctoLpi       = BIT(11),
      kLcTopotoLpi   = BIT(12),
      kLcPIDtoLpi    = BIT(13),
    };

    AliHFTreeHandlerLc2V0bachelor();
    AliHFTreeHandlerLc2V0bachelor(int PIDopt);

    virtual ~AliHFTreeHandlerLc2V0bachelor();

    virtual TTree* BuildTree(TString name = "tree", TString title = "tree");
    virtual bool SetVariables(int runnumber, unsigned int eventID, AliAODRecoDecayHF* cand, float bfield, int masshypo = 0, AliPIDResponse* pidrespo = 0x0);
    virtual void FillTree();

    void SetCalcSecoVtx(int opt) {fCalcSecoVtx=opt;}

    void SetIsLctoLpi(int isSeltoLpi, int isSelTopotoLpi, int isSelPIDtoLpi) {
      if(isSeltoLpi) fCandTypeMap |= kLctoLpi;
      else fCandTypeMap &= ~kLctoLpi;
      if(isSelTopotoLpi) fCandTypeMap |= kLcTopotoLpi;
      else fCandTypeMap &= ~kLcTopotoLpi;
      if(isSelPIDtoLpi) fCandTypeMap |= kLcPIDtoLpi;
      else fCandTypeMap &= ~kLcPIDtoLpi;
    }

  private:

    vector<float> fImpParProng[knMaxProngs]; ///vectors of prong impact parameter
    vector<float> fImpParK0s; /// vector of impact parameter K0s
    vector<float> fDecayLengthK0s; /// vector of decay length K0s
    vector<float> fInvMassK0s; /// invariant mass of K0s
    vector<float> fDCAK0s; ///vector of DCA K0s prongs
    vector<float> fPtK0s; ///vector of K0s pt
    vector<float> fEtaK0s; ///vector of K0s pseudorapidity
    vector<float> fPhiK0s; ///vector of K0s azimuthal angle
    vector<float> fcTauK0s; /// vector of cTau of the K0s
    vector<float> fV0PointingAngle; ///vector of K0s pointing angle
    vector<float> fCosThetaStar; ///vector of cos theta star (proton - Lc)
    vector<float> fsignd0; //vector of signed d0 proton (different from standard d0)
    vector<float> fArmqTOverAlpha; ///vector of Armenteros qT/|alpha| of the K0s
    int fCalcSecoVtx; /// flag to calculate secondary vertex for Lc (if false, CommonDmesonVarBranches are not filled)

    /// \cond CLASSIMP
    ClassDef(AliHFTreeHandlerLc2V0bachelor, 2); ///
    /// \endcond
};

#endif
