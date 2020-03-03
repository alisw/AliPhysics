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
    virtual bool SetVariables(int runnumber, int eventID, int eventID_Ext, Long64_t eventID_Long, float ptgen, AliAODRecoDecayHF* cand, float bfield, int masshypo = 0, AliPIDResponse* pidrespo = nullptr);

    void SetCalcSecoVtx(int opt) {fCalcSecoVtx=opt;}

    void SetIsLctoLpi(int isSeltoLpi, int isSelTopotoLpi, int isSelPIDtoLpi) {
      if(isSeltoLpi) fCandType |= kLctoLpi;
      else fCandType &= ~kLctoLpi;
      if(isSelTopotoLpi) fCandType |= kLcTopotoLpi;
      else fCandType &= ~kLcTopotoLpi;
      if(isSelPIDtoLpi) fCandType |= kLcPIDtoLpi;
      else fCandType &= ~kLcPIDtoLpi;
    }

  private:

    float fImpParProng[knMaxProngs]; ///prong impact parameter
    float fImpParK0s; /// impact parameter K0s
    float fDecayLengthK0s; /// decay length K0s
    float fInvMassK0s; /// invariant mass of K0s
    float fDCAK0s; ///DCA K0s prongs
    float fPtK0s; ///K0s pt
    float fEtaK0s; ///K0s pseudorapidity
    float fPhiK0s; ///K0s azimuthal angle
    float fcTauK0s; /// cTau of the K0s
    float fV0PointingAngle; ///K0s pointing angle
    float fCosThetaStar; ///cos theta star (proton - Lc)
    float fsignd0; //signed d0 proton (different from standard d0)
    float fArmqTOverAlpha; ///Armenteros qT/|alpha| of the K0s
    int fCalcSecoVtx; /// flag to calculate secondary vertex for Lc (if false, CommonDmesonVarBranches are not filled)

    /// \cond CLASSIMP
    ClassDef(AliHFTreeHandlerLc2V0bachelor, 3); ///
    /// \endcond
};
#endif
