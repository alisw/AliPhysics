#ifndef ALIHFTREEHANDLERBPLUSTOD0PI_H
#define ALIHFTREEHANDLERBPLUSTOD0PI_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliHFTreeHandlerBplustoD0pi
// \brief helper class to handle a tree for B+ cut optimisation and MVA analyses
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

class AliHFTreeHandlerBplustoD0pi : public AliHFTreeHandler
{
  public:
    AliHFTreeHandlerBplustoD0pi();
    AliHFTreeHandlerBplustoD0pi(int PIDopt);

    virtual ~AliHFTreeHandlerBplustoD0pi();

    virtual TTree* BuildTree(TString name="tree", TString title="tree");
    virtual bool SetVariables(int runnumber, unsigned int eventID, AliAODRecoDecayHF* cand, float bfield, int masshypo=0, AliPIDResponse *pidrespo=0x0);

  private:

    float fImpParProng[knMaxProngs]; ///vectors of prong impact parameter
    float fCosThetaStar; ///vector of candidate cos theta star
    float fImpParProd; ///vector of candidate product of impact parameter
    float fNormd0MeasMinusExp; ///vector of candidate topomatic variable
    float fAngleProngs; ///vector of angle between candidates prongs
    
    float fInvMass_D0; ///vector of candidate invariant mass D0
    float fPt_D0; ///vector of D0 pt
    //We save Y, Eta and Phi in common framework only for the Bplus and the 3 prongs, so missing D0 candidate. Added as Bplus specific variables. (Could be deleted if output size grows to big?)
    float fY_D0; ///vector of D0 rapidity
    float fEta_D0; ///vector of D0 pseudorapidity
    float fPhi_D0; ///vector of D0 azimuthal angle
    float fDecayLength_D0; ///vector of D0 decay length
    float fDecayLengthXY_D0; ///vector of D0 decay length in the transverse plane
    float fNormDecayLengthXY_D0; ///vector of D0 normalised decay length in the transverse plane
    float fCosP_D0; ///vector of D0 cosine of pointing angle
    float fCosPXY_D0; ///vector of D0 cosine of pointing angle in the transcverse plane
    float fImpParXY_D0; ///vector of D0 impact parameter in the transverse plane
    float fCosThetaStar_D0; ///vector of D0 cos theta star
    float fImpParProd_D0; ///vector of D0 product of impact parameter
    float fNormd0MeasMinusExp_D0; ///vector of D0 topomatic variable
    float fDCA_D0; ///vector of D0 DCA variable
    float fAngleProngs_D0; ///vector of angle between D0's prongs
    
    /// \cond CLASSIMP
    ClassDef(AliHFTreeHandlerBplustoD0pi,3); /// 
    /// \endcond
};
#endif
