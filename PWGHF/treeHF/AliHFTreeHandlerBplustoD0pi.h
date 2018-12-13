#ifndef ALIHFTREEHANDLERBPLUSTOD0PI_H
#define ALIHFTREEHANDLERBPLUSTOD0PI_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliHFTreeHandlerBplustoD0pi
// \brief helper class to handle a tree for D+ cut optimisation and MVA analyses
// \authors:
// F. Catalano, fabio.catalano@cern.ch
// A. Festanti, andrea.festanti@cern.ch
// F. Grosa, fabrizio.grosa@cern.ch
// G. Innocenti, gian.michele.innocenti@cern.ch
// F. Prino, prino@to.infn.it
// L. Vermunt, luuk.vermunt@cern.ch
// L. van Doremalen, lennart.van.doremalen@cern.ch
/////////////////////////////////////////////////////////////

#include "AliHFTreeHandler.h"

using std::vector;

class AliHFTreeHandlerBplustoD0pi : public AliHFTreeHandler
{
  public:
    AliHFTreeHandlerBplustoD0pi();
    AliHFTreeHandlerBplustoD0pi(int PIDopt);

    virtual ~AliHFTreeHandlerBplustoD0pi();

    virtual TTree* BuildTree(TString name="tree", TString title="tree");
    virtual bool SetVariables(AliAODRecoDecayHF* cand, float bfield, int masshypo=0, AliAODPidHF* pidHF=0x0);
    virtual void FillTree(); //to be called for each event, not each candidate!

  private:

    vector<float> fImpParProng[knMaxProngs]; ///vectors of prong impact parameter
    vector<float> fCosThetaStar; ///vector of candidate cos theta star
    vector<float> fImpParProd; ///vector of candidate product of impact parameter
    vector<float> fNormd0MeasMinusExp; ///vector of candidate topomatic variable
    vector<float> fDCA; ///vector of candidate DCA variable
    vector<float> fAngleProngs; ///vector of angle between candidates prongs
    
    vector<float> fInvMass_D0; ///vector of candidate invariant mass D0
    vector<float> fPt_D0; ///vector of D0 pt
    //We save Y, Eta and Phi in common framework only for the Bplus and the 3 prongs, so missing D0 candidate. Added as Bplus specific variables. (Could be deleted if output size grows to big?)
    vector<float> fY_D0; ///vector of D0 rapidity
    vector<float> fEta_D0; ///vector of D0 pseudorapidity
    vector<float> fPhi_D0; ///vector of D0 azimuthal angle
    vector<float> fDecayLength_D0; ///vector of D0 decay length
    vector<float> fDecayLengthXY_D0; ///vector of D0 decay length in the transverse plane
    vector<float> fNormDecayLengthXY_D0; ///vector of D0 normalised decay length in the transverse plane
    vector<float> fCosP_D0; ///vector of D0 cosine of pointing angle
    vector<float> fCosPXY_D0; ///vector of D0 cosine of pointing angle in the transcverse plane
    vector<float> fImpParXY_D0; ///vector of D0 impact parameter in the transverse plane
    vector<float> fCosThetaStar_D0; ///vector of D0 cos theta star
    vector<float> fImpParProd_D0; ///vector of D0 product of impact parameter
    vector<float> fNormd0MeasMinusExp_D0; ///vector of D0 topomatic variable
    vector<float> fDCA_D0; ///vector of D0 DCA variable
    vector<float> fAngleProngs_D0; ///vector of angle between D0's prongs
    
    /// \cond CLASSIMP
    ClassDef(AliHFTreeHandlerBplustoD0pi,1); /// 
    /// \endcond
};
#endif