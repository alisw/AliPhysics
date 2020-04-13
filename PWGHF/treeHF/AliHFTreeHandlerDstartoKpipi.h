#ifndef ALIHFTREEHANDLERDSTARTOKPIPI_H
#define ALIHFTREEHANDLERDSTARTOKPIPI_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliHFTreeHandlerDstartoKpipi
// \brief helper class to handle a tree for Dstar cut optimisation and MVA analyses
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

class AliHFTreeHandlerDstartoKpipi : public AliHFTreeHandler
{
  public:
    AliHFTreeHandlerDstartoKpipi();
    AliHFTreeHandlerDstartoKpipi(int PIDopt);

    virtual ~AliHFTreeHandlerDstartoKpipi();

    virtual TTree* BuildTree(TString name="tree", TString title="tree");
    virtual bool SetVariables(int runnumber, int eventID, int eventID_Ext, Long64_t eventID_Long, float ptgen, AliAODRecoDecayHF* cand, float bfield, int masshypo=0, AliPIDResponse *pidrespo=nullptr);

  private:

    float fImpParProng[knMaxProngs]; ///prong impact parameter
    float fCosThetaStar; ///candidate cos theta star
    float fImpParProd; ///D0 product of impact parameter
    float fNormd0MeasMinusExp; ///candidate topomatic variable
    float fAngleD0dkpPisoft; ///angle between D0 decay plane and soft pion

    float fInvMass_D0; ///candidate invariant mass D0
    float fPt_D0; ///D0 pt
    float fY_D0; ///D0 rapidity
    float fEta_D0; ///D0 pseudorapidity
    float fPhi_D0; ///D0 azimuthal angle

    /// \cond CLASSIMP
    ClassDef(AliHFTreeHandlerDstartoKpipi,4); ///
    /// \endcond
};
#endif
