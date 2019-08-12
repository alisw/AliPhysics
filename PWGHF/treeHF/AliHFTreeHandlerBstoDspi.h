#ifndef ALIHFTREEHANDLERBSTODSPI_H
#define ALIHFTREEHANDLERBSTODSPI_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliHFTreeHandlerBstoDspi
// \brief helper class to handle a tree for B+ cut optimisation and MVA analyses
// \authors:
// L. Vermunt, luuk.vermunt@cern.ch
// G. Innocenti, gian.michele.innocenti@cern.ch
/////////////////////////////////////////////////////////////

#include "AliHFTreeHandler.h"
#include "AliRDHFCutsDstoKKpi.h"

class AliHFTreeHandlerBstoDspi : public AliHFTreeHandler
{
  public:
  
    //Need a bit for reflections/different decay? To check
    //static const int kDplustoKKpi = BIT(11);

    AliHFTreeHandlerBstoDspi();
    AliHFTreeHandlerBstoDspi(int PIDopt);

    virtual ~AliHFTreeHandlerBstoDspi();

    virtual TTree* BuildTree(TString name="tree", TString title="tree");
    virtual bool SetVariables(int runnumber, unsigned int eventID, float ptgen, AliAODRecoDecayHF* cand, float bfield, int masshypo=0, AliPIDResponse *pidrespo=0x0);
    Int_t IsBsPionSelected(TObject* obj, AliRDHFCutsDstoKKpi* cutsDs, AliAODPidHF* fPidHFDs, AliAODEvent* aod, AliAODVertex *vtx);

  private:

    float fImpParProng[knMaxProngs]; ///prong impact parameter
    float fCosThetaStar; ///candidate cos theta star
    float fImpParProd; ///candidate product of impact parameter
    float fNormd0MeasMinusExp; ///candidate topomatic variable
  
    float fInvMass_Ds; ///candidate invariant mass Ds
    float fPt_Ds; ///Ds pt
    float fY_Ds; ///Ds rapidity
    float fEta_Ds; ///Ds pseudorapidity
    float fPhi_Ds; ///Ds azimuthal angle
    float fDecayLength_Ds; ///Ds decay length
    float fDecayLengthXY_Ds; ///Ds decay length in the transverse plane
    float fNormDecayLengthXY_Ds; ///Ds normalised decay length in the transverse plane
    float fCosP_Ds; ///Ds cosine of pointing angle
    float fCosPXY_Ds; ///Ds cosine of pointing angle in the transcverse plane
    float fImpParXY_Ds; ///Ds impact parameter in the transverse plane
    float fDCA_Ds; ///Ds DCA variable
    float fSigmaVertex_Ds; /// Ds sigma vertex
    float fMassKK_Ds; /// Ds massKK
    float fCosPiDs_Ds; /// Ds cos3piDs
    float fCosPiKPhi_Ds; /// Ds cospiKphi
    float fNormd0MeasMinusExp_Ds; ///Ds topomatic variable
  
    /// \cond CLASSIMP
    ClassDef(AliHFTreeHandlerBstoDspi,1); ///
    /// \endcond
};
#endif
