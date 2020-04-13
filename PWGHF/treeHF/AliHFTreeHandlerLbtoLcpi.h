#ifndef ALIHFTREEHANDLERLBTOLCPI_H
#define ALIHFTREEHANDLERLBTOLCPI_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliHFTreeHandlerLbtoLcpi
// \brief helper class to handle a tree for Lb cut optimisation and MVA analyses
// \authors:
// D. Andreou, dimitra.andreou@cern.ch
// L. Vermunt, luuk.vermunt@cern.ch
// G. Innocenti, gian.michele.innocenti@cern.ch
/////////////////////////////////////////////////////////////

#include "AliHFTreeHandler.h"
#include "AliRDHFCutsLctopKpi.h"

class AliHFTreeHandlerLbtoLcpi : public AliHFTreeHandler
{
  public:
  
    //Need a bit for reflections/different decay? or an enum for resdecaytype? To check
    //enum resdecaytype {kNonResonant = 1, kL1520 = 2, kKstar = 3, kDelta = 4};
  
    AliHFTreeHandlerLbtoLcpi();
    AliHFTreeHandlerLbtoLcpi(int PIDopt);

    virtual ~AliHFTreeHandlerLbtoLcpi();

    virtual TTree* BuildTree(TString name="tree", TString title="tree");
    virtual bool SetVariables(int runnumber, int eventID, int eventID_Ext, Long64_t eventID_Long, float ptgen, AliAODRecoDecayHF* cand, float bfield, int masshypo=0, AliPIDResponse *pidrespo=0x0);
    Int_t IsLbPionSelected(TObject* obj, AliRDHFCutsLctopKpi* cutsLc, AliAODPidHF* fPidHFLc, AliAODEvent* aod, AliAODVertex *vtx);

  private:

    //Variables for Lb->Lcpi
    float fImpParProng[knMaxProngs]; ///prong impact parameter
    float fCosThetaStar; /// candidate costhetastar
    float fImpParProd; /// daughter impact-parameter product
    float fcTau; /// cTau of the Lb
    float fChi2OverNDF; /// chi2 over NDF of secondary vertex

    //Variables for Lc->pKpi
    float fInvMass_Lc; ///Lc invariant mass
    float fImpPar_Lc; /// impact parameter Lc
    float fPt_Lc; ///Lc pt
    float fY_Lc; ///Lc rapidity
    float fEta_Lc; ///Lc pseudorapidity
    float fPhi_Lc; ///Lc azimuthal angle
    float fDecayLength_Lc; ///Lc decay length
    float fDecayLengthXY_Lc; ///Lc decay length in the transverse plane
    float fNormDecayLengthXY_Lc; ///Lc normalised decay length in the transverse plane
    float fCosP_Lc; ///Lc cosine of pointing angle
    float fCosPXY_Lc; ///Lc cosine of pointing angle in the transcverse plane
    float fImpParXY_Lc; ///Lc impact parameter in the transverse plane
    float fDCA_Lc; ///Lc DCA variable
    float fDCAProng_Lc[knMaxProngs]; ///Lc prong DCA (pr0pr1, pr0pr2, pr1pr2)
    float fSigmaVertex_Lc; ///Lc sigma vertex
    float fDist12toPrim_Lc; ///Lc distance between track 1-2 vertex to primary vertex
    float fDist23toPrim_Lc; ///Lc distance between track 2-3 vertex to primary vertex
    float fNormd0MeasMinusExp_Lc; ///Lc topomatic variable
    float fSumImpParProngs_Lc; ///sum of Lc prong impact parameter squared
  
    /// \cond CLASSIMP
    ClassDef(AliHFTreeHandlerLbtoLcpi,4); ///
    /// \endcond
};
#endif
