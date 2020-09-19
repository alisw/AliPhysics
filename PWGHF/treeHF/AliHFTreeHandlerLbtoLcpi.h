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
#include "AliAODRecoDecayHF2Prong.h"

class AliHFTreeHandlerLbtoLcpi : public AliHFTreeHandler
{
  public:
  
    //Need a bit for reflections/different decay? or an enum for resdecaytype? To check
    //enum resdecaytype {kNonResonant = 1, kL1520 = 2, kKstar = 3, kDelta = 4};
  
    AliHFTreeHandlerLbtoLcpi();
    AliHFTreeHandlerLbtoLcpi(int PIDopt);

    virtual ~AliHFTreeHandlerLbtoLcpi();

    virtual TTree* BuildTree(TString name="tree", TString title="tree");
    virtual bool SetVariables(int runnumber, int eventID, int eventID_Ext, Long64_t eventID_Long, Float_t ptgen, AliAODRecoDecayHF* cand, Float_t bfield, int masshypo=0, AliPIDResponse *pidrespo=0x0);
    Int_t IsLbPionSelected(TObject* obj, AliRDHFCutsLctopKpi* cutsLc, AliAODPidHF* fPidHFLc, AliAODEvent* aod, AliAODVertex *vtx);
    Int_t IsLbSelected(AliAODRecoDecayHF2Prong* lb);

  private:

    //Variables for Lb->Lcpi
    Float_t fImpParProng[knMaxProngs]; ///prong impact parameter
    Float_t fCosThetaStar; /// candidate costhetastar
    Float_t fImpParProd; /// daughter impact-parameter product
    Float_t fcTau; /// cTau of the Lb
    Float_t fChi2OverNDF; /// chi2 over NDF of secondary vertex

    //Variables for Lc->pKpi
    Float_t fInvMass_Lc; ///Lc invariant mass
    Float_t fImpPar_Lc; /// impact parameter Lc
    Float_t fPt_Lc; ///Lc pt
    Float_t fY_Lc; ///Lc rapidity
    Float_t fEta_Lc; ///Lc pseudorapidity
    Float_t fPhi_Lc; ///Lc azimuthal angle
    Float_t fDecayLength_Lc; ///Lc decay length
    Float_t fDecayLengthXY_Lc; ///Lc decay length in the transverse plane
    Float_t fNormDecayLengthXY_Lc; ///Lc normalised decay length in the transverse plane
    Float_t fCosP_Lc; ///Lc cosine of pointing angle
    Float_t fCosPXY_Lc; ///Lc cosine of pointing angle in the transcverse plane
    Float_t fImpParXY_Lc; ///Lc impact parameter in the transverse plane
    Float_t fDCA_Lc; ///Lc DCA variable
    Float_t fDCAProng_Lc[knMaxProngs]; ///Lc prong DCA (pr0pr1, pr0pr2, pr1pr2)
    Float_t fSigmaVertex_Lc; ///Lc sigma vertex
    Float_t fDist12toPrim_Lc; ///Lc distance between track 1-2 vertex to primary vertex
    Float_t fDist23toPrim_Lc; ///Lc distance between track 2-3 vertex to primary vertex
    Float_t fNormd0MeasMinusExp_Lc; ///Lc topomatic variable
    Float_t fSumImpParProngs_Lc; ///sum of Lc prong impact parameter squared
    

    Float_t fInvMassLbCut; ///Cut on invariant mass for Lb selection
    Float_t fPtLbCut; ///Cut on pT for Lb selection
    Float_t fImpParProdLbCut; ///Cut on d0xd0 for Lb selection
    Float_t fCosPLbCut; ///Cut on cos pointing angle for Lb selection
    Float_t fCosPXYLbCut; ///Cut on cos pointing angle xy for Lb selection

    /// \cond CLASSIMP
    ClassDef(AliHFTreeHandlerLbtoLcpi,4); ///
    /// \endcond
};
#endif
