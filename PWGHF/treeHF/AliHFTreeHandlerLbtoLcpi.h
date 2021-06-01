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
  
    //Flag to store injected candidate + HIJING track for background shape studies
    enum lcmesontype {
      kLcPrompt  = BIT(11), //includes all charmed hadron decays to Lc
      kLcFDBplus = BIT(12),
      kLcFDB0    = BIT(13),
      kLcFDLb0   = BIT(14),
      kLcFDBs0   = BIT(15)
    };

    //Need a bit for reflections/different decay? or an enum for resdecaytype? To check
    //enum resdecaytype {kNonResonant = 1, kL1520 = 2, kKstar = 3, kDelta = 4};
  
    AliHFTreeHandlerLbtoLcpi();
    AliHFTreeHandlerLbtoLcpi(int PIDopt);

    virtual ~AliHFTreeHandlerLbtoLcpi();

    virtual TTree* BuildTree(TString name="tree", TString title="tree");
    virtual bool SetVariables(int runnumber, int eventID, int eventID_Ext, Long64_t eventID_Long, float ptgen, AliAODRecoDecayHF* cand, float bfield, int masshypo=0, AliPIDResponse *pidrespo=0x0, AliAODPidHF *pidhf=0x0);
    Int_t IsLbPionSelected(TObject* obj, AliRDHFCutsLctopKpi* cutsLc, AliAODPidHF* fPidHFLc, AliAODEvent* aod, AliAODVertex *vtx);
    Int_t IsLbSelected(AliAODRecoDecayHF2Prong* lb);

    void SetLcBackgroundShapeType(bool isPr, bool isFDBplus, bool isFDB0, bool isFDLb0, bool isFDBs0);
    void SetLbSelectionValues(float invmass, float pt, float impparprod, float cosp, float cospxy){
      fInvMassLbCut = invmass;
      fPtLbCut = pt;
      fImpParProdLbCut = impparprod;
      fCosPLbCut = cosp;
      fCosPXYLbCut = cospxy;
    }

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
  
    float fInvMassLbCut; ///Cut on invariant mass for Lb selection
    float fPtLbCut; ///Cut on pT for Lb selection
    float fImpParProdLbCut; ///Cut on d0xd0 for Lb selection
    float fCosPLbCut; ///Cut on cos pointing angle for Lb selection
    float fCosPXYLbCut; ///Cut on cos pointing angle xy for Lb selection

    /// \cond CLASSIMP
    ClassDef(AliHFTreeHandlerLbtoLcpi,5); ///
    /// \endcond
};
#endif
