#ifndef ALIHFTREEHANDLERBSTODSPI_H
#define ALIHFTREEHANDLERBSTODSPI_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//*************************************************************************
// \class AliHFTreeHandlerBstoDspi
// \brief helper class to handle a tree for Bs cut optimisation and MVA analyses
// \authors:
// L. Vermunt, luuk.vermunt@cern.ch
// G. Innocenti, gian.michele.innocenti@cern.ch
/////////////////////////////////////////////////////////////

#include "AliHFTreeHandler.h"
#include "AliRDHFCutsDstoKKpi.h"
#include "AliAODRecoDecayHF2Prong.h"

class AliHFTreeHandlerBstoDspi : public AliHFTreeHandler
{
  public:
  
    //Flag to store injected candidate + HIJING track for background shape studies
    enum dsmesontype {
      kDsPrompt  = BIT(11),
      kDsFDBplus = BIT(12),
      kDsFDB0    = BIT(13),
      kDsFDLb0   = BIT(14),
      kDsFDBs0   = BIT(15)
    };

    //Need a bit for reflections/different decay? To check
    //static const int kDplustoKKpi = BIT(11);

    AliHFTreeHandlerBstoDspi();
    AliHFTreeHandlerBstoDspi(int PIDopt);

    virtual ~AliHFTreeHandlerBstoDspi();

    virtual TTree* BuildTree(TString name="tree", TString title="tree");
    virtual bool SetVariables(int runnumber, int eventID, int eventID_Ext, Long64_t eventID_Long, float ptgen, AliAODRecoDecayHF* cand, float bfield, int masshypo=0, AliPIDResponse *pidrespo=0x0);
    Int_t IsBsPionSelected(TObject* obj, AliRDHFCutsDstoKKpi* cutsDs, AliAODPidHF* fPidHFDs, AliAODEvent* aod, AliAODVertex *vtx);
    Int_t IsBsSelected(AliAODRecoDecayHF2Prong* bs);
    void SetDsBackgroundShapeType(bool isPr, bool isFDBplus, bool isFDB0, bool isFDLb0, bool isFDBs0);

    void SetBsSelectionValues(float invmass, float pt, float impparprod, float cosp, float cospxy){
      fInvMassBsCut = invmass;
      fPtBsCut = pt;
      fImpParProdBsCut = impparprod;
      fCosPBsCut = cosp;
      fCosPXYBsCut = cospxy;
    }

  private:

    //Variables for Bs->Dspi
    float fImpParProng[knMaxProngs]; ///prong impact parameter
    float fCosThetaStar; ///candidate cos theta star
    float fImpParProd; ///candidate product of impact parameter
    float fNormd0MeasMinusExp; ///candidate topomatic variable
  
    //Variables for Ds->KKpi
    float fInvMass_Ds; ///Ds invariant mass
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
  
    float fInvMassBsCut; ///Cut on invariant mass for Bs selection
    float fPtBsCut; ///Cut on pT for Bs selection
    float fImpParProdBsCut; ///Cut on d0xd0 for Bs selection
    float fCosPBsCut; ///Cut on cos pointing angle for Bs selection
    float fCosPXYBsCut; ///Cut on cos pointing angle xy for Bs selection

    /// \cond CLASSIMP
    ClassDef(AliHFTreeHandlerBstoDspi,2); ///
    /// \endcond
};
#endif
