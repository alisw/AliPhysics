// Copyright CERN. This software is distributed under the terms of the GNU
// General Public License v3 (GPL Version 3).
//
// See http://www.gnu.org/licenses/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

//**************************************************************************************
// \class AliGammaMLResponseGammatoEPlusEMinus
/////////////////////////////////////////////////////////////////////////////////////////

#include <TDatabasePDG.h>

#include "AliGammaMLResponseGammatoEPlusEMinus.h"


/// \cond CLASSIMP
ClassImp(AliGammaMLResponseGammatoEPlusEMinus);
/// \endcond

//________________________________________________________________
AliGammaMLResponseGammatoEPlusEMinus::AliGammaMLResponseGammatoEPlusEMinus() : AliGammaMLResponse()
{
    //
    // Default constructor
    //
}

//________________________________________________________________
AliGammaMLResponseGammatoEPlusEMinus::AliGammaMLResponseGammatoEPlusEMinus(const Char_t *name, const Char_t *title, 
                                               const std::string configfilepath) : AliGammaMLResponse(name, title, configfilepath)
{
    //
    // Standard constructor
}

//________________________________________________________________
AliGammaMLResponseGammatoEPlusEMinus::~AliGammaMLResponseGammatoEPlusEMinus()
{
    //
    // Destructor
    //
}

//--------------------------------------------------------------------------
AliGammaMLResponseGammatoEPlusEMinus::AliGammaMLResponseGammatoEPlusEMinus(const AliGammaMLResponseGammatoEPlusEMinus &source) : AliGammaMLResponse(source)
{
    //
    // Copy constructor
    //
}

AliGammaMLResponseGammatoEPlusEMinus &AliGammaMLResponseGammatoEPlusEMinus::operator=(const AliGammaMLResponseGammatoEPlusEMinus &source)
{
    //
    // assignment operator
    //
    if (&source == this)
        return *this;

    AliGammaMLResponse::operator=(source);

    return *this;
}

//________________________________________________________________
void AliGammaMLResponseGammatoEPlusEMinus::SetMapOfVariables(AliAODConversionPhoton *cand, AliVEvent* fInputEvent, AliConversionPhotonCuts* fiPhotonCut, AliV0ReaderV1* fV0Reader)
{
   
  AliPIDResponse* pidResonse = ((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetPIDResponse();

  AliVTrack * ENegTrack = fiPhotonCut->GetTrack(fInputEvent, cand->GetTrackLabelNegative());
  AliVTrack * EPosTrack = fiPhotonCut->GetTrack(fInputEvent, cand->GetTrackLabelPositive());
  
  if(ENegTrack && EPosTrack ) {
      fVars["fInvMass_PhotonML"]        = cand->M();
      fVars["fPt_PhotonML"]             = cand->Pt();
      fVars["fAlpha_PhotonML"]          = cand->GetArmenterosAlpha();
      fVars["fQt_PhotonML"]             = cand->GetArmenterosQt();
      fVars["fPsiPair_PhotonML"]        = cand->GetPsiPair();
      fVars["fEta_PhotonML"]            = cand->Eta();
      fVars["fPhi_PhotonML"]            = cand->GetPhotonPhi();
      fVars["fR_PhotonML"]              = cand->GetConversionRadius();
      fVars["fCosThetaPA_PhotonML"]     = fiPhotonCut->GetCosineOfPointingAngle(cand,fInputEvent);
      fVars["fChi2PerNDF_PhotonML"]     = cand->GetChi2perNDF();
      fVars["fPhotonQuality_PhotonML"]  = cand->GetPhotonQuality();
      fVars["fDCAr_PhotonML"]           = cand->GetDCArToPrimVtx();

      if (fiPhotonCut->CorrectedTPCClusterCut(cand, fInputEvent)){ // This condition is used to maintain uniformity with manual cut
        fVars["fNegClusterTPCToF_PhotonML"] = ENegTrack->GetTPCClusterInfo(2,0,fiPhotonCut->GetFirstTPCRow(cand->GetConversionRadius()));
        fVars["fPosClusterTPCToF_PhotonML"] = EPosTrack->GetTPCClusterInfo(2,0,fiPhotonCut->GetFirstTPCRow(cand->GetConversionRadius()));
      }

      fVars["fDeDx_ITS_ENeg_PhotonML"]   = pidResonse->NumberOfSigmasITS(ENegTrack,AliPID::kElectron);
      fVars["fTOF_ENeg_PhotonML"]        = pidResonse->NumberOfSigmasTOF(ENegTrack,AliPID::kElectron);
      fVars["fDeDx_TPC_ENeg_PhotonML"]   = fiPhotonCut->GetCorrectedElectronTPCResponse(ENegTrack->Charge(),pidResonse->NumberOfSigmasTPC(ENegTrack,AliPID::kElectron), ENegTrack->P(), ENegTrack->Eta(), ENegTrack->GetTPCNcls(),cand->GetConversionRadius() ); 
      fVars["fDeDx_PiM_ITS_PhotonML"]    = pidResonse->NumberOfSigmasITS(ENegTrack,AliPID::kPion);
      fVars["fDeDx_PiM_TPC_PhotonML"]    = pidResonse->NumberOfSigmasTPC(ENegTrack,AliPID::kPion);
      fVars["fTOF_PiM_PhotonML"]         = pidResonse->NumberOfSigmasTOF(ENegTrack,AliPID::kPion);

      fVars["fPt_ENeg_PhotonML"]         = ENegTrack->Pt();
      fVars["fEta_ENeg_PhotonML"]        = ENegTrack->Eta();
           
      fVars["fDeDx_ITS_EPos_PhotonML"]   = pidResonse->NumberOfSigmasITS(EPosTrack,AliPID::kElectron);
      fVars["fTOF_EPos_PhotonML"]        = pidResonse->NumberOfSigmasTOF(EPosTrack,AliPID::kElectron);
      fVars["fDeDx_TPC_EPos_PhotonML"]   =  fiPhotonCut->GetCorrectedElectronTPCResponse(EPosTrack->Charge(),pidResonse->NumberOfSigmasTPC(ENegTrack,AliPID::kElectron), EPosTrack->P(), EPosTrack->Eta(), EPosTrack->GetTPCNcls(),cand->GetConversionRadius() ); 

      fVars["fDeDx_PiP_ITS_PhotonML"]    =  pidResonse->NumberOfSigmasITS(EPosTrack,AliPID::kPion);
      fVars["fDeDx_PiP_TPC_PhotonML"]    =  pidResonse->NumberOfSigmasTPC(EPosTrack,AliPID::kPion);
      fVars["fTOF_PiP_PhotonML"]         = pidResonse->NumberOfSigmasTOF(EPosTrack,AliPID::kPion);
      
      fVars["fPt_EPos_PhotonML"]         = EPosTrack->Pt();
      fVars["fEta_EPos_PhotonML"]        = EPosTrack->Eta();
      fVars["fPAsymmetry_EPos_PhotonML"] = EPosTrack->P()/cand->GetPhotonP();
      fVars["fKind"]                     = 9;
    }

}
