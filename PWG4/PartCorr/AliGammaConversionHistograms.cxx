/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Ana Marin, Kathrin Koch, Kenneth Aamodt                        *
 * Version 1.0                                                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
/**
 * Class containing histograms 
 //Change here
 // here we need a description of the naming scheme of the histograms.

*/

#include "AliGammaConversionHistograms.h"
#include "TMath.h"

using namespace std;

ClassImp(AliGammaConversionHistograms)


AliGammaConversionHistograms::AliGammaConversionHistograms() :
  fOutputContainer(NULL),
  fNPhiIndex(0),
  fNRIndex(0),
  fMinRadius(0.),
  fMaxRadius(0.),
  fDeltaR(0.),
  fMinPhi(0.),
  fMaxPhi(0.),
  fDeltaPhi(0.),
  fMC_EP_R(NULL),
  fMC_EP_Z_R(NULL),
  fMC_EP_X_Y(NULL),
  fMC_EP_OpeningAngle(NULL),
  fMC_E_Energy(NULL),
  fMC_E_Pt(NULL),
  fMC_E_Eta(NULL),
  fMC_E_Phi(NULL),
  fMC_P_Energy(NULL),
  fMC_P_Pt(NULL),
  fMC_P_Eta(NULL),
  fMC_P_Phi(NULL),
  fMC_Gamma_Energy(NULL),
  fMC_Gamma_Pt(NULL),
  fMC_Gamma_Eta(NULL),
  fMC_Gamma_Phi(NULL),
  fMC_DirectGamma_Energy(NULL),
  fMC_DirectGamma_Pt(NULL),
  fMC_DirectGamma_Eta(NULL),
  fMC_DirectGamma_Phi(NULL),
  fMC_Mapping(),
  fMC_Mapping_Phi(),
  fMC_Mapping_R(),
  fMC_Match_Gamma_Eta(NULL),
  fMC_Match_Gamma_Phi(NULL),
  fMC_Match_Gamma_Pt(NULL),
  fMC_Match_Gamma_Energy(NULL),
  fMC_Match_Gamma_Mass(NULL),
  fMC_Match_Gamma_OpeningAngle(NULL),
  fMC_Match_Gamma_R(NULL),
  fMC_Match_Gamma_Z_R(NULL),
  fMC_Match_Gamma_X_Y(NULL),
  fMC_Pi0_Eta(NULL),
  fMC_Pi0_Phi(NULL),
  fMC_Pi0_Pt(NULL),
  fMC_Pi0_Energy(NULL),
  fMC_Pi0_Mass(NULL),
  fMC_Pi0_OpeningAngleGamma(NULL),
  fMC_Pi0_R(NULL),
  fMC_Pi0_Z_R(NULL),
  fMC_Pi0_X_Y(NULL),
  fMC_Pi0Secondaries_Eta(NULL),
  fMC_Pi0Secondaries_Phi(NULL),
  fMC_Pi0Secondaries_Pt(NULL),
  fMC_Pi0Secondaries_Energy(NULL),
  fMC_Pi0Secondaries_Mass(NULL),
  fMC_Pi0Secondaries_OpeningAngleGamma(NULL),
  fMC_Pi0Secondaries_R(NULL),
  fMC_Pi0Secondaries_Z_R(NULL),
  fMC_Pi0Secondaries_X_Y(NULL),
  fMC_Eta_Eta(NULL),
  fMC_Eta_Phi(NULL),
  fMC_Eta_Pt(NULL),
  fMC_Eta_Energy(NULL),
  fMC_Eta_Mass(NULL),
  fMC_Eta_OpeningAngleGamma(NULL),
  fMC_Eta_R(NULL),
  fMC_Eta_Z_R(NULL),
  fMC_Eta_X_Y(NULL),
  fESD_EP_R(NULL),
  fESD_EP_Z_R(NULL),
  fESD_EP_X_Y(NULL),
  fESD_EP_OpeningAngle(NULL),
  fESD_E_Energy(NULL),
  fESD_E_Pt(NULL),
  fESD_E_Eta(NULL),
  fESD_E_Phi(NULL),
  fESD_P_Energy(NULL),
  fESD_P_Pt(NULL),
  fESD_P_Eta(NULL),
  fESD_P_Phi(NULL),
  fESD_Gamma_Energy(NULL),
  fESD_Gamma_Pt(NULL),
  fESD_Gamma_Eta(NULL),
  fESD_Gamma_Phi(NULL),
  fESD_Mapping(),
  fESD_Mapping_Phi(),
  fESD_Mapping_R(),
  fESD_Match_Gamma_OpeningAngle(NULL),
  fESD_Match_Gamma_Energy(NULL),
  fESD_Match_Gamma_Pt(NULL),
  fESD_Match_Gamma_Eta(NULL),
  fESD_Match_Gamma_Phi(NULL),
  fESD_Match_Gamma_Mass(NULL),
  fESD_Match_Gamma_Width(NULL),
  fESD_Match_Gamma_Chi2(NULL),
  fESD_Match_Gamma_NDF(NULL),
  fESD_Match_Gamma_R(NULL),
  fESD_Match_Gamma_Z_R(NULL),
  fESD_Match_Gamma_X_Y(NULL),
  fESD_Pi0_OpeningAngleGamma(NULL),
  fESD_Pi0_Energy(NULL),
  fESD_Pi0_Pt(NULL),
  fESD_Pi0_Eta(NULL),
  fESD_Pi0_Phi(NULL),
  fESD_Pi0_Mass(NULL),
  fESD_Pi0_R(NULL),
  fESD_Pi0_Z_R(NULL),
  fESD_Pi0_X_Y(NULL),
  fESD_Eta_OpeningAngleGamma(NULL),
  fESD_Eta_Energy(NULL),
  fESD_Eta_Pt(NULL),
  fESD_Eta_Eta(NULL),
  fESD_Eta_Phi(NULL),
  fESD_Eta_Mass(NULL),
  fESD_Eta_R(NULL),
  fESD_Eta_Z_R(NULL),
  fESD_Eta_X_Y(NULL),
  fESD_Background_OpeningAngleGamma(NULL),
  fESD_Background_Energy(NULL),
  fESD_Background_Pt(NULL),
  fESD_Background_Eta(NULL),
  fESD_Background_Phi(NULL),
  fESD_Background_Mass(NULL),
  fESD_Background_R(NULL),
  fESD_Background_Z_R(NULL),
  fESD_Background_X_Y(NULL),
  fResolution_dPt(NULL),
  fResolution_dR(NULL),
  fResolution_dZ(NULL),
  fResolution_dR_dPt(NULL),
  fResolution_MC_Pt(NULL),
  fResolution_MC_R(NULL),
  fResolution_MC_Z(NULL),
  fResolution_ESD_Pt(NULL),
  fResolution_ESD_R(NULL),
  fResolution_ESD_Z(NULL),
  fNumberOfV0s(NULL),
  fNumberOfSurvivingV0s(NULL),
  fV0MassDebugCut1(NULL),
  fV0MassDebugCut2(NULL),
  fV0MassDebugCut3(NULL),
  fV0MassDebugCut4(NULL),
  fV0MassDebugCut5(NULL),
  fV0MassDebugCut6(NULL),
  fV0MassDebugCut7(NULL),
  fV0MassDebugCut8(NULL)
{

}


AliGammaConversionHistograms::AliGammaConversionHistograms(const AliGammaConversionHistograms & original) :
  fOutputContainer(original.fOutputContainer),
  fNPhiIndex(0),
  fNRIndex(0),
  fMinRadius(0.),
  fMaxRadius(0.),
  fDeltaR(0.),
  fMinPhi(0.),
  fMaxPhi(0.),
  fDeltaPhi(0.),
  fMC_EP_R(original.fMC_EP_R),
  fMC_EP_Z_R(original.fMC_EP_Z_R),
  fMC_EP_X_Y(original.fMC_EP_X_Y),
  fMC_EP_OpeningAngle(original.fMC_EP_OpeningAngle),
  fMC_E_Energy(original.fMC_E_Energy),
  fMC_E_Pt(original.fMC_E_Pt),
  fMC_E_Eta(original.fMC_E_Eta),
  fMC_E_Phi(original.fMC_E_Phi),
  fMC_P_Energy(original.fMC_P_Energy),
  fMC_P_Pt(original.fMC_P_Pt),
  fMC_P_Eta(original.fMC_P_Eta),
  fMC_P_Phi(original.fMC_P_Phi),
  fMC_Gamma_Energy(original.fMC_Gamma_Energy),
  fMC_Gamma_Pt(original.fMC_Gamma_Pt),
  fMC_Gamma_Eta(original.fMC_Gamma_Eta),
  fMC_Gamma_Phi(original.fMC_Gamma_Phi),
  fMC_DirectGamma_Energy(original.fMC_Gamma_Energy),
  fMC_DirectGamma_Pt(original.fMC_Gamma_Pt),
  fMC_DirectGamma_Eta(original.fMC_Gamma_Eta),
  fMC_DirectGamma_Phi(original.fMC_Gamma_Phi),
  fMC_Mapping(),
  fMC_Mapping_Phi(),
  fMC_Mapping_R(),
  fMC_Match_Gamma_Eta(original.fMC_Match_Gamma_Eta),
  fMC_Match_Gamma_Phi(original.fMC_Match_Gamma_Phi),
  fMC_Match_Gamma_Pt(original.fMC_Match_Gamma_Pt),
  fMC_Match_Gamma_Energy(original.fMC_Match_Gamma_Energy),
  fMC_Match_Gamma_Mass(original.fMC_Match_Gamma_Mass),
  fMC_Match_Gamma_OpeningAngle(original.fMC_Match_Gamma_OpeningAngle),
  fMC_Match_Gamma_R(original.fMC_Match_Gamma_R),
  fMC_Match_Gamma_Z_R(original.fMC_Match_Gamma_Z_R),
  fMC_Match_Gamma_X_Y(original.fMC_Match_Gamma_X_Y),
  fMC_Pi0_Eta(original.fMC_Pi0_Eta),
  fMC_Pi0_Phi(original.fMC_Pi0_Phi),
  fMC_Pi0_Pt(original.fMC_Pi0_Pt),
  fMC_Pi0_Energy(original.fMC_Pi0_Energy),
  fMC_Pi0_Mass(original.fMC_Pi0_Mass),
  fMC_Pi0_OpeningAngleGamma(original.fMC_Pi0_OpeningAngleGamma),
  fMC_Pi0_R(original.fMC_Pi0_R),
  fMC_Pi0_Z_R(original.fMC_Pi0_Z_R),
  fMC_Pi0_X_Y(original.fMC_Pi0_X_Y),
  fMC_Pi0Secondaries_Eta(original.fMC_Pi0_Eta),
  fMC_Pi0Secondaries_Phi(original.fMC_Pi0_Phi),
  fMC_Pi0Secondaries_Pt(original.fMC_Pi0_Pt),
  fMC_Pi0Secondaries_Energy(original.fMC_Pi0_Energy),
  fMC_Pi0Secondaries_Mass(original.fMC_Pi0_Mass),
  fMC_Pi0Secondaries_OpeningAngleGamma(original.fMC_Pi0_OpeningAngleGamma),
  fMC_Pi0Secondaries_R(original.fMC_Pi0_R),
  fMC_Pi0Secondaries_Z_R(original.fMC_Pi0_Z_R),
  fMC_Pi0Secondaries_X_Y(original.fMC_Pi0Secondaries_X_Y),
  fMC_Eta_Eta(original.fMC_Eta_Eta),
  fMC_Eta_Phi(original.fMC_Eta_Phi),
  fMC_Eta_Pt(original.fMC_Eta_Pt),
  fMC_Eta_Energy(original.fMC_Eta_Energy),
  fMC_Eta_Mass(original.fMC_Eta_Mass),
  fMC_Eta_OpeningAngleGamma(original.fMC_Eta_OpeningAngleGamma),
  fMC_Eta_R(original.fMC_Eta_R),
  fMC_Eta_Z_R(original.fMC_Eta_Z_R),
  fMC_Eta_X_Y(original.fMC_Eta_X_Y),
  fESD_EP_R(original.fESD_EP_R),
  fESD_EP_Z_R(original.fESD_EP_Z_R),
  fESD_EP_X_Y(original.fESD_EP_X_Y),
  fESD_EP_OpeningAngle(original.fESD_EP_OpeningAngle),
  fESD_E_Energy(original.fESD_E_Energy),
  fESD_E_Pt(original.fESD_E_Pt),
  fESD_E_Eta(original.fESD_E_Eta),
  fESD_E_Phi(original.fESD_E_Phi),
  fESD_P_Energy(original.fESD_P_Energy),
  fESD_P_Pt(original.fESD_P_Pt),
  fESD_P_Eta(original.fESD_P_Eta),
  fESD_P_Phi(original.fESD_P_Phi),
  fESD_Gamma_Energy(original.fESD_Gamma_Energy),
  fESD_Gamma_Pt(original.fESD_Gamma_Pt),
  fESD_Gamma_Eta(original.fESD_Gamma_Eta),
  fESD_Gamma_Phi(original.fESD_Gamma_Phi),
  fESD_Mapping(original.fESD_Mapping),
  fESD_Mapping_Phi(original.fESD_Mapping_Phi),
  fESD_Mapping_R(original.fESD_Mapping_R),
  fESD_Match_Gamma_OpeningAngle(original.fESD_Match_Gamma_OpeningAngle),
  fESD_Match_Gamma_Energy(original.fESD_Match_Gamma_Energy),
  fESD_Match_Gamma_Pt(original.fESD_Match_Gamma_Pt),
  fESD_Match_Gamma_Eta(original.fESD_Match_Gamma_Eta),
  fESD_Match_Gamma_Phi(original.fESD_Match_Gamma_Phi),
  fESD_Match_Gamma_Mass(original.fESD_Match_Gamma_Mass),
  fESD_Match_Gamma_Width(original.fESD_Match_Gamma_Width),
  fESD_Match_Gamma_Chi2(original.fESD_Match_Gamma_Chi2),
  fESD_Match_Gamma_NDF(original.fESD_Match_Gamma_NDF),
  fESD_Match_Gamma_R(original.fESD_Match_Gamma_R),
  fESD_Match_Gamma_Z_R(original.fESD_Match_Gamma_Z_R),
  fESD_Match_Gamma_X_Y(original.fESD_Match_Gamma_X_Y),
  fESD_Pi0_OpeningAngleGamma(original.fESD_Pi0_OpeningAngleGamma),
  fESD_Pi0_Energy(original.fESD_Pi0_Energy),
  fESD_Pi0_Pt(original.fESD_Pi0_Pt),
  fESD_Pi0_Eta(original.fESD_Pi0_Eta),
  fESD_Pi0_Phi(original.fESD_Pi0_Phi),
  fESD_Pi0_Mass(original.fESD_Pi0_Mass),
  fESD_Pi0_R(original.fESD_Pi0_R),
  fESD_Pi0_Z_R(original.fESD_Pi0_Z_R),
  fESD_Pi0_X_Y(original.fESD_Pi0_X_Y),
  fESD_Eta_OpeningAngleGamma(original.fESD_Eta_OpeningAngleGamma),
  fESD_Eta_Energy(original.fESD_Eta_Energy),
  fESD_Eta_Pt(original.fESD_Eta_Pt),
  fESD_Eta_Eta(original.fESD_Eta_Eta),
  fESD_Eta_Phi(original.fESD_Eta_Phi),
  fESD_Eta_Mass(original.fESD_Eta_Mass),
  fESD_Eta_R(original.fESD_Eta_R),
  fESD_Eta_Z_R(original.fESD_Eta_Z_R),
  fESD_Eta_X_Y(original.fESD_Eta_X_Y),
  fESD_Background_OpeningAngleGamma(original.fESD_Background_OpeningAngleGamma),
  fESD_Background_Energy(original.fESD_Background_Energy),
  fESD_Background_Pt(original.fESD_Background_Pt),
  fESD_Background_Eta(original.fESD_Background_Eta),
  fESD_Background_Phi(original.fESD_Background_Phi),
  fESD_Background_Mass(original.fESD_Background_Mass),
  fESD_Background_R(original.fESD_Background_R),
  fESD_Background_Z_R(original.fESD_Background_Z_R),
  fESD_Background_X_Y(original.fESD_Background_X_Y),
  fResolution_dPt(original.fResolution_dPt),
  fResolution_dR(original.fResolution_dR),
  fResolution_dZ(original.fResolution_dZ),
  fResolution_dR_dPt(original.fResolution_dR_dPt),
  fResolution_MC_Pt(original.fResolution_MC_Pt),
  fResolution_MC_R(original.fResolution_MC_R),
  fResolution_MC_Z(original.fResolution_MC_Z),
  fResolution_ESD_Pt(original.fResolution_ESD_Pt),
  fResolution_ESD_R(original.fResolution_ESD_R),
  fResolution_ESD_Z(original.fResolution_ESD_Z),
  fNumberOfV0s(original.fNumberOfV0s),
  fNumberOfSurvivingV0s(original.fNumberOfSurvivingV0s),
  fV0MassDebugCut1(original.fV0MassDebugCut1),
  fV0MassDebugCut2(original.fV0MassDebugCut2),
  fV0MassDebugCut3(original.fV0MassDebugCut3),
  fV0MassDebugCut4(original.fV0MassDebugCut4),
  fV0MassDebugCut5(original.fV0MassDebugCut5),
  fV0MassDebugCut6(original.fV0MassDebugCut6),
  fV0MassDebugCut7(original.fV0MassDebugCut7),
  fV0MassDebugCut8(original.fV0MassDebugCut8)
{    

}


AliGammaConversionHistograms & AliGammaConversionHistograms::operator = (const AliGammaConversionHistograms & /*source*/)
{
  // assignment operator
  return *this;
}


AliGammaConversionHistograms::~AliGammaConversionHistograms() {
  if(fOutputContainer != NULL){ delete fOutputContainer;}

  if(fMC_EP_R != NULL){ delete fMC_EP_R;}
  if(fMC_EP_Z_R != NULL){ delete fMC_EP_Z_R;}
  if(fMC_EP_X_Y != NULL){ delete fMC_EP_X_Y;}
  if(fMC_EP_OpeningAngle != NULL){ delete fMC_EP_OpeningAngle;}

  if(fMC_E_Energy != NULL){ delete fMC_E_Energy;}
  if(fMC_E_Pt != NULL){ delete fMC_E_Pt;}
  if(fMC_E_Eta != NULL){ delete fMC_E_Eta;}
  if(fMC_E_Phi != NULL){ delete fMC_E_Phi;}

  if(fMC_P_Energy != NULL){ delete fMC_P_Energy;}
  if(fMC_P_Pt != NULL){ delete fMC_P_Pt;}
  if(fMC_P_Eta != NULL){ delete fMC_P_Eta;}
  if(fMC_P_Phi != NULL){ delete fMC_P_Phi;}

  if(fMC_Gamma_Energy != NULL){ delete fMC_Gamma_Energy;}
  if(fMC_Gamma_Pt != NULL){ delete fMC_Gamma_Pt;}
  if(fMC_Gamma_Eta != NULL){ delete fMC_Gamma_Eta;}
  if(fMC_Gamma_Phi != NULL){ delete fMC_Gamma_Phi;}

  if(fMC_DirectGamma_Energy != NULL){ delete fMC_DirectGamma_Energy;}
  if(fMC_DirectGamma_Pt != NULL){ delete fMC_DirectGamma_Pt;}
  if(fMC_DirectGamma_Eta != NULL){ delete fMC_DirectGamma_Eta;}
  if(fMC_DirectGamma_Phi != NULL){ delete fMC_DirectGamma_Phi;}

  //mapping
  for(UInt_t i=0;i<fMC_Mapping.size();i++){
    for(UInt_t j=0;j<fMC_Mapping[i].size();j++){
      if(fMC_Mapping[i][j] != NULL){delete fMC_Mapping[i][j];}
      fMC_Mapping[i][j]=NULL;
    }
    fMC_Mapping[i].clear();
  }
  fMC_Mapping.clear();

  for(UInt_t i=0;i<fMC_Mapping_Phi.size();i++){
    if(fMC_Mapping_Phi[i] != NULL){delete fMC_Mapping_Phi[i];}
    fMC_Mapping_Phi[i]=NULL;
  }
  fMC_Mapping_Phi.clear();

  for(UInt_t i=0;i<fMC_Mapping_R.size();i++){
    if(fMC_Mapping_R[i] != NULL){delete fMC_Mapping_R[i];}
    fMC_Mapping_R[i]=NULL;
  }
  fMC_Mapping_R.clear();


  if(fMC_Match_Gamma_Eta != NULL){ delete fMC_Match_Gamma_Eta;}
  if(fMC_Match_Gamma_Phi != NULL){ delete fMC_Match_Gamma_Phi;}
  if(fMC_Match_Gamma_Pt != NULL){ delete fMC_Match_Gamma_Pt;}
  if(fMC_Match_Gamma_Energy != NULL){ delete fMC_Match_Gamma_Energy;}
  if(fMC_Match_Gamma_Mass != NULL){ delete fMC_Match_Gamma_Mass;}
  if(fMC_Match_Gamma_OpeningAngle != NULL){ delete fMC_Match_Gamma_OpeningAngle;}
  if(fMC_Match_Gamma_R != NULL){ delete fMC_Match_Gamma_R;}
  if(fMC_Match_Gamma_Z_R != NULL){ delete fMC_Match_Gamma_Z_R;}
  if(fMC_Match_Gamma_X_Y != NULL){ delete fMC_Match_Gamma_X_Y;}

  if(fMC_Pi0_Eta != NULL){ delete fMC_Pi0_Eta;}
  if(fMC_Pi0_Phi != NULL){ delete fMC_Pi0_Phi;}
  if(fMC_Pi0_Pt != NULL){ delete fMC_Pi0_Pt;}
  if(fMC_Pi0_Energy != NULL){ delete fMC_Pi0_Energy;}
  if(fMC_Pi0_Mass != NULL){ delete fMC_Pi0_Mass;}
  if(fMC_Pi0_OpeningAngleGamma != NULL){ delete fMC_Pi0_OpeningAngleGamma;}
  if(fMC_Pi0_R != NULL){ delete fMC_Pi0_R;}
  if(fMC_Pi0_Z_R != NULL){ delete fMC_Pi0_Z_R;}
  if(fMC_Pi0_X_Y != NULL){ delete fMC_Pi0_X_Y;}
  if(fMC_Pi0Secondaries_X_Y != NULL){ delete fMC_Pi0Secondaries_X_Y;}

  if(fMC_Eta_Eta != NULL){ delete fMC_Eta_Eta;}
  if(fMC_Eta_Phi != NULL){ delete fMC_Eta_Phi;}
  if(fMC_Eta_Pt != NULL){ delete fMC_Eta_Pt;}
  if(fMC_Eta_Energy != NULL){ delete fMC_Eta_Energy;}
  if(fMC_Eta_Mass != NULL){ delete fMC_Eta_Mass;}
  if(fMC_Eta_OpeningAngleGamma != NULL){ delete fMC_Eta_OpeningAngleGamma;}
  if(fMC_Eta_R != NULL){ delete fMC_Eta_R;}
  if(fMC_Eta_Z_R != NULL){ delete fMC_Eta_Z_R;}
  if(fMC_Eta_X_Y != NULL){ delete fMC_Eta_X_Y;}
    
  // Histograms from esd tracks
  if(fESD_EP_R != NULL){ delete fESD_EP_R;}
  if(fESD_EP_Z_R != NULL){ delete fESD_EP_Z_R;}
  if(fESD_EP_X_Y != NULL){ delete fESD_EP_X_Y;}
  if(fESD_EP_OpeningAngle != NULL){ delete fESD_EP_OpeningAngle;}

  if(fESD_E_Energy != NULL){ delete fESD_E_Energy;}
  if(fESD_E_Pt != NULL){ delete fESD_E_Pt;}
  if(fESD_E_Eta != NULL){ delete fESD_E_Eta;}
  if(fESD_E_Phi != NULL){ delete fESD_E_Phi;}

  if(fESD_P_Energy != NULL){ delete fESD_P_Energy;}
  if(fESD_P_Pt != NULL){ delete fESD_P_Pt;}
  if(fESD_P_Eta != NULL){ delete fESD_P_Eta;}
  if(fESD_P_Phi != NULL){ delete fESD_P_Phi;}

  if(fESD_Gamma_Energy != NULL){ delete fESD_Gamma_Energy;}
  if(fESD_Gamma_Pt != NULL){ delete fESD_Gamma_Pt;}
  if(fESD_Gamma_Eta != NULL){ delete fESD_Gamma_Eta;}
  if(fESD_Gamma_Phi != NULL){ delete fESD_Gamma_Phi;}

  //mapping
  for(UInt_t i=0;i<fESD_Mapping.size();i++){
    for(UInt_t j=0;j<fESD_Mapping[i].size();j++){
      if(fESD_Mapping[i][j] != NULL){delete fESD_Mapping[i][j];}
      fESD_Mapping[i][j]=NULL;
    }
    fESD_Mapping[i].clear();
  }
  fESD_Mapping.clear();

  for(UInt_t i=0;i<fESD_Mapping_Phi.size();i++){
    if(fESD_Mapping_Phi[i] != NULL){delete fESD_Mapping_Phi[i];}
    fESD_Mapping_Phi[i]=NULL;
  }
  fESD_Mapping_Phi.clear();

  for(UInt_t i=0;i<fESD_Mapping_R.size();i++){
    if(fESD_Mapping_R[i] != NULL){delete fESD_Mapping_R[i];}
    fESD_Mapping_R[i]=NULL;
  }
  fESD_Mapping_R.clear();

  if(fESD_Match_Gamma_OpeningAngle != NULL){ delete fESD_Match_Gamma_OpeningAngle;}
  if(fESD_Match_Gamma_Energy != NULL){ delete fESD_Match_Gamma_Energy;}
  if(fESD_Match_Gamma_Pt != NULL){ delete fESD_Match_Gamma_Pt;}
  if(fESD_Match_Gamma_Eta != NULL){ delete fESD_Match_Gamma_Eta;}
  if(fESD_Match_Gamma_Phi != NULL){ delete fESD_Match_Gamma_Phi;}
  if(fESD_Match_Gamma_Mass != NULL){ delete fESD_Match_Gamma_Mass;}
  if(fESD_Match_Gamma_Width != NULL){ delete fESD_Match_Gamma_Width;}
  if(fESD_Match_Gamma_Chi2 != NULL){ delete fESD_Match_Gamma_Chi2;}
  if(fESD_Match_Gamma_NDF != NULL){ delete fESD_Match_Gamma_NDF;}
  if(fESD_Match_Gamma_R != NULL){ delete fESD_Match_Gamma_R;}
  if(fESD_Match_Gamma_Z_R != NULL){ delete fESD_Match_Gamma_Z_R;}
  if(fESD_Match_Gamma_X_Y != NULL){ delete fESD_Match_Gamma_X_Y;}

  if(fESD_Pi0_OpeningAngleGamma != NULL){ delete fESD_Pi0_OpeningAngleGamma;}
  if(fESD_Pi0_Energy != NULL){ delete fESD_Pi0_Energy;}
  if(fESD_Pi0_Pt != NULL){ delete fESD_Pi0_Pt;}
  if(fESD_Pi0_Eta != NULL){ delete fESD_Pi0_Eta;}
  if(fESD_Pi0_Phi != NULL){ delete fESD_Pi0_Phi;}
  if(fESD_Pi0_Mass != NULL){ delete fESD_Pi0_Mass;}
  if(fESD_Pi0_R != NULL){ delete fESD_Pi0_R;}
  if(fESD_Pi0_Z_R != NULL){ delete fESD_Pi0_Z_R;}
  if(fESD_Pi0_X_Y != NULL){ delete fESD_Pi0_X_Y;}

  if(fESD_Eta_OpeningAngleGamma != NULL){ delete fESD_Eta_OpeningAngleGamma;}
  if(fESD_Eta_Energy != NULL){ delete fESD_Eta_Energy;}
  if(fESD_Eta_Pt != NULL){ delete fESD_Eta_Pt;}
  if(fESD_Eta_Eta != NULL){ delete fESD_Eta_Eta;}
  if(fESD_Eta_Phi != NULL){ delete fESD_Eta_Phi;}
  if(fESD_Eta_Mass != NULL){ delete fESD_Eta_Mass;}
  if(fESD_Eta_R != NULL){ delete fESD_Eta_R;}
  if(fESD_Eta_Z_R != NULL){ delete fESD_Eta_Z_R;}
  if(fESD_Eta_X_Y != NULL){ delete fESD_Eta_X_Y;}

  if(fESD_Background_OpeningAngleGamma != NULL){ delete fESD_Background_OpeningAngleGamma;}
  if(fESD_Background_Energy != NULL){ delete fESD_Background_Energy;}
  if(fESD_Background_Pt != NULL){ delete fESD_Background_Pt;}
  if(fESD_Background_Eta != NULL){ delete fESD_Background_Eta;}
  if(fESD_Background_Phi != NULL){ delete fESD_Background_Phi;}
  if(fESD_Background_Mass != NULL){ delete fESD_Background_Mass;}
  if(fESD_Background_R != NULL){ delete fESD_Background_R;}
  if(fESD_Background_Z_R != NULL){ delete fESD_Background_Z_R;}
  if(fESD_Background_X_Y != NULL){ delete fESD_Background_X_Y;}

  if(fResolution_dPt != NULL){ delete fResolution_dPt;}
  if(fResolution_dR != NULL){ delete fResolution_dR;}
  if(fResolution_dZ != NULL){ delete fResolution_dZ;}
  if(fResolution_dR_dPt != NULL){ delete fResolution_dR_dPt;}
  if(fResolution_MC_Pt != NULL){ delete fResolution_MC_Pt;}
  if(fResolution_MC_R != NULL){ delete fResolution_MC_R;}
  if(fResolution_MC_Z != NULL){ delete fResolution_MC_Z;}
  if(fResolution_ESD_Pt != NULL){ delete fResolution_ESD_Pt;}
  if(fResolution_ESD_R != NULL){ delete fResolution_ESD_R;}
  if(fResolution_ESD_Z != NULL){ delete fResolution_ESD_Z;}
  
  if(fNumberOfV0s != NULL){delete fNumberOfV0s;}
  if(fNumberOfSurvivingV0s != NULL){delete fNumberOfSurvivingV0s;}
  if(fV0MassDebugCut1 != NULL){delete fV0MassDebugCut1;}
  if(fV0MassDebugCut2 != NULL){delete fV0MassDebugCut2;}
  if(fV0MassDebugCut3 != NULL){delete fV0MassDebugCut3;}
  if(fV0MassDebugCut4 != NULL){delete fV0MassDebugCut4;}
  if(fV0MassDebugCut5 != NULL){delete fV0MassDebugCut5;}
  if(fV0MassDebugCut6 != NULL){delete fV0MassDebugCut6;}
  if(fV0MassDebugCut7 != NULL){delete fV0MassDebugCut7;}
  if(fV0MassDebugCut8 != NULL){delete fV0MassDebugCut8;}

}


TList * AliGammaConversionHistograms::GetOutputContainer(){
  //checking if the container is alrerady created
  if(fOutputContainer != NULL){
    delete fOutputContainer;
    fOutputContainer=NULL;
  }
  fOutputContainer = new TList();
  TList*  fMappingContainer = new TList();
  fMappingContainer->SetName("Mapping Histograms");

  if(fMC_EP_R != NULL){ fOutputContainer->Add(fMC_EP_R);}
  if(fMC_EP_Z_R != NULL){ fOutputContainer->Add(fMC_EP_Z_R);}
  if(fMC_EP_X_Y != NULL){ fOutputContainer->Add(fMC_EP_X_Y);}
  if(fMC_EP_OpeningAngle != NULL){ fOutputContainer->Add(fMC_EP_OpeningAngle);}

  if(fMC_E_Energy != NULL){ fOutputContainer->Add(fMC_E_Energy);}
  if(fMC_E_Pt != NULL){ fOutputContainer->Add(fMC_E_Pt);}
  if(fMC_E_Eta != NULL){ fOutputContainer->Add(fMC_E_Eta);}
  if(fMC_E_Phi != NULL){ fOutputContainer->Add(fMC_E_Phi);}

  if(fMC_P_Energy != NULL){ fOutputContainer->Add(fMC_P_Energy);}
  if(fMC_P_Pt != NULL){ fOutputContainer->Add(fMC_P_Pt);}
  if(fMC_P_Eta != NULL){ fOutputContainer->Add(fMC_P_Eta);}
  if(fMC_P_Phi != NULL){ fOutputContainer->Add(fMC_P_Phi);}

  if(fMC_Gamma_Energy != NULL){ fOutputContainer->Add(fMC_Gamma_Energy);}
  if(fMC_Gamma_Pt != NULL){ fOutputContainer->Add(fMC_Gamma_Pt);}
  if(fMC_Gamma_Eta != NULL){ fOutputContainer->Add(fMC_Gamma_Eta);}
  if(fMC_Gamma_Phi != NULL){ fOutputContainer->Add(fMC_Gamma_Phi);}

  if(fMC_DirectGamma_Energy != NULL){ fOutputContainer->Add(fMC_DirectGamma_Energy);}
  if(fMC_DirectGamma_Pt != NULL){ fOutputContainer->Add(fMC_DirectGamma_Pt);}
  if(fMC_DirectGamma_Eta != NULL){ fOutputContainer->Add(fMC_DirectGamma_Eta);}
  if(fMC_DirectGamma_Phi != NULL){ fOutputContainer->Add(fMC_DirectGamma_Phi);}

  //mapping
  for(UInt_t i=0;i<fMC_Mapping.size();i++){
    for(UInt_t j=0;j<fMC_Mapping[i].size();j++){
      if(fMC_Mapping[i][j] != NULL){fMappingContainer->Add(fMC_Mapping[i][j]);}
    }
  }
  for(UInt_t i=0;i<fMC_Mapping_Phi.size();i++){
    if(fMC_Mapping_Phi[i] != NULL){fMappingContainer->Add(fMC_Mapping_Phi[i]);}
  }
  for(UInt_t i=0;i<fMC_Mapping_R.size();i++){
    if(fMC_Mapping_R[i] != NULL){fMappingContainer->Add(fMC_Mapping_R[i]);}
  }
  if(fMC_Match_Gamma_Eta != NULL){ fOutputContainer->Add(fMC_Match_Gamma_Eta);}
  if(fMC_Match_Gamma_Phi != NULL){ fOutputContainer->Add(fMC_Match_Gamma_Phi);}
  if(fMC_Match_Gamma_Pt != NULL){ fOutputContainer->Add(fMC_Match_Gamma_Pt);}
  if(fMC_Match_Gamma_Energy != NULL){ fOutputContainer->Add(fMC_Match_Gamma_Energy);}
  if(fMC_Match_Gamma_Mass != NULL){ fOutputContainer->Add(fMC_Match_Gamma_Mass);}
  if(fMC_Match_Gamma_OpeningAngle != NULL){ fOutputContainer->Add(fMC_Match_Gamma_OpeningAngle);}
  if(fMC_Match_Gamma_R != NULL){ fOutputContainer->Add(fMC_Match_Gamma_R);}
  if(fMC_Match_Gamma_Z_R != NULL){ fOutputContainer->Add(fMC_Match_Gamma_Z_R);}
  if(fMC_Match_Gamma_X_Y != NULL){ fOutputContainer->Add(fMC_Match_Gamma_X_Y);}

  if(fMC_Pi0_Eta != NULL){ fOutputContainer->Add(fMC_Pi0_Eta);}
  if(fMC_Pi0_Phi != NULL){ fOutputContainer->Add(fMC_Pi0_Phi);}
  if(fMC_Pi0_Pt != NULL){ fOutputContainer->Add(fMC_Pi0_Pt);}
  if(fMC_Pi0_Energy != NULL){ fOutputContainer->Add(fMC_Pi0_Energy);}
  if(fMC_Pi0_Mass != NULL){ fOutputContainer->Add(fMC_Pi0_Mass);}
  if(fMC_Pi0_OpeningAngleGamma != NULL){ fOutputContainer->Add(fMC_Pi0_OpeningAngleGamma);}
  if(fMC_Pi0_R != NULL){ fOutputContainer->Add(fMC_Pi0_R);}
  if(fMC_Pi0_Z_R != NULL){ fOutputContainer->Add(fMC_Pi0_Z_R);}
  if(fMC_Pi0_X_Y != NULL){ fOutputContainer->Add(fMC_Pi0_X_Y);}
  if(fMC_Pi0Secondaries_X_Y != NULL){ fOutputContainer->Add(fMC_Pi0Secondaries_X_Y);}

  if(fMC_Eta_Eta != NULL){ fOutputContainer->Add(fMC_Eta_Eta);}
  if(fMC_Eta_Phi != NULL){ fOutputContainer->Add(fMC_Eta_Phi);}
  if(fMC_Eta_Pt != NULL){ fOutputContainer->Add(fMC_Eta_Pt);}
  if(fMC_Eta_Energy != NULL){ fOutputContainer->Add(fMC_Eta_Energy);}
  if(fMC_Eta_Mass != NULL){ fOutputContainer->Add(fMC_Eta_Mass);}
  if(fMC_Eta_OpeningAngleGamma != NULL){ fOutputContainer->Add(fMC_Eta_OpeningAngleGamma);}
  if(fMC_Eta_R != NULL){ fOutputContainer->Add(fMC_Eta_R);}
  if(fMC_Eta_Z_R != NULL){ fOutputContainer->Add(fMC_Eta_Z_R);}
  if(fMC_Eta_X_Y != NULL){ fOutputContainer->Add(fMC_Eta_X_Y);}
    
  // Histograms from esd tracks
  if(fESD_EP_R != NULL){ fOutputContainer->Add(fESD_EP_R);}
  if(fESD_EP_Z_R != NULL){ fOutputContainer->Add(fESD_EP_Z_R);}
  if(fESD_EP_X_Y != NULL){ fOutputContainer->Add(fESD_EP_X_Y);}
  if(fESD_EP_OpeningAngle != NULL){ fOutputContainer->Add(fESD_EP_OpeningAngle);}

  if(fESD_E_Energy != NULL){ fOutputContainer->Add(fESD_E_Energy);}
  if(fESD_E_Pt != NULL){ fOutputContainer->Add(fESD_E_Pt);}
  if(fESD_E_Eta != NULL){ fOutputContainer->Add(fESD_E_Eta);}
  if(fESD_E_Phi != NULL){ fOutputContainer->Add(fESD_E_Phi);}

  if(fESD_P_Energy != NULL){ fOutputContainer->Add(fESD_P_Energy);}
  if(fESD_P_Pt != NULL){ fOutputContainer->Add(fESD_P_Pt);}
  if(fESD_P_Eta != NULL){ fOutputContainer->Add(fESD_P_Eta);}
  if(fESD_P_Phi != NULL){ fOutputContainer->Add(fESD_P_Phi);}

  if(fESD_Gamma_Energy != NULL){ fOutputContainer->Add(fESD_Gamma_Energy);}
  if(fESD_Gamma_Pt != NULL){ fOutputContainer->Add(fESD_Gamma_Pt);}
  if(fESD_Gamma_Eta != NULL){ fOutputContainer->Add(fESD_Gamma_Eta);}
  if(fESD_Gamma_Phi != NULL){ fOutputContainer->Add(fESD_Gamma_Phi);}

  //mapping
  for(UInt_t i=0;i<fESD_Mapping.size();i++){
    for(UInt_t j=0;j<fESD_Mapping[i].size();j++){
      if(fESD_Mapping[i][j] != NULL){fMappingContainer->Add(fESD_Mapping[i][j]);}
    }
  }
  for(UInt_t i=0;i<fESD_Mapping_Phi.size();i++){
    if(fESD_Mapping_Phi[i] != NULL){fMappingContainer->Add(fESD_Mapping_Phi[i]);}
  }
  for(UInt_t i=0;i<fESD_Mapping_R.size();i++){
    if(fESD_Mapping_R[i] != NULL){fMappingContainer->Add(fESD_Mapping_R[i]);}
  }

  fOutputContainer->Add(fMappingContainer);

  if(fESD_Match_Gamma_OpeningAngle != NULL){ fOutputContainer->Add(fESD_Match_Gamma_OpeningAngle);}
  if(fESD_Match_Gamma_Energy != NULL){ fOutputContainer->Add(fESD_Match_Gamma_Energy);}
  if(fESD_Match_Gamma_Pt != NULL){ fOutputContainer->Add(fESD_Match_Gamma_Pt);}
  if(fESD_Match_Gamma_Eta != NULL){ fOutputContainer->Add(fESD_Match_Gamma_Eta);}
  if(fESD_Match_Gamma_Phi != NULL){ fOutputContainer->Add(fESD_Match_Gamma_Phi);}
  if(fESD_Match_Gamma_Mass != NULL){ fOutputContainer->Add(fESD_Match_Gamma_Mass);}
  if(fESD_Match_Gamma_Width != NULL){ fOutputContainer->Add(fESD_Match_Gamma_Width);}
  if(fESD_Match_Gamma_Chi2 != NULL){ fOutputContainer->Add(fESD_Match_Gamma_Chi2);}
  if(fESD_Match_Gamma_NDF != NULL){ fOutputContainer->Add(fESD_Match_Gamma_NDF);}
  if(fESD_Match_Gamma_R != NULL){ fOutputContainer->Add(fESD_Match_Gamma_R);}
  if(fESD_Match_Gamma_Z_R != NULL){ fOutputContainer->Add(fESD_Match_Gamma_Z_R);}
  if(fESD_Match_Gamma_X_Y != NULL){ fOutputContainer->Add(fESD_Match_Gamma_X_Y);}

  if(fESD_Pi0_OpeningAngleGamma != NULL){ fOutputContainer->Add(fESD_Pi0_OpeningAngleGamma);}
  if(fESD_Pi0_Energy != NULL){ fOutputContainer->Add(fESD_Pi0_Energy);}
  if(fESD_Pi0_Pt != NULL){ fOutputContainer->Add(fESD_Pi0_Pt);}
  if(fESD_Pi0_Eta != NULL){ fOutputContainer->Add(fESD_Pi0_Eta);}
  if(fESD_Pi0_Phi != NULL){ fOutputContainer->Add(fESD_Pi0_Phi);}
  if(fESD_Pi0_Mass != NULL){ fOutputContainer->Add(fESD_Pi0_Mass);}
  if(fESD_Pi0_R != NULL){ fOutputContainer->Add(fESD_Pi0_R);}
  if(fESD_Pi0_Z_R != NULL){ fOutputContainer->Add(fESD_Pi0_Z_R);}
  if(fESD_Pi0_X_Y != NULL){ fOutputContainer->Add(fESD_Pi0_X_Y);}

  if(fESD_Eta_OpeningAngleGamma != NULL){ fOutputContainer->Add(fESD_Eta_OpeningAngleGamma);}
  if(fESD_Eta_Energy != NULL){ fOutputContainer->Add(fESD_Eta_Energy);}
  if(fESD_Eta_Pt != NULL){ fOutputContainer->Add(fESD_Eta_Pt);}
  if(fESD_Eta_Eta != NULL){ fOutputContainer->Add(fESD_Eta_Eta);}
  if(fESD_Eta_Phi != NULL){ fOutputContainer->Add(fESD_Eta_Phi);}
  if(fESD_Eta_Mass != NULL){ fOutputContainer->Add(fESD_Eta_Mass);}
  if(fESD_Eta_R != NULL){ fOutputContainer->Add(fESD_Eta_R);}
  if(fESD_Eta_Z_R != NULL){ fOutputContainer->Add(fESD_Eta_Z_R);}
  if(fESD_Eta_X_Y != NULL){ fOutputContainer->Add(fESD_Eta_X_Y);}

  if(fESD_Background_OpeningAngleGamma != NULL){ fOutputContainer->Add(fESD_Background_OpeningAngleGamma);}
  if(fESD_Background_Energy != NULL){ fOutputContainer->Add(fESD_Background_Energy);}
  if(fESD_Background_Pt != NULL){ fOutputContainer->Add(fESD_Background_Pt);}
  if(fESD_Background_Eta != NULL){ fOutputContainer->Add(fESD_Background_Eta);}
  if(fESD_Background_Phi != NULL){ fOutputContainer->Add(fESD_Background_Phi);}
  if(fESD_Background_Mass != NULL){ fOutputContainer->Add(fESD_Background_Mass);}
  if(fESD_Background_R != NULL){ fOutputContainer->Add(fESD_Background_R);}
  if(fESD_Background_Z_R != NULL){ fOutputContainer->Add(fESD_Background_Z_R);}
  if(fESD_Background_X_Y != NULL){ fOutputContainer->Add(fESD_Background_X_Y);}

  if(fResolution_dPt != NULL){ fOutputContainer->Add(fResolution_dPt);}
  if(fResolution_dR != NULL){ fOutputContainer->Add(fResolution_dR);}
  if(fResolution_dZ != NULL){ fOutputContainer->Add(fResolution_dZ);}
  if(fResolution_dR_dPt != NULL){ fOutputContainer->Add(fResolution_dR_dPt);}
  if(fResolution_MC_Pt != NULL){ fOutputContainer->Add(fResolution_MC_Pt);}
  if(fResolution_MC_R != NULL){ fOutputContainer->Add(fResolution_MC_R);}
  if(fResolution_MC_Z != NULL){ fOutputContainer->Add(fResolution_MC_Z);}
  if(fResolution_ESD_Pt != NULL){ fOutputContainer->Add(fResolution_ESD_Pt);}
  if(fResolution_ESD_R != NULL){ fOutputContainer->Add(fResolution_ESD_R);}
  if(fResolution_ESD_Z != NULL){ fOutputContainer->Add(fResolution_ESD_Z);}

  if(fNumberOfV0s != NULL){fOutputContainer->Add(fNumberOfV0s);}
  if(fNumberOfSurvivingV0s != NULL){fOutputContainer->Add(fNumberOfSurvivingV0s);}
  if(fV0MassDebugCut1 != NULL){fOutputContainer->Add(fV0MassDebugCut1);}
  if(fV0MassDebugCut2 != NULL){fOutputContainer->Add(fV0MassDebugCut2);}
  if(fV0MassDebugCut3 != NULL){fOutputContainer->Add(fV0MassDebugCut3);}
  if(fV0MassDebugCut4 != NULL){fOutputContainer->Add(fV0MassDebugCut4);}
  if(fV0MassDebugCut5 != NULL){fOutputContainer->Add(fV0MassDebugCut5);}
  if(fV0MassDebugCut6 != NULL){fOutputContainer->Add(fV0MassDebugCut6);}
  if(fV0MassDebugCut7 != NULL){fOutputContainer->Add(fV0MassDebugCut7);}
  if(fV0MassDebugCut8 != NULL){fOutputContainer->Add(fV0MassDebugCut8);}

  return fOutputContainer;
}

Int_t AliGammaConversionHistograms::GetRBin(Double_t radius){
  Int_t iResult=0;
  if(fDeltaR>0){
    iResult = (Int_t)((radius - fMinRadius)/fDeltaR);
  }
  return iResult;
}

Int_t AliGammaConversionHistograms::GetPhiBin(Double_t phi){
  Int_t iResult=0;
  if(fDeltaPhi>0){
    if(phi>TMath::Pi()){
      phi-=2*TMath::Pi();
    }
    iResult = (Int_t)((phi - fMinPhi)/fDeltaPhi);
  }
  return iResult;
}



// Initializing the valuse for the mapping
void AliGammaConversionHistograms::Initialize_MappingValues(Int_t nPhiIndex, Int_t nRIndex, Int_t nBinsR, Double_t minRadius, Double_t maxRadius,Int_t nBinsPhi, Double_t minPhi, Double_t maxPhi){
  fNPhiIndex = nPhiIndex;
  fNRIndex   = nRIndex;
  fMinRadius      = minRadius;
  fMaxRadius      = maxRadius;
  if(nBinsR>0){
    fDeltaR       = (fMaxRadius - fMinRadius)/nRIndex;
  }
  fMinPhi         = minPhi;
  fMaxPhi         = maxPhi;
  if(nBinsPhi>0){
    fDeltaPhi     = (fMaxPhi-fMinPhi)/nPhiIndex;
  }
}

//Initializing functions for the histograms
void AliGammaConversionHistograms::Initialize_MC_EP_R(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fMC_EP_R = new TH1F("MC_EP_R", "", nXBins, firstX, lastX);
  fMC_EP_R->GetXaxis()->SetTitle(xAxisTitle);
  fMC_EP_R->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_MC_EP_Z_R(Int_t nXBins, Double_t firstX, Double_t lastX, Int_t nYBins, Double_t firstY, Double_t lastY, TString xAxisTitle, TString yAxisTitle){
  fMC_EP_Z_R = new TH2F("MC_EP_Z_R", "", nXBins, firstX, lastX, nYBins, firstY, lastY);
  fMC_EP_Z_R->GetXaxis()->SetTitle(xAxisTitle);
  fMC_EP_Z_R->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_MC_EP_X_Y(Int_t nXBins, Double_t firstX, Double_t lastX, Int_t nYBins, Double_t firstY, Double_t lastY, TString xAxisTitle, TString yAxisTitle){
  fMC_EP_X_Y = new TH2F("MC_EP_X_Y", "", nXBins, firstX, lastX, nYBins, firstY, lastY);
  fMC_EP_X_Y->GetXaxis()->SetTitle(xAxisTitle);
  fMC_EP_X_Y->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_MC_EP_OpeningAngle(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fMC_EP_OpeningAngle = new TH1F("MC_EP_OpeningAngle", "", nXBins, firstX, lastX);
  fMC_EP_OpeningAngle->GetXaxis()->SetTitle(xAxisTitle);
  fMC_EP_OpeningAngle->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_MC_E_Energy(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fMC_E_Energy = new TH1F("MC_E_Energy", "", nXBins, firstX, lastX);
  fMC_E_Energy->GetXaxis()->SetTitle(xAxisTitle);
  fMC_E_Energy->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_MC_E_Pt(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fMC_E_Pt = new TH1F("MC_E_Pt", "", nXBins, firstX, lastX);
  fMC_E_Pt->GetXaxis()->SetTitle(xAxisTitle);
  fMC_E_Pt->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_MC_E_Eta(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fMC_E_Eta = new TH1F("MC_E_Eta", "", nXBins, firstX, lastX);
  fMC_E_Eta->GetXaxis()->SetTitle(xAxisTitle);
  fMC_E_Eta->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_MC_E_Phi(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fMC_E_Phi = new TH1F("MC_E_Phi", "", nXBins, firstX, lastX);
  fMC_E_Phi->GetXaxis()->SetTitle(xAxisTitle);
  fMC_E_Phi->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_MC_P_Energy(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fMC_P_Energy = new TH1F("MC_P_Energy", "", nXBins, firstX, lastX);
  fMC_P_Energy->GetXaxis()->SetTitle(xAxisTitle);
  fMC_P_Energy->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_MC_P_Pt(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fMC_P_Pt = new TH1F("MC_P_Pt", "", nXBins, firstX, lastX);
  fMC_P_Pt->GetXaxis()->SetTitle(xAxisTitle);
  fMC_P_Pt->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_MC_P_Eta(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fMC_P_Eta = new TH1F("MC_P_Eta", "", nXBins, firstX, lastX);
  fMC_P_Eta->GetXaxis()->SetTitle(xAxisTitle);
  fMC_P_Eta->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_MC_P_Phi(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fMC_P_Phi = new TH1F("MC_P_Phi", "", nXBins, firstX, lastX);
  fMC_P_Phi->GetXaxis()->SetTitle(xAxisTitle);
  fMC_P_Phi->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_MC_Gamma_Energy(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fMC_Gamma_Energy = new TH1F("MC_Gamma_Energy", "", nXBins, firstX, lastX);
  fMC_Gamma_Energy->GetXaxis()->SetTitle(xAxisTitle);
  fMC_Gamma_Energy->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_MC_Gamma_Pt(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fMC_Gamma_Pt = new TH1F("MC_Gamma_Pt", "", nXBins, firstX, lastX);
  fMC_Gamma_Pt->GetXaxis()->SetTitle(xAxisTitle);
  fMC_Gamma_Pt->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_MC_Gamma_Eta(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fMC_Gamma_Eta = new TH1F("MC_Gamma_Eta", "", nXBins, firstX, lastX);
  fMC_Gamma_Eta->GetXaxis()->SetTitle(xAxisTitle);
  fMC_Gamma_Eta->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_MC_Gamma_Phi(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fMC_Gamma_Phi = new TH1F("MC_Gamma_Phi", "", nXBins, firstX, lastX);
  fMC_Gamma_Phi->GetXaxis()->SetTitle(xAxisTitle);
  fMC_Gamma_Phi->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_MC_DirectGamma_Energy(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fMC_DirectGamma_Energy = new TH1F("MC_DirectGamma_Energy", "", nXBins, firstX, lastX);
  fMC_DirectGamma_Energy->GetXaxis()->SetTitle(xAxisTitle);
  fMC_DirectGamma_Energy->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_MC_DirectGamma_Pt(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fMC_DirectGamma_Pt = new TH1F("MC_DirectGamma_Pt", "", nXBins, firstX, lastX);
  fMC_DirectGamma_Pt->GetXaxis()->SetTitle(xAxisTitle);
  fMC_DirectGamma_Pt->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_MC_DirectGamma_Eta(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fMC_DirectGamma_Eta = new TH1F("MC_DirectGamma_Eta", "", nXBins, firstX, lastX);
  fMC_DirectGamma_Eta->GetXaxis()->SetTitle(xAxisTitle);
  fMC_DirectGamma_Eta->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_MC_DirectGamma_Phi(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fMC_DirectGamma_Phi = new TH1F("MC_DirectGamma_Phi", "", nXBins, firstX, lastX);
  fMC_DirectGamma_Phi->GetXaxis()->SetTitle(xAxisTitle);
  fMC_DirectGamma_Phi->GetYaxis()->SetTitle(yAxisTitle);
}

//mapping


void AliGammaConversionHistograms::Initialize_MappingHistograms(Int_t nPhiIndex, Int_t nRIndex,Int_t nXBins, Double_t firstX, Double_t lastX, Int_t nYBins, Double_t firstY, Double_t lastY, TString xAxisTitle, TString yAxisTitle){
  
  for(Int_t phi =0; phi<=fNPhiIndex;phi++){
    AliConversionMappingVector tmpRowMC;
    AliConversionMappingVector tmpRowESD;
    for(Int_t r =0; r<fNRIndex;r++){
      //MC
      Char_t nameMC[100];
      sprintf(nameMC,"MC_EP_Mapping-Phi%02d-R%02d",phi,r);
      TH2F * tmpMC = new TH2F(nameMC,nameMC,nXBins,firstX,lastX,nYBins,firstY,lastY);
      tmpMC->SetXTitle(xAxisTitle);
      tmpMC->SetYTitle(yAxisTitle);
      tmpRowMC.push_back(tmpMC);
      //ESD
      Char_t nameESD[100];
      sprintf(nameESD,"ESD_EP_Mapping-Phi%02d-R%02d",phi,r);
      TH2F * tmpESD = new TH2F(nameESD,nameESD,nXBins,firstX,lastX,nYBins,firstY,lastY);
      tmpESD->SetXTitle(xAxisTitle);
      tmpESD->SetYTitle(yAxisTitle);
      tmpRowESD.push_back(tmpESD);
    }
    fMC_Mapping.push_back(tmpRowMC);
    fESD_Mapping.push_back(tmpRowESD);
  }

  for(Int_t phi =0; phi<=nPhiIndex;phi++){ 
    //MC
    Char_t nameMC[100];
    sprintf(nameMC,"MC_EP_Mapping-Phi%02d",phi);
    TH2F * tmpMC = new TH2F(nameMC,nameMC,nXBins,firstX,lastX,nYBins,firstY,lastY);
    tmpMC->SetXTitle(xAxisTitle);
    tmpMC->SetYTitle(yAxisTitle);
    fMC_Mapping_Phi.push_back(tmpMC);
    //ESD
    Char_t nameESD[100];
    sprintf(nameESD,"ESD_EP_Mapping-Phi%02d",phi);
    TH2F * tmpESD = new TH2F(nameESD,nameESD,nXBins,firstX,lastX,nYBins,firstY,lastY);
    tmpESD->SetXTitle(xAxisTitle);
    tmpESD->SetYTitle(yAxisTitle);
    fESD_Mapping_Phi.push_back(tmpESD);   
  }
  for(Int_t r =0; r<=nRIndex;r++){
    //MC
    Char_t nameMC[100];
    sprintf(nameMC,"MC_EP_Mapping-R%02d",r);
    TH2F * tmpMC = new TH2F(nameMC,nameMC,nXBins,firstX,lastX,nYBins,firstY,lastY);
    tmpMC->SetXTitle(xAxisTitle);
    tmpMC->SetYTitle(yAxisTitle);
    fMC_Mapping_R.push_back(tmpMC);
    //ESD
    Char_t nameESD[100];
    sprintf(nameESD,"ESD_EP_Mapping-R%02d",r);
    TH2F * tmpESD = new TH2F(nameESD,nameESD,nXBins,firstX,lastX,nYBins,firstY,lastY);
    tmpESD->SetXTitle(xAxisTitle);
    tmpESD->SetYTitle(yAxisTitle);
    fESD_Mapping_R.push_back(tmpESD);
  }
}

void AliGammaConversionHistograms::Initialize_MC_Match_Gamma_Eta(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fMC_Match_Gamma_Eta = new TH1F("MC_Match_Gamma_Eta", "", nXBins, firstX, lastX);
  fMC_Match_Gamma_Eta->GetXaxis()->SetTitle(xAxisTitle);
  fMC_Match_Gamma_Eta->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_MC_Match_Gamma_Phi(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fMC_Match_Gamma_Phi = new TH1F("MC_Match_Gamma_Phi", "", nXBins, firstX, lastX);
  fMC_Match_Gamma_Phi->GetXaxis()->SetTitle(xAxisTitle);
  fMC_Match_Gamma_Phi->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_MC_Match_Gamma_Pt(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fMC_Match_Gamma_Pt = new TH1F("MC_Match_Gamma_Pt", "", nXBins, firstX, lastX);
  fMC_Match_Gamma_Pt->GetXaxis()->SetTitle(xAxisTitle);
  fMC_Match_Gamma_Pt->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_MC_Match_Gamma_Energy(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fMC_Match_Gamma_Energy = new TH1F("MC_Match_Gamma_Energy", "", nXBins, firstX, lastX);
  fMC_Match_Gamma_Energy->GetXaxis()->SetTitle(xAxisTitle);
  fMC_Match_Gamma_Energy->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_MC_Match_Gamma_Mass(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fMC_Match_Gamma_Mass = new TH1F("MC_Match_Gamma_Mass", "", nXBins, firstX, lastX);
  fMC_Match_Gamma_Mass->GetXaxis()->SetTitle(xAxisTitle);
  fMC_Match_Gamma_Mass->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_MC_Match_Gamma_OpeningAngle(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fMC_Match_Gamma_OpeningAngle = new TH1F("MC_Match_Gamma_OpeningAngle", "", nXBins, firstX, lastX);
  fMC_Match_Gamma_OpeningAngle->GetXaxis()->SetTitle(xAxisTitle);
  fMC_Match_Gamma_OpeningAngle->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_MC_Match_Gamma_R(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fMC_Match_Gamma_R = new TH1F("MC_Match_Gamma_R", "", nXBins, firstX, lastX);
  fMC_Match_Gamma_R->GetXaxis()->SetTitle(xAxisTitle);
  fMC_Match_Gamma_R->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_MC_Match_Gamma_Z_R(Int_t nXBins, Double_t firstX, Double_t lastX, Int_t nYBins, Double_t firstY, Double_t lastY, TString xAxisTitle, TString yAxisTitle){
  fMC_Match_Gamma_Z_R = new TH2F("MC_Match_Gamma_Z_R", "", nXBins, firstX, lastX, nYBins, firstY, lastY);
  fMC_Match_Gamma_Z_R->GetXaxis()->SetTitle(xAxisTitle);
  fMC_Match_Gamma_Z_R->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_MC_Match_Gamma_X_Y(Int_t nXBins, Double_t firstX, Double_t lastX, Int_t nYBins, Double_t firstY, Double_t lastY, TString xAxisTitle, TString yAxisTitle){
  fMC_Match_Gamma_X_Y = new TH2F("MC_Match_Gamma_X_Y", "", nXBins, firstX, lastX, nYBins, firstY, lastY);
  fMC_Match_Gamma_X_Y->GetXaxis()->SetTitle(xAxisTitle);
  fMC_Match_Gamma_X_Y->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_MC_Pi0_Eta(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fMC_Pi0_Eta = new TH1F("MC_Pi0_Eta", "", nXBins, firstX, lastX);
  fMC_Pi0_Eta->GetXaxis()->SetTitle(xAxisTitle);
  fMC_Pi0_Eta->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_MC_Pi0_Phi(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fMC_Pi0_Phi = new TH1F("MC_Pi0_Phi", "", nXBins, firstX, lastX);
  fMC_Pi0_Phi->GetXaxis()->SetTitle(xAxisTitle);
  fMC_Pi0_Phi->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_MC_Pi0_Pt(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fMC_Pi0_Pt = new TH1F("MC_Pi0_Pt", "", nXBins, firstX, lastX);
  fMC_Pi0_Pt->GetXaxis()->SetTitle(xAxisTitle);
  fMC_Pi0_Pt->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_MC_Pi0_Energy(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fMC_Pi0_Energy = new TH1F("MC_Pi0_Energy", "", nXBins, firstX, lastX);
  fMC_Pi0_Energy->GetXaxis()->SetTitle(xAxisTitle);
  fMC_Pi0_Energy->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_MC_Pi0_Mass(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fMC_Pi0_Mass = new TH1F("MC_Pi0_Mass", "", nXBins, firstX, lastX);
  fMC_Pi0_Mass->GetXaxis()->SetTitle(xAxisTitle);
  fMC_Pi0_Mass->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_MC_Pi0_OpeningAngleGamma(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fMC_Pi0_OpeningAngleGamma = new TH1F("MC_Pi0_OpeningAngleGamma", "", nXBins, firstX, lastX);
  fMC_Pi0_OpeningAngleGamma->GetXaxis()->SetTitle(xAxisTitle);
  fMC_Pi0_OpeningAngleGamma->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_MC_Pi0_R(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fMC_Pi0_R = new TH1F("MC_Pi0_R", "", nXBins, firstX, lastX);
  fMC_Pi0_R->GetXaxis()->SetTitle(xAxisTitle);
  fMC_Pi0_R->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_MC_Pi0_Z_R(Int_t nXBins, Double_t firstX, Double_t lastX, Int_t nYBins, Double_t firstY, Double_t lastY, TString xAxisTitle, TString yAxisTitle){
  fMC_Pi0_Z_R = new TH2F("MC_Pi0_Z_R", "", nXBins, firstX, lastX, nYBins, firstY, lastY);
  fMC_Pi0_Z_R->GetXaxis()->SetTitle(xAxisTitle);
  fMC_Pi0_Z_R->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_MC_Pi0_X_Y(Int_t nXBins, Double_t firstX, Double_t lastX, Int_t nYBins, Double_t firstY, Double_t lastY, TString xAxisTitle, TString yAxisTitle){
  fMC_Pi0_X_Y = new TH2F("MC_Pi0_X_Y", "", nXBins, firstX, lastX, nYBins, firstY, lastY);
  fMC_Pi0_X_Y->GetXaxis()->SetTitle(xAxisTitle);
  fMC_Pi0_X_Y->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_MC_Pi0Secondaries_X_Y(Int_t nXBins, Double_t firstX, Double_t lastX, Int_t nYBins, Double_t firstY, Double_t lastY, TString xAxisTitle, TString yAxisTitle){
  fMC_Pi0Secondaries_X_Y = new TH2F("MC_Pi0Secondaries_X_Y", "", nXBins, firstX, lastX, nYBins, firstY, lastY);
  fMC_Pi0Secondaries_X_Y->GetXaxis()->SetTitle(xAxisTitle);
  fMC_Pi0Secondaries_X_Y->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_MC_Eta_Eta(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fMC_Eta_Eta = new TH1F("MC_Eta_Eta", "", nXBins, firstX, lastX);
  fMC_Eta_Eta->GetXaxis()->SetTitle(xAxisTitle);
  fMC_Eta_Eta->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_MC_Eta_Phi(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fMC_Eta_Phi = new TH1F("MC_Eta_Phi", "", nXBins, firstX, lastX);
  fMC_Eta_Phi->GetXaxis()->SetTitle(xAxisTitle);
  fMC_Eta_Phi->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_MC_Eta_Pt(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fMC_Eta_Pt = new TH1F("MC_Eta_Pt", "", nXBins, firstX, lastX);
  fMC_Eta_Pt->GetXaxis()->SetTitle(xAxisTitle);
  fMC_Eta_Pt->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_MC_Eta_Energy(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fMC_Eta_Energy = new TH1F("MC_Eta_Energy", "", nXBins, firstX, lastX);
  fMC_Eta_Energy->GetXaxis()->SetTitle(xAxisTitle);
  fMC_Eta_Energy->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_MC_Eta_Mass(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fMC_Eta_Mass = new TH1F("MC_Eta_Mass", "", nXBins, firstX, lastX);
  fMC_Eta_Mass->GetXaxis()->SetTitle(xAxisTitle);
  fMC_Eta_Mass->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_MC_Eta_OpeningAngleGamma(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fMC_Eta_OpeningAngleGamma = new TH1F("MC_Eta_OpeningAngleGamma", "", nXBins, firstX, lastX);
  fMC_Eta_OpeningAngleGamma->GetXaxis()->SetTitle(xAxisTitle);
  fMC_Eta_OpeningAngleGamma->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_MC_Eta_R(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fMC_Eta_R = new TH1F("MC_Eta_R", "", nXBins, firstX, lastX);
  fMC_Eta_R->GetXaxis()->SetTitle(xAxisTitle);
  fMC_Eta_R->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_MC_Eta_Z_R(Int_t nXBins, Double_t firstX, Double_t lastX, Int_t nYBins, Double_t firstY, Double_t lastY, TString xAxisTitle, TString yAxisTitle){
  fMC_Eta_Z_R = new TH2F("MC_Eta_Z_R", "", nXBins, firstX, lastX, nYBins, firstY, lastY);
  fMC_Eta_Z_R->GetXaxis()->SetTitle(xAxisTitle);
  fMC_Eta_Z_R->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_MC_Eta_X_Y(Int_t nXBins, Double_t firstX, Double_t lastX, Int_t nYBins, Double_t firstY, Double_t lastY, TString xAxisTitle, TString yAxisTitle){
  fMC_Eta_X_Y = new TH2F("MC_Eta_X_Y", "", nXBins, firstX, lastX, nYBins, firstY, lastY);
  fMC_Eta_X_Y->GetXaxis()->SetTitle(xAxisTitle);
  fMC_Eta_X_Y->GetYaxis()->SetTitle(yAxisTitle);
}   
// esd

void AliGammaConversionHistograms::Initialize_ESD_EP_R(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fESD_EP_R = new TH1F("ESD_EP_R", "", nXBins, firstX, lastX);
  fESD_EP_R->GetXaxis()->SetTitle(xAxisTitle);
  fESD_EP_R->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_EP_Z_R(Int_t nXBins, Double_t firstX, Double_t lastX, Int_t nYBins, Double_t firstY, Double_t lastY, TString xAxisTitle, TString yAxisTitle){
  fESD_EP_Z_R = new TH2F("ESD_EP_Z_R", "", nXBins, firstX, lastX, nYBins, firstY, lastY);
  fESD_EP_Z_R->GetXaxis()->SetTitle(xAxisTitle);
  fESD_EP_Z_R->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_EP_X_Y(Int_t nXBins, Double_t firstX, Double_t lastX, Int_t nYBins, Double_t firstY, Double_t lastY, TString xAxisTitle, TString yAxisTitle){
  fESD_EP_X_Y = new TH2F("ESD_EP_X_Y", "", nXBins, firstX, lastX, nYBins, firstY, lastY);
  fESD_EP_X_Y->GetXaxis()->SetTitle(xAxisTitle);
  fESD_EP_X_Y->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_EP_OpeningAngle(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fESD_EP_OpeningAngle = new TH1F("ESD_EP_OpeningAngle", "", nXBins, firstX, lastX);
  fESD_EP_OpeningAngle->GetXaxis()->SetTitle(xAxisTitle);
  fESD_EP_OpeningAngle->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_E_Energy(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fESD_E_Energy = new TH1F("ESD_E_Energy", "", nXBins, firstX, lastX);
  fESD_E_Energy->GetXaxis()->SetTitle(xAxisTitle);
  fESD_E_Energy->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_E_Pt(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fESD_E_Pt = new TH1F("ESD_E_Pt", "", nXBins, firstX, lastX);
  fESD_E_Pt->GetXaxis()->SetTitle(xAxisTitle);
  fESD_E_Pt->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_E_Eta(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fESD_E_Eta = new TH1F("ESD_E_Eta", "", nXBins, firstX, lastX);
  fESD_E_Eta->GetXaxis()->SetTitle(xAxisTitle);
  fESD_E_Eta->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_E_Phi(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fESD_E_Phi = new TH1F("ESD_E_Phi", "", nXBins, firstX, lastX);
  fESD_E_Phi->GetXaxis()->SetTitle(xAxisTitle);
  fESD_E_Phi->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_P_Energy(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fESD_P_Energy = new TH1F("ESD_P_Energy", "", nXBins, firstX, lastX);
  fESD_P_Energy->GetXaxis()->SetTitle(xAxisTitle);
  fESD_P_Energy->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_P_Pt(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fESD_P_Pt = new TH1F("ESD_P_Pt", "", nXBins, firstX, lastX);
  fESD_P_Pt->GetXaxis()->SetTitle(xAxisTitle);
  fESD_P_Pt->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_P_Eta(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fESD_P_Eta = new TH1F("ESD_P_Eta", "", nXBins, firstX, lastX);
  fESD_P_Eta->GetXaxis()->SetTitle(xAxisTitle);
  fESD_P_Eta->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_P_Phi(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fESD_P_Phi = new TH1F("ESD_P_Phi", "", nXBins, firstX, lastX);
  fESD_P_Phi->GetXaxis()->SetTitle(xAxisTitle);
  fESD_P_Phi->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_Gamma_Energy(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fESD_Gamma_Energy = new TH1F("ESD_Gamma_Energy", "", nXBins, firstX, lastX);
  fESD_Gamma_Energy->GetXaxis()->SetTitle(xAxisTitle);
  fESD_Gamma_Energy->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_Gamma_Pt(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fESD_Gamma_Pt = new TH1F("ESD_Gamma_Pt", "", nXBins, firstX, lastX);
  fESD_Gamma_Pt->GetXaxis()->SetTitle(xAxisTitle);
  fESD_Gamma_Pt->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_Gamma_Eta(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fESD_Gamma_Eta = new TH1F("ESD_Gamma_Eta", "", nXBins, firstX, lastX);
  fESD_Gamma_Eta->GetXaxis()->SetTitle(xAxisTitle);
  fESD_Gamma_Eta->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_Gamma_Phi(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fESD_Gamma_Phi = new TH1F("ESD_Gamma_Phi", "", nXBins, firstX, lastX);
  fESD_Gamma_Phi->GetXaxis()->SetTitle(xAxisTitle);
  fESD_Gamma_Phi->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_Match_Gamma_OpeningAngle(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fESD_Match_Gamma_OpeningAngle = new TH1F("ESD_Match_Gamma_OpeningAngle", "", nXBins, firstX, lastX);
  fESD_Match_Gamma_OpeningAngle->GetXaxis()->SetTitle(xAxisTitle);
  fESD_Match_Gamma_OpeningAngle->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_Match_Gamma_Energy(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fESD_Match_Gamma_Energy = new TH1F("ESD_Match_Gamma_Energy", "", nXBins, firstX, lastX);
  fESD_Match_Gamma_Energy->GetXaxis()->SetTitle(xAxisTitle);
  fESD_Match_Gamma_Energy->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_Match_Gamma_Pt(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fESD_Match_Gamma_Pt = new TH1F("ESD_Match_Gamma_Pt", "", nXBins, firstX, lastX);
  fESD_Match_Gamma_Pt->GetXaxis()->SetTitle(xAxisTitle);
  fESD_Match_Gamma_Pt->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_Match_Gamma_Eta(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fESD_Match_Gamma_Eta = new TH1F("ESD_Match_Gamma_Eta", "", nXBins, firstX, lastX);
  fESD_Match_Gamma_Eta->GetXaxis()->SetTitle(xAxisTitle);
  fESD_Match_Gamma_Eta->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_Match_Gamma_Phi(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fESD_Match_Gamma_Phi = new TH1F("ESD_Match_Gamma_Phi", "", nXBins, firstX, lastX);
  fESD_Match_Gamma_Phi->GetXaxis()->SetTitle(xAxisTitle);
  fESD_Match_Gamma_Phi->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_Match_Gamma_Mass(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fESD_Match_Gamma_Mass = new TH1F("ESD_Match_Gamma_Mass", "", nXBins, firstX, lastX);
  fESD_Match_Gamma_Mass->GetXaxis()->SetTitle(xAxisTitle);
  fESD_Match_Gamma_Mass->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_Match_Gamma_Width(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fESD_Match_Gamma_Width = new TH1F("ESD_Match_Gamma_Width", "", nXBins, firstX, lastX);
  fESD_Match_Gamma_Width->GetXaxis()->SetTitle(xAxisTitle);
  fESD_Match_Gamma_Width->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_Match_Gamma_Chi2(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fESD_Match_Gamma_Chi2 = new TH1F("ESD_Match_Gamma_Chi2", "", nXBins, firstX, lastX);
  fESD_Match_Gamma_Chi2->GetXaxis()->SetTitle(xAxisTitle);
  fESD_Match_Gamma_Chi2->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_Match_Gamma_NDF(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fESD_Match_Gamma_NDF = new TH1F("ESD_Match_Gamma_NDF", "", nXBins, firstX, lastX);
  fESD_Match_Gamma_NDF->GetXaxis()->SetTitle(xAxisTitle);
  fESD_Match_Gamma_NDF->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_Match_Gamma_R(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fESD_Match_Gamma_R = new TH1F("ESD_Match_Gamma_R", "", nXBins, firstX, lastX);
  fESD_Match_Gamma_R->GetXaxis()->SetTitle(xAxisTitle);
  fESD_Match_Gamma_R->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_Match_Gamma_Z_R(Int_t nXBins, Double_t firstX, Double_t lastX, Int_t nYBins, Double_t firstY, Double_t lastY, TString xAxisTitle, TString yAxisTitle){
  fESD_Match_Gamma_Z_R = new TH2F("ESD_Match_Gamma_Z_R", "", nXBins, firstX, lastX, nYBins, firstY, lastY);
  fESD_Match_Gamma_Z_R->GetXaxis()->SetTitle(xAxisTitle);
  fESD_Match_Gamma_Z_R->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_Match_Gamma_X_Y(Int_t nXBins, Double_t firstX, Double_t lastX, Int_t nYBins, Double_t firstY, Double_t lastY, TString xAxisTitle, TString yAxisTitle){
  fESD_Match_Gamma_X_Y = new TH2F("ESD_Match_Gamma_X_Y", "", nXBins, firstX, lastX, nYBins, firstY, lastY);
  fESD_Match_Gamma_X_Y->GetXaxis()->SetTitle(xAxisTitle);
  fESD_Match_Gamma_X_Y->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_Pi0_OpeningAngleGamma(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fESD_Pi0_OpeningAngleGamma = new TH1F("ESD_Pi0_OpeningAngleGamma", "", nXBins, firstX, lastX);
  fESD_Pi0_OpeningAngleGamma->GetXaxis()->SetTitle(xAxisTitle);
  fESD_Pi0_OpeningAngleGamma->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_Pi0_Energy(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fESD_Pi0_Energy = new TH1F("ESD_Pi0_Energy", "", nXBins, firstX, lastX);
  fESD_Pi0_Energy->GetXaxis()->SetTitle(xAxisTitle);
  fESD_Pi0_Energy->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_Pi0_Pt(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fESD_Pi0_Pt = new TH1F("ESD_Pi0_Pt", "", nXBins, firstX, lastX);
  fESD_Pi0_Pt->GetXaxis()->SetTitle(xAxisTitle);
  fESD_Pi0_Pt->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_Pi0_Eta(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fESD_Pi0_Eta = new TH1F("ESD_Pi0_Eta", "", nXBins, firstX, lastX);
  fESD_Pi0_Eta->GetXaxis()->SetTitle(xAxisTitle);
  fESD_Pi0_Eta->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_Pi0_Phi(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fESD_Pi0_Phi = new TH1F("ESD_Pi0_Phi", "", nXBins, firstX, lastX);
  fESD_Pi0_Phi->GetXaxis()->SetTitle(xAxisTitle);
  fESD_Pi0_Phi->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_Pi0_Mass(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fESD_Pi0_Mass = new TH1F("ESD_Pi0_Mass", "", nXBins, firstX, lastX);
  fESD_Pi0_Mass->GetXaxis()->SetTitle(xAxisTitle);
  fESD_Pi0_Mass->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_Pi0_R(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fESD_Pi0_R = new TH1F("ESD_Pi0_R", "", nXBins, firstX, lastX);
  fESD_Pi0_R->GetXaxis()->SetTitle(xAxisTitle);
  fESD_Pi0_R->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_Pi0_Z_R(Int_t nXBins, Double_t firstX, Double_t lastX, Int_t nYBins, Double_t firstY, Double_t lastY, TString xAxisTitle, TString yAxisTitle){
  fESD_Pi0_Z_R = new TH2F("ESD_Pi0_Z_R", "", nXBins, firstX, lastX, nYBins, firstY, lastY);
  fESD_Pi0_Z_R->GetXaxis()->SetTitle(xAxisTitle);
  fESD_Pi0_Z_R->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_Pi0_X_Y(Int_t nXBins, Double_t firstX, Double_t lastX, Int_t nYBins, Double_t firstY, Double_t lastY, TString xAxisTitle, TString yAxisTitle){
  fESD_Pi0_X_Y = new TH2F("ESD_Pi0_X_Y", "", nXBins, firstX, lastX, nYBins, firstY, lastY);
  fESD_Pi0_X_Y->GetXaxis()->SetTitle(xAxisTitle);
  fESD_Pi0_X_Y->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_Eta_OpeningAngleGamma(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fESD_Eta_OpeningAngleGamma = new TH1F("ESD_Eta_OpeningAngleGamma", "", nXBins, firstX, lastX);
  fESD_Eta_OpeningAngleGamma->GetXaxis()->SetTitle(xAxisTitle);
  fESD_Eta_OpeningAngleGamma->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_Eta_Energy(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fESD_Eta_Energy = new TH1F("ESD_Eta_Energy", "", nXBins, firstX, lastX);
  fESD_Eta_Energy->GetXaxis()->SetTitle(xAxisTitle);
  fESD_Eta_Energy->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_Eta_Pt(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fESD_Eta_Pt = new TH1F("ESD_Eta_Pt", "", nXBins, firstX, lastX);
  fESD_Eta_Pt->GetXaxis()->SetTitle(xAxisTitle);
  fESD_Eta_Pt->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_Eta_Eta(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fESD_Eta_Eta = new TH1F("ESD_Eta_Eta", "", nXBins, firstX, lastX);
  fESD_Eta_Eta->GetXaxis()->SetTitle(xAxisTitle);
  fESD_Eta_Eta->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_Eta_Phi(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fESD_Eta_Phi = new TH1F("ESD_Eta_Phi", "", nXBins, firstX, lastX);
  fESD_Eta_Phi->GetXaxis()->SetTitle(xAxisTitle);
  fESD_Eta_Phi->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_Eta_Mass(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fESD_Eta_Mass = new TH1F("ESD_Eta_Mass", "", nXBins, firstX, lastX);
  fESD_Eta_Mass->GetXaxis()->SetTitle(xAxisTitle);
  fESD_Eta_Mass->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_Eta_R(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fESD_Eta_R = new TH1F("ESD_Eta_R", "", nXBins, firstX, lastX);
  fESD_Eta_R->GetXaxis()->SetTitle(xAxisTitle);
  fESD_Eta_R->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_Eta_Z_R(Int_t nXBins, Double_t firstX, Double_t lastX, Int_t nYBins, Double_t firstY, Double_t lastY, TString xAxisTitle, TString yAxisTitle){
  fESD_Eta_Z_R = new TH2F("ESD_Eta_Z_R", "", nXBins, firstX, lastX, nYBins, firstY, lastY);
  fESD_Eta_Z_R->GetXaxis()->SetTitle(xAxisTitle);
  fESD_Eta_Z_R->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_Eta_X_Y(Int_t nXBins, Double_t firstX, Double_t lastX, Int_t nYBins, Double_t firstY, Double_t lastY, TString xAxisTitle, TString yAxisTitle){
  fESD_Eta_X_Y = new TH2F("ESD_Eta_X_Y", "", nXBins, firstX, lastX, nYBins, firstY, lastY);
  fESD_Eta_X_Y->GetXaxis()->SetTitle(xAxisTitle);
  fESD_Eta_X_Y->GetYaxis()->SetTitle(yAxisTitle);
}
void AliGammaConversionHistograms::Initialize_ESD_Background_OpeningAngleGamma(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fESD_Background_OpeningAngleGamma = new TH1F("ESD_Background_OpeningAngleGamma", "", nXBins, firstX, lastX);
  fESD_Background_OpeningAngleGamma->GetXaxis()->SetTitle(xAxisTitle);
  fESD_Background_OpeningAngleGamma->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_Background_Energy(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fESD_Background_Energy = new TH1F("ESD_Background_Energy", "", nXBins, firstX, lastX);
  fESD_Background_Energy->GetXaxis()->SetTitle(xAxisTitle);
  fESD_Background_Energy->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_Background_Pt(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fESD_Background_Pt = new TH1F("ESD_Background_Pt", "", nXBins, firstX, lastX);
  fESD_Background_Pt->GetXaxis()->SetTitle(xAxisTitle);
  fESD_Background_Pt->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_Background_Eta(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fESD_Background_Eta = new TH1F("ESD_Background_Eta", "", nXBins, firstX, lastX);
  fESD_Background_Eta->GetXaxis()->SetTitle(xAxisTitle);
  fESD_Background_Eta->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_Background_Phi(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fESD_Background_Phi = new TH1F("ESD_Background_Phi", "", nXBins, firstX, lastX);
  fESD_Background_Phi->GetXaxis()->SetTitle(xAxisTitle);
  fESD_Background_Phi->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_Background_Mass(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fESD_Background_Mass = new TH1F("ESD_Background_Mass", "", nXBins, firstX, lastX);
  fESD_Background_Mass->GetXaxis()->SetTitle(xAxisTitle);
  fESD_Background_Mass->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_Background_R(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fESD_Background_R = new TH1F("ESD_Background_R", "", nXBins, firstX, lastX);
  fESD_Background_R->GetXaxis()->SetTitle(xAxisTitle);
  fESD_Background_R->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_Background_Z_R(Int_t nXBins, Double_t firstX, Double_t lastX, Int_t nYBins, Double_t firstY, Double_t lastY, TString xAxisTitle, TString yAxisTitle){
  fESD_Background_Z_R = new TH2F("ESD_Background_Z_R", "", nXBins, firstX, lastX, nYBins, firstY, lastY);
  fESD_Background_Z_R->GetXaxis()->SetTitle(xAxisTitle);
  fESD_Background_Z_R->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_ESD_Background_X_Y(Int_t nXBins, Double_t firstX, Double_t lastX, Int_t nYBins, Double_t firstY, Double_t lastY, TString xAxisTitle, TString yAxisTitle){
  fESD_Background_X_Y = new TH2F("ESD_Background_X_Y", "", nXBins, firstX, lastX, nYBins, firstY, lastY);
  fESD_Background_X_Y->GetXaxis()->SetTitle(xAxisTitle);
  fESD_Background_X_Y->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_Resolution_dPt(Int_t nXBins, Double_t firstX, Double_t lastX, Int_t nYBins, Double_t firstY, Double_t lastY, TString xAxisTitle, TString yAxisTitle){
  fResolution_dPt = new TH2F("Resolution_dPt", "", nXBins, firstX, lastX, nYBins, firstY, lastY);
  fResolution_dPt->GetXaxis()->SetTitle(xAxisTitle);
  fResolution_dPt->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_Resolution_dR(Int_t nXBins, Double_t firstX, Double_t lastX, Int_t nYBins, Double_t firstY, Double_t lastY, TString xAxisTitle, TString yAxisTitle){
  fResolution_dR = new TH2F("Resolution_dR", "", nXBins, firstX, lastX, nYBins, firstY, lastY);
  fResolution_dR->GetXaxis()->SetTitle(xAxisTitle);
  fResolution_dR->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_Resolution_dZ(Int_t nXBins, Double_t firstX, Double_t lastX, Int_t nYBins, Double_t firstY, Double_t lastY, TString xAxisTitle, TString yAxisTitle){
  fResolution_dZ = new TH2F("Resolution_dZ", "", nXBins, firstX, lastX, nYBins, firstY, lastY);
  fResolution_dZ->GetXaxis()->SetTitle(xAxisTitle);
  fResolution_dZ->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_Resolution_dR_dPt(Int_t nXBins, Double_t firstX, Double_t lastX, Int_t nYBins, Double_t firstY, Double_t lastY, TString xAxisTitle, TString yAxisTitle){
  fResolution_dR_dPt = new TH2F("Resolution_dR_dPt", "", nXBins, firstX, lastX, nYBins, firstY, lastY);
  fResolution_dR_dPt->GetXaxis()->SetTitle(xAxisTitle);
  fResolution_dR_dPt->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_Resolution_MC_Pt(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fResolution_MC_Pt = new TH1F("Resolution_MC_Pt", "", nXBins, firstX, lastX);
  fResolution_MC_Pt->GetXaxis()->SetTitle(xAxisTitle);
  fResolution_MC_Pt->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_Resolution_MC_R(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fResolution_MC_R = new TH1F("Resolution_MC_R", "", nXBins, firstX, lastX);
  fResolution_MC_R->GetXaxis()->SetTitle(xAxisTitle);
  fResolution_MC_R->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_Resolution_MC_Z(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fResolution_MC_Z = new TH1F("Resolution_MC_Z", "", nXBins, firstX, lastX);
  fResolution_MC_Z->GetXaxis()->SetTitle(xAxisTitle);
  fResolution_MC_Z->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_Resolution_ESD_Pt(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fResolution_ESD_Pt = new TH1F("Resolution_ESD_Pt", "", nXBins, firstX, lastX);
  fResolution_ESD_Pt->GetXaxis()->SetTitle(xAxisTitle);
  fResolution_ESD_Pt->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_Resolution_ESD_R(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fResolution_ESD_R = new TH1F("Resolution_ESD_R", "", nXBins, firstX, lastX);
  fResolution_ESD_R->GetXaxis()->SetTitle(xAxisTitle);
  fResolution_ESD_R->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_Resolution_ESD_Z(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fResolution_ESD_Z = new TH1F("Resolution_ESD_Z", "", nXBins, firstX, lastX);
  fResolution_ESD_Z->GetXaxis()->SetTitle(xAxisTitle);
  fResolution_ESD_Z->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_NumberOfV0s(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fNumberOfV0s = new TH1F("NumberOfV0s", "", nXBins, firstX, lastX);
  fNumberOfV0s->GetXaxis()->SetTitle(xAxisTitle);
  fNumberOfV0s->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_NumberOfSurvivingV0s(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fNumberOfSurvivingV0s = new TH1F("NumberOfSurvivingV0s", "", nXBins, firstX, lastX);
  fNumberOfSurvivingV0s->GetXaxis()->SetTitle(xAxisTitle);
  fNumberOfSurvivingV0s->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_V0MassDebugCut1(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fV0MassDebugCut1 = new TH1F("V0MassDebugCut1", "", nXBins, firstX, lastX);
  fV0MassDebugCut1->GetXaxis()->SetTitle(xAxisTitle);
  fV0MassDebugCut1->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_V0MassDebugCut2(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fV0MassDebugCut2 = new TH1F("V0MassDebugCut2", "", nXBins, firstX, lastX);
  fV0MassDebugCut2->GetXaxis()->SetTitle(xAxisTitle);
  fV0MassDebugCut2->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_V0MassDebugCut3(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fV0MassDebugCut3 = new TH1F("V0MassDebugCut3", "", nXBins, firstX, lastX);
  fV0MassDebugCut3->GetXaxis()->SetTitle(xAxisTitle);
  fV0MassDebugCut3->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_V0MassDebugCut4(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fV0MassDebugCut4 = new TH1F("V0MassDebugCut4", "", nXBins, firstX, lastX);
  fV0MassDebugCut4->GetXaxis()->SetTitle(xAxisTitle);
  fV0MassDebugCut4->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_V0MassDebugCut5(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fV0MassDebugCut5 = new TH1F("V0MassDebugCut5", "", nXBins, firstX, lastX);
  fV0MassDebugCut5->GetXaxis()->SetTitle(xAxisTitle);
  fV0MassDebugCut5->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_V0MassDebugCut6(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fV0MassDebugCut6 = new TH1F("V0MassDebugCut6", "", nXBins, firstX, lastX);
  fV0MassDebugCut6->GetXaxis()->SetTitle(xAxisTitle);
  fV0MassDebugCut6->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_V0MassDebugCut7(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fV0MassDebugCut7 = new TH1F("V0MassDebugCut7", "", nXBins, firstX, lastX);
  fV0MassDebugCut7->GetXaxis()->SetTitle(xAxisTitle);
  fV0MassDebugCut7->GetYaxis()->SetTitle(yAxisTitle);
}

void AliGammaConversionHistograms::Initialize_V0MassDebugCut8(Int_t nXBins, Double_t firstX, Double_t lastX, TString xAxisTitle, TString yAxisTitle){
  fV0MassDebugCut8 = new TH1F("V0MassDebugCut8", "", nXBins, firstX, lastX);
  fV0MassDebugCut8->GetXaxis()->SetTitle(xAxisTitle);
  fV0MassDebugCut8->GetYaxis()->SetTitle(yAxisTitle);
}
