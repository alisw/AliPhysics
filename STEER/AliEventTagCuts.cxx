/**************************************************************************
 * Author: Panos Christakoglou.                                           *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//-----------------------------------------------------------------
//           AliEventTagCuts class
//   This is the class to deal with the event tag level cuts
//   Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-----------------------------------------------------------------

class AliLog;
class AliESD;

#include "AliEventTag.h"
#include "AliEventTagCuts.h"

ClassImp(AliEventTagCuts)


//----------------------------------------//
AliEventTagCuts::AliEventTagCuts()
{
  //Default constructor which calls the Reset method.
  Reset();
}

//----------------------------------------//
AliEventTagCuts::~AliEventTagCuts()
{  
  //Defaut destructor.
}

//----------------------------------------//
void AliEventTagCuts::Reset()
{
  //Sets dummy values to every private member.
  fVxFlag = kFALSE;
  fVyFlag = kFALSE;
  fVzFlag = kFALSE;
  fParticipantsFlag = kFALSE;
  fImpactParamFlag = kFALSE;
  fPVFlag = kFALSE;
  fZDCNeutronEnergyFlag = kFALSE;
  fZDCProtonEnergyFlag = kFALSE;
  fZDCEMEnergyFlag = kFALSE;
  fT0VertexZFlag = kFALSE;
  fMultFlag = kFALSE;
  fMultPosFlag = kFALSE;
  fMultNegFlag = kFALSE;
  fMultNeutrFlag = kFALSE;
  fV0sFlag = kFALSE;
  fCascadesFlag = kFALSE;
  fkinksFlag = kFALSE;
  fMaxJetEnergyFlag = kFALSE;
  fNHardPhotonsCandidatesFlag = kFALSE;
  fMaxNeutralFlag = kFALSE;
  fChargedAbove1GeVFlag = kFALSE;
  fChargedAbove3GeVFlag = kFALSE;
  fChargedAbove10GeVFlag = kFALSE;
  fMuonsAbove1GeVFlag = kFALSE;
  fMuonsAbove3GeVFlag = kFALSE;
  fMuonsAbove10GeVFlag = kFALSE;
  fElectronsAbove1GeVFlag = kFALSE;
  fElectronsAbove3GeVFlag = kFALSE;
  fElectronsAbove10GeVFlag = kFALSE;
  fElectronsFlag = kFALSE;
  fMuonsFlag = kFALSE;
  fPionsFlag = kFALSE;
  fKaonsFlag = kFALSE;
  fProtonsFlag = kFALSE;
  fLambdasFlag = kFALSE;
  fPhotonFlag = kFALSE;
  fPi0sFlag = kFALSE;
  fNeutronsFlag = kFALSE;
  fKaon0sFlag = kFALSE;
  fTotalPFlag = kFALSE;
  fMeanPtFlag = kFALSE;
  fMaxPtFlag = kFALSE;
  fTotalNeutralPFlag = kFALSE;
  fMeanNeutralPtFlag = kFALSE;
  fMaxNeutralPtFlag = kFALSE;
  fEventPlaneAngleFlag = kFALSE;
  fHBTRadiiFlag = kFALSE;
  
  fVxMin = -1000.0;
  fVxMax = 1000.0; 
  fVyMin = -1000.0;
  fVyMax = 1000.0;  
  fVzMin = -1000.0;
  fVzMax = 1000.0;

  fMultMin = 0;
  fMultMax = 100000;  

  fParticipantsMin = -1;
  fParticipantMax = 10000;
  fImpactParamMin = -1.0;
  fImpactParamMax = 1000.0;
  fPrimaryVertexFlag = 1;
 
  fZDCNeutronEnergyMin = -1.0;
  fZDCNeutronEnergyMax = 100000.0;
  fZDCProtonEnergyMin = -1.0;
  fZDCProtonEnergyMax = 100000.0;
  fZDCEMEnergyMin = -1.0;
  fZDCEMEnergyMax = 100000.0;
  fT0VertexZMin = -10000.0;
  fT0VertexZMax = 10000.0;
  
  fMultPosMin = -1;
  fMultPosMax = 100000;
  fMultNegMin = -1;
  fMultNegMax = 100000;
  fMultNeutrMin = -1;
  fMultNeutrMax = 100000;
  fV0sMin = -1;
  fV0sMax = 1000000;
  fCascadesMin = -1;
  fCascadesMax = 100000;
  fkinksMin = -1;
  fkinksMax = 1000000;

  fMaxJetEnergy = -1.0; 

  fNHardPhotonsCandidatesMin = -1;
  fNHardPhotonsCandidatesMax = 100000;
  fMaxNeutralEnergy = -1.0; 
  
  fChargedAbove1GeVMin = -1;
  fChargedAbove1GeVMax = 100000;
  fChargedAbove3GeVMin = -1;
  fChargedAbove3GeVMax = 100000;
  fChargedAbove10GeVMin = -1;
  fChargedAbove10GeVMax = 100000;
  fMuonsAbove1GeVMin = -1;
  fMuonsAbove1GeVMax = 100000;
  fMuonsAbove3GeVMin = -1;
  fMuonsAbove3GeVMax = 100000;
  fMuonsAbove10GeVMin = -1;
  fMuonsAbove10GeVMax = 100000; 
  fElectronsAbove1GeVMin = -1;
  fElectronsAbove1GeVMax = 100000;
  fElectronsAbove3GeVMin = -1;
  fElectronsAbove3GeVMax = 100000;
  fElectronsAbove10GeVMin = -1;
  fElectronsAbove10GeVMax = 100000;

  fElectronsMin = -1;
  fElectronsMax = 100000;
  fMuonsMin = -1;
  fMuonsMax = 100000;
  fPionsMin = -1;
  fPionsMax = 100000;
  fKaonsMin = -1;
  fKaonsMax = 100000;
  fProtonsMin = -1;
  fProtonsMax = 100000;
  fLambdasMin = -1;
  fLambdasMax = 100000;
  fPhotonsMin = -1;
  fPhotonsMax = 100000;
  fPi0sMin = -1;
  fPi0sMax = 100000; 
  fNeutronsMin = -1;
  fNeutronsMax = 100000; 
  fKaon0sMin = -1;
  fKaon0sMax = 100000; 

  fTotalPMin = -1.0;
  fTotalPMax = 1000000.0;
  fMeanPtMin = -1.0;
  fMeanPtMax = 100000.0;
  fMaxPt = -1.0; 
  fTotalNeutralPMin = -1.0;
  fTotalNeutralPMax = 1000000.0;  
  fMeanNeutralPtMin = -1.0;
  fMeanNeutralPtMax = 1000000.0; 
  fMaxNeutralPt = -1.0; 
  fEventPlaneAngleMin = -10000000.0;
  fEventPlaneAngleMax = 10000000.0; 
  fHBTRadiiMin = -1.0;
  fHBTRadiiMax = 100000.0; 
}

//----------------------------------------//
void AliEventTagCuts::SetPrimaryVertexXRange(Float_t r1, Float_t r2)
{
  //Sets the primary vertex x range 
  //and the corresponding flag to kTRUE if the cut is used.
  fVxMin = r1;
  fVxMax = r2; 
  fVxFlag = kTRUE;
}

//----------------------------------------//
void AliEventTagCuts::SetPrimaryVertexYRange(Float_t r1, Float_t r2)
{
  //Sets the primary vertex y range 
  //and the corresponding flag to kTRUE if the cut is used.
  fVyMin = r1;
  fVyMax = r2; 
  fVyFlag = kTRUE;
}

//----------------------------------------//
void AliEventTagCuts::SetPrimaryVertexZRange(Float_t r1, Float_t r2)
{
  //Sets the primary vertex z range 
  //and the corresponding flag to kTRUE if the cut is used.
  fVzMin = r1;
  fVzMax = r2; 
  fVzFlag = kTRUE;
}

//----------------------------------------//
void AliEventTagCuts::SetMultiplicityRange(Int_t n1, Int_t n2)
{
  //Sets the primary multiplicity range 
  //and the corresponding flag to kTRUE if the cut is used.
  fMultMin = n1;
  fMultMax = n2;
  fMultFlag = kTRUE;
}

//----------------------------------------//
void AliEventTagCuts::SetParticipantsRange(Int_t i1, Int_t i2)
{
  //Sets the number of participants range 
  //and the corresponding flag to kTRUE if the cut is used.
  fParticipantsMin = i1;
  fParticipantMax = i2;
  fParticipantsFlag = kTRUE;
}

//----------------------------------------//
void AliEventTagCuts::SetImpactParamRange(Float_t r1, Float_t r2)
{
  //Sets the impact parameter range 
  //and the corresponding flag to kTRUE if the cut is used.
  fImpactParamMin = r1;
  fImpactParamMax = r2;
  fImpactParamFlag = kTRUE;
}
 

//----------------------------------------//
void AliEventTagCuts::SetPrimaryVertexFlag(Int_t i)
{
  //Sets the primary vertex flag cut 
  //and the corresponding flag to kTRUE if the cut is used.
  fPrimaryVertexFlag = i;
  fPVFlag = kTRUE;
}

//----------------------------------------//
void AliEventTagCuts::SetZDCNeutrRange(Float_t r1, Float_t r2)
{
  //Sets the ZDC's neutron energy range 
  //and the corresponding flag to kTRUE if the cut is used.
  fZDCNeutronEnergyMin = r1;
  fZDCNeutronEnergyMax = r2;
  fZDCNeutronEnergyFlag = kTRUE;
}
//----------------------------------------//
void AliEventTagCuts::SetZDCProtRange(Float_t r1, Float_t r2)
{
  //Sets the ZDC's proton energy range 
  //and the corresponding flag to kTRUE if the cut is used.
  fZDCProtonEnergyMin = r1;
  fZDCProtonEnergyMax = r2;
  fZDCProtonEnergyFlag = kTRUE;
}
//----------------------------------------//
void AliEventTagCuts::SetZDCEMRange(Float_t r1, Float_t r2)
{
  //Sets the ZDC's e/m energy range 
  //and the corresponding flag to kTRUE if the cut is used.
  fZDCEMEnergyMin = r1;
  fZDCEMEnergyMax = r2;
  fZDCEMEnergyFlag = kTRUE;
}

//----------------------------------------//
void AliEventTagCuts::SetT0VertexZRange(Float_t r1, Float_t r2)
{
  //Sets the T0's Vz range 
  //and the corresponding flag to kTRUE if the cut is used.
  fT0VertexZMin = r1;
  fT0VertexZMax = r2;
  fT0VertexZFlag = kTRUE;
}

//----------------------------------------//
void AliEventTagCuts::SetPosMultiplicityRange(Int_t n1, Int_t n2)
{
  //Sets the positive multiplicity range 
  //and the corresponding flag to kTRUE if the cut is used.
  fMultPosMin = n1;
  fMultPosMax = n2;
  fMultPosFlag = kTRUE;
}


//----------------------------------------//
void AliEventTagCuts::SetNegMultiplicityRange(Int_t n1, Int_t n2)
{
  //Sets the negative multiplicity range 
  //and the corresponding flag to kTRUE if the cut is used.
  fMultNegMin = n1;
  fMultNegMax = n2;
  fMultNegFlag = kTRUE;
}


//----------------------------------------//
void AliEventTagCuts::SetNeutrMultiplicityRange(Int_t n1, Int_t n2)
{
  //Sets the neutral particle multiplicity range 
  //and the corresponding flag to kTRUE if the cut is used.
  fMultNeutrMin = n1;
  fMultNeutrMax = n2;
  fMultNeutrFlag = kTRUE;
}

//----------------------------------------//
void AliEventTagCuts::SetV0sRange(Int_t n1, Int_t n2)
{
  //Sets the v0s multiplicity range 
  //and the corresponding flag to kTRUE if the cut is used.
  fV0sMin = n1;
  fV0sMax = n2;
  fV0sFlag = kTRUE;
}

//----------------------------------------//
void AliEventTagCuts::SetCascadesRange(Int_t n1, Int_t n2)
{
  //Sets the cascades multiplicity range 
  //and the corresponding flag to kTRUE if the cut is used.
  fCascadesMin = n1;
  fCascadesMax = n2;
  fCascadesFlag = kTRUE;
}

//----------------------------------------//
void AliEventTagCuts::SetKinksRange(Int_t n1, Int_t n2)
{
  //Sets the kinks multipliicity range 
  //and the corresponding flag to kTRUE if the cut is used.
  fkinksMin = n1;
  fkinksMax = n2;
  fkinksFlag = kTRUE;
}

//----------------------------------------//
void AliEventTagCuts::SetMaxJetEnergy(Float_t r1)
{
  //Sets the lower limit of the maximum jet energy
  //and the corresponding flag to kTRUE if the cut is used.
  fMaxJetEnergy = r1; 
  fMaxJetEnergyFlag = kTRUE;
}
//----------------------------------------//
void AliEventTagCuts::SetMaxNeutralEnergy(Float_t r1)
{
  //Sets the lower limit of the maximum neutral jet energy
  //and the corresponding flag to kTRUE if the cut is used.
  fMaxNeutralEnergy = r1; 
  fMaxNeutralFlag = kTRUE;
}
//----------------------------------------//
void AliEventTagCuts::SetHardPhotonsRange(Int_t i1, Int_t i2)
{
  //Sets the hard photons multiplicity range
  //and the corresponding flag to kTRUE if the cut is used.
  fNHardPhotonsCandidatesMin = i1;
  fNHardPhotonsCandidatesMax = i2;
  fNHardPhotonsCandidatesFlag = kTRUE;
} 

//----------------------------------------//
void AliEventTagCuts::SetNChargedAbove1GeVRange(Int_t i1, Int_t i2)
{
  //Sets the number of charged above 1GeV range
  //and the corresponding flag to kTRUE if the cut is used.
  fChargedAbove1GeVMin = i1;
  fChargedAbove1GeVMax = i2;
  fChargedAbove1GeVFlag = kTRUE;
}

//----------------------------------------//
 void AliEventTagCuts::SetNChargedAbove3GeVRange(Int_t i1, Int_t i2)
{
  //Sets the number of charged above 3GeV range
  //and the corresponding flag to kTRUE if the cut is used.
  fChargedAbove3GeVMin = i1;
  fChargedAbove3GeVMax = i2;
  fChargedAbove3GeVFlag = kTRUE;
}


//----------------------------------------//
void AliEventTagCuts::SetNChargedAbove10GeVRange(Int_t i1, Int_t i2)
{
  //Sets the number of charged above 10GeV range
  //and the corresponding flag to kTRUE if the cut is used.
  fChargedAbove10GeVMin = i1;
  fChargedAbove10GeVMax = i2;
  fChargedAbove10GeVFlag = kTRUE;
}


//----------------------------------------//
void AliEventTagCuts::SetNMuonsAbove1GeVRange(Int_t i1, Int_t i2)
{
  //Sets the number of muons above 1GeV range
  //and the corresponding flag to kTRUE if the cut is used.
  fMuonsAbove1GeVMin = i1;
  fMuonsAbove1GeVMax = i2;
  fMuonsAbove1GeVFlag = kTRUE;
}


//----------------------------------------//
void AliEventTagCuts::SetNMuonsAbove3GeVRange(Int_t i1, Int_t i2)
{
  //Sets the number of muons above 3GeV range
  //and the corresponding flag to kTRUE if the cut is used.
  fMuonsAbove3GeVMin = i1;
  fMuonsAbove3GeVMax = i2;
  fMuonsAbove3GeVFlag = kTRUE;
} 

//----------------------------------------//
void AliEventTagCuts::SetNMuonsAbove10GeVRange(Int_t i1, Int_t i2)
{
  //Sets the number of muons above 10GeV range
  //and the corresponding flag to kTRUE if the cut is used.
  fMuonsAbove10GeVMin = i1;
  fMuonsAbove10GeVMax = i2; 
  fMuonsAbove10GeVFlag = kTRUE;
}


//----------------------------------------//
void AliEventTagCuts::SetNElectronsAbove1GeVRange(Int_t i1, Int_t i2)
{
  //Sets the number of electrons above 1GeV range
  //and the corresponding flag to kTRUE if the cut is used.
  fElectronsAbove1GeVMin = i1;
  fElectronsAbove1GeVMax = i2;
  fElectronsAbove1GeVFlag = kTRUE;
}

//----------------------------------------//
void AliEventTagCuts::SetNElectronsAbove3GeVRange(Int_t i1, Int_t i2)
{
  //Sets the number of electrons above 3GeV range
  //and the corresponding flag to kTRUE if the cut is used.
  fElectronsAbove3GeVMin = i1;
  fElectronsAbove3GeVMax = i2;
  fElectronsAbove3GeVFlag = kTRUE;
}

//----------------------------------------//
void AliEventTagCuts::SetNElectronsAbove10GeVRange(Int_t i1, Int_t i2)
{  
  //Sets the number of electrons above 10GeV range
  //and the corresponding flag to kTRUE if the cut is used.
  fElectronsAbove10GeVMin = i1;
  fElectronsAbove10GeVMax = i2;
  fElectronsAbove10GeVFlag = kTRUE;
}
//----------------------------------------//
void AliEventTagCuts::SetNElectronRange(Int_t n1, Int_t n2)
{
  //Sets the electron multiplicity range
  //and the corresponding flag to kTRUE if the cut is used.
  fElectronsMin = n1;
  fElectronsMax = n2;
  fElectronsFlag = kTRUE;
}
//----------------------------------------//
void AliEventTagCuts::SetNMuonRange(Int_t n1, Int_t n2)
{
  //Sets the muon multiplicity range
  //and the corresponding flag to kTRUE if the cut is used.
  fMuonsMin = n1;
  fMuonsMax = n2;
  fMuonsFlag = kTRUE;
}

//----------------------------------------//
void AliEventTagCuts::SetNPionRange(Int_t n1, Int_t n2)
{
  //Sets the pion multiplicity range
  //and the corresponding flag to kTRUE if the cut is used.
  fPionsMin = n1;
  fPionsMax = n2;
  fPionsFlag = kTRUE;
} 

//----------------------------------------//
void AliEventTagCuts::SetNKaonRange(Int_t n1, Int_t n2)
{
  //Sets the kaon multiplicity range
  //and the corresponding flag to kTRUE if the cut is used.
  fKaonsMin = n1;
  fKaonsMax = n2;
  fKaonsFlag = kTRUE;
}

//----------------------------------------//
void AliEventTagCuts::SetNProtonRange(Int_t n1, Int_t n2)
{
  //Sets the proton multiplicity range
  //and the corresponding flag to kTRUE if the cut is used.
  fProtonsMin = n1;
  fProtonsMax = n2;
  fProtonsFlag = kTRUE;
} 

//----------------------------------------//
void AliEventTagCuts::SetNLambdaRange(Int_t n1, Int_t n2)
{
  //Sets the lambda multiplicity range
  //and the corresponding flag to kTRUE if the cut is used.
  fLambdasMin = n1;
  fLambdasMax = n2;
  fLambdasFlag = kTRUE;
} 
//----------------------------------------//
void AliEventTagCuts::SetNPhotonRange(Int_t n1, Int_t n2)
{
  //Sets the photon multiplicity range
  //and the corresponding flag to kTRUE if the cut is used.
  fPhotonsMin = n1;
  fPhotonsMax = n2;
  fPhotonFlag = kTRUE;
} 
//----------------------------------------//
void AliEventTagCuts::SetNPi0Range(Int_t n1, Int_t n2)
{
  //Sets the pi0 multiplicity range
  //and the corresponding flag to kTRUE if the cut is used.
  fPi0sMin = n1;
  fPi0sMax = n2; 
  fPi0sFlag = kTRUE;
}  

//----------------------------------------//
void AliEventTagCuts::SetNNeutronRange(Int_t n1, Int_t n2)
{
  //Sets the neutron multiplicity range
  //and the corresponding flag to kTRUE if the cut is used.
  fNeutronsMin = n1;
  fNeutronsMax = n2; 
  fNeutronsFlag = kTRUE;
}

//----------------------------------------//
void AliEventTagCuts::SetNKaon0Range(Int_t n1, Int_t n2)
{  
  //Sets the K0s multiplicity range
  //and the corresponding flag to kTRUE if the cut is used.
  fKaon0sMin = n1;
  fKaon0sMax = n2; 
  fKaon0sFlag = kTRUE;
}

//----------------------------------------//
void AliEventTagCuts::SetTotalPRange(Float_t r1, Float_t r2)
{
  //Sets the total momentum range
  //and the corresponding flag to kTRUE if the cut is used.
  fTotalPMin = r1;
  fTotalPMax = r2;
  fTotalPFlag = kTRUE;
}

//----------------------------------------//
void AliEventTagCuts::SetMeanPtRange(Float_t r1, Float_t r2)
{
  //Sets the mean Pt range
  //and the corresponding flag to kTRUE if the cut is used.
  fMeanPtMin = r1;
  fMeanPtMax = r2;
  fMeanPtFlag = kTRUE;
}  

//----------------------------------------//
void AliEventTagCuts::SetMaxPt(Float_t r1)
{
  //Sets the lower limit of the max Pt value
  //and the corresponding flag to kTRUE if the cut is used.
  fMaxPt = r1; 
  fMaxPtFlag = kTRUE;
}

//----------------------------------------//
void AliEventTagCuts::SetTotalNeutralPRange(Float_t r1, Float_t r2)
{  
  //Sets the total momentum of neutral particles range
  //and the corresponding flag to kTRUE if the cut is used.
  fTotalNeutralPMin =r1 ;
  fTotalNeutralPMax = r2;  
  fTotalNeutralPFlag = kTRUE;
}
//----------------------------------------//
void AliEventTagCuts::SetMeanNeutralPtPRange(Float_t r1, Float_t r2)
{  
  //Sets the mean Pt of neutral particles range
  //and the corresponding flag to kTRUE if the cut is used.
  fMeanNeutralPtMin = r1;
  fMeanNeutralPtMax = r2; 
  fMeanNeutralPtFlag = kTRUE;
} 
//----------------------------------------//
void AliEventTagCuts::SetMaxNeutralPt(Float_t r1)
{  
  //Sets the lower limit of the maximum Pt of neutral particles
  //and the corresponding flag to kTRUE if the cut is used.
  fMaxNeutralPt = r1; 
  fMaxNeutralPtFlag = kTRUE;
}

//----------------------------------------//
void AliEventTagCuts::SetEvPlaneAngleRange(Float_t r1, Float_t r2)
{
  //Sets the event plane range
  //and the corresponding flag to kTRUE if the cut is used.
  fEventPlaneAngleMin = r1;
  fEventPlaneAngleMax = r2; 
  fEventPlaneAngleFlag = kTRUE;
}

//----------------------------------------//
void AliEventTagCuts::SetHBTRadiiRange(Float_t r1, Float_t r2)
{
  //Sets the HBT radii range
  //and the corresponding flag to kTRUE if the cut is used.
  fHBTRadiiMin = r1;
  fHBTRadiiMax = r2; 
  fHBTRadiiFlag = kTRUE;
}

//----------------------------------------//
Bool_t AliEventTagCuts::IsAccepted(AliEventTag *EvTag) const
{
  //Returns true if the event is accepted otherwise false.
  if(fMultFlag)
    if((EvTag->GetNumOfTracks() < fMultMin) || (EvTag->GetNumOfTracks() > fMultMax))
      return kFALSE; 
  
  if(fVzFlag)
    if((EvTag->GetVertexZ() < fVzMin) || (EvTag->GetVertexZ() > fVzMax))
      return kFALSE;
  
  if(fVyFlag)
    if((EvTag->GetVertexY() < fVyMin) || (EvTag->GetVertexY() > fVyMax))
      return kFALSE;
  
  if(fVxFlag)
    if((EvTag->GetVertexX() < fVxMin) || (EvTag->GetVertexX() > fVxMax))
      return kFALSE;
  
  if(fParticipantsFlag)
    if((EvTag->GetNumOfParticipants() < fParticipantsMin) || (EvTag->GetNumOfParticipants() > fParticipantMax))
      return kFALSE; 
  
  if(fImpactParamFlag)
    if((EvTag->GetImpactParameter() < fImpactParamMin) || (EvTag->GetImpactParameter() > fImpactParamMax))
      return kFALSE; 
  
  if(fPVFlag)
    if((EvTag->GetVertexFlag() != fPrimaryVertexFlag))
      return kFALSE; 
  
  if(fZDCNeutronEnergyFlag)
    if((EvTag->GetZDCNeutronEnergy() < fZDCNeutronEnergyMin) || (EvTag->GetZDCNeutronEnergy() > fZDCNeutronEnergyMax))
      return kFALSE; 
  
  if(fZDCProtonEnergyFlag)
    if((EvTag->GetZDCProtonEnergy() < fZDCProtonEnergyMin) || (EvTag->GetZDCProtonEnergy() > fZDCProtonEnergyMax))
      return kFALSE; 
  
  if(fZDCEMEnergyFlag)
    if((EvTag->GetZDCEMEnergy() < fZDCEMEnergyMin) || (EvTag->GetZDCEMEnergy() > fZDCEMEnergyMax))
      return kFALSE; 
  
  if(fT0VertexZFlag)
    if((EvTag->GetT0VertexZ() < fT0VertexZMin) || (EvTag->GetT0VertexZ() > fT0VertexZMax))
      return kFALSE; 
  
  if(fMultPosFlag)
    if((EvTag->GetNumOfPosTracks() < fMultPosMin) || (EvTag->GetNumOfPosTracks() > fMultPosMax))
      return kFALSE; 
  
  if(fMultNegFlag)
    if((EvTag->GetNumOfNegTracks() < fMultNegMin) || (EvTag->GetNumOfNegTracks() > fMultNegMax))
      return kFALSE; 
  
  if(fMultNeutrFlag)
    if((EvTag->GetNumOfNeutrTracks() < fMultNeutrMin) || (EvTag->GetNumOfNeutrTracks() > fMultNeutrMax))
      return kFALSE; 
  
  if(fV0sFlag)
    if((EvTag->GetNumOfV0s() < fV0sMin) || (EvTag->GetNumOfV0s() > fV0sMax))
      return kFALSE; 
  
  if(fCascadesFlag)
    if((EvTag->GetNumOfCascades() < fCascadesMin) || (EvTag->GetNumOfCascades() > fCascadesMax))
      return kFALSE; 
  
  if(fkinksFlag)
    if((EvTag->GetNumOfKinks() < fkinksMin) || (EvTag->GetNumOfKinks() > fkinksMax))
      return kFALSE; 
  
  if(fMaxJetEnergyFlag)
    if((EvTag->GetMaxJetEnergy() < fMaxJetEnergy))
      return kFALSE; 
  
  if(fNHardPhotonsCandidatesFlag)
    if((EvTag->GetNumOfHardPhotonsCandidates() < fNHardPhotonsCandidatesMin) || (EvTag->GetNumOfHardPhotonsCandidates() > fNHardPhotonsCandidatesMax))
      return kFALSE; 
  
  if(fMaxNeutralFlag)
    if((EvTag->GetMaxNeutralEnergy() < fMaxNeutralEnergy))
      return kFALSE; 
  
  if(fChargedAbove1GeVFlag)
    if((EvTag->GetNumOfChargedAbove1GeV() < fChargedAbove1GeVMin) || (EvTag->GetNumOfChargedAbove1GeV() > fChargedAbove1GeVMax))
      return kFALSE; 
  
  if(fChargedAbove3GeVFlag)
    if((EvTag->GetNumOfChargedAbove3GeV() < fChargedAbove3GeVMin) || (EvTag->GetNumOfChargedAbove3GeV() > fChargedAbove3GeVMax))
      return kFALSE; 
  
  if(fChargedAbove10GeVFlag)
    if((EvTag->GetNumOfChargedAbove10GeV() < fChargedAbove10GeVMin) || (EvTag->GetNumOfChargedAbove10GeV() > fChargedAbove10GeVMax))
      return kFALSE; 
  
  if(fMuonsAbove1GeVFlag)
    if((EvTag->GetNumOfMuonsAbove1GeV() < fMuonsAbove1GeVMin) || (EvTag->GetNumOfMuonsAbove1GeV() > fMuonsAbove1GeVMax))
      return kFALSE; 
  
  if(fMuonsAbove3GeVFlag)
    if((EvTag->GetNumOfMuonsAbove3GeV() < fMuonsAbove3GeVMin) || (EvTag->GetNumOfMuonsAbove3GeV() > fMuonsAbove3GeVMax))
      return kFALSE; 
  
  if(fMuonsAbove10GeVFlag)
    if((EvTag->GetNumOfMuonsAbove10GeV() < fMuonsAbove10GeVMin) || (EvTag->GetNumOfMuonsAbove10GeV() > fMuonsAbove10GeVMax))
      return kFALSE; 
  
  if(fElectronsAbove1GeVFlag)
    if((EvTag->GetNumOfElectronsAbove1GeV()  < fElectronsAbove1GeVMin) || (EvTag->GetNumOfElectronsAbove1GeV()  > fElectronsAbove1GeVMax))
      return kFALSE; 
  
  if(fElectronsAbove3GeVFlag)
    if((EvTag->GetNumOfElectronsAbove3GeV() < fElectronsAbove3GeVMin) || (EvTag->GetNumOfElectronsAbove3GeV() > fElectronsAbove3GeVMax))
      return kFALSE; 
  
  if(fElectronsAbove10GeVFlag)
    if((EvTag->GetNumOfElectronsAbove10GeV() < fElectronsAbove10GeVMin) || (EvTag->GetNumOfElectronsAbove10GeV() > fElectronsAbove10GeVMax))
      return kFALSE; 
  
  if(fElectronsFlag)
    if((EvTag->GetNumOfElectrons() < fElectronsMin) || (EvTag->GetNumOfElectrons() > fElectronsMax))
      return kFALSE; 
  
  if(fMuonsFlag)
    if((EvTag->GetNumOfMuons() < fMuonsMin) || (EvTag->GetNumOfMuons() > fMuonsMax))
      return kFALSE; 
  
  if(fPionsFlag)
    if((EvTag->GetNumOfPions() < fPionsMin) || (EvTag->GetNumOfPions() > fPionsMax))
      return kFALSE; 
  
  if(fKaonsFlag)
    if((EvTag->GetNumOfKaons() < fKaonsMin) || (EvTag->GetNumOfKaons() > fKaonsMax))
      return kFALSE; 
  
  if(fProtonsFlag)
    if((EvTag->GetNumOfProtons() < fProtonsMin) || (EvTag->GetNumOfProtons() > fProtonsMax))
      return kFALSE; 
  
  if(fLambdasFlag)
    if((EvTag->GetNumOfLambdas() < fLambdasMin) || (EvTag->GetNumOfLambdas() > fLambdasMax))
      return kFALSE; 
  
  if(fPhotonFlag)
    if((EvTag->GetNumOfPhotons() < fPhotonsMin) || (EvTag->GetNumOfPhotons() > fPhotonsMax))
      return kFALSE; 
  
  if(fPi0sFlag)
    if((EvTag->GetNumOfPi0s() < fPi0sMin) || (EvTag->GetNumOfPi0s() > fPi0sMax))
      return kFALSE; 
  
  if(fNeutronsFlag)
    if((EvTag->GetNumOfNeutrons() < fNeutronsMin) || (EvTag->GetNumOfNeutrons() > fNeutronsMax))
      return kFALSE; 
  
  if(fKaon0sFlag)
    if((EvTag->GetNumOfKaon0s() < fKaon0sMin) || (EvTag->GetNumOfKaon0s() > fKaon0sMax))
      return kFALSE; 
  
  if(fTotalPFlag)
    if((EvTag->GetTotalMomentum() < fTotalPMin) || (EvTag->GetTotalMomentum() > fTotalPMax))
      return kFALSE; 
  
  if(fMeanPtFlag)
    if((EvTag->GetMeanPt() < fMeanPtMin) || (EvTag->GetMeanPt() > fMeanPtMax))
      return kFALSE; 
  
  if(fMaxPtFlag)
    if((EvTag->GetMaxPt() < fMaxPt))
      return kFALSE; 
  
  if(fTotalNeutralPFlag)
    if((EvTag->GetNeutralTotalMomentum() < fTotalNeutralPMin) || (EvTag->GetNeutralTotalMomentum() > fTotalNeutralPMax))
      return kFALSE; 
  
  if(fMeanNeutralPtFlag)
    if((EvTag->GetNeutralMeanPt() < fMeanNeutralPtMin) || (EvTag->GetNeutralMeanPt() >fMeanNeutralPtMax ))
      return kFALSE; 
  
  if(fMaxNeutralPtFlag)
    if((EvTag->GetNeutralMaxPt() < fMaxNeutralPt))
      return kFALSE; 
  
  if(fEventPlaneAngleFlag)
    if((EvTag->GetEventPlaneAngle() < fEventPlaneAngleMin) || (EvTag->GetEventPlaneAngle() > fEventPlaneAngleMax))
      return kFALSE; 
  
  if(fHBTRadiiFlag)
    if((EvTag->GetHBTRadii() < fHBTRadiiMin) || (EvTag->GetHBTRadii() > fHBTRadiiMax))
      return kFALSE; 
  
  return kTRUE;
}
