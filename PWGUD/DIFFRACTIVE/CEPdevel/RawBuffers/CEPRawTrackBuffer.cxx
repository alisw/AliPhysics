// ////////////////////////////////////////////////////////////////////////////
//
// CEP track buffer
//
// structure to hold track information
//
// Authors:
// P. Buehler, paul.buehler@oeaw.ac.at       27.06.2016
//
//
// ____________________________________________________________________________
#include "CEPTrackBuffer.h"
#include "CEPRawTrackBuffer.h"

ClassImp(CEPRawTrackBuffer)

// ____________________________________________________________________________
CEPRawTrackBuffer::CEPRawTrackBuffer()
  : TObject()
  , fTrackLength(CEPTrackBuffer::kdumval)
  , fGlobalChi2(CEPTrackBuffer::kdumval)
  , fPIDITSsignal(CEPTrackBuffer::kdumval)
  , fPIDHMPIDsignal(CEPTrackBuffer::kdumval)
  , fPIDTRDsignal(CEPTrackBuffer::kdumval)
  , fPIDTOFsignalRaw(CEPTrackBuffer::kdumval)
  , fXY(CEPTrackBuffer::kdumval)
  , fZ(CEPTrackBuffer::kdumval)
  , fDx(CEPTrackBuffer::kdumval)
  , fDz(CEPTrackBuffer::kdumval)
  , fTrkPhiOnEMCal(CEPTrackBuffer::kdumval)
  , fTrkEtaOnEMCal(CEPTrackBuffer::kdumval)
  , fTrkPtOnEMCal(CEPTrackBuffer::kdumval)
  , fTrkPOnEMCal(CEPTrackBuffer::kdumval)
{

}

// ____________________________________________________________________________
void CEPRawTrackBuffer::Reset()
{
    fTrackLength     = CEPTrackBuffer::kdumval;
    fGlobalChi2      = CEPTrackBuffer::kdumval;
    fPIDITSsignal    = CEPTrackBuffer::kdumval;
    fPIDHMPIDsignal  = CEPTrackBuffer::kdumval;
    fPIDTRDsignal    = CEPTrackBuffer::kdumval;
    fPIDTOFsignalRaw = CEPTrackBuffer::kdumval;
    fXY              = CEPTrackBuffer::kdumval;
    fZ               = CEPTrackBuffer::kdumval;
    fDx              = CEPTrackBuffer::kdumval;
    fDz              = CEPTrackBuffer::kdumval;
    fTrkPhiOnEMCal   = CEPTrackBuffer::kdumval;
    fTrkEtaOnEMCal   = CEPTrackBuffer::kdumval;
    fTrkPtOnEMCal    = CEPTrackBuffer::kdumval;
    fTrkPOnEMCal     = CEPTrackBuffer::kdumval;
}

// ____________________________________________________________________________
void CEPRawTrackBuffer::SetTrackVariables(AliESDtrack* track)
{
    this->SetTrackLenght( track->GetIntegratedLength() );
    this->SetGlobalChi2( track->GetGlobalChi2() );

    this->SetPIDITSsignal( track->GetITSsignal() );
    this->SetPIDHMPsignal( track->GetHMPIDsignal() );
    this->SetPIDTRDsignal( track->GetTRDsignal() );
    this->SetPIDTOFsignalRaw( track->GetTOFsignalRaw() );

    Float_t impactXY(0.0), impactZ(0.0);
    track->GetImpactParameters(impactXY, impactZ);
    this->SetImpactXY( impactXY );
    this->SetImpactZ( impactZ  );
    this->SetLocalTOFImpX( track->GetTOFsignalDx() );
    this->SetLocalTOFImpZ( track->GetTOFsignalDz() );

    this->SetTrkPhiOnEMC( track->GetTrackPhiOnEMCal() );
    this->SetTrkEtaOnEMC( track->GetTrackEtaOnEMCal() );
    this->SetTrkPtOnEMC( track->GetTrackPtOnEMCal() );
    this->SetTrkPOnEMC( track->GetTrackPOnEMCal() );
}

// ____________________________________________________________________________

