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
  , fPIDITSsignalTuned(CEPTrackBuffer::kdumval)
  , fPIDTOFsignalTuned(CEPTrackBuffer::kdumval)
  , fPIDTPCsignalTuned(CEPTrackBuffer::kdumval)
  , fITSChi2(CEPTrackBuffer::kdumval)
  , fTPCChi2(CEPTrackBuffer::kdumval)
  , fXY(CEPTrackBuffer::kdumval)
  , fZ(CEPTrackBuffer::kdumval)
  , fDCAx(CEPTrackBuffer::kdumval)
  , fDCAy(CEPTrackBuffer::kdumval)
  , fDCAz(CEPTrackBuffer::kdumval)
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
    fTrackLength = CEPTrackBuffer::kdumval;
    fGlobalChi2  = CEPTrackBuffer::kdumval;

    fPIDITSsignal    = CEPTrackBuffer::kdumval;
    fPIDHMPIDsignal  = CEPTrackBuffer::kdumval;
    fPIDTRDsignal    = CEPTrackBuffer::kdumval;
    fPIDTOFsignalRaw = CEPTrackBuffer::kdumval;

    fPIDITSsignalTuned = CEPTrackBuffer::kdumval;
    fPIDTOFsignalTuned = CEPTrackBuffer::kdumval;
    fPIDTPCsignalTuned = CEPTrackBuffer::kdumval;

    fITSChi2 = CEPTrackBuffer::kdumval;
    fTPCChi2 = CEPTrackBuffer::kdumval;

    fXY = CEPTrackBuffer::kdumval;
    fZ  = CEPTrackBuffer::kdumval;

    fDCAx = CEPTrackBuffer::kdumval;
    fDCAy = CEPTrackBuffer::kdumval;
    fDCAz = CEPTrackBuffer::kdumval;

    fDx = CEPTrackBuffer::kdumval;
    fDz = CEPTrackBuffer::kdumval;

    fTrkPhiOnEMCal = CEPTrackBuffer::kdumval;
    fTrkEtaOnEMCal = CEPTrackBuffer::kdumval;
    fTrkPtOnEMCal  = CEPTrackBuffer::kdumval;
    fTrkPOnEMCal   = CEPTrackBuffer::kdumval;
}

// ____________________________________________________________________________
void CEPRawTrackBuffer::SetTrackVariables(AliESDtrack* track, AliESDVertex* vertex)
{
    this->SetTrackLenght( track->GetIntegratedLength() );
    this->SetGlobalChi2( track->GetGlobalChi2() );

    this->SetPIDITSsignal( track->GetITSsignal() );
    this->SetPIDHMPsignal( track->GetHMPIDsignal() );
    this->SetPIDTRDsignal( track->GetTRDsignal() );
    this->SetPIDTOFsignalRaw( track->GetTOFsignalRaw() );

    this->SetPIDITSsigTunedOnData( track->GetITSsignalTunedOnData() );
    this->SetPIDTPCsigTunedOnData( track->GetTPCsignalTunedOnData() );
    this->SetPIDTOFsigTunedOnData( track->GetTOFsignalTunedOnData() );

    this->SetITSChi2( track->GetITSchi2() );
    this->SetTPCChi2( track->GetTPCchi2() );

    Float_t impactXY(0.0), impactZ(0.0);
    track->GetImpactParameters(impactXY, impactZ);
    this->SetImpactXY( impactXY );
    this->SetImpactZ( impactZ  );

    Double_t primaryVertex[3];
    primaryVertex[0] = vertex->GetX();
    primaryVertex[1] = vertex->GetY();
    primaryVertex[2] = vertex->GetZ();

    Double_t v[3], dca[3];
    track->GetXYZ(v);
    for(Int_t i = 0; i < 3; i++){
        dca[i] = v[i] - primaryVertex[i];
    }
    this->SetDCAx( (dca[0]) );
    this->SetDCAy( (dca[1]) );
    this->SetDCAz( (dca[2]) );

    this->SetLocalTOFImpX( track->GetTOFsignalDx() );
    this->SetLocalTOFImpZ( track->GetTOFsignalDz() );

    this->SetTrkPhiOnEMC( track->GetTrackPhiOnEMCal() );
    this->SetTrkEtaOnEMC( track->GetTrackEtaOnEMCal() );
    this->SetTrkPtOnEMC( track->GetTrackPtOnEMCal() );
    this->SetTrkPOnEMC( track->GetTrackPOnEMCal() );
}
// ____________________________________________________________________________

