#include <stdlib.h>
#include <iostream.h>

#include "AliTag.h"

ClassImp(AliRunTag)
ClassImp(AliLHCTag)
ClassImp(AliDetectorTag)
ClassImp(AliEventTag)

TClonesArray *AliRunTag::fgEvents = 0;
TClonesArray *AliRunTag::fgDetectors = 0;

//______________________________________________________________________________
AliRunTag::AliRunTag()
{
	if (!fgEvents) fgEvents = new TClonesArray("AliEventTag", 1000);
	fEventTag = fgEvents;
	fNumEvents = 0;

	if (!fgDetectors) fgDetectors = new TClonesArray("AliDetectorTag", 1000);
	fDetectorTag = fgDetectors;
	fNumDetectors = 0;

	fAliceMagneticField = 0.0;
	fAliceRunStartTime = 0;
	fAliceRunStopTime = 0;
	fAliceReconstructionVersion = 0;
	fAliceRunQuality = 0;
	fAliceBeamEnergy = 0.0;
	fAliceCalibrationVersion = 0;
}

//______________________________________________________________________________
AliRunTag::~AliRunTag()
{
}

//______________________________________________________________________________
void AliRunTag::SetLHCTag(Float_t lumin, char *type)
{
	fLHCTag.SetLHCTag(lumin,type);
}

//______________________________________________________________________________
void AliRunTag::SetDetectorTag(AliDetectorTag *DetTag)
{
	TClonesArray &detectors = *fDetectorTag;
	new(detectors[fNumDetectors++]) AliDetectorTag(DetTag);
}

//______________________________________________________________________________
void AliRunTag::AddEventTag(AliEventTag *EvTag)
{
	TClonesArray &events = *fEventTag;
	new(events[fNumEvents++]) AliEventTag(EvTag);
}

//______________________________________________________________________________
void AliRunTag::Clear()
{
	fNumEvents = 0;
	fNumDetectors = 0;
}


//______________________________________________________________________________
//______________________________________________________________________________
AliLHCTag::AliLHCTag()
{
	fLHCLuminosity = -1.0;
}

//______________________________________________________________________________
AliLHCTag::~AliLHCTag()
{
}

//______________________________________________________________________________
//______________________________________________________________________________
AliDetectorTag::AliDetectorTag()
{
	fITS = 0;
	fTPC = 0;
	fTRD = 0;
	fTOF = 0;
	fHMPID = 0;
	fPHOS = 0;
	fZDC = 0;
	fMUON = 0;
	fABSORBER = 0;
	fPMD = 0;
	fRICH = 0;
	fEMCAL = 0;
	fVZERO = 0;
	fTZERO = 0;
}

//______________________________________________________________________________
AliDetectorTag::AliDetectorTag(AliDetectorTag *DetTag)
{
	// DetectorTag copy constructor
	CopyTag(DetTag);
}

//______________________________________________________________________________
AliDetectorTag::~AliDetectorTag()
{
}

//______________________________________________________________________________
void AliDetectorTag::CopyTag(AliDetectorTag *DetTag)
{
	SetITS(DetTag->GetITS());
	SetTPC(DetTag->GetTPC());
	SetTRD(DetTag->GetTRD());
	SetTOF(DetTag->GetTOF());
	SetHMPID(DetTag->GetHMPID());
	SetPHOS(DetTag->GetPHOS());
	SetZDC(DetTag->GetZDC());
	SetMUON(DetTag->GetMUON());
	SetABSORBER(DetTag->GetABSORBER());
	SetPMD(DetTag->GetPMD());
	SetRICH(DetTag->GetRICH());
	SetEMCAL(DetTag->GetEMCAL());
	SetVZERO(DetTag->GetVZERO());
	SetTZERO(DetTag->GetTZERO());
}


//______________________________________________________________________________
//______________________________________________________________________________
AliEventTag::AliEventTag()
{
	fAliceEventId = 0;
	fGUID = 0;

	fNumberOfParticipants = -10;
	fImpactParameter = -10.0;

	fPrimaryVertexX = -100.0;
	fPrimaryVertexY = -100.0;
	fPrimaryVertexZ = -100.0;

	fTriggerInfo = -10;

	fZDCNeutronEnergy = -10.0;
	fZDCProtonEnergy = -10.0;
	fZDCEMEnergy = -10.0;

	fT0VertexZ = -10.0;

	fNumberOfTracks = -10;
	fNumberOfPositiveTracks = -10;
	fNumberOfNegativeTracks = -10;
	fNumberOfNeutralTracks = -10;

	fNumberOfV0s = -10;
	fNumberOfCascades = -10;
	fNumberOfKinks = -10;

	fNumberOfPMDTracks = -10;
	fNumberOfPHOSTracks = -10;
	fNumberOfEMCALTracks = -10;
	fNumberOfFMDTracks = -10;

	fNumberOfJetCandidates = -10;
	fNumberOfHardPhotonsCandidates = -10;

	fNumberOfElectrons = -10;
	fNumberOfMuons = -10;
	fNumberOfPions = -10;
	fNumberOfKaons = -10;
	fNumberOfProtons = -10;
	fNumberOfLambdas = -10;

	fNumberOfJPsiCandidates = -10;
	fNumberOfPsiPrimeCandidates = -10;
	fNumberOfUpsilonCandidates = -10;
	fNumberOfUpsilonPrimeCandidates = -10;
	fNumberOfUpsilonDoublePrimeCandidates = -10;
	fNumberOfCharmParticleCandidates = -10;
	fNumberOfBeautyParticleCandidates = -10;

	fK0PeakPosition = -10.0;
	fK0PeakWidth = -10.0;

	fTotalP = -10.0;
	fMeanPt = -10.0;
	fMaxPt = -10.0;

	fFlowV1 = -10.0;
	fFlowV2 = -10.0;
}


//______________________________________________________________________________
AliEventTag::AliEventTag(AliEventTag *EvTag)
{
	// EventTag copy constructor
	CopyTag(EvTag);
}
//______________________________________________________________________________
AliEventTag::~AliEventTag()
{
}

//______________________________________________________________________________
void AliEventTag::CopyTag(AliEventTag *EvTag)
{
	SetEventId(EvTag->GetEventId());
	SetGUID(EvTag->GetGUID());

	SetNumOfParticipants(EvTag->GetNumOfParticipants());
	SetImpactParameter(EvTag->GetImpactParameter());

	SetVertexX(EvTag->GetVertexX());
	SetVertexY(EvTag->GetVertexY());
	SetVertexZ(EvTag->GetVertexZ());

	SetTrigger(EvTag->GetTrigger());

	SetZDCNeutronEnergy(EvTag->GetZDCNeutronEnergy());
	SetZDCProtonEnergy(EvTag->GetZDCProtonEnergy());
	SetZDCEMEnergy(EvTag->GetZDCEMEnergy());

	SetT0VertexZ(EvTag->GetT0VertexZ());

	SetNumOfTracks(EvTag->GetNumOfTracks());
	SetNumOfPosTracks(EvTag->GetNumOfPosTracks());
	SetNumOfNegTracks(EvTag->GetNumOfNegTracks());
	SetNumOfNeutrTracks(EvTag->GetNumOfNeutrTracks());

	SetNumOfV0s(EvTag->GetNumOfV0s());
	SetNumOfCascades(EvTag->GetNumOfCascades());
	SetNumOfKinks(EvTag->GetNumOfKinks());

	SetNumOfPMDTracks(EvTag->GetNumOfPMDTracks());
	SetNumOfPHOSTracks(EvTag->GetNumOfPHOSTracks());
	SetNumOfEMCALTracks(EvTag->GetNumOfEMCALTracks());
	SetNumOfFMDTracks(EvTag->GetNumOfFMDTracks());

	SetNumOfJetCandidates(EvTag->GetNumOfJetCandidates());
	SetNumOfHardPhotonsCandidates(EvTag->GetNumOfHardPhotonsCandidates());
	SetNumOfJPsiCandidates(EvTag->GetNumOfJPsiCandidates());
	SetNumOfPsiPrimeCandidates(EvTag->GetNumOfPsiPrimeCandidates());
	SetNumOfUpsilonCandidates(EvTag->GetNumOfUpsilonCandidates());
	SetNumOfUpsilonPrimeCandidates(EvTag->GetNumOfUpsilonPrimeCandidates());
	SetNumOfUpsilonDoublePrimeCandidates(EvTag->GetNumOfUpsilonDoublePrimeCandidates());
	SetNumOfCharmCandidates(EvTag->GetNumOfCharmCandidates());
	SetNumOfBeautyCandidates(EvTag->GetNumOfBeautyCandidates());

	SetNumOfElectrons(EvTag->GetNumOfElectrons());
	SetNumOfMuons(EvTag->GetNumOfMuons());
	SetNumOfPions(EvTag->GetNumOfPions());
	SetNumOfKaons(EvTag->GetNumOfKaons());
	SetNumOfProtons(EvTag->GetNumOfProtons());
	SetNumOfLambdas(EvTag->GetNumOfLambdas());

	SetK0Peak(EvTag->GetK0Peak());
	SetK0Width(EvTag->GetK0Width());

	SetTotalMomentum(EvTag->GetTotalMomentum());
	SetMeanPt(EvTag->GetMeanPt());
	SetMaxPt(EvTag->GetMaxPt());

	SetFlowV1(EvTag->GetFlowV1());
	SetFlowV2(EvTag->GetFlowV2());
}
