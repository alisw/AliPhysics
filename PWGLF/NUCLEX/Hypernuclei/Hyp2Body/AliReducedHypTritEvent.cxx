
#include "AliReducedHypTritEvent.h"
#include "TLorentzVector.h"
#include <TClonesArray.h>

ClassImp(AliReducedHypTritTrack)
ClassImp(AliReducedHypTritV0)
ClassImp(AliReducedHypTritEvent)

AliReducedHypTritTrack::AliReducedHypTritTrack() :
  TObject(),
  fP(),
  fPtrack(-1),
  fDca(-1),
  fDcaSigned(0),
  fDedx(-1),
  fDedxSigma(-9999),
  fDedxSigmaTriton(-9999),
  fEta(-999),
  fPhi(-999),
  fTpcNClusters(-1),
  fITSNClusters(-1),
  fTpcChi2(-1),
  fTPCrefit(-1),
  fITSrefit(-1),
  fKink(-1),
  fGeoLength(-1),
  fTrkCutsPassed(0),
	fTRDvalid(0),
	fTRDtrigHNU(0),	 
	fTRDtrigHQU(0),
	fTRDPid(0),
	fTRDnTracklets(0),	
	fTRDPt(0),
	fTRDLayerMask(0),	
	fTRDStack(0),
	fTRDSector(0),
	fTRDPID0(0),
	fTRDPID1(0),
	fTRDPID2(0),
	fTRDPID3(0),
	fTRDPID4(0),
	fTRDPID5(0),
	fTRDSagitta(-1) {
}

AliReducedHypTritTrack::~AliReducedHypTritTrack() {
}

AliReducedHypTritV0::AliReducedHypTritV0() :
  TObject(),
  fPosition(),
  fPvect(),
  fPiTrack(0x0),
  fHeTrack(0x0),
  fP(-1),
  fPt(-1),
  fM(-1),
  fDcaV0(-999),
  fCosPointingAngle(-999),
  fDecayLength(-1),
  fMcTruth(-1),
  fRapidity(-999),
  fParticleSpecies(0),
  fCharge(-999),
  fOnFlyStatus() {
  if (!fPiTrack) fPiTrack = new AliReducedHypTritTrack();
  if (!fHeTrack) fHeTrack = new AliReducedHypTritTrack();
}

AliReducedHypTritV0::~AliReducedHypTritV0() {
}

AliReducedHypTritEvent::AliReducedHypTritEvent() :
  TObject(),
  fVertexPosition(),
  fV0s(0x0),
  fNumberV0s(0),
  fCentrality(-1),
  fRunNumber(0),
  fEvCutsPassed(0),
  fEventID(0), 
	fMagField(0),
  fTrigger(0),
  fTrigMB(0),
  fTrigHNU(0),      
  fTrigHQU(0),    
  fTrigHJT(0),    
  fTrigHSE(0),    
  fTrigV0(0),   
  fTrigSPD(0),     
  fTriggerClasses(),
  fSPDFiredChips0(0),
	fSPDFiredChips1(0),
	fSPDTracklets(0),
	fSPDCluster(0),
	fV0Multiplicity(0),
	fMultV0M(0),		
	fMultOfV0M(0),			
	fMultSPDTracklet(0),	
	fMultSPDCluster(0),	
	fMultRef05(0),			
	fMultRef08(0){
  if (!fV0s) fV0s = new TClonesArray("AliReducedHypTritV0", 100000);
}

AliReducedHypTritEvent::~AliReducedHypTritEvent() {
}

void AliReducedHypTritEvent::ClearEvent() {
  if (fV0s) fV0s->Delete();
  fNumberV0s = 0;
}
