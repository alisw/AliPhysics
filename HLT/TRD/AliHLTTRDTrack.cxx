#include "AliHLTTRDTrack.h"
#include "AliHLTTRDTracklet.h"


/**
 * Default Constructor
 */
//============================================================================
AliHLTTRDTrack::AliHLTTRDTrack():
  fDE(-1),
  fFakeRatio(-1),
  fChi2(-1),
  fN(-1),
  fIntegratedLength(-1),
  fX(-1),
  fAlpha(-1),
  fSize(sizeof(AliHLTTRDTrack))
{
  InitArrays();
  // not to be used
}

/**
 * Constructor
 * Creates hltTrack from TRDtrackV1
 */
//============================================================================
AliHLTTRDTrack::AliHLTTRDTrack(const AliTRDtrackV1* const inTrack):
  fDE(inTrack->fDE),
  fFakeRatio(inTrack->fFakeRatio),
  fChi2(inTrack->fChi2),
  fN(inTrack->fN),
  fIntegratedLength(inTrack->GetIntegratedLength()),
  fX(inTrack->GetX()),
  fAlpha(inTrack->GetAlpha()),
  fSize(sizeof(AliHLTTRDTrack))
{
  CopyDataMembers(inTrack);
}

/**
 * Default Destructor
 * In principle should not be empty, but... we do not use it
 */
//============================================================================
AliHLTTRDTrack::~AliHLTTRDTrack()
{
  
}

/**
 * Copy data members (except tracklets) from TRDtrackV1 to HLTTRDTrack.
 */
//============================================================================
void AliHLTTRDTrack::CopyDataMembers(const AliTRDtrackV1* const inTrack)
{
  for(Int_t i = 0; i < AliPID::kSPECIES; i++)
    {
      fPID[i] = inTrack->fPID[i];
    }
  
  for (Int_t i = 0; i < 3; i++)
    {
      fBudget[i] = inTrack->GetBudget(i);
    }
  
  const Double_t* const Ptemp = inTrack->GetParameter();
  for (Int_t i = 0; i < 5; i++)
    {
      fP[i] = Ptemp[i];
    }
  const Double_t* const Ctemp = inTrack->GetCovariance();
  for (Int_t i = 0; i < 15; i++)
    {
      fC[i] = Ctemp[i];
    }

  for (Int_t iTracklet = 0; iTracklet < AliTRDtrackV1::kNplane; iTracklet++)
    {
      AliTRDseedV1* trdTracklet = inTrack->GetTracklet(iTracklet);
      if (trdTracklet){
	AliHLTTRDTracklet* hltTracklet = new (GetEndPointer()) AliHLTTRDTracklet(trdTracklet);
	fSize += hltTracklet->GetSize();
	fTrackletAtPlane[iTracklet] = kTRUE;
      }
      else fTrackletAtPlane[iTracklet] = kFALSE;
    }
}

/**
 * Copy data to the output TRDtrackV1
 */
//============================================================================
void AliHLTTRDTrack::ExportTRDTrack(AliTRDtrackV1* const outTrack) const
{
  //outTrack->Reset(); we always use a new fresh trdtrack as input, so this is useless
  outTrack->SetBit(AliTRDtrackV1::kOwner);  

  outTrack->fDE=fDE;
  outTrack->fFakeRatio=fFakeRatio;
  outTrack->fChi2=fChi2;
  outTrack->fN=fN;
  outTrack->SetIntegratedLength(fIntegratedLength);
  outTrack->Set(fX, fAlpha, fP, fC);

  for(Int_t i = 0; i < AliPID::kSPECIES; i++)
    {
      outTrack->fPID[i] = fPID[i];
    }
  for (Int_t i = 0; i < 3; i++)
    {
      outTrack->SetBudget(i, fBudget[i]);
    }

  AliHLTUInt8_t *iterPtr = (AliHLTUInt8_t*)this+sizeof(*this);
  AliHLTTRDTracklet* hltTracklet;
  
  for (Int_t iTracklet = 0; iTracklet < AliTRDtrackV1::kNplane; iTracklet++){
    if (fTrackletAtPlane[iTracklet]){
      AliTRDseedV1* trdTracklet = new AliTRDseedV1();
      hltTracklet = (AliHLTTRDTracklet*) iterPtr;
      hltTracklet->ExportTRDTracklet(trdTracklet);
      outTrack->SetTracklet(trdTracklet,iTracklet);
      iterPtr += hltTracklet->GetSize();
    }
  }

}
  

/**
 * Init of arrays
 */
//============================================================================
void AliHLTTRDTrack::InitArrays()
{
  for(Int_t i = 0; i < AliTRDtrackV1::kNplane; i++){
    fTrackletAtPlane[i]=kFALSE;
  }

  for(Int_t i = 0; i < AliPID::kSPECIES; i++)
    {
      fPID[i] = -1;
    }
  
  for (Int_t i = 0; i < 3; i++)
    {
      fBudget[i] = -1;
    }
  for (Int_t i = 0; i < 5; i++)
    {
      fP[i] = -1;
    }
  for (Int_t i = 0; i < 15; i++)
    {
      fC[i] = -1;
    }
}

/**
 * Print main values for HLTTrack
 */
//============================================================================
void AliHLTTRDTrack::Print(Bool_t printTracklets) const
{
  printf("--hltTrack-- addr %p; fSize %i\n", this, fSize);
  printf("   fX = %f; fAlpha = %f\n", fX, fAlpha);
  printf("   ");
  
  for(Int_t i = 0; i < AliPID::kSPECIES; i++)
    {
      printf("fPID[%i] = %f; ",i, fPID[i]);
    }
  printf("\n   ");
  
  for (Int_t i = 0; i < 3; i++)
    {
      printf("fBudget[%i] = %f; ",i, fBudget[i]);
    }
  printf("\n");

  if (printTracklets)
    {
      AliHLTUInt8_t *iterPtr = (AliHLTUInt8_t*)this+sizeof(*this);
      AliHLTTRDTracklet* hltTracklet;

      for (Int_t i = 0; i < AliTRDtrackV1::kNplane; i++){
	if (fTrackletAtPlane[i]){
	  printf("[%i]",i);
	  hltTracklet = (AliHLTTRDTracklet*) iterPtr;
	  hltTracklet->Print();
	  iterPtr += hltTracklet->GetSize();
	}
	else
	  printf(" NULL ");
      }
    }

  printf("\n");
}

/**
 * Read tracklets from the memory. 
 * Number of tracklets should be already known
 */
//============================================================================
// void AliHLTTRDTrack::ReadTrackletsFromMemory(void* input)
// {
//   AliHLTUInt8_t *iterPtr = (AliHLTUInt8_t*) input;
//   AliHLTTRDTracklet* hltTracklet = NULL;
  
//   for (Int_t iTracklet = 0; iTracklet < AliTRDtrackV1::kNplane; iTracklet++){
//     if (fTracklet[iTracklet]){
//       hltTracklet = (AliHLTTRDTracklet*) iterPtr;
//       //hltTracklet->ReadClustersFromMemory(iterPtr+sizeof(AliHLTTRDTracklet));
//       fTracklet[iTracklet] = hltTracklet;
//       iterPtr += hltTracklet->GetSize();
//       //hltTracklet->Print(kFALSE);
//     }
//   }
// }

