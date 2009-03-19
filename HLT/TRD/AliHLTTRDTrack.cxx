#include "AliHLTTRDTrack.h"
#include "AliHLTTRDTracklet.h"


/**
 * Default Constructor
 */
//============================================================================
AliHLTTRDTrack::AliHLTTRDTrack():
  fSize(sizeof(AliHLTTRDTrack)),
  fTRDtrack(NULL),
  fPIDquality(0),
  fDE(-1),
  fFakeRatio(-1),
  fChi2(-1),
  fMass(-1),
  fLab(-1),
  fN(-1),
  fIntegratedLength(-1),
  fX(-1),
  fAlpha(-1)
{
  InitArrays();
  // not to be used
}

/**
 * Constructor
 * Creates hltTrack from TRDtrackV1
 */
//============================================================================
AliHLTTRDTrack::AliHLTTRDTrack(AliTRDtrackV1 * inTrack):
  fSize(sizeof(AliHLTTRDTrack)),
  fTRDtrack(NULL),
  fPIDquality(0),
  fDE(-1),
  fFakeRatio(-1),
  fChi2(-1),
  fMass(-1),
  fLab(-1),
  fN(-1),
  fIntegratedLength(-1),
  fX(-1),
  fAlpha(-1)
{
  InitArrays();
  
  fTRDtrack = inTrack;

  CopyDataMembers();
  
  AddTracklets();
  
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
 * Copy tracklets to the HLTTRDTrack
 */
//============================================================================
void AliHLTTRDTrack::AddTracklets()
{
  for (Int_t iTracklet = 0; iTracklet < kNplane; iTracklet++)
    {
      //      if (fTracklet[iTracklet])
	// 	HLTWarning("Trying to rewrite tracklets in the Track. Not good.");
      AliTRDseedV1 * trdTracklet = fTRDtrack->GetTracklet(iTracklet);
      if (trdTracklet){
	AliHLTTRDTracklet * hltTracklet = new (GetEndPointer()) AliHLTTRDTracklet(trdTracklet);
	fSize += hltTracklet->GetSize();
	//HLTDebug("tracklet %i; adr 0x%x; endPointer 0x%x; fSize %i", iTracklet, hltTracklet, hltTracklet->GetEndPointer(), fSize);
	fTracklet[iTracklet] = hltTracklet;
      }
      else 
	fTracklet[iTracklet] = NULL;
    }
}


/**
 * Copy data members (except tracklets) from TRDtrackV1 to HLTTRDTrack.
 */
//============================================================================
void AliHLTTRDTrack::CopyDataMembers()
{

  fPIDquality = fTRDtrack->GetNumberOfTrackletsPID();
  
  for(Int_t i = 0; i < AliPID::kSPECIES; i++)
    {
      fPID[i] = fTRDtrack->GetPID(i);
    }
  
  for (Int_t i = 0; i < 3; i++)
    {
      fBudget[i] = fTRDtrack->GetBudget(i);
    }
  fDE = fTRDtrack->GetEdep();
  fFakeRatio = fTRDtrack->GetFakeRatio();
  fChi2 = fTRDtrack->GetChi2();
  fMass = fTRDtrack->GetMass();
  fLab = fTRDtrack->GetLabel();
  fN = fTRDtrack->GetNumberOfClusters();
  fIntegratedLength = fTRDtrack->GetIntegratedLength();
  
  fX = fTRDtrack->GetX();
  fAlpha = fTRDtrack->GetAlpha();
  const Double_t *Ptemp = fTRDtrack->GetParameter();
  for (Int_t i = 0; i < 5; i++)
    {
      fP[i] = Ptemp[i];
    }
  const Double_t *Ctemp = fTRDtrack->GetCovariance();
  for (Int_t i = 0; i < 15; i++)
    {
      fC[i] = Ctemp[i];
    }
  
}

/**
 * Copy data to the output TRDtrackV1
 */
//============================================================================
void AliHLTTRDTrack::ExportTRDTrack(AliTRDtrackV1 *outTrack)
{
  outTrack->Reset();
  
  //Set members from AliTRDtrackV1
  outTrack->SetPIDquality(fPIDquality);
  outTrack->SetEdep(fDE);
  for(Int_t i = 0; i < AliPID::kSPECIES; i++)
    {
      outTrack->SetPID(i,fPID[i]);
    }
  for (Int_t i = 0; i < 3; i++)
    {
      outTrack->SetBudget(i, fBudget[i]);
    }
  for (Int_t iTracklet = 0; iTracklet < kNplane; iTracklet++){
    if (fTracklet[iTracklet]){
      AliTRDseedV1* trdTracklet = new AliTRDseedV1();
      fTracklet[iTracklet]->ExportTRDTracklet(trdTracklet);
      outTrack->SetTracklet(trdTracklet,iTracklet);
    }
  }

  //Set members from AliKalmanTrack
  outTrack->SetFakeRatio(fFakeRatio);
  //outTrack->SetChi2(fChi2);
  outTrack->SetMass(fMass);
  outTrack->SetLabel(fLab);
  (dynamic_cast<AliKalmanTrack*>(outTrack))->SetNumberOfClusters(fN);
  outTrack->SetIntegratedLength(fIntegratedLength);
  
  //Set members from AliExternalTrackParam
  outTrack->Set(fX, fAlpha, fP, fC);
}
  

/**
 * Init of arrays
 */
//============================================================================
void AliHLTTRDTrack::InitArrays()
{
  for(Int_t i = 0; i < kNplane; i++){
    fTracklet[i]=NULL;
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
void AliHLTTRDTrack::Print(Bool_t printTracklets)
{
  //printf("--hltTrack-- addr 0x%p(%i); fSize %i\n", this, (int)this, fSize);
  //printf("   fPIDquality = %s; fX = %f; fAlpha = %f\n", fPIDquality, fX, fAlpha);
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
      for (Int_t i = 0; i < kNplane; i++){
	if (fTracklet[i]){
	  printf("[%i]",i);
	  fTracklet[i]->Print();
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
void AliHLTTRDTrack::ReadTrackletsFromMemory(void* input)
{
  AliHLTUInt8_t *iterPtr = (AliHLTUInt8_t*) input;
  AliHLTTRDTracklet* hltTracklet = NULL;
  
  for (Int_t iTracklet = 0; iTracklet < kNplane; iTracklet++){
    if (fTracklet[iTracklet]){
      hltTracklet = (AliHLTTRDTracklet*) iterPtr;
      hltTracklet->ReadClustersFromMemory(iterPtr+sizeof(AliHLTTRDTracklet));
      fTracklet[iTracklet] = hltTracklet;
      iterPtr += hltTracklet->GetSize();
      //hltTracklet->Print(kFALSE);
    }
  }
}
