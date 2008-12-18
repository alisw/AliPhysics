#include "AliHLTTRDTracklet.h"
#include "AliHLTTRDCluster.h"

/**
 * Default Constructor
 */
//============================================================================
AliHLTTRDTracklet::AliHLTTRDTracklet():
  fTRDtracklet(NULL),
  fSize(sizeof(AliHLTTRDTracklet)),
  fSigmaY(-1),
  fSigmaY2(-1),
  fTilt(-1),
  fPadLength(-1),
  fX0(-1),
  fMeanz(-1),
  fZProb(-1),
  fN(-1),
  fN2(-1),
  fNUsed(-1),
  fFreq(-1),
  fNChange(-1),
  fMPads(-1),
  fC(-1),
  fCC(-1),
  fChi2(-1),
  fChi2Z(-1),
  fDet(-1),
  fMom(-1),
  fdX(-1)
{
  InitArrays();
}

/**
 * Main Constructor
 */
//============================================================================
AliHLTTRDTracklet::AliHLTTRDTracklet(AliTRDseedV1 * inTracklet):
  fTRDtracklet(NULL),
  fSize(sizeof(AliHLTTRDTracklet)),
  fSigmaY(-1),
  fSigmaY2(-1),
  fTilt(-1),
  fPadLength(-1),
  fX0(-1),
  fMeanz(-1),
  fZProb(-1),
  fN(-1),
  fN2(-1),
  fNUsed(-1),
  fFreq(-1),
  fNChange(-1),
  fMPads(-1),
  fC(-1),
  fCC(-1),
  fChi2(-1),
  fChi2Z(-1),
  fDet(-1),
  fMom(-1),
  fdX(-1)
{
  InitArrays();
 
  fTRDtracklet = inTracklet;
  CopyDataMembers();
  //  Print(kFALSE);
  AddClusters();
}

/**
 * Add clusters to the HLTTRDTracklet
 */
//============================================================================
void AliHLTTRDTracklet::AddClusters()
{
  for (Int_t iTimeBin = 0; iTimeBin < knTimebins; iTimeBin++)
    {
//       if (fClusters[iTimeBin])
// 	HLTWarning("Trying to rewrite cluster in tracklet. Not good.");
      AliTRDcluster* trdCluster = fTRDtracklet->GetClusters(iTimeBin);
      if (trdCluster){
	AliHLTTRDCluster * hltCluster = new (GetEndPointer()) AliHLTTRDCluster(trdCluster);
	fSize += hltCluster->GetSize();
	//HLTInfo("cluster %i; adr 0x%x; endPointer 0x%x; fSize %i", iTimeBin, hltCluster, GetEndPointer(), fSize);
	fClusters[iTimeBin] = hltCluster;
      }
      else 
	fClusters[iTimeBin] = NULL;
      
    }
}

/**
 * Copy simple (non-pointer) data members from TRDTracklet to HLTTRDTracklet 
 */
//============================================================================  
void AliHLTTRDTracklet::CopyDataMembers()
{
  for (Int_t i=0; i < 2; i++){
    fYref[i] = fTRDtracklet->GetYref(i);
    fZref[i] = fTRDtracklet->GetZref(i);
  }
  fSigmaY = fTRDtracklet->GetSigmaY();
  fSigmaY2 = fTRDtracklet->GetSigmaY2();
  fTilt = fTRDtracklet->GetTilt();
  fPadLength = fTRDtracklet->GetPadLength();
  
  fX0 = fTRDtracklet->GetX0();
  for (Int_t i = 0; i < knTimebins; i++){
    fX[i] = fTRDtracklet->GetX(i);
    fY[i] = fTRDtracklet->GetY(i);
    fZ[i] = fTRDtracklet->GetZ(i);
    fIndexes[i] = fTRDtracklet->GetIndexes(i);
    fUsable[i] = fTRDtracklet->IsUsable(i);
  }

  for (Int_t i=0; i < 2; i++){
    fYfit[i] = fTRDtracklet->GetYfit(i);
    fZfit[i] = fTRDtracklet->GetZfit(i);
    fYfitR[i] = fTRDtracklet->GetYfitR(i);
    fZfitR[i] = fTRDtracklet->GetZfitR(i);
    fLabels[i] = fTRDtracklet->GetLabels(i);
  }
  
  fMeanz   = fTRDtracklet->GetMeanz();
  fZProb   = fTRDtracklet->GetZProb();
  fN       = fTRDtracklet->GetNbClusters();
  fN2      = fTRDtracklet->GetN2();
  fNUsed   = fTRDtracklet->GetNUsed();
  fFreq    = fTRDtracklet->GetFreq();
  fNChange = fTRDtracklet->GetNChange();
  fMPads = fTRDtracklet->GetMPads();
   
  fC = fTRDtracklet->GetC();
  fCC = fTRDtracklet->GetCC();
  fChi2 = fTRDtracklet->GetChi2();
  fChi2Z = fTRDtracklet->GetChi2Z();
  
  fDet = fTRDtracklet->GetDetector();
  fMom = fTRDtracklet->GetMomentum();
  fdX = fTRDtracklet->GetdX();
  
}

/**
 * Copy data to the output TRDseedV1
 */
//============================================================================
void AliHLTTRDTracklet::ExportTRDTracklet(AliTRDseedV1 *outTracklet)
{
  outTracklet->Reset();
  /* ======= From AliTRDseedV1 ======== */
  outTracklet->SetDetector(fDet);
  outTracklet->SetMomentum(fMom);
  outTracklet->SetDX(fdX);
  
  /* ======= From AliTRDseed ======== */
  for (Int_t i=0; i < 2; i++){
    outTracklet->SetYref(i, fYref[i]);
    outTracklet->SetZref(i, fZref[i]);
  }
  
  outTracklet->SetSigmaY(fSigmaY);
  outTracklet->SetSigmaY2(fSigmaY2);
  outTracklet->SetTilt(fTilt);
  outTracklet->SetPadLength(fPadLength);
  outTracklet->SetX0(fX0);

  for (Int_t i=0; i < knTimebins; i++){
    outTracklet->SetX(i,fX[i]);
    outTracklet->SetX(i,fY[i]);
    outTracklet->SetX(i,fZ[i]);
    outTracklet->SetIndexes(i, fIndexes[i]);
    outTracklet->SetUsable(i, fUsable[i]);
  }
  for (Int_t i=0; i < 2; i++){
    outTracklet->SetYfit(i,fYfit[i]);
    outTracklet->SetYfitR(i,fYfitR[i]);
    outTracklet->SetZfit(i,fZfit[i]);
    outTracklet->SetZfitR(i,fZfitR[i]);
    outTracklet->SetLabels(i,fLabels[i]);
  }
  
  outTracklet->SetMeanz(fMeanz);
  outTracklet->SetZProb(fZProb);
  outTracklet->SetN(fN);
  outTracklet->SetN2(fN2);
  outTracklet->SetNUsed(fNUsed);
  outTracklet->SetFreq(fFreq);
  outTracklet->SetNChange(fNChange);
  outTracklet->SetMPads(fMPads);
  outTracklet->SetC(fC);
  outTracklet->SetCC(fCC);
  outTracklet->SetChi2(fChi2);
  outTracklet->SetChi2Z(fChi2Z);

  for (Int_t iCluster = 0; iCluster < knTimebins; iCluster++){
    if (fClusters[iCluster]){
      AliTRDcluster *trdCluster = new AliTRDcluster();
      fClusters[iCluster]->ExportTRDCluster(trdCluster);
      outTracklet->SetClusters(iCluster, trdCluster);
    }
  }
  
}


/**
 * Init arrays
 */
//============================================================================
void AliHLTTRDTracklet::InitArrays()
{
  for (Int_t i=0; i < knTimebins; i++){
    fClusters[i] = NULL;
  }

  for (Int_t i=0; i < 2; i++){
    fYref[i] = -1;
    fZref[i] = -1;
  }
  for (Int_t i = 0; i < knTimebins; i++){
    fX[i] = -1;
    fY[i] = -1;
    fZ[i] = -1;
    fIndexes[i] = -1;
    fUsable[i] = -1;
  }

  for (Int_t i=0; i < 2; i++){
    fYfit[i] = -1;
    fZfit[i] = -1;
    fYfitR[i] = -1;
    fZfitR[i] = -1;
    fLabels[i] = -1;
  }
  
}

/**
 * Prints main info about tracklet
 */
//============================================================================
void AliHLTTRDTracklet::Print(Bool_t printClusters)
{
  //printf("--hltTracklet-- addr 0x%p(%i); fSize %i\n", this, (int)this, fSize);
  printf("      fDet %i; dMom %f; fdX %f fN %i\n", fDet, fMom, fdX, fN);

  if (printClusters){
    
    for (Int_t iCluster = 0; iCluster < knTimebins; iCluster++){
      printf(" [%i] ",iCluster);
      if (fClusters[iCluster]){
	fClusters[iCluster]->Print();
      }
      else 
	printf("      NULL\n");
    }
  }
  
}

/**
 * Read clusters to TRDtracklet from the memory
 */
//============================================================================
void AliHLTTRDTracklet::ReadClustersFromMemory(void *input)
{
  AliHLTUInt8_t *iterPtr = (AliHLTUInt8_t*) input;
  AliHLTTRDCluster* hltCluster = NULL;
  
  for (Int_t iCluster = 0; iCluster < knTimebins; iCluster++){
    // if we had something in the fClusters[iCluster] before copying,
    // then this entry in the array should not be empty. Fill it.
    if (fClusters[iCluster]){
      hltCluster = (AliHLTTRDCluster*) iterPtr;
      fClusters[iCluster] = hltCluster;
      iterPtr += hltCluster->GetSize();
      //hltCluster->Print();
    }
    
  }
}
