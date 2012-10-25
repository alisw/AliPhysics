#include "AliITSUClusterizer.h"


ClassImp(AliITSUClusterizer)

//______________________________________________________________________________
AliITSUClusterizer::AliITSUClusterizer(Int_t initNRow) 
:  fVolID(-1)
  ,fSegm(0)
  ,fInputDigits(0)
  ,fInputDigitsReadIndex(0)
  ,fOutputClusters(0)
  ,fDigitFreelist(0)
  ,fPartFreelist(0)
  ,fCandFreelist(0)
  ,fDigitFreelistBptrFirst(0)
  ,fDigitFreelistBptrLast(0)
  ,fPartFreelistBptr(0)
  ,fCandFreelistBptr(0)
{
  SetUniqueID(0);
  // c-tor
  SetNRow(initNRow);
}

//______________________________________________________________________________
void AliITSUClusterizer::SetSegmentation(const AliITSUSegmentationPix *segm) 
{
  // attach segmentation, if needed, reinitialize array
  fSegm = segm;
  SetNRow(fSegm->GetNRow()); // reinitialize if needed

}

//______________________________________________________________________________
void AliITSUClusterizer::SetNRow(Int_t nr)
{
  // update buffers
  int nrOld = GetUniqueID();
  if (nrOld>=nr) return;
  SetUniqueID(nr);
  while (fDigitFreelistBptrFirst) {
    AliITSUClusterizerClusterDigit *next = fDigitFreelistBptrFirst[kDigitChunkSize-1].fNext;
    delete[] fDigitFreelistBptrFirst;
    fDigitFreelistBptrFirst=next;
  }
  delete[] fPartFreelistBptr;
  delete[] fCandFreelistBptr;
  //
  fPartFreelist=fPartFreelistBptr = new AliITSUClusterizerClusterPart[nr+1];
  fCandFreelist=fCandFreelistBptr = new AliITSUClusterizerClusterCand[nr+1];
  for (int i=0;i<nr;++i) {
    fPartFreelistBptr[i].fNextInRow = &fPartFreelistBptr[i+1];
    fCandFreelistBptr[i].fNext      = &fCandFreelistBptr[i+1];
  }  
}

//______________________________________________________________________________
AliITSUClusterizer::~AliITSUClusterizer() 
{
  // d-tor
  while (fDigitFreelistBptrFirst) {
    AliITSUClusterizerClusterDigit *next = fDigitFreelistBptrFirst[kDigitChunkSize-1].fNext;
    delete[] fDigitFreelistBptrFirst;
    fDigitFreelistBptrFirst=next;
  }
  delete[] fPartFreelistBptr;
  delete[] fCandFreelistBptr;
}

//______________________________________________________________________________
AliITSUClusterizer::AliITSUClusterizerClusterDigit* AliITSUClusterizer::AllocDigitFreelist() 
{
  // allocate aux space
  AliITSUClusterizerClusterDigit *tmp = new AliITSUClusterizerClusterDigit[kDigitChunkSize];
  for (int i=0;i<kDigitChunkSize-2;++i) tmp[i].fNext=&tmp[i+1];
  tmp[kDigitChunkSize-2].fNext=0;
  tmp[kDigitChunkSize-1].fNext=0;
  if (!fDigitFreelistBptrFirst) fDigitFreelistBptrFirst=tmp;
  else                          fDigitFreelistBptrLast[kDigitChunkSize-1].fNext=tmp;
  fDigitFreelistBptrLast=tmp;
  return tmp;
}

//______________________________________________________________________________
AliITSUClusterizer::AliITSUClusterizerClusterDigit* AliITSUClusterizer::NextDigit() 
{
  // get next digit
  if (fInputDigitsReadIndex<fInputDigits->GetEntriesFast()) {
    AliITSdigit *tmp=static_cast<AliITSdigit*>(fInputDigits->UncheckedAt(fInputDigitsReadIndex++));
    AliITSUClusterizerClusterDigit *digit=AllocDigit();
    digit->fDigit=tmp;
    // IMPORTANT: A lexiographical order (fV,fU) is assumed
    digit->fU=tmp->GetCoord1();
    digit->fV=tmp->GetCoord2();
    return digit;
  }
  else
    return 0;
}

//______________________________________________________________________________
void AliITSUClusterizer::AttachPartToCand(AliITSUClusterizerClusterCand *cand,AliITSUClusterizerClusterPart *part) 
{
  // attach part
  part->fParent = cand;
  part->fPrevInCluster = 0;
  part->fNextInCluster = cand->fFirstPart;
  if (cand->fFirstPart) cand->fFirstPart->fPrevInCluster = part;
  cand->fFirstPart=part;
}

//______________________________________________________________________________
void AliITSUClusterizer::MergeCands(AliITSUClusterizerClusterCand *a,AliITSUClusterizerClusterCand *b) 
{
  // merge cluster parts
  AliITSUClusterizerClusterPart *ipart=b->fFirstPart; 
  AliITSUClusterizerClusterPart *jpart; 
  do {
    jpart=ipart;
    jpart->fParent=a;
  } while ((ipart=ipart->fNextInCluster));
  jpart->fNextInCluster=a->fFirstPart;
  jpart->fNextInCluster->fPrevInCluster=jpart;
  a->fFirstPart=b->fFirstPart;
  // merge digits
  b->fLastDigit->fNext=a->fFirstDigit;
  a->fFirstDigit=b->fFirstDigit;
  //  DeallocCand(b);
}

//______________________________________________________________________________
void AliITSUClusterizer::Transform(AliCluster *cluster,AliITSUClusterizerClusterCand *cand) 
{
  // convert set of digits to clusted data
  Double_t su=0.,sv=0.;
  Int_t n=0;
  cand->fLastDigit->fNext=0;
  for (AliITSUClusterizerClusterDigit *idigit=cand->fFirstDigit;idigit;idigit=idigit->fNext) {
    su+=idigit->fU;
    sv+=idigit->fV;
    ++n;
  }
  Double_t fac=1./n; // Todo: weighting by signal
  Double_t detX = fac*sv, detZ = fac*su;
  
  // Set local coordinates
  Float_t x = detX, z = detZ;
  if (fSegm) { // local coordinates in cm
    x = (-0.5*fSegm->Dx() + detX*fSegm->Dpx() + 0.5*fSegm->Dpx());
    z = (-0.5*fSegm->Dz() + detZ*fSegm->Dpz(0)+ 0.5*fSegm->Dpz(0));
  }
  cluster->SetX(x);
  cluster->SetZ(z);
  // Set Volume id
  cluster->SetVolumeId(fVolID);
  //    printf("mod %d: (%.4lf,%.4lf)cm\n",fVolID,x,z);
}

//______________________________________________________________________________
void AliITSUClusterizer::CloseCand(AliITSUClusterizerClusterCand *cand) 
{
  // finish cluster
  AliCluster *cluster=NextCluster();
  Transform(cluster,cand);
  DeallocDigits(cand->fFirstDigit,cand->fLastDigit);
  DeallocCand(cand);
}

//______________________________________________________________________________
void AliITSUClusterizer::ClosePart(AliITSUClusterizerClusterPart *part) 
{
  // finish cluster part
  AliITSUClusterizerClusterCand *cand=part->fParent;
  DetachPartFromCand(cand,part);
  DeallocPart(part);
  if (!cand->fFirstPart) CloseCand(cand); 
}

//______________________________________________________________________________
void AliITSUClusterizer::CloseRemainingParts(AliITSUClusterizerClusterPart *part) 
{
  // finish what is left
  while (part) {
    AliITSUClusterizerClusterPart *next=part->fNextInRow;
    ClosePart(part);
    part=next;
  } 
}

//______________________________________________________________________________
void AliITSUClusterizer::Clusterize() 
{
  // main algo
  AliITSUClusterizerClusterDigit *iDigit=NextDigit();
  AliITSUClusterizerClusterPart *iPrevRowBegin=0;
  AliITSUClusterizerClusterPart *iNextRowBegin=0;
  AliITSUClusterizerClusterPart *iPrevRow=0;
  AliITSUClusterizerClusterPart *iNextRow=0;
  Int_t lastV=0;
  while (iDigit) {
    if (iDigit->fV!=lastV) {
      // NEW ROW
      if (iNextRow) iNextRow->fNextInRow=0;
      if (iPrevRowBegin) CloseRemainingParts(iPrevRowBegin);
      if (iDigit->fV==lastV+1) {
	iPrevRowBegin=iNextRowBegin;
	iPrevRow     =iNextRowBegin;
      }
      else {
	// there was an empty row
	CloseRemainingParts(iNextRowBegin);
	iPrevRowBegin=0;
	iPrevRow     =0;
      }
      iNextRowBegin=0;
      iNextRow     =0;
      lastV=iDigit->fV; 
    }
    // skip cluster parts before this digit
    while (iPrevRow && iPrevRow->fUEnd<iDigit->fU) {
      iPrevRow=iPrevRow->fNextInRow;
    }
    // find the longest continous line of digits [iDigit,pDigit]=[iDigit,jDigit)
    AliITSUClusterizerClusterCand *cand=AllocCand(); 
    AliITSUClusterizerClusterDigit *pDigit=iDigit;
    AliITSUClusterizerClusterDigit *jDigit=NextDigit();
    cand->fFirstPart=0;
    cand->fFirstDigit=cand->fLastDigit=iDigit; // NB: first diggit is attached differently
    iDigit->fNext=0;
    Int_t lastU =iDigit->fU;
    Int_t lastU1=lastU+1;
    while (jDigit && jDigit->fU==lastU1 && jDigit->fV==lastV) {
      pDigit=jDigit;
      jDigit=NextDigit();
      AttachDigitToCand(cand,pDigit);
      ++lastU1;
    }
    --lastU1;
    AliITSUClusterizerClusterPart *part=AllocPart();
    part->fUBegin=lastU ;
    part->fUEnd  =lastU1;
    AttachPartToCand(cand,part);
    // merge all cluster candidates of the previous line touching this one,
    // advance to the last one, but keep that one the next active one
    AliITSUClusterizerClusterPart *jPrevRow=iPrevRow;
    while (jPrevRow && jPrevRow->fUBegin<=lastU1) {
      if (jPrevRow->fParent!=cand) {
	MergeCands(jPrevRow->fParent,cand);
	DeallocCand(cand);
	cand=jPrevRow->fParent;
      }
      iPrevRow=jPrevRow;
      jPrevRow=jPrevRow->fNextInRow;
    }
    if (iNextRow)
      iNextRow->fNextInRow=part;
    else
      iNextRowBegin=part;
    iNextRow=part;
    iDigit=jDigit;
  }
  // remove remaining cluster parts
  CloseRemainingParts(iPrevRowBegin);
  if (iNextRow) iNextRow->fNextInRow=0;
  CloseRemainingParts(iNextRowBegin);
  return;
}
