#include <TTree.h>
#include <TObjArray.h>
#include <TParticle.h>
#include <TMath.h>
#include "AliRun.h"
#include "AliMC.h"
#include <AliITSUSegmentationPix.h>
#include "AliITSUClusterizer.h"
#include "AliITSUClusterPix.h"
#include "AliITSUGeomTGeo.h"
#include "AliITSUSegmentationPix.h"
#include "AliITSdigit.h"
#include "AliITSURecoParam.h"
#include "AliITSUAux.h"
using namespace TMath;
using namespace AliITSUAux;

ClassImp(AliITSUClusterizer)

//______________________________________________________________________________
AliITSUClusterizer::AliITSUClusterizer(Int_t initNRow) 
:  fVolID(-1)
  ,fAllowDiagonalClusterization(kFALSE)
  ,fSegm(0)
  ,fRecoParam(0)
  ,fInputDigits(0)
  ,fInputDigitsReadIndex(0)
  ,fLayerID(0)
  ,fNLabels(0)
  ,fRawData(kFALSE)
  ,fLorAngCorrection(0)
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
void AliITSUClusterizer::Transform(AliITSUClusterPix *cluster,AliITSUClusterizerClusterCand *cand) 
{
  // convert set of digits to cluster data in LOCAL frame
  const double k1to12 = 1./12;
  static int maxLbinDigit = AliITSdigit::GetNTracks();
  //
  fNLabels = 0;
  Int_t n=0;
  cand->fLastDigit->fNext=0;
  double x=0,z=0,xmn=1e9,xmx=-1e9,zmn=1e9,zmx=-1e9,px=0,pz=0;
  float  cx,cz;
  int charge=0;
  for (AliITSUClusterizerClusterDigit *idigit=cand->fFirstDigit;idigit;idigit=idigit->fNext) {
    AliITSdigit* dig = idigit->fDigit;
    fSegm->DetToLocal(dig->GetCoord2(),dig->GetCoord1(),cx,cz); // center of pixel
    //
    // account for possible diod shift
    double ddx,ddz, dx=fSegm->Dpx(dig->GetCoord2()), dz=fSegm->Dpz(dig->GetCoord1());
    fSegm->GetDiodShift(dig->GetCoord2(),dig->GetCoord1(),ddx,ddz);
    //
    charge += dig->GetSignal();
    x += cx;
    z += cz;
    if (cx<xmn) xmn=cx;
    if (cx>xmx) xmx=cx;
    if (cz<zmn) zmn=cz;
    if (cz>zmx) zmx=cz;
    x += ddx*dx;
    z += ddz*dz;
    px += dx;
    pz += dz;
    //
    if (!fRawData) {
      for(Int_t dlab=0;dlab<maxLbinDigit;dlab++){
	Int_t digitlab = (dig->GetTracks())[dlab];
	if(digitlab<0) continue;
	AddLabel(digitlab);
      }
    }
    //
    ++n;
  }
  UChar_t nx=1,nz=1;
  double dx = xmx-xmn, dz = zmx-zmn;
  if (n>1) {
    double fac=1./n;
    x  *= fac;  // mean coordinates
    z  *= fac;
    px *= fac;  // mean pitch
    pz *= fac; 
    nx = 1+Nint(dx/px);
    nz = 1+Nint(dz/pz);
  }
  x -= fLorAngCorrection;  // LorentzAngle correction
  cluster->SetX(x);
  cluster->SetZ(z);
  cluster->SetY(0);
  cluster->SetSigmaZ2(nz>1 ? dz*dz*k1to12 : pz*pz*k1to12);
  cluster->SetSigmaY2(nx>1 ? dx*dx*k1to12 : px*px*k1to12);
  cluster->SetSigmaYZ(0);
  cluster->SetFrameLoc();
  cluster->SetNxNzN(nx,nz,n);
  cluster->SetQ(charge); // note: this is MC info
  //
  if (!fRawData) {
    CheckLabels();
    int nl = Min(kMaxLabInCluster,fNLabels);
    for (int i=nl;i--;) cluster->SetLabel(fCurrLabels[i],i);
  }
  //
  // Set Volume id
  cluster->SetVolumeId(fVolID);
  //    printf("mod %d: (%.4lf,%.4lf)cm\n",fVolID,x,z);
}

//______________________________________________________________________________
void AliITSUClusterizer::CloseCand(AliITSUClusterizerClusterCand *cand) 
{
  // finish cluster
  AliITSUClusterPix *cluster = (AliITSUClusterPix*)NextCluster();
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
  // main algo for MC clustererization
  SetRawData(kFALSE);
  //
  AliITSUClusterizerClusterDigit *iDigit=NextDigit();
  AliITSUClusterizerClusterPart *iPrevRowBegin=0;
  AliITSUClusterizerClusterPart *iNextRowBegin=0;
  AliITSUClusterizerClusterPart *iPrevRow=0;
  AliITSUClusterizerClusterPart *iNextRow=0;
  Int_t lastV=0;
  while (iDigit) {
    if (iDigit->fDigit->GetCoord2()!=lastV) {
      // NEW ROW
      if (iNextRow) iNextRow->fNextInRow=0;
      if (iPrevRowBegin) CloseRemainingParts(iPrevRowBegin);
      if (iDigit->fDigit->GetCoord2()==lastV+1) {
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
      lastV=iDigit->fDigit->GetCoord2(); 
    }
    // skip cluster parts before this digit
    int limCol = iDigit->fDigit->GetCoord1()-fAllowDiagonalClusterization;
    while (iPrevRow && iPrevRow->fUEnd < limCol) {
      iPrevRow=iPrevRow->fNextInRow;
    }
    // find the longest continous line of digits [iDigit,pDigit]=[iDigit,jDigit)
    AliITSUClusterizerClusterCand *cand=AllocCand(); 
    AliITSUClusterizerClusterDigit *pDigit=iDigit;
    AliITSUClusterizerClusterDigit *jDigit=NextDigit();
    cand->fFirstPart=0;
    cand->fFirstDigit=cand->fLastDigit=iDigit; // NB: first diggit is attached differently
    iDigit->fNext=0;
    Int_t lastU =iDigit->fDigit->GetCoord1();
    Int_t lastU1=lastU+1;
    while (jDigit && jDigit->fDigit->GetCoord1()==lastU1 && jDigit->fDigit->GetCoord2()==lastV) {
      pDigit=jDigit;
      jDigit=NextDigit();
      AttachDigitToCand(cand,pDigit);
      ++lastU1;
    }
    if (!fAllowDiagonalClusterization) --lastU1;
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

//______________________________________________________________________________
void AliITSUClusterizer::PrepareLorentzAngleCorrection(Double_t bz)
{
  // calculate parameters for Lorentz Angle correction. Must be called 
  // after setting segmentation and recoparams
  fLorAngCorrection = 0.5*fRecoParam->GetTanLorentzAngle(fLayerID)*bz/kNominalBz*fSegm->Dy();
}

//______________________________________________________________________
void AliITSUClusterizer::CheckLabels() 
{
  // Tries to find mother's labels
  //
  if (fNLabels<1) return;
  AliRunLoader *rl = AliRunLoader::Instance();
  if(!rl) return;
  TTree *trK=(TTree*)rl->TreeK();
  if (!trK) return;
  //
  static int   labS[kMaxLabels];
  static float kine[kMaxLabels];
  Int_t nlabels = fNLabels; 
  Int_t ntracks = gAlice->GetMCApp()->GetNtrack();
  for (Int_t i=fNLabels;i--;) labS[i] = fCurrLabels[i];
  //
  for (Int_t i=0;i<nlabels;i++) {
    Int_t label = labS[i];
    if (label>=ntracks) continue;
    TParticle *part=(TParticle*)gAlice->GetMCApp()->Particle(label);
    kine[i] = part->Energy() - part->GetCalcMass(); // kinetic energy 
    if (kine[i] < 0.02) {    // reduce soft particles from the same cluster
      Int_t m=part->GetFirstMother();
      if (m<0) continue; // primary
      //
      if (part->GetStatusCode()>0) continue;
      //
      // if the parent is within the same cluster, assign parent's label
      for (int j=0;j<nlabels;j++) if (labS[j]==m) { labS[i] = m; break;}
    }
  } 
  //
  if (nlabels>kMaxLabInCluster) { // only 3 labels are stored in cluster, sort in decreasing momentum
    static int ind[kMaxLabels],labSS[kMaxLabels];
    TMath::Sort(nlabels,kine,ind);
    for (int i=nlabels;i--;) labSS[i] = labS[i];
    for (int i=nlabels;i--;) labS[i] = labSS[ind[i]]; 
  }
  //
  //compress labels -- if multi-times the same
  for (Int_t i=0;i<nlabels;i++) fCurrLabels[i]=-2;
  fNLabels = 0;
  int j=0;
  for (int i=0;i<nlabels;i++) {
    for (j=fNLabels;j--;) if (labS[i]==fCurrLabels[j]) break; // the label already there
    if (j<0) fCurrLabels[fNLabels++] = labS[i];
  }
  //
}
