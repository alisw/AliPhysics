#include "AliITSUTrackCond.h"
#include "AliITSMFTAux.h"
#include "AliLog.h"
#include <TMath.h>

using namespace AliITSMFTAux;
using namespace TMath;

Char_t    AliITSUTrackCond::fgkClSharing = 0;
Int_t     AliITSUTrackCond::fgkMaxBranches = 50;
Int_t     AliITSUTrackCond::fgkMaxCandidates = 500;
Float_t   AliITSUTrackCond::fgkMaxTr2ClChi2  = 50.;
Float_t   AliITSUTrackCond::fgkMaxChi2GloNrm = 50.;
Float_t   AliITSUTrackCond::fgkMissPenalty  = 2.;
Float_t   AliITSUTrackCond::fgkMaxMatchChi2 = 15.;
Float_t   AliITSUTrackCond::fgkMaxITSSAChi2 = 15;
//______________________________________________________________
AliITSUTrackCond::AliITSUTrackCond(int nLayers)
  :fInitDone(kFALSE)
  ,fActiveLrInner(0)
  ,fActiveLrOuter(0)
  ,fAllowLayers(0)
  ,fNLayers(0)
  ,fMaxClus(0)
  ,fMaxITSTPCMatchChi2(fgkMaxMatchChi2)
  ,fClSharing(0)
  ,fMaxBranches(0)
  ,fMaxCandidates(0)
  ,fMaxITSSAChi2(0)
  ,fMaxTr2ClChi2(0)
  ,fMaxChi2GloNrm(0)
  ,fMissPenalty(0)
  ,fNSigmaRoadY(0)
  ,fNSigmaRoadZ(0)
  ,fNConditions(0)
  ,fConditions(0)
  ,fAuxData(0)
{
  // def c-tor
  if (nLayers) SetNLayers(nLayers);
}

//______________________________________________________________
AliITSUTrackCond::AliITSUTrackCond(const AliITSUTrackCond& src)
  :TObject(src)
  ,fInitDone(src.fInitDone)
  ,fActiveLrInner(src.fActiveLrInner)
  ,fActiveLrOuter(src.fActiveLrOuter)
  ,fAllowLayers(src.fAllowLayers)
  ,fNLayers(0)
  ,fMaxClus(0)
  ,fMaxITSTPCMatchChi2(src.fMaxITSTPCMatchChi2)
  ,fClSharing(0)
  ,fMaxBranches(0)
  ,fMaxCandidates(0)
  ,fMaxITSSAChi2(0)
  ,fMaxTr2ClChi2(0)
  ,fMaxChi2GloNrm(0)
  ,fMissPenalty(0)
  ,fNSigmaRoadY(0)
  ,fNSigmaRoadZ(0)
  ,fNConditions(src.fNConditions)
  ,fConditions(src.fConditions)
  ,fAuxData(src.fAuxData)
{
  // copy c-tor
  SetNLayers(src.fNLayers);
  fMaxClus = src.fMaxClus;
  for (int i=fNLayers;i--;) {
    SetMaxBranches(i,src.GetMaxBranches(i));
    SetMaxCandidates(i,src.GetMaxCandidates(i));
    SetMaxTr2ClChi2(i,src.GetMaxTr2ClChi2(i));
    SetMaxChi2GloNrm(i,src.GetMaxChi2GloNrm(i));
    SetMissPenalty(i,src.GetMissPenalty(i));
    SetNSigmaRoadY(i,src.GetNSigmaRoadY(i));
    SetNSigmaRoadZ(i,src.GetNSigmaRoadZ(i));
    SetClSharing(i,src.GetClSharing(i));
  }
  //
  for (int i=fMaxClus;i--;) SetMaxITSSAChi2(1+i,src.GetMaxITSSAChi2(1+i));
  //
}

//______________________________________________________________
AliITSUTrackCond& AliITSUTrackCond::operator=(const AliITSUTrackCond& src)
{
  // copy op.
  if (this!=&src) {
    fInitDone = src.fInitDone;
    fActiveLrInner = src.fActiveLrInner;
    fActiveLrOuter = src.fActiveLrOuter;
    //
    fAllowLayers = src.fAllowLayers;
    fNConditions = src.fNConditions;
    fConditions  = src.fConditions;
    fMaxITSTPCMatchChi2 = src.fMaxITSTPCMatchChi2;
    //
    SetNLayers(src.fNLayers);
    //
    for (int i=fNLayers;i--;) {
      SetMaxBranches(i,src.GetMaxBranches(i));
      SetMaxCandidates(i,src.GetMaxCandidates(i));
      SetMaxTr2ClChi2(i,src.GetMaxTr2ClChi2(i));
      SetMaxChi2GloNrm(i,src.GetMaxChi2GloNrm(i));
      SetMissPenalty(i,src.GetMissPenalty(i));
      SetNSigmaRoadY(i,src.GetNSigmaRoadY(i));
      SetNSigmaRoadZ(i,src.GetNSigmaRoadZ(i));
      SetClSharing(i,src.GetClSharing(i));
    }
    for (int i=fMaxClus;i--;) SetMaxITSSAChi2(1+i,src.GetMaxITSSAChi2(1+i));
    //
    fAuxData = src.fAuxData;
  }
  return *this;
}

//______________________________________________________________
void AliITSUTrackCond::SetNLayers(int nLayers)
{
  // set number of layers
  fInitDone = kFALSE;
  if (fNLayers) {
    delete[] fClSharing;
    delete[] fMaxBranches;
    delete[] fMaxCandidates;
    delete[] fMaxTr2ClChi2;
    delete[] fMaxChi2GloNrm;
    delete[] fMissPenalty;
    delete[] fNSigmaRoadY;
    delete[] fNSigmaRoadZ;
    delete[] fMaxITSSAChi2;
  }
  fNLayers = nLayers;
  fMaxClus = 2*fNLayers;
  fAllowLayers = 0;
  //
  if (fNLayers>0) {
    fActiveLrInner = 0;
    fActiveLrOuter = fNLayers-1;
    fClSharing     = new Char_t[fNLayers];
    fMaxBranches   = new Short_t[fNLayers];
    fMaxCandidates = new Short_t[fNLayers];
    fMaxTr2ClChi2  = new Float_t[fNLayers];
    fMaxChi2GloNrm = new Float_t[fNLayers];
    fMissPenalty   = new Float_t[fNLayers];
    fNSigmaRoadY   = new Float_t[fNLayers];
    fNSigmaRoadZ   = new Float_t[fNLayers];
    fMaxITSSAChi2  = new Float_t[fMaxClus];
    for (int i=fNLayers;i--;) {
      fAllowLayers |= 0x1<<i;
      SetClSharing(i,fgkClSharing);
      SetMaxBranches(i,fgkMaxBranches);
      SetMaxCandidates(i,fgkMaxCandidates);
      SetMaxTr2ClChi2(i,fgkMaxTr2ClChi2);
      SetMaxChi2GloNrm(i,fgkMaxChi2GloNrm);
      SetMissPenalty(i,fgkMissPenalty);
      SetNSigmaRoadY(i,-1); // force recalculation
      SetNSigmaRoadZ(i,-1); // force recalculation
    }
    for (int i=fMaxClus;i--;) SetMaxITSSAChi2(1+i,fgkMaxITSSAChi2);
  }
  else {
    fClSharing     = 0;
    fMaxBranches   = 0;
    fMaxCandidates = 0;
    fMaxTr2ClChi2  = 0;
    fMaxChi2GloNrm = 0;
    fMissPenalty   = 0;
    fNSigmaRoadY   = 0;
    fNSigmaRoadZ   = 0;
  }
  //
}

//______________________________________________________________
void AliITSUTrackCond::AddGroupPattern(UShort_t patt,Int_t minCl)
{
  // add new group pattern to last condition: the track should have at least minCl clusters at layers given by patt
  if (fNConditions<1) AliFatal("Can be called only after AddCondition");
  if (minCl>int(AliITSMFTAux::kMaxLayers)) AliFatal(Form("Requested Nlayers=%d exceeds max alowed %d",minCl,AliITSMFTAux::kMaxLayers));
  if (minCl<1)                           AliFatal(Form("Requested Nlayers=%d for pattern %x",minCl,patt));
  int ind = fConditions.GetSize();
  fConditions.Set(ind+1);
  fConditions[ind] = (patt&AliITSMFTAux::kMaxLrMask) | (minCl<<kShiftNcl);
  fAuxData[(fNConditions-1)*kNAuxSz + kNGroups]++;
}

//______________________________________________________________
void AliITSUTrackCond::AddNewCondition(Int_t minClusters)
{
  // add new track condition
  fAuxData.Set( (1+fNConditions)*kNAuxSz );
  fAuxData[fNConditions*kNAuxSz+kCondStart] = fConditions.GetSize();
  fAuxData[fNConditions*kNAuxSz+kNGroups]   = 0;
  fAuxData[fNConditions*kNAuxSz+kMinClus]   = minClusters;
  fNConditions++;
  //
}

//______________________________________________________________
Bool_t AliITSUTrackCond::CheckPattern(UShort_t patt) const
{
  // check if the pattern matches to some condition
  Short_t *arrAux = (Short_t*)fAuxData.GetArray();
  UInt_t  *arrGrp = (UInt_t*)fConditions.GetArray();  
  int ncl = NumberOfBitsSet(patt);
  int cntCond = 0;
  for (int ic=0;ic<fNConditions;ic++) {
    if (arrAux[cntCond+kMinClus]>ncl) {cntCond+=kNAuxSz; continue;} // check number of clusters
    int grAddr = arrAux[cntCond+kCondStart]; // 1st group pattern address in the condition
    Bool_t ok = kTRUE;
    // if every group of the condition does not match, check next contition
    for (int ig=arrAux[cntCond+kNGroups];ig--;) {
      UInt_t pattReq = arrGrp[grAddr++];
      UShort_t actLr = (pattReq&AliITSMFTAux::kMaxLrMask)&patt;  // patter of active layers satisfying to mask
      if (!actLr || NumberOfBitsSet(actLr)<int(pattReq>>kShiftNcl)) {ok = kFALSE; break;}
    }
    if (ok) return kTRUE;
    cntCond += kNAuxSz;
  }
  return kFALSE;
}

//______________________________________________________________
void AliITSUTrackCond::Print(Option_t*) const
{
  // print conditions
  int nc = GetNConditions();  
  Short_t *arrAux = (Short_t*)fAuxData.GetArray();
  UInt_t  *arrGrp = (UInt_t*)fConditions.GetArray();  
  int cntCond = 0;
  printf("Conditions set ID=%d : %d entries\n",GetID(),nc);
  for (int i=0;i<nc;i++) {
    printf("#%2d: MinCl:%2d | %d groups :",i,arrAux[cntCond+kMinClus],arrAux[cntCond+kNGroups]);
    int grAddr = arrAux[cntCond+kCondStart];
    for (int ig=arrAux[cntCond+kNGroups];ig--;) {
      UInt_t patt = arrGrp[grAddr];
      printf("{");
      PrintBits(patt, fNLayers);      
      printf("|%d}",patt>>kShiftNcl);
      grAddr++;
    }
    printf("\n");
    cntCond += kNAuxSz;
  }
  if (fAllowLayers) {
    printf("Allowed Layers: ");
    for (int i=0;i<fNLayers;i++) if (!IsLayerExcluded(i)) printf(" %d",i); printf("\n");
  }
  printf("Cuts:\t%6s\t%6s\t%4s\t%8s\t%8s\t%8s\t%8s\t%8s\n", "MaxBrn","MaxCnd","ClSh","Chi2Cl","Chi2Glo","Mis.Pen.","NSig.Y","NSig.Z");
  for (int i=0;i<fNLayers;i++) {
    printf("Lr%2d:\t%6d\t%6d\t%4d\t%8.1f\t%8.2f\t%8.2f\t%8.2f\t%8.2f\n",i,
	   fMaxBranches[i],fMaxCandidates[i],fClSharing[i],fMaxTr2ClChi2[i],fMaxChi2GloNrm[i],fMissPenalty[i],fNSigmaRoadY[i],fNSigmaRoadZ[i]);
  }
  //
  printf("ITS/TPC matching MaxChi2: %.3f\n",fMaxITSTPCMatchChi2);
  printf("ITS_SA BWD fit   MaxChi2 vs Ncl :");
  for (int i=1;i<=fMaxClus;i++) if (GetMaxITSSAChi2(i)>1e-6) printf("\t%d: %.2f",i,GetMaxITSSAChi2(i)); printf("\n");
  //
}

//______________________________________________________________
void AliITSUTrackCond::Init()
{
  // finalize and check consistency
  if (fInitDone) return;
  //
  fActiveLrInner = -1;
  for (int ilr=0;ilr<fNLayers;ilr++) {
    if (IsLayerExcluded(ilr)) continue;
    if (fActiveLrInner<0) fActiveLrInner = ilr;
    fActiveLrOuter = ilr;
    float nsig = Sqrt(2*GetMaxTr2ClChi2(ilr));
    if (GetNSigmaRoadY(ilr)<0) SetNSigmaRoadY(ilr,nsig);
    if (GetNSigmaRoadZ(ilr)<0) SetNSigmaRoadZ(ilr,nsig);
    //
  }
  for (int i=fMaxClus;i--;) if (GetMaxITSSAChi2(1+1)<1e-6)  SetMaxITSSAChi2(1+i,fgkMaxMatchChi2);
  if (fMaxITSTPCMatchChi2<1e-6) SetMaxITSTPCMatchChi2(fgkMaxMatchChi2);
  //
  fInitDone = kTRUE;
}
