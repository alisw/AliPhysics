#include "AliITSUTrackCond.h"
#include "AliITSUAux.h"
#include "AliLog.h"
#include <TMath.h>

using namespace AliITSUAux;
using namespace TMath;

Char_t    AliITSUTrackCond::fgkClSharing = 0;
Int_t     AliITSUTrackCond::fgkMaxBranches = 50;
Int_t     AliITSUTrackCond::fgkMaxCandidates = 500;
Float_t   AliITSUTrackCond::fgkMaxTr2ClChi2  = 50.;
Float_t   AliITSUTrackCond::fgkMissPenalty  = 2.;
Float_t   AliITSUTrackCond::fgkMaxMatchChi2 = 15.;
Float_t   AliITSUTrackCond::fgkMaxITSSAChi2 = 15;
//______________________________________________________________
AliITSUTrackCond::AliITSUTrackCond(int nLayers)
  :fInitDone(kFALSE)
  ,fNLayers(0)
  ,fMaxITSTPCMatchChi2(fgkMaxMatchChi2)
  ,fMaxITSSAChi2(fgkMaxITSSAChi2)
  ,fClSharing(0)
  ,fMaxBranches(0)
  ,fMaxCandidates(0)
  ,fMaxTr2ClChi2(0)
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
  ,fNLayers(0)
  ,fMaxITSTPCMatchChi2(src.fMaxITSTPCMatchChi2)
  ,fMaxITSSAChi2(src.fMaxITSSAChi2)
  ,fClSharing(0)
  ,fMaxBranches(0)
  ,fMaxCandidates(0)
  ,fMaxTr2ClChi2(0)
  ,fMissPenalty(0)
  ,fNSigmaRoadY(0)
  ,fNSigmaRoadZ(0)
  ,fNConditions(src.fNConditions)
  ,fConditions(src.fConditions)
  ,fAuxData(src.fAuxData)
{
  // copy c-tor
  SetNLayers(src.fNLayers);
  for (int i=fNLayers;i--;) {
    SetMaxBranches(i,src.GetMaxBranches(i));
    SetMaxCandidates(i,src.GetMaxCandidates(i));
    SetMaxTr2ClChi2(i,src.GetMaxTr2ClChi2(i));
    SetMissPenalty(i,src.GetMissPenalty(i));
    SetNSigmaRoadY(i,src.GetNSigmaRoadY(i));
    SetNSigmaRoadZ(i,src.GetNSigmaRoadZ(i));
    SetClSharing(i,src.GetClSharing(i));
  }
}

//______________________________________________________________
AliITSUTrackCond& AliITSUTrackCond::operator=(const AliITSUTrackCond& src)
{
  // copy op.
  if (this!=&src) {
    fInitDone = src.fInitDone;
    fNConditions = src.fNConditions;
    fConditions  = src.fConditions;
    fMaxITSTPCMatchChi2 = src.fMaxITSTPCMatchChi2;
    fMaxITSSAChi2       = src.fMaxITSSAChi2;
    //
    SetNLayers(src.fNLayers);
    //
    for (int i=fNLayers;i--;) {
      SetMaxBranches(i,src.GetMaxBranches(i));
      SetMaxCandidates(i,src.GetMaxCandidates(i));
      SetMaxTr2ClChi2(i,src.GetMaxTr2ClChi2(i));
      SetMissPenalty(i,src.GetMissPenalty(i));
      SetNSigmaRoadY(i,src.GetNSigmaRoadY(i));
      SetNSigmaRoadZ(i,src.GetNSigmaRoadZ(i));
      SetClSharing(i,src.GetClSharing(i));
    }
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
    delete[] fMissPenalty;
    delete[] fNSigmaRoadY;
    delete[] fNSigmaRoadZ;
  }
  fNLayers = nLayers;
  //
  if (fNLayers>0) {
    fClSharing     = new Char_t[fNLayers];
    fMaxBranches   = new Short_t[fNLayers];
    fMaxCandidates = new Short_t[fNLayers];
    fMaxTr2ClChi2  = new Float_t[fNLayers];
    fMissPenalty   = new Float_t[fNLayers];
    fNSigmaRoadY   = new Float_t[fNLayers];
    fNSigmaRoadZ   = new Float_t[fNLayers];
    for (int i=fNLayers;i--;) {
      SetClSharing(i,fgkClSharing);
      SetMaxBranches(i,fgkMaxBranches);
      SetMaxCandidates(i,fgkMaxCandidates);
      SetMaxTr2ClChi2(i,fgkMaxTr2ClChi2);
      SetMissPenalty(i,fgkMissPenalty);
      SetNSigmaRoadY(i,-1); // force recalculation
      SetNSigmaRoadZ(i,-1); // force recalculation
    }
  }
  else {
    fClSharing     = 0;
    fMaxBranches   = 0;
    fMaxCandidates = 0;
    fMaxTr2ClChi2  = 0;
    fMissPenalty   = 0;
    fNSigmaRoadY   = 0;
    fNSigmaRoadZ   = 0;
  }
  //
}

//______________________________________________________________
void AliITSUTrackCond::AddGroupPattern(UShort_t patt)
{
  // add new group pattern to last condition
  if (fNConditions<1) AliFatal("Can be called only after AddCondition");
  int ind = fConditions.GetSize();
  fConditions.Set(ind+1);
  fConditions[ind] = patt;
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
  Short_t *arrGrp = (Short_t*)fConditions.GetArray();  
  int ncl = NumberOfBitsSet(patt);
  int cntCond = 0;
  for (int ic=0;ic<fNConditions;ic++) {
    if (arrAux[cntCond+kMinClus]>ncl) {cntCond+=kNAuxSz; continue;} // check number of clusters
    int grAddr = arrAux[cntCond+kCondStart]; // 1st group pattern address in the condition
    Bool_t ok = kTRUE;
    // if every group of the condition does not match, check next contition
    for (int ig=arrAux[cntCond+kNGroups];ig--;) if ( !(patt&arrGrp[grAddr++]) ) {ok = kFALSE; break;}
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
  Short_t *arrGrp = (Short_t*)fConditions.GetArray();  
  int cntCond = 0;
  printf("Conditions set ID=%d : %d entries\n",GetID(),nc);
  for (int i=0;i<nc;i++) {
    printf("#%2d: MinCl:%2d | %d groups :",i,arrAux[cntCond+kMinClus],arrAux[cntCond+kNGroups]);
    int grAddr = arrAux[cntCond+kCondStart];
    for (int ig=arrAux[cntCond+kNGroups];ig--;) {
      printf("{");
      PrintBits(arrGrp[grAddr++], fNLayers);
      printf("}");
    }
    printf("\n");
    cntCond += kNAuxSz;
  }
  printf("Cuts:\t%6s\t%6s\t%4s\t%8s\t%8s\t%8s\t%8s\n", "MaxBrn","MaxCnd","ClSh","Chi2Cl","Mis.Pen.","NSig.Y","NSig.Z");
  for (int i=0;i<fNLayers;i++) {
    printf("Lr%2d:\t%6d\t%6d\t%4d\t%8.1f\t%8.2f\t%8.2f\t%8.2f\n",i,
	   fMaxBranches[i],fMaxCandidates[i],fClSharing[i],fMaxTr2ClChi2[i],fMissPenalty[i],fNSigmaRoadY[i],fNSigmaRoadZ[i]);
  }
  //
  printf("ITS/TPC matching MaxChi2: %.3f\n",fMaxITSTPCMatchChi2);
  printf("ITS_SA BWD fit   MaxChi2: %.3f\n",fMaxITSSAChi2);
  //
}

//______________________________________________________________
void AliITSUTrackCond::Init()
{
  // finalize and check consistency
  if (fInitDone) return;
  //
  for (int ilr=0;ilr<fNLayers;ilr++) {
    if (IsLayerExcluded(ilr)) continue;
    float nsig = Sqrt(GetMaxTr2ClChi2(ilr));
    if (GetNSigmaRoadY(ilr)<0) SetNSigmaRoadY(ilr,nsig);
    if (GetNSigmaRoadZ(ilr)<0) SetNSigmaRoadZ(ilr,nsig);
    //
  }
  if (fMaxITSTPCMatchChi2<1e-6) SetMaxITSTPCMatchChi2(fgkMaxMatchChi2);
  if (fMaxITSSAChi2<1e-6)       SetMaxITSSAChi2(fgkMaxITSSAChi2);
  //
  fInitDone = kTRUE;
}
