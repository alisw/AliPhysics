#include "AliITSUTrackCond.h"
#include "AliITSUAux.h"
#include "AliLog.h"

using namespace AliITSUAux;

//______________________________________________________________
AliITSUTrackCond::AliITSUTrackCond(int nLayers)
  :fNLayers(0)
  ,fMaxBranches(0)
  ,fMaxCandidates(0)
  ,fNConditions(0)
  ,fConditions(0)
  ,fAuxData(0)
{
  // def c-tor
  if (nLayers) SetNLayers(nLayers);
}

//______________________________________________________________
AliITSUTrackCond::AliITSUTrackCond(const AliITSUTrackCond& src)
  :TObject(src),
   fNLayers(src.fNLayers)
  ,fMaxBranches(0)
  ,fMaxCandidates(0)
  ,fNConditions(src.fNConditions)
  ,fConditions(src.fConditions)
  ,fAuxData(src.fAuxData)
{
  // copy c-tor
  if (fNLayers>0) {
    fMaxBranches = new Short_t[fNLayers];
    fMaxCandidates = new Short_t[fNLayers];
    for (int i=fNLayers;i--;) {
      SetMaxBranches(i,src.GetMaxBranches(i));
      SetMaxCandidates(i,src.GetMaxCandidates(i));
    }
  }
}

//______________________________________________________________
AliITSUTrackCond& AliITSUTrackCond::operator=(const AliITSUTrackCond& src)
{
  // copy op.
  if (this!=&src) {
    fNLayers = src.fNLayers;
    fNConditions = src.fNConditions;
    fConditions  = src.fConditions;
    if (fNLayers) {
      delete fMaxBranches;
      delete fMaxCandidates;
      fMaxBranches = new Short_t[fNLayers];
      fMaxCandidates = new Short_t[fNLayers];
      for (int i=fNLayers;i--;) {
	SetMaxBranches(i,src.GetMaxBranches(i));
	SetMaxCandidates(i,src.GetMaxCandidates(i));
      }
    }
    fAuxData = src.fAuxData;
  }
  return *this;
}

//______________________________________________________________
void AliITSUTrackCond::SetNLayers(int nLayers)
{
  // set number of layers
  fNLayers = nLayers;
  if (fNLayers>0) {
    fMaxBranches = new Short_t[fNLayers];
    fMaxCandidates = new Short_t[fNLayers];
    for (int i=fNLayers;i--;) {
      SetMaxBranches(i,kMaxBranches);
      SetMaxCandidates(i,kMaxCandidates);
    }
  }
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
  printf("Max allowed branches/candidates per seed: ");
  for (int i=0;i<fNLayers;i++) printf("L%d: %d/%d ",i,fMaxBranches[i],fMaxCandidates[i]); printf("\n");
}
