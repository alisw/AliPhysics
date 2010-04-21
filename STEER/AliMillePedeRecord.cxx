#include "AliMillePedeRecord.h"
#include <TMath.h>
#include "AliLog.h"

/**********************************************************************************************/
/* AliMillePedeRecords: class to store the data of single track processing                    */
/* Format: for each measured point the data is stored consequtively                           */
/* INDEX                                                      VALUE                           */
/* -1                                                         residual                        */
/* Local_param_id                                             dResidual/dLocal_param          */
/* ...                                                        ...                             */
/* -2                                                         weight of the measurement       */
/* Global_param_od                                            dResidual/dGlobal_param         */
/* ...                                                        ...                             */
/*                                                                                            */
/* The records for all processed tracks are stored in the temporary tree in orgder to be      */
/* reused for multiple iterations of MillePede                                                */
/*                                                                                            */ 
/* Author: ruben.shahoyan@cern.ch                                                             */
/*                                                                                            */ 
/**********************************************************************************************/

ClassImp(AliMillePedeRecord)

//_____________________________________________________________________________________________
AliMillePedeRecord::AliMillePedeRecord() : 
fSize(0),fNGroups(0),fGroupID(0),fIndex(0),fValue(0),fWeight(1) {SetUniqueID(0);}

//_____________________________________________________________________________________________
AliMillePedeRecord::AliMillePedeRecord(const AliMillePedeRecord& src) : 
  TObject(src),fSize(src.fSize),fNGroups(src.fNGroups),fGroupID(0),fIndex(0),fValue(0),fWeight(src.fWeight)
{
  fIndex = new Int_t[GetDtBufferSize()];
  memcpy(fIndex,src.fIndex,fSize*sizeof(Int_t));
  fValue = new Double_t[GetDtBufferSize()];
  memcpy(fValue,src.fValue,fSize*sizeof(Double_t));
  fGroupID = new UShort_t[GetGrBufferSize()];
  memcpy(fGroupID,src.fGroupID,GetGrBufferSize()*sizeof(UShort_t));
}

//_____________________________________________________________________________________________
AliMillePedeRecord& AliMillePedeRecord::operator=(const AliMillePedeRecord& rhs)
{
  if (this!=&rhs) {
    Reset();
    for (int i=0;i<rhs.GetSize();i++) {
      Double_t val;
      Int_t    ind;
      rhs.GetIndexValue(i,ind,val);
      AddIndexValue(ind,val);
    }
    fWeight = rhs.fWeight;
    for (int i=0;i<rhs.GetNGroups();i++) MarkGroup(rhs.GetGroupID(i));
  }
  return *this;
}

//_____________________________________________________________________________________________
AliMillePedeRecord::~AliMillePedeRecord() {delete[] fIndex; delete[] fValue; delete[] fGroupID;}

//_____________________________________________________________________________________________
void AliMillePedeRecord::Reset()
{
  fSize = 0;
  for (int i=fNGroups;i--;) fGroupID[i] = 0;
  fNGroups = 0;
  fWeight = 1.;
}

//_____________________________________________________________________________________________
void AliMillePedeRecord::Print(const Option_t *) const
{
  if (!fSize) {AliInfo("No data"); return;}
  int cnt=0,point=0;
  //  
  if (fNGroups) printf("Groups: ");
  for (int i=0;i<fNGroups;i++) printf("%4d |",GetGroupID(i)); 
  printf(" Weight: %+.2e\n",fWeight);
  while(cnt<fSize) {
    //
    Double_t resid = fValue[cnt++];
    Double_t *derLoc = GetValue()+cnt;
    int    *indLoc = GetIndex()+cnt;
    int     nLoc = 0;
    while(!IsWeight(cnt)) {nLoc++;cnt++;}
    Double_t weight = GetValue(cnt++);
    Double_t *derGlo = GetValue()+cnt;
    int    *indGlo = GetIndex()+cnt;
    int     nGlo = 0;
    while(!IsResidual(cnt) && cnt<fSize) {nGlo++; cnt++;} 
    //
    printf("\n*** Point#%2d | Residual = %+.4e | Weight = %+.4e\n",point++,resid,weight);
    printf("Locals : "); 
    for (int i=0;i<nLoc;i++) printf("[%5d] %+.4e|",indLoc[i],derLoc[i]); printf("\n");
    printf("Globals: "); 
    for (int i=0;i<nGlo;i++) printf("[%5d] %+.4e|",indGlo[i],derGlo[i]); printf("\n");
    //
  }
  //
}

//_____________________________________________________________________________________________
void AliMillePedeRecord::ExpandDtBuffer(Int_t bfsize)
{
  // add extra space for derivatives data
  bfsize = TMath::Max(bfsize,GetDtBufferSize());
  Int_t *tmpI = new Int_t[bfsize];
  memcpy(tmpI,fIndex, fSize*sizeof(Int_t));
  delete fIndex;
  fIndex = tmpI;
  //
  Double_t *tmpD = new Double_t[bfsize];
  memcpy(tmpD,fValue, fSize*sizeof(Double_t));
  delete fValue;
  fValue = tmpD;
  //
  SetDtBufferSize(bfsize);
}

//_____________________________________________________________________________________________
void AliMillePedeRecord::ExpandGrBuffer(Int_t bfsize)
{
  // add extra space for groupID data 
  bfsize = TMath::Max(bfsize,GetGrBufferSize());
  UShort_t *tmpI = new UShort_t[bfsize];
  memcpy(tmpI,fGroupID, fNGroups*sizeof(UShort_t));
  delete[] fGroupID;
  fGroupID = tmpI;
  for (int i=fNGroups;i<bfsize;i++) fGroupID[i] = 0;
  //
  SetGrBufferSize(bfsize);
}

//_____________________________________________________________________________________________
void AliMillePedeRecord::MarkGroup(Int_t id)
{
  // mark the presence of the detector group
  id++; // groupID is stored as realID+1
  if (fNGroups>0 && fGroupID[fNGroups-1]==id) return; // already there
  if (fNGroups>=GetGrBufferSize()) ExpandGrBuffer(2*(fNGroups+1));
  fGroupID[fNGroups++] = id;  
}

