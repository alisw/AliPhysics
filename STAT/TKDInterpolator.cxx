#include "TKDInterpolator.h"
#include "TKDNodeInfo.h"

#include "TError.h"
#include "TClonesArray.h"

ClassImp(TKDInterpolator)



//_________________________________________________________________
TKDInterpolator::TKDInterpolator() :
  TKDInterpolatorBase()
{
// Default constructor. To be used with care since in this case building
// of data structure is completly left to the user responsability.
}

//_________________________________________________________________
TKDInterpolator::TKDInterpolator(Int_t ndim, Int_t npoints) :
  TKDInterpolatorBase(ndim)
{
// Wrapper constructor for the TKDTree.

  if(npoints) TKDInterpolatorBase::Build(npoints);
}


//_________________________________________________________________
TKDInterpolator::~TKDInterpolator()
{
}

//_________________________________________________________________
void TKDInterpolator::AddNode(const TKDNodeInfo &node)
{
  if(!fTNodes){
    Warning("TKDInterpolator::SetNode()", "Node array not defined.");
    return;
  }

  new((*fTNodes)[fNTNodes++]) TKDNodeInfo(node);
}

//_________________________________________________________________
void TKDInterpolator::Build(Int_t npoints, Int_t ndim)
{
  fNSize = ndim;
  TKDInterpolatorBase::Build(npoints);
}

//_________________________________________________________________
Int_t TKDInterpolator::GetNodeIndex(const Float_t *p)
{
// 	printf("TKDInterpolator::GetNodeIndex() ...\n");
//   printf("Looking for p[");
//   for(int i=0; i<fNSize; i++) printf("%f ", p[i]);
//   printf("] ...\n");

  for(Int_t i=fNTNodes; i--;)
    if(((TKDNodeInfo*)(*fTNodes)[i])->Has(p)) return i;
  
  printf("Point p[");
  for(int i=0; i<fNSize; i++) printf("%f ", p[i]);
  printf("] outside range.\n");
  return -1;
}


//_________________________________________________________________
Bool_t TKDInterpolator::SetNode(Int_t inode, const TKDNodeInfo &ref)
{
  if(!fTNodes){
    Warning("TKDInterpolator::SetNode()", "Node array not defined.");
    return kFALSE;
  }
  if(inode >= fNTNodes){
    Warning("TKDInterpolator::SetNode()", Form("Node array defined up to %d.", fNTNodes));
    return kFALSE;
  }
  TKDNodeInfo *node = (TKDNodeInfo*)(*fTNodes)[inode];
  (*node) = ref;
  return kTRUE;
}

