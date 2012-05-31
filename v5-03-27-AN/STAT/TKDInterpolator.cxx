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
  if(!fNodes){
    Warning("TKDInterpolator::SetNode()", "Node array not defined.");
    return;
  }

  Int_t n(GetNTNodes());
  new((*fNodes)[n++]) TKDNodeInfo(node);
}

//_________________________________________________________________
Bool_t TKDInterpolator::Build(Int_t npoints, Int_t ndim)
{
  fNSize = ndim;
  return TKDInterpolatorBase::Build(npoints);
}

//_________________________________________________________________
Int_t TKDInterpolator::GetNodeIndex(const Float_t *p)
{
// 	printf("TKDInterpolator::GetNodeIndex() ...\n");
//   printf("Looking for p[");
//   for(int i=0; i<fNSize; i++) printf("%f ", p[i]);
//   printf("] ...\n");

  for(Int_t i=GetNTNodes(); i--;)
    if(((TKDNodeInfo*)(*fNodes)[i])->Has(p)) return i;
  
  printf("Point p[");
  for(int i=0; i<fNSize; i++) printf("%f ", p[i]);
  printf("] outside range.\n");
  return -1;
}


//_________________________________________________________________
Bool_t TKDInterpolator::SetNode(Int_t inode, const TKDNodeInfo &ref)
{
  if(!fNodes){
    Warning("TKDInterpolator::SetNode()", "Node array not defined.");
    return kFALSE;
  }
  if(inode >= GetNTNodes()){
    Warning("TKDInterpolator::SetNode()", "Node array defined up to %d.", GetNTNodes());
    return kFALSE;
  }
  TKDNodeInfo *node = (TKDNodeInfo*)(*fNodes)[inode];
  (*node) = ref;
  return kTRUE;
}

