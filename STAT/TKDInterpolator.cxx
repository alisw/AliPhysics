#include "TKDInterpolator.h"
#include "TKDNodeInfo.h"

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
		printf("W - TKDInterpolator::SetNode() : Node array not defined.\n");
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
/*	printf("TKDInterpolator::GetNodeIndex() ...\n");
	printf("Looking for p[");
	for(int i=0; i<fNSize; i++) printf("%f ", p[i]);
	printf("] ...\n");*/
	TKDNodeInfo *node;
	for(int inode=0; inode<fNTNodes; inode++){
		node = (TKDNodeInfo*)(*fTNodes)[inode];
		//node->Print();
		if(node->Has(p)) return inode;
	}
	return -1;
}


//_________________________________________________________________
Bool_t TKDInterpolator::SetNode(const Int_t inode, const TKDNodeInfo &ref)
{
	if(!fTNodes){
		printf("W - TKDInterpolator::SetNode() : Node array not defined.\n");
		return kFALSE;
	}
	if(inode >= fNTNodes){
		printf("W - TKDInterpolator::SetNode() : Node array defined up to %d.\n", fNTNodes);
		return kFALSE;
	}
	TKDNodeInfo *node = (TKDNodeInfo*)(*fTNodes)[inode];
	(*node) = ref;
}

