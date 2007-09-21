#include "TKDTree.h"

#include "TString.h"
#include <string.h>

#ifndef R__ALPHA
templateClassImp(TKDTree)
templateClassImp(TKDNode)
#endif
//////////////////////////////////////////////////////////////////////
// Memory setup of protected data members:
// 
// 	kDataOwner;  // Toggle ownership of the data
// 	fNnodes:
//  size of branch node array, and index of first terminal node
// 	fNDim;       // number of variables
// 	fNpoints;    // number of multidimensional points
// 	fBucketSize; // number of data points per terminal node
// 
// 	fIndPoints : array containing rearanged indexes according to nodes:
// 	| point indexes from 1st terminal node (fBucketSize elements) | ... | point indexes from the last terminal node (<= fBucketSize elements)
// 
// 	Value   **fData;
// 	fRange : array containing the boundaries of the domain:
// 	| 1st dimension (min + max) | 2nd dimension (min + max) | ...
// 	fBoundaries : nodes boundaries
// 	| 1st node {1st dim * 2 elements | 2nd dim * 2 elements | ...} | 2nd node {...} | ...
// 	the nodes are arranged in the following order :
// 	 - first fNnodes are primary nodes
// 	 - after are the terminal nodes
// 	fNodes : array of primary nodes
///////////////////////////////////////////////////////////////////////
//_________________________________________________________________
template <typename  Index, typename Value>
TKDTree<Index, Value>::TKDTree() :
	TObject()
	,kDataOwner(kFALSE)
	,fNnodes(0)
	,fNDim(0)
	,fNDimm(0)
	,fNpoints(0)
	,fBucketSize(0)
	,fNodes(0x0)
	,fRange(0x0)
	,fData(0x0)
	,fBoundaries(0x0)
	,fkNNdim(0)
	,fkNN(0x0)
	,fkNNdist(0x0)
	,fDistBuffer(0x0)
	,fIndBuffer(0x0)
	,fIndPoints(0x0)
	,fRowT0(0)
	,fCrossNode(0)
	,fOffset(0)
{
// Default constructor. Nothing is built
}


//_________________________________________________________________
template <typename  Index, typename Value>
TKDTree<Index, Value>::TKDTree(Index npoints, Index ndim, UInt_t bsize, Value **data) :
	TObject()
	,kDataOwner(kFALSE)
	,fNnodes(0)
	,fNDim(ndim)
	,fNDimm(2*ndim)
	,fNpoints(npoints)
	,fBucketSize(bsize)
	,fNodes(0x0)
	,fRange(0x0)
	,fData(data) //Columnwise!!!!!
	,fBoundaries(0x0)
	,fkNNdim(0)
	,fkNN(0x0)
	,fkNNdist(0x0)
	,fDistBuffer(0x0)
	,fIndBuffer(0x0)
	,fIndPoints(0x0)
	,fRowT0(0)
	,fCrossNode(0)
	,fOffset(0)
{
// Allocate data by hand. See TKDTree(TTree*, const Char_t*) for an automatic method.

	Build();
}

//_________________________________________________________________
template <typename  Index, typename Value>
TKDTree<Index, Value>::~TKDTree()
{
	if (fIndBuffer) delete [] fIndBuffer;
	if (fDistBuffer) delete [] fDistBuffer;
	if (fkNNdist) delete [] fkNNdist;
	if (fkNN) delete [] fkNN;
	if (fNodes) delete [] fNodes;
	if (fIndPoints) delete [] fIndPoints;
	if (fRange) delete [] fRange;
	if (fBoundaries) delete [] fBoundaries;
	if (kDataOwner && fData){
		for(int idim=0; idim<fNDim; idim++) delete [] fData[idim];
		delete [] fData;
	}
}


//_________________________________________________________________
template <typename  Index, typename Value>
void TKDTree<Index, Value>::Build(){
	//
	// Build binning
	//
	// 1. calculate number of nodes
	// 2. calculate first terminal row
	// 3. initialize index array
	// 4. non recursive building of the binary tree
	//
	//1.
	fNnodes = fNpoints/fBucketSize-1;
	if (fNpoints%fBucketSize) fNnodes++;
	
	//2.
	fRowT0=0;
	for (;(fNnodes+1)>(1<<fRowT0);fRowT0++);
	fRowT0-=1;
	//         2 = 2**0 + 1
	//         3 = 2**1 + 1
	//         4 = 2**1 + 2
	//         5 = 2**2 + 1
	//         6 = 2**2 + 2
	//         7 = 2**2 + 3
	//         8 = 2**2 + 4
	
	//3.
	// allocate space for boundaries
	fRange = new Value[2*fNDim];
	fIndPoints= new Index[fNpoints];
	for (Index i=0; i<fNpoints; i++) fIndPoints[i] = i;
	fNodes  = new TKDNode<Index, Value>[fNnodes];
	//
	fCrossNode = (1<<(fRowT0+1))-1;
	if (fCrossNode<fNnodes) fCrossNode = 2*fCrossNode+1;
	//
	//  fOffset = (((fNnodes+1)-(1<<fRowT0)))*2;
	Int_t   over   = (fNnodes+1)-(1<<fRowT0);
	Int_t   filled = ((1<<fRowT0)-over)*fBucketSize;
	fOffset = fNpoints-filled;
	//
// 	printf("Row0      %d\n", fRowT0);
// 	printf("CrossNode %d\n", fCrossNode);
// 	printf("Offset    %d\n", fOffset);
	//
	//
	//4.
	//    stack for non recursive build - size 128 bytes enough
	Int_t rowStack[128];
	Int_t nodeStack[128];
	Int_t npointStack[128];
	Int_t posStack[128];
	Int_t currentIndex = 0;
	Int_t iter =0;
	rowStack[0]    = 0;
	nodeStack[0]   = 0;
	npointStack[0] = fNpoints;
	posStack[0]   = 0;
	//
	Int_t nbucketsall =0;
	while (currentIndex>=0){
		iter++;
		//
		Int_t npoints  = npointStack[currentIndex];
		if (npoints<=fBucketSize) {
			//printf("terminal node : index %d iter %d\n", currentIndex, iter);
			currentIndex--;
			nbucketsall++;
			continue; // terminal node
		}
		Int_t crow     = rowStack[currentIndex];
		Int_t cpos     = posStack[currentIndex];
		Int_t cnode    = nodeStack[currentIndex];		
		TKDNode<Index, Value> * node = &(fNodes[cnode]);
		//printf("currentIndex %d npoints %d node %d\n", currentIndex, npoints, cnode);
		//
		// divide points
		Int_t nbuckets0 = npoints/fBucketSize;           //current number of  buckets
		if (npoints%fBucketSize) nbuckets0++;            //
		Int_t restRows = fRowT0-rowStack[currentIndex];  // rest of fully occupied node row
		if (restRows<0) restRows =0;
		for (;nbuckets0>(2<<restRows); restRows++);
		Int_t nfull = 1<<restRows;
		Int_t nrest = nbuckets0-nfull;
		Int_t nleft =0, nright =0;
		//
		if (nrest>(nfull/2)){
			nleft  = nfull*fBucketSize;
			nright = npoints-nleft;
		}else{
			nright = nfull*fBucketSize/2;
			nleft  = npoints-nright;
		}
	
		//
		//find the axis with biggest spread
		Value maxspread=0;
		Value tempspread, min, max;
		Index axspread=0;
		Value *array;
		for (Int_t idim=0; idim<fNDim; idim++){
			array = fData[idim];
			Spread(npoints, array, fIndPoints+cpos, min, max);
			tempspread = max - min;
			if (maxspread < tempspread) {
				maxspread=tempspread;
				axspread = idim;
			}
			if(cnode==0) {fRange[2*idim] = min; fRange[2*idim+1] = max;}
		}
		array = fData[axspread];
		KOrdStat(npoints, array, nleft, fIndPoints+cpos);
		node->fAxis  = axspread;
		node->fValue = array[fIndPoints[cpos+nleft]];
		//printf("Set node %d : ax %d val %f\n", cnode, node->fAxis, node->fValue);
		//
		//
		npointStack[currentIndex] = nleft;
		rowStack[currentIndex]    = crow+1;
		posStack[currentIndex]    = cpos;
		nodeStack[currentIndex]   = cnode*2+1;
		currentIndex++;
		npointStack[currentIndex] = nright;
		rowStack[currentIndex]    = crow+1;
		posStack[currentIndex]    = cpos+nleft;
		nodeStack[currentIndex]   = (cnode*2)+2;
		//
		if (0){
			// consistency check
			Info("Build()", Form("points %d left %d right %d", npoints, nleft, nright));
			if (nleft<nright) Warning("Build", "Problem Left-Right");
			if (nleft<0 || nright<0) Warning("Build()", "Problem Negative number");
		}
	}
	
	//printf("NBuckets\t%d\n", nbucketsall);
	//fData = 0;
}


//_________________________________________________________________
template <typename  Index, typename Value>
Int_t	TKDTree<Index, Value>::FindNearestNeighbors(const Value *point, const Int_t k, Index *&in, Value *&d)
{
// Find "k" nearest neighbors to "point".
//
// Return true in case of success and false in case of failure.
// The indexes of the nearest k points are stored in the array "in" in
// increasing order of the distance to "point" and the maxim distance
// in "d".
//
// The array "in" is managed internally by the TKDTree.

	Index inode = FindNode(point);
	if(inode < fNnodes){
		Error("FindNearestNeighbors()", "Point outside data range.");
		return kFALSE;
	}

	UInt_t debug = 0;
	Int_t nCalculations = 0;
	
	// allocate working memory
	if(!fDistBuffer){
		fDistBuffer = new Value[fBucketSize];
		fIndBuffer  = new Index[fBucketSize];
	}
	if(fkNNdim < k){
		//if(debug>=1){
			if(fkNN) printf("Reallocate memory %d -> %d \n", fkNNdim, 2*k);
			else printf("Allocate %d memory\n", 2*k);
		//}
		fkNNdim  = 2*k;
		if(fkNN){
			delete [] fkNN; 
			delete [] fkNNdist;
		}
		fkNN     = new Index[fkNNdim];
		fkNNdist = new Value[fkNNdim];
	}
	memset(fkNN, -1, k*sizeof(Index));
	for(int i=0; i<k; i++) fkNNdist[i] = 9999.;
	Index itmp, jtmp; Value ftmp, gtmp;
	
	// calculate number of boundaries with the data domain.	
	Index nBounds = 0;
	if(!fBoundaries) MakeBoundaries();
	Value *bounds = &fBoundaries[inode*2*fNDim];
	for(int idim=0; idim<fNDim; idim++){
		if(bounds[2*idim] == fRange[2*idim]) nBounds++;
		if(bounds[2*idim+1] == fRange[2*idim+1]) nBounds++;
	}
	if(debug>=1) printf("Calculate boundaries [nBounds = %d]\n", nBounds);
	
	// traverse tree
	TKDNode<Index, Value> *node;
	Int_t nodeStack[128], nodeIn[128];
	Index currentIndex = 0;
	nodeStack[0]   = inode;
	nodeIn[0] = inode;
	while(currentIndex>=0){
		Int_t tnode = nodeStack[currentIndex];
		Int_t entry = nodeIn[currentIndex];
		currentIndex--;
		if(debug>=1) printf("Analyse tnode %d entry %d\n", tnode, entry);
		
		// check if node is still eligible
		node = &fNodes[tnode/2 + (tnode%2) - 1];
		if((TMath::Abs(point[node->fAxis] - node->fValue) > fkNNdist[k-1])
											&&
			((point[node->fAxis] > node->fValue && tnode%2) ||
			 (point[node->fAxis] < node->fValue && !tnode%2))) {
			//printf("\tREMOVE NODE %d\n", tnode/2 + (tnode%2) - 1);

			// mark bound
			nBounds++;
			// end all recursions
			if(nBounds==2 * fNDim) break;
			continue;
		}
		
		if(IsTerminal(tnode)){
			if(debug>=2) printf("\tProcess terminal node %d\n", tnode);
			// Link data to terminal node
			Int_t offset = (tnode >= fCrossNode) ? (tnode-fCrossNode)*fBucketSize : fOffset+(tnode-fNnodes)*fBucketSize;
			Index *indexPoints = &fIndPoints[offset];
			Int_t nbs = (tnode == 2*fNnodes) ? fNpoints%fBucketSize : fBucketSize;
			nbs = nbs ? nbs : fBucketSize;
			nCalculations += nbs;
			
			Int_t npoints = 0;
			for(int idx=0; idx<nbs; idx++){
				// calculate distance in the L1 metric
				fDistBuffer[npoints] = 0.;
				for(int idim=0; idim<fNDim; idim++) fDistBuffer[npoints] += TMath::Abs(point[idim] - fData[idim][indexPoints[idx]]);
				// register candidate neighbor
				if(fDistBuffer[npoints] < fkNNdist[k-1]){
					fIndBuffer[npoints] = indexPoints[idx];
					npoints++;
				}
			}
			for(int ibrowse=0; ibrowse<npoints; ibrowse++){
				if(fDistBuffer[ibrowse] >= fkNNdist[k-1]) continue;
				//insert neighbor
				int iNN=0;
				while(fDistBuffer[ibrowse] >= fkNNdist[iNN]) iNN++;
				if(debug>=2) printf("\t\tinsert data %d @ %d distance %f\n", fIndBuffer[ibrowse], iNN, fDistBuffer[ibrowse]);
				
				itmp = fkNN[iNN]; ftmp = fkNNdist[iNN];
				fkNN[iNN]     = fIndBuffer[ibrowse];
				fkNNdist[iNN] = fDistBuffer[ibrowse];
				for(int ireplace=iNN+1; ireplace<k; ireplace++){
					jtmp = fkNN[ireplace]; gtmp = fkNNdist[ireplace];
					fkNN[ireplace]     = itmp; fkNNdist[ireplace] = ftmp;
					itmp = jtmp; ftmp = gtmp;
					if(ftmp == 9999.) break;
				}
				if(debug>=3){
					for(int i=0; i<k; i++){
						if(fkNNdist[i] == 9999.) break;
						printf("\t\tfkNNdist[%d] = %f\n", i, fkNNdist[i]);
					}
				}				
			}
		}
		
		// register parent
		Int_t parent_node = tnode/2 + (tnode%2) - 1;
		if(parent_node >= 0 && entry != parent_node){
			// check if parent node is eligible at all
			node = &fNodes[parent_node];
			if(TMath::Abs(point[node->fAxis] - node->fValue) > fkNNdist[k-1]){
				// mark bound
				nBounds++;
				// end all recursions
				if(nBounds==2 * fNDim) break;
			}
			currentIndex++;
			nodeStack[currentIndex]=tnode/2 + (tnode%2) - 1;
			nodeIn[currentIndex]=tnode;
			if(debug>=2) printf("\tregister %d\n", nodeStack[currentIndex]);
		}
		
		// register children nodes
		Int_t child_node;
		Bool_t kAllow[] = {kTRUE, kTRUE};
		node = &fNodes[tnode];
		if(TMath::Abs(point[node->fAxis] - node->fValue) > fkNNdist[k-1]){
			if(point[node->fAxis] > node->fValue) kAllow[0] = kFALSE;
			else kAllow[1] = kFALSE;
		}
		for(int ic=1; ic<=2; ic++){
			if(!kAllow[ic-1]) continue;
			child_node = (tnode*2)+ic;
			if(child_node < fNnodes + GetNTerminalNodes() && entry != child_node){
				currentIndex++;
				nodeStack[currentIndex] = child_node;
				nodeIn[currentIndex]=tnode;
				if(debug>=2) printf("\tregister %d\n", nodeStack[currentIndex]);
			}
		}		
	}
	// save results
	in = fkNN;
	d  = fkNNdist;
	
	return nCalculations; //kTRUE;
}



//_________________________________________________________________
template <typename  Index, typename Value>
Index TKDTree<Index, Value>::FindNode(const Value * point){
  //
  // find the terminal node to which point belongs
  
	Index stackNode[128], inode;
	Int_t currentIndex =0;
	stackNode[0] = 0;
	//
	while (currentIndex>=0){
		inode    = stackNode[currentIndex];
		currentIndex--;
		if (IsTerminal(inode)) return inode;
		
		TKDNode<Index, Value> & node = fNodes[inode];
		if (point[node.fAxis]<=node.fValue){
			currentIndex++;
			stackNode[currentIndex]=(inode*2)+1;
		}
		if (point[node.fAxis]>=node.fValue){
			currentIndex++;
			stackNode[currentIndex]=(inode*2)+2;
		}
	}
	
	return -1;
}



//_________________________________________________________________
template <typename  Index, typename Value>
void TKDTree<Index, Value>::FindPoint(Value * point, Index &index, Int_t &iter){
  //
  // find the index of point
  // works only if we keep fData pointers
  
	Int_t stackNode[128];
	Int_t currentIndex =0;
	stackNode[0] = 0;
	iter =0;
	//
	while (currentIndex>=0){
		iter++;
		Int_t inode    = stackNode[currentIndex];
		currentIndex--;
		if (IsTerminal(inode)){
			// investigate terminal node
			Int_t indexIP  = (inode >= fCrossNode) ? (inode-fCrossNode)*fBucketSize : (inode-fNnodes)*fBucketSize+fOffset;
			printf("terminal %d indexP %d\n", inode, indexIP);
			for (Int_t ibucket=0;ibucket<fBucketSize;ibucket++){
				Bool_t isOK    = kTRUE;
				indexIP+=ibucket;
				printf("ibucket %d index %d\n", ibucket, indexIP);
				if (indexIP>=fNpoints) continue;
				Int_t index0   = fIndPoints[indexIP];
				for (Int_t idim=0;idim<fNDim;idim++) if (fData[idim][index0]!=point[idim]) isOK = kFALSE;
				if (isOK) index = index0;
			}
			continue;
		}
		
		TKDNode<Index, Value> & node = fNodes[inode];
		if (point[node.fAxis]<=node.fValue){
			currentIndex++;
			stackNode[currentIndex]=(inode*2)+1;
		}
		if (point[node.fAxis]>=node.fValue){
			currentIndex++;
			stackNode[currentIndex]=(inode*2)+2;
		}
	}
  //
  //  printf("Iter\t%d\n",iter);
}

//_________________________________________________________________
template <typename  Index, typename Value>
void TKDTree<Index, Value>::FindInRangeA(Value * point, Value * delta, Index *res , Index &npoints, Index & iter, Int_t bnode)
{
 //
  // Find all points in the range specified by    (point +- range)
  // res     - Resulting indexes are stored in res array 
  // npoints - Number of selected indexes in range
  // NOTICE:
  // For some cases it is better to don't keep data - because of memory consumption
  // If the data are not kept - only boundary conditions are investigated
  // some of the data can be outside of the range
  // What is guranteed in this mode: All of the points in the range are selected + some fraction of others (but close)

	Index stackNode[128];
  iter=0;
  Index currentIndex = 0;
  stackNode[currentIndex] = bnode; 
  while (currentIndex>=0){
    iter++;
    Int_t inode    = stackNode[currentIndex];
    //
    currentIndex--;
    if (!IsTerminal(inode)){
      // not terminal
      TKDNode<Index, Value> * node = &(fNodes[inode]);
      if (point[node->fAxis] - delta[node->fAxis] <  node->fValue) {
	currentIndex++; 
	stackNode[currentIndex]= (inode*2)+1;
      }
      if (point[node->fAxis] + delta[node->fAxis] >= node->fValue){
	currentIndex++; 
	stackNode[currentIndex]= (inode*2)+2;
      }
    }else{
      Int_t indexIP  = 
	(inode >= fCrossNode) ? (inode-fCrossNode)*fBucketSize : (inode-fNnodes)*fBucketSize+fOffset;
      for (Int_t ibucket=0;ibucket<fBucketSize;ibucket++){
	if (indexIP+ibucket>=fNpoints) break;
	res[npoints]   = fIndPoints[indexIP+ibucket];
	npoints++;
      }
    }
  }
  if (fData){
    //
    //  compress rest if data still accesible
    //
    Index npoints2 = npoints;
    npoints=0;
    for (Index i=0; i<npoints2;i++){
      Bool_t isOK = kTRUE;
      for (Index idim = 0; idim< fNDim; idim++){
	if (TMath::Abs(fData[idim][res[i]]- point[idim])>delta[idim]) 
	  isOK = kFALSE;	
      }
      if (isOK){
	res[npoints] = res[i];
	npoints++;
      }
    }    
  }
}


//_________________________________________________________________
template <typename  Index, typename Value>
void TKDTree<Index, Value>::FindInRangeB(Value * point, Value * delta, Index *res , Index &npoints,Index & iter, Int_t bnode)
{
  //
  Long64_t  goldStatus = (1<<(2*fNDim))-1;  // gold status
  Index stackNode[128];
  Long64_t   stackStatus[128]; 
  iter=0;
  Index currentIndex   = 0;
  stackNode[currentIndex]   = bnode; 
  stackStatus[currentIndex] = 0;
  while (currentIndex>=0){
    Int_t inode     = stackNode[currentIndex];
    Long64_t status = stackStatus[currentIndex];
    currentIndex--;
    iter++;
    if (IsTerminal(inode)){
      Int_t indexIP  = 
	(inode >= fCrossNode) ? (inode-fCrossNode)*fBucketSize : (inode-fNnodes)*fBucketSize+fOffset;
      for (Int_t ibucket=0;ibucket<fBucketSize;ibucket++){
	if (indexIP+ibucket>=fNpoints) break;
	res[npoints]   = fIndPoints[indexIP+ibucket];
	npoints++;
      }
      continue;
    }      
    // not terminal    
    if (status == goldStatus){
      Int_t ileft  = inode;
      Int_t iright = inode;
      for (;ileft<fNnodes; ileft   = (ileft<<1)+1);
      for (;iright<fNnodes; iright = (iright<<1)+2);
      Int_t indexL  = 
	(ileft >= fCrossNode) ? (ileft-fCrossNode)*fBucketSize : (ileft-fNnodes)*fBucketSize+fOffset;
      Int_t indexR  = 
	(iright >= fCrossNode) ? (iright-fCrossNode)*fBucketSize : (iright-fNnodes)*fBucketSize+fOffset;
      if (indexL<=indexR){
	Int_t endpoint = indexR+fBucketSize;
	if (endpoint>fNpoints) endpoint=fNpoints;
	for (Int_t ipoint=indexL;ipoint<endpoint;ipoint++){
	  res[npoints]   = fIndPoints[ipoint];
	  npoints++;
	}	  
	continue;
      }
    }
    //
    TKDNode<Index, Value> * node = &(fNodes[inode]);
    if (point[node->fAxis] - delta[node->fAxis] <  node->fValue) {
      currentIndex++; 
      stackNode[currentIndex]= (inode<<1)+1;
      if (point[node->fAxis] + delta[node->fAxis] > node->fValue) 
	stackStatus[currentIndex]= status | (1<<(2*node->fAxis));
    }
    if (point[node->fAxis] + delta[node->fAxis] >= node->fValue){
      currentIndex++; 
      stackNode[currentIndex]= (inode<<1)+2;
      if (point[node->fAxis] - delta[node->fAxis]<node->fValue) 
	stackStatus[currentIndex]= status | (1<<(2*node->fAxis+1));
    }
  }
  if (fData){
    // compress rest
    Index npoints2 = npoints;
    npoints=0;
    for (Index i=0; i<npoints2;i++){
      Bool_t isOK = kTRUE;
      for (Index idim = 0; idim< fNDim; idim++){
	if (TMath::Abs(fData[idim][res[i]]- point[idim])>delta[idim]) 
	  isOK = kFALSE;	
      }
      if (isOK){
	res[npoints] = res[i];
	npoints++;
      }
    }    
  }
}



//_________________________________________________________________
template <typename  Index, typename Value>
void TKDTree<Index, Value>::SetData(Index npoints, Index ndim, UInt_t bsize, Value **data)
{
// TO DO
// 
// Check reconstruction/reallocation of memory of data. Maybe it is not
// necessary to delete and realocate space but only to use the same space

	Clear();
  
	//Columnwise!!!!
   fData = data;
   fNpoints = npoints;
   fNDim = ndim;
   fBucketSize = bsize;

	Build();
}



//_________________________________________________________________
template <typename  Index, typename Value>
void TKDTree<Index, Value>::Spread(Index ntotal, Value *a, Index *index, Value &min, Value &max)
{
  //Value min, max;
  Index i;
  min = a[index[0]];
  max = a[index[0]];
  for (i=0; i<ntotal; i++){
    if (a[index[i]]<min) min = a[index[i]];
    if (a[index[i]]>max) max = a[index[i]];
  }
}


//_________________________________________________________________
template <typename  Index, typename Value>
Value TKDTree<Index, Value>::KOrdStat(Index ntotal, Value *a, Index k, Index *index)
{
  //
   //copy of the TMath::KOrdStat because I need an Index work array

   Index i, ir, j, l, mid;
   Index arr;
   Index temp;

   Index rk = k;
   l=0;
   ir = ntotal-1;
  for(;;) {
      if (ir<=l+1) { //active partition contains 1 or 2 elements
         if (ir == l+1 && a[index[ir]]<a[index[l]])
	    {temp = index[l]; index[l]=index[ir]; index[ir]=temp;}
         Value tmp = a[index[rk]];
         return tmp;
      } else {
         mid = (l+ir) >> 1; //choose median of left, center and right
         {temp = index[mid]; index[mid]=index[l+1]; index[l+1]=temp;}//elements as partitioning element arr.
         if (a[index[l]]>a[index[ir]])  //also rearrange so that a[l]<=a[l+1]
	    {temp = index[l]; index[l]=index[ir]; index[ir]=temp;}

         if (a[index[l+1]]>a[index[ir]])
	    {temp=index[l+1]; index[l+1]=index[ir]; index[ir]=temp;}

         if (a[index[l]]>a[index[l+1]])
    	    {temp = index[l]; index[l]=index[l+1]; index[l+1]=temp;}

         i=l+1;        //initialize pointers for partitioning
         j=ir;
         arr = index[l+1];
         for (;;) {
	    do i++; while (a[index[i]]<a[arr]);
	    do j--; while (a[index[j]]>a[arr]);
	    if (j<i) break;  //pointers crossed, partitioning complete
	       {temp=index[i]; index[i]=index[j]; index[j]=temp;}
         }
         index[l+1]=index[j];
         index[j]=arr;
         if (j>=rk) ir = j-1; //keep active the partition that
         if (j<=rk) l=i;      //contains the k_th element
      }
   }
}

//_________________________________________________________________
template <typename Index, typename Value>
void TKDTree<Index, Value>::MakeBoundaries(Value *range)
{
// Build boundaries for each node.


	if(range) memcpy(fRange, range, 2*fNDim*sizeof(Value));
	// total number of nodes including terminal nodes
	Int_t totNodes = fNnodes + GetNTerminalNodes();
	fBoundaries = new Value[totNodes*2*fNDim];
	Info("MakeBoundaries(Value*)", Form("Allocate boundaries for %d nodes", totNodes));
	
	// loop
	Int_t child_index;
	for(int inode=fNnodes-1; inode>=0; inode--){
		memcpy(&fBoundaries[inode*2*fNDim], fRange, 2*fNDim*sizeof(Value));
				
		// check left child node
		child_index = 2*inode+1;
		if(IsTerminal(child_index)) CookBoundaries(inode, kTRUE);
		for(int idim=0; idim<fNDim; idim++) fBoundaries[2*fNDim*inode+2*idim] = fBoundaries[2*fNDim*child_index+2*idim];
		
		// check right child node
		child_index = 2*inode+2;
		if(IsTerminal(child_index)) CookBoundaries(inode, kFALSE);
		for(int idim=0; idim<fNDim; idim++) fBoundaries[2*fNDim*inode+2*idim+1] = fBoundaries[2*fNDim*child_index+2*idim+1];
	}
}

//_________________________________________________________________
template <typename Index, typename Value>
void TKDTree<Index, Value>::CookBoundaries(Int_t parent_node, Bool_t LEFT)
{
	// define index of this terminal node
	Int_t index = 2 * parent_node + (LEFT ? 1 : 2);
	
	// build and initialize boundaries for this node
	//printf("\tAllocate terminal node %d\n", index);
	memcpy(&fBoundaries[2*fNDim*index], fRange, 2*fNDim*sizeof(Value));
	Bool_t flag[256];  // cope with up to 128 dimensions 
	memset(flag, kFALSE, 2*fNDim);
	Int_t nvals = 0;
		
	// recurse parent nodes
	TKDNode<Index, Value> *node = 0x0;
	while(parent_node >= 0 && nvals < 2 * fNDim){
		node = &(fNodes[parent_node]);
		if(LEFT){
			if(!flag[2*node->fAxis+1]) {
				fBoundaries[2*fNDim*index+2*node->fAxis+1] = node->fValue;
				flag[2*node->fAxis+1] = kTRUE;
				nvals++;
			}
		} else {
			if(!flag[2*node->fAxis]) {
				fBoundaries[2*fNDim*index+2*node->fAxis] = node->fValue;
				flag[2*node->fAxis] = kTRUE;
				nvals++;
			}
		}
		LEFT = parent_node%2;
		parent_node =  parent_node/2 + (parent_node%2) - 1;
	}
}

template class TKDTree<Int_t, Float_t>;
template class TKDTree<Int_t, Double_t>;

