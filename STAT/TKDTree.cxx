#include "TKDTree.h"
#include "TRandom.h"

#include "TString.h"
#include <string.h>

#ifndef R__ALPHA
templateClassImp(TKDTree)
#endif

/*



Content:
1. What is kd-tree
2. How to cosntruct kdtree - Pseudo code
3. TKDTree implementation
4. User interface



1. What is kdtree ? ( http://en.wikipedia.org/wiki/Kd-tree )

In computer science, a kd-tree (short for k-dimensional tree) is a space-partitioning data structure for organizing points in a k-dimensional space. kd-trees are a useful data structure for several applications, such as searches involving a multidimensional search key (e.g. range searches and nearest neighbour searches). kd-trees are a special case of BSP trees.

A kd-tree uses only splitting planes that are perpendicular to one of the coordinate system axes. This differs from BSP trees, in which arbitrary splitting planes can be used. In addition, in the typical definition every node of a kd-tree, from the root to the leaves, stores a point.[1] This differs from BSP trees, in which leaves are typically the only nodes that contain points (or other geometric primitives). As a consequence, each splitting plane must go through one of the points in the kd-tree. kd-tries are a variant that store data only in leaf nodes. It is worth noting that in an alternative definition of kd-tree the points are stored in its leaf nodes only, although each splitting plane still goes through one of the points.[2]
<img src=".....">

2. Constructing a kd-tree ( Pseudo code)

Since there are many possible ways to choose axis-aligned splitting planes, there are many different ways to construct kd-trees. The canonical method of kd-tree construction has the following constraints:

    * As one moves down the tree, one cycles through the axes used to select the splitting planes. (For example, the root would have an x-aligned plane, the root's children would both have y-aligned planes, the root's grandchildren would all have z-aligned planes, and so on.)
    * At each step, the point selected to create the splitting plane is the median of the points being put into the kd-tree, with respect to their coordinates in the axis being used. (Note the assumption that we feed the entire set of points into the algorithm up-front.)

This method leads to a balanced kd-tree, in which each leaf node is about the same distance from the root. However, balanced trees are not necessarily optimal for all applications.

Note also that it is not required to select the median point. In that case, the result is simply that there is no guarantee that the tree will be balanced. A simple heuristic to avoid coding a complex linear-time median-finding algorithm nor using an O(n log n) sort is to use sort to find the median of a fixed number of randomly selected points to serve as the cut line. Practically this technique results in nicely balanced trees.

Given a list of n points, the following algorithm will construct a balanced kd-tree containing those points.
<img src=".....">

function kdtree (list of points pointList, int depth)
{
    if pointList is empty
        return nil;
    else
    {
        // Select axis based on depth so that axis cycles through all valid values
        var int axis := depth mod k;

        // Sort point list and choose median as pivot element
        select median from pointList;

        // Create node and construct subtrees
        var tree_node node;
        node.location := median;
        node.leftChild := kdtree(points in pointList before median, depth+1);
        node.rightChild := kdtree(points in pointList after median, depth+1);
        return node;
    }
}


3. TKDtree implementation. 

// TKDtree is optimized to minimize memory consumption.
// TKDTree is by default balanced. Nodes are mapped to the 1 dimensional arrays of size
// fNnodes.
// 
// fAxix[fNnodes]  - Division axis (0,1,2,3 ...)
// fValue[fNnodes] - Division value 
//
// Navigation to the Left and right doughter node is expressed by  analytical formula
// Let's parent node id  is inode
// Then: 
// Left   = inode*2+1
// Right  =  (inode+1)*2  
// Let's  daughter node id is inode
// Then:
// Parent = inode/2
// 
// 
// The mapping of the nodes to the 1D array enables us to use the indexing mechanism.
// This is important if some additional kind of information 
// (not directly connected to kd-tree, e.g function fit parameters in the node)
// follow binary tree structure. The mapping can be reused.
// 
// Drawback:   Insertion to the TKDtree is not supported.  
// Advantage:  Random access is supported -  simple formulas  
// 
// Number of division nodes and number of terminals :
// fNnodes = (fNPoints/fBucketSize)
// 
// TKDTree is mapped to one dimensional array.
// Mapping: (see for example the TKDTree::FindNode) 
// 
// The nodes are filled always from left side to the right side:
//
//
// TERMINOLOGY: 
// - inode  -  index of node 
// - irow   -  index of row

// Ideal case:
// Number of terminal nodes = 2^N - N=3
//
//            INode
// irow 0     0                                                                   -  1 inode 
// irow 1     1                              2                                    -  2 inodes
// irow 2     3              4               5               6                    -  4 inodes
// irow 3     7       8      9      10       11     12       13      14           -  8 inodes
//
//
// Non ideal case:
// Number of terminal nodes = 2^N+k  - N=3  k=1

//           INode
// irow 0     0                                                                   - 1 inode
// irow 1     1                              2                                    - 2 inodes
// irow 2     3              4               5               6                    - 3 inodes
// irow 3     7       8      9      10       11     12       13      14           - 8 inodes
// irow 4     15  16                                                              - 2 inodes 
//
//
// The division algorithm:
// The  n points are divided to tho groups.
//
// n=2^k+rest
//
// Left  = 2^k-1 +(rest>2^k-2) ?  2^k-2      : rest
// Right = 2^k-1 +(rest>2^k-2) ?  rest-2^k-2 : 0
//
//
 



*/



//////////////////////////////////////////////////////////////////////
// Memory setup of protected data members:
// 
// 	fDataOwner;  // Toggle ownership of the data
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
  ,fDataOwner(kFALSE)
  ,fNnodes(0)
  ,fNDim(0)
  ,fNDimm(0)
  ,fNpoints(0)
  ,fBucketSize(0)
  ,fAxis(0x0)
  ,fValue(0x0)
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
  ,fDataOwner(kFALSE)
  ,fNnodes(0)
  ,fNDim(ndim)
  ,fNDimm(2*ndim)
  ,fNpoints(npoints)
  ,fBucketSize(bsize)
  ,fAxis(0x0)
  ,fValue(0x0)
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
  // Destructor
  if (fIndBuffer) delete [] fIndBuffer;
  if (fDistBuffer) delete [] fDistBuffer;
  if (fkNNdist) delete [] fkNNdist;
  if (fkNN) delete [] fkNN;
  if (fAxis) delete [] fAxis;
  if (fValue) delete [] fValue;
  if (fIndPoints) delete [] fIndPoints;
  if (fRange) delete [] fRange;
  if (fBoundaries) delete [] fBoundaries;
  if (fDataOwner && fData){
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
  // Teh tree is recursivaly divided. The division point is choosen in such a way, that the left side 
  // alway has 2^k points, while at the same time trying 
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
  fAxis  = new UChar_t[fNnodes];
  fValue = new Value[fNnodes];
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
      if(cnode) continue;
      //printf("set %d %6.3f %6.3f\n", idim, min, max);
      fRange[2*idim] = min; fRange[2*idim+1] = max;
    }
    array = fData[axspread];
    KOrdStat(npoints, array, nleft, fIndPoints+cpos);
    fAxis[cnode]  = axspread;
    fValue[cnode] = array[fIndPoints[cpos+nleft]];
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
}


//_________________________________________________________________
template <typename  Index, typename Value>
Bool_t	TKDTree<Index, Value>::FindNearestNeighbors(const Value *point, const Int_t k, Index *&in, Value *&d)
{
  //
  // NOT MI - to be checked
  //
  //
  // Find "k" nearest neighbors to "point".
  //
  // Return true in case of success and false in case of failure.
  // The indexes of the nearest k points are stored in the array "in" in
  // increasing order of the distance to "point" and the maxim distance
  // in "d".
  //
  // The arrays "in" and "d" are managed internally by the TKDTree.
  
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
    if(debug>=1){
      if(fkNN) printf("Reallocate memory %d -> %d \n", fkNNdim, 2*k);
      else printf("Allocate %d memory\n", 2*k);
    }
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
  if(!fBoundaries) MakeBoundaries();
  Value *bounds = &fBoundaries[inode*2*fNDim];
  Index nBounds = 0;
  for(int idim=0; idim<fNDim; idim++){
    if((bounds[2*idim] - fRange[2*idim]) < 1.E-10) nBounds++;
    if((bounds[2*idim+1] - fRange[2*idim+1]) < 1.E-10) nBounds++;
  }
  if(debug>=1) printf("Calculate boundaries [nBounds = %d]\n", nBounds);
  
  
  // traverse tree
  UChar_t ax;
  Value   val, dif;
  Int_t nAllNodes = fNnodes + fNpoints/fBucketSize + ((fNpoints%fBucketSize)?1:0);
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
    Int_t inode = (tnode-1)>>1; //tnode/2 + (tnode%2) - 1;
    ax = fAxis[inode];
    val = fValue[inode];
    dif = point[ax] - val;
    if((TMath::Abs(dif) > fkNNdist[k-1]) &&
       ((dif > 0. && (tnode&1)) || (dif < 0. && !(tnode&1)))) {
      if(debug>=1) printf("\tremove %d\n", (tnode-1)>>1);
      
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
    Int_t pn = (tnode-1)>>1;//tnode/2 + (tnode%2) - 1;
    if(pn >= 0 && entry != pn){
      // check if parent node is eligible at all
      if(TMath::Abs(point[fAxis[pn]] - fValue[pn]) > fkNNdist[k-1]){
	// mark bound
	nBounds++;
	// end all recursions
	if(nBounds==2 * fNDim) break;
      }
      currentIndex++;
      nodeStack[currentIndex]=pn;
      nodeIn[currentIndex]=tnode;
      if(debug>=2) printf("\tregister %d\n", nodeStack[currentIndex]);
    }
    if(IsTerminal(tnode)) continue;
    
    // register children nodes
    Int_t cn;
    Bool_t kAllow[] = {kTRUE, kTRUE};
    ax = fAxis[tnode];
    val = fValue[tnode];
    if(TMath::Abs(point[ax] - val) > fkNNdist[k-1]){
      if(point[ax] > val) kAllow[0] = kFALSE;
      else kAllow[1] = kFALSE;
    }
    for(int ic=1; ic<=2; ic++){
      if(!kAllow[ic-1]) continue;
      cn = (tnode<<1)+ic;
      if(cn < nAllNodes && entry != cn){
	currentIndex++;
	nodeStack[currentIndex] = cn;
	nodeIn[currentIndex]=tnode;
	if(debug>=2) printf("\tregister %d\n", nodeStack[currentIndex]);
      }
    }		
  }
  // save results
  in = fkNN;
  d  = fkNNdist;  
  return kTRUE;
}



//_________________________________________________________________
template <typename  Index, typename Value>
Index TKDTree<Index, Value>::FindNode(const Value * point){
  //
  // find the terminal node to which point belongs

  Index stackNode[128], inode;
  Int_t currentIndex =0;
  stackNode[0] = 0;
  while (currentIndex>=0){
    inode    = stackNode[currentIndex];
    if (IsTerminal(inode)) return inode;
    
    currentIndex--;
    if (point[fAxis[inode]]<=fValue[inode]){
      currentIndex++;
      stackNode[currentIndex]=(inode<<1)+1;
    }
    if (point[fAxis[inode]]>=fValue[inode]){
      currentIndex++;
      stackNode[currentIndex]=(inode+1)<<1;
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
    
    if (point[fAxis[inode]]<=fValue[inode]){
      currentIndex++;
      stackNode[currentIndex]=(inode*2)+1;
    }
    if (point[fAxis[inode]]>=fValue[inode]){
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
      //TKDNode<Index, Value> * node = &(fNodes[inode]);
      if (point[fAxis[inode]] - delta[fAxis[inode]] <  fValue[inode]) {
	currentIndex++; 
	stackNode[currentIndex]= (inode*2)+1;
      }
      if (point[fAxis[inode]] + delta[fAxis[inode]] >= fValue[inode]){
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
    // TKDNode<Index, Value> * node = &(fNodes[inode]);
    if (point[fAxis[inode]] - delta[fAxis[inode]] <  fValue[inode]) {
      currentIndex++; 
      stackNode[currentIndex]= (inode<<1)+1;
      if (point[fAxis[inode]] + delta[fAxis[inode]] > fValue[inode])
	stackStatus[currentIndex]= status | (1<<(2*fAxis[inode]));
    }
    if (point[fAxis[inode]] + delta[fAxis[inode]] >= fValue[inode]){
      currentIndex++; 
      stackNode[currentIndex]= (inode<<1)+2;
      if (point[fAxis[inode]] - delta[fAxis[inode]]<fValue[inode])
	stackStatus[currentIndex]= status | (1<<(2*fAxis[inode]+1));
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
void TKDTree<Index, Value>::Spread(Index ntotal, Value *a, Index *index, Value &min, Value &max) const
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
Value TKDTree<Index, Value>::KOrdStat(Index ntotal, Value *a, Index k, Index *index) const
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


	if(range) memcpy(fRange, range, fNDimm*sizeof(Value));
	// total number of nodes including terminal nodes
	Int_t totNodes = fNnodes + fNpoints/fBucketSize + ((fNpoints%fBucketSize)?1:0);
	fBoundaries = new Value[totNodes*fNDimm];
	//Info("MakeBoundaries(Value*)", Form("Allocate boundaries for %d nodes", totNodes));

	
	// loop
	Value *tbounds = 0x0, *cbounds = 0x0;
	Int_t cn;
	for(int inode=fNnodes-1; inode>=0; inode--){
		tbounds = &fBoundaries[inode*fNDimm];
		memcpy(tbounds, fRange, fNDimm*sizeof(Value));
				
		// check left child node
		cn = (inode<<1)+1;
		if(IsTerminal(cn)) CookBoundaries(inode, kTRUE);
		cbounds = &fBoundaries[fNDimm*cn];
		for(int idim=0; idim<fNDim; idim++) tbounds[idim<<1] = cbounds[idim<<1];
		
		// check right child node
		cn = (inode+1)<<1;
		if(IsTerminal(cn)) CookBoundaries(inode, kFALSE);
		cbounds = &fBoundaries[fNDimm*cn];
		for(int idim=0; idim<fNDim; idim++) tbounds[(idim<<1)+1] = cbounds[(idim<<1)+1];
	}
}

//_________________________________________________________________
template <typename Index, typename Value>
void TKDTree<Index, Value>::CookBoundaries(const Int_t node, Bool_t LEFT)
{
	// define index of this terminal node
	Int_t index = (node<<1) + (LEFT ? 1 : 2);
	//Info("CookBoundaries()", Form("Node %d", index));
	
	// build and initialize boundaries for this node
	Value *tbounds = &fBoundaries[fNDimm*index];
	memcpy(tbounds, fRange, fNDimm*sizeof(Value));
	Bool_t flag[256];  // cope with up to 128 dimensions 
	memset(flag, kFALSE, fNDimm);
	Int_t nvals = 0;
		
	// recurse parent nodes
	Int_t pn = node;
	while(pn >= 0 && nvals < fNDimm){
		if(LEFT){
			index = (fAxis[pn]<<1)+1;
			if(!flag[index]) {
				tbounds[index] = fValue[pn];
				flag[index] = kTRUE;
				nvals++;
			}
		} else {
			index = fAxis[pn]<<1;
			if(!flag[index]) {
				tbounds[index] = fValue[pn];
				flag[index] = kTRUE;
				nvals++;
			}
		}
		LEFT = pn&1;
		pn =  (pn - 1)>>1;
	}
}



//______________________________________________________________________
TKDTreeIF *TKDTreeTestBuild(const Int_t npoints, const Int_t bsize){  
  //
  // Example function to  
  //
  Float_t *data0 =  new Float_t[npoints*2];
  Float_t *data[2];
  data[0] = &data0[0];
  data[1] = &data0[npoints];
  for (Int_t i=0;i<npoints;i++) {
    data[1][i]= gRandom->Rndm();
    data[0][i]= gRandom->Rndm();
  }
  TKDTree<Int_t, Float_t> *kdtree = new TKDTreeIF(npoints, 2, bsize, data);
  return kdtree;	
}



template class TKDTree<Int_t, Float_t>;
template class TKDTree<Int_t, Double_t>;

