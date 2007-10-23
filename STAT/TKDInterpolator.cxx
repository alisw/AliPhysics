#include "TKDInterpolator.h"

#include "TLinearFitter.h"
#include "TTree.h"
#include "TH2.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TBox.h"
#include "TGraph.h"
#include "TMarker.h"
#include "TVectorD.h"
#include "TMatrixD.h"

ClassImp(TKDInterpolator)
ClassImp(TKDInterpolator::TKDNodeInfo)

/////////////////////////////////////////////////////////////////////
// Memory setup of protected data members
// fRefPoints : evaluation point of PDF for each terminal node of underlying KD Tree.
// | 1st terminal node (fNDim point coordinates) | 2nd terminal node (fNDim point coordinates) | ...
//
// fRefValues : evaluation value/error of PDF for each terminal node of underlying KD Tree.
// | 1st terminal node (value) | 2nd terminal node (value) | ... | 1st terminal node (error) | 2nd terminal node (error) | ...
//
// status = |0|0|0|0|0|1(tri-cubic weights)|1(STORE)|1 INT(0 COG )|
/////////////////////////////////////////////////////////////////////

//_________________________________________________________________
TKDInterpolator::TKDNodeInfo::TKDNodeInfo(const Int_t dim): 
	fNDim(dim)
	,fRefPoint(0x0)
	,fRefValue(0.)
	,fCov(0x0)
	,fPar(0x0)
{
	if(fNDim) Build(dim);
}

//_________________________________________________________________
TKDInterpolator::TKDNodeInfo::~TKDNodeInfo()
{
	if(fRefPoint) delete [] fRefPoint;
	if(fCov){
		delete fPar;
		delete fCov;
	}
}

//_________________________________________________________________
void TKDInterpolator::TKDNodeInfo::Build(const Int_t dim)
{
// Allocate/Reallocate space for this node.

	if(!dim) return;

	fNDim = dim;
	Int_t lambda = Int_t(1 + fNDim + .5*fNDim*(fNDim+1));
	if(fRefPoint) delete [] fRefPoint;
	fRefPoint = new Float_t[fNDim];
	if(fCov){
		fCov->ResizeTo(lambda, lambda);
		fPar->ResizeTo(lambda);
	}
	return;
}

//_________________________________________________________________
void TKDInterpolator::TKDNodeInfo::Store(const TVectorD &par, const TMatrixD &cov)
{
	if(!fCov){
		fCov = new TMatrixD(cov.GetNrows(), cov.GetNrows());
		fPar = new TVectorD(par.GetNrows());
	}
	(*fPar) = par;
	(*fCov) = cov;
}

//_________________________________________________________________
Double_t TKDInterpolator::TKDNodeInfo::CookPDF(const Double_t *point, Double_t &result, Double_t &error)
{
// Recalculate the PDF for one node from the results of interpolation (parameters and covariance matrix)

	if(fNDim>10) return 0.; // support only up to 10 dimensions
	
	Int_t lambda = 1 + fNDim + (fNDim*(fNDim+1)>>1);
	Double_t fdfdp[66];
	Int_t ipar = 0;
	fdfdp[ipar++] = 1.;
	for(int idim=0; idim<fNDim; idim++){
		fdfdp[ipar++] = point[idim];
		for(int jdim=idim; jdim<fNDim; jdim++) fdfdp[ipar++] = point[idim]*point[jdim];
	}

	// calculate estimation
	result =0.; error = 0.;
	for(int i=0; i<lambda; i++){
		result += fdfdp[i]*(*fPar)(i);
		for(int j=0; j<lambda; j++) error += fdfdp[i]*fdfdp[j]*(*fCov)(i,j);
	}	
	error = TMath::Sqrt(error);
 	
	//printf("TKDNodeInfo::CookPDF() : %6.3f +- %6.3f\n", result, error);

	return 0.;
}


//_________________________________________________________________
TKDInterpolator::TKDInterpolator() : TKDTreeIF()
	,fNTNodes(0)
	,fTNodes(0x0)
	,fStatus(4)
	,fLambda(0)
	,fDepth(-1)
	,fAlpha(.5)
	,fRefPoints(0x0)
	,fBuffer(0x0)
	,fKDhelper(0x0)
	,fFitter(0x0)
{
// Default constructor. To be used with care since in this case building
// of data structure is completly left to the user responsability.
}

//_________________________________________________________________
TKDInterpolator::TKDInterpolator(Int_t npoints, Int_t ndim, UInt_t bsize, Float_t **data) : TKDTreeIF(npoints, ndim, bsize, data)
	,fNTNodes(GetNTNodes())
	,fTNodes(0x0)
	,fStatus(4)
	,fLambda(0)
	,fDepth(-1)
	,fAlpha(.5)
	,fRefPoints(0x0)
	,fBuffer(0x0)
	,fKDhelper(0x0)
	,fFitter(0x0)
{
// Wrapper constructor for the TKDTree.

	Build();
}


//_________________________________________________________________
TKDInterpolator::TKDInterpolator(TTree *t, const Char_t *var, const Char_t *cut, UInt_t bsize, Long64_t nentries, Long64_t firstentry) : TKDTreeIF()
	,fNTNodes(0)
	,fTNodes(0x0)
	,fStatus(4)
	,fLambda(0)
	,fDepth(-1)
	,fAlpha(.5)
	,fRefPoints(0x0)
	,fBuffer(0x0)
	,fKDhelper(0x0)
	,fFitter(0x0)
{
// Alocate data from a tree. The variables which have to be analysed are
// defined in the "var" parameter as a colon separated list. The format should
// be identical to that used by TTree::Draw().
//
// 

	TObjArray *vars = TString(var).Tokenize(":");
	fNDim = vars->GetEntriesFast(); fNDimm = 2*fNDim;
	if(fNDim > 6/*kDimMax*/) Warning("TKDInterpolator(TTree*, const Char_t, const Char_t, UInt_t)", Form("Variable number exceed maximum dimension %d. Results are unpredictable.", 6/*kDimMax*/));
	fBucketSize = bsize;

	Int_t np;
	Double_t *v;
	for(int idim=0; idim<fNDim; idim++){
		if(!(np = t->Draw(((TObjString*)(*vars)[idim])->GetName(), cut, "goff", nentries, firstentry))){
			Warning("TKDInterpolator(TTree*, const Char_t, const Char_t, UInt_t)", Form("Can not access data for keys %s. Key defined on tree :", ((TObjString*)(*vars)[idim])->GetName() ));
			TIterator *it = (t->GetListOfLeaves())->MakeIterator();
			TObject *o;
			while((o = (*it)())) printf("\t%s\n", o->GetName());
			continue;
		}
		if(!fNpoints){
			fNpoints = np;
			//Info("TKDInterpolator(TTree*, const Char_t, const Char_t, UInt_t)", Form("Allocating %d data points in %d dimensions.", fNpoints, fNDim));
			fData = new Float_t*[fNDim];
			for(int idim=0; idim<fNDim; idim++) fData[idim] = new Float_t[fNpoints];
			kDataOwner = kTRUE;
		}
		v = t->GetV1();
		for(int ip=0; ip<fNpoints; ip++) fData[idim][ip] = (Float_t)v[ip];
	}
	TKDTreeIF::Build();
	Build();
}

//_________________________________________________________________
TKDInterpolator::~TKDInterpolator()
{
	if(fFitter) delete fFitter;
	if(fKDhelper) delete fKDhelper;
	if(fBuffer) delete [] fBuffer;
	
	if(fRefPoints){
		for(int idim=0; idim<fNDim; idim++) delete [] fRefPoints[idim] ;
		delete [] fRefPoints;
	}
	if(fTNodes) delete [] fTNodes;
}

//_________________________________________________________________
void TKDInterpolator::Build()
{
// Fill interpolator's data array i.e.
//  - estimation points 
//  - corresponding PDF values

	fNTNodes = TKDTreeIF::GetNTNodes();
	if(!fBoundaries) MakeBoundaries();
	fLambda = 1 + fNDim + (fNDim*(fNDim+1)>>1);
	//printf("after MakeBoundaries() %d\n", memory());
	
	// allocate memory for data
	fTNodes = new TKDNodeInfo[fNTNodes];
	for(int in=0; in<fNTNodes; in++) fTNodes[in].Build(fNDim);
	//printf("after BuildNodes() %d\n", memory());

	Float_t *bounds = 0x0;
	Int_t *indexPoints;
	for(int inode=0, tnode = fNnodes; inode<fNTNodes-1; inode++, tnode++){
		fTNodes[inode].fRefValue =  Float_t(fBucketSize)/fNpoints;
		bounds = GetBoundary(tnode);
		for(int idim=0; idim<fNDim; idim++) fTNodes[inode].fRefValue /= (bounds[2*idim+1] - bounds[2*idim]);
		
		indexPoints = GetPointsIndexes(tnode);
		// loop points in this terminal node
		for(int idim=0; idim<fNDim; idim++){
			fTNodes[inode].fRefPoint[idim] = 0.;
			for(int ip = 0; ip<fBucketSize; ip++){
/*				printf("\t\tindex[%d] = %d %f\n", ip, indexPoints[ip], fData[idim][indexPoints[ip]]);*/
				fTNodes[inode].fRefPoint[idim] += fData[idim][indexPoints[ip]];
			}
			fTNodes[inode].fRefPoint[idim] /= fBucketSize;
		}
	}

	// analyze last (incomplete) terminal node
	Int_t counts = fNpoints%fBucketSize;
	counts = counts ? counts : fBucketSize;
	Int_t inode = fNTNodes - 1, tnode = inode + fNnodes;
	fTNodes[inode].fRefValue =  Float_t(counts)/fNpoints;
	bounds = GetBoundary(tnode);
	for(int idim=0; idim<fNDim; idim++) fTNodes[inode].fRefValue /= (bounds[2*idim+1] - bounds[2*idim]);

	// loop points in this terminal node
	indexPoints = GetPointsIndexes(tnode);
	for(int idim=0; idim<fNDim; idim++){
		fTNodes[inode].fRefPoint[idim] = 0.;
		for(int ip = 0; ip<counts; ip++) fTNodes[inode].fRefPoint[idim] += fData[idim][indexPoints[ip]];
		fTNodes[inode].fRefPoint[idim] /= counts;
	}
}

//__________________________________________________________________
void TKDInterpolator::GetStatus()
{
// Prints the status of the interpolator

	printf("Interpolator Status :\n");
	printf("  Method : %s\n", fStatus&1 ? "INT" : "COG");
	printf("  Store  : %s\n", fStatus&2 ? "YES" : "NO");
	printf("  Weights: %s\n", fStatus&4 ? "YES" : "NO");
	return;
	
	printf("fNTNodes %d\n", fNTNodes);        //Number of evaluation data points
	for(int i=0; i<fNTNodes; i++){
		printf("%d ", i);
		for(int idim=0; idim<fNDim; idim++) printf("%f ", fTNodes[i].fRefPoint[idim]);
		printf("[%f] %s\n", fTNodes[i].fRefValue, fTNodes[i].fCov ? "true" : "false");
		printf("Fit parameters : ");
		if(!fTNodes[i].fPar){
			printf("Not defined.\n");
			continue;
		}
		for(int ip=0; ip<3; ip++) printf("p%d[%f] ", ip, (*fTNodes[i].fPar)(ip));
		printf("\n");
	}
}

//_________________________________________________________________
Double_t TKDInterpolator::Eval(const Double_t *point, Double_t &result, Double_t &error, Bool_t force)
{
// Evaluate PDF for "point". The result is returned in "result" and error in "error". The function returns the chi2 of the fit.
//
// Observations:
//
// 1. The default method used for interpolation is kCOG.
// 2. The initial number of neighbors used for the estimation is set to Int(alpha*fLambda) (alpha = 1.5)
			
	Float_t pointF[50]; // local Float_t conversion for "point"
	for(int idim=0; idim<fNDim; idim++) pointF[idim] = (Float_t)point[idim];
	Int_t node = FindNode(pointF) - fNnodes;
	if((fStatus&1) && fTNodes[node].fCov && !force) return fTNodes[node].CookPDF(point, result, error);

	// Allocate memory
	if(!fBuffer) fBuffer = new Double_t[2*fLambda];
	if(!fKDhelper){ 
		fRefPoints = new Float_t*[fNDim];
		for(int id=0; id<fNDim; id++){ 
			fRefPoints[id] = new Float_t[fNTNodes];
			for(int in=0; in<fNTNodes; in++) fRefPoints[id][in] = fTNodes[in].fRefPoint[id];
		}
// 		for(int in=0; in<fNTNodes; in++){
// 			printf("%3d ", in);
// 			for(int id=0; id<fNDim; id++) printf("%6.3f ", fTNodes[in].fRefPoint[id]/*fRefPoints[id][in]*/);
// 			printf("\n");
// 		}
		fKDhelper = new TKDTreeIF(fNTNodes, fNDim, 30, fRefPoints);
		fKDhelper->MakeBoundaries();
	}
	if(!fFitter) fFitter = new TLinearFitter(fLambda, Form("hyp%d", fLambda-1));
	
	// generate parabolic for nD
	//Float_t alpha = Float_t(2*lambda + 1) / fNTNodes; // the bandwidth or smoothing parameter
	//Int_t npoints = Int_t(alpha * fNTNodes);
	//printf("Params : %d NPoints %d\n", lambda, npoints);
	// prepare workers

	Int_t *index,  // indexes of NN 
	      ipar,    // local looping variable
				npoints = Int_t((1.+fAlpha)*fLambda); // number of data points used for interpolation
	Float_t *dist, // distances of NN
					d,     // NN normalized distance
					w0,    // work
					w;     // tri-cubic weight function
	Double_t sig   // bucket error 
	        = TMath::Sqrt(1./fBucketSize);

	do{
		// find nearest neighbors
		for(int idim=0; idim<fNDim; idim++) pointF[idim] = (Float_t)point[idim];
		if(!fKDhelper->FindNearestNeighbors(pointF, npoints+1, index, dist)){
			Error("Eval()", Form("Failed retriving %d neighbours for point:", npoints));
			for(int idim=0; idim<fNDim; idim++) printf("%f ", point[idim]);
			printf("\n");
			return -1;
		}
		// add points to fitter
		fFitter->ClearPoints();
		TKDNodeInfo *node = 0x0;
		for(int in=0; in<npoints; in++){
			node = &fTNodes[index[in]];
			if(fStatus&1){ // INT
				Float_t *bounds = GetBoundary(FindNode(node->fRefPoint));
				ipar = 0;
				for(int idim=0; idim<fNDim; idim++){
					fBuffer[ipar++] = .5*(bounds[2*idim] + bounds[2*idim+1]);
					fBuffer[ipar++] = (bounds[2*idim]*bounds[2*idim] + bounds[2*idim] * bounds[2*idim+1] + bounds[2*idim+1] * bounds[2*idim+1])/3.;
					for(int jdim=idim+1; jdim<fNDim; jdim++) fBuffer[ipar++] = (bounds[2*idim] + bounds[2*idim+1]) * (bounds[2*jdim] + bounds[2*jdim+1]) * .25;
				}
			} else { // COG
				Float_t *p = node->fRefPoint;
				ipar = 0;
				for(int idim=0; idim<fNDim; idim++){
					fBuffer[ipar++] = p[idim];
					for(int jdim=idim; jdim<fNDim; jdim++) fBuffer[ipar++] = p[idim]*p[jdim];
				}
			}

			// calculate tri-cubic weighting function
			if(fStatus&4){
				d = dist[in]/ dist[npoints];
				w0 = (1. - d*d*d); w = w0*w0*w0;
			} else w = 1.;
			 
			//for(int idim=0; idim<fNDim; idim++) printf("%f ", fBuffer[idim]);
			//printf("\nd[%f] w[%f] sig[%f]\n", d, w, sig);
			fFitter->AddPoint(fBuffer, node->fRefValue, node->fRefValue*sig/w);
		}
		npoints += 4;
	} while(fFitter->Eval());

	// retrive fitter results
	TMatrixD cov(fLambda, fLambda);
	TVectorD par(fLambda);
	fFitter->GetCovarianceMatrix(cov);
	fFitter->GetParameters(par);
	Double_t chi2 = fFitter->GetChisquare()/(npoints - 4 - fLambda);

	// store results
	if(fStatus&2 && fStatus&1) fTNodes[node].Store(par, cov);
		
	// Build df/dpi|x values
	Double_t *fdfdp = &fBuffer[fLambda];
	ipar = 0;
	fdfdp[ipar++] = 1.;
	for(int idim=0; idim<fNDim; idim++){
		fdfdp[ipar++] = point[idim];
		for(int jdim=idim; jdim<fNDim; jdim++) fdfdp[ipar++] = point[idim]*point[jdim];
	}

	// calculate estimation
	result =0.; error = 0.;
	for(int i=0; i<fLambda; i++){
		result += fdfdp[i]*par(i);
		for(int j=0; j<fLambda; j++) error += fdfdp[i]*fdfdp[j]*cov(i,j);
	}	
	error = TMath::Sqrt(error);

	return chi2;
}

//_________________________________________________________________
void TKDInterpolator::DrawNodes(UInt_t ax1, UInt_t ax2, Int_t depth)
{
// Draw nodes structure projected on plane "ax1:ax2". The parameter
// "depth" specifies the bucket size per node. If depth == -1 draw only
// terminal nodes and evaluation points (default -1 i.e. bucket size per node equal bucket size specified by the user)
//
// Observation:
// This function creates the nodes (TBox) array for the specified depth
// but don't delete it. Abusing this function may cause memory leaks !


	if(!fBoundaries) MakeBoundaries();

	// Count nodes in specific view
	Int_t nnodes = 0;
	for(int inode = 0; inode <= 2*fNnodes; inode++){
		if(depth == -1){
			if(!IsTerminal(inode)) continue;
		} else if((inode+1) >> depth != 1) continue;
		nnodes++;
	}

	//printf("depth %d nodes %d\n", depth, nnodes);
	
	TH2 *h2 = new TH2S("hNodes", "", 100, fRange[2*ax1], fRange[2*ax1+1], 100, fRange[2*ax2], fRange[2*ax2+1]);
	h2->GetXaxis()->SetTitle(Form("x_{%d}", ax1));
	h2->GetYaxis()->SetTitle(Form("x_{%d}", ax2));
	h2->Draw();
	
	const Float_t kBorder = 0.;//1.E-4;
	TBox *nodeArray = new TBox[nnodes], *node;
	Float_t *bounds = 0x0;
	nnodes = 0;
	for(int inode = 0; inode <= 2*fNnodes; inode++){
		if(depth == -1){
			if(!IsTerminal(inode)) continue;
		} else if((inode+1) >> depth != 1) continue;

		node = &nodeArray[nnodes++];
		//node = new TBox(bounds[2*ax1]+border, bounds[2*ax2]+border, bounds[2*ax1+1]-border, bounds[2*ax2+1]-border);
		node->SetFillStyle(3002);	
		node->SetFillColor(50+inode/*Int_t(gRandom->Uniform()*50.)*/);
		bounds = GetBoundary(inode);
		node->DrawBox(bounds[2*ax1]+kBorder, bounds[2*ax2]+kBorder, bounds[2*ax1+1]-kBorder, bounds[2*ax2+1]-kBorder);
	}
	if(depth != -1) return;

	// Draw reference points
	TGraph *ref = new TGraph(fNTNodes);
	ref->SetMarkerStyle(3);
	ref->SetMarkerSize(.7);
	ref->SetMarkerColor(2);
	for(int inode = 0; inode < fNTNodes; inode++) ref->SetPoint(inode, fTNodes[inode].fRefPoint[ax1], fTNodes[inode].fRefPoint[ax2]);
	ref->Draw("p");
	return;
}

//_________________________________________________________________
void TKDInterpolator::DrawNode(Int_t tnode, UInt_t ax1, UInt_t ax2)
{
// Draw node "node" and the data points within.
//
// Observation:
// This function creates some graphical objects
// but don't delete it. Abusing this function may cause memory leaks !

	if(tnode < 0 || tnode >= fNTNodes){
		Warning("DrawNode()", Form("Terminal node %d outside defined range.", tnode));
		return;
	}

	Int_t inode = tnode;
	tnode += fNnodes;
	// select zone of interest in the indexes array
	Int_t *index = GetPointsIndexes(tnode);
	Int_t nPoints = (tnode == 2*fNnodes) ? fNpoints%fBucketSize : fBucketSize;

	// draw data points
	TGraph *g = new TGraph(nPoints);
	g->SetMarkerStyle(7);
	for(int ip = 0; ip<nPoints; ip++) g->SetPoint(ip, fData[ax1][index[ip]], fData[ax2][index[ip]]);

	// draw estimation point
	TMarker *m=new TMarker(fTNodes[inode].fRefPoint[ax1], fTNodes[inode].fRefPoint[ax2], 20);
	m->SetMarkerColor(2);
	m->SetMarkerSize(1.7);
	
	// draw node contour
	Float_t *bounds = GetBoundary(tnode);
	TBox *n = new TBox(bounds[2*ax1], bounds[2*ax2], bounds[2*ax1+1], bounds[2*ax2+1]);
	n->SetFillStyle(0);

	g->Draw("ap");
	m->Draw();
	n->Draw();
	
	return;
}


//__________________________________________________________________
void TKDInterpolator::SetInterpolationMethod(const Bool_t on)
{
// Set interpolation bit to "on".
	
	if(on) fStatus += fStatus&1 ? 0 : 1;
	else fStatus += fStatus&1 ? -1 : 0;
}


//_________________________________________________________________
void TKDInterpolator::SetStore(const Bool_t on)
{
// Set store bit to "on"
	
	if(on) fStatus += fStatus&2 ? 0 : 2;
	else fStatus += fStatus&2 ? -2 : 0;
}

//_________________________________________________________________
void TKDInterpolator::SetWeights(const Bool_t on)
{
// Set weights bit to "on"
	
	if(on) fStatus += fStatus&4 ? 0 : 4;
	else fStatus += fStatus&4 ? -4 : 0;
}
