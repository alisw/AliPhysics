#include "TKDInterpolator.h"

#include "TLinearFitter.h"
#include "TVector.h"
#include "TTree.h"
#include "TH2.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TPad.h"
#include "TBox.h"
#include "TGraph.h"
#include "TMarker.h"
#include "TRandom.h"
#include "TROOT.h"

ClassImp(TKDInterpolator)
ClassImp(TKDInterpolator::TKDNodeInfo)

/////////////////////////////////////////////////////////////////////
// Memory setup of protected data memebers
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
	,fCov()
	,fPar()
	,fPDFstatus(kFALSE)
{
	if(fNDim) Build(dim);
}

//_________________________________________________________________
TKDInterpolator::TKDNodeInfo::~TKDNodeInfo()
{
	if(fRefPoint) delete [] fRefPoint;
}

//_________________________________________________________________
void TKDInterpolator::TKDNodeInfo::Build(const Int_t dim)
{
	if(!dim) return;

	fNDim = dim;
	Int_t lambda = Int_t(1 + fNDim + .5*fNDim*(fNDim+1));
	if(fRefPoint) delete [] fRefPoint;
	fRefPoint = new Float_t[fNDim];
	fCov.ResizeTo(lambda, lambda);
	fPar.ResizeTo(lambda);
	return;
}


//_________________________________________________________________
TKDInterpolator::TKDInterpolator() : TKDTreeIF()
	,fNTNodes(0)
	,fTNodes(0x0)
	,fStatus(4)
	,fLambda(0)
	,fDepth(-1)
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
	,fNTNodes(GetNTerminalNodes())
	,fTNodes(0x0)
	,fStatus(4)
	,fLambda(0)
	,fDepth(-1)
	,fRefPoints(0x0)
	,fBuffer(0x0)
	,fKDhelper(0x0)
	,fFitter(0x0)
{
// Wrapper constructor for the similar TKDTree one.
	
	Build();
}


//_________________________________________________________________
TKDInterpolator::TKDInterpolator(TTree *t, const Char_t *var, const Char_t *cut, UInt_t bsize, Long64_t nentries, Long64_t firstentry) : TKDTreeIF()
	,fNTNodes(0)
	,fTNodes(0x0)
	,fStatus(4)
	,fLambda(0)
	,fDepth(-1)
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
			while(o = (*it)()) printf("\t%s\n", o->GetName());
			continue;
		}
		if(!fNpoints){
			fNpoints = np;
			Info("TKDInterpolator(TTree*, const Char_t, const Char_t, UInt_t)", Form("Allocating %d data points in %d dimensions.", fNpoints, fNDim));
			fData = new Float_t*[fNDim];
			for(int idim=0; idim<fNDim; idim++) fData[idim] = new Float_t[fNpoints];
			kDataOwner = kTRUE;
		}
		v = t->GetV1();
		for(int ip=0; ip<fNpoints; ip++) fData[idim][ip] = (Float_t)v[ip];
	}
	TKDTreeIF::Build();
	fNTNodes = GetNTerminalNodes();
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

	if(!fBoundaries) MakeBoundaries();
	fLambda = 1 + fNDim + fNDim*(fNDim+1)/2;

	// allocate memory for data
	fTNodes = new TKDNodeInfo[fNTNodes];
	for(int in=0; in<fNTNodes; in++) fTNodes[in].Build(fNDim);

	Float_t *bounds = 0x0;
	Int_t *indexPoints;
	for(int inode=0, tnode = fNnodes; inode<fNTNodes-1; inode++, tnode++){
		fTNodes[inode].fRefValue =  Float_t(fBucketSize)/fNpoints;
		bounds = GetBoundary(tnode);
		for(int idim=0; idim<fNDim; idim++) fTNodes[inode].fRefValue /= (bounds[2*idim+1] - bounds[2*idim]);

		indexPoints = GetPointsIndexes(tnode);
		// loop points in this terminal node
		for(int idim=0; idim<fNDim; idim++){
			for(int ip = 0; ip<fBucketSize; ip++) fTNodes[inode].fRefPoint[idim] += fData[idim][indexPoints[ip]];
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

	indexPoints = GetPointsIndexes(tnode);
	// loop points in this terminal node
	for(int idim=0; idim<fNDim; idim++){
		for(int ip = 0; ip<counts; ip++) fTNodes[inode].fRefPoint[idim] += fData[idim][indexPoints[ip]];
		fTNodes[inode].fRefPoint[idim] /= counts;
	}

	//GetStatus();
}

//__________________________________________________________________
void TKDInterpolator::GetStatus()
{
	printf("Interpolator Status :\n");
	printf("  Method : %s\n", fStatus&1 ? "INT" : "COG");
	printf("  Store  : %s\n", fStatus&2 ? "YES" : "NO");
	printf("  Weights: %s\n", fStatus&4 ? "YES" : "NO");

	printf("nnodes %d\n", fNTNodes);        //Number of evaluation data points
	printf("nodes 0x%x\n", fTNodes);    //[fNTNodes]
	for(int i=0; i<fNTNodes; i++){
		printf("\t%d ", i);
		for(int idim=0; idim<fNDim; idim++) printf("%f ", fTNodes[i].fRefPoint[idim]);
		printf("[%f] %s\n", fTNodes[i].fRefValue, fTNodes[i].fPDFstatus ? "true" : "false");
		for(int ip=0; ip<3; ip++) printf("p%d[%f] ", ip, fTNodes[i].fPar(ip));
		printf("\n");
	}
}

//_________________________________________________________________
Double_t TKDInterpolator::Eval(const Double_t *point, Double_t &result, Double_t &error)
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
	if((fStatus&1) && fTNodes[node].fPDFstatus) return CookPDF(point, node, result, error); // maybe move to TKDNodeInfo

	// Allocate memory
	if(!fBuffer) fBuffer = new Double_t[2*fLambda];
	if(!fKDhelper){ 
		fRefPoints = new Float_t*[fNDim];
		for(int id=0; id<fNDim; id++){ 
			fRefPoints[id] = new Float_t[fNTNodes];
			for(int in=0; in<fNTNodes; in++) fRefPoints[id][in] = fTNodes[in].fRefPoint[id];
		}
		fKDhelper = new TKDTreeIF(fNTNodes, fNDim, 30, fRefPoints);
	}
	if(!fFitter) SetIntInterpolation(kFALSE);
	
	// generate parabolic for nD
	//Float_t alpha = Float_t(2*lambda + 1) / fNTNodes; // the bandwidth or smoothing parameter
	//Int_t npoints = Int_t(alpha * fNTNodes);
	//printf("Params : %d NPoints %d\n", lambda, npoints);
	// prepare workers

	Int_t *index,  // indexes of NN 
	      ipar,    // local looping variable
				npoints = Int_t(1.5*fLambda); // number of data points used for interpolation
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
		for(int in=0; in<npoints; in++){
			if(fStatus&1){ // INT
				//for(int idim=0; idim<fNDim; idim++) pointF[idim] = fRefPoints[idim][index[in]];
				Float_t *bounds = GetBoundary(FindNode(fTNodes[index[in]].fRefPoint/*pointF*/));
				
				ipar = 0;
				for(int idim=0; idim<fNDim; idim++){
					fBuffer[ipar++] = .5*(bounds[2*idim] + bounds[2*idim+1]);
					fBuffer[ipar++] = (bounds[2*idim]*bounds[2*idim] + bounds[2*idim] * bounds[2*idim+1] + bounds[2*idim+1] * bounds[2*idim+1])/3.;
					for(int jdim=idim+1; jdim<fNDim; jdim++) fBuffer[ipar++] = (bounds[2*idim] + bounds[2*idim+1]) * (bounds[2*jdim] + bounds[2*jdim+1]) * .25;
				}
			} else { // COG
				for(int idim=0; idim<fNDim; idim++) fBuffer[idim] = fTNodes[index[in]].fRefPoint[idim];
			}

			// calculate tri-cubic weighting function
			if(fStatus&4){
				d = dist[in]/ dist[npoints];
				w0 = (1. - d*d*d); w = w0*w0*w0;
			} else w = 1.;
			 
			//for(int idim=0; idim<fNDim; idim++) printf("%f ", fBuffer[idim]);
			//printf("\nd[%f] w[%f] sig[%f]\n", d, w, sig);
			fFitter->AddPoint(fBuffer, fTNodes[index[in]].fRefValue, fTNodes[index[in]].fRefValue*sig/w);
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
	if(fStatus&2 && fStatus&1){
		fTNodes[node].fPar = par;
		fTNodes[node].fCov = cov;
		fTNodes[node].fPDFstatus = kTRUE;
	}
		
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

// //_________________________________________________________________
// Double_t TKDInterpolator::Eval1(const Double_t *point, Int_t npoints, Double_t &result, Double_t &error)
// {
// // Evaluate PDF at k-dimensional position "point". The initial number of
// // neighbour estimation points is set to "npoints". The default method
// // used for interpolation is kCOG.
// 
// 	// calculate number of parameters in the parabolic expresion
// 	Int_t lambda = 1 + fNDim + fNDim*(fNDim+1)/2;
// 
// 	if(!fBuffer) fBuffer = new Double_t[lambda-1];
// 	if(!fKDhelper) fKDhelper = new TKDTreeIF(GetNTerminalNodes(), fNDim, npoints, fRefPoints);
// 
// 	if(!fFitter) fFitter = new TLinearFitter(lambda, Form("hyp%d", fNDim+1));
// 	else fFitter->SetFormula(Form("hyp%d", fNDim+1));
// 
// 
// 	Float_t pointF[50];
// 	for(int idim=0; idim<fNDim; idim++) pointF[idim] = point[idim];
// 	Int_t istart = 0;
// 	Int_t *index, ipar;
// 	Float_t *bounds, *dist, *w = new Float_t[fNDim];
// 	Double_t uncertainty = TMath::Sqrt(1./fBucketSize);
// 	fFitter->ClearPoints();
// 	do{
// 		if(!fKDhelper->FindNearestNeighbors(pointF, npoints+1, index, dist)){
// 			Error("Eval()", Form("Failed retriving %d neighbours for point:", npoints));
// 			for(int idim=0; idim<fNDim; idim++) printf("%f ", point[idim]);
// 			printf("\n");
// 			return -1;
// 		}
// 		for(int in=istart; in<npoints; in++){
// 			for(int idim=0; idim<fNDim; idim++) w[idim] = fRefPoints[idim][index[in]];
// 			bounds = GetBoundary(FindNode(w));
// 
// 			ipar = 0;
// 			for(int idim=0; idim<fNDim; idim++){
// 				fBuffer[ipar++] = .5*(bounds[2*idim] + bounds[2*idim+1]);
// 				fBuffer[ipar++] = (bounds[2*idim]*bounds[2*idim] + bounds[2*idim] * bounds[2*idim+1] + bounds[2*idim+1] * bounds[2*idim+1])/3.;
// 				for(int jdim=idim+1; jdim<fNDim; jdim++) fBuffer[ipar++] = (bounds[2*idim] + bounds[2*idim+1]) * (bounds[2*jdim] + bounds[2*jdim+1]) * .25;
// 			}
// 
// 			fFitter->AddPoint(fBuffer, fRefValues[index[in]], fRefValues[index[in]]*uncertainty);
// 		}
// 		istart = npoints;
// 		npoints += 4;
// 	} while(fFitter->Eval());
// 	delete [] w;
// 
// 	// calculate evaluation
// //	fFitter->PrintResults(3);
// 	TMatrixD cov(lambda, lambda);
// 	TVectorD par(lambda);
// 	fFitter->GetCovarianceMatrix(cov);
// 	fFitter->GetParameters(par);
// 
// 	// Build temporary array to keep values df/dpi|x
// 	Double_t f[100];
// 	ipar = 0;
// 	f[ipar++] = 1.;
// 	for(int idim=0; idim<fNDim; idim++){
// 		f[ipar++] = point[idim];
// 		for(int jdim=idim; jdim<fNDim; jdim++) f[ipar++] = point[idim]*point[jdim];
// 	}
// 	result =0.; error = 0.;
// 	for(int i=0; i<lambda; i++){
// 		result += f[i]*par[i];
// 		for(int j=0; j<lambda; j++) error += f[i]*f[j]*cov(i,j);
// 	}
// 	error = TMath::Sqrt(error);
// 	Double_t chi2 = fFitter->GetChisquare()/(npoints - 4 - lambda);
// 
//  	for(int ipar=0; ipar<lambda; ipar++) printf("%d %8.6e %8.6e\n", ipar, par[ipar], TMath::Sqrt(cov(ipar, ipar)));
//  	printf("result %6.3f +- %6.3f [%f]\n", result, error, chi2);
// 	return chi2;
// }


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
	
	TH2 *h2 = 0x0;
	if(!(h2 = (TH2S*)gROOT->FindObject("hNodes"))) h2 = new TH2S("hNodes", "", 100, fRange[2*ax1], fRange[2*ax1+1], 100, fRange[2*ax2], fRange[2*ax2+1]);
	h2->GetXaxis()->SetTitle(Form("x_{%d}", ax1));
	h2->GetYaxis()->SetTitle(Form("x_{%d}", ax2));
	h2->Draw();
	
	const Float_t border = 0.;//1.E-4;
	TBox *node_array = new TBox[nnodes], *node;
	Float_t *bounds = 0x0;
	nnodes = 0;
	for(int inode = 0; inode <= 2*fNnodes; inode++){
		if(depth == -1){
			if(!IsTerminal(inode)) continue;
		} else if((inode+1) >> depth != 1) continue;

		node = &node_array[nnodes++];
		//node = new TBox(bounds[2*ax1]+border, bounds[2*ax2]+border, bounds[2*ax1+1]-border, bounds[2*ax2+1]-border);
		node->SetFillStyle(3002);	
		node->SetFillColor(50+Int_t(gRandom->Uniform()*50.));
		bounds = GetBoundary(inode);
		node->DrawBox(bounds[2*ax1]+border, bounds[2*ax2]+border, bounds[2*ax1+1]-border, bounds[2*ax2+1]-border);
	}
	if(depth != -1) return;

	// Draw reference points
	TGraph *ref = new TGraph(GetNTerminalNodes());
	ref->SetMarkerStyle(3);
	ref->SetMarkerSize(.7);
	ref->SetMarkerColor(2);
	for(int inode = 0; inode < GetNTerminalNodes(); inode++) ref->SetPoint(inode, fTNodes[inode].fRefPoint[ax1], fTNodes[inode].fRefPoint[ax2]);
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

	if(tnode < 0 || tnode >= GetNTerminalNodes()){
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

	if(gPad) gPad->Clear();	
	g->Draw("ap");
	m->Draw();
	n->Draw();
	
	return;
}


//__________________________________________________________________
void TKDInterpolator::SetIntInterpolation(const Bool_t on)
{
// Set interpolation bit to "on" and build/delete memory
	
	if(on) fStatus += fStatus&1 ? 0 : 1;
	else fStatus += fStatus&1 ? -1 : 0;
	TString formula;
	if(on) formula = Form("hyp%d", fLambda-1);
	else {
		formula = "1";
		for(int idim=0; idim<fNDim; idim++){
			formula += Form("++x[%d]", idim);
			for(int jdim=idim; jdim<fNDim; jdim++) formula += Form("++x[%d]*x[%d]", idim, jdim);
		}
	}
	if(!fFitter) fFitter = new TLinearFitter(fLambda, formula.Data());
	else fFitter->SetFormula(formula.Data());
}


//_________________________________________________________________
void TKDInterpolator::SetSetStore(const Bool_t on)
{
// Set store bit to "on" and build/delete memory
	
	if(on){
		fStatus += fStatus&2 ? 0 : 2;
/*		if(!fCov){
			fPDFstatus = new Bool_t[fNTNodes];
			fCov = new TMatrixD[fNTNodes];
			fPar = new TVectorD[fNTNodes];
			for(int i=0; i<fNTNodes; i++){
				fPDFstatus[i] = kFALSE;
				fCov[i].ResizeTo(fLambda, fLambda);
				fPar[i].ResizeTo(fLambda);
			}
		}*/
	} else {
		fStatus += fStatus&2 ? -2 : 0;
/*		if(fCov){
			delete [] fPar;
			delete [] fCov;
			delete [] fPDFstatus;
		}*/
	}
}

//_________________________________________________________________
void TKDInterpolator::SetUseWeights(const Bool_t on)
{
	if(on) fStatus += fStatus&4 ? 0 : 4;
	else fStatus += fStatus&4 ? -4 : 0;
}


//_________________________________________________________________
Double_t TKDInterpolator::CookPDF(const Double_t *point, const Int_t node, Double_t &result, Double_t &error)
{
// Recalculate the PDF for one node from the results of interpolation (parameters and covariance matrix)

	Info("CookPDF()", Form("Called for node %d", node));

	if(!fBuffer) fBuffer = new Double_t[2*fLambda];
	Double_t *fdfdp = &fBuffer[fLambda];
	Int_t ipar = 0;
	fdfdp[ipar++] = 1.;
	for(int idim=0; idim<fNDim; idim++){
		fdfdp[ipar++] = point[idim];
		for(int jdim=idim; jdim<fNDim; jdim++) fdfdp[ipar++] = point[idim]*point[jdim];
	}

	// calculate estimation
	result =0.; error = 0.;
	for(int i=0; i<fLambda; i++){
		result += fdfdp[i]*fTNodes[node].fPar(i);
		for(int j=0; j<fLambda; j++) error += fdfdp[i]*fdfdp[j]*fTNodes[node].fCov(i,j);
	}	
	error = TMath::Sqrt(error);
 	printf("result[CookPDF] %6.3f +- %6.3f\n", result, error);

	return 0.;
}

