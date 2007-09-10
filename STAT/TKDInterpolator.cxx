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

/////////////////////////////////////////////////////////////////////
// Memory setup of protected data memebers
// fRefPoints : evaluation point of PDF for each terminal node of underlying KD Tree.
// | 1st terminal node (fNDim point coordinates) | 2nd terminal node (fNDim point coordinates) | ...
//
// fRefValues : evaluation value/error of PDF for each terminal node of underlying KD Tree.
// | 1st terminal node (value) | 2nd terminal node (value) | ... | 1st terminal node (error) | 2nd terminal node (error) | ...
/////////////////////////////////////////////////////////////////////

//_________________________________________________________________
TKDInterpolator::TKDInterpolator() : TKDTreeIF()
	,fNTNodes(0)
	,fRefPoints(0x0)
	,fRefValues(0x0)
	,fDepth(-1)
	,fTmpPoint(0x0)
	,fKDhelper(0x0)
	,fFitter(0x0)
{
// Default constructor. To be used with care since in this case building
// of data structure is completly left to the user responsability.
}

//_________________________________________________________________
TKDInterpolator::TKDInterpolator(Int_t npoints, Int_t ndim, UInt_t bsize, Float_t **data) : TKDTreeIF(npoints, ndim, bsize, data)
	,fNTNodes(GetNTerminalNodes())
	,fRefPoints(0x0)
	,fRefValues(0x0)
	,fDepth(-1)
	,fTmpPoint(0x0)
	,fKDhelper(0x0)
	,fFitter(0x0)
{
// Wrapper constructor for the similar TKDTree one.
	
	Build();
}


//_________________________________________________________________
TKDInterpolator::TKDInterpolator(TTree *t, const Char_t *var, const Char_t *cut, UInt_t bsize) : TKDTreeIF()
	,fNTNodes(0)
	,fRefPoints(0x0)
	,fRefValues(0x0)
	,fDepth(-1)
	,fTmpPoint(0x0)
	,fKDhelper(0x0)
	,fFitter(0x0)
{
// Alocate data from a tree. The variables which have to be analysed are
// defined in the "var" parameter as a colon separated list. The format should
// be identical to that used by TTree::Draw().
//
// 

	TObjArray *vars = TString(var).Tokenize(":");
	fNDim = vars->GetEntriesFast();
	if(fNDim > 6/*kDimMax*/) Warning("TKDInterpolator(TTree*, const Char_t, const Char_t, UInt_t)", Form("Variable number exceed maximum dimension %d. Results are unpredictable.", 6/*kDimMax*/));
	fBucketSize = bsize;

	Int_t np;
	Double_t *v;
	for(int idim=0; idim<fNDim; idim++){
		if(!(np = t->Draw(((TObjString*)(*vars)[idim])->GetName(), cut, "goff"))){
			Warning("TKDInterpolator(TTree*, const Char_t, const Char_t, UInt_t)", Form("Can not access data for %s", ((TObjString*)(*vars)[idim])->GetName() ));
			continue;
		}
		if(!fNpoints){
			fNpoints = np;
			Info("TKDInterpolator(TTree*, const Char_t, const Char_t, UInt_t)", Form("Allocating %d data points in %d dimensions.", fNpoints, fNDim));
			//Float_t *mem = new Float_t[fNDim*fNpoints];
			fData = new Float_t*[fNDim];
			for(int idim=0; idim<fNDim; idim++) fData[idim] = new Float_t[fNpoints]; //&mem[idim*fNpoints];
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
	if(fTmpPoint) delete [] fTmpPoint;
	
	if(fRefPoints){
		for(int idim=0; idim<fNDim; idim++) delete [] fRefPoints[idim] ;
		delete [] fRefPoints;
	}
	if(fRefValues) delete [] fRefValues;
}

//_________________________________________________________________
void TKDInterpolator::Build()
{
// Fill interpolator's data array i.e.
//  - estimation points 
//  - corresponding PDF values

	if(!fBoundaries) MakeBoundaries();
	
	// allocate memory for data
	fRefValues = new Float_t[fNTNodes];
	fRefPoints = new Float_t*[fNDim];
	for(int id=0; id<fNDim; id++){
		fRefPoints[id] = new Float_t[fNTNodes];
		for(int in=0; in<fNTNodes; in++) fRefPoints[id][in] = 0.;
	}

	Float_t *bounds = 0x0;
	Int_t *indexPoints;
	for(int inode=0, tnode = fNnodes; inode<fNTNodes-1; inode++, tnode++){
		fRefValues[inode] =  Float_t(fBucketSize)/fNpoints;
		bounds = GetBoundary(tnode);
		for(int idim=0; idim<fNDim; idim++) fRefValues[inode] /= (bounds[2*idim+1] - bounds[2*idim]);

		indexPoints = GetPointsIndexes(tnode);
		// loop points in this terminal node
		for(int idim=0; idim<fNDim; idim++){
			for(int ip = 0; ip<fBucketSize; ip++) fRefPoints[idim][inode] += fData[idim][indexPoints[ip]];
			fRefPoints[idim][inode] /= fBucketSize;
		}
	}

	// analyze last (incomplete) terminal node
	Int_t counts = fNpoints%fBucketSize;
	counts = counts ? counts : fBucketSize;
	Int_t inode = fNTNodes - 1, tnode = inode + fNnodes;
	fRefValues[inode] =  Float_t(counts)/fNpoints;
	bounds = GetBoundary(tnode);
	for(int idim=0; idim<fNDim; idim++) fRefValues[inode] /= (bounds[2*idim+1] - bounds[2*idim]);

	indexPoints = GetPointsIndexes(tnode);
	// loop points in this terminal node
	for(int idim=0; idim<fNDim; idim++){
		for(int ip = 0; ip<counts; ip++) fRefPoints[idim][inode] += fData[idim][indexPoints[ip]];
		fRefPoints[idim][inode] /= counts;
	}
}

//_________________________________________________________________
Double_t TKDInterpolator::Eval(const Double_t *point, Int_t npoints)
{
// Evaluate PDF at k-dimensional position "point". The initial number of
// neighbour estimation points is set to "npoints"
	
	//Int_t npoints = Int_t(alpha * fNTNodes);
	//printf("Params : %d NPoints %d\n", lambda, npoints);
	// prepare workers
	if(!fTmpPoint) fTmpPoint = new Double_t[fNDim];
	if(!fKDhelper) fKDhelper = new TKDTreeIF(GetNTerminalNodes(), fNDim, npoints, fRefPoints);
	if(!fFitter){
		// generate parabolic for nD
		
		// calculate number of parameters in the parabolic expresion
		Int_t lambda = 1 + fNDim + fNDim*(fNDim+1)/2;
		//Float_t alpha = Float_t(2*lambda + 1) / fNTNodes; // the bandwidth or smoothing parameter
		TString formula("1");
		for(int idim=0; idim<fNDim; idim++){
			formula += Form("++x[%d]", idim);
			for(int jdim=idim; jdim<fNDim; jdim++) formula += Form("++x[%d]*x[%d]", idim, jdim);
		}
		fFitter = new TLinearFitter(lambda, formula.Data());
		Info("Eval(const Double_t*, Int_t)", Form("Using %s for local interpolation.", formula.Data()));
	}

	Float_t pointF[50];
	for(int idim=0; idim<fNDim; idim++) pointF[idim] = point[idim];
	Int_t istart = 0;
	Int_t *index;
	Float_t dist, d0, w0, w;
	Double_t uncertainty = TMath::Sqrt(1./fBucketSize); 
	fFitter->ClearPoints();
	do{
		if(!fKDhelper->FindNearestNeighbors(pointF, npoints+1, index, dist)){
			Error("Eval()", Form("Failed retriving %d neighbours for point:", npoints));
			for(int idim=0; idim<fNDim; idim++) printf("%f ", point[idim]);
			printf("\n");
			return -1;
		}
		for(int in=istart; in<npoints; in++){
			//printf("%d index[%2d] x(", in, index[in]);
			d0 = 0.;
			for(int idim=0; idim<fNDim; idim++){
				fTmpPoint[idim] = fRefPoints[idim][index[in]];
				//printf("%6.4f ", fTmpPoint[idim]);
				d0 += TMath::Abs(fTmpPoint[idim] - point[idim]);
			}
			d0 /= dist;
			w0 = (1. - d0*d0*d0);
			w = w0*w0*w0;

			//printf(") f = %f [%f] d0 = %6.4f   w = %6.4f  \n", fRefValues[index[in]], TMath::Log(fRefValues[index[in]]), d0, w);
			fFitter->AddPoint(fTmpPoint, TMath::Log(fRefValues[index[in]]), uncertainty/w);
		}
		istart = npoints;
		npoints += 4;
	} while(fFitter->Eval());

	// calculate evaluation
	Int_t ipar = 0;
	Double_t result = fFitter->GetParameter(ipar++);
	for(int idim=0; idim<fNDim; idim++){
		result += fFitter->GetParameter(ipar++)*point[idim];
		for(int jdim=idim; jdim<fNDim; jdim++) result += fFitter->GetParameter(ipar++)*point[idim]*point[jdim];
	}
	//printf("\tResult : %f [%f]\n", TMath::Exp(result), result);
	return TMath::Exp(result);
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
	for(int inode = 0; inode < GetNTerminalNodes(); inode++) ref->SetPoint(inode, fRefPoints[ax1][inode], fRefPoints[ax2][inode]);
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
	TMarker *m=new TMarker(fRefPoints[ax1][inode], fRefPoints[ax2][inode], 20);
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

