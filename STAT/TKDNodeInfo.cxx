#include "TKDNodeInfo.h"

#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMath.h"

ClassImp(TKDNodeInfo)


//_________________________________________________________________
TKDNodeInfo::TKDNodeInfo(Int_t dim):
	TObject()
	,fNDim(3*dim)
	,fData(0x0)
	,fCov(0x0)
	,fPar(0x0)
{
  // Default constructor
	fVal[0] = 0.; fVal[1] = 0.;
	Build(dim);
}

//_________________________________________________________________
TKDNodeInfo::TKDNodeInfo(const TKDNodeInfo &ref):
	TObject((TObject&) ref)
	,fNDim(fNDim)
	,fData(0x0)
	,fCov(0x0)
	,fPar(0x0)
{
  // Copy constructor
	Build(fNDim/3);

	memcpy(fData, ref.fData, fNDim*sizeof(Float_t));
	fVal[0] = ref.fVal[0];
	fVal[1] = ref.fVal[1];
	if(ref.fCov) (*fCov) = (*ref.fCov);
	if(ref.fPar) (*fPar) = (*ref.fPar);
}


//_________________________________________________________________
TKDNodeInfo::~TKDNodeInfo()
{
  // Destructor
	if(fData) delete [] fData;
	if(fCov){
		delete fPar;
		delete fCov;
	}
}

//_________________________________________________________________
TKDNodeInfo& TKDNodeInfo::operator=(const TKDNodeInfo & ref)
{
//	Info("operator==()", "...");
	
	Int_t ndim = fNDim/3;
	if(fNDim != ref.fNDim){
		fNDim = ref.fNDim;
		Build(ndim);
	}
	memcpy(fData, ref.fData, fNDim*sizeof(Float_t));
	fVal[0] = ref.fVal[0];
	fVal[1] = ref.fVal[1];
	if(ref.fCov) (*fCov) = (*ref.fCov);
	if(ref.fPar) (*fPar) = (*ref.fPar);
	
	return *this;
}

//_________________________________________________________________
void TKDNodeInfo::Build(Int_t dim)
{
// Allocate/Reallocate space for this node.

//	Info("Build()", "...");

	if(!dim) return;
	
	Int_t lambda = Int_t(1 + dim + .5*dim*(dim+1));
	if(fData) delete [] fData;
	fData = new Float_t[fNDim];
	if(fCov){
		fCov->ResizeTo(lambda, lambda);
		fPar->ResizeTo(lambda);
	}
	return;
}

//_________________________________________________________________
void TKDNodeInfo::Print(const Option_t *) const
{
  // Print the content of the node
	Int_t dim = fNDim/3;
	printf("x[");
	for(int idim=0; idim<dim; idim++) printf("%f ", fData[idim]);
	printf("] f = [%f +- %f]\n", fVal[0], fVal[1]);

	//	Float_t *bounds = &fData[dim];
	printf("range[");
	for(int idim=0; idim<dim; idim++) printf("(%f %f) ", fData[2*idim], fData[2*idim+1]);
	printf("]\n");
	
	printf("Fit parameters : ");
	if(!fCov){
		printf("Not defined.\n");
		return;
	}
	
	//	Int_t lambda = Int_t(1 + dim + .5*dim*(dim+1));
	for(int ip=0; ip<3; ip++) printf("p%d[%f] ", ip, (*fPar)(ip));
	printf("\n");
}

//_________________________________________________________________
void TKDNodeInfo::Store(const TVectorD &par, const TMatrixD &cov)
{
  // Store the parameters and the covariance in the node
	if(!fCov){
		fCov = new TMatrixD(cov.GetNrows(), cov.GetNrows());
		fPar = new TVectorD(par.GetNrows());
	}
	(*fPar) = par;
	(*fCov) = cov;
}

//_________________________________________________________________
Double_t TKDNodeInfo::CookPDF(const Double_t *point, Double_t &result, Double_t &error)
{
// Recalculate the PDF for one node from the results of interpolation (parameters and covariance matrix)

	Int_t ndim = fNDim/3;
	if(ndim>10) return 0.; // support only up to 10 dimensions

	Int_t lambda = 1 + ndim + (ndim*(ndim+1)>>1);
	Double_t fdfdp[66];
	Int_t ipar = 0;
	fdfdp[ipar++] = 1.;
	for(int idim=0; idim<ndim; idim++){
		fdfdp[ipar++] = point[idim];
		for(int jdim=idim; jdim<ndim; jdim++) fdfdp[ipar++] = point[idim]*point[jdim];
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

