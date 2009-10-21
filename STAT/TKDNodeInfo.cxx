////////////////////////////////////////////////////////
//
// Bucket representation for TKDInterpolator(Base)
//
// The class store data and provides the interface to draw and print.
// The bucket - generalized histogram bin in N dimensions is represented by unprocessed data like
//   - experimental PDF value and statistical error 
//   - COG position (n-tuplu)
//   - boundaries
// and interpolated data like
//   - parameters of the local parabolic fit
//   - their covariance matrix
//  
// For drawing 2D projections the helper class TKDNodeInfo::TKDNodeDraw is used.
//
// Author Alexandru Bercuci <A.Bercuci@gsi.de>
//
////////////////////////////////////////////////////////

#include "TKDNodeInfo.h"

#include "TRandom.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMath.h"

ClassImp(TKDNodeInfo)
ClassImp(TKDNodeInfo::TKDNodeDraw)


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
  fNDim = 3*dim;  
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
void TKDNodeInfo::SetNode(Int_t ndim, Float_t *data, Float_t *pdf)
{
  Build(ndim);
  memcpy(fData, data, fNDim*sizeof(Float_t));
  fVal[0]=pdf[0]; fVal[1]=pdf[1];
}


//_________________________________________________________________
void TKDNodeInfo::Print(const Option_t *) const
{
  // Print the content of the node
  Int_t dim = fNDim/3;
  printf("x[");
  for(int idim=0; idim<dim; idim++) printf("%f ", fData[idim]);
  printf("] f = [%f +- %f]\n", fVal[0], fVal[1]);

  Float_t *bounds = &fData[dim];
  printf("range[");
  for(int idim=0; idim<dim; idim++) printf("(%f %f) ", bounds[2*idim], bounds[2*idim+1]);
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



//_________________________________________________________________
TKDNodeInfo::TKDNodeDraw::TKDNodeDraw() 
  :TBox()
  ,fCOG()
  ,fNode(0x0)
{
  SetFillStyle(3002);
  SetFillColor(50+Int_t(gRandom->Uniform()*50.));

  fCOG.SetMarkerStyle(3);
  fCOG.SetMarkerSize(.7);
  fCOG.SetMarkerColor(2);
}


//_________________________________________________________________
void TKDNodeInfo::TKDNodeDraw::Draw(Option_t* option)
{
  TBox::Draw(option);
  fCOG.Draw("p");
}

//_________________________________________________________________
void TKDNodeInfo::TKDNodeDraw::SetNode(TKDNodeInfo *node, UChar_t size, UChar_t ax1, UChar_t ax2)
{
  fNode=node;
  const Float_t kBorder = 0.;//1.E-4;
  Float_t *bounds = &(node->Data()[size]);
  fX1=bounds[2*ax1]+kBorder;
  fX2=bounds[2*ax1+1]-kBorder;
  fY1=bounds[2*ax2]+kBorder;
  fY2=bounds[2*ax2+1]-kBorder;
  
  Float_t x(node->Data()[ax1]), y(node->Data()[ax2]);
  fCOG.SetX(x); fCOG.SetY(y);
}


//_________________________________________________________________
void TKDNodeInfo::TKDNodeDraw::Print(const Option_t* option) const
{
  if(!fNode) return;
  fNode->Print(option);
}
