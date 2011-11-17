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

#include "TVectorD.h"
#include "TMatrixD.h"
#include "TRandom.h"
#include "TMath.h"

ClassImp(TKDNodeInfo)
ClassImp(TKDNodeInfo::TKDNodeDraw)


//_________________________________________________________________
TKDNodeInfo::TKDNodeInfo(Int_t dim):
  TObject()
  ,fNDim(3*dim)
  ,fData(NULL)
  ,fNpar(0)
  ,fNcov(0)
  ,fPar(NULL)
  ,fCov(NULL)
{
  // Default constructor
  fVal[0] = 0.; fVal[1] = 0.;
  Build(dim);
}

//_________________________________________________________________
TKDNodeInfo::TKDNodeInfo(const TKDNodeInfo &ref):
  TObject((TObject&) ref)
  ,fNDim(ref.fNDim)
  ,fData(NULL)
  ,fNpar(0)
  ,fNcov(0)
  ,fPar(NULL)
  ,fCov(NULL)
{
  // Copy constructor
  Build(fNDim/3);

  fData = new Float_t[fNDim];
  memcpy(fData, ref.fData, fNDim*sizeof(Float_t));
  fVal[0] = ref.fVal[0];
  fVal[1] = ref.fVal[1];
  if(ref.fNpar&&ref.fPar){ 
    fNpar = ref.fNpar;
    fPar=new Double_t[fNpar];
    memcpy(fPar, ref.fPar, fNpar*sizeof(Double_t));
  }
  if(ref.fNcov && ref.fCov){ 
    fNcov = ref.fNcov;
    fCov=new Double_t[fNcov];
    memcpy(fCov, ref.fCov, fNcov*sizeof(Double_t));
  }
}


//_________________________________________________________________
TKDNodeInfo::~TKDNodeInfo()
{
  // Destructor
  if(fData) delete [] fData;
  if(fPar) delete [] fPar;
  if(fCov) delete [] fCov;
}

//_________________________________________________________________
TKDNodeInfo& TKDNodeInfo::operator=(const TKDNodeInfo & ref)
{
//	Info("operator==()", "...");
  
  if(this == &ref) return *this;
  Int_t ndim = fNDim/3;
  if(fNDim != ref.fNDim){
    fNDim = ref.fNDim;
    Build(ndim);
  }
  memcpy(fData, ref.fData, fNDim*sizeof(Float_t));
  fVal[0] = ref.fVal[0];
  fVal[1] = ref.fVal[1];
  if(ref.fNpar&&ref.fPar){ 
    fNpar = ref.fNpar;
    fPar=new Double_t[fNpar];
    memcpy(fPar, ref.fPar, fNpar*sizeof(Double_t));
  }
  if(ref.fNcov && ref.fCov){ 
    fNcov = ref.fNcov;
    fCov=new Double_t[fNcov];
    memcpy(fCov, ref.fCov, fNcov*sizeof(Double_t));
  }
  return *this;
}

//_________________________________________________________________
void TKDNodeInfo::Build(Int_t dim)
{
// Allocate/Reallocate space for this node.

//	Info("Build()", "...");

  if(!dim) return;
  fNDim = 3*dim;
  if(fData) delete [] fData;
  fData = new Float_t[fNDim];
  return;
}

//_________________________________________________________________
void TKDNodeInfo::Bootstrap()
{
  if(!fNpar || !fPar) return;

  Int_t ndim = Int_t(.5*(TMath::Sqrt(1.+8.*fNpar)-1.))-1;
  fNDim = 3*ndim;
}

//_________________________________________________________________
void TKDNodeInfo::SetNode(Int_t ndim, Float_t *data, Float_t *pdf)
{
  Build(ndim);
  memcpy(fData, data, fNDim*sizeof(Float_t));
  fVal[0]=pdf[0]; fVal[1]=pdf[1];
}


//_________________________________________________________________
void TKDNodeInfo::Print(const Option_t *opt) const
{
  // Print the content of the node
  Int_t dim = Int_t(fNDim/3.);
  printf("x[");
  for(int idim=0; idim<dim; idim++) printf("%f ", fData?fData[idim]:0.);
  printf("] f = [%f +- %f]\n", fVal[0], fVal[1]);

  if(fData){
    Float_t *bounds = &fData[dim];
    printf("range[");
    for(int idim=0; idim<dim; idim++) printf("(%f %f) ", bounds[2*idim], bounds[2*idim+1]);
    printf("]\n");
  }
  if(strcmp(opt, "a")!=0) return;

  if(fNpar){ 
    printf("Fit parameters : \n");
    for(int ip=0; ip<fNpar; ip++) printf("p%d[%f] ", ip, fPar[ip]);
    printf("\n");
  }
  if(!fNcov) return;
  for(int ip(0), n(0); ip<fNpar; ip++){
    for(int jp(ip); jp<fNpar; jp++) printf("c(%d %d)[%f] ", ip, jp, fCov[n++]);
    printf("\n");
  }
}

//_________________________________________________________________
void TKDNodeInfo::Store(TVectorD const *par, TMatrixD const *cov)
{
// Store the parameters and the covariance in the node

  if(!fPar){SetNpar(); fPar = new Double_t[fNpar];}
  for(int ip=0; ip<fNpar; ip++) fPar[ip] = (*par)[ip];

  if(!cov) return;
  if(!fCov){SetNcov(); fCov = new Double_t[fNcov];}
  for(int ip(0), np(0); ip<fNpar; ip++)
    for(int jp=ip; jp<fNpar; jp++) fCov[np++] = (*cov)(ip,jp);
}

//_________________________________________________________________
Bool_t TKDNodeInfo::CookPDF(const Double_t *point, Double_t &result, Double_t &error) const
{
// Recalculate the PDF for one node from the results of interpolation (parameters and covariance matrix)

  Int_t ndim = Int_t(fNDim/3.);
  if(ndim>10) return kFALSE; // support only up to 10 dimensions
  //printf("ndim[%d] npar[%d] ncov[%d]\n", ndim, fNpar, fNcov);

  Double_t fdfdp[66]; memset(fdfdp, 0, ndim*sizeof(Double_t));
  Int_t ipar = 0;
  fdfdp[ipar++] = 1.;
  for(int idim=0; idim<ndim; idim++){
    fdfdp[ipar++] = point[idim];
    for(int jdim=idim; jdim<ndim; jdim++) fdfdp[ipar++] = point[idim]*point[jdim];
  }

  // calculate estimation
  result =0.; error = 0.;
  for(int i=0; i<fNpar; i++) result += fdfdp[i]*fPar[i];
  if(!fNcov) return kTRUE;

  for(int i(0), n(0); i<fNpar; i++){
    error += fdfdp[i]*fdfdp[i]*fCov[n++];
    for(int j(i+1); j<fNpar; j++) error += 2.*fdfdp[i]*fdfdp[j]*fCov[n++];
  }	
  error = TMath::Sqrt(error);
  
  //printf("TKDNodeInfo::CookPDF() : %6.3f +- %6.3f\n", result, error);

  return kTRUE;
}



//_________________________________________________________________
TKDNodeInfo::TKDNodeDraw::TKDNodeDraw() 
  :TBox()
  ,fCOG()
  ,fNode(NULL)
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
