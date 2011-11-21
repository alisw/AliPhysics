#include "TKDInterpolatorBase.h"
#include "TKDNodeInfo.h"
#include "TKDTree.h"

#include "TROOT.h"
#include "TClonesArray.h"
#include "TLinearFitter.h"
#include "TTree.h"
#include "TH2.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TMath.h"
#include "TVectorD.h"
#include "TMatrixD.h"

ClassImp(TKDInterpolatorBase)

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
TKDInterpolatorBase::TKDInterpolatorBase(Int_t dim) :
  fNSize(dim)
  ,fNodes(NULL)
  ,fNodesDraw(NULL)
  ,fStatus(0)
  ,fLambda(1 + dim + (dim*(dim+1)>>1))
  ,fDepth(-1)
  ,fAlpha(.5)
  ,fRefPoints(NULL)
  ,fBuffer(NULL)
  ,fKDhelper(NULL)
  ,fFitter(NULL)
{
// Default constructor. To be used with care since in this case building
// of data structure is completly left to the user responsability.
  UseWeights();
}

//_________________________________________________________________
Bool_t TKDInterpolatorBase::Build(Int_t n)
{
  // allocate memory for data
  if(Int_t((1.+fAlpha)*fLambda) > n){ // check granularity
    Error("TKDInterpolatorBase::Build()", "Minimum number of points [%d] needed for interpolation exceeds number of evaluation points [%d]. Please increase granularity.", Int_t((1.+fAlpha)*fLambda), n);
    return kFALSE;
  }

  if(fNodes){
    Warning("TKDInterpolatorBase::Build()", "Data already allocated.");
    fNodes->Delete();
  } else {
    fNodes = new TClonesArray("TKDNodeInfo", n); 
    fNodes->SetOwner();
  }

  for(int in=0; in<n; in++) new ((*fNodes)[in]) TKDNodeInfo(fNSize);

  return kTRUE;
}

//_________________________________________________________________
Bool_t TKDInterpolatorBase::Bootstrap()
{
  if(!fNodes){
    Error("TKDInterpolatorBase::Bootstrap()", "Nodes missing. Nothing to bootstrap from.");
    return kFALSE;
  }
  Int_t in = GetNTNodes(); TKDNodeInfo *n(NULL);
  while(in--){ 
    if(!(n=(TKDNodeInfo*)(*fNodes)[in])){
      Error("TKDInterpolatorBase::Bootstrap()", "Node @ %d missing.", in);
      return kFALSE;
    }
    n->Bootstrap();
    if(!fNSize) fNSize  = n->GetDimension();
    //n->SetNode(fNSize, ...);
  }
  fLambda = n->GetNpar();
  return kTRUE;
}

//_________________________________________________________________
TKDInterpolatorBase::~TKDInterpolatorBase()
{
  if(fFitter) delete fFitter;
  if(fKDhelper) delete fKDhelper;
  if(fBuffer) delete [] fBuffer;
  
  if(fRefPoints){
    for(int idim=0; idim<fNSize; idim++) delete [] fRefPoints[idim] ;
    delete [] fRefPoints;
  }
  if(fNodes){ 
    fNodes->Delete();
    delete fNodes;
  }
  if(fNodesDraw) delete [] fNodesDraw;

  TH2 *h2=NULL;
  if((h2 = (TH2S*)gROOT->FindObject("hKDnodes"))) delete h2;
}


//__________________________________________________________________
Bool_t	TKDInterpolatorBase::GetCOGPoint(Int_t inode, Float_t *&coord, Float_t &val, Float_t &err) const
{
  if(inode < 0 || inode > GetNTNodes()) return kFALSE;

  TKDNodeInfo *node = (TKDNodeInfo*)(*fNodes)[inode];
  coord = &(node->Data()[0]);
  val = node->Val()[0];
  err = node->Val()[1];
  return kTRUE;
}

//_________________________________________________________________
TKDNodeInfo* TKDInterpolatorBase::GetNodeInfo(Int_t inode) const
{
  if(!fNodes || inode >= GetNTNodes()) return NULL;
  return (TKDNodeInfo*)(*fNodes)[inode];
}

//_________________________________________________________________
Int_t TKDInterpolatorBase::GetNTNodes() const 
{
  return fNodes?fNodes->GetEntriesFast():0;
}

//_________________________________________________________________
Bool_t TKDInterpolatorBase::GetRange(Int_t ax, Float_t &min, Float_t &max) const
{
  if(!fNodes) return kFALSE;
  Int_t ndim = ((TKDNodeInfo*)(*fNodes)[0])->GetDimension();
  if(ax<0 || ax>=ndim){
    min=0.; max=0.;
    return kFALSE;
  }
  min=1.e10; max=-1.e10;
  Float_t axmin, axmax;
  for(Int_t in=GetNTNodes(); in--; ){ 
    TKDNodeInfo *node = (TKDNodeInfo*)((*fNodes)[in]);
    node->GetBoundary(ax, axmin, axmax);
    if(axmin<min) min = axmin;
    if(axmax>max) max = axmax;
  }
  
  return kTRUE;
}

//__________________________________________________________________
void TKDInterpolatorBase::GetStatus(Option_t *opt)
{
// Prints the status of the interpolator

  printf("Interpolator Status[%d] :\n", fStatus);
  printf("  Dim    : %d [%d]\n", fNSize, fLambda);
  printf("  Method : %s\n", UseCOG() ? "COG" : "INT");
  printf("  Store  : %s\n", HasStore() ? "YES" : "NO");
  printf("  Weights: %s\n", UseWeights() ? "YES" : "NO");
  
  if(strcmp(opt, "all") != 0 ) return;
  printf("GetNTNodes() %d\n", GetNTNodes());        //Number of evaluation data points
  for(int i=0; i<GetNTNodes(); i++){
    TKDNodeInfo *node = (TKDNodeInfo*)(*fNodes)[i]; 
    printf("%d ", i); node->Print();
  }
}

//_________________________________________________________________
Double_t TKDInterpolatorBase::Eval(const Double_t *point, Double_t &result, Double_t &error, Bool_t force)
{
// Evaluate PDF for "point". The result is returned in "result" and error in "error". The function returns the chi2 of the fit.
//
// Observations:
//
// 1. The default method used for interpolation is kCOG.
// 2. The initial number of neighbors used for the estimation is set to Int(alpha*fLambda) (alpha = 1.5)
                   
  Float_t pointF[50]; // local Float_t conversion for "point"
  for(int idim=0; idim<fNSize; idim++) pointF[idim] = (Float_t)point[idim];
  Int_t nodeIndex = GetNodeIndex(pointF);
  if(nodeIndex<0){ 
    Error("TKDInterpolatorBase::Eval()", "Can not retrive node for data point.");      
    result = 0.;   
    error = 1.E10;
    return 0.;
  }
  TKDNodeInfo *node = (TKDNodeInfo*)(*fNodes)[nodeIndex];
  if(node->Par() && !force){ 
    //printf("Node @ %d\n", nodeIndex); node->Print("a");
    return node->CookPDF(point, result, error);
  }

  // Allocate memory
  if(!fBuffer) fBuffer = new Double_t[2*fLambda];
  if(!fKDhelper){ 
    fRefPoints = new Float_t*[fNSize];
    for(int id=0; id<fNSize; id++){
      fRefPoints[id] = new Float_t[GetNTNodes()];
      for(int in=0; in<GetNTNodes(); in++) fRefPoints[id][in] = ((TKDNodeInfo*)(*fNodes)[in])->Data()[id];
    }
    Info("TKDInterpolatorBase::Eval()", "Build TKDTree(%d, %d, %d)", GetNTNodes(), fNSize, kNhelper);
    fKDhelper = new TKDTreeIF(GetNTNodes(), fNSize, kNhelper, fRefPoints);
    fKDhelper->Build();
    fKDhelper->MakeBoundariesExact();
  }
  if(!fFitter) fFitter = new TLinearFitter(fLambda, Form("hyp%d", fLambda-1));
  
  // generate parabolic for nD
  //Float_t alpha = Float_t(2*lambda + 1) / GetNTNodes(); // the bandwidth or smoothing parameter
  //Int_t npoints = Int_t(alpha * GetNTNodes());
  //printf("Params : %d NPoints %d\n", lambda, npoints);
  // prepare workers

  Int_t ipar,    // local looping variable
        npoints_new = Int_t((1.+fAlpha)*fLambda),
        npoints(0); // number of data points used for interpolation
  Int_t *index = new Int_t[2*npoints_new];  // indexes of NN 
  Float_t *dist = new Float_t[2*npoints_new], // distances of NN
          d,     // NN normalized distance
          w0,    // work
          w;     // tri-cubic weight function

  Bool_t kDOWN = kFALSE;
  do{
    if(npoints){
      Info("TKDInterpolatorBase::Eval()", "Interpolation failed. Trying to increase the number of interpolation points from %d to %d.", npoints, npoints_new);
    }
    if(npoints == npoints_new){
      Error("TKDInterpolatorBase::Eval()", "Interpolation failed and number of interpolation points (%d) Can not be increased further.", npoints);
      result = 0.;
      error = 1.E10;
      // clean memory
      delete [] dist; delete [] index;
      return 0.;
    } else npoints = npoints_new;
    if(npoints > GetNTNodes()){
      Warning("TKDInterpolatorBase::Eval()", "The number of interpolation points requested (%d) exceeds number of PDF values (%d). Downscale.", npoints, GetNTNodes());
      npoints = GetNTNodes();
      kDOWN = kTRUE;
    }

    // find nearest neighbors
    for(int idim=0; idim<fNSize; idim++) pointF[idim] = (Float_t)point[idim];
    fKDhelper->FindNearestNeighbors(pointF, npoints+1, index, dist);

    // add points to fitter
    fFitter->ClearPoints();
    TKDNodeInfo *tnode = NULL;
    for(int in=0; in<npoints; in++){
      tnode = (TKDNodeInfo*)(*fNodes)[index[in]];
      //tnode->Print();
      if(UseCOG()){ // COG
        Float_t *p = &(tnode->Data()[0]);
        ipar = 0;
        for(int idim=0; idim<fNSize; idim++){
          fBuffer[ipar++] = p[idim];
          for(int jdim=idim; jdim<fNSize; jdim++) fBuffer[ipar++] = p[idim]*p[jdim];
        }
      } else { // INT
        Float_t *bounds = &(tnode->Data()[fNSize]);
        ipar = 0;
        for(int idim=0; idim<fNSize; idim++){
          fBuffer[ipar++] = .5*(bounds[2*idim] + bounds[2*idim+1]);
          fBuffer[ipar++] = (bounds[2*idim]*bounds[2*idim] + bounds[2*idim] * bounds[2*idim+1] + bounds[2*idim+1] * bounds[2*idim+1])/3.;
          for(int jdim=idim+1; jdim<fNSize; jdim++) fBuffer[ipar++] = (bounds[2*idim] + bounds[2*idim+1]) * (bounds[2*jdim] + bounds[2*jdim+1]) * .25;
        }
      }

      // calculate tri-cubic weighting function
      if(UseWeights()){
        d = dist[in]/dist[npoints];
        w0 = (1. - d*d*d); w = w0*w0*w0;
        if(w<1.e-30) continue;
      } else w = 1.;
      
//       printf("%2d d[%f] w[%f] x[", index[in], d, w);
//       for(int idim=0; idim<fLambda-1; idim++) printf("%f ", fBuffer[idim]);
//       printf("]\n");  printf("v[%f +- %f] (%f, %f)\n", tnode->Val()[0], tnode->Val()[1]/w, tnode->Val()[1], w);
      fFitter->AddPoint(fBuffer, tnode->Val()[0], tnode->Val()[1]/w);
    }
    npoints_new = npoints+ (kDOWN ? 0 : kdN);
  } while(fFitter->Eval());
  delete [] index;
  delete [] dist;

  // retrive fitter results
  TMatrixD cov(fLambda, fLambda);
  TVectorD par(fLambda);
  fFitter->GetCovarianceMatrix(cov);
  fFitter->GetParameters(par);
  Double_t chi2 = fFitter->GetChisquare()/(npoints - 4 - fLambda);

  // store results
  node->Store(&par, HasStore()?&cov:NULL);
    
  // Build df/dpi|x values
  Double_t *fdfdp = &fBuffer[fLambda];
  ipar = 0;
  fdfdp[ipar++] = 1.;
  for(int idim=0; idim<fNSize; idim++){
    fdfdp[ipar++] = point[idim];
    for(int jdim=idim; jdim<fNSize; jdim++) fdfdp[ipar++] = point[idim]*point[jdim];
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
void TKDInterpolatorBase::DrawProjection(UInt_t ax1, UInt_t ax2)
{
// Draw nodes structure projected on plane "ax1:ax2". The parameter
// "depth" specifies the bucket size per node. If depth == -1 draw only
// terminal nodes and evaluation points (default -1 i.e. bucket size per node equal bucket size specified by the user)
//
  
  Float_t ax1min, ax1max, ax2min, ax2max;
  GetRange(ax1, ax1min, ax1max);
  GetRange(ax2, ax2min, ax2max);
  TH2 *h2 = NULL;
  if(!(h2 = (TH2S*)gROOT->FindObject("hKDnodes"))){
    h2 = new TH2S("hKDnodes", "", 100, ax1min, ax1max, 100, ax2min, ax2max);
  }
  h2->GetXaxis()->SetRangeUser(ax1min, ax1max);
  h2->GetXaxis()->SetTitle(Form("x_{%d}", ax1));
  h2->GetYaxis()->SetRangeUser(ax2min, ax2max);
  h2->GetYaxis()->SetTitle(Form("x_{%d}", ax2));
  h2->Draw();


  if(!fNodesDraw) fNodesDraw = new TKDNodeInfo::TKDNodeDraw[GetNTNodes()]; 
  TKDNodeInfo::TKDNodeDraw *box = NULL;
  for(Int_t in=GetNTNodes(); in--; ){ 
    box = &(fNodesDraw[in]);
    box->SetNode((TKDNodeInfo*)((*fNodes)[in]), fNSize, ax1, ax2);
    box->Draw();
  }

  return;
}

//_________________________________________________________________
void TKDInterpolatorBase::SetAlpha(Float_t a)
{
  if(a<0.5){ 
    Warning("TKDInterpolatorBase::SetAlpha()", "The scale parameter has to be larger than 0.5");
    fAlpha = 0.5;
    return;
  }
  // check value
  if(Int_t((a+1.)*fLambda) > GetNTNodes()){
    fAlpha = TMath::Max(0.5, Float_t(GetNTNodes())/fLambda-1.);
    Warning("TKDInterpolatorBase::SetAlpha()", "Interpolation neighborhood  exceeds number of evaluation points. Downscale alpha to %f", fAlpha);
    //printf("n[%d] nodes[%d]\n", Int_t((fAlpha+1.)*fLambda), GetNTNodes());
    return;
  }
  fAlpha = a;
  return;
}

