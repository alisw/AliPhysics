#include "TKDInterpolatorBase.h"
#include "TKDNodeInfo.h"
#include "TKDTree.h"

#include "TClonesArray.h"
#include "TLinearFitter.h"
#include "TTree.h"
#include "TH2.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TBox.h"
#include "TGraph.h"
#include "TMarker.h"
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
  ,fNTNodes(0)
  ,fTNodes(0x0)
  ,fStatus(4)
  ,fLambda(1 + dim + (dim*(dim+1)>>1))
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
void	TKDInterpolatorBase::Build(Int_t n)
{
  // allocate memory for data

  if(fTNodes) delete fTNodes;
  fNTNodes = n;
  fTNodes = new TClonesArray("TKDNodeInfo", fNTNodes);
  for(int in=0; in<fNTNodes; in++) new ((*fTNodes)[in]) TKDNodeInfo(fNSize);
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
  if(fTNodes) delete fTNodes;
}


//__________________________________________________________________
Bool_t	TKDInterpolatorBase::GetCOGPoint(Int_t inode, Float_t *&coord, Float_t &val, Float_t &err) const
{
  if(inode < 0 || inode > fNTNodes) return kFALSE;

  TKDNodeInfo *node = (TKDNodeInfo*)(*fTNodes)[inode];
  coord = &(node->Data()[0]);
  val = node->Val()[0];
  err = node->Val()[1];
  return kTRUE;
}

//_________________________________________________________________
TKDNodeInfo* TKDInterpolatorBase::GetNodeInfo(Int_t inode) const
{
  if(!fTNodes || inode >= fNTNodes) return 0x0;
  return (TKDNodeInfo*)(*fTNodes)[inode];
}


//__________________________________________________________________
void TKDInterpolatorBase::GetStatus()
{
// Prints the status of the interpolator

  printf("Interpolator Status :\n");
  printf("  Dim    : %d [%d]\n", fNSize, fLambda);
  printf("  Method : %s\n", fStatus&1 ? "INT" : "COG");
  printf("  Store  : %s\n", fStatus&2 ? "YES" : "NO");
  printf("  Weights: %s\n", fStatus&4 ? "YES" : "NO");
  
  printf("fNTNodes %d\n", fNTNodes);        //Number of evaluation data points
  for(int i=0; i<fNTNodes; i++){
    TKDNodeInfo *node = (TKDNodeInfo*)(*fTNodes)[i]; 
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
    result = 0.;
    error = 1.E10;
    return 0.;
  }
  TKDNodeInfo *node = (TKDNodeInfo*)(*fTNodes)[nodeIndex];
  if((fStatus&1) && node->Cov() && !force) return node->CookPDF(point, result, error);

  // Allocate memory
  if(!fBuffer) fBuffer = new Double_t[2*fLambda];
  if(!fKDhelper){ 
    fRefPoints = new Float_t*[fNSize];
    for(int id=0; id<fNSize; id++){
      fRefPoints[id] = new Float_t[fNTNodes];
      for(int in=0; in<fNTNodes; in++) fRefPoints[id][in] = ((TKDNodeInfo*)(*fTNodes)[in])->Data()[id];
    }
    fKDhelper = new TKDTreeIF(fNTNodes, fNSize, 30, fRefPoints);
    fKDhelper->MakeBoundaries();
  }
  if(!fFitter) fFitter = new TLinearFitter(fLambda, Form("hyp%d", fLambda-1));
  
  // generate parabolic for nD
  //Float_t alpha = Float_t(2*lambda + 1) / fNTNodes; // the bandwidth or smoothing parameter
  //Int_t npoints = Int_t(alpha * fNTNodes);
  //printf("Params : %d NPoints %d\n", lambda, npoints);
  // prepare workers

  Int_t ipar,    // local looping variable
        npoints = Int_t((1.+fAlpha)*fLambda); // number of data points used for interpolation
  Int_t *index = new Int_t[2*npoints];  // indexes of NN 
  Float_t *dist = new Float_t[2*npoints], // distances of NN
          d,     // NN normalized distance
          w0,    // work
          w;     // tri-cubic weight function

  do{
    // find nearest neighbors
    for(int idim=0; idim<fNSize; idim++) pointF[idim] = (Float_t)point[idim];
    fKDhelper->FindNearestNeighbors(pointF, npoints+1, index, dist);
    // add points to fitter
    fFitter->ClearPoints();
    TKDNodeInfo *tnode = 0x0;
    for(int in=0; in<npoints; in++){
      tnode = (TKDNodeInfo*)(*fTNodes)[index[in]];
      //tnode->Print();
      if(fStatus&1){ // INT
        Float_t *bounds = &(tnode->Data()[fNSize]);
        ipar = 0;
        for(int idim=0; idim<fNSize; idim++){
          fBuffer[ipar++] = .5*(bounds[2*idim] + bounds[2*idim+1]);
          fBuffer[ipar++] = (bounds[2*idim]*bounds[2*idim] + bounds[2*idim] * bounds[2*idim+1] + bounds[2*idim+1] * bounds[2*idim+1])/3.;
          for(int jdim=idim+1; jdim<fNSize; jdim++) fBuffer[ipar++] = (bounds[2*idim] + bounds[2*idim+1]) * (bounds[2*jdim] + bounds[2*jdim+1]) * .25;
        }
      } else { // COG
        Float_t *p = &(tnode->Data()[0]);
        ipar = 0;
        for(int idim=0; idim<fNSize; idim++){
          fBuffer[ipar++] = p[idim];
          for(int jdim=idim; jdim<fNSize; jdim++) fBuffer[ipar++] = p[idim]*p[jdim];
        }
      }

      // calculate tri-cubic weighting function
      if(fStatus&4){
        d = dist[in]/ dist[npoints];
        w0 = (1. - d*d*d); w = w0*w0*w0;
      } else w = 1.;
      
// 			printf("x[");
// 			for(int idim=0; idim<fLambda-1; idim++) printf("%f ", fBuffer[idim]);
// 			printf("]  v[%f +- %f] (%f, %f)\n", tnode->Val()[0], tnode->Val()[1]/w, tnode->Val()[1], w);
      fFitter->AddPoint(fBuffer, tnode->Val()[0], tnode->Val()[1]/w);
    }
    npoints += 4;
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
  if(fStatus&2 && fStatus&1) node->Store(par, cov);
    
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
void TKDInterpolatorBase::DrawBins(UInt_t ax1, UInt_t ax2, Float_t ax1min, Float_t ax1max, Float_t ax2min, Float_t ax2max)
{
// Draw nodes structure projected on plane "ax1:ax2". The parameter
// "depth" specifies the bucket size per node. If depth == -1 draw only
// terminal nodes and evaluation points (default -1 i.e. bucket size per node equal bucket size specified by the user)
//
// Observation:
// This function creates the nodes (TBox) array for the specified depth
// but don't delete it. Abusing this function may cause memory leaks !


  
  TH2 *h2 = new TH2S("hNodes", "", 100, ax1min, ax1max, 100, ax2min, ax2max);
  h2->GetXaxis()->SetTitle(Form("x_{%d}", ax1));
  h2->GetYaxis()->SetTitle(Form("x_{%d}", ax2));
  h2->Draw();
  
  const Float_t kBorder = 0.;//1.E-4;
  TBox *boxArray = new TBox[fNTNodes], *box;
  Float_t *bounds = 0x0;
  for(int inode = 0; inode < fNTNodes; inode++){
    box = &boxArray[inode];
    box->SetFillStyle(3002);
    box->SetFillColor(50+inode/*Int_t(gRandom->Uniform()*50.)*/);
    
    bounds = &(((TKDNodeInfo*)(*fTNodes)[inode])->Data()[fNSize]);
    box->DrawBox(bounds[2*ax1]+kBorder, bounds[2*ax2]+kBorder, bounds[2*ax1+1]-kBorder, bounds[2*ax2+1]-kBorder);
  }

  // Draw reference points
  TGraph *ref = new TGraph(fNTNodes);
  ref->SetMarkerStyle(3);
  ref->SetMarkerSize(.7);
  ref->SetMarkerColor(2);
  for(int inode = 0; inode < fNTNodes; inode++){
    TKDNodeInfo *node = (TKDNodeInfo*)(*fTNodes)[inode];
    ref->SetPoint(inode, node->Data()[ax1], node->Data()[ax2]);
  }
  ref->Draw("p");
  return;
}

//__________________________________________________________________
void TKDInterpolatorBase::SetInterpolationMethod(Bool_t on)
{
// Set interpolation bit to "on".
  
  if(on) fStatus += fStatus&1 ? 0 : 1;
  else fStatus += fStatus&1 ? -1 : 0;
}


//_________________________________________________________________
void TKDInterpolatorBase::SetStore(Bool_t on)
{
// Set store bit to "on"
  
  if(on) fStatus += fStatus&2 ? 0 : 2;
  else fStatus += fStatus&2 ? -2 : 0;
}

//_________________________________________________________________
void TKDInterpolatorBase::SetWeights(Bool_t on)
{
// Set weights bit to "on"
  
  if(on) fStatus += fStatus&4 ? 0 : 4;
  else fStatus += fStatus&4 ? -4 : 0;
}
