#include "AliITSUClusterLines.h"
#include <TMath.h>
#include <iostream>
#include "AliStrLine.h"

using std::cout;
using std::endl;

ClassImp(AliITSUClusterLines);

//////////////////////////////////////////////////////////////////////
// This class is used by the AliITSUVertexer to compute and store   //
// information about vertex candidates.                             //           
// The cluster is seeded starting from two lines, then it is        //
// possible to attach other lines. Whenever a new line is attached  //
// the weight and coefficient matrices are computed.                //
// Origin puccio@to.infn.it  Feb. 20 2014                           //
//////////////////////////////////////////////////////////////////////


//____________________________________________________________________________________________________________
AliITSUClusterLines::AliITSUClusterLines() : TObject(),
					     fA(),
					     fB(),
					     fLabels(),
					     fV() {
  // Default Constructor

  for(Int_t i=0;i<9;++i) fW[i] = 0.;
}

//____________________________________________________________________________________________________________
AliITSUClusterLines::AliITSUClusterLines(UInt_t first, AliStrLine *line1, UInt_t second, AliStrLine *line2,Bool_t weight) : TObject(),
															    fA(),
															    fB(),
															    fLabels(),
															    fV() {
  // Standard constructor
  fLabels.push_back(first);
  fLabels.push_back(second);
  Double_t wmat1[9],wmat2[9],p1[3],p2[3],cd1[3],cd2[3],sq1[3]={1.,1.,1.},sq2[3]={1.,1.,1.};
  line1->GetWMatrix(wmat1);
  line2->GetWMatrix(wmat2);
  for(Int_t i=0;i<9;++i) fW[i]=wmat1[i]+wmat2[i];
  line1->GetP0(p1);
  line2->GetP0(p2);
  line1->GetCd(cd1);
  line2->GetCd(cd2);
  if(weight) {
    line1->GetSigma2P0(sq1);
    line2->GetSigma2P0(sq2);
  }
  
  Double_t det1=cd1[2]*cd1[2]*sq1[0]*sq1[1]+cd1[1]*cd1[1]*sq1[0]*sq1[2]+cd1[0]*cd1[0]*sq1[1]*sq1[2];
  Double_t det2=cd2[2]*cd2[2]*sq2[0]*sq2[1]+cd2[1]*cd2[1]*sq2[0]*sq2[2]+cd2[0]*cd2[0]*sq2[1]*sq2[2];

  fA[0]=(cd1[2]*cd1[2]*sq1[1]+cd1[1]*cd1[1]*sq1[2])/det1+(cd2[2]*cd2[2]*sq2[1]+cd2[1]*cd2[1]*sq2[2])/det2;
  fA[1]=-cd1[0]*cd1[1]*sq1[2]/det1-cd2[0]*cd2[1]*sq2[2]/det2;
  fA[2]=-cd1[0]*cd1[2]*sq1[1]/det1-cd2[0]*cd2[2]*sq2[1]/det2;
  fA[3]=(cd1[2]*cd1[2]*sq1[0]+cd1[0]*cd1[0]*sq1[2])/det1+(cd2[2]*cd2[2]*sq2[0]+cd2[0]*cd2[0]*sq2[2])/det2;
  fA[4]=-cd1[1]*cd1[2]*sq1[0]/det1-cd2[1]*cd2[2]*sq2[0]/det2;
  fA[5]=(cd1[1]*cd1[1]*sq1[0]+cd1[0]*cd1[0]*sq1[1])/det1+(cd2[1]*cd2[1]*sq2[0]+cd2[0]*cd2[0]*sq2[1])/det2;
  
  fB[0]=(cd1[1]*sq1[2]*(-cd1[1]*p1[0]+cd1[0]*p1[1])+cd1[2]*sq1[1]*(-cd1[2]*p1[0]+cd1[0]*p1[2]))/det1;
  fB[0]+=(cd2[1]*sq2[2]*(-cd2[1]*p2[0]+cd2[0]*p2[1])+cd2[2]*sq2[1]*(-cd2[2]*p2[0]+cd2[0]*p2[2]))/det2;
  fB[1]=(cd1[0]*sq1[2]*(-cd1[0]*p1[1]+cd1[1]*p1[0])+cd1[2]*sq1[0]*(-cd1[2]*p1[1]+cd1[1]*p1[2]))/det1;
  fB[1]+=(cd2[0]*sq2[2]*(-cd2[0]*p2[1]+cd2[1]*p2[0])+cd2[2]*sq2[0]*(-cd2[2]*p2[1]+cd2[1]*p2[2]))/det2;
  fB[2]=(cd1[0]*sq1[1]*(-cd1[0]*p1[2]+cd1[2]*p1[0])+cd1[1]*sq1[0]*(-cd1[1]*p1[2]+cd1[2]*p1[1]))/det1;
  fB[2]+=(cd2[0]*sq2[1]*(-cd2[0]*p2[2]+cd2[2]*p2[0])+cd2[1]*sq2[0]*(-cd2[1]*p2[2]+cd2[2]*p2[1]))/det2;

  this->ComputeClusterCentroid();

}

//____________________________________________________________________________________________________________
AliITSUClusterLines::~AliITSUClusterLines() {
  // Destructor
}

//____________________________________________________________________________________________________________
void AliITSUClusterLines::Add(UInt_t label, AliStrLine *line, Bool_t weight) {
  // Add a line to the cluster. It changes the weight matrix of the cluster and its parameters
  fLabels.push_back(label);
  Double_t wmat[9],p[3],cd[3],sq[3]={1.,1.,1.};
  line->GetWMatrix(wmat);
  line->GetP0(p);
  line->GetCd(cd);

  for(Int_t i=0;i<9;++i) fW[i]+=wmat[i];
  if(weight) line->GetSigma2P0(sq);
    
  Double_t det=cd[2]*cd[2]*sq[0]*sq[1]+cd[1]*cd[1]*sq[0]*sq[2]+cd[0]*cd[0]*sq[1]*sq[2];
  fA[0]+=(cd[2]*cd[2]*sq[1]+cd[1]*cd[1]*sq[2])/det;
  fA[1]+=-cd[0]*cd[1]*sq[2]/det;
  fA[2]+=-cd[0]*cd[2]*sq[1]/det;
  fA[3]+=(cd[2]*cd[2]*sq[0]+cd[0]*cd[0]*sq[2])/det;
  fA[4]+=-cd[1]*cd[2]*sq[0]/det;
  fA[5]+=(cd[1]*cd[1]*sq[0]+cd[0]*cd[0]*sq[1])/det;

  fB[0]+=(cd[1]*sq[2]*(-cd[1]*p[0]+cd[0]*p[1])+cd[2]*sq[1]*(-cd[2]*p[0]+cd[0]*p[2]))/det;
  fB[1]+=(cd[0]*sq[2]*(-cd[0]*p[1]+cd[1]*p[0])+cd[2]*sq[0]*(-cd[2]*p[1]+cd[1]*p[2]))/det; 
  fB[2]+=(cd[0]*sq[1]*(-cd[0]*p[2]+cd[2]*p[0])+cd[1]*sq[0]*(-cd[1]*p[2]+cd[2]*p[1]))/det;

  this->ComputeClusterCentroid();
}

//____________________________________________________________________________________________________________
Int_t AliITSUClusterLines::Compare(const TObject* obj) const {
  // Comparison criteria between two clusters
  const AliITSUClusterLines *cl=(const AliITSUClusterLines*)obj;
  if(fLabels.size()<cl->GetSize()) return 1;
  if(fLabels.size()>cl->GetSize()) return -1;
  return 0;
}

//____________________________________________________________________________________________________________
void AliITSUClusterLines::ComputeClusterCentroid() {
  // Calculation of the centroid
  Double_t *a=fA;
  Double_t *b=fB;
  Double_t *v=fV;

  Double_t det=a[0]*(a[3]*a[5]-a[4]*a[4])-a[1]*(a[1]*a[5]-a[4]*a[2])+a[2]*(a[1]*a[4]-a[2]*a[3]);

  if(det==0) {
    cout << "Could not invert weight matrix" << endl;
    return;
  }
  v[0]=-(b[0]*(a[3]*a[5]-a[4]*a[4])-a[1]*(b[1]*a[5]-a[4]*b[2])+a[2]*(b[1]*a[4]-b[2]*a[3]))/det;
  v[1]=-(a[0]*(b[1]*a[5]-b[2]*a[4])-b[0]*(a[1]*a[5]-a[4]*a[2])+a[2]*(a[1]*b[2]-a[2]*b[1]))/det;
  v[2]=-(a[0]*(a[3]*b[2]-b[1]*a[4])-a[1]*(a[1]*b[2]-b[1]*a[2])+b[0]*(a[1]*a[4]-a[2]*a[3]))/det;
}

//____________________________________________________________________________________________________________
void AliITSUClusterLines::GetCovMatrix(Float_t cov[6]) {
  // Returns the covariance matrix (single precision)
  Double_t *w=fW;
  Double_t den=w[0]*(w[4]*w[8]-w[5]*w[7])-w[1]*(w[3]*w[8]-w[5]*w[6])+w[2]*(w[3]*w[7]-w[4]*w[6]);
  if(den==0) {
    cout << "Could not invert weight matrix" << endl;
    return;
  }
  cov[0]=(w[4]*w[8]-w[5]*w[7])/den;
  cov[1]=-(w[1]*w[8]-w[7]*w[2])/den;
  cov[2]=(w[1]*w[5]-w[4]*w[2])/den;
  cov[3]=(w[0]*w[8]-w[6]*w[2])/den;
  cov[4]=-(w[0]*w[5]-w[3]*w[2])/den;
  cov[5]=(w[0]*w[4]-w[1]*w[3])/den;
}

//____________________________________________________________________________________________________________
void AliITSUClusterLines::GetCovMatrix(Double_t cov[6]) {
  // Returns the covariance matrix (double precision)
  Double_t *w=fW;
  Double_t den=w[0]*(w[4]*w[8]-w[5]*w[7])-w[1]*(w[3]*w[8]-w[5]*w[6])+w[2]*(w[3]*w[7]-w[4]*w[6]);
  if(den==0) {
    cout << "Could not invert weight matrix" << endl;
    return;
  }
  cov[0]=(w[4]*w[8]-w[5]*w[7])/den;
  cov[1]=-(w[1]*w[8]-w[7]*w[2])/den;
  cov[2]=(w[1]*w[5]-w[4]*w[2])/den;
  cov[3]=(w[0]*w[8]-w[6]*w[2])/den;
  cov[4]=-(w[0]*w[5]-w[3]*w[2])/den;
  cov[5]=(w[0]*w[4]-w[1]*w[3])/den;

}

//____________________________________________________________________________________________________________
Bool_t AliITSUClusterLines::IsEqual(const TObject* obj) const {
  // Comparison criteria between two clusters
  const AliITSUClusterLines *cl=(const AliITSUClusterLines*)obj;
  if(fLabels.size()==cl->GetSize()) return kTRUE;
  return kFALSE;
}

