#include <AliITSUVertexCandidate.h>

AliITSUVertexCandidate::AliITSUVertexCandidate()
:fA()
,fB()
,fV()
,fSize(0) {
  // Default Constructor

  for(int i = 0; i < 9; ++i)
    fW[i] = 0.;
  for(int i = 0; i < 6; ++i)
    fA[i] = 0.;
  for(int i = 0; i < 3; ++i)
    fB[i] = 0.;

}


void AliITSUVertexCandidate::Add(Line &line) {
  // Add a line to the cluster. It changes the weight matrix of the cluster and its parameters
  /* W (weight matrix): invert the covariance matrix  */
  double wmat[9] = {1.,0.,0.,0.,1.,0.,0.,0.,1.};
  double *p = line.x;          // point coordinates
  double *cd = line.c;         // director cosines
  double sq[3] = {1.,1.,1.};   // squared sigma of vertex.

  fSize++;

  for(int i = 0; i < 9; ++i)
    fW[i]+=wmat[i];
  // if(weight) line->GetSigma2P0(sq);
  double det = 1. / (cd[2]*cd[2]*sq[0]*sq[1]+cd[1]*cd[1]*sq[0]*sq[2]+cd[0]*cd[0]*sq[1]*sq[2]); 
  fA[0] += (cd[2]*cd[2]*sq[1]+cd[1]*cd[1]*sq[2])*det;
  fA[1] += -cd[0]*cd[1]*sq[2]*det;
  fA[2] += -cd[0]*cd[2]*sq[1]*det;
  fA[3] += (cd[2]*cd[2]*sq[0]+cd[0]*cd[0]*sq[2])*det;
  fA[4] += -cd[1]*cd[2]*sq[0]*det;
  fA[5] += (cd[1]*cd[1]*sq[0]+cd[0]*cd[0]*sq[1])*det;

  fB[0] += (cd[1]*sq[2]*(-cd[1]*p[0]+cd[0]*p[1])+cd[2]*sq[1]*(-cd[2]*p[0]+cd[0]*p[2]))*det;
  fB[1] += (cd[0]*sq[2]*(-cd[0]*p[1]+cd[1]*p[0])+cd[2]*sq[0]*(-cd[2]*p[1]+cd[1]*p[2]))*det;
  fB[2] += (cd[0]*sq[1]*(-cd[0]*p[2]+cd[2]*p[0])+cd[1]*sq[0]*(-cd[1]*p[2]+cd[2]*p[1]))*det;
}

void AliITSUVertexCandidate::ComputeClusterCentroid() {
  // Calculation of the centroid
  if (fSize < 2) {
    fV[0] = fV[1] = fV[2] = 0.;
    return;
  }
  double *a = fA;
  double *b = fB;
  double *v = fV;

  double det = a[0]*(a[3]*a[5]-a[4]*a[4])-a[1]*(a[1]*a[5]-a[4]*a[2])+a[2]*(a[1]*a[4]-a[2]*a[3]);

  if(det*det < 1.e-16) {
    DEBUG_MSG("Could not invert weight matrix");
    return;
  }
  det = 1. / det;

  v[0] = -(b[0]*(a[3]*a[5]-a[4]*a[4])-a[1]*(b[1]*a[5]-a[4]*b[2])+a[2]*(b[1]*a[4]-b[2]*a[3])) * det;
  v[1] = -(a[0]*(b[1]*a[5]-b[2]*a[4])-b[0]*(a[1]*a[5]-a[4]*a[2])+a[2]*(a[1]*b[2]-a[2]*b[1])) * det;
  v[2] = -(a[0]*(a[3]*b[2]-b[1]*a[4])-a[1]*(a[1]*b[2]-b[1]*a[2])+b[0]*(a[1]*a[4]-a[2]*a[3])) * det;
}

void AliITSUVertexCandidate::GetCovMatrix(double cov[6]) {
  // Returns the covariance matrix as the invers of the weight matrix (double precision)
  double *w = fW;
  double den = w[0]*(w[4]*w[8]-w[5]*w[7])-w[1]*(w[3]*w[8]-w[5]*w[6])+w[2]*(w[3]*w[7]-w[4]*w[6]);
  if(den == 0) {
    DEBUG_MSG("Could not invert weight matrix");
    return;
  }
  den = 1. / den;
  cov[0] =  (w[4] * w[8] - w[5] * w[7]) * den;
  cov[1] = -(w[1] * w[8] - w[7] * w[2]) * den;
  cov[2] =  (w[1] * w[5] - w[4] * w[2]) * den;
  cov[3] =  (w[0] * w[8] - w[6] * w[2]) * den;
  cov[4] = -(w[0] * w[5] - w[3] * w[2]) * den;
  cov[5] =  (w[0] * w[4] - w[1] * w[3]) * den;

}

void AliITSUVertexCandidate::Add(AliITSUVertexCandidate &cand) {
  double a[6],b[3],w[9];
  cand.GetA(a);
  cand.GetB(b);
  cand.GetWeight(w);
  fSize += cand.GetSize();
  for (int i = 0; i < 9; ++i) {
    fW[i] += w[i];
  }
  for (int i = 0; i < 6; ++i) {
    fA[i] += a[i];
  }
  for (int i = 0; i < 3; ++i) {
    fB[i] += b[i];
  }
}
