#ifndef ALIITSUTRACKERSAAUX_H
#define ALIITSUTRACKERSAAUX_H

#include <vector>
using std::vector;

/* struct CompareAsc { // Adapted from TMath ROOT code */
/* CompareAsc(Float_t *d,Float_t offset) : fData(d) {} */
  
/*   bool operator()(int i1, int i2) { */
/*     return *(fData + i1) < *(fData + i2); */
/*   } */
  
/*   Float_t fData; */
/* }; */

struct itsCluster {
itsCluster():isUsed(false),x(0.f),y(0.f),z(0.f),varx(0.f),covxy(0.f),vary(0.f),phi(0.f),phiM(0.f) {}
itsCluster(const float &X,const float &Y, const float &Z, const float &varX, const float &covXY, const float &varY,const float &Phi, const float &PhiM) : 
  isUsed(false),x(X),y(Y),z(Z),varx(varX),covxy(covXY),vary(varY),phi(Phi),phiM(PhiM) {}
  bool isUsed;
  float x,y,z;            // Global coordinates
  float varx,covxy,vary;  // Local covariance matrix
  float phi,phiM;         // phi of the cluster and phi angle of the module containing the cluster
};

struct nPlets {
nPlets() : id0(-1),id1(-1),id2(-1),level(0),tanPhi(),tanLambda(),neighbours() {}
nPlets(int arg0,int arg1) : id0(arg0),id1(arg1),id2(-1),level(0),tanPhi(),tanLambda(),neighbours()  {}
  int id0,id1,id2;
  int level;
  float tanPhi,tanLambda;
  vector<int> neighbours;
};

#endif
