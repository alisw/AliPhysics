#ifndef ALIITSUTRACKERSAAUX_H
#define ALIITSUTRACKERSAAUX_H

#ifdef __DEBUG__ 
#include <iostream>
using std::ostream;
using std::endl;
#endif

#include <vector>
using std::vector;
#include "AliExternalTrackParam.h"
#include "AliITSUAux.h"

struct itsCluster {
itsCluster():isUsed(false),x(0.f),y(0.f),z(0.f),varx(0.f),covxy(0.f),vary(0.f),phi(0.f),phiM(0.f) {}
itsCluster(const float &X,const float &Y, const float &Z, const float &varX, const float &covXY, const float &varY,const float &Phi, const float &PhiM) : 
  isUsed(false),x(X),y(Y),z(Z),varx(varX),covxy(covXY),vary(varY),phi(Phi),phiM(PhiM) {}
  bool isUsed;
  float x,y,z;            // Global coordinates
  float varx,covxy,vary;  // Local covariance matrix
  float phi,phiM;         // phi of the cluster and phi angle of the module containing the cluster
  
#ifdef __DEBUG__
  friend ostream& operator<<(ostream& out, const itsCluster& cl) {
    out << "pos = (" << cl.x << ", " << cl.y << ", "<< cl.z <<")"<<" phi="<<cl.phi <<endl;
    return out;
  }
#endif
};

struct nPlets {
nPlets() : id0(-1),id1(-1),id2(-1),level(1),tanPhi(),tanLambda(),neighbours() {}
nPlets(int arg0,int arg1) : id0(arg0),id1(arg1),id2(-1),level(1),tanPhi(),tanLambda(),neighbours()  {}
  int id0,id1,id2;
  int level;
  float tanPhi,tanLambda;
  vector<int> neighbours;
#ifdef __DEBUG__
  friend ostream& operator<<(ostream& out, const nPlets& cl) {
    out << "id = (" << cl.id0 << ", " << cl.id1 << ", "<< cl.id2 <<")"<< endl;
    out << "tanPhi="<< cl.tanPhi <<" tanLambda="<<cl.tanLambda << " level=" << cl.level <<endl;
    out << "neighbours= ";
    for( unsigned int i = 0; i< cl.neighbours.size(); ++i ) out << cl.neighbours[i];
    out << endl;
    return out;
  }
#endif
};

class trackC : public AliExternalTrackParam {
 public : 
 trackC() : AliExternalTrackParam(), fChi2( 0. ), fPoints() {
    for ( unsigned int i = 0; i < 2*AliITSUAux::kMaxLayers; ++i ) fPoints[i]=-1;
  }
#ifdef __DEBUG__
  friend ostream& operator<<(ostream& out, const trackC& cl) {
    out << "points = (";
    for( unsigned int i = 0; i < 2*AliITSUAux::kMaxLayers; ++i ) out << cl.fPoints[i] << " ";
    out << "), chi2: " << cl.fChi2 <<endl;
    const double* par = cl.GetParameter();
    const double* cov = cl.GetCovariance();    
    out << "X: " << cl.GetX() << " Alpha: " << cl.GetAlpha() << endl;
    out << "Param: \n";
    for (int i=0;i<5;i++) out << par[i] << " "; out << endl;
    out << "Covar: \n";
    int cnt = 0;
    for (int i=0;i<5;i++) {for (int j=i+1;j--;) out << cov[cnt++] << " "; out << endl;}
    return out;
  }
#endif
  
  /* bool operator < ( const trackC &t ) { return fChi2 < t.fChi2; } */
  /* bool operator <= ( const trackC &t ) { return fChi2 <= t.fChi2; } */
  /* bool operator > ( const trackC &t ) { return fChi2 > t.fChi2; } */
  /* bool operator >= ( const trackC &t ) { return fChi2 >= t.fChi2; } */
  /* bool operator == ( const trackC &t ) { return fChi2 == t.fChi2; } */
  /* bool operator != ( const trackC &t ) { return fChi2 != t.fChi2; } */
  Double_t fChi2;
  Int_t fPoints[2*AliITSUAux::kMaxLayers]; 
  
};


struct CompDesc { //Adapted from TMath ROOT code 
CompDesc(vector<trackC> *d) : fData(d) {} 
  
  bool operator()(int i1, int i2) { 
    return fData->at(i1).fChi2 > fData->at(i2).fChi2; 
  } 
  
  vector<trackC> *fData; 
}; 

#endif
