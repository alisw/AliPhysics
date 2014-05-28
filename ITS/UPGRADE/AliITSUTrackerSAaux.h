#ifndef ALIITSUTRACKERSAAUX_H
#define ALIITSUTRACKERSAAUX_H

#ifdef __DEBUG__
#include <iostream>
using std::ostream;
using std::endl;
using std::cout;
#endif

#include <vector>
using std::vector;
#include "AliExternalTrackParam.h"
#include "AliITSUAux.h"

struct itsCluster {
itsCluster():isUsed(false),x(0.f),y(0.f),z(0.f),varx(0.f),covxy(0.f),vary(0.f),phi(0.f),phiM(0.f) 
#ifdef __DEBUG__ 
    ,pid()
#endif
  {}
itsCluster(const float &X,const float &Y, const float &Z, const float &varX, const float &covXY, const float &varY,const float &Phi, const float &PhiM) :
  isUsed(false),x(X),y(Y),z(Z),varx(varX),covxy(covXY),vary(varY),phi(Phi),phiM(PhiM) 
#ifdef __DEBUG__ 
    ,pid()
#endif
  {}
#ifdef __DEBUG__
itsCluster(const float &X,const float &Y, const float &Z, const float &varX, const float &covXY, const float &varY,const float &Phi, const float &PhiM,int &Pid) :
  isUsed(false),x(X),y(Y),z(Z),varx(varX),covxy(covXY),vary(varY),phi(Phi),phiM(PhiM),pid(Pid) {}
#endif
  bool isUsed;
  float x,y,z;            // Global coordinates
  float varx,covxy,vary;  // Local covariance matrix
  float phi,phiM;         // phi of the cluster and phi angle of the module containing the cluster

#ifdef __DEBUG__
  int pid;
  friend ostream& operator<<(ostream& out, const itsCluster& cl) {
    out << "pos = (" << cl.x << ", " << cl.y << ", "<< cl.z <<")"<<" phi="<<cl.phi <<endl;
    return out;
  }
#endif
};

struct Road {
  Road() : fElements(), fNElements(0) {
    ResetElements(); 
  }
  Road(const Road& copy) : fElements(), fNElements(copy.fNElements) {
    for ( int i=0; i<6; ++i ) {
      fElements[i] = copy.fElements[i];
    }
  }

  void ResetElements() {
    for ( int i=0; i<6; ++i ) {
      fElements[i] = -1;
    }
  }

  void AddElement(int i, int el) {
    fNElements++;
    fElements[i] = el;
  }
  
  int fElements[6];
  int fNElements;

  #ifdef __DEBUG__
  friend ostream& operator<<(ostream& out, const Road& cl) {
    out << "Elements ("<< cl.fNElements <<"): ";
    for ( int i=0; i<6; ++i ) out << cl.fElements[i]; 
    return out;
  }
  #endif
};

struct nPlets {
nPlets() : id0(-1),id1(-1),id2(-1),level(1),tanPhi(),tanLambda(),neighbours() 
#ifdef __DEBUG__ 
    ,pid0()
    ,pid1()
    ,pid2()
#endif
  {}
nPlets(int arg0,int arg1) : id0(arg0),id1(arg1),id2(-1),level(1),tanPhi(),tanLambda(),neighbours()
#ifdef __DEBUG__ 
    ,pid0()
    ,pid1()
    ,pid2()
#endif
  {}
  int id0,id1,id2;
  int level;
  float tanPhi,tanLambda;
  vector<int> neighbours;
#ifdef __DEBUG__
nPlets(int arg0,int arg1,int pd0,int pd1) : id0(arg0),id1(arg1),id2(-1),level(1),tanPhi(),tanLambda(),neighbours(),pid0(pd0),pid1(pd1),pid2(-1)  {}
  int pid0,pid1,pid2;
  friend ostream& operator<<(ostream& out, const nPlets& cl) {
    out << "id = (" << cl.id0 << ", " << cl.id1 << ", "<< cl.id2 <<")"<< endl;
    out << "pid = (" << cl.pid0 << ", " << cl.pid1 << ", "<< cl.pid2 <<")"<< endl;
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

 trackC() : AliExternalTrackParam(), 
 fChi2( 0. ), 
 fPoints(), 
 fNPoints(0), 
 fInnermostLayer(-1), 
 fOutermostLayer(-1) 
 {
   for ( unsigned int i = 0; i < 2*AliITSUAux::kMaxLayers; ++i ) fPoints[i]=-1;
 }

trackC(const trackC &copy) : AliExternalTrackParam(), 
fChi2(copy.fChi2), 
fPoints(),
fNPoints(copy.fNPoints), 
fInnermostLayer(copy.fInnermostLayer), 
fOutermostLayer(copy.fOutermostLayer)
{
 for ( unsigned int i = 0; i < 2*AliITSUAux::kMaxLayers; ++i ) fPoints[i]=copy.fPoints[i];
}

 trackC(int points[7]) : AliExternalTrackParam(), 
 fChi2( 0. ), 
 fPoints(),
 fNPoints( 0 ), 
 fInnermostLayer( -1 ), 
 fOutermostLayer( -1 )
 {
   bool outer=false;
   for ( unsigned int i = 0; i < 2*AliITSUAux::kMaxLayers; ++i ) fPoints[i]=-1;
   for ( int i=6; i--;  ) {
     if (points[i]!=-1) {
       if (! outer ) { 
         outer = true;
         fOutermostLayer = points[i];
       }
       fInnermostLayer = points[i];
       ++fNPoints;
     }
     fPoints[i<<0x1] = points[i];
   }
 }

  #ifdef __DEBUG__
 friend ostream& operator<<(ostream& out, const trackC& cl) {
  out << "points = (";
    for( unsigned int i = 0; i < 2*AliITSUAux::kMaxLayers; ++i ) 
      out << cl.fPoints[i] << " ";
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
  void ResetPoints() { for(unsigned int i = 0; i < 2*AliITSUAux::kMaxLayers; ++i) fPoints[i]=-1; }
  //
  Double_t fChi2;
  Int_t fPoints[2*AliITSUAux::kMaxLayers];
  Int_t fNPoints;
  Int_t fInnermostLayer;
  Int_t fOutermostLayer;

};

struct CompDesc { //Adapted from TMath ROOT code
CompDesc(vector<trackC> *d) : fData(d) {}

  bool operator()(int i1, int i2) {
    return fData->at(i1).fChi2 > fData->at(i2).fChi2;
  }

  vector<trackC> *fData;
};

#endif
