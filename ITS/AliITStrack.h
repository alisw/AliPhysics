#ifndef ALIITSTRACK_H
#define ALIITSTRACK_H

#include <TObject.h>
#include <TMatrix.h>
#include <TVector.h>

#include "../TPC/AliTPCtrack.h"

class TObjArray;
//   ITS Track Class
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it
//
class AliITStrack : public TObject { 

public:

  AliITStrack() ;
  AliITStrack(AliTPCtrack &obj);
  AliITStrack(const AliITStrack &cobj);
  AliITStrack &operator=(AliITStrack obj);
  ~AliITStrack();
  Float_t &operator()(Int_t i) const {  return fvTrack(i);}
  Float_t &operator()(Int_t r, Int_t c) const {return (*fmCovariance)(r,c);}
  Int_t GetNumClust() { return fNumClustInTrack;}
  void AddClustInTrack() { fNumClustInTrack++;}
  TObjArray *GetListOfCluster() { return flistCluster;}
  void SetChi2(Double_t chi2) { fChi2 = chi2;}
  Double_t GetChi2() { return fChi2;}
  Double_t GetZ() const {return fvTrack(1);}
  Double_t GetTgl() const {return fvTrack(3);}   
  Double_t Getrtrack() const{return rtrack;}
  Double_t Getphi()  const{return fvTrack(0);}
  Double_t GetC() const {return fvTrack(4);}
  Double_t GetD() const{return fvTrack(2);} 
  Double_t GetPt() const {return 0.3*0.2/(fvTrack(4)*100.);}            
  void SetVertex(TVector &vert) { for(Int_t i=0;i<3;i++) fVertex(i) = vert(i);}
  void SetErrorVertex(TVector &evert) {for(Int_t i=0;i<3;i++) fErrorVertex(i) = evert(i);}

  void LmTPC(); // trasform state vector and covariance matrix from local TPC to master
  TVector GetVertex() { return fVertex;}
  TVector GetErrorVertex() { return fErrorVertex;}
  Long_t  GetLabel() { return flabel;}
  void SetLabel(Long_t label) { flabel = label;}
  Int_t  GetLayer() { return fLayer;}
  TMatrix &GetCMatrix() { return *fmCovariance;}
  TVector &GetVector() { return fvTrack;}
  void SetLayer(Int_t layer) { fLayer = layer;}
  AliTPCtrack *GetTPCtrack() { return fTPCtrack;}

  void PutCluster(Int_t layerc, TVector vecclust);
  TVector GetLabTrack(Int_t lay);
  void Search(TVector VecTotLabref, Long_t &labref, Int_t &freq); 
  Float_t GetZclusterTrack(Int_t lay) {return ((Float_t) (*ClusterInTrack)(lay,2));}  
  void GetClusters();
  Int_t GetLabTPC() {return (*fTPCtrack).GetLabel();}
  Int_t GetIdPoint(Int_t lay) {return ((Int_t) (*ClusterInTrack)(lay,4));}
  Int_t GetIdModule(Int_t lay) {return ((Int_t) (*ClusterInTrack)(lay,5));}
      

  Int_t DoNotCross(Double_t rk) const;
  Double_t argA(Double_t rk) const;
  Double_t arga(Double_t rk) const;
  Double_t argB(Double_t rk) const;
  Double_t argC(Double_t rk) const;             
  void  Propagation(Double_t rk) ;

  Double_t GetSigmaphi() const{return ((Double_t)(*fmCovariance)(0,0));}
  Double_t GetSigmaZ() const{return ((Double_t)(*fmCovariance)(1,1));}
  void AddEL(Double_t signdE,  Bool_t flagtot, Double_t mass=0.1396); 
  void AddMS();
  void Correct(Double_t rk); 
  void SetDv(Double_t x) {Dv=x;}
  void SetZv(Double_t x) {Zv=x;}
  Double_t GetDv() {return Dv;}
  Double_t GetZv() {return Zv;}
  void SetsigmaDv( Double_t x) {sigmaDv=x;}
  void SetsigmaZv( Double_t x) {sigmaZv=x;} 
  Double_t GetsigmaDv() {return sigmaDv;}
  Double_t GetsigmaZv() {return sigmaZv;} 
  void PrimaryTrack();
  void Setd2(TVector &x) {for(Int_t i=0; i<6; i++){d2(i)=x(i);}}
  void Settgl2(TVector &x) {for(Int_t i=0; i<6; i++){tgl2(i)=x(i);}}
  void Setdtgl(TVector &x) {for(Int_t i=0; i<6; i++){dtgl(i)=x(i);}} 
  TVector Getd2() { return d2;}
  TVector Gettgl2() { return tgl2;}
  TVector Getdtgl() { return dtgl;}
  Double_t Getd2(Int_t i){return (Double_t)d2(i);}  
  Double_t Gettgl2(Int_t i){return (Double_t)tgl2(i);}
  Double_t Getdtgl(Int_t i){return (Double_t)dtgl(i);}
  Double_t GetxoTPC() {return xoTPC;} 
 // Double_t PhiDef(Double_t x, Double_t y); 
 // Double_t  Getalphaprov(){return alphaprov;}  //provvirio      	        
//////////////////////////////////////////////////////////////////////////////////////// 

 private:  
   
  AliTPCtrack     *fTPCtrack;           // reference to TPC track
   
  TVector         fvTrack;              // state vector: |phi/z/D/tgl/C 
  Double_t        rtrack;               // radius of courrent layer     
  TMatrix         *fmCovariance;        // Covariance Matrix
      
  Double_t        fChi2;                // fChi^2 of track         
  TObjArray       *flistCluster;        // list of clusters of the track
  Int_t           fNumClustInTrack;     // total number of clusters
  Long_t          flabel;               // label of the track
  TVector         fVertex;              // vertex coordinates of the track
  TVector         fErrorVertex;         // error on the vertex coordinates
  Int_t           fLayer;               // current Layer of the track
  TMatrix        *ClusterInTrack;       // matrix of clusters belonging to the  track
                                        // row index = layer-1; 
                                        // cols index = master coordinates of the clusters
  
  
  Double_t          Dv;                 // radial impact parameter for vertex  constraint
  Double_t          Zv;                 // longitudinal impact parameter for vertex constraint
  Double_t          sigmaDv;            // sigma for Dv extraction
  Double_t          sigmaZv;            // sigma for Zv extraction
  TVector           d2;                 // C(2,2)  per primary track
  TVector           tgl2;               // C(3,3)   per primary track
  TVector           dtgl;               // C(2,3)     per primary track

  Double_t          xoTPC;
  		   
  //Double_t alphaprov;  //provviosorio		   

 
  ClassDef(AliITStrack, 1)
   
};

#endif

