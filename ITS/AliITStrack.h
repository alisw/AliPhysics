#ifndef ALIITSTRACK_H
#define ALIITSTRACK_H

#include <TObject.h>
#include <TMatrix.h>
#include <TVector.h>

#include "../TPC/AliTPCtrack.h"

class TObjArray;
class AliITSRad;

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
  Int_t GetNumClust() { return fNumClustInTrack;}
  void AddClustInTrack() { fNumClustInTrack++;}
  TObjArray *GetListOfCluster() { return flistCluster;}
  void SetChi2(Double_t chi2) { fChi2 = chi2;}
  Double_t GetChi2() { return fChi2;}
  Double_t GetZ() const {return fX1;}
  Double_t GetTgl() const {return fX3;}   
  Double_t Getrtrack() const{return rtrack;}
  Double_t Getphi()  const{return fX0;}
  Double_t GetC() const {return fX4;}
  Double_t GetD() const{return fX2;} 
  Double_t GetPt() const {return 0.299792458*0.2/(fX4*100.);}            
  void SetVertex(TVector &vert) { for(Int_t i=0;i<3;i++) fVertex(i) = vert(i);}
  void SetErrorVertex(TVector &evert) {for(Int_t i=0;i<3;i++) fErrorVertex(i) = evert(i);}

  void LmTPC(); // trasform state vector and covariance matrix from local TPC to master
  TVector GetVertex() { return fVertex;}
  TVector GetErrorVertex() { return fErrorVertex;}
  Long_t  GetLabel() { return flabel;}
  void SetLabel(Long_t label) { flabel = label;}
  Int_t  GetLayer() { return fLayer;}
  

  void PutCElements(Double_t C00, Double_t C10, Double_t C11, Double_t C20, Double_t C21, 
  Double_t C22, Double_t C30, Double_t C31, Double_t C32, Double_t C33, Double_t C40, 
  Double_t C41, Double_t C42, Double_t C43, Double_t C44);
  
  void GetCElements(Double_t &C00, Double_t &C10, Double_t &C11, Double_t &C20, Double_t &C21, 
  Double_t &C22, Double_t &C30, Double_t &C31, Double_t &C32, Double_t &C33, Double_t &C40, 
  Double_t &C41, Double_t &C42, Double_t &C43, Double_t &C44);
   
  void GetXElements(Double_t &X0, Double_t &X1, Double_t &X2, Double_t &X3, Double_t &X4);
  void PutXElements(Double_t X0, Double_t X1, Double_t X2, Double_t X3, Double_t X4);
    
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
  Float_t GetIdParticle(Int_t lay) {return (*ClusterInTrack)(lay,3);}
      

  Int_t DoNotCross(Double_t rk) const;
  Double_t argA(Double_t rk) const;
  Double_t arga(Double_t rk) const;
  Double_t argB(Double_t rk) const;
  Double_t argC(Double_t rk) const;             
  void  Propagation(Double_t rk) ;

  Double_t GetSigmaphi() const{return fC00;}
  Double_t GetSigmaZ() const{return  fC11;}
  void AddEL(AliITSRad *rl,Double_t signdE,  Bool_t flagtot, Double_t mass=0.1396); 
  void AddMS(AliITSRad *rl);
  void Correct(Double_t rk); 
  void SetDv(Double_t x) {Dv=x;}
  void SetZv(Double_t x) {Zv=x;}
  Double_t GetDv() {return Dv;}
  Double_t GetZv() {return Zv;}
  void SetsigmaDv( Double_t x) {sigmaDv=x;}
  void SetsigmaZv( Double_t x) {sigmaZv=x;} 
  Double_t GetsigmaDv() {return sigmaDv;}
  Double_t GetsigmaZv() {return sigmaZv;} 
  void PrimaryTrack(AliITSRad *rl);
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
  Double_t PhiDef(Double_t x, Double_t y); 
  //Int_t  Getfreq(){return freq;}  //provvisorio 
  // void  Setfreq(Int_t xfreq){freq=xfreq;}  //provvisorio        	        
//////////////////////////////////////////////////////////////////////////////////////// 

 private:  
   
  AliTPCtrack     *fTPCtrack;           // reference to TPC track

  Double_t        fX0,fX1,fX2,fX3,fX4;  // state vector: |phi/z/D/tgl/C 
  Double_t        rtrack;               // radius of courrent layer     
  
  Double_t        fC00, fC10, fC11,     // Covariance Matrix
                  fC20, fC21, fC22,     //      "       "
						fC30, fC31, fC32,     //      "       " 
						fC33, fC40, fC41,     //      "       " 
						fC42, fC43, fC44;     //      "       "
      
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
  		   
 // Int_t freq; //provvisorio	   

 
  ClassDef(AliITStrack, 1)
   
};

#endif

