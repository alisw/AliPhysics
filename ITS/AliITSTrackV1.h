#ifndef ALIITSTRACKV1_H
#define ALIITSTRACKV1_H
//   ITS Track Class
//Origin  A. Badala' and G.S. Pappalardo:  e-mail Angela.Badala@ct.infn.it, Giuseppe.S.Pappalardo@ct.infn.it
// It contain all the usefull information for the track and the method to calculate, modify or extract them
#include <TObject.h>
#include <TMatrix.h>
#include <TVector.h>

#include "../TPC/AliTPCtrack.h"

class TObjArray;
class AliITSRad;

class AliITSTrackV1 : public TObject { 

public:

  AliITSTrackV1();  // default constructor
  //AliITSTrackV1(const char *opt);  // Standard constructor
  AliITSTrackV1(Double_t fieldfactor);  // Standard constructor
  //AliITSTrackV1(AliTPCtrack &obj);  // Standard constructor
  AliITSTrackV1(AliTPCtrack &obj, Double_t fieldfactor);  // Standard constructor
  AliITSTrackV1(const AliITSTrackV1 &cobj);  // copy constructor
  AliITSTrackV1 &operator=(AliITSTrackV1 obj);  // operator =
  ~AliITSTrackV1(); // default destructor
  Int_t GetNumClust() const { return fNumClustInTrack;}  // gets the num of cluster in a track
  void AddClustInTrack() { fNumClustInTrack++;}  // adds a cluster in track
  TObjArray *GetListOfCluster() const { return flistCluster;}  // gets the list of cluster in the track
  void SetChi2(Double_t chi2) { fChi2 = chi2;}   // sets the chisquare value for the track
  Double_t GetChi2() const { return fChi2;}            // gets the chisquare value for the track
  Double_t GetZ() const {return fX1;}            // gets the Z value for the track
  Double_t GetTgl() const {return fX3;}          // gets the tgl value for the track
  Double_t Getrtrack() const{return frtrack;}    // gets the raius value for the current track
  Double_t Getphi()  const{return fX0;}          // gets the phi value for the track
  Double_t GetC() const {return fX4;}            // gets the curvature value for the track
  Double_t GetD() const{return fX2;}             // gets the radial impact parameter for the track
  Double_t GetPt() const {return 0.299792458*0.2*fFieldFactor/(fX4*100.);} // gets the transvers momentum 
  Float_t GetdEdx() const {return fdEdx;}        //gets fdEdx    // oggi

                                                                           // value for the track           
  void SetVertex(TVector &vert) { for(Int_t i=0;i<3;i++) fVertex(i) = vert(i);} // sets the vertex
                                                                                // cohordinates
  void SetErrorVertex(TVector &evert) {for(Int_t i=0;i<3;i++) fErrorVertex(i) = evert(i);} // sets the errors
                                                                                    //for vertex cohordinates

  void LmTPC(); // trasform state vector and covariance matrix from local TPC to master
  TVector GetVertex() const { return fVertex;}
  TVector GetErrorVertex() const { return fErrorVertex;}
  Long_t  GetLabel() const { return flabel;}
  void SetLabel(Long_t label) { flabel = label;}
  Int_t  GetLayer() const { return fLayer;}
  

  void PutCElements(Double_t C00, Double_t C10, Double_t C11, Double_t C20, Double_t C21, 
  Double_t C22, Double_t C30, Double_t C31, Double_t C32, Double_t C33, Double_t C40, 
  Double_t C41, Double_t C42, Double_t C43, Double_t C44);  // put elements of covariance matrix
  
  void GetCElements(Double_t &C00, Double_t &C10, Double_t &C11, Double_t &C20, Double_t &C21, 
  Double_t &C22, Double_t &C30, Double_t &C31, Double_t &C32, Double_t &C33, Double_t &C40, 
  Double_t &C41, Double_t &C42, Double_t &C43, Double_t &C44) const;  // get elements of covariance matrix
   
  void GetXElements(Double_t &X0, Double_t &X1, Double_t &X2, Double_t &X3, Double_t &X4) const; // get elements
                                                                                           // of state vector
  void PutXElements(Double_t X0, Double_t X1, Double_t X2, Double_t X3, Double_t X4);  // put elements
 
  void PutMass(Double_t mass) {fMass=mass;} // put the particle mass
  Double_t GetMass() const {return fMass;}  // get the particle mass   // oggi                                                                                      // of state vector

    
  void SetLayer(Int_t layer) { fLayer = layer;}      // set current layer
  AliTPCtrack *GetTPCtrack() const { return fTPCtrack;}    // get hte TPC track

  void PutCluster(Int_t layerc, TVector vecclust);   // put information for clusters
  TVector GetLabTrack(Int_t lay);                    // get the label of the track
  void Search(TVector VecTotLabref, Long_t &labref, Int_t &freq); // determine the label and the frequency of
                                                                  // it for the current track
  Float_t GetZclusterTrack(Int_t lay) {return ((Float_t) (*fClusterInTrack)(lay,2));} // get the Z
                             //cohordinate of the cluster gelonging to the track for a given layer
  void GetClusters();     // prints the clusters belonging to the current track
  Int_t GetLabTPC() const {return (*fTPCtrack).GetLabel();}  // get the TPC label for the current track
  Int_t GetIdPoint(Int_t lay) {return ((Int_t) (*fClusterInTrack)(lay,4));}  // get the label identifiing the
                                                                             //point of the track
  Int_t GetIdModule(Int_t lay) {return ((Int_t) (*fClusterInTrack)(lay,5));} // get the label identifiing the
                                                                             // module of the track
  Float_t GetIdParticle(Int_t lay) {return (*fClusterInTrack)(lay,3);}       // get the label to identify
                                                                             // the particle      
  Int_t DoNotCross(Double_t rk) const;  // determine if the track cross a layer 
  Double_t ArgA(Double_t rk) const;     // quantity usefull in propagation
  Double_t Arga(Double_t rk) const;     // quantity usefull in propagation
  Double_t ArgB(Double_t rk) const;     // quantity usefull in propagation
  Double_t ArgC(Double_t rk) const;     // quantity usefull in propagation             
  void  Propagation(Double_t rk) ;      // propagation of the track to a layer of radius rk

  Double_t GetSigmaphi() const{return fC00;}    // gets the phi variance
  Double_t GetSigmaZ() const{return  fC11;}     // gets the Z variance
  void AddEL(Double_t signdE,  Bool_t flagtot, Double_t mass=0.1396);  // adds the energy loss
  void AddMS(AliITSRad *rl,Double_t mass=0.1396);  // modify the covariance matrix to take into account the multiple scattering
  void Correct(Double_t rk);  // correct the track to take into account the real detector geometry
  void SetDv(Double_t x) {fDv=x;}  // sets the radial impact parameter for vertex constraint 
  void SetZv(Double_t x) {fZv=x;}  // sets longitudinal impact parameter for vertex constraint
  Double_t GetDv() const {return fDv;}   // gets the radial impact parameter for vertex constraint 
  Double_t GetZv() const {return fZv;}   // gets longitudinal impact parameter for vertex constraint
  void SetsigmaDv( Double_t x) {fsigmaDv=x;} // sets sigma for Dv extraction
  void SetsigmaZv( Double_t x) {fsigmaZv=x;} // sets sigma for Zv extraction
  void Setfnoclust() {fnoclust++;}          //modify fnoclust 
  Double_t GetsigmaDv() const {return fsigmaDv;}   // gets sigma for Dv extraction
  Double_t GetsigmaZv() const {return fsigmaZv;}   // gets sigma for Zv extraction
  void PrimaryTrack(AliITSRad *rl);   // calculation of part of covariance matrix for vertex constraint
  void Setd2(TVector &x) {for(Int_t i=0; i<6; i++){fd2(i)=x(i);}} // sets the vector fd2
  void Settgl2(TVector &x) {for(Int_t i=0; i<6; i++){ftgl2(i)=x(i);}}  // sets the vector ftgl2
  void Setdtgl(TVector &x) {for(Int_t i=0; i<6; i++){fdtgl(i)=x(i);}}  // sets the vector fdtgl
  TVector Getd2() const { return fd2;}  // gets the vector fd2
  TVector Gettgl2() const { return ftgl2;}  // gets the vector ftgl2
  TVector Getdtgl() const { return fdtgl;}  // gets the vector dtgl
  Double_t Getd2(Int_t i){return (Double_t)fd2(i);}     // gets the i element of the vector fd2
  Double_t Gettgl2(Int_t i){return (Double_t)ftgl2(i);} // gets the i element of the vector tgl2
  Double_t Getdtgl(Int_t i){return (Double_t)fdtgl(i);} // gets the i element of the vector fdtgl
  //Double_t GetxoTPC() const {return fxoTPC;}  // gets fxoTPC
  Int_t  Getfnoclust() const {return fnoclust;}  //gets fnoclust 
  Double_t GetPredChi2(Double_t m[2], Double_t sigma[2]) const; //get predicted chi2
  void Setfcor()                  //set correction for layer                // oggi
   {if(fLayer>=3) fcor[fLayer-3] = 1./TMath::Sqrt(1.+ fX3*fX3);}            // oggi
  Float_t Getfcor(Int_t i) {return fcor[i];}  //return correction for layer // oggi

        	        
//////////////////////////////////////////////////////////////////////////////////////// 

 private:  
   
  AliTPCtrack     *fTPCtrack;           // reference to TPC track

  Double_t        fX0,fX1,fX2,fX3,fX4;  // state vector: |phi/z/D/tgl/C 
  Double_t        frtrack;               // radius of courrent layer     
  
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
  TMatrix        *fClusterInTrack;      // matrix of clusters belonging to the  track
                                        // row index = layer-1; 
                                        // cols index = master coordinates of the clusters

  Double_t        fFieldFactor;         // magnetic field factor
  Double_t          fDv;                // radial impact parameter for vertex  constraint
  Double_t          fZv;                // longitudinal impact parameter for vertex constraint
  Double_t          fsigmaDv;           // sigma for Dv extraction
  Double_t          fsigmaZv;           // sigma for Zv extraction
  TVector           fd2;                // C(2,2)  for primary track
  TVector           ftgl2;              // C(3,3)  for primary track
  TVector           fdtgl;              // C(2,3)  for primary track

  Double_t          fMass;         //  tracking particle mass
  
  Float_t           fdEdx ;         // energy loss                  // oggi
  Float_t           fcor[4];          // corrections for dE/dx      // oggi
 
  
  Int_t   fnoclust;  //nm of layers in which tracking doesn't add a cluster to the track
  		   
  

 
  ClassDef(AliITSTrackV1, 1)
   
};

#endif

