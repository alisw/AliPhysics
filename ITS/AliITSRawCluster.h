#ifndef ALIITSRAWCLUSTER_H
#define ALIITSRAWCLUSTER_H 


////////////////////////////////////////////////////
//  Cluster classes for set:ITS                   //
////////////////////////////////////////////////////

#include <TObject.h>


class AliITSRawCluster : public TObject {
  
  // this class is subject to changes ! - info used for declustering 
  // and  eventual debugging
  
public:
  
  AliITSRawCluster() {
    /*
      for (int k=0;k<100;k++) {
      fIndexMap[k]=-1;
      }
      fNcluster[0]=fNcluster[1]=-1;
      fChi2=-1; 
    */
    fMultiplicity=0;
  }
  
  virtual ~AliITSRawCluster() {
    // destructor
  }
  virtual Bool_t IsSortable() const {
    // is sortable
    return kTRUE;
  }
  
public:

  Int_t       fMultiplicity;      // cluster multiplicity
  //Int_t       fIndexMap[100];     // indices of digits
  //Int_t       fNcluster[2];
  //Float_t     fChi2;
  
  ClassDef(AliITSRawCluster,1)  // AliITSRawCluster class
    };

//---------------------------------------------
class AliITSRawClusterSPD : public AliITSRawCluster {
  
  // these classes are subject to changes - keep them temporarily for
  // compatibility !!!
  
public:
  
  AliITSRawClusterSPD() {
    // constructor
    fX=fZ=fQ;
    fZStart=fZStop;
    fNClZ=fNClX=fXStart=fXStop=fXStartf=fXStopf=fZend;
       }
  
  AliITSRawClusterSPD(Float_t clz,Float_t clx,Float_t Charge,
		      Int_t ClusterSizeZ,Int_t ClusterSizeX,Int_t xstart,Int_t xstop,
		      Int_t xstartf,Int_t xstopf,Float_t zstart,Float_t zstop,Int_t zend);
  virtual ~AliITSRawClusterSPD() {
    // destructor
  }
  
  void Add(AliITSRawClusterSPD* clJ); 
  Bool_t Brother(AliITSRawClusterSPD* cluster,Float_t dz,Float_t dx);
  void Print();
  // Getters
  Float_t Q() const {
    // Q
    return fQ ;
  }
  Float_t Z() const {
    // Z
    return fZ ;
  }
  Float_t X() const {
    // X
    return fX ;
  }
  Float_t NclZ() const {
    // NclZ
    return fNClZ ;
  }
  Float_t NclX() const {
    // NclX
    return fNClX ;
  }
  Int_t   XStart() const {
    //XStart
    return fXStart;
  }
  Int_t   XStop() const {
    //XStop
    return fXStop;
  }
  Int_t   XStartf() const {
    //XStartf
    return fXStartf;
  }
  Int_t   XStopf() const {
    //XStopf
    return fXStopf;
  }
  Float_t ZStart() const {
    //ZStart
    return fZStart;
  }
  Float_t ZStop() const {
    //ZStop
    return fZStop;
  }
  Int_t   Zend() const {
    //Zend
    return fZend;
  }
  
protected:
  
  Float_t   fX;                 // X of cluster
  Float_t   fZ;                 // Z of cluster
  Float_t   fQ;                 // Q of cluster
  Int_t     fNClZ;        // Cluster size in Z direction
  Int_t     fNClX;        // Cluster size in X direction
  Int_t     fXStart;      // number of first pixel in cluster
  Int_t     fXStop;       // number of last pixel in cluster
  Int_t     fXStartf;     // number of first pixel in full cluster
  Int_t     fXStopf;      // number of last pixel in full cluster
  Float_t   fZStart;      // number of first pixel in cluster
  Float_t   fZStop;       // number of last pixel in cluster
  Int_t     fZend;        // Zend
  
  ClassDef(AliITSRawClusterSPD,1)  // AliITSRawCluster class for SPD

    };

//---------------------------------------------
class AliITSRawClusterSDD : public AliITSRawCluster {
  
public:
  
  AliITSRawClusterSDD() {
    // constructor
    fX=fZ=fQ;
    fWing=fNsamples=0;
    fNanodes=1;
    fAnode=fTime=fPeakAmplitude=0;
  }
  
  AliITSRawClusterSDD(Int_t wing, Float_t Anode,Float_t Time,Float_t Charge,
		      Float_t PeakAmplitude,Float_t Asigma, Float_t Tsigma,Float_t DriftPath, Float_t AnodeOffset,Int_t Samples);
  virtual ~AliITSRawClusterSDD() {
    // destructor
  }
  
  void Add(AliITSRawClusterSDD* clJ); 
  Bool_t Brother(AliITSRawClusterSDD* cluster,Float_t dz,Float_t dx);
  void Print();
  // Getters
  Float_t X() const {
    //X
    return fX ;
  }
  Float_t Z() const {
    //Z
    return fZ ;
  }
  Float_t Q() const {
    //Q
    return fQ ;
  }
  Float_t A() const {
    //A
    return fAnode ;
  }
  Float_t T() const {
    //T
    return fTime ;
  }
  Float_t W() const {
    //W
    return fWing ;
  }
  Int_t Anodes() const {
    //Anodes
    return fNanodes ;
  }
  Int_t Samples() const {
    //Samples
    return fNsamples ;
  }
  Float_t PeakAmpl() const {
    //PeakAmpl
    return fPeakAmplitude ;
  }
  
protected:
  
  Float_t   fX;                 // X of cluster
  Float_t   fZ;                 // Z of cluster
  Float_t   fQ;                 // Q of cluster
  Int_t     fWing;              // Wing number
  Float_t   fAnode;             // Anode number
  Float_t   fTime;              // Drift Time
  Float_t   fPeakAmplitude;     // Peak Amplitude
  Int_t     fNanodes;           // N of anodes used for the cluster
  Int_t     fNsamples;          // N of samples used for the cluster
  
  ClassDef(AliITSRawClusterSDD,1)  // AliITSRawCluster class for SDD
    };

//-----------------------------------------
class AliITSRawClusterSSD : public AliITSRawCluster {
    
public:
  
  AliITSRawClusterSSD() {
    /*
      for (int k=0;k<100;k++) {
      fIndexMapN[k]=-1;
      }
    */
    fMultiplicityN=0;
    // fProbability=0;
    fQErr=0; 
    //fChi2N=-1;
    fStatus=-1;
  }
  AliITSRawClusterSSD(Float_t Prob,Int_t Sp,Int_t Sn);
  virtual ~AliITSRawClusterSSD() {
    // destructor
  }
  
  Int_t  GetStatus() const {
    // get status
    return fStatus;
  }
  void   SetStatus(Int_t status) {
    // set status
    fStatus=status;
  }

   
public:
  // Float_t fProbability;    // The probability that this is a "real" point
  //Int_t   fIndexMapN[100];  // indices of digits for Nside - the corresponding
                                // info for P side is carried in the base class
  //Float_t fChi2N;
  
  Int_t   fMultiplicityN;  // The number of N side strips involved 
  // in this point calculations
  Int_t   fStatus;         // Flag status : 0 - real point
  //               1 - ghost 
  //               2 - EIC ? 
  //               3 - single side 
  
  Float_t fQErr;           // Total charge error
 
  ClassDef(AliITSRawClusterSSD,1)  // AliITSRawCluster class for SSD

};


#endif
