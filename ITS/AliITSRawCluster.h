#ifndef ALIITSRAWCLUSTER_H
#define ALIITSRAWCLUSTER_H

////////////////////////////////////////////////////
//  Cluster classes for set:ITS                   //
////////////////////////////////////////////////////

#include <TObject.h>
// this class is subject to changes !!! - info used for declustering
// and  eventual debugging

class AliITSRawCluster : public TObject {  
 public:
    AliITSRawCluster() { // default constructor
	fMultiplicity=0;}
    virtual ~AliITSRawCluster() {// destructor
    }
    virtual Bool_t IsSortable() const {// is sortable
	return kTRUE;}
 public:
    Int_t       fMultiplicity;        // cluster multiplicity
  
    ClassDef(AliITSRawCluster,1)  // AliITSRawCluster class
};
//----------------------------------------------------------------------
// these classes are subject to changes - keep them temporarily for
// compatibility !!!
class AliITSRawClusterSPD : public AliITSRawCluster {
 public:
    AliITSRawClusterSPD() {// constructor
	fX=fZ=fQ=0.0;
	fZStart=fZStop=0.0;
	fNClZ=fNClX=fXStart=fXStop=fZend=fNTracks=0;
	fTracks[0]=fTracks[1]=fTracks[2]=-3,fModule=0;
    }
    AliITSRawClusterSPD(Float_t clz,Float_t clx,Float_t Charge,
			Int_t ClusterSizeZ,Int_t ClusterSizeX,
			Int_t xstart,Int_t xstop,Float_t zstart,
			Float_t zstop,Int_t zend,Int_t module);
    virtual ~AliITSRawClusterSPD() {// destructor
    }
    void Add(AliITSRawClusterSPD* clJ); 
    Bool_t Brother(AliITSRawClusterSPD* cluster,Float_t dz,Float_t dx);
    void PrintInfo();
    // Getters
    Float_t Q() const {// Q
	return fQ ;}
    Float_t Z() const {// Z
	return fZ ;}
    Float_t X() const {// X
	return fX ;}
    Int_t NclZ() const {// NclZ
	return fNClZ ;}
    Int_t NclX() const {// NclX
	return fNClX ;}
    Int_t   XStart() const {//XStart
	return fXStart;}
    Int_t   XStop() const {//XStop
	return fXStop;}
    Int_t   XStartf() const {//XStartf
	return fXStart;}
    Int_t   XStopf() const {//XStopf
	return fXStop;}
    Float_t ZStart() const {//ZStart
	return fZStart;}
    Float_t ZStop() const {//ZStop
	return fZStop;}
    Int_t   Zend() const {//Zend
	return fZend;}
    Int_t   NTracks() const {//NTracks
	return fNTracks;}
    Int_t Module() const {//Returns module where this cluster came from
	return fModule;}
    void GetTracks(Int_t &track0,Int_t &track1,Int_t &track2) const {
	// returns tracks created a cluster
	track0=fTracks[0];
	track1=fTracks[1];
	track2=fTracks[2];
	return;
  };
    void   SetTracks(Int_t track0, Int_t track1, Int_t track2) {
	// set tracks in cluster (not more than three ones)
	fTracks[0]=track0;
	fTracks[1]=track1;
	fTracks[2]=track2;
    }
    void   SetNTracks(Int_t ntracks) {
	// set ntracks
	fNTracks=ntracks;
    }
 protected:
    Float_t   fX;           // X of cluster
    Float_t   fZ;           // Z of cluster
    Float_t   fQ;           // Q of cluster
    Int_t     fNClZ;        // Cluster size in Z direction
    Int_t     fNClX;        // Cluster size in X direction
    Int_t     fXStart;      // number of first pixel in cluster
    Int_t     fXStop;       // number of last pixel in cluster
    Float_t   fZStart;      // number of first pixel in cluster
    Float_t   fZStop;       // number of last pixel in cluster
    Int_t     fZend;        // Zend
    Int_t     fNTracks;     // number of tracks created a cluster
    Int_t     fTracks[3];   // tracks created a cluster
    Int_t     fModule;      // Module number for this culuster
  
  ClassDef(AliITSRawClusterSPD,2)  // AliITSRawCluster class for SPD
};
//----------------------------------------------------------------------
class AliITSRawClusterSDD : public AliITSRawCluster {
 public:
    AliITSRawClusterSDD() {
	// constructor
	fX=fZ=fQ=0;
	fWing=0;
	fNanodes=1;
	fAnode=fTime=fPeakAmplitude=0;
	fPeakPosition=-1;
	fMultiplicity=0;
	fTstart=fTstop=fTstartf=fTstopf=0;
	fAstart=fAstop=0;
  }
    AliITSRawClusterSDD(Int_t wing,Float_t Anode,Float_t Time,Float_t Charge, 
			Float_t PeakAmplitude,Int_t PeakPosition,
			Float_t Asigma,Float_t Tsigma,Float_t DriftPath,
			Float_t AnodeOffset,  Int_t Samples, 
			Int_t Tstart, Int_t Tstop,Int_t Tstartf,Int_t Tstopf,
			Int_t Anodes, Int_t Astart, Int_t Astop);
    AliITSRawClusterSDD( const AliITSRawClusterSDD & source);
    virtual ~AliITSRawClusterSDD() {// destructor
    }
    void Add(AliITSRawClusterSDD* clJ); 
    Bool_t Brother(AliITSRawClusterSDD* cluster,Float_t dz,Float_t dx);
    void PrintInfo();
    // Setters
    void SetX(Float_t x) {fX=x;}
    void SetZ(Float_t z) {fZ=z;}
    void SetQ(Float_t q) {fQ=q;}
    void SetAnode(Float_t anode) {fAnode=anode;}
    void SetTime(Float_t time) {fTime=time;}
    void SetAsigma(Float_t asigma) {fAsigma=asigma;}
    void SetTsigma(Float_t tsigma) {fTsigma=tsigma;}
    void SetWing(Int_t wing) {fWing=wing;}
    void SetNanodes(Int_t na) {fNanodes=na;}
    void SetNsamples(Int_t ns) {fMultiplicity=ns;}
    void SetPeakAmpl(Float_t ampl) {fPeakAmplitude=ampl;}
    void SetPeakPos(Int_t pos) {fPeakPosition=pos;}
    // Getters
    Float_t X() const {//X
	return fX ;}
    Float_t Z() const {//Z
	return fZ ;}
    Float_t Q() const {//Q
	return fQ ;}
    Float_t A() const {//A
	return fAnode ;}
    Float_t T() const {//T
	return fTime ;}
    Float_t Asigma() const {//Asigma
	return fAsigma ;}
    Float_t Tsigma() const {//Tsigma
	return fTsigma ;}
    Float_t W() const {//W
	return fWing ;}
    Int_t Anodes() const {//Anodes
	return fNanodes ;}
    Int_t Samples() const {//Samples
	return fMultiplicity;}
    Float_t PeakAmpl() const {//PeakAmpl
	return fPeakAmplitude ;}
    Float_t SumAmpl() const {//PeakAmpl
	return fSumAmplitude ;}
    Int_t PeakPos() {return fPeakPosition;}
    Int_t Tstart() const {//Tstart
	return fTstart ;}
    Int_t Tstartf() const {//Tstartf
	return fTstartf ;}
    Int_t Tstop() const {//Tstop
	return fTstop;}
    Int_t Tstopf() const {//Tstopf
	return fTstopf ;}
    Int_t Astart() const {//Astart
	return fAstart ;}
    Int_t Astop() const {//Astop
	return fAstop ;}
 public:
    Float_t   fX;                 // X of cluster
    Float_t   fZ;                 // Z of cluster
    Float_t   fQ;                 // Q of cluster
    Int_t     fWing;              // Wing number
    Float_t   fAnode;             // Anode number
    Float_t   fTime;              // Drift Time
    Float_t   fAsigma;            //
    Float_t   fTsigma;            //
    Float_t   fPeakAmplitude;     // Peak Amplitude
    Float_t   fSumAmplitude;      // Total Amplitude (for weighting)  
    Int_t     fPeakPosition;      // index of digit corresponding to peak
    Int_t     fNanodes;           // N of anodes used for the cluster
    Int_t     fTstart;            // First sample in 1D cluster   
    Int_t     fTstop;             // Last sample in 1D cluster 
    Int_t     fTstartf;           // First sample in the full 2D cluster
    Int_t     fTstopf;            // Last sample in the full 2D cluster  
    Int_t     fAstart;            // First anode in the 2D cluster
    Int_t     fAstop;             // last anode in the 2D cluster
  
    ClassDef(AliITSRawClusterSDD,1)  // AliITSRawCluster class for SDD
};
//----------------------------------------------------------------------
class AliITSRawClusterSSD : public AliITSRawCluster {
 public:
    AliITSRawClusterSSD() {
	fMultiplicityN=0;
	fQErr=0; 
	fSignalP=0;
	fSignalN=0;
	fStatus=-1;
	fNtracks=0;
    }
    AliITSRawClusterSSD(Float_t Prob,Int_t Sp,Int_t Sn);
    virtual ~AliITSRawClusterSSD() {// destructor
    }
    Int_t  GetStatus() const {// get status
	return fStatus;}
    void   SetStatus(Int_t status) {// set status
	fStatus=status;}
 public:
    Int_t   fMultiplicityN;  // The number of N side strips involved 
                             // in this point calculations
    Float_t fQErr;           // Total charge error
    Float_t fSignalP;        // Signal of P side cluster
    Float_t fSignalN;        // Signal of N side cluster
    Int_t   fStatus;         // Flag status : 0 - real point
                             //               1 - ghost 
                             //               2 - EIC ? 
                             //               3 - single side 
    Int_t fNtracks;          // Number of tracks created the cluster
    ClassDef(AliITSRawClusterSSD,1)  // AliITSRawCluster class for SSD
};
#endif
