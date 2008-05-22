#ifndef ALIITSRAWCLUSTERSDD_H
#define ALIITSRAWCLUSTERSDD_H

#include "AliITSRawCluster.h"

////////////////////////////////////////////////////
//  Cluster classes for set:ITS                   //
//  Raw Clusters for SDD                          //
//                                                //
////////////////////////////////////////////////////

class AliITSRawClusterSDD : public AliITSRawCluster {
 public:
    AliITSRawClusterSDD(); 
    AliITSRawClusterSDD(Int_t wing, Float_t anode, Float_t time, Float_t charge, 
			Float_t peakAmplitude,Int_t peakPosition,
			Float_t asigma, Float_t tsigma, Float_t driftPath,
			Float_t anodeOffset,  Int_t samples, 
			Int_t tstart, Int_t tstop, Int_t tstartf, Int_t tstopf,
			Int_t anodes, Int_t astart, Int_t astop);
    AliITSRawClusterSDD( const AliITSRawClusterSDD & source);
    virtual ~AliITSRawClusterSDD() {// destructor
    }
    void Add(AliITSRawClusterSDD* clJ); 
    Bool_t Brother(AliITSRawClusterSDD* cluster,Float_t dz,Float_t dx);
    void PrintInfo() const;
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
    Int_t PeakPos() const {return fPeakPosition;}
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
 protected:
    Float_t   fX;                 // X of cluster
    Float_t   fZ;                 // Z of cluster
    Float_t   fQ;                 // Q of cluster
    Int_t     fWing;              // Wing number
    Float_t   fAnode;             // Anode number
    Float_t   fTime;              // Drift Time
    Float_t   fAsigma;            // sigma (anode)
    Float_t   fTsigma;            // sigma (time)
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
#endif
