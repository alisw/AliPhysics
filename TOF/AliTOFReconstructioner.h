#ifndef ALITOFRECONSTRUCTIONER_H
#define ALITOFRECONSTRUCTIONER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


//_________________________________________________________________________
//  Task Class for Reconstruction in TOF      
//                  
//-- Authors: Bologna-ITEP-Salerno Group


#include "TTask.h"
#include "TString.h"
#include <TF1.h>
#include <TF2.h>
#include <TH1.h>
#include <TH2.h>
#include <TRandom.h>
#include <TNtuple.h>
#include "TParticle.h"
#include "AliTOF.h"
#include "AliDetector.h"

class AliTOFPad;
class AliTOFRecHit;
class AliTOFTrack;
class TVirtualMC;

class AliTOFReconstructioner: public TTask {

public:
  AliTOFReconstructioner() ;          // ctor
  AliTOFReconstructioner(char* headerFile, Option_t* opt, char *RecFile = 0) ; 
  virtual ~AliTOFReconstructioner() ; // dtor   
  virtual void  InitArray(Float_t array[], Int_t nlocations);
  virtual void  InitArray(Int_t   array[], Int_t nlocations);
  virtual void  Exec(const char* datafile, Option_t* option); // do the main work
  virtual void  ReadTOFHits(Int_t ntracks, TTree* treehits, TClonesArray* tofhits, Int_t *** MapPixels, Int_t* kTOFhitFirst, AliTOFPad* pixelArray ,  Int_t* iTOFpixel, Float_t* toftime, AliTOFRecHit* hitArray, Int_t& isHitOnFiredPad, Int_t& ipixel);
  virtual void  AddNoiseFromOuter(Option_t *option, Int_t *** MapPixels, AliTOFPad* pixelArray , AliTOFRecHit* hitArray, Int_t& isHitOnFiredPad, Int_t& ipixel);
  virtual void  SetMinDistance(AliTOFRecHit* hitArray, Int_t ilastEntry);
  // this line has to be commented till TPC will provide fPx fPy fPz 
  // and fL in AliTPChit class
  //virtual void  ReadTPCHits(Int_t ntracks, TTree* treehits, TClonesArray* tpchits, Int_t* iTrackPt, Int_t* iparticle, Float_t* ptTrack, AliTOFTrack* trackArray, Int_t& itrack);
  virtual void  ReadTPCTracks(TFile* /*tpcReconFile*/){};
  virtual void  Matching(AliTOFTrack* trackArray, AliTOFRecHit* hitArray, Int_t *** mapPixel, AliTOFPad* pixelArray, Int_t* kTOFhitFirst, Int_t& ipixel, Int_t* iTrackPt, Int_t* iTOFpixel, Int_t ntotTpcTracks);

  virtual void FillNtuple(Int_t ntracks, AliTOFTrack* trackArray, AliTOFRecHit* hitArray, AliTOFPad* pixelArray, Int_t* iTOFpixel, Int_t* iparticle, Float_t* toftime, Int_t& ipixelLastEntry, Int_t itrack); 
  void          Init(Option_t* opt);
  void          CreateNTuple();
  void          SetNEvents(Int_t Nevents) {fNevents = Nevents;}
  void          SetFirstEvent(Int_t firstevent) {fFirstEvent = firstevent;}
  void          SetLastEvent(Int_t lastevent) {fLastEvent = lastevent;}
  Int_t         GetNEvents() const {return fNevents;}
  const char*   GetRecFile() const {return fRecFile.Data();}
  Int_t         PDGtoGeantCode(Int_t pdgcode); 
  virtual void  IsInsideThePad(TVirtualMC* vmc, Float_t x, Float_t y, Float_t z, Int_t *nGeom, Float_t& zPad, Float_t& xPad);
  virtual void  BorderEffect(Float_t z0, Float_t x0, Float_t geantTime, Int_t& nActivatedPads, Int_t& nFiredPads, Bool_t* isFired, Int_t* nPlace, Float_t* qInduced, Float_t* tofTime, Float_t& averageTime);
  virtual void  EpMulScatt(Float_t& px, Float_t& py, Float_t& pz, Float_t& p, Float_t& theta);
  void  SetDbg(Int_t dbgflag)                        {fdbg=dbgflag;}
  void  SetRecFile(char * file ) ;
  void  SetTimeResolution(Float_t timeResolution)    {fTimeResolution=timeResolution;}
  void  SetPadefficiency(Float_t padefficiency)      {fpadefficiency=padefficiency;}
  void  SetEdgeEffect(Int_t   edgeEffect)            {fEdgeEffect=edgeEffect;}
  void  SetEdgeTails(Int_t   edgeTails)              {fEdgeTails=edgeTails;}
  void  SetHparameter(Float_t hparameter)            {fHparameter=hparameter;}
  void  SetH2parameter(Float_t h2parameter)          {fH2parameter=h2parameter;}
  void  SetKparameter(Float_t kparameter)            {fKparameter=kparameter;}
  void  SetK2parameter(Float_t k2parameter)          {fK2parameter=k2parameter;}
  void  SetEffCenter(Float_t effCenter)              {fEffCenter=effCenter;}
  void  SetEffBoundary(Float_t effBoundary)          {fEffBoundary=effBoundary;}
  void  SetEff2Boundary(Float_t eff2Boundary)        {fEff2Boundary=eff2Boundary;}
  void  SetEff3Boundary(Float_t eff3Boundary)        {fEff3Boundary=eff3Boundary;}
  void  SetResCenter (Float_t resCenter)             {fResCenter=resCenter;}
  void  SetResBoundary(Float_t resBoundary)          {fResBoundary=resBoundary;}
  void  SetResSlope(Float_t resSlope)                {fResSlope=resSlope;}
  void  SetTimeWalkCenter(Float_t timeWalkCenter)    {fTimeWalkCenter=timeWalkCenter;}
  void  SetTimeWalkBoundary(Float_t timeWalkBoundary){fTimeWalkBoundary=timeWalkBoundary;}
  void  SetTimeWalkSlope(Float_t timeWalkSlope)      {fTimeWalkSlope=timeWalkSlope;}

  void  SetTimeDelayFlag(Int_t timeDelayFlag)        {fTimeDelayFlag=timeDelayFlag;}
  void  SetPulseHeightSlope(Float_t pulseHeightSlope){fPulseHeightSlope=pulseHeightSlope;}
  void  SetTimeDelaySlope(Float_t timeDelaySlope)    {fTimeDelaySlope=timeDelaySlope;}
  void  SetMinimumCharge(Float_t minimumCharge)      {fMinimumCharge=minimumCharge;}
  void  SetChargeSmearing(Float_t chargeSmearing)    {fChargeSmearing=chargeSmearing;}
  void  SetLogChargeSmearing(Float_t logChargeSmearing){fLogChargeSmearing=logChargeSmearing;}
  void  SetTimeSmearing(Float_t timeSmearing)        {fTimeSmearing=timeSmearing;}
  void  SetAverageTimeFlag(Int_t averageTimeFlag)    {fAverageTimeFlag=averageTimeFlag;}
  void  SetChargeFactorForMatching(Int_t chargeFactorForMatching){fChargeFactorForMatching=chargeFactorForMatching;}
  void  SetMatchingStyle(Int_t matchingStyle)            {fMatchingStyle=matchingStyle;}
  void  SetTrackingEfficiency(Float_t trackingEfficiency){fTrackingEfficiency=trackingEfficiency;}
  void  SetSigmavsp(Float_t sigmavsp)                    {fSigmavsp=sigmavsp;}
  void  SetSigmaZ(Float_t sigmaZ)                        {fSigmaZ=sigmaZ;}
  void  SetSigmarphi(Float_t sigmarphi)                  {fSigmarphi=sigmarphi;}
  void  SetSigmap(Float_t sigmap)                        {fSigmap=sigmap;}
  void  SetSigmaPhi(Float_t sigmaPhi)                    {fSigmaPhi=sigmaPhi;}
  void  SetSigmaTheta(Float_t sigmaTheta)                {fSigmaTheta=sigmaTheta;}
  void  SetNoise(Float_t noise)                          {fNoise=noise;}
  void  SetNoiseSlope(Float_t noiseSlope)                {fNoiseSlope=noiseSlope;}
  void  SetNoiseMeanTof(Float_t noiseTof)                    {fNoiseMeanTof=noiseTof;}
  void  SetField(Float_t field)                          {fField=field;}
  void  SetRadLenTPC(Float_t radLenTPC)                  {fRadLenTPC=radLenTPC;}
  void  SetCorrectionTRD(Float_t correctionTRD)          {fCorrectionTRD=correctionTRD;}
  void  SetLastTPCRow(Int_t lastTPCRow)                  {fLastTPCRow=lastTPCRow;}
  void  SetRadiusvtxBound(Float_t radiusvtxBound)        {fRadiusvtxBound=radiusvtxBound;}
  void  SetMaxTestTracks(Int_t maxTestTracks)            {fMaxTestTracks=maxTestTracks;}
  void  SetStep(Float_t step)                            {fStep=step;}
  void  SetMaxPixels(Int_t maxPixels)                    {fMaxPixels=maxPixels;}
  void  SetMaxAllTracks(Int_t maxAllTracks)              {fMaxAllTracks=maxAllTracks;}
  void  SetMaxTracks(Int_t maxTracks)                    {fMaxTracks=maxTracks;}
  void  SetMaxTOFHits(Int_t maxTOFHits)                  {fMaxTOFHits=maxTOFHits;}
  void  SetPBound(Float_t pBound)                        {fPBound=pBound;}

  Int_t    GetDbgFlag()          const {return fdbg;}
  Float_t  GetTimeResolution()   const {return fTimeResolution;}
  Float_t  GetPadefficiency()    const {return fpadefficiency;}
  Int_t    GetEdgeEffect()       const {return fEdgeEffect;}
  Int_t    GetEdgeTails()        const {return fEdgeTails;}
  Float_t  GetHparameter()       const {return fHparameter;}
  Float_t  GetH2parameter()      const {return fH2parameter;}
  Float_t  GetKparameter()       const {return fKparameter;}
  Float_t  GetK2parameter()      const {return fK2parameter;}
  Float_t  GetEffCenter()        const {return fEffCenter;}
  Float_t  GetEffBoundary()      const {return fEffBoundary;}
  Float_t  GetEff2Boundary()     const {return fEff2Boundary;}
  Float_t  GetEff3Boundary()     const {return fEff3Boundary;}
  Float_t  GetResCenter ()       const {return fResCenter;}
  Float_t  GetResBoundary()      const {return fResBoundary;}
  Float_t  GetResSlope()         const {return fResSlope;}
  Float_t  GetTimeWalkCenter()   const {return fTimeWalkCenter;}
  Float_t  GetTimeWalkBoundary() const {return fTimeWalkBoundary;}
  Float_t  GetTimeWalkSlope()    const {return fTimeWalkSlope;}
  Int_t    GetTimeDelayFlag()    const {return fTimeDelayFlag;}
  Float_t  GetPulseHeightSlope() const {return fPulseHeightSlope;}
  Float_t  GetTimeDelaySlope()   const {return fTimeDelaySlope;}
  Float_t  GetMinimumCharge()    const {return fMinimumCharge;}
  Float_t  GetChargeSmearing()   const {return fChargeSmearing;}
  Float_t  GetLogChargeSmearing()const {return fLogChargeSmearing;}
  Float_t  GetTimeSmearing()     const {return fTimeSmearing;}
  Int_t    GetAverageTimeFlag()  const {return fAverageTimeFlag;}
  Int_t    GetChargeFactorForMatching()const{return fChargeFactorForMatching;}
  Int_t    GetMatchingStyle()          const{return fMatchingStyle;}
  Float_t  GetTrackingEfficiency() const{return fTrackingEfficiency;}
  Float_t  GetSigmavsp()           const{return fSigmavsp;}
  Float_t  GetSigmaZ()             const{return fSigmaZ;}
  Float_t  GetSigmarphi()          const{return fSigmarphi;}
  Float_t  GetSigmap()             const{return fSigmap;}
  Float_t  GetSigmaPhi()           const{return fSigmaPhi;}
  Float_t  GetSigmaTheta()         const{return fSigmaTheta;}
  Float_t  GetNoise()              const{return fNoise;}
  Float_t  GetNoiseSlope()         const{return fNoiseSlope;}
  Float_t  GetNoiseMeanTof()       const{return fNoiseMeanTof;}
  Float_t  GetField()              const{return fField;}
  Float_t  GetRadLenTPC()          const{return fRadLenTPC;}
  Float_t  GetCorrectionTRD()      const{return fCorrectionTRD;}
  Int_t    GetLastTPCRow()         const{return fLastTPCRow;}
  Float_t  GetRadiusvtxBound()     const{return fRadiusvtxBound;}
  Int_t    GetMaxTestTracks()      const{return fMaxTestTracks;}
  Float_t  GetStep()               const{return fStep;}
  Int_t    GetMaxPixels()          const{return fMaxPixels;}
  Int_t    GetMaxAllTracks()       const{return fMaxAllTracks;}
  Int_t    GetMaxTracks()          const{return fMaxTracks;}
  Int_t    GetMaxTOFHits()         const{return fMaxTOFHits;}
  Float_t  GetPBound()             const{return fPBound;}

  virtual void PrintParameters() const ;
  virtual void Print(Option_t* option) const ;
  void  UseHitsFrom(const char * filename) ;
  Bool_t   operator == (const AliTOFReconstructioner & tofrec) const ;

private:
  TFile   *foutputfile;     //! pointer to output file
  TNtuple *foutputntuple;   //! pointer to output ntuple
  TF1     *fZnoise;         // pointer to formula giving the noise along z direction
  TF1     *ftail;           // pointer to formula for time with tail
  Int_t   fdbg;             // Flag for debug, 0 no debug, 1 debug
  Int_t   fNevents;         // Number of events to reconstruct  
  Int_t   fFirstEvent;      // First event to reconstruct
  Int_t   fLastEvent;       // Last event to reconstruct
  TString fRecFile;         // output file 
  TString fHeadersFile;     // input file
  // Intrisic MRPC time resolution and pad edge effect parameters
  Float_t fTimeResolution;  // time resolution of the MRPC (ns)
  Float_t fpadefficiency;   // intrinsic pad efficiency, used if fEdgeEffect==0
  Int_t   fEdgeEffect;      // edge effects option
  Int_t   fEdgeTails;       // edge tails option
  Float_t fHparameter;      // sensitive edge (to produce hits on the
                            // neighbouring pads) =0.7, new = 0.4 cm
  Float_t fH2parameter;     // parameter to fit the efficiency
  Float_t fKparameter;      // sensitive edge (going ahead towards the
                            // center no delay effects are suffered) =1.0, new = 0.5 cm
  Float_t fK2parameter;     // parameter to fit the efficiency
  // Pad Efficiency and Resolution parameters
  Float_t fEffCenter;       // efficiency in the central region of the pad
  Float_t fEffBoundary;     // efficiency at the boundary of the pad
  Float_t fEff2Boundary;    // efficiency value at H2parameter
  Float_t fEff3Boundary;    // efficiency value at K2parameter
  Float_t fResCenter;       // resolution (ps) in the central region of the pad
  Float_t fResBoundary;     // resolution (ps)  at the boundary of the pad
  Float_t fResSlope;        // slope (ps/K) for neighbouring pad
  // Time Walk parameters
  Float_t fTimeWalkCenter;  // time walk (ps) in the central region of the pad
  Float_t fTimeWalkBoundary;// time walk (ps) at the boundary of the pad
  Float_t fTimeWalkSlope;   // slope (ps/K) for neighbouring pad
  Int_t   fTimeDelayFlag;   // flag for delay due to the PulseHeightEffect
  Float_t fPulseHeightSlope;// It determines the charge amount induced
                            // due to edge effect, using the formula
                            // qInduced=exp(-PulseHeightSlope*x)
  Float_t fTimeDelaySlope;  // It determines the time delay. This is the slope
                            // in the T1-T2 vs log(q1/q2) plot
  // ADC-TDC correlation parameters
  Float_t fMinimumCharge;   // Minimum charge amount which could be induced
  Float_t fChargeSmearing;  // Smearing in charge in (q1/q2) vs x plot
  Float_t fLogChargeSmearing;// Smearing in log of charge ratio
  Float_t fTimeSmearing;    // Smearing in time in time vs log(q1/q2) plot
  Int_t   fAverageTimeFlag; // flag (see the setter for details)
  Int_t   fChargeFactorForMatching; // if set to 1, during matching procedure
                                     //  probe hits are weighted according to
                                     // the pulse height
  Int_t   fMatchingStyle;   // Matching style option (see setter for details)

  // TPC tracking parameters
  Float_t fTrackingEfficiency; //tracking efficiency in the TPC     0.88
  Float_t fSigmavsp;           //!=0 - sigmas depend on momentum, SIGMA VS P
                               // =0 - sigmas do not depend on momentum
  Float_t fSigmaZ;    //sigma(z) (cm)         0.044,  0.03/P(GeV/c) -> av.0.083, see AN-97-39 table 2 
  Float_t fSigmarphi; //sigma(R(phi)) (cm)    0.023,  0.015/P(GeV/c) -> av.0.041, see AN-97-39 table 2 
  Float_t fSigmap;    //sigma(delta(P)/P)     0.019,  0.01*(fabs((logP(GeV/c)+0.5)/0.7)**3+1.5) -> av.0.017
  Float_t fSigmaPhi;  //sigma(phi) (rad)      0.0050, 0.001*((1-logP(GeV/c))**3+0.3) for P<10 -> av.0.003 
  Float_t fSigmaTheta;//sigma(theta) (rad)    0.0035, 0.001*((1-logP(GeV/c))**3+0.3) for P<10 -> av.0.003

  // Parameters for additional noise hits
  Float_t fNoise;          //number of noise hits  6600/7800 with/without the holes,
                           //7800-holes*1200 for V3.02(TDR), 11000-holes*1200 for V3.04
                           // 10000-holes*? for V3.05 with the half z-length 370 cm, for 350 cm 
                           //                         it should be 5% less, 9500
                           //for the field 0.4 T:    ? /6333 with/without the holes
                           //for the field 0.4 T:  8400 for V3.05
                           //for PYTHIA p+p at 14 TeV: 26 for V3.05
  Float_t fNoiseSlope;     //slope parameter (ns) in the time distribution of the add. noise hits: ~exp(-tau/NOISETAU)
  Float_t fNoiseMeanTof;       // mean value of the time of flight for noise from outer regions (ns)
  Float_t fField;          //magnetic field (tesla), 0.2 tesla = 2 kilogauss
  Float_t fRadLenTPC;      //radiation length of the outer wall of TPC, 0.03+0.14/0.03 with/out TRD
  Float_t fCorrectionTRD;  //!=0 px, py on the last row of TPC are corrected
                           //    using x, y position on the last layer of TRD, see void correctionTRD()
                           //=0 without the correction
  Int_t   fLastTPCRow;     //the number of the last TPC row 111 for V3.05
  Float_t fRadiusvtxBound; //vertex radius (cm) for selected tracks
  Int_t   fMaxTestTracks;  //max.number of test tracks  20
                           //for PYTHIA p+p at 14 TeV: 500
  Float_t fStep;           //space step (cm) along circle, shuold be < the pixel width YP=ZAZOR+2*DY=12.3*0.05=0.615
  Int_t   fMaxPixels;      //max.number of pixels involved in the matching procedure, 70000
  Int_t   fMaxAllTracks;   //max.number of all tracks including the neutral ones, 65000
  Int_t   fMaxTracks;      //max.number of tracks selected for matching (hit on TPC with Rvtx<RVTXBOUND), 15000
  Int_t   fMaxTOFHits;     //max.number of TOF hits per event, 35000
  Float_t fPBound;         //tracks/hits with P(GeV/c)<PBOUND do not take into account (kinematical cut)

 protected:

  ClassDef(AliTOFReconstructioner,1)  // Task class for TOF reconstruction

};

#endif // AliTOFRECONSTRUCTIONER_H
