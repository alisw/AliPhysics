#ifndef AliRICHRecon_h
#define AliRICHRecon_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TTask.h>
#include "AliRICH.h"
#include <AliRun.h>

class AliRICH;

class AliRICHRecon : public TTask {
   public :
 
   AliRICHRecon(const char*, const char*);
   ~AliRICHRecon(){EndProcessEvent();}

   AliRICH* Rich() {return fRich;}
   void StartProcessEvent();
   void EndProcessEvent();
     
   void InitRecon();

   void PatRec();

   void Minimization();

   void FillHistograms();

   void FindThetaPhotonCerenkov();
   void FindAreaAndPortionOfRing();

   void FindEmissionPoint();

   Int_t PhotonInBand();
   Float_t PhotonPositionOnCathode();

   void FindPhotonAnglesInDRS();
   void FindPhiPoint();
   void FindThetaAtQuartz(Float_t ThetaCerenkov);

   Float_t FindMassOfParticle();

   Float_t Cerenkovangle(Float_t n, Float_t beta);

   void HoughResponse();
   void HoughFiltering(float HCS[]);
   void FlagPhotons();
   void FindWeightThetaCerenkov();

   void EstimationOfTheta();
   void FindIntersectionWithDetector();
   Int_t CheckDetectorAcceptance();

   void DrawRing();
   void DrawEvent(Int_t flag);

   void waiting();

   //////////////////////////////////////
   Float_t GetEventVertexZ() { return fEventVertZ;}
   Float_t GetEventMultiplicity() { return fEventMultiplicity;}
   Float_t GetPhotonEnergy() { return fPhotonEnergy;}
   Float_t GetFreonRefractiveIndex() { return fFreonRefractiveIndex;}
   Float_t GetQuartzRefractiveIndex() { return fQuartzRefractiveIndex;}
   Float_t GetGasRefractiveIndex() { return fGasRefractiveIndex;}

   Float_t GetFreonScaleFactor() { return fFreonScaleFactor;}

   Float_t GetEmissionPoint() { return fLengthEmissionPoint;}
   Float_t GetMassHypotesis() { return fMassHypotesis;}
   Float_t GetBetaOfParticle() { return fTrackBeta;}
   Float_t GetEntranceX() { return fXtoentr;}
   Float_t GetEntranceY() { return fYtoentr;}
   Float_t GetThetaCerenkov() { return fThetaCerenkov;}
   Float_t GetThetaPhotonCerenkov() { return fThetaPhotonCerenkov;}
   Float_t GetTrackMomentum() { return fTrackMomentum;}
   Float_t GetTrackEta() { return fTrackEta;}
   Float_t GetTrackTheta() { return fTrackTheta;}
   Float_t GetTrackPhi() { return fTrackPhi;}
   Float_t GetTrackPt() { return fTrackPt;}
   Int_t   GetTrackCharge() { return fTrackCharge;}
   Float_t GetTrackTPCLastZ() { return fTrackTPCLastZ;}
   Float_t GetMinDist() { return fMinDist;}
   Float_t GetXPointOnCathode() { return fPhotonLimitX;}
   Float_t GetYPointOnCathode() { return fPhotonLimitY;}


   Float_t GetThetaPhotonInDRS() { return fThetaPhotonInDRS;}
   Float_t GetPhiPhotonInDRS() { return fPhiPhotonInDRS;}
   Float_t GetThetaPhotonInTRS() { return fThetaPhotonInTRS;}
   Float_t GetPhiPhotonInTRS() { return fPhiPhotonInTRS;}

   Float_t GetThetaAtQuartz() { return fThetaAtQuartz;}

   Float_t GetPhiPoint(){ return fPhiPoint;}
   Float_t GetXCoordOfEmission() {return fXEmiss;}
   Float_t GetYCoordOfEmission() {return fYEmiss;}

   Float_t GetXInnerRing() {return fXInner;}
   Float_t GetYInnerRing() {return fYInner;}
   Float_t GetRadiusInnerRing() {return fInnerRadius;}

   Float_t GetXOuterRing() {return fXOuter;}
   Float_t GetYOuterRing() {return fYOuter;}
   Float_t GetRadiusOuterRing() {return fOuterRadius;}
   Float_t GetShiftX() { return fShiftX;}
   Float_t GetShiftY() { return fShiftY;}
   Float_t GetDetectorWhereX() { return fXcoord;}
   Float_t GetDetectorWhereY() { return fYcoord;}
   Float_t GetIntersectionX() {return fIntersectionX;}
   Float_t GetIntersectionY() {return fIntersectionY;}

   Float_t GetThetaOfRing() { return fThetaOfRing;}
   Float_t GetAreaOfRing() {return fAreaOfRing;}
   Float_t GetPortionOfRing() {return fPortionOfRing;}
   Float_t GetHoughArea() { return fHoughArea;}

   Int_t   GetPhotonsNumber() { return fPhotonsNumber;}
   Int_t   GetPhotonIndex() { return fPhotonIndex;}
   Float_t GetPhotonEta() { return fPhotonEta[fPhotonIndex];}
   Int_t   GetPhotonFlag() { return fPhotonFlag[fPhotonIndex];}
   Float_t GetPhotonWeight() { return fPhotonWeight[fPhotonIndex];}

   Float_t GetHoughRMS() { return fHoughRMS;}

   Int_t GetMipIndex() { return fMipIndex;}
   Int_t GetTrackIndex() { return fTrackIndex;}
   Float_t* GetCandidatePhotonX() { return fCandidatePhotonX;}
   Float_t* GetCandidatePhotonY() { return fCandidatePhotonY;}
   Int_t GetCandidatePhotonsNumber() { return fCandidatePhotonsNumber;}

   Int_t GetHoughPhotons() { return fHoughPhotons;}
   Float_t GetHoughPhotonsNorm() { return fHoughPhotonsNorm;}

   Float_t GetFittedTrackTheta() { return fFittedTrackTheta;}
   Float_t GetFittedTrackPhi() { return fFittedTrackPhi;}
   Float_t GetFittedThetaCerenkov() { return fFittedThetaCerenkov;}
   Int_t   GetFittedHoughPhotons() { return fFittedHoughPhotons;}
   Float_t GetEstimationOfTheta() { return fEstimationOfTheta;}
   Float_t GetEstimationOfThetaRMS() { return fEstimationOfThetaRMS;}

   void SetEventVertexZ(Float_t EventVertZ) { fEventVertZ = EventVertZ;}
   void SetEventMultiplicity(Float_t EventMultiplicity) { fEventMultiplicity = EventMultiplicity;}
   void SetPhotonEnergy(Float_t PhotonEnergy) { fPhotonEnergy = PhotonEnergy;} 
   void SetFreonRefractiveIndex() {fFreonRefractiveIndex = fFreonScaleFactor*(1.177+0.0172*fPhotonEnergy);}
   void SetQuartzRefractiveIndex() {fQuartzRefractiveIndex = sqrt(1+(46.411/(113.763556-TMath::Power(fPhotonEnergy,2)))+(228.71/(328.51563-TMath::Power(fPhotonEnergy,2))));}
   void SetGasRefractiveIndex() { fGasRefractiveIndex = 1.;}

   void SetFreonScaleFactor(Float_t FreonScaleFactor) {fFreonScaleFactor = FreonScaleFactor;}

   void SetEmissionPoint(Float_t LengthEmissionPoint) { fLengthEmissionPoint = LengthEmissionPoint;}
   void SetMassHypotesis(Float_t mass) {fMassHypotesis = mass;}

   void SetBetaOfParticle() { fTrackBeta = fTrackMomentum/sqrt(TMath::Power(fTrackMomentum,2)+TMath::Power(fMassHypotesis,2));}

   void SetEntranceX(Float_t Xtoentr) { fXtoentr = Xtoentr;}
   void SetEntranceY(Float_t Ytoentr) { fYtoentr = Ytoentr;}

   void SetThetaPhotonInTRS(Float_t Theta) {fThetaPhotonInTRS = Theta;}
   void SetPhiPhotonInTRS(Float_t Phi) {fPhiPhotonInTRS = Phi;}
   void SetThetaPhotonInDRS(Float_t Theta) {fThetaPhotonInDRS = Theta;}
   void SetPhiPhotonInDRS(Float_t Phi) {fPhiPhotonInDRS = Phi;}

   void SetThetaAtQuartz(Float_t ThetaAtQuartz) {fThetaAtQuartz = ThetaAtQuartz;}

   void SetPhiPoint(Float_t PhiPoint){ fPhiPoint = PhiPoint;}

   void SetXCoordOfEmission(Float_t XEmiss) {fXEmiss = XEmiss;}
   void SetYCoordOfEmission(Float_t YEmiss) {fYEmiss = YEmiss;}


   void SetXPointOnCathode(Float_t PhotonLimitX) { fPhotonLimitX = PhotonLimitX;}
   void SetYPointOnCathode(Float_t PhotonLimitY) { fPhotonLimitY = PhotonLimitY;}

   void SetXInnerRing(Float_t XInner) {fXInner = XInner;}
   void SetYInnerRing(Float_t YInner) {fYInner = YInner;}
   void SetRadiusInnerRing(Float_t InnerRadius) {fInnerRadius = InnerRadius;}

   void SetXOuterRing(Float_t XOuter) {fXOuter = XOuter;}
   void SetYOuterRing(Float_t YOuter) {fYOuter = YOuter;}
   void SetRadiusOuterRing(Float_t OuterRadius) {fOuterRadius = OuterRadius;}

   void SetThetaCerenkov(Float_t ThetaCer) {fThetaCerenkov = ThetaCer;}
   void SetThetaPhotonCerenkov(Float_t ThetaPhotCer) {fThetaPhotonCerenkov = ThetaPhotCer;}

   void SetTrackMomentum(Float_t TrackMomentum) {fTrackMomentum = TrackMomentum;}
   void SetTrackEta(Float_t TrackEta) {fTrackEta = TrackEta;}
   void SetTrackTheta(Float_t TrackTheta) { fTrackTheta = TrackTheta;}
   void SetTrackPhi(Float_t TrackPhi) { fTrackPhi = TrackPhi;}
   void SetTrackPt(Float_t TrackPt) { fTrackPt = TrackPt;}
   void SetTrackCharge(Int_t TrackCharge) { fTrackCharge = TrackCharge;}
   void SetTrackTPCLastZ(Float_t TrackTPCLastZ) { fTrackTPCLastZ = TrackTPCLastZ;}
   void SetMinDist(Float_t MinDist) { fMinDist = MinDist;}
   void SetShiftX(Float_t ShiftX) { fShiftX = ShiftX;}
   void SetShiftY(Float_t ShiftY) { fShiftY = ShiftY;}

   void SetDetectorWhereX(Float_t Xcoord) { fXcoord = Xcoord;}
   void SetDetectorWhereY(Float_t Ycoord) { fYcoord = Ycoord;}

   void SetIntersectionX(Float_t IntersectionX) { fIntersectionX = IntersectionX;}
   void SetIntersectionY(Float_t IntersectionY) { fIntersectionY = IntersectionY;}

   void SetThetaOfRing(Float_t ThetaOfRing) { fThetaOfRing = ThetaOfRing;}
   void SetAreaOfRing(Float_t AreaOfRing) { fAreaOfRing = AreaOfRing;}
   void SetPortionOfRing(Float_t PortionOfRing) { fPortionOfRing = PortionOfRing;}
   void SetHoughArea(Float_t HoughArea) { fHoughArea = HoughArea;}


   void SetPhotonsNumber(Int_t PhotonsNumber) { fPhotonsNumber = PhotonsNumber;}
   void SetPhotonIndex(Int_t PhotonIndex) { fPhotonIndex = PhotonIndex;}
   void SetPhotonEta(Float_t PhotonEta) { fPhotonEta[fPhotonIndex] = PhotonEta;}
   void SetPhotonFlag(Int_t PhotonFlag) { fPhotonFlag[fPhotonIndex] = PhotonFlag;}
   void SetPhotonWeight(Float_t PhotonWeight) { fPhotonWeight[fPhotonIndex] = PhotonWeight;}

   void SetHoughRMS(Float_t HoughRMS) { fHoughRMS = HoughRMS;}
   void SetMipIndex(Int_t MipIndex) { fMipIndex = MipIndex;}
   void SetTrackIndex(Int_t TrackIndex) { fTrackIndex = TrackIndex;}

   void SetCandidatePhotonX(Float_t *CandidatePhotonX) { fCandidatePhotonX = CandidatePhotonX;}
   void SetCandidatePhotonY(Float_t *CandidatePhotonY) { fCandidatePhotonY = CandidatePhotonY;}
   void SetCandidatePhotonsNumber(Int_t CandidatePhotonsNumber) { fCandidatePhotonsNumber = CandidatePhotonsNumber;}
   void SetHoughPhotons(Int_t HoughPhotons) { fHoughPhotons = HoughPhotons;}
   void SetHoughPhotonsNorm(Float_t HoughPhotonsNorm) { fHoughPhotonsNorm = HoughPhotonsNorm;}

   void SetFittedTrackTheta(Float_t FittedTrackTheta)    { fFittedTrackTheta = FittedTrackTheta;}
   void SetFittedTrackPhi(Float_t FittedTrackPhi)    { fFittedTrackPhi = FittedTrackPhi;}
   void SetFittedThetaCerenkov(Float_t FittedThetaCerenkov) { fFittedThetaCerenkov = FittedThetaCerenkov;}
   void SetFittedHoughPhotons(Int_t FittedHoughPhotons) { fFittedHoughPhotons = FittedHoughPhotons;}
   void SetEstimationOfTheta(Float_t EstimationOfTheta) { fEstimationOfTheta = EstimationOfTheta;}
   void SetEstimationOfThetaRMS(Float_t EstimationOfThetaRMS) { fEstimationOfThetaRMS = EstimationOfThetaRMS;}

   void FindBetaFromTheta(Float_t ThetaCerenkov) {fTrackBeta = 1/(fFreonRefractiveIndex*cos(ThetaCerenkov));}

   Float_t SnellAngle(Float_t n1, Float_t n2, Float_t theta1);

   Float_t FromEmissionToCathode();
   //////////////////////////////////////
   private:

   AliRICH* fRich;
       
   Float_t fEventVertZ;
   Float_t fEventMultiplicity;

   Float_t fTrackTheta;
   Float_t fTrackPhi;
   Float_t fTrackMomentum;
   Float_t fTrackEta;
   Float_t fTrackPt;
   Int_t   fTrackCharge;
   Float_t fTrackTPCLastZ;
   Float_t fMinDist;
   Float_t fTrackBeta;

   Float_t fXtoentr;
   Float_t fYtoentr;

   Float_t fThetaPhotonInTRS;
   Float_t fPhiPhotonInTRS;

   Float_t fThetaPhotonInDRS;
   Float_t fPhiPhotonInDRS;

   Float_t fThetaAtQuartz;
   Float_t fPhiPoint;

   Float_t fXEmiss;
   Float_t fYEmiss;

   Float_t fXInner;
   Float_t fYInner;
   Float_t fXOuter;
   Float_t fYOuter;
   Float_t fInnerRadius;
   Float_t fOuterRadius;

   Float_t fPhotonEnergy;
   Float_t fFreonRefractiveIndex;
   Float_t fQuartzRefractiveIndex;
   Float_t fGasRefractiveIndex;

   Float_t fFreonScaleFactor;

   Float_t fLengthEmissionPoint;

   Float_t fPhotonLimitX;
   Float_t fPhotonLimitY;
   Float_t fDistanceFromCluster;

   Float_t fMassHypotesis;

   Float_t fCerenkovAnglePad;

   Float_t fThetaPhotonCerenkov;

   Float_t fShiftX;
   Float_t fShiftY;

   Float_t fXcoord;
   Float_t fYcoord;

   Float_t fIntersectionX;
   Float_t fIntersectionY;

   Float_t fThetaOfRing;
   Float_t fAreaOfRing;
   Float_t fPortionOfRing;
   Float_t fHoughArea;

   Int_t fPhotonsNumber;
   Int_t fPhotonIndex;
   Float_t fPhotonEta[3000];
   Int_t fPhotonFlag[3000];
   Float_t fPhotonWeight[3000];

   Float_t fHoughRMS;

   Int_t fMipIndex;
   Int_t fTrackIndex;

   Float_t* fCandidatePhotonX;
   Float_t* fCandidatePhotonY;
   Int_t fCandidatePhotonsNumber;

   Int_t fHoughPhotons;
   Float_t fHoughPhotonsNorm;

   Float_t fFittedTrackTheta;
   Float_t fFittedTrackPhi;
   Float_t fFittedThetaCerenkov;
   Int_t   fFittedHoughPhotons;

   Float_t fEstimationOfTheta;
   Float_t fEstimationOfThetaRMS;

   public:

   Int_t   fNumEtaPhotons;                 // Number of photons
   Int_t   fEtaFlag[3000];                 // flag for good photons
   Float_t fEtaPhotons[3000];              // Cerenkov angle each photon
   Float_t fWeightPhotons[3000];           // weight for each photon
   Float_t fThetaCerenkov;                 // Theta angle for Hough
   Float_t fWeightThetaCerenkov;           // Theta Cerenkov angle weighted
   Float_t fThetaPeakPos;                  // Peak position


ClassDef(AliRICHRecon,0)

};
    
#endif // #ifdef AliRICHRecon_cxx

