#ifndef AliRICHRecon_h
#define AliRICHRecon_h

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// AliRICHRecon                                                         //
//                                                                      //
// RICH class to perfom pattern recognition based on Hough transfrom    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include <TTask.h>

class AliRICHHelix;

class AliRICHRecon : public TTask 
{
public : 
    AliRICHRecon(AliRICHHelix *pHelix,TClonesArray *pClusters,Int_t iMipId);
    virtual ~AliRICHRecon(){;}

  Double_t ThetaCerenkov();                   // it returns reconstructed Theta Cerenkov
  void FindThetaPhotonCerenkov();             //
  void FindAreaAndPortionOfRing();            //
  void FindEmissionPoint();                   //
  void FindPhotonAnglesInDRS();               //
  void FindPhiPoint();                        //
  void FindThetaAtQuartz(Float_t ThetaCer);   //
  void HoughResponse();                       //
  void HoughFiltering(float HCS[]);           //
  void FlagPhotons();                         //
  void FindWeightThetaCerenkov();             //
  void FindIntersectionWithDetector();        //
  Float_t Cerenkovangle(Float_t n, Float_t b);//
  Int_t   PhotonInBand();                     //
  Int_t   CheckDetectorAcceptance() const;    //
  Int_t   GetFittedHoughPhotons()                   const{ return fFittedHoughPhotons;}             //
  Int_t   GetPhotonFlag()                           const{ return fPhotonFlag[fPhotonIndex];}       //
  Int_t   GetTrackCharge()                          const{ return fTrackCharge;}                    //
  Int_t   GetPhotonsNumber()                        const{ return fPhotonsNumber;}                  //
  Int_t   GetPhotonIndex()                          const{ return fPhotonIndex;}                    //
  Int_t   GetMipIndex()                             const{ return fMipIndex;}                       //
  Int_t   GetTrackIndex()                           const{ return fTrackIndex;}                     //
  Int_t   GetCandidatePhotonsNumber()               const{ return fCandidatePhotonsNumber;}         //
  Int_t   GetHoughPhotons()                         const{ return fHoughPhotons;}                   //
  Float_t GetPhotonEnergy()                         const{ return fPhotonEnergy;}                   //
  Float_t GetFreonRefractiveIndex()                 const{ return fFreonRefractiveIndex;}           //
  Float_t GetQuartzRefractiveIndex()                const{ return fQuartzRefractiveIndex;}          //
  Float_t GetGasRefractiveIndex()                   const{ return fGasRefractiveIndex;}             //
  Float_t GetEmissionPoint()                        const{ return fLengthEmissionPoint;}            //
  Float_t GetMassHypotesis()                        const{ return fMassHypotesis;}                  //
  Float_t GetBetaOfParticle()                       const{ return fTrackBeta;}                      //
  Float_t GetEntranceX()                            const{ return fXtoentr;}                        //
  Float_t GetEntranceY()                            const{ return fYtoentr;}                        //
  Float_t GetThetaCerenkov()                        const{ return fThetaCerenkov;}                  //
  Float_t GetThetaPhotonCerenkov()                  const{ return fThetaPhotonCerenkov;}            //
  Float_t GetTrackTheta()                           const{ return fTrackTheta;}                     //
  Float_t GetTrackPhi()                             const{ return fTrackPhi;}                       //
  Float_t GetXPointOnCathode()                      const{ return fPhotonLimitX;}                   //
  Float_t GetYPointOnCathode()                      const{ return fPhotonLimitY;}                   //
  Float_t GetThetaPhotonInDRS()                     const{ return fThetaPhotonInDRS;}               //
  Float_t GetPhiPhotonInDRS()                       const{ return fPhiPhotonInDRS;}                 //
  Float_t GetThetaPhotonInTRS()                     const{ return fThetaPhotonInTRS;}               //
  Float_t GetPhiPhotonInTRS()                       const{ return fPhiPhotonInTRS;}                 //
  Float_t GetThetaAtQuartz()                        const{ return fThetaAtQuartz;}                  //
  Float_t GetPhiPoint()                             const{ return fPhiPoint;}                       //
  Float_t GetXCoordOfEmission()                     const{ return fXEmiss;}                         //
  Float_t GetYCoordOfEmission()                     const{ return fYEmiss;}                         //
  Float_t GetXInnerRing()                           const{ return fXInner;}                         //
  Float_t GetYInnerRing()                           const{ return fYInner;}                         //
  Float_t GetRadiusInnerRing()                      const{ return fInnerRadius;}                    //
  Float_t GetXOuterRing()                           const{ return fXOuter;}                         //
  Float_t GetYOuterRing()                           const{ return fYOuter;}                         //
  Float_t GetRadiusOuterRing()                      const{ return fOuterRadius;}                    //
  Float_t GetShiftX()                               const{ return fShiftX;}                         //
  Float_t GetShiftY()                               const{ return fShiftY;}                         //
  Float_t GetDetectorWhereX()                       const{ return fXcoord;}                         //
  Float_t GetDetectorWhereY()                       const{ return fYcoord;}                         //
  Float_t GetIntersectionX()                        const{ return fIntersectionX;}                  //
  Float_t GetIntersectionY()                        const{ return fIntersectionY;}                  //
  Float_t GetThetaOfRing()                          const{ return fThetaOfRing;}                    //
  Float_t GetAreaOfRing()                           const{ return fAreaOfRing;}                     //
  Float_t GetPortionOfRing()                        const{ return fPortionOfRing;}                  //
  Float_t GetHoughArea()                            const{ return fHoughArea;}                      //
  Float_t GetPhotonEta()                            const{ return fPhotonEta[fPhotonIndex];}        //
  Float_t GetPhotonWeight()                         const{ return fPhotonWeight[fPhotonIndex];}     //
  Float_t GetHoughRMS()                             const{ return fHoughRMS;}                       //
  Float_t GetFittedTrackTheta()                     const{ return fFittedTrackTheta;}               //
  Float_t GetFittedTrackPhi()                       const{ return fFittedTrackPhi;}                 //
  Float_t GetFittedThetaCerenkov()                  const{ return fFittedThetaCerenkov;}            //
  void SetPhotonEnergy(Float_t PhotonEnergy) { fPhotonEnergy = PhotonEnergy;}                       //
  void SetFreonRefractiveIndex() {fFreonRefractiveIndex = fFreonScaleFactor*(1.177+0.0172*fPhotonEnergy);}//
  void SetQuartzRefractiveIndex() {fQuartzRefractiveIndex = sqrt(1+(46.411/(113.763556-TMath::Power(fPhotonEnergy,2)))+(228.71/(328.51563-TMath::Power(fPhotonEnergy,2))));}//
  void SetGasRefractiveIndex() { fGasRefractiveIndex = 1.;}                                         //
  void SetFreonScaleFactor(Float_t FreonScaleFactor) {fFreonScaleFactor = FreonScaleFactor;}        //
  void SetEmissionPoint(Float_t LengthEmissionPoint) { fLengthEmissionPoint = LengthEmissionPoint;} //
  void SetEntranceX(Float_t Xtoentr) { fXtoentr = Xtoentr;}                                         //
  void SetEntranceY(Float_t Ytoentr) { fYtoentr = Ytoentr;}                                         //
  void SetThetaPhotonInTRS(Float_t Theta) {fThetaPhotonInTRS = Theta;}                              //
  void SetPhiPhotonInTRS(Float_t Phi) {fPhiPhotonInTRS = Phi;}                                      //
  void SetThetaPhotonInDRS(Float_t Theta) {fThetaPhotonInDRS = Theta;}                              //
  void SetPhiPhotonInDRS(Float_t Phi) {fPhiPhotonInDRS = Phi;}                                      //
  void SetThetaAtQuartz(Float_t ThetaAtQuartz) {fThetaAtQuartz = ThetaAtQuartz;}                    //
  void SetPhiPoint(Float_t PhiPoint){ fPhiPoint = PhiPoint;}                                        //
  void SetXCoordOfEmission(Float_t XEmiss) {fXEmiss = XEmiss;}                                      //
  void SetYCoordOfEmission(Float_t YEmiss) {fYEmiss = YEmiss;}                                      //
  void SetXPointOnCathode(Float_t PhotonLimitX) { fPhotonLimitX = PhotonLimitX;}                    //
  void SetYPointOnCathode(Float_t PhotonLimitY) { fPhotonLimitY = PhotonLimitY;}                    //
  void SetXInnerRing(Float_t XInner) {fXInner = XInner;}                                            //
  void SetYInnerRing(Float_t YInner) {fYInner = YInner;}                                            //
  void SetRadiusInnerRing(Float_t InnerRadius) {fInnerRadius = InnerRadius;}                        //
  void SetXOuterRing(Float_t XOuter) {fXOuter = XOuter;}                                            //
  void SetYOuterRing(Float_t YOuter) {fYOuter = YOuter;}                                            //
  void SetRadiusOuterRing(Float_t OuterRadius) {fOuterRadius = OuterRadius;}                        //
  void SetThetaCerenkov(Float_t ThetaCer) {fThetaCerenkov = ThetaCer;}                              //
  void SetThetaPhotonCerenkov(Float_t ThetaPhotCer) {fThetaPhotonCerenkov = ThetaPhotCer;}          //
  void SetTrackTheta(Float_t TrackTheta) { fTrackTheta = TrackTheta;}                               //
  void SetTrackPhi(Float_t TrackPhi) { fTrackPhi = TrackPhi;}                                       //
  void SetTrackCharge(Int_t TrackCharge) { fTrackCharge = TrackCharge;}                             //
  void SetShiftX(Float_t ShiftX) { fShiftX = ShiftX;}                                               //
  void SetShiftY(Float_t ShiftY) { fShiftY = ShiftY;}                                               //
  void SetDetectorWhereX(Float_t Xcoord) { fXcoord = Xcoord;}                                       //
  void SetDetectorWhereY(Float_t Ycoord) { fYcoord = Ycoord;}                                       //
  void SetIntersectionX(Float_t IntersectionX) { fIntersectionX = IntersectionX;}                   //
  void SetIntersectionY(Float_t IntersectionY) { fIntersectionY = IntersectionY;}                   //
  void SetThetaOfRing(Float_t ThetaOfRing) { fThetaOfRing = ThetaOfRing;}                           //
  void SetAreaOfRing(Float_t AreaOfRing) { fAreaOfRing = AreaOfRing;}                               //
  void SetPortionOfRing(Float_t PortionOfRing) { fPortionOfRing = PortionOfRing;}                   //
  void SetHoughArea(Float_t HoughArea) { fHoughArea = HoughArea;}                                   //
  void SetPhotonsNumber(Int_t PhotonsNumber) { fPhotonsNumber = PhotonsNumber;}                     //
  void SetPhotonIndex(Int_t PhotonIndex) { fPhotonIndex = PhotonIndex;}                             //
  void SetPhotonEta(Float_t PhotonEta) { fPhotonEta[fPhotonIndex] = PhotonEta;}                     //
  void SetPhotonFlag(Int_t PhotonFlag) { fPhotonFlag[fPhotonIndex] = PhotonFlag;}                   //
  void SetPhotonWeight(Float_t PhotonWeight) { fPhotonWeight[fPhotonIndex] = PhotonWeight;}         //
  void SetHoughRMS(Float_t HoughRMS) { fHoughRMS = HoughRMS;}                                       //
  void SetMipIndex(Int_t MipIndex) { fMipIndex = MipIndex;}                                         //
  void SetTrackIndex(Int_t TrackIndex) { fTrackIndex = TrackIndex;}                                 //
  void SetHoughPhotons(Int_t HoughPhotons) { fHoughPhotons = HoughPhotons;}                         //
  void SetHoughPhotonsNorm(Float_t HoughPhotonsNorm) { fHoughPhotonsNorm = HoughPhotonsNorm;}       //
  void SetFittedTrackTheta(Float_t FittedTrackTheta)    { fFittedTrackTheta = FittedTrackTheta;}    //
  void SetFittedTrackPhi(Float_t FittedTrackPhi)    { fFittedTrackPhi = FittedTrackPhi;}            //
  void SetFittedThetaCerenkov(Float_t FittedThetaCerenkov) { fFittedThetaCerenkov = FittedThetaCerenkov;}//
  void SetFittedHoughPhotons(Int_t FittedHoughPhotons) { fFittedHoughPhotons = FittedHoughPhotons;} //
  void FindBetaFromTheta(Float_t ThetaCerenkov) {fTrackBeta = 1/(fFreonRefractiveIndex*cos(ThetaCerenkov));}//
  Float_t SnellAngle(Float_t n1, Float_t n2, Float_t theta1);                                       //
  Float_t FromEmissionToCathode();                                                                  //

protected:
  TClonesArray *fpClusters;                   // poiter to clusters
  Int_t   fTrackCharge;                       // charge track
  Int_t fMipIndex;                            // index for Mip
  Int_t fTrackIndex;                          // index for track
  Int_t fPhotonsNumber;                       // Number of photons candidate
  Int_t fPhotonIndex;                         // index of photons
  Int_t fPhotonFlag[3000];                    // flag for photons
  Int_t fCandidatePhotonsNumber;              // number of candidate photons
  Int_t fHoughPhotons;                        // n. photons after Hough
  Int_t   fFittedHoughPhotons;                // n. photons after Hough and after minimization

  Float_t fTrackTheta;                        // Theta of track at RICH
  Float_t fTrackPhi;                          // Phi of track at RICH
  Float_t fMinDist;                           // min distance between extrapolated track and MIP
  Float_t fTrackBeta;                         // beta of the track
  Float_t fXtoentr;                           // X entrance to RICH
  Float_t fYtoentr;                           // Y entrance to RICH
  Float_t fThetaPhotonInTRS;                  // Theta of photon in the Track Reference System (TRS)
  Float_t fPhiPhotonInTRS;                    // Phi of photon in TRS
  Float_t fThetaPhotonInDRS;                  // Theta of photon in Detector Reference System (DRS)
  Float_t fPhiPhotonInDRS;                    // Phi of photon in DRS
  Float_t fThetaAtQuartz;                     // Theta at the quartz entrance
  Float_t fPhiPoint;                          // phi of ring point
  Float_t fXEmiss;                            //  x emission
  Float_t fYEmiss;                            //  y emission
  Float_t fXInner;                            // X inner ring
  Float_t fYInner;                            // Y inner ring
  Float_t fXOuter;                            // X outer ring
  Float_t fYOuter;                            // Y outer ring
  Float_t fInnerRadius;                       // inner radius
  Float_t fOuterRadius;                       // outer radius
  Float_t fPhotonEnergy;                      // photon energy
  Float_t fFreonRefractiveIndex;              // n freon
  Float_t fQuartzRefractiveIndex;             // n quartz
  Float_t fGasRefractiveIndex;                // n gas
  Float_t fFreonScaleFactor;                  // scale factor for n freon
  Float_t fLengthEmissionPoint;               // lenght of emmission point
  Float_t fPhotonLimitX;                      // X phys limit for photon
  Float_t fPhotonLimitY;                      // Y phys limit for photon 
  Float_t fDistanceFromCluster;               // distance from cluster
  Float_t fCerenkovAnglePad;                  // cherenkov angle of pad
  Float_t fThetaPhotonCerenkov;               // theta cerenkov for photon
  Float_t fShiftX;                            // x shift to entrance in radiator
  Float_t fShiftY;                            // y shift to entrance in radiator
  Float_t fXcoord;                            // ..
  Float_t fYcoord;                            // ..
  Float_t fIntersectionX;                     // ..
  Float_t fIntersectionY;                     // ..
  Float_t fMassHypotesis;                     //
  Float_t fThetaOfRing;                       // theta of ring
  Float_t fAreaOfRing;                        // area of the ring
  Float_t fPortionOfRing;                     // fraction of the accepted ring
  Float_t fHoughArea;                         // area Hough
  Float_t fPhotonEta[3000];                   // theta cerenkov of photon candidates
  Float_t fPhotonWeight[3000];                // weigth
  Float_t fHoughRMS;                          // rms Hough
  Float_t* fCandidatePhotonX;                 // x photon candidates
  Float_t* fCandidatePhotonY;                 // y photon candidates
  Float_t fHoughPhotonsNorm;                  // n. photons norm.
  Float_t fFittedTrackTheta;                  // theta track after minim.
  Float_t fFittedTrackPhi;                    // phi track after minim.
  Float_t fFittedThetaCerenkov;               // thetacerenkov after minim.
  Int_t   fThetaBin;                          // bin in theta
  Float_t fThetaMin,fThetaMax;                // min max
  Float_t fXmin,fXmax,fYmin,fYmax;            // xy min max
  Int_t   fNrings;                            //current number of reconstructed rings
  Bool_t  fIsWEIGHT;                          // flag to consider weight procedure
  Bool_t  fIsBACKGROUND;                      // flag to simulate bkg
  Float_t fRadiatorWidth;                     // radiator width
  Float_t fQuartzWidth;                       // quartz width
  Float_t fGapWidth;                          // gap width
  Float_t fDTheta;                            // Step for sliding window
  Float_t fWindowWidth;                       // Hough width of sliding window
  
  Int_t   fNumEtaPhotons;                     // Number of photons
  Int_t   fEtaFlag[3000];                     // flag for good photons
  Float_t fEtaPhotons[3000];                  // Cerenkov angle each photon
  Float_t fWeightPhotons[3000];               // weight for each photon
  Float_t fThetaCerenkov;                     // Theta angle for Hough
  Float_t fWeightThetaCerenkov;               // Theta Cerenkov angle weighted
  Float_t fThetaPeakPos;                      // Peak position


  ClassDef(AliRICHRecon,0)
};
    
#endif // #ifdef AliRICHRecon_cxx

