/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// Implementation for correlation analysis

#include "AliJJtCorrelations.h"

#include "../AliJBaseTrack.h"
#include "../AliJCard.h"


AliJJtCorrelations::AliJJtCorrelations( AliJCard *cardIn, AliJJtHistograms *histosIn) :
  AliJCorrelationInterface(),
  fcard(cardIn),
  fhistos(histosIn),
  fAcceptanceCorrection(0x0),
  fnReal(0),
  fnMix(0),
  fsamplingMethod(0),    // flat by default
  fmaxEtaRange(0),
  frandom(0x0),
  fptt(0),
  fpta(0),
  fTrackPairEfficiency(0),
  fpttBin(0),
  fptaBin(0),
  fPhiTrigger(0),
  fPhiAssoc(0),
  fEtaTrigger(0),
  fEtaAssoc(0),
  fDeltaPhi(0),
  fDeltaPhiPiPi(0),
  fDeltaEta(0),
  fXlong(0),
  fZBin(0),
  fNearSide(true),
  fNearSide3D(true),
  fEtaGapBin(0),
  fPhiGapBinNear(0),
  fRGapBinNear(0),
  fCentralityBin(0),
  fXlongBin(0),
  fUseKlongBins(false),
  fIsLikeSign(false),
  fUseZVertexBinsAcceptance(false)
{
  // constructor
  fmaxEtaRange = fcard->Get("EtaRange");
  if(fcard->Get("EnableKlongBins") == 1) fUseKlongBins = true;
  
  frandom = new TRandom3(); // Random number generator for background randomization
  frandom->SetSeed(0);
  
}

AliJJtCorrelations::AliJJtCorrelations() :
  AliJCorrelationInterface(),
  fcard(0x0),
  fhistos(0x0),
  fAcceptanceCorrection(0x0),
  fnReal(0),
  fnMix(0),
  fsamplingMethod(0),    // flat by default
  fmaxEtaRange(0),
  frandom(0x0),
  fptt(0),
  fpta(0),
  fTrackPairEfficiency(0),
  fpttBin(0),
  fptaBin(0),
  fPhiTrigger(0),
  fPhiAssoc(0),
  fEtaTrigger(0),
  fEtaAssoc(0),
  fDeltaPhi(0),
  fDeltaPhiPiPi(0),
  fDeltaEta(0),
  fXlong(0),
  fZBin(0),
  fNearSide(true),
  fNearSide3D(true),
  fEtaGapBin(0),
  fPhiGapBinNear(0),
  fRGapBinNear(0),
  fCentralityBin(0),
  fXlongBin(0),
  fUseKlongBins(false),
  fIsLikeSign(false),
  fUseZVertexBinsAcceptance(false)
{
  // default constructor
}

AliJJtCorrelations::AliJJtCorrelations(const AliJJtCorrelations& in) :
  fcard(in.fcard),
  fhistos(in.fhistos),
  fAcceptanceCorrection(in.fAcceptanceCorrection),
  fnReal(in.fnReal),
  fnMix(in.fnMix),
  fsamplingMethod(in.fsamplingMethod),
  fmaxEtaRange(in.fmaxEtaRange),
  frandom(in.frandom),
  fptt(in.fptt),
  fpta(in.fpta),
  fTrackPairEfficiency(in.fTrackPairEfficiency),
  fpttBin(in.fpttBin),
  fptaBin(in.fptaBin),
  fPhiTrigger(in.fPhiTrigger),
  fPhiAssoc(in.fPhiAssoc),
  fEtaTrigger(in.fEtaTrigger),
  fEtaAssoc(in.fEtaAssoc),
  fDeltaPhi(in.fDeltaPhi),
  fDeltaPhiPiPi(in.fDeltaPhiPiPi),
  fDeltaEta(in.fDeltaEta),
  fXlong(in.fXlong),
  fZBin(in.fZBin),
  fNearSide(in.fNearSide),
  fNearSide3D(in.fNearSide3D),
  fEtaGapBin(in.fEtaGapBin),
  fPhiGapBinNear(in.fPhiGapBinNear),
  fRGapBinNear(in.fRGapBinNear),
  fCentralityBin(in.fCentralityBin),
  fXlongBin(in.fXlongBin),
  fUseKlongBins(in.fUseKlongBins),
  fIsLikeSign(in.fIsLikeSign),
  fUseZVertexBinsAcceptance(in.fUseZVertexBinsAcceptance)
{
  // The pointers to card and histos are just copied. I think this is safe, since they are not created by
  // AliJJtCorrelations and thus should not disappear if the AliJCorrelation managing them is destroyed.
  
  frandom = new TRandom3(); // // Random number generator for background randomization
  frandom->SetSeed(0);
}

AliJJtCorrelations& AliJJtCorrelations::operator=(const AliJJtCorrelations& in){
  // Assingment operator
  
  if (&in==this) return *this;
  
  fptt = in.fptt;
  fpta = in.fpta;
  fTrackPairEfficiency = in.fTrackPairEfficiency;
  fpttBin = in.fpttBin;
  fptaBin = in.fptaBin;
  fPhiTrigger = in.fPhiTrigger;
  fPhiAssoc = in.fPhiAssoc;
  fEtaTrigger = in.fEtaTrigger;
  fEtaAssoc = in.fEtaAssoc;
  fDeltaPhi = in.fDeltaPhi;
  fDeltaPhiPiPi = in.fDeltaPhiPiPi;
  fDeltaEta = in.fDeltaEta;
  fXlong = in.fXlong;
  fZBin = in.fZBin;
  fNearSide = in.fNearSide;
  fNearSide3D = in.fNearSide3D;
  fEtaGapBin = in.fEtaGapBin;
  fPhiGapBinNear = in.fPhiGapBinNear;
  fRGapBinNear = in.fRGapBinNear;
  fCentralityBin = in.fCentralityBin;
  fXlongBin = in.fXlongBin;
  fUseKlongBins = in.fUseKlongBins;
  fIsLikeSign = in.fIsLikeSign;
  fnReal = in.fnReal;
  fnMix = in.fnMix;
  fsamplingMethod = in.fsamplingMethod;
  fmaxEtaRange = in.fmaxEtaRange;
  fUseZVertexBinsAcceptance = in.fUseZVertexBinsAcceptance;
  
  // The pointers to card and histos are just copied. I think this is safe, since they are not created by
  // AliJJtCorrelations and thus should not disappear if the AliJCorrelation managing them is destroyed.
  fcard = in.fcard;
  fhistos = in.fhistos;
  fAcceptanceCorrection = in.fAcceptanceCorrection;

  frandom = new TRandom3(); // Random number generator for background randomization
  frandom->SetSeed(0);
  
  return *this;
  // copy constructor
}

AliJJtCorrelations::~AliJJtCorrelations(){
  // destructor
  
  delete frandom;

}


/*
 * Histogram filled based on correlation type. Only calls the main histogram filler in case of correct correlation type.
 *
 *  corrFillType cFTyp = corraletion type. Histograms only filled for kAzimuthFill
 *  fillType fTyp = real or mixed events
 *  int CentBin = centrality bin index
 *  int zBin = z-vertex bin index
 *  AliJBaseTrack *ftk1 = Track for trigger particle
 *  AliJBaseTrack *ftk2 = Track for associated particle
 */
void AliJJtCorrelations::FillHisto(corrFillType cFTyp, fillType fTyp, int cBin, int zBin, AliJBaseTrack *ftk1, AliJBaseTrack *ftk2){
  // histo filler
  if( cFTyp == kAzimuthFill ){
    FillCorrelationHistograms( fTyp, cBin, zBin, ftk1, ftk2);
  }
  
}

/*
 * Main correlation histogram filler. Reads basic info from tracks and calls other methods that fill the histograms.
 *
 *  fillType fTyp = real or mixed events
 *  int CentBin = centrality bin index
 *  int zBin = z-vertex bin index
 *  AliJBaseTrack *ftk1 = Track for trigger particle
 *  AliJBaseTrack *ftk2 = Track for associated particle
 */
void AliJJtCorrelations::FillCorrelationHistograms(fillType fTyp, int CentBin, int zBin, AliJBaseTrack *ftk1, AliJBaseTrack *ftk2)
{
  // histo filler
  bool twoTracks = false;
  if(ftk1->GetParticleType()==kJHadron && ftk2->GetParticleType()==kJHadron) twoTracks =true;
  
  //double-counting check
  if(fTyp == kReal && twoTracks && ftk1->GetID()==ftk2->GetID()) return;
  
  // Check the signs of the paired particles
  fIsLikeSign = false;
  if(ftk1->GetCharge() > 0 && ftk2->GetCharge() > 0) fIsLikeSign = true;
  if(ftk1->GetCharge() < 0 && ftk2->GetCharge() < 0) fIsLikeSign = true;
  
  //----------------------------------------------------------------
  fptt = ftk1->Pt();
  fpta = ftk2->Pt();
  
  fTrackPairEfficiency = 1./( ftk1->GetTrackEff() * ftk2->GetTrackEff() );
  
  TLorentzVector vTrigger = ftk1->GetLorentzVector(), vAssoc = ftk2->GetLorentzVector();    // Lorentz vectors for trigger and associated particles
  
  fpttBin       = ftk1->GetTriggBin();
  fptaBin       = ftk2->GetAssocBin();
  fPhiTrigger   = ftk1->Phi();
  fPhiAssoc     = ftk2->Phi();
  fDeltaPhiPiPi = atan2(sin(fPhiTrigger-fPhiAssoc), cos(fPhiTrigger-fPhiAssoc));
  fDeltaEta     = ftk1->Eta() - ftk2->Eta();
  fEtaTrigger   = ftk1->Eta();
  fEtaAssoc     = ftk2->Eta();
  
  fNearSide     = cos(fPhiTrigger-fPhiAssoc) > 0 ? true : false;  // Traditional near side definition using deltaPhi
  fNearSide3D   = vTrigger.Vect().Dot(vAssoc.Vect()) > 0 ? true : false; // Near side definition using half ball around the trigger

  fEtaGapBin = fcard->GetBin( kEtaGapType, fabs(fDeltaEta));
  fPhiGapBinNear = fcard->GetBin( kEtaGapType, fabs(fDeltaPhiPiPi) );
  fRGapBinNear   = fcard->GetBin( kRGapType, fabs(ftk1->DeltaR(*ftk2)));
  fCentralityBin = CentBin;
  fZBin = zBin;
  
  fXlong = vTrigger.Vect().Dot(vAssoc.Vect())/pow(vTrigger.P(),2);
  fXlongBin = fcard->GetBin(kXeType, TMath::Abs(fXlong));
  
  // Check that the bins make sense. If not, do not fill anything.
  if(fpttBin<0 || fptaBin<0 || fEtaGapBin<0 ){
    cout<<"Error in FillJtHistograms: some pT or eta out of bin. pttBin="<<fpttBin<<" pTaBin="<<fptaBin <<" etaGapBin="<< fEtaGapBin << endl;
    ftk1->Print();
    ftk2->Print();
    exit(-1);
    //return;
  }
  
  
  // ===================================================================
  // =====================  Fill Histograms  ===========================
  // ===================================================================
  
  bool fill2DBackgroundQualityControlHistograms = false;  // Choose whether to fill the DeltaPhi DeltaEta histograms for jT background
  
  // If quality control level is high enough, fill 2D quality control histograms
  if(fcard->Get("QualityControlLevel")>1) fill2DBackgroundQualityControlHistograms = true;
  
  FillJtHistograms(fTyp, ftk1, ftk2, fill2DBackgroundQualityControlHistograms);  // Fill the jT and pout histograms togerher with some background quality assurance histograms
  FillDeltaEtaHistograms(fTyp);  // Fill all the delta eta histograms
  if(fcard->Get("EnableDeltaEtaDeltaPhiHistograms")==1) FillDeltaEtaDeltaPhiHistograms(fTyp);  // Fill the 2D correlation functions
  FillPtaHistograms(fTyp);  // Fill various pTa histograms
  
}



/*
 *  Find out the bin index based on the associated binning type and trigger and asociated vectors
 *
 *  int assocType = Associated particle binning type. 0 = klong, 1 = xlong, 2 = pta
 *  TLorentzVector *vTrigger = vector for trigger particle
 *  TLorentzVector *vAssoc = vector for associated particle
 *
 */
int AliJJtCorrelations::GetBinIndex(int assocType, TLorentzVector *vTrigger, TLorentzVector *vAssoc)
{
 
  // Find out the correct bin to fill the value. Require that we are in the near side.
  // Notice, that the 3D near side requirement can be written as a dot product between trigger and associated particles
  // If the dot product is negative, the associated particle must be on the away side
  // For the pTa bins a traditional near side definition is used
  int iBin = -1;
  if(assocType == 0){ // klong bin
    iBin = fcard->GetBin(kLongType, vTrigger->Vect().Dot(vAssoc->Vect())/vTrigger->P());
  } else if (assocType == 1){ // xlong bin
    iBin = fcard->GetBin(kXeType, vTrigger->Vect().Dot(vAssoc->Vect())/pow(vTrigger->P(),2) );
  } else { // pta bin
    if(cos(vTrigger->Phi()-vAssoc->Phi()) > 0) iBin = fptaBin;
  }
  
  return iBin;
}

/*
 * Method for filling the jT and invariant mass histograms
 *
 *  fillType fTyp = Tell what kind of particles are filled (real events vs. mixed events)
 *  int assocType = Associated particle binning type. 0 = klong, 1 = xlong, 2 = pta
 *  TLorentzVector *vTrigger = vector for trigger particle
 *  TLorentzVector *vAssoc = vector for associated particle
 *  AliJTH1D &hDistribution = main jT histogram
 *  AliJTH1D &hDistributionLikeSign = jT histogram for the like sign pairs
 *  AliJTH1D &hDistributionUnlikeSign = jT histogram for the unlike sign pairs
 *  AliJTH1D &hInvariantMass = invariant mass histogram
 *  AliJTH1D &hInvariantMassLikeSign = invariant mass histogram for the like sign pairs
 *  AliJTH1D &hInvariantMassUnlikeSign = invariant mass histogram for the unlike sign pairs
 *
 */
void AliJJtCorrelations::FillJtDistributionHistograms(fillType fTyp, int assocType, TLorentzVector *vTrigger, TLorentzVector *vAssoc, AliJTH1D &hDistribution, AliJTH1D &hDistributionLikeSign, AliJTH1D &hDistributionUnlikeSign, AliJTH1D &hInvariantMass, AliJTH1D &hInvariantMassLikeSign, AliJTH1D &hInvariantMassUnlikeSign)
{
  
  // Find out the correct bin to fill the value
  int iBin = GetBinIndex(assocType, vTrigger, vAssoc);
  
  // Calculate jT, invariant mass and the correction factor for jT
  double jt = vAssoc->Perp(vTrigger->Vect());
  double geometricAcceptanceCorrection = 0;
  if(fUseZVertexBinsAcceptance){
    geometricAcceptanceCorrection = fAcceptanceCorrection->GetAcceptanceCorrection(assocType,fsamplingMethod,fDeltaEta,fDeltaPhiPiPi,fCentralityBin,fZBin,fpttBin);
  } else {
    geometricAcceptanceCorrection = fAcceptanceCorrection->GetAcceptanceCorrection(assocType,fsamplingMethod,fDeltaEta,fDeltaPhiPiPi,fCentralityBin,fpttBin);
  }
  double weight = jt > 1e-3 ? geometricAcceptanceCorrection * fTrackPairEfficiency/jt : 0;
  double invariantMass = sqrt(2*(vTrigger->P()*vAssoc->P()-vTrigger->Vect().Dot(vAssoc->Vect())));
  
  // Fill the calculated values to histograms
  if(iBin>=0){
    hDistribution[fTyp][fCentralityBin][fpttBin][iBin]->Fill(jt, weight);
    hInvariantMass[fTyp][fCentralityBin][fpttBin][iBin]->Fill(invariantMass);
    
    // Fill also the signed distributions
    if(fIsLikeSign){
      hDistributionLikeSign[fTyp][fCentralityBin][fpttBin][iBin]->Fill(jt, weight);
      hInvariantMassLikeSign[fTyp][fCentralityBin][fpttBin][iBin]->Fill(invariantMass);
    } else {
      hDistributionUnlikeSign[fTyp][fCentralityBin][fpttBin][iBin]->Fill(jt, weight);
      hInvariantMassUnlikeSign[fTyp][fCentralityBin][fpttBin][iBin]->Fill(invariantMass);
    }
  }
  
}

/*
 * Method for filling all the necessary background histograms for a background estimation method
 *
 *  int assocType = Associated particle binning type. 0 = klong, 1 = xlong, 2 = pta
 *  int gapType = Used gap type for background. 0 = eta gap, 1 = R gap, 2 = phi gap
 *  TLorentzVector *vTrigger = trigger particle used for the background estimation
 *  TLorentzVector *vAssoc = associated particle used for the background estimation
 *  AliJTH1D &hBackground = main background histogram
 *  AliJTH1D &hBackgroundLikeSign = background histogram for the like sign pairs
 *  AliJTH1D &hBackgroundUnlikeSign = background histogram for the unlike sign pairs
 *  AliJTH1D &hPtAssoc = associated particle pT distribution
 *  AliJTH2D &hBackground2D = two dimensional quality control histogram
 *  bool fill2DBackground = Switch for quality control histograms. False = do not fill, true = fill
 *
 */
void AliJJtCorrelations::FillJtBackgroundHistograms(int assocType, int gapType, TLorentzVector *vTrigger, TLorentzVector *vAssoc, AliJTH1D &hBackground, AliJTH1D &hBackgroundLikeSign, AliJTH1D &hBackgroundUnlikeSign, AliJTH1D &hPtAssoc, AliJTH2D &hBackground2D, bool fill2DBackground)
{
  
  // Choose the correct bin based on associated binning type
  int iBin = GetBinIndex(assocType, vTrigger, vAssoc);
  
  // Choose the correct gap based on the gap type
  int maxGap = -1;
  if(gapType == 0){ // eta gap
    maxGap = fEtaGapBin;
  } else if (gapType == 1){ // R gap
    maxGap = fRGapBinNear;
  } else { // phi gap
    maxGap = fPhiGapBinNear;
  }
  
  // Calculate the jt for the background pair
  double jt = vAssoc->Perp(vTrigger->Vect());
  
  // Find the acceptance correction for the pair
  double dEtaRndm = vTrigger->Eta() - vAssoc->Eta();
  
  // Find the phi difference in the interval ]-Pi,Pi]
  double dPhiRndm = atan2(sin(vTrigger->Phi()-vAssoc->Phi()), cos(vTrigger->Phi()-vAssoc->Phi()));
  
  // Find the acceptance correction for the found particle pair
  double geoAccCorrRndm = 0;
  if(fUseZVertexBinsAcceptance){
    geoAccCorrRndm = fAcceptanceCorrection->GetAcceptanceCorrection(assocType,fsamplingMethod,dEtaRndm,dPhiRndm,fCentralityBin,fZBin,fpttBin);
  } else {
    geoAccCorrRndm = fAcceptanceCorrection->GetAcceptanceCorrection(assocType,fsamplingMethod,dEtaRndm,dPhiRndm,fCentralityBin,fpttBin);
  }
  
  // Fill the background histograms
  if(iBin>=0){ //
    if(jt>1e-3) { //here we used and BgFidCut
      for(int iGap=0; iGap<=maxGap; iGap++){
        hBackground[fCentralityBin][iGap][fpttBin][iBin]->Fill(jt, geoAccCorrRndm * fTrackPairEfficiency/jt);
        if(fill2DBackground) hBackground2D[fCentralityBin][iGap][fpttBin][iBin]->Fill(dEtaRndm, dPhiRndm, fTrackPairEfficiency);
        hPtAssoc[fCentralityBin][iGap][fpttBin][iBin]->Fill(fpta, fTrackPairEfficiency);
        
        // Fill the background histogram for signed pairs
        if(fIsLikeSign){
          hBackgroundLikeSign[fCentralityBin][iGap][fpttBin][iBin]->Fill(jt, geoAccCorrRndm * fTrackPairEfficiency/jt);
        } else {
          hBackgroundUnlikeSign[fCentralityBin][iGap][fpttBin][iBin]->Fill(jt, geoAccCorrRndm * fTrackPairEfficiency/jt);
        }
      }
    }
  }
}

/*
 *  Method for filling all the histograms needed in the jT analysis
 *
 *  fillType ftyp = Tell what kind of particles are filled (real events vs. mixed events)
 *  AliJBaseTrack *ftk1 = Track information for the trigger particle
 *  AliJBaseTrack *ftk2 = Track information for the associated particle
 *  bool fill2DBackground = Switch for quality control histograms. False = do not fill, true = fill
 *
 */
void AliJJtCorrelations::FillJtHistograms(fillType fTyp, AliJBaseTrack *ftk1, AliJBaseTrack *ftk2, bool fill2DBackground)
{
  // This method fills the jT and pout histograms
  
  double pout = fabs(fpta*sin(fPhiTrigger-fPhiAssoc));
  
  // Fill the distributions only if the near side definition is fullfilled
  
  TLorentzVector vTrigger = ftk1->GetLorentzVector(), vAssoc = ftk2->GetLorentzVector();   // Lorentz vectors for tracks
  
  /*
   *  Arguments for the distribution histogram filler
   *
   *  1)  fillType fTyp = Tell what kind of particles are filled (real events vs. mixed events)
   *  2)  int assocType = Associated particle binning type. 0 = klong, 1 = xlong, 2 = pta
   *  3)  TLorentzVector *vTrigger = vector for trigger particle
   *  4)  TLorentzVector *vAssoc = vector for associated particle
   *  5)  AliJTH1D &hDistribution = main jT histogram
   *  6)  AliJTH1D &hDistributionLikeSign = jT histogram for the like sign pairs
   *  7)  AliJTH1D &hDistributionUnlikeSign = jT histogram for the unlike sign pairs
   *  8)  AliJTH1D &hInvariantMass = invariant mass histogram
   *  9)  AliJTH1D &hInvariantMassLikeSign = invariant mass histogram for the like sign pairs
   *  10) AliJTH1D &hInvariantMassUnlikeSign = invariant mass histogram for the unlike sign pairs
   *
   */
  
  // Fill the jT histograms for klong binning. Require 3D near side
  if(fNearSide3D && fUseKlongBins) FillJtDistributionHistograms(fTyp,0,&vTrigger,&vAssoc,fhistos->fhJTKlong,fhistos->fhJTKlongLikeSign,fhistos->fhJTKlongUnlikeSign,fhistos->fhInvariantMassKlong,fhistos->fhInvariantMassKlongLikeSign,fhistos->fhInvariantMassKlongUnlikeSign);
  
  // Fill the jT histograms for xlong binning. Require 3D near side
  if(fNearSide3D) FillJtDistributionHistograms(fTyp,1,&vTrigger,&vAssoc,fhistos->fhJT,fhistos->fhJTLikeSign,fhistos->fhJTUnlikeSign,fhistos->fhInvariantMassXe,fhistos->fhInvariantMassXeLikeSign,fhistos->fhInvariantMassXeUnlikeSign);
  
  // Fill the jT histograms for pta binning. Require traditional near side
  if(fNearSide) FillJtDistributionHistograms(fTyp,2,&vTrigger,&vAssoc,fhistos->fhJTPta,fhistos->fhJTPtaLikeSign,fhistos->fhJTPtaUnlikeSign,fhistos->fhInvariantMassPta,fhistos->fhInvariantMassPtaLikeSign,fhistos->fhInvariantMassPtaUnlikeSign);
  
  //-------------------------------------
  // Rapidity/R/Phi Gap background
  // ------------------------------------
  
  // Vectors for randomized background particles
  TLorentzVector vAssocRndm;
  TLorentzVector vAssocRndmRgap;
  TLorentzVector vAssocRndmPhiGap;
  TLorentzVector vThrustRndm;
  TLorentzVector vThrustRndmRgap;
  TLorentzVector vThrustRndmPhiGap;
  
  // Variables for randomized eta and phi
  double etaTriggRndm, etaAssocRndm;
  double phiTriggRndm, phiAssocRndm;
  
  // Only do the randomization, if there is a background pair according to at least one of the background methods
  if(fTyp==kReal && ( fEtaGapBin >=0 || fRGapBinNear >=0 || fPhiGapBinNear >=0  ) ){

    for(int itrial=0; itrial<20; itrial++){  //For each high gap track generate twenty others
      
      // Randomize the vectors for the background pair
      
      if(fsamplingMethod == 0){ // Flat sampling
        
        // Randomize eta
        etaTriggRndm  = frandom->Uniform(-fmaxEtaRange, fmaxEtaRange);
        etaAssocRndm  = frandom->Uniform(-fmaxEtaRange, fmaxEtaRange);
        
        // Randomize phi
        phiTriggRndm = kJPi * frandom->Uniform(-1, 1);
        phiAssocRndm = kJPi * frandom->Uniform(-1, 1);
        
      } else { // Inclusive sampling
        
        // Randomize eta
        etaTriggRndm  = fhistos->fhIetaTriggFromFile[fCentralityBin][fpttBin]->GetRandom();
        etaAssocRndm  = fhistos->fhIetaAssocFromFile[fCentralityBin][fptaBin]->GetRandom();
        
        // Randomize phi
        phiTriggRndm = fhistos->fhIphiTriggFromFile[fCentralityBin][fpttBin]->GetRandom();
        phiAssocRndm = fhistos->fhIphiAssocFromFile[fCentralityBin][fptaBin]->GetRandom();
      }
      
      // Set up randomized trigger and associated vectors for different background scenarios
      vThrustRndm.SetPtEtaPhiM(vTrigger.Pt(),etaTriggRndm, fPhiTrigger, 0);  // Trigger with only eta randomized
      vAssocRndm.SetPtEtaPhiM(fpta, etaAssocRndm, fPhiAssoc, 0);             // Associated with only eta randomized
      
      vThrustRndmRgap.SetPtEtaPhiM(vTrigger.Pt(),etaTriggRndm, phiTriggRndm, 0);  // Trigger with eta and phi randomized
      vAssocRndmRgap.SetPtEtaPhiM(fpta, etaAssocRndm, phiAssocRndm, 0);           // Associated with eta and phi randomized
      
      vThrustRndmPhiGap.SetPtEtaPhiM(vTrigger.Pt(),fEtaTrigger, phiTriggRndm, 0);  // Trigger with only phi randomized
      vAssocRndmPhiGap.SetPtEtaPhiM(fpta, fEtaAssoc, phiAssocRndm, 0);             // Associated with only phi randomized
      
      /*
       * Arguments for the background histogram filler
       *
       *  1)  int assocType = Associated particle binning type. 0 = klong, 1 = xlong, 2 = pta
       *  2)  int gapType = Used gap type for background. 0 = eta gap, 1 = R gap, 2 = phi gap
       *  3)  TLorentzVector *vTrigger = trigger particle used for the background estimation
       *  4)  TLorentzVector *vAssoc = associated particle used for the background estimation
       *  5)  AliJTH1D &hBackground = main background histogram
       *  6)  AliJTH1D &hBackgroundLikeSign = background histogram for the like sign pairs
       *  7)  AliJTH1D &hBackgroundUnlikeSign = background histogram for the unlike sign pairs
       *  8)  AliJTH1D &hPtAssoc = associated particle pT distribution
       *  9)  AliJTH2D &hBackground2D = two dimensional quality control histogram
       *  10) bool fill2DBackground = Switch for quality control histograms. False = do not fill, true = fill
       *
       */
      
      /*
       * Note: We must accept all the pairs for 3D near side calculation.
       * The reason for this is that we do not want to bias the DeltaEta and DeltaPhi distributions.
       * For example requiring large eta gap biases DeltaPhi distribution towards small values.
       * Thus more background is seen in the small jT region.
       * In reality, we could also have large DeltaPhi pairs with small DeltaEta in the background.
       * For this reason, no 3D near side requirement is required for the pairs entering the background calculation.
       * FillJtBackgroundHistograms method ensures that the generated pair is only filled if it happens
       * to be in the 3D near side.
       */
      
      // Only fill klong bins is enabled from JCard
      if(fUseKlongBins){
        
        // Fill background histograms for klong binning, eta gap
        FillJtBackgroundHistograms(0,0,&vThrustRndm,&vAssocRndm,fhistos->fhJTKlongBg,fhistos->fhJTKlongBgLikeSign,fhistos->fhJTKlongBgUnlikeSign,fhistos->fhBgAssocKlongEta,fhistos->fhDphiDetaBgKlongEta,fill2DBackground);
        
        // Fill background histograms for klong binning, R gap
        FillJtBackgroundHistograms(0,1,&vThrustRndmRgap,&vAssocRndmRgap,fhistos->fhJTKlongBgR,fhistos->fhJTKlongBgRLikeSign,fhistos->fhJTKlongBgRUnlikeSign,fhistos->fhBgAssocKlongR,fhistos->fhDphiDetaBgKlongR,fill2DBackground);
        
        // Fill background histograms for klong binning, phi gap
        FillJtBackgroundHistograms(0,2,&vThrustRndmPhiGap,&vAssocRndmPhiGap,fhistos->fhJTKlongBgPhi,fhistos->fhJTKlongBgPhiLikeSign,fhistos->fhJTKlongBgPhiUnlikeSign,fhistos->fhBgAssocKlongPhi,fhistos->fhDphiDetaBgKlongPhi,fill2DBackground);
        
      }
      
      // Fill background histograms for xlong binning, eta gap
      FillJtBackgroundHistograms(1,0,&vThrustRndm,&vAssocRndm,fhistos->fhJTBg,fhistos->fhJTBgLikeSign,fhistos->fhJTBgUnlikeSign,fhistos->fhBgAssocXlongEta,fhistos->fhDphiDetaBgXlongEta,fill2DBackground);
      
      // Fill background histograms for xlong binning, R gap
      FillJtBackgroundHistograms(1,1,&vThrustRndmRgap,&vAssocRndmRgap,fhistos->fhJTBgR,fhistos->fhJTBgRLikeSign,fhistos->fhJTBgRUnlikeSign,fhistos->fhBgAssocXlongR,fhistos->fhDphiDetaBgXlongR,fill2DBackground);
      
      // Fill background histograms for xlong binning, phi gap
      FillJtBackgroundHistograms(1,2,&vThrustRndmPhiGap,&vAssocRndmPhiGap,fhistos->fhJTBgPhi,fhistos->fhJTBgPhiLikeSign,fhistos->fhJTBgPhiUnlikeSign,fhistos->fhBgAssocXlongPhi,fhistos->fhDphiDetaBgXlongPhi,fill2DBackground);
      
      
      // Only fill the pTa background histograms, if the original pair was in traditional near side. Note that the background histogram filler method also checks, that the newly generated pair is on the traditional near side.
      if(fNearSide){
        
        // Fill background histograms for pta binning, eta gap
        FillJtBackgroundHistograms(2,0,&vThrustRndm,&vAssocRndm,fhistos->fhJTPtaBg,fhistos->fhJTPtaBgLikeSign,fhistos->fhJTPtaBgUnlikeSign,fhistos->fhBgAssocPtaEta,fhistos->fhDphiDetaBgPtaEta,fill2DBackground);
        
        // Fill background histograms for pta binning, R gap
        FillJtBackgroundHistograms(2,1,&vThrustRndmRgap,&vAssocRndmRgap,fhistos->fhJTPtaBgR,fhistos->fhJTPtaBgRLikeSign,fhistos->fhJTPtaBgRUnlikeSign,fhistos->fhBgAssocPtaR,fhistos->fhDphiDetaBgPtaR,fill2DBackground);
        
        // Fill background histograms for pta binning, phi gap
        FillJtBackgroundHistograms(2,2,&vThrustRndmPhiGap,&vAssocRndmPhiGap,fhistos->fhJTPtaBgPhi,fhistos->fhJTPtaBgPhiLikeSign,fhistos->fhJTPtaBgPhiUnlikeSign,fhistos->fhBgAssocPtaPhi,fhistos->fhDphiDetaBgPtaPhi,fill2DBackground);
        
      } // Traditional near side check
      
    } // trials
  } // background randomization
}

/*
 * Fill the deltaEta histogram for simplistic acceptance correction
 *
 *  fillType fTyp = real or mixed events
 */
void AliJJtCorrelations::FillDeltaEtaHistograms(fillType fTyp)
{
  // This method fills the mixed event DeltaEta histograms
  
  if( fNearSide ){ //one could check the phiGapBin, but in the pi/2 <1.6 and thus phiGap is always>-1
    if( fTyp == 1 ) {
      fhistos->fhDetaNearMixAcceptance[fCentralityBin][fpttBin][fptaBin]->Fill( fDeltaEta, fTrackPairEfficiency);
    }
  }

}

/*
 * Fill deltaEta deltaPhi histograms to be used for acceptance correction
 *
 *  fillType fTyp = real or mixed events
 */
void AliJJtCorrelations::FillDeltaEtaDeltaPhiHistograms(fillType fTyp)
{
  // This method fills the two dimensional DeltaEta,DeltaPhi histograms
  // No acceptance correction here, since we want to see tha structure caused by acceptance effects
  
  // Fill the histogram in pTa bins
  if(fNearSide){
    fhistos->fhDphiDetaPta[fTyp][fCentralityBin][fZBin][fpttBin][fptaBin]->Fill(fDeltaEta, fDeltaPhiPiPi, fTrackPairEfficiency);
  }
  
  // Fill the histogram in xlong bins
  if(fNearSide3D && fXlongBin >= 0){
    fhistos->fhDphiDetaXlong[fTyp][fCentralityBin][fZBin][fpttBin][fXlongBin]->Fill(fDeltaEta, fDeltaPhiPiPi, fTrackPairEfficiency);
  }
  
}

/*
 * Fill associated particle pT distribution in pTt and pTa bins
 *
 *  fillType fTyp = real or mixed events
 */
void AliJJtCorrelations::FillPtaHistograms(fillType fTyp)
{
  // This method fills various pta histograms
  
  if ( fTyp == kReal ) {
    // pTa distribution in pTt and pTa bins
    fhistos->fhAssocPtBin[fCentralityBin][fpttBin][fptaBin]->Fill(fpta); // Not weighted by efficiency?
    
    fnReal++;
  } else { // only mix
    fnMix++;
  }
}





















