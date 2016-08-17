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

#include "AliJIaaCorrelations.h"

#include "../AliJBaseTrack.h"
#include "../AliJCard.h"

double DeltaPhi(double phi1, double phi2) {
	// dphi
	double res =  atan2(sin(phi1-phi2), cos(phi1-phi2));
	return res > -kJPi*9./20. ? res : kJTwoPi+res ;
}

AliJIaaCorrelations::AliJIaaCorrelations( AliJCard *cardIn, AliJIaaHistograms *histosIn) :
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
  fIsLikeSign(false),
  fUseZVertexBinsAcceptance(false)
{
  // constructor
  //fcard->PrintOut();
  fmaxEtaRange    = fcard->Get("EtaRange");
  
  frandom = new TRandom3(); // Random number generator for background randomization
  frandom->SetSeed(0);
  
}

AliJIaaCorrelations::AliJIaaCorrelations() :
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
  fIsLikeSign(false),
  fUseZVertexBinsAcceptance(false)
{
  // default constructor
}

AliJIaaCorrelations::AliJIaaCorrelations(const AliJIaaCorrelations& in) :
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
  fIsLikeSign(in.fIsLikeSign),
  fUseZVertexBinsAcceptance(in.fUseZVertexBinsAcceptance)
{
  // The pointers to card and histos are just copied. I think this is safe, since they are not created by
  // AliJIaaCorrelations and thus should not disappear if the AliJCorrelation managing them is destroyed.
  
  frandom = new TRandom3(); // // Random number generator for background randomization
  frandom->SetSeed(0);
}

AliJIaaCorrelations& AliJIaaCorrelations::operator=(const AliJIaaCorrelations& in){
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
  fIsLikeSign = in.fIsLikeSign;
  fnReal = in.fnReal;
  fnMix = in.fnMix;
  fsamplingMethod = in.fsamplingMethod;
  fmaxEtaRange = in.fmaxEtaRange;
  fUseZVertexBinsAcceptance = in.fUseZVertexBinsAcceptance;
  
  // The pointers to card and histos are just copied. I think this is safe, since they are not created by
  // AliJIaaCorrelations and thus should not disappear if the AliJCorrelation managing them is destroyed.
  fcard = in.fcard;
  fhistos = in.fhistos;
  fAcceptanceCorrection = in.fAcceptanceCorrection;

  frandom = new TRandom3(); // Random number generator for background randomization
  frandom->SetSeed(0);
  
  return *this;
  // copy constructor
}

AliJIaaCorrelations::~AliJIaaCorrelations(){
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
void AliJIaaCorrelations::FillHisto(corrFillType cFTyp, fillType fTyp, int cBin, int zBin, AliJBaseTrack *ftk1, AliJBaseTrack *ftk2){
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
void AliJIaaCorrelations::FillCorrelationHistograms(fillType fTyp, int CentBin, int zBin, AliJBaseTrack *ftk1, AliJBaseTrack *ftk2)
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
  fDeltaPhi     = DeltaPhi(fPhiTrigger, fPhiAssoc);
  //fDeltaPhi     = ftk1->DeltaPhi(*ftk2);
  fDeltaPhiPiPi = atan2(sin(fPhiTrigger-fPhiAssoc), cos(fPhiTrigger-fPhiAssoc));
  fDeltaEta     = ftk1->Eta() - ftk2->Eta();
  fEtaTrigger   = ftk1->Eta();
  fEtaAssoc     = ftk2->Eta();
  
  fNearSide     = cos(fPhiTrigger-fPhiAssoc) > 0 ? true : false;  // Traditional near side definition using deltaPhi
  fNearSide3D   = vTrigger.Vect().Dot(vAssoc.Vect()) > 0 ? true : false; // Near side definition using half ball around the trigger

  fEtaGapBin     = fcard->GetBinFast( kEtaGapType, fabs(fDeltaEta));
  fPhiGapBinNear = fcard->GetBin( kEtaGapType, fabs(fDeltaPhiPiPi) );
  fRGapBinNear   = fcard->GetBin( kRGapType, fabs(ftk1->DeltaR(*ftk2)));
  fCentralityBin = CentBin;
  fZBin = zBin;
  
  fXlong = vTrigger.Vect().Dot(vAssoc.Vect())/pow(vTrigger.P(),2);
  fXlongBin = fcard->GetBin(kXeType, TMath::Abs(fXlong));
  
  // Check that the bins make sense. If not, do not fill anything.
  if(fpttBin<0 || fptaBin<0 || fEtaGapBin<0 ){
	cout<<"Error in FillCorrelationHistograms: some pT or eta out of bin. pttBin="<<fpttBin<<" pTaBin="<<fptaBin <<" etaGapBin="<< fEtaGapBin <<" deltaEta=" << fDeltaEta << endl;
	//ftk1->Print();
	//ftk2->Print();
    exit(-1);
    //return;
  }
  
  
  // ===================================================================
  // =====================  Fill Histograms  ===========================
  // ===================================================================
  
  bool fill2DBackgroundQualityControlHistograms = false;  // Choose whether to fill the DeltaPhi DeltaEta histograms for jT background
  
  // If quality control level is high enough, fill 2D quality control histograms
  if(fcard->Get("QualityControlLevel")>1) fill2DBackgroundQualityControlHistograms = true;
  
  FillDeltaEtaHistograms(fTyp);  // Fill all the delta eta histograms
  FillDeltaEtaDeltaPhiHistograms(fTyp);  // Fill the 2D correlation functions
  FillDeltaPhiHistograms(fTyp); // Fill all the delta phi histograms
  
}



/*
 *  Find out the bin index based on the associated binning type and trigger and asociated vectors
 *
 *  int assocType = Associated particle binning type. 0 = klong, 1 = xlong, 2 = pta
 *  TLorentzVector *vTrigger = vector for trigger particle
 *  TLorentzVector *vAssoc = vector for associated particle
 *
 */
int AliJIaaCorrelations::GetBinIndex(int assocType, TLorentzVector *vTrigger, TLorentzVector *vAssoc)
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
 * Fill the deltaEta histogram for simplistic acceptance correction
 *
 *  fillType fTyp = real or mixed events
 */
void AliJIaaCorrelations::FillDeltaEtaHistograms(fillType fTyp)
{
  // This method fills the DeltaEta histograms

  if( fNearSide ){ //one could check the phiGapBin, but in the pi/2 <1.6 and thus phiGap is always>-1
	if( fTyp == 0 ) {
	  fhistos->fhDEtaNear[fCentralityBin][fPhiGapBinNear][fpttBin][fptaBin]->Fill( fDeltaEta ,  fTrackPairEfficiency ); //fGeometricAcceptanceCorrection?
	} else {
	  fhistos->fhDEtaNearM[fCentralityBin][fPhiGapBinNear][fpttBin][fptaBin]->Fill( fDeltaEta ,  fTrackPairEfficiency );//fGeometricAcceptanceCorrection?
	  fhistos->fhDetaNearMixAcceptance[fCentralityBin][fpttBin][fptaBin]->Fill( fDeltaEta, fTrackPairEfficiency);
	}
  }
}

/*
 * Fill deltaEta deltaPhi histograms to be used for acceptance correction
 *
 *  fillType fTyp = real or mixed events
 */
void AliJIaaCorrelations::FillDeltaEtaDeltaPhiHistograms(fillType fTyp)
{
  // This method fills the two dimensional DeltaEta,DeltaPhi histograms
  // No acceptance correction here, since we want to see the structure caused by acceptance effects
  

	fhistos->fhDphiDetaPta[fTyp][fCentralityBin][fpttBin][fptaBin]->Fill(fDeltaEta, fDeltaPhi/kJPi, fTrackPairEfficiency);
//		fhistos->fhDphiDetaPtaRgap[fTyp][fRGapBinNear][fCentralityBin][fpttBin][fptaBin]->Fill(fDeltaEta, fDeltaPhi/kJPi, fTrackPairEfficiency);
}

void AliJIaaCorrelations::FillDeltaPhiHistograms(fillType fTyp)
{
  // This method fills the DeltaPhi histograms

  // When hists are filled for thresholds they are not properly normalized and need to be subtracted
  // This induced improper errors - subtraction of not-independent entries

  fhistos->fhDphiAssoc[fTyp][fCentralityBin][fEtaGapBin][fpttBin][fptaBin]->Fill( fDeltaPhi/kJPi , fTrackPairEfficiency); //fGeometricAcceptanceCorrection?
}






















