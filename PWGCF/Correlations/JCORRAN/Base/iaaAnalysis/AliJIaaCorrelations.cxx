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
	return res > -kJPi*0.5 ? res : kJTwoPi+res ;
}

AliJIaaCorrelations::AliJIaaCorrelations( AliJCard *cardIn, AliJIaaHistograms *histosIn) :
	AliJCorrelationInterface(),
	fcard(cardIn),
	fhistos(histosIn),
//	fAcceptanceCorrection(0x0),
	fnReal(0),
	fnMix(0),
	fsamplingMethod(0),    // flat by default
	fmaxEtaRange(0),
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
	fMagneticFieldPolarity(0),
	fTrackMergeCut(0),
    fNearSide(true),
    fNearSide3D(true),
    fEtaGapBin(0),
	fPhiGapBinNear(0),
	fRGapBinNear(0),
	fCentralityBin(0),
	fXlongBin(0),
	fUseZVertexBinsAcceptance(false),
    fRequireLikeSign(0),
    fIsLikeSign(false),
    fUseTrackMergingCorr(0)
{
	// constructor
    fmaxEtaRange = fcard->Get("EtaRange");
    fRequireLikeSign = fcard->Get("RequireLikeSign"); // 0: all, 1: like-sign, -1: opposite-sign
    fUseTrackMergingCorr = fcard->Get("UseTrackMergingCorr"); // correct for track merging (1) or not (0)
    fTrackMergeCut = fcard->Get("TrackMergeCut"); // minimum of acceptable track separation
}

AliJIaaCorrelations::AliJIaaCorrelations() :
	AliJCorrelationInterface(),
	fcard(0x0),
	fhistos(0x0),
	//fAcceptanceCorrection(0x0),
	fnReal(0),
	fnMix(0),
	fsamplingMethod(0),    // flat by default
	fmaxEtaRange(0),
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
	fMagneticFieldPolarity(0),
	fTrackMergeCut(0),
	fNearSide(true),
	fNearSide3D(true),
    fEtaGapBin(0),
	fPhiGapBinNear(0),
	fRGapBinNear(0),
	fCentralityBin(0),
	fXlongBin(0),
	fUseZVertexBinsAcceptance(false),
    fRequireLikeSign(0),
    fIsLikeSign(false),
    fUseTrackMergingCorr(0)
{
	// default constructor
}

AliJIaaCorrelations::AliJIaaCorrelations(const AliJIaaCorrelations& in) :
	fcard(in.fcard),
	fhistos(in.fhistos),
	//fAcceptanceCorrection(in.fAcceptanceCorrection),
	fnReal(in.fnReal),
	fnMix(in.fnMix),
	fsamplingMethod(in.fsamplingMethod),
	fmaxEtaRange(in.fmaxEtaRange),
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
	fMagneticFieldPolarity(in.fMagneticFieldPolarity),
	fTrackMergeCut(in.fTrackMergeCut),
	fNearSide(in.fNearSide),
	fNearSide3D(in.fNearSide3D),
    fEtaGapBin(in.fEtaGapBin),
	fPhiGapBinNear(in.fPhiGapBinNear),
	fRGapBinNear(in.fRGapBinNear),
	fCentralityBin(in.fCentralityBin),
	fXlongBin(in.fXlongBin),
	fUseZVertexBinsAcceptance(in.fUseZVertexBinsAcceptance),
    fRequireLikeSign(in.fRequireLikeSign),
    fIsLikeSign(in.fIsLikeSign),
    fUseTrackMergingCorr(in.fUseTrackMergingCorr)
{
	// The pointers to card and histos are just copied. I think this is safe, since they are not created by
	// AliJIaaCorrelations and thus should not disappear if the AliJCorrelation managing them is destroyed.
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
	fMagneticFieldPolarity = in.fMagneticFieldPolarity;
	fTrackMergeCut = in.fTrackMergeCut;
	fNearSide = in.fNearSide;
	fNearSide3D = in.fNearSide3D;
    fEtaGapBin = in.fEtaGapBin;
	fPhiGapBinNear = in.fPhiGapBinNear;
	fRGapBinNear = in.fRGapBinNear;
	fCentralityBin = in.fCentralityBin;
	fXlongBin = in.fXlongBin;
	fnReal = in.fnReal;
	fnMix = in.fnMix;
	fsamplingMethod = in.fsamplingMethod;
	fmaxEtaRange = in.fmaxEtaRange;
	fUseZVertexBinsAcceptance = in.fUseZVertexBinsAcceptance;
    fRequireLikeSign = in.fRequireLikeSign;
    fIsLikeSign = in.fIsLikeSign;
    fUseTrackMergingCorr = in.fUseTrackMergingCorr;

	// The pointers to card and histos are just copied. I think this is safe, since they are not created by
	// AliJIaaCorrelations and thus should not disappear if the AliJCorrelation managing them is destroyed.
	fcard = in.fcard;
	fhistos = in.fhistos;
	//fAcceptanceCorrection = in.fAcceptanceCorrection;

	return *this;
	// copy constructor
}

AliJIaaCorrelations::~AliJIaaCorrelations()
{
	// destructor
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
	if(ftk1->GetParticleType()==kJHadronMC && ftk2->GetParticleType()==kJHadronMC) twoTracks =true;

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
	fDeltaPhiPiPi = atan2(sin(fPhiTrigger-fPhiAssoc), cos(fPhiTrigger-fPhiAssoc));
	fDeltaEta     = ftk1->Eta() - ftk2->Eta();
	fEtaTrigger   = ftk1->Eta();
	fEtaAssoc     = ftk2->Eta();

	fNearSide     = cos(fPhiTrigger-fPhiAssoc) > 0 ? true : false;  // Traditional near side definition using deltaPhi
	fNearSide3D   = vTrigger.Vect().Dot(vAssoc.Vect()) > 0 ? true : false; // Near side definition using half ball around the trigger

	fEtaGapBin     = fcard->GetBinFast( kEtaGapType, fabs(fDeltaEta));
  fPhiGapBinNear = fcard->GetBin( kPhiGapType, fabs(fDeltaPhiPiPi) );
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

	// Track merging correction
    if(fUseTrackMergingCorr)
    {
        if (TMath::Abs(fDeltaEta) < fTrackMergeCut * 2.5 * 3)
        {
            double triggerCharge = ftk1->GetCharge();
            double associatedCharge = ftk2->GetCharge();
            double twoTrackCutMinRadius = 0.039; // Value from Monika Kofarago

            // check first boundaries to see if is worth to loop and find the minimum
            Float_t dphistar1 = GetDPhiStar(fPhiTrigger, fptt, triggerCharge, fPhiAssoc, fpta, associatedCharge, twoTrackCutMinRadius, fMagneticFieldPolarity);
            Float_t dphistar2 = GetDPhiStar(fPhiTrigger, fptt, triggerCharge, fPhiAssoc, fpta, associatedCharge, 2.5, fMagneticFieldPolarity);

            const Float_t kLimit = fTrackMergeCut * 3;

            Float_t dphistarminabs = 1e5;
            Float_t dphistarmin = 1e5;
            if (TMath::Abs(dphistar1) < kLimit || TMath::Abs(dphistar2) < kLimit || dphistar1 * dphistar2 < 0){
                for (Double_t radius=twoTrackCutMinRadius; radius<2.51; radius+=0.01){
                    Float_t dphistar = GetDPhiStar(fPhiTrigger, fptt, triggerCharge, fPhiAssoc, fpta, associatedCharge, radius, fMagneticFieldPolarity);

                    Float_t dphistarabs = TMath::Abs(dphistar);

                    if (dphistarabs < dphistarminabs){
                        dphistarmin = dphistar;
                        dphistarminabs = dphistarabs;
                    }
                }

                // Remove the pair if deltaPhi* and deltaEta values are too close
                if (dphistarminabs < fTrackMergeCut && TMath::Abs(fDeltaEta) < fTrackMergeCut){
                    return;
                }
            }
        }
    }

    // check pair sign only if +1:same or -1:opposite check is required
    if(fRequireLikeSign!=0) {
        // break if requirement is not met
        if(fIsLikeSign==true && fRequireLikeSign!=+1) return;
        if(fIsLikeSign==false && fRequireLikeSign!=-1 ) return;
    }
	// ===================================================================
	// =====================  Fill Histograms  ===========================
	// ===================================================================

	if(ResonanceCut(ftk1,ftk2))
	{
        //FillDeltaEtaHistograms(fTyp);  // Fill all the delta eta histograms
		FillDeltaEtaDeltaPhiHistograms(fTyp);  // Fill the 2D correlation functions
	} else {
        FillResonanceHistograms(fTyp);
	}
}

/*
 * Fill the deltaEta histogram for simplistic acceptance correction
 *
 *  fillType fTyp = real or mixed events
 */
void AliJIaaCorrelations::FillDeltaEtaHistograms(fillType fTyp)
{
	// This method fills the DeltaEta histograms

	if( fNearSide ){
		if( fTyp == 0 ) {
            //fhistos->fhDEtaNear[fCentralityBin][fZBin][fPhiGapBinNear][fpttBin][fptaBin]->Fill( fDeltaEta ,  fTrackPairEfficiency );
		} else {
            //fhistos->fhDEtaNearM[fCentralityBin][fZBin][fPhiGapBinNear][fpttBin][fptaBin]->Fill( fDeltaEta ,  fTrackPairEfficiency );
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

	fhistos->fhDphiDetaPta[fTyp][fCentralityBin][fZBin][fpttBin][fptaBin]->Fill(fDeltaEta, fDeltaPhi, fTrackPairEfficiency);
}

void AliJIaaCorrelations::FillResonanceHistograms(fillType fTyp)
{
	//fhistos->fhResonanceCut[fTyp][fCentralityBin][fZBin][fpttBin][fptaBin]->Fill(fDeltaEta, fDeltaPhi, fTrackPairEfficiency);
}

bool AliJIaaCorrelations::ResonanceCut(AliJBaseTrack *ftk1, AliJBaseTrack *ftk2)
{
	//Conversions
	if( ftk1->GetCharge() * ftk2->GetCharge() < 0)
	{
		Float_t mass = GetInvMassSquaredCheap(fptt, ftk1->Eta(), fPhiTrigger, fpta, ftk2->Eta(), fPhiAssoc, 0.510e-3, 0.510e-3);
		if (mass < 0.1){
			mass = GetInvMassSquared(fptt, ftk1->Eta(), fPhiTrigger, fpta, ftk2->Eta(), fPhiAssoc, 0.510e-3, 0.510e-3);
			if(mass < 0.04 * 0.04)
				return kFALSE;
		}
	}
	//K0s
	if ( ftk1->GetCharge() * ftk2->GetCharge() < 0){
		Float_t mass = GetInvMassSquaredCheap(fptt, ftk1->Eta(), fPhiTrigger, fpta, ftk2->Eta(), fPhiAssoc, 0.1396, 0.1396);
		const Float_t kK0smass = 0.4976;
		if (TMath::Abs(mass -kK0smass * kK0smass)  < 0.1){
			mass = GetInvMassSquared(fptt, ftk1->Eta(), fPhiTrigger, fpta, ftk2->Eta(), fPhiAssoc, 0.1396, 0.1396);
			if(mass > (kK0smass-0.02)*(kK0smass-0.02) && mass < (kK0smass+0.02)*(kK0smass+0.02))
				return kFALSE;
		}
	}
	//Lambda
	if (ftk1->GetCharge() * ftk2->GetCharge() < 0){
		Float_t mass1 = GetInvMassSquaredCheap(fptt, ftk1->Eta(), fPhiTrigger, fpta, ftk2->Eta(), fPhiAssoc, 0.1396, 0.9383);
		Float_t mass2 = GetInvMassSquaredCheap(fptt, ftk1->Eta(), fPhiTrigger, fpta, ftk2->Eta(), fPhiAssoc, 0.9383, 0.1396);
		const Float_t kLambdaMass = 1.115;
		if (TMath::Abs(mass1 -kLambdaMass * kLambdaMass)  < 0.1){
			mass1 = GetInvMassSquared(fptt, ftk1->Eta(), fPhiTrigger, fpta, ftk2->Eta(), fPhiAssoc, 0.1396, 0.9383);
			if(mass1 > (kLambdaMass-0.02)*(kLambdaMass-0.02) && mass1 < (kLambdaMass+0.02)*(kLambdaMass+0.02))
				return kFALSE;
		}
		if (TMath::Abs(mass2 -kLambdaMass * kLambdaMass)  < 0.1){
			mass2 = GetInvMassSquared(fptt, ftk1->Eta(), fPhiTrigger, fpta, ftk2->Eta(), fPhiAssoc, 0.1396, 0.9383);
			if(mass2 > (kLambdaMass-0.02)*(kLambdaMass-0.02) && mass2 < (kLambdaMass+0.02)*(kLambdaMass+0.02))
				return kFALSE;
		}
	}

	return kTRUE;
}




/*
 * Calculate deltaPhi*
 *
 *  Float_t phi1 = azimuthal angle of the trigger particle
 *  Float_t pt1 = transverse momentum of the trigger particle
 *  Float_t charge1 = charge of the trigger particle
 *  Float_t phi2 = azimuthal angle of the associated particle
 *  Float_t pt2 = transverse momentum of the associated particle
 *  Float_t charge2 = charge of the associated particle
 *  Float_t radius = radius inside TPC
 *  Float_t bSign = magnetic field polarisation
 */
Float_t AliJIaaCorrelations::GetDPhiStar(Float_t phi1, Float_t pt1, Float_t charge1, Float_t phi2, Float_t pt2, Float_t charge2, Float_t radius, Float_t bSign){
	//
	// calculates dphistar
	//

	Float_t dphistar = phi1 - phi2 - charge1 * bSign * TMath::ASin(0.075 * radius / pt1) + charge2 * bSign * TMath::ASin(0.075 * radius / pt2);

	if (dphistar > kJPi)
		dphistar = kJTwoPi - dphistar;
	if (dphistar < -1.*kJPi)
		dphistar = -kJTwoPi - dphistar;
	if (dphistar > kJPi) // might look funny but is needed
		dphistar = kJTwoPi - dphistar;

	return dphistar;
}















