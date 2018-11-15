/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

// Implementation for correlation analysis

//===========================================================
// AliJIaaCorrelations.h
//
// Author: Marton Vargyas, Jussi Viinikainen
//===========================================================

#ifndef AliJIaaCORRELATIONS_H
#define AliJIaaCORRELATIONS_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <TRandom3.h>

#include "AliJIaaHistograms.h"
#include "../AliJAcceptanceCorrection.h"
#include "../AliJCorrelationInterface.h"

using namespace std;

class AliJIaaCorrelations;
class AliJIaaHistograms;
class AliJBaseTrack;
class AliJCard;

class AliJIaaCorrelations : public AliJCorrelationInterface{

public:

    AliJIaaCorrelations(); // default constructor
    AliJIaaCorrelations( AliJCard *cardIn, AliJIaaHistograms *histosIn); // constructor
    AliJIaaCorrelations(const AliJIaaCorrelations& in); // copy constructor
    virtual ~AliJIaaCorrelations(); //destructor
    AliJIaaCorrelations& operator=(const AliJIaaCorrelations& obj); // equal sign operator

    void PrintOut(){cout<<"Number of events = "<<fnReal<<"  Mixed events = "<<fnMix<<endl;} // Event count print

    void FillHisto(corrFillType cFTyp, fillType fTyp, int cBin, int zBin, AliJBaseTrack *ftk1, AliJBaseTrack *ftk2); // correlation histogram filler based on correlation type
    void FillCorrelationHistograms (fillType fTyp, int cBin, int zBin, AliJBaseTrack *ftk1, AliJBaseTrack *ftk2); // correlation histogram filler

    void SetSamplingInclusive(){fsamplingMethod = 1;} // Setter for inclusive sampling
    void UseZVertexAcceptance(bool useZ){fUseZVertexBinsAcceptance = useZ;} // Setter for fUseZVertexBinsAcceptance
    //void SetAcceptanceCorrection(AliJAcceptanceCorrection *accCorr){fAcceptanceCorrection = accCorr;} // Setter for acceptance correction
    void SetMagneticFieldPolarity(double magneticFieldPolarity){fMagneticFieldPolarity = magneticFieldPolarity;}

protected:

    AliJCard* fcard; // Card with binning information etc.
    AliJIaaHistograms* fhistos;  // Histograms needed in the analysis
    //AliJAcceptanceCorrection *fAcceptanceCorrection;  // Acceptance correction container
    int fnReal; // Number of events
    int fnMix; // Number of mixed events
    int fsamplingMethod; // Sampling method flag (flat or inclusive)
    double fmaxEtaRange; // Eta range used in the analysis

    double fptt;  // pT of the trigger particle in the correlation loop
    double fpta;  // pT ot the asociated particle in the correlation loop
    double fTrackPairEfficiency;  // pair efficiency for the tracks in the correlation loop
    int fpttBin;  // Bin index for the trigger pT bin
    int fptaBin;  // Bin index for the associated pT bin
    double fPhiTrigger;  // Azimuthal angle of the trigger particle
    double fPhiAssoc;  // Asimuthal angle of the associated particle
    double fEtaTrigger;  // Pseodurapidity of the trigger particle
    double fEtaAssoc;  // Pseudorapidity of the associated particle
    double fDeltaPhi;  // Difference of the azimuthal angles of trigger and associated particles
    double fDeltaPhiPiPi;  // The same as above but measured from -pi to pi
    double fDeltaEta;  // Difference of the pseudorapidities of the trigger and associated particles
    double fXlong;  // The xlong value of the trigger and associated particles
    int fZBin; // z-vertex bin
    double fMagneticFieldPolarity;  // polarity of the magnetic field
    double fTrackMergeCut; // dphi/deta smaller than this will be cut out (disable with 0)

    bool fNearSide;  // true if near side correlation, false if away side
    bool fNearSide3D; // near side defined by the half ball around the trigger
    int fEtaGapBin;  // Bin index for the current eta gap
    int fPhiGapBinNear;  // Bin index for the phi gap in the near side
    int fRGapBinNear;  // Bin index for the R gap in the near side
    int fCentralityBin;  // Bin index for the centrality bin
    int fXlongBin;  // Bin index for xlong bin

    bool fUseZVertexBinsAcceptance; // false = integrate over z-vertex bins, true = Do acceptance correction in z-vertex bins
    int  fRequireLikeSign; // Requirement from card (+1: only like-sign, -1: only opposite sign, 0: all)
    bool fIsLikeSign; // Track property (true: like-sign, false: opposite-sign)
    bool fUseTrackMergingCorr; // true: correct, false: don't correct


private:
    bool ResonanceCut(AliJBaseTrack *ftk1, AliJBaseTrack *ftk2); // cut on long-lived resonances (K0, Lambda), and on conversion photons. True if no cut.
    bool ResonanceCutCheap(AliJBaseTrack *ftk1, AliJBaseTrack *ftk2); // Approximate but quicker solution for the upper function.
    void FillDeltaEtaHistograms(fillType fTyp); // deltaEta histogram filler
    void FillDeltaEtaDeltaPhiHistograms(fillType fTyp); // deltaEta deltaPhi histogram filler
    void FillResonanceHistograms(fillType fTyp); // histogram filler for cut-out resonances
    Float_t GetDPhiStar(Float_t phi1, Float_t pt1, Float_t charge1, Float_t phi2, Float_t pt2, Float_t charge2, Float_t radius, Float_t bSign);

    Float_t GetInvMassSquared(Float_t pt1, Float_t eta1, Float_t phi1, Float_t pt2, Float_t eta2, Float_t phi2, Float_t m0_1, Float_t m0_2)
    {
        // calculate inv mass squared
        Float_t tantheta1 = 1e10;
        if (eta1 < -1e-10 || eta1 > 1e-10){
            Float_t expTmp = TMath::Exp(-eta1);
            tantheta1 = 2.0 * expTmp/(1.0 -expTmp*expTmp);
        }

        Float_t tantheta2 = 1e10;
        if (eta2 < -1e-10 || eta2 > 1e-10){
            Float_t expTmp = TMath::Exp(-eta2);
            tantheta2 = 2* expTmp/(1.0 - expTmp*expTmp);
        }

        Float_t e1squ = m0_1 * m0_1 + pt1 * pt1 * (1.0 + 1.0 / tantheta1 / tantheta1);
        Float_t e2squ = m0_2 * m0_2 + pt2 * pt2 * (1.0 + 1.0 / tantheta2 / tantheta2);
        Float_t mass2 = m0_1 * m0_1 + m0_2 * m0_2 + 2 * (TMath::Sqrt(e1squ * e2squ) - ( pt1 * pt2 * ( TMath::Cos(phi1 - phi2) + 1.0 / tantheta1 / tantheta2) ) );
        return mass2;
    }
    Float_t GetInvMassSquaredCheap(Float_t pt1, Float_t eta1, Float_t phi1, Float_t pt2, Float_t eta2, Float_t phi2, Float_t m0_1, Float_t m0_2)
    {
        // calculate inv mass squared approximately
        Float_t tantheta1 = 1e10;
        if (eta1 < -1e-10 || eta1 > 1e-10){
            Float_t expTmp = 1.0-eta1+eta1*eta1/2-eta1*eta1*eta1/6+eta1*eta1*eta1*eta1/24;
            tantheta1 = 2.0 * expTmp / ( 1.0 - expTmp*expTmp);
        }

        Float_t tantheta2 = 1e10;
        if (eta2 < -1e-10 || eta2 > 1e-10){
            Float_t expTmp = 1.0-eta2+eta2*eta2/2-eta2*eta2*eta2/6+eta2*eta2*eta2*eta2/24;
            tantheta2 = 2.0 * expTmp / ( 1.0 - expTmp*expTmp);
        }

        Float_t e1squ = m0_1 * m0_1 + pt1 * pt1 * (1.0 + 1.0 / tantheta1 / tantheta1);
        Float_t e2squ = m0_2 * m0_2 + pt2 * pt2 * (1.0 + 1.0 / tantheta2 / tantheta2);

        // fold onto 0...pi
        Float_t deltaPhi = TMath::Abs(phi1 - phi2);
        while (deltaPhi > TMath::TwoPi())
            deltaPhi -= TMath::TwoPi();
        if (deltaPhi > TMath::Pi())
            deltaPhi = TMath::TwoPi() - deltaPhi;

        Float_t cosDeltaPhi = 0;
        if (deltaPhi < TMath::Pi()/3)
            cosDeltaPhi = 1.0 - deltaPhi*deltaPhi/2 + deltaPhi*deltaPhi*deltaPhi*deltaPhi/24;
        else if (deltaPhi < 2*TMath::Pi()/3)
            cosDeltaPhi = -(deltaPhi - TMath::Pi()/2) + 1.0/6 * TMath::Power((deltaPhi - TMath::Pi()/2), 3);
        else
            cosDeltaPhi = -1.0 + 1.0/2.0*(deltaPhi - TMath::Pi())*(deltaPhi - TMath::Pi()) - 1.0/24.0 * TMath::Power(deltaPhi - TMath::Pi(), 4);

        Float_t mass2 = m0_1 * m0_1 + m0_2 * m0_2 + 2 * ( TMath::Sqrt(e1squ * e2squ) - ( pt1 * pt2 * ( cosDeltaPhi + 1.0 / tantheta1 / tantheta2 ) ) );
        //   Printf(Form("%f %f %f %f %f %f %f %f %f", pt1, eta1, phi1, pt2, eta2, phi2, m0_1, m0_2, mass2));
        return mass2;
    }

};

#endif






















