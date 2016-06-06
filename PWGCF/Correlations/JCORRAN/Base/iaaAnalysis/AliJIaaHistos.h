/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

// Interface for all analysis histograms.
// This is the minimum amount of histograms all analysis must have defined

//===========================================================
// AliJIAAHistos.h
//===========================================================

#ifndef ALIJIAAHISTOS_H
#define ALIJIAAHISTOS_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TFile.h>
#include <TList.h>
#include <TLorentzVector.h>

#include "../AliJConst.h"
#include "../AliJHistogramInterface.h"
#include "../AliJHistManager.h"

class AliJCard;
class AliJBaseTrack;
class AliJTrack;

using namespace std;

class AliJIaaHistos : public AliJHistogramInterface {
  
public:
    AliJIaaHistos(AliJCard* cardP); //constructor
    virtual ~AliJIaaHistos();    //destructor
    AliJIaaHistos(const AliJIaaHistos& obj);
    AliJIaaHistos& operator=(const AliJIaaHistos& obj);

    // create histograms
    void CreateCorrelationHistos();
    void CreateEventTrackHistos();

    void CreatePtCorrHistos();
    void CreateRunByRunHistos(int runID, int runcounter) const;

    void ReadInclusiveHistos(const char *inclusFileName);

    bool Is2DHistosEnabled(){ return fenable2DHistos; }
    void Set2DHistoCreate(bool isenable) { fenable2DHistos = isenable; }
    void SetAcceptanceCorrectionQA(bool isenable) { fEnableAcceptanceQAHistos = isenable; }

    AliJTProfile fhMixStat; // comment me

    //==Pt stat. fcorrelations ===============================================
    AliJTH1D fhPtNear; // comment me
    AliJTH1D fhPtFar ; // comment me

    //assorted
    AliJTH1D fhPhi; // comment me
    AliJTH1D fhDphiAssoc; // comment me


    AliJTH2D fhDphiDetaPta;      // 2D histogram from deltaPhi-deltaEta plane in pta bins

    AliJTH1D fhDetaNearMixAcceptance;   // Mixed event uncorrected deltaEta histogram for acceptance correction
    AliJTH1D fhDeta3DNearMixAcceptance; // Mixed event uncorrected deltaEta histogram in 3D near side for acceptance corfection

    AliJTH1D fhDEtaNear; // comment me
    AliJTH1D fhDEtaNearM; // comment me

    AliJTH1D fhDEtaFar; // comment me
    AliJTH1D fhIphiTrigg; // comment me
    AliJTH1D fhIetaTrigg; // comment me
    AliJTH1D fhIphiAssoc; // comment me
    AliJTH1D fhIetaAssoc; // comment me
    AliJTH1D fhTriggPtBin; // comment me
    AliJTH1D fhTriggMult; // comment me


    //===================================================
    // Event/Track histograms
    //===================================================
    AliJTH1D fhLPpt; // comment me
    AliJTH1D fhLPpairPt; // comment me
    AliJTH1D fhChargedPt, fhChargedPtNoCorr; // comment me
    AliJTH1D fhChargedPtJacek; // comment me
    AliJTH1D fhChargedPtJacekEta; // 3 bins in eta
    AliJTH1D fhChargedPtFiete; // comment me
    AliJTH1D fhVdelta2,fhVdelta3, fhVN; // comment me
    AliJTProfile fhTrackingEfficiency; // comment needed
    AliJTProfile fpV2, fpV3, fpVdeltaNorm; // comment me
    AliJTH1D fhChargedEta; // comment me
    AliJTH1D fhLPeta; // comment me
    AliJTH1D fhAssocMult; // comment me
    AliJTH1D fhChargedMult, fhChargedMultCut; // comment me
    AliJTH2D fhChargedMultCent; // comment me

    AliJTH1D  fhV0AMult; // comment needed

    AliJTH1D fhZVertRawErr; // comment me
    AliJTH1D fhZVert; // comment me
    AliJTH1D fhCentr; // comment me
    AliJTH1D fhiCentr; // comment me
    AliJTH1D fhEventPerRun; // comment me


  
protected:
    double fmaxEtaRange;                       // maximum eta range
    double fmaxTriggEtaRange;                  // should be the same as above. Use for GeoAccCorr
    double fLowRange, fHighRange;              // lower and upper range for dphi histos
    bool fenable2DHistos;                      // enable the filling of two dimensional histograms
    bool fEnableAcceptanceQAHistos;            // enable the filling of acceptance correction quality control histograms

    int fNJacek;        // Number of bins in Jacek binning
    double *fPttJacek;  // Bin borders in Jacek binning
    int fNEta;          // Number of bins in eta binning
    double *fEta;       // Bin borders in eta binning

};

#endif






















