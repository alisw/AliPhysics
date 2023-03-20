/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//______________________________________________________________________________
// Basic histogram implimentation via AliJHistogramInterface.
// author: O. Saarimaki, D.J. Kim (dong.jo.kim@cern.ch)
// ALICE Group University of Jyvaskyla, Finland
// Note: Adapted for AliAnalysisTaskSE
//////////////////////////////////////////////////////////////////////////////
#ifndef ALIJCDIJETHISTOS_H
#define ALIJCDIJETHISTOS_H

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

#include "AliJHistogramInterface.h"


using namespace std;

class AliJCDijetHistos : public AliJHistogramInterface
{

    public:

        AliJCDijetHistos(); //constructor
        AliJCDijetHistos(const AliJCDijetHistos& obj); // copy constructor
        virtual ~AliJCDijetHistos();    //destructor
        AliJCDijetHistos& operator=(const AliJCDijetHistos& obj); // equal sign operator

        void CreateEventTrackHistos();
        void SetName(TString slMngrName) { sMngrName = slMngrName.Data(); }
        TString GetName() { return sMngrName; }

        static int GetCentralityClass(Double_t);
        static int GetDijetMClass(Double_t);
        void SetCentralityBinsHistos( vector<double> centralityClasses ) {
            CentBin=centralityClasses;
            fNCentBin = CentBin.size();
        }
        void SetDijetMBinsHistos( TString sMBins ) {
            fSMBins = sMBins;
            std::vector<TString> ar = Tokenize( sMBins, "\t ,");
            fnNewBinsDijet1=ar.size()-1;
            for (unsigned i=0; i<ar.size(); i++) {
                dijetMBin.push_back(ar.at(i).Atof());
            }
        }
        void SetNJetClasses( int ljetclasses ) {
            fNJetClasses = ljetclasses;
        }
        static int fNCentBin;
        static vector<double> CentBin;//[NCENT+1]; //8
        TString fSMBins;//
        static int fNJetClasses;
        static int fnNewBinsDijet1;//
        static vector<double> dijetMBin;//

        AliJHistManager * fHMG; //! Histogram manager
        TString sMngrName = "jcdijet"; // Histogram manager name tag
        AliJBin fHistCentBin;   //! Centrality bin
        AliJBin fJetBin;        //! Jet bin
        AliJBin fMBin;        //! Dijet mass bin
        //===================================================
        // Event/Track histograms
        //===================================================
        AliJTH1D fh_events;       //! // for counting events, jets, dijets and so on.
        AliJTH1D fh_eventSel;     //! // for studying event selection.
        AliJTH1D fh_info;         //! // General information about the run.
        AliJTH1D fh_centrality;   //! // centrality histogram
        AliJTH1D fh_zvtx;         //! // z-vertex histogram
        AliJTH1D fh_nch;          //! // Number of charged tracks
        AliJTH1D fh_hisWeight;    //! // Histogram weights

        AliJTH1D fh_pt;     //! // for pt dist of tracks
        AliJTH1D fh_ptPosEta; //! // for pt dist of tracks with positive eta
        AliJTH1D fh_ptNegEta; //! // for pt dist of tracks with negative eta
        AliJTH1D fh_eta;    //! // for eta dist of tracks
        AliJTH2D fh_etaPhi; //! // for (eta,phi) dist of tracks
        AliJTH1D fh_phi;    //! // for phi dist of tracks

        AliJTH1D fh_rho;        //! // for event pt density
        AliJTH1D fh_rhoHighPt;  //! // for event pt density when high pt jets in event
        AliJTH1D fh_rhom;       //! // for event mt density
        AliJTH1D fh_rhomHighPt; //! // for event mt density when high pt jets in event
        AliJTH1D fh_rhoLin;        //! // for event pt density. Linear bins
        AliJTH1D fh_rhoLinHighPt;  //! // for event pt density when high pt jets in event. Linear bins
        AliJTH1D fh_rhomLin;       //! // for event mt density. Linear bins
        AliJTH1D fh_rhomLinHighPt; //! // for event mt density when high pt jets in event. Linear bins
        AliJTH1D fh_coveredRatio;        //! // for event pt density

        AliJTH1D fh_jetPt;      //! // for pt dist of jets
        AliJTH1D fh_jetPt_ALICE;//! // for pt dist of jets, with ALICE bins
        AliJTH1D fh_jetPtLeadSublead_ALICE;//! // for pt dist of leading and subleading jets, with ALICE bins
        AliJTH1D fh_jetPtLeadSubleadDeltaPhi_ALICE;//! // for pt dist of leading and subleading jets, with ALICE bins
        AliJTH1D fh_jetPtLeadSubleadMBin_ALICE;//! // for pt dist of leading and subleading jets, with ALICE bins
        AliJTH1D fh_jetPtLeadSubleadDeltaPhiMBin_ALICE;//! // for pt dist of leading and subleading jets, with ALICE bins
        
        AliJTH1D fh_jetPtTransBGSub; //! // for pt dist of jets with BG subtraction only for transverse momentum
        AliJTH1D fh_jetN;      //! // for number of jets
        AliJTH1D fh_jetEta;     //! // for eta dist of jets
        AliJTH1D fh_jetPhi;     //! // for phi dist of jets
        AliJTH2D fh_jetEtaPhi;  //! // for (eta,phi) dist of jets
        AliJTH2D fh_randConeEtaPhi; //! // for (eta,phi) dist of random cones
        AliJTH1D fh_jetArea;    //! // for jet area spectrum
        AliJTH1D fh_jetAreaRho; //! // for jet area*pt spectrum
        AliJTH1D fh_deltaPt; //! // for delta-pt distribution
        AliJTH1D fh_maxJetptOverPtHard; //! // for jet pt / pt_hard bin
        AliJTH1D fh_ptHard; //! // for pt_hard
        AliJTH1D fh_pythiaSigma; //! // for pythia sigma_gen
        AliJTH1D fh_pythiaTrial; //! // for pythia trial

        AliJTH1D fh_dijetInvM;                     //! // for dijet invariant mass
        AliJTH1D fh_dijetInvMLin;                  //! // for dijet invariant mass
        AliJTH1D fh_dijetInvMTrunc;                //! // for dijet invariant mass truncated above and below
        AliJTH1D fh_dijetInvMTrunc2;               //! // for dijet invariant mass truncated above and below
        AliJTH1D fh_dijetPtPair;                   //! // for dijet pt
        AliJTH1D fh_dijetDeltaPhi;                 //! // for dijet deltaPhi
        AliJTH1D fh_dijetCosDeltaPhi;              //! // for cos(dijet deltaPhi)
        AliJTH1D fh_dijetDeltaEta;                 //! // for dijet deltaPhi
        AliJTH1D fh_dijetCoshDeltaEta;             //! // for Cosh(dijet deltaPhi)
        AliJTH1D fh_dijetSqrtGeometry;             //! // for sqrt(Cosh(dijet deltaPhi) - cos(dijet deltaPhi))
        AliJTH1D fh_dijetSqrt2pt12;                //! // for sqrt(2*pt1*pt2)
        AliJTH1D fh_dijetPtPairDeltaPhiCut;        //! // for dijet pt after deltaPhi cut
        AliJTH1D fh_dijetInvMDeltaPhiCut;          //! // for dijet invariant mass after deltaPhi cut
        AliJTH1D fh_dijetInvMDeltaPhiCutLin;       //! // for dijet invariant mass after deltaPhi cut
        AliJTH1D fh_dijetInvMDeltaPhiCutTrunc;     //! // for dijet invariant mass after deltaPhi cut truncated above and below
        AliJTH1D fh_dijetInvMDeltaPhiCutTrunc2;    //! // for dijet invariant mass after deltaPhi cut truncated above and below
        AliJTH1D fh_dijetPtPairDeltaPhiCutWithCut; //! // for dijet pt after deltaPhi cut
        AliJTH1D fh_dijetDeltaPhiWithCut;          //! // for dijet delta phi after deltaPhi cut
        AliJTH1D fh_dijetCosDeltaPhiWithCut;       //! // for cos(dijet deltaPhi) (with deltaPhi cut)
        AliJTH1D fh_dijetDeltaEtaWithCut;          //! // for dijet deltaPhi (with deltaPhi cut)
        AliJTH1D fh_dijetCoshDeltaEtaWithCut;      //! // for Cosh(dijet deltaPhi) (with deltaPhi cut)
        AliJTH1D fh_dijetSqrtGeometryWithCut;      //! // for sqrt(Cosh(dijet deltaPhi) - cos(dijet deltaPhi)) (with deltaPhi cut)
        AliJTH1D fh_dijetSqrt2pt12WithCut;         //! // for sqrt(2*pt1*pt2) (with deltaPhi cut)

        AliJTH1D fh_responseInfo;                  //! // for counting response related things.
        AliJTH1D fh_jetResponseDeltaR;             //! // true jet vs. detector jet deltaR
        AliJTH1D fh_jetResponseDeltaRClosest;      //! // true jet vs. detector jet deltaR no limits
        AliJTH1D fh_jetResponseDeltaPt;            //! // true jet vs. detector jet pt normed with true jet pt
        AliJTH1D fh_jetDeltaRMin;                  //! // Minimum deltaR between jets
        AliJTH1D fh_jetBGSubtrDeltaR;              //! // DeltaR between BG subtr jet and raw jet.
        AliJTH2D fh_jetResponse;                   //! // Jet response matrix
        AliJTH2D fh_jetResponse_ALICE;             //! // Jet response matrix with ALICE bins
        AliJTH2D fh_jetResponse_ALICEScalar;       //! // Jet response matrix with ALICE bins (bg subtr with scalar rho*A quantity)
        AliJTH2D fh_deltaPtResponse;               //! // delta-pt response matrix 
        AliJTH2D fh_deltaPtResponse_ALICE;         //! // delta-pt response matrix with ALICE bins
        AliJTH2D fh_deltaPtResponseEvery;          //! // delta-pt response matrix, filled for every true bin
        AliJTH2D fh_deltaPtResponseEvery_ALICE;    //! // delta-pt response matrix with ALICE bin, filled for every true bin
        AliJTH2D fh_dijetResponse;                 //! // Dijet response matrix
        AliJTH2D fh_dijetResponseLin;              //! // Dijet response matrix
        AliJTH2D fh_dijetResponseLinNoMatching;    //! // Dijet response matrix

        AliJTH1D fh_doubleConeM;              //! // Double cone invariant mass
        AliJTH1D fh_doubleConeMAlt;              //! // Double cone invariant mass
        AliJTH1D fh_jet2Cone1Dist;            //! // Distance between jet 2 and cone 1
        AliJTH1D fh_jet1Cone2Dist;            //! // Distance between jet 1 and cone 2
        AliJTH1D fh_jet1Cone2AltDist;            //! // Distance between jet 1 and cone 2
        AliJTH1D fh_localRho1;                //! // Localrho1
        AliJTH1D fh_localRho2;                //! // Localrho2
        AliJTH1D fh_localRho2Alt;                //! // Localrho2
        AliJTH1D fh_deltaRho1;                //! // Difference of rho and localrho1
        AliJTH1D fh_deltaRho2;                //! // Difference of rho and localrho2
        AliJTH1D fh_deltaRho2Alt;                //! // Difference of rho and localrho2
        AliJTH1D fh_deltaLocalRho;            //! // Difference of localrho1 and localrho2
        AliJTH1D fh_deltaLocalRhoAlt;            //! // Difference of localrho1 and localrho2
        AliJTH1D fh_dijetdeltaM5;             //! // Dijet deltaM
        AliJTH1D fh_dijetdeltaM5Alt;             //! // Dijet deltaM
        AliJTH1D fh_dijetdeltaM5Binned;             //! // Dijet deltaM
        AliJTH1D fh_dijetdeltaMScaled;             //! // Dijet deltaM
        AliJTH1D fh_dijetdeltaMScaledBinned;             //! // Dijet deltaM
        AliJTH1D fh_dijetdeltaM5NearCone;             //! // Dijet deltaM
        AliJTH1D fh_dijetdeltaM5NearConeAlt;             //! // Dijet deltaM
        AliJTH1D fh_dijetMLocalRho;            //! // Dijet deltaM with localrho
        AliJTH1D fh_dijetMLocalRhoAlt;            //! // Dijet deltaM with localrho
        AliJTH1D fh_dijetMLocalRhoAltMjjBins;            //! // Dijet deltaM with localrho
        AliJTH2D fh_deltaMResponse;           //! // deltaM response matrix 
        AliJTH2D fh_dijetResponseTrunc;            //! // Dijet response matrix truncated from above and below
        AliJTH2D fh_dijetResponseTrunc2;            //! // Dijet response matrix truncated from above and below
        AliJTH2D fh_dijetResponseDeltaPhiCut;      //! // Dijet response matrix with deltaPhi cut
        AliJTH2D fh_dijetResponseDeltaPhiCutLin;      //! // Dijet response matrix with deltaPhi cut
        AliJTH2D fh_dijetResponseDeltaPhiCutLinNoMatching;      //! // Dijet response matrix with deltaPhi cut
        AliJTH2D fh_dijetResponseDeltaPhiCutTrunc; //! // Dijet response matrix with deltaPhi cut truncated from above and below
        AliJTH2D fh_dijetResponseDeltaPhiCutTrunc2; //! // Dijet response matrix with deltaPhi cut truncated from above and below
};

#endif //ALIJCDIJETHISTOS_H
