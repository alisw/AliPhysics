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
        void SetCentralityBinsHistos( vector<double> centralityClasses ) {
            CentBin=centralityClasses;
            fNCentBin = CentBin.size();
        }
        static int fNCentBin;
        static vector<double> CentBin;//[NCENT+1]; //8

        AliJHistManager * fHMG; //! Histogram manager
        TString sMngrName = "jcdijet"; // Histogram manager name tag
        AliJBin fHistCentBin;   //! Centrality bin
        AliJBin fJetBin;        //! Jet bin
        //===================================================
        // Event/Track histograms
        //===================================================
        AliJTH1D fh_events;       //! // for counting events, jets, dijets and so on.
        AliJTH1D fh_info;         //! // General information about the run.
        AliJTH1D fh_centrality;   //! // centrality histogram
        AliJTH1D fh_zvtx;         //! // z-vertex histogram

        AliJTH1D fh_pt;     //! // for pt dist of tracks
        AliJTH1D fh_eta;    //! // for eta dist of tracks
        AliJTH2D fh_etaPhi; //! // for (eta,phi) dist of tracks
        AliJTH1D fh_phi;    //! // for phi dist of tracks

        AliJTH1D fh_rho;        //! // for event pt density
        AliJTH1D fh_rhom;       //! // for event mt density

        AliJTH1D fh_jetPt;      //! // for pt dist of jets
        AliJTH1D fh_jetEta;     //! // for eta dist of jets
        AliJTH1D fh_jetPhi;     //! // for phi dist of jets
        AliJTH2D fh_jetEtaPhi;  //! // for (eta,phi) dist of jets
        AliJTH1D fh_jetArea;    //! // for jet area spectrum
        AliJTH1D fh_jetAreaRho; //! // for jet area*pt spectrum

        AliJTH1D fh_dijetInvM;                //! // for dijet invariant mass
        AliJTH1D fh_dijetPtPair;              //! // for dijet pt
        AliJTH1D fh_dijetDeltaPhi;            //! // for dijet deltaPhi
        AliJTH1D fh_dijetPtPairDeltaPhiCut;   //! // for dijet pt after deltaPhi cut
        AliJTH1D fh_dijetInvMDeltaPhiCut;     //! // for dijet invariant mass after deltaPhi cut
        AliJTH1D fh_dijetDeltaPhiWithCut;     //! // for dijet delta phi after deltaPhi cut

        AliJTH1D fh_responseInfo;             //! // for counting response related things.
        AliJTH1D fh_jetResponseDeltaR;        //! // true jet vs. detector jet deltaR
        AliJTH1D fh_jetDeltaRMin;             //! // Minimum deltaR between jets
        AliJTH1D fh_jetBGSubtrDeltaR;         //! // DeltaR between BG subtr jet and raw jet.
        AliJTH2D fh_jetResponse;              //! // Jet response matrix
        AliJTH2D fh_dijetResponse;            //! // Dijet response matrix
        AliJTH2D fh_dijetResponseDeltaPhiCut; //! // Dijet response matrix with deltaPhi cut
};

#endif //ALIJCDIJETHISTOS_H
