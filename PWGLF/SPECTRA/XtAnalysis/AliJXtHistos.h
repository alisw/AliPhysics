/**************************************************************************
 * AliJXtHistos.h
 * This class encapsulated all histograms that the analysis provides
 *
 * contact: Sami Räsänen
 *          University of Jyväskylä, Finland
 *          sami.s.rasanen@jyu.fi
 **************************************************************************/

#ifndef ALIJXTHISTOS_H
#define ALIJXTHISTOS_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TFile.h>
#include <TList.h>

#define kMaxNoCentrBin 10   // Maximum no of centrality bins defined in JCard.h

class AliJCard;

using namespace std;

class AliJXtHistos {

    public:
        AliJXtHistos(AliJCard* cardP); //constructor
        virtual ~AliJXtHistos(){delete fhistoList;}    //destructor
        AliJXtHistos(const AliJXtHistos& obj);
        AliJXtHistos& operator=(const AliJXtHistos& obj);

        // create histograms 
        void CreateXtHistos();

        TList *GetHistoList() { return fhistoList; } //return the list of histograms

        bool UseDirectory() const { return fUseDirectory; } // Are there directories in the final root file, yes/no
        void UseDirectory(bool b) { fUseDirectory=b; } // Decide weather to use directories

        TDirectory * MakeDirectory(TString name); // Make new directory into final root file
        TDirectory * JumpToDefaultDirectory(); // move into default directory at the final root file
    
        // Fill various histograms
        void FillCentralityHistos(double fcent, int fCentBin){
            // Fill centrality histograms
            fhCentr->Fill(fcent);
            fhiCentr->Fill(1.0*fCentBin);
        }
        void FillRawVertexHisto(double zVert){fhZVertRaw->Fill(zVert);} // fill raw z-vertex histogram
        void FillAcceptedVertexHisto(double zVert, int centBin){fhZVert[centBin]->Fill(zVert);} // fill accepted z-vertex histograms
        void FillInclusiveHistograms(double pT, double xT, double eta, double phi, double effCorr, int centBin); // Fill histograms
        void FillIsolatedHistograms(double pT, double xT, double eta, double phi, double effCorr, int centBin);  // Fill histograms
        void FillInclusiveConeActivities(double pT, double sumPt){fhConeActivity->Fill(pT,sumPt);} // Fill TProfile
        void FillIsolatedConeActivities(double pT, double sumPt){fhConeActivityIsolated->Fill(pT,sumPt);}  // Fill TProfile
    
    protected:
        bool   fUseDirectory; // to create sub-directories in the final results, used more in case of JCORRAN
	TDirectory * fTopDirectory; // top directory, different analysis (JCORRAN) in sub-directories.
	
	AliJCard  *fCard; // parameters in card
    
        double fmaxEtaRange; // charged track eta acceptance
    
        TList *fhistoList; // list of histograms

    private:
        char  fhname[40], fhtit[40]; // dummy variables to create histogram names and titles

        TH1D *fhChargedPt[kMaxNoCentrBin], *fhInvariantChargedPt[kMaxNoCentrBin], *fhChargedPtNoCorr[kMaxNoCentrBin]; // inclusive pT distribution - dN/dpT and 1/pT dN/dpT, last without efficiency correction
        TH1D *fhIsolatedChargedPt[kMaxNoCentrBin], *fhInvariantIsolatedChargedPt[kMaxNoCentrBin]; // the same as above but for isolated pT
        TH1D *fhChargedXt[kMaxNoCentrBin], *fhInvariantChargedXt[kMaxNoCentrBin], *fhChargedXtNoCorr[kMaxNoCentrBin]; // The same as inclusive, but for xT
        TH1D *fhIsolatedChargedXt[kMaxNoCentrBin], *fhInvariantIsolatedChargedXt[kMaxNoCentrBin]; // the same, but for isolated xT
        TH1D *fhChargedPhi[kMaxNoCentrBin], *fhChargedPhiNoCorr[kMaxNoCentrBin], *fhIsolatedChargedPhi[kMaxNoCentrBin]; // phi -distributions - inclusive, without efficiency correction and for isolated
        TH1D *fhChargedEta[kMaxNoCentrBin], *fhChargedEtaNoCorr[kMaxNoCentrBin], *fhIsolatedChargedEta[kMaxNoCentrBin]; // the same as above but for eta

        TProfile *fhConeActivity;          // pT sum in cone, to be compared to the ALICE UE results
        TProfile *fhConeActivityIsolated;  // activity for isolated particles

		TH1D *fhZVertRaw; // raw z-vertex distribution
		TH1D *fhZVert[kMaxNoCentrBin]; // accepted z-vertex as a function of centrality
		TH1D *fhCentr; // centrality
		TH1D *fhiCentr; // centrality, accepted
		TH1D *fhEventPerRun; // number of events per run
		TH2D* fhVertexZTriggVtx; // z-vertex distribution of triggered events, needed somewhere?

};

#endif
