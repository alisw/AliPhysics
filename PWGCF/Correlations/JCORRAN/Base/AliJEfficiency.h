#ifndef ALIJEFFICIENCY_H
#define ALIJEFFICIENCY_H
// Class to Load and Get efficiency inforamtion
// ..
// TODO

#include <TString.h>
#include <TFile.h>
#include "AliJTrackCut.h"
#include "AliJRunTable.h"
#include "AliJConst.h"
#include <TGraphErrors.h>
#include <TAxis.h>
#include <iostream>
using namespace std;

class AliJEfficiency{
    public:
        enum Mode { kNotUse, kPeriod, kRunNumber, kAuto };
        enum Type { kRE, kMC, kMerge };
        AliJEfficiency();
        AliJEfficiency(const AliJEfficiency& obj);
        AliJEfficiency& operator=(const AliJEfficiency& obj);

        void SetMode( int i ){ fMode = i; }

        void SetDataPath(TString s ){ fDataPath=s; }
        void SetEffFile(TString s ){ fInputRootName=s; }
        void SetName(TString s ){ fName=s; }
        void SetPeriod(int period){ fPeriod = period; }
        void SetPeriod(TString s){ fPeriodStr = s; }
        void SetMCPeriod(TString s){ fMCPeriodStr = s; }
        void SetRunNumber( Long64_t runnum ){ fRunNumber=runnum; }
        void SetTag(TString s){ fTag=s; }

        TString GetName() const { return fName; }
        double GetCorrection( double pt, int icut, double cent ) const ;
        TString GetEffName() ;
        TString GetEffFullName() ;
        bool   Load();
        void   PrintOut() const {
            cout<<fInputRootName<<endl;
        }
        void Write();

    private:
        int      fMode;             // Mode. see enum Mode
        int      fPeriod;           // Data Period index
        AliJTrackCut fTrackCut;     // Track Cut Object. TODO:why not pointer?
        AliJRunTable fRunTable;     // run Table. TODO:why not pointer?

        TString fDataPath;          // locaction of eff files
        TString fName;              // name of efficiency. usually empty
        TString fPeriodStr;         // DATA period
        TString fMCPeriodStr;       // MC period
        Long64_t fRunNumber;        // Runnumber
        TString fTag;               // Tags to distinguish special eff file
        TString fInputRootName;     // name of input

        TFile * fInputRoot;         // input file  
        TDirectory * fEffDir[3];    // root directory of efficiency. only second item of fEffDir with "Efficiency" is being used.
        TGraphErrors * fCorrection[20][20][20]; // Storage of Correction factor 
        TAxis * fCentBin;     // Bin of Centrality. replace with AliJBin?
};
#endif
