//===========================================================
// AliJEbeHistos.h
//  
//   J
//===========================================================

#ifndef ALIJEBEHISTOS_H
#define ALIJEBEHISTOS_H

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

#include  "AliJConst.h"


class AliJCard;

using namespace std;

class AliJEbeHistos {

    public:
        AliJEbeHistos(AliJCard* cardP); //constructor
        virtual ~AliJEbeHistos(){delete fhistoList;}    //destructor
      AliJEbeHistos(const AliJEbeHistos& obj);
      AliJEbeHistos& operator=(const AliJEbeHistos& obj);

        // create histograms 
	void CreateUnfoldingHistos();
        TList *GetHistoList() { return fhistoList; }

        void UseDirectory(bool b) { fUseDirectory=b; }
        bool UseDrectory(){ return fUseDirectory; }

        TDirectory * MakeDirectory(TString name){
            JumpToDefalutDirectory();
            TDirectory * dir = gDirectory->GetDirectory(name);
            if( !dir ) dir = gDirectory->mkdir(name);
            dir->cd();
            return dir;
        }
        TDirectory * JumpToDefalutDirectory(){
            fTopDirectory->cd();
            return gDirectory;
        }

    public:
        AliJCard  *fcard; // comment me
        char  fhname[50], fhtit[50]; // comment me
        TString fhtyp[3]; // comment me


		//Unfolding
		TH1D* fhVnObsVector[kMaxNoCentrBin][kNHarmonics];
		TH2D* fhResponseDist[kMaxNoCentrBin][kNHarmonics];
		TH1D* fhMultiCount[kMaxNoCentrBin];
		TH1D* fhVnObsEP[kMaxNoCentrBin][kNHarmonics];
		TH1D* fhCosndPhiPt[kMaxNoCentrBin][kNHarmonics];
		TH1D* fheCosndPhiPt[kMaxNoCentrBin][kNHarmonics];
		TH1D* fhCounter[kMaxNoCentrBin][kNHarmonics];
		TH1D* fhEventPlane[kMaxNoCentrBin][kNHarmonics];
		TH1D *fhEPCosndPhi[kNHarmonics];
		TH1D *fhEPCosndPhi2[kNHarmonics];

		TH1D *fhQvectorV0[kMaxNoCentrBin][kNHarmonics];
		TH1D *fhQvectorV0A[kMaxNoCentrBin][kNHarmonics];
		TH1D *fhQvectorV0C[kMaxNoCentrBin][kNHarmonics];
		TH2D *fhQvectorCorrelation[kMaxNoCentrBin][kNHarmonics];
		TH2D *fhVnObsVsQvectorCorrelation[kMaxNoCentrBin][kNHarmonics];

		TH1D *fhResolution[kMaxNoCentrBin][kNHarmonics];
		


	protected:
		double fmaxEtaRange;                       // evident
		double fmaxTriggEtaRange;                  // should be the same as above. Use for GeoAccCorr
		double ftriggFiducCut;                     // fiducial cut for the trigger part in eta. Not in use I think (Jan) 

		TList *fhistoList; // comment me
		bool   fUseDirectory;
		TDirectory * fTopDirectory;

};

#endif






















