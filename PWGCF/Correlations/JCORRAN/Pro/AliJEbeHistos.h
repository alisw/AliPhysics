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
#include  "AliJHistManager.h"


class AliJCard;

using namespace std;

class AliJEbeHistos {

    public:
        AliJEbeHistos(); //constructor
        AliJEbeHistos(AliJCard* cardP); //constructor
        virtual ~AliJEbeHistos();    //destructor
      AliJEbeHistos(const AliJEbeHistos& obj);
      AliJEbeHistos& operator=(const AliJEbeHistos& obj);

        // create histograms 
      void CreateUnfoldingHistos();
      //TList *GetHistoList() { return fhistoList; }

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

        
        AliJHistManager * fHMG; // Histogram manager
        AliJBin  fCentBin;
        AliJBin  fHarmonicBin;


		//Unfolding
		AliJTH1D  fhVnObsVector;
		AliJTH1D  fhVnObsVectorAfterSelection;
		AliJTH2D  fhResponseDist;
		AliJTH1D  fhMultiCount;
		AliJTH1D  fhVnObsEP;
		AliJTH1D  fhCosndPhiPt;
		AliJTH1D  fheCosndPhiPt;
		AliJTH1D  fhCounter;
		AliJTH1D  fhEventPlane;
		AliJTH1D  fhEPCosndPhi;
		AliJTH1D  fhEPCosndPhi2;

		AliJTH1D  fhQvectorV0;
		AliJTH1D  fhQvectorV0A;
		AliJTH1D  fhQvectorV0C;
		AliJTH2D  fhQvectorCorrelation;
		AliJTH2D  fhVnObsVsQvectorCorrelation;

		AliJTH1D  fhResolution;
		


	protected:
		double fmaxEtaRange;                       // evident
		double fmaxTriggEtaRange;                  // should be the same as above. Use for GeoAccCorr
		double ftriggFiducCut;                     // fiducial cut for the trigger part in eta. Not in use I think (Jan) 

		//TList *fhistoList; // comment me
		bool   fUseDirectory;
		TDirectory * fTopDirectory;

};

#endif






















