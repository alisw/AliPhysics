/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* Add a description of your MPI analysis */

#ifndef ALIANALYSISTASKGENMCKNOUE_H
#define ALIANALYSISTASKGENMCKNOUE_H
// ROOT includes

#include <TObject.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <TH1F.h>
#include <THnSparse.h>
#include <TArrayD.h>
#include <TParticle.h>
#include <TClonesArray.h>
#include <TList.h>
#include <TObjArray.h>

// AliRoot includes

#include <AliAnalysisTaskSE.h>
#include <AliMCEvent.h>
#include <AliStack.h>
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliInputEventHandler.h"

class TList;
class TH1F;
class TH1I;
class TH2D;
class TH3D;
class TParticle;
class AliStack;
class AliVVertex;
class AliVParticle;


class AliAnalysisTaskGenMcKnoUe : public AliAnalysisTaskSE
{
	public:

		AliAnalysisTaskGenMcKnoUe();
		AliAnalysisTaskGenMcKnoUe(const char *name);

		virtual      ~AliAnalysisTaskGenMcKnoUe();
		virtual void UserCreateOutputObjects();
		virtual void UserExec(Option_t* option);
		virtual void Terminate(Option_t* option);
		virtual void SetPtMin(Double_t val){fPtMin = val;}
		virtual void SetEtaMax(Double_t val){fEtaCut = val;}
		virtual void SetEtaRhoMax(Double_t val){fEtaCutRho = val;} 
		virtual void SetnEtaBinsRho(Double_t val){fnEtaBinsRho = val;}
		virtual void SetnPhiBinsRho(Double_t val){fnPhiBinsRho = val;}
		virtual void SetIsPP(Bool_t val){fIsPP = val;}
		virtual void SetGenerator(Int_t val){fGenerator = val;}
		virtual void SetIsFirstPart(Bool_t val){fFirstPart = val;}
		virtual void SetParametersRho(Double_t etarange); 

	private:

		void       GetGenLeadingObject();
		void       GetGenUEObservables();
		Int_t      GetPidCode(Int_t pdgCode);
		void       MakeALICE3Analysis();
		void       MakeALICE3AnalysisP2();
		Double_t   DeltaPhi(Double_t phia, Double_t phib,Double_t rangeMin = -TMath::Pi()/2.0, Double_t rangeMax = 3.0*TMath::Pi()/2.0 );

	protected:

		Bool_t IsMCEventSelected(TObject* obj);
		AliMCEvent*  fMC;                   //!<! MC event
		AliInputEventHandler*  fMcHandler;  //!<! MCEventHandler
		AliStack*  fMCStack;

		Bool_t fFirstPart;
		Int_t fGenerator;
		Double_t fEtaCut;
		Double_t fEtaCutRho;
		Int_t fnEtaBinsRho;
		Int_t fnPhiBinsRho;
		Bool_t   fIsPP;
		Double_t fPtMin;

		Double_t fGenLeadPhi; 
		Double_t fGenLeadEta;
		Double_t fGenLeadPt;
		Int_t    fGenLeadIn;
		// UE 
		TH1D * hPtLeadingTrue;
		TH2D * hPtLVsV0A;
		TH2D * hnchmpi;
		TH3F * hnchmpirho;
		TH2D * hnchrho;
		TH2D * hmpirho;
		TH3D * hPtVsPtLeadingTrue[3];
		TH1D * hphiKNO;
		TH1D * hphiKNO1;
		TH1D * hphiKNO2;
		TH1D * hNchTforKNOana;
		TH1D * hNchTforKNOanaMin;
		TH1D * hNchTforKNOanaMax;

		TH1D * hPtLeadingRho[3];// 0: all, 1: low rho, 2: high rho
		TH2D * hNchRho[3];
		TH1D * hEtaLeadingRho[3];
		TH2D * hDetaDphiRho[3];
		TH2D * hDetaDphiRhoWideEta[3];
		TH3D * hNchPtPidRho[4][3]; // region, rho 

		TList*  fOutputList;    //!<! Output list of objects

		AliAnalysisTaskGenMcKnoUe(const AliAnalysisTaskGenMcKnoUe&);
		AliAnalysisTaskGenMcKnoUe& operator=(const AliAnalysisTaskGenMcKnoUe&);

		ClassDef(AliAnalysisTaskGenMcKnoUe, 1);
};

#endif
