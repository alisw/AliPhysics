/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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
//
// Task for the analysis of the production of 3He as a function of multiplicity
//
// Author:
//  Sebastian Hornung <Sebastian.Hornung@cern.ch>

#ifndef AliAnalysisTask3HeMC_h
#define AliAnalysisTask3HeMC_h

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

class AliAnalysisUtils;
class AliEventCuts;
class AliPIDResponse;
class AliAODMCHeader;
class TClonesArray;
class TH2D;
class THnSparse;
class TProfile2D;

class AliAnalysisTask3HeMC: public AliAnalysisTaskSE {
public:

   //class constructors
   AliAnalysisTask3HeMC();
   AliAnalysisTask3HeMC(const char *name);
   // class destructor
   virtual ~AliAnalysisTask3HeMC();

   // called once at beginning or runtime
   virtual void UserCreateOutputObjects(); // called for each event
   virtual void UserExec(Option_t* option); // called at end of analysis
   virtual void Terminate(Option_t* option);

   // Setter
   void SetAODAnalysis(){kAODanalysis =kTRUE;};
   void SetESDAnalysis(){kAODanalysis =kFALSE;};
   void SetHasMCData(Bool_t UseMC) {kHasMCdata=UseMC;};
   void SetppAnalysis() {kpp = kTRUE; kpPb =kFALSE; kPbPb = kFALSE;};
   void SetpPbAnalysis() {kpp = kFALSE; kpPb =kTRUE; kPbPb = kFALSE;};
   void SetPbPbAnalysis() {kpp = kFALSE; kpPb =kFALSE; kPbPb = kTRUE;};
   void SetEtaCut(Float_t Min, Float_t Max){fEtaMin = Min; fEtaMax = Max;};
   void SetDCACut(Double_t XY, Double_t Z){fMaxDCAxy = XY; fMaxDCAz = Z;};
   void RejectKinkMothers(){fRejectKinkMother=kTRUE;};

   // Getter
   Bool_t IsAODanalysis() const { return kAODanalysis; };
   Bool_t IsESDanalysis() const { return !kAODanalysis; };
   Bool_t HasMCData() const { return kHasMCdata; };
   Bool_t Ispp() const { return kpp; };
   Bool_t IspPb() const { return kpPb; };
   Bool_t IsPbPb() const { return kPbPb; };

private:
   void ProcessAOD();						// Loop over tracks and perform analysis
   void BinLogAxis(TAxis* axis);			// Method for the correct logarithmic binning of an axis
   Double_t GetDCAxy (AliAODTrack *track); // At the moment only for AODs
   Double_t GetDCAz (AliAODTrack *track); // At the moment only for AODs
   Double_t GetDCAxy (AliAODMCParticle *part); // At the moment only for AODs
   Double_t GetDCAz (AliAODMCParticle *part); // At the moment only for AODs
   Double_t CalculateRapidity(Double_t px, Double_t py, Double_t pz, TString type);

   // Variables
   Bool_t kAODanalysis;
   Bool_t kHasMCdata;
   Bool_t kpp;									// Collision system = pp?
   Bool_t kpPb;								// Collision system = pPb?
   Bool_t kPbPb;								// Collision system = PbPb?
   Float_t fEtaMin;							// Minimum pseudorapidity
   Float_t fEtaMax;							// Maximum pseudorapidity
   Double_t fMaxDCAxy;						// Maximum DCA in radial (xy) direction
   Double_t fMaxDCAz;               	// Maximum DCA in z direction
   Bool_t fRejectKinkMother;				// Reject Kink mothers?
   Double_t Multiplicity;					// Multiplicity
   Double_t MultiplicityPercentile;      // Multiplicity percentile
   Double_t kTPCnSigmaCut;					// loose TPC PID cut

   // QA histograms
   TH2D* fHistMultEstZDep;					//! QA histogram for the dependence of the multiplicity value on the z position of the primary vertex
   TH2D* fHistMultValueVsPercientile;	//! Q Ahistogram for the dependence of the multiplicity value and the percentile
   TH2D* fHistTPCdEdxRigidity; 			//! QA histogram for TPC dE/dx
   TProfile2D* fHistTPCsigHe3;			//! To check TPC spline quality He3
   TProfile2D* fHistTPCsigHe4;			//! To check TPC spline quality He4

   // Results
   THnSparse* fHistTPCdEdxSigmaHe3;
   THnSparse* fHistTPCdEdxSigmaHe4;

   // MC Results
   TH2D* fHistPtTrueRecHe3;	         //! MC pT true vs pT reconstructed
   TH2D* fHistPtTrueRecHe4;				//! MC pT true vs pT reconstructed
   THnSparse* fHistTrueHe3;				//! MC primary He3 paritcle
   THnSparse* fHistTrueHe4;				//! MC primary He4 paritcle
   
   // Objects
   AliEventCuts fEventCuts;            // Event cuts
   AliPIDResponse* fPIDResponse;			//! PID response
   AliAnalysisUtils* fUtils;				//
   TList *fOutput;							//! Container for Task Output
   TList *fQAList;							//! Container for QA Output
   AliAODMCHeader *fAODMCHeader;       // ! MC info AOD
   TClonesArray *fAODArrayMCParticles;      // ! MC info particle AOD


   ClassDef(AliAnalysisTask3HeMC , 1);
};
#endif /* AliAnalysisTask3HeMC_h */
