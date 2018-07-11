/* 
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 
 * See cxx source for full Copyright notice 
 * $Id$ 
 */

/*******************************************************************************
* Analysis task for anisotropic flow analysis of data taken by ALICE 				   *
* with different methods for two- and multiparticle correlations     				   *
*																			                                         *
* Author: Cindy Mordasini (cindy.mordasini@cern.ch)		            					   *
*******************************************************************************/

#ifndef ALIANALYSISTASKTWOMULTICORRELATIONS_H
#define ALIANALYSISTASKTWOMULTICORRELATIONS_H

#include "AliAnalysisTaskSE.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliVEvent.h"
#include "TProfile.h"
#include "TComplex.h"
#include "TH1F.h"

//==============================================================================

class AliAnalysisTaskTwoMultiCorrelations : public AliAnalysisTaskSE{
  public:
  // The six following lines are mandatory for AliAnalysisTaskSE to use the class properly
    AliAnalysisTaskTwoMultiCorrelations();
    AliAnalysisTaskTwoMultiCorrelations(const char *name, Bool_t useParticleWeights=kFALSE);
    virtual ~AliAnalysisTaskTwoMultiCorrelations();

    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *);
    virtual void Terminate(Option_t *);


// Organisation of the methods
  // Setters and getters for the data members
  // Methods called in UserCreateOutputObjects()
  // Methods called in UserExec(Option_t *)
  // Methods called in Terminate(Option_t *)

  // Setters and getters for the data members
    void SetParameters(Int_t nCorr, Bool_t usePartweights, Bool_t checkNestedLoops)
    {
      this->fNparticlesCorrelations = nCorr;
      this->fUseParticleWeights = usePartweights;
      this->fDoNestedLoops = checkNestedLoops;
    };  // End of void SetParameters(Int_t, Bool_t, Bool_t)

    void SetHarmonics(Int_t nOne, Int_t nTwo, Int_t nThree, Int_t nFour, Int_t nFive, Int_t nSix, Int_t nSeven, Int_t nEight, Int_t nNine, Int_t nTen, Int_t nEleven, Int_t nTwelve, Int_t nThirteen, Int_t nFourteen)
    {
      this->fHarmonicOne = nOne;
      this->fHarmonicTwo = nTwo;
      this->fHarmonicThree = nThree;
      this->fHarmonicFour = nFour;
      this->fHarmonicFive = nFive;
      this->fHarmonicSix = nSix;
      this->fHarmonicSeven = nSeven;
      this->fHarmonicEight = nEight;
      this->fHarmonicNine = nNine;
      this->fHarmonicTen = nTen;
      this->fHarmonicEleven = nEleven;
      this->fHarmonicTwelve = nTwelve;
      this->fHarmonicThirteen = nThirteen;
      this->fHarmonicFourteen = nFourteen;
    };  // End of void SetHarmonics(Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t)

    void SetCentralityBinning(Int_t const nBins, Float_t minBin, Float_t maxBin)
    {
      this->fMinCentrality = minBin;
      this->fMaxCentrality = maxBin;
    };  // End of void SetCentralityBinning(Int_t const, Float_t, Float_t)

    void SetControlOutputList(TList* const col) {this->fControlOutputList = col;};
    TList* GetControlOutputList() const {return this->fControlOutputList;}
    void SetDraftOutputList(TList* const dol) {this->fDraftOutputList = dol;};
    TList* GetDraftOutputList() const {return this->fDraftOutputList;}
    void SetFinalOutputList(TList* const fol) {this->fFinalOutputList = fol;};
    TList* GetFinalOutputList() const {return this->fFinalOutputList;}

  // Methods called in UserCreateOutputObjects()
    virtual void BookAndNestAllLists();
    virtual void BookControlList();
    virtual void BookDraftList();
    virtual void BookFinalList();

  // Methods called in UserExec(Option_t *)
    TComplex CalculateQvector(Int_t n, Int_t p, Int_t nParticles, Double_t phi[], Double_t particleWeight[]);
    TComplex CalculateRecursionWithQvectors(Int_t nParticles, Double_t phi[], Double_t particleWeight[], Int_t nCorr, Int_t harmonics[], Int_t p = 1, Int_t skip = 0);
    virtual void ComputeCorrelationsWithQvectors(Int_t nParticles, Double_t phi[], Double_t particleWeight[], Int_t harmonics[], Int_t nCorr);
    virtual void ComputeCorrelationsWithTwoNestedLoops(Int_t nOne, Int_t nTwo, Int_t nParticles, Double_t phi[], Double_t particleWeight[]);
    virtual void ComputeCorrelationsWithFourNestedLoops(Int_t nOne, Int_t nTwo, Int_t nThree, Int_t nFour, Int_t nParticles, Double_t phi[], Double_t particleWeight[]);
    virtual void ComputeCorrelationsWithStandAloneQvectors(Int_t n, Int_t p, Int_t nParticles, Double_t phi[], Double_t particleWeight[]);

  // Methods called in Terminate(Option_t *)
    //virtual void EstimateFlowWithCumulants(Int_t maxMcorr);


  private:
    AliAnalysisTaskTwoMultiCorrelations(const AliAnalysisTaskTwoMultiCorrelations& aattmc);
    AliAnalysisTaskTwoMultiCorrelations& operator=(const AliAnalysisTaskTwoMultiCorrelations& aattmc);

// Organisation of the data members
  // Initial general parameters
  // Selection of the methods to compute
  // Structure of the output
  // Histograms and TProfiles

  // Initial general parameters
    Int_t fNparticlesCorrelations;  // Number of m-particle correlations and harmonics (2-14)
    Int_t fHarmonicOne; // Harmonic n_1
    Int_t fHarmonicTwo; // Harmonic n_2
    Int_t fHarmonicThree; // Harmonic n_3
    Int_t fHarmonicFour;  // Harmonic n_4
    Int_t fHarmonicFive;  // Harmonic n_5
    Int_t fHarmonicSix; // Harmonic n_6
    Int_t fHarmonicSeven; // Harmonic n_7
    Int_t fHarmonicEight; // Harmonic n_8
    Int_t fHarmonicNine;  // Harmonic n_9
    Int_t fHarmonicTen; // Harmonic n_10
    Int_t fHarmonicEleven;  // Harmonic n_11
    Int_t fHarmonicTwelve;  // Harmonic n_12
    Int_t fHarmonicThirteen;  // Harmonic n_13
    Int_t fHarmonicFourteen;  // Harmonic n_14
    Float_t fMinCentrality; // Minimum of the centrality
    Float_t fMaxCentrality; // Maximum of the centrality

  // Selection of the methods to compute (unit-weighted Q-vectors are default)
    Bool_t fUseParticleWeights; // Use non-unit particle weights
    Bool_t fDoNestedLoops;  // Cross-check the results with nested loops

  // Structure of the output
    TList* fOutputList; // Main list holding all the output objects
    TList* fControlOutputList;  // List holding all the control objects
    TList* fDraftOutputList;  // List holding all the intermediate objects
    TList* fFinalOutputList;  // List holding all the final results

  // Histograms and TProfiles
    TH1F *fPtControlHisto;  //! Control histogram for the transverse momentum
    TH1F *fEtaControlHisto; //! Control histogram for the pseudorapidity
    TH1F *fControlPhiHisto; //! Control histogram for the azimuthal angles
    TH1F *fCentralityHisto; //! Control histogram for the centrality
    TH1F *fMultiplicityDist;  //! Control histogram for the multiplicity distribution
    TProfile *fCorrelationWithQvectorsProfile;  //! m-p correlation estimated with Q-vectors
    TProfile *fCorrelationWithNestedLoopsProfile; //! m-p correlation estimated with nested loops
    TProfile *fCorrelationWithQvectorsSaProfile;  //! 2-p correlation estimated with stand-alone Q-vectors
    //TProfile *fEstimatedFlowWithQcProfile;  //! Anisotropic flow estimated with Q-cumulants

  // Version counter for the submissions on Grid
    // Increase the counter by one when the latest version changes the structure
    // of the output file
    ClassDef(AliAnalysisTaskTwoMultiCorrelations,3);

}; // End of the definition of the class AliAnalysisTaskTwoMultiCorrelations

#endif

