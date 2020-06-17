/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskUPCforwardMC_H
#define AliAnalysisTaskUPCforwardMC_H

/**
 * @file   AliAnalysisTaskUPCforwardMC.h
 * @author Simone Ragoni <Simone.Ragoni@cern.ch>
 * @date   February 2019
 */

#include "AliAnalysisTaskSE.h"
#include "TBits.h"

class AliMuonTrackCuts; 					// Include class for standard muon tack cuts

/**
 * \file AliAnalysisTaskUPCforwardMC.h
 * \brief Contains the declaration of the AliAnalysisTaskUPCforwardMC class
 */

/**
 * \class AliAnalysisTaskUPCforwardMC
 * \brief Contains the tools to every part of my analysis
 */
class AliAnalysisTaskUPCforwardMC : public AliAnalysisTaskSE
{
    public:
                                /**
                                 * Create a new AliAnalysisTaskUPCforwardMC with
                                 * default constructor. Based on my understanding
                                 * this is mostly related to input/output
                                 * processes.
                                 */
                                AliAnalysisTaskUPCforwardMC();

                                /**
                                 * Create a new AliAnalysisTaskUPCforwardMC with
                                 * a specific name. This serves to instantiate
                                 * the AliAnalysisTaskSE object built with the
                                 * constructor, that will take this specific
                                 * name.
                                 *
                                 * \param name , the name taken by the AliAnalysisTaskSE object.
                                 * \param isMC , 0 if Data, 1 if MC(look at the AddTask please).
                                 */
                                AliAnalysisTaskUPCforwardMC(const char *name);

                                /**
                                 * Virtual destructor for the class. It will
                                 * take care of all the particular objects to be
                                 * destroyed for the class.
                                 */
        virtual                 ~AliAnalysisTaskUPCforwardMC();

                                /**
                                 * The function related to the instantiation of
                                 * all the histograms and the output list.
                                 */
        virtual void            UserCreateOutputObjects();

                                /**
                                 * Everything happens here. Here, the cuts are
                                 * applied, the histograms are filled and the
                                 * J/Psi peak is manifested inside the histograms.
                                 *
                                 * \param option , actually it is not used for now...
                                 */
        virtual void            UserExec(Option_t* option);

                                /**
                                 * Called at the END of the analysis (when all
                                 * events are processed). But it is not actually
                                 * doing anything! I guess it is mostly needed
                                 * for I/O purposes and GRID interfacing...
                                 */
        virtual void            Terminate(Option_t* option);

                                /**
                                 * Implement the NotifyRun to search for the new
                                 * parameters at each new runs. Sets the run
                                 * number for the successive cuts.
                                 */
        virtual void            NotifyRun();

                                /**
                                 * Set the cap for the luminosity.
                                 * For 2018 MC this is 40k*lumi per run...
                                 * For 2015 it is different because the
                                 * amount of kCohJPsiToMu is double the
                                 * other processes!
                                 */
        virtual void            SetLuminosityCap();


                                /**
                                 * This will fill the vector containing the good
                                 * run numbers. For now this function will be
                                 * inside the constructor of the class.
                                 */
        // void                    FillGoodRunVector(std::vector<Int_t> &fVectorGoodRunNumbers);

                                /**
                                 * This function substitutes the roel of the
                                 * triggers for the MC.
                                 */
        Bool_t                  IsTriggered();

                                /**
                                 * This function is needed to record the MC
                                 * truth information inside the nanoAODs.
                                 * Everything MC related should go in here.
                                 * This should help for a subsequent merging
                                 * of the MC and DATA classes.
                                 */
        void                    ProcessMCParticles(AliMCEvent* fMCEvent);

                                /**
                                 * This function computes the cos(theta) in the
                                 * Collins-Soper frame for the helicity
                                 * analysis.
                                 */
        Double_t                CosThetaCollinsSoper( TLorentzVector muonPositive,
                                                      TLorentzVector muonNegative,
                                                      TLorentzVector possibleJPsi );

                                /**
                                 * This function computes the cos(theta) in the
                                 * Helicity frame for the helicity
                                 * analysis.
                                 */
        Double_t                CosThetaHelicityFrame( TLorentzVector muonPositive,
                                                       TLorentzVector muonNegative,
                                                       TLorentzVector possibleJPsi );

                                /**
                                 * This function computes the cos(phi) in the
                                 * Collins-Soper frame for the helicity
                                 * analysis.
                                 */
        Double_t                CosPhiCollinsSoper( TLorentzVector muonPositive,
                                                    TLorentzVector muonNegative,
                                                    TLorentzVector possibleJPsi );

                                /**
                                 * This function computes the cos(phi) in the
                                 * Helicity frame for the helicity
                                 * analysis.
                                 */
        Double_t                CosPhiHelicityFrame( TLorentzVector muonPositive,
                                                     TLorentzVector muonNegative,
                                                     TLorentzVector possibleJPsi );

                                /**
                                 * Use the class as a data member. It contains
                                 * the cuts for the muon track.
                                 */
        AliMuonTrackCuts*       fMuonTrackCuts;


    private:

                                /// The input events for the analysis.
        AliAODEvent*            fAOD;               //!

                                /**
                                 * The output list containing all the histograms
                                 * required for the analysis. In a second time
                                 * I will probably make it so to include every
                                 * possible cut variation to better compute the
                                 * systematics.
                                 */
        TList*                  fOutputList;        //!

                                /**
                                 * The corresponding MC event as explained in
                                 * the Analysis Tutorial...
                                 */
        AliMCEvent*             fMCEvent;       //!

                                /**
                                 * Utility type histo. It counts the GOOD muons
                                 * per event.
                                 */
        TH1F*                   fNumberMuonsH;      //!

                                /**
                                 * In this histogram the number of events passing
                                 * each cut is recorded. So it is highest at the
                                 * 0-th cut, and lowest, even possibly null, at
                                 * the last possible cut.
                                 */
        TH1F*                   fCounterH;          //!

                                /**
                                 * As far as I understand, it should be the
                                 * pseudorapidity distribution of the many muons.
                                 */
        TH1F*                   fEtaMuonH;          //!

                                /**
                                 *
                                 */
        TH1F*                   fRAbsMuonH;         //!

                                /**
                                 * This histogram records the invariant mass
                                 * distribution of the dimuon system.
                                 */
        TH1F*                   fInvariantMassDistributionH;                        //!
        TH1F*                   fInvariantMassDistributionRapidityBinsH[6];         //!

                                /**
                                 * This histogram records the number of entries
                                 * against the respective run number.
                                 */
        TH1F*                   fEntriesAgainstRunNumberH;         //!

                                /**
                                 * This histogram records the number of entries
                                 * against the respective run number, but the
                                 * binning is done properly, like normally
                                 * found for this type of histograms.
                                 */
        TH1F*                   fEntriesAgainstRunNumberProperlyH;         //!

                                /**
                                 * This histogram records the invariant mass
                                 * distribution of the dimuon system, only
                                 * coherent component, so as to say, only
                                 * pt < 0.25 GeV/c for pt of the dimuon pair.
                                 */
        TH1F*                   fInvariantMassDistributionCoherentH;                      //!
        TH1F*                   fInvariantMassDistributionCoherentRapidityBinsH[6];       //!

                                /**
                                 * This histogram records the invariant mass
                                 * distribution of the dimuon system, only
                                 * coherent component, so as to say, only
                                 * pt > 0.25 GeV/c for pt of the dimuon pair.
                                 */
        TH1F*                   fInvariantMassDistributionIncoherentH;                    //!
        TH1F*                   fInvariantMassDistributionIncoherentRapidityBinsH[6];     //!

                                /**
                                 * This histogram records the pt-ditribution
                                 * of the dimuon pairs.
                                 */
        TH1F*                   fDimuonPtDistributionH;         //!

                                /**
                                 * This histogram records the pt-ditribution
                                 * of the dimuon pairs. This is the template
                                 * to be used for the Pt-distribution fit.
                                 *
                                 * Rapidity bin study too.
                                 */
        TH1F*                   fTemplatePtDistributionH;                       //!
        TH1F*                   fTemplatePtDistributionRapidityH[3];            //!
        TH1F*                   fTemplatePtDistributionRapidityMoreH[6];        //!
        TH1F*                   fTemplatePtDistributionHLowerSide;              //!
        TH1F*                   fTemplatePtDistributionRapidityHLowerSide[3];   //!
        TH1F*                   fTemplatePtDistributionHHigherSide;             //!
        TH1F*                   fTemplatePtDistributionRapidityHHigherSide[3];  //!

                                /**
                                 * This histogram records the invariant mass
                                 * distribution of the dimuon system, but the
                                 * novelty is the use of a different method to
                                 * extract the data. While before we were always
                                 * using the fMomentum, this histogram will be
                                 * filled wth the information from the
                                 * fMomentumAtDCA, hopefully this may lead to a
                                 * change in the sigma of the peaks. If not, it
                                 * will just do nothing!
                                 */
        TH2F*                   fDcaAgainstInvariantMassH;         //!


        //_______________________________
        // Cloned histograms with EXTENDED Range (0,20)->(0,40).

                                /**
                                 * This histogram records the invariant mass
                                 * distribution of the dimuon system.
                                 */
        TH1F*                   fInvariantMassDistributionExtendedH;         //!

                                /**
                                 * This histogram records the invariant mass
                                 * distribution of the dimuon system, only
                                 * coherent component, so as to say, only
                                 * pt < 0.25 GeV/c for pt of the dimuon pair.
                                 */
        TH1F*                   fInvariantMassDistributionCoherentExtendedH;       //!

                                /**
                                 * This histogram records the invariant mass
                                 * distribution of the dimuon system, only
                                 * coherent component, so as to say, only
                                 * pt > 0.25 GeV/c for pt of the dimuon pair.
                                 */
        TH1F*                   fInvariantMassDistributionIncoherentExtendedH;     //!

                                /**
                                 * This histogram shows the angular distribution
                                 * of the positive muon in the rest frame of the
                                 * J/Psi. This histogram is needed to evaluate
                                 * the polarization of the J/Psi!
                                 */
        TH1F*                   fAngularDistribOfPositiveMuonRestFrameJPsiH;      //!

                                /**
                                 * This histogram shows the angular distribution
                                 * of the negative muon in the rest frame of the
                                 * J/Psi. This histogram is needed to evaluate
                                 * the polarization of the J/Psi!
                                 */
        TH1F*                   fAngularDistribOfNegativeMuonRestFrameJPsiH;      //!

                                /**
                                 * This histogram represents the check over
                                 * the helicity of the J/Psi. It should be flat
                                 * if I remember well enough!
                                 */
        TH1F*                   fCheckHelicityRestFrameJPsiH;

                                /**
                                 * This histogram shows the angular distribution
                                 * of the positive muon in the rest frame of the
                                 * J/Psi. This histogram is needed to evaluate
                                 * the polarization of the J/Psi! Divided in
                                 * psudorapidity bins... 8??
                                 */
        TH1F*                   fThetaDistribOfPositiveMuonRestFrameJPsiRapidityBinH[8];


        //_______________________________
        // MC TRUTH PLOTS
                                /**
                                 * This histogram records the pdgCodes of all
                                 * the MC GENERATED particles...
                                 */
        TH1F*                   fMCpdgCodesH;           //!

                                /**
                                 * This histogram records the pdgCodes of all
                                 * the MC GENERATED PRIMARY particles...
                                 */
        TH1F*                   fMCpdgCodesOnlyPrimaryH;           //!

        // MC TRUTH GENERATED FEATURES
                                /**
                                 * This histogram records the PHI of all
                                 * the MC GENERATED particles...
                                 */
        TH1F*                   fMCphiGeneratedTruthH;           //!

                                /**
                                 * This histogram records the ETA of all
                                 * the MC GENERATED PRIMARY particles...
                                 */
        TH1F*                   fMCetaGeneratedTruthH;           //!

                                /**
                                 * This histogram records the PSEUDORAPIDITY of
                                 * all the MC GENERATED PRIMARY particles...
                                 */
        TH1F*                   fMCpseudorapidityGeneratedTruthH;           //!

                                /**
                                 * This histogram records the Pt of all
                                 * the MC GENERATED particles...
                                 */
        TH1F*                   fMCptGeneratedTruthH;           //!

        // MC TRUTH GENERATED DIMUONS FEATURES
                                /**
                                 * This histogram records the PHI of all
                                 * the MC GENERATED muons of the J/Psi...
                                 */
        TH1F*                   fMCphiDimuonGeneratedTruthH;           //!

                                /**
                                 * This histogram records the ETA of all
                                 * the MC GENERATED muons of the J/Psi...
                                 */
        TH1F*                   fMCetaDimuonGeneratedTruthH;           //!

                                /**
                                 * This histogram records the PSEUDORAPIDITY of
                                 * all the MC GENERATED muons of the J/Psi...
                                 */
        TH1F*                   fMCpseudorapidityDimuonGeneratedTruthH;           //!

                                /**
                                 * This histogram records the Pt of all
                                 * the MC GENERATED muons of the J/Psi...
                                 */
        TH1F*                   fMCptDimuonGeneratedTruthSingleMuonsH;           //!

                                /**
                                 * This histogram records the Pt of all
                                 * the MC GENERATED J/Psi.
                                 */
        TH1F*                   fMCptDimuonGeneratedTruthH;           //!






                                /**
                                 * Invariant Mass Distribution of the J/Psi
                                 * before EVENT and TRACK selection for the MC
                                 * (aka GENERATED).
                                 */
        TH1F*                   fMCinvariantMassDistrJPsiGeneratedTruthH;           //!

                                /**
                                 * Invariant Mass Distribution of the J/Psi
                                 * after EVENT and TRACK selection for the MC
                                 * WITH MC TRUTH.
                                 */
        TH1F*                   fMCinvariantMassDistrJPsiAfterEvtAndTrkSelectionTruthH;           //!

                                /**
                                 * This histogram shows the angular distribution
                                 * of the positive muon in the rest frame of the
                                 * J/Psi. This histogram is needed to evaluate
                                 * the polarization of the J/Psi!
                                 */
        TH1F*                   fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthH;           //!

                                /**
                                 * This histogram shows the angular distribution
                                 * of the negative muon in the rest frame of the
                                 * J/Psi. This histogram is needed to evaluate
                                 * the polarization of the J/Psi!
                                 */
        TH1F*                   fMCthetaDistribOfNegativeMuonRestFrameJPsiGeneratedTruthH;           //!

                                /**
                                 * This histogram shows the angular distribution
                                 * of the positive muon in the rest frame of the
                                 * J/Psi. This histogram is needed to evaluate
                                 * the polarization of the J/Psi! Divided in
                                 * psudorapidity bins... 8??
                                 */
        TH1F*                   fMCthetaDistribOfPositiveMuonRestFrameJPsiGeneratedTruthRapidityBinH[8];

                                /**
                                 * This is the GENERATED
                                 * value of the cos(theta) in the HELICITY
                                 * frame. It should be
                                 * needed for the "bin migration" study...
                                 * As far as I can imagine, in each event there
                                 * should only be a single J/Psi because these
                                 * are all UPC events, as such I should record
                                 * all possible values for cos(theta) and plot
                                 * only when the corresponding J/Psi falls
                                 * inside the detector acceptance...
                                 * The resulting plot should be TH2F and be like
                                 * a lego plot. See for an example:
                                 *https://www.researchgate.net/figure/The-migration-matrix-for-leading-p-jet-T-Element-i-j-is-the-probability-for-a-particle_fig1_222896619
                                 */
        Double_t                fCosThetaGeneratedHelicityFrame;           //!
        Double_t                fPhiGeneratedHelicityFrame;           //!

                                /**
                                 * This is the RECONSTRUCTED
                                 * value of the cos(theta) in the HELICITY
                                 * frame. It should be
                                 * needed for the "bin migration" study...
                                 * As far as I can imagine, in each event there
                                 * should only be a single J/Psi because these
                                 * are all UPC events, as such I should record
                                 * all possible values for cos(theta) and plot
                                 * only when the corresponding J/Psi falls
                                 * inside the detector acceptance...
                                 * The resulting plot should be TH2F and be like
                                 * a lego plot. See for an example:
                                 *https://www.researchgate.net/figure/The-migration-matrix-for-leading-p-jet-T-Element-i-j-is-the-probability-for-a-particle_fig1_222896619
                                 */
        Double_t                fCosThetaReconHelicityFrame;           //!
        Double_t                fPhiReconHelicityFrame;           //!

                                /**
                                 * Counter for the UPC events to access vectors.
                                 * The only thing that matters is a counter to
                                 * give me the proper value of the vector for
                                 * the generated value I should consider inside
                                 * the vector of the GENERATED level.
                                 * It is important to remember there is at MOST
                                 * a single J/Psi for UPC event!!!
                                 */
        Int_t                   fCounterUPCevent;           //!

                                /**
                                 * This histogram shows the bin migration for
                                 * the helicity Analysis...
                                 * NOTE: this is only in CosTheta.
                                 */
        TH2F*                   fBinMigrationHelicityH;           //!

                                /**
                                 * This histogram shows the bin migration for
                                 * the helicity Analysis...
                                 * NOTE: this is only in Phi.
                                 */
        TH2F*                   fBinMigrationForPhiHelicityH;           //!

        //_______________________________
        // HELICITY AND COLLINS-SOPER

                                /**
                                 * This histogram shows the angular distribution
                                 * of the positive muon in the HELICITY frame.
                                 * COS(THETA) distribution.
                                 */
        TH1F*                   fCosThetaHelicityFrameJPsiH;

                                /**
                                 * This histogram shows the angular distribution
                                 * of the positive muon in the HELICITY frame.
                                 * PHI distribution.
                                 */
        TH1F*                   fPhiHelicityFrameJPsiH;

                                /**
                                 * This histogram shows the angular distribution
                                 * of the positive muon in the COLLINS-SOPER
                                 * frame.  COS(THETA) distribution.
                                 */
        TH1F*                   fCosThetaCollinsSoperFrameJPsiH;

                                /**
                                 * This histogram shows the angular distribution
                                 * of the positive muon in the COLLINS-SOPER
                                 * frame. PHI distribution.
                                 */
        TH1F*                   fPhiCollinsSoperFrameJPsiH;

                                /**
                                 * This histogram shows the angular distribution
                                 * of the positive muon in the HELICITY frame.
                                 * COS(THETA) distribution. Divided per
                                 * rapidity bins. RECONSTRUCTED.
                                 */
        TH1F*                   fCosThetaHelicityFrameJPsiRapidityBinsH[8];

                                /**
                                 * This histogram shows the angular distribution
                                 * of the positive muon in the HELICITY frame.
                                 * PHI distribution. Divided per
                                 * rapidity bins. RECONSTRUCTED.
                                 */
        TH1F*                   fPhiHelicityFrameJPsiRapidityBinsH[8];

                                /**
                                 * This histogram shows the angular distribution
                                 * of the positive muon in the COLLINS-SOPER
                                 * frame.  COS(THETA) distribution. Divided per
                                 * rapidity bins. RECONSTRUCTED.
                                 */
        TH1F*                   fCosThetaCollinsSoperFrameJPsiRapidityBinsH[8];

                                /**
                                 * This histogram shows the angular distribution
                                 * of the positive muon in the COLLINS-SOPER
                                 * frame. PHI distribution. Divided per
                                 * rapidity bins. RECONSTRUCTED.
                                 */
        TH1F*                   fPhiCollinsSoperFrameJPsiRapidityBinsH[8];

                                /**
                                 * This histogram shows the angular distribution
                                 * of the positive muon in the HELICITY frame.
                                 * COS(THETA) distribution. Divided per
                                 * rapidity bins. RECONSTRUCTED. 10 rapidity bins.
                                 */
        TH1F*                   fCosThetaHelicityFrameJPsiTenRapidityBinsH[10];

                                /**
                                 * This histogram shows the angular distribution
                                 * of the positive muon in the HELICITY frame.
                                 * PHI distribution. Divided per
                                 * rapidity bins. RECONSTRUCTED. 10 rapidity bins.
                                 */
        TH1F*                   fPhiHelicityFrameJPsiTenRapidityBinsH[10];

                                /**
                                 * This histogram shows the angular distribution
                                 * of the positive muon in the COLLINS-SOPER
                                 * frame.  COS(THETA) distribution. Divided per
                                 * rapidity bins. RECONSTRUCTED. 10 rapidity bins.
                                 */
        TH1F*                   fCosThetaCollinsSoperFrameJPsiTenRapidityBinsH[10];

                                /**
                                 * This histogram shows the angular distribution
                                 * of the positive muon in the COLLINS-SOPER
                                 * frame. PHI distribution. Divided per
                                 * rapidity bins. RECONSTRUCTED. 10 rapidity bins.
                                 */
        TH1F*                   fPhiCollinsSoperFrameJPsiTenRapidityBinsH[10];

        /*________________________________________
         *
         * FINAL INCARNATION OF THE
         * POLARISATION ANALYSIS
         *
         */
                                /**
                                 * This histogram shows the angular distribution
                                 * of the positive muon in the HELICITY frame.
                                 * COS(THETA) distribution. GENERATED
                                 */
        TH1F*                   fMCCosThetaHelicityFrameJPsiH;

                                /**
                                 * This histogram shows the angular distribution
                                 * of the positive muon in the HELICITY frame.
                                 * PHI distribution. GENERATED
                                 */
        TH1F*                   fMCPhiHelicityFrameJPsiH;

                                /**
                                 * This histogram shows the angular distribution
                                 * of the positive muon in the COLLINS-SOPER
                                 * frame.  COS(THETA) distribution. GENERATED
                                 */
        TH1F*                   fMCCosThetaCollinsSoperFrameJPsiH;

                                /**
                                 * This histogram shows the angular distribution
                                 * of the positive muon in the COLLINS-SOPER
                                 * frame. PHI distribution. GENERATED
                                 */
        TH1F*                   fMCPhiCollinsSoperFrameJPsiH;

                                /**
                                 * TEMPLATES NEEDED for the SIGNAL EXTRACTION
                                 *
                                 * Signal extraction in Phi, CosTheta, and
                                 * TildePhi, with 25 bins only...
                                 *
                                 * NB: He => helicity frame
                                 * NB: Cs => CS       frame
                                 *
                                 * NB: TH1F*** is the 2D ANALYSIS
                                 */
        TH1F*                   fInvariantMassDistributionOnlyPhiHeFrameTwentyfiveBinsH[25];       //!
        TH1F*                   fInvariantMassDistributionOnlyCosThetaHeFrameTwentyfiveBinsH[25];  //!
        TH1F*                   fInvariantMassDistributionOnlyTildePhiHeFrameTwentyfiveBinsH[25];  //!
        TH1F*                   fInvariantMassDistributionOnlyPhiCsFrameTwentyfiveBinsH[25];       //!
        TH1F*                   fInvariantMassDistributionOnlyCosThetaCsFrameTwentyfiveBinsH[25];  //!
        TH1F*                   fInvariantMassDistributionOnlyTildePhiCsFrameTwentyfiveBinsH[25];  //!

                                /**
                                 * This histogram shows the invariant mass
                                 * distribution of the dimuon pairs in terms
                                 * of bins of cos theta of the positive muon
                                 * in the helicity frame of the J/Psi.
                                 *
                                 * My variable binning: [-0.65, -0.35, -0.15,
                                 * -0.05, 0.05, 0.15, 0.35, 0.65] vs
                                 * 3.14*[1, 19/20, 18/20, 17/20, 13/20, 9/20,
                                 * 6/20, 4/20, 2/20, 1/10, 0, negative side].
                                 *
                                 * NOTE: the first is in the helicity frame,
                                 *       the second with Collins-Soper.
                                 */
        TH1F***                 fInvariantMassDistributionForSignalExtractionHelicityFrameMyBinningH;  //!
        TH1F***                 fInvariantMassDistributionForSignalExtractionCsFrameMyBinningH;        //!

        //_____________________________________

                                /**
                                 * This histogram shows the angular distribution
                                 * of the positive muon in the HELICITY frame.
                                 * COS(THETA) distribution. Divided per
                                 * rapidity bins. GENERATED.
                                 */
        TH1F*                   fMCCosThetaHelicityFrameJPsiRapidityBinsH[8];

                                /**
                                 * This histogram shows the angular distribution
                                 * of the positive muon in the HELICITY frame.
                                 * PHI distribution. Divided per
                                 * rapidity bins. GENERATED.
                                 */
        TH1F*                   fMCPhiHelicityFrameJPsiRapidityBinsH[8];

                                /**
                                 * This histogram shows the angular distribution
                                 * of the positive muon in the COLLINS-SOPER
                                 * frame.  COS(THETA) distribution. Divided per
                                 * rapidity bins. GENERATED.
                                 */
        TH1F*                   fMCCosThetaCollinsSoperFrameJPsiRapidityBinsH[8];

                                /**
                                 * This histogram shows the angular distribution
                                 * of the positive muon in the COLLINS-SOPER
                                 * frame. PHI distribution. Divided per
                                 * rapidity bins. GENERATED.
                                 */
        TH1F*                   fMCPhiCollinsSoperFrameJPsiRapidityBinsH[8];

                                /**
                                 * This histogram shows the angular distribution
                                 * of the positive muon in the HELICITY frame.
                                 * COS(THETA) distribution. Divided per
                                 * rapidity bins. GENERATED. 10 rapidity bins.
                                 */
        TH1F*                   fMCCosThetaHelicityFrameJPsiTenRapidityBinsH[10];

                                /**
                                 * This histogram shows the angular distribution
                                 * of the positive muon in the HELICITY frame.
                                 * PHI distribution. Divided per
                                 * rapidity bins. GENERATED. 10 rapidity bins.
                                 */
        TH1F*                   fMCPhiHelicityFrameJPsiTenRapidityBinsH[10];

                                /**
                                 * This histogram shows the angular distribution
                                 * of the positive muon in the COLLINS-SOPER
                                 * frame.  COS(THETA) distribution. Divided per
                                 * rapidity bins. GENERATED. 10 rapidity bins.
                                 */
        TH1F*                   fMCCosThetaCollinsSoperFrameJPsiTenRapidityBinsH[10];

                                /**
                                 * This histogram shows the angular distribution
                                 * of the positive muon in the COLLINS-SOPER
                                 * frame. PHI distribution. Divided per
                                 * rapidity bins. GENERATED. 10 rapidity bins.
                                 */
        TH1F*                   fMCPhiCollinsSoperFrameJPsiTenRapidityBinsH[10];

                                /**
                                 * This histogram shows the angular distribution
                                 * of the positive muon in the HELICITY
                                 * frame.  COS(THETA) distribution. RECON.
                                 */
        TH1F*                   fCosThetaHeFrameForSignalExH;

                                /**
                                 * This histogram shows the angular distribution
                                 * of the positive muon in the HELICITY
                                 * frame.  PHI distribution. RECON.
                                 */
        TH1F*                   fPhiHeFrameForSignalExH;

                                /**
                                 * This histogram shows the angular distribution
                                 * of the positive muon in the HELICITY
                                 * frame.  COS(THETA) distribution. GEN.
                                 */
        TH1F*                   fMCCosThetaHeFrameForSignalExH;

                                /**
                                 * This histogram shows the angular distribution
                                 * of the positive muon in the HELICITY
                                 * frame.  PHI distribution. GEN.
                                 */
        TH1F*                   fMCPhiHeFrameForSignalExH;

                                /**
                                 * This histogram shows the invariant mass
                                 * distribution of the dimuon pairs in terms
                                 * of bins of cos theta of the positive muon
                                 * in the helicity frame of the J/Psi.
                                 *
                                 * What it means is that we divide in 5 bins of
                                 * possible CosTheta of the decaying J/Psi,
                                 * meaning  (-1,-0.8), (-0.8,-0.6), (-0.6,-0.4),
                                 * (-0.4,-0.2) and so on until (0.8,1). We fill
                                 * the invariant mass distribution of the
                                 * dimuons in this many bins.
                                 *
                                 * The next step is to fit this invariant mass
                                 * distributions, so as to obtain the relative
                                 * contribution of J/Psi and GammaGamma to the
                                 * angular distributions. This should help in
                                 * validating our results...
                                 *
                                 * NEW: This histogram how this thing in 2D.
                                 * IT is shown the same distribution in terms of
                                 * CosTheta and Phi bins. Let's see the results!
                                 *
                                 * IMPORTANT: RECONSTRUCTED
                                 */
        TH2F*                   fInvariantMassDistributionBinsOfCosThetaAndPhiHelicityFrameH;  //!

                                /**
                                 * This histogram shows the invariant mass
                                 * distribution of the dimuon pairs in terms
                                 * of bins of cos theta of the positive muon
                                 * in the helicity frame of the J/Psi.
                                 *
                                 * What it means is that we divide in 5 bins of
                                 * possible CosTheta of the decaying J/Psi,
                                 * meaning  (-1,-0.8), (-0.8,-0.6), (-0.6,-0.4),
                                 * (-0.4,-0.2) and so on until (0.8,1). We fill
                                 * the invariant mass distribution of the
                                 * dimuons in this many bins.
                                 *
                                 * The next step is to fit this invariant mass
                                 * distributions, so as to obtain the relative
                                 * contribution of J/Psi and GammaGamma to the
                                 * angular distributions. This should help in
                                 * validating our results...
                                 *
                                 * NEW: This histogram how this thing in 2D.
                                 * IT is shown the same distribution in terms of
                                 * CosTheta and Phi bins. Let's see the results!
                                 *
                                 * IMPORTANT: GENERATED
                                 *
                                 * This is useful because we can divide the
                                 * Reconstructed for the Generated to obtain
                                 * a generalised ACCxEFF!
                                 */
        TH2F*                   fMCInvariantMassDistributionBinsOfCosThetaAndPhiHelicityFrameH;  //!

                                /**
                                 * This histogram shows the invariant mass
                                 * distribution of the dimuon pairs in terms
                                 * of bins of cos theta of the positive muon
                                 * in the helicity frame of the J/Psi.
                                 *
                                 * What it means is that we divide in 5 bins of
                                 * possible CosTheta of the decaying J/Psi,
                                 * meaning  (-1,-0.8), (-0.8,-0.6), (-0.6,-0.4),
                                 * (-0.4,-0.2) and so on until (0.8,1). We fill
                                 * the invariant mass distribution of the
                                 * dimuons in this many bins.
                                 *
                                 * The next step is to fit this invariant mass
                                 * distributions, so as to obtain the relative
                                 * contribution of J/Psi and GammaGamma to the
                                 * angular distributions. This should help in
                                 * validating our results...
                                 *
                                 * NEW: This histogram shows this thing in 2D.
                                 * IT is shown the same distribution in terms of
                                 * CosTheta and Phi bins. Let's see the results!
                                 *
                                 * NEW: this specific histogram follows the
                                 * binning of the inclusive people...
                                 * Hopefully this would yield better results!
                                 *
                                 * RECONSTRUCTED level.
                                 */
        TH2F*                   fCosThetaAndPhiHelicityFrameInclusivePeopleBinningH;  //!

                                /**
                                 * This histogram shows the invariant mass
                                 * distribution of the dimuon pairs in terms
                                 * of bins of cos theta of the positive muon
                                 * in the helicity frame of the J/Psi.
                                 *
                                 * What it means is that we divide in 5 bins of
                                 * possible CosTheta of the decaying J/Psi,
                                 * meaning  (-1,-0.8), (-0.8,-0.6), (-0.6,-0.4),
                                 * (-0.4,-0.2) and so on until (0.8,1). We fill
                                 * the invariant mass distribution of the
                                 * dimuons in this many bins.
                                 *
                                 * The next step is to fit this invariant mass
                                 * distributions, so as to obtain the relative
                                 * contribution of J/Psi and GammaGamma to the
                                 * angular distributions. This should help in
                                 * validating our results...
                                 *
                                 * NEW: This histogram shows this thing in 2D.
                                 * IT is shown the same distribution in terms of
                                 * CosTheta and Phi bins. Let's see the results!
                                 *
                                 * NEW: this specific histogram follows the
                                 * binning of the inclusive people...
                                 * Hopefully this would yield better results!
                                 *
                                 * GENERATED level.
                                 */
        TH2F*                   fMCCosThetaAndPhiHelicityFrameInclusivePeopleBinningH;  //!

                                /**
                                 * This histogram shows the CosTheta vs Phi
                                 * distribution for the J/Psi.
                                 * The binning here reflects the binning from
                                 * the signal extraction, meaning that this time
                                 * we are operating in 10x10 bins, from
                                 * -1.00 < CosTheta < 1.00 and
                                 * -3.14 < Phi      < 3.14 .
                                 * The name is the same just to remember its
                                 * purpose...
                                 * RECON level.
                                 */
        TH2F*                   fInvariantMassDistributionForSignalExtractionHelicityFrameH;  //!

                                /**
                                 * This histogram shows the CosTheta vs Phi
                                 * distribution for the J/Psi.
                                 * The binning here reflects the binning from
                                 * the signal extraction, meaning that this time
                                 * we are operating in 10x10 bins, from
                                 * -1.00 < CosTheta < 1.00 and
                                 * -3.14 < Phi      < 3.14 .
                                 * The name is the same just to remember its
                                 * purpose...
                                 * GENERATED level.
                                 */
        TH2F*                   fMCInvariantMassDistributionForSignalExtractionHelicityFrameH;  //!

                                /**
                                 * This histogram shows the invariant mass
                                 * distribution of the dimuon pairs in terms
                                 * of bins of cos theta of the positive muon
                                 * in the helicity frame of the J/Psi.
                                 *
                                 * What it means is that we divide in 40 bins of
                                 * possible CosTheta of the decaying J/Psi,
                                 * meaning  (-1,-0.8), (-0.8,-0.6), (-0.6,-0.4),
                                 * (-0.4,-0.2) and so on until (0.8,1). We fill
                                 * the invariant mass distribution of the
                                 * dimuons in this many bins.
                                 *
                                 * The next step is to fit this invariant mass
                                 * distributions, so as to obtain the relative
                                 * contribution of J/Psi and GammaGamma to the
                                 * angular distributions. This should help in
                                 * validating our results...
                                 *
                                 * NEW: the mass range has been extended.
                                 * Meaning there are no bounds on the invariant
                                 * mass range BEFORE signal extraction.
                                 */
        TH1F*                   fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameH[40];  //!

                                /**
                                 * This histogram shows the invariant mass
                                 * distribution of the dimuon pairs in terms
                                 * of bins of cos theta of the positive muon
                                 * in the helicity frame of the J/Psi.
                                 *
                                 * What it means is that we divide in 40 bins of
                                 * possible CosTheta of the decaying J/Psi,
                                 * meaning  (-1,-0.8), (-0.8,-0.6), (-0.6,-0.4),
                                 * (-0.4,-0.2) and so on until (0.8,1). We fill
                                 * the invariant mass distribution of the
                                 * dimuons in this many bins.
                                 *
                                 * The next step is to fit this invariant mass
                                 * distributions, so as to obtain the relative
                                 * contribution of J/Psi and GammaGamma to the
                                 * angular distributions. This should help in
                                 * validating our results...
                                 *
                                 * NEW: the mass range has been extended.
                                 * Meaning there are no bounds on the invariant
                                 * mass range BEFORE signal extraction.
                                 */
        TH1F*                   fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameH[50];  //!

        //_______________________________
        // Efficiency plots.
                                /**
                                 * This histogram shows the entries distribution
                                 * per run of the RECONSTRUCTED level.
                                 * This has to be divided by the corresponding
                                 * fMCEfficiencyPerRunH to extract the
                                 * efficiency on a run-by-run basis.
                                 */
        TH1F*                   fEfficiencyPerRunH;             //!
        TH1F*                   fEfficiencyPerRunRapidityH[6];  //!

                                /**
                                 * This histogram shows the entries distribution
                                 * per run of the GENERATED level.
                                 * I believe that the errors in each bin
                                 * have to be set to 0...
                                 * If this doesn't appear on the histograms,
                                 * it is because this doesn't apply anymore
                                 * and I forgot to eliminate  the comment
                                 * from here....
                                 */
        TH1F*                   fMCEfficiencyPerRunH;             //!
        TH1F*                   fMCEfficiencyPerRunRapidityH[6];  //!

                                /**
                                 * This array of histograms shows the
                                 * distribution of the dead zones of the
                                 * detector on a run-by-run basis.
                                 * To avoid a complicated logic, I have
                                 * declared 60k histograms.
                                 * This is because there are runs from about
                                 * 240k to 300k.
                                 * However, this is too heavy for ROOT.
                                 * Hence, I have declared them as
                                 * SPARSE. Most of these 60k histograms
                                 * are empty anyway, so it is best to save
                                 * on space...
                                 */
        // THnSparseF*             fDeadZoneEtaVsPhiPerRunH[60000 + 1]; //!
        // THnSparseF*             fDeadZoneEtaVsPhiPerRunH[364]; //!    number of total runs
        TH2F*                   fDeadZoneEtaVsPhiPerRunH[364]; //!    number of total runs

                                /**
                                 * This histogram shows CosTheta and Phi
                                 * distribution with my variable binning.
                                 * RECON level.
                                 *
                                 * NB: both HELICITY and COLLINS-SOPER
                                 */
        TH2F*                   fCosThetaAndPhiHelicityFrameMyBinningH;              //!
        TH2F*                   fCosThetaAndPhiCsFrameMyBinningH;                    //!
        TH2F*                   fCosThetaAndPhiHelicityFrameMyBinningReweightingH;   //!
        TH2F*                   fCosThetaAndPhiCsFrameMyBinningReweightingH;         //!
        TH2F*                   fCosThetaAndPhiHelicityFrameMyBinningTriggerH[7];    //!
        TH2F*                   fCosThetaAndPhiCsFrameMyBinningTriggerH[7];          //!


                                /**
                                 * This histogram shows CosTheta and Phi
                                 * distribution with my variable binning.
                                 * GENERATED level.
                                 *
                                 * NB: both HELICITY and COLLINS-SOPER
                                 */
        TH2F*                   fMCCosThetaAndPhiHelicityFrameMyBinningH;              //!
        TH2F*                   fMCCosThetaAndPhiCsFrameMyBinningH;                    //!
        TH2F*                   fMCCosThetaAndPhiHelicityFrameMyBinningReweightingH;   //!
        TH2F*                   fMCCosThetaAndPhiCsFrameMyBinningReweightingH;         //!
        TH2F*                   fMCCosThetaAndPhiHelicityFrameReweightingH;            //!
        TH2F*                   fMCCosThetaAndPhiCsFrameReweightingH;                  //!
        TH2F*                   fMCCosThetaAndPhiHelicityFrameMyBinningTriggerH[7];    //!
        TH2F*                   fMCCosThetaAndPhiCsFrameMyBinningTriggerH[7];          //!

        /* - 1D analysis.
         * - My variable binning.
         */
                                /**
                                 * This histogram shows CosTheta
                                 * distribution with my variable binning.
                                 * RECON level.
                                 */
        TH1F*                   fCosThetaHelicityFrameMyBinningH;   //!

                                /**
                                 * This histogram shows CosTheta
                                 * distribution with my variable binning.
                                 * GENERATED level.
                                 */
        TH1F*                   fMCCosThetaHelicityFrameMyBinningH;   //!

                                /**
                                 * This histogram shows CosTheta
                                 * distribution with my variable 17 binning.
                                 * RECON level.
                                 */
        TH1F*                   fCosThetaHelicityFrameMySeventeenBinningH;   //!

                                /**
                                 * This histogram shows CosTheta
                                 * distribution with my variable 17 binning.
                                 * GENERATED level.
                                 */
        TH1F*                   fMCCosThetaHelicityFrameMySeventeenBinningH;   //!

                                /**
                                 * This histogram shows CosTheta
                                 * distribution with my variable SMALL binning.
                                 * RECON level.
                                 */
        TH1F*                   fCosThetaHelicityFrameMyBinningSmallH;   //!

                                /**
                                 * This histogram shows CosTheta
                                 * distribution with my variable SMALL binning.
                                 * GENERATED level.
                                 */
        TH1F*                   fMCCosThetaHelicityFrameMyBinningSmallH;   //!

                                /**
                                 * This histogram shows CosTheta
                                 * distribution with 25 bins.
                                 * This implies a purity well above 80%.
                                 * RECON level.
                                 *
                                 * NB: CS => Collins-Soper Frame
                                 */
        TH1F*                   fCosThetaHelicityFrameTwentyfiveBinsH;             //!
        TH1F*                   fCosThetaCsFrameTwentyfiveBinsH;                   //!
        TH1F*                   fCosThetaHelicityFrameTwentyfiveBinsTriggerH[7];   //!
        TH1F*                   fCosThetaCsFrameTwentyfiveBinsTriggerH[7];         //!
        TH1F*                   fCosThetaHelicityFrameTwentyfiveBinsSidebandsH;    //!
        TH1F*                   fCosThetaCsFrameTwentyfiveBinsSidebandsH;          //!

                                /**
                                 * This histogram shows CosTheta
                                 * distribution with 25 bins.
                                 * This implies a purity well above 80%.
                                 * GENERATED level..
                                 *
                                 * NB: CS => Collins-Soper Frame
                                 */
        TH1F*                   fMCCosThetaHelicityFrameTwentyfiveBinsH;             //!
        TH1F*                   fMCCosThetaCsFrameTwentyfiveBinsH;                   //!
        TH1F*                   fMCCosThetaHelicityFrameTwentyfiveBinsTriggerH[7];   //!
        TH1F*                   fMCCosThetaCsFrameTwentyfiveBinsTriggerH[7];         //!
        TH2F*                   fMCCosThetaHeVsCsH;                                  //!
        TH2F*                   fMCCosThetaHeVsCsFlatH;                              //!
        TH1F*                   fMCCosThetaHelicityFrameTwentyfiveBinsSidebandsH;    //!
        TH1F*                   fMCCosThetaCsFrameTwentyfiveBinsSidebandsH;          //!

                                /**
                                 * This histogram shows Phi
                                 * distribution with 25 bins.
                                 * This implies a purity well above 80% (?).
                                 * RECON level..
                                 *
                                 * NB: CS => Collins-Soper Frame
                                 */
        TH1F*                   fPhiHelicityFrameTwentyfiveBinsH;             //!
        TH1F*                   fPhiCsFrameTwentyfiveBinsH;                   //!
        TH1F*                   fPhiHelicityFrameTwentyfiveBinsTriggerH[7];   //!
        TH1F*                   fPhiCsFrameTwentyfiveBinsTriggerH[7];         //!
        TH1F*                   fPhiHelicityFrameTwentyfiveBinsSidebandsH;    //!
        TH1F*                   fPhiCsFrameTwentyfiveBinsSidebandsH;          //!

                                /**
                                 * This histogram shows Phi
                                 * distribution with 25 bins.
                                 * This implies a purity well above 80% (?).
                                 * GENERATED level..
                                 *
                                 * NB: CS => Collins-Soper Frame
                                 */
        TH1F*                   fMCPhiHelicityFrameTwentyfiveBinsH;             //!
        TH1F*                   fMCPhiCsFrameTwentyfiveBinsH;                   //!
        TH1F*                   fMCPhiHelicityFrameTwentyfiveBinsTriggerH[7];   //!
        TH1F*                   fMCPhiCsFrameTwentyfiveBinsTriggerH[7];         //!
        TH1F*                   fMCPhiHelicityFrameTwentyfiveBinsSidebandsH;    //!
        TH1F*                   fMCPhiCsFrameTwentyfiveBinsSidebandsH;          //!

                                /**
                                 * This histogram shows the TildePhi
                                 * distribution with 25 bins.
                                 * This implies a purity well above 80% (?).
                                 * RECON level..
                                 *
                                 * NB: CS => Collins-Soper Frame
                                 */
        TH1F*                   fTildePhiHelicityFrameTwentyfiveBinsH;             //!
        TH1F*                   fTildePhiCsFrameTwentyfiveBinsH;                   //!
        TH1F*                   fTildePhiHelicityFrameTwentyfiveBinsTriggerH[7];   //!
        TH1F*                   fTildePhiCsFrameTwentyfiveBinsTriggerH[7];         //!
        TH1F*                   fTildePhiHelicityFrameTwentyfiveBinsSidebandsH;    //!
        TH1F*                   fTildePhiCsFrameTwentyfiveBinsSidebandsH;          //!

                                /**
                                 * This histogram shows the TildePhi
                                 * distribution with 25 bins.
                                 * This implies a purity well above 80% (?).
                                 * GENERATED level..
                                 *
                                 * NB: CS => Collins-Soper Frame
                                 */
        TH1F*                   fMCTildePhiHelicityFrameTwentyfiveBinsH;             //!
        TH1F*                   fMCTildePhiCsFrameTwentyfiveBinsH;                   //!
        TH1F*                   fMCTildePhiHelicityFrameTwentyfiveBinsTriggerH[7];   //!
        TH1F*                   fMCTildePhiCsFrameTwentyfiveBinsTriggerH[7];         //!
        TH1F*                   fMCTildePhiHelicityFrameTwentyfiveBinsSidebandsH;    //!
        TH1F*                   fMCTildePhiCsFrameTwentyfiveBinsSidebandsH;          //!

                                /**
                                 * This histogram shows  Phi
                                 * distribution with my variable binning.
                                 * RECON level.
                                 */
        TH1F*                   fPhiHelicityFrameMyBinningH;   //!

                                /**
                                 * This histogram shows  Phi
                                 * distribution with my variable binning.
                                 * GENERATED level.
                                 */
        TH1F*                   fMCPhiHelicityFrameMyBinningH;   //!

                                /**
                                 * This histogram shows the invariant mass
                                 * distribution of the dimuon pairs in terms
                                 * of bins of cos theta of the positive muon
                                 * in the helicity frame of the J/Psi.
                                 * My variable binning.
                                 */
        TH1F*                   fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameMyVariableBinningH[26];  //!

                                /**
                                 * This histogram shows the invariant mass
                                 * distribution of the dimuon pairs in terms
                                 * of bins of cos theta of the positive muon
                                 * in the helicity frame of the J/Psi.
                                 * My variable binning.
                                 */
        TH1F*                   fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameMyVariableBinningH[30];  //!


        //_______________________________
        // CUTS
        /*
         * The following is all the possible checks for the event selections
         * and the track selection as well. Enjoy.
         */
        Int_t                   fRunNum;        //!
        Int_t                   fTracklets;     //!
        Double_t                fLumiPerRun;    //!
        Double_t*               fEtaAndPhi;     //!

        UInt_t                  fL0inputs;      //!
      	UInt_t                  fL1inputs;      //!

      	Double_t                fZem1Energy;    //!
      	Double_t                fZem2Energy;    //!

      	Double_t                fZNCEnergy;     //!
      	Double_t                fZNAEnergy;     //!
      	Double_t                fZPCEnergy;     //!
      	Double_t                fZPAEnergy;     //!
      	Double_t                fZNATDC[4];     //!
      	Double_t                fZNCTDC[4];     //!
      	Double_t                fZPATDC[4];     //!
      	Double_t                fZPCTDC[4];     //!
      	Double_t                fZNATime;       //!
      	Double_t                fZNCTime;       //!
      	Int_t                   fV0ADecision;   //!
      	Int_t                   fV0CDecision;   //!
      	Int_t                   fADADecision;   //!
      	Int_t                   fADCDecision;   //!
        TBits                   fIR1Map;        //!
        TBits                   fIR2Map;        //!


        Bool_t                  fV0Hits[64];    //!
        Int_t                   fV0TotalNCells; //!
        //_______________________________


        //_______________________________
        // TRIGGER INPUTS for MC

        // V0 inputs
        Bool_t                  fBBFlag[64];    //!
        Bool_t                  fBGFlag[64];    //!
        UInt_t                  fBBAFlags;      //!
        UInt_t                  fBBCFlags;      //!
        UInt_t                  fBGAFlags;      //!
        UInt_t                  fBGCFlags;      //!

        // AD inputs
        Bool_t                  fBBFlagAD[16];  //!
        Bool_t                  fBGFlagAD[16];  //!
        UInt_t                  fBBAFlagsAD;    //!
        UInt_t                  fBBCFlagsAD;    //!
        UInt_t                  fBGAFlagsAD;    //!
        UInt_t                  fBGCFlagsAD;    //!
        Double_t                fCounterGeneratedLevel[60000];    //!

        // FINISHED TRIGGER INPUTS for MC
        //_______________________________


        /**
         * This is the vector containing the GOOD RunNumbers.
         */
        // std::vector<Int_t> fVectorGoodRunNumbers;

        /**
         * Not implemented yet...
         */
        AliAnalysisTaskUPCforwardMC(const AliAnalysisTaskUPCforwardMC&);

        /**
         * Not implemented yet...
         */
        AliAnalysisTaskUPCforwardMC& operator=(const AliAnalysisTaskUPCforwardMC&);

        /**
         * This is important for ROOT only. I do not remember the reason anymore.
         * If I happen to encounter it again in the future, I will make sure to
         * record it!
         */
        ClassDef(AliAnalysisTaskUPCforwardMC, 39);
};

#endif
