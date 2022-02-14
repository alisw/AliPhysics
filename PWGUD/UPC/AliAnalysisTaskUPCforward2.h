/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskUPCforward2_H
#define AliAnalysisTaskUPCforward2_H

/**
 * @file   AliAnalysisTaskUPCforward2.h
 * @author Simone Ragoni <Simone.Ragoni@cern.ch>
 * @date   February 2019
 */

#include "AliAnalysisTaskSE.h"
#include "TBits.h"

class AliMuonTrackCuts; 											// Include class for standard muon tack cuts

/**
 * \class MatrixTH1F
 * \brief Alias for a 2D vector of TH1F, roughly speaking an arary of TH1F.
 */
typedef std::vector< std::vector< TH1F* > >  MatrixTH1F;

/**
 * \file AliAnalysisTaskUPCforward2.h
 * \brief Contains the declaration of the AliAnalysisTaskUPCforward2 class
 */

/**
 * \class AliAnalysisTaskUPCforward2
 * \brief Contains the tools to every part of my analysis
 */
class AliAnalysisTaskUPCforward2 : public AliAnalysisTaskSE
{
    public:
                                /**
                                 * Create a new AliAnalysisTaskUPCforward2 with
                                 * default constructor. Based on my understanding
                                 * this is mostly related to input/output
                                 * processes.
                                 */
                                AliAnalysisTaskUPCforward2();

                                /**
                                 * Create a new AliAnalysisTaskUPCforward2 with
                                 * a specific name. This serves to instantiate
                                 * the AliAnalysisTaskSE object built with the
                                 * constructor, that will take this specific
                                 * name.
                                 *
                                 * \param name , the name taken by the AliAnalysisTaskSE object.
                                 */
                                AliAnalysisTaskUPCforward2(const char *name);

                                /**
                                 * Virtual destructor for the class. It will
                                 * take care of all the particular objects to be
                                 * destroyed for the class.
                                 */
        virtual                 ~AliAnalysisTaskUPCforward2();

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
                                 * This will fill the vector containing the good
                                 * run numbers. For now this function will be
                                 * inside the constructor of the class.
                                 */
        void                    FillGoodRunVector(std::vector<Int_t> &fVectorGoodRunNumbers);

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


    protected:

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
        TH1F*                   fInvariantMassDistributionH;         //!

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
        TH1F*                   fInvariantMassDistributionCoherentH;       //!

                                /**
                                 * This histogram records the invariant mass
                                 * distribution of the dimuon system, only
                                 * coherent component, so as to say, only
                                 * pt > 0.25 GeV/c for pt of the dimuon pair.
                                 */
        TH1F*                   fInvariantMassDistributionIncoherentH;     //!

                                /**
                                 * This histogram records the pt-ditribution
                                 * of the dimuon pairs.
                                 */
        TH1F*                   fDimuonPtDistributionH;         //!

                                /**
                                 * This histogram records the energy distri-
                                 * bution of the neutron ZDC. This plot should
                                 * show the relative components of the 0 neutron
                                 * peak, the 1 neutron peak and possibly the
                                 * 2 neutrons peak. Anything higher than that,
                                 * requires help from the user and is more like
                                 * a guess...
                                 */
        TH1F*                   fZNCEnergyAgainstEntriesH;         //!

                                /**
                                 * This histogram records the energy distri-
                                 * bution of the neutron ZDC. This plot should
                                 * show the relative components of the 0 neutron
                                 * peak, the 1 neutron peak and possibly the
                                 * 2 neutrons peak. Anything higher than that,
                                 * requires help from the user and is more like
                                 * a guess...
                                 */
        TH1F*                   fZNAEnergyAgainstEntriesH;         //!

                                /**
                                 * This histogram records the energy distri-
                                 * bution of the neutron ZDC. BEFORE timing
                                 * selection.
                                 */
        TH1F*                   fZNCEnergyBeforeTimingSelectionH;         //!

                                /**
                                 * This histogram records the energy distri-
                                 * bution of the neutron ZDC. BEFORE timing
                                 * selection.
                                 */
        TH1F*                   fZNAEnergyBeforeTimingSelectionH;         //!

                                /**
                                 * This histogram records the energy distri-
                                 * bution of the neutron ZDC. CALIBRATED.
                                 */
        TH1F*                   fZNCEnergyCalibratedH;         //!

                                /**
                                 * This histogram records the energy distri-
                                 * bution of the neutron ZDC. CALIBRATED.
                                 */
        TH1F*                   fZNAEnergyCalibratedH;         //!

                                /**
                                 * This histogram records the energy distri-
                                 * bution of the neutron ZDC. UNCALIBRATED.
                                 */
        TH1F*                   fZNCEnergyUncalibratedH;         //!

                                /**
                                 * This histogram records the energy distri-
                                 * bution of the neutron ZDC. UNCALIBRATED.
                                 */
        TH1F*                   fZNAEnergyUncalibratedH;         //!

                                /**
                                 * This histogram records the energy distri-
                                 * bution of the neutron ZDC. CALIBRATED.
                                 * This is only a trial version for a possible
                                 * future unified plot with the supposedly
                                 * not calibrated runs. What happens here is
                                 * that this plots are filled with the value
                                 * obtained from the LOW RESOLUTION (HIGH GAIN)
                                 * getter. This multiplies for 8 times the
                                 * value obtained fro mnormal getters, and it
                                 * should provide a better description of the
                                 * lower part of the ZDC energy spectrum.
                                 */
        TH1F*                   fZNCEnergyCalibratedHigherGainH;         //!

                                /**
                                 * This histogram records the energy distri-
                                 * bution of the neutron ZDC. CALIBRATED.
                                 * This is only a trial version for a possible
                                 * future unified plot with the supposedly
                                 * not calibrated runs. What happens here is
                                 * that this plots are filled with the value
                                 * obtained from the LOW RESOLUTION (HIGH GAIN)
                                 * getter. This multiplies for 8 times the
                                 * value obtained fro mnormal getters, and it
                                 * should provide a better description of the
                                 * lower part of the ZDC energy spectrum.
                                 */
        TH1F*                   fZNAEnergyCalibratedHigherGainH;         //!







                                /**
                                 * This histogram records the invariant mass
                                 * distribution pt-integrated
                                 * of the dimuon pairs when the
                                 * neutron ZNC has not seen any neutrons.
                                 */
        TH1F*                   fInvariantMassDistributionNoNeutronsH;//!

                                /**
                                 * This histogram records the invariant mass
                                 * distribution pt-integrated
                                 * of the dimuon pairs when the
                                 * neutron ZNC has seen a SINGLE neutron.
                                 */
        TH1F*                   fInvariantMassDistributionOneNeutronH;//!

                                /**
                                 * This histogram records the invariant mass
                                 * distribution pt-integrated
                                 * of the dimuon pairs when the
                                 * neutron ZNC has seen at LEAST a single
                                 * neutron.
                                 */
        TH1F*                   fInvariantMassDistributionAtLeastOneNeutronH;//!


                                /**
                                 * This histogram records the invariant mass
                                 * distribution for the COHERENT component
                                 * of the dimuon pairs when the
                                 * neutron ZNC has not seen any neutrons.
                                 */
        TH1F*                   fInvariantMassDistributionCoherentNoNeutronsH;//!

                                /**
                                 * This histogram records the invariant mass
                                 * distribution for the COHERENT component
                                 * of the dimuon pairs when the
                                 * neutron ZNC has seen a SINGLE neutron.
                                 */
        TH1F*                   fInvariantMassDistributionCoherentOneNeutronH;//!

                                /**
                                 * This histogram records the invariant mass
                                 * distribution for the COHERENT component
                                 * of the dimuon pairs when the
                                 * neutron ZNC has seen at LEAST a single
                                 * neutron.
                                 */
        TH1F*                   fInvariantMassDistributionCoherentAtLeastOneNeutronH;//!

                                /**
                                 * This histogram records the invariant mass
                                 * distribution for the INCOHERENT component
                                 * of the dimuon pairs when the
                                 * neutron ZNC has not seen any neutrons.
                                 */
        TH1F*                   fInvariantMassDistributionIncoherentNoNeutronsH;//!

                                /**
                                 * This histogram records the invariant mass
                                 * distribution for the INCOHERENT component
                                 * of the dimuon pairs when the
                                 * neutron ZNC has seen a SINGLE neutron.
                                 */
        TH1F*                   fInvariantMassDistributionIncoherentOneNeutronH;//!

                                /**
                                 * This histogram records the invariant mass
                                 * distribution for the INCOHERENT component
                                 * of the dimuon pairs when the
                                 * neutron ZNC has seen at LEAST a single
                                 * neutron.
                                 */
        TH1F*                   fInvariantMassDistributionIncoherentAtLeastOneNeutronH;//!

                                /**
                                 * This histogram records the time distribution
                                 * of neutron ZDC.
                                 */
        TH1F*                   fZNCTimeAgainstEntriesH;         //!

                                /**
                                 * This histogram records the time distribution
                                 * of neutron ZDC.
                                 */
        TH1F*                   fZNATimeAgainstEntriesH;         //!

                                /**
                                 * This histogram records the time distribution
                                 * of neutron ZDC. The timing window here is
                                 * a bit stricter at -1.<t<1.
                                 */
        TH1F*                   fZNCTimeStrictTimeWindowH;         //!

                                /**
                                 * This histogram records the time distribution
                                 * of neutron ZDC.The timing window here is
                                 * a bit stricter at -1.<t<1.
                                 */
        TH1F*                   fZNATimeStrictTimeWindowH;         //!

                                /**
                                 * This histogram records the time distribution
                                 * of neutron ZDC. No timing selection to check
                                 * if bunch crossings influence the output.
                                 */
        TH1F*                   fZNCTimeWithoutTimingH[4];         //!

                                /**
                                 * This histogram records the time distribution
                                 * of neutron ZDC. No timing selection to check
                                 * if bunch crossings influence the output.
                                 */
        TH1F*                   fZNATimeWithoutTimingH[4];         //!

                                /**
                                 * This histogram records the time distribution
                                 * of neutron ZDC. I am filling it with all the
                                 * possible TDC info. Basically it is the sum
                                 * of the all the 4 histograms.
                                 */
        TH1F*                   fZNCTime4FillingH;         //!

                                /**
                                 * This histogram records the time distribution
                                 * of neutron ZDC. I am filling it with all the
                                 * possible TDC info. Basically it is the sum
                                 * of the all the 4 histograms.
                                 */
        TH1F*                   fZNATime4FillingH;         //!

                                /**
                                 * This histogram records the time distribution
                                 * of neutron ZDC. What we are
                                 * plotting here is the distribution of:
                                 * (ZNC-ZNA time) vs (ZNC+ZNA time).
                                 * Such a plot should show up all the different
                                 * collision possibilities (main-main,
                                 * main-satellite/satellite-main and
                                 * satellite-satellite).
                                 */
        TH2F*                   fZNCminusZNAtimeVsZNCplusZNAtimeH[4];  //!

                                /**
                                 * This histogram records the time distribution
                                 * of neutron ZDC. I am filling it with all the
                                 * possible TDC info. Basically it is the sum
                                 * of the all the 4 histograms. What we are
                                 * plotting here is the distribution of:
                                 * (ZNC-ZNA time) vs (ZNC+ZNA time).
                                 * Such a plot should show up all the different
                                 * collision possibilities (main-main,
                                 * main-satellite/satellite-main and
                                 * satellite-satellite).
                                 */
        TH2F*                   fZNCminusZNAtimeVsZNCplusZNAtime4FillingH;  //!

                                /**
                                 * This histogram records which information
                                 * I am looking at when I fill the time
                                 * histogram, and this effect should
                                 * reflect itself into the energy plots
                                 * because it shoudl give an indication of
                                 * what could possibly go wrong in the ZDC
                                 * analysis. ZNA case.
                                 */
        TH1F*                   fCounterZNAH;         //!

                                /**
                                 * This histogram records which information
                                 * I am looking at when I fill the time
                                 * histogram, and this effect should
                                 * reflect itself into the energy plots
                                 * because it shoudl give an indication of
                                 * what could possibly go wrong in the ZDC
                                 * analysis. ZNC case.
                                 */
        TH1F*                   fCounterZNCH;         //!

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
        TH1F*                   fInvariantMassDistributionAtDcaH;         //!

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
        // Finished cloned histograms.
        //_______________________________


        //_______________________________
        // DIFFERENTIAL NEUTRON EMISSION PLOTS

                                /**
                                 * This histogram records the invariant mass
                                 * distribution of the dimuon system, only
                                 * coherent component, so as to say, only
                                 * pt < 0.25 GeV/c for pt of the dimuon pair.
                                 * The ZNC has 0 neutrons.
                                 * The ZNA has 0 neutrons.
                                 */
        TH1F*                   fInvariantMassDistributionCoherentZNCzeroZNAzeroH;       //!

                                /**
                                 * This histogram records the invariant mass
                                 * distribution of the dimuon system, only
                                 * coherent component, so as to say, only
                                 * pt < 0.25 GeV/c for pt of the dimuon pair.
                                 * The ZNC has 0 neutrons.
                                 * The ZNA has ANY neutrons.
                                 */
        TH1F*                   fInvariantMassDistributionCoherentZNCzeroZNAanyH;       //!

                                /**
                                 * This histogram records the invariant mass
                                 * distribution of the dimuon system, only
                                 * coherent component, so as to say, only
                                 * pt < 0.25 GeV/c for pt of the dimuon pair.
                                 * The ZNC has ANY neutrons.
                                 * The ZNA has 0 neutrons.
                                 */
        TH1F*                   fInvariantMassDistributionCoherentZNCanyZNAzeroH;       //!

                                /**
                                 * This histogram records the invariant mass
                                 * distribution of the dimuon system, only
                                 * coherent component, so as to say, only
                                 * pt < 0.25 GeV/c for pt of the dimuon pair.
                                 * The ZNC has ANY neutrons.
                                 * The ZNA has ANY neutrons.
                                 */
        TH1F*                   fInvariantMassDistributionCoherentZNCanyZNAanyH;       //!

                                /**
                                 * This histogram records the invariant mass
                                 * distribution of the dimuon system, only
                                 * coherent component, so as to say, only
                                 * pt > 0.25 GeV/c for pt of the dimuon pair.
                                 * The ZNC has 0 neutrons.
                                 * The ZNA has 0 neutrons.
                                 */
        TH1F*                   fInvariantMassDistributionIncoherentZNCzeroZNAzeroH;       //!

                                /**
                                 * This histogram records the invariant mass
                                 * distribution of the dimuon system, only
                                 * coherent component, so as to say, only
                                 * pt > 0.25 GeV/c for pt of the dimuon pair.
                                 * The ZNC has 0 neutrons.
                                 * The ZNA has ANY neutrons.
                                 */
        TH1F*                   fInvariantMassDistributionIncoherentZNCzeroZNAanyH;       //!

                                /**
                                 * This histogram records the invariant mass
                                 * distribution of the dimuon system, only
                                 * coherent component, so as to say, only
                                 * pt > 0.25 GeV/c for pt of the dimuon pair.
                                 * The ZNC has ANY neutrons.
                                 * The ZNA has 0 neutrons.
                                 */
        TH1F*                   fInvariantMassDistributionIncoherentZNCanyZNAzeroH;       //!

                                /**
                                 * This histogram records the invariant mass
                                 * distribution of the dimuon system, only
                                 * coherent component, so as to say, only
                                 * pt > 0.25 GeV/c for pt of the dimuon pair.
                                 * The ZNC has ANY neutrons.
                                 * The ZNA has ANY neutrons.
                                 */
        TH1F*                   fInvariantMassDistributionIncoherentZNCanyZNAanyH;       //!

                                /**
                                 * This histogram shows the angular distribution
                                 * of the positive muon in the rest frame of the
                                 * J/Psi. This histogram is needed to evaluate
                                 * the polarization of the J/Psi!
                                 */
        TH1F*                   fAngularDistribOfPositiveMuonRestFrameJPsiH;       //!

                                /**
                                 * This histogram shows the angular distribution
                                 * of the negative muon in the rest frame of the
                                 * J/Psi. This histogram is needed to evaluate
                                 * the polarization of the J/Psi!
                                 */
        TH1F*                   fAngularDistribOfNegativeMuonRestFrameJPsiH;       //!

                                /**
                                 * This histogram represents the check over
                                 * the helicity of the J/Psi. It should be flat
                                 * if I remember well enough!
                                 */
        TH1F*                   fCheckHelicityRestFrameJPsiH;       //!

                                /**
                                 * This histogram shows the angular distribution
                                 * of the positive muon in the rest frame of the
                                 * J/Psi. This histogram is needed to evaluate
                                 * the polarization of the J/Psi! Divided in
                                 * psudorapidity bins... 8??
                                 */
        TH1F*                   fThetaDistribOfPositiveMuonRestFrameJPsiRapidityBinH[8];       //!



        // END DIFFERENTIAL NEUTRON EMISSION PLOTS
        //_______________________________

        //_______________________________
        // HELICITY AND COLLINS-SOPER

                                /**
                                 * This histogram shows the angular distribution
                                 * of the positive muon in the HELICITY frame.
                                 * COS(THETA) distribution.
                                 */
        TH1F*                   fCosThetaHelicityFrameJPsiH;       //!

                                /**
                                 * This histogram shows the angular distribution
                                 * of the positive muon in the HELICITY frame.
                                 * PHI distribution.
                                 */
        TH1F*                   fPhiHelicityFrameJPsiH;       //!

                                /**
                                 * This histogram shows the angular distribution
                                 * of the positive muon in the COLLINS-SOPER
                                 * frame.  COS(THETA) distribution.
                                 */
        TH1F*                   fCosThetaCollinsSoperFrameJPsiH;       //!

                                /**
                                 * This histogram shows the angular distribution
                                 * of the positive muon in the COLLINS-SOPER
                                 * frame. PHI distribution.
                                 */
        TH1F*                   fPhiCollinsSoperFrameJPsiH;       //!

                                /**
                                 * This histogram shows the angular distribution
                                 * of the positive muon in the HELICITY frame.
                                 * COS(THETA) distribution. Divided per
                                 * rapidity bins.
                                 */
        TH1F*                   fCosThetaHelicityFrameJPsiRapidityBinsH[8];       //!

                                /**
                                 * This histogram shows the angular distribution
                                 * of the positive muon in the HELICITY frame.
                                 * PHI distribution. Divided per
                                 * rapidity bins.
                                 */
        TH1F*                   fPhiHelicityFrameJPsiRapidityBinsH[8];       //!

                                /**
                                 * This histogram shows the angular distribution
                                 * of the positive muon in the COLLINS-SOPER
                                 * frame.  COS(THETA) distribution. Divided per
                                 * rapidity bins.
                                 */
        TH1F*                   fCosThetaCollinsSoperFrameJPsiRapidityBinsH[8];       //!

                                /**
                                 * This histogram shows the angular distribution
                                 * of the positive muon in the COLLINS-SOPER
                                 * frame. PHI distribution. Divided per
                                 * rapidity bins.
                                 */
        TH1F*                   fPhiCollinsSoperFrameJPsiRapidityBinsH[8];       //!

                                /**
                                 * This histogram shows the angular distribution
                                 * of the positive muon in the HELICITY frame.
                                 * COS(THETA) distribution. Divided per
                                 * rapidity bins. 10 of them.
                                 */
        TH1F*                   fCosThetaHelicityFrameJPsiTenRapidityBinsH[10];       //!

                                /**
                                 * This histogram shows the angular distribution
                                 * of the positive muon in the HELICITY frame.
                                 * PHI distribution. Divided per
                                 * rapidity bins. 10 of them.
                                 */
        TH1F*                   fPhiHelicityFrameJPsiTenRapidityBinsH[10];       //!

                                /**
                                 * This histogram shows the angular distribution
                                 * of the positive muon in the COLLINS-SOPER
                                 * frame.  COS(THETA) distribution. Divided per
                                 * rapidity bins. 10 of them.
                                 */
        TH1F*                   fCosThetaCollinsSoperFrameJPsiTenRapidityBinsH[10];       //!

                                /**
                                 * This histogram shows the angular distribution
                                 * of the positive muon in the COLLINS-SOPER
                                 * frame. PHI distribution. Divided per
                                 * rapidity bins. 10 of them.
                                 */
        TH1F*                   fPhiCollinsSoperFrameJPsiTenRapidityBinsH[10];       //!

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
                                 */
        TH1F*                   fInvariantMassDistributionInBinsOfCosThetaHelicityFrameH[10];       //!

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
                                 * NEW: This histogram shows this thing in 2D.
                                 * IT is shown the same distribution in terms of
                                 * CosTheta and Phi bins. Let's see the results!
                                 *
                                 * NEW: this specific histogram follows the
                                 * binning of the inclusive people...
                                 * Hopefully this would yield better results!
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
                                 * NEW: the mass range has been extended.
                                 * Meaning there are no bounds on the invariant
                                 * mass range BEFORE signal extraction.
                                 */
        TH1F***                 fInvariantMassDistributionForSignalExtractionHelicityFrameH;  //!

        //                         /**
        //                          * This histogram shows the invariant mass
        //                          * distribution of the dimuon pairs in terms
        //                          * of bins of cos theta of the positive muon
        //                          * in the helicity frame of the J/Psi.
        //                          *
        //                          * What it means is that we divide in 5 bins of
        //                          * possible CosTheta of the decaying J/Psi,
        //                          * meaning  (-1,-0.8), (-0.8,-0.6), (-0.6,-0.4),
        //                          * (-0.4,-0.2) and so on until (0.8,1). We fill
        //                          * the invariant mass distribution of the
        //                          * dimuons in this many bins.
        //                          *
        //                          * The next step is to fit this invariant mass
        //                          * distributions, so as to obtain the relative
        //                          * contribution of J/Psi and GammaGamma to the
        //                          * angular distributions. This should help in
        //                          * validating our results...
        //                          *
        //                          * NEW: the mass range has been extended.
        //                          * Meaning there are no bounds on the invariant
        //                          * mass range BEFORE signal extraction.
        //                          */
        // MatrixTH1F              fInvariantMassDistributionForSignalExtractionHelicityFrameH;  //!


        //_______________________________
        // CUTS
        /*
         * The following is all the possible checks for the event selections
         * and the track selection as well. Enjoy.
         */
        Int_t                   fRunNum;        //!
        Int_t                   fTracklets;     //!

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

        /**
         * This is the vector containing the GOOD RunNumbers.
         */
        std::vector<Int_t> fVectorGoodRunNumbers;       //!



        /**
         * Not implemented yet...
         */
        AliAnalysisTaskUPCforward2(const AliAnalysisTaskUPCforward2&);

        /**
         * Not implemented yet...
         */
        AliAnalysisTaskUPCforward2& operator=(const AliAnalysisTaskUPCforward2&);


        /**
         * This is important for ROOT only. I do not remember the reason anymore.
         * If I happen to encounter it again in the future, I will make sure to
         * record it!
         */
        ClassDef(AliAnalysisTaskUPCforward2, 1);
};

#endif
