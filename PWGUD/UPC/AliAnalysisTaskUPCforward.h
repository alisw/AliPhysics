/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskUPCforward_H
#define AliAnalysisTaskUPCforward_H

/**
 * @file   AliAnalysisTaskUPCforward.h
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
// typedef std::vector< std::vector< TH1F* > >  MatrixTH1F;

/**
 * \file AliAnalysisTaskUPCforward.h
 * \brief Contains the declaration of the AliAnalysisTaskUPCforward class
 */

/**
 * \class AliAnalysisTaskUPCforward
 * \brief Contains the tools to every part of my analysis
 */
class AliAnalysisTaskUPCforward : public AliAnalysisTaskSE
{
    public:
                                /**
                                 * Create a new AliAnalysisTaskUPCforward with
                                 * default constructor. Based on my understanding
                                 * this is mostly related to input/output
                                 * processes.
                                 */
                                AliAnalysisTaskUPCforward();

                                /**
                                 * Create a new AliAnalysisTaskUPCforward with
                                 * a specific name. This serves to instantiate
                                 * the AliAnalysisTaskSE object built with the
                                 * constructor, that will take this specific
                                 * name.
                                 *
                                 * \param name , the name taken by the AliAnalysisTaskSE object.
                                 */
                                AliAnalysisTaskUPCforward( const char *name, Int_t _fSetSingleMuonPt );

                                /**
                                 * Virtual destructor for the class. It will
                                 * take care of all the particular objects to be
                                 * destroyed for the class.
                                 */
        virtual                 ~AliAnalysisTaskUPCforward();

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
        // void                    FillGoodRunVector(std::vector<Int_t> &fVectorGoodRunNumbers);

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
                                 * This function computes the acceptance
                                 * value to fix the weight of the entry
                                 * to the CosTheta in the helicity frame.
                                 */
        Double_t                AccEffCorrection( Double_t CosThetaToBeWeighted );

                                /**
                                 * This function sets the muon pt for the
                                 * Trigger Efficiency computation in the
                                 * range (0.85->1.15) GeV with 0.05 Gev steps.
                                 *
                                 * The aim is to check how the different
                                 * selections change the polarisation
                                 * parameters for now, maybe in the future
                                 * the Neutron Emission Yield.
                                 */
        void                    SetSingleMuonPt( Int_t _fSetSingleMuonPt ){ fSetSingleMuonPt = _fSetSingleMuonPt; }

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
                                 * The flag to select the single muons
                                 * Pt for the Trigger Efficiency
                                 * systematic source.
                                 */
        Int_t                   fSetSingleMuonPt;

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
                                 * Pseudorapidity distribution of the many muons.
                                 */
        TH1F*                   fEtaMuonH;          //!

                                /**
                                 * Pt distribution of the many muons.
                                 */
        TH1F*                   fPtSingleMuonH;          //!

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
                                 * This histogram records the number of entries
                                 * against the respective run number.
                                 * But only for the CMUP11 Trigger Class.
                                 */
        TH1F*                   fRunNumberTriggerCMUP11ClassH;         //!

                                /**
                                 * This histogram records the number of entries
                                 * against the respective run number, but the
                                 * binning is done properly, like normally
                                 * found for this type of histograms.
                                 * But only for the CMUP11 Trigger Class.
                                 */
        TH1F*                   fRunNumberTriggerCMUP11ClassProperlyH;         //!

                                /**
                                 * This histogram records the number of entries
                                 * against the respective run number.
                                 * But only for the CMUP26 Trigger Class.
                                 */
        TH1F*                   fRunNumberTriggerCMUP26ClassH;         //!

                                /**
                                 * This histogram records the number of entries
                                 * against the respective run number, but the
                                 * binning is done properly, like normally
                                 * found for this type of histograms.
                                 * But only for the CMUP26 Trigger Class.
                                 */
        TH1F*                   fRunNumberTriggerCMUP26ClassProperlyH;         //!

                                /**
                                 * This histogram records the number of entries
                                 * against the respective run number.
                                 * But only for the CMUP6 Trigger Class.
                                 */
        TH1F*                   fRunNumberTriggerCMUP6ClassH;         //!

                                /**
                                 * This histogram records the number of entries
                                 * against the respective run number, but the
                                 * binning is done properly, like normally
                                 * found for this type of histograms.
                                 * But only for the CMUP6 Trigger Class.
                                 */
        TH1F*                   fRunNumberTriggerCMUP6ClassProperlyH;         //!

                                /**
                                 * This histogram records the number of entries
                                 * against the respective run number.
                                 * But only for the CMUP10 Trigger Class.
                                 */
        TH1F*                   fRunNumberTriggerCMUP10ClassH;         //!

                                /**
                                 * This histogram records the number of entries
                                 * against the respective run number, but the
                                 * binning is done properly, like normally
                                 * found for this type of histograms.
                                 * But only for the CMUP10 Trigger Class.
                                 */
        TH1F*                   fRunNumberTriggerCMUP10ClassProperlyH;         //!

                                /**
                                 * This histogram records the number of entries
                                 * against the respective run number.
                                 * But only for the CMUP13 Trigger Class.
                                 */
        TH1F*                   fRunNumberTriggerCMUP13ClassH;         //!

                                /**
                                 * This histogram records the number of entries
                                 * against the respective run number, but the
                                 * binning is done properly, like normally
                                 * found for this type of histograms.
                                 * But only for the CMUP13 Trigger Class.
                                 */
        TH1F*                   fRunNumberTriggerCMUP13ClassProperlyH;         //!

                                /**
                                 * This histogram records the number of entries
                                 * against the respective run number.
                                 * This is Evgeny's style, so it is a
                                 * TH2F where the x-axis represents the
                                 * trigger class...
                                 */
        TH2F*                   fTriggersVsRunH;         //!

                                /**
                                 * This histogram records the invariant mass
                                 * distribution of the dimuon system, only
                                 * coherent component, so as to say, only
                                 * pt < 0.25 GeV/c for pt of the dimuon pair.
                                 *
                                 * The array is simply the same plot but
                                 * divided in rapidity bins.
                                 * From 0 up to 6 it is Y = -4 to -2.5.
                                 */
        TH1F*                   fInvariantMassDistributionCoherentH;                      //!
        TH1F*                   fInvariantMassDistributionCoherentRapidityBinsH[6];       //!

                                /**
                                 * This histogram records the invariant mass
                                 * distribution of the dimuon system, only
                                 * coherent component, so as to say, only
                                 * pt < 0.2 (MinusTwoShift) GeV/c for pt
                                 * of the dimuon pair.
                                 */
        TH1F*                   fInvariantMassDistributionCoherentShiftMinusTwoH;       //!

                                /**
                                 * This histogram records the invariant mass
                                 * distribution of the dimuon system, only
                                 * coherent component, so as to say, only
                                 * pt < 0.225 (MinusOneShift) GeV/c for
                                 * pt of the dimuon pair.
                                 */
        TH1F*                   fInvariantMassDistributionCoherentShiftMinusOneH;       //!

                                /**
                                 * This histogram records the invariant mass
                                 * distribution of the dimuon system, only
                                 * coherent component, so as to say, only
                                 * pt < 0.275 (PlusOneShift) GeV/c for
                                 * pt of the dimuon pair.
                                 */
        TH1F*                   fInvariantMassDistributionCoherentShiftPlusOneH;       //!

                                /**
                                 * This histogram records the invariant mass
                                 * distribution of the dimuon system, only
                                 * coherent component, so as to say, only
                                 * pt < 0.3 (PlusTwoShift) GeV/c for
                                 * pt of the dimuon pair.
                                 */
        TH1F*                   fInvariantMassDistributionCoherentShiftPlusTwoH;       //!

                                /**
                                 * This histogram records the invariant mass
                                 * distribution of the dimuon system, only
                                 * coherent component, so as to say, only
                                 * pt > 0.25 GeV/c for pt of the dimuon pair.
                                 */
        TH1F*                   fInvariantMassDistributionIncoherentH;                   //!
        TH1F*                   fInvariantMassDistributionIncoherentShiftMinusTwoH;      //!
        TH1F*                   fInvariantMassDistributionIncoherentShiftMinusOneH;      //!
        TH1F*                   fInvariantMassDistributionIncoherentShiftPlusOneH;       //!
        TH1F*                   fInvariantMassDistributionIncoherentShiftPlusTwoH;       //!

                                /**
                                 * This histogram records the pt-ditribution
                                 * of the dimuon pairs.
                                 *
                                 * Shift +1 => 20 Mev/c shift
                                 */
        TH1F*                   fDimuonPtDistributionH;                     //!
        TH1F*                   fDimuonPtDistributionShiftPlusOneH;         //!

                                /**
                                 * This histogram records the energy distri-
                                 * bution of the neutron ZDC. This plot should
                                 * show the relative components of the 0 neutron
                                 * peak, the 1 neutron peak and possibly the
                                 * 2 neutrons peak. Anything higher than that,
                                 * requires help from the user and is more like
                                 * a guess...
                                 */
        TH1F*                   fZNCEnergyAgainstEntriesH;                 //!
        TH1F*                   fZNCEnergyAgainstEntriesExtendedH;         //!
        TH1F*                   fZNCEnergyAgainstEntriesExtendedHv2;       //!

                                /**
                                 * This histogram records the energy distri-
                                 * bution of the neutron ZDC. This plot should
                                 * show the relative components of the 0 neutron
                                 * peak, the 1 neutron peak and possibly the
                                 * 2 neutrons peak. Anything higher than that,
                                 * requires help from the user and is more like
                                 * a guess...
                                 */
        TH1F*                   fZNAEnergyAgainstEntriesH;                 //!
        TH1F*                   fZNAEnergyAgainstEntriesExtendedH;         //!
        TH1F*                   fZNAEnergyAgainstEntriesExtendedHv2;       //!

                                /**
                                 * This histogram records the energy distri-
                                 * bution of the neutron ZDC. BEFORE timing
                                 * selection.
                                 */
        TH1F*                   fZNCEnergyBeforeTimingSelectionH;                 //!
        TH1F*                   fZNCEnergyBeforeTimingSelectionExtendedH;         //!

                                /**
                                 * This histogram records the energy distri-
                                 * bution of the neutron ZDC. BEFORE timing
                                 * selection.
                                 */
        TH1F*                   fZNAEnergyBeforeTimingSelectionH;                 //!
        TH1F*                   fZNAEnergyBeforeTimingSelectionExtendedH;         //!

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

                                /**
                                 * This histogram records the DCA vs the Pt
                                 * of the vector meson.
                                 */
        TH2F*                   fDcaAgainstPtOfVectorMesonH;         //!

                                /**
                                 * This histogram records the invariant mass
                                 * distribution of the dimuon system with a
                                 * strict Pt and a strics DCA.
                                 */
        TH1F*                   fInvariantMassDistributionStrictPtStrictDcaH;         //!

                                /**
                                 * This histogram records the invariant mass
                                 * distribution of the dimuon system VS Pt
                                 * and a strics DCA.
                                 */
        TH2F*                   fInvariantMassDistributionVsPtStrictDcaH;         //!

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
                                 *
                                 * Shift -2 => pt < 0.200 GeV/c
                                 * Shift -1 => pt < 0.225 GeV/c
                                 * Shift +1 => pt < 0.275 GeV/c
                                 * Shift +2 => pt < 0.300 GeV/c
                                 *
                                 *
                                 * NEW: v2 uses the isZNAfired() methods
                                 * instead of the threshold to extract the
                                 * numbers.
                                 */
        TH1F*                   fInvariantMassDistributionCoherentZNCzeroZNAzeroH;                   //!
        TH1F*                   fInvariantMassDistributionCoherentZNCzeroZNAzeroShiftMinusTwoH;      //!
        TH1F*                   fInvariantMassDistributionCoherentZNCzeroZNAzeroShiftMinusOneH;      //!
        TH1F*                   fInvariantMassDistributionCoherentZNCzeroZNAzeroShiftPlusOneH;       //!
        TH1F*                   fInvariantMassDistributionCoherentZNCzeroZNAzeroShiftPlusTwoH;       //!
        TH1F*                   fInvariantMassDistributionCoherentZNCzeroZNAzeroHv2;                 //!
        TH1F*                   fInvariantMassDistributionCoherentZNCzeroZNAzeroRapidityHv2[3];      //!


                                /**
                                 * This histogram records the invariant mass
                                 * distribution of the dimuon system, only
                                 * coherent component, so as to say, only
                                 * pt < 0.25 GeV/c for pt of the dimuon pair.
                                 * The ZNC has 0 neutrons.
                                 * The ZNA has ANY neutrons.
                                 *
                                 * Shift -2 => pt < 0.200 GeV/c
                                 * Shift -1 => pt < 0.225 GeV/c
                                 * Shift +1 => pt < 0.275 GeV/c
                                 * Shift +2 => pt < 0.300 GeV/c
                                 *
                                 *
                                 * NEW: v2 uses the isZNAfired() methods
                                 * instead of the threshold to extract the
                                 * numbers.
                                 */
        TH1F*                   fInvariantMassDistributionCoherentZNCzeroZNAanyH;                   //!
        TH1F*                   fInvariantMassDistributionCoherentZNCzeroZNAanyShiftMinusTwoH;      //!
        TH1F*                   fInvariantMassDistributionCoherentZNCzeroZNAanyShiftMinusOneH;      //!
        TH1F*                   fInvariantMassDistributionCoherentZNCzeroZNAanyShiftPlusOneH;       //!
        TH1F*                   fInvariantMassDistributionCoherentZNCzeroZNAanyShiftPlusTwoH;       //!
        TH1F*                   fInvariantMassDistributionCoherentZNCzeroZNAanyHv2;                 //!
        TH1F*                   fInvariantMassDistributionCoherentZNCzeroZNAanyRapidityHv2[3];      //!

                                /**
                                 * This histogram records the invariant mass
                                 * distribution of the dimuon system, only
                                 * coherent component, so as to say, only
                                 * pt < 0.25 GeV/c for pt of the dimuon pair.
                                 * The ZNC has ANY neutrons.
                                 * The ZNA has 0 neutrons.
                                 *
                                 * Shift -2 => pt < 0.200 GeV/c
                                 * Shift -1 => pt < 0.225 GeV/c
                                 * Shift +1 => pt < 0.275 GeV/c
                                 * Shift +2 => pt < 0.300 GeV/c
                                 *
                                 *
                                 * NEW: v2 uses the isZNAfired() methods
                                 * instead of the threshold to extract the
                                 * numbers.
                                 */
        TH1F*                   fInvariantMassDistributionCoherentZNCanyZNAzeroH;                   //!
        TH1F*                   fInvariantMassDistributionCoherentZNCanyZNAzeroShiftMinusTwoH;      //!
        TH1F*                   fInvariantMassDistributionCoherentZNCanyZNAzeroShiftMinusOneH;      //!
        TH1F*                   fInvariantMassDistributionCoherentZNCanyZNAzeroShiftPlusOneH;       //!
        TH1F*                   fInvariantMassDistributionCoherentZNCanyZNAzeroShiftPlusTwoH;       //!
        TH1F*                   fInvariantMassDistributionCoherentZNCanyZNAzeroHv2;                 //!
        TH1F*                   fInvariantMassDistributionCoherentZNCanyZNAzeroRapidityHv2[3];      //!

                                /**
                                 * This histogram records the invariant mass
                                 * distribution of the dimuon system, only
                                 * coherent component, so as to say, only
                                 * pt < 0.25 GeV/c for pt of the dimuon pair.
                                 * The ZNC has ANY neutrons.
                                 * The ZNA has ANY neutrons.
                                 *
                                 * Shift -2 => pt < 0.200 GeV/c
                                 * Shift -1 => pt < 0.225 GeV/c
                                 * Shift +1 => pt < 0.275 GeV/c
                                 * Shift +2 => pt < 0.300 GeV/c
                                 *
                                 *
                                 * NEW: v2 uses the isZNAfired() methods
                                 * instead of the threshold to extract the
                                 * numbers.
                                 */
        TH1F*                   fInvariantMassDistributionCoherentZNCanyZNAanyH;                   //!
        TH1F*                   fInvariantMassDistributionCoherentZNCanyZNAanyShiftMinusTwoH;      //!
        TH1F*                   fInvariantMassDistributionCoherentZNCanyZNAanyShiftMinusOneH;      //!
        TH1F*                   fInvariantMassDistributionCoherentZNCanyZNAanyShiftPlusOneH;       //!
        TH1F*                   fInvariantMassDistributionCoherentZNCanyZNAanyShiftPlusTwoH;       //!
        TH1F*                   fInvariantMassDistributionCoherentZNCanyZNAanyHv2;                 //!
        TH1F*                   fInvariantMassDistributionCoherentZNCanyZNAanyRapidityHv2[3];      //!

                                /**
                                 * This histogram records the invariant mass
                                 * distribution of the dimuon system, only
                                 * coherent component, so as to say, only
                                 * pt > 0.25 GeV/c for pt of the dimuon pair.
                                 * The ZNC has 0 neutrons.
                                 * The ZNA has 0 neutrons.
                                 *
                                 * Shift -2 => pt > 0.200 GeV/c
                                 * Shift -1 => pt > 0.225 GeV/c
                                 * Shift +1 => pt > 0.275 GeV/c
                                 * Shift +2 => pt > 0.300 GeV/c
                                 *
                                 *
                                 * NEW: v2 uses the isZNAfired() methods
                                 * instead of the threshold to extract the
                                 * numbers.
                                 */
        TH1F*                   fInvariantMassDistributionIncoherentZNCzeroZNAzeroH;                   //!
        TH1F*                   fInvariantMassDistributionIncoherentZNCzeroZNAzeroShiftMinusTwoH;      //!
        TH1F*                   fInvariantMassDistributionIncoherentZNCzeroZNAzeroShiftMinusOneH;      //!
        TH1F*                   fInvariantMassDistributionIncoherentZNCzeroZNAzeroShiftPlusOneH;       //!
        TH1F*                   fInvariantMassDistributionIncoherentZNCzeroZNAzeroShiftPlusTwoH;       //!
        TH1F*                   fInvariantMassDistributionIncoherentZNCzeroZNAzeroHv2;                 //!
        TH1F*                   fInvariantMassDistributionIncoherentZNCzeroZNAzeroRapidityHv2[3];      //!

                                /**
                                 * This histogram records the invariant mass
                                 * distribution of the dimuon system, only
                                 * coherent component, so as to say, only
                                 * pt > 0.25 GeV/c for pt of the dimuon pair.
                                 * The ZNC has 0 neutrons.
                                 * The ZNA has ANY neutrons.
                                 *
                                 * Shift -2 => pt > 0.200 GeV/c
                                 * Shift -1 => pt > 0.225 GeV/c
                                 * Shift +1 => pt > 0.275 GeV/c
                                 * Shift +2 => pt > 0.300 GeV/c
                                 *
                                 *
                                 * NEW: v2 uses the isZNAfired() methods
                                 * instead of the threshold to extract the
                                 * numbers.
                                 */
        TH1F*                   fInvariantMassDistributionIncoherentZNCzeroZNAanyH;                   //!
        TH1F*                   fInvariantMassDistributionIncoherentZNCzeroZNAanyShiftMinusTwoH;      //!
        TH1F*                   fInvariantMassDistributionIncoherentZNCzeroZNAanyShiftMinusOneH;      //!
        TH1F*                   fInvariantMassDistributionIncoherentZNCzeroZNAanyShiftPlusOneH;       //!
        TH1F*                   fInvariantMassDistributionIncoherentZNCzeroZNAanyShiftPlusTwoH;       //!
        TH1F*                   fInvariantMassDistributionIncoherentZNCzeroZNAanyHv2;                 //!
        TH1F*                   fInvariantMassDistributionIncoherentZNCzeroZNAanyRapidityHv2[3];      //!

                                /**
                                 * This histogram records the invariant mass
                                 * distribution of the dimuon system, only
                                 * coherent component, so as to say, only
                                 * pt > 0.25 GeV/c for pt of the dimuon pair.
                                 * The ZNC has ANY neutrons.
                                 * The ZNA has 0 neutrons.
                                 *
                                 * Shift -2 => pt > 0.200 GeV/c
                                 * Shift -1 => pt > 0.225 GeV/c
                                 * Shift +1 => pt > 0.275 GeV/c
                                 * Shift +2 => pt > 0.300 GeV/c
                                 *
                                 *
                                 * NEW: v2 uses the isZNAfired() methods
                                 * instead of the threshold to extract the
                                 * numbers.
                                 */
        TH1F*                   fInvariantMassDistributionIncoherentZNCanyZNAzeroH;                   //!
        TH1F*                   fInvariantMassDistributionIncoherentZNCanyZNAzeroShiftMinusTwoH;      //!
        TH1F*                   fInvariantMassDistributionIncoherentZNCanyZNAzeroShiftMinusOneH;      //!
        TH1F*                   fInvariantMassDistributionIncoherentZNCanyZNAzeroShiftPlusOneH;       //!
        TH1F*                   fInvariantMassDistributionIncoherentZNCanyZNAzeroShiftPlusTwoH;       //!
        TH1F*                   fInvariantMassDistributionIncoherentZNCanyZNAzeroHv2;                 //!
        TH1F*                   fInvariantMassDistributionIncoherentZNCanyZNAzeroRapidityHv2[3];      //!

                                /**
                                 * This histogram records the invariant mass
                                 * distribution of the dimuon system, only
                                 * coherent component, so as to say, only
                                 * pt > 0.25 GeV/c for pt of the dimuon pair.
                                 * The ZNC has ANY neutrons.
                                 * The ZNA has ANY neutrons.
                                 *
                                 * Shift -2 => pt > 0.200 GeV/c
                                 * Shift -1 => pt > 0.225 GeV/c
                                 * Shift +1 => pt > 0.275 GeV/c
                                 * Shift +2 => pt > 0.300 GeV/c
                                 *
                                 *
                                 * NEW: v2 uses the isZNAfired() methods
                                 * instead of the threshold to extract the
                                 * numbers.
                                 */
        TH1F*                   fInvariantMassDistributionIncoherentZNCanyZNAanyH;                   //!
        TH1F*                   fInvariantMassDistributionIncoherentZNCanyZNAanyShiftMinusTwoH;      //!
        TH1F*                   fInvariantMassDistributionIncoherentZNCanyZNAanyShiftMinusOneH;      //!
        TH1F*                   fInvariantMassDistributionIncoherentZNCanyZNAanyShiftPlusOneH;       //!
        TH1F*                   fInvariantMassDistributionIncoherentZNCanyZNAanyShiftPlusTwoH;       //!
        TH1F*                   fInvariantMassDistributionIncoherentZNCanyZNAanyHv2;                 //!
        TH1F*                   fInvariantMassDistributionIncoherentZNCanyZNAanyRapidityHv2[3];      //!

                                /**
                                 * This histogram records the pt-ditribution
                                 * of the dimuon pairs. This is divided in
                                 * neutron emission classes.
                                 * ZNC=0n, ZNA=0n.
                                 *
                                 * Shift +1 => 20 MeV/c shift
                                 *
                                 *
                                 * NEW: v2 uses the isZNAfired() methods
                                 * instead of the threshold to extract the
                                 * numbers.
                                 */
        TH1F*                   fDimuonPtDistributionZNCzeroZNAzeroH;                     //!
        TH1F*                   fDimuonPtDistributionZNCzeroZNAzeroShiftPlusOneH;         //!
        TH1F*                   fDimuonPtDistributionZNCzeroZNAzeroHv2;                   //!
        TH1F*                   fDimuonPtDistributionZNCzeroZNAzeroRapidityHv2[3];        //!

                                /**
                                 * This histogram records the pt-ditribution
                                 * of the dimuon pairs. This is divided in
                                 * neutron emission classes.
                                 * ZNC=0n, ZNA=Xn.
                                 *
                                 * Shift +1 => 20 MeV/c shift
                                 *
                                 *
                                 * NEW: v2 uses the isZNAfired() methods
                                 * instead of the threshold to extract the
                                 * numbers.
                                 */
        TH1F*                   fDimuonPtDistributionZNCzeroZNAanyH;                     //!
        TH1F*                   fDimuonPtDistributionZNCzeroZNAanyShiftPlusOneH;         //!
        TH1F*                   fDimuonPtDistributionZNCzeroZNAanyHv2;                   //!
        TH1F*                   fDimuonPtDistributionZNCzeroZNAanyRapidityHv2[3];        //!

                                /**
                                 * This histogram records the pt-ditribution
                                 * of the dimuon pairs. This is divided in
                                 * neutron emission classes.
                                 * ZNC=Xn, ZNA=0n.
                                 *
                                 * Shift +1 => 20 MeV/c shift
                                 *
                                 *
                                 * NEW: v2 uses the isZNAfired() methods
                                 * instead of the threshold to extract the
                                 * numbers.
                                 */
        TH1F*                   fDimuonPtDistributionZNCanyZNAzeroH;                     //!
        TH1F*                   fDimuonPtDistributionZNCanyZNAzeroShiftPlusOneH;         //!
        TH1F*                   fDimuonPtDistributionZNCanyZNAzeroHv2;                   //!
        TH1F*                   fDimuonPtDistributionZNCanyZNAzeroRapidityHv2[3];        //!

                                /**
                                 * This histogram records the pt-ditribution
                                 * of the dimuon pairs. This is divided in
                                 * neutron emission classes and pt.
                                 * ZNC=Xn, ZNA=Xn.
                                 *
                                 * Shift +1 => 20 MeV/c shift
                                 *
                                 *
                                 * NEW: v2 uses the isZNAfired() methods
                                 * instead of the threshold to extract the
                                 * numbers.
                                 */
        TH1F*                   fDimuonPtDistributionZNCanyZNAanyH;                     //!
        TH1F*                   fDimuonPtDistributionZNCanyZNAanyShiftPlusOneH;         //!
        TH1F*                   fDimuonPtDistributionZNCanyZNAanyHv2;                   //!
        TH1F*                   fDimuonPtDistributionZNCanyZNAanyRapidityHv2[3];        //!

                                /**
                                 * This histogram records the pt-ditribution
                                 * of the dimuon pairs. This is divided in
                                 * neutron emission classes and pt.
                                 * COHERENT, ZNC=0n, ZNA=0n.
                                 */
        TH1F*                   fDimuonPtDistributionCoherentZNCzeroZNAzeroH;         //!

                                /**
                                 * This histogram records the pt-ditribution
                                 * of the dimuon pairs. This is divided in
                                 * neutron emission classes and pt.
                                 * COHERENT, ZNC=0n, ZNA=Xn.
                                 */
        TH1F*                   fDimuonPtDistributionCoherentZNCzeroZNAanyH;         //!

                                /**
                                 * This histogram records the pt-ditribution
                                 * of the dimuon pairs. This is divided in
                                 * neutron emission classes and pt.
                                 * COHERENT, ZNC=Xn, ZNA=0n.
                                 */
        TH1F*                   fDimuonPtDistributionCoherentZNCanyZNAzeroH;         //!

                                /**
                                 * This histogram records the pt-ditribution
                                 * of the dimuon pairs. This is divided in
                                 * neutron emission classes and pt.
                                 * COHERENT, ZNC=Xn, ZNA=Xn.
                                 */
        TH1F*                   fDimuonPtDistributionCoherentZNCanyZNAanyH;         //!

                                /**
                                 * This histogram records the pt-ditribution
                                 * of the dimuon pairs. This is divided in
                                 * neutron emission classes and pt.
                                 * INCOHERENT, ZNC=0n, ZNA=0n.
                                 */
        TH1F*                   fDimuonPtDistributionIncoherentZNCzeroZNAzeroH;         //!

                                /**
                                 * This histogram records the pt-ditribution
                                 * of the dimuon pairs. This is divided in
                                 * neutron emission classes and pt.
                                 * INCOHERENT, ZNC=0n, ZNA=Xn.
                                 */
        TH1F*                   fDimuonPtDistributionIncoherentZNCzeroZNAanyH;         //!

                                /**
                                 * This histogram records the pt-ditribution
                                 * of the dimuon pairs. This is divided in
                                 * neutron emission classes and pt.
                                 * INCOHERENT, ZNC=Xn, ZNA=0n.
                                 */
        TH1F*                   fDimuonPtDistributionIncoherentZNCanyZNAzeroH;         //!

                                /**
                                 * This histogram records the pt-ditribution
                                 * of the dimuon pairs. This is divided in
                                 * neutron emission classes and pt.
                                 * INCOHERENT, ZNC=Xn, ZNA=Xn.
                                 */
        TH1F*                   fDimuonPtDistributionIncoherentZNCanyZNAanyH;         //!


        //_______________________________

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
        TH1F*                   fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameTryingH[40];  //!

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
        TH1F*                   fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameTryingH[50];  //!

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
                                 *
                                 *
                                 * NB: 0N0N implies the use of the neutron
                                 *     emission class to suppress the
                                 *     feed-down contribution and the
                                 *     incoherent component too!!
                                 */
        TH1F***                 fInvariantMassDistributionForSignalExtractionHelicityFrameMyBinningH;      //!
        TH1F***                 fInvariantMassDistributionForSignalExtractionCsFrameMyBinningH;            //!

        TH1F***                 fInvariantMassDistributionForSignalExtractionHelicityFrameMyBinning0N0NH;  //!
        TH1F***                 fInvariantMassDistributionForSignalExtractionCsFrameMyBinning0N0NH;        //!


                                /**
                                 * This histogram shows the invariant mass
                                 * distribution of the dimuon pairs in a
                                 * restricted bin of pt: 0.4 < pt < 0.6 GeV/c.
                                 */
        TH1F*                   fInvariantMassDistributionStrictPtH;  //!

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

                                /**
                                 * This histogram shows the invariant mass
                                 * distribution of the dimuon pairs in terms
                                 * of bins of cos theta of the positive muon
                                 * in the helicity frame of the J/Psi.
                                 * My STRICT variable binning.
                                 */
        TH1F*                   fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameMyStrictVariableBinningH[9];  //!

                                /**
                                 * This histogram shows the invariant mass
                                 * distribution of the dimuon pairs in terms
                                 * of bins of cos theta of the positive muon
                                 * in the helicity frame of the J/Psi.
                                 * My STRICT variable binning.
                                 */
        TH1F*                   fInvariantMassDistributionOnlyPhiForSignalExtractionHelicityFrameMyStrictVariableBinningH[15];  //!

                                /**
                                 * This histogram shows the invariant mass
                                 * distribution of the dimuon pairs in terms
                                 * of bins of cos theta of the positive muon
                                 * in the helicity frame of the J/Psi.
                                 * My variable binning 17 bins.
                                 */
        TH1F*                   fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameMySeventeenBinsVariableBinningH[17];  //!

                                /**
                                 * This histogram shows the angular distribution
                                 * of the positive muon in the HELICITY frame.
                                 * COS(THETA) distribution.
                                 * Already corrected by ACCxEFF thanks to
                                 * a weight...
                                 *
                                 * Plus the distributions in invariant mass
                                 * for signal extraction.
                                 */
        TH1F*                   fCosThetaHelicityFrameJPsiAlreadyCorrectedH;       //!
        TH1F*                   fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameForAlreadyCorrectedFiftyH[50];  //!
        TH1F*                   fInvariantMassDistributionOnlyCosThetaForSignalExtractionHelicityFrameForAlreadyCorrectedHundredH[100];  //!

                                /**
                                 * Signal extraction in Phi, CosTheta, and
                                 * TildePhi, with 25 bins only...
                                 *
                                 * NB: He => helicity frame
                                 * NB: Cs => CS       frame
                                 *
                                 * NB: 0N0N implies the use of the neutron
                                 *     emission class to suppress the
                                 *     feed-down contribution and the
                                 *     incoherent component too!!
                                 */
        TH1F*                   fInvariantMassDistributionOnlyPhiHeFrameTwentyfiveBinsH[25];           //!
        TH1F*                   fInvariantMassDistributionOnlyCosThetaHeFrameTwentyfiveBinsH[25];      //!
        TH1F*                   fInvariantMassDistributionOnlyTildePhiHeFrameTwentyfiveBinsH[25];      //!

        TH1F*                   fInvariantMassDistributionOnlyPhiCsFrameTwentyfiveBinsH[25];           //!
        TH1F*                   fInvariantMassDistributionOnlyCosThetaCsFrameTwentyfiveBinsH[25];      //!
        TH1F*                   fInvariantMassDistributionOnlyTildePhiCsFrameTwentyfiveBinsH[25];      //!

        TH1F*                   fInvariantMassDistributionOnlyPhiHeFrameTwentyfiveBins0N0NH[25];       //!
        TH1F*                   fInvariantMassDistributionOnlyCosThetaHeFrameTwentyfiveBins0N0NH[25];  //!
        TH1F*                   fInvariantMassDistributionOnlyTildePhiHeFrameTwentyfiveBins0N0NH[25];  //!

        TH1F*                   fInvariantMassDistributionOnlyPhiCsFrameTwentyfiveBins0N0NH[25];       //!
        TH1F*                   fInvariantMassDistributionOnlyCosThetaCsFrameTwentyfiveBins0N0NH[25];  //!
        TH1F*                   fInvariantMassDistributionOnlyTildePhiCsFrameTwentyfiveBins0N0NH[25];  //!

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
        // std::vector<Int_t> fVectorGoodRunNumbers;       //!



        /**
         * Not implemented yet...
         */
        AliAnalysisTaskUPCforward(const AliAnalysisTaskUPCforward&);

        /**
         * Not implemented yet...
         */
        AliAnalysisTaskUPCforward& operator=(const AliAnalysisTaskUPCforward&);


        /**
         * This is important for ROOT only. I do not remember the reason anymore.
         * If I happen to encounter it again in the future, I will make sure to
         * record it!
         */
        ClassDef(AliAnalysisTaskUPCforward, 45);
};

#endif
