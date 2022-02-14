/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskADin2018_H
#define AliAnalysisTaskADin2018_H

/**
 * @file   AliAnalysisTaskADin2018.h
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
 * \file AliAnalysisTaskADin2018.h
 * \brief Contains the declaration of the AliAnalysisTaskADin2018 class
 */

/**
 * \class AliAnalysisTaskADin2018
 * \brief Contains the tools to every part of my analysis
 */
class AliAnalysisTaskADin2018 : public AliAnalysisTaskSE
{
    public:
                                /**
                                 * Create a new AliAnalysisTaskADin2018 with
                                 * default constructor. Based on my understanding
                                 * this is mostly related to input/output
                                 * processes.
                                 */
                                AliAnalysisTaskADin2018();

                                /**
                                 * Create a new AliAnalysisTaskADin2018 with
                                 * a specific name. This serves to instantiate
                                 * the AliAnalysisTaskSE object built with the
                                 * constructor, that will take this specific
                                 * name.
                                 *
                                 * \param name , the name taken by the AliAnalysisTaskSE object.
                                 */
                                AliAnalysisTaskADin2018(const char *name);

                                /**
                                 * Virtual destructor for the class. It will
                                 * take care of all the particular objects to be
                                 * destroyed for the class.
                                 */
        virtual                 ~AliAnalysisTaskADin2018();

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
                                 * Implement the check for the AD.
                                 * 0 = no ADC at all;
                                 * 1 = ADC decision
                                 */
        virtual void            SetCheckAD( Int_t fFlagValue ){ fADcheck = fFlagValue; }

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
                                 * Switch for the AD analysis.
                                 */
        Int_t                   fADcheck;           //  NB: without streamer,
                                                    //      otherwise you can't
                                                    //      set it from outside
                                                    //      in the wagon configs...

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
        TH1F*                   fDimuonPtDistributionRapidityHv3[6];        //!
        TH1F*                   fDimuonPtDistributionRapidityH[6];          //!


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
                                 *
                                 * NEW: v3 uses a variable binning...
                                 */
        TH1F*                   fDimuonPtDistributionZNCzeroZNAzeroH;                              //!
        TH1F*                   fDimuonPtDistributionZNCzeroZNAzeroShiftPlusOneH;                  //!
        TH1F*                   fDimuonPtDistributionZNCzeroZNAzeroHv2;                            //!
        TH1F*                   fDimuonPtDistributionZNCzeroZNAzeroRapidityHv2[3];                 //!
        TH1F*                   fDimuonPtDistributionZNCzeroZNAzeroHv3;                            //!
        TH1F*                   fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3[3];                 //!
        TH1F*                   fDimuonPtDistributionZNCzeroZNAzeroHv2LowerSide;                   //!
        TH1F*                   fDimuonPtDistributionZNCzeroZNAzeroRapidityHv2LowerSide[3];        //!
        TH1F*                   fDimuonPtDistributionZNCzeroZNAzeroHv3LowerSide;                   //!
        TH1F*                   fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3LowerSide[3];        //!
        TH1F*                   fDimuonPtDistributionZNCzeroZNAzeroHv2HigherSide;                  //!
        TH1F*                   fDimuonPtDistributionZNCzeroZNAzeroRapidityHv2HigherSide[3];       //!
        TH1F*                   fDimuonPtDistributionZNCzeroZNAzeroHv3HigherSide;                  //!
        TH1F*                   fDimuonPtDistributionZNCzeroZNAzeroRapidityHv3HigherSide[3];       //!

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
                                 *
                                 * NEW: v3 uses a variable binning...
                                 */
        TH1F*                   fDimuonPtDistributionZNCzeroZNAanyH;                              //!
        TH1F*                   fDimuonPtDistributionZNCzeroZNAanyShiftPlusOneH;                  //!
        TH1F*                   fDimuonPtDistributionZNCzeroZNAanyHv2;                            //!
        TH1F*                   fDimuonPtDistributionZNCzeroZNAanyRapidityHv2[3];                 //!
        TH1F*                   fDimuonPtDistributionZNCzeroZNAanyHv3;                            //!
        TH1F*                   fDimuonPtDistributionZNCzeroZNAanyRapidityHv3[3];                 //!
        TH1F*                   fDimuonPtDistributionZNCzeroZNAanyHv2LowerSide;                   //!
        TH1F*                   fDimuonPtDistributionZNCzeroZNAanyRapidityHv2LowerSide[3];        //!
        TH1F*                   fDimuonPtDistributionZNCzeroZNAanyHv3LowerSide;                   //!
        TH1F*                   fDimuonPtDistributionZNCzeroZNAanyRapidityHv3LowerSide[3];        //!
        TH1F*                   fDimuonPtDistributionZNCzeroZNAanyHv2HigherSide;                  //!
        TH1F*                   fDimuonPtDistributionZNCzeroZNAanyRapidityHv2HigherSide[3];       //!
        TH1F*                   fDimuonPtDistributionZNCzeroZNAanyHv3HigherSide;                  //!
        TH1F*                   fDimuonPtDistributionZNCzeroZNAanyRapidityHv3HigherSide[3];       //!

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
                                 *
                                 * NEW: v3 uses a variable binning...
                                 */
        TH1F*                   fDimuonPtDistributionZNCanyZNAzeroH;                              //!
        TH1F*                   fDimuonPtDistributionZNCanyZNAzeroShiftPlusOneH;                  //!
        TH1F*                   fDimuonPtDistributionZNCanyZNAzeroHv2;                            //!
        TH1F*                   fDimuonPtDistributionZNCanyZNAzeroRapidityHv2[3];                 //!
        TH1F*                   fDimuonPtDistributionZNCanyZNAzeroHv3;                            //!
        TH1F*                   fDimuonPtDistributionZNCanyZNAzeroRapidityHv3[3];                 //!
        TH1F*                   fDimuonPtDistributionZNCanyZNAzeroHv2LowerSide;                   //!
        TH1F*                   fDimuonPtDistributionZNCanyZNAzeroRapidityHv2LowerSide[3];        //!
        TH1F*                   fDimuonPtDistributionZNCanyZNAzeroHv3LowerSide;                   //!
        TH1F*                   fDimuonPtDistributionZNCanyZNAzeroRapidityHv3LowerSide[3];        //!
        TH1F*                   fDimuonPtDistributionZNCanyZNAzeroHv2HigherSide;                  //!
        TH1F*                   fDimuonPtDistributionZNCanyZNAzeroRapidityHv2HigherSide[3];       //!
        TH1F*                   fDimuonPtDistributionZNCanyZNAzeroHv3HigherSide;                  //!
        TH1F*                   fDimuonPtDistributionZNCanyZNAzeroRapidityHv3HigherSide[3];       //!

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
                                 *
                                 * NEW: v3 uses a variable binning...
                                 */
        TH1F*                   fDimuonPtDistributionZNCanyZNAanyH;                              //!
        TH1F*                   fDimuonPtDistributionZNCanyZNAanyShiftPlusOneH;                  //!
        TH1F*                   fDimuonPtDistributionZNCanyZNAanyHv2;                            //!
        TH1F*                   fDimuonPtDistributionZNCanyZNAanyRapidityHv2[3];                 //!
        TH1F*                   fDimuonPtDistributionZNCanyZNAanyHv3;                            //!
        TH1F*                   fDimuonPtDistributionZNCanyZNAanyRapidityHv3[3];                 //!
        TH1F*                   fDimuonPtDistributionZNCanyZNAanyHv2LowerSide;                   //!
        TH1F*                   fDimuonPtDistributionZNCanyZNAanyRapidityHv2LowerSide[3];        //!
        TH1F*                   fDimuonPtDistributionZNCanyZNAanyHv3LowerSide;                   //!
        TH1F*                   fDimuonPtDistributionZNCanyZNAanyRapidityHv3LowerSide[3];        //!
        TH1F*                   fDimuonPtDistributionZNCanyZNAanyHv2HigherSide;                  //!
        TH1F*                   fDimuonPtDistributionZNCanyZNAanyRapidityHv2HigherSide[3];       //!
        TH1F*                   fDimuonPtDistributionZNCanyZNAanyHv3HigherSide;                  //!
        TH1F*                   fDimuonPtDistributionZNCanyZNAanyRapidityHv3HigherSide[3];       //!

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
        // - AD charge / multiplicity
        /* -
         * - Checks needed to make sure
         * - that the bkg is under control.
         */
                                /**
                                 * Collected multiplicity in the AD
                                 * per channel and neutron emission
                                 * class, plus total.
                                 */
        TH1F*                   fADmultiplicityH[16];                  //!
        TH1F*                   fADmultiplicity0N0NclassH[16];         //!
        TH1F*                   fADmultiplicity0NXNclassH[16];         //!
        TH1F*                   fADmultiplicityXN0NclassH[16];         //!
        TH1F*                   fADmultiplicityXNXNclassH[16];         //!
        TH1F*                   fADmultiplicity0N0NclassRapidityH[48]; //!
        TH1F*                   fADmultiplicity0NXNclassRapidityH[48]; //!
        TH1F*                   fADmultiplicityXN0NclassRapidityH[48]; //!
        TH1F*                   fADmultiplicityXNXNclassRapidityH[48]; //!
        TH1F*                   fADmultiplicityTotalH;                 //!
        TH1F*                   fADmultiplicity0N0NclassTotalH;        //!
        TH1F*                   fADmultiplicity0NXNclassTotalH;        //!
        TH1F*                   fADmultiplicityXN0NclassTotalH;        //!
        TH1F*                   fADmultiplicityXNXNclassTotalH;        //!
        TH1F*                   fADAmultiplicityTotalH;                //!
        TH1F*                   fADAmultiplicity0N0NclassTotalH;       //!
        TH1F*                   fADAmultiplicity0NXNclassTotalH;       //!
        TH1F*                   fADAmultiplicityXN0NclassTotalH;       //!
        TH1F*                   fADAmultiplicityXNXNclassTotalH;       //!
        TH1F*                   fADCmultiplicityTotalH;                //!
        TH1F*                   fADCmultiplicity0N0NclassTotalH;       //!
        TH1F*                   fADCmultiplicity0NXNclassTotalH;       //!
        TH1F*                   fADCmultiplicityXN0NclassTotalH;       //!
        TH1F*                   fADCmultiplicityXNXNclassTotalH;       //!
        TH2F*                   fADCmultiplicityTotalVsZNCenergyH;     //!
        TH2F*                   fADAmultiplicityTotalVsZNAenergyH;     //!

        TH2F*                   fADCmultiplicityTotalVsZNCenergyH_ADCno;     //!
        TH2F*                   fADAmultiplicityTotalVsZNAenergyH_ADAno;     //!


        //_______________________________
        // - VZERO charge / multiplicity
        /* -
         * - Checks needed to make sure
         * - that the bkg is under control.
         */
                                /**
                                 * Collected multiplicity in the AD
                                 * per channel and neutron emission
                                 * class, plus total.
                                 */
        TH1F*                   fVZEROmultiplicityTotalH;                 //!
        TH1F*                   fVZEROmultiplicity0N0NclassTotalH;        //!
        TH1F*                   fVZEROmultiplicity0NXNclassTotalH;        //!
        TH1F*                   fVZEROmultiplicityXN0NclassTotalH;        //!
        TH1F*                   fVZEROmultiplicityXNXNclassTotalH;        //!
        TH1F*                   fVZEROAmultiplicityTotalH;                //!
        TH1F*                   fVZEROAmultiplicity0N0NclassTotalH;       //!
        TH1F*                   fVZEROAmultiplicity0NXNclassTotalH;       //!
        TH1F*                   fVZEROAmultiplicityXN0NclassTotalH;       //!
        TH1F*                   fVZEROAmultiplicityXNXNclassTotalH;       //!
        TH1F*                   fVZEROCmultiplicityTotalH;                //!
        TH1F*                   fVZEROCmultiplicity0N0NclassTotalH;       //!
        TH1F*                   fVZEROCmultiplicity0NXNclassTotalH;       //!
        TH1F*                   fVZEROCmultiplicityXN0NclassTotalH;       //!
        TH1F*                   fVZEROCmultiplicityXNXNclassTotalH;       //!
        TH2F*                   fVZEROCmultiplicityTotalVsZNCenergyH;     //!
        TH2F*                   fVZEROAmultiplicityTotalVsZNAenergyH;     //!

        TH2F*                   fADCmultiplicityVsVZEROCmultiplicityH;    //!
        TH2F*                   fADAmultiplicityVsVZEROAmultiplicityH;    //!

        TH2F*                   fVZEROCmultiplicityTotalVsZNCenergyH_VZEROCno;     //!
        TH2F*                   fVZEROAmultiplicityTotalVsZNAenergyH_VZEROAno;     //!

        TH2F*                   fADCmultiplicityVsVZEROCmultiplicityH_VZEROCno;    //!
        TH2F*                   fADCmultiplicityVsVZEROCmultiplicityH_ADCno_VZEROCno;    //!
        TH2F*                   fADAmultiplicityVsVZEROAmultiplicityH_VZEROAno;    //!
        TH2F*                   fADAmultiplicityVsVZEROAmultiplicityH_ADAno_VZEROAno;    //!


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
        AliAnalysisTaskADin2018(const AliAnalysisTaskADin2018&);

        /**
         * Not implemented yet...
         */
        AliAnalysisTaskADin2018& operator=(const AliAnalysisTaskADin2018&);


        /**
         * This is important for ROOT only. I do not remember the reason anymore.
         * If I happen to encounter it again in the future, I will make sure to
         * record it!
         */
        ClassDef(AliAnalysisTaskADin2018, 12);
};

#endif
