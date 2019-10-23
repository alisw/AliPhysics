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
        Int_t                   fADcheck;           //!

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
        ClassDef(AliAnalysisTaskADin2018, 1);
};

#endif
