/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskMatchTriggerForward_H
#define AliAnalysisTaskMatchTriggerForward_H

/**
 * @file   AliAnalysisTaskMatchTriggerForward.h
 * @author Simone Ragoni <Simone.Ragoni@cern.ch>
 * @date   February 2019
 */

#include "AliAnalysisTaskSE.h"
#include "TBits.h"

class AliMuonTrackCuts; 					// Include class for standard muon tack cuts

/**
 * \file AliAnalysisTaskMatchTriggerForward.h
 * \brief Contains the declaration of the AliAnalysisTaskMatchTriggerForward class
 */

/**
 * \class AliAnalysisTaskMatchTriggerForward
 * \brief Contains the tools to every part of my analysis
 */
class AliAnalysisTaskMatchTriggerForward : public AliAnalysisTaskSE
{
    public:
                                /**
                                 * Create a new AliAnalysisTaskMatchTriggerForward with
                                 * default constructor. Based on my understanding
                                 * this is mostly related to input/output
                                 * processes.
                                 */
                                AliAnalysisTaskMatchTriggerForward();

                                /**
                                 * Create a new AliAnalysisTaskMatchTriggerForward with
                                 * a specific name. This serves to instantiate
                                 * the AliAnalysisTaskSE object built with the
                                 * constructor, that will take this specific
                                 * name.
                                 *
                                 * \param name , the name taken by the AliAnalysisTaskSE object.
                                 * \param isMC , 0 if Data, 1 if MC(look at the AddTask please).
                                 */
                                AliAnalysisTaskMatchTriggerForward(const char *name);

                                /**
                                 * Virtual destructor for the class. It will
                                 * take care of all the particular objects to be
                                 * destroyed for the class.
                                 */
        virtual                 ~AliAnalysisTaskMatchTriggerForward();

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
                                 * Substitutes the trigger in MC.
                                 */
        Bool_t                  IsTriggered();


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
                                 * This function is needed to record the MC
                                 * truth information inside the nanoAODs.
                                 * Everything MC related should go in here.
                                 * This should help for a subsequent merging
                                 * of the MC and DATA classes.
                                 */
        void                    ProcessMCParticles(AliMCEvent* fMCEvent);

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
                                 * As far as I understand, it should be the
                                 * pseudorapidity distribution of the many muons.
                                 */
        TH1F*                   fEtaMuonH;          //!

                                /**
                                 * Counter to check everything goes alright.
                                 */
        TH1F*                   fCounterH;          //!

                                /**
                                 *
                                 */
        TH1F*                   fRAbsMuonH;         //!

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


        //_______________________________
        // Efficiency plots.
                                /**
                                 * This histogram shows the entries distribution
                                 * per run of the RECONSTRUCTED level.
                                 * This has to be divided by the corresponding
                                 * fMCEfficiencyPerRunH to extract the
                                 * efficiency on a run-by-run basis.
                                 */
        TH1F*                   fEfficiencyPerRunH;  //!

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
        TH1F*                   fMCEfficiencyPerRunH;  //!

                                /**
                                 * This histogram shows the entries distribution
                                 * per run of the RECONSTRUCTED level.
                                 * This has to be divided by the corresponding
                                 * fMCEfficiencyPerRunH to extract the
                                 * efficiency on a run-by-run basis.
                                 * NOTE: this is with the GetTrigger() request.
                                 */
        TH1F*                   fEfficiencyPerRunWithTriggeringH;  //!

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
        TH2F*                   fDeadZoneEtaVsPhiPerRunH[364];               //!    number of total runs
        TH2F*                   fDeadZoneEtaVsPhiPerRunWithTriggeringH[364]; //!    number of total runs

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
                                 * This array of histograms records the
                                 * energy spectra
                                 * of the neutron ZDC.
                                 * On a run-by-run basis.
                                 */
        TH1F*                   fZNCEnergyPerRunH[364];         //!
        TH1F*                   fZNAEnergyPerRunH[364];         //!

                                /**
                                 * This histogram records the pt-ditribution
                                 * for the single muons.
                                 */
        TH1F*                   fSingleMuonPtDistributionH;         //!


        //_______________________________
        // CUTS
        /*
         * The following is all the possible checks for the event selections
         * and the track selection as well. Enjoy.
         */
        Int_t                   fRunNum;        //!
        Int_t                   fTracklets;     //!
        Double_t                fLumiPerRun;    //!

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
        AliAnalysisTaskMatchTriggerForward(const AliAnalysisTaskMatchTriggerForward&);

        /**
         * Not implemented yet...
         */
        AliAnalysisTaskMatchTriggerForward& operator=(const AliAnalysisTaskMatchTriggerForward&);

        /**
         * This is important for ROOT only. I do not remember the reason anymore.
         * If I happen to encounter it again in the future, I will make sure to
         * record it!
         */
        ClassDef(AliAnalysisTaskMatchTriggerForward, 5);
};

#endif
