/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskCTrue_H
#define AliAnalysisTaskCTrue_H

/**
 * @file   AliAnalysisTaskCTrue.h
 * @author Simone Ragoni <Simone.Ragoni@cern.ch>
 * @date   April 2019
 */


#include "AliAnalysisTaskSE.h"

class TTree;
class TH1I;
// class TGraphErrors;
class TList;
class AliAODEvent;


typedef std::map< Int_t, std::pair<Double_t, Double_t> > AliRunWithMuAndWeight;

/**
 * \file AliAnalysisTaskCTrue.h
 * \brief Contains the declaration of the AliAnalysisTaskCTrue class
 */

/**
 * \class AliAnalysisTaskCTrue
 * \brief This analysis task is entirely derived from the files that were
 *        passed onto me for the CTRUE analysis. hopefully I should implement
 *        the full analysis for all the CTRUE classes, classes being:
 *        -  B : collision conditions;
 *        -  E : empty;
 *        -  A : beam on A side;
 *        -  C : beam on C side.
 */
class AliAnalysisTaskCTrue : public AliAnalysisTaskSE
{
public:
                                /**
                                 * Create a new AliAnalysisTaskCTrue with
                                 * default constructor. Based on my understanding
                                 * this is mostly related to input/output
                                 * processes.
                                 */
                                AliAnalysisTaskCTrue();

                                /**
                                 * Create a new AliAnalysisTaskCTrue with
                                 * a specific name. This serves to instantiate
                                 * the AliAnalysisTaskSE object built with the
                                 * constructor, that will take this specific
                                 * name.
                                 *
                                 * \param name , the name taken by the AliAnalysisTaskSE object.
                                 */
                                AliAnalysisTaskCTrue(const char *name);

                                /**
                                 * Virtual destructor for the class. It will
                                 * take care of all the particular objects to be
                                 * destroyed for the class.
                                 */
    virtual                     ~AliAnalysisTaskCTrue();

                                /**
                                 * The function related to the instantiation of
                                 * all the histograms and the output list.
                                 */
    virtual void                UserCreateOutputObjects();

                                /**
                                 * The full analysis happens here.
                                 * What will happen is that the code will
                                 * divide in the many different cases the
                                 * results, so B/E/A/C. At least, this
                                 * should be its final incarnation.
                                 *
                                 * \param option , actually it is not used for now...
                                 */
    virtual void                UserExec(Option_t* option);

                                /**
                                 * This will fill the map containing the good
                                 * run numbers. For now this function will be
                                 * inside the constructor of the class.
                                 */
    void                        FillGoodRunMapInfo(AliRunWithMuAndWeight &fMapGoodRunsToMuAndWeight);

                                /**
                                 * Set trigger inputs for 2018 Pb--Pb data...
                                 */
    void                        Set2018PbPb();

                                /**
                                 * Set trigger inputs for 2015 Pb--Pb data...
                                 */
    void                        Set2015PbPb();

                                /**
                                 * Set trigger inputs for 2017 Xe--Xe data...
                                 */
    void                        SetXeXe();

    //                             /**
    //                              * Function computing the binomial error.
    //                              * It has two inputs.
    //                              *
    //                              * \param cut , the threshold value
    //                              * \param CTRUE, the number of CTRUE events
    //                              *
    //                              * If I understand well enough the efficiency
    //                              * would be something like cut/CTRUE.
    //                              * Hence the binomial error would be
    //                              * something similar to the
    //                              * \sqrt[eff*(1-eff)/CTRUE].
    //                              */
    // Double_t                    BinomialError(Double_t cut, Double_t CTRUE);
    //
    //                             /**
    //                              * This function executes the polinomial
    //                              * fit to the wanted distribution.
    //                              */
    // Double_t                    FitPolinomial(Double_t *x, Double_t *p);
    //
    //                             /**
    //                              * Compute the total efficiency.
    //                              */
    // void                        ComputeEfficiency(  Int_t n,
    //                                                 Double_t *weight,
    //                                                 Double_t *mu,
    // 	                                              Double_t p0,
    //                                                 Double_t p0e,
    //                                                 Double_t p1,
    //                                                 Double_t p1e,
    // 	                                              Double_t *eff
    //                                                 // std::vector<Double_t>* fVectorEfficiencyC
    //                                                 );
    //
    //                             /**
    //                              * Normalize the total efficiency and
    //                              * retrieve the value with the corresponding
    //                              * error.
    //                              */
    // Double_t*                   NormalizeEfficiency();
    //
    //                             /**
    //                              * This function takes as argument 2 TH1F*
    //                              * and one TGraph. Basically it is like a fit.
    //                              * What happens here is that we obtain the final
    //                              * results we wanted for the efficiency.
    //                              * This needs all the other histograms to
    //                              * have already been filled however...
    //                              * This is the reason why it has to be called
    //                              * at the end of the analysis, inside the
    //                              * Terminate() function!!
    //                              *
    //                              * \param signal: the signal histo
    //                              * \param bkg:    the bkg    histo
    //                              * \param graphToBeFilled: self-explanatory
    //                              */
    // void                        DoPlot( TH1F*         signal,
    //                                     TH1F*         background,
    //                                     TGraphErrors* graphToBeFilled
    //                                     );

                                /**
                                 * Called at the END of the analysis (when all
                                 * events are processed). But it is not actually
                                 * doing anything! I guess it is mostly needed
                                 * for I/O purposes and GRID interfacing...
                                 * Infact, in the LEGO trains it is easier
                                 * to set up a wagon if this doesn't require
                                 * histograms to be prepared in the terminate
                                 * function. So it is best to go with empty or
                                 * almost empty Terminate() function...
                                 */
    virtual void                Terminate(Option_t* option);

private:

                                /**
                                 * Input event for the analysis.
                                 */
    AliAODEvent*                fAOD;           //!

                                /**
                                 * The output list containing all the histograms
                                 * required for the analysis. In a second time
                                 * I will probably make it so to include every
                                 * possible cut variation to better compute the
                                 * systematics.
                                 */
    TList*                      fOutputList;    //!

                                /**
                                 * The analysis tree created by the original
                                 * version present in the UPC directory.
                                 * It should be still kept to check for
                                 * consistency with the original code.
                                 */
    TTree*                      fAnaTree;       //!

                                /**
                                 * Histogram containing the counter information.
                                 */
    TH1I*                       fCounterH;       //!

                                /**
                                 * This histogram is my tentative to instantiate
                                 * histograms inside the terminate function.
                                 * This could be useful when I want to do the
                                 * fits inside the AliRoot code!!!
                                 */
    TH1I*                       fCounterTryingH;       //!

                                /**
                                 * This is the IMPORTANT flag needed to
                                 * switch between the many different types of
                                 * CTRUE classes.
                                 */
    Int_t                       fCtrue;           //!

                                /**
                                 * Same as above, for other trigger classes.
                                 */
    Int_t                       fC1zed;           //!

                                /**
                                 * This is the map containing the GOOD RunNumbers.
                                 * Furthermore, I am using a map this time in
                                 * order to retrieve the information
                                 * specific to that run...
                                 * Basically, the format I would like to apply
                                 * for now is that I use the run number as
                                 * the key to access a Double_t[2] array.
                                 * The basic idea is that.
                                 * What I will do is to create instead a
                                 * std::map< Int_t, std::pair<Double_t, Double_t> >.
                                 * This way I am accessing the element
                                 * corresponding to the right run number,
                                 * and then using the right information
                                 * for the weight and the mu value.
                                 */
    AliRunWithMuAndWeight       fMapGoodRunsToMuAndWeight;

                                /**
                                 * This histogram records the number of entries
                                 * against the respective run number, but the
                                 * binning is done properly, like normally
                                 * found for this type of histograms.
                                 * What is different is the usage of proper
                                 * bin labels...
                                 */
    TH1F*                       fEntriesAgainstRunNumberProperlyH;         //!

    /*_________________________________________________________________________
    ***************************************************************************
    *                                                                         *
    *                                                                         *
    *                             CTRUE == 1 (B)                              *
    *                                                                         *
    *                                                                         *
    ***************************************************************************
    __________________________________________________________________________*/

                                /**
                                 * This histogram records the number of entries
                                 * against the respective run number, but the
                                 * binning is done properly, like normally
                                 * found for this type of histograms.
                                 * What is different is the usage of proper
                                 * bin labels...
                                 *
                                 * Counting the number of CTRUE==1 events.
                                 */
    TH1F*                       fCTrueBEventsPerRunNumberH;         //!

                                /**
                                 * This histogram records the number of entries
                                 * against the respective run number, but the
                                 * binning is done properly, like normally
                                 * found for this type of histograms.
                                 * What is different is the usage of proper
                                 * bin labels...
                                 *
                                 * Counting the number of CTRUE==1 events
                                 * with the following conditions:
                                 * !f0UBA
                                 * !f0UBC
                                 * fADCDecision == 0
                                 * fADCDecision == 0
                                 * !f0VBA
                                 * fV0ADecision == 0
                                 * !f0VBC
                                 * fV0CDecision == 0
                                 */
    TH1F*                       fCTrueBEventsPerRunNumberConditionsH;         //!

                                /**
                                 * This histogram records the number of entries
                                 * against the respective run number, but the
                                 * binning is done properly, like normally
                                 * found for this type of histograms.
                                 * What is different is the usage of proper
                                 * bin labels...
                                 *
                                 * Recording the number of tracklets per
                                 * run number.
                                 */
    TH1F*                       fTrackletsPerRunNumberH;         //!

    //_______________________________
    // - V0A PLOTS
                                /**
                                 * This histogram records the number of entries
                                 * against the respective run number, but the
                                 * binning is done properly, like normally
                                 * found for this type of histograms.
                                 * What is different is the usage of proper
                                 * bin labels...
                                 *
                                 * Recording the number of VBA firing per
                                 * run number.
                                 */
    TH1F*                       fRecurringVetoAH;         //!

                                /**
                                 * This histogram records the number of entries
                                 * against the respective run number, but the
                                 * binning is done properly, like normally
                                 * found for this type of histograms.
                                 * What is different is the usage of proper
                                 * bin labels...
                                 *
                                 * Recording the number of VBA firing per
                                 * run number.
                                 */
    TH1F*                       fRecurringVetoAtrackletsH;         //!

                                /**
                                 * This histogram records the number of entries
                                 * against the respective run number, but the
                                 * binning is done properly, like normally
                                 * found for this type of histograms.
                                 * What is different is the usage of proper
                                 * bin labels...
                                 *
                                 * Recording the number of VBA firing per
                                 * run number.
                                 */
    TH1F*                       fVBAforRunNumberH;         //!

                                /**
                                 * This histogram records the number of entries
                                 * against the respective run number, but the
                                 * binning is done properly, like normally
                                 * found for this type of histograms.
                                 * What is different is the usage of proper
                                 * bin labels...
                                 *
                                 * Recording the number of VBA firing
                                 * tracklets per
                                 * run number.
                                 */
    TH1F*                       fVBATrackletsForRunNumberH;         //!

                                /**
                                 * This histogram records the number of entries
                                 * against the respective run number, but the
                                 * binning is done properly, like normally
                                 * found for this type of histograms.
                                 * What is different is the usage of proper
                                 * bin labels...
                                 *
                                 * Recording the number of VDA firing per
                                 * run number.
                                 */
    TH1F*                       fVDAforRunNumberH;         //!

                                /**
                                 * This histogram records the number of entries
                                 * against the respective run number, but the
                                 * binning is done properly, like normally
                                 * found for this type of histograms.
                                 * What is different is the usage of proper
                                 * bin labels...
                                 *
                                 * Recording the number of VDA firing
                                 * tracklets per
                                 * run number.
                                 */
    TH1F*                       fVDATrackletsForRunNumberH;         //!

                                /**
                                 * This histogram records the number of entries
                                 * against the respective run number, but the
                                 * binning is done properly, like normally
                                 * found for this type of histograms.
                                 * What is different is the usage of proper
                                 * bin labels...
                                 *
                                 * Recording the number of VBA and VDA
                                 * firing per
                                 * run number.
                                 */
    TH1F*                       fVBAandVDAforRunNumberH;         //!

                                /**
                                 * This histogram records the number of entries
                                 * against the respective run number, but the
                                 * binning is done properly, like normally
                                 * found for this type of histograms.
                                 * What is different is the usage of proper
                                 * bin labels...
                                 *
                                 * Recording the number of VBA and VDA
                                 * firing tracklets per
                                 * run number.
                                 */
    TH1F*                       fVBAandVDATrackletsForRunNumberH;         //!

    //_______________________________
    // - V0C PLOTS
                                /**
                                 * This histogram records the number of entries
                                 * against the respective run number, but the
                                 * binning is done properly, like normally
                                 * found for this type of histograms.
                                 * What is different is the usage of proper
                                 * bin labels...
                                 *
                                 * Recording the number of VBA firing per
                                 * run number.
                                 */
    TH1F*                       fRecurringVetoCH;         //!

                                /**
                                 * This histogram records the number of entries
                                 * against the respective run number, but the
                                 * binning is done properly, like normally
                                 * found for this type of histograms.
                                 * What is different is the usage of proper
                                 * bin labels...
                                 *
                                 * Recording the number of VBA firing per
                                 * run number.
                                 */
    TH1F*                       fRecurringVetoCtrackletsH;         //!

                                /**
                                 * This histogram records the number of entries
                                 * against the respective run number, but the
                                 * binning is done properly, like normally
                                 * found for this type of histograms.
                                 * What is different is the usage of proper
                                 * bin labels...
                                 *
                                 * Recording the number of VBC firing per
                                 * run number.
                                 */
    TH1F*                       fVBCforRunNumberH;         //!

                                /**
                                 * This histogram records the number of entries
                                 * against the respective run number, but the
                                 * binning is done properly, like normally
                                 * found for this type of histograms.
                                 * What is different is the usage of proper
                                 * bin labels...
                                 *
                                 * Recording the number of VBC firing
                                 * tracklets per
                                 * run number.
                                 */
    TH1F*                       fVBCTrackletsForRunNumberH;         //!

                                /**
                                 * This histogram records the number of entries
                                 * against the respective run number, but the
                                 * binning is done properly, like normally
                                 * found for this type of histograms.
                                 *
                                 * Recording the number of VDC firing per
                                 * run number.
                                 */
    TH1F*                       fVDCforRunNumberH;         //!

                                /**
                                 * This histogram records the number of entries
                                 * against the respective run number, but the
                                 * binning is done properly, like normally
                                 * found for this type of histograms.
                                 * What is different is the usage of proper
                                 * bin labels...
                                 *
                                 * Recording the number of VDC firing
                                 * tracklets per
                                 * run number.
                                 */
    TH1F*                       fVDCTrackletsForRunNumberH;         //!

                                /**
                                 * This histogram records the number of entries
                                 * against the respective run number, but the
                                 * binning is done properly, like normally
                                 * found for this type of histograms.
                                 * What is different is the usage of proper
                                 * bin labels...
                                 *
                                 * Recording the number of VBC and VDC
                                 * firing per
                                 * run number.
                                 */
    TH1F*                       fVBCandVDCforRunNumberH;         //!

                                /**
                                 * This histogram records the number of entries
                                 * against the respective run number, but the
                                 * binning is done properly, like normally
                                 * found for this type of histograms.
                                 * What is different is the usage of proper
                                 * bin labels...
                                 *
                                 * Recording the number of VBC and
                                 * VDC firing tracklets per
                                 * run number.
                                 */
    TH1F*                       fVBCandVDCTrackletsForRunNumberH;         //!

                                /**
                                 * This histogram shows the distribution of
                                 * fV0CDecision VS fADCDecision. This should
                                 * behave like some kind of Truth Table, with
                                 * most data points focused on the True-True
                                 * and False-False cells.
                                 */
    TH2F*                       fCTRUEBfV0CDecisionVSfADCDecisionH;

                                /**
                                 * This histogram shows the distribution of
                                 * fV0ADecision VS fADADecision. This should
                                 * behave like some kind of Truth Table, with
                                 * most data points focused on the True-True
                                 * and False-False cells.
                                 */
    TH2F*                       fCTRUEBfV0ADecisionVSfADADecisionH;



    //__________________________________________________________________________
    // RESULTS TGRAPHERRORS
    //                             /**
    //                              * TGraphErrors to be filled in the
    //                              * Terminate() function from which we could
    //                              * possibly retrieve the efficiency...
    //                              *
    //                              * Only tracklets without V0A or V0C requests.
    //                              */
    // TGraphErrors*               fOnlyTrackletsGE;
    //
    //                             /**
    //                              * TGraphErrors to be filled in the
    //                              * Terminate() function from which we could
    //                              * possibly retrieve the efficiency...
    //                              *
    //                              * Only VBA requested.
    //                              */
    // TGraphErrors*               fVbaGE;
    //
    //                             /**
    //                              * TGraphErrors to be filled in the
    //                              * Terminate() function from which we could
    //                              * possibly retrieve the efficiency...
    //                              *
    //                              * Only VBA with tracklets requested.
    //                              */
    // TGraphErrors*               fVbaTrackletsGE;
    //
    //                             /**
    //                              * TGraphErrors to be filled in the
    //                              * Terminate() function from which we could
    //                              * possibly retrieve the efficiency...
    //                              *
    //                              * Only VDA requested.
    //                              */
    // TGraphErrors*               fVdaGE;
    //
    //                             /**
    //                              * TGraphErrors to be filled in the
    //                              * Terminate() function from which we could
    //                              * possibly retrieve the efficiency...
    //                              *
    //                              * Only VDA with tracklets requested.
    //                              */
    // TGraphErrors*               fVdaTrackletsGE;
    //
    //                             /**
    //                              * TGraphErrors to be filled in the
    //                              * Terminate() function from which we could
    //                              * possibly retrieve the efficiency...
    //                              *
    //                              * VBA and VDA requested.
    //                              */
    // TGraphErrors*               fVbaAndVdaGE;
    //
    //                             /**
    //                              * TGraphErrors to be filled in the
    //                              * Terminate() function from which we could
    //                              * possibly retrieve the efficiency...
    //                              *
    //                              * VBA and VDA and tracklets requested.
    //                              */
    // TGraphErrors*               fVbaAndVdaTrackletsGE;
    //
    //                             /**
    //                              * TGraphErrors to be filled in the
    //                              * Terminate() function from which we could
    //                              * possibly retrieve the efficiency...
    //                              *
    //                              * VBC requested.
    //                              */
    // TGraphErrors*               fVbcGE;
    //
    //                             /**
    //                              * TGraphErrors to be filled in the
    //                              * Terminate() function from which we could
    //                              * possibly retrieve the efficiency...
    //                              *
    //                              * VBC and tracklets requested.
    //                              */
    // TGraphErrors*               fVbcTrackletsGE;
    //
    //                             /**
    //                              * TGraphErrors to be filled in the
    //                              * Terminate() function from which we could
    //                              * possibly retrieve the efficiency...
    //                              *
    //                              * VDC requested.
    //                              */
    // TGraphErrors*               fVdcGE;
    //
    //                             /**
    //                              * TGraphErrors to be filled in the
    //                              * Terminate() function from which we could
    //                              * possibly retrieve the efficiency...
    //                              *
    //                              * VDC and tracklets requested.
    //                              */
    // TGraphErrors*               fVdcTrackletsGE;
    //
    //                             /**
    //                              * TGraphErrors to be filled in the
    //                              * Terminate() function from which we could
    //                              * possibly retrieve the efficiency...
    //                              *
    //                              * VBC and VDC requested.
    //                              */
    // TGraphErrors*               fVbcAndVdcGE;
    //
    //                             /**
    //                              * TGraphErrors to be filled in the
    //                              * Terminate() function from which we could
    //                              * possibly retrieve the efficiency...
    //                              *
    //                              * VBC and VDC and tracklets requested.
    //                              */
    // TGraphErrors*               fVbcAndVdcTrackletsGE;


    // _______________________________
    // CUTS
    /**
     * The following is all the possible checks for the event selections
     * and the track selection as well. Enjoy.
     */
    Int_t                       fRunNum;          //!
    Int_t                       fTracklets;       //!

    UInt_t                      fL0inputs;        //!
    UInt_t                      fL1inputs;        //!

    Double_t                    fZem1Energy;      //!
    Double_t                    fZem2Energy;      //!

                                /**
                                 * ZDC stuff:
                                 * Neutron ZDC energy, time, and TDC values.
                                 * Proton  ZDC same.
                                 */
    Double_t                    fZNCEnergy;       //!
    Double_t                    fZNAEnergy;       //!
    Double_t                    fZPCEnergy;       //!
    Double_t                    fZPAEnergy;       //!
    Double_t                    fZNATDC[4];       //!
    Double_t                    fZNCTDC[4];       //!
    Double_t                    fZPATDC[4];       //!
    Double_t                    fZPCTDC[4];       //!
    Double_t                    fZNATime;         //!
    Double_t                    fZNCTime;         //!

                                /**
                                 * V0 stuff:
                                 * V0A decision.
                                 * V0C decision.
                                 */
    Int_t                       fV0ADecision;     //!
    Int_t                       fV0CDecision;     //!

                                /**
                                 * AD stuff:
                                 * ADA decision.
                                 * ADC decision.
                                 */
    Int_t                       fADADecision;     //!
    Int_t                       fADCDecision;     //!

    TBits                       fIR1Map;          //!
    TBits                       fIR2Map;          //!
    UShort_t                    fBCrossNum;       //!


    //_______________________________
    /* - Trigger Input Definitions
       -
     */
    UChar_t inputId_0VBA;     //! // 0VBA: >=1 V0A cell fired in BB timing gate
    UChar_t inputId_0VBC;     //! // 0VBC: >=1 V0C cell fired in BB timing gate
    UChar_t inputId_0UBA;     //! // 0UBA: >=1 ADA cell fired in BB timing gate
    UChar_t inputId_0UBC;     //! // 0UBC: >=1 ADC cell fired in BB timing gate
    UChar_t inputId_0SH1;     //! // 0SH1: FO fired in SPD
    UChar_t inputId_0STG;     //!
    UChar_t inputId_1ZED;     //!
    UChar_t inputId_0MUL;     //!
    UChar_t inputId_0OM2;     //!
    UChar_t inputId_0VOM;     //!

                                /**
                                 * This is the total weight needed for
                                 * computing efficiencies and similar stuff.
                                 */
    Double_t                    TotalWeightForEfficiency;     //!

                                /**
                                 * Vector containing the efficiency
                                 * values. There will be three of them, one
                                 * for each of the efficiencies...
                                 */
    std::vector<Double_t>*      fVectorEfficiency;     //!

    /**
     * Not implemented yet...
     */
    AliAnalysisTaskCTrue(const AliAnalysisTaskCTrue&);

    /**
     * Not implemented yet...
     */
    AliAnalysisTaskCTrue& operator=(const AliAnalysisTaskCTrue&);


    /**
     * This is important for ROOT only. I do not remember the reason anymore.
     * If I happen to encounter it again in the future, I will make sure to
     * record it!
     */
    ClassDef(AliAnalysisTaskCTrue, 6);
};

#endif
