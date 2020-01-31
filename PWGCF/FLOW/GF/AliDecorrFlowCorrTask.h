/* Copyright(c) 2019, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef UNIFLOWCORRTASK_H
#define UNIFLOWCORRTASK_H

class AliDecorrFlowCorrTask : public TObject
{
    public:
        AliDecorrFlowCorrTask(); // default ctor
        AliDecorrFlowCorrTask( // actual ctor
            Bool_t doRef,
            Bool_t doDiff,
            Bool_t doPtA,
            Bool_t doPtRef,
            Bool_t doPtB,
            std::vector<Int_t> harms,
            std::vector<Double_t> gaps = std::vector<Double_t>());
        virtual ~AliDecorrFlowCorrTask() { fiHarm.clear(); fdGaps.clear(); }

        Bool_t      HasGap() const { return (Bool_t) fiNumGaps; }; // check if Gap
        void        PrintTask() const; // print AliDecorrFlowCorrTask properties

        Bool_t                fbDoRef; // which particles are procesed (RFPs / POIs / both )
        Bool_t                fbDoDiff; // which particles are procesed (RFPs / POIs / both )
        Bool_t                fbDoPtA;  //Do Reference particles in pt bins
        Bool_t                fbDoPtRef;  //Single-differential
        Bool_t                fbDoPtB;  //Double-Differential
        Int_t                 fiNumHarm; // correlation order <M>
        Int_t                 fiNumGaps; // number of subevents
        Int_t                 fMaxWeightPower; //max power in q vector filling
        Int_t                 fMaxHarm; //max harm in q vector filling
        TString               fsName; // automatically generated name: see Init() for format
        TString               fsLabel; // automatically generated label see Init() for format
        std::vector<Int_t>    fiHarm; // harmonics n1,n2,...,nM
        std::vector<Double_t> fdGaps; // gaps between subevents (standard GF notation)
    protected:
    private:

    ClassDef(AliDecorrFlowCorrTask, 2);
};


#endif
