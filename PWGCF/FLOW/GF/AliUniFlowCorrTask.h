/* Copyright(c) 2019, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef UNIFLOWCORRTASK_H
#define UNIFLOWCORRTASK_H

class AliUniFlowCorrTask : public TObject
{
    public:
        AliUniFlowCorrTask(); // default ctor
        AliUniFlowCorrTask( // actual ctor
            Bool_t doRFPs,
            Bool_t doPOIs,
            std::vector<Int_t> harms,
            std::vector<Double_t> gaps = std::vector<Double_t>());
        virtual ~AliUniFlowCorrTask() { fiHarm.clear(); fdGaps.clear(); }

        Bool_t      HasGap() const { return (Bool_t) fiNumGaps; }; // check if Gap
        void        PrintTask() const; // print AliUniFlowCorrTask properties

        Bool_t                fbDoRefs; // which particles are procesed (RFPs / POIs / both )
        Bool_t                fbDoPOIs; // which particles are procesed (RFPs / POIs / both )
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

    ClassDef(AliUniFlowCorrTask, 2);
};


#endif
