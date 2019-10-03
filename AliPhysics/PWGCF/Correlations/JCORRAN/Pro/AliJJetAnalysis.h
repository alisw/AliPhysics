/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

// Short comment describing what this class does needed!

//===========================================================
// Dummy comment, should be replaced by a real one
//===========================================================

#ifndef ALIJJETANALYSIS_H
#define ALIJJETANALYSIS_H

#include <TObjArray.h>
#include "AliJJet.h"


class AliJJetAnalysis{
    public:


        AliJJetAnalysis();
        AliJJetAnalysis(const AliJJetAnalysis& ap);
        AliJJetAnalysis& operator = (const AliJJetAnalysis& ap);
        ~AliJJetAnalysis();


        void AddJets(TObjArray * jets ){
            if( !jets ) {
                //return;
            }
            fJetListOfList.Add( (TObject*)jets );
            //if( !jets ) return;
            for( int i=0;i<jets->GetEntriesFast();i++ ){
                //((AliJJet*)jets->At(i))->ReSum();
            }
        } 

        void Run();
        void CompareTwoJets(AliJJet *j1, AliJJet *j2, double & dE ,int & dN);



    private:
        TObjArray fJetListOfList; // !comment needed
};

#endif


