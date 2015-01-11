/**************************************************************************
 * Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// Comment describing what this class does needed!

//===========================================================
// Dummy comment, should be replaced by a real one
// comment
// comment
// Simple class for the jt anlyais by Beomkyu Kim and Dongjo Kim
//===========================================================

#include <TMath.h>
#include "AliJJet.h"
#include "AliJJetAnalysis.h"
#include "TClonesArray.h"

AliJJetAnalysis::AliJJetAnalysis():
	 fJetListOfList()
{

}

AliJJetAnalysis::AliJJetAnalysis(const AliJJetAnalysis& ap) :
	fJetListOfList(ap.fJetListOfList)
{

}

AliJJetAnalysis& AliJJetAnalysis::operator = (const AliJJetAnalysis& ap)
{
	// assignment operator

	this->~AliJJetAnalysis();
	new(this) AliJJetAnalysis(ap);
	return *this;
}


AliJJetAnalysis::~AliJJetAnalysis(){
}


void AliJJetAnalysis::Run(){

    int iS1 = 0;
    int iS2 = 3;
    TObjArray * jetfinder1 = (TObjArray*) fJetListOfList[iS1];
    TObjArray * jetfinder2 = (TObjArray*) fJetListOfList[iS2];
    AliJJet *jet1 = NULL;
    AliJJet *jet2 = NULL;
    for (int ijet = 0; ijet<jetfinder1->GetEntriesFast(); ijet++){
        jet1 = dynamic_cast<AliJJet*>( jetfinder1->At(ijet) );
        if (!jet1) continue;
        for (int jjet = 0; jjet<jetfinder2->GetEntriesFast(); jjet++){
            jet2 = dynamic_cast<AliJJet*>( jetfinder2->At(jjet) );
            if (!jet2) continue;
            if (jet2->E() < 5 ) continue;
            if( TMath::Abs(jet1->Eta()-jet2->Eta()) < 0.4 ) { 
                //CompareTwoJets(jet1, jet2, -1000., -1000);
            }
        }
    }

}

void AliJJetAnalysis::CompareTwoJets(AliJJet *jet1, AliJJet *jet2, double & dE , int &dN  ){

    if (!jet1 || !jet2) return;
    double commEsum=0;
    double chargedEsum=0;
    int commN=0;
    for (int icon = 0; icon<jet1->GetConstituents()->GetEntriesFast(); icon++){
        AliJBaseTrack *con1 = jet1->GetConstituent(icon);
        if (!con1) continue;
        chargedEsum = 0;
        for (int jcon = 0; jcon<jet2->GetConstituents()->GetEntriesFast(); jcon++){
            AliJBaseTrack *con2 = jet2->GetConstituent(jcon);
            if (!con2) continue;
            chargedEsum += con2->E();
            if(con1->GetID() == con2->GetID()){ 

                commEsum += con2->E();     
                commN ++;       
            }
            
            
        }
    }
    //deltaE = jet2->E()- commEsum;
    dE = chargedEsum - commEsum;
    dN = jet2->GetConstituents()->GetEntriesFast() - commN;

    

}
