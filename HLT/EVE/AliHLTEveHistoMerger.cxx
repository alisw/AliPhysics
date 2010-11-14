
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        *
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Kalliopi Kanaki <Kalliopi.Kanaki@ift.uib.no>          *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** @file   AliHLTEveHistoMerger.cxx
    @author Kalliopi Kanaki
    @date
    @brief  The Histogram Handler component
*/

#if __GNUC__>= 3
using namespace std;
#endif
#include "AliHLTEveHistoMerger.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TH1.h"
#include "TTimeStamp.h"
#include "TSystem.h"
#include <iostream>

AliHLTEveHistoMerger::AliHLTEveHistoMerger()
        :
        fStore()
{
}

AliHLTEveHistoMerger::~AliHLTEveHistoMerger() {
    // see header file for class documentation
    Clear();
}



void AliHLTEveHistoMerger::Clear()
{
    // reset the store

    for ( unsigned int i=0; i<fStore.size(); i++ ) {
        for ( unsigned int j=0; j<fStore[i].fInstances.size(); j++ ) {
            delete fStore[i].fInstances[j].fObject;
        }
        delete fStore[i].fMergedObject;
    }
    fStore.clear();
}




TObject* AliHLTEveHistoMerger::Process(const TObject * evtData, AliHLTUInt32_t spec)
{
    // see header file for class documentation
    //cout<<"\n\nDoEvent called"<<endl;


    if ( !evtData ) return 0;

    if ( !evtData->InheritsFrom(TH1::Class())
            && !evtData->InheritsFrom(TSeqCollection::Class()) ) return 0;

    std::cout<<"received object "<<evtData->GetName()<<" with id="<< spec << std::endl;

    //search for the base entry, if not exist then create a new entry

    int iColl = -1;
    for ( unsigned int i=0; i<fStore.size(); i++ ) {
        if ( fStore[i].fInstances.size()<1 ) continue;
        TObject * obj = fStore[i].fInstances[0].fObject;
        if ( !obj ) continue;
        if ( TString(obj->GetName()).CompareTo(evtData->GetName())==0) {
            iColl = i;
            break;
        }
    }
    cout<<"Collection found: "<<iColl<<endl;
    if ( iColl<0 ) {
        AliHLTGlobalHCCollection c;
        c.fMergedObject = 0;
        c.fNeedToMerge = 1;
        fStore.push_back(c);
        iColl = fStore.size()-1;
    } else {
        fStore[iColl].fNeedToMerge = 1;
    }

    // search for the specific entry, if not exist then create a new one
    {
        AliHLTGlobalHCCollection &c = fStore[iColl];

        int iSpec=-1;
        for ( unsigned int i=0; i<c.fInstances.size(); i++ ) {
            AliHLTGlobalHCInstance &inst = c.fInstances[i];
            if ( inst.fHLTSpecification == spec ) {
                iSpec = i;
                break;
            }
        }
        cout<<"Instance found:"<<iSpec<<endl;
        if ( iSpec<0 ) {
            AliHLTGlobalHCInstance inst;
            inst.fHLTSpecification = spec;
            inst.fObject = 0;
            c.fInstances.push_back(inst);
            iSpec = c.fInstances.size()-1;
        } else {
            delete c.fInstances[iSpec].fObject;
        }

        c.fInstances[iSpec].fObject = evtData->Clone();

        cout<<"index = "<<iColl<<","<<iSpec<<endl;
    }

// merge histos

    for ( unsigned int jColl = 0; jColl<fStore.size(); jColl++) {
        AliHLTGlobalHCCollection &c = fStore[jColl];
        if ( !c.fNeedToMerge && c.fMergedObject ) continue;
        if ( c.fInstances.size() <1 ) continue;
        delete c.fMergedObject;
        c.fMergedObject = c.fInstances[0].fObject->Clone();
        TList l;
        for ( unsigned int i=1; i<c.fInstances.size(); i++ ) {
            l.Add(c.fInstances[i].fObject);
        }

        if ( c.fMergedObject->InheritsFrom(TH1::Class()) ) {
            TH1 *histo = dynamic_cast<TH1*>(c.fMergedObject);
            if ( histo ) histo->Merge(&l);
        }
        else if ( c.fMergedObject->InheritsFrom(TSeqCollection::Class()) ) {
            TSeqCollection *list = dynamic_cast<TSeqCollection*>(c.fMergedObject);
            if ( list ) list->Merge(&l);
        }
        c.fNeedToMerge = 0;

    }

    return fStore[iColl].fMergedObject;

}

