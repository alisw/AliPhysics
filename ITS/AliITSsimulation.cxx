/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

#include "TClonesArray.h"

#include "AliITSsimulation.h"
#include "AliITSpList.h"

ClassImp(AliITSsimulation)

AliITSsimulation::AliITSsimulation(){
    // constructor
    fSegmentation = 0;
    fResponse     = 0;
    fpList        = 0;
}
//__________________________________________________________________________
AliITSsimulation::~AliITSsimulation(){
    // destructor
    fSegmentation = 0; // local copies of pointer, do not delete
    fResponse     = 0; // local copies of pointer, do not delete
    delete fpList;
}
//__________________________________________________________________________
AliITSsimulation::AliITSsimulation(const AliITSsimulation &source) : 
    TObject(source){
    //     Copy Constructor 
 
    if(&source == this) return;
    this->fResponse     = source.fResponse;
    this->fSegmentation = source.fSegmentation;
    this->fModule       = source.fModule;
    this->fEvent        = source.fEvent;
    this->fpList        = source.fpList;
    return;
}

//_________________________________________________________________________
AliITSsimulation&  AliITSsimulation::operator=(const AliITSsimulation &source){
    //    Assignment operator

    if(&source == this) return *this;
    this->fResponse     = source.fResponse; 
    this->fSegmentation = source.fSegmentation;
    this->fModule       = source.fModule;
    this->fEvent        = source.fEvent;
    this->fpList        = source.fpList;
    return *this;
}
//______________________________________________________________________
Bool_t AliITSsimulation::AddSDigitsToModule(TClonesArray *pItemA,Int_t mask ){
    // Add Summable digits to module maps.
    Int_t nItems = pItemA->GetEntries();
    Bool_t sig = kFALSE;
 
    // cout << "Adding "<< nItems <<" SDigits to module " << fModule << endl;
    for( Int_t i=0; i<nItems; i++ ) {
        AliITSpListItem * pItem = (AliITSpListItem *)(pItemA->At( i ));
        if( pItem->GetModule() != fModule ) {
            Error( "AddSDigitsToModule","Error reading, SDigits module %d "
		   "!= current module %d: exit",
		   pItem->GetModule(), fModule );
            return sig;
        } // end if
	if(pItem->GetSignal()>0.0 ) sig = kTRUE;
        fpList->AddItemTo( mask, pItem );
    } // end for i
    return sig;
}
