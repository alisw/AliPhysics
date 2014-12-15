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
//////////////////////////////////////////////////////////////////////////////
// This is the base class for ITS detector signal simulations. Data members //
// include are a pointer to the AliITSDetTypeSim clas in order to access    //
// segmentation and response objects                                        // 
// classes. See the detector specific implementations for the propper code. //
//////////////////////////////////////////////////////////////////////////////
#include "TClonesArray.h"

#include "AliITSsimulation.h"
#include "AliITSpList.h"

ClassImp(AliITSsimulation)

//______________________________________________________________________
AliITSsimulation::AliITSsimulation(): TObject(),
fDetType(0),
fpList(0),
fModule(0),
fEvent(0),
fDebug(0){
    // Default constructor
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    a default constructed AliITSsimulation class
}
//______________________________________________________________________
AliITSsimulation::AliITSsimulation(AliITSDetTypeSim *dettyp): TObject(),
fDetType(dettyp),
fpList(0),
fModule(0),
fEvent(0),
fDebug(0){
    // Default constructor
    // Inputs:
    //    AliITSDetTypeSim * : object used to access segmentation and response
    // Outputs:
    //    none.
    // Return:
    //    a default constructed AliITSsimulation class
}
//__________________________________________________________________________
AliITSsimulation::~AliITSsimulation(){
    // destructor
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    none.

    if(fpList){
      delete fpList;
      fpList = 0;
    }
   }
//__________________________________________________________________________
AliITSsimulation::AliITSsimulation(const AliITSsimulation &s) : TObject(s),
fDetType(s.fDetType),
fpList(s.fpList),
fModule(s.fModule),
fEvent(s.fEvent),
fDebug(s.fDebug){
    //     Copy Constructor
    // Inputs:
    //    const AliITSsimulation &s  simulation class to copy from
    // Outputs:
    //    none.
    // Return:
    //    a standard constructed AliITSsimulation class with values the same
    //    as that of s.
 
}

//_________________________________________________________________________
AliITSsimulation&  AliITSsimulation::operator=(const AliITSsimulation &s){
    //    Assignment operator
    // Inputs:
    //    const AliITSsimulation &s  simulation class to copy from
    // Outputs:
    //    none.
    // Return:
    //    a standard constructed AliITSsimulation class with values the same
    //    as that of s.

    if(&s == this) return *this;
    this->fModule       = s.fModule;
    this->fEvent        = s.fEvent;
    this->fpList        = s.fpList;
    return *this;
}
//______________________________________________________________________
Bool_t AliITSsimulation::AddSDigitsToModule(TClonesArray *pItemA,Int_t mask ){
    // Add Summable digits to module maps.
    // Inputs:
    //    TClonesArray *pItemA  Array of AliITSpListItems (SDigits).
    //    Int_t         mask    Track number off set value (see 
    //                          AliITSpList::AddItemTo).
    // Outputs:
    //    none.
    // Return:
    //    kTRUE if there is a signal >0 else kFALSE
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
