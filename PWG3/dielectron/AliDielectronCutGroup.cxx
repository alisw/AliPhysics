/*************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

//////////////////////////////////////////////////////////////////////////
//                           CutGroup                                   //
//                                                                      //
//                                                                      //
//   Allow to define groups of cut conditions which are tested with     //
//      an OR condition between groups and an AND within groups         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "AliDielectronCutGroup.h"

ClassImp(AliDielectronCutGroup)

AliDielectronCutGroup::AliDielectronCutGroup(Bool_t compOperator /*=kCompOR*/) :
  AliAnalysisCuts(),
  fCutGroupList(0x0),
  fCompOperator(compOperator)
{
  //
  // Default constructor
  //
}

//_____________________________________________________________________
AliDielectronCutGroup::AliDielectronCutGroup(const char* name, const char* title, Bool_t compOperator /*=kCompOR*/) :
  AliAnalysisCuts(name, title),
  fCutGroupList(0x0),
  fCompOperator(compOperator)
{
  //
  // Named Constructor
  //
}

//_____________________________________________________________________
AliDielectronCutGroup::~AliDielectronCutGroup() 
{
  //
  //Default Destructor
  //
}

//_____________________________________________________________________
Bool_t AliDielectronCutGroup::IsSelected(TObject* track) 
{
  //
  // Selection-finder handling different comparison operations
  //
  
  
  //Different init for and/or makes code shorter
  Bool_t selectionResult=fCompOperator;
  
  TIter listIterator(&fCutGroupList);
  while (AliAnalysisCuts *thisCut = (AliAnalysisCuts*) listIterator()) {
    if (fCompOperator == kCompOR) {
      selectionResult = (selectionResult || thisCut->IsSelected(track));
    }
    else { //kCompAND
      selectionResult = (selectionResult && thisCut->IsSelected(track));
      if (selectionResult==kFALSE) break; //Save loops vs. additional check?
    }
    
  }
  return selectionResult;
}

//_____________________________________________________________________

void AliDielectronCutGroup::AddCut(AliAnalysisCuts* fCut) 
{
  //
  // Add a defined cut to the list
  //
  
  fCutGroupList.Add(fCut);
}

//_____________________________________________________________________
void AliDielectronCutGroup::SetCompOperator(Bool_t compOperator) 
{
  //
  // Switch between AND/OR
  //
  
  fCompOperator = compOperator;
}
