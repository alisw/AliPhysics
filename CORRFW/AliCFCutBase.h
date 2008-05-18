#ifndef ALICFCUTBASE_H
#define ALICFCUTBASE_H
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
// Base class for selecton classes for the correction framework 
// Inherits from AliAnalysisCuts. It includes additional methods to handle QA 
// histograms and if needed, study the cut statistics & correlations 
// Author S.Arcelli
// silvia.Arcelli@cern.ch

#include <AliAnalysisCuts.h>
class TBits;
class TList;
//___________________________________________________________________________
class AliCFCutBase : public AliAnalysisCuts
{
 public:
  AliCFCutBase(); //default ctor
  AliCFCutBase(const char* name, const char* title); //ctor
  AliCFCutBase(const AliCFCutBase& obj); //copy ctor  
  virtual ~AliCFCutBase() {;} //dtor
  virtual Bool_t IsQAOn() const {return fIsQAOn;}; //QA flag getter
  virtual void SetQAOn(TList* list) {fIsQAOn=kTRUE; AddQAHistograms(list);} //QA flag setter
  virtual void SetEvtInfo(TObject *) {;}; //Pass pointer to event-level info
  
 protected:
  Bool_t fIsQAOn;//qa checking on/off
  virtual void AddQAHistograms(TList*) {;}; //QA Histos

  ClassDef(AliCFCutBase, 1); // Base class for Correction Framework Cuts
};
 
#endif
