/**
* @Author: Pascal Dillenseger <pascaldillenseger>
* @Date:   2016-11-11, 11:35:56
* @Email:  pdillens@cern.ch
* @Last modified by:   pascaldillenseger
* @Last modified time: 2016-11-11, 17:38:34
*/



#ifndef ALIDIELECTRONQNEPCORRECTION_H
#define ALIDIELECTRONQNEPCORRECTION_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//#############################################################
//#                                                           #
//#         Class AliDielectronQnEPcorrection                 #
//#         Provides needed function for the 2016 est.        #
//#         QnFramework and can be used to set the cuts       #
//#         for the auto-correlation removal in the           #
//#         AliDielectronVarManager                           #
//#                                                           #
//#  Authors:                                                 #
//#   Anton     Andronic, GSI / A.Andronic@gsi.de             #
//#   Ionut C.  Arsene,   GSI / I.C.Arsene@gsi.de             #
//#   Julian    Book,     Uni Ffm / Julian.Book@cern.ch       #
//#   Frederick Kramer,   Uni Ffm, / Frederick.Kramer@cern.ch #
//#   Magnus    Mager,    CERN / Magnus.Mager@cern.ch         #
//#   WooJin J. Park,     GSI / W.J.Park@gsi.de               #
//#   Jens      Wiechula, Uni HD / Jens.Wiechula@cern.ch      #
//#                                                           #
//#############################################################

#include <Rtypes.h>
#include <TBits.h>
#include <TList.h>

#include <AliAnalysisCuts.h>
#include <AliVTrack.h>
#include <AliESDtrack.h>
#include <AliAODTrack.h>

#include "AliQnCorrectionsQnVector.h"
#include "AliDielectronPair.h"




class AliDielectronQnEPcorrection : public TNamed{
public:
  enum EPairCutActive{
    kPairMass = 0,
    kPairPt,
    kPairCutActiveMax
  };
  enum ECutValues{
    kLow = 0,
    kHigh
  };

  AliDielectronQnEPcorrection();
  virtual ~AliDielectronQnEPcorrection();

  // setters
  void SetActivatePairCut(EPairCutActive activateCut) {fPairCutActive[activateCut] = kTRUE;}
  void SetActivatePairCut(EPairCutActive activateCut, Double_t min, Double_t max) {fPairCutActive[activateCut] = kTRUE; fPairCutValues[activateCut][kLow] = min; fPairCutValues[activateCut][kHigh] = max;}
  // getters
  static AliQnCorrectionsQnVector* GetQnVectorFromList( TList *list,            const char *subdetector,
                      const char *expectedstep,
                      const char *altstep);
  Double_t GetACcorrectedQnTPCEventplane(const AliDielectronPair *pair, TList *qnlist);
  Bool_t IsSelected(const AliDielectronPair *pair);
  //
  //Analysis cuts interface
  //

  //
  // Cut information
  //





 private:

  Bool_t IsSelected(const AliVTrack *track);
  Bool_t IsSelectedESDtrack(const AliESDtrack *track); // NOTE Not implemented
  Bool_t IsSelectedAODtrack(const AliAODTrack *track);

  Bool_t fPairCutActive[kPairCutActiveMax];
  Double_t fPairCutValues[kPairCutActiveMax][kHigh+1];
  ClassDef(AliDielectronQnEPcorrection,1)
};


//
//Inline functions
//

#endif
