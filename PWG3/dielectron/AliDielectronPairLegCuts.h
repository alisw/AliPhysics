#ifndef ALIDIELECTRONPAIRLEGCUTS_H
#define ALIDIELECTRONPAIRLEGCUTS_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//#############################################################
//#                                                           # 
//#         Class AliDielectronPairLegCuts                    #
//#         Manage Cuts on the legs of the pair               #
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

#include <AliAnalysisFilter.h>

#include <AliAnalysisCuts.h>

class AliDielectronPairLegCuts : public AliAnalysisCuts {
public:
  enum CutType { kBothLegs=0, kAnyLeg, kMixLegs };

  AliDielectronPairLegCuts();
  AliDielectronPairLegCuts(const char* name, const char* title);
  virtual ~AliDielectronPairLegCuts() {;}
  //TODO: make copy constructor and assignment operator public
  //      and implement them
  
  //
  //AliAnalysisCuts interface
  //
  virtual Bool_t IsSelected(TObject* track);
  virtual Bool_t IsSelected(TList*   /* list */ ) {return kFALSE;}
//   virtual Long64_t Merge(TCollection* /* list */)      { return 0; }

  AliAnalysisFilter& GetLeg1Filter() { return fFilterLeg1; }
  AliAnalysisFilter& GetLeg2Filter() { return fFilterLeg2; }

  void SetCutType(CutType type) {fCutType=type;}
private:
  AliAnalysisFilter fFilterLeg1;     // Analysis Filter for leg1
  AliAnalysisFilter fFilterLeg2;     // Analysis Filter for leg2

  CutType fCutType;                  // Type of the cut

  AliDielectronPairLegCuts(const AliDielectronPairLegCuts &c);
  AliDielectronPairLegCuts &operator=(const AliDielectronPairLegCuts &c);
  
  ClassDef(AliDielectronPairLegCuts,1)         //Cut class providing cuts for both legs of a pair
};

#endif
