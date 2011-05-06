#ifndef ALIDIELECTRONDEBUGTREE_H
#define ALIDIELECTRONDEBUGTREE_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//#############################################################
//#                                                           # 
//#         Class AliDielectronDebugTree                     #
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

#include <TNamed.h>
#include <TString.h>

#include "AliDielectronVarManager.h"

class TTreeSRedirector;
class AliDielectronPair;
class AliDielectron;

class AliDielectronDebugTree : public TNamed {
public:
  AliDielectronDebugTree();
  AliDielectronDebugTree(const char*name, const char* title);

  virtual ~AliDielectronDebugTree();

  void SetOutputFileName(const char* file) { fFileName=file; }
  
  void AddPairVariable(AliDielectronVarManager::ValueTypes type) { fVariables[fNVars++]=(Int_t)type; }
  void AddLegVariable(AliDielectronVarManager::ValueTypes type)  { fVariablesLeg[fNVarsLeg++]=(Int_t)type; }
  
  void Fill(AliDielectronPair *pair);

  void SetDielectron(AliDielectron * const dielectron) { fDielectron=dielectron; }
  
  void DeleteStreamer();
  void WriteTree();
private:
  TString fFileName;                                          //output file name
  
  Int_t  fNVars;                                              //number of configured variables
  Int_t  fVariables[AliDielectronVarManager::kNMaxValues];    //configured variables
  Int_t  fNVarsLeg;                                           //number of configured variables
  Int_t  fVariablesLeg[AliDielectronVarManager::kNMaxValues]; //configured variables for the legs

  TTreeSRedirector *fStreamer;     //! Tree Redirector
  AliDielectron *fDielectron;      //! pointer to mother dielectron manager

  AliDielectronDebugTree(const AliDielectronDebugTree &c);
  AliDielectronDebugTree &operator=(const AliDielectronDebugTree &c);

  
  ClassDef(AliDielectronDebugTree,1)         // Dielectron DebugTree
};



#endif
