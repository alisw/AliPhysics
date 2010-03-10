/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/// @file   AliHLTEveBase.h
/// @author Svein Lindal
/// @brief  Base class for the HLT eve display elements

#ifndef ALIHLTEVEBASE_H
#define ALIHLTEVEBASE_H

#include "Rtypes.h"
class AliHLTHOMERBlockDesc;
class AliEveHOMERManager;
class TCanvas;
class TString;

class AliHLTEveBase{

public:
  
  /** Default constructor prohibited **/
  AliHLTEveBase();

  /** Destructor **/
  virtual ~AliHLTEveBase();

  /** Process the incoming blocks, must be implemented by children */
  virtual void ProcessBlock(AliHLTHOMERBlockDesc * block) = 0;

  /** Update the elements after new event loaded, to be implemented by children */
  virtual void UpdateElements() = 0;

  /** Reset the elements before reading in new event, to be implemented by children */
  virtual void ResetElements() = 0;

  /** Set the parent AliEveHOMERManager instance */
  void SetEventManager(AliEveHOMERManager * em) { fEventManager = em; };

  

protected:

  /** Create a new canvas tab */
  TCanvas * CreateCanvas(TString  tabTitle, TString  canvasTitle );
 
  /** Addhistograms to the canvas */
  virtual void AddHistogramsToCanvas(AliHLTHOMERBlockDesc * block, TCanvas * canvas, Int_t &cdCount );

  
  AliEveHOMERManager * fEventManager; //Pointer to AliEveHOMERManager instance
  TCanvas * fCanvas;                  //Canvas for histograms
  Int_t fHistoCount;                  //Counter for histograms, to track where to draw the next one


private:

  /** copy constructor prohibited */
  AliHLTEveBase(const AliHLTEveBase&);
  /** assignment operator prohibited */
  AliHLTEveBase& operator=(const AliHLTEveBase&);

  ClassDef(AliHLTEveBase, 0);
};

#endif
