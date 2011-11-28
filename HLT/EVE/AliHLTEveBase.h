//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTEVEBASE_H
#define ALIHLTEVEBASE_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/// @file   AliHLTEveBase.h
/// @author Svein Lindal <slindal@fys.uio.no>
/// @brief  Base class for the HLT eve display elements

#include "Rtypes.h"
#include "TString.h"
#include "TEveElement.h"
class AliHLTHOMERBlockDesc;
class AliEveHLTEventManager;
class TCanvas;


class AliHLTEveBase : public TEveElementList {

public:
  
  /** Default constructor prohibited **/
  AliHLTEveBase(const char * name);

  /** Destructor **/
  virtual ~AliHLTEveBase();

  /** Process the incoming blocks, must be implemented by children */
  virtual void ProcessBlock(AliHLTHOMERBlockDesc * block) = 0;

  /** Update the elements after new event loaded, to be implemented by children */
  virtual void UpdateElements() = 0;

  /** Reset the elements before reading in new event, to be implemented by children */
  virtual void ResetElements() = 0;

  /** Set the parent AliEveHLTEventManager instance */
  void SetEventManager(AliEveHLTEventManager * em) { fEventManager = em; };

  

protected:

  /** Create a new canvas tab */
  TCanvas * CreateCanvas(TString  tabTitle, TString  canvasTitle );
 
  /** Addhistograms to the canvas */
  virtual void AddHistogramsToCanvas(AliHLTHOMERBlockDesc * block, TCanvas * canvas, Int_t &cdCount );


	/// Getters and setters for the max number of histograms
  void SetMaxHistograms(Int_t mh) {fMaxHistos = mh;}
  Int_t GetMaxHistograms() const {return fMaxHistos;}
  
	///Getter and setter for the detector string
  void SetDetector(TString det) {fDetector = det;}
  TString GetDetector() const {return fDetector;}

  
  AliEveHLTEventManager * fEventManager; //Pointer to AliEveHLTEventManager instance
  TCanvas * fCanvas;                  //Canvas for histograms
  Int_t fHistoCount;                  //Counter for histograms, to track where to draw the next one


private:

  /** Default constructor prohibited **/
  AliHLTEveBase();
  /** copy constructor prohibited */
  AliHLTEveBase(const AliHLTEveBase&);
  /** assignment operator prohibited */
  AliHLTEveBase& operator=(const AliHLTEveBase&);

  Int_t fMaxHistos;  // Maximum number histograms there is room for. 

  TString fDetector;  //String denoting the detector

  ClassDef(AliHLTEveBase, 0);
};

#endif
