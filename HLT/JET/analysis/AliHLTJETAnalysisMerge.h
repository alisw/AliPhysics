//-*- Mode: C++ -*-

// $Id: AliHLTJETAnalysisMerge.h  $

#ifndef ALIHLTJETANALYSISMERGE_H
#define ALIHLTJETANALYSISMERGE_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTJETAnalysisMerge.h
    @author Jochen Thaeder
    @date   
    @brief  Container merging analysis objects
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt


#include "AliHLTLogging.h"

#include "AliHLTJETAnalysisBase.h"
#include "AliHLTJETAnalysisJets.h"

class TH1;
class TCanvas;
class TObjArray;

class AliHLTJETJets;


/**
 * @class AliHLTJETAnalysisMerge
 * This class merges analysis objects oft the type
 * AliHLTJETAnalysisJets. It takes care of the plotting.
 *
 * @ingroup alihlt_jet
 * @ingroup alihlt_jet_analysis
 */

class AliHLTJETAnalysisMerge : public TObject, public AliHLTLogging {
public:
  
  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */

  /** Constructor */
  AliHLTJETAnalysisMerge();
  
  /** Destructor */
  ~AliHLTJETAnalysisMerge();

  /*
   * ---------------------------------------------------------------------------------
   *                         Initialize / Setup / Reset - public
   * ---------------------------------------------------------------------------------
   */

  /** Initialize class and members */
  Int_t Initialize();

  /*
   * ---------------------------------------------------------------------------------
   *                                Setter - public
   * ---------------------------------------------------------------------------------
   */

  /** Add analysis histogram object to list of AliHLTJETAnalysisJets
   * param jets  Ptr to AliHLTJETAnalysisJets
   */
  void AddJets( AliHLTJETAnalysisJets* jets );

  /*
   * ---------------------------------------------------------------------------------
   *                                 Getter - public
   * ---------------------------------------------------------------------------------
   */

  /** Get List of canvases */
  TObjArray* GetCanvasList() const { return fCanvasArray; }

  /*
   * ---------------------------------------------------------------------------------
   *                             Output - public
   * ---------------------------------------------------------------------------------
   */

  /** Create all Canvases */
  void CreateCanvas();
  
  ///////////////////////////////////////////////////////////////////////////////////
  
 private:
 
  /** copy constructor prohibited */
  AliHLTJETAnalysisMerge(const AliHLTJETAnalysisMerge&);
  
  /** assignment operator prohibited */
  AliHLTJETAnalysisMerge& operator=(const AliHLTJETAnalysisMerge&);

  /*
   * ---------------------------------------------------------------------------------
   *                             Output - private
   * ---------------------------------------------------------------------------------
   */

  /** Create canvas for spectra histograms */
  void CreateCanvasSpectra();

  /** Create canvas for delta histograms */
  void CreateCanvasDelta();

  /** Create canvas for matched histograms */
  void CreateCanvasMatched();
  
  /*
   * ---------------------------------------------------------------------------------
   *                               Helper - private
   * ---------------------------------------------------------------------------------
   */

  /** Add a new canvas to the canvas list
   *  @param name     Name of the canvas
   *  @param divideX  Number of columns of histograms
   *  @param divideY  Number of lines of histograms
   *  @return         Ptr to new created canvas
   */
  TCanvas* AddCanvas( TString name, Int_t divideX, Int_t divideY ); 

  /** Draw a histogram
   *  @param canvas  Ptr to canvas, to be drawn in.
   *  @param idx     Pad index in the canvas
   *  @param hist    Ptr to histogram
   *  @param bScale  If hist should be scaled 
   *  @param bLogY   If hist Y-axis should be logarithmic 
   */
  void DrawHistogram( TCanvas* canvas, Int_t idx, TH1* hist, Bool_t bScale, Bool_t bLogY);

  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */

  /** List of canvases */ 
  TObjArray *fCanvasArray;            //! transient

  /** List of AliHLTJETAnalysisJets objects */
  TObjArray *fAnalysisJetsArray;      //! transient

  ClassDef(AliHLTJETAnalysisMerge, 2);
};
#endif
