// -*- mode: C++ -*- 

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
 */
// This class computes background corrections for the FMD. The
// correction is computed in eta,phi cells and the objects stored can
// be put into alien to use with analysis.
//
// Author: Hans Hjersing Dalsgaard, NBI, hans.dalsgaard@cern.ch
//
// 

#ifndef ALIFMDBACKGROUNDCORRECTION_H
#define ALIFMDBACKGROUNDCORRECTION_H

#include "AliFMDInput.h"
#include <TObjArray.h>
#include <AliRunLoader.h>
#include <AliFMDFloatMap.h>
#include <TH2F.h>
class AliTrackReference;

/**
 * Background correction map.
 * 
 */
class AliFMDBackgroundCorrection : public TNamed 
{
public:
  /** 
   * Constructor
   */
  AliFMDBackgroundCorrection() ;
  /** 
   * Destructor
   */
  virtual ~AliFMDBackgroundCorrection() {};
  /** 
   * Generate the background correction
   * 
   * @param from_hits     Whether we run over hits
   * @param nvtxbins      Number of vertex bins
   * @param zvtxcut       Outer cut on interaction point Z coordinate
   * @param nBinsEta      Number of psuedo-rapidity bins
   * @param storeInAlien  Whether to store the result in AliEn
   * @param runNo         Run number 
   * @param endRunNo      End validity 
   * @param filename      The output file name 
   * @param simulate      Whether to run the simulation or not
   * @param nEvents       Number of events to analyse 
   * @param inFile        Whether an input file is specified
   * @param infilename    Input file name
   */
  void GenerateBackgroundCorrection(Bool_t from_hits=kFALSE,
				    Int_t nvtxbins=10,
				    Float_t zvtxcut=10,
				    Int_t nBinsEta=100, 
				    Bool_t storeInAlien = kFALSE, 
				    Int_t runNo =0, 
				    Int_t endRunNo=999999999, 
				    const Char_t* filename="background.root", 
				    Bool_t simulate = kFALSE, 
				    Int_t nEvents=10,
				    Bool_t inFile = kFALSE,
				    const Char_t* infilename="");
  /** 
   * Nested class that derives from AliFMDInput to do the actual
   * processing  
   */  
  class AliFMDInputBG : public AliFMDInput 
  {
  
  public :
    //AliFMDInputBG() ; 
    /** 
     * Constructo
     * 
     * @param hits_not_trackrefs Use hits rather than track references
     */
    AliFMDInputBG(Bool_t hits_not_trackrefs);
    /** 
     * Initialize the object
     * 
     * 
     * @return @c true on success.
     */
    Bool_t Init();
    
    /** 
     * Get number of primaries seen
     * 
     * 
     * @return Number of primaries seen 
     */    
    Int_t GetNprim() {return fPrim;}
    /** 
     * Get number of hits (total number of particles seen)
     * 
     * 
     * @return Total number of particles seen
     */
    Int_t GetNhits() {return fHits;}
    /** 
     * Set the interaction point Z coordinate cut 
     * 
     * @param vtxCut The Z coordinate cut 
     */
    void  SetVtxCutZ(Double_t vtxCut) { fZvtxCut = vtxCut;}
    /** 
     * Set the number of vertex bins
     * 
     * @param nBins Number of vertex bins 
     */
    void  SetNvtxBins(Int_t nBins) { fNvtxBins = nBins;}
    /** 
     * Set the number of eta bins
     * 
     * @param nBins Number of eta bins
     */
    void  SetNbinsEta(Int_t nBins) { fNbinsEta = nBins;}
    /** 
     * Get a list of hits i.e., the particles that hit the FMD
     * 
     * @return List of particles
     */
    TObjArray*  GetHits() {return &fHitArray;}
    /** 
     * Get a list of primaries i.e., the primary particles that hit the FMD
     * 
     * @return List of particles
     */
    TObjArray*  GetPrimaries() {return &fPrimaryArray;}
    /** 
     * Get the run loader used
     * 
     * @return Run loader used 
     */
    AliRunLoader* GetRunLoader() {return fLoader; }
  private:
    /** 
     * Process a single hit
     * 
     * @param h Hit data
     * @param p Particle that made the hit
     * 
     * @return @c false on failure
     */    
    Bool_t ProcessHit(AliFMDHit* h, TParticle* p );
    /** 
     * Process a track reference 
     * 
     * @param tr Track reference
     * @param p  Particle that made the track reference
     * 
     * @return  @c false on failure
     */
    Bool_t ProcessTrackRef(AliTrackReference* tr, TParticle* p );
    /** 
     * Process a single event
     * 
     * @param det     Detector
     * @param ring    Ring
     * @param sector  Sector 
     * @param strip   Strip
     * @param nTrack  Track number
     * @param charge  Charge 
     * 
     * @return @c false on failure
     */
    Bool_t ProcessEvent(UShort_t det,
			Char_t ring, 
			UShort_t sector, 
			UShort_t strip,
			Int_t nTrack,
			Float_t charge);
    /** 
     * Called at start of event
     * 
     * @param event Event number
     * 
     * @return @c false on failure
     */		
    Bool_t Begin(Int_t event );
    /** 
     * Called at end of event 
     * 
     * @return @c false on failure
     */
    Bool_t End();
    TObjArray      fPrimaryArray;     // List of primaries
    TObjArray      fHitArray;         // List of all particles
    TH2F           fPrimaryMapInner;  // Histogram for inners
    TH2F           fPrimaryMapOuter;  // Histogram for outers
    // FIXME: Consider using AliFMDUShortMap, or maybe new class AliFMDIntMap
    AliFMDFloatMap fHitMap;           // A map of hits
    AliFMDFloatMap fLastTrackByStrip; // A map of last track 
    Int_t          fPrim;             // Number of primaries
    Int_t          fHits;             // Number of hits
    Double_t       fZvtxCut;          // Vertex cut
    Int_t          fNvtxBins;         // Number of vertex bins
    Int_t          fPrevTrack;        // Previous track number
    Int_t          fPrevDetector;     // Previous detector number
    Char_t         fPrevRing;         // Previous ring number
    Int_t          fPrevSec;          // Previous sector number
    Int_t          fNbinsEta;         // Number of eta bins
  };
  
private:
  /** 
   * Run a simulation
   * 
   * @param e_t Not used
   */  
  void Simulate(Int_t e);
  /** 
   * Process all primaries from the run loader
   * 
   * @param rl Run loader
   */
  void ProcessPrimaries(AliRunLoader* rl);
  TObjArray fCorrectionArray; // Array of corrections
  TList     fPrimaryList;     // List of primaries
  ClassDef(AliFMDBackgroundCorrection,0)
  
};
#endif
// EOF
