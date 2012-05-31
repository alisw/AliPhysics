#ifndef AliFMDFANCY_H
#define AliFMDFANCY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
 */
/** @file    AliFMDFancy.h
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 12:39:09 2006
    @brief   FMD Event display (as fancys)
*/
//___________________________________________________________________
//
// The classes defined here, are utility classes for reading in data
// for the FMD.  They are  put in a seperate library to not polute the
// normal libraries.  The classes are intended to be used as base
// classes for customized class that do some sort of analysis on the
// various types of data produced by the FMD. 
//
#include "AliFMDDisplay.h"
#include <TObjArray.h>
#include <TGraph2D.h>
#include <TLatex.h>
#include <TLine.h>
class TCanvas;
class TPad;
class TH1;
class TH2;
// class TH3;


//___________________________________________________________________
/** @class AliFMDFancy 
    @brief Utility class to visualize FMD data in 2D. 
    @ingroup FMD_util
 */
class AliFMDFancy : public AliFMDDisplay
{
public:
  /** Sigh! the code checker thinks that all structs are POD - morron
      */ 
  class AliFancyDetector 
  {
  public:
    /** CTOR */
    AliFancyDetector(UShort_t id);
    /** DTOR */
    ~AliFancyDetector();
    /** Initialise */
    void Init();
    /** Called at the beginning of an event */
    void Begin(Int_t event=0);
    /** Clear display */
    void Clear(Int_t event=0);
    /** Called that the end of an event */
    void End();
    /** Add a marker */
    void AddMarker(Char_t rng, UShort_t sec, UShort_t str, 
		   Float_t v, Float_t max);

    /**
     * Get  Mother frame 
     * 
     * @return  Mother frame 
     */
    TH1*& GetFrame() { return fFrame; }
    /**
     * Set  Mother frame 
     *
     * @param p
     */
    void SetFrame(TH1* p) { fFrame = p; }
    /**
     * Get  Identifier 
     * 
     * @return  Identifier 
     */
    Int_t GetId() const { return fId; }
    /**
     * Set  Identifier 
     *
     * @param p
     */
    void SetId(Int_t p) { fId = p; }
    /**
     * Get  Array of shapes 
     * 
     * @return  Array of shapes 
     */
    TObjArray& GetShapes() { return fShapes; }
    /**
     * Get  Number of inner hits 
     * 
     * @return  Number of inner hits 
     */
    Int_t& GetNInnerHits() { return fNInnerHits; }
    /**
     * Get  Graph of inner hits 
     * 
     * @return  Graph of inner hits 
     */
    TGraph2D& GetInnerHits() { return fInnerHits; }
    /**
     * Get  Number of outer hits 
     * 
     * @return  Number of outer hits 
     */
    Int_t& GetNOuterHits() { return fNOuterHits; }
    /**
     * Get  Graph  of outer hits 
     * 
     * @return  Graph  of outer hits 
     */
    TGraph2D& GetOuterHits() { return fOuterHits; }
    /**
     * Get  Maximum radius 
     * 
     * @return  Maximum radius 
     */
    Double_t GetMaxR() const { return fMaxR; }
    /**
     * Set  Maximum radius 
     *
     * @param p
     */
    void SetMaxR(Double_t p) { fMaxR = p; }
    /**
     * Get  Minimum Z 
     * 
     * @return  Minimum Z 
     */
    Double_t GetMinZ() const { return fMinZ; }
    /**
     * Set  Minimum Z 
     *
     * @param p
     */
    void SetMinZ(Double_t p) { fMinZ = p; }
    /**
     * Get  Maximum Z 
     * 
     * @return  Maximum Z 
     */
    Double_t GetMaxZ() const { return fMaxZ; }
    /**
     * Set  Maximum Z 
     *
     * @param p
     */
    void SetMaxZ(Double_t p) { fMaxZ = p; }
  protected:
    /** Mother frame */
    TH1*      fFrame;	// Mother frame
    /** Identifier */
    Int_t     fId;	// Identifier
    /** Array of shapes */
    TObjArray fShapes;	// Array of shapes
    /** Number of inner hits */
    Int_t     fNInnerHits;	// Number of inner hits
    /** Graph of inner hits */
    TGraph2D  fInnerHits;	// Graph of inner hits
    /** Number of outer hits */
    Int_t     fNOuterHits;	// Number of outer hits
    /** Graph  of outer hits */
    TGraph2D  fOuterHits;	// Graph  of outer hits
    /** Maximum radius */
    Double_t  fMaxR;	// Maximum radius
    /** Minimum Z */
    Double_t  fMinZ;	// Minimum Z
    /** Maximum Z */
    Double_t  fMaxZ;	// Maximum Z
  private:
    /** Add a histogram to a 2D graph.  For some reason the code
	checker thinks that this function can be made const - well, it
	cannot, since the graph passed down is a member of this
	object, and would be const in the this context if the member
	function is const.   Since we modify the graph, we cannot make
	it a const reference, no matter how much we'd like to. */
    void AddHistogram(TGraph2D& g, const char* toopt="");
    /** Copy ctor */
    AliFancyDetector(const AliFancyDetector& );
    /** Assignement operator */
    AliFancyDetector& operator=(const AliFancyDetector& ) { return *this; }
  };
    
  /** Constructor
      @param gAliceFile galice file*/
  AliFMDFancy(const char* gAliceFile="galice.root");
  /** DTOR */
  virtual ~AliFMDFancy();

  /** Initialize
      @return  @c false on error */
  virtual Bool_t Init();
  /** Called at beginning of an event 
      @param event Event number
      @return @c false on error  */
  virtual Bool_t Begin(Int_t event);
  /** Called at end of an event 
      @return @c false on error  */
  virtual Bool_t End();
protected:
  /** Copy ctor */
  AliFMDFancy(const AliFMDFancy& );
  /** Assignement operator */
  AliFMDFancy& operator=(const AliFMDFancy& ) { return *this; }
  virtual void AddMarker(Float_t x, Float_t y, Float_t z, 
			 TObject* o, Float_t s, Float_t min, Float_t max)
  {
    AliFMDDisplay::AddMarker(x, y, z, o, s, min, max);
  }
  /** Add a marker to the display
      @param det Detector
      @param rng Ring
      @param sec Sector 
      @param str Strip
      @param o   Object to refer to
      @param s   Signal 
      @param max Maximum of signal */
  virtual void AddMarker(UShort_t det, Char_t rng, UShort_t sec,
			 UShort_t str, TObject* o, Float_t s, 
			 Float_t min, Float_t max);
  /** Process a hit 
      @param hit hit to process */
  virtual Bool_t ProcessHit(AliFMDHit* hit, TParticle*);

  /** FMD1 Pad */
  TPad*  fFMD1Pad;		// FMD1 Pad
  /** FMD1 Frame */
  AliFancyDetector fFMD1;	// FMD1 Frame
  /** FMD2 Pad  */
  TPad*  fFMD2Pad;		// FMD2 Pad 
  /** FMD2 Frame */
  AliFancyDetector fFMD2;	// FMD2 Frame
  /** FMD3 Pad */
  TPad*  fFMD3Pad;		// FMD3 Pad
  /** FMD3 Frame */
  AliFancyDetector fFMD3;	// FMD3 Frame
  /** Summary pad */
  TPad*    fSummary;		// Summary pad
  /** Text fields */
  TLatex fEvent;		// Text fields
  /** Number of hits in FMD1i */
  TLatex fFMD1IHits;		// Number of hits in FMD1i
  /** Number of hits in FMD2i */
  TLatex fFMD2IHits;		// Number of hits in FMD2i
  /** Number of hits in FMD2o */
  TLatex fFMD2OHits;		// Number of hits in FMD2o
  /** Number of hits in FMD3i */
  TLatex fFMD3IHits;		// Number of hits in FMD3i
  /** Number of hits in FMD3o */
  TLatex fFMD3OHits;		// Number of hits in FMD3o
  /** Just a line */
  TLine  fLine;			// Just a line
  /** Number of hits in FMD */
  TLatex fTotal;		// Number of hits in FMD
  
  ClassDef(AliFMDFancy,0)
};


#endif
//____________________________________________________________________
//
// Local Variables:
//   mode: C++
// End:
//
// EOF
//
