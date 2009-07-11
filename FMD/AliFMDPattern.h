#ifndef AliFMDPATTERN_H
#define AliFMDPATTERN_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
 */
/** @file    AliFMDPattern.h
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 12:39:09 2006
    @brief   FMD Event display (as patterns)
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
#include <TArrayI.h>
#include <TObjArray.h>
// #include <TGraph.h>
#include <TLatex.h>
#include <TLine.h>
class TCanvas;
class TPad;
class TH2;


//___________________________________________________________________
/** @class AliFMDPattern 
    @brief Utility class to visualize FMD data in 2D. 
    @ingroup FMD_util
 */
class AliFMDPattern : public AliFMDDisplay
{
public:
  /** @class AliFMDPatternDetector 
      @brief Utility for the pattern display. 
      The name is this long, because the code-checker even cares about
      nested names, even though it seems a bit nit-picking and
      retareded.   It's a class rather than a structure, because the
      code checker thinks that structs are only for POD - sigh! */
  class AliFMDPatternDetector 
  {
  public:
    /** Constructor */
    AliFMDPatternDetector(UShort_t id);
    /** Destructor */
    ~AliFMDPatternDetector();
    /** Clear this detector */
    void  Clear();
    /** Get the total number of hits */
    Int_t Total() const { return Int_t(fCounts.GetSum()); }
    /** Called at the  end of things */ 
    void  End();
    /** Initiate. 
	@param nlvl Number of levels. 
	@param r    Radius 
	@param inners Array of inner shapes  
	@param outers Array of outer shapes */
    void  Begin(Int_t nlvl, Double_t r, TObjArray& inners, 
		TObjArray&  outers);
    /** Draw everything 
	@param a Array of shapes to draw */
    void  DrawShape(TObjArray& a);
    /** Add a marker at specified coordinates 
	@param X X coordinate 
	@param Y Y coordinate 
	@param Z Z coordinate 
	@param max The maximum value to scale to */
    void  AddMarker(Double_t x, Double_t y, Float_t s, Float_t max);
    Int_t     fId;     // Identifier # 
    TArrayI   fCounts; // Number of counts at each level 
    TObjArray fGraphs; // Array of graphs - one for each level
    TH2*      fFrame;  // The mother frame 
  private:
    /** Copy constructor 
	- Not implemented. */ 
    AliFMDPatternDetector(const AliFMDPatternDetector&);
    /** Assignement operator 
	-- Not implemented */
    AliFMDPatternDetector& operator=(const AliFMDPatternDetector&);
    void CopyShapes(TObjArray& input, TObjArray& own, 
		    Double_t ang=0, Double_t fx=1, Double_t fy=1);
    /** Our own cache of shapes */
    TObjArray fInners;
    /** Our own cache of shapes */
    TObjArray fOuters;
  };
  
  
  /** Constructor
      @param gAliceFile galice file*/
  AliFMDPattern(const char* gAliceFile="galice.root");
  /** DTOR */
  virtual ~AliFMDPattern();

  /** Initialize
      @return  @c false on error */
  virtual Bool_t Init();
  /** Called at beginning of an event 
      @param event Event number
      @return @c false on error  */
  virtual Bool_t Begin(Int_t event);
 protected:
  /** Copy constructor 
      - Not implemented. */ 
  AliFMDPattern(const AliFMDPattern&);
  /** Assignement operator 
      -- Not implemented */
  AliFMDPattern& operator=(const AliFMDPattern&);
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
  /** @param hit Hit to process */
  virtual Bool_t ProcessHit(AliFMDHit* hit, TParticle*);
  /** Re-draw the patterns */
  virtual void Redisplay();
  /** Called at the end. */
  virtual void AtEnd();
  /** Graph to show shape of inner sensor */
  TObjArray fInners;	// Graph to show shape of inner sensor
  /** Graph to show shape of outer sensor */
  TObjArray fOuters;	// Graph to show shape of outer sensor
  /** Max inner radius */
  Float_t fInnerMax;	// Max inner radius
  /** Max outer radius */
  Float_t fOuterMax;	// Max outer radius
  /** FMD1 Pad */
  TPad*  fFMD1Pad;	// FMD1 Pad
  /** FMD1 Frame */
  AliFMDPatternDetector fFMD1;	// FMD1 Frame
  /** FMD2 Pad  */
  TPad*  fFMD2Pad;	// FMD2 Pad 
  /** FMD2 Frame */
  AliFMDPatternDetector fFMD2;	// FMD2 Frame
  /** FMD3 Pad */
  TPad*  fFMD3Pad;	// FMD3 Pad
  /** FMD3 Frame */
  AliFMDPatternDetector fFMD3;	// FMD3 Frame
  /** Summary pad */
  TPad* fSummary;	// Summary pad
  /** Text fields */
  TLatex fEvent;	// Text fields
  TLatex fFMD1Sum;	// Total in FMD1
  TLatex fFMD2Sum;	// Total in FMD1
  TLatex fFMD3Sum;	// Total in FMD1
  TLine  fLine;		// Just a line 
  TLatex fTotal;	// Total in FMD

  Double_t fFMD1Area;   // 
  Double_t fFMD2Area;   // 
  Double_t fFMD3Area;   // 

  ClassDef(AliFMDPattern,0) // Display FMD data as hit-patterns. 
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
