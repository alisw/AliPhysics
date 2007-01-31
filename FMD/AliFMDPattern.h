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
#include <TGraph.h>
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
  struct Detector 
  {
    Detector(UShort_t id);
    ~Detector();
    void  Clear();
    Int_t Total() const { return Int_t(fCounts.GetSum()); }
    void  End();
    void  Begin(Int_t nlvl, Double_t r, TObjArray& inners, TObjArray& outers);
    void  DrawShape(TObjArray& a);
    void  AddMarker(Double_t x, Double_t y, Float_t s, Float_t max);
    Int_t     fId;
    TArrayI   fCounts;
    TObjArray fGraphs;
    TH2*      fFrame;
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
  /** Add a marker to the display
      @param det Detector
      @param rng Ring
      @param sec Sector 
      @param str Strip
      @param o   Object to refer to
      @param s   Signal 
      @param max Maximum of signal */
  virtual void AddMarker(UShort_t det, Char_t rng, UShort_t sec, UShort_t str, 
			 TObject* o, Float_t s, Float_t max);
  virtual Bool_t ProcessHit(AliFMDHit* hit, TParticle*);
  virtual void Redisplay();
  virtual void AtEnd();
  /** Graph to show shape of inner sensor */
  TObjArray fInners;
  /** Graph to show shape of outer sensor */
  TObjArray fOuters;
  /** Max inner radius */
  Float_t fInnerMax;
  /** Max outer radius */
  Float_t fOuterMax;
  /** FMD1 Pad */
  TPad*  fFMD1Pad;
  /** FMD1 Frame */ 
  Detector fFMD1;
  /** FMD2 Pad  */
  TPad*  fFMD2Pad;
  /** FMD2 Frame */ 
  Detector fFMD2;
  /** FMD3 Pad */
  TPad*  fFMD3Pad;
  /** FMD3 Frame */ 
  Detector fFMD3;
  /** Summary pad */
  TPad* fSummary;
  /** Text fields */
  TLatex fEvent;
  TLatex fFMD1Sum;
  TLatex fFMD2Sum;
  TLatex fFMD3Sum;
  TLine  fLine;
  TLatex fTotal;

  ClassDef(AliFMDPattern,0)
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
