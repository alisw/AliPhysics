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
class TH3;


//___________________________________________________________________
/** @class AliFMDFancy 
    @brief Utility class to visualize FMD data in 2D. 
    @ingroup FMD_util
 */
class AliFMDFancy : public AliFMDDisplay
{
public:
  struct Detector 
  {
    Detector(UShort_t id);
    ~Detector();
    
    void Init();
    void Begin(Int_t event=0);
    void Clear(Int_t event=0);
    void End();
    void AddMarker(Char_t rng, UShort_t sec, UShort_t str, 
		   Float_t v, Float_t max);
    TH1*      fFrame;
    Int_t     fId;
    TObjArray fShapes;
    Int_t     fNInnerHits;
    TGraph2D  fInnerHits;
    Int_t     fNOuterHits;
    TGraph2D  fOuterHits;
    Double_t  fMaxR;
    Double_t  fMinZ;
    Double_t  fMaxZ;
  private:
    void AddHistogram(TGraph2D& g, const char* opt="");
    Detector(const Detector& );
    Detector& operator=(const Detector& ) { return *this; }
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
  AliFMDFancy(const AliFMDFancy& );
  AliFMDFancy& operator=(const AliFMDFancy& ) { return *this; }
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
  TPad*    fSummary;
  /** Text fields */
  TLatex fEvent;
  TLatex fFMD1IHits;
  TLatex fFMD2IHits;
  TLatex fFMD2OHits;
  TLatex fFMD3IHits;
  TLatex fFMD3OHits;
  TLine  fLine;
  TLatex fTotal;
  
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
