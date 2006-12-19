#ifndef AliFMDDISPLAY_H
#define AliFMDDISPLAY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
 */
/** @file    AliFMDDisplay.h
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Mon Mar 27 12:39:09 2006
    @brief   FMD Event display 
*/
//___________________________________________________________________
//
// The classes defined here, are utility classes for reading in data
// for the FMD.  They are  put in a seperate library to not polute the
// normal libraries.  The classes are intended to be used as base
// classes for customized class that do some sort of analysis on the
// various types of data produced by the FMD. 
//
#include "AliFMDInput.h"
#include <TObjArray.h>
class TCanvas;
class TPad;
class TButton;
class TSlider;
class TH1;

//___________________________________________________________________
/** @class AliFMDDisplay 
    @brief Utility class to visualize FMD data in geometry. 
    @ingroup FMD_util
 */
class AliFMDDisplay : public AliFMDInput
{
public:
  /** Constructor
      @param onlyFMD Only show the FMD
      @param gAliceFile galice file*/
  AliFMDDisplay(Bool_t onlyFMD=kTRUE, const char* gAliceFile="galice.root");
  /** DTOR */
  virtual ~AliFMDDisplay();
  /** Singleton access function
      @return Singleton object. */
  static AliFMDDisplay* Instance();

  /** Continue to next event */
  void  Continue() { fWait = kFALSE; }
  /** Zoom mode */
  void  Zoom() { fZoomMode = kTRUE; }
  /** Pick mode */
  void  Pick() { fZoomMode = kFALSE; }
  /** Redisplay the event */ 
  virtual void Redisplay(); // *MENU*
  /** Change cut */
  virtual void ChangeCut();
  /** Called when a mouse or similar event happens in the display. 
      @param event Event type
      @param px    where the event happened in pixels along X
      @param py    where the event happened in pixels along Y */
  virtual void  ExecuteEvent(Int_t event, Int_t px, Int_t py);
  /** Calculate distance from point @f$ (p_x,p_y)@f$ to this object. 
      @param px Pixel X coordinate 
      @param py Pixel Y coordinate 
      @return distance. */
  virtual Int_t DistancetoPrimitive(Int_t px, Int_t py);
  /** Paint into canvas 
      @param option Not used */
  virtual void  Paint(Option_t* option="") { (void)option; }

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
  /** Visualize a hit
      @param hit Hit
      @param p   Track
      @return @c false on error  */
  virtual Bool_t ProcessHit(AliFMDHit* hit, TParticle* p);
  /** Visualize a digit
      @param digit Digit to draw
      @return @c false on error  */
  virtual Bool_t ProcessDigit(AliFMDDigit* digit);
  /** Visualize a raw digit
      @param digit Raw digit.
      @return @c false on error  */
  virtual Bool_t ProcessRaw(AliFMDDigit* digit);
  /** Visualize a reconstructed point.
      @param recpoint Reconstructed point
      @return @c false on error  */
  virtual Bool_t ProcessRecPoint(AliFMDRecPoint* recpoint);
  /** Process ESD data for the FMD.  Users should overload this to
      deal with ESD data. 
      @param d    Detector number (1-3)
      @param r    Ring identifier ('I' or 'O')
      @param s    Sector number (0-19, or 0-39)
      @param t    Strip number (0-511, or 0-255)
      @param eta  Psuedo-rapidity 
      @param mult Psuedo-multiplicity 
      @return  @c false on error  */
  virtual Bool_t ProcessESD(UShort_t d, Char_t r, UShort_t s, UShort_t t, 
			    Float_t eta, Float_t mult);
  /** Look up a color index, based on the value @a x and the maximum
      value of @a x
      @param x   Value 
      @param max Maximum (for example 1023 for digits)
      @return @c false on error  */
  virtual Int_t  LookupColor(Float_t x, Float_t max)  const;
  /** Set multiplicity cut 
      @param cut Cut-off in multiplicity */
  virtual void SetMultiplicityCut(Float_t c=.01) { fMultCut = c; }//*MENU*
  /** Set pedestal width factor 
      @param fac Factor */
  virtual void SetPedestalFactor(Float_t f=3) { fPedestalFactor = f; }//*MENU*
protected:
  /** Copy constructor 
      @param o Object to copy from  */
  AliFMDDisplay(const AliFMDDisplay& o) 
    : AliFMDInput(o),
      fWait(kFALSE),
      fMarkers(0),
      fHits(0),
      fCanvas(0),
      fPad(0),
      fSlider(0),
      fZoomMode(0),
      fX0(0),
      fY0(0),
      fX1(0),
      fY1(0),
      fMultCut(0),
      fPedestalFactor(0),
      fXPixel(0),
      fYPixel(0),
      fOldXPixel(0),
      fOldYPixel(0),
      fLineDrawn(0)
  { } 
  /** Assignment operator 
      @return Reference to this object */
  AliFMDDisplay& operator=(const AliFMDDisplay&) { return *this; } 
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
  
  /** Show only the FMD detectors. */
  void ShowOnlyFMD();
  /** Make base canvas */ 
  virtual void MakeCanvas(const char** which);
  virtual void MakeAux();
  virtual void DrawAux();
  virtual void Idle();
  virtual void AtEnd();
  
  static AliFMDDisplay* fgInstance; // Static instance 
  Bool_t                fWait;      // Wait until user presses `Continue'
  TObjArray*            fMarkers;   // Cache of markers
  TObjArray*            fHits;      // Cache of `hits'
  TCanvas*              fCanvas;    // Canvas to draw in 
  TPad*                 fPad;       // View pad. 
  TObjArray             fButtons;   // Continue button
  TSlider*              fSlider;    // Cut slider
  Bool_t                fZoomMode;  // Whether we're in Zoom mode
  Float_t               fX0;        // X at lower left corner or range 
  Float_t               fY0;        // Y at lower left corner or range 
  Float_t               fX1;        // X at upper right corner or range 
  Float_t               fY1;        // Y at upper right corner or range 
  Float_t               fMultCut;   // Multiplicity cut  
  Float_t               fPedestalFactor; // ADC acceptance factor 
  Int_t                 fXPixel;    // X pixel of mark
  Int_t                 fYPixel;    // Y pixel of mark
  Int_t                 fOldXPixel; // Old x pixel of mark
  Int_t                 fOldYPixel; // Old y pixel of mark
  Bool_t                fLineDrawn; // Whether we're drawing a box
  Bool_t                fOnlyFMD;
  TH1*                  fSpec;
  TH1*                  fSpecCut;
  TCanvas*              fAux;
  ClassDef(AliFMDDisplay,0)  // FMD specialised event display
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
