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
class TObjArray;
class TCanvas;
class TPad;
class TButton;

//___________________________________________________________________
/** @class AliFMDDisplay 
    @brief Utility class to visualize FMD data in geometry. 
    @ingroup FMD_util
 */
class AliFMDDisplay : public AliFMDInput
{
public:
  /** Constructor
      @param gAliceFile galice file*/
  AliFMDDisplay(const char* gAliceFile="galice.root");
  /** DTOR */
  virtual ~AliFMDDisplay() {}
  /** Singleton access function
      @return Singleton object. */
  static AliFMDDisplay* Instance();

  /** Continue to next event */
  void  Continue() { fWait = kFALSE; }
  /** Zoom mode */
  void  Zoom() { fZoomMode = kTRUE; }
  /** Pick mode */
  void  Pick() { fZoomMode = kFALSE; }
  /** Called when a mouse or similar event happens in the display. 
      @param event Event type
      @param px    where the event happened in pixels along X
      @param py    where the event happened in pixels along Y */
  void  ExecuteEvent(Int_t event, Int_t px, Int_t py);
  /** Calculate distance from point @f$ (p_x,p_y)@f$ to this object. 
      @param px Pixel X coordinate 
      @param py Pixel Y coordinate 
      @return distance. */
  Int_t DistancetoPrimitive(Int_t px, Int_t py);
  /** Paint into canvas 
      @param option Not used */
  void  Paint(Option_t* option="") { (void)option; }

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
  /** Visualize data in ESD 
      @param esd FMD ESD data 
      @return  Always @c true */
  virtual Bool_t ProcessESD(AliESDFMD* esd);
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
      fButton(0),
      fZoom(0),
      fPick(0),
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
  void AddMarker(UShort_t det, Char_t rng, UShort_t sec, UShort_t str, 
		 TObject* o, Float_t s, Float_t max);
  
  static AliFMDDisplay* fgInstance; // Static instance 
  Bool_t                fWait;      // Wait until user presses `Continue'
  TObjArray*            fMarkers;   // Cache of markers
  TObjArray*            fHits;      // Cache of `hits'
  TCanvas*              fCanvas;    // Canvas to draw in 
  TPad*                 fPad;       // View pad. 
  TButton*              fButton;    // Continue button
  TButton*              fZoom;      // Zoom button
  TButton*              fPick;      // Pick button
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
