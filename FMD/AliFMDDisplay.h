#ifndef AliFMDDISPLAY_H
#define AliFMDDISPLAY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
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
  /** Look up a color index, based on the value @a x and the maximum
      value of @a x
      @param x   Value 
      @param max Maximum (for example 1023 for digits)
      @return @c false on error  */
  virtual Int_t  LookupColor(Float_t x, Float_t max)  const;
protected:
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
