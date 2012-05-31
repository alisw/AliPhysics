#ifndef AliFMDSPECTRADISPLAY_H
#define AliFMDSPECTRADISPLAY_H
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
#include "AliFMDPattern.h"
#include <TGFrame.h>
#include <TGListTree.h>
#include <TObjArray.h>
#include <TQObject.h>
#include <RQ_OBJECT.h>

class TGPicture;
class TH1;
class TAxis;

//____________________________________________________________________
/**
 * FMD event and spectra display
 * 
 */
class AliFMDSpectraDisplay : public AliFMDPattern
{
public:
  class AliFMDSpectraDisplayTop;
  class AliFMDSpectraDisplayDetector;
  class AliFMDSpectraDisplayRing;
  class AliFMDSpectraDisplaySector;
  class AliFMDSpectraDisplayStrip;

  //__________________________________________________________________
  /**
   * Base class for elements in the tree.
   * 
   */
  class AliFMDSpectraDisplayElement : public TNamed 
  {
  public:
    /** 
     * Destructor.
     * 
     * 
     */  
    virtual ~AliFMDSpectraDisplayElement() {}
    /** 
     * Draw it
     * 
     * @param option Draw option
     * @param l      Low cut 
     * @param h      High cut
     */
    virtual void Show(Option_t* option, Double_t l, Double_t h);
    /** 
     * Get the top of the tree
     * 
     * @return Top element
     */
    virtual AliFMDSpectraDisplayTop& GetTop() = 0;
    /** 
     * Make histograms for this element.
     * 
     * @param axis Axis specs
     */
    virtual void MakeHistograms(const TAxis* axis);
    /** 
     * Compare to object
     * 
     * @param o 
     * 
     * @return 
     */  
    virtual Int_t Compare(const TObject* o) const;
  protected:
    /** 
     * Constructor
     * 
     * @param name  Name 
     * @param title Title
     */  
    AliFMDSpectraDisplayElement(const char* name, const char* title) 
      : TNamed(name, title), fFull(0), fCut(0)
    {}
    AliFMDSpectraDisplayElement(const AliFMDSpectraDisplayElement& o) 
      : TNamed(o), 
	fFull(o.fFull), 
	fCut(o.fCut) 
    {}
    AliFMDSpectraDisplayElement&
    operator=(const AliFMDSpectraDisplayElement&) { return *this; }
    /** 
     * Fill histogram
     * 
     * @param v Value to fill
     */
    void DoFill(Double_t v);
    TH1* fFull; // Full spectra
    TH1* fCut;  // Spectra after cut
    ClassDef(AliFMDSpectraDisplayElement,0); // Element in FMD spectra display
  };

  //__________________________________________________________________
  /**
   * Top element in FMD spectra display
   * 
   */
  class AliFMDSpectraDisplayTop : public AliFMDSpectraDisplayElement, 
				  public TQObject
  {
  public:
    /** 
     * Constructor
     * 
     * @param frame   Parent frame
     * @param canvas  Canvas
     */
    AliFMDSpectraDisplayTop(TGCompositeFrame& frame, TCanvas* canvas);
    /** 
     * Get the list tree
     * 
     * @return List tree
     */
    TGListTree& GetList() { return fList; }
    /** 
     * Get or add a detector element
     * 
     * @param id Element id
     * 
     * @return element
     */
    AliFMDSpectraDisplayDetector& GetOrAdd(UShort_t id);
    /** 
     * Fill histograms
     * 
     * @param d    Detector
     * @param ring Ring 
     * @param sec  Sector
     * @param str  Strip
     * @param v    value
     */
    void      Fill(UShort_t d, Char_t ring, 
		   UShort_t sec, UShort_t str, Double_t v);
    /** 
     * Get the axis spec.
     *
     * @return Axis spec,
     */
    TAxis* GetAxis() const { return fAxis; }
    /** 
     * Set the axis spec.
     * 
     * @param a Axis spec.
     */
    void SetAxis(TAxis* a);
    /** 
     * Get the top.
     * 
     * @return top element
     */
    AliFMDSpectraDisplayTop&   GetTop() { return *this; }
       
    /** 
     * Handle entries 
     *
     * @param e  selected entry, if any 
     * @param id Id of entry 
     */
    virtual void HandleEntry(TGListTreeItem* e, Int_t id);
    /** 
     * Handle key strokes 
     * @param f      Item selected, if any 
     * @param keysym Key symbol 
     * @param mask   Modifier mask 
     */
    virtual void HandleKey(TGListTreeItem* f, UInt_t keysym, UInt_t mask);
    /**
     * Handle Return 
     * @param f Selected item, if any 
     */
    virtual void HandleReturn(TGListTreeItem* f);
    /** 
     * Clear the list 
     */
    virtual void ClearList();
    /** 
     * Clear the canvas 
     */ 
    virtual void ClearCanvas();
    /** 
     * Update canvas 
     */ 
    virtual void UpdateCanvas();
    /** 
     * Update canvas 
     */ 
    virtual void UpdateList();
    /** 
     * Return the currently selected entry 
     */ 
    TGListTreeItem* CurrentEntry() const { return fCurrentEntry; }
    /** 
     * @return the currently selected user data (possibly 0) 
     */
    TObject* Current() const;
    /** 
     * Selection changed signal 
     */
    void SelectionChanged() { Emit("SelectionChanged()"); }//*SIGNAL*
    /** 
     * Get Picture for 1D histogram 
     */
    const TGPicture* GetH1Pic() const { return fkHist1DIcon; }
    /** 
     * 2D Histogram Icon 
     */
    const TGPicture* GetH2Pic() const { return fkHist2DIcon; }
    /** 
     * 3D Histogram Icon 
     */
    const TGPicture* GetH3Pic() const { return fkHist3DIcon; }
    /** 
     * Graph Icon 
     */
    const TGPicture* GetGPic() const { return fkGraphIcon; }
    /** 
     * Get the entry 
     *
     * @return 
     */
    TGListTreeItem&  GetEntry() const { return fEntry; }
    /** 
     * Compare to object
     * 
     * @param o 
     * 
     * @return 
     */  
    virtual Int_t Compare(const TObject* o) const;
  protected:
    AliFMDSpectraDisplayTop(const AliFMDSpectraDisplayTop& o) 
      : AliFMDSpectraDisplayElement(o), 
	TQObject(),
	fHints(o.fHints),
	fContainer(), // o.fContainer), 
	fList(), // o.fList), 
	fChildren(o.fChildren), 
	fCurrentEntry(o.fCurrentEntry),
	fCanvas(o.fCanvas),
	fkHist1DIcon(o.fkHist1DIcon),
	fkHist2DIcon(o.fkHist2DIcon),
	fkHist3DIcon(o.fkHist3DIcon),
	fkGraphIcon(o.fkGraphIcon),
	fAxis(o.fAxis),
	fEntry(o.fEntry)
    {}
  
    AliFMDSpectraDisplayTop& operator=(const AliFMDSpectraDisplayTop&) 
    {
      return *this;
    }
    TGLayoutHints fHints;                // Layout hints
    TGCanvas      fContainer;            // Container
    TGListTree    fList;                 // List
    TObjArray     fChildren;             // Children
    TGListTreeItem* fCurrentEntry;       // Current entry
    TCanvas*        fCanvas;             // Canvas
    const TGPicture*    fkHist1DIcon;    // 1D Histogram Icon
    const TGPicture*    fkHist2DIcon;    // 2D Histogram Icon
    const TGPicture*    fkHist3DIcon;    // 3D Histogram Icon
    const TGPicture*    fkGraphIcon;     // Graph Icon 
    TAxis* fAxis;                        // The axis to use 
    TGListTreeItem& fEntry;              // Entry

    ClassDef(AliFMDSpectraDisplayTop,0);
  };

  //__________________________________________________________________
  /**
   * Detector element in FMD spectra display
   * 
   */
  class AliFMDSpectraDisplayDetector : public AliFMDSpectraDisplayElement 
  {
  public:
    /** 
     * Constructor
     * 
     * @param det 
     * @param top 
     */
    AliFMDSpectraDisplayDetector(UShort_t det, AliFMDSpectraDisplayTop& top);
    /** 
     * Get identifier
     * 
     * @return 
     */
    UShort_t Id() const { return fId; }
    /** 
     * Get top of hierarcy
     * 
     * @return 
     */
    AliFMDSpectraDisplayTop&     GetTop() { return fParent; }
    /** 
     * Get parent
     * 
     * @return 
     */
    AliFMDSpectraDisplayTop&     GetParent() const { return fParent; }
    /** 
     * Get or add a ring element
     * 
     * @param id 
     * 
     * @return 
     */
    AliFMDSpectraDisplayRing&    GetOrAdd(Char_t id);
    /** 
     * Fill histograms
     * 
     * @param ring 
     * @param sec 
     * @param str 
     * @param v 
     */
    void     Fill(Char_t ring, UShort_t sec, UShort_t str, Double_t v);
    /** 
     * Get the entry
     * 
     * 
     * @return 
     */
    TGListTreeItem&  GetEntry() const { return fEntry; }
    /** 
     * Compare to object
     * 
     * @param o 
     * 
     * @return 
     */  
    virtual Int_t Compare(const TObject* o) const;
  protected:
    UShort_t                            fId; // Identifier
    AliFMDSpectraDisplayTop&            fParent; // Parent
    TObjArray                           fChildren; // Children
    TGListTreeItem&                     fEntry; // The entry
    ClassDef(AliFMDSpectraDisplayDetector,0); // Detector element in FMD spectra
  };

  //__________________________________________________________________
  /**
   * Ring element in FMD spectra display
   * 
   */
  class AliFMDSpectraDisplayRing : public AliFMDSpectraDisplayElement 
  {
  public:
    /** 
     * Constructor
     * 
     * @param id 
     * @param d 
     */
    AliFMDSpectraDisplayRing(Char_t id, AliFMDSpectraDisplayDetector& d);
    /** 
     * Get identifier
     * 
     * @return 
     */
    Char_t    Id() const { return fId; }
    /** 
     * Get detector identifier
     * 
     * @return 
     */
    UShort_t  DetectorId() const { return fParent.Id(); }
    /** 
     * Get top of hierarcy
     * 
     * @return 
     */  
    AliFMDSpectraDisplayTop&      GetTop() { return fParent.GetTop(); }
    /** 
     * Get parent detector
     * 
     * @return 
     */
    AliFMDSpectraDisplayDetector& GetDetector() const { return GetParent(); }
    /** 
     * Get parent detector
     * 
     * @return 
     */
    AliFMDSpectraDisplayDetector& GetParent() const { return fParent; }
    /** 
     * Get or add a sector element
     * 
     * @param id 
     * 
     * @return 
     */
    AliFMDSpectraDisplaySector&   GetOrAdd(UShort_t id);
    /** 
     * Fill histograms
     * 
     * @param sec 
     * @param str 
     * @param v 
     */
    void      Fill(UShort_t sec, UShort_t str, Double_t v);
    /** 
     * Get the entry
     * 
     * @return 
     */
    TGListTreeItem&  GetEntry() const { return fEntry; }
    /** 
     * Compare to object
     * 
     * @param o 
     * 
     * @return 
     */  
    virtual Int_t Compare(const TObject* o) const;
  protected:
    AliFMDSpectraDisplayDetector&       fParent; // Parent
    Char_t                              fId;        // Identifier
    TObjArray                           fChildren;  // Children
    TGListTreeItem&                     fEntry;     // Entry
    ClassDef(AliFMDSpectraDisplayRing,0); // Ring element in FMD spectra
  };
  
  //__________________________________________________________________
  /**
   * Sector element in FMD spectra display
   * 
   */
  class AliFMDSpectraDisplaySector : public AliFMDSpectraDisplayElement 
  {
  public:
    /** 
     * Constructor
     * 
     * @param id 
     * @param r 
     */
    AliFMDSpectraDisplaySector(UShort_t id, AliFMDSpectraDisplayRing& r);
    /** 
     * Get identifier
     * 
     * @return 
     */
    UShort_t  Id() const    { return fId; }
    /** 
     * Get detector identifier
     * 
     * @return 
     */
    UShort_t  DetectorId()  const { return fParent.DetectorId(); }
    /** 
     * Get ring identifier
     * 
     * @return 
     */
    Char_t    RingId()      const { return fParent.Id(); }
    /** 
     * Get top of hierarcy
     * 
     * @return 
     */
    AliFMDSpectraDisplayTop&      GetTop()    { return fParent.GetTop(); }
    /** 
     * Get parent detector
     * 
     * @return 
     */
    AliFMDSpectraDisplayDetector& GetDetector() { return fParent.GetDetector(); }
    /** 
     * Get parent ring.
     * 
     * @return 
     */
    AliFMDSpectraDisplayRing&     GetRing()     { return fParent; }
    /** 
     * Get parent ring.
     * 
     * @return 
     */
    AliFMDSpectraDisplayRing&     GetParent()   { return fParent; }
    /** 
     * Get or add a strip element
     * 
     * @param id 
     * 
     * @return 
     */
    AliFMDSpectraDisplayStrip&    GetOrAdd(UShort_t id);
    /** 
     * Fill histograms
     * 
     * @param str 
     * @param v 
     */  
    void      Fill(UShort_t str, Double_t v);
    /** 
     * Get the entry
     * 
     * @return 
     */  
    TGListTreeItem&  GetEntry() const { return fEntry; }
    /** 
     * Compare to object
     * 
     * @param o 
     * 
     * @return 
     */  
    virtual Int_t Compare(const TObject* o) const;
  protected:
    AliFMDSpectraDisplayRing&           fParent;  // PArent
    UShort_t                            fId;      // Identifier
    TObjArray                           fChildren;// Children
    TGListTreeItem&                     fEntry;   // Entry
    ClassDef(AliFMDSpectraDisplaySector,0); // Sector element in FMD spectra
  };

  //__________________________________________________________________
  /**
   * Strip element in FMD spectra display
   * 
   */
  class AliFMDSpectraDisplayStrip : public AliFMDSpectraDisplayElement 
  {
  public:
    /** 
     * Constructor
     * 
     * @param id 
     * @param s 
     */
    AliFMDSpectraDisplayStrip(UShort_t id, AliFMDSpectraDisplaySector& s);
    /** 
     * Get identifier
     * 
     * @return 
     */
    UShort_t  Id() const    { return fId; }
    /** 
     * Get detector identifier
     * 
     * @return 
     */
    UShort_t  DetectorId()  const { return fParent.DetectorId(); }
    /** 
     * Get ring identifier
     * 
     * @return 
     */
    Char_t    RingId()  const    { return fParent.RingId(); }
    /** 
     * Get sector identifier
     * 
     * @return 
     */
    UShort_t  SectorId()  const  { return fParent.Id(); }
    /** 
     * Get top of hierarcy
     * 
     * @return 
     */
    AliFMDSpectraDisplayTop&      GetTop()    { return fParent.GetTop(); }
    /** 
     * Get parent detector
     * 
     * @return 
     */
    AliFMDSpectraDisplayDetector& GetDetector()  { return fParent.GetDetector(); }
    /** 
     * Get parent ring
     * 
     * @return 
     */
    AliFMDSpectraDisplayRing&     GetRing()    { return fParent.GetRing(); }
    /** 
     * Get parent sector
     * 
     * @return 
     */
    AliFMDSpectraDisplaySector&   GetSector()  { return fParent; }
    /** 
     * Get parent sector
     * 
     * @return 
     */
    AliFMDSpectraDisplaySector&   GetParent() { return fParent; }
    /** 
     * Fill histograms
     * 
     * @param v 
     */  
    void      Fill(Double_t v);
    /** 
     * Get entry
     * 
     * @return 
     */
    TGListTreeItem&  GetEntry() const { return fEntry; }
    /** 
     * Compare to object
     * 
     * @param o 
     * 
     * @return 
     */  
    virtual Int_t Compare(const TObject* o) const;
  protected:
    AliFMDSpectraDisplaySector&         fParent; // Parent
    UShort_t                            fId;     // Identifier 
    TGListTreeItem&                     fEntry;  // Entry
    ClassDef(AliFMDSpectraDisplayStrip,0);
  };

  /** 
   * Constructor
   * 
   */
  AliFMDSpectraDisplay();
  /** 
   * Handle draw
   * 
   * @return 
   */
  Bool_t HandleDraw(); 
  /** 
   * Make AUX canvas
   * 
   */
  void MakeAux();
  /** 
   * Draw spectra
   * 
   */
  void DrawAux();
  /** 
   * Process a hit
   * 
   * @param hit 
   * @param p 
   * 
   * @return 
   */
  Bool_t ProcessHit(AliFMDHit* hit, TParticle* p);
  /** 
   * Process a digit
   * 
   * @param digit 
   * 
   * @return 
   */
  Bool_t ProcessDigit(AliFMDDigit* digit);
  /** 
   * Process a summable digit
   * 
   * @param sdigit 
   * 
   * @return 
   */
  Bool_t ProcessSDigit(AliFMDSDigit* sdigit);
  /** 
   * Process a raw digit
   * 
   * @param digit 
   * 
   * @return 
   */
  Bool_t ProcessRawDigit(AliFMDDigit* digit);
  /** 
   * Process a reconstruction point
   * 
   * @param recpoint 
   * 
   * @return 
   */  
  Bool_t ProcessRecPoint(AliFMDRecPoint* recpoint);
  /** 
   * Process and ESD entry
   * 
   * @param det 
   * @param rng 
   * @param sec 
   * @param str 
   * @param x 
   * @param mult 
   * 
   * @return 
   */
  Bool_t ProcessESD(UShort_t det, Char_t rng, UShort_t sec, UShort_t str,
		    Float_t x, Float_t mult);
protected:
  TGMainFrame  fSelector; // GUI to select spectra
  AliFMDSpectraDisplayTop fTop;     // Top of hierarcy
  ClassDef(AliFMDSpectraDisplay,0); // FMD event and spectra display
}
  ;

#endif
//____________________________________________________________________
//
// Local Variables:
//   mode: C++
// End:
//
// EOF
//

