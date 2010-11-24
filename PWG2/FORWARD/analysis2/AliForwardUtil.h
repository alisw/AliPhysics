#ifndef ALIROOT_PWG2_FORWARD_ALIFORWARDUTIL_H
#define ALIROOT_PWG2_FORWARD_ALIFORWARDUTIL_H
#include <TObject.h>
#include <TString.h>
class TH2D;
class TH1I;
class TH1;
class TAxis;
class AliESDEvent;

/** 
 * Utilities used in the forward multiplcity analysis 
 * 
 * @ingroup pwg2_forward_analysis 
 */
class AliForwardUtil : public TObject
{
public:
  /** 
   * Read the trigger information from the ESD event 
   * 
   * @param esd        ESD event 
   * @param triggers   On return, contains the trigger bits 
   * @param hTriggers  Histogram to fill 
   * 
   * @return @c true on success, @c false otherwise 
   */
  static Bool_t ReadTriggers(AliESDEvent* esd, UInt_t& triggers, 
			     TH1I* hTriggers);
  /** 
   * Read the vertex information from the ESD event 
   * 
   * @param esd  ESD event 
   * @param vz   On return, the vertex Z position 
   * 
   * @return @c true on success, @c false otherwise 
   */
  static Bool_t ReadVertex(AliESDEvent* esd, Double_t& vz,
			   Double_t maxErr=0.1);
  //__________________________________________________________________
  /** 
   * Structure to hold histograms 
   *
   * @ingroup pwg2_forward_analysis 
   */
  struct Histos : public TObject
  {	
    /** 
     * Constructor 
     * 
     * 
     */
    Histos() : fFMD1i(0), fFMD2i(0), fFMD2o(0), fFMD3i(0), fFMD3o(0) {}
    /** 
     * Copy constructor 
     * 
     * @param o Object to copy from 
     */
    Histos(const Histos& o) 
      : TObject(o), 
	fFMD1i(o.fFMD1i), 
	fFMD2i(o.fFMD2i), 
	fFMD2o(o.fFMD2o), 
	fFMD3i(o.fFMD3i), 
	fFMD3o(o.fFMD3o)
    {}
    /** 
     * Assignement operator 
     * 
     * @return Reference to this 
     */
    Histos& operator=(const Histos&) { return *this;}
    /** 
     * Destructor
     */
    ~Histos();
    /** 
     * Initialize the object 
     * 
     * @param etaAxis Eta axis to use 
     */
    void Init(const TAxis& etaAxis);
    /** 
     * Make a histogram 
     * 
     * @param d        Detector
     * @param r        Ring 
     * @param etaAxis  Eta axis to use
     * 
     * @return Newly allocated histogram 
     */
    TH2D* Make(UShort_t d, Char_t r, const TAxis& etaAxis) const;
    /** 
     * Clear data 
     * 
     * @param option Not used 
     */
    void  Clear(Option_t* option="");
    // const TH2D* Get(UShort_t d, Char_t r) const;
    /** 
     * Get the histogram for a particular detector,ring
     * 
     * @param d Detector 
     * @param r Ring 
     * 
     * @return Histogram for detector,ring or nul 
     */
    TH2D* Get(UShort_t d, Char_t r) const;
    TH2D* fFMD1i; // Histogram for FMD1i
    TH2D* fFMD2i; // Histogram for FMD2i
    TH2D* fFMD2o; // Histogram for FMD2o
    TH2D* fFMD3i; // Histogram for FMD3i
    TH2D* fFMD3o; // Histogram for FMD3o

    ClassDef(Histos,1) 
  };

  //__________________________________________________________________
  struct RingHistos : public TObject
  {
    RingHistos() : fDet(0), fRing('\0'), fName("") {}
    RingHistos(UShort_t d, Char_t r) 
      : fDet(d), fRing(r), fName(TString::Format("FMD%d%c", d, r)) 
    {}
    RingHistos(const RingHistos& o) 
      : TObject(o), fDet(o.fDet), fRing(o.fRing), fName(o.fName)
    {}
    virtual ~RingHistos() {}
    RingHistos& operator=(const RingHistos& o) 
    {
      TObject::operator=(o);
      fDet  = o.fDet;
      fRing = o.fRing;
      fName = o.fName;
      return *this;
    }
    TList* DefineOutputList(TList* d) const;
    TList* GetOutputList(TList* d) const;
    TH1* GetOutputHist(TList* d, const char* name) const;
    UShort_t fDet; 
    Char_t   fRing;
    TString  fName;

    ClassDef(RingHistos,1) 
  };
    
};

#endif
// Local Variables:
//  mode: C++
// End:

