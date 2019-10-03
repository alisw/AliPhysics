#ifndef ALIFMDENCODEDEDX_H
#define ALIFMDENCODEDEDX_H
/**
 * @file   AliFMDEncodedEdx.h
 * @author Christian Holm Christensen <cholm@nbi.dk>
 * @date   Thu Oct 10 11:50:59 2013
 * 
 * @brief  Class to encode the energy loss @f$d\Delta@f$ and path length
 * @f$dx@f$ of a particle into track reference bits.
 */
#ifndef __CINT__
# include <TMath.h>
# include <TArrayD.h>
# include <TH1D.h>
# include <TH2D.h>
#else
class TArrayD;
class TH1;
class TH2;
#endif

/**
 * Class to encode the energy loss @f$d\Delta@f$ and path length
 * @f$dx@f$ of a particle into track reference bits.  A total of 13
 * bits are available.  Of these 8 are used for @f$d\Delta@f$ and 5
 * for @f$dx@f$.
 *
 * The encoding is done by binning.  That is, for a given value of
 * @f$d\Delta@f$ or @f$dx@f$ we calculate a bin number and store that.
 * The reverse decoding is done by looking up the bin center of the
 * bin number stored.  Note, that the bin numbers go from 0 to 255
 * (for @f$d\Delta@f$) and 31 (for @f$dx@f$).
 *
 * The bins become progressively wider.  That is, we define 3 regions. 
 * 
 * - Lower region which spans the lower 10% of the distribution has
 *   75% of all avaible bins.
 * - Middle region which spans from the 10% to 20% of the distribtion
 *    and has 20% of the available bins.
 * - Upper region which covers the rest of the distribtion in 5% of
 *   the available bins
 *
 *  | Type          | N bins | Range | bin | value | bin | value |
 *  | ------------- | ------ | ----- | --- | ----- | --- | ----- |
 *  | @f$d\Delta@f$ |  256   | 0-11  | 192 |  1.1  | 243 |  2.2  |
 *  | @f$dx@f$      |   32   | 0-0.7 |  24 |  0.07 |  30 |  1.4  |
 *      
 * Currently there's no support of other schemas 
 * 
 * 
 */
class AliFMDEncodedEdx
{
public:
  /**
   * Specification of the bins 
   * 
   */
  struct Spec {
    UShort_t nBins;     // Total number of bins
    Double_t min;       // Least value 
    Double_t max;       // Largest value 
    UShort_t cutBin1;   // Last bin (plus one) of lower region
    Double_t cutVal1;   // Upper cut on lower region
    UShort_t cutBin2;   // Last bin (plus one) of middle region
    Double_t cutVal2;	// Upper cut on middle region	      
    /** 
     * Constructor 
     * 
     * @param nb Total number of bins 
     * @param l  Lower value 
     * @param h  Upper value 
     */
    Spec(UShort_t nb, Double_t l, Double_t h)
      : nBins(nb), min(l), max(h), 
	cutBin1(UShort_t(nb * 0.75 + .5)), 
	cutVal1((max-min) / 10 + min),
	cutBin2(UShort_t(nb * 0.95 + .5)),
	cutVal2((max-min) / 5  + min)
    {}
    /** 
     * Encode a value 
     * 
     * @param v Value to encode 
     * 
     * @return Encoded (bin number) value 
     */
    UInt_t Encode(Double_t v) const
    {
      UInt_t   off  = 0;
      UInt_t   n    = cutBin1;
      Double_t low  = min;
      Double_t high = cutVal1;
      if (v > cutVal2) { 
	// Upper part of the plot
	off  = cutBin2;
	n    = nBins - cutBin2;
	low  = cutVal2;
	high = max;
      }
      else if (v > cutVal1) {
	// Middle part 
	off  = cutBin1;
	n    = cutBin2 - cutBin1;
	low  = cutVal1;
	high = cutVal2;
      }
      return off + UInt_t(n*(v-low)/(high-low));
    }
    /** 
     * Decode a bin into a value 
     * 
     * @param b Encoded (bin number) value 
     * @param w On return, the width of the bin
     * 
     * @return Decoded value (center of bin)
     */
    Double_t Decode(UInt_t b, Double_t& w) const
    {
      Double_t off  = min;
      UInt_t   n    = cutBin1;
      Double_t high = cutVal1;
      if (b >= cutBin2) {
	// Upper part
	off  = cutVal2;
	n    = nBins - cutBin2;
	high = max;
	b    -= cutBin2;
      }
      else if (b >= cutBin1) {
	// Middle part
	off  = cutVal1;
	n    = cutBin2 - cutBin1;
	high = cutVal2;
	b    -= cutBin1;
      }
      w = (high-off)/n;
      return off + w * (b+.5);
    }
    /** 
     * Decode a bin into a value 
     * 
     * @param b Encoded (bin number) value 
     * 
     * @return Decoded value (center of bin)
     */
    Double_t Decode(UInt_t b) const 
    {
      Double_t w;
      return Decode(b, w);
    }
    /** 
     * Fill an array with values appropriate for defining a histogram 
     * axis with the _natural_ binning of the encoding 
     * 
     * @param a On return, the modified bin-border array 
     */
    void FillBinArray(TArrayD& a) const
    {
      a.Set(nBins+1);
      a[0] = min;
      Double_t w0 = (cutVal1 - min)     / cutBin1;
      Double_t w1 = (cutVal2 - cutVal1) / (cutBin2 - cutBin1);
      Double_t w2 = (max     - cutVal2) / (nBins   - cutBin2);
      for (UInt_t i = 1;         i <= cutBin1; i++) a[i] = a[i-1] + w0;
      for (UInt_t i = cutBin1+1; i <= cutBin2; i++) a[i] = a[i-1] + w1;
      for (UInt_t i = cutBin2+1; i <= nBins;   i++) a[i] = a[i-1] + w2;
    }
    /** 
     * Print information.
     * 
     * @param opt If this starts with `T' also run a test 
     */
    void Print(Option_t* opt="") const 
    {
      Printf("Spec: [%8.4f,%8.4f] in %3d bins, cuts %8.4f (%3d) %8.4f (%3d)",
	     min, max, nBins, cutVal1, cutBin1, cutVal2, cutBin2);
      if (opt[0] == 'T' || opt[0] == 't') {
	for (Int_t i = 0; i < nBins; i++) { 
	  Double_t w = 0;
	  Double_t x = Decode(i,w );
	  UInt_t   j = Encode(x);
	  Printf("%3d -> %8.4f (%7.4f) -> %3d", i, x, w, j);
	}
      }
    }
    /** 
     * Run a test 
     * 
     */
    static void Test() 
    {
      Spec s(125, 0, 125);
      s.Print("T");
    }
  };
  /**
   * How the 13 bits are distributed 
   */
  enum { 
    kNEBits = 8, 
    kNLBits = 5,
    kEMask  = 0xFF, // (1 << kNEBits) - 1
    kLMask  = 0x1F  // (1 << kNLBits) - 1
  };
  /** 
   * Constructor - a no-op
   */
  AliFMDEncodedEdx() {}
  /** 
   * Destructor - a no-op 
   */
  virtual ~AliFMDEncodedEdx() {}
  /** 
   * Get the @f$d\Delta@f$ bin specification. If not initialized
   * already, do so .
   * 
   * @return Constant reference to @f$d\Delta@f$ bin specification
   */
  static const Spec& GetdEAxis() 
  { 
    static Spec* dEAxis = 0;
    if (!dEAxis) dEAxis = new Spec((1<<kNEBits), 0 /*0.0000025*/, 11);
    return *dEAxis;
  }
  /** 
   * Get the @f$dx@f$ bin specification. If not initialized
   * already, do so .
   * 
   * @return Constant reference to @f$dx@f$ bin specification
   */
  static const Spec& GetdLAxis() 
  { 
    static Spec* dLAxis = 0;
    if (!dLAxis) dLAxis = new Spec((1<<kNLBits), 0 /*0.00014*/,   0.7);
    return *dLAxis;
  }
  /** 
   * Encode @f$d\Delta@f$ and @f$dx@f$ into a 13bit number.  
   * 
   * @param edep   @f$d\Delta@f$
   * @param length @f$dx@f$ 
   * 
   * @return 13-bit (lower) encoded value 
   */
  static UInt_t Encode(Double_t edep, Double_t length)
  {
    UInt_t uE = EncodeOne(edep,   GetdEAxis());
    UInt_t uL = EncodeOne(length, GetdLAxis());
    return (((uE & kEMask) << 0) | ((uL & kLMask) << kNEBits));
  }
  /** 
   * Decode the lower 13-bit of the input into @f$d\Delta@f$ and @f$dx@f$ 
   * 
   * @param bits   Encoded 13-bit word (lower 13 bit)
   * @param edep   On return, the @f$d\Delta@f$
   * @param length On return, the @f$dx@f$ 
   */
  static void Decode(UInt_t bits, Double_t& edep, Double_t& length)
  {
    edep   = DecodeOne((bits >> 0)       & kEMask, GetdEAxis());
    length = DecodeOne((bits >> kNEBits) & kLMask, GetdLAxis());
  }
  /** 
   * Decode the lower 13-bit of the input into @f$d\Delta@f$ and @f$dx@f$ 
   * 
   * @param bits    Encoded 13-bit word (lower 13 bit)
   * @param edep    On return, the @f$d\Delta@f$
   * @param length  On return, the @f$dx@f$ 
   * @param wEdep   On return, the width of the corresponding @f$d\Delta@f$ bin
   * @param wLength On return, the width of the corresponding @f$dx@f$ bin
   */
  static void Decode(UInt_t bits, Double_t& edep, Double_t& length, 
		     Double_t& wEdep, Double_t& wLength)
  {
    edep   = DecodeOne((bits >> 0)       & kEMask, wEdep,   GetdEAxis());
    length = DecodeOne((bits >> kNEBits) & kLMask, wLength, GetdLAxis());
  }    
  /** 
   * Make a 1-dimension histogram with the natural binning for the
   * encoding for either @f$d\Delta@f$ or @f$dx@f$
   * 
   * @param name   Name of produced histogram
   * @param title  Title of produced histogram 
   * @param mode   If 0, make histogram for @f$d\Delta@f$. If 1
   * for @f$dx@f$, if 2 for @f$d\Delta/dx@f$
   * 
   * @return Newly allocated histogram
   */
  static TH1* Make1D(const char* name, const char* title, UShort_t mode=true)
  {
    const Spec& a = (mode==0 || mode==2 ? GetdEAxis() : GetdLAxis());
    TArrayD     aa; a.FillBinArray(aa);

    if (mode == 2) 
      // In case we need to do dE/dx, extend the range by a factor 100
      for (Int_t i = 0; i < aa.GetSize(); i++) aa[i] *= 100;

    // Make the histogram 
    TH1* h = new TH1D(name, title, aa.GetSize()-1, aa.GetArray());
    h->SetXTitle(mode == 0 ? "d#Delta [MeV]" : 
		 mode == 1 ? "dx [cm]" : 
		 mode == 2 ? "d#Delta/dx [MeV/cm]" : "?");
    h->SetFillStyle(3001);
    h->SetMarkerStyle(20);
    h->Sumw2();
    
    return h;
  }
  /** 
   * Make a 2-dimension histogram with the natural binning for the
   * encoding of @f$d\Delta@f$ versus @f$dx@f$ (or vice versa)
   * 
   * @param name   Name of produced histogram
   * @param title  Title of produced histogram 
   * @param xedep  If true, put @f$d\Delta@f$ on the X-axis, otherwise
   * for @f$dx@f$
   * 
   * @return Newly allocated histogram
   */
  static TH2* Make2D(const char* name, const char* title, Bool_t xedep=true)
  {
    const Spec& a1 = (xedep ? GetdEAxis() : GetdLAxis());
    const Spec& a2 = (xedep ? GetdLAxis() : GetdEAxis());
    TArrayD     aa1; a1.FillBinArray(aa1);
    TArrayD     aa2; a2.FillBinArray(aa2);
    TH2* h = new TH2D(name, title, 
		      aa1.GetSize()-1, aa1.GetArray(),
		      aa2.GetSize()-1, aa2.GetArray());
    h->SetXTitle(xedep ? "d#Delta [MeV]" : "dx [cm]");
    h->SetYTitle(xedep ? "dx [cm]"       : "d#Delta [MeV]");
    return h;
  }
  /** 
   * Check if the encoded @f$d\Delta@f$ and @f$dx@f$ are available in
   * the upper 13 bits of the unique ID field of track references.
   * 
   * @param alirootRev AliROOT revision of the code that _produced_
   * the track references.  Note, it cannot be the revision of AliROOT
   * running, since that can be much newer (or older) than the code
   * that made the track references.  One should get this information
   * from the object stored in the ESD tree's user objects.
   * 
   * @return true if @a alirootRev is larger than some fixed number
   * set when this class was committed to SVN.
   */
  static Bool_t IsAvailable(UInt_t alirootRev) 
  {
    const UInt_t target = 60000;
    return alirootRev >= target;
  }
private:
  /** 
   * Encode one value 
   * 
   * @param v  Value 
   * @param a  Bin specification
   * 
   * @return Encoded value 
   */
  static UInt_t EncodeOne(Double_t v, const Spec& a) 
  {
    if (v < a.min) return 0;
    if (v > a.max) return a.nBins;
    return a.Encode(v);
  }
  /** 
   * Decode a value 
   * 
   * @param b   Encoded value 
   * @param a   Bin specification
   * 
   * @return Decoded value 
   */  
  static Double_t DecodeOne(UInt_t b, const Spec& a) 
  {
    if (b >= a.nBins) b = a.nBins-1;
    return a.Decode(b);
  }
  /** 
   * Decode a value 
   * 
   * @param b   Encoded value 
   * @param a   Bin specification
   * @param w   On return, the bin width 
   *
   * @return Decoded value 
   */  
  static Double_t DecodeOne(UInt_t b, Double_t& w, const Spec& a) 
  {
    if (b >= a.nBins) b = a.nBins-1;
    return a.Decode(b, w);
  }

  ClassDef(AliFMDEncodedEdx,1); // En-/Decode dE/dx for/from track references
};
#endif
// Local Variables:
//  mode: C++
// End:
