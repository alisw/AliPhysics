#ifndef ALIFMDCALIBDRAWER_H
#define ALIFMDCALIBDRAWER_H
#include <TObject.h>
class TH2D;
class TH1D;
class TH1;

class AliFMDCalibDrawer : public TObject 
{
public:
  /**
   * Things we can draw 
   */
  enum EWhat { 
    kPedestal, 
    kNoise, 
    kGain, 
    kDead,
    kRate, 
    kRange, 
    kZeroSuppression 
  };

  AliFMDCalibDrawer() {}
  /** 
   * Initialize the drawer
   * 
   * @param runNo  Run number
   * @param ocdb   Source of parameters (if 0, do not set)
   */
  void Init(Int_t runNo, const char* ocdb=0);
  /** 
   * Draw pedestals for the selected region.
   * 
]   * @param d Detector number (if 0 or negative, draw all)
   * @param r Ring id (if null, draw all rings)
   * @param s Sector number (if negative, draw all sectors) 
   * @param t Strip number (if negative, draw all strips)
   */
  void DrawPedestals(Short_t d=-1, Char_t r='\0', 
		     Short_t s=-1, Short_t t=-1) const 
  { 
    DrawOne(kPedestal, d, r, s, t);
  }
  /** 
   * Draw noise for the selected region.
   * 
   * @param d Detector number (if 0 or negative, draw all)
   * @param r Ring id (if null, draw all rings)
   * @param s Sector number (if negative, draw all sectors) 
   * @param t Strip number (if negative, draw all strips)
   */
  void DrawNoise(Short_t d=-1, Char_t r='\0', 
		 Short_t s=-1, Short_t t=-1) const 
  {
    DrawOne(kNoise, d, r, s, t);
  }
  /** 
   * Draw gains for the selected region.
   * 
   * @param d Detector number (if 0 or negative, draw all)
   * @param r Ring id (if null, draw all rings)
   * @param s Sector number (if negative, draw all sectors) 
   * @param t Strip number (if negative, draw all strips)
   */
  void DrawGains(Short_t d=-1, Char_t r='\0', 
		     Short_t s=-1, Short_t t=-1) const
  {
    DrawOne(kGain, d, r, s, t);
  }
  /** 
   * Draw dead for the selected region.
   * 
   * @param d Detector number (if 0 or negative, draw all)
   * @param r Ring id (if null, draw all rings)
   * @param s Sector number (if negative, draw all sectors) 
   * @param t Strip number (if negative, draw all strips)
   */
  void DrawDead(Short_t d=-1, Char_t r='\0', Short_t s=-1, Short_t t=-1) const
  {
    DrawOne(kDead, d, r, s, t);
  }
  
  void DrawRates(Short_t d=-1, Char_t r='\0', Short_t s=-1, Short_t t=-1) const
  {
    DrawOne(kRate, d, r, s, t);
  }
  void DrawRanges(Short_t d=-1, Char_t r='\0', Short_t s=-1, Short_t t=-1) const
  {
    DrawOne(kRange, d, r, s, t);
  }
  void DrawThresholds(Short_t d=-1, Char_t r='\0', Short_t s=-1, 
		     Short_t t=-1) const
  {
    DrawOne(kZeroSuppression, d, r, s, t);
  }
  /** 
   * Draw one thing
   * 
   * @param what  What to draw 
   * @param d Detector number (if 0 or negative, draw all)
   * @param r Ring id (if null, draw all rings)
   * @param s Sector number (if negative, draw all sectors) 
   * @param t Strip number (if negative, draw all strips)
   */
  void  DrawOne(EWhat what, Short_t d=-1, Char_t r='\0', 
		Short_t s=-1, Short_t t=-1) const;
protected:
  void SetAttributes(TH1* ret, EWhat what, UShort_t d, Char_t r) const;
  Double_t GetHistMax(EWhat what) const;
  Double_t GetHistMin(EWhat what) const;
  /** 
   * Get the base histogram name
   * 
   * @param what What is drawn
   * 
   * @return Histogram base name
   */
  const char* GetHistName(EWhat what) const;
  /** 
   * Get the base histogram title
   * 
   * @param what What is drawn
   * 
   * @return Histogram base title
   */
  const char* GetHistTitle(EWhat what) const;
  /** 
   * Get the numbers 
   * 
   * @param what What to get
   * @param d    Detector number
   * @param r    Ring id
   * @param s    Sector number
   * @param t    Strip number
   * @param val  On return, the value 
   * @param err  On return, the error on value (or negative if no error)
   */
  void GetNumber(EWhat what, UShort_t d, Char_t r, UShort_t s, UShort_t t, 
		 Double_t& val, Double_t& err) const;
  /** 
   * Make a 1D histogram
   *  
   * @param what What to make the histogram for 
   * @param d    Detector 
   * @param r    Ring 
   * @param s    Sector
   * 
   * @return Newly allocated histogram
   */
  TH1D* Make1D(EWhat what, UShort_t d, Char_t r, UShort_t s) const;
  /** 
   * Make a 2D histogram
   *  
   * @param what What to make the histogram for 
   * @param d    Detector 
   * @param r    Ring 
   * 
   * @return Newly allocated histogram
   */
  TH2D* Make2D(EWhat what, UShort_t d, Char_t r) const;
  /** 
   * Fill a per-ring 2D histogram 
   * 
   * @param what What to fill in
   * @param d    Detector 
   * @param r    Ring 
   * 
   * @return Filled histogram 
   */
  TH1* FillRing(EWhat what, UShort_t d, Char_t r) const;
  /** 
   * Fill a per-sector 1D histogram 
   * 
   * @param what What to fill in
   * @param d    Detector
   * @param r    Ring 
   * @param s    Sector 
   * 
   * @return Filled histogram 
   */
  TH1* FillSector(EWhat what, UShort_t d, Char_t r, UShort_t s) const;

  Int_t GetRingColor(UShort_t d, Char_t r) const;
  ClassDef(AliFMDCalibDrawer,0); // Draw calibrations 
};
#endif
// Local Variables: 
//  mode: C++
// End:

