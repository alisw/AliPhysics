/** 
 * 
 * 
 * @ingroup pwg2_forward_analysis_scripts_tests
 */
namespace {
  enum { 
    kSolid        = 0x000, 
    kHollow       = 0x001, 
    kCircle       = 0x002,
    kSquare       = 0x004, 
    kUpTriangle   = 0x006, 
    kDownTriangle = 0x008, 
    kDiamond      = 0x00a,
    kCross        = 0x00c,
    kStar         = 0x00e
  };
  /** 
   * 
   * 
   * @param bits 
   * 
   * @return 
   * @ingroup pwg2_forward_analysis_scripts_tests
   */
  Int_t MarkerStyle(UInt_t bits)
  {
    Int_t  base   = bits & (0xFE);
    Bool_t hollow = bits & kHollow;
    switch (base) { 
    case kCircle:       return (hollow ? 24 : 20);
    case kSquare:       return (hollow ? 25 : 21);
    case kUpTriangle:   return (hollow ? 26 : 22);
    case kDownTriangle: return (hollow ? 32 : 23);
    case kDiamond:      return (hollow ? 27 : 33); 
    case kCross:        return (hollow ? 28 : 34); 
    case kStar:         return (hollow ? 30 : 29); 
    }
    return 1;
  }
  /** 
   * 
   * 
   * @param style 
   * 
   * @return 
   * @ingroup pwg2_forward_analysis_scripts_tests
   */
  UShort_t MarkerBits(Int_t style) 
  { 
    UShort_t bits = 0;
    switch (style) { 
    case 24: case 25: case 26: case 27: case 28: case 30: case 32: 
      bits |= kHollow; break;
    }
    switch (style) { 
    case 20: case 24: bits |= kCircle;       break;
    case 21: case 25: bits |= kSquare;       break;
    case 22: case 26: bits |= kUpTriangle;   break;
    case 23: case 32: bits |= kDownTriangle; break;
    case 27: case 33: bits |= kDiamond;      break;
    case 28: case 34: bits |= kCross;        break;
    case 29: case 30: bits |= kStar;         break;
    }
    return bits;
  }
  /** 
   * 
   * 
   * @param style 
   * 
   * @return 
   * @ingroup pwg2_forward_analysis_scripts_tests
   */
  Int_t FlipHollow(Int_t style) 
  {
    UShort_t bits = MarkerBits(style);
    Int_t ret = MarkerStyle(bits ^ kHollow);
    Info("FlipHollow", "style=%2d -> bits=0x%02x -> mask=0x%02x -> ret=%02d", 
	 style, bits, (bits ^ kHollow), ret);
    return ret;
  }
}

/** 
 * 
 * 
 * @param what 
 * @param base 
 * @param y 
 * @ingroup pwg2_forward_analysis_scripts_tests
 */
void DrawOne(const char* what, UShort_t base, Double_t y)
{
  TLatex* l = new TLatex(.07, y, what);
  l->SetTextSize(0.03);
  l->Draw();

  Int_t filled = MarkerStyle(base);
  // Info("DrawOne", "%2d (%16s) -> %d", base, what, style);
  TMarker* p = new TMarker(.35, y, filled);
  p->SetMarkerSize(1.5);
  p->Draw();

  Int_t hollow = MarkerStyle(base|kHollow);
  p = new TMarker(.60, y, hollow);
  p->SetMarkerSize(1.5);
  p->Draw();

  p = new TMarker(.75, y, FlipHollow(filled));
  p->SetMarkerSize(1.5);
  p->Draw();

  p = new TMarker(.85, y, FlipHollow(hollow));
  p->SetMarkerSize(1.5);
  p->Draw();
    
}

//
// EOF
//
