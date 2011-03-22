//____________________________________________________________________
//
// $id$
//
// Check integrety of Hardware2Detector and Detector2Hardware
//
/** @defgroup HW_Test Hardware map test
    @ingroup FMD_script 
*/
//____________________________________________________________________
/** 
 * Convert hardware address to board,altro,channel
 * 
 * @param w Hardware address
 * @param b out, Board
 * @param a out, Altro 
 * @param c out, Channel
 */
void
HwAddr2Channel(UShort_t w, UShort_t& b, UShort_t& a, UShort_t& c)
{
  b = (w >> 7) & 0x1F;
  a = (w >> 4) & 0x7;
  c = w & 0xF;
}

//____________________________________________________________________
/** 
 * Convert board,altro,channel to hardware address
 * 
 * @param b Board
 * @param a Altro
 * @param c Channel
 * 
 * @return Hardware address
 */
UShort_t
Channel2HwAddr(UShort_t b, UShort_t a, UShort_t c)
{
  return (((b & 0x1F) << 7) | ((a & 0x7) << 4) | (c & 0xF));
}
//____________________________________________________________________
/** 
 * Format a hardware addresss
 * 
 * @ingroup HW_test
 * 
 * @param l  DDL ID
 * @param b  Board
 * @param a  Altro
 * @param c  Channel
 * @param tb Timebin
 * 
 * @return Formatted hardware address
 */
Char_t* 
Addr2Str(UInt_t l, UShort_t b, UShort_t a, UShort_t c, UShort_t tb)
{
  static TString s;
  s = Form("(0x%01X,0x%02X,0x%1X,0x%1X)-%4d", l, b, a, c, tb);
  return s.Data();
}

//____________________________________________________________________
/** 
 * Format a hardware addresss
 * 
 * @ingroup HW_test
 * 
 * @param l  DDL ID
 * @param w  Hardware address
 * @param tb Timebin
 * 
 * @return Formatted hardware address
 */
Char_t* 
Addr2Str(UShort_t l, UShort_t w, UShort_t tb)
{
  UShort_t b, a, c;
  HwAddr2Channel(w, b, a, c);
  return Addr2Str(l, b, a, c, tb);
}

//____________________________________________________________________
/** 
 * Format a detector address 
 * 
 * @ingroup HW_test
 * 
 * @param d Detector
 * @param r Ring 
 * @param s Sector  
 * @param t Strip
 * @param n Sample
 * 
 * @return 
 */
Char_t* 
Det2Str(UShort_t d, Char_t r, UShort_t s, UShort_t t, UShort_t n)
{
  static TString str;
  str = Form("FMD%d%c[%2d,%3d]-%d", d, r, s, t, n);
  return str.Data();
}

//____________________________________________________________________
/** 
 * Print full transformation
 * 
 * @ingroup HW_test
 *
 * @param d    Detector
 * @param r    Ring
 * @param s    Sector 
 * @param t    Strip
 * @param n    Sample
 * @param l    DDL ID
 * @param w    Hardware address
 * @param tb   Timebin
 * @param od   Detector output
 * @param oR   Ring output    
 * @param os   Sector output  
 * @param ot   Strip output   
 * @param on   Sample output  
 */
void 
PrintTrans(UShort_t d,   Char_t r,   UShort_t s,  UShort_t t,  UShort_t n, 
	   UShort_t l,   UShort_t w, UShort_t tb,
	   UShort_t od,  Char_t oR,  UShort_t os, Short_t  ot, UShort_t on)
{
  static TString s1, s2, s3;
  s1 = Det2Str(d, r, s, t, n);
  s2 = Addr2Str(l,w,tb);
  s3 = Det2Str(od, oR, os, ot, on);
  Info("TestHWMap","%s -> %s -> %s", s1.Data(), s2.Data(), s3.Data());
}
//____________________________________________________________________
/** 
 * Print full transformation
 * 
 * @ingroup HW_test
 *
 * @param d    Detector
 * @param r    Ring
 * @param s    Sector 
 * @param t    Strip
 * @param n    Sample
 * @param l    DDL ID
 * @param b    Board
 * @param a    ALTRO
 * @param c    Channel
 * @param tb   Timebin
 * @param od   Detector output
 * @param oR   Ring output    
 * @param os   Sector output  
 * @param ot   Strip output   
 * @param on   Sample output  
 */
void 
PrintTrans(UShort_t d,   Char_t   r,  UShort_t s,  UShort_t t,  UShort_t n, 
	   UShort_t l,   UShort_t b,  UShort_t a,  UShort_t c,  UShort_t tb,
	   UShort_t od,  Char_t   oR, UShort_t os, Short_t  ot, UShort_t on)
{
  static TString s1, s2, s3;
  s1 = Det2Str(d, r, s, t, n);
  s2 = Addr2Str(l,b,a,c,tb);
  s3 = Det2Str(od, oR, os, ot, on);
  Info("TestHWMap","%s -> %s -> %s", s1.Data(), s2.Data(), s3.Data());
}

//____________________________________________________________________
/** 
 * Check transformation
 * 
 * @ingroup HW_test
 *
 * @param d    Detector
 * @param r    Ring
 * @param s    Sector 
 * @param t    Strip
 * @param n    Sample
 * @param od   Detector output
 * @param oR   Ring output    
 * @param os   Sector output  
 * @param ot   Strip output   
 * @param on   Sample output  
 * 
 * @return true on success 
 */
Bool_t
CheckTrans(UShort_t d,  Char_t r,  UShort_t s,  UShort_t t,  UShort_t n, 
	   UShort_t od, Char_t oR, UShort_t os, Short_t  ot, UShort_t on)
{
  bool ok = true;
  TString what;
  if (d != od) {
    what.Append(Form("\n\tDetector # differ %d != %d", d, od));
    ok = false;
  }
  if (r != oR) {
    what.Append(Form("\n\tRing Id differ %c != %c", r, oR));  
    ok = false; 
  }
  if (s != os) {
    what.Append(Form("\n\tSector # differ %d != %d", s, os));
    ok = false;
  }
  if (t != ot) {
    what.Append(Form("\n\tStrip # differ %d != %d", (t / 128) * 128, ot));
    ok = false;
  }
  if (!ok) {
    static TString s1, s3;
    s1 = Det2Str(d, r, s, t);
    s3 = Det2Str(od, oR, os, ot);
    Warning("TestHWMap", "%s -> %s %s", s1.Data(), s3.Data(), what.Data());
  }
  return ok;
}

//____________________________________________________________________
/** 
 * Test hardware address map by converting from detector coordinates
 * to hardware address and then back again.
 * 
 * @ingroup HW_test
 *
 * @param useHwAddr Whether to use 12 bit hardware address, or
 *                  board,altro,channel 
 * @param few       Only do a few - 1 detector, one sector 
 */
void
TestHWMap(bool useHwAddr=false, bool few=false)
{
  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://${ALICE_ROOT}/OCDB");
  cdb->SetRun(0);

  AliFMDParameters* param = AliFMDParameters::Instance();
  param->Init(true, AliFMDParameters::kAltroMap);
  param->SetSampleRate(2);
  AliLog::SetModuleDebugLevel("FMD", 51);

  UInt_t ol = 0, ow = 0xFFFFFFFF;
  UShort_t nd = (few ? 1 : 3); 
  for (UShort_t d = 1; d <= nd; d++) {
    UShort_t nr = (d == 1 ? 1 : 2);
    for (UShort_t q = 0; q < nr; q++) {
      Char_t r  = (q ==   0 ? 'I' : 'O');
      Int_t  ns = (few ? 1 : (r == 'I' ?  20 :  40));
      Int_t  nt =            (r == 'I' ? 512 : 256);
      for (UShort_t s = 0; s < ns; s++) {
	for (UShort_t t = 0; t < nt; t++) {
	  Int_t nsam = param->GetSampleRate(d, r, s, t);
	  for (UShort_t n=0; n < nsam; n++) {
	    UShort_t l, b, a, c, w, tb;
	    bool ret = true;
	    if (useHwAddr) 
	      ret = param->Detector2Hardware(d, r, s, t, n, l, w, tb);
	    else {
	      ret = param->Detector2Hardware(d, r, s, t, n, l, b, a, c, tb);
	      w   = Channel2HwAddr(b, a, c);
	    }
	    if (!ret) {
	      Warning("TestHWMap", "detector to hardware failed on %s", 
		      Det2Str(d, r, s, t, n));
	      continue;
	    }
	    UShort_t od, os, on;
	    Short_t  ot;
	    Char_t   oR;
	    if (useHwAddr) 
	      ret = param->Hardware2Detector(l,w,tb,od,oR,os,ot,on);
	    else 
	      ret = param->Hardware2Detector(l,b,a,c,tb,od,oR,os,ot,on);
	    if (!ret) {
	      Warning("TestHWMap", "hardware to detector failed on %s", 
		      Addr2Str(l, w, tb));
	      continue;
	    }
	    bool print = few;
	    if (ol != l || ow != w) {
	      print = true;
	      ol    = l;
	      ow    = w;
	    }
	    if (!CheckTrans(d,r,s,t,n,od,oR,os,ot,on)) print = true;
	    if (print) {
	      if (useHwAddr)  
		PrintTrans(d,r,s,t,n,l,w,tb,od,oR,os,ot,on);
	      else
		PrintTrans(d,r,s,t,n,l,b,a,c,tb,od,oR,os,ot,on);
	    }
	  } // Loop over samples
	}// Loop over strips
      } // Loop over sectors
    } // Loop over rings
  } // Loop over detectors
}

  
//____________________________________________________________________
//
// EOF
//
    
	

	  
