//____________________________________________________________________
//
// $id$
//
// Check integrety of Hardware2Detector and Detector2Hardware
//
/** @defgroup ALTRO_test ALTRO test
    @ingroup FMD_script 
*/
//____________________________________________________________________
/** @ingroup ALTRO_test
    @param ddl 
    @param hwaddr 
    @return  */
Char_t* 
Addr2Str(UInt_t ddl, UInt_t hwaddr)
{
  static TString s;
  UInt_t board = (hwaddr >> 7) & 0x1F;
  UInt_t chip  = (hwaddr >> 4) & 0x7;
  UInt_t chan  = hwaddr & 0xF;
  s = Form("(0x%05X,0x%02X,0x%1X,0x%1X)", ddl, board, chip, chan);
  return s.Data();
}

//____________________________________________________________________
/** @ingroup ALTRO_test
    @param det 
    @param ring 
    @param sec 
    @param str 
    @return  */
Char_t* 
Det2Str(UShort_t det, Char_t ring, UShort_t sec, UShort_t str)
{
  static TString s;
  s = Form("FMD%d%c[%2d,%3d]", det, ring, sec, str);
  return s.Data();
}

//____________________________________________________________________
/** @ingroup ALTRO_test
    @param det 
    @param ring 
    @param sec 
    @param str 
    @param ddl 
    @param hwaddr 
    @param odet 
    @param oring 
    @param osec 
    @param ostr 
*/
void 
PrintTrans(UShort_t det, Char_t ring, UShort_t sec, UShort_t str, 
	   UInt_t ddl, UInt_t hwaddr,
	   UShort_t odet, Char_t oring, UShort_t osec, UShort_t ostr)
{
  static TString s1, s2, s3;
  s1 = Det2Str(det, ring, sec, str);
  s2 = Addr2Str(ddl,hwaddr);
  s3 = Det2Str(odet, oring, osec, ostr);
  Info("TestHWMap","%s -> %s -> %s", s1.Data(), s2.Data(), s3.Data());
}

//____________________________________________________________________
/** @ingroup ALTRO_test
    @param det 
    @param ring 
    @param sec 
    @param str 
    @param odet 
    @param oring 
    @param osec 
    @param ostr 
*/
void
CheckTrans(UShort_t det, Char_t ring, UShort_t sec, UShort_t str, 
	   UShort_t odet, Char_t oring, UShort_t osec, UShort_t ostr)
{
  if (det != odet) 
    Warning("TestHWMap", "Detector # differ %d != %d", det, odet);
  if (ring != oring) 
    Warning("TestHWMap", "Ring Id differ %c != %c", ring, oring);  
  if (sec != osec) 
    Warning("TestHWMap", "Sector # differ %d != %d", sec, osec);
  if (str != ostr) 
    Warning("TestHWMap", "Strip # differ %d != %d", str, ostr);
}

//____________________________________________________________________
/** @ingroup ALTRO_test
 */
void
TestAltroMapping(Int_t min=2, Int_t max=0)
{
  // if (min < 1 || min > 3) min = 1;
  if (max < min)          max = min;
  AliFMDParameters* param = AliFMDParameters::Instance();
  AliFMDAltroMapping m;
  
  for (UShort_t det = min; det <= max; det++) {
    for (UShort_t rng = 0; rng < 2; rng++) {
      Char_t ring = (rng == 0 ? 'I' : 'O');
      Int_t  nsec = (ring == 'I' ?  20 :  40);
      Int_t  nstr = (ring == 'I' ? 512 : 256);
      for (UShort_t sec = 0; sec < nsec; sec++) {
	for (UShort_t str = 0; str < nstr; str += 128) {
	  UInt_t ddl, hwaddr;
	  if (!m.Detector2Hardware(det, ring, sec, str, ddl, hwaddr)) {
	    Warning("TestHWMap", "detector to hardware failed on %s", 
		    Det2Str(det, ring, sec, str));
	    continue;
	  }
	  UShort_t odet, osec, ostr;
	  Char_t   oring;
	  if (!m.Hardware2Detector(ddl, hwaddr, odet, oring, osec, ostr)){
	    Warning("TestHWMap", "hardware to detector failed on %s", 
		    Addr2Str(ddl, hwaddr));
	    continue;
	  }
	  PrintTrans(det,ring,sec,str,ddl,hwaddr,odet,oring,osec,ostr);
	  CheckTrans(det,ring,sec,str,odet,oring,osec,ostr);
	}// Loop over strips
      } // Loop over sectors
    } // Loop over rings
  } // Loop over detectors
}

  
//____________________________________________________________________
//
// EOF
//
    
	

	  
