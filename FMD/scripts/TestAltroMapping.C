//____________________________________________________________________
//
// $id$
//
// Check integrety of Hardware2Detector and Detector2Hardware
//
/** @defgroup FMD_ALTRO_test ALTRO test
    @ingroup FMD_script 
*/
#ifndef __CINT__
# include <TString.h>
// # include <FMD/AliFMDParameters.h>
# include <FMD/AliFMDAltroMapping.h>
// # include <FMD/AliFMDUShortMap.h>
// # include <FMD/AliFMDBoolMap.h>
# include <AliLog.h>
# include <TError.h>
# include <iostream>
#endif
bool show_all=false;

//____________________________________________________________________
/** @ingroup FMD_ALTRO_test
    @param ddl 
    @param hwaddr 
    @return  */
const Char_t* 
Addr2Str(UInt_t ddl, UInt_t hwaddr, UShort_t timebin)
{
  static TString s;
  UInt_t board = (hwaddr >> 7) & 0x1F;
  UInt_t chip  = (hwaddr >> 4) & 0x7;
  UInt_t chan  = hwaddr & 0xF;
  s = Form("(0x%05X,0x%02X,0x%1X,0x%1X,%04d)", ddl, board, chip, chan, timebin);
  return s.Data();
}

//____________________________________________________________________
/** @ingroup FMD_ALTRO_test
    @param det 
    @param ring 
    @param sec 
    @param str 
    @return  */
const Char_t* 
Det2Str(UShort_t det, Char_t ring, UShort_t sec, UShort_t str, UShort_t sam)
{
  static TString s;
  s = Form("FMD%d%c[%2d,%3d]-%d", det, ring, sec, str, sam);
  return s.Data();
}

//____________________________________________________________________
/** @ingroup FMD_ALTRO_test
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
PrintTrans(UShort_t det, Char_t ring, UShort_t sec, Short_t str, UShort_t sam,
	   UInt_t ddl, UInt_t hwaddr, UShort_t timebin,
	   UShort_t odet, Char_t oring, UShort_t osec, Short_t ostr, 
	   UShort_t osam)
{
  static TString s1, s2, s3;
  s1 = Det2Str(det, ring, sec, str, sam);
  s2 = Addr2Str(ddl,hwaddr,timebin);
  s3 = Det2Str(odet, oring, osec, ostr, osam);
  Info("TestHWMap","%s -> %s -> %s", s1.Data(), s2.Data(), s3.Data());
}

//____________________________________________________________________
/** @ingroup FMD_ALTRO_test
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
CheckTrans(UShort_t det, Char_t ring, UShort_t sec, UShort_t str, UShort_t sam,
	   UInt_t ddl, UInt_t hwaddr, UShort_t timebin,
	   UShort_t odet, Char_t oring, UShort_t osec, UShort_t ostr,
	   UShort_t osam)
{
  bool ok = true;
  if (det != odet) {
    Warning("TestHWMap", "Detector # differ %d != %d", det, odet);
    ok = false;
  }
  if (ring != oring) {
    Warning("TestHWMap", "Ring Id differ %c != %c", ring, oring);  
    ok = false;
  }
  if (sec != osec) {
    Warning("TestHWMap", "Sector # differ %d != %d", sec, osec);
    ok = false;
  }
  if (str != ostr) { 
    ok = false;
    Warning("TestHWMap", "Strip # differ %d != %d", str, ostr);
  }
  if (sam != osam) { 
    ok = false;
    Warning("TestHWMap", "Sample # differ %d != %d", sam, osam);
  }

  if (!show_all) { 
    if (!ok) 
      PrintTrans(det,ring,sec,str,sam,
		 ddl,hwaddr,timebin,
		 odet,oring,osec,ostr,osam);
  }
}

//____________________________________________________________________
/** @ingroup FMD_ALTRO_test
 */
void
TestAltroMapping(bool sa=false, Int_t min=1, Int_t max=3)
{
  show_all = sa;
  // AliLog::SetModuleDebugLevel("FMD", 1);
  // if (min < 1 || min > 3) min = 1;
  if (max < min)          max = min;
  // AliFMDParameters* param = AliFMDParameters::Instance();
  AliFMDAltroMapping m;
  UShort_t presamp  = 19;
  UShort_t oversamp = 4;

  for (UShort_t det = min; det <= max; det++) {
    for (UShort_t rng = 0; rng < 2; rng++) {
      Char_t ring = (rng == 0 ? 'I' : 'O');
      Int_t  nsec = (ring == 'I' ?  20 :  40);
      Int_t  nstr = (ring == 'I' ? 512 : 256);
      for (UShort_t sec = 0; sec < nsec; sec++) {
	for (Short_t str = 0; str < nstr; str ++ /*= 128*/) {
	  for(UShort_t sam = 0; sam < oversamp; sam++) {
	    UShort_t ddl, hwaddr;
	    UShort_t timebin;
	    if (!m.Detector2Hardware(det, ring, sec, str, sam, 
				     presamp, oversamp, 
				     ddl, hwaddr, timebin)) {
	      Warning("TestHWMap", "detector to hardware failed on %s", 
		      Det2Str(det, ring, sec, str, sam));
	      continue;
	    }
	    UShort_t odet, osec, osam;
	    Short_t ostr;
	    Char_t   oring;
	    if (!m.Hardware2Detector(ddl, hwaddr, timebin, 
				     presamp, oversamp, 
				     odet, oring, osec, ostr, osam)){
	      Warning("TestHWMap", "hardware to detector failed on %s", 
		      Addr2Str(ddl, hwaddr, timebin));
	      continue;
	    }
	    if (show_all) 
	      PrintTrans(det,ring,sec,str,sam,
			 ddl,hwaddr,timebin,
			 odet,oring,osec,ostr,osam);
	    CheckTrans(det,ring,sec,str,sam,
		       ddl,hwaddr,timebin,
		       odet,oring,osec,ostr,osam);
	  }
	}// Loop over strips
      } // Loop over sectors
    } // Loop over rings
  } // Loop over detectors
}

  
//____________________________________________________________________
//
// EOF
//
    
	

	  
