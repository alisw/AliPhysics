/**
 * A script that will run the most simple chain possible. 
 *
 * What is done, depends on a bit mask:
 *
 *  Bit | Job
 *  ----+-----------------------
 *   0  | Make fake hits
 *   1  | Make digits from hits
 *   2  | Make sdigits from hits
 *   3  | Make raw DDL files from digits
 *   4  | Make ESD from raw ddl files
 *   5  | Display the data
 *
 *  So a mask of 0x3f means do everything, while 0x2 means
 *  only make digits from hits, and 0x1F is all but display.
 * 
 */

const char* path = "$ALICE_ROOT/FMD/scripts/";

void
RunAndDisplay(const char* script, bool show, const char* display)
{
  gROOT->Macro(Form("%s/%s", path, script));
  if (show) 
    gSystem->Exec(Form("aliroot -q -l %s/%s &", path, display));
}

void
RunSimpleChain(Int_t what=0x3b)
{
  if (what == 0) { 
    std::cout << "Usage: RunSimpleChain(MASK)\n\n"
	      << "MASK is a bit mask of things to do:\n\n"
	      << "\tBit | Job\n"
	      << "\t----+-----------------------\n"
	      << "\t 0  | Make fake hits\n"
	      << "\t 1  | Make digits from hits\n"
	      << "\t 2  | Make sdigits from hits\n"
	      << "\t 3  | Make raw DDL files from digits\n"
	      << "\t 4  | Make ESD from raw ddl files\n"
	      << "\t 5  | Display the data\n\n"
	      << "So a mask of 0x3f means do everything, while 0x2 means\n"
	      << "only make digits from hits, and 0x1F is all but display.\n"
	      << "\n" << std::endl;
    return;
  }
  Bool_t display = (what & (1 << 5));
  
  if (what & (1 << 0)) { 
    std::cout << "Making fake hits" << std::endl;
    RunAndDisplay("MakeFakeHits.C", display, "PatternHits.C");
  }
  
  if (what & (1 << 1)) { 
    std::cout << "Converting hits to digits" << std::endl;
    RunAndDisplay("Hits2Digits.C", display, "PatternDigits.C");
  }
  
  if (what & (1 << 2)) { 
    std::cout << "Converting hits to summable digits" << std::endl;
    RunAndDisplay("Hits2SDigits.C", display, "PatternSDigits.C");
  }
  
  if (what & (1 << 3)) { 
    std::cout << "Converting digits to raw data" << std::endl;
    gSystem->Exec("rm -rf raw0");
    RunAndDisplay("Digits2Raw.C", false, "");

    std::cout << "Moving ddl files" << std::endl;
    gSystem->mkdir("raw0");
    gSystem->Exec("mv FMD_*.ddl raw0/");
    gSystem->Exec("touch raw0/run0");
    if (display)
      gSystem->Exec(Form("aliroot -q %s/PatternRaw.C\\(\\\"\\\"\\) &", path));
  }
  
  if (what & (1 << 4)) {
    std::cout << "Making ESD from raw data" << std::endl;
    RunAndDisplay("Raw2ESD.C", display, "PatternESD.C");
  }
  
  std::cout << "All done" << std::endl;
}

  
//
// EOF
//
