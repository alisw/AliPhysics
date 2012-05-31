// $Id$
/**
 * @file changeRunValidity.C
 * @brief Change the run validity for OCDB entries. It increments the version of the input by 1,
 * in case there exists another object in the same folder with the same run validity. 
 * Please BE CAREFUL where you run it and where you copy the result. Check that the version of the
 * produced file is correct to enable its loading during the reconstruction.
 *
 * THIS IS A TEMPORARY SOLUTION UNTIL THE PENDOLINO STARTS RUNNING. WE SHOULD NOT RELY ON THIS MACRO.
 *
 * <pre>
 * Usage: aliroot -b -q -l changeRunValidity.C'("RunXXX_YYY_v?_s?.root",runmin,runmax)'
 * Default runmin = 0, runmax = 999999999
 *
 * Used for OCDB objects whose run validity is restrictic (run specific),
 * e.g. TPC/Calib/AltroConfig
 *      TPC/Calib/Temperature
 *      TPC/Calib/HighVoltage
 *      GRP/CTP/CTPtiming
 *
 *
 * @author Kalliopi.Kanaki@ift.uib.no
 * @ingroup alihlt_tutorial
 */
 
 void changeRunValidity(const char* file, int runmin=0, int runmax=999999999){

  TFile *f_in = TFile::Open(file);
  AliCDBEntry* entry = AliCDBEntry;

  AliCDBId id = entry->GetId();  
  printf("existing runmin: %d, runmax: %d, version: %d, subversion: %d\n", id.GetFirstRun(), id.GetLastRun(), id.GetVersion(), id.GetSubVersion() );  
  
  id.SetFirstRun(runmin);
  id.SetLastRun(runmax);  
  //id.SetRunRange(nnn,mmm); 
  //id.SetLastRun(AliCDBRunRange::Infinity()); 
 
  id.SetVersion(id.GetVersion()+1);
  entry->SetId(id);
  
  printf("change to runmin: %d, runmax: %d, version: %d, subversion: %d\n", id.GetFirstRun(), id.GetLastRun(), id.GetVersion(), id.GetSubVersion() );  
  
  
  TString out;
  out.Form("Run%i_%i_v%i_s%i.root",runmin, runmax, id.GetVersion(), id.GetSubVersion());
     
  TFile f_out(out, "RECREATE");
  entry->Write();
  f_out.Close(); 

}
