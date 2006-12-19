//____________________________________________________________________
//
// $Id$
//
// Print calibration constants
//
/** @file    PrintCalibration.C
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Sun Mar 26 17:58:27 2006
    @brief   Make fake calibration data.
*/
/** Make fake calibration data 
    @ingroup simple_script
 */
void
PrintCalibration(Int_t r=0, const char* what="gain")
{
  AliCDBManager* cdb   = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://$ALICE_ROOT");
  cdb->SetRun(r);
  AliFMDParameters* p = AliFMDParameters::Instance();
  p->Init(kTRUE);
  p->Print("fmd3*[8,0]");
  // p->Draw(what);
}
//____________________________________________________________________
//
// EOF
//
