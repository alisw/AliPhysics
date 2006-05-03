//____________________________________________________________________
//
// $Id$
//
// Make fake alignment data.
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
PrintCalibration()
{
  AliCDBManager* cdb   = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://$ALICE_ROOT");
  cdb->SetRun(0);
  AliFMDParameters* p = AliFMDParameters::Instance();
  p->Init();
  p->Print("fmd1I[0,0]");
  p->Draw("pedestal");
}
//____________________________________________________________________
//
// EOF
//
