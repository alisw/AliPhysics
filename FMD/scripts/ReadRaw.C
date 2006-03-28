/** @file    ReadRaw.C
    @author  Christian Holm Christensen <cholm@nbi.dk>
    @date    Tue Mar 28 12:39:08 2006
    @brief   Script to read raw data 
*/
/** @ingroup FMD_script
    @brief Read raw data into a TClonesArray - for testing 
 */
void
ReadRaw()
{
  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://$ALICE_ROOT");
  AliLog::SetModuleDebugLevel("FMD", 10);
  AliFMDParameters::Instance()->Init();
  AliRawReader* r = new AliRawReaderFile(0);
  AliFMDRawReader* fr = new AliFMDRawReader(r, 0);
  TClonesArray* a = new TClonesArray("AliFMDDigit");
  fr->ReadAdcs(a);
}
//____________________________________________________________________
//
// EOF
//
