// $Id$
/**
 * @file rec-sdd-from-hltout.C
 * @brief SDD reconstruction from HLTOUT data
 *
 * Example macro to run SDD reconstruction from the HLTOUT data
 * instead of the original detector data.
 *
 * Usage: aliroot -b -q rec-sdd-from-hltout.C | tee rec-sdd-from-hltout.log
 *
 * The macro asumes SDD ddl raw data to be available in the HLTOUT,
 * e.g. simulated by sim-sddraw-hltout.C this publishes ITSSDD raw data
 * into HLTOUT. Input replacement is run for the ITS.
 *
 * The HLTOUT handler implemented in libAliHLTITS.lib is used to
 * redirect data blocks of type {DDLRAW :ITS } from the HLTOUT to the
 * reconstruction input.
 *
 * @author Matthias.Richter@ift.uib.no
 * @ingroup alihlt_its
 */
void rec_sdd_from_hltout()
{
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  //
  // setup of the HLT system
  AliHLTSystem* pHLT=AliHLTPluginBase::GetInstance();
  if (!pHLT) {
    cerr << "fatal error: can not get HLT instance" << endl;
  }
  // disable chains from the library agent, not of interest for
  // this example of data redirection
  pHLT->ScanOptions("chains=");

  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  //
  // setup of the reconstruction
  AliReconstruction rec;
  rec.SetInput("./");
  rec.SetRunLocalReconstruction("ITS HLT");
  rec.SetRunTracking("ITS");
  rec.SetFillESD("ITS HLT");
  rec.SetRunQA(":");
  rec.SetFillTriggerESD(kFALSE);

  // set the redirection for ITS SDD
  rec.SetUseHLTData("ITSSDD");

  rec.Run();
}
