// $Id$
/**
 * @file rec-from-hltout.C
 * @brief Detector reconstruction from HLTOUT data
 *
 * Example macro to run detector reconstruction from the HLTOUT data
 * instead of the original detector data.
 *
 * Usage: aliroot -b -q rec-from-hltout.C | tee rec-from-hltout.log
 *
 * The macro asumes detector ddl raw data to be available in the HLTOUT,
 * e.g. simulated by sim-hlt-rawddl.C
 * Input replacement is run for the ITS (Note that only the ITSSDD was
 * published to HLTOUT)
 *
 * The HLTOUT handler implemented in libAliHLTSample.lib is used to
 * redirect data blocks of type {DDLRAW :SMPL} from the HLTOUT to the
 * reconstruction input. \b Note: the origin \em SMPL is a special
 * id just for the purpose of examples.
 *
 * @author Matthias.Richter@ift.uib.no
 * @ingroup alihlt_tutorial
 */
void rec_from_hltout()
{
  /////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////
  //
  // setup of the HLT system
  AliHLTSystem* pHLT=AliHLTPluginBase::GetInstance();
  if (!pHLT) {
    cerr << "fatal error: can not get HLT instance" << endl;
  }
  // load libAliHLTSample in order to make HLTOUT handler available 
  pHLT->LoadComponentLibraries("libAliHLTSample.so");
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
