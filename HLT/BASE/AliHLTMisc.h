// @(#) $Id$

#ifndef ALIHLTMISC_H
#define ALIHLTMISC_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               */

/** @file   AliHLTMisc.h
    @author Matthias Richter
    @date   
    @brief  Definition of various glue functions implemented in dynamically
            loaded libraries
*/

#define ALIHLTMISC_LIBRARY "libHLTrec.so"
#define ALIHLTMISC_INIT_CDB "AliHLTMiscInitCDB"
#define ALIHLTMISC_SET_CDB_RUNNO "AliHLTMiscSetCDBRunNo"

#ifdef __cplusplus
extern "C" {
#endif

  /**
   * Init the CDB access for the running instance.
   * The method is used from the C wrapper interface utilized by the  on-line
   * framework. The path of the (H)CDB is set to the specified path.<br>
   * When running from AliRoot, the CDB path is set in the startup of the
   * reconstruction.<br>
   * If cdbpath is nil or empty and the CDB is not already initialized, the
   * CDB storage is set to local://$ALICE_ROOT/OCDB and the run no to 0.
   * @param cdbpath     path to the CDB
   * @return neg. error code if failed
   * @note function implemented in libHLTrec
   */
  int AliHLTMiscInitCDB(const char* cdbpath);
  typedef int (*AliHLTMiscInitCDB_t)(const char* cdbpath);

  /**
   * Init the Run no for the CDB access.
   * @param runNo       the run no
   * @return neg. error code if failed
   * @note function implemented in libHLTrec
   */
  int AliHLTMiscSetCDBRunNo(int runNo);
  typedef int (*AliHLTMiscSetCDBRunNo_t)(int runNo);

#ifdef __cplusplus
}
#endif
#endif //ALIHLTMISC_H
