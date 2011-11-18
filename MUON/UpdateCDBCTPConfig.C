/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

/// \ingroup macros
/// \file UpdateCDBCTPConfig.C
/// \brief The macro writes a new GRP/CTP/Config trigger configuration file
/// from the GRP/CTP/MUON.cfg trigger descriptor, corresponding to a
/// run where only MUON is used as a trigger detector.
/// 
/// The macro writes a new GRP/CTP/Config trigger configuration file
/// from the GRP/CTP/MUON.cfg trigger descriptor, corresponding to a
/// run where only MUON is used as a trigger detector. The
/// compatibility is check against the Config.C used in the simulation
/// (which means that the MUON detector must be there). When "check"
/// is "true", the macro only shows the last CTP configuration stored in
/// the GRP.
/// This is necessary because at the first step of the simulation (digits)
/// the trigger configuration is guessed from the detectors defined in the
/// Config.C, while the reconstruction is performed in a second separate
/// step, having no more knowledge of the Config.C file.
/// This has to be done before starting the simulations, only once after
/// the installation of AliRoot:
/// 
/// <pre> 
///.L $ALICE_ROOT/MUON/UpdateCDBCTPConfig.C+
/// UpdateCDBCTPConfig(1);    // just checking
/// UpdateCDBCTPConfig();     // update the GRP/CDB entry
/// </pre>
///
/// AliRoot comes with a default entry corresponding to a pp trigger. In
/// this case, at the reconstruction phase error messages will appear (without
/// breaking the reconstruction): \n
///
/// E-AliCentralTrigger::CheckTriggeredDetectors: Wrong cluster mask from trigger
/// classes (7ffff), expecting (20c00)! Loaded trigger configuration is possibly wrong!
///
/// \author B. Vulpescu, LPC Clermont-Ferrand

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "ARVersion.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliCDBId.h"
#include "AliCDBMetaData.h"
#include "AliTriggerConfiguration.h"
#include "AliTriggerUtils.h"
#include "AliSimulation.h"
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#endif

void UpdateCDBCTPConfig(Bool_t check = false) {
  
  // AliSimulation object must exist, as it is used via AliMC
  // which is used in AliTriggerUtils::CheckConfiguration()
  AliSimulation sim;

  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdb->SetRun(0);

  if (check) {
    // current entry in GRP/CTP/Config trigger configuration
    AliCDBEntry *entry;
    entry = cdb->Get("GRP/CTP/Config");
    AliCDBMetaData *md = entry->GetMetaData();
    Printf("AliRoot version: %s",md->GetAliRootVersion());
    Printf("Comment: %s ",md->GetComment());
    Printf("Responsible: %s ",md->GetResponsible());
    return;
  }

  const Char_t* alice = gSystem->Getenv("ALICE_ROOT");

  // construct the CTP configuration starting from GRP/CTP/<CTPcfg>.cfg file

  // Config.C detector configuration
  TString cfgFile(Form("%s/MUON/Config.C",alice));

  // MUON.cfg trigger configuration
  TString cfgCTP(Form("%s/GRP/CTP/MUON.cfg",alice));

  AliTriggerConfiguration *trconfig = AliTriggerConfiguration::LoadConfiguration(cfgCTP);
  if (!trconfig) {
    Printf("Invalid cfg file! Exiting...");
    return;
  }

  // check if Config.C is compatible with the trigger configuration requested
  AliTriggerUtils tru;
  if (!tru.CheckConfiguration(cfgFile,trconfig)) {
    Printf("CTP configuration is incompatible with the specified Config.C and AliRoot version! Exiting...");
    return;
  }

  // put the new trigger configuration "trconfig" in the GRP/CTP/Config

  AliCDBId id("GRP/CTP/Config",0,AliCDBRunRange::Infinity());
  AliCDBMetaData *md= new AliCDBMetaData();

  // ROOT and AliRoot versions
  const char* rootv = gROOT->GetVersion();
  TString av(ALIROOT_SVN_BRANCH);
  Int_t revnum = ALIROOT_SVN_REVISION;

  Printf("root version: %s.  AliRoot %s, revision number %d",rootv,av.Data(),revnum);

  md->SetAliRootVersion(av.Data());
  md->SetComment(Form("Default CTP configuration for MUON mode produced with root version %s and AliRoot version %s revision %d ",rootv,av.Data(),revnum));

  AliCDBStorage* storage = cdb->GetStorage("local://$ALICE_ROOT/OCDB");
  storage->Put(trconfig,id,md);
  
}
