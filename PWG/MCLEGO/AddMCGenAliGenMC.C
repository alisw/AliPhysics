/************************************************************************************
 * Copyright (C) 2020, Copyright Holders of the ALICE Collaboration                 *
 * All rights reserved.                                                             *
 *                                                                                  *
 * Redistribution and use in source and binary forms, with or without               *
 * modification, are permitted provided that the following conditions are met:      *
 *     * Redistributions of source code must retain the above copyright             *
 *       notice, this list of conditions and the following disclaimer.              *
 *     * Redistributions in binary form must reproduce the above copyright          *
 *       notice, this list of conditions and the following disclaimer in the        *
 *       documentation and/or other materials provided with the distribution.       *
 *     * Neither the name of the <organization> nor the                             *
 *       names of its contributors may be used to endorse or promote products       *
 *       derived from this software without specific prior written permission.      *
 *                                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY              *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     *
 ************************************************************************************/
#include "AliGenerator.h"
#include "AliGenExtExec.h"

/**
 * @brief Macro adding event generation via aligenmc
 * 
 * Configuring the event generation process handled by aligenmc, where
 * aligenmc steers external event generators (pythia, herwig, sherpa, ...)
 * and events are passed to the aliroot process as HepMC events which are
 * read in via AliGenExtExec and distributed for analysis. Parameters to
 * be specified match the command line arguments which need to be provided
 * to aligenmc.
 * 
 * @param generator event generator called by aligenmc
 * @param package package needed to be loaded in order to run the event generator
 * @param aligenmc_version aligenmc package
 * @param tune generator tune
 * @param energy centre-of-mass energy
 * @param pthardmin min. pt hard (for hard QCD processes)
 * @param pthardmax max. pt hard (for hard QCD processes)
 * @param nevents Number of events (should only be set if the macro has no acces to the grid plugin)
 * @return AliGenExtExec steering aligenmc process, nullptr if configuration failed 
 */
AliGenerator *AddMCGenAliGenMC(const char *generator,
                               const char *package = "",
                               const char *aligenmc_version = "",
                               const char *tune = "",
                               Int_t energy = 13000.,
                               Int_t pthardmin = -1.,
                               Int_t pthardmax = -1.,
                               Int_t nevents = -1
                               ) {
  if(!strlen(generator)) {
    std::cerr << "AddMCGenAliGenMC: Generator needs to be specified" << std::endl;
    return nullptr;
  }

  // Handling of the number of events in case of aligenmc
  // aligenmc simulates the number of events defined when calling. In order
  // to have the number consistent between aligenmc and AliRoot the same number of
  // events must be set. In the MCgen train the number of events per job is defined
  // in the variable "SPLIT_MAX_INPUT_FILE_NUMBER". The number is specified in the train
  // env and accessible at train generation time (when the add macro is called). Therefore 
  // the number is passed to the gen script via the add macro. The macro is intended to be
  // used in a LEGO train environment. Therefore we assume that the number of events is 
  // already configured in the analysis manager, so we can take it from there. In case the macro
  // is used without the LEGO train mode the number of events to be generated must be set 
  // explicitly by the user via the argument nevents.  
  int nmcevents = 0;
  if(nevents > 0) {
    nmcevents = nevents;
  } else {
    auto analysismgr = AliAnalysisManager::GetAnalysisManager();
    nmcevents = analysismgr->GetNMCevents(); 
  }

  if(!nmcevents) {
    std::cerr << "Number of events 0, please set it via the AliAnalysisAlien plugin or the argument nevents" << std::endl;
    std::cerr << "Generator cannot be created" << std::endl;
    return nullptr;
  }

  AliGenExtExec *genhandler = new AliGenExtExec(Form("$ALICE_PHYSICS/PWG/MCLEGO/ALIGENMC/gen.sh"));
  TString argstring = Form("-g=%s -n=%d", generator, nmcevents);
  if(strlen(package)) {
    argstring += Form(" -p=%s -e=%d", package, energy);
  }
  if(strlen(aligenmc_version)) {
    argstring += Form(" -a=%s", aligenmc_version);
  }
  if(strlen(tune)) {
    argstring += Form(" -t=%s", tune);
  }
  if(pthardmin >= 0) {
    argstring += Form(" --pthardmin=%d", pthardmin);
  }
  if(pthardmax > 0) {
    argstring += Form(" --pthardmax=%d", pthardmax);
  }
  genhandler->SetGeneratorOptionalArguments(argstring);
  return genhandler;
}
