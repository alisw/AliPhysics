//************************************************************************
// Copyright(c) 2004, ALICE Experiment at CERN, All rights reserved. *
//                                                                        *
// Author: The ALICE Off-line Project.                                    *
// Contributors are mentioned in the code where appropriate.              *
//                                                                        *
// Permission to use, copy, modify and distribute this software and its   *
// documentation strictly for non-commercial purposes is hereby granted   *
// without fee, provided that the above copyright notice appears in all   *
// copies and that both the copyright notice and this permission notice   *
// appear in the supporting documentation. The authors make no claims     *
// about the suitability of this software for any purpose. It is          *
// provided "as is" without express or implied warranty.                  *
//************************************************************************/
// $Id$ */
/**
 * @file    AliFMDDigitizer.cxx
 * 
 * @author  Christian Holm Christensen <cholm@nbi.dk>
 * @date    Mon Mar 27 12:38:26 2006
 * @brief   FMD Digitizers implementation
 *
 * @ingroup FMD_sim
 */
//////////////////////////////////////////////////////////////////////////////
//
//  This class contains the procedures simulation ADC  signal for the
//  Forward Multiplicity detector  : SDigits->Digits
// 
//  Digits consists of
//   - Detector #
//   - Ring ID                                             
//   - Sector #     
//   - Strip #
//   - ADC count in this channel                                  
//
//  Digits consists of
//   - Detector #
//   - Ring ID                                             
//   - Sector #     
//   - Strip #
//   - Total energy deposited in the strip
//   - ADC count in this channel                                  
//
// As the Digits and SDigits have so much in common, the classes
// AliFMDDigitizer and AliFMDSDigitizer are implemented via a base
// class AliFMDBaseDigitizer.
//
//                 +---------------------+
//                 | AliFMDBaseDigitizer |
//                 +---------------------+
//                           ^
//                           |
//                +----------+---------+
//                |                    |
//      +-----------------+     +------------------+
//      | AliFMDDigitizer |	| AliFMDSDigitizer |
//      +-----------------+	+------------------+
//                |
//     +-------------------+
//     | AliFMDSSDigitizer |
//     +-------------------+
//
// These classes has several paramters: 
//
//     fPedestal
//     fPedestalWidth
//         (Only AliFMDDigitizer)
//         Mean and width of the pedestal.  The pedestal is simulated
//         by a Guassian, but derived classes my override MakePedestal
//         to simulate it differently (or pick it up from a database).
//
//     fVA1MipRange
//         The dymamic MIP range of the VA1_ALICE pre-amplifier chip 
//
//     fAltroChannelSize
//         The largest number plus one that can be stored in one
//         channel in one time step in the ALTRO ADC chip. 
//
//     fSampleRate
//         How many times the ALTRO ADC chip samples the VA1_ALICE
//         pre-amplifier signal.   The VA1_ALICE chip is read-out at
//         10MHz, while it's possible to drive the ALTRO chip at
//         25MHz.  That means, that the ALTRO chip can have time to
//         sample each VA1_ALICE signal up to 2 times.  Although it's
//         not certain this feature will be used in the production,
//         we'd like have the option, and so it should be reflected in
//         the code.
//
//
// The shaping function of the VA1_ALICE is generally given by 
//
//      f(x) = A(1 - exp(-Bx))
//
// where A is the total charge collected in the pre-amp., and B is a
// paramter that depends on the shaping time of the VA1_ALICE circut.
// 
// When simulating the shaping function of the VA1_ALICe
// pre-amp. chip, we have to take into account, that the shaping
// function depends on the previous value of read from the pre-amp. 
//
// That results in the following algorithm:
//
//    last = 0;
//    FOR charge IN pre-amp. charge train DO 
//      IF last < charge THEN 
//        f(t) = (charge - last) * (1 - exp(-B * t)) + last
//      ELSE
//        f(t) = (last - charge) * exp(-B * t) + charge)
//      ENDIF
//      FOR i IN # samples DO 
//        adc_i = f(i / (# samples))
//      DONE
//      last = charge
//   DONE
//
// Here, 
//
//   pre-amp. charge train 
//       is a series of 128 charges read from the VA1_ALICE chip
//
//   # samples
//       is the number of times the ALTRO ADC samples each of the 128
//       charges from the pre-amp. 
//
// Where Q is the total charge collected by the VA1_ALICE
// pre-amplifier.   Q is then given by 
//
//           E S 
//      Q =  - -
//           e R
//
// where E is the total energy deposited in a silicon strip, R is the
// dynamic range of the VA1_ALICE pre-amp (fVA1MipRange), e is the
// energy deposited by a single MIP, and S ALTRO channel size in each
// time step (fAltroChannelSize).  
//
// The energy deposited per MIP is given by 
//
//      e = M * rho * w 
//
// where M is the universal number 1.664, rho is the density of
// silicon, and w is the depth of the silicon sensor. 
//
// The final ADC count is given by 
//
//      C' = C + P
//
// where P is the (randomized) pedestal (see MakePedestal)
//
// This class uses the class template AliFMDMap<Type> to make an
// internal cache of the energy deposted of the hits.  The class
// template is instantasized as 
//
//  typedef AliFMDMap<std::pair<Float_t, UShort_t> > AliFMDEdepMap;
//
// The first member of the values is the summed energy deposition in a
// given strip, while the second member of the values is the number of
// hits in a given strip.  Using the second member, it's possible to
// do some checks on just how many times a strip got hit, and what
// kind of error we get in our reconstructed hits.  Note, that this
// information is currently not written to the digits tree.  I think a
// QA (Quality Assurance) digit tree is better suited for that task.
// However, the information is there to be used in the future. 
//
//
// Latest changes by Christian Holm Christensen
//
//////////////////////////////////////////////////////////////////////////////

//      /1
//      |           A(-1 + B + exp(-B))
//      | f(x) dx = ------------------- = 1
//      |                    B
//      / 0
//
// and B is the a parameter defined by the shaping time (fShapingTime).  
//
// Solving the above equation, for A gives
//
//                 B
//      A = ----------------
//          -1 + B + exp(-B)
//
// So, if we define the function g: [0,1] -> [0:1] by 
//
//               / v
//               |              Bu + exp(-Bu) - Bv - exp(-Bv) 
//      g(u,v) = | f(x) dx = -A -----------------------------
//               |                            B
//               / u
//
// we can evaluate the ALTRO sample of the VA1_ALICE pre-amp between
// any two times (u, v), by 
//       
//
//                                B	    Bu + exp(-Bu) - Bv - exp(-Bv)
//      C = Q g(u,v) = - Q ---------------- -----------------------------
//		           -1 + B + exp(-B)              B	            
//
//               Bu + exp(-Bu) - Bv - exp(-Bv) 
//        = -  Q -----------------------------
//                    -1 + B + exp(-B)
//

#include <TTree.h>		// ROOT_TTree
#include <TFile.h>
#include "AliFMDDebug.h" 	// Better debug macros
#include "AliFMDDigitizer.h"	// ALIFMDSSDIGITIZER_H
#include "AliFMD.h"		// ALIFMD_H
#include "AliFMDSDigit.h"	// ALIFMDDIGIT_H
#include "AliFMDDigit.h"	// ALIFMDDIGIT_H
#include "AliFMDParameters.h"   // ALIFMDPARAMETERS_H
#include <AliDigitizationInput.h>	// ALIRUNDIGITIZER_H
#include <AliRun.h>		// ALIRUN_H
#include <AliLoader.h>		// ALILOADER_H
#include <AliRunLoader.h>	// ALIRUNLOADER_H
    
//====================================================================
ClassImp(AliFMDDigitizer)
#if 0
;
#endif

//____________________________________________________________________
Bool_t
AliFMDDigitizer::Init()
{
  // 
  // Initialisation
  // 
  if (!AliFMDBaseDigitizer::Init()) return kFALSE;
  
#if 0
  // Get the AliRun object 
  AliRun* run = fRunLoader->GetAliRun();
  if (!run) { 
    AliWarning("Loading gAlice");
    fRunLoader->LoadgAlice();
    if (!run) { 
      AliError("Can not get Run from Run Loader");
      return kFALSE;
    }
  }
  
  // Get the AliFMD object 
  fFMD = static_cast<AliFMD*>(run->GetDetector("FMD"));
  if (!fFMD) {
    AliError("Can not get FMD from gAlice");
    return kFALSE;
  }  
#endif
  return kTRUE;
}


//____________________________________________________________________
void
AliFMDDigitizer::Digitize(Option_t*)
{
  // 
  // Execute this digitizer.  
  // This member function will be called once per event by the passed
  // AliDigitizationInput* digInput object. 
  // 
  // Parameters:
  //    options Not used 
  //
  if (!fDigInput) { 
    AliError("No digitisation input defined");
    return;
  }

  // Clear array of deposited energies 
  fEdep.Reset();

  AliRunLoader* runLoader = 0;
  if (!gAlice) { 
    TString folderName(fDigInput->GetInputFolderName(0));
    runLoader = AliRunLoader::GetRunLoader(folderName.Data());
    if (!runLoader) { 
      AliError(Form("Failed at getting run loader from %s",
		    folderName.Data()));
      return;
    }
    if (!runLoader->GetAliRun()) runLoader->LoadgAlice();
    runLoader->GetAliRun();
  }
  if (!gAlice) { 
    AliError("Can not get Run from Run Loader");
    return;
  }
  
  // Get the AliFMD object 
  fFMD = static_cast<AliFMD*>(gAlice->GetDetector("FMD"));
  if (!fFMD) {
    AliError("Can not get FMD from gAlice");
    return;
  }  


  // Loop over input files
  Int_t nFiles= fDigInput->GetNinputs();
  AliFMDDebug(1, (" Digitizing event number %d, got %d inputs",
		  fDigInput->GetOutputEventNr(), nFiles));
  for (Int_t inputFile = 0; inputFile < nFiles; inputFile++) {
    AliFMDDebug(5, ("Now reading input # %d", inputFile));
    // Get the current loader 
    AliRunLoader* currentLoader = 
      AliRunLoader::GetRunLoader(fDigInput->GetInputFolderName(inputFile));
    if (!currentLoader) { 
      Error("Exec", "no run loader for input file # %d", inputFile);
      continue;
    }

    // Cache contriutions 
    AliFMDDebug(5, ("Now summing the contributions from input # %d",inputFile));

    // Get the FMD loader 
    AliLoader* inFMD = currentLoader->GetLoader("FMDLoader");
    // And load the summable digits
    inFMD->LoadSDigits("READ");
  
    // Get the tree of summable digits
    TTree* sdigitsTree = inFMD->TreeS();
    if (!sdigitsTree)  {
      AliError("No sdigit tree from input");
      continue;
    }
    if (AliLog::GetDebugLevel("FMD","") >= 10) {
      TFile* file = sdigitsTree->GetCurrentFile();
      if (!file) {
	AliWarning("Input tree has no file!");
      }
      else { 
	AliFMDDebug(10, ("Input tree file %s content:", file->GetName()));
	file->ls();
      }
      // AliFMDDebug(5, ("Input tree %s file structure:", 
      //                 sdigitsTree->GetName()));
      // sdigitsTree->Print();
    }

    // Get the FMD branch 
    TBranch* sdigitsBranch = sdigitsTree->GetBranch("FMD");
    if (!sdigitsBranch) {
      AliError("Failed to get sdigit branch");
      return;
    }

    // Set the branch addresses 
    fFMD->SetTreeAddress();

    // Sum contributions from the sdigits
    AliFMDDebug(3, ("Will now sum contributions from SDigits"));
    SumContributions(sdigitsBranch);

    // Unload the sdigits
    inFMD->UnloadSDigits();
  }  

  TString       outFolder(fDigInput->GetOutputFolderName());
  AliRunLoader* out    = AliRunLoader::GetRunLoader(outFolder.Data());
  AliLoader*    outFMD = out->GetLoader("FMDLoader");
  if (!outFMD) { 
    AliError("Cannot get the FMDLoader output folder");
    return;
  }
  TTree* outTree = MakeOutputTree(outFMD);
  if (!outTree) { 
    AliError("Failed to get output tree");
    return;
  }
  // Set the branch address 
  fFMD->SetTreeAddress();
  
  // And digitize the cached data 
  DigitizeHits();
  
  // Fill the tree
  Int_t write = outTree->Fill();
  AliFMDDebug(5, ("Wrote %d bytes to digit tree", write));
  
  // Store the digits
  StoreDigits(outFMD);
}

//____________________________________________________________________
void
AliFMDDigitizer::SumContributions(TBranch* sdigitsBranch) 
{
  // 
  // Sum contributions from SDigits 
  // 
  // Parameters:
  //    sdigitsBranch Branch of SDigit data 
  //
  AliFMDDebug(3, ("Runnin our version of SumContributions"));

  // Get a list of hits from the FMD manager 
  TClonesArray *fmdSDigits = fFMD->SDigits();
  
  // Get number of entries in the tree 
  Int_t nevents  = Int_t(sdigitsBranch->GetEntries());
  
  Int_t read = 0;
  // Loop over the events in the 
  for (Int_t event = 0; event < nevents; event++)  {
    // Read in entry number `event' 
    read += sdigitsBranch->GetEntry(event);
    
    // Get the number of sdigits 
    Int_t nsdigits = fmdSDigits->GetEntries ();
    AliFMDDebug(3, ("Got %5d SDigits", nsdigits));
    for (Int_t sdigit = 0; sdigit < nsdigits; sdigit++) {
      // Get the sdigit number `sdigit'
      AliFMDSDigit* fmdSDigit = 
	static_cast<AliFMDSDigit*>(fmdSDigits->UncheckedAt(sdigit));

      AliFMDDebug(5, ("Adding contribution of %d tracks", 
		      fmdSDigit->GetNTrack()));
      AliFMDDebug(15, ("Contrib from FMD%d%c[%2d,%3d] (%s) from track %d", 
		       fmdSDigit->Detector(),
		       fmdSDigit->Ring(),
		       fmdSDigit->Sector(),
		       fmdSDigit->Strip(),
		       fmdSDigit->GetName(), 
		       fmdSDigit->GetTrack(0)));
      
      // Extract parameters 
      AddContribution(fmdSDigit->Detector(),
		      fmdSDigit->Ring(),
		      fmdSDigit->Sector(),
		      fmdSDigit->Strip(),
		      fmdSDigit->Edep(), 
		      kTRUE,
		      fmdSDigit->GetNTrack(),
		      fmdSDigit->GetTracks());
    }  // sdigit loop
  } // event loop


  AliFMDDebug(3, ("Size of cache: %d bytes, read %d bytes", 
		  int(sizeof(fEdep)), read));
}

//____________________________________________________________________
//
// EOF
// 




