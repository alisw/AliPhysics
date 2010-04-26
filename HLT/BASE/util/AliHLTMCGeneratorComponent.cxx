// -*- Mode: C++ -*-
// $Id: AliHLTMCGeneratorComponent.cxx  $
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Jochen Thaeder <thaeder@kip.uni-heidelberg.de>        *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTMCGeneratorComponent.cxx
    @author Jochen Thaeder
    @date   
    @brief  Component for generating MC events
    @note   The class is used in Offline (AliRoot) context
*/

#include "AliHLTMCGeneratorComponent.h"

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TH1F.h>
#include <TStopwatch.h>
#include <TDatime.h>
#include <TRandom.h>
#include <TDatabasePDG.h>
#include <TParticle.h>
#include <TArrayI.h>

#include "AliGenerator.h"
#include "AliPDG.h"
#include "AliRunLoader.h"
#include "AliRun.h"
#include "AliStack.h"
#include "AliHeader.h"
#include "AliGenPythia.h"
#include "AliPythia.h"
#endif

#include "TPythia6Calls.h"
#include "TPythia6.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTMCGeneratorComponent)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
AliHLTMCGeneratorComponent::AliHLTMCGeneratorComponent() :  
  AliHLTDataSource(),
  fpMC(NULL),
  fpHLTMC(NULL),
  fNEvents(0),
  fCurrentEvent(0),
  fEventNumber(0),
  fRunNumber(0),
  fRunLoader(NULL),
  fGenerator(NULL),
  fComment(""),
  fSeed(0),
  fRunType(kPythia6Jets104_125),
  fEcms(14000.),
  fJetEtaMax(0.2), 
  fJetEtaMin(-0.2),
  fJetEtMax(1000.),
  fJetEtMin(10.0),
  fJetConeRadius(0.4),
  fPtHardMin(10.0),
  fPtHardMax(-1.0),
  fQuenching(0),                   
  fQhat(20.),
  fApplyParticleCuts(kFALSE) {
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

// #################################################################################
AliHLTMCGeneratorComponent::~AliHLTMCGeneratorComponent() {
  // see header file for class documentation

  // file list and file name list are owner of their objects and
  // delete all the objects
}

/*
 * ---------------------------------------------------------------------------------
 *                              Initialize static const
 * ---------------------------------------------------------------------------------
 */

//##################################################################################
const Char_t *AliHLTMCGeneratorComponent::fgkPprRunName[] = {
  "kPythia6Jets20_24",   "kPythia6Jets24_29",   "kPythia6Jets29_35",
  "kPythia6Jets35_42",   "kPythia6Jets42_50",   "kPythia6Jets50_60",
  "kPythia6Jets60_72",   "kPythia6Jets72_86",   "kPythia6Jets86_104",
  "kPythia6Jets104_125", "kPythia6Jets125_150", "kPythia6Jets150_180",
  "kPyJetJet", "kPyGammaJetPHOS"
};

/*
 * ---------------------------------------------------------------------------------
 * Public functions to implement AliHLTComponent's interface.
 * These functions are required for the registration process
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
const Char_t* AliHLTMCGeneratorComponent::GetComponentID() {
  // see header file for class documentation
  return "MCGenerator";
}

// #################################################################################
AliHLTComponentDataType AliHLTMCGeneratorComponent::GetOutputDataType() {
  // see header file for class documentation

  return kAliHLTDataTypeMCObject|kAliHLTDataOriginHLT;
}

// #################################################################################
void AliHLTMCGeneratorComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ) {
  // see header file for class documentation

  constBase = 20000;
  inputMultiplier = 1.0;
}

// #################################################################################
AliHLTComponent* AliHLTMCGeneratorComponent::Spawn() {
  // see header file for class documentation
  return new AliHLTMCGeneratorComponent;
}

/*
 * ---------------------------------------------------------------------------------
 * Protected functions to implement AliHLTComponent's interface.
 * These functions provide initialization as well as the actual processing
 * capabilities of the component. 
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
Int_t AliHLTMCGeneratorComponent::DoInit(int argc, const char** argv) {
  // see header file for class documentation

  Int_t iResult = 0;
  Int_t bMissingParam=0;
  
  TString argument="";

  // -- Loop over all arguments
  for ( Int_t iter = 0; iter<argc && iResult>=0; iter++) {
    argument=argv[iter];
    if (argument.IsNull()) continue;
    
    // -- seed
    if ( !argument.CompareTo("-seed") ) {
      if ((bMissingParam=(++iter>=argc))) break;
      
      TString parameter(argv[iter]);
      parameter.Remove(TString::kLeading, ' ');

      if ( parameter.IsDigit() ) {
	fSeed = parameter.Atoi();
      }
      
      else {
	HLTError("Wrong parameter %s for argument %s.", parameter.Data(), argument.Data());
	iResult=-EINVAL;
      }
    } 

    // -- runtype
    else if ( !argument.CompareTo("-runtype") ) {
      if ((bMissingParam=(++iter>=argc))) break;
      
      TString parameter(argv[iter]);
      parameter.Remove(TString::kLeading, ' ');

      Bool_t found = kFALSE;

      if ( parameter.IsDigit() ) {
	Int_t idx = parameter.Atoi();
	if ( idx < kRunMax ) {
	  fRunType = (PprRun_t) idx;
	  found = kTRUE;
	}
      }
      else {
	for ( Int_t iRun = 0; iRun < kRunMax && !found; iRun++ ) {
	  if ( ! parameter.CompareTo(fgkPprRunName[iRun]) ) {
	    fRunType = (PprRun_t) iRun;  
	    found = kTRUE;
	  }
	}      
      }
      if ( !found ) {
	HLTError("Wrong parameter %s for argument %s.", parameter.Data(), argument.Data());
	iResult=-EINVAL;
      }
    }

    // -- nevents
    else if ( !argument.CompareTo("-nevents") ) {
      if ((bMissingParam=(++iter>=argc))) break;
      
      TString parameter(argv[iter]);
      parameter.Remove(TString::kLeading, ' ');

      if ( parameter.IsDigit() ) {
	fNEvents = parameter.Atoi();
      }
      
      else {
	HLTError("Wrong parameter %s for argument %s.", parameter.Data(), argument.Data());
	iResult=-EINVAL;
      }
    } 

    // -- coneRadius
    else if ( !argument.CompareTo("-coneRadius") ) {
      if ((bMissingParam=(++iter>=argc))) break;
      
      TString parameter(argv[iter]);
      parameter.Remove(TString::kLeading, ' ');
      
      if ( parameter.IsFloat() ) {
	fJetConeRadius = parameter.Atof();
      }
      else {
	HLTError("Wrong parameter %s for argument %s.", parameter.Data(), argument.Data());
	iResult=-EINVAL;
      }
    } 

    // -- jetCutMinEt
    else if ( !argument.CompareTo("-jetCutMinEt") ) {
      if ((bMissingParam=(++iter>=argc))) break;
      
      TString parameter(argv[iter]);
      parameter.Remove(TString::kLeading, ' ');
      
      if ( parameter.IsFloat() ) {
	fJetEtMin = parameter.Atof();
      }
      else {
	HLTError("Wrong parameter %s for argument %s.", parameter.Data(), argument.Data());
	iResult=-EINVAL;
      }
    } 

    // -- applyParticleCuts
    else if ( !argument.CompareTo("-applyParticleCuts") ) {
      fApplyParticleCuts = kTRUE;
    } 

    // -- Argument not known
    else {
      HLTError("Unknown argument %s.", argument.Data());
      iResult = -EINVAL;
    }
    
  } // for ( Int iter = 0; iter<argc && iResult>=0; iter++) {
  
  // -- Check if parameter is missing
  if ( bMissingParam ) {
    HLTError("Missing parameter for argument %s.", argument.Data());
     iResult=-EINVAL;
  }

  //  -- Check if seed and nevents has been specified
  if ( !fSeed || !fNEvents ) {
    HLTError("-seed and -nevents need to be specified.");
    iResult=-EINVAL;
  }

  // -- Set Random Number seed
  gRandom->SetSeed( fSeed );

  // -- Run loader 
  // ---------------
  AliPDG::AddParticlesToPdgDataBase();

  fRunLoader = AliRunLoader::Open("galice.root","FASTRUN","recreate");
  if (fRunLoader == NULL ) {
    //    gAlice->Fatal("AliHLTMCGeneratorComponent","Can not instatiate the Run Loader");
    HLTFatal("Can not instatiate the Run Loader");
    return -EINPROGRESS;
  }

  fRunLoader->SetCompressionLevel(2);
  fRunLoader->SetNumberOfEventsPerFile(fNEvents);
  fRunLoader->LoadKinematics("RECREATE");
  fRunLoader->MakeTree("E");
  gAlice->SetRunLoader(fRunLoader); 

  // -- Create stack
  fRunLoader->MakeStack();
 
  // -- Create and Initialize Generator
  fGenerator = GeneratorFactory();
  fGenerator->SetStack(fRunLoader->Stack());
  fGenerator->Init();

  HLTInfo("MC Generator setup with : %s", fComment.Data() );

  return iResult;
}

// #################################################################################
Int_t AliHLTMCGeneratorComponent::DoDeinit() {
  // see header file for class documentation

  if ( fpMC ) 
    delete fpMC;
  fpMC = NULL;

  if ( fpHLTMC ) 
    delete fpHLTMC;
  fpHLTMC = NULL;

  return 0;
}

// #################################################################################
Int_t AliHLTMCGeneratorComponent::GetEvent( const AliHLTComponentEventData& /*evtData*/,
						AliHLTComponentTriggerData& /*trigData*/,
						AliHLTUInt8_t* /*outputPtr*/, 
						AliHLTUInt32_t& size,
						vector<AliHLTComponentBlockData>& /*outputBlocks*/ ) {
  // see header file for class documentation

  if ( !IsDataEvent() ) 
    return 0;

  Int_t iResult=0;
  size=0;

  AliStack*  stack  = fRunLoader->Stack();
  AliHeader* header = fRunLoader->GetHeader();
  
  // -- Initialize event
  header->Reset(0,fCurrentEvent);

  // -- Reset stack 
  stack->Reset();

  // -- Generate event
  fGenerator->Generate();
    
  // -- Finish event
  header->SetNprimary(stack->GetNprimary());
  header->SetNtrack(stack->GetNtrack());  

  // -- Finish Event	
  stack->FinishEvent();
  header->SetStack(stack);

  // -- Create HLT MC Event
  if ( fpHLTMC ) 
    delete fpHLTMC;
  fpHLTMC = new AliHLTMCEvent( fApplyParticleCuts);

  // -- Fill HLT MC event
  iResult = fpHLTMC->FillMCEvent( stack, header );

  // -- Send HLT MC event
  if ( fpHLTMC && !iResult )
    PushBack( fpHLTMC, kAliHLTDataTypeMCObject|kAliHLTDataOriginHLT, 0 ); 

  fCurrentEvent++;

  return iResult;
}

// #################################################################################
AliGenerator* AliHLTMCGeneratorComponent::GeneratorFactory() {
  
  // -- Configuration for hard QCD processes generation with PYTHIA 
  // ----------------------------------------------------------------
  AliGenPythia * gener = new AliGenPythia(-1);

  gener->SetEnergyCMS(fEcms); // Centre of mass energy
  gener->SetGluonRadiation(1,1);
  gener->SetForceDecay(kAll); //  Decay type

  // -- Final state kinematic cuts
  gener->SetJetEtaRange(fJetEtaMin, fJetEtaMax);
  gener->SetJetPhiRange(0., 360.);
  gener->SetJetEtRange(fJetEtMin, fJetEtMax);
  // gener->SetMomentumRange(0,99999999);

  // -- Structure function
  gener->SetStrucFunc(kCTEQ4L);

  /************************************************************************
   * void AliGenPythia::SetPycellParameters()
   *  Float_t etamax   -- fPycellEtaMax     -- PARU(51) -- +/- eta max of detector
   *  Int_t neta       -- fPycellNEta       -- MSTU(51) -- N cells in eta
   *  Int_t nphi       -- fPycellNPhi       -- MSTU(52) -- N cells in phi
   *  Float_t thresh   -- fPycellThreshold  -- PARU(58) -- cells with Et below treshold are discared
   *  Float_t etseed   -- fPycellEtSeed     -- PARU(52) -- Et threshold for seed cells
   *  Float_t minet    -- fPycellMinEtJet   -- PARU(53) -- Min Et for jets
   *  Float_t r        -- fPycellMaxRadius  -- PARU(54) -- cone radius
   ************************************************************************/ 
  gener->SetPycellParameters(fJetEtaMax, 274, 432, 0., 4., fJetEtMin, fJetConeRadius);

  // gener->SetPtKick(5.);      // set the intrinsic kt to 5 GeV/c
  
  // -- Jet Quenching
  gener->SetQuench(fQuenching);
  if( fQuenching == 1){
    Float_t k = 6e5*(fQhat/1.7) ; // qhat=1.7, k = 6e5, default  value 
    AliPythia::Instance()->InitQuenching(0.,0.1,k,0,0.95,6);
  }

  switch (fRunType) {
  case kPythia6Jets20_24:
    {
      fComment = fComment.Append(":Pythia jets 20-24 GeV @ 14 TeV");
      gener->SetProcess(kPyJets);//        Process type
      gener->SetPtHard(20., 24.);// Pt transfer of the hard scattering
      fPtHardMin=20.;
    }
    break;
  case kPythia6Jets24_29:
    {
      fComment = fComment.Append(":Pythia jets 24-29 GeV @ 14 TeV");
      gener->SetProcess(kPyJets);//        Process type
      gener->SetPtHard(24., 29.);// Pt transfer of the hard scattering
      fPtHardMin=24.;
    }
    break;
  case kPythia6Jets29_35:
    {
      fComment = fComment.Append(":Pythia jets 29-35 GeV @ 14 TeV");
      gener->SetProcess(kPyJets);//        Process type
      gener->SetPtHard(29., 35.);// Pt transfer of the hard scattering
      fPtHardMin=29.;
    }
    break;
  case kPythia6Jets35_42:
    {
      fComment = fComment.Append(":Pythia jets 35-42 GeV @ 14 TeV");
      gener->SetProcess(kPyJets);//        Process type
      gener->SetPtHard(35., 42.);// Pt transfer of the hard scattering
      fPtHardMin=35.;
    }
    break;
  case kPythia6Jets42_50:
    {
      fComment = fComment.Append(":Pythia jets 42-50 GeV @ 14 TeV");
      gener->SetProcess(kPyJets);//        Process type
      gener->SetPtHard(42., 50.);// Pt transfer of the hard scattering
      fPtHardMin=42.;
    }
    break;
  case kPythia6Jets50_60:
    {
      fComment = fComment.Append(":Pythia jets 50-60 GeV @ 14 TeV");
      gener->SetProcess(kPyJets);//        Process type
      gener->SetPtHard(50., 60.);// Pt transfer of the hard scattering
      fPtHardMin=50.;
    }
    break;
  case kPythia6Jets60_72:
    {
      fComment = fComment.Append(":Pythia jets 60-72 GeV @ 14 TeV");
      gener->SetProcess(kPyJets);//        Process type
      gener->SetPtHard(60., 72.);// Pt transfer of the hard scattering
      fPtHardMin=60.;
    }
    break;
  case kPythia6Jets72_86:
    {
      fComment = fComment.Append(":Pythia jets 72-86 GeV @ 14 TeV");
      gener->SetProcess(kPyJets);//        Process type
      gener->SetPtHard(72., 86.);// Pt transfer of the hard scattering
      fPtHardMin=72.;
    }
    break;
  case kPythia6Jets86_104:
    {
      fComment = fComment.Append(":Pythia jets 86-104 GeV @ 14 TeV");
      gener->SetProcess(kPyJets);//        Process type
      gener->SetPtHard(86., 104.);// Pt transfer of the hard scattering
      fPtHardMin=86.;
    }
    break;
  case kPythia6Jets104_125:
    {
      fComment = fComment.Append(":Pythia jets 104-125 GeV @ 14 TeV");
      gener->SetProcess(kPyJets);//        Process type
      gener->SetPtHard(104., 125.);// Pt transfer of the hard scattering
      fPtHardMin=104.;
    }
    break;
  case kPythia6Jets125_150:
    {
      fComment = fComment.Append(":Pythia jets 125-150 GeV @ 14 TeV");
      gener->SetProcess(kPyJets);//        Process type
      gener->SetPtHard(125., 150.);// Pt transfer of the hard scattering
      fPtHardMin=125.;
    }
    break;
  case kPythia6Jets150_180:
    {
      fComment = fComment.Append(":Pythia jets 150-180 GeV @ 14 TeV");
      gener->SetProcess(kPyJets);//        Process type
      gener->SetPtHard(150., 180.);// Pt transfer of the hard scattering
      fPtHardMin=150.;
    }
    break;
  case kPyJetJet:
    {
      fComment = fComment.Append(" Jet-jet @ 14 TeV");
      gener->SetProcess(kPyJets); 
      gener->SetPtHard(fPtHardMin,fPtHardMax);
      gener->SetEventListRange(0,1);  // XXXX
    }
    break;
  case kPyGammaJetPHOS:
    {
      fComment = fComment.Append(" pp->jet + Gamma over PHOS @ 14 TeV");
      gener->SetProcess(kPyDirectGamma);
      gener->SetPtHard(fPtHardMin,fPtHardMax);
      //gener->SetYHard(-1.0,1.0);
      gener->SetGammaEtaRange(-0.13,0.13);
      gener->SetGammaPhiRange(218.,322.); //Over 5 modules +-2 deg  
      gener->SetEventListRange(0,1); // XXXX
    }
    break;
    
  default: break;
  }
  
  return gener;
}
