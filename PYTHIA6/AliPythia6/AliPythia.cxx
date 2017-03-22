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

// Pythia 6 interface used by AliGenPythia
// Some settings are done by AliGenPythia, others here :)
//
/* $Id$ */

#include "AliPythia.h"
#include "AliPythiaRndm.h"
#include "AliFastGlauber.h"
#include "AliQuenchingWeights.h"
#include "AliOmegaDalitz.h"
#include "AliDecayerExodus.h"
#include "AliLog.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "PyquenCommon.h"

ClassImp(AliPythia)

#ifndef WIN32
# define pyclus pyclus_
# define pycell pycell_
# define pyshow pyshow_
# define pyrobo pyrobo_
# define pyquen pyquen_
# define pyevnw pyevnw_
# define pyshowq pyshowq_
# define qpygin0 qpygin0_
# define pytune  pytune_
# define py2ent  py2ent_
# define setpowwght setpowwght_
# define type_of_call
#else
# define pyclus PYCLUS
# define pycell PYCELL
# define pyrobo PYROBO
# define pyquen PYQUEN
# define pyevnw PYEVNW
# define pyshowq PYSHOWQ
# define qpygin0 QPYGIN0
# define pytune  PYTUNE
# define py2ent  PY2ENT
# define setpowwght SETPOWWGHT
# define type_of_call _stdcall
#endif

extern "C" void type_of_call pyclus(Int_t & );
extern "C" void type_of_call pycell(Int_t & );
extern "C" void type_of_call pyshow(Int_t &, Int_t &, Double_t &);
extern "C" void type_of_call pyrobo(Int_t &, Int_t &, Double_t &, Double_t &, Double_t &, Double_t &, Double_t &);
extern "C" void type_of_call pyquen(Double_t &, Int_t &, Double_t &);
extern "C" void type_of_call pyevnw();
extern "C" void type_of_call pyshowq(Int_t &, Int_t &, Double_t &);
extern "C" void type_of_call pytune(Int_t &);
extern "C" void type_of_call py2ent(Int_t &, Int_t&, Int_t&, Double_t&);
extern "C" void type_of_call qpygin0();
extern "C" void type_of_call setpowwght(Double_t &);
//_____________________________________________________________________________

AliPythia* AliPythia::fgAliPythia=NULL;

AliPythia::AliPythia():
    fProcess(kPyMb),
    fEcms(0.),
    fStrucFunc(kCTEQ5L),
    fProjectile("p"),
    fTarget("p"),
    fXJet(0.),
    fYJet(0.),
    fNGmax(30),
    fZmax(0.97),
    fGlauber(0),
    fQuenchingWeights(0),
    fItune(-1),
    fOmegaDalitz(),
    fExodus()
{
// Default Constructor
//
//  Set random number
    if (!AliPythiaRndm::GetPythiaRandom()) 
      AliPythiaRndm::SetPythiaRandom(GetRandom());
    fGlauber          = 0;
    fQuenchingWeights = 0;
    Int_t i = 0;
    for (i = 0; i <  501; i++) fDefMDCY[i] = 0;
    for (i = 0; i < 2001; i++) fDefMDME[i] = 0;
    for (i = 0; i <    4; i++) fZQuench[i] = 0;
}

AliPythia::AliPythia(const AliPythia& pythia):
    TPythia6(pythia), 
    AliRndm(pythia),
    fProcess(kPyMb),
    fEcms(0.),
    fStrucFunc(kCTEQ5L),
    fProjectile("p"),
    fTarget("p"),
    fXJet(0.),
    fYJet(0.),
    fNGmax(30),
    fZmax(0.97),
    fGlauber(0),
    fQuenchingWeights(0),
    fItune(-1),
    fOmegaDalitz(),
    fExodus()
{
    // Copy Constructor
    Int_t i;
    for (i = 0; i <  501; i++) fDefMDCY[i] = 0;
    for (i = 0; i < 2001; i++) fDefMDME[i] = 0;
    for (i = 0; i <    4; i++) fZQuench[i] = 0;
    pythia.Copy(*this);
}

void AliPythia::ProcInit(Process_t process, Float_t energy, StrucFunc_t strucfunc, Int_t itune)
{
// Initialise the process to generate 
    if (!AliPythiaRndm::GetPythiaRandom()) 
      AliPythiaRndm::SetPythiaRandom(GetRandom());
    
    fItune = itune;
    
    fProcess = process;
    fEcms = energy;
    fStrucFunc = strucfunc;
//...Switch off decay of pi0, K0S, Lambda, Sigma+-, Xi0-, Omega-.
    SetMDCY(Pycomp(111) ,1,0); // pi0
    SetMDCY(Pycomp(310) ,1,0); // K0S
    SetMDCY(Pycomp(3122),1,0); // kLambda
    SetMDCY(Pycomp(3112),1,0); // sigma -
    SetMDCY(Pycomp(3222),1,0); // sigma +
    SetMDCY(Pycomp(3312),1,0); // xi - 
    SetMDCY(Pycomp(3322),1,0); // xi 0
    SetMDCY(Pycomp(3334),1,0); // omega-
    // Select structure function 
    SetMSTP(52,2);
    SetMSTP(51, AliStructFuncType::PDFsetIndex(strucfunc));
    // Particles produced in string fragmentation point directly to either of the two endpoints
    // of the string (depending in the side they were generated from).
    SetMSTU(16,2);

//
// Pythia initialisation for selected processes//
//
// Make MSEL clean
//
    for (Int_t i=1; i<= 200; i++) {
	SetMSUB(i,0);
    }
//  select charm production
    switch (process) 
    {
    case kPyOldUEQ2ordered:  //Old underlying events with Q2 ordered QCD processes
//        Multiple interactions on.
	SetMSTP(81,1);
// Double Gaussian matter distribution.
	SetMSTP(82,4);
	SetPARP(83,0.5);
	SetPARP(84,0.4);
//  pT0.
	SetPARP(82,2.0);
//  Reference energy for pT0 and energy rescaling pace.
	SetPARP(89,1800);
	SetPARP(90,0.25);
//  String drawing almost completely minimizes string length.
	SetPARP(85,0.9);
	SetPARP(86,0.95);
// ISR and FSR activity.
	SetPARP(67,4);
	SetPARP(71,4);
// Lambda_FSR scale.
	SetPARJ(81,0.29);
	break;
    case kPyOldUEQ2ordered2:   
// Old underlying events with Q2 ordered QCD processes
// Multiple interactions on.
	SetMSTP(81,1);
// Double Gaussian matter distribution.
	SetMSTP(82,4);
	SetPARP(83,0.5);
	SetPARP(84,0.4);
// pT0.
	SetPARP(82,2.0);
// Reference energy for pT0 and energy rescaling pace.
	SetPARP(89,1800);
	SetPARP(90,0.16);  // here is the difference with  kPyOldUEQ2ordered
// String drawing almost completely minimizes string length.
	SetPARP(85,0.9);
	SetPARP(86,0.95);
// ISR and FSR activity.
	SetPARP(67,4);
	SetPARP(71,4);
// Lambda_FSR scale.
	SetPARJ(81,0.29);	
	break;
    case kPyOldPopcorn:  
// Old production mechanism: Old Popcorn
	SetMSEL(1);
	SetMSTJ(12,3); 
// (D=2) Like MSTJ(12)=2 but added prod ofthe 1er rank baryon
	SetMSTP(88,2); 
// (D=1)see can be used to form  baryons (BARYON JUNCTION)
	SetMSTJ(1,1);  
	AtlasTuning();
	break;
    case kPyCharm:
	SetMSEL(4);
//  heavy quark masses

	SetPMAS(4,1,1.2);
//
//    primordial pT
	SetMSTP(91,1);
	SetPARP(91,1.);
	SetPARP(93,5.);
//
	break;
    case kPyBeauty:
	SetMSEL(5);
	SetPMAS(5,1,4.75);
	break;
    case kPyJpsi:
	SetMSEL(0);
// gg->J/Psi g
	SetMSUB(86,1);
	break;
    case kPyJpsiChi:
	SetMSEL(0);
// gg->J/Psi g
	SetMSUB(86,1);
// gg-> chi_0c g
	SetMSUB(87,1);
// gg-> chi_1c g
	SetMSUB(88,1);
// gg-> chi_2c g
	SetMSUB(89,1);	
	break;
    case kPyCharmUnforced:
	SetMSEL(0);
// gq->qg   
	SetMSUB(28,1);
// gg->qq
	SetMSUB(53,1);
// gg->gg
	SetMSUB(68,1);
	break;
    case kPyBeautyUnforced:
	SetMSEL(0);
// gq->qg   
	SetMSUB(28,1);
// gg->qq
	SetMSUB(53,1);
// gg->gg
	SetMSUB(68,1);
	break;
    case kPyMb:
// Minimum Bias pp-Collisions
//
//   
//      select Pythia min. bias model
	SetMSEL(0);
	SetMSUB(92,1);             // single diffraction AB-->XB
	SetMSUB(93,1);             // single diffraction AB-->AX
	SetMSUB(94,1);             // double diffraction
	SetMSUB(95,1);	           // low pt production

	AtlasTuning();
	break;
	
    case kPyMbAtlasTuneMC09:
// Minimum Bias pp-Collisions
//
//   
//      select Pythia min. bias model
	SetMSEL(0);
	SetMSUB(92,1);             // single diffraction AB-->XB
	SetMSUB(93,1);             // single diffraction AB-->AX
	SetMSUB(94,1);             // double diffraction
	SetMSUB(95,1);	           // low pt production

	AtlasTuningMC09();
	break;

    case kPyMbWithDirectPhoton:
// Minimum Bias pp-Collisions with direct photon processes added 
//
//   
//      select Pythia min. bias model
	SetMSEL(0);
	SetMSUB(92,1);             // single diffraction AB-->XB
	SetMSUB(93,1);             // single diffraction AB-->AX
	SetMSUB(94,1);             // double diffraction
	SetMSUB(95,1);	           // low pt production

	SetMSUB(14,1);             //
	SetMSUB(18,1);             //
	SetMSUB(29,1);             //
	SetMSUB(114,1);            //
	SetMSUB(115,1);            //


	AtlasTuning();
	break;

    case kPyMbDefault:
// Minimum Bias pp-Collisions
//
//   
//      select Pythia min. bias model
	SetMSEL(0);
	SetMSUB(92,1);             // single diffraction AB-->XB
	SetMSUB(93,1);             // single diffraction AB-->AX
	SetMSUB(94,1);             // double diffraction
	SetMSUB(95,1);	           // low pt production
	break;
    case kPyLhwgMb:
// Les Houches Working Group 05 Minimum Bias pp-Collisions: hep-ph/0604120
//  -> Pythia 6.3 or above is needed
//   
	SetMSEL(0);
	SetMSUB(92,1);             // single diffraction AB-->XB
	SetMSUB(93,1);             // single diffraction AB-->AX
	SetMSUB(94,1);             // double diffraction
	SetMSUB(95,1);	           // low pt production

        SetMSTP(51,AliStructFuncType::PDFsetIndex(kCTEQ6ll));      // CTEQ6ll pdf
	SetMSTP(52,2);
	SetMSTP(68,1);
	SetMSTP(70,2);
	SetMSTP(81,1);             // Multiple Interactions ON
	SetMSTP(82,4);             // Double Gaussian Model
	SetMSTP(88,1);

	SetPARP(82,2.3);           // [GeV]    PT_min at Ref. energy
	SetPARP(83,0.5);           // Core density in proton matter distribution (def.value)
	SetPARP(84,0.5);           // Core radius
	SetPARP(85,0.9);           // Regulates gluon prod. mechanism
	SetPARP(90,0.2);           // 2*epsilon (exponent in power law)

	break;
    case kPyMbNonDiffr:
// Minimum Bias pp-Collisions
//
//   
//      select Pythia min. bias model
	SetMSEL(0);
	SetMSUB(95,1);	           // low pt production

	AtlasTuning();
	break;
    case kPyMbMSEL1:
	ConfigHeavyFlavor();
// Intrinsic <kT^2>
        SetMSTP(91,1);// Width (1=gaussian) primordial kT dist. inside hadrons
        SetPARP(91,1.);     // <kT^2> = PARP(91,1.)^2
        SetPARP(93,5.);     // Upper cut-off
// Set Q-quark mass
        SetPMAS(4,1,1.2);   // Charm quark mass
        SetPMAS(5,1,4.78);  // Beauty quark mass
	SetPARP(71,4.);     // Defaut value
// Atlas Tuning
	AtlasTuning();
	break;
    case kPyJets:
//
//  QCD Jets
//
	SetMSEL(1);

 // Pythia Tune A (CDF)
 //
	if (fItune < 0) {
	  SetPARP(67,2.5);           // Regulates Initial State Radiation (value from best fit to D0 dijet analysis)
	  SetMSTP(82,4);             // Double Gaussian Model
	  SetPARP(82,2.0);           // [GeV]    PT_min at Ref. energy
	  SetPARP(84,0.4);           // Core radius
	  SetPARP(85,0.90) ;         // Regulates gluon prod. mechanism
	  SetPARP(86,0.95);          // Regulates gluon prod. mechanism
	  SetPARP(89,1800.);         // [GeV]   Ref. energy
	  SetPARP(90,0.25);	     // 2*epsilon (exponent in power law)
	}
	  break;
    case kPyDirectGamma:
	SetMSEL(10);
	break;
    case kPyCharmPbPbMNR:
    case kPyD0PbPbMNR:
    case kPyDPlusPbPbMNR:
    case kPyDPlusStrangePbPbMNR:
      // Tuning of Pythia parameters aimed to get a resonable agreement
      // between with the NLO calculation by Mangano, Nason, Ridolfi for the
      // c-cbar single inclusive and double differential distributions.
      // This parameter settings are meant to work with Pb-Pb collisions
      // (AliGenPythia::SetNuclei) and with kCTEQ4L PDFs.
      // To get a good agreement the minimum ptHard (AliGenPythia::SetPtHard)
      // has to be set to 2.1GeV. Example in ConfigCharmPPR.C.
	ConfigHeavyFlavor();
      // Intrinsic <kT>
      SetMSTP(91,1);
      SetPARP(91,1.304);
      SetPARP(93,6.52);
      // Set c-quark mass
      SetPMAS(4,1,1.2);
      break;
    case kPyCharmpPbMNR:
    case kPyD0pPbMNR:
    case kPyDPluspPbMNR:
    case kPyDPlusStrangepPbMNR:
      // Tuning of Pythia parameters aimed to get a resonable agreement
      // between with the NLO calculation by Mangano, Nason, Ridolfi for the
      // c-cbar single inclusive and double differential distributions.
      // This parameter settings are meant to work with p-Pb collisions
      // (AliGenPythia::SetNuclei) and with kCTEQ4L PDFs.
      // To get a good agreement the minimum ptHard (AliGenPythia::SetPtHard)
      // has to be set to 2.1GeV. Example in ConfigCharmPPR.C.
	ConfigHeavyFlavor();
      // Intrinsic <kT>
	SetMSTP(91,1);
	SetPARP(91,1.16);
	SetPARP(93,5.8);
	
      // Set c-quark mass
	SetPMAS(4,1,1.2);
      break;
    case kPyCharmppMNR:
    case kPyD0ppMNR:
    case kPyDPlusppMNR:
    case kPyDPlusStrangeppMNR:
    case kPyLambdacppMNR:
      // Tuning of Pythia parameters aimed to get a resonable agreement
      // between with the NLO calculation by Mangano, Nason, Ridolfi for the
      // c-cbar single inclusive and double differential distributions.
      // This parameter settings are meant to work with pp collisions
      // (AliGenPythia::SetNuclei) and with kCTEQ4L PDFs.
      // To get a good agreement the minimum ptHard (AliGenPythia::SetPtHard)
      // has to be set to 2.1GeV. Example in ConfigCharmPPR.C.
	ConfigHeavyFlavor();
      // Intrinsic <kT^2>
	SetMSTP(91,1);
	SetPARP(91,1.);
	SetPARP(93,5.);
	
      // Set c-quark mass
	SetPMAS(4,1,1.2);
      break;
    case kPyCharmppMNRwmi:
      // Tuning of Pythia parameters aimed to get a resonable agreement
      // between with the NLO calculation by Mangano, Nason, Ridolfi for the
      // c-cbar single inclusive and double differential distributions.
      // This parameter settings are meant to work with pp collisions
      // and with kCTEQ5L PDFs.
      // Added multiple interactions according to ATLAS tune settings.
      // To get a "reasonable" agreement with MNR results, events have to be 
      // generated with the minimum ptHard (AliGenPythia::SetPtHard)
      // set to 2.76 GeV.
      // To get a "perfect" agreement with MNR results, events have to be 
      // generated in four ptHard bins with the following relative 
      // normalizations:
      // 2.76-3 GeV: 25%
      //    3-4 GeV: 40%
      //    4-8 GeV: 29%
      //     >8 GeV:  6%
	ConfigHeavyFlavor();
      // Intrinsic <kT^2>
	SetMSTP(91,1);
	SetPARP(91,1.);
	SetPARP(93,5.);

      // Set c-quark mass
	SetPMAS(4,1,1.2);
	AtlasTuning();
	break;
    case kPyBeautyPbPbMNR:
      // Tuning of Pythia parameters aimed to get a resonable agreement
      // between with the NLO calculation by Mangano, Nason, Ridolfi for the
      // b-bbar single inclusive and double differential distributions.
      // This parameter settings are meant to work with Pb-Pb collisions
      // (AliGenPythia::SetNuclei) and with kCTEQ4L PDFs.
      // To get a good agreement the minimum ptHard (AliGenPythia::SetPtHard)
      // has to be set to 2.75GeV. Example in ConfigBeautyPPR.C.
	ConfigHeavyFlavor();
      // QCD scales
	SetPARP(67,1.0);
	SetPARP(71,1.0);
      // Intrinsic <kT>
	SetMSTP(91,1);
	SetPARP(91,2.035);
	SetPARP(93,10.17);
      // Set b-quark mass
	SetPMAS(5,1,4.75);
      break;
    case kPyBeautypPbMNR:
      // Tuning of Pythia parameters aimed to get a resonable agreement
      // between with the NLO calculation by Mangano, Nason, Ridolfi for the
      // b-bbar single inclusive and double differential distributions.
      // This parameter settings are meant to work with p-Pb collisions
      // (AliGenPythia::SetNuclei) and with kCTEQ4L PDFs.
      // To get a good agreement the minimum ptHard (AliGenPythia::SetPtHard)
      // has to be set to 2.75GeV. Example in ConfigBeautyPPR.C.
	ConfigHeavyFlavor();
      // QCD scales
	SetPARP(67,1.0);
	SetPARP(71,1.0);
      // Intrinsic <kT>
	SetMSTP(91,1);
	SetPARP(91,1.60);
	SetPARP(93,8.00);
      // Set b-quark mass
	SetPMAS(5,1,4.75);
      break;
    case kPyBeautyppMNR:
      // Tuning of Pythia parameters aimed to get a resonable agreement
      // between with the NLO calculation by Mangano, Nason, Ridolfi for the
      // b-bbar single inclusive and double differential distributions.
      // This parameter settings are meant to work with pp collisions
      // (AliGenPythia::SetNuclei) and with kCTEQ4L PDFs.
      // To get a good agreement the minimum ptHard (AliGenPythia::SetPtHard)
      // has to be set to 2.75GeV. Example in ConfigBeautyPPR.C.
	ConfigHeavyFlavor();
      // QCD scales
	SetPARP(67,1.0);
	SetPARP(71,1.0);
	
	// Intrinsic <kT>
	SetMSTP(91,1);
	SetPARP(91,1.);
	SetPARP(93,5.);
	
	// Set b-quark mass
	SetPMAS(5,1,4.75);
      break;
     case kPyBeautyJets:
     case kPyBeautyppMNRwmi:
      // Tuning of Pythia parameters aimed to get a resonable agreement
      // between with the NLO calculation by Mangano, Nason, Ridolfi for the
      // b-bbar single inclusive and double differential distributions.
      // This parameter settings are meant to work with pp collisions
      // and with kCTEQ5L PDFs.
      // Added multiple interactions according to ATLAS tune settings.
      // To get a "reasonable" agreement with MNR results, events have to be 
      // generated with the minimum ptHard (AliGenPythia::SetPtHard)
      // set to 2.76 GeV.
      // To get a "perfect" agreement with MNR results, events have to be 
      // generated in four ptHard bins with the following relative 
      // normalizations:
      // 2.76-4 GeV:  5% 
      //    4-6 GeV: 31%
      //    6-8 GeV: 28%
      //     >8 GeV: 36%
	 ConfigHeavyFlavor();
      // QCD scales
	 SetPARP(67,1.0);
	 SetPARP(71,1.0);
	 
	 // Intrinsic <kT>
	 SetMSTP(91,1);
	 SetPARP(91,1.);
	 SetPARP(93,5.);

      // Set b-quark mass
	 SetPMAS(5,1,4.75);

	 AtlasTuning();
	 break; 
    case kPyW:

      //Inclusive production of W+/-
      SetMSEL(0);
      //f fbar -> W+ 
      SetMSUB(2,1);
      // 	//f fbar -> g W+
      // 	SetMSUB(16,1);
      // 	//f fbar -> gamma W+
      // 	SetMSUB(20,1);
      // 	//f g -> f W+  
      // 	SetMSUB(31,1);
      // 	//f gamma -> f W+
      // 	SetMSUB(36,1);
      
      // Initial/final parton shower on (Pythia default)
      // With parton showers on we are generating "W inclusive process"
      SetMSTP(61,1); //Initial QCD & QED showers on
      SetMSTP(71,1); //Final QCD & QED showers on
      
      break;  

    case kPyZ:

      //Inclusive production of Z
      SetMSEL(0);
      //f fbar -> Z/gamma
      SetMSUB(1,1);
      
      //       // f fbar -> g Z/gamma
      //       SetMSUB(15,1);
      //       // f fbar -> gamma Z/gamma
      //       SetMSUB(19,1);
      //       // f g -> f Z/gamma
      //       SetMSUB(30,1);
      //       // f gamma -> f Z/gamma
      //       SetMSUB(35,1);
      
      //only Z included, not gamma
      SetMSTP(43,2);
      
      // Initial/final parton shower on (Pythia default)
      // With parton showers on we are generating "Z inclusive process"
      SetMSTP(61,1); //Initial QCD & QED showers on
      SetMSTP(71,1); //Final QCD & QED showers on
      
      break;  
    case kPyZgamma:
        //Inclusive production of Z
        SetMSEL(0);
        //f fbar -> Z/gamma
        SetMSUB(1,1);
        // Initial/final parton shower on (Pythia default)
        // With parton showers on we are generating "Z inclusive process"
        SetMSTP(61,1); //Initial QCD & QED showers on
        SetMSTP(71,1); //Final QCD & QED showers on
      break;  
      case kPyMBRSingleDiffraction:
      case kPyMBRDoubleDiffraction:
      case kPyMBRCentralDiffraction:
      break;
      case kPyJetsPWHG:
      //    N.B.
      //    ====
      //    For the case of jet production the following parameter setting
      //    limits the transverse momentum of secondary scatterings, due
      //    to multiple parton interactions, to be less than that of the
      //    primary interaction (see POWHEG Dijet paper arXiv:1012.3380
      //    [hep-ph] sec. 4.1 and also the PYTHIA Manual).
      SetMSTP(86,1);
      
      //    maximum number of errors before pythia aborts (def=10)
      SetMSTU(22,10);
      //    number of warnings printed on the shell
      SetMSTU(26,20);
      break;
    
      case kPyCharmPWHG:
      case kPyBeautyPWHG:
      case kPyWPWHG:
      //    number of warnings printed on the shell
      SetMSTU(26,20);
            
      break;
    }
//
//  Initialize PYTHIA
//
//  Select the tune
    if (itune > -1) {
      Pytune(itune);
      if (GetMSTP(192) > 1 || GetMSTP(193) > 1) {
	AliWarning(Form("Structure function for tune %5d set to %5s\n", 
			itune,  AliStructFuncType::PDFsetName(strucfunc).Data()));
	SetMSTP(52,2);
	SetMSTP(51, AliStructFuncType::PDFsetIndex(strucfunc));
      }
    }
//  
    SetMSTP(41,1);   // all resonance decays switched on
    if (process == kPyJetsPWHG || process == kPyCharmPWHG || process == kPyBeautyPWHG || process == kPyWPWHG) {
      Initialize("USER","","",0.);
    } else {	
      Initialize("CMS",fProjectile,fTarget,fEcms);
    }
    fOmegaDalitz.Init();
    fExodus.Init();
}

Int_t AliPythia::CheckedLuComp(Int_t kf)
{
// Check Lund particle code (for debugging)
    Int_t kc=Pycomp(kf);
    printf("\n Lucomp kf,kc %d %d",kf,kc);
    return kc;
}

void AliPythia::SetNuclei(Int_t a1, Int_t a2, Int_t pdf)
{
// Treat protons as inside nuclei with mass numbers a1 and a2  
//    The MSTP array in the PYPARS common block is used to enable and 
//    select the nuclear structure functions. 
//    MSTP(52)  : (D=1) choice of proton and nuclear structure-function library
//            =1: internal PYTHIA acording to MSTP(51) 
//            =2: PDFLIB proton  s.f., with MSTP(51)  = 1000xNGROUP+NSET
//    If the following mass number both not equal zero, nuclear corrections of the stf are used.
//    MSTP(192) : Mass number of nucleus side 1
//    MSTP(193) : Mass number of nucleus side 2
//    MSTP(194) : Nuclear structure function: 0: EKS98 8:EPS08 9:EPS09LO 19:EPS09NLO
    SetMSTP(52,2);
    SetMSTP(192, a1);
    SetMSTP(193, a2); 
    SetMSTP(194, pdf);
}
	

AliPythia* AliPythia::Instance()
{ 
// Set random number generator 
    if (fgAliPythia) {
	return fgAliPythia;
    } else {
	fgAliPythia = new AliPythia();
	return fgAliPythia;
    }
}

void AliPythia::PrintParticles()
{ 
// Print list of particl properties
    Int_t np = 0;
    char*   name = new char[16];    
    for (Int_t kf=0; kf<1000000; kf++) {
	for (Int_t c = 1;  c > -2; c-=2) {
	    Int_t kc = Pycomp(c*kf);
	    if (kc) {
		Float_t mass  = GetPMAS(kc,1);
		Float_t width = GetPMAS(kc,2);	
		Float_t tau   = GetPMAS(kc,4);

		Pyname(kf,name);
	
		np++;
		
		printf("\n mass, width, tau: %6d %s %10.3f %10.3e %10.3e", 
		       c*kf, name, mass, width, tau);
	    }
	}
    }
    printf("\n Number of particles %d \n \n", np);
}

void  AliPythia::ResetDecayTable()
{
//  Set default values for pythia decay switches
    Int_t i;
    for (i = 1; i <  501; i++) SetMDCY(i,1,fDefMDCY[i]);
    for (i = 1; i < 2001; i++) SetMDME(i,1,fDefMDME[i]);
}

void  AliPythia::SetDecayTable()
{
//  Set default values for pythia decay switches
//
    Int_t i;
    for (i = 1; i <  501; i++) fDefMDCY[i] = GetMDCY(i,1);
    for (i = 1; i < 2001; i++) fDefMDME[i] = GetMDME(i,1);
}

void  AliPythia::Pyclus(Int_t& njet)
{
//  Call Pythia clustering algorithm
//
    pyclus(njet);
}

void  AliPythia::Pycell(Int_t& njet)
{
//  Call Pythia jet reconstruction algorithm
//
    pycell(njet);
}

void  AliPythia::Pyshow(Int_t ip1, Int_t ip2, Double_t qmax)
{
//  Call Pythia jet reconstruction algorithm
//
    pyshow(ip1, ip2, qmax);
}

void AliPythia::Pyrobo(Int_t imi, Int_t ima, Double_t the, Double_t phi, Double_t bex, Double_t bey, Double_t bez)
{
    pyrobo(imi, ima, the, phi, bex, bey, bez);
}

void AliPythia::Pytune(Int_t itune)
{
/*
C
C ITUNE    NAME (detailed descriptions below)
C     0 Default : No settings changed => linked Pythia version's defaults.
C ====== Old UE, Q2-ordered showers ==========================================
C   100       A : Rick Field's CDF Tune A 
C   101      AW : Rick Field's CDF Tune AW
C   102      BW : Rick Field's CDF Tune BW
C   103      DW : Rick Field's CDF Tune DW
C   104     DWT : Rick Field's CDF Tune DW with slower UE energy scaling
C   105      QW : Rick Field's CDF Tune QW (NB: needs CTEQ6.1M pdfs externally)
C   106 ATLAS-DC2: Arthur Moraes' (old) ATLAS tune (ATLAS DC2 / Rome)
C   107     ACR : Tune A modified with annealing CR
C   108      D6 : Rick Field's CDF Tune D6 (NB: needs CTEQ6L pdfs externally)
C   109     D6T : Rick Field's CDF Tune D6T (NB: needs CTEQ6L pdfs externally)
C ====== Intermediate Models =================================================
C   200    IM 1 : Intermediate model: new UE, Q2-ordered showers, annealing CR
C   201     APT : Tune A modified to use pT-ordered final-state showers
C ====== New UE, interleaved pT-ordered showers, annealing CR ================
C   300      S0 : Sandhoff-Skands Tune 0 
C   301      S1 : Sandhoff-Skands Tune 1
C   302      S2 : Sandhoff-Skands Tune 2
C   303     S0A : S0 with "Tune A" UE energy scaling
C   304    NOCR : New UE "best try" without colour reconnections
C   305     Old : New UE, original (primitive) colour reconnections
C   306 ATLAS-CSC: Arthur Moraes' (new) ATLAS tune (needs CTEQ6L externally)
C ======= The Uppsala models =================================================
C   ( NB! must be run with special modified Pythia 6.215 version )
C   ( available from http://www.isv.uu.se/thep/MC/scigal/        )
C   400   GAL 0 : Generalized area-law model. Old parameters
C   401   SCI 0 : Soft-Colour-Interaction model. Old parameters
C   402   GAL 1 : Generalized area-law model. Tevatron MB retuned (Skands)
*/
    pytune(itune);
}

void AliPythia::Py2ent(Int_t idx, Int_t pdg1, Int_t pdg2, Double_t p){
  // Inset 2-parton system at line idx
  py2ent(idx, pdg1, pdg2, p);
}

void AliPythia::SetWeightPower(Double_t pow)
{
    setpowwght(pow);
    SetMSTP(142, 1); // Tell Pythia to use pyevwt to calculate event wghts
    if (GetCKIN(3) <= 0) 
      AliWarning("Need to set minimum p_T,hard to nonzero value for weighted event generation");
}

void AliPythia::InitQuenching(Float_t cMin, Float_t cMax, Float_t k, Int_t iECMethod, Float_t zmax, Int_t ngmax)
{
// Initializes 
// (1) The quenching model using quenching weights according to C. Salgado and U. Wiedemann
// (2) The nuclear geometry using the Glauber Model
//     
    
    fGlauber = AliFastGlauber::Instance();
    fGlauber->Init(2);
    fGlauber->SetCentralityClass(cMin, cMax); 

    fQuenchingWeights = new AliQuenchingWeights();
    fQuenchingWeights->InitMult();
    fQuenchingWeights->SetK(k);
    fQuenchingWeights->SetECMethod(AliQuenchingWeights::kECMethod(iECMethod));
    fNGmax = ngmax;
    fZmax  = zmax;
    
}


void  AliPythia::Quench()
{
//
//
//  Simple Jet Quenching routine:
//  =============================
//  The jet formed by all final state partons radiated by the parton created 
//  in the hard collisions is quenched by a factor (1-z) using light cone variables in 
//  the initial parton reference frame:
//  (E + p_z)new = (1-z) (E + p_z)old
//
//
//
//
//   The lost momentum is first balanced by one gluon with virtuality > 0.   
//   Subsequently the gluon splits to yield two gluons with E = p.
//
//
// 
    static Float_t eMean = 0.;
    static Int_t   icall = 0;
    
    Double_t p0[4][5];
    Double_t p1[4][5];
    Double_t p2[4][5];
    Int_t   klast[4] = {-1, -1, -1, -1};

    Int_t numpart   = fPyjets->N;
    Double_t px = 0., py = 0., pz = 0., e = 0., m = 0., p = 0., pt = 0., theta = 0., phi = 0.;
    Double_t pxq[4], pyq[4], pzq[4], eq[4], yq[4], mq[4], pq[4], phiq[4], thetaq[4], ptq[4];
    Bool_t  quenched[4];
    Double_t wjtKick[4] = {0., 0., 0., 0.};
    Int_t nGluon[4];
    Int_t qPdg[4];
    Int_t   imo, kst, pdg;
    
//
//  Sore information about Primary partons
//
//  j =
//  0, 1 partons from hard scattering
//  2, 3 partons from initial state radiation
// 
    for (Int_t i = 2; i <= 7; i++) {
	Int_t j = 0;
	// Skip gluons that participate in hard scattering
	if (i == 4 || i == 5) continue;
	// Gluons from hard Scattering
	if (i == 6 || i == 7) {
	    j = i - 6;
	    pxq[j]    = fPyjets->P[0][i];
	    pyq[j]    = fPyjets->P[1][i];
	    pzq[j]    = fPyjets->P[2][i];
	    eq[j]     = fPyjets->P[3][i];
	    mq[j]     = fPyjets->P[4][i];
	} else {
	    // Gluons from initial state radiation
	    //
	    // Obtain 4-momentum vector from difference between original parton and parton after gluon 
	    // radiation. Energy is calculated independently because initial state radition does not 
	    // conserve strictly momentum and energy for each partonic system independently.
	    //
	    // Not very clean. Should be improved !
	    //
	    //
	    j = i;
	    pxq[j]    = fPyjets->P[0][i] - fPyjets->P[0][i+2];
	    pyq[j]    = fPyjets->P[1][i] - fPyjets->P[1][i+2];
	    pzq[j]    = fPyjets->P[2][i] - fPyjets->P[2][i+2];
	    mq[j]     = fPyjets->P[4][i];
	    eq[j]     = TMath::Sqrt(pxq[j] * pxq[j] + pyq[j] * pyq[j] + pzq[j] * pzq[j] + mq[j] * mq[j]);
	}
//
//  Calculate some kinematic variables
//
	yq[j]     = 0.5 * TMath::Log((eq[j] + pzq[j] + 1.e-14) / (eq[j] - pzq[j] + 1.e-14));
	pq[j]     = TMath::Sqrt(pxq[j] * pxq[j] + pyq[j] * pyq[j] + pzq[j] * pzq[j]);
	phiq[j]   = TMath::Pi()+TMath::ATan2(-pyq[j], -pxq[j]);
	ptq[j]    = TMath::Sqrt(pxq[j] * pxq[j] + pyq[j] * pyq[j]);
	thetaq[j] = TMath::ATan2(ptq[j], pzq[j]);
	qPdg[j]   =  fPyjets->K[1][i];
    }
  
    Double_t int0[4];
    Double_t int1[4];
    
    fGlauber->GetI0I1ForPythiaAndXY(4, phiq, int0, int1, fXJet, fYJet, 15.);

    for (Int_t j = 0; j < 4; j++) {
	//
	// Quench only central jets and with E > 10.
	//


	Int_t itype = (qPdg[j] == 21) ? 2 : 1;
	Double_t eloss = fQuenchingWeights->GetELossRandomKFast(itype, int0[j], int1[j], eq[j]);

	if (TMath::Abs(yq[j]) > 2.5 || eq[j] < 10.) {
	    fZQuench[j] = 0.;
	} else {
	    if (eq[j] > 40. && TMath::Abs(yq[j]) < 0.5) {
		icall ++;
		eMean += eloss;
	    }
	    //
	    // Extra pt
	    Double_t l =   fQuenchingWeights->CalcLk(int0[j], int1[j]);	    
	    wjtKick[j] = TMath::Sqrt(l *  fQuenchingWeights->CalcQk(int0[j], int1[j]));
	    //
	    // Fractional energy loss
	    fZQuench[j] = eloss / eq[j];
	    //
	    // Avoid complete loss
	    //
	    if (fZQuench[j] > fZmax) fZQuench[j] = fZmax;
	    //
	    // Some debug printing

	    
//	    printf("Initial parton # %3d, Type %3d Energy %10.3f Phi %10.3f Length %10.3f Loss %10.3f Kick %10.3f Mean: %10.3f %10.3f\n", 
//		   j, itype, eq[j], phiq[j], l, eloss, wjtKick[j], eMean / Float_t(icall+1), yq[j]);
	    
//	    fZQuench[j] = 0.8;
//	    while (fZQuench[j] >= 0.95)  fZQuench[j] = gRandom->Exp(0.2);
	}
	
	quenched[j] = (fZQuench[j] > 0.01);
    } // primary partons
    
    

    Double_t pNew[1000][4];
    Int_t    kNew[1000];
    Int_t icount = 0;
    Double_t zquench[4];
    
//
//  System Loop    
    for (Int_t isys = 0; isys < 4; isys++) {
//      Skip to next system if not quenched.
	if (!quenched[isys]) continue;
	
	nGluon[isys]   = 1 + Int_t(fZQuench[isys] / (1. - fZQuench[isys]));
	if (nGluon[isys] > fNGmax) nGluon[isys] = fNGmax;
	zquench[isys] = 1. - TMath::Power(1. - fZQuench[isys], 1./Double_t(nGluon[isys]));
	wjtKick[isys]  = wjtKick[isys] / TMath::Sqrt(Double_t(nGluon[isys]));


	
	Int_t igMin = -1;
	Int_t igMax = -1;
	Double_t pg[4] = {0., 0., 0., 0.};
	
//
// Loop on radiation events

	for (Int_t iglu = 0; iglu < nGluon[isys]; iglu++) {
	    while (1) {
		icount = 0;
		for (Int_t k = 0; k < 4; k++)
		{
		    p0[isys][k] = 0.;
		    p1[isys][k] = 0.;
		    p2[isys][k] = 0.;
		}
//      Loop over partons
		for (Int_t i = 0; i < numpart; i++)
		{
		    imo =  fPyjets->K[2][i];
		    kst =  fPyjets->K[0][i];
		    pdg =  fPyjets->K[1][i];
		    
		
		
//      Quarks and gluons only
		    if (pdg != 21 && TMath::Abs(pdg) > 6) continue;
//      Particles from hard scattering only
		    
		    if (imo > 8 && imo < 1000) imo = fPyjets->K[2][imo - 1];
		    Int_t imom = imo % 1000;
		    if ((isys == 0 || isys == 1) && ((imom != (isys + 7)))) continue;
		    if ((isys == 2 || isys == 3) && ((imom != (isys + 1)))) continue;		    
		    
		    
//      Skip comment lines
		    if (kst != 1 && kst != 2) continue;
//
//      Parton kinematic
		    px    = fPyjets->P[0][i];
		    py    = fPyjets->P[1][i];
		    pz    = fPyjets->P[2][i];
		    e     = fPyjets->P[3][i];
		    m     = fPyjets->P[4][i];
		    pt    = TMath::Sqrt(px * px + py * py);
		    p     = TMath::Sqrt(px * px + py * py + pz * pz); 
		    phi   = TMath::Pi() + TMath::ATan2(-py, -px);
		    theta = TMath::ATan2(pt, pz);
		
//
//      Save 4-momentum sum for balancing
		    Int_t index = isys;
		    
		    p0[index][0] += px;
		    p0[index][1] += py;
		    p0[index][2] += pz;
		    p0[index][3] += e;
		
		    klast[index] = i;
		    
//
//      Fractional energy loss
		    Double_t z = zquench[index];
		    
		    
//      Don't fully quench radiated gluons
//
		    if (imo > 1000) {
//      This small factor makes sure that the gluons are not too close in phase space to avoid recombination
//

			z = 0.02;
		    }
//		    printf("z: %d %f\n", imo, z);
		    

//
		    
		    //
		    //
		    //      Transform into frame in which initial parton is along z-axis
		    //
		    TVector3 v(px, py, pz);
		    v.RotateZ(-phiq[index]);  v.RotateY(-thetaq[index]);
		    Double_t pxs = v.X(); Double_t pys = v.Y(); Double_t pl  = v.Z();

		    Double_t jt  = TMath::Sqrt(pxs * pxs + pys * pys);
		    Double_t mt2 = jt * jt + m * m;
		    Double_t zmax = 1.;	    
		    //
		    // Kinematic limit on z
		    //
		    if (m > 0.) zmax = 1. - m / TMath::Sqrt(m * m + jt * jt);
		    //
		    // Change light-cone kinematics rel. to initial parton
		    //	
		    Double_t eppzOld = e + pl;
		    Double_t empzOld = e - pl;
		    
		    Double_t eppzNew = (1. - z) * eppzOld;
		    Double_t empzNew = empzOld - mt2 * z / eppzOld;
		    Double_t eNew    = 0.5 * (eppzNew + empzNew);
		    Double_t plNew   = 0.5 * (eppzNew - empzNew);
		    
		    Double_t jtNew;
		    //
		    // if mt very small (or sometimes even < 0 for numerical reasons) set it to 0
		    Double_t mt2New = eppzNew * empzNew;
		    if (mt2New < 1.e-8) mt2New = 0.;
		    if (z < zmax) {
			if (m * m > mt2New) {
			    //
			    // This should not happen 
			    //
			    Fatal("Quench()", "This should never happen %e %e %e!", m, eppzNew, empzNew);
			    jtNew = 0;
			} else {
			    jtNew    = TMath::Sqrt(mt2New - m * m);
			}
		    } else {
			// If pT is to small (probably a leading massive particle) we scale only the energy
			// This can cause negative masses of the radiated gluon
			// Let's hope for the best ...
			jtNew = jt;
			eNew  = TMath::Sqrt(plNew * plNew + mt2);
			
		    }
		    //
		    //     Calculate new px, py
		    //
		    Double_t pxNew   = 0;
		    Double_t pyNew   = 0;
		    
		    if (jt>0) {
		      pxNew = jtNew / jt * pxs;
		      pyNew = jtNew / jt * pys;
		    }	
//		    Double_t dpx = pxs - pxNew;
//		    Double_t dpy = pys - pyNew;
//		    Double_t dpz = pl  - plNew;
//		    Double_t de  = e   - eNew;
//		    Double_t dmass2 = de * de  - dpx * dpx - dpy * dpy - dpz * dpz;
//		    printf("New mass (1) %e %e %e %e %e %e %e \n", dmass2, jt, jtNew, pl, plNew, e, eNew);
//		    printf("New mass (2) %e %e \n", pxNew, pyNew);
		    //
		    //      Rotate back
		    //	
		    TVector3 w(pxNew, pyNew, plNew);
		    w.RotateY(thetaq[index]); w.RotateZ(phiq[index]);
		    pxNew = w.X(); pyNew = w.Y(); plNew = w.Z();
		
		    p1[index][0] += pxNew;
		    p1[index][1] += pyNew;
		    p1[index][2] += plNew;
		    p1[index][3] += eNew;	
		    //
		    // Updated 4-momentum vectors
		    //
		    pNew[icount][0]  = pxNew;
		    pNew[icount][1]  = pyNew;
		    pNew[icount][2]  = plNew;
		    pNew[icount][3]  = eNew;
		    kNew[icount]     = i;
		    icount++;
		} // parton loop
		//
		// Check if there was phase-space for quenching
		//

		if (icount == 0) quenched[isys] = kFALSE;
		if (!quenched[isys]) break;
		
		for (Int_t j = 0; j < 4; j++) 
		{
		    p2[isys][j] = p0[isys][j] - p1[isys][j];
		}
		p2[isys][4] = p2[isys][3] * p2[isys][3] - p2[isys][0] * p2[isys][0] - p2[isys][1] * p2[isys][1] - p2[isys][2] * p2[isys][2];
		if (p2[isys][4] > 0.) {
		    p2[isys][4] = TMath::Sqrt(p2[isys][4]);
		    break;
		} else {
		    printf("Warning negative mass squared in system %d %f ! \n", isys, zquench[isys]);
		    printf("4-Momentum: %10.3e %10.3e %10.3e %10.3e %10.3e \n", p2[isys][0], p2[isys][1], p2[isys][2], p2[isys][3], p2[isys][4]);
		    if (p2[isys][4] < -0.01) {
			printf("Negative mass squared !\n");
			// Here we have to put the gluon back to mass shell
			// This will lead to a small energy imbalance
			p2[isys][4]  = 0.;
			p2[isys][3]  = TMath::Sqrt(p2[isys][0] * p2[isys][0] + p2[isys][1] * p2[isys][1] + p2[isys][2] * p2[isys][2]);
			break;
		    } else {
			p2[isys][4] = 0.;
			break;
		    }
		}
		/*
		zHeavy *= 0.98;
		printf("zHeavy lowered to %f\n", zHeavy);
		if (zHeavy < 0.01) {
		    printf("No success ! \n");
		    icount = 0;
		    quenched[isys] = kFALSE;
		    break;
		}
		*/
	    } // iteration on z (while)
 	    
//	    Update  event record
	    for (Int_t k = 0; k < icount; k++) {
//		printf("%6d %6d %10.3e %10.3e %10.3e %10.3e\n", k, kNew[k], pNew[k][0],pNew[k][1], pNew[k][2], pNew[k][3] );
		fPyjets->P[0][kNew[k]] = pNew[k][0];
		fPyjets->P[1][kNew[k]] = pNew[k][1];
		fPyjets->P[2][kNew[k]] = pNew[k][2];
		fPyjets->P[3][kNew[k]] = pNew[k][3];
	    }
	    //
	    // Add the gluons
	    //
	    Int_t ish = 0;    
	    Int_t iGlu;
	    if (!quenched[isys]) continue;
//
//      Last parton from shower i
	    Int_t in = klast[isys];
//
//      Continue if no parton in shower i selected
	    if (in == -1) continue;
//  
//      If this is the second initial parton and it is behind the first move pointer by previous ish
	    if (isys == 1 && klast[1] > klast[0]) in += ish;
//
//      Starting index
	    
//	    jmin = in - 1;
// How many additional gluons will be generated
	    ish  = 1;
	    if (p2[isys][4] > 0.05) ish = 2;
//
//      Position of gluons
	    iGlu = numpart;
	    if (iglu == 0) igMin = iGlu;
	    igMax = iGlu;
	    numpart += ish;
	    (fPyjets->N) += ish;
	    
	    if (ish == 1) {
		fPyjets->P[0][iGlu] = p2[isys][0];
		fPyjets->P[1][iGlu] = p2[isys][1];
		fPyjets->P[2][iGlu] = p2[isys][2];
		fPyjets->P[3][iGlu] = p2[isys][3];
		fPyjets->P[4][iGlu] = p2[isys][4];
		
		fPyjets->K[0][iGlu] = 1;
		if (iglu == nGluon[isys] - 1) fPyjets->K[0][iGlu] = 1;
		fPyjets->K[1][iGlu] = 21;	
		fPyjets->K[2][iGlu] = fPyjets->K[2][in] + 1000;
		fPyjets->K[3][iGlu] = -1;	
		fPyjets->K[4][iGlu] = -1;
		
		pg[0] += p2[isys][0];
		pg[1] += p2[isys][1];
		pg[2] += p2[isys][2];
		pg[3] += p2[isys][3];
	    } else {
		//
		// Split gluon in rest frame.
		//
		Double_t bx   =  p2[isys][0] / p2[isys][3];
		Double_t by   =  p2[isys][1] / p2[isys][3];
		Double_t bz   =  p2[isys][2] / p2[isys][3];
		Double_t pst  =  p2[isys][4] / 2.;
		//
		// Isotropic decay ????
		Double_t cost = 2. * gRandom->Rndm() - 1.;
		Double_t sint = TMath::Sqrt((1.-cost)*(1.+cost));
		Double_t phis =  2. * TMath::Pi() * gRandom->Rndm();
		
		Double_t pz1 =   pst * cost;
		Double_t pz2 =  -pst * cost;
		Double_t pt1 =   pst * sint;
		Double_t pt2 =  -pst * sint;
		Double_t px1 =   pt1 * TMath::Cos(phis);
		Double_t py1 =   pt1 * TMath::Sin(phis);	    
		Double_t px2 =   pt2 * TMath::Cos(phis);
		Double_t py2 =   pt2 * TMath::Sin(phis);	    
		
		fPyjets->P[0][iGlu] = px1;
		fPyjets->P[1][iGlu] = py1;
		fPyjets->P[2][iGlu] = pz1;
		fPyjets->P[3][iGlu] = pst;
		fPyjets->P[4][iGlu] = 0.;
		
		fPyjets->K[0][iGlu] = 1 ;
		fPyjets->K[1][iGlu] = 21;	
		fPyjets->K[2][iGlu] = fPyjets->K[2][in] + 1000;
		fPyjets->K[3][iGlu] = -1;	
		fPyjets->K[4][iGlu] = -1;
		
		fPyjets->P[0][iGlu+1] = px2;
		fPyjets->P[1][iGlu+1] = py2;
		fPyjets->P[2][iGlu+1] = pz2;
		fPyjets->P[3][iGlu+1] = pst;
		fPyjets->P[4][iGlu+1] = 0.;
		
		fPyjets->K[0][iGlu+1] = 1;
		if (iglu == nGluon[isys] - 1) fPyjets->K[0][iGlu+1] = 1;
		fPyjets->K[1][iGlu+1] = 21;	
		fPyjets->K[2][iGlu+1] = fPyjets->K[2][in] + 1000;
		fPyjets->K[3][iGlu+1] = -1;	
		fPyjets->K[4][iGlu+1] = -1;
		SetMSTU(1,0);
		SetMSTU(2,0);
		//
		// Boost back
		//
		Pyrobo(iGlu + 1, iGlu + 2, 0., 0., bx, by, bz);
	    }
/*    
	    for (Int_t ig = iGlu; ig < iGlu+ish; ig++) {
		Double_t px, py, pz;
		px = fPyjets->P[0][ig]; 
		py = fPyjets->P[1][ig]; 
		pz = fPyjets->P[2][ig]; 
		TVector3 v(px, py, pz);
		v.RotateZ(-phiq[isys]);
		v.RotateY(-thetaq[isys]);
		Double_t pxs     = v.X(); Double_t pys = v.Y(); Double_t pzs  = v.Z();	   
		Double_t r       = AliPythiaRndm::GetPythiaRandom()->Rndm();
		Double_t jtKick  = 0.3 * TMath::Sqrt(-TMath::Log(r));
		if (ish == 2)   jtKick  = wjtKick[i] * TMath::Sqrt(-TMath::Log(r)) / TMath::Sqrt(2.);
		Double_t phiKick = 2. * TMath::Pi() * AliPythiaRndm::GetPythiaRandom()->Rndm();
		pxs += jtKick * TMath::Cos(phiKick);
		pys += jtKick * TMath::Sin(phiKick);
		TVector3 w(pxs, pys, pzs);
		w.RotateY(thetaq[isys]);
		w.RotateZ(phiq[isys]);
		fPyjets->P[0][ig] = w.X(); 
		fPyjets->P[1][ig] = w.Y(); 
		fPyjets->P[2][ig] = w.Z(); 
		fPyjets->P[2][ig] = w.Mag();
	    }
*/
	} // kGluon	    
	
	
    // Check energy conservation
	Double_t pxs = 0.;
	Double_t pys = 0.;
	Double_t pzs = 0.;	
	Double_t es  = 14000.;
	
	for (Int_t i = 0; i < numpart; i++)
	{
	    kst =  fPyjets->K[0][i];
	    if (kst != 1 && kst != 2) continue;
	    pxs += fPyjets->P[0][i];
	    pys += fPyjets->P[1][i];
	    pzs += fPyjets->P[2][i];	    
	    es  -= fPyjets->P[3][i];	    
	}
	if (TMath::Abs(pxs) > 1.e-2 ||
	    TMath::Abs(pys) > 1.e-2 ||
	    TMath::Abs(pzs) > 1.e-1) {
	    printf("%e %e %e %e\n", pxs, pys, pzs, es);
//		Fatal("Quench()", "4-Momentum non-conservation");
	}
	
    } // end quenching loop (systems)
// Clean-up
    for (Int_t i = 0; i < numpart; i++)
    {
	imo =  fPyjets->K[2][i];
	if (imo > 1000) {
	    fPyjets->K[2][i] = fPyjets->K[2][i] % 1000;
	}
    }
//	this->Pylist(1);
} // end quench


void AliPythia::Pyquen(Double_t a, Int_t ibf, Double_t b)
{
    // Igor Lokthine's quenching routine
    // http://lokhtin.web.cern.ch/lokhtin/pyquen/pyquen.txt

    pyquen(a, ibf, b);
}

void AliPythia::SetPyquenParameters(Double_t t0, Double_t tau0, Int_t nf, Int_t iengl, Int_t iangl)
{
    // Set the parameters for the PYQUEN package.
    // See comments in PyquenCommon.h
    
    
    PYQPAR.t0    = t0;
    PYQPAR.tau0  = tau0;
    PYQPAR.nf    = nf;
    PYQPAR.iengl = iengl;
    PYQPAR.iangl = iangl;
}


void AliPythia::Pyevnw()
{
    // New multiple interaction scenario
    pyevnw();
}

void  AliPythia::Pyshowq(Int_t ip1, Int_t ip2, Double_t qmax)
{
    //  Call medium-modified Pythia jet reconstruction algorithm
    //
    pyshowq(ip1, ip2, qmax);
}
 void AliPythia::Qpygin0()                                                                      
 {                                                                                              
     // New multiple interaction scenario                                                       
     qpygin0();                                                                                 
 }

void AliPythia::GetQuenchingParameters(Double_t& xp, Double_t& yp, Double_t z[4])
{
    // Return event specific quenching parameters
    xp = fXJet;
    yp = fYJet;
    for (Int_t i = 0; i < 4; i++) z[i] = fZQuench[i];

}

void AliPythia::ConfigHeavyFlavor()
{
    //
    // Default configuration for Heavy Flavor production
    //
    // All QCD processes
    //
    SetMSEL(1);
    

    if (fItune < 0) {
      // No multiple interactions
      SetMSTP(81,0);
      SetPARP(81, 0.);
      SetPARP(82, 0.);    
    }
    // Initial/final parton shower on (Pythia default)
    SetMSTP(61,1);
    SetMSTP(71,1);
    
    // 2nd order alpha_s
    SetMSTP(2,2);
    
    // QCD scales
    SetMSTP(32,2);
    SetPARP(34,1.0);
}

void AliPythia::AtlasTuning()
{
    //
    // Configuration for the ATLAS tuning
    if (fItune > -1) return;
    printf("ATLAS TUNE \n");
    
    SetMSTP(51, AliStructFuncType::PDFsetIndex(kCTEQ5L));      // CTEQ5L pdf
    SetMSTP(81,1);             // Multiple Interactions ON
    SetMSTP(82,4);             // Double Gaussian Model
    SetPARP(81,1.9);           // Min. pt for multiple interactions (default in 6.2-14) 
    SetPARP(82,1.8);           // [GeV]    PT_min at Ref. energy
    SetPARP(89,1000.);         // [GeV]   Ref. energy
    SetPARP(90,0.16);          // 2*epsilon (exponent in power law)
    SetPARP(83,0.5);           // Core density in proton matter distribution (def.value)
    SetPARP(84,0.5);           // Core radius
    SetPARP(85,0.33);          // Regulates gluon prod. mechanism
    SetPARP(86,0.66);          // Regulates gluon prod. mechanism
    SetPARP(67,1);             // Regulates Initial State Radiation
}

void AliPythia::AtlasTuningMC09()
{
    //
    // Configuration for the ATLAS tuning
    if (fItune > -1) return;
    printf("ATLAS New TUNE MC09\n");
    SetMSTP(81,21);             // treatment for MI, ISR, FSR and beam remnants: MI on, new model
    SetMSTP(82, 4);             // Double Gaussian Model
    SetMSTP(52, 2);             // External PDF
    SetMSTP(51, 20650);         // MRST LO*
  
    
    SetMSTP(70, 0);             // (was 2: def manual 1, def code 0) virtuality scale for ISR 
    SetMSTP(72, 1);             // (was 0: def 1) maximum scale for FSR
    SetMSTP(88, 1);             // (was 0: def 1) strategy for qq junction to di-quark or baryon in beam remnant
    SetMSTP(90, 0);             // (was 1: def 0) strategy of compensate the primordial kT

    SetPARP(78, 0.3);           // the amount of color reconnection in the final state
    SetPARP(80, 0.1);           // probability of color partons kicked out from beam remnant
    SetPARP(82, 2.3);           // [GeV]    PT_min at Ref. energy    
    SetPARP(83, 0.8);           // Core density in proton matter distribution (def.value)    
    SetPARP(84, 0.7);           // Core radius
    SetPARP(90, 0.25);          //  2*epsilon (exponent in power law)
    SetPARJ(81, 0.29);          // (was 0.14: def 0.29) Labmda value in running alpha_s for parton showers

    SetMSTP(95, 6);
    SetPARJ(41, 0.3);           // a and b parameters of the symmm. Lund FF
    SetPARJ(42, 0.58);
    SetPARJ(46, 0.75);          // mod. of the Lund FF for heavy end-point quarks
    SetPARP(89,1800.);         // [GeV]   Ref. energy
}

AliPythia& AliPythia::operator=(const  AliPythia& rhs)
{
// Assignment operator
    rhs.Copy(*this);
    return *this;
}

 void AliPythia::Copy(TObject&) const
{
    //
    // Copy 
    //
    Fatal("Copy","Not implemented!\n");
}

void AliPythia::DalitzDecays()
{

  //
  // Replace all omega dalitz decays with the correct matrix element decays
  //
  Int_t nt = fPyjets->N;
  for (Int_t i = 0; i < nt; i++) {
    if (fPyjets->K[1][i] != 223)                       continue;	
    Int_t fd = fPyjets->K[3][i] - 1;
    Int_t ld = fPyjets->K[4][i] - 1;
    if (fd < 0)                                        continue;
    if ((ld - fd) != 2)                                continue;
    if ((fPyjets->K[1][fd] != 111) ||
	((TMath::Abs(fPyjets->K[1][fd+1]) != 11) && (TMath::Abs(fPyjets->K[1][fd+1]) != 13)))
	continue;
    TLorentzVector omega(fPyjets->P[0][i], fPyjets->P[1][i], fPyjets->P[2][i], fPyjets->P[3][i]);
    Int_t pdg = TMath::Abs(fPyjets->K[1][fd+1]);
    fOmegaDalitz.Decay(pdg, &omega);
    for (Int_t j = 0; j < 3; j++) {
      for (Int_t k = 0; k < 4; k++) {
	TLorentzVector vec = (fOmegaDalitz.Products())[2-j];
	fPyjets->P[k][fd+j] = vec[k];
      }
    }
  }
}


//
// Replace all dalitz(pi0,eta,omega,eta',phi) and resonance(rho,omega,phi,jpsi) decays with the correct matrix element decays
// for di-electron cocktail calculations
//


void AliPythia::PizeroDalitz()
{

  Int_t nt = fPyjets->N;
  for (Int_t i = 0; i < nt; i++) {
  if (fPyjets->K[1][i] != 111)                       continue;
  Int_t fd = fPyjets->K[3][i] - 1;
  Int_t ld = fPyjets->K[4][i] - 1;
  if (fd < 0)                                        continue;
  if ((ld - fd) != 2)                                continue;
  if ((fPyjets->K[1][fd] != 22) || (TMath::Abs(fPyjets->K[1][fd+1]) != 11) )
  continue;
  TLorentzVector pizero(fPyjets->P[0][i], fPyjets->P[1][i], fPyjets->P[2][i], fPyjets->P[3][i]);
  Int_t pdg = TMath::Abs(fPyjets->K[1][i]);
  fExodus.Decay(pdg, &pizero);
  for (Int_t j = 0; j < 3; j++) {
  for (Int_t k = 0; k < 4; k++) {
  TLorentzVector vec = (fExodus.Products_pion())[2-j];
  fPyjets->P[k][fd+j] = vec[k];
        }
      }
  }
}


void AliPythia::EtaDalitz()
{

  Int_t nt = fPyjets->N;
  for (Int_t i = 0; i < nt; i++) {
  if (fPyjets->K[1][i] != 221)                       continue;
  Int_t fd = fPyjets->K[3][i] - 1;
  Int_t ld = fPyjets->K[4][i] - 1;
  if (fd < 0)                                        continue;
  if ((ld - fd) != 2)                                continue;
  if ((fPyjets->K[1][fd] != 22) || (TMath::Abs(fPyjets->K[1][fd+1]) != 11))
  continue;
  TLorentzVector eta(fPyjets->P[0][i], fPyjets->P[1][i], fPyjets->P[2][i], fPyjets->P[3][i]);
  Int_t pdg = TMath::Abs(fPyjets->K[1][i]);
  fExodus.Decay(pdg, &eta);
  for (Int_t j = 0; j < 3; j++) {
  for (Int_t k = 0; k < 4; k++) {
  TLorentzVector vec = (fExodus.Products_eta())[2-j];
  fPyjets->P[k][fd+j] = vec[k];
        }
     }
  }
}


void AliPythia::RhoDirect()
{

  Int_t nt = fPyjets->N;
  for (Int_t i = 0; i < nt; i++) {
  if (fPyjets->K[1][i] != 113)                       continue;
  Int_t fd = fPyjets->K[3][i] - 1;
  Int_t ld = fPyjets->K[4][i] - 1;
  if (fd < 0)                                        continue;
  if ((ld - fd) != 1)                                continue;
  if ((TMath::Abs(fPyjets->K[1][fd]) != 11))
  continue;
  TLorentzVector rho(fPyjets->P[0][i], fPyjets->P[1][i], fPyjets->P[2][i], fPyjets->P[3][i]);
  Int_t pdg = TMath::Abs(fPyjets->K[1][i]);
  fExodus.Decay(pdg, &rho);
  for (Int_t j = 0; j < 2; j++) {
  for (Int_t k = 0; k < 4; k++) {
  TLorentzVector vec = (fExodus.Products_rho())[1-j];
  fPyjets->P[k][fd+j] = vec[k];
        }
     }
  }
}


void AliPythia::OmegaDalitz()
{

  Int_t nt = fPyjets->N;
  for (Int_t i = 0; i < nt; i++) {
  if (fPyjets->K[1][i] != 223)                       continue;
  Int_t fd = fPyjets->K[3][i] - 1;
  Int_t ld = fPyjets->K[4][i] - 1;
  if (fd < 0)                                        continue;
  if ((ld - fd) != 2)                                continue;
  if ((fPyjets->K[1][fd] != 111) || (TMath::Abs(fPyjets->K[1][fd+1]) != 11))
  continue;
  TLorentzVector omegadalitz(fPyjets->P[0][i], fPyjets->P[1][i], fPyjets->P[2][i], fPyjets->P[3][i]);
  Int_t pdg = TMath::Abs(fPyjets->K[1][i]);
  fExodus.Decay(pdg, &omegadalitz);
  for (Int_t j = 0; j < 3; j++) {
  for (Int_t k = 0; k < 4; k++) {
  TLorentzVector vec = (fExodus.Products_omega_dalitz())[2-j];
  fPyjets->P[k][fd+j] = vec[k];
        }
     }
  }
}


void AliPythia::OmegaDirect()
{

  Int_t nt = fPyjets->N;
  for (Int_t i = 0; i < nt; i++) {
  if (fPyjets->K[1][i] != 223)                       continue;
  Int_t fd = fPyjets->K[3][i] - 1;
  Int_t ld = fPyjets->K[4][i] - 1;
  if (fd < 0)                                        continue;
  if ((ld - fd) != 1)                                continue;
  if ((TMath::Abs(fPyjets->K[1][fd]) != 11))
  continue;
  TLorentzVector omegadirect(fPyjets->P[0][i], fPyjets->P[1][i], fPyjets->P[2][i], fPyjets->P[3][i]);
  Int_t pdg = TMath::Abs(fPyjets->K[1][i]);
  fExodus.Decay(pdg, &omegadirect);
  for (Int_t j = 0; j < 2; j++) {
  for (Int_t k = 0; k < 4; k++) {
  TLorentzVector vec = (fExodus.Products_omega())[1-j];
  fPyjets->P[k][fd+j] = vec[k];
        }
     }
  }
}


void AliPythia::EtaprimeDalitz()
{

  Int_t nt = fPyjets->N;
  for (Int_t i = 0; i < nt; i++) {
  if (fPyjets->K[1][i] != 331)                       continue;
  Int_t fd = fPyjets->K[3][i] - 1;
  Int_t ld = fPyjets->K[4][i] - 1;
  if (fd < 0)                                        continue;
  if ((ld - fd) != 2)                                continue;
  if ((fPyjets->K[1][fd] != 22) || (TMath::Abs(fPyjets->K[1][fd+1]) != 11))
  continue;
  TLorentzVector etaprime(fPyjets->P[0][i], fPyjets->P[1][i], fPyjets->P[2][i], fPyjets->P[3][i]);
  Int_t pdg = TMath::Abs(fPyjets->K[1][i]);
  fExodus.Decay(pdg, &etaprime);
  for (Int_t j = 0; j < 3; j++) {
  for (Int_t k = 0; k < 4; k++) {
  TLorentzVector vec = (fExodus.Products_etaprime())[2-j];
  fPyjets->P[k][fd+j] = vec[k];
        }
     }
  }
}


void AliPythia::PhiDalitz()
{

  Int_t nt = fPyjets->N;
  for (Int_t i = 0; i < nt; i++) {
  if (fPyjets->K[1][i] != 333)                       continue;
  Int_t fd = fPyjets->K[3][i] - 1;
  Int_t ld = fPyjets->K[4][i] - 1;
  if (fd < 0)                                        continue;
  if ((ld - fd) != 2)                                continue;
  if ((fPyjets->K[1][fd] != 221) || (TMath::Abs(fPyjets->K[1][fd+1]) != 11))
  continue;
  TLorentzVector phidalitz(fPyjets->P[0][i], fPyjets->P[1][i], fPyjets->P[2][i], fPyjets->P[3][i]);
  Int_t pdg = TMath::Abs(fPyjets->K[1][i]);
  fExodus.Decay(pdg, &phidalitz);
  for (Int_t j = 0; j < 3; j++) {
  for (Int_t k = 0; k < 4; k++) {
  TLorentzVector vec = (fExodus.Products_phi_dalitz())[2-j];
  fPyjets->P[k][fd+j] = vec[k];
        }
     }
  }
}


void AliPythia::PhiDirect()
{

  Int_t nt = fPyjets->N;
  for (Int_t i = 0; i < nt; i++) {
  if (fPyjets->K[1][i] != 333)                       continue;
  Int_t fd = fPyjets->K[3][i] - 1;
  Int_t ld = fPyjets->K[4][i] - 1;
  if (fd < 0)                                        continue;
  if ((ld - fd) != 1)                                continue;
  if ((TMath::Abs(fPyjets->K[1][fd]) != 11))
  continue;
  TLorentzVector phi(fPyjets->P[0][i], fPyjets->P[1][i], fPyjets->P[2][i], fPyjets->P[3][i]);
  Int_t pdg = TMath::Abs(fPyjets->K[1][i]);
  fExodus.Decay(pdg, &phi);
  for (Int_t j = 0; j < 2; j++) {
  for (Int_t k = 0; k < 4; k++) {
  TLorentzVector vec = (fExodus.Products_phi())[1-j];
  fPyjets->P[k][fd+j] = vec[k];
        }
     }
  }
}

void AliPythia::JPsiDirect()
{

  Int_t nt = fPyjets->N;
  for (Int_t i = 0; i < nt; i++) {
  if (fPyjets->K[1][i] != 443)                       continue;
  Int_t fd = fPyjets->K[3][i] - 1;
  Int_t ld = fPyjets->K[4][i] - 1;
  if (fd < 0)                                        continue;
  if ((ld - fd) != 1)                                continue;
  if ((TMath::Abs(fPyjets->K[1][fd]) != 11))
  continue;
  TLorentzVector jpsi(fPyjets->P[0][i], fPyjets->P[1][i], fPyjets->P[2][i], fPyjets->P[3][i]);
  Int_t pdg = TMath::Abs(fPyjets->K[1][i]);
  fExodus.Decay(pdg, &jpsi);
  for (Int_t j = 0; j < 2; j++) {
  for (Int_t k = 0; k < 4; k++) {
  TLorentzVector vec = (fExodus.Products_jpsi())[1-j];
  fPyjets->P[k][fd+j] = vec[k];
        }
     }
  }
}

Int_t AliPythia::GetNMPI()
{
  // returns the number of parton-parton interactions
  return (GetMSTI(31));
}
