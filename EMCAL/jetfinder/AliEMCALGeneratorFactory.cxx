/**************************************************************************
 * Copyright(c) 1998-2002, ALICE Experiment at CERN, All rights reserved. *
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

/*
$Log$
Revision 1.1  2005/08/01 16:11:18  pavlinov
"Version

*/
#include "assert.h"
#include "TString.h"
#include "AliEMCALGeneratorFactory.h"
#include "AliGenerator.h"
#include "AliGenFixed.h"
#include "AliGenBox.h"
#include "AliGenHIJINGpara.h"
#include "AliGenHIJINGparaBa.h"
#include "AliGenHijing.h"
#include "AliGenCocktail.h"
#include "AliGenPythia.h"

//*-- Authors: Aleksei Pavlinov (WSU)
//*   Initial variant is in Config.C for EMCAL production.
ClassImp(AliEMCALGeneratorFactory)

AliEMCALGeneratorFactory::AliEMCALGeneratorFactory(PprRunFact_t run, PprRadFact_t rad) 
{
    fGenerator = fBgGenerator = fSignalGenerator = 0;
    fRunType   = run;
    fRadiation = rad;
    fMomentum  = 0;

    Int_t isw = 3;
    if (rad == kNoGluonRadiation) isw = 0;    
    Float_t thmin=0, thmax=0;

    AliGenFixed *genGG=0;
    AliGenBox   *genB=0;
    AliGenHIJINGparaBa *genHijParaB=0, *bg=0;
    AliGenHIJINGpara   *genHijPara=0;
    AliGenHijing       *genHij=0;
    AliGenCocktail     *genCoct=0;
    AliGenPythia       *genPy=0, *jets=0;
    //    AliPythia          *aliPy = 0;

    fComment = new TString;

    switch (run) {
    case kGammaGun:
        genGG = new AliGenFixed(1);
	genGG->SetMomentum(100.);
	genGG->SetPhi(240.);
	genGG->SetTheta(91.);
	genGG->SetPart(kGamma);
	fGenerator = (AliGenerator*)genGG;
     break;

    case kGammaBox:
	genB = new AliGenBox(100);
	genB->SetMomentumRange(100., 100.1);
	genB->SetPhiRange(0,360);
	genB->SetThetaRange(45., 135.);
	genB->SetPart(kGamma);
        fGenerator = (AliGenerator*)genB;
	break;
	
    case kTest50:
	genHijParaB = new AliGenHIJINGparaBa(50);
	genHijParaB->SetMomentumRange(0, 999999.);
	genHijParaB->SetPhiRange(-180., 180.);
	// Set pseudorapidity range from -8 to 8.
	thmin = EtaToTheta(8);   // theta min. <---> eta max
	thmax = EtaToTheta(-8);  // theta max. <---> eta min 
	genHijParaB->SetThetaRange(thmin,thmax);
        fGenerator = (AliGenerator*)genHijParaB;
	break;

    case kParam_8000:
	//coment= fComment.Append(":HIJINGparam N=8000");
	genHijPara = new AliGenHIJINGpara(86030);
	genHijPara->SetMomentumRange(0, 999999.);
	genHijPara->SetPhiRange(-180., 180.);
	// Set pseudorapidity range from -8 to 8.
	thmin = EtaToTheta(8);   // theta min. <---> eta max
	thmax = EtaToTheta(-8);  // theta max. <---> eta min 
	genHijPara->SetThetaRange(thmin,thmax);
        fGenerator = (AliGenerator*)genHijPara;
	break;
    case kParam_4000:
	genHijPara = new AliGenHIJINGpara(43015);
	genHijPara->SetMomentumRange(0, 999999.);
	genHijPara->SetPhiRange(-180., 180.);
	// Set pseudorapidity range from -8 to 8.
	thmin = EtaToTheta(8);   // theta min. <---> eta max
	thmax = EtaToTheta(-8);  // theta max. <---> eta min 
	genHijPara->SetThetaRange(thmin,thmax);
        fGenerator = (AliGenerator*)genHijPara;
	break;
    case kParam_2000:
        (*fComment) = "HIJINGparam N=2000";
	genHijPara = new AliGenHIJINGpara(21507);
	genHijPara->SetMomentumRange(0, 999999.);
	genHijPara->SetPhiRange(-180., 180.);
	// Set pseudorapidity range from -8 to 8.
	thmin = EtaToTheta(8);   // theta min. <---> eta max
	thmax = EtaToTheta(-8);  // theta max. <---> eta min 
	genHijPara->SetThetaRange(thmin,thmax);
        fGenerator = (AliGenerator*)genHijPara;
	break;

    case kParam_8000_Ecal:
	genHijParaB = new AliGenHIJINGparaBa(82534);
	genHijParaB->SetMomentumRange(0, 999999.);
	genHijParaB->SetPhiRange(-180., 180.);
	// Set pseudorapidity range from -8 to 8.
	thmin = EtaToTheta( 5);   // theta min. <---> eta max
	thmax = EtaToTheta(-5);  // theta max. <---> eta min 
	genHijParaB->SetThetaRange(thmin,thmax);
        fGenerator = (AliGenerator*)genHijParaB;
	break;

    case kParam_4000_Ecal:
	genHijParaB = new AliGenHIJINGparaBa(82534/2);
	genHijParaB->SetMomentumRange(0, 999999.);
	genHijParaB->SetPhiRange(-180., 180.);
	// Set pseudorapidity range from -8 to 8.
	thmin = EtaToTheta( 5);   // theta min. <---> eta max
	thmax = EtaToTheta(-5);  // theta max. <---> eta min 
	genHijParaB->SetThetaRange(thmin,thmax);
        fGenerator = (AliGenerator*)genHijParaB;
	break;
//
//  Hijing Central
//
    case kHijing_cent1:
	genHij = HijingStandard();
// impact parameter range
	genHij->SetImpactParameterRange(0., 5.);
        fGenerator = (AliGenerator*)genHij;
	break;
    case kHijing_cent2:
	genHij = HijingStandard();
// impact parameter range
	genHij->SetImpactParameterRange(0., 2.);
        fGenerator = (AliGenerator*)genHij;
	break;
//
// Hijing Peripheral 
//
    case kHijing_per1:
	genHij = HijingStandard();
// impact parameter range
	genHij->SetImpactParameterRange(5., 8.6);
        fGenerator = (AliGenerator*)genHij;
	break;
    case kHijing_per2:
	//coment= comment.Append("HIJING per2");
	genHij = HijingStandard();
// impact parameter range
	genHij->SetImpactParameterRange(8.6, 11.2);
        fGenerator = (AliGenerator*)genHij;
	break;
    case kHijing_per3:
	genHij = HijingStandard();
// impact parameter range
	genHij->SetImpactParameterRange(11.2, 13.2);
        fGenerator = (AliGenerator*)genHij;
	break;
    case kHijing_per4:
	genHij = HijingStandard();
// impact parameter range
	genHij->SetImpactParameterRange(13.2, 15.);
        fGenerator = (AliGenerator*)genHij;
	break;
    case kHijing_per5:
	//coment= comment.Append("HIJING per5");
	genHij = HijingStandard();
// impact parameter range
	genHij->SetImpactParameterRange(15., 100.);
        fGenerator = (AliGenerator*)genHij;
	break;
//
//  Jet-Jet
//
    case kHijing_jj25:
	//coment= comment.Append("HIJING Jet 25 GeV");
	genHij = HijingStandard();
// impact parameter range
	genHij->SetImpactParameterRange(0., 5.);
	// trigger
	genHij->SetTrigger(1);
	genHij->SetPtJet(25.);
	genHij->SetSimpleJets(1);
	genHij->SetRadiation(isw);
	genHij->SetJetEtaRange(-0.3,0.3);
	genHij->SetJetPhiRange(15.,105.);   
        fGenerator = (AliGenerator*)genHij;
	break;

    case kHijing_jj50:
	//coment= comment.Append("HIJING Jet 50 GeV");
	genHij = HijingStandard();
// impact parameter range
	genHij->SetImpactParameterRange(0., 5.);
	// trigger
	genHij->SetTrigger(1);
	genHij->SetPtJet(50.);
	genHij->SetSimpleJets(1);
	genHij->SetRadiation(isw);
	genHij->SetJetEtaRange(-0.3,0.3);
	genHij->SetJetPhiRange(15.,105.);   
        fGenerator = (AliGenerator*)genHij;
	break;

    case kHijing_jj75:
	//coment= comment.Append("HIJING Jet 75 GeV");
	genHij = HijingStandard();
// impact parameter range
	genHij->SetImpactParameterRange(0., 5.);
	// trigger
	genHij->SetTrigger(1);
	genHij->SetPtJet(75.);
	genHij->SetSimpleJets(1);
	genHij->SetRadiation(isw);
	genHij->SetJetEtaRange(-0.3,0.3);
	genHij->SetJetPhiRange(15.,105.);   
        fGenerator = (AliGenerator*)genHij;
	break;

    case kHijing_jj100:
	//coment= comment.Append("HIJING Jet 100 GeV");
	genHij = HijingStandard();
// impact parameter range
	genHij->SetImpactParameterRange(0., 5.);
	// trigger
	genHij->SetTrigger(1);
	genHij->SetPtJet(100.);
	genHij->SetSimpleJets(1);
	genHij->SetRadiation(isw);
	genHij->SetJetEtaRange(-0.3,0.3);
	genHij->SetJetPhiRange(15.,105.);   
        fGenerator = (AliGenerator*)genHij;
	break;

    case kHijing_jj125:
	//coment= comment.Append("HIJING Jet 125 GeV");
	genHij = HijingStandard();
// impact parameter range
	genHij->SetImpactParameterRange(0., 5.);
	// trigger
	genHij->SetTrigger(1);
	genHij->SetPtJet(125.);
	genHij->SetSimpleJets(1);
	genHij->SetRadiation(isw);
	genHij->SetJetEtaRange(-0.3,0.3);
	genHij->SetJetPhiRange(15.,105.);   
        fGenerator = (AliGenerator*)genHij;
	break;
//
// Gamma-Jet
//
    case kHijing_gj25:
	//coment= comment.Append("HIJING Gamma 25 GeV");
	genHij = HijingStandard();
// impact parameter range
	genHij->SetImpactParameterRange(0., 5.);
	// trigger
	genHij->SetTrigger(2);
	genHij->SetPtJet(25.);
	genHij->SetSimpleJets(1);
	genHij->SetRadiation(isw);
	genHij->SetJetEtaRange(-0.3,0.3);
	genHij->SetJetPhiRange(15.,105.);   
        fGenerator = (AliGenerator*)genHij;
	break;

    case kHijing_gj50:
	//coment= comment.Append("HIJING Gamma 50 GeV");
	genHij = HijingStandard();
// impact parameter range
	genHij->SetImpactParameterRange(0., 5.);
	// trigger
	genHij->SetTrigger(2);
	genHij->SetPtJet(50.);
	genHij->SetSimpleJets(1);
	genHij->SetRadiation(isw);
	genHij->SetJetEtaRange(-0.3,0.3);
	genHij->SetJetPhiRange(15.,105.);   
        fGenerator = (AliGenerator*)genHij;
	break;

    case kHijing_gj75:
	//coment= comment.Append("HIJING Gamma 75 GeV");
	genHij = HijingStandard();
// impact parameter range
	genHij->SetImpactParameterRange(0., 5.);
	// trigger
	genHij->SetTrigger(2);
	genHij->SetPtJet(75.);
	genHij->SetSimpleJets(1);
	genHij->SetRadiation(isw);
	genHij->SetJetEtaRange(-0.3,0.3);
	genHij->SetJetPhiRange(15.,105.);   
        fGenerator = (AliGenerator*)genHij;
	break;

    case kHijing_gj100:
	//coment= comment.Append("HIJING Gamma 100 GeV");
	genHij = HijingStandard();
// impact parameter range
	genHij->SetImpactParameterRange(0., 5.);
	// trigger
	genHij->SetTrigger(2);
	genHij->SetPtJet(100.);
	genHij->SetSimpleJets(1);
	genHij->SetRadiation(isw);
	genHij->SetJetEtaRange(-0.3,0.3);
	genHij->SetJetPhiRange(15.,105.);   
        fGenerator = (AliGenerator*)genHij;
	break;

    case kHijing_gj125:
	//coment= comment.Append("HIJING Gamma 125 GeV");
	genHij = HijingStandard();
// impact parameter range
	genHij->SetImpactParameterRange(0., 5.);
	// trigger
	genHij->SetTrigger(2);
	genHij->SetPtJet(125.);
	genHij->SetSimpleJets(1);
	genHij->SetRadiation(isw);
	genHij->SetJetEtaRange(-0.3,0.3);
	genHij->SetJetPhiRange(15.,105.);   
        fGenerator = (AliGenerator*)genHij;
	break;
    case kJetPlusBg:
	genCoct = new AliGenCocktail();
	genCoct->SetMomentumRange(0, 999999.);
	genCoct->SetPhiRange(-180., 180.);
	// Set pseudorapidity range from -8 to 8.
	thmin = EtaToTheta( 5.);   // theta min. <---> eta max
	thmax = EtaToTheta(-5.);  // theta max. <---> eta min 
	genCoct->SetThetaRange(thmin,thmax);

//
//      Underlying Event
//
//	AliGenHIJINGparaBa *bg = new AliGenHIJINGparaBa(82534);
	bg = new AliGenHIJINGparaBa(10);
        fBgGenerator = (AliGenerator*)bg;
//
//      Jets from Pythia
//
	jets = new AliGenPythia(-1);
        fSignalGenerator = (AliGenerator*)jets; 
//   Centre of mass energy 
	jets->SetEnergyCMS(5500.);
//   Process type
	jets->SetProcess(kPyJets);
//   final state kinematic cuts
	jets->SetJetEtaRange(-0.3, 0.3);
	jets->SetJetPhiRange(15., 105.);
//   Structure function
	jets->SetStrucFunc(kGRVLO98);
//   
//   Pt transfer of the hard scattering
	jets->SetPtHard(100.,100.1);
//   Decay type (semielectronic, semimuonic, nodecay)
	jets->SetForceDecay(kAll);
//
//      Add all to cockail ...    
//
	genCoct->AddGenerator(jets,"Jets",1);
	genCoct->AddGenerator(bg,"Underlying Event", 1);
        fGenerator = (AliGenerator*)genCoct;

     break;	
    case kGammaPlusBg:
	genCoct = new AliGenCocktail();
	genCoct->SetMomentumRange(0, 999999.);
	genCoct->SetPhiRange(-180., 180.);
	// Set pseudorapidity range from -8 to 8.
	thmin = EtaToTheta( 5.);   // theta min. <---> eta max
	thmax = EtaToTheta(-5.);  // theta max. <---> eta min 
	genCoct->SetThetaRange(thmin,thmax);
//
//      Underlying Event
//
	bg = new AliGenHIJINGparaBa(82534);
        fBgGenerator = (AliGenerator*)bg;
//
//      Jets from Pythia
//
	jets = new AliGenPythia(-1);
        fSignalGenerator = (AliGenerator*)jets; 
//   Centre of mass energy 
	jets->SetEnergyCMS(5500.);
//   Process type
	jets->SetProcess(kPyDirectGamma);
//   final state kinematic cuts
	jets->SetJetEtaRange(-0.3, 0.3);
	jets->SetJetPhiRange(15., 105.);
	jets->SetGammaEtaRange(-0.12, 0.12);
	jets->SetGammaPhiRange(220., 320.);
//   Structure function
	jets->SetStrucFunc(kGRVLO98);
//   
//   Pt transfer of the hard scattering
	jets->SetPtHard(100.,100.1);
//   Decay type (semielectronic, semimuonic, nodecay)
	jets->SetForceDecay(kAll);
//
//      Add all to cockail ...    
//
	genCoct->AddGenerator(jets,"Jets",1);
	genCoct->AddGenerator(bg,"Underlying Event", 1);
        fGenerator = (AliGenerator*)genCoct;

	break;
    case kJets_50:
//  50 GeV Jets
        genPy = PythiaJets(50.);
        fGenerator = (AliGenerator*)genPy;
        break; 
    case kJets_75:
//  75 GeV Jets
        genPy = PythiaJets(75.);
        fGenerator = (AliGenerator*)genPy;
        break; 
    case kJets_100:
//  100 GeV Jets  
	genPy = PythiaJets(100.);
        fGenerator = (AliGenerator*)genPy;
	break; 
    case kJets_200:
//  200 GeV Jets
        genPy = PythiaJets(200.);
        fGenerator = (AliGenerator*)genPy;
        break;

    case kJets_100RadOn:
//  100 GeV Jets with radiation on - 22-mar-2002 
//  See AliPythia.cxx for default
	genPy = PythiaJets(100.);
	//        genPy->SetKeyPartonJets(1);   // for jet partons
	//        genPy->DefineParametersForPartonsJets(); 

        fGenerator = (AliGenerator*)genPy;
	break; 

    case kGammaJets_50:
//  50 GeV Jets + Gamma
        genPy = PythiaJets(-1);
        genPy->SetEnergyCMS(5500.);
        genPy->SetProcess(kPyDirectGamma);
        genPy->SetJetEtaRange(-0.3,+0.3);
        genPy->SetJetPhiRange(15.,105.);
        genPy->SetGammaEtaRange(-0.12, 0.12);
        genPy->SetGammaPhiRange(220., 320.);
        genPy->SetStrucFunc(kGRVLO98);
        genPy->SetPtHard(50.,50.001);
        genPy->SetForceDecay(kAll);
        fGenerator = (AliGenerator*)genPy;
        break;
    case kGammaJets_75:
//  75 GeV Jets + Gamma 
        genPy = PythiaJets(-1);
        genPy->SetEnergyCMS(5500.);
        genPy->SetProcess(kPyDirectGamma);
        genPy->SetJetEtaRange(-0.3,+0.3);
        genPy->SetJetPhiRange(15.,105.);
        genPy->SetGammaEtaRange(-0.12, 0.12);
        genPy->SetGammaPhiRange(220., 320.);
        genPy->SetStrucFunc(kGRVLO98);
        genPy->SetPtHard(75.,75.001);
        genPy->SetForceDecay(kAll);
        fGenerator = (AliGenerator*)genPy;
        break; 
    case kGammaJets_100:
// 100 GeV Jets + Gamma
        genPy = PythiaJets(-1);
        genPy->SetEnergyCMS(5500.);
        genPy->SetProcess(kPyDirectGamma);
        genPy->SetJetEtaRange(-0.3,+0.3);
        genPy->SetJetPhiRange(15.,105.);
        genPy->SetGammaEtaRange(-0.12, 0.12);
        genPy->SetGammaPhiRange(220., 320.);
        genPy->SetStrucFunc(kGRVLO98);
        genPy->SetPtHard(100.,100.001);
        genPy->SetForceDecay(kAll);
        fGenerator = (AliGenerator*)genPy;
        break; 
    case kGammaJets_200:
//  200 GeV Jets + Gamma
        genPy = PythiaJets(-1);
        genPy->SetEnergyCMS(5500.);
        genPy->SetProcess(kPyDirectGamma);
        genPy->SetJetEtaRange(-0.3,+0.3);
        genPy->SetJetPhiRange(15.,105.);
        genPy->SetGammaEtaRange(-0.12, 0.12);
        genPy->SetGammaPhiRange(220., 320.);
        genPy->SetStrucFunc(kGRVLO98);
        genPy->SetPtHard(200.,200.001);
        genPy->SetForceDecay(kAll);
        fGenerator = (AliGenPythia*)genPy;
        break; 
    case kGammaJets_250:
//  250 GeV Jets + Gamma
        genPy = PythiaJets(-1);
        genPy->SetEnergyCMS(5500.);
        genPy->SetProcess(kPyDirectGamma);
        genPy->SetJetEtaRange(-0.3,+0.3);
        genPy->SetJetPhiRange(15.,105.);
        genPy->SetGammaEtaRange(-0.12, 0.12);
        genPy->SetGammaPhiRange(220., 320.);
        genPy->SetStrucFunc(kGRVLO98);
        genPy->SetPtHard(250.,250.001);
        genPy->SetForceDecay(kAll);
        fGenerator = (AliGenerator*)genPy;
        break; 
    case kGammaJets_300:
//  300 GeV Jets + Gamma
        genPy = PythiaJets(-1);
        genPy->SetEnergyCMS(5500.);
        genPy->SetProcess(kPyDirectGamma);
        genPy->SetJetEtaRange(-0.3,+0.3);
        genPy->SetJetPhiRange(15.,105.);
        genPy->SetGammaEtaRange(-0.12, 0.12);
        genPy->SetGammaPhiRange(220., 320.);
        genPy->SetStrucFunc(kGRVLO98);
        genPy->SetPtHard(300.,300.001);
        genPy->SetForceDecay(kAll);
        fGenerator = (AliGenerator*)genPy;
        break;
    default:
        printf("<I> wrong parameter for generator run %i rad %i\n", run, rad);
        assert(0);
    }
    if(fGenerator) fGenerator->SetPtRange(0.,1.e10); // discard the limit on pT

}

AliEMCALGeneratorFactory::AliEMCALGeneratorFactory(PprRunFact_t run, Float_t p) 
{
    fGenerator = fBgGenerator = fSignalGenerator = 0;
    fRunType   = run;
    fMomentum  = p;

    AliGenBox   *genB=0;

    switch (run) {
    case kGammaBoxOne:
        genB = OneParticleWithFixedEnergy(kGamma, p);  
    break;
    case kPi0BoxOne:
        genB = OneParticleWithFixedEnergy(kPi0, p);  
    break;
    default:
        printf("<I> wrong parameter for generator run %i \n",run);
        assert(0);
    }
    if(genB) fGenerator = (AliGenerator*)genB;
    //    if(fGenerator) fGenerator->SetPtRange(0.,1.e10); // discard the limit on pT - 23-aug-04
}

AliGenHijing* AliEMCALGeneratorFactory::HijingStandard()
{
    AliGenHijing *gener = new AliGenHijing(-1);
// centre of mass energy 
    gener->SetEnergyCMS(5500.);
// reference frame
    gener->SetReferenceFrame("CMS");
// projectile
    gener->SetProjectile("A", 208, 82);
    gener->SetTarget    ("A", 208, 82);
// tell hijing to keep the full parent child chain
    gener->KeepFullEvent();
// enable jet quenching
    gener->SetJetQuenching(1);
// enable shadowing
    gener->SetShadowing(1);
// neutral pion and heavy particle decays switched off
    gener->SetDecaysOff(1);
// Don't track spectators
    gener->SetSpectators(0);
// kinematic selection
    gener->SetSelectAll(0);
    return gener;
}

AliGenPythia* AliEMCALGeneratorFactory::PythiaJets(Float_t energy)
{
    AliGenPythia *gener = new AliGenPythia(-1);
//   Centre of mass energy 
    gener->SetEnergyCMS(5500.);
//   Process type
    gener->SetProcess(kPyJets);
//   final state kinematic cuts
    gener->SetJetEtaRange(-0.3, 0.3);
    gener->SetJetPhiRange(15., 105.);
//   Structure function
    gener->SetStrucFunc(kGRVLO98);
//   
//   Pt transfer of the hard scattering
    gener->SetPtHard(energy, energy+0.1);
//   Decay type (semielectronic, semimuonic, nodecay)
    gener->SetForceDecay(kAll);
//
    return gener;
}


AliGenPythia* AliEMCALGeneratorFactory::PythiaGamma(Float_t energy)
{
    AliGenPythia *gener = new AliGenPythia(-1);
//   Centre of mass energy 
    gener->SetEnergyCMS(5500.);
//   Process type
    gener->SetProcess(kPyDirectGamma);
//   final state kinematic cuts
    gener->SetJetEtaRange(-0.3, 0.3);
    gener->SetJetPhiRange(15., 105.);
    gener->SetGammaEtaRange(-0.12, 0.12);
    gener->SetGammaPhiRange(220., 320.);
//   Structure function
    gener->SetStrucFunc(kGRVLO98);
//   
//   Pt transfer of the hard scattering
    gener->SetPtHard(energy, energy+0.1);
//   Decay type (semielectronic, semimuonic, nodecay)
    gener->SetForceDecay(kAll);
//
    return gener;
}

//
// Staff of Aleksei Pavlinov.
//
AliGenBox* AliEMCALGeneratorFactory::OneParticleWithFixedEnergy(Int_t type, Float_t p)
{// one particle in EMCAL acceptance
   Float_t thmin  = EtaToTheta(0.7), thmax  = EtaToTheta(-0.7); 
   Float_t phimin=0, phimax=120;
   Float_t pmin=p, pmax=p+0.01;

   AliGenBox *gen = new AliGenBox(1);
   gen->SetPart(type);
   gen->SetNumberParticles(1); 
   gen->SetThetaRange(thmin, thmax);
   gen->SetPhiRange(phimin, phimax);
   gen->SetMomentumRange(pmin, pmax);

   printf("<I> AliEMCALGeneratorFactory::OneParticleWithFixedEnergy \n");
   printf("       type of particle -> %i \n", type);
   printf("    %6.4f <    eta   < %6.4f\n", thmin,  thmax);
   printf("    %6.4f <    phi   < %6.4f\n", phimin, phimax);
   printf("    %7.2f < Momentum < %7.2f\n", pmin, pmax);
   printf("     TestBit(kPtRange) %i | TestBit(kMomentumRange) %i\n", 
   gen->TestBit(BIT(17)), gen->TestBit(BIT(19)));
   return gen;
}
