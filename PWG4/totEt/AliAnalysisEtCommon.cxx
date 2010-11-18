//Created by Christine Nattrass, Rebecca Scott, Irakli Martashvili
//University of Tennessee at Knoxville
//
// This class is designed for the analysis of the hadronic component of 
// transverse energy.  It is used by AliAnalysisTaskHadEt.
// This gets information about the hadronic component of the transverse energy 
// from tracks reconstructed in an event
// it has daughters, AliAnalysisEtCommonMonteCarlo and 
// AliAnalysisEtCommonReconstructed which loop over either Monte Carlo data or 
// real data to get Et

#include "AliAnalysisEtCommon.h"
#include "TMath.h"
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include <iostream>
#include "AliAnalysisEtCuts.h"
#include "AliMCEvent.h"
#include "AliVEvent.h"
#include "AliStack.h"
#include "AliESDtrackCuts.h"
#include "TDatabasePDG.h"
#include "TParticle.h"
#include "Rtypes.h"
#include "AliPDG.h"

using namespace std;

ClassImp(AliAnalysisEtCommon);
//These are from the PDG database but by making them static the code is a bit more efficient and has no problems running with the plugin
//Cuts are defined in $ROOTSYS/etc/pdg_table.txt

Float_t AliAnalysisEtCommon::fgPionMass = 0.13957;
Float_t AliAnalysisEtCommon::fgKaonMass = 0.493677;
Float_t AliAnalysisEtCommon::fgProtonMass = 0.938272;
Float_t AliAnalysisEtCommon::fgElectronMass = 0.000510999;
Int_t AliAnalysisEtCommon::fgPiPlusCode = 211;
Int_t AliAnalysisEtCommon::fgPiMinusCode = -211;
Int_t AliAnalysisEtCommon::fgKPlusCode = 321;
Int_t AliAnalysisEtCommon::fgKMinusCode = -321;
Int_t AliAnalysisEtCommon::fgProtonCode = 2212;
Int_t AliAnalysisEtCommon::fgAntiProtonCode = -2212;
Int_t AliAnalysisEtCommon::fgLambdaCode = 3122;
Int_t AliAnalysisEtCommon::fgAntiLambdaCode = -3122;
Int_t AliAnalysisEtCommon::fgK0SCode = 310;
Int_t AliAnalysisEtCommon::fgOmegaCode = 3334;
Int_t AliAnalysisEtCommon::fgAntiOmegaCode = -3334;
Int_t AliAnalysisEtCommon::fgXi0Code = 3322;
Int_t AliAnalysisEtCommon::fgAntiXi0Code = -3322;
Int_t AliAnalysisEtCommon::fgXiCode = 3312;
Int_t AliAnalysisEtCommon::fgAntiXiCode = -3312;
Int_t AliAnalysisEtCommon::fgSigmaCode = 3112;
Int_t AliAnalysisEtCommon::fgAntiSigmaCode = -3112;
Int_t AliAnalysisEtCommon::fgK0LCode = 130;
Int_t AliAnalysisEtCommon::fgNeutronCode = 2112;
Int_t AliAnalysisEtCommon::fgAntiNeutronCode = -2112;
Int_t AliAnalysisEtCommon::fgEPlusCode = -11;
Int_t AliAnalysisEtCommon::fgEMinusCode = 11;
Int_t AliAnalysisEtCommon::fgMuPlusCode = -13;
Int_t AliAnalysisEtCommon::fgMuMinusCode = 13;
Int_t AliAnalysisEtCommon::fgGammaCode = 22;
Int_t AliAnalysisEtCommon::fgPi0Code = 111;
Int_t AliAnalysisEtCommon::fgEtaCode = 221;
Int_t AliAnalysisEtCommon::fgOmega0Code = 223;



Float_t AliAnalysisEtCommon::fgPtTPCCutOff = 0.15;
Float_t AliAnalysisEtCommon::fgPtITSCutOff = 0.10;

AliAnalysisEtCommon::AliAnalysisEtCommon() :
  fHistogramNameSuffix("")
  ,fCuts(0)
  ,fEsdtrackCutsITSTPC(0)
  ,fEsdtrackCutsTPC(0)
  ,fEsdtrackCutsITS(0)
{//default constructor

}

AliAnalysisEtCommon::~AliAnalysisEtCommon()
{//destructor
  delete fCuts;
  delete fEsdtrackCutsITSTPC;
  delete fEsdtrackCutsITS;
  delete fEsdtrackCutsTPC;
}

Int_t AliAnalysisEtCommon::AnalyseEvent(AliVEvent *event)
{ //this line is basically here to eliminate a compiler warning that event is not used.  Making it a virtual function did not work with the plugin.
  cout << "This event has " << event->GetNumberOfTracks() << " tracks" << endl;
  ResetEventValues();
  return 0;
}


void AliAnalysisEtCommon::Init()
{// clear variables, set up cuts and PDG info

}

void AliAnalysisEtCommon::ResetEventValues()
{//Resets event values of et to zero
  
  if (!fCuts) { // some Init's needed
    cout << __FILE__ << ":" << __LINE__ << " : Init " << endl;
    if (!fCuts) {
      cout << " setting up Cuts " << endl;
      fCuts = new AliAnalysisEtCuts();
    }
  }
}


Float_t AliAnalysisEtCommon::Et(TParticle *part, float mass){//function to calculate et in the same way as it would be calculated in a calorimeter
  if(mass+1000<0.01){//if no mass given return default.  The default argument is -1000
    if(TMath::Abs(part->GetPDG(0)->PdgCode())==2212 || TMath::Abs(part->GetPDG(0)->PdgCode())==2112){
      if(part->GetPDG(0)->PdgCode()==-2212 || part->GetPDG(0)->PdgCode()==-2112){//antiproton or antineutron
	//for antinucleons we specifically want to return the kinetic energy plus twice the rest mass
	return (part->Energy()+part->GetMass())*TMath::Sin(part->Theta());
      }
      if(part->GetPDG(0)->PdgCode()==2212 || part->GetPDG(0)->PdgCode()==2112){//proton or neutron
	//for nucleons we specifically want to return the kinetic energy only
	return (part->Energy()-part->GetMass())*TMath::Sin(part->Theta());
      }
    }
    else{//otherwise go to the default
      return part->Energy()*TMath::Sin(part->Theta());
    }
  }
  else{//otherwise use the mass that was given
    return (TMath::Sqrt(TMath::Power(part->P(),2.0)+TMath::Power(mass,2.0)))*TMath::Sin(part->Theta());
  }
  return 0.0;
}
Float_t AliAnalysisEtCommon::Et(Float_t p, Float_t theta, Int_t pid, Short_t charge) const {//function to calculate et in the same way as it would be calculated in a calorimeter
  if(pid==fgPiPlusCode || pid==fgPiMinusCode){//Nothing special for pions
    return TMath::Sqrt(p*p + fgPionMass*fgPionMass) * TMath::Sin(theta);
  }
  if(pid==fgKPlusCode || pid==fgKMinusCode){//Nothing special for kaons
    return TMath::Sqrt(p*p + fgKaonMass*fgKaonMass) * TMath::Sin(theta);
  }
  if(pid==fgEPlusCode || pid==fgEMinusCode){//Nothing special for electrons
    return TMath::Sqrt(p*p + fgElectronMass*fgElectronMass) * TMath::Sin(theta);
  }
  if(pid==fgProtonCode || pid==fgAntiProtonCode){//But for protons we must be careful...
    if(charge<0.0){//antiprotns: kinetic energy plus twice the rest mass
      return (TMath::Sqrt(p*p + fgProtonMass*fgProtonMass) + fgProtonMass) * TMath::Sin(theta);
    }
    if(charge>0.0){//antiprotns: kinetic energy only
      return (TMath::Sqrt(p*p + fgProtonMass*fgProtonMass) - fgProtonMass) * TMath::Sin(theta);
    }
  }
  cerr<<"Uh-oh!  Et not set properly!"<<endl;
  return 0.0;
}
