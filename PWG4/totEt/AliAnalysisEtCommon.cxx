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
#include "TF1.h"
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

AliAnalysisEtCommon::AliAnalysisEtCommon() : TObject()
					   ,fHistogramNameSuffix("")
					   ,fCuts(0)
					   ,fDataSet(2010)
					   ,fV0ScaleDataSet(2010)
					   ,fEsdtrackCutsITSTPC(0)
					   ,fEsdtrackCutsTPC(0)
					   ,fEsdtrackCutsITS(0)
					   ,fK0PythiaD6T(0)
					   ,fLambdaPythiaD6T(0)
					   ,fAntiLambdaPythiaD6T(0)
					   ,fK0Data(0)
					   ,fLambdaData(0)
					   ,fAntiLambdaData(0)
					   ,fLambdaEnhancement(0)
					   ,fProtonEnhancement(0)
					   ,fCentralityMethod("V0M")
					   ,fNCentBins(21)
					   ,fCentBin(-1)
{//default constructor

}

AliAnalysisEtCommon::~AliAnalysisEtCommon()
{//destructor
  delete fCuts;
  delete fEsdtrackCutsITSTPC;
  delete fEsdtrackCutsITS;
  delete fEsdtrackCutsTPC;
  delete fK0PythiaD6T;
  delete fLambdaPythiaD6T;
  delete fAntiLambdaPythiaD6T;
  delete fK0Data;
  delete fLambdaData;
  delete fAntiLambdaData;
  delete fLambdaEnhancement;
  delete fProtonEnhancement;
}

Int_t AliAnalysisEtCommon::AnalyseEvent(AliVEvent */*event*/)
{ //this line is basically here to eliminate a compiler warning that event is not used.  Making it a virtual function did not work with the plugin.
//  cout << "This event has " << event->GetNumberOfTracks() << " tracks" << endl;
  ResetEventValues();
  return 0;
}


void AliAnalysisEtCommon::Init()
{// clear variables, set up cuts and PDG info
  // LevyPt function described in LevyFitEvaluate below
  //parameter 0 = dNdy
  //parameter 1 = temp
  //parameter 2 = power
  if(fK0PythiaD6T) delete fK0PythiaD6T;
  if(fLambdaPythiaD6T) delete fLambdaPythiaD6T;
  if(fAntiLambdaPythiaD6T) delete fAntiLambdaPythiaD6T;
  if(fK0Data) delete fK0Data;
  if(fLambdaData) delete fLambdaData;
  if(fAntiLambdaData) delete fAntiLambdaData;

  fK0PythiaD6T = new TF1("K0PythiaD6T", LevyPtEvaluate, 0,50,4);
  fLambdaPythiaD6T = new TF1("LambdaPythiaD6T", LevyPtEvaluate,0,50,4);
  fAntiLambdaPythiaD6T = new TF1("AntiLambdaPythiaD6T", LevyPtEvaluate,0,50,4);
  fK0Data = new TF1("K0Data", LevyPtEvaluate,0,50,4);
  fLambdaData = new TF1("LambdaData", LevyPtEvaluate,0,50,4);
  fAntiLambdaData = new TF1("AntiLambdaData", LevyPtEvaluate,0,50,4);

  fK0PythiaD6T->FixParameter(3,0.493677);
  fK0Data->FixParameter(3,0.493677);
  fLambdaPythiaD6T->FixParameter(3,1.115683);
  fAntiLambdaPythiaD6T->FixParameter(3,1.115683);
  fLambdaData->FixParameter(3,1.115683);
  fAntiLambdaData->FixParameter(3,1.115683);
  if(fV0ScaleDataSet==2009){
    //These data are from the ALICE 900 GeV p+p paper
    fK0PythiaD6T->SetParameter(0,0.1437);
    fK0PythiaD6T->SetParameter(1,0.1497);
    fK0PythiaD6T->SetParameter(2,6.94);
    fLambdaPythiaD6T->SetParameter(0,0.0213);
    fLambdaPythiaD6T->SetParameter(1,0.1315);
    fLambdaPythiaD6T->SetParameter(2,4.60);
    fAntiLambdaPythiaD6T->SetParameter(0,0.0213);
    fAntiLambdaPythiaD6T->SetParameter(1,0.1315);
    fAntiLambdaPythiaD6T->SetParameter(2,4.60);
    fK0Data->SetParameter(0,0.184);
    fK0Data->SetParameter(1,0.168);
    fK0Data->SetParameter(2,6.6);
    fLambdaData->SetParameter(0,0.048);
    fLambdaData->SetParameter(1,0.229);
    fLambdaData->SetParameter(2,10.8);
    fAntiLambdaData->SetParameter(0,0.047);
    fAntiLambdaData->SetParameter(1,0.210);
    fAntiLambdaData->SetParameter(2,9.2);
  }
  if(fV0ScaleDataSet==2010 ||fV0ScaleDataSet==20100 ){
    //These data are from the CMS analysis note on 7 TeV spectra
    //http://cdsweb.cern.ch/record/1279344/files/QCD-10-007-pas.pdf
    //Note the CMS parameterization of the Levy function differs from the ALICE parameterization by a constant.
    //CMS does not list the overall constant in their fit, the ratios of the dN/dy(y=0) is used.
    fK0PythiaD6T->SetParameter(0,0.72);
    fK0PythiaD6T->SetParameter(1,0.183);
    fK0PythiaD6T->SetParameter(2,7.41);
    fLambdaPythiaD6T->SetParameter(0,0.54);
    fLambdaPythiaD6T->SetParameter(1,0.216);
    fLambdaPythiaD6T->SetParameter(2,5.71);
    fAntiLambdaPythiaD6T->SetParameter(0,0.54);
    fAntiLambdaPythiaD6T->SetParameter(1,0.216);
    fAntiLambdaPythiaD6T->SetParameter(2,5.71);
    fK0Data->SetParameter(0,1.0);
    fK0Data->SetParameter(1,0.215);
    fK0Data->SetParameter(2,6.79);
    fLambdaData->SetParameter(0,1.0);
    fLambdaData->SetParameter(1,0.290);
    fLambdaData->SetParameter(2,9.28);
    fAntiLambdaData->SetParameter(0,1.0);
    fAntiLambdaData->SetParameter(1,0.290);
    fAntiLambdaData->SetParameter(2,9.28);
  }
  if(fLambdaEnhancement) delete fLambdaEnhancement;
  fLambdaEnhancement = new TF1("fLambdaEnhancement","([0]*pow(x,[1])*exp(-pow(x/[2],[3])))/([4]*exp(-pow([5]/x,[6]))+[7]*x)",0,50);
   fLambdaEnhancement->SetParameter(0,0.5630487);
   fLambdaEnhancement->SetParameter(1,1.388818);
   fLambdaEnhancement->SetParameter(2,3.954147);
   fLambdaEnhancement->SetParameter(3,3.443772);
   fLambdaEnhancement->SetParameter(4,2.844288);
   fLambdaEnhancement->SetParameter(5,2);
   fLambdaEnhancement->SetParameter(6,0.4747893);
   fLambdaEnhancement->SetParameter(7,-0.2250856);
   if(fProtonEnhancement) delete fProtonEnhancement;
   fProtonEnhancement = new TF1("fProtonEnhancement","[0]*pow(x,[1])*exp(-pow(x/[2],[3]))/([4]+[5]*x)",0,50);
   fProtonEnhancement->SetParameter(0,0.5630487*1.6);
   fProtonEnhancement->SetParameter(1,1.388818);
   fProtonEnhancement->SetParameter(2,3.954147/1.5);
   fProtonEnhancement->SetParameter(3,3.443772/2.5);
   fProtonEnhancement->SetParameter(4,0.5);
   fProtonEnhancement->SetParameter(5,-.03);
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

Double_t AliAnalysisEtCommon::LevyPtEvaluate(const Double_t *pt, 
					     const Double_t *par) 
{//LevyPt function for TF1's
  
  Double_t lMass  = par[3];
  Double_t ldNdy  = par[0];
  Double_t l2pi   = 2*TMath::Pi();
  Double_t lTemp = par[1];
  Double_t lPower = par[2];
  
  Double_t lBigCoef = ((lPower-1)*(lPower-2)) / (l2pi*lPower*lTemp*(lPower*lTemp+lMass*(lPower-2)));
  Double_t lInPower = 1 + (TMath::Sqrt(pt[0]*pt[0]+lMass*lMass)-lMass) / (lPower*lTemp);
  
  return ldNdy * pt[0] * lBigCoef * TMath::Power(lInPower,(-1)*lPower);
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

Float_t AliAnalysisEtCommon::K0Weight(Float_t pt){
  Float_t data = fK0Data->Eval(pt);
  Float_t mc = fK0PythiaD6T->Eval(pt);
  return data/mc;
}
Float_t AliAnalysisEtCommon::LambdaWeight(Float_t pt){
  Float_t data = fLambdaData->Eval(pt);
  Float_t mc = fLambdaPythiaD6T->Eval(pt);
  return data/mc;
}
Float_t AliAnalysisEtCommon::AntiLambdaWeight(Float_t pt){
  Float_t data = fAntiLambdaData->Eval(pt);
  Float_t mc = fAntiLambdaPythiaD6T->Eval(pt);
  return data/mc;
}

Float_t AliAnalysisEtCommon::LambdaBaryonEnhancement(Float_t pt){
  if(pt<0.8) return 1.0;
  return fLambdaEnhancement->Eval(pt);
};//Function which gives the factor to reweigh a lambda or antilambda so it roughly matches baryon enhancement seen at RHIC
Float_t AliAnalysisEtCommon::ProtonBaryonEnhancement(Float_t pt){
  if(pt<0.8) return 1.0;
  return fProtonEnhancement->Eval(pt);
}//Function which gives the factor to reweigh a lambda or antilambda so it roughly matches baryon enhancement seen at RHIC
