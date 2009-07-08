 /**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes hereby granted      *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
/* $Id: $ */

//_________________________________________________________________________
//
// Class for the electron identification.
// Clusters from EMCAL matched to tracks
// and kept in the AOD. Few histograms produced.
//
// -- Author: J.L. Klay (Cal Poly)
//////////////////////////////////////////////////////////////////////////////
  
  
// --- ROOT system --- 
#include <TH2F.h>
#include <TClonesArray.h>
#include <TObjString.h>
#include <Riostream.h>

// --- Analysis system --- 
#include "AliAnaElectron.h" 
#include "AliCaloTrackReader.h"
#include "AliMCAnalysisUtils.h"
#include "AliFidutialCut.h"
#include "AliESDCaloCluster.h"
//#include "AliESDCaloCells.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliCaloPID.h"
#include "AliVEvent.h"

ClassImp(AliAnaElectron)
  
//____________________________________________________________________________
AliAnaElectron::AliAnaElectron() 
 : AliAnaPartCorrBaseClass(), fCalorimeter(""),
  fpOverEmin(0.),fpOverEmax(0.),fResidualCut(0.),
  //matching checks
  fh1pOverE(0),fh1dR(0),fh2EledEdx(0),fh2MatchdEdx(0),fh2dEtadPhi(0),
  fh2dEtadPhiMatched(0),fh2dEtadPhiUnmatched(0),
  fh2OuterPtVsExtrapPt(0),fh2OuterPhiVsExtrapPhi(0),fh2OuterEtaVsExtrapEta(0),
  fh2TrackPVsClusterE(0),fh2TrackPtVsClusterE(0),fh2TrackPhiVsClusterPhi(0),fh2TrackEtaVsClusterEta(0),
  //reco
  fhPtElectron(0),fhPhiElectron(0),fhEtaElectron(0),
  //MC
  fhPtConversion(0),fhPhiConversion(0),fhEtaConversion(0),
  fhPtBottom(0),fhPhiBottom(0),fhEtaBottom(0),
  fhPtCharm(0),fhPhiCharm(0),fhEtaCharm(0),
  fhPtCFromB(0),fhPhiCFromB(0),fhEtaCFromB(0),
  fhPtDalitz(0),fhPhiDalitz(0),fhEtaDalitz(0),
  fhPtWDecay(0),fhPhiWDecay(0),fhEtaWDecay(0),
  fhPtZDecay(0),fhPhiZDecay(0),fhEtaZDecay(0),
  fhPtPrompt(0),fhPhiPrompt(0),fhEtaPrompt(0),
  fhPtUnknown(0),fhPhiUnknown(0),fhEtaUnknown(0)
//  fhMCElePt(0),fhMCElePhi(0),fhMCEleEta(0)
{
  //default ctor
  
  //Initialize parameters
  InitParameters();

}

//____________________________________________________________________________
AliAnaElectron::AliAnaElectron(const AliAnaElectron & g) 
 : AliAnaPartCorrBaseClass(g), fCalorimeter(g.fCalorimeter),
   fpOverEmin(g.fpOverEmin),fpOverEmax(g.fpOverEmax),fResidualCut(g.fResidualCut),
   //matching checks
   fh1pOverE(g.fh1pOverE),fh1dR(g.fh1dR),
   fh2EledEdx(g.fh2EledEdx),fh2MatchdEdx(g.fh2MatchdEdx),fh2dEtadPhi(g.fh2dEtadPhi),
   fh2dEtadPhiMatched(g.fh2dEtadPhiMatched),fh2dEtadPhiUnmatched(g.fh2dEtadPhiUnmatched),
   fh2OuterPtVsExtrapPt(g.fh2OuterPtVsExtrapPt),fh2OuterPhiVsExtrapPhi(g.fh2OuterPhiVsExtrapPhi),
   fh2OuterEtaVsExtrapEta(g.fh2OuterEtaVsExtrapEta),
   fh2TrackPVsClusterE(g.fh2TrackPVsClusterE),fh2TrackPtVsClusterE(g.fh2TrackPtVsClusterE),
   fh2TrackPhiVsClusterPhi(g.fh2TrackPhiVsClusterPhi),fh2TrackEtaVsClusterEta(g.fh2TrackEtaVsClusterEta),   
   //reco
   fhPtElectron(g.fhPtElectron),fhPhiElectron(g.fhPhiElectron),fhEtaElectron(g.fhEtaElectron),
   //MC
   fhPtConversion(g.fhPtConversion),fhPhiConversion(g.fhPhiConversion),fhEtaConversion(g.fhEtaConversion),
   fhPtBottom(g.fhPtBottom),fhPhiBottom(g.fhPhiBottom),fhEtaBottom(g.fhEtaBottom),
   fhPtCharm(g.fhPtCharm),fhPhiCharm(g.fhPhiCharm),fhEtaCharm(g.fhEtaCharm),
   fhPtCFromB(g.fhPtCFromB),fhPhiCFromB(g.fhPhiCFromB),fhEtaCFromB(g.fhEtaCFromB),
   fhPtDalitz(g.fhPtDalitz),fhPhiDalitz(g.fhPhiDalitz),fhEtaDalitz(g.fhEtaDalitz),
   fhPtWDecay(g.fhPtWDecay),fhPhiWDecay(g.fhPhiWDecay),fhEtaWDecay(g.fhEtaWDecay),
   fhPtZDecay(g.fhPtZDecay),fhPhiZDecay(g.fhPhiZDecay),fhEtaZDecay(g.fhEtaZDecay),
   fhPtPrompt(g.fhPtPrompt),fhPhiPrompt(g.fhPhiPrompt),fhEtaPrompt(g.fhEtaPrompt),
   fhPtUnknown(g.fhPtUnknown),fhPhiUnknown(g.fhPhiUnknown),fhEtaUnknown(g.fhEtaUnknown)
   //   fhMCElePt(g.fhMCElePt),fhMCElePhi(g.fhMCElePhi),fhMCEleEta(g.fhMCEleEta)
{
  // cpy ctor
  
}

//_________________________________________________________________________
AliAnaElectron & AliAnaElectron::operator = (const AliAnaElectron & g)
{
  // assignment operator
  
  if(&g == this) return *this;
  fCalorimeter = g.fCalorimeter;
  fpOverEmin = g.fpOverEmin;
  fpOverEmax = g.fpOverEmax;
  fResidualCut = g.fResidualCut;
  //fhEnergy = g.fhEnergy;
  //fhClusMult = g.fhClusMult;
  //fhClusters = g.fhClusters;
  //fhDigitsEvent = g.fhDigitsEvent;
  fh1pOverE = g.fh1pOverE;
  fh1dR = g.fh1dR;
  fh2EledEdx = g.fh2EledEdx;
  fh2MatchdEdx = g.fh2MatchdEdx;
  fh2dEtadPhi = g.fh2dEtadPhi;
  fh2dEtadPhiMatched = g.fh2dEtadPhiMatched;
  fh2dEtadPhiUnmatched = g.fh2dEtadPhiUnmatched;
  fh2OuterPtVsExtrapPt = g.fh2OuterPtVsExtrapPt;
  fh2OuterPhiVsExtrapPhi = g.fh2OuterPhiVsExtrapPhi;
  fh2OuterEtaVsExtrapEta = g.fh2OuterEtaVsExtrapEta;
  fh2TrackPVsClusterE = g.fh2TrackPVsClusterE;
  fh2TrackPtVsClusterE = g.fh2TrackPtVsClusterE;
  fh2TrackPhiVsClusterPhi = g.fh2TrackPhiVsClusterPhi;
  fh2TrackEtaVsClusterEta = g.fh2TrackEtaVsClusterEta;   
  fhPtElectron = g.fhPtElectron;
  fhPhiElectron = g.fhPhiElectron;
  fhEtaElectron = g.fhEtaElectron;
  fhPtConversion = g.fhPtConversion;
  fhPhiConversion = g.fhPhiConversion;
  fhEtaConversion = g.fhEtaConversion;
  fhPtBottom = g.fhPtBottom;
  fhPhiBottom = g.fhPhiBottom;
  fhEtaBottom = g.fhEtaBottom;
  fhPtCharm = g.fhPtCharm;
  fhPhiCharm = g.fhPhiCharm;
  fhEtaCharm = g.fhEtaCharm;
  fhPtCFromB = g.fhPtCFromB;
  fhPhiCFromB = g.fhPhiCFromB;
  fhEtaCFromB = g.fhEtaCFromB;
  fhPtDalitz = g.fhPtDalitz;
  fhPhiDalitz = g.fhPhiDalitz;
  fhEtaDalitz = g.fhEtaDalitz;
  fhPtWDecay = g.fhPtWDecay;
  fhPhiWDecay = g.fhPhiWDecay;
  fhEtaWDecay = g.fhEtaWDecay;
  fhPtZDecay = g.fhPtZDecay;
  fhPhiZDecay = g.fhPhiZDecay;
  fhEtaZDecay = g.fhEtaZDecay;
  fhPtPrompt = g.fhPtPrompt;
  fhPhiPrompt = g.fhPhiPrompt;
  fhEtaPrompt = g.fhEtaPrompt;
  fhPtUnknown = g.fhPtUnknown;
  fhPhiUnknown = g.fhPhiUnknown;
  fhEtaUnknown = g.fhEtaUnknown;

  /*
  fhMCElePt = g.fhMCElePt;
  fhMCElePhi = g.fhMCElePhi;
  fhMCEleEta = g.fhMCEleEta;
  */

  return *this;
  
}

//____________________________________________________________________________
AliAnaElectron::~AliAnaElectron() 
{
  //dtor

}


//________________________________________________________________________
TList *  AliAnaElectron::GetCreateOutputObjects()
{  
  // Create histograms to be saved in output file and 
  // store them in outputContainer
  TList * outputContainer = new TList() ; 
  outputContainer->SetName("ElectronHistos") ; 
  
  Int_t nptbins  = GetHistoNPtBins();
  Int_t nphibins = GetHistoNPhiBins();
  Int_t netabins = GetHistoNEtaBins();
  Float_t ptmax  = GetHistoPtMax();
  Float_t phimax = GetHistoPhiMax();
  Float_t etamax = GetHistoEtaMax();
  Float_t ptmin  = GetHistoPtMin();
  Float_t phimin = GetHistoPhiMin();
  Float_t etamin = GetHistoEtaMin();	

  fh1pOverE = new TH1F("h1pOverE","EMCAL-TRACK matches p/E",100,0.,10.);
  fh1dR = new TH1F("h1dR","EMCAL-TRACK matches dR",300, 0.,TMath::Pi());
  fh2EledEdx = new TH2F("h2EledEdx","dE/dx vs. p for electrons",200,0.,50.,200,0.,400.);
  fh2MatchdEdx = new TH2F("h2MatchdEdx","dE/dx vs. p for all matches",200,0.,50.,200,0.,400.);
  fh2dEtadPhi = new TH2F("h2dEtadPhi","#Delta#eta vs. #Delta#phi for all track-cluster pairs",200,0.,1.4,300,0.,TMath::Pi());
  fh2dEtadPhiMatched = new TH2F("h2dEtadPhiMatched","#Delta#eta vs. #Delta#phi for matched track-cluster pairs",200,0.,1.4,300,0.,TMath::Pi());
  fh2dEtadPhiUnmatched = new TH2F("h2dEtadPhiUnmatched","#Delta#eta vs. #Delta#phi for unmatched track-cluster pairs",200,0.,1.4,300,0.,TMath::Pi());
  fh2OuterPtVsExtrapPt = new TH2F("h2OuterPtVsExtrapPt","h2OuterPtVsExtrapPt",nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
  fh2OuterPhiVsExtrapPhi = new TH2F("h2OuterPhiVsExtrapPhi","h2OuterPhiVsExtrapPhi",nphibins,phimin,phimax,nphibins,phimin,phimax);
  fh2OuterEtaVsExtrapEta = new TH2F("h2OuterEtaVsExtrapEta","h2OuterEtaVsExtrapEta",netabins,etamin,etamax,netabins,etamin,etamax);

  fh2TrackPVsClusterE = new TH2F("h2TrackPVsClusterE","h2TrackPVsClusterE",nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
  fh2TrackPtVsClusterE = new TH2F("h2TrackPtVsClusterE","h2TrackPtVsClusterE",nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
  fh2TrackPhiVsClusterPhi = new TH2F("h2TrackPhiVsClusterPhi","h2TrackPhiVsClusterPhi",nphibins,phimin,phimax,nphibins,phimin,phimax);
  fh2TrackEtaVsClusterEta = new TH2F("h2TrackEtaVsClusterEta","h2TrackEtaVsClusterEta",netabins,etamin,etamax,netabins,etamin,etamax);

  outputContainer->Add(fh1pOverE) ; 
  outputContainer->Add(fh1dR) ; 
  outputContainer->Add(fh2EledEdx) ;
  outputContainer->Add(fh2MatchdEdx) ;
  outputContainer->Add(fh2dEtadPhi) ;
  outputContainer->Add(fh2dEtadPhiMatched) ;
  outputContainer->Add(fh2dEtadPhiUnmatched) ;
  outputContainer->Add(fh2OuterPtVsExtrapPt) ;
  outputContainer->Add(fh2OuterPhiVsExtrapPhi) ;
  outputContainer->Add(fh2OuterEtaVsExtrapEta) ;
  outputContainer->Add(fh2TrackPVsClusterE) ;
  outputContainer->Add(fh2TrackPtVsClusterE) ;
  outputContainer->Add(fh2TrackPhiVsClusterPhi) ;
  outputContainer->Add(fh2TrackEtaVsClusterEta) ;
  
  fhPtElectron = new TH1F("hPtElectron","Electron pT",nptbins,ptmin,ptmax);
  fhPhiElectron = new TH2F("hPhiElectron","Electron phi vs pT",nptbins,ptmin,ptmax,nphibins,phimin,phimax);
  fhEtaElectron = new TH2F("hEtaElectron","Electron eta vs. eta",nptbins,ptmin,ptmax,netabins,etamin,etamax);

  outputContainer->Add(fhPtElectron) ; 
  outputContainer->Add(fhPhiElectron) ; 
  outputContainer->Add(fhEtaElectron) ;
  
  if(IsDataMC()){
    
    fhPtConversion = new TH1F("hPtConversion","Conversion electron pT",nptbins,ptmin,ptmax);
    fhPhiConversion = new TH2F("hPhiConversion","Conversion Electron phi vs pT",nptbins,ptmin,ptmax,nphibins,phimin,phimax);
    fhEtaConversion = new TH2F("hEtaConversion","Conversion Electron eta vs. eta",nptbins,ptmin,ptmax,netabins,etamin,etamax);
    fhPtBottom = new TH1F("hPtBottom","Bottom electron pT",nptbins,ptmin,ptmax);
    fhPhiBottom = new TH2F("hPhiBottom","Bottom Electron phi vs pT",nptbins,ptmin,ptmax,nphibins,phimin,phimax);
    fhEtaBottom = new TH2F("hEtaBottom","Bottom Electron eta vs. eta",nptbins,ptmin,ptmax,netabins,etamin,etamax);
    fhPtCharm = new TH1F("hPtCharm","Charm electron pT",nptbins,ptmin,ptmax);
    fhPhiCharm = new TH2F("hPhiCharm","Charm Electron phi vs pT",nptbins,ptmin,ptmax,nphibins,phimin,phimax);
    fhEtaCharm = new TH2F("hEtaCharm","Charm Electron eta vs. eta",nptbins,ptmin,ptmax,netabins,etamin,etamax);
    fhPtCFromB = new TH1F("hPtCFromB","Charm from Bottom electron pT",nptbins,ptmin,ptmax);
    fhPhiCFromB = new TH2F("hPhiCFromB","Charm from Bottom Electron phi vs pT",nptbins,ptmin,ptmax,nphibins,phimin,phimax);
    fhEtaCFromB = new TH2F("hEtaCFromB","Charm from Bottom Electron eta vs. eta",nptbins,ptmin,ptmax,netabins,etamin,etamax);
    fhPtDalitz = new TH1F("hPtDalitz","Dalitz electron pT",nptbins,ptmin,ptmax);
    fhPhiDalitz = new TH2F("hPhiDalitz","Dalitz Electron phi vs pT",nptbins,ptmin,ptmax,nphibins,phimin,phimax);
    fhEtaDalitz = new TH2F("hEtaDalitz","Dalitz Electron eta vs. eta",nptbins,ptmin,ptmax,netabins,etamin,etamax);
    fhPtWDecay = new TH1F("hPtWDecay","W-boson Electron pT",nptbins,ptmin,ptmax);
    fhPhiWDecay = new TH2F("hPhiWDecay","W-boson electron phi vs pT",nptbins,ptmin,ptmax,nphibins,phimin,phimax);
    fhEtaWDecay = new TH2F("hEtaWDecay","W-boson Electron eta vs. eta",nptbins,ptmin,ptmax,netabins,etamin,etamax);
    fhPtZDecay = new TH1F("hPtZDecay","Z-boson electron pT",nptbins,ptmin,ptmax);
    fhPhiZDecay = new TH2F("hPhiZDecay","Z-boson Electron phi vs pT",nptbins,ptmin,ptmax,nphibins,phimin,phimax);
    fhEtaZDecay = new TH2F("hEtaZDecay","Z-boson Electron eta vs. eta",nptbins,ptmin,ptmax,netabins,etamin,etamax);
    fhPtPrompt = new TH1F("hPtPrompt","Prompt electron pT",nptbins,ptmin,ptmax);
    fhPhiPrompt = new TH2F("hPhiPrompt","Prompt Electron phi vs pT",nptbins,ptmin,ptmax,nphibins,phimin,phimax);
    fhEtaPrompt = new TH2F("hEtaPrompt","Prompt Electron eta vs. eta",nptbins,ptmin,ptmax,netabins,etamin,etamax);
    fhPtUnknown = new TH1F("hPtUnknown","Unknown electron pT",nptbins,ptmin,ptmax);
    fhPhiUnknown = new TH2F("hPhiUnknown","Unknown Electron phi vs pT",nptbins,ptmin,ptmax,nphibins,phimin,phimax);
    fhEtaUnknown = new TH2F("hEtaUnknown","Unknown Electron eta vs. eta",nptbins,ptmin,ptmax,netabins,etamin,etamax);

    outputContainer->Add(fhPtConversion);
    outputContainer->Add(fhPhiConversion);
    outputContainer->Add(fhEtaConversion);
    outputContainer->Add(fhPtBottom);
    outputContainer->Add(fhPhiBottom);
    outputContainer->Add(fhEtaBottom);
    outputContainer->Add(fhPtCharm);
    outputContainer->Add(fhPhiCharm);
    outputContainer->Add(fhEtaCharm);
    outputContainer->Add(fhPtCFromB);
    outputContainer->Add(fhPhiCFromB);
    outputContainer->Add(fhEtaCFromB);
    outputContainer->Add(fhPtDalitz);
    outputContainer->Add(fhPhiDalitz);
    outputContainer->Add(fhEtaDalitz);
    outputContainer->Add(fhPtWDecay);
    outputContainer->Add(fhPhiWDecay);
    outputContainer->Add(fhEtaWDecay);
    outputContainer->Add(fhPtZDecay);
    outputContainer->Add(fhPhiZDecay);
    outputContainer->Add(fhEtaZDecay);
    outputContainer->Add(fhPtPrompt);
    outputContainer->Add(fhPhiPrompt);
    outputContainer->Add(fhEtaPrompt);
    outputContainer->Add(fhPtUnknown);
    outputContainer->Add(fhPhiUnknown);
    outputContainer->Add(fhEtaUnknown);

  /*
    fhMCElePt = new TH1F("hMCElePt","MC Electron pT",nptbins,ptmin,ptmax);
    fhMCElePhi = new TH2F("hMCElePhi","MC Electron phi",nptbins,ptmin,ptmax,nphibins,phimin,phimax);
    fhMCEleEta = new TH2F("hMCEleEta","MC Electron eta",nptbins,ptmin,ptmax,netabins,etamin,etamax);
    
    outputContainer->Add(fhMCElePt) ; 
    outputContainer->Add(fhMCElePhi) ; 
    outputContainer->Add(fhMCEleEta) ;
  */

  }//Histos with MC
  
  //Save parameters used for analysis
  TString parList ; //this will be list of parameters used for this analysis.
  char onePar[255] ;
  
  sprintf(onePar,"--- AliAnaElectron ---\n") ;
  parList+=onePar ;	
  sprintf(onePar,"fCalorimeter: %s\n",fCalorimeter.Data()) ;
  parList+=onePar ;  
  sprintf(onePar,"fpOverEmin: %f\n",fpOverEmin) ;
  parList+=onePar ;  
  sprintf(onePar,"fpOverEmax: %f\n",fpOverEmax) ;
  parList+=onePar ;  
  sprintf(onePar,"fResidualCut: %f\n",fResidualCut) ;
  parList+=onePar ;  
  
  //Get parameters set in base class.
  parList += GetBaseParametersList() ;
  
  //Get parameters set in FidutialCut class (not available yet)
  //parlist += GetFidCut()->GetFidCutParametersList() 
  
  TObjString *oString= new TObjString(parList) ;
  outputContainer->Add(oString);
  
  return outputContainer ;
  
}

//____________________________________________________________________________
void AliAnaElectron::Init()
{
  
  //Init
  //Do some checks
//  if(fCalorimeter == "PHOS" && !GetReader()->IsPHOSSwitchedOn()){
//    printf("AliAnaElectron::Init() - !!ABORT: You want to use PHOS in analysis but it is not read!! \n!!Check the configuration file!!\n");
//    abort();
//  }
//  else  if(fCalorimeter == "EMCAL" && !GetReader()->IsEMCALSwitchedOn()){
//    printf("AliAnaElectron::Init() - !!ABORT: You want to use EMCAL in analysis but it is not read!! \n!!Check the configuration file!!\n");
//    abort();
//  }
  
}


//____________________________________________________________________________
void AliAnaElectron::InitParameters()
{
  
  //Initialize the parameters of the analysis.
  SetOutputAODClassName("AliAODPWG4Particle");
  SetOutputAODName("PWG4Particle");

  AddToHistogramsName("AnaElectron_");

  fCalorimeter = "EMCAL" ;
  fpOverEmin = 0.5;
  fpOverEmax = 1.5;
  fResidualCut = 0.02;

}

//__________________________________________________________________
void  AliAnaElectron::MakeAnalysisFillAOD() 
{
  //
  // Do analysis and fill aods with electron candidates
  // These AODs will be used to do subsequent histogram filling
  //
  // Also fill some QA histograms
  //

  //Search for electrons in fCalorimeter 
  //TRefArray * caloData = new TRefArray(); 
  //TRefArray * ctsData = new TRefArray();
	cout<<"Event type "<<GetReader()->GetInputEvent()->GetName()<<endl;
	if((strcmp(GetReader()->GetInputEvent()->GetName(),"AliESDEvent"))) {
		printf("AliAnaElectron::MakeAnalysisFillAOD() - !!ABORT: Analysis working only with ESDs!!\n");
		abort();
	}

  //Get vertex for cluster momentum calculation
  Double_t vertex[]={0,0,0} ; //vertex ;
  if(GetReader()->GetDataType() != AliCaloTrackReader::kMC) GetReader()->GetVertex(vertex);
  //Get bfield for track extrapolation to Calorimeter
  Double_t bfield = GetReader()->GetInputEvent()->GetMagneticField();  //from V0 finder
  Double_t radius = 0.;

  //Get the CTS tracks
  //ctsData = GetAODCTS();
  //Select the Calorimeter of the electron
  if(fCalorimeter == "PHOS") {
    //caloData = GetAODPHOS();
    radius = 425.0;  //FIXME
  } else if (fCalorimeter == "EMCAL") {
    //caloData = GetAODEMCAL();
    radius = 441.0;  //[cm] EMCAL radius +13cm FIXME
  }

  //if(!(ctsData && caloData) || (ctsData->GetEntriesFast() == 0 || caloData->GetEntriesFast() == 0)) return;
  //if(!caloData ||  caloData->GetEntriesFast() == 0) return;

  ////////////////////////////////////////////////
  //Start from tracks and get associated clusters 
  //
  //Note: an alternative method would be to start from clusters and get associated tracks -
  //which is better?  For electrons, probably tracks-->clusters
  ////////////////////////////////////////////////
//  for(Int_t itrk = 0; itrk < ctsData->GetEntriesFast(); itrk++){
//    AliAODTrack *track = (AliAODTrack*)ctsData->At(itrk);

    AliESDEvent *esd = (AliESDEvent*) GetReader()->GetInputEvent();
    Int_t nTracks   = esd->GetNumberOfTracks() ;
    for (Int_t itrk =  0; itrk <  nTracks; itrk++) {////////////// track loop
      AliESDtrack * track = (AliESDtrack*) esd->GetTrack(itrk) ;
      //extrapolate track to Calorimeter
      Double_t emcmom[3] = {0.,0.,0.};
      Double_t emcpos[3] = {0.,0.,0.};
      AliExternalTrackParam *outerparam = (AliExternalTrackParam*)track->GetOuterParam();
      if(!outerparam) continue;
      
      Bool_t okpos = outerparam->GetXYZAt(radius,bfield,emcpos);
      Bool_t okmom = outerparam->GetPxPyPzAt(radius,bfield,emcmom);
      if(!(okpos && okmom)) continue;
      
      TVector3 pos(emcpos[0],emcpos[1],emcpos[2]);
      TVector3 mom(emcmom[0],emcmom[1],emcmom[2]);
      Double_t tphi = pos.Phi();
      Double_t teta = pos.Eta();
      Double_t tmom = mom.Mag();
      
      TLorentzVector mom2(mom,0.);
      Bool_t in =  GetFidutialCut()->IsInFidutialCut(mom2,fCalorimeter) ;
      if(GetDebug() > 1) printf("AliAnaElectron::MakeAnalysisFillAOD() - Track pt %2.2f, phi %2.2f, eta %2.2f in fidutial cut %d\n",track->Pt(), track->Phi(), track->Eta(), in);
      if(mom.Pt() > GetMinPt() && in) {

        printf("\tExtrapolated pt %2.2f, phi %2.2f, eta %2.2f \n",mom.Pt(),mom.Phi(),mom.Eta());
	fh2OuterPtVsExtrapPt->Fill(mom.Pt(),track->Pt());
	fh2OuterPhiVsExtrapPhi->Fill(mom.Phi(),track->Phi());
	fh2OuterEtaVsExtrapEta->Fill(mom.Eta(),track->Eta());

	Int_t ntot = esd->GetNumberOfCaloClusters();//caloData->GetEntriesFast();
	Double_t res = 999.;
	Double_t pOverE = -999.;
	
	//For tracks in EMCAL acceptance, pair them with all clusters
	//and fill the dEta vs dPhi for these pairs:
	for(Int_t iclus = 0; iclus < esd->GetNumberOfCaloClusters(); iclus++) {
	  AliESDCaloCluster * clus = (AliESDCaloCluster*) esd->GetCaloCluster(iclus);
	  if(!clus) continue;
	  if(fCalorimeter == "PHOS"  && !clus->IsPHOS())  continue;
	  if(fCalorimeter == "EMCAL" && !clus->IsEMCAL()) continue;
	  
	  Float_t x[3];
	  clus->GetPosition(x);
	  TVector3 cluspos(x[0],x[1],x[2]);
	  Double_t deta = teta - cluspos.Eta();
	  Double_t dphi = tphi - cluspos.Phi();
	  if(dphi > TMath::Pi()) dphi -= 2*TMath::Pi();
	  if(dphi < -TMath::Pi()) dphi += 2*TMath::Pi();
	  if(track->GetEMCALcluster() < -9000) fh2dEtadPhiUnmatched->Fill(deta,dphi);
	  fh2dEtadPhi->Fill(deta,dphi);
	  fh2TrackPVsClusterE->Fill(clus->E(),track->P());
	  fh2TrackPtVsClusterE->Fill(clus->E(),track->Pt());
	  fh2TrackPhiVsClusterPhi->Fill(cluspos.Phi(),mom.Phi());
	  fh2TrackEtaVsClusterEta->Fill(cluspos.Eta(),mom.Eta());
	}
	
	/////////////////////////////////
	//Perform electron cut analysis//
	/////////////////////////////////
	Bool_t isElectron = kFALSE;
	
	Int_t iCluster = track->GetEMCALcluster();
	if(iCluster < -9000) {printf("NOT MATCHED"); continue; }//no match
	if(iCluster > ntot) continue; //index out of range; shouldn't happen
	if(iCluster < 0 && iCluster > -9000) { //this should only happen in MC events
	  printf("AliAnaElectron::MakeAnalysisFillAOD() - Track has a fake match: %d\n",iCluster);
	  continue;
	}
	AliESDCaloCluster * clus = (AliESDCaloCluster*) esd->GetCaloCluster(iCluster);
	if(!clus) continue;
	if(fCalorimeter == "PHOS"  && !clus->IsPHOS())  continue;
	if(fCalorimeter == "EMCAL" && !clus->IsEMCAL()) continue;
	
	Float_t x[3];
	clus->GetPosition(x);
	TVector3 cluspos(x[0],x[1],x[2]);
	Double_t deta = teta - cluspos.Eta();
	Double_t dphi = tphi - cluspos.Phi();
	if(dphi > TMath::Pi()) dphi -= 2*TMath::Pi();
	if(dphi < -TMath::Pi()) dphi += 2*TMath::Pi();
	res = sqrt(dphi*dphi + deta*deta);
	fh1dR->Fill(res);
	fh2dEtadPhiMatched->Fill(deta,dphi);
	
	if(res < fResidualCut) {
	  //Good match
	  fh2MatchdEdx->Fill(track->P(),track->GetTPCsignal());
	  
	  Double_t energy = clus->E(); 
	  if(energy > 0) pOverE = tmom/energy;
	  fh1pOverE->Fill(pOverE);
	  
	  Int_t mult = clus->GetNumberOfDigits();
	  //	Int_t mcClus = clus->GetLabel();
	  AliESDCaloCluster * esdcalo = (AliESDCaloCluster*) esd->GetCaloCluster(clus->GetID());
	  Int_t matchIndex = esdcalo->GetTrackMatched();
	  
	  if(matchIndex != itrk) printf("Track and cluster don't agree! track %d, cluster %d",itrk,matchIndex);
	  if(mult < 2) printf("Single digit cluster.");
	  
	  //////////////////////////////
	  //Electron cuts happen here!//
	  //////////////////////////////
	  if(pOverE > fpOverEmin && pOverE < fpOverEmax) isElectron = kTRUE;
	  
	} //good matching residual
	
	///////////////////////////
	//Fill AOD with electrons//
	///////////////////////////
	if(isElectron) {
	  
	  fh2EledEdx->Fill(track->P(),track->GetTPCsignal());

	Double_t eMass = 0.511/1000; //mass in GeV
	Double_t eleE = sqrt(track->P()*track->P() + eMass*eMass);
	AliAODPWG4Particle tr = AliAODPWG4Particle(track->Px(),track->Py(),track->Pz(),eleE);
	tr.SetLabel(tr.GetLabel());
	tr.SetCaloLabel(clus->GetID(),-1); //sets the indices of the original caloclusters
	tr.SetDetector(fCalorimeter);
	tr.SetPdg(11);
	
	//Play with the MC stack if available
	//Check origin of the candidates
	if(IsDataMC()){
	  
	  //FIXME:  Need to re-think this for track-oriented analysis
	  tr.SetTag(GetMCAnalysisUtils()->CheckOrigin(clus->GetLabel(),GetMCStack()));
	  
	  if(GetDebug() > 0) printf("AliAnaElectron::MakeAnalysisFillAOD() - Origin of candidate %d\n",tr.GetTag());
	}//Work with stack also   

	AddAODParticle(tr);

	if(GetDebug() > 1) printf("AliAnaElectron::MakeAnalysisFillAOD() - Electron selection cuts passed: pT %3.2f, pdg %d\n",tr.Pt(),tr.GetPdg());	
      }//electron
    }//pt, fiducial selection                                                                                  
  }//track loop                         
  
  //FIXME:  Should we also check from the calocluster side, just in
  //case?

  if(GetDebug() > 1) printf("AliAnaElectron::MakeAnalysisFillAOD()  End fill AODs \n");  

}

//__________________________________________________________________
void  AliAnaElectron::MakeAnalysisFillHistograms() 
{
  //Do analysis and fill histograms
  
  //Loop on stored AOD electrons
  Int_t naod = GetOutputAODBranch()->GetEntriesFast();
  if(GetDebug() > 0) printf("AliAnaElectron::MakeAnalysisFillHistograms() - aod branch entries %d\n", naod);
  
  for(Int_t iaod = 0; iaod < naod ; iaod++){
    AliAODPWG4Particle* ele =  (AliAODPWG4Particle*) (GetOutputAODBranch()->At(iaod));
    Int_t pdg = ele->GetPdg();
    
    if(GetDebug() > 3) 
      printf("AliAnaElectron::MakeAnalysisFillHistograms() - PDG %d, MC TAG %d, Calorimeter %s\n", ele->GetPdg(),ele->GetTag(), (ele->GetDetector()).Data()) ;
    
    if(pdg != AliCaloPID::kElectron) continue; 
    if(ele->GetDetector() != fCalorimeter) continue;
    
    if(GetDebug() > 1) 
      printf("AliAnaElectron::MakeAnalysisFillHistograms() - ID Electron: pt %f, phi %f, eta %f\n", ele->Pt(),ele->Phi(),ele->Eta()) ;
    
    //Fill electron histograms 
    Float_t ptele = ele->Pt();
    Float_t phiele = ele->Phi();
    Float_t etaele = ele->Eta();
    
    fhPtElectron  ->Fill(ptele);
    fhPhiElectron ->Fill(ptele,phiele);
    fhEtaElectron ->Fill(ptele,etaele);
    
    if(IsDataMC()){
      Int_t tag = ele->GetTag();
      if(tag == AliMCAnalysisUtils::kMCConversion){
	fhPtConversion  ->Fill(ptele);
	fhPhiConversion ->Fill(ptele,phiele);
	fhEtaConversion ->Fill(ptele,etaele);
      }
      else if(tag==AliMCAnalysisUtils::kMCEFromB)
	{
	  fhPtBottom  ->Fill(ptele);
	  fhPhiBottom ->Fill(ptele,phiele);
	  fhEtaBottom ->Fill(ptele,etaele);
	}
      else if(tag==AliMCAnalysisUtils::kMCEFromC)
	{
	  fhPtCharm  ->Fill(ptele);
	  fhPhiCharm ->Fill(ptele,phiele);
	  fhEtaCharm ->Fill(ptele,etaele);
	}
      else if(tag==AliMCAnalysisUtils::kMCEFromCFromB)
	{
	  fhPtCFromB  ->Fill(ptele);
	  fhPhiCFromB ->Fill(ptele,phiele);
	  fhEtaCFromB ->Fill(ptele,etaele);
	}
      else if(tag==AliMCAnalysisUtils::kMCPi0Decay || tag==AliMCAnalysisUtils::kMCEtaDecay || tag==AliMCAnalysisUtils::kMCOtherDecay)
	{
	  fhPtDalitz  ->Fill(ptele);
	  fhPhiDalitz ->Fill(ptele,phiele);
	  fhEtaDalitz ->Fill(ptele,etaele);
	}
      else if(tag==AliMCAnalysisUtils::kMCWDecay)
	{
	  fhPtWDecay  ->Fill(ptele);
	  fhPhiWDecay ->Fill(ptele,phiele);
	  fhEtaWDecay ->Fill(ptele,etaele);
	}
      else if(tag==AliMCAnalysisUtils::kMCZDecay)
	{
	  fhPtZDecay  ->Fill(ptele);
	  fhPhiZDecay ->Fill(ptele,phiele);
	  fhEtaZDecay ->Fill(ptele,etaele);
	}
      else if(tag==AliMCAnalysisUtils::kMCElectron)
	{
	  fhPtPrompt  ->Fill(ptele);
	  fhPhiPrompt ->Fill(ptele,phiele);
	  fhEtaPrompt ->Fill(ptele,etaele);	  
	}
      else{
	fhPtUnknown  ->Fill(ptele);
	fhPhiUnknown ->Fill(ptele,phiele);
	fhEtaUnknown ->Fill(ptele,etaele);
      }
    }//Histograms with MC
    
  }// aod loop

  ////////////////////////////////////////////////////////
  //Fill histograms of pure MC kinematics from the stack//
  ////////////////////////////////////////////////////////
  if(IsDataMC()) {
    AliStack * stack =  GetMCStack() ;

    if(!stack)
      printf("AliAnaElectron::MakeAnalysisFillHistograms() *** no stack ***: \n");

    //FIXME:  Fill pure kine histograms here

  }

}


//__________________________________________________________________
void AliAnaElectron::Print(const Option_t * opt) const
{
  //Print some relevant parameters set for the analysis
  
  if(! opt)
    return;
  
  printf("**** Print %s %s ****\n", GetName(), GetTitle() ) ;
  AliAnaPartCorrBaseClass::Print(" ");

  printf("Calorimeter            =     %s\n", fCalorimeter.Data()) ;
  printf("pOverE range           =     %f - %f\n",fpOverEmin,fpOverEmax);
  printf("residual cut           =     %f\n",fResidualCut);
  printf("    \n") ;
	
} 

//________________________________________________________________________                          
void AliAnaElectron::ReadHistograms(TList* outputList)
{
  // Needed when Terminate is executed in distributed environment                             
  // Refill analysis histograms of this class with corresponding
  // histograms in output list.   

  // Histograms of this analsys are kept in the same list as other
  // analysis, recover the position of
  // the first one and then add the next                                                      
  Int_t index = outputList->IndexOf(outputList->FindObject(GetAddedHistogramsStringToName()+"fh1pOverE"));

  //Read histograms, must be in the same order as in
  //GetCreateOutputObject.                   
  fh1pOverE     = (TH1F *) outputList->At(index);
  fh1dR         = (TH1F *) outputList->At(index++);
  fh2EledEdx    = (TH2F *) outputList->At(index++);
  fh2MatchdEdx  = (TH2F *) outputList->At(index++);
  
}

//__________________________________________________________________                                
void  AliAnaElectron::Terminate(TList* outputList)
{

  //Do some plots to end
  //Recover histograms from output histograms list, needed for
  //distributed analysis.                
  //ReadHistograms(outputList);

  printf(" AliAnaElectron::Terminate()  *** %s Report: %d outputs", GetName(), outputList->GetEntries()) ;

}

