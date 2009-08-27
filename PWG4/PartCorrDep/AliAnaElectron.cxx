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
// -- Author: J.L. Klay (Cal Poly), M. Heinz (Yale)
//////////////////////////////////////////////////////////////////////////////
  
// --- ROOT system --- 
#include <TH2F.h>
#include <TParticle.h>
#include <TNtuple.h>
#include <TClonesArray.h>
#include <TObjString.h>
//#include <Riostream.h>

// --- Analysis system --- 
#include "AliAnaElectron.h" 
#include "AliCaloTrackReader.h"
#include "AliMCAnalysisUtils.h"
#include "AliAODCaloCluster.h"
#include "AliFidutialCut.h"
#include "AliAODTrack.h"
#include "AliAODPid.h"
#include "AliCaloPID.h"
#include "AliAODMCParticle.h"
#include "AliStack.h"

ClassImp(AliAnaElectron)
  
//____________________________________________________________________________
AliAnaElectron::AliAnaElectron() 
 : AliAnaPartCorrBaseClass(), fCalorimeter(""),
  fpOverEmin(0.),fpOverEmax(0.),fResidualCut(0.),
  fDrCut(0.),fPairDcaCut(0.),fDecayLenCut(0.),fImpactCut(0.),
  fAssocPtCut(0.),fMassCut(0.),fSdcaCut(0.),fITSCut(0),
  fWriteNtuple(kFALSE),
  //matching checks
  fEleNtuple(0),
  fh1pOverE(0),fh1dR(0),fh2EledEdx(0),fh2MatchdEdx(0),fh2dEtadPhi(0),
  fh2dEtadPhiMatched(0),fh2dEtadPhiUnmatched(0),
  fh2TrackPVsClusterE(0),fh2TrackPtVsClusterE(0),fh2TrackPhiVsClusterPhi(0),fh2TrackEtaVsClusterEta(0),
  //reco
  fhPtElectron(0),fhPhiElectron(0),fhEtaElectron(0),
  fhPtConversion(0),fhPhiConversion(0),fhEtaConversion(0),
  fhPtBottom(0),fhPhiBottom(0),fhEtaBottom(0),
  fhPtCharm(0),fhPhiCharm(0),fhEtaCharm(0),
  fhPtCFromB(0),fhPhiCFromB(0),fhEtaCFromB(0),
  fhPtDalitz(0),fhPhiDalitz(0),fhEtaDalitz(0),
  fhPtWDecay(0),fhPhiWDecay(0),fhEtaWDecay(0),
  fhPtZDecay(0),fhPhiZDecay(0),fhEtaZDecay(0),
  fhPtPrompt(0),fhPhiPrompt(0),fhEtaPrompt(0),
  fhPtUnknown(0),fhPhiUnknown(0),fhEtaUnknown(0),
  //B-tagging
  fhBtagCut1(0),fhBtagCut2(0),fhBtagCut3(0),
  //MC
  fMCEleNtuple(0)
{
  //default ctor
  
  //Initialize parameters
  InitParameters();

}

//____________________________________________________________________________
AliAnaElectron::AliAnaElectron(const AliAnaElectron & g) 
 : AliAnaPartCorrBaseClass(g), fCalorimeter(g.fCalorimeter),
   fpOverEmin(g.fpOverEmin),fpOverEmax(g.fpOverEmax),fResidualCut(g.fResidualCut),
   fDrCut(g.fDrCut),fPairDcaCut(g.fPairDcaCut),fDecayLenCut(g.fDecayLenCut),fImpactCut(g.fImpactCut),
   fAssocPtCut(g.fAssocPtCut),fMassCut(g.fMassCut),fSdcaCut(g.fSdcaCut),fITSCut(g.fITSCut),
   fWriteNtuple(g.fWriteNtuple),
   //matching checks
   fEleNtuple(g.fEleNtuple),
   fh1pOverE(g.fh1pOverE),fh1dR(g.fh1dR),
   fh2EledEdx(g.fh2EledEdx),fh2MatchdEdx(g.fh2MatchdEdx),fh2dEtadPhi(g.fh2dEtadPhi),
   fh2dEtadPhiMatched(g.fh2dEtadPhiMatched),fh2dEtadPhiUnmatched(g.fh2dEtadPhiUnmatched),
   fh2TrackPVsClusterE(g.fh2TrackPVsClusterE),fh2TrackPtVsClusterE(g.fh2TrackPtVsClusterE),
   fh2TrackPhiVsClusterPhi(g.fh2TrackPhiVsClusterPhi),fh2TrackEtaVsClusterEta(g.fh2TrackEtaVsClusterEta),   
   //reco
   fhPtElectron(g.fhPtElectron),fhPhiElectron(g.fhPhiElectron),fhEtaElectron(g.fhEtaElectron),
   fhPtConversion(g.fhPtConversion),fhPhiConversion(g.fhPhiConversion),fhEtaConversion(g.fhEtaConversion),
   fhPtBottom(g.fhPtBottom),fhPhiBottom(g.fhPhiBottom),fhEtaBottom(g.fhEtaBottom),
   fhPtCharm(g.fhPtCharm),fhPhiCharm(g.fhPhiCharm),fhEtaCharm(g.fhEtaCharm),
   fhPtCFromB(g.fhPtCFromB),fhPhiCFromB(g.fhPhiCFromB),fhEtaCFromB(g.fhEtaCFromB),
   fhPtDalitz(g.fhPtDalitz),fhPhiDalitz(g.fhPhiDalitz),fhEtaDalitz(g.fhEtaDalitz),
   fhPtWDecay(g.fhPtWDecay),fhPhiWDecay(g.fhPhiWDecay),fhEtaWDecay(g.fhEtaWDecay),
   fhPtZDecay(g.fhPtZDecay),fhPhiZDecay(g.fhPhiZDecay),fhEtaZDecay(g.fhEtaZDecay),
   fhPtPrompt(g.fhPtPrompt),fhPhiPrompt(g.fhPhiPrompt),fhEtaPrompt(g.fhEtaPrompt),
   fhPtUnknown(g.fhPtUnknown),fhPhiUnknown(g.fhPhiUnknown),fhEtaUnknown(g.fhEtaUnknown),
   //B-tagging
   fhBtagCut1(g.fhBtagCut1),fhBtagCut2(g.fhBtagCut2),fhBtagCut3(g.fhBtagCut3),
   //MC
   fMCEleNtuple(g.fMCEleNtuple)
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
  fDrCut = g.fDrCut;
  fPairDcaCut = g.fPairDcaCut;
  fDecayLenCut = g.fDecayLenCut;
  fImpactCut = g.fImpactCut;
  fAssocPtCut = g.fAssocPtCut;
  fMassCut = g.fMassCut;
  fSdcaCut = g.fSdcaCut;
  fITSCut = g.fITSCut;
  fWriteNtuple = g.fWriteNtuple;
  fEleNtuple = g.fEleNtuple;
  fh1pOverE = g.fh1pOverE;
  fh1dR = g.fh1dR;
  fh2EledEdx = g.fh2EledEdx;
  fh2MatchdEdx = g.fh2MatchdEdx;
  fh2dEtadPhi = g.fh2dEtadPhi;
  fh2dEtadPhiMatched = g.fh2dEtadPhiMatched;
  fh2dEtadPhiUnmatched = g.fh2dEtadPhiUnmatched;
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
  fMCEleNtuple = g.fMCEleNtuple;

  //B-tagging
  fhBtagCut1 = g.fhBtagCut1;
  fhBtagCut2 = g.fhBtagCut2;
  fhBtagCut3 = g.fhBtagCut3;

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
  
  //created ele ntuple for further analysis
  if(fWriteNtuple) {
    fEleNtuple = new TNtuple("EleNtuple","Electron Ntuple","mclabel:pt:phi:eta:p:E:deta:dphi:nCells:dEdx:eProb:impXY:impZ");//14 vars
    outputContainer->Add(fEleNtuple) ;
  }

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

  //B-tagging
  fhBtagCut1 = new TH2F("hbtag_cut1","B-tag result cut1", 10,0,10 ,nptbins,ptmin,ptmax);
  fhBtagCut2 = new TH2F("hbtag_cut2","B-tag result cut2", 10,0,10 ,nptbins,ptmin,ptmax);
  fhBtagCut3 = new TH2F("hbtag_cut3","B-tag result cut3", 10,0,10 ,nptbins,ptmin,ptmax);
  
  outputContainer->Add(fhBtagCut1) ;
  outputContainer->Add(fhBtagCut2) ;
  outputContainer->Add(fhBtagCut3) ;
  
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
    
    //created ele ntuple for further analysis
    if(fWriteNtuple) {
      fMCEleNtuple = new TNtuple("MCEleNtuple","MC Electron Ntuple","mclabel:pt:phi:eta:x:y:z");
      outputContainer->Add(fMCEleNtuple) ;
    }

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

  //do some initialization
  if(fCalorimeter == "PHOS") {
    printf("AliAnaElectron::Init() - !!STOP: You want to use PHOS in analysis but this is not (yet) supported!!\n!!Check the configuration file!!\n");
    fCalorimeter = "EMCAL";
  }
  if(fCalorimeter == "EMCAL" && !GetReader()->IsEMCALSwitchedOn()){
    printf("AliAnaElectron::Init() - !!STOP: You want to use EMCAL in analysis but it is not read!!\n!!Check the configuration file!!\n");
    abort();
  }

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
  //B-tagging
  fDrCut       = 1.0; 
  fPairDcaCut  = 0.02;
  fDecayLenCut = 1.0;
  fImpactCut   = 0.5;
  fAssocPtCut  = 1.0;
  fMassCut     = 1.5;
  fSdcaCut     = 0.1;
  fITSCut      = 4;

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

  TObjArray *cl = new TObjArray();

  //Get vertex for cluster momentum calculation
  Double_t vertex[]={0,0,0} ; //vertex ;
  if(GetReader()->GetDataType() != AliCaloTrackReader::kMC) GetReader()->GetVertex(vertex);

  //Select the calorimeter of the electron
  if(fCalorimeter == "PHOS") {
    cl = GetAODPHOS();
  }
  else if (fCalorimeter == "EMCAL") {
    cl = GetAODEMCAL();
  }

  ////////////////////////////////////////////////
  //Start from tracks and get associated clusters 
  ////////////////////////////////////////////////
  if(!GetAODCTS() || GetAODCTS()->GetEntriesFast() == 0) return ;
  Int_t ntracks = GetAODCTS()->GetEntriesFast();
  if(GetDebug() > 0)
    printf("AliAnaElectron::MakeAnalysisFillAOD() - In CTS aod entries %d\n", ntracks);

  //Unfortunately, AliAODTracks don't have associated EMCAL clusters.
  //we have to redo track-matching, I guess
  Int_t iCluster = -999;

  for (Int_t itrk =  0; itrk <  ntracks; itrk++) {////////////// track loop
    iCluster = -999; //start with no match
    AliAODTrack * track = (AliAODTrack*) (GetAODCTS()->At(itrk)) ;
    AliAODPid* pid = (AliAODPid*) track->GetDetPid();

    Double_t emcpos[3];
    pid->GetEMCALPosition(emcpos);
    Double_t emcmom[3];
    pid->GetEMCALMomentum(emcmom);
    
    TVector3 pos(emcpos[0],emcpos[1],emcpos[2]);
    TVector3 mom(emcmom[0],emcmom[1],emcmom[2]);
    Double_t tphi = pos.Phi();
    Double_t teta = pos.Eta();
    Double_t tmom = mom.Mag();

    TLorentzVector mom2(mom,0.);
    Bool_t in =  GetFidutialCut()->IsInFidutialCut(mom2,fCalorimeter) ;
    if(GetDebug() > 1) printf("AliAnaElectron::MakeAnalysisFillAOD() - Track pt %2.2f, phi %2.2f, eta %2.2f in fidutial cut %d\n",track->Pt(), track->Phi(), track->Eta(), in);
    if(mom.Pt() > GetMinPt() && in) {
            
      Int_t ntot = cl->GetEntriesFast();
      Double_t res = 999.;
      Double_t pOverE = -999.;

      Bool_t isElectron = kFALSE;      
      //For tracks in EMCAL acceptance, pair them with all clusters
      //and fill the dEta vs dPhi for these pairs:
      for(Int_t iclus = 0; iclus < ntot; iclus++) {
	AliAODCaloCluster * clus = (AliAODCaloCluster*) (cl->At(iclus));
	if(!clus) continue;
	
	Float_t x[3];
	clus->GetPosition(x);
	TVector3 cluspos(x[0],x[1],x[2]);
	Double_t deta = teta - cluspos.Eta();
	Double_t dphi = tphi - cluspos.Phi();
	if(dphi > TMath::Pi()) dphi -= 2*TMath::Pi();
	if(dphi < -TMath::Pi()) dphi += 2*TMath::Pi();
	fh2dEtadPhi->Fill(deta,dphi);
	fh2TrackPVsClusterE->Fill(clus->E(),track->P());
	fh2TrackPtVsClusterE->Fill(clus->E(),track->Pt());
	fh2TrackPhiVsClusterPhi->Fill(cluspos.Phi(),mom.Phi());
	fh2TrackEtaVsClusterEta->Fill(cluspos.Eta(),mom.Eta());

	res = sqrt(dphi*dphi + deta*deta);
	fh1dR->Fill(res);
	fh2dEtadPhiMatched->Fill(deta,dphi);

	/////////////////////////////////
	//Perform electron cut analysis//
	/////////////////////////////////
	//Good match
	if(res < fResidualCut) {
	  iCluster = iclus;

	  Double_t xyz[3];
	  track->XYZAtDCA(xyz);
	  Float_t xy = TMath::Sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
	  Float_t z = xyz[2];

	  //Double_t tpcProb[5]; //e,pi,mu,k,p
	  //pid->GetTPCpid(tpcProb);
	  
	  Float_t mclabel = -1;

	  if(IsDataMC()) {
	    //Input from second AOD?
		Int_t input = 0;
	    if(GetReader()->GetAODCTSNormalInputEntries() <= itrk) input = 1;
	    mclabel = (Float_t)GetMCAnalysisUtils()->CheckOrigin(track->GetLabel(),GetReader(),input);
		//Do you want the cluster or the track label?
//		if(GetReader()->GetAODEMCALNormalInputEntries() <= iclus) input = 1;
//		mclabel = (Float_t)GetMCAnalysisUtils()->CheckOrigin(clus->GetLabel(0),GetReader(),input);
		  
	  }

	  if(fWriteNtuple) {
	    fEleNtuple->Fill(mclabel,track->Pt(),track->Phi(),track->Eta(),track->P(),clus->E(),deta,dphi,clus->GetNCells(),pid->GetTPCsignal(),0.,xy,z);
	  }
	  
	  fh2MatchdEdx->Fill(track->P(),pid->GetTPCsignal());
	  
	  Double_t energy = clus->E(); 
	  if(energy > 0) pOverE = tmom/energy;
	  fh1pOverE->Fill(pOverE);
	  
	  Int_t mult = clus->GetNCells();
	  if(mult < 2 &&  GetDebug() > 0) printf("Single digit cluster.\n");
	  
	  //////////////////////////////
	  //Electron cuts happen here!//
	  //////////////////////////////
	  if(pOverE > fpOverEmin && pOverE < fpOverEmax) isElectron = kTRUE;
	} else {
	  fh2dEtadPhiUnmatched->Fill(deta,dphi);
	}
	  
      } //calocluster loop

      ///////////////////////////
      //Fill AOD with electrons//
      ///////////////////////////
      if(isElectron) {

	//B-tagging
	Int_t btag = GetBtag(track);
	
	fh2EledEdx->Fill(track->P(),pid->GetTPCsignal());
	
	Double_t eMass = 0.511/1000; //mass in GeV
	Double_t eleE = sqrt(track->P()*track->P() + eMass*eMass);
	AliAODPWG4Particle tr = AliAODPWG4Particle(track->Px(),track->Py(),track->Pz(),eleE);
	tr.SetLabel(track->GetLabel());
	tr.SetCaloLabel(iCluster,-1); //sets the indices of the original caloclusters
	tr.SetTrackLabel(track->GetID(),-1); //sets the indices of the original tracks
	tr.SetDetector(fCalorimeter);
	if(GetReader()->GetAODCTSNormalInputEntries() <= itrk) tr.SetInputFileIndex(1);
	tr.SetPdg(11);
	tr.SetBtag(btag);

	//Play with the MC stack if available
	//Check origin of the candidates
	if(IsDataMC()){
	  
	  //FIXME:  Need to re-think this for track-oriented analysis
	  tr.SetTag(GetMCAnalysisUtils()->CheckOrigin(tr.GetLabel(),GetReader(),tr.GetInputFileIndex()));
	  
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

  AliStack * stack = 0x0;
  TParticle * primary = 0x0;
  TClonesArray * mcparticles0 = 0x0;
  TClonesArray * mcparticles1 = 0x0;
  AliAODMCParticle * aodprimary = 0x0;

  if(IsDataMC()) {
    if(GetReader()->ReadStack()){
      stack =  GetMCStack() ;
      
      if(!stack)
	printf("AliAnaElectron::MakeAnalysisFillHistograms() *** no stack ***: \n");
      
    }
    else if(GetReader()->ReadAODMCParticles()){
      //Get the list of MC particles
      mcparticles0 = GetReader()->GetAODMCParticles(0);
      if(!mcparticles0 && GetDebug() > 0)     {
	printf("AliAnaElectron::MakeAnalysisFillHistograms() -  Standard MCParticles not available!\n");
      }
      if(GetReader()->GetSecondInputAODTree()){
	mcparticles1 = GetReader()->GetAODMCParticles(1);
	if(!mcparticles1 && GetDebug() > 0)     {
	  printf("AliAnaElectron::MakeAnalysisFillHistograms() -  Second input MCParticles not available!\n");
	}
      }
      
    }
  }// is data and MC
  
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
      if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion)){
	fhPtConversion  ->Fill(ptele);
	fhPhiConversion ->Fill(ptele,phiele);
	fhEtaConversion ->Fill(ptele,etaele);
      }
      else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEFromB))
	{
	  fhPtBottom  ->Fill(ptele);
	  fhPhiBottom ->Fill(ptele,phiele);
	  fhEtaBottom ->Fill(ptele,etaele);
	}
      else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEFromC))
	{
	  fhPtCharm  ->Fill(ptele);
	  fhPhiCharm ->Fill(ptele,phiele);
	  fhEtaCharm ->Fill(ptele,etaele);
	}
      else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEFromCFromB))
	{
	  fhPtCFromB  ->Fill(ptele);
	  fhPhiCFromB ->Fill(ptele,phiele);
	  fhEtaCFromB ->Fill(ptele,etaele);
	}
      else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0Decay) || GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEtaDecay) || GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCOtherDecay))
	{
	  fhPtDalitz  ->Fill(ptele);
	  fhPhiDalitz ->Fill(ptele,phiele);
	  fhEtaDalitz ->Fill(ptele,etaele);
	}
      else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCWDecay))
	{
	  fhPtWDecay  ->Fill(ptele);
	  fhPhiWDecay ->Fill(ptele,phiele);
	  fhEtaWDecay ->Fill(ptele,etaele);
	}
      else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCZDecay))
	{
	  fhPtZDecay  ->Fill(ptele);
	  fhPhiZDecay ->Fill(ptele,phiele);
	  fhEtaZDecay ->Fill(ptele,etaele);
	}
      else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCElectron))
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
    if(GetReader()->ReadStack()) {
      for(Int_t ipart = 0; ipart < stack->GetNtrack(); ipart++) {
	primary = stack->Particle(ipart);
	Int_t pdgcode = primary->GetPdgCode();
	//we only care about electrons
	if(abs(pdgcode) != 11) continue;
	//we only want TRACKABLE electrons (TPC 85-250cm)
	if(primary->R() > 200.) continue;
	
	//find out what the ancestry of this electron is
	Float_t mclabel = -1;
	Int_t input = 0;
	mclabel = (Float_t)GetMCAnalysisUtils()->CheckOrigin(ipart,GetReader(),input);

	//fill ntuple
	if(fWriteNtuple) {
	  fMCEleNtuple->Fill(mclabel,primary->Pt(),primary->Phi(),primary->Eta(),primary->Vx(),primary->Vy(),primary->Vz());
	}
	
      }
      
    } else if(GetReader()->ReadAODMCParticles()) {
      Int_t npart0 = mcparticles0->GetEntriesFast();
	  Int_t npart1 = 0;
	  if(mcparticles1) npart1 = mcparticles1->GetEntriesFast();
      Int_t npart = npart0+npart1;
      for(Int_t ipart = 0; ipart < npart; ipart++) {
	if(ipart < npart0) aodprimary = (AliAODMCParticle*)mcparticles0->At(ipart);
	else aodprimary = (AliAODMCParticle*)mcparticles1->At(ipart-npart0);
	if(!aodprimary) {
	  printf("AliAnaElectron::MakeAnalysisFillHistograms() *** no primary ***:  label %d \n", ipart);
	  continue;
	}
	Int_t pdgcode = aodprimary->GetPdgCode();
	//we only care about electrons
	if(abs(pdgcode) != 11) continue;
	//we only want TRACKABLE electrons (TPC 85-250cm)
	Double_t radius = TMath::Sqrt(aodprimary->Xv()*aodprimary->Xv() + aodprimary->Yv()*aodprimary->Yv());
	if(radius > 200.) continue;
	
	//find out what the ancestry of this electron is
	Float_t mclabel = -1;
	Int_t input = 0;
	Int_t ival = ipart;
	if(ipart > npart0) { ival -= npart0; input = 1;}
	mclabel = (Float_t)GetMCAnalysisUtils()->CheckOrigin(ival,GetReader(),input);
	
	//fill ntuple
	if(fWriteNtuple) {
	  fMCEleNtuple->Fill(mclabel,aodprimary->Pt(),aodprimary->Phi(),aodprimary->Eta(),aodprimary->Xv(),aodprimary->Yv(),aodprimary->Zv());
	}
	
      }
    }
  } //pure MC kine histos
    
}

//__________________________________________________________________
Int_t AliAnaElectron::GetBtag(AliAODTrack * tr )
{
  
  //UChar_t itsmap = tr->GetITSClusterMap();
  //JLK 
  //DON'T KNOW HOW TO USE THIS???
  //if (itsmap < fITSCut) return 0;
  //JLK
  Double_t x[3];
  Bool_t gotit = tr->XYZAtDCA(x);
  if(!gotit) { printf("##Problem getting impact parameter!"); return 0; }

  Double_t d1 = TMath::Sqrt(x[0]*x[0] + x[1]*x[1]);
  if (TMath::Abs(d1) > fImpactCut ) return 0;
  if (TMath::Abs(x[2]) > fImpactCut ) return 0;

  Double_t p1[3]={0,0,0};
  gotit = tr->PxPyPzAtDCA(p1);
  if(!gotit) { printf("##Problem getting inner P"); return 0;}

  TVector3 mom1(p1[0],p1[1],p1[2]) ;
  int nvtx1 = 0;
  int nvtx2 = 0;
  int nvtx3 = 0;

  for (int k2 =0; k2 < GetAODCTS()->GetEntriesFast() ; k2++) {
    //loop over assoc
    AliAODTrack* track2 = (AliAODTrack*) (GetAODCTS()->At(k2));
    int id1 = tr->GetID();
    int id2 = track2->GetID();
    if(id1 == id2) continue;

    gotit = tr->XYZAtDCA(x);
    if(!gotit) { printf("##Problem getting impact parameter!"); continue; }

    d1 = TMath::Sqrt(x[0]*x[0] + x[1]*x[1]);
    if (TMath::Abs(d1) > fImpactCut ) continue;
    if (TMath::Abs(x[2]) > fImpactCut ) continue;

    //JLK
    //HOW TO IMPLEMENT?
    //if (track2->GetITSclusters(0) < fITSCut) continue;
    //JLK

    Double_t p2[3]={0,0,0};
    gotit = track2->PxPyPzAtDCA(p2);
    if(!gotit) { printf("###Problem getting inner P"); continue; }

    TVector3 mom2(p2[0],p2[1],p2[2]) ;
    if (mom2.Perp() <fAssocPtCut) continue;
    
    float dphi = mom2.Phi() - mom1.Phi();
    float deta = mom2.Eta() - mom1.Eta();
    float dr = sqrt(deta*deta + dphi*dphi);

    if (dr> fDrCut) continue;
    
    Double_t sDca1 = ComputeSignDca(tr, track2, 1.0);
    printf("### sec vtx-signdca :%f",sDca1);
    if (sDca1 > 0.1) nvtx1++;
    Double_t sDca2 = ComputeSignDca(tr, track2, 1.5);
    printf("### sec vtx-signdca :%f",sDca2);
    if (sDca2 > 0.1) nvtx2++;
    Double_t sDca3 = ComputeSignDca(tr, track2, 1.8);
    printf("### sec vtx-signdca :%f",sDca3);
    if (sDca3 > 0.1) nvtx3++;
  } //loop over hadrons

  if (nvtx1>0) printf("result1 of btagging: %d \n",nvtx1);
  if (nvtx2>0) printf("result2 of btagging: %d \n",nvtx2);
  if (nvtx3>0) printf("result3 of btagging: %d \n",nvtx3);

  //fill QA histograms
  fhBtagCut1->Fill(nvtx1,mom1.Perp());
  fhBtagCut2->Fill(nvtx2,mom1.Perp());
  fhBtagCut3->Fill(nvtx3,mom1.Perp());

  return nvtx2;
}

//__________________________________________________________________
Double_t AliAnaElectron::ComputeSignDca(AliAODTrack *tr, AliAODTrack *tr2 , float masscut)
{
  //Compute the signed dca between two tracks
  //and return the result
  //This method needs to be fixed to work with AODS

  Double_t signDca=-999;
  //Do something to avoid compiler warning
  if(GetDebug() > 1 ) printf("Labels: track1 %d, track2 %d, masscut %f", tr->GetLabel(), tr2->GetLabel(), masscut);
  /*
  //=====Now calculate DCA between both tracks=======  
  float massE = 0.000511;
  float massK = 0.493677;

  //No bfield member in AliCaloTrackReader...
  float bfield = ev->GetMagneticField();

  Double_t vertex[3] ;
  if(!GetReader()->GetDataType()== AliCaloTrackReader::kMC) GetReader()->GetVertex(vertex);
  TVector3 PV(vertex[0],vertex[1],vertex[2]) ;

  double xplane1=0; double xplane2=0;
  double pairdca = tr->GetDCA(tr2,bfield,xplane1,xplane2);

  AliExternalTrackParam innerparam1(*tr);
  AliExternalTrackParam innerparam2(*track2);
  AliESDv0 bvertex(innerparam1,id1,innerparam2,id2);
  double vx,vy,vz;
  bvertex.GetXYZ(vx,vy,vz);
  double emom[3];
  double hmom[3];
  innerparam1.GetPxPyPz(emom);
  innerparam2.GetPxPyPz(hmom);
  TVector3 emomAtB(emom[0],emom[1],emom[2]);
  TVector3 hmomAtB(hmom[0],hmom[1],hmom[2]);
  TVector3 secvtxpt(vx,vy,vz);
  TVector3 decayvector(0,0,0);
  decayvector = secvtxpt - PV; //decay vector from PrimVtx
  float decaylength = decayvector.Mag();

  if (emomAtB.Mag()>0 && pairdca < fPairDcaCut && decaylength < fDecayLenCut ) {
    TVector3 sumMom = emomAtB+hmomAtB;
    float ener1 = sqrt(pow(emomAtB.Mag(),2) + massE*massE);
    float ener2 = sqrt(pow(hmomAtB.Mag(),2) + massK*massK);
    float mass = sqrt(pow((ener1+ener2),2) - pow(sumMom.Mag(),2));
    if (mass > masscut) signDca = decayvector.Dot(emomAtB)/emomAtB.Mag();
    printf("### sec vtx-signdca :%f",signDca);
  }
  */
  
  return signDca;
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

  printf(" AliAnaElectron::Terminate()  *** %s Report: %d outputs\n", GetName(), outputList->GetEntries()) ;

}

