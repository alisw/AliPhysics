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
#include <TH3F.h>
#include <TParticle.h>
#include <TNtuple.h>
#include <TClonesArray.h>
#include <TObjString.h>
//#include <Riostream.h>

// --- Analysis system --- 
#include "AliAnaElectron.h" 
#include "AliCaloTrackReader.h"
#include "AliMCAnalysisUtils.h"
#include "AliVCluster.h"
#include "AliFiducialCut.h"
#include "AliAODTrack.h"
#include "AliAODPid.h"
#include "AliCaloPID.h"
#include "AliAODMCParticle.h"
#include "AliStack.h"
#include "AliExternalTrackParam.h"
#include "AliESDv0.h"
#include "AliESDtrack.h"
#include "AliAODJet.h"
#include "AliAODEvent.h"
#include "AliGenPythiaEventHeader.h"

ClassImp(AliAnaElectron)
  
//____________________________________________________________________________
AliAnaElectron::AliAnaElectron() 
: AliAnaPartCorrBaseClass(),fCalorimeter(""),
  fpOverEmin(0.),fpOverEmax(0.),fResidualCut(0.),fMinClusEne(0.),
  fDrCut(0.),fPairDcaCut(0.),fDecayLenCut(0.),fImpactCut(0.),
  fAssocPtCut(0.),fMassCut(0.),fSdcaCut(0.),fITSCut(0),
  fNTagTrkCut(0),fIPSigCut(0.),fJetEtaCut(0.3),fJetPhiMin(1.8),fJetPhiMax(2.9),
  fWriteNtuple(kFALSE),
  //event QA histos
  fhImpactXY(0),fhRefMult(0),fhRefMult2(0),
  //matching checks
  fh3pOverE(0),fh3EOverp(0),fh3pOverE2(0),fh3EOverp2(0),fh3pOverE3(0),fh3EOverp3(0),
  fh2pOverE(0),fh2EOverp(0),fh2pOverE2(0),fh2EOverp2(0),
  fh1dR(0),fh2EledEdx(0),fh2MatchdEdx(0),fh2dEtadPhi(0),
  fh2dEtadPhiMatched(0),fh2dEtadPhiUnmatched(0),fh2TrackPVsClusterE(0),
  fh2TrackPtVsClusterE(0),fh2TrackPhiVsClusterPhi(0),fh2TrackEtaVsClusterEta(0),
  //Photonic electron checks
  fh1OpeningAngle(0),fh1MinvPhoton(0),
  //Reconstructed electrons
  fhPtElectron(0),fhPhiElectron(0),fhEtaElectron(0),
  fhPtNPE(0),fhPhiNPE(0),fhEtaNPE(0),
  fhPtPE(0),fhPhiPE(0),fhEtaPE(0),
  //for comparisons with tracking detectors
  fhPtHadron(0),fhPtNPEleTPC(0),fhPtNPEleTPCTRD(0),fhPtNPEleTTE(0),
  fhPtNPEleEMCAL(0),
  //DVM B-tagging
  fhDVMBtagCut1(0),fhDVMBtagCut2(0),fhDVMBtagCut3(0),fhDVMBtagQA1(0),fhDVMBtagQA2(0),
  fhDVMBtagQA3(0),fhDVMBtagQA4(0),fhDVMBtagQA5(0),
  //IPSig B-tagging
  fhIPSigBtagQA1(0),fhIPSigBtagQA2(0),fhTagJetPt1x4(0),fhTagJetPt2x3(0),fhTagJetPt3x2(0),
  fhePlusTagJetPt1x4(0),fhePlusTagJetPt2x3(0),fhePlusTagJetPt3x2(0),
  //B-Jet histograms
  fhJetType(0),fhLeadJetType(0),fhBJetXsiFF(0),fhBJetPtFF(0),fhBJetEtaPhi(0),
  fhNonBJetXsiFF(0),fhNonBJetPtFF(0),fhNonBJetEtaPhi(0),
  /////////////////////////////////////////////////////////////
  //Histograms that rely on MC info (not filled for real data)
  fEleNtuple(0),
  //reco electrons from various sources
  fhPhiConversion(0),fhEtaConversion(0),
  //for comparisons with tracking detectors
  fhPtTrack(0),
  fhPtNPEBHadron(0),
  //for computing efficiency of B-jet tags
  fhBJetPt1x4(0),fhBJetPt2x3(0),fhBJetPt3x2(0),
  fhFakeJetPt1x4(0),fhFakeJetPt2x3(0),fhFakeJetPt3x2(0),fhDVMJet(0),
  //MC rate histograms/ntuple
  fMCEleNtuple(0),fhMCBJetElePt(0),fhMCBHadronElePt(0),fhPtMCHadron(0),fhPtMCElectron(0),
  fhMCXYConversion(0),fhMCRadPtConversion(0)
{
  //default ctor
  
  //Initialize parameters
  InitParameters();

}
/*
//____________________________________________________________________________
AliAnaElectron::AliAnaElectron(const AliAnaElectron & g) 
  : AliAnaPartCorrBaseClass(g),fCalorimeter(g.fCalorimeter),
    fpOverEmin(g.fpOverEmin),fpOverEmax(g.fpOverEmax),
    fResidualCut(g.fResidualCut),fMinClusEne(g.fMinClusEne),
    fDrCut(g.fDrCut),fPairDcaCut(g.fPairDcaCut),fDecayLenCut(g.fDecayLenCut),fImpactCut(g.fImpactCut),
    fAssocPtCut(g.fAssocPtCut),fMassCut(g.fMassCut),fSdcaCut(g.fSdcaCut),fITSCut(g.fITSCut),
    fNTagTrkCut(g.fNTagTrkCut),fIPSigCut(g.fIPSigCut),
    fJetEtaCut(g.fJetEtaCut),fJetPhiMin(g.fJetPhiMin),fJetPhiMax(g.fJetPhiMax),
    fWriteNtuple(g.fWriteNtuple),
    //event QA histos
    fhImpactXY(g.fhImpactXY),fhRefMult(g.fhRefMult),fhRefMult2(g.fhRefMult2),
    //matching checks
    fh3pOverE(g.fh3pOverE),fh3EOverp(g.fh3EOverp),
    fh3pOverE2(g.fh3pOverE2),fh3EOverp2(g.fh3EOverp2),
    fh3pOverE3(g.fh3pOverE3),fh3EOverp3(g.fh3EOverp3),
    fh2pOverE(g.fh2pOverE),fh2EOverp(g.fh2EOverp),
    fh2pOverE2(g.fh2pOverE2),fh2EOverp2(g.fh2EOverp2),
    fh1dR(g.fh1dR),fh2EledEdx(g.fh2EledEdx),
    fh2MatchdEdx(g.fh2MatchdEdx),fh2dEtadPhi(g.fh2dEtadPhi),
    fh2dEtadPhiMatched(g.fh2dEtadPhiMatched),fh2dEtadPhiUnmatched(g.fh2dEtadPhiUnmatched),
    fh2TrackPVsClusterE(g.fh2TrackPVsClusterE),fh2TrackPtVsClusterE(g.fh2TrackPtVsClusterE),
    fh2TrackPhiVsClusterPhi(g.fh2TrackPhiVsClusterPhi),fh2TrackEtaVsClusterEta(g.fh2TrackEtaVsClusterEta),
    //Photonic electron checks
    fh1OpeningAngle(g.fh1OpeningAngle),fh1MinvPhoton(g.fh1MinvPhoton),
    //Reconstructed electrons
    fhPtElectron(g.fhPtElectron),fhPhiElectron(g.fhPhiElectron),fhEtaElectron(g.fhEtaElectron),
    fhPtNPE(g.fhPtNPE),fhPhiNPE(g.fhPhiNPE),fhEtaNPE(g.fhEtaNPE),
    fhPtPE(g.fhPtPE),fhPhiPE(g.fhPhiPE),fhEtaPE(g.fhEtaPE),
    //for comparisons with tracking detectors
    fhPtHadron(g.fhPtHadron),fhPtNPEleTPC(g.fhPtNPEleTPC),
    fhPtNPEleTPCTRD(g.fhPtNPEleTPCTRD),fhPtNPEleTTE(g.fhPtNPEleTTE),
    fhPtNPEleEMCAL(g.fhPtNPEleEMCAL),
    //DVM B-tagging
    fhDVMBtagCut1(g.fhDVMBtagCut1),fhDVMBtagCut2(g.fhDVMBtagCut2),fhDVMBtagCut3(g.fhDVMBtagCut3),
    fhDVMBtagQA1(g.fhDVMBtagQA1),fhDVMBtagQA2(g.fhDVMBtagQA2),
    fhDVMBtagQA3(g.fhDVMBtagQA3),fhDVMBtagQA4(g.fhDVMBtagQA4),fhDVMBtagQA5(g.fhDVMBtagQA5),
    //IPSig B-tagging
    fhIPSigBtagQA1(g.fhIPSigBtagQA1),fhIPSigBtagQA2(g.fhIPSigBtagQA2),
    fhTagJetPt1x4(g.fhTagJetPt1x4),fhTagJetPt2x3(g.fhTagJetPt2x3),fhTagJetPt3x2(g.fhTagJetPt3x2),
    fhePlusTagJetPt1x4(g.fhePlusTagJetPt1x4),fhePlusTagJetPt2x3(g.fhePlusTagJetPt2x3),
    fhePlusTagJetPt3x2(g.fhePlusTagJetPt3x2),
    //B-Jet histograms
    fhJetType(g.fhJetType),fhLeadJetType(g.fhLeadJetType),fhBJetXsiFF(g.fhBJetXsiFF),
    fhBJetPtFF(g.fhBJetPtFF),fhBJetEtaPhi(g.fhBJetEtaPhi),fhNonBJetXsiFF(g.fhNonBJetXsiFF),
    fhNonBJetPtFF(g.fhNonBJetPtFF),fhNonBJetEtaPhi(g.fhNonBJetEtaPhi),
    /////////////////////////////////////////////////////////////
    //Histograms that rely on MC info (not filled for real data)
    fEleNtuple(g.fEleNtuple),
    //reco electrons from various sources
    fhPhiConversion(g.fhPhiConversion),fhEtaConversion(g.fhEtaConversion),
    //for comparisons with tracking detectors
    fhPtTrack(g.fhPtTrack),
    fhPtNPEBHadron(g.fhPtNPEBHadron),
    //for computing efficiency of B-jet tags
    fhBJetPt1x4(g.fhBJetPt1x4),fhBJetPt2x3(g.fhBJetPt2x3),
    fhBJetPt3x2(g.fhBJetPt3x2),
    fhFakeJetPt1x4(g.fhFakeJetPt1x4),fhFakeJetPt2x3(g.fhBJetPt2x3),
    fhFakeJetPt3x2(g.fhFakeJetPt3x2),fhDVMJet(g.fhDVMJet),
    //MC rate histograms/ntuple
    fMCEleNtuple(g.fMCEleNtuple),fhMCBJetElePt(g.fhMCBJetElePt),
    fhMCBHadronElePt(g.fhMCBHadronElePt),
    fhPtMCHadron(g.fhPtMCHadron),fhPtMCElectron(g.fhPtMCElectron),
    fhMCXYConversion(g.fhMCXYConversion),fhMCRadPtConversion(g.fhMCRadPtConversion)
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
  fMinClusEne = g.fMinClusEne;
  fDrCut = g.fDrCut;
  fPairDcaCut = g.fPairDcaCut;
  fDecayLenCut = g.fDecayLenCut;
  fImpactCut = g.fImpactCut;
  fAssocPtCut = g.fAssocPtCut;
  fMassCut = g.fMassCut;
  fSdcaCut = g.fSdcaCut;
  fITSCut = g.fITSCut;
  fNTagTrkCut = g.fNTagTrkCut;
  fIPSigCut = g.fIPSigCut;
  fJetEtaCut = g.fJetEtaCut;
  fJetPhiMin = g.fJetPhiMin;
  fJetPhiMax = g.fJetPhiMax;
  fWriteNtuple = g.fWriteNtuple;
  //event QA histos
  fhImpactXY = g.fhImpactXY;
  fhRefMult  = g.fhRefMult;
  fhRefMult2 = g.fhRefMult2;
  //matching checks
  fh3pOverE  = g.fh3pOverE;
  fh3EOverp  = g.fh3EOverp;
  fh3pOverE2 = g.fh3pOverE2;
  fh3EOverp2 = g.fh3EOverp2;
  fh3pOverE3 = g.fh3pOverE3;
  fh3EOverp3 = g.fh3EOverp3;
  fh2pOverE  = g.fh2pOverE;
  fh2EOverp  = g.fh2EOverp;
  fh2pOverE2 = g.fh2pOverE2;
  fh2EOverp2 = g.fh2EOverp2;
  fh1dR     = g.fh1dR;
  fh2EledEdx = g.fh2EledEdx;
  fh2MatchdEdx = g.fh2MatchdEdx;
  fh2dEtadPhi = g.fh2dEtadPhi;
  fh2dEtadPhiMatched = g.fh2dEtadPhiMatched;
  fh2dEtadPhiUnmatched = g.fh2dEtadPhiUnmatched;
  fh2TrackPVsClusterE = g.fh2TrackPVsClusterE;
  fh2TrackPtVsClusterE = g.fh2TrackPtVsClusterE;
  fh2TrackPhiVsClusterPhi = g.fh2TrackPhiVsClusterPhi;
  fh2TrackEtaVsClusterEta = g.fh2TrackEtaVsClusterEta;
  //Photonic electron checks
  fh1OpeningAngle = g.fh1OpeningAngle;
  fh1MinvPhoton = g.fh1MinvPhoton;
  //Reconstructed electrons
  fhPtElectron = g.fhPtElectron; 
  fhPhiElectron = g.fhPhiElectron; 
  fhEtaElectron = g.fhEtaElectron; 
  fhPtNPE = g.fhPtNPE;
  fhPhiNPE = g.fhPhiNPE;
  fhEtaNPE = g.fhEtaNPE; 
  fhPtPE = g.fhPtPE;
  fhPhiPE = g.fhPhiPE;
  fhEtaPE = g.fhEtaPE; 
  //for comparisons with tracking detectors
  fhPtHadron = g.fhPtHadron; fhPtNPEleTPC = g.fhPtNPEleTPC; 
  fhPtNPEleTPCTRD = g.fhPtNPEleTPCTRD; fhPtNPEleTTE = g.fhPtNPEleTTE; 
  fhPtNPEleEMCAL = g.fhPtNPEleEMCAL; 
  //DVM B-tagging
  fhDVMBtagCut1 = g.fhDVMBtagCut1;
  fhDVMBtagCut2 = g.fhDVMBtagCut2; 
  fhDVMBtagCut3 = g.fhDVMBtagCut3; 
  fhDVMBtagQA1 = g.fhDVMBtagQA1; 
  fhDVMBtagQA2 = g.fhDVMBtagQA2; 
  fhDVMBtagQA3 = g.fhDVMBtagQA3; 
  fhDVMBtagQA4 = g.fhDVMBtagQA4; 
  fhDVMBtagQA5 = g.fhDVMBtagQA5; 
  //IPSig B-tagging
  fhIPSigBtagQA1 = g.fhIPSigBtagQA1; 
  fhIPSigBtagQA2 = g.fhIPSigBtagQA2; 
  fhTagJetPt1x4 = g.fhTagJetPt1x4; 
  fhTagJetPt2x3 = g.fhTagJetPt2x3; 
  fhTagJetPt3x2 = g.fhTagJetPt3x2; 
  fhePlusTagJetPt1x4 = g.fhePlusTagJetPt1x4; 
  fhePlusTagJetPt2x3 = g.fhePlusTagJetPt2x3; 
  fhePlusTagJetPt3x2 = g.fhePlusTagJetPt3x2; 
  //B-Jet histograms
  fhJetType = g.fhJetType; 
  fhLeadJetType = g.fhLeadJetType; 
  fhBJetXsiFF = g.fhBJetXsiFF; 
  fhBJetPtFF = g.fhBJetPtFF; 
  fhBJetEtaPhi = g.fhBJetEtaPhi; 
  fhNonBJetXsiFF = g.fhNonBJetXsiFF; 
  fhNonBJetPtFF = g.fhNonBJetPtFF; 
  fhNonBJetEtaPhi = g.fhNonBJetEtaPhi; 
  /////////////////////////////////////////////////////////////
  //Histograms that rely on MC info (not filled for real data)
  fEleNtuple = g.fEleNtuple; 
  //reco electrons from various sources
  fhPhiConversion = g.fhPhiConversion; 
  fhEtaConversion = g.fhEtaConversion;
  //for comparisons with tracking detectors
  fhPtTrack = g.fhPtTrack;
  fhPtNPEBHadron = g.fhPtNPEBHadron;
  //for computing efficiency of B-jet tags
  fhBJetPt1x4 = g.fhBJetPt1x4; fhBJetPt2x3 = g.fhBJetPt2x3; 
  fhBJetPt3x2 = g.fhBJetPt3x2;
  fhFakeJetPt1x4 = g.fhFakeJetPt1x4; fhFakeJetPt2x3 = g.fhFakeJetPt2x3; 
  fhFakeJetPt3x2 = g.fhFakeJetPt3x2; fhDVMJet = g.fhDVMJet;
  //MC rate histograms/ntuple
  fMCEleNtuple = g.fMCEleNtuple; fhMCBJetElePt = g.fhMCBJetElePt; 
  fhMCBHadronElePt = g.fhMCBHadronElePt;
  fhPtMCHadron = g.fhPtMCHadron; fhPtMCElectron = g.fhPtMCElectron; 
  fhMCXYConversion = g.fhMCXYConversion;
  fhMCRadPtConversion = g.fhMCRadPtConversion;

  return *this;
  
}
*/

//____________________________________________________________________________
AliAnaElectron::~AliAnaElectron() 
{
  //dtor

}

//________________________________________________________________________
TObjString *  AliAnaElectron::GetAnalysisCuts()
{
	//Save parameters used for analysis
	 TString parList ; //this will be list of parameters used for this analysis.
  const Int_t buffersize = 255; 
	 char onePar[buffersize] ;
	 
	 snprintf(onePar,buffersize,"--- AliAnaElectron ---\n") ;
	 parList+=onePar ;	
	 snprintf(onePar,buffersize,"fCalorimeter: %s\n",fCalorimeter.Data()) ;
	 parList+=onePar ;  
	 snprintf(onePar,buffersize,"fpOverEmin: %f\n",fpOverEmin) ;
	 parList+=onePar ;  
	 snprintf(onePar,buffersize,"fpOverEmax: %f\n",fpOverEmax) ;
	 parList+=onePar ;  
	 snprintf(onePar,buffersize,"fResidualCut: %f\n",fResidualCut) ;
	 parList+=onePar ;  
	 snprintf(onePar,buffersize,"fMinClusEne: %f\n",fMinClusEne) ;
	 parList+=onePar ;
	 snprintf(onePar,buffersize,"---DVM Btagging\n");
	 parList+=onePar ;
	 snprintf(onePar,buffersize,"max IP-cut (e,h): %f\n",fImpactCut);
	 parList+=onePar ;
	 snprintf(onePar,buffersize,"min ITS-hits: %d\n",fITSCut);
	 parList+=onePar ;
	 snprintf(onePar,buffersize,"max dR (e,h): %f\n",fDrCut);
	 parList+=onePar ;
	 snprintf(onePar,buffersize,"max pairDCA: %f\n",fPairDcaCut);
	 parList+=onePar ;
	 snprintf(onePar,buffersize,"max decaylength: %f\n",fDecayLenCut);
	 parList+=onePar ;
	 snprintf(onePar,buffersize,"min Associated Pt: %f\n",fAssocPtCut);
	 parList+=onePar ;
	 snprintf(onePar,buffersize,"---IPSig Btagging\n");
	 parList+=onePar ;
	 snprintf(onePar,buffersize,"min tag track: %d\n",fNTagTrkCut);
	 parList+=onePar ;
	 snprintf(onePar,buffersize,"min IP significance: %f\n",fIPSigCut);
	 parList+=onePar ;
	//
	 //Get parameters set in base class.
	 parList += GetBaseParametersList() ;
	 
	 //Get parameters set in FiducialCut class (not available yet)
	 //parlist += GetFidCut()->GetFidCutParametersList() 
	 
	 return new TObjString(parList) ;
  	
}

//________________________________________________________________________
TList *  AliAnaElectron::GetCreateOutputObjects()
{  
  // Create histograms to be saved in output file and 
  // store them in outputContainer
  TList * outputContainer = new TList() ; 
  outputContainer->SetName("ElectronHistos") ; 

  Int_t nptbins  = GetHistoPtBins();
  Int_t nphibins = GetHistoPhiBins();
  Int_t netabins = GetHistoEtaBins();
  Float_t ptmax  = GetHistoPtMax();
  Float_t phimax = GetHistoPhiMax();
  Float_t etamax = GetHistoEtaMax();
  Float_t ptmin  = GetHistoPtMin();
  Float_t phimin = GetHistoPhiMin();
  Float_t etamin = GetHistoEtaMin();	

  //event QA
  fhImpactXY = new TH1F("hImpactXY","Impact parameter for all tracks",200,-10,10.);
  fhRefMult = new TH1F("hRefMult" ,"refmult QA: " ,5000,0,5000);
  fhRefMult2  = new TH1F("hRefMult2" ,"refmult2 QA: " ,5000,0,5000);

  outputContainer->Add(fhImpactXY);
  outputContainer->Add(fhRefMult);
  outputContainer->Add(fhRefMult2);
  
  //matching checks
  fh3pOverE  = new TH3F("h3pOverE"    ,"EMCAL-TRACK matches p/E",nptbins,ptmin,ptmax,200,0.,10.,30,0,30);
  fh3EOverp  = new TH3F("h3EOverp"    ,"EMCAL-TRACK matches E/p",nptbins,ptmin,ptmax,200,0.,5. ,30,0,30);
  fh3pOverE2 = new TH3F("h3pOverE_Trk","EMCAL-TRACK matches p/E",nptbins,ptmin,ptmax,200,0.,10.,30,0,30);
  fh3EOverp2 = new TH3F("h3EOverp_Trk","EMCAL-TRACK matches E/p",nptbins,ptmin,ptmax,200,0.,5. ,30,0,30);
  fh3pOverE3 = new TH3F("h3pOverE_Tpc","EMCAL-TRACK matches p/E",nptbins,ptmin,ptmax,200,0.,10.,30,0,30);
  fh3EOverp3 = new TH3F("h3EOverp_Tpc","EMCAL-TRACK matches E/p",nptbins,ptmin,ptmax,200,0.,5. ,30,0,30);
  fh2pOverE  = new TH2F("h2pOverE"    ,"EMCAL-TRACK matches p/E",nptbins,ptmin,ptmax,200,0.,10.);
  fh2EOverp  = new TH2F("h2EOverp"    ,"EMCAL-TRACK matches E/p",nptbins,ptmin,ptmax,200,0.,5. );
  fh2pOverE2 = new TH2F("h2pOverE_Trk","EMCAL-TRACK matches p/E",nptbins,ptmin,ptmax,200,0.,10.);
  fh2EOverp2 = new TH2F("h2EOverp_Trk","EMCAL-TRACK matches E/p",nptbins,ptmin,ptmax,200,0.,5. );

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

  outputContainer->Add(fh3pOverE) ;
  outputContainer->Add(fh3EOverp) ;
  outputContainer->Add(fh3pOverE2) ;
  outputContainer->Add(fh3EOverp2) ;
  outputContainer->Add(fh3pOverE3) ;
  outputContainer->Add(fh3EOverp3) ;
  outputContainer->Add(fh2pOverE) ;
  outputContainer->Add(fh2EOverp) ;
  outputContainer->Add(fh2pOverE2) ;
  outputContainer->Add(fh2EOverp2) ;
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
  
  //photonic electron checks
  fh1OpeningAngle = new TH1F("hOpeningAngle","Opening angle between e+e- pairs",100,0.,TMath::Pi());
  fh1MinvPhoton = new TH1F("hMinvPhoton","Invariant mass of e+e- pairs",200,0.,2.);

  outputContainer->Add(fh1OpeningAngle);
  outputContainer->Add(fh1MinvPhoton);

  //Reconstructed electrons
  fhPtElectron = new TH1F("hPtElectron","Electron pT",nptbins,ptmin,ptmax);
  fhPhiElectron = new TH2F("hPhiElectron","Electron phi vs pT",nptbins,ptmin,ptmax,nphibins,phimin,phimax);
  fhEtaElectron = new TH2F("hEtaElectron","Electron eta vs. eta",nptbins,ptmin,ptmax,netabins,etamin,etamax);
  fhPtNPE = new TH1F("hPtNPE","Non-photonic Electron pT",nptbins,ptmin,ptmax);
  fhPhiNPE = new TH2F("hPhiNPE","Non-photonic Electron phi vs pT",nptbins,ptmin,ptmax,nphibins,phimin,phimax);
  fhEtaNPE = new TH2F("hEtaNPE","Non-photonic Electron eta vs. eta",nptbins,ptmin,ptmax,netabins,etamin,etamax);
  fhPtPE = new TH1F("hPtPE","Photonic Electron pT",nptbins,ptmin,ptmax);
  fhPhiPE = new TH2F("hPhiPE","Photonic Electron phi vs pT",nptbins,ptmin,ptmax,nphibins,phimin,phimax);
  fhEtaPE = new TH2F("hEtaPE","Photonic Electron eta vs. eta",nptbins,ptmin,ptmax,netabins,etamin,etamax);

  outputContainer->Add(fhPtElectron) ; 
  outputContainer->Add(fhPhiElectron) ; 
  outputContainer->Add(fhEtaElectron) ;
  outputContainer->Add(fhPtNPE) ; 
  outputContainer->Add(fhPhiNPE) ; 
  outputContainer->Add(fhEtaNPE) ;
  outputContainer->Add(fhPtPE) ; 
  outputContainer->Add(fhPhiPE) ; 
  outputContainer->Add(fhEtaPE) ;

  //These histograms are mixed REAL/MC:
  //Bins along y-axis are:  
  //0 - unfiltered (filled for both real and MC data) 
  //1 - bottom, 2 - charm, 3 - charm from bottom  (MC only)
  //4 - conversion, 5 - Dalitz, 6 - W and Z, 7 - junk/unknown (MC only)
  //8 - misidentified (MC only)

  //histograms for comparison to tracking detectors
  fhPtHadron = new TH2F("hPtHadron","Charged hadrons w/in EMCAL acceptance",nptbins,ptmin,ptmax,10,0,10);
  fhPtNPEleTPC = new TH2F("hPtNPEleTPC","Non-phot. Electrons identified by TPC w/in EMCAL acceptance",nptbins,ptmin,ptmax,10,0,10);
  fhPtNPEleTPCTRD = new TH2F("hPtNPEleTPCTRD","Non-phot. Electrons identified by TPC+TRD w/in EMCAL acceptance",nptbins,ptmin,ptmax,10,0,10);
  fhPtNPEleTTE = new TH2F("hPtNPEleTTE","Non-phot. Electrons identified by TPC+TRD+EMCAL w/in EMCAL acceptance",nptbins,ptmin,ptmax,10,0,10);    
  fhPtNPEleEMCAL = new TH2F("hPtNPEleEMCAL","Non-phot. Electrons identified by EMCAL w/in EMCAL acceptance",nptbins,ptmin,ptmax,10,0,10);
  
  outputContainer->Add(fhPtHadron);
  outputContainer->Add(fhPtNPEleTPC);
  outputContainer->Add(fhPtNPEleTPCTRD);
  outputContainer->Add(fhPtNPEleTTE);
  outputContainer->Add(fhPtNPEleEMCAL);

  //B-tagging
  fhDVMBtagCut1 = new TH2F("hdvmbtag_cut1","DVM B-tag result cut1", 10,0,10 ,nptbins,ptmin,ptmax);
  fhDVMBtagCut2 = new TH2F("hdvmbtag_cut2","DVM B-tag result cut2", 10,0,10 ,nptbins,ptmin,ptmax);
  fhDVMBtagCut3 = new TH2F("hdvmbtag_cut3","DVM B-tag result cut3", 10,0,10 ,nptbins,ptmin,ptmax);
  fhDVMBtagQA1  = new TH2F("hdvmbtag_qa1" ,"DVM B-tag QA: pairDCA vs length", 100,0,0.2 ,100,0,1.0);
  fhDVMBtagQA2  = new TH2F("hdvmbtag_qa2" ,"DVM B-tag QA: signDCA vs mass"  , 200,-0.5,0.5 ,100,0,10);
  fhDVMBtagQA3  = new TH1F("hdvmbtag_qa3" ,"DVM B-tag QA: ITS-Hits electron" ,7,0,7);
  fhDVMBtagQA4  = new TH1F("hdvmbtag_qa4" ,"DVM B-tag QA: IP d electron" ,200,-3,3);
  fhDVMBtagQA5  = new TH1F("hdvmbtag_qa5" ,"DVM B-tag QA: IP z electron" ,200,-3,3);

  outputContainer->Add(fhDVMBtagCut1) ;
  outputContainer->Add(fhDVMBtagCut2) ;
  outputContainer->Add(fhDVMBtagCut3) ;
  outputContainer->Add(fhDVMBtagQA1) ;
  outputContainer->Add(fhDVMBtagQA2) ;
  outputContainer->Add(fhDVMBtagQA3) ;
  outputContainer->Add(fhDVMBtagQA4) ;
  outputContainer->Add(fhDVMBtagQA5) ;

  //IPSig B-tagging
  fhIPSigBtagQA1  = new TH1F("hipsigbtag_qa1" ,"IPSig B-tag QA: # tag tracks", 20,0,20);
  fhIPSigBtagQA2  = new TH1F("hipsigbtag_qa2" ,"IPSig B-tag QA: IP significance", 200,-10.,10.);
  fhTagJetPt1x4 = new TH1F("hTagJetPt1x4","tagged jet pT (1 track, ipSignif>4);p_{T}",300,0.,300.);
  fhTagJetPt2x3 = new TH1F("hTagJetPt2x3","tagged jet pT (2 track, ipSignif>3);p_{T}",300,0.,300.);
  fhTagJetPt3x2 = new TH1F("hTagJetPt3x2","tagged jet pT (3 track, ipSignif>2);p_{T}",300,0.,300.);
  fhePlusTagJetPt1x4 = new TH1F("hePlusTagJetPt1x4","tagged eJet pT (1 track, ipSignif>4);p_{T}",300,0.,300.);
  fhePlusTagJetPt2x3 = new TH1F("hePlusTagJetPt2x3","tagged eJet pT (2 track, ipSignif>3);p_{T}",300,0.,300.);
  fhePlusTagJetPt3x2 = new TH1F("hePlusTagJetPt3x2","tagged eJet pT (3 track, ipSignif>2);p_{T}",300,0.,300.);

  outputContainer->Add(fhIPSigBtagQA1) ;
  outputContainer->Add(fhIPSigBtagQA2) ;
  outputContainer->Add(fhTagJetPt1x4);
  outputContainer->Add(fhTagJetPt2x3);
  outputContainer->Add(fhTagJetPt3x2);
  outputContainer->Add(fhePlusTagJetPt1x4);
  outputContainer->Add(fhePlusTagJetPt2x3);
  outputContainer->Add(fhePlusTagJetPt3x2);

  //B-Jet histograms
  fhJetType = new TH2F("hJetType","# jets passing each tag method vs jet pt",15,0,15,300,0.,300.);
  fhLeadJetType = new TH2F("hLeadJetType","# leading jets passing each tag method vs jet pt",15,0,15,300,0.,300.);
  fhBJetXsiFF = new TH2F("hBJetXsiFF","B-jet #Xsi Frag. Fn.",100,0.,10.,300,0.,300.);
  fhBJetPtFF = new TH2F("hBJetPtFF","B-jet p_{T} Frag. Fn.",nptbins,ptmin,ptmax,300,0.,300.);
  fhBJetEtaPhi = new TH2F("hBJetEtaPhi","B-jet eta-phi distribution",netabins,etamin,etamax,nphibins,phimin,phimax);
  fhNonBJetXsiFF = new TH2F("hNonBJetXsiFF","Non B-jet #Xsi Frag. Fn.",100,0.,10.,300,0.,300.);
  fhNonBJetPtFF = new TH2F("hNonBJetPtFF","Non B-jet p_{T} Frag. Fn.",nptbins,ptmin,ptmax,300,0.,300.);
  fhNonBJetEtaPhi = new TH2F("hNonBJetEtaPhi","Non B-jet eta-phi distribution",netabins,etamin,etamax,nphibins,phimin,phimax);

  outputContainer->Add(fhJetType);
  outputContainer->Add(fhLeadJetType);
  outputContainer->Add(fhBJetXsiFF);
  outputContainer->Add(fhBJetPtFF);
  outputContainer->Add(fhBJetEtaPhi);
  outputContainer->Add(fhNonBJetXsiFF);
  outputContainer->Add(fhNonBJetPtFF);
  outputContainer->Add(fhNonBJetEtaPhi);

  //Histograms that use MC information
  if(IsDataMC()){

    //electron ntuple for further analysis
    if(fWriteNtuple) {
      fEleNtuple = new TNtuple("EleNtuple","Electron Ntuple","tmctag:cmctag:pt:phi:eta:p:E:deta:dphi:nCells:dEdx:pidProb:impXY:impZ");
      outputContainer->Add(fEleNtuple) ;
    }

    //electrons from various MC sources
    fhPhiConversion = new TH2F("hPhiConversion","Conversion Electron phi vs pT",nptbins,ptmin,ptmax,nphibins,phimin,phimax);
    fhEtaConversion = new TH2F("hEtaConversion","Conversion Electron eta vs. pT",nptbins,ptmin,ptmax,netabins,etamin,etamax);

    outputContainer->Add(fhPhiConversion);
    outputContainer->Add(fhEtaConversion);

    //Bins along y-axis are:  0 - unfiltered, 1 - bottom, 2 - charm, 3 - charm from bottom,
    //4 - conversion, 5 - Dalitz, 6 - W and Z, 7 - junk/unknown, 8 - misidentified

    //histograms for comparison to tracking detectors
    fhPtTrack  = new TH2F("hPtTrack","Track w/in EMCAL acceptance",nptbins,ptmin,ptmax,10,0,10);
    fhPtNPEBHadron = new TH2F("hPtNPEBHadron","Non-phot. b-electrons (TPC+TRD+EMCAL) vs B-hadron pt w/in EMCAL acceptance",nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);

    outputContainer->Add(fhPtTrack);
    outputContainer->Add(fhPtNPEBHadron);

    //for computing efficiency of IPSig tag
    fhBJetPt1x4 = new TH1F("hBJetPt1x4","tagged B-jet pT (1 track, ipSignif>4);p_{T}",300,0.,300.);
    fhBJetPt2x3 = new TH1F("hBJetPt2x3","tagged B-jet pT (2 track, ipSignif>3);p_{T}",300,0.,300.);
    fhBJetPt3x2 = new TH1F("hBJetPt3x2","tagged B-jet pT (3 track, ipSignif>2);p_{T}",300,0.,300.);
    fhFakeJetPt1x4 = new TH1F("hFakeJetPt1x4","fake tagged B-jet pT (1 track, ipSignif>4);p_{T}",300,0.,300.);
    fhFakeJetPt2x3 = new TH1F("hFakeJetPt2x3","fake tagged B-jet pT (2 track, ipSignif>3);p_{T}",300,0.,300.);
    fhFakeJetPt3x2 = new TH1F("hFakeJetPt3x2","fake tagged B-jet pT (3 track, ipSignif>2);p_{T}",300,0.,300.);
    fhDVMJet = new TH2F("hDVM_algo","# DVM jets passing vs Mc-Bjet",10,0,10,300,0.,300.);

    outputContainer->Add(fhBJetPt1x4);
    outputContainer->Add(fhBJetPt2x3);
    outputContainer->Add(fhBJetPt3x2);
    outputContainer->Add(fhFakeJetPt1x4);
    outputContainer->Add(fhFakeJetPt2x3);
    outputContainer->Add(fhFakeJetPt3x2);
    outputContainer->Add(fhDVMJet);

    //MC Only histograms
    
    //MC ele ntuple for further analysis
    if(fWriteNtuple) {
      fMCEleNtuple = new TNtuple("MCEleNtuple","MC Electron Ntuple","mctag:pt:phi:eta:x:y:z");
      outputContainer->Add(fMCEleNtuple) ;
    }

    fhMCBJetElePt = new TH2F("hMCBJetElePt","MC B-jet pT vs. electron pT",300,0.,300.,300,0.,300.);
    fhMCBHadronElePt = new TH2F("hMCBHadronElePt","MC B-hadron pT vs. electron pT",nptbins,ptmin,ptmax,nptbins,ptmin,ptmax);
    fhPtMCHadron = new TH1F("hPtMCHadron","MC Charged hadrons w/in EMCAL acceptance",nptbins,ptmin,ptmax);

    //Bins along y-axis are:  0 - unfiltered, 1 - bottom, 2 - charm, 3 - charm from bottom,
    //4 - conversion, 5 - Dalitz, 6 - W and Z, 7 - junk/unknown
    fhPtMCElectron = new TH2F("hPtMCElectron","MC electrons from various sources w/in EMCAL acceptance",nptbins,ptmin,ptmax,10,0,10);

    fhMCXYConversion = new TH2F("hMCXYConversion","XvsY of conversion electrons",400,-400.,400.,400,-400.,400.);
    fhMCRadPtConversion = new TH2F("hMCRadPtConversion","Radius vs pT of conversion electrons",200,0.,400.,nptbins,ptmin,ptmax);

    outputContainer->Add(fhMCBJetElePt);
    outputContainer->Add(fhMCBHadronElePt);
    outputContainer->Add(fhPtMCHadron);
    outputContainer->Add(fhPtMCElectron);
    outputContainer->Add(fhMCXYConversion);
    outputContainer->Add(fhMCRadPtConversion);

  }//Histos with MC
  

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
  fpOverEmax = 1.2;
  fResidualCut = 0.02;
  fMinClusEne = 4.0;
  //DVM B-tagging
  fDrCut       = 1.0; 
  fPairDcaCut  = 0.02;
  fDecayLenCut = 1.0;
  fImpactCut   = 0.5;
  fAssocPtCut  = 1.0;
  fMassCut     = 1.5;
  fSdcaCut     = 0.1;
  fITSCut      = 4;
  //IPSig B-tagging
  fNTagTrkCut  = 2;
  fIPSigCut    = 3.0;
  //Jet fiducial cuts
  fJetEtaCut = 0.3;
  fJetPhiMin = 1.8;
  fJetPhiMax = 2.9;
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

  Double_t bfield = 0.;
  if(GetReader()->GetDataType() != AliCaloTrackReader::kMC) bfield = GetReader()->GetBField();

  //Select the calorimeter of the electron
  if(fCalorimeter != "EMCAL") {
    printf("This class not yet implemented for PHOS\n");
    abort();
  }
  
  TObjArray *cl = GetEMCALClusters();
  
  ////////////////////////////////////////////////
  //Start from tracks and get associated clusters 
  ////////////////////////////////////////////////
  if(!GetCTSTracks() || GetCTSTracks()->GetEntriesFast() == 0) return ;
  Int_t ntracks = GetCTSTracks()->GetEntriesFast();
  Int_t refmult = 0; Int_t refmult2 = 0;
  if(GetDebug() > 0)
    printf("AliAnaElectron::MakeAnalysisFillAOD() - In CTS aod entries %d\n", ntracks);

  //Unfortunately, AliAODTracks don't have associated EMCAL clusters.
  //we have to redo track-matching, I guess
  Int_t iCluster = -999;
  Int_t bt = 0; //counter for event b-tags

  for (Int_t itrk =  0; itrk <  ntracks; itrk++) {////////////// track loop
    iCluster = -999; //start with no match
    AliAODTrack * track = (AliAODTrack*) (GetCTSTracks()->At(itrk)) ;
    //Added negative track condition. Seems that about 3-4% of the tracks have negative label. Causes crashes sometimes.
    if(track->GetLabel()<0){
      printf("Negative track label! Not sure what it  means, abort track. \n");
      continue;
    }

    if (TMath::Abs(track->Eta())< 0.5) refmult++;
    Double_t imp[2] = {-999.,-999.}; Double_t cov[3] = {-999.,-999.,-999.};
    Bool_t dcaOkay = GetDCA(track,imp,cov);  //homegrown dca calculation until AOD is fixed
    if(!dcaOkay) printf("AliAnaElectron::Problem computing DCA to primary vertex for track %d.  Skipping it...\n",itrk);
    if(TMath::Abs(track->Eta())< 0.5 && TMath::Abs(imp[0])<1.0 && TMath::Abs(imp[1])<1.0) refmult2++;
    fhImpactXY->Fill(imp[0]);

    //JLK CHECK
    //AliESDtrack esdTrack(track);
    //Double_t tpcpid[AliPID::kSPECIES];
    //esdTrack.GetTPCpid(tpcpid);
    //Double_t eProb = tpcpid[AliPID::kElectron];
    //if(eProb > 0) printf("<%d> ESD eProb = %2.2f\n",itrk,eProb);

    AliAODPid* pid = (AliAODPid*) track->GetDetPid();
    if(pid == 0) {
      if(GetDebug() > 0) printf("AliAnaElectron::MakeAnalysisFillAOD() - No PID object - skipping track %d",itrk);
      continue;
    } else {
      Double_t emcpos[3];
      pid->GetEMCALPosition(emcpos);
      Double_t emcmom[3];
      pid->GetEMCALMomentum(emcmom);
      
      TVector3 pos(emcpos[0],emcpos[1],emcpos[2]);
      TVector3 mom(emcmom[0],emcmom[1],emcmom[2]);
      Double_t tphi = pos.Phi();
      Double_t teta = pos.Eta();
      Double_t tmom = mom.Mag();
      
      //TLorentzVector mom2(mom,0.);
      Bool_t in = kFALSE;
      //Removed two lines below. Seems wrong to demand that the momentum vector should point to the EMCAL.
      //if(mom.Phi()*180./TMath::Pi() > 80. && mom.Phi()*180./TMath::Pi() < 190. &&
      // mom.Eta() > -0.7 && mom.Eta() < 0.7) in = kTRUE;
      //Also check the track
      if(track->Phi()*180./TMath::Pi() > 80. && track->Phi()*180./TMath::Pi() < 190. &&
	 track->Eta() > -0.7 && track->Eta() < 0.7) in = kTRUE;
      ////////////////////////////
      //THIS HAS A MEM LEAK JLK 24-Oct-09
      //Bool_t in =  GetFiducialCut()->IsInFiducialCut(mom2,fCalorimeter) ;
      ///////////////////////////
      if(GetDebug() > 1) printf("AliAnaElectron::MakeAnalysisFillAOD() - Track(Extrap) pt %2.2f(%2.2f), phi %2.2f(%2.2f), eta %2.2f(%2.2f) in fiducial cut %d\n",track->Pt(), mom.Pt(), track->Phi(), mom.Phi(), track->Eta(),mom.Eta(), in);

      if(mom.Pt() > GetMinPt() && in) {
	
	Double_t dEdx = pid->GetTPCsignal();
		
	Int_t ntot = cl->GetEntriesFast();
	Double_t res = 999.;
	Double_t pOverE = -999.;
	
	Int_t pidProb = track->GetMostProbablePID();
	Bool_t tpcEle = kFALSE; if(dEdx > 70.) tpcEle = kTRUE;
	Bool_t trkEle = kFALSE; if(pidProb == AliAODTrack::kElectron) trkEle = kTRUE;
	Bool_t trkChgHad = kFALSE; if(pidProb == AliAODTrack::kPion || pidProb == AliAODTrack::kKaon || pidProb == AliAODTrack::kProton) trkChgHad = kTRUE;

	Int_t tmctag = -1;

	//Check against V0 for conversion, only if it is flagged as electron
	Bool_t photonic = kFALSE;
	if(tpcEle || trkEle) photonic = PhotonicV0(itrk);
	if(trkEle && !photonic) fhPtNPEleTPCTRD->Fill(track->Pt(),0); //0 = no MC info
	if(tpcEle && !photonic) fhPtNPEleTPC->Fill(track->Pt(),0); //0 = no MC info

	if(trkChgHad) fhPtHadron->Fill(track->Pt(),0); //0 = no MC info
	if(IsDataMC()) {
	  //Input from second AOD?
	  Int_t input = 0;
	  //if(GetReader()->GetCTSTracksNormalInputEntries() <= itrk) input = 1;
	  tmctag = GetMCAnalysisUtils()->CheckOrigin(track->GetLabel(),GetReader(),input);

	  if(trkChgHad) fhPtHadron->Fill(track->Pt(),GetMCSource(tmctag));
	  if(tpcEle && !photonic) fhPtNPEleTPC->Fill(track->Pt(),GetMCSource(tmctag));
	  if(trkEle && !photonic) fhPtNPEleTPCTRD->Fill(track->Pt(),GetMCSource(tmctag));
	  fhPtTrack->Fill(track->Pt(),GetMCSource(tmctag));
	}

	Bool_t emcEle = kFALSE;      
	//For tracks in EMCAL acceptance, pair them with all clusters
	//and fill the dEta vs dPhi for these pairs:

	Double_t minR  = 99;
        Double_t minPe =-1;
        Double_t minEp =-1;
        Double_t minMult = -1;
        Double_t minPt   = -1;

	for(Int_t iclus = 0; iclus < ntot; iclus++) {
	  AliVCluster * clus = (AliVCluster*) (cl->At(iclus));
	  if(!clus) continue;

	  //As of 11-Oct-2009
	  //only select "good" clusters	  
	  if (clus->GetNCells()       < 2    ) continue;
          if (clus->GetNCells()       > 30   ) continue;
          if (clus->E()               < fMinClusEne ) continue;
          if (clus->GetDispersion()   > 1    ) continue;
          if (clus->GetM20()          > 0.4  ) continue;
          if (clus->GetM02()          > 0.4  ) continue;
          if (clus->GetM20()          < 0.03 ) continue;
          if (clus->GetM02()          < 0.03 ) continue;
	  
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
	  
	  /////////////////////////////////
	  //Perform electron cut analysis//
	  /////////////////////////////////
	  //Good match
	  if(res < fResidualCut) {
	    fh2dEtadPhiMatched->Fill(deta,dphi);
	    fh2MatchdEdx->Fill(track->P(),dEdx);
	    iCluster = iclus;

	    Double_t energy = clus->E(); 
	    if(energy > 0) pOverE = tmom/energy;

	    if (res< minR) {
              minR  = res;
              minPe = pOverE;
              minEp = energy/tmom;
              minMult = clus->GetNCells() ;
              minPt = track->Pt();
            }

	    Int_t cmctag = -1;	    
	    if(IsDataMC()) {  
	      //Do you want the cluster or the track label?
	      Int_t input = 0;
	      //if(GetReader()->GetEMCALClustersNormalInputEntries() <= iclus) input = 1;
	      cmctag = GetMCAnalysisUtils()->CheckOrigin(clus->GetLabel(),GetReader(),input);
	    }
	    
	    if(fWriteNtuple) {
	      fEleNtuple->Fill(tmctag,cmctag,track->Pt(),track->Phi(),track->Eta(),track->P(),clus->E(),deta,dphi,clus->GetNCells(),dEdx,pidProb,imp[0],imp[1]);
	    }

	  } else {
	    fh2dEtadPhiUnmatched->Fill(deta,dphi);
	  }//res cut
	}//calo cluster loop

	fh3pOverE->Fill(minPt,minPe ,minMult);
	fh3EOverp->Fill(minPt,minEp ,minMult);
	if (trkEle) {
	  fh3pOverE2->Fill(minPt,minPe ,minMult);
	  fh3EOverp2->Fill(minPt,minEp ,minMult);
	}
	if (tpcEle) {
	  fh3pOverE3->Fill(minPt,minPe ,minMult);
	  fh3EOverp3->Fill(minPt,minEp ,minMult);
	}
	//new
	if (tmctag>-1 && GetMCSource(tmctag)<8 ) {
          fh2pOverE->Fill(minPt,minPe );
          fh2EOverp->Fill(minPt,minEp );
          if (trkEle) {
            fh2pOverE2->Fill(minPt,minPe );
            fh2EOverp2->Fill(minPt,minEp );
          }
        }

	//////////////////////////////
	//Electron cuts happen here!//
	//////////////////////////////
	if(minPe > fpOverEmin && minPe < fpOverEmax) emcEle = kTRUE;

	///////////////////////////
	//Fill AOD with electrons//
	///////////////////////////
	//Take all emcal electrons, but the others only if pT < 10 GeV
	if(emcEle || ( (tpcEle || trkEle) && track->Pt() < 10.) ) {

	  //B-tagging
	  if(GetDebug() > 1) printf("Found Electron - do b-tagging\n");
	  Int_t dvmbtag = GetDVMBtag(track); bt += dvmbtag;
	  fh2EledEdx->Fill(track->P(),dEdx);
	  
	  Double_t eMass = 0.511/1000; //mass in GeV
	  Double_t eleE = sqrt(track->P()*track->P() + eMass*eMass);
	  AliAODPWG4Particle tr = AliAODPWG4Particle(track->Px(),track->Py(),track->Pz(),eleE);
	  tr.SetLabel(track->GetLabel());
	  tr.SetCaloLabel(iCluster,-1); //sets the indices of the original caloclusters
	  tr.SetTrackLabel(itrk,-1); //sets the indices of the original tracks
	  if(emcEle) {//PID determined by EMCAL
	    tr.SetDetector(fCalorimeter);
	  } else {
	    tr.SetDetector("CTS"); //PID determined by CTS
	  }

	  //if(GetReader()->GetCTSTracksNormalInputEntries() <= itrk) tr.SetInputFileIndex(1);
	  //Make this preserve sign of particle
	  if(track->Charge() < 0) tr.SetPdg(11); //electron is 11
	  else  tr.SetPdg(-11); //positron is -11
	  Int_t btag = 0;
	  if(dvmbtag > 0) tr.SetBTagBit(btag,tr.kDVMTag0);
	  if(dvmbtag > 1) tr.SetBTagBit(btag,tr.kDVMTag1);
	  if(dvmbtag > 2) tr.SetBTagBit(btag,tr.kDVMTag2);
	  tr.SetBtag(btag);
	  
	  //Play with the MC stack if available
	  //Check origin of the candidates
	  if(IsDataMC()){
	    
	    //FIXME:  Need to re-think this for track-oriented analysis
	    //JLK DO WE WANT TRACK TAG OR CLUSTER TAG?
	    tr.SetTag(GetMCAnalysisUtils()->CheckOrigin(tr.GetLabel(),GetReader(),tr.GetInputFileIndex()));
	    
	    if(GetDebug() > 0) printf("AliAnaElectron::MakeAnalysisFillAOD() - Origin of candidate %d\n",tr.GetTag());
	  }//Work with stack also   
	  
	  AddAODParticle(tr);
	  
	  if(GetDebug() > 1) printf("AliAnaElectron::MakeAnalysisFillAOD() - Electron selection cuts passed: pT %3.2f, pdg %d\n",tr.Pt(),tr.GetPdg());	
	}//electron
      }//pt, fiducial selection
    }//pid check
  }//track loop                         
  
  fhRefMult->Fill(refmult);
  fhRefMult2->Fill(refmult2);

  if(GetDebug() > 1 && bt > 0) printf("AliAnaElectron::MakeAnalysisFillAOD() *** Event Btagged *** \n");
  if(GetDebug() > 1) printf("AliAnaElectron::MakeAnalysisFillAOD()  End fill AODs \n");  
  
}

//__________________________________________________________________
void  AliAnaElectron::MakeAnalysisFillHistograms() 
{
  //Do analysis and fill histograms
  
  AliStack * stack = 0x0;
  TParticle * primary = 0x0;
  AliAODMCParticle * aodprimary = 0x0;
  
  Int_t ph1 = 0;  //photonic 1 count
  Int_t ph2 = 0;  //photonic 2 count
  Int_t phB = 0;  //both count
  
  if(IsDataMC()) {
    if(GetReader()->ReadStack()){
      stack =  GetMCStack() ;      
      if(!stack){
        printf("AliAnaElectron::MakeAnalysisFillHistograms() *** no stack ***: \n");
      }
      
    }
  }// is data and MC
  
  ////////////////////////////////////
  //Loop over jets and check for b-tag
  ////////////////////////////////////
  Int_t njets = (GetReader()->GetOutputEvent())->GetNJets();
  if(njets > 0) {
    if(GetDebug() > 0) printf("AliAnaElectron::MakeAnalysisFillHistograms() - Jet AOD branch has %d jets.  Performing b-jet tag analysis\n",njets);
    
    for(Int_t ijet = 0; ijet < njets ; ijet++) {
      AliAODJet * jet = (AliAODJet*)(GetReader()->GetOutputEvent())->GetJet(ijet) ;
      //Only consider jets with pt > 10 GeV (the rest have to be junk)
      //printf("AODJet<%d> pt = %2.2f\n",ijet,jet->Pt());
      if(jet->Pt() < 10.) continue;
      
      if(GetDebug() > 3) {
        printf("AliAODJet ijet = %d\n",ijet);
        jet->Print("");
      }
      //Skip jets not inside a smaller fiducial volume to ensure that
      //they are completely contained in the EMCAL
      if(TMath::Abs(jet->Eta()) > fJetEtaCut) continue;
      if(jet->Phi() < fJetPhiMin || jet->Phi() > fJetPhiMax) continue;
      
      //To "tag" the jet, we will look for it to pass our various criteria
      //For e jet tag, we just look to see which ones have NPEs
      //For DVM jet tag, we will look for DVM electrons
      //For IPSig, we compute the IPSig for all tracks and if the
      //number passing is above the cut, it passes
      Bool_t leadJet  = kFALSE;
      if (ijet==0) leadJet= kTRUE;
      Bool_t eJet = kFALSE;  
      Bool_t eJet2    = kFALSE; //electron triggered
      Bool_t hadJet   = kFALSE; //hadron triggered 
      Bool_t dvmJet = kFALSE;  
      Bool_t ipsigJet = kFALSE;
      TRefArray* rt = jet->GetRefTracks();
      Int_t ntrk = rt->GetEntries();
      Int_t trackCounter[4] = {0,0,0,0}; //for ipsig
      for(Int_t itrk = 0; itrk < ntrk; itrk++) {
      	AliAODTrack* jetTrack = (AliAODTrack*)jet->GetTrack(itrk);
        if( GetIPSignificance(jetTrack, jet->Phi()) > fIPSigCut) trackCounter[0]++;
        if( GetIPSignificance(jetTrack, jet->Phi()) > 4.) trackCounter[1]++;
        if( GetIPSignificance(jetTrack, jet->Phi()) > 3.) trackCounter[2]++;
        if( GetIPSignificance(jetTrack, jet->Phi()) > 2.) trackCounter[3]++;
        Bool_t isNPE = CheckTrack(jetTrack,"NPE");
        if(isNPE) eJet = kTRUE;
        if ( isNPE && jetTrack->Pt()>10.0 ) eJet2  =kTRUE;
        if (!isNPE && jetTrack->Pt()>10.0) hadJet =kTRUE;
        Bool_t isDVM = CheckTrack(jetTrack,"DVM");
        if(isDVM) dvmJet = kTRUE;
      }
      fhIPSigBtagQA1->Fill(trackCounter[0]);
      if(trackCounter[1]>0) fhTagJetPt1x4->Fill(jet->Pt());
      if(trackCounter[2]>1) fhTagJetPt2x3->Fill(jet->Pt());
      if(trackCounter[3]>2) fhTagJetPt3x2->Fill(jet->Pt());
      
      if(trackCounter[1]>0 && eJet) fhePlusTagJetPt1x4->Fill(jet->Pt());
      if(trackCounter[2]>1 && eJet) fhePlusTagJetPt2x3->Fill(jet->Pt());
      if(trackCounter[3]>2 && eJet) fhePlusTagJetPt3x2->Fill(jet->Pt());
      
      if(trackCounter[0] > fNTagTrkCut) ipsigJet = kTRUE;
      
      if(IsDataMC()) {
        //determine tagging efficiency & mis-tagging rate
        //using b-quarks from stack
        Bool_t isTrueBjet = IsMcBJet(jet->Eta(), jet->Phi());
        Bool_t isTrueDjet = IsMcDJet(jet->Eta(), jet->Phi());
        if (isTrueBjet && GetDebug() > 0) printf("== True Bjet==\n");
        if (isTrueDjet && GetDebug() > 0) printf("== True Charm-jet==\n");
        if (dvmJet && GetDebug() > 0)     printf("== found DVM jet==\n");
        
        if(isTrueBjet && dvmJet) fhDVMJet->Fill(0.,jet->Pt()); // good tagged
        if(isTrueBjet && !dvmJet) fhDVMJet->Fill(1.,jet->Pt()); // missed tagged
        if(!isTrueBjet && dvmJet) fhDVMJet->Fill(2.,jet->Pt());  // fake tagged
        if(!isTrueBjet && !dvmJet) fhDVMJet->Fill(3.,jet->Pt());  // others
        if(isTrueDjet && !isTrueBjet &&   dvmJet) fhDVMJet->Fill(4.,jet->Pt()); // charm-tagged
        if(isTrueDjet && !isTrueBjet &&  !dvmJet) fhDVMJet->Fill(5.,jet->Pt()); // charm -not tagged
        if(!(isTrueDjet||isTrueBjet ) &&  dvmJet) fhDVMJet->Fill(6.,jet->Pt()); // light flavor -tagged
        if(!(isTrueDjet||isTrueBjet ) && !dvmJet) fhDVMJet->Fill(7.,jet->Pt()); // light flavor -not tagged
        if(isTrueBjet &&  eJet && dvmJet) fhDVMJet->Fill(8.,jet->Pt()); // bjet with electron
        if(isTrueBjet && !eJet && dvmJet) fhDVMJet->Fill(9.,jet->Pt()); // needs more thought
        
        if(isTrueBjet) {
          if(trackCounter[1]>0) fhBJetPt1x4->Fill(jet->Pt());
          if(trackCounter[2]>1) fhBJetPt2x3->Fill(jet->Pt());
          if(trackCounter[3]>2) fhBJetPt3x2->Fill(jet->Pt());
        } else {
          if(trackCounter[1]>0) fhFakeJetPt1x4->Fill(jet->Pt());
          if(trackCounter[2]>1) fhFakeJetPt2x3->Fill(jet->Pt());
          if(trackCounter[3]>2) fhFakeJetPt3x2->Fill(jet->Pt());
        }
      }
      
      //Fill bjet histograms here
      if(!(eJet || ipsigJet || dvmJet)) fhJetType->Fill(0.,jet->Pt()); //none
      if(eJet && !(ipsigJet || dvmJet)) fhJetType->Fill(1.,jet->Pt()); //only ejet
      if(dvmJet && !(eJet || ipsigJet)) fhJetType->Fill(2.,jet->Pt()); //only dvm
      if(ipsigJet && !(eJet || dvmJet)) fhJetType->Fill(3.,jet->Pt()); //only ipsig
      if(eJet && dvmJet && !ipsigJet)   fhJetType->Fill(4.,jet->Pt()); //ejet & dvm
      if(eJet && ipsigJet && !dvmJet)   fhJetType->Fill(5.,jet->Pt()); //ejet & ipsig
      if(dvmJet && ipsigJet && !eJet)   fhJetType->Fill(6.,jet->Pt()); //dvm & ipsig
      if(dvmJet && ipsigJet && eJet)    fhJetType->Fill(7.,jet->Pt()); //all
      if(dvmJet || ipsigJet || eJet)    fhJetType->Fill(8.,jet->Pt()); //any of them
      if(eJet  )    fhJetType->Fill(9.,jet->Pt()); //any of them
      if(dvmJet)    fhJetType->Fill(10.,jet->Pt()); //any of them
      if(eJet2 )    fhJetType->Fill(11.,jet->Pt()); //any of them
      if(hadJet)    fhJetType->Fill(12.,jet->Pt()); //any of them
      
      if(eJet || ipsigJet || dvmJet) fhBJetEtaPhi->Fill(jet->Eta(),jet->Phi());
      else fhNonBJetEtaPhi->Fill(jet->Eta(),jet->Phi());
      
      //leading jets
      if (leadJet) {
        fhLeadJetType->Fill(0.,jet->Pt()); //all
        if(eJet   )    fhLeadJetType->Fill(1.,jet->Pt());
        if(eJet2  )    fhLeadJetType->Fill(2.,jet->Pt());
        if(hadJet )    fhLeadJetType->Fill(3.,jet->Pt());
        if(eJet   && (dvmJet || ipsigJet) )    fhLeadJetType->Fill(4.,jet->Pt());
        if(eJet2  && (dvmJet || ipsigJet) )    fhLeadJetType->Fill(5.,jet->Pt());
        if(hadJet && (dvmJet || ipsigJet) )    fhLeadJetType->Fill(6.,jet->Pt());
      }
      
      for(Int_t itrk = 0; itrk < ntrk; itrk++) {
        AliAODTrack* jetTrack = (AliAODTrack*)jet->GetTrack(itrk);
        Double_t xsi = TMath::Log(jet->Pt()/jetTrack->Pt());
        if(eJet || ipsigJet || dvmJet) {
          if(GetDebug() > 0) printf("AliAnaElectron::MakeAnalysisFillHistograms - We have a bjet!\n");
          fhBJetXsiFF->Fill(xsi,jet->Pt());
          fhBJetPtFF->Fill(jetTrack->Pt(),jet->Pt());
        } else {
          //Fill non-bjet histograms here
          fhNonBJetXsiFF->Fill(xsi,jet->Pt());
          fhNonBJetPtFF->Fill(jetTrack->Pt(),jet->Pt());
        }
      }
      
    } //jet loop
  } //jets exist
  
  //////////////////////////////
  //Loop on stored AOD electrons
  //////////////////////////////
  Int_t naod = GetOutputAODBranch()->GetEntriesFast();
  if(GetDebug() > 0) printf("AliAnaElectron::MakeAnalysisFillHistograms() - aod branch entries %d\n", naod);
  
  for(Int_t iaod = 0; iaod < naod ; iaod++){
    AliAODPWG4Particle* ele =  (AliAODPWG4Particle*) (GetOutputAODBranch()->At(iaod));
    Int_t pdg = ele->GetPdg();
    
    if(GetDebug() > 3) 
      printf("AliAnaElectron::MakeAnalysisFillHistograms() - PDG %d, MC TAG %d, Calorimeter %s\n", ele->GetPdg(),ele->GetTag(), (ele->GetDetector()).Data()) ;
    
    if(TMath::Abs(pdg) != AliCaloPID::kElectron) continue; 
    
    if(GetDebug() > 1) 
      printf("AliAnaElectron::MakeAnalysisFillHistograms() - ID Electron: pt %f, phi %f, eta %f\n", ele->Pt(),ele->Phi(),ele->Eta()) ;
    
    //MC tag of this electron
    Int_t mctag = ele->GetTag();
    
    //Filter for photonic electrons based on opening angle and Minv
    //cuts, also fill histograms
    Bool_t photonic = kFALSE;
    Bool_t photonic1 = kFALSE;
    photonic1 = PhotonicPrim(ele); //check against primaries
    if(photonic1) ph1++;
    Bool_t photonic2 = kFALSE;
    photonic2 = PhotonicV0(ele->GetTrackLabel(0)); //check against V0s
    if(photonic2) ph2++;
    if(photonic1 && photonic2) phB++;
    if(photonic1 || photonic2) photonic = kTRUE;
    
    //Fill electron histograms 
    Float_t ptele = ele->Pt();
    Float_t phiele = ele->Phi();
    Float_t etaele = ele->Eta();
    
    //"Best reconstructed electron spectrum" = EMCAL or tracking
    //detectors say it is an electron and it does not form a V0
    //with Minv near a relevant resonance
    if(!photonic) {
      fhPtNPEleTTE->Fill(ptele,0); //0 = no MC info
      if(ele->GetDetector() == fCalorimeter) fhPtNPEleEMCAL->Fill(ptele,0);
      if(IsDataMC()) {
        fhPtNPEleTTE->Fill(ptele,GetMCSource(mctag));
        if(ele->GetDetector() == "EMCAL") fhPtNPEleEMCAL->Fill(ptele,GetMCSource(mctag));
        if(GetMCSource(mctag) == 1) { //it's a bottom electron, now
          //get the parent's pt
          Double_t ptbHadron = GetBParentPt(ele->GetLabel());
          fhPtNPEBHadron->Fill(ptele,ptbHadron);
        } //mctag
      } //isdatamc
    } //!photonic
    
    //kept for historical reasons?
    fhPtElectron  ->Fill(ptele);
    fhPhiElectron ->Fill(ptele,phiele);
    fhEtaElectron ->Fill(ptele,etaele);
    
    if(photonic) {
      fhPtPE->Fill(ptele);
      fhPhiPE->Fill(ptele,phiele);
      fhEtaPE->Fill(ptele,etaele);
    } else {
      fhPtNPE->Fill(ptele);
      fhPhiNPE->Fill(ptele,phiele);
      fhEtaNPE->Fill(ptele,etaele);
    }
    
    if(IsDataMC()){
      if(GetMCAnalysisUtils()->CheckTagBit(mctag,AliMCAnalysisUtils::kMCConversion)){
        fhPhiConversion ->Fill(ptele,phiele);
        fhEtaConversion ->Fill(ptele,etaele);
      }
    }//Histograms with MC
    
  }// aod loop
  
  ////////////////////////////////////////////////////////
  //Fill histograms of pure MC kinematics from the stack//
  ////////////////////////////////////////////////////////
  if(IsDataMC() && stack) {
    
    //MC Jets
    TVector3 bjetVect[4];
    Int_t nPythiaGenJets = 0;
    AliGenPythiaEventHeader*  pythiaGenHeader = (AliGenPythiaEventHeader*)GetReader()->GetGenEventHeader();
    if(pythiaGenHeader){
      //Get Jets from MC header
      nPythiaGenJets = pythiaGenHeader->NTriggerJets();
      Int_t iCount = 0;
      for(int ip = 0;ip < nPythiaGenJets;++ip){
        if (iCount>3) break;
        Float_t p[4];
        pythiaGenHeader->TriggerJet(ip,p);
        TVector3 tempVect(p[0],p[1],p[2]);
        if ( TMath::Abs(tempVect.Eta())>fJetEtaCut || tempVect.Phi() < fJetPhiMin || tempVect.Phi() > fJetPhiMax) continue;
        //Only store it if it has a b-quark within dR < 0.2 of jet axis ?
        if(IsMcBJet(tempVect.Eta(),tempVect.Phi())) {
          bjetVect[iCount].SetXYZ(p[0], p[1], p[2]);
          iCount++;
        }
      }
    }
    
    Int_t nPart = GetNumAODMCParticles();
    if(GetReader()->ReadStack()) nPart = stack->GetNtrack();
    
    for(Int_t ipart = 0; ipart < nPart; ipart++) {
      
      //All the variables we want from MC particles
      Double_t px = 0.; Double_t py = 0.; Double_t pz = 0.; Double_t e = 0.;
      Double_t vx = -999.; Double_t vy = -999.; Double_t vz = -999.; Double_t vt = -999.;
      Int_t pdg = 0; Int_t mpdg = 0; Double_t mpt = 0.;
      
      if(GetReader()->ReadStack()) {
        primary = stack->Particle(ipart);
        pdg = primary->GetPdgCode();
        px = primary->Px(); py = primary->Py(); pz = primary->Pz(); e = primary->Energy();
        vx = primary->Vx(); vy = primary->Vy(); vz = primary->Vz(); vt = primary->T();
        if(primary->GetMother(0)>=0) {
          TParticle *parent = stack->Particle(primary->GetMother(0));
          if (parent) {
            mpdg = parent->GetPdgCode();
            mpt = parent->Pt();
          }
        }
      } else if(GetReader()->ReadAODMCParticles()) {
        aodprimary =  (AliAODMCParticle*)GetMCParticle(ipart);
        pdg = aodprimary->GetPdgCode();
        px = aodprimary->Px(); py = aodprimary->Py(); pz = aodprimary->Pz(); e = aodprimary->E();
        vx = aodprimary->Xv(); vy = aodprimary->Yv(); vz = aodprimary->Zv(); vt = aodprimary->T();
        Int_t parentId = aodprimary->GetMother();
        if(parentId>=0) {
          AliAODMCParticle *parent = (AliAODMCParticle*)GetMCParticle(parentId);
          if (parent) {
            mpdg = parent->GetPdgCode();
            mpt = parent->Pt();
          }
        }	
      }
      
      TLorentzVector mom(px,py,pz,e);
      TLorentzVector pos(vx,vy,vz,vt);
      Bool_t in = kFALSE;
      if(mom.Phi()*180./TMath::Pi() > 80. && mom.Phi()*180./TMath::Pi() < 190. &&
         mom.Eta() > -0.7 && mom.Eta() < 0.7) in = kTRUE;
      /////////////////////////////////
      //THIS HAS A MEM LEAK JLK 24-Oct-09
      //Bool_t in = GetFiducialCut()->IsInFiducialCut(mom,fCalorimeter);
      ////////////////////////////////
      if(mom.Pt() < GetMinPt()) continue;
      if(!in) continue;
      
      if(TMath::Abs(pdg) == 211 || TMath::Abs(pdg) == 321 || TMath::Abs(pdg) == 2212)
        fhPtMCHadron->Fill(mom.Pt());
      
      //we only care about electrons
      if(TMath::Abs(pdg) != 11) continue;
      //we only want TRACKABLE electrons (TPC 85-250cm)
      if(pos.Rho() > 200.) continue;
      //Ignore low pt electrons
      if(mom.Pt() < 0.2) continue;
      
      //find out what the ancestry of this electron is
      Int_t mctag = -1;
      Int_t input = 0;
      mctag = GetMCAnalysisUtils()->CheckOrigin(ipart,GetReader(),input);
      
      if(GetMCSource(mctag)==1) { //bottom electron
        //See if it is within dR < 0.4 of a bjet
        for(Int_t ij = 0; ij < nPythiaGenJets; ij++) {
          Double_t deta = primary->Eta() - bjetVect[ij].Eta();
          Double_t dphi = primary->Phi() - bjetVect[ij].Phi();
          Double_t dR = TMath::Sqrt(deta*deta + dphi*dphi);
          if(dR < 0.4) {
            fhMCBJetElePt->Fill(primary->Pt(),bjetVect[ij].Pt());
          }
        }
      }
      
      if ((TMath::Abs(mpdg) >500  && TMath::Abs(mpdg) <600 ) ||
          (TMath::Abs(mpdg) >5000 && TMath::Abs(mpdg) <6000 ) )
      {
        fhMCBHadronElePt->Fill(mom.Pt(), mpt); 
      }
      //CHECK THAT THIS IS CORRECTLY FILLED - SHOULD WE USE MCSOURCE HERE?
      fhPtMCElectron->Fill(mom.Pt(),0);  //0 = unfiltered
      fhPtMCElectron->Fill(mom.Pt(),GetMCSource(mctag));
      
      if(GetMCSource(mctag) == 4) {//conversion
        fhMCXYConversion->Fill(vx,vy);
        fhMCRadPtConversion->Fill(TMath::Sqrt(vx*vx+vy*vy),mom.Pt());
      }
      
      //fill ntuple
      if(fWriteNtuple) {
        fMCEleNtuple->Fill(mctag,mom.Pt(),mom.Phi(),mom.Eta(),vx,vy,vz);
      }
    }
  } //MC loop
  
  //if(GetDebug() > 0) 
  printf("\tAliAnaElectron::Photonic electron counts: ph1 %d, ph2 %d, Both %d\n",ph1,ph2,phB);
}

//__________________________________________________________________
Int_t AliAnaElectron::GetDVMBtag(AliAODTrack * tr )
{
  //This method uses the Displaced Vertex between electron-hadron
  //pairs and the primary vertex to determine whether an electron is
  //likely from a B hadron.

  Int_t ncls1 = 0;
  for(Int_t l = 0; l < 6; l++) if(TESTBIT(tr->GetITSClusterMap(),l)) ncls1++;

  fhDVMBtagQA3->Fill(ncls1);
  if (ncls1 < fITSCut) return 0;
  Double_t imp[2] = {-999.,-999.}; Double_t cov[3] = {-999.,-999.,-999.};
  Bool_t dcaOkay = GetDCA(tr,imp,cov);  //homegrown dca calculation until AOD is fixed                  
  if(!dcaOkay) {
    printf("AliAnaElectron::Problem computing DCA to primary vertex for track %d",tr->GetID());
    return 0;
  }

  fhDVMBtagQA4->Fill(imp[0]);
  if (TMath::Abs(imp[0])   > fImpactCut ) return 0;
  fhDVMBtagQA5->Fill(imp[1]);
  if (TMath::Abs(imp[1])   > fImpactCut ) return 0;

  Int_t nvtx1 = 0;
  Int_t nvtx2 = 0;
  Int_t nvtx3 = 0;

  for (Int_t k2 =0; k2 < GetCTSTracks()->GetEntriesFast() ; k2++) {
    //loop over assoc
    AliAODTrack* track2 = (AliAODTrack*)GetCTSTracks()->At(k2);
    Int_t id1 = tr->GetID();
    Int_t id2 = track2->GetID();
    if(id1 == id2) continue;

    Int_t ncls2 = 0;
    for(Int_t l = 0; l < 6; l++) if(TESTBIT(track2->GetITSClusterMap(),l)) ncls2++;
    if (ncls2 < fITSCut) continue;

    if(track2->Pt() < fAssocPtCut) continue;

    Double_t dphi = tr->Phi() - track2->Phi();
    if(dphi > TMath::Pi()) dphi -= 2*TMath::Pi();
    if(dphi < -TMath::Pi()) dphi += 2*TMath::Pi();
    Double_t deta = tr->Eta() - track2->Eta();
    Double_t dr = sqrt(deta*deta + dphi*dphi);

    if(dr > fDrCut) continue;
    Double_t sDca1 = ComputeSignDca(tr, track2, 1.0);
    if (sDca1 > fSdcaCut) nvtx1++;
    Double_t sDca2 = ComputeSignDca(tr, track2, 1.5);
    if (sDca2 > fSdcaCut) nvtx2++;
    Double_t sDca3 = ComputeSignDca(tr, track2, 1.8);
    if (sDca3 > fSdcaCut) nvtx3++;

  } //loop over hadrons

  if(GetDebug() > 0) {
    if (nvtx1>0) printf("result1 of btagging: %d \n",nvtx1);
    if (nvtx2>0) printf("result2 of btagging: %d \n",nvtx2);
    if (nvtx3>0) printf("result3 of btagging: %d \n",nvtx3);
  }

  //fill QA histograms
  fhDVMBtagCut1->Fill(nvtx1,tr->Pt());
  fhDVMBtagCut2->Fill(nvtx2,tr->Pt());
  fhDVMBtagCut3->Fill(nvtx3,tr->Pt());

  return nvtx2;

}

//__________________________________________________________________
Double_t AliAnaElectron::ComputeSignDca(AliAODTrack *tr, AliAODTrack *tr2 , float masscut)
{
  //Compute the signed dca between two tracks
  //and return the result

  Double_t signDca=-999.;
  if(GetDebug() > 2 ) printf(">>ComputeSdca:: track1 %d, track2 %d, masscut %f \n", tr->GetLabel(), tr2->GetLabel(), masscut);

  //=====Now calculate DCA between both tracks=======  
  Double_t massE = 0.000511;
  Double_t massK = 0.493677;

  Double_t vertex[3] = {-999.,-999.,-999}; //vertex
  if(GetReader()->GetDataType() != AliCaloTrackReader::kMC) {
    GetVertex(vertex); //If only one file, get the vertex from there
    //FIXME:  Add a check for whether file 2 is PYTHIA or HIJING
    //If PYTHIA, then set the vertex from file 2, if not, use the
    //vertex from file 1
    //if(GetReader()->GetSecondInputAODTree()) GetReader()->GetSecondInputAODVertex(vertex);
  }
  
  TVector3 primV(vertex[0],vertex[1],vertex[2]) ;

  if(GetDebug() > 5) printf(">>ComputeSdca:: primary vertex = %2.2f,%2.2f,%2.2f \n",vertex[0],vertex[1],vertex[2]) ;

  AliExternalTrackParam *param1 = new AliExternalTrackParam(tr);
  AliExternalTrackParam *param2 = new AliExternalTrackParam(tr2);

  //Replaced functions to get the B-field. They do not work.
  Double_t bfield[3];
  param1->GetBxByBz(bfield);
  bfield[0]=0;
  bfield[1]=0;
  bfield[2]=5.0;
  Double_t bz = 5.0; //param1->GetBz();
  printf("ZOMG1 bfield is: %f, %f",bz,bfield[2]);
  Double_t xplane1 = 0.; Double_t xplane2 = 0.;
  Double_t pairdca = param1->GetDCA(param2,bz,xplane1,xplane2);

  param1->PropagateToBxByBz(xplane1,bfield);
  param2->PropagateToBxByBz(xplane2,bfield);

  Int_t id1 = 0, id2 = 0;
  AliESDv0 bvertex(*param1,id1,*param2,id2);
  Double_t vx,vy,vz;
  bvertex.GetXYZ(vx,vy,vz);

  Double_t emom[3];
  Double_t hmom[3];
  param1->PxPyPz(emom);
  param2->PxPyPz(hmom);
  TVector3 emomAtB(emom[0],emom[1],emom[2]);
  TVector3 hmomAtB(hmom[0],hmom[1],hmom[2]);
  TVector3 secvtxpt(vx,vy,vz);
  TVector3 decayvector(0,0,0);
  decayvector = secvtxpt - primV; //decay vector from PrimVtx
  Double_t decaylength = decayvector.Mag();

  //printf("\t JLK pairDCA = %2.2f\n",pairdca);

  if(GetDebug() > 0) {
    printf(">>ComputeSdca:: mom1=%f, mom2=%f \n", emomAtB.Perp(), hmomAtB.Perp() );
    printf(">>ComputeSdca:: pairDCA=%f, length=%f \n", pairdca,decaylength );
  }

  if (masscut<1.1) fhDVMBtagQA1->Fill(pairdca,decaylength);

  if (emomAtB.Mag()>0 && pairdca < fPairDcaCut && decaylength < fDecayLenCut ) {
    TVector3 sumMom = emomAtB+hmomAtB;
    Double_t ener1 = sqrt(pow(emomAtB.Mag(),2) + massE*massE);
    Double_t ener2 = sqrt(pow(hmomAtB.Mag(),2) + massK*massK);
    Double_t ener3 = sqrt(pow(hmomAtB.Mag(),2) + massE*massE);
    Double_t mass = sqrt(pow((ener1+ener2),2) - pow(sumMom.Mag(),2));
    Double_t massPhot = sqrt(pow((ener1+ener3),2) - pow(sumMom.Mag(),2));
    Double_t sDca = decayvector.Dot(emomAtB)/emomAtB.Mag();

    if (masscut<1.1) fhDVMBtagQA2->Fill(sDca, mass);

    if (mass > masscut && massPhot > 0.1) signDca = sDca;
    
    if(GetDebug() > 0) printf("\t>>ComputeSdca:: mass=%f \n", mass);
    if(GetDebug() > 0) printf("\t>>ComputeSdca:: sec vtx-signdca :%f\n",signDca);
  }

  //clean up
  delete param1;
  delete param2;

  return signDca;
}

//__________________________________________________________________
Double_t AliAnaElectron::GetIPSignificance(AliAODTrack *tr, Double_t jetPhi)
{
  //get signed impact parameter significance of the given AOD track
  //for the given jet

  Int_t trackIndex = 0;
  Int_t ntrk = GetCTSTracks()->GetEntriesFast();
  for (Int_t k2 =0; k2 < ntrk ; k2++) {
    //loop over assoc
    AliAODTrack* track2 = (AliAODTrack*)GetCTSTracks()->At(k2);
    int id1 = tr->GetID();
    int id2 = track2->GetID();
    if(id1 == id2) {
      trackIndex = k2;//FIXME: check if GetCTSTracks stores tracks in the
		      //same order of the event
      break;
    }
  }

  Double_t significance=0;
  Double_t maxD = 10000.;
  Double_t impPar[] = {0,0};
  Double_t ipCov[]={0,0,0};
  Double_t ipVec2D[] = {0,0};

  AliVEvent* vEvent = (AliVEvent*)GetReader()->GetInputEvent();
  if(!vEvent) return -97;
  AliVVertex* vv = (AliVVertex*)vEvent->GetPrimaryVertex();
  if(!vv) return -98;
  AliVTrack* vTrack = (AliVTrack*)vEvent->GetTrack(trackIndex);
  if(!vTrack) return -99;
  AliESDtrack esdTrack(vTrack);
  //Replaced function to get B-field, it does not work
  Double_t bfield[3];
  //esdTrack.GetBxByBz(bfield);
  bfield[0]=0;
  bfield[1]=0;
  bfield[2]=5.;
  //printf("ZOMG2 ESD bfield is: %f",bfield[2]);
  if(!esdTrack.PropagateToDCABxByBz(vv, bfield, maxD, impPar, ipCov)) return -100;
  if(ipCov[0]<0) return -101;

  Double_t Pxy[] = {esdTrack.Px(), esdTrack.Py()};
  Double_t Txy[] = {esdTrack.Xv(), esdTrack.Yv()};
  Double_t Vxy[] = {vv->GetX(),  vv->GetY()};
  GetImpactParamVect(Pxy, Txy, Vxy, ipVec2D);
	Double_t phiIP = TMath::ATan2(ipVec2D[1], ipVec2D[0]) + (TMath::Abs(ipVec2D[1])-ipVec2D[1])/TMath::Abs(ipVec2D[1])*TMath::Pi();
  Double_t cosTheta = TMath::Cos(jetPhi - phiIP);
  Double_t sign = cosTheta/TMath::Abs(cosTheta);
  significance = TMath::Abs(impPar[0])/TMath::Sqrt(ipCov[0])*sign;
  //printf("\t JLK significance = %2.2f\n",significance);
  //ip = fabs(impPar[0]);
  fhIPSigBtagQA2->Fill(significance);
  return significance;
}

//__________________________________________________________________
void AliAnaElectron::GetImpactParamVect(Double_t Pxy[2], Double_t Txy[2], Double_t Vxy[2], Double_t IPxy[2])
{
  //px,py: momentum components at the origin of the track; tx, ty:
  //origin (x,y) of track; vx, vy: coordinates of primary vertex
  // analytical geometry auxiliary variables
  Double_t mr = Pxy[1]/Pxy[0]; //angular coeficient of the straight
			      //line that lies on top of track
			      //momentum
  Double_t br = Txy[1] - mr*Txy[0]; //linear coeficient of the straight
				   //line that lies on top of track
				   //momentum
  Double_t ms = -1./mr; //angular coeficient of the straight line that
		       //lies on top of the impact parameter line
  //  Double_t bs = Vxy[1] - ms*Vxy[0]; //linear coeficient of the straight
				   //line that lies on top of the
				   //impact parameter line 
  Double_t xIntersection = (mr*Txy[0] - ms*Vxy[0] + Vxy[1] - Txy[1])/(mr - ms);
  Double_t yIntersection = mr*xIntersection + br;
  //if(ceil(10000*yIntersection) - ceil(10000*(ms*xIntersection + bs))
  //!= 0 )cout<<yIntersection<<", "<<ms*xIntersection + bs<<endl;
  IPxy[0] = xIntersection - Vxy[0];
  IPxy[1] = yIntersection - Vxy[1];
  return;
}

//__________________________________________________________________
Bool_t AliAnaElectron::PhotonicPrim(const AliAODPWG4Particle* part) 
{
  //This method checks the opening angle and invariant mass of
  //electron pairs within the AliAODPWG4Particle list to see if 
  //they are likely to be photonic electrons

  Bool_t itIS = kFALSE;

  Double_t massE = 0.000511;
  Double_t massEta = 0.547;
  Double_t massRho0 = 0.770;
  Double_t massOmega = 0.782;
  Double_t massPhi = 1.020;

  Int_t pdg1 = part->GetPdg();
  Int_t trackId = part->GetTrackLabel(0);
  AliAODTrack* track = (AliAODTrack*)GetCTSTracks()->At(trackId);
  if(!track) {
    if(GetDebug() > 0) printf("AliAnaElectron::PhotonicPrim - can't get the AOD Track from the particle!  Skipping the photonic check");
    return kFALSE; //Don't proceed because we can't get the track
  }

  AliExternalTrackParam *param1 = new AliExternalTrackParam(track);

  //Loop on stored AOD electrons and compute the angle differences and Minv
  for (Int_t k2 =0; k2 < GetOutputAODBranch()->GetEntriesFast() ; k2++) {
    AliAODPWG4Particle* part2 = (AliAODPWG4Particle*) GetOutputAODBranch()->At(k2);
    Int_t track2Id = part2->GetTrackLabel(0);
    if(trackId == track2Id) continue;
    Int_t pdg2 = part2->GetPdg();
    if(TMath::Abs(pdg2) != AliCaloPID::kElectron) continue;
    if(part2->GetDetector() != fCalorimeter) continue;

    //JLK: Check opp. sign pairs only
    if(pdg1*pdg2 > 0) continue; //skip same-sign pairs

    //propagate to common vertex and check opening angle
    AliAODTrack* track2 = (AliAODTrack*)GetCTSTracks()->At(track2Id);
    if(!track2) {
      if(GetDebug() >0) printf("AliAnaElectron::PhotonicPrim - problem getting the partner track.  Continuing on to the next one");
      continue;
    }
    AliExternalTrackParam *param2 = new AliExternalTrackParam(track2);
    Int_t id1 = 0, id2 = 0;
    AliESDv0 photonVtx(*param1,id1,*param2,id2);
    Double_t vx,vy,vz;
    photonVtx.GetXYZ(vx,vy,vz);

    Double_t p1mom[3];
    Double_t p2mom[3];
    param1->PxPyPz(p1mom);
    param2->PxPyPz(p2mom);

    TVector3 p1momAtB(p1mom[0],p1mom[1],p1mom[2]);
    TVector3 p2momAtB(p2mom[0],p2mom[1],p2mom[2]);
    TVector3 sumMom = p1momAtB+p2momAtB;

    Double_t ener1 = sqrt(pow(p1momAtB.Mag(),2) + massE*massE);
    Double_t ener2 = sqrt(pow(p2momAtB.Mag(),2) + massE*massE);
    Double_t mass = sqrt(pow((ener1+ener2),2) - pow(sumMom.Mag(),2));

    Double_t dphi = p1momAtB.DeltaPhi(p2momAtB);
    fh1OpeningAngle->Fill(dphi);
    fh1MinvPhoton->Fill(mass);

    if(mass < 0.1 ||
       (mass > massEta-0.05 || mass < massEta+0.05) ||
       (mass > massRho0-0.05 || mass < massRho0+0.05) ||
       (mass > massOmega-0.05 || mass < massOmega+0.05) ||
       (mass > massPhi-0.05 || mass < massPhi+0.05)) 
      {
      
	if(GetDebug() > 0) printf("######PROBABLY A PHOTON\n");
	itIS = kTRUE;
      }
    
    //clean up
    delete param2;
    
  }

  delete param1;
  return itIS;

}

//__________________________________________________________________
Bool_t AliAnaElectron::PhotonicV0(Int_t id) 
{
  //This method checks to see whether a track that has been flagged as
  //an electron was determined to match to a V0 candidate with
  //invariant mass consistent with photon conversion

  Bool_t itIS = kFALSE;

  Double_t massEta = 0.547;
  Double_t massRho0 = 0.770;
  Double_t massOmega = 0.782;
  Double_t massPhi = 1.020;
  
  //---Get V0s---
  AliAODEvent *aod = (AliAODEvent*) GetReader()->GetInputEvent();
  int nv0s = aod->GetNumberOfV0s();
  for (Int_t iV0 = 0; iV0 < nv0s; iV0++) {
    AliAODv0 *v0 = aod->GetV0(iV0);
    if (!v0) continue;
    double radius = v0->RadiusV0();
    double mass = v0->InvMass2Prongs(0,1,11,11);
    if(GetDebug() > 0) {
      printf("## PhotonicV0() :: v0: %d, radius: %f \n", iV0 , radius );
      printf("## PhotonicV0() :: neg-id: %d, pos-id: %d, THIS id: %d\n", v0->GetNegID(), v0->GetPosID(), id);
      printf("## PhotonicV0() :: Minv(e,e): %f \n", v0->InvMass2Prongs(0,1,11,11) );
    }
    if (mass < 0.100 ||
	(mass > massEta-0.05 || mass < massEta+0.05) ||
	(mass > massRho0-0.05 || mass < massRho0+0.05) ||
	(mass > massOmega-0.05 || mass < massOmega+0.05) ||
	(mass > massPhi-0.05 || mass < massPhi+0.05)) {
      if ( id == v0->GetNegID() || id == v0->GetPosID()) {
	itIS=kTRUE;
	if(GetDebug() > 0) printf("## PhotonicV0() :: It's a conversion electron!!! \n" );
      }
    } }
  return itIS;

}

//__________________________________________________________________
Bool_t AliAnaElectron::GetDCA(const AliAODTrack* track,Double_t impPar[2], Double_t cov[3]) 
{
  //Use the Event vertex and AOD track information to get
  //a real impact parameter for the track
  //Once alice-off gets its act together and fixes the AOD, this
  //should become obsolete.

  Double_t maxD = 100000.; //max transverse IP
  if(GetReader()->GetDataType() != AliCaloTrackReader::kMC) {
    AliVEvent* ve = (AliVEvent*)GetReader()->GetInputEvent();
    AliVVertex *vv = (AliVVertex*)ve->GetPrimaryVertex();
    AliESDtrack esdTrack(track);
    Double_t bfield[3];
    //This GetBxByBz function does not work properly.
    //esdTrack.GetBxByBz(bfield);
    bfield[0]=0;
    bfield[1]=0;
    bfield[2]=5.;
    //printf("ZOMG2 ESD bfield is: %f",bfield[2]);
    Bool_t gotit = esdTrack.PropagateToDCABxByBz(vv,bfield,maxD,impPar,cov);
    //printf("\t JLK impPar = %2.2f\n",impPar[0]);
    return gotit;
  }

  return kFALSE;

}

//__________________________________________________________________
Bool_t AliAnaElectron::CheckTrack(const AliAODTrack* track, const char* type) 
{
  //Check this track to see if it is also tagged as an electron in the
  //AliAODPWG4Particle list and if it is non-photonic
  
  Bool_t pass = kFALSE;
  
  Int_t trackId = track->GetID(); //get the index in the reader
  
  Int_t naod = GetOutputAODBranch()->GetEntriesFast();
  if(GetDebug() > 3) printf("AliAnaElectron::CheckTrack() - aod branch entries %d\n", naod);
  for(Int_t iaod = 0; iaod < naod ; iaod++){
    AliAODPWG4Particle* ele =  (AliAODPWG4Particle*) (GetOutputAODBranch()->At(iaod));
    Int_t label = ele->GetTrackLabel(0);
    if(label != trackId) continue;  //skip to the next one if they don't match
    
    if(strcmp(type,"DVM")==0) { 
      if(ele->CheckBTagBit(ele->GetBtag(),AliAODPWG4Particle::kDVMTag1) ||
         ele->CheckBTagBit(ele->GetBtag(),AliAODPWG4Particle::kDVMTag2))
        pass = kTRUE;
      
    } else if (strcmp(type,"NPE")==0) {
      
      Bool_t photonic = kFALSE;
      Bool_t photonic1 = kFALSE;
      photonic1 = PhotonicPrim(ele); //check against primaries
      Bool_t photonic2 = kFALSE;
      photonic2 = PhotonicV0(ele->GetTrackLabel(0)); //check against V0s
      if(photonic1 || photonic2) photonic = kTRUE;
      
      if(!photonic) pass = kTRUE;
      
    } else {
      return kFALSE;
    }
  }
  
  return pass;
  
}

//__________________________________________________________________
Double_t AliAnaElectron::GetBParentPt(Int_t ipart)
{
  //return MC B parent pt
  if(GetReader()->ReadStack()) { //only done if we have the stack                                                                                               
    AliStack* stack = GetMCStack();
    if(!stack) {
      printf("Problem getting stack\n");
      return 0.;
    }
    TParticle* prim = stack->Particle(ipart);
    if(prim->GetMother(0)>=0) {
      Int_t mpdg = 0;
      TParticle *parent = stack->Particle(prim->GetMother(0));
      if(parent){
        mpdg = parent->GetPdgCode();
      
        if ((TMath::Abs(mpdg) >500  && TMath::Abs(mpdg) <600 ) ||
            (TMath::Abs(mpdg) >5000 && TMath::Abs(mpdg) <6000 ) )
          return parent->Pt();
      }
    }
  } else if(GetReader()->ReadAODMCParticles()){
    AliAODMCParticle* prim = (AliAODMCParticle*)GetMCParticle(ipart);
    if(prim->GetMother()>=0) {
      Int_t mpdg = 0;
      AliAODMCParticle* parent = (AliAODMCParticle*)GetMCParticle(prim->GetMother());
      if(parent){
        mpdg = parent->GetPdgCode();
        if ((TMath::Abs(mpdg) >500  && TMath::Abs(mpdg) <600 ) ||
          (TMath::Abs(mpdg) >5000 && TMath::Abs(mpdg) <6000 ) )
          return parent->Pt();
      }
    }
  }
  return 0.;
}

//__________________________________________________________________
Int_t AliAnaElectron::GetMCSource(Int_t tag)
{
  //For determining how to classify electrons using MC info
  //the number returned is the bin along one axis of 2-d histograms in
  //which to fill this electron

  //Do this first
  if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion)) return 4;

  if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCElectron)) {
    //Bottom
    if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEFromB)) return 1;
    //Charm only
    else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEFromC)
	    && !GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEFromB)) return 2;
    //Charm from bottom
    else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEFromCFromB)) return 3;
    //    //Conversion
    //else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCConversion)) return 4;
    //Dalitz
    else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCPi0Decay) 
       || GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCEtaDecay) 
       || GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCOtherDecay)) return 5; 
    //W,Z
    else if(GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCWDecay)
	    || GetMCAnalysisUtils()->CheckTagBit(tag,AliMCAnalysisUtils::kMCZDecay)) return 6;
    //Everything else
    else 
      return 7;
  } else {
    //Misidentified electron
    return 8;
  }

}

//__________________________________________________________________
Int_t AliAnaElectron::GetNumAODMCParticles() 
{
  //Get the number of AliAODMCParticles, if any
  Int_t num = 0;
  Int_t npart0 = 0;
  TClonesArray * mcparticles0 = 0x0;
//  TClonesArray * mcparticles1 = 0x0;

  if(GetReader()->ReadAODMCParticles()){
    //Get the list of MC particles
    //                                                                                                 
    mcparticles0 = GetReader()->GetAODMCParticles(0);
    if(!mcparticles0) {
     if(GetDebug() > 0) printf("AliAnaElectron::MakeAnalysisFillHistograms() -  Standard MCParticles not available!\n");
    }
//    if(GetReader()->GetSecondInputAODTree()){
//      mcparticles1 = GetReader()->GetAODMCParticles(1);
//      if(!mcparticles1 && GetDebug() > 0) {
//        printf("AliAnaElectron::MakeAnalysisFillHistograms() -  Second input MCParticles not available!\n");
//      }
//    }
    else{
      npart0 = mcparticles0->GetEntriesFast();
    }
    //Int_t npart1 = 0;
    //if(mcparticles1) npart1 = mcparticles1->GetEntriesFast();
    //Int_t npart = npart0;//+npart1;
    return npart0;

  }

  return num;
}
//__________________________________________________________________
AliAODMCParticle* AliAnaElectron::GetMCParticle(Int_t ipart) 
{
  //Get the MC particle at position ipart
  
  AliAODMCParticle* aodprimary = 0x0;
  TClonesArray * mcparticles0 = 0x0;
  //TClonesArray * mcparticles1 = 0x0;
  
  if(GetReader()->ReadAODMCParticles()){
    //Get the list of MC particles                                                                                                                           
    mcparticles0 = GetReader()->GetAODMCParticles(0);
    if(!mcparticles0) {
      if (GetDebug() > 0) printf("AliAnaElectron::MakeAnalysisFillHistograms() -  Standard MCParticles not available!\n");
    }
    //    if(GetReader()->GetSecondInputAODTree()){
    //      mcparticles1 = GetReader()->GetAODMCParticles(1);
    //      if(!mcparticles1 && GetDebug() > 0) {
    //	printf("AliAnaElectron::MakeAnalysisFillHistograms() -  Second input MCParticles not available!\n");
    //      }
    //    }
    else{
      Int_t npart0 = mcparticles0->GetEntriesFast();
      //Int_t npart1 = 0;
      //if(mcparticles1) npart1 = mcparticles1->GetEntriesFast();
      if(ipart < npart0) aodprimary = (AliAODMCParticle*)mcparticles0->At(ipart);
      //else aodprimary = (AliAODMCParticle*)mcparticles1->At(ipart-npart0);
      if(!aodprimary) {
        printf("AliAnaElectron::GetMCParticle() *** no primary ***:  label %d \n", ipart);
        return 0x0;
      }
    }
  } else {
    printf("AliAnaElectron::GetMCParticle() - Asked for AliAODMCParticle but we have a stack reader.\n");
  }
  return aodprimary;
  
}

//__________________________________________________________________
Bool_t  AliAnaElectron::IsMcBJet(Double_t jeta, Double_t jphi)
{
  //Check the jet eta,phi against that of the b-quark
  //to decide whether it is an MC B-jet
  Bool_t bjet=kFALSE;

  //      printf("MTH: McStack ,nparticles=%d \n", stack->GetNtrack() );

  AliStack* stack = 0x0;
  
  for(Int_t ipart = 0; ipart < 100; ipart++) {

    Double_t pphi = -999.;
    Double_t peta = -999.;
    Int_t pdg = 0;
    if(GetReader()->ReadStack()) {
      stack = GetMCStack();
      if(!stack) {
	printf("AliAnaElectron::IsMCBJet() *** no stack ***: \n");
	return kFALSE;
      }
      TParticle* primary = stack->Particle(ipart);
      if (!primary) continue;
      pdg = primary->GetPdgCode();
      pphi = primary->Phi();
      peta = primary->Eta();
    } else if(GetReader()->ReadAODMCParticles()) {
      AliAODMCParticle* aodprimary = GetMCParticle(ipart);
      if(!aodprimary) continue;
      pdg = aodprimary->GetPdgCode();
      pphi = aodprimary->Phi();
      peta = aodprimary->Eta();
    }
    if ( TMath::Abs(pdg) != 5) continue;
    
    //      printf("MTH: IsMcBJet : %d, pdg=%d : pt=%f \n", ipart, pdgcode, primary->Pt());
    Double_t dphi = jphi - pphi;
    Double_t deta = jeta - peta;
    Double_t dr = sqrt(deta*deta + dphi*dphi);
    
    if (dr < 0.2) {
      bjet=kTRUE;
      //printf("MTH: **** found matching MC-Bjet: PDG=%d, pt=%f,dr=%f \n", pdgcode, primary->Pt(),dr );
      break;
    }
  }
  return bjet;

}

//__________________________________________________________________
Bool_t  AliAnaElectron::IsMcDJet(Double_t jeta, Double_t jphi)
{
  //Check if this jet is a charm jet
  Bool_t cjet=kFALSE;

  AliStack* stack = 0x0;

  for(Int_t ipart = 0; ipart < 100; ipart++) {
    
    Double_t pphi = -999.;
    Double_t peta = -999.;
    Int_t pdg = 0;
    if(GetReader()->ReadStack()) {
      stack = GetMCStack();
      if(!stack) {
	printf("AliAnaElectron::IsMCDJet() *** no stack ***: \n");
	return kFALSE;
      }
      TParticle* primary = stack->Particle(ipart);
      if (!primary) continue;
      pdg = primary->GetPdgCode();
      pphi = primary->Phi();
      peta = primary->Eta();
    } else if(GetReader()->ReadAODMCParticles()) {
      AliAODMCParticle* aodprimary = GetMCParticle(ipart);
      if(!aodprimary) continue;
      pdg = aodprimary->GetPdgCode();
      pphi = aodprimary->Phi();
      peta = aodprimary->Eta();
    }

    if ( TMath::Abs(pdg) != 4) continue;

    Double_t dphi = jphi - pphi;
    Double_t deta = jeta - peta;
    Double_t dr = sqrt(deta*deta + dphi*dphi);
    
    if (dr < 0.2) {
      cjet=kTRUE;
      break;
    }
  }

  return cjet;

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
  printf("---DVM Btagging\n");
  printf("max IP-cut (e,h)       =     %f\n",fImpactCut);
  printf("min ITS-hits           =     %d\n",fITSCut);
  printf("max dR (e,h)           =     %f\n",fDrCut);
  printf("max pairDCA            =     %f\n",fPairDcaCut);
  printf("max decaylength        =     %f\n",fDecayLenCut);
  printf("min Associated Pt      =     %f\n",fAssocPtCut);
  printf("---IPSig Btagging\n");
  printf("min tag track          =     %d\n",fNTagTrkCut);
  printf("min IP significance    =     %f\n",fIPSigCut);
  printf("    \n") ;
	
} 

//________________________________________________________________________
void AliAnaElectron::ReadHistograms(TList* /* outputList */)
{
  // Needed when Terminate is executed in distributed environment                             
  // Refill analysis histograms of this class with corresponding
  // histograms in output list.   

  // Histograms of this analsys are kept in the same list as other
  // analysis, recover the position of
  // the first one and then add the next                                                      
  //Int_t index = outputList->IndexOf(outputList->FindObject(GetAddedHistogramsStringToName()+"fh1pOverE"));

  //Read histograms, must be in the same order as in
  //GetCreateOutputObject.                   
  //fh1pOverE     = (TH1F *) outputList->At(index);
  //fh1dR         = (TH1F *) outputList->At(index++);
  //fh2EledEdx    = (TH2F *) outputList->At(index++);
  //fh2MatchdEdx  = (TH2F *) outputList->At(index++);
  
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

