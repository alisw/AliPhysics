/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

/////////////////////////////////////////////////////////////
//
// AliAnalysisTaskSE for D0 candidates invariant mass histogram
// and comparison with the MC truth and cut variables distributions.
//
// Authors: A.Dainese, andrea.dainese@lnl.infn.it
// Chiara Bianchin, chiara.bianchin@pd.infn.it (invariant mass)
// Carmelo Di Giglio, carmelo.digiglio@ba.infn.it (like sign)
/////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TClonesArray.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <TList.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TDatabasePDG.h>

#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include "AliAnalysisManager.h"
#include "AliESDtrack.h"
#include "AliVertexerTracks.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskSED0Mass.h"
#include "AliNormalizationCounter.h"

ClassImp(AliAnalysisTaskSED0Mass)


//________________________________________________________________________
AliAnalysisTaskSED0Mass::AliAnalysisTaskSED0Mass():
AliAnalysisTaskSE(),
  fOutputMass(0),
  fDistr(0),
  fNentries(0), 
  fCuts(0),
  fArray(0),
  fReadMC(0),
  fCutOnDistr(0),
  fUsePid4Distr(0),
  fCounter(0),
  fNPtBins(1),
  fLsNormalization(1.),
  fFillOnlyD0D0bar(0),
  fDaughterTracks(),
  fIsSelectedCandidate(0),
  fFillVarHists(kTRUE),
  fSys(0),
  fIsRejectSDDClusters(0)
{
  // Default constructor
}

//________________________________________________________________________
AliAnalysisTaskSED0Mass::AliAnalysisTaskSED0Mass(const char *name,AliRDHFCutsD0toKpi* cuts):
  AliAnalysisTaskSE(name),
  fOutputMass(0), 
  fDistr(0),
  fNentries(0),
  fCuts(0),
  fArray(0),
  fReadMC(0),
  fCutOnDistr(0),
  fUsePid4Distr(0),
  fCounter(0),
  fNPtBins(1),
  fLsNormalization(1.),
  fFillOnlyD0D0bar(0),
  fDaughterTracks(),
  fIsSelectedCandidate(0),
  fFillVarHists(kTRUE),
  fSys(0),
  fIsRejectSDDClusters(0)
{
  // Default constructor

  fNPtBins=cuts->GetNPtBins();
    
  fCuts=cuts;

  // Output slot #1 writes into a TList container (mass with cuts)
  DefineOutput(1,TList::Class());  //My private output
  // Output slot #2 writes into a TList container (distributions)
  if(fFillVarHists) DefineOutput(2,TList::Class());  //My private output
  // Output slot #3 writes into a TH1F container (number of events)
  DefineOutput(3,TH1F::Class());  //My private output
  // Output slot #4 writes into a TList container (cuts)
  DefineOutput(4,AliRDHFCutsD0toKpi::Class());  //My private output
  // Output slot #5 writes Normalization Counter 
  DefineOutput(5,AliNormalizationCounter::Class());
}

//________________________________________________________________________
AliAnalysisTaskSED0Mass::~AliAnalysisTaskSED0Mass()
{
  if (fOutputMass) {
    delete fOutputMass;
    fOutputMass = 0;
  }
  if (fDistr) {
    delete fDistr;
    fDistr = 0;
  }
  if (fCuts) {
    delete fCuts;
    fCuts = 0;
  }
  if (fNentries){
    delete fNentries;
    fNentries = 0;
  }
  if(fCounter){
    delete fCounter;
    fCounter=0;
  }
 
}  

//________________________________________________________________________
void AliAnalysisTaskSED0Mass::Init()
{
  // Initialization

  if(fDebug > 1) printf("AnalysisTaskSED0Mass::Init() \n");

  
  AliRDHFCutsD0toKpi* copyfCuts=new AliRDHFCutsD0toKpi(*fCuts);
  const char* nameoutput=GetOutputSlot(4)->GetContainer()->GetName();
  copyfCuts->SetName(nameoutput);
  // Post the data
  PostData(4,copyfCuts);


  return;
}

//________________________________________________________________________
void AliAnalysisTaskSED0Mass::UserCreateOutputObjects()
{

  // Create the output container
  //
  if(fDebug > 1) printf("AnalysisTaskSED0Mass::UserCreateOutputObjects() \n");

  // Several histograms are more conveniently managed in a TList
  fOutputMass = new TList();
  fOutputMass->SetOwner();
  fOutputMass->SetName("listMass");

  fDistr = new TList();
  fDistr->SetOwner();
  fDistr->SetName("distributionslist");

  TString nameMass=" ",nameSgn27=" ",nameSgn=" ", nameBkg=" ", nameRfl=" ",nameMassNocutsS =" ",nameMassNocutsB =" ", namedistr=" ";

  for(Int_t i=0;i<fCuts->GetNPtBins();i++){

    nameMass="histMass_";
    nameMass+=i;
    nameSgn27="histSgn27_";
    nameSgn27+=i;
    nameSgn="histSgn_";
    nameSgn+=i;
    nameBkg="histBkg_";
    nameBkg+=i;
    nameRfl="histRfl_";
    nameRfl+=i;
    nameMassNocutsS="hMassS_";
    nameMassNocutsS+=i;
    nameMassNocutsB="hMassB_";
    nameMassNocutsB+=i;

    //histograms of cut variable distributions
    if(fFillVarHists){
      if(fReadMC){
	//  dca
	namedistr="hdcaS_";
	namedistr+=i;
	TH1F *hdcaS = new TH1F(namedistr.Data(), "DCA distribution;dca [cm]",200,0.,0.1);
	// impact parameter
	namedistr="hd0piS_";
	namedistr+=i;
	TH1F *hd0piS = new TH1F(namedistr.Data(), "Impact parameter distribution (pions);d0(#pi) [cm]",200,-0.1,0.1);

	namedistr="hd0KS_";
	namedistr+=i;
	TH1F *hd0KS = new TH1F(namedistr.Data(), "Impact parameter distribution (kaons);d0(K) [cm]",200,-0.1,0.1);
	namedistr="hd0d0S_";
	namedistr+=i;
	TH1F *hd0d0S = new TH1F(namedistr.Data(), "d_{0}#timesd_{0} distribution;d_{0}#timesd_{0} [cm^{2}]",200,-0.001,0.001);
 
	//decay lenght
	namedistr="hdeclS_";
	namedistr+=i;
	TH1F *hdeclengthS = new TH1F(namedistr.Data(), "Decay Length^{2} distribution;Decay Length^{2} [cm]",200,0,0.015);

	namedistr="hnormdeclS_";
	namedistr+=i;
	TH1F *hnormdeclengthS = new TH1F(namedistr.Data(), "Normalized Decay Length^{2} distribution;(Decay Length/Err)^{2} ",200,0,12.);

	namedistr="hdeclxyS_";
	namedistr+=i;
	TH1F* hdeclxyS=new TH1F(namedistr.Data(),"Decay Length XY distribution;Decay Length XY [cm]",200,0,0.15);
	namedistr="hnormdeclxyS_";
	namedistr+=i;
	TH1F* hnormdeclxyS=new TH1F(namedistr.Data(),"Normalized decay Length XY distribution;Decay Length XY/Err",200,0,6.); 

	namedistr="hdeclxyd0d0S_";
	namedistr+=i;
	TH2F* hdeclxyd0d0S=new TH2F(namedistr.Data(),"Correlation decay Length XY - d_{0}#times d_{0};Decay Length XY [cm];d_{0}#times d_{0}[cm^{2}]",200,0,0.15,200,-0.001,0.001);

	namedistr="hnormdeclxyd0d0S_";
	namedistr+=i;
	TH2F* hnormdeclxyd0d0S=new TH2F(namedistr.Data(),"Correlation normalized decay Length XY - d_{0}#times d_{0};Decay Length XY/Err;d_{0}#times d_{0}[cm^{2}]",200,0,6,200,-0.001,0.001);

	//  costhetapoint
	namedistr="hcosthetapointS_";
	namedistr+=i;
	TH1F *hcosthetapointS = new TH1F(namedistr.Data(), "cos#theta_{Point} distribution;cos#theta_{Point}",200,0,1.);

	namedistr="hcosthetapointxyS_";
	namedistr+=i;
	TH1F *hcosthetapointxyS = new TH1F(namedistr.Data(), "cos#theta_{Point} XYdistribution;cos#theta_{Point}",300,0.,1.);

	TH1F* tmpMS = new TH1F(nameMassNocutsS.Data(),"D^{0} invariant mass; M [GeV]; Entries",300,1.5648,2.1648); //range (MD0-300MeV, mD0 + 300MeV)
	tmpMS->Sumw2();


	fDistr->Add(hdcaS);

	fDistr->Add(hd0piS);
	fDistr->Add(hd0KS);

	fDistr->Add(hd0d0S);

	fDistr->Add(hcosthetapointS);

	fDistr->Add(hcosthetapointxyS);

	fDistr->Add(hdeclengthS);

	fDistr->Add(hnormdeclengthS);

	fDistr->Add(hdeclxyS);

	fDistr->Add(hnormdeclxyS);

	fDistr->Add(hdeclxyd0d0S);
	fDistr->Add(hnormdeclxyd0d0S);

	fDistr->Add(tmpMS);
      }

      //  dca
      namedistr="hdcaB_";
      namedistr+=i;
      TH1F *hdcaB = new TH1F(namedistr.Data(), "DCA distribution;dca [cm]",200,0.,0.1);

      // impact parameter
      namedistr="hd0B_";
      namedistr+=i;
      TH1F *hd0B = new TH1F(namedistr.Data(), "Impact parameter distribution (both);d0 [cm]",200,-0.1,0.1);
 
      namedistr="hd0d0B_";
      namedistr+=i;
      TH1F *hd0d0B = new TH1F(namedistr.Data(), "d_{0}#timesd_{0} distribution;d_{0}#timesd_{0} [cm^{2}]",200,-0.001,0.001);

      //decay lenght
      namedistr="hdeclB_";
      namedistr+=i;
      TH1F *hdeclengthB = new TH1F(namedistr.Data(), "Decay Length^{2} distribution;Decay Length^{2} [cm^{2}]",200,0,0.015);

      namedistr="hnormdeclB_";
      namedistr+=i;
      TH1F *hnormdeclengthB = new TH1F(namedistr.Data(), "Normalized Decay Length distribution;(Decay Length/Err)^{2} ",200,0,12.);

      namedistr="hdeclxyB_";
      namedistr+=i;
      TH1F* hdeclxyB=new TH1F(namedistr.Data(),"Decay Length XY distribution;Decay Length XY [cm]",200,0,0.15);
      namedistr="hnormdeclxyB_";
      namedistr+=i;
      TH1F* hnormdeclxyB=new TH1F(namedistr.Data(),"Normalized decay Length XY distribution;Decay Length XY/Err",200,0,6.); 

      namedistr="hdeclxyd0d0B_";
      namedistr+=i;
      TH2F* hdeclxyd0d0B=new TH2F(namedistr.Data(),"Correlation decay Length XY - d_{0}#times d_{0};Decay Length XY [cm];d_{0}#times d_{0}[cm^{2}]",200,0,0.15,200,-0.001,0.001);

      namedistr="hnormdeclxyd0d0B_";
      namedistr+=i;
      TH2F* hnormdeclxyd0d0B=new TH2F(namedistr.Data(),"Correlation normalized decay Length XY - d_{0}#times d_{0};Decay Length XY/Err;d_{0}#times d_{0}[cm^{2}]",200,0,6,200,-0.001,0.001);

      //  costhetapoint
      namedistr="hcosthetapointB_";
      namedistr+=i;
      TH1F *hcosthetapointB = new TH1F(namedistr.Data(), "cos#theta_{Point} distribution;cos#theta_{Point}",200,0,1.);

      namedistr="hcosthetapointxyB_";
      namedistr+=i;
      TH1F *hcosthetapointxyB = new TH1F(namedistr.Data(), "cos#theta_{Point} XY distribution;cos#theta_{Point} XY",300,0.,1.);

      TH1F* tmpMB = new TH1F(nameMassNocutsB.Data(),"D^{0} invariant mass; M [GeV]; Entries",300,1.5648,2.1648); //range (MD0-300MeV, mD0 + 300MeV)
      tmpMB->Sumw2();



      fDistr->Add(hdcaB);

      fDistr->Add(hd0B);

      fDistr->Add(hd0d0B);

      fDistr->Add(hcosthetapointB);

      fDistr->Add(hcosthetapointxyB);

      fDistr->Add(hdeclengthB);

      fDistr->Add(hnormdeclengthB);

      fDistr->Add(hdeclxyB);

      fDistr->Add(hnormdeclxyB);

      fDistr->Add(hdeclxyd0d0B);
      fDistr->Add(hnormdeclxyd0d0B);

      fDistr->Add(tmpMB);

      //histograms filled only when the secondary vertex is recalculated w/o the daughter tracks (as requested in the cut object)

      if(fCuts->GetIsPrimaryWithoutDaughters()){
	if(fReadMC){
	  namedistr="hd0vpiS_";
	  namedistr+=i;
	  TH1F *hd0vpiS = new TH1F(namedistr.Data(), "Impact parameter distribution (pions)(vtx w/o these tracks);d0(#pi) [cm]",200,-0.1,0.1);
	  
	  namedistr="hd0vKS_";
	  namedistr+=i;
	  TH1F *hd0vKS = new TH1F(namedistr.Data(), "Impact parameter distribution (kaons) (vtx w/o these tracks);d0(K) [cm]",200,-0.1,0.1);

	  namedistr="hd0d0vS_";
	  namedistr+=i;
	  TH1F *hd0d0vS = new TH1F(namedistr.Data(), "d_{0}#timesd_{0} distribution (vtx w/o these tracks);d_{0}#timesd_{0} [cm^{2}]",200,-0.001,0.001);
 
	  namedistr="hdeclvS_";
	  namedistr+=i;
	  TH1F *hdeclengthvS = new TH1F(namedistr.Data(), "Decay Length distribution (vtx w/o tracks);Decay Length [cm]",200,0,0.6);

	  namedistr="hnormdeclvS_";
	  namedistr+=i;
	  TH1F *hnormdeclengthvS = new TH1F(namedistr.Data(), "Normalized Decay Length distribution (vtx w/o tracks);Decay Length/Err ",200,0,10.);

	  fDistr->Add(hd0vpiS);
	  fDistr->Add(hd0vKS);
	  fDistr->Add(hd0d0vS);
	  fDistr->Add(hdeclengthvS);
	  fDistr->Add(hnormdeclengthvS);

	}

	namedistr="hd0vmoresB_";
	namedistr+=i;
	TH1F *hd0vmoresB = new TH1F(namedistr.Data(), "Impact parameter distribution (both);d0 [cm]",200,-0.1,0.1);

	namedistr="hd0d0vmoresB_";
	namedistr+=i;
	TH1F *hd0d0vmoresB = new TH1F(namedistr.Data(), "Impact parameter distribution (prong +);d0 [cm]",200,-0.001,0.001);


	namedistr="hd0vB_";
	namedistr+=i;
	TH1F *hd0vB = new TH1F(namedistr.Data(), "Impact parameter distribution (vtx w/o these tracks);d0 [cm]",200,-0.1,0.1);

	namedistr="hd0vp0B_";
	namedistr+=i;
	TH1F *hd0vp0B = new TH1F(namedistr.Data(), "Impact parameter distribution (prong + ** vtx w/o these tracks);d0 [cm]",200,-0.1,0.1);

	namedistr="hd0vp1B_";
	namedistr+=i;
	TH1F *hd0vp1B = new TH1F(namedistr.Data(), "Impact parameter distribution (prong - ** vtx w/o these tracks);d0 [cm]",200,-0.1,0.1);

	namedistr="hd0d0vB_";
	namedistr+=i;
	TH1F *hd0d0vB = new TH1F(namedistr.Data(), "d_{0}#timesd_{0} distribution (vtx w/o these tracks);d_{0}#timesd_{0} [cm^{2}]",200,-0.001,0.001);

	namedistr="hdeclvB_";
	namedistr+=i;
	TH1F *hdeclengthvB = new TH1F(namedistr.Data(), "Decay Length distribution (vtx w/o tracks);Decay Length [cm]",200,0,0.6);

	namedistr="hnormdeclvB_";
	namedistr+=i;
	TH1F *hnormdeclengthvB = new TH1F(namedistr.Data(), "Normalized Decay Length distribution (vtx w/o tracks);Decay Length/Err ",200,0,10.);

	fDistr->Add(hd0vB);
	fDistr->Add(hd0vp0B);
	fDistr->Add(hd0vp1B);
	fDistr->Add(hd0vmoresB);

	fDistr->Add(hd0d0vB);
	fDistr->Add(hd0d0vmoresB);

	fDistr->Add(hdeclengthvB);

	fDistr->Add(hnormdeclengthvB);
      }


    }
    //histograms of invariant mass distributions


    //MC signal
    if(fReadMC){
      TH1F* tmpSt = new TH1F(nameSgn.Data(), "D^{0} invariant mass - MC; M [GeV]; Entries",200,1.5648,2.1648);

      TH1F *tmpSl=(TH1F*)tmpSt->Clone();
      tmpSt->Sumw2();
      tmpSl->Sumw2();

      //Reflection: histo filled with D0Mass which pass the cut (also) as D0bar and with D0bar which pass (also) the cut as D0
      TH1F* tmpRt = new TH1F(nameRfl.Data(), "Reflected signal invariant mass - MC; M [GeV]; Entries",200,1.5648,2.1648);
      //TH1F *tmpRl=(TH1F*)tmpRt->Clone();
      TH1F* tmpBt = new TH1F(nameBkg.Data(), "Background invariant mass - MC; M [GeV]; Entries",200,1.5648,2.1648);
      //TH1F *tmpBl=(TH1F*)tmpBt->Clone();
      tmpBt->Sumw2();
      //tmpBl->Sumw2();
      tmpRt->Sumw2();
      //tmpRl->Sumw2();

      fOutputMass->Add(tmpSt);
      fOutputMass->Add(tmpRt);
      fOutputMass->Add(tmpBt);

    }

    //mass
    TH1F* tmpMt = new TH1F(nameMass.Data(),"D^{0} invariant mass; M [GeV]; Entries",200,1.5648,2.1648);
    //TH1F *tmpMl=(TH1F*)tmpMt->Clone();
    tmpMt->Sumw2();
    //tmpMl->Sumw2();
    //distribution w/o cuts large range
    // TH1F* tmpMS = new TH1F(nameMassNocutsS.Data(),"D^{0} invariant mass; M [GeV]; Entries",300,0.7,3.);

    fOutputMass->Add(tmpMt);

    if(fSys==0){ //histograms filled only in pp to save time in PbPb
      if(fFillVarHists){

	if(fReadMC){
	  //  pT
	  namedistr="hptpiS_";
	  namedistr+=i;
	  TH1F *hptpiS = new TH1F(namedistr.Data(), "P_{T} distribution (pions);p_{T} [GeV/c]",200,0.,8.);
	
	  namedistr="hptKS_";
	  namedistr+=i;
	  TH1F *hptKS = new TH1F(namedistr.Data(), "P_{T} distribution (kaons);p_{T} [GeV/c]",200,0.,8.);

	  //  costhetastar
	  namedistr="hcosthetastarS_";
	  namedistr+=i;
	  TH1F *hcosthetastarS = new TH1F(namedistr.Data(), "cos#theta* distribution;cos#theta*",200,-1.,1.);

	  //pT no mass cut

	  namedistr="hptpiSnoMcut_";
	  namedistr+=i;
	  TH1F *hptpiSnoMcut = new TH1F(namedistr.Data(), "P_{T} distribution (pions);p_{T} [GeV/c]",200,0.,8.);

	  namedistr="hptKSnoMcut_";
	  namedistr+=i;
	  TH1F *hptKSnoMcut = new TH1F(namedistr.Data(), "P_{T} distribution (kaons);p_{T} [GeV/c]",200,0.,8.);

	  fDistr->Add(hptpiS);
	  fDistr->Add(hptKS);
	  fDistr->Add(hcosthetastarS);

	  fDistr->Add(hptpiSnoMcut);
	  fDistr->Add(hptKSnoMcut);

	  //  costhetapoint vs d0 or d0d0
	  namedistr="hcosthpointd0S_";
	  namedistr+=i;
	  TH2F *hcosthpointd0S= new TH2F(namedistr.Data(),"Correlation cos#theta_{Point}-d_{0};cos#theta_{Point};d_{0} [cm^{2}]",200,0,1.,200,-0.001,0.001);
	  namedistr="hcosthpointd0d0S_";
	  namedistr+=i;
	  TH2F *hcosthpointd0d0S= new TH2F(namedistr.Data(),"Correlation cos#theta_{Point}-d_{0}#timesd_{0};cos#theta_{Point};d_{0}#timesd_{0} [cm^{2}]",200,0,1.,200,-0.001,0.001);

	  fDistr->Add(hcosthpointd0S);
	  fDistr->Add(hcosthpointd0d0S);

	  //to compare with AliAnalysisTaskCharmFraction
	  TH1F* tmpS27t = new TH1F(nameSgn27.Data(),"D^{0} invariant mass in M(D^{0}) +/- 27 MeV - MC; M [GeV]; Entries",200,1.5648,2.1648);
	  TH1F *tmpS27l=(TH1F*)tmpS27t->Clone();
	  tmpS27t->Sumw2();
	  tmpS27l->Sumw2();
 
	  fOutputMass->Add(tmpS27t);
	  fOutputMass->Add(tmpS27l);

	}

	//  pT
	namedistr="hptB_";
	namedistr+=i;
	TH1F *hptB = new TH1F(namedistr.Data(), "P_{T} distribution;p_{T} [GeV/c]",200,0.,8.);

	//  costhetastar
	namedistr="hcosthetastarB_";
	namedistr+=i;
	TH1F *hcosthetastarB = new TH1F(namedistr.Data(), "cos#theta* distribution;cos#theta*",200,-1.,1.);

	//pT no mass cut
	namedistr="hptB1prongnoMcut_";
	namedistr+=i;
	TH1F *hptB1pnoMcut = new TH1F(namedistr.Data(), "P_{T} distribution;p_{T} [GeV/c]",200,0.,8.);

	namedistr="hptB2prongsnoMcut_";
	namedistr+=i;
	TH1F *hptB2pnoMcut = new TH1F(namedistr.Data(), "P_{T} distribution;p_{T} [GeV/c]",200,0.,8.);
    
	fDistr->Add(hptB);
	fDistr->Add(hcosthetastarB);

	fDistr->Add(hptB1pnoMcut);
	fDistr->Add(hptB2pnoMcut);

	//impact parameter of negative/positive track
	namedistr="hd0p0B_";
	namedistr+=i;
	TH1F *hd0p0B = new TH1F(namedistr.Data(), "Impact parameter distribution (prong +);d0 [cm]",200,-0.1,0.1);

	namedistr="hd0p1B_";
	namedistr+=i;
	TH1F *hd0p1B = new TH1F(namedistr.Data(), "Impact parameter distribution (prong -);d0 [cm]",200,-0.1,0.1);

	//impact parameter corrected for strangeness
	namedistr="hd0moresB_";
	namedistr+=i;
	TH1F *hd0moresB = new TH1F(namedistr.Data(), "Impact parameter distribution (both);d0 [cm]",200,-0.1,0.1);

	namedistr="hd0d0moresB_";
	namedistr+=i;
	TH1F *hd0d0moresB = new TH1F(namedistr.Data(), "Impact parameter distribution (prong +);d0 [cm]",200,-0.001,0.001);

    
	namedistr="hcosthetapointmoresB_";
	namedistr+=i;
	TH1F *hcosthetapointmoresB = new TH1F(namedistr.Data(), "cos#theta_{Point} distribution;cos#theta_{Point}",200,0,1.);

	//  costhetapoint vs d0 or d0d0
	namedistr="hcosthpointd0B_";
	namedistr+=i;
	TH2F *hcosthpointd0B= new TH2F(namedistr.Data(),"Correlation cos#theta_{Point}-d_{0};cos#theta_{Point};d_{0} [cm^{2}]",200,0,1.,200,-0.001,0.001);

	namedistr="hcosthpointd0d0B_";
	namedistr+=i;
	TH2F *hcosthpointd0d0B= new TH2F(namedistr.Data(),"Correlation cos#theta_{Point}-d_{0}#timesd_{0};cos#theta_{Point};d_{0}#timesd_{0} [cm^{2}]",200,0,1.,200,-0.001,0.001);

	fDistr->Add(hd0p0B);
	fDistr->Add(hd0p1B);

	fDistr->Add(hd0moresB);
	fDistr->Add(hd0d0moresB);
	fDistr->Add(hcosthetapointmoresB);


	fDistr->Add(hcosthpointd0B);


	fDistr->Add(hcosthpointd0d0B);
      }

    } //end pp histo
  }


  //for Like sign analysis

  if(fArray==1){
    namedistr="hpospair";
    TH1F* hpospair=new TH1F(namedistr.Data(),"Number of positive pairs",fCuts->GetNPtBins(),-0.5,fCuts->GetNPtBins()-0.5);
    namedistr="hnegpair";
    TH1F* hnegpair=new TH1F(namedistr.Data(),"Number of negative pairs",fCuts->GetNPtBins(),-0.5,fCuts->GetNPtBins()-0.5);
    fOutputMass->Add(hpospair);
    fOutputMass->Add(hnegpair);
  }

  const char* nameoutput=GetOutputSlot(3)->GetContainer()->GetName();

  fNentries=new TH1F(nameoutput, "Integral(1,2) = number of AODs *** Integral(2,3) = number of candidates selected with cuts *** Integral(3,4) = number of D0 selected with cuts *** Integral(4,5) = events with good vertex ***  Integral(5,6) = pt out of bounds", 17,-0.5,16.5);

  fNentries->GetXaxis()->SetBinLabel(1,"nEventsAnal");
  fNentries->GetXaxis()->SetBinLabel(2,"nCandSel(Cuts)");
  if(fReadMC) fNentries->GetXaxis()->SetBinLabel(3,"nD0Selected");
  else fNentries->GetXaxis()->SetBinLabel(3,"Dstar<-D0");
  fNentries->GetXaxis()->SetBinLabel(4,"nEventsGoodVtxS");
  fNentries->GetXaxis()->SetBinLabel(5,"ptbin = -1");
  fNentries->GetXaxis()->SetBinLabel(6,"no daughter");
  if(fSys==0) fNentries->GetXaxis()->SetBinLabel(7,"nCandSel(Tr)");
  if(fFillVarHists){
    fNentries->GetXaxis()->SetBinLabel(8,"PID=0");
    fNentries->GetXaxis()->SetBinLabel(9,"PID=1");
    fNentries->GetXaxis()->SetBinLabel(10,"PID=2");
    fNentries->GetXaxis()->SetBinLabel(11,"PID=3");
  }
  if(fReadMC && fSys==0){
    fNentries->GetXaxis()->SetBinLabel(12,"K");
    fNentries->GetXaxis()->SetBinLabel(13,"Lambda");
  }
  fNentries->GetXaxis()->SetBinLabel(14,"Pile-up Rej");
  fNentries->GetXaxis()->SetBinLabel(15,"N. of 0SMH");
  if(fSys==1) fNentries->GetXaxis()->SetBinLabel(16,"Nev in centr");
  if(fIsRejectSDDClusters) fNentries->GetXaxis()->SetBinLabel(17,"SDD-Cls Rej");
  fNentries->GetXaxis()->SetNdivisions(1,kFALSE);

  fCounter = new AliNormalizationCounter(Form("%s",GetOutputSlot(5)->GetContainer()->GetName()));
  fCounter->Init();

  // Post the data
  PostData(1,fOutputMass);
  if(fFillVarHists) PostData(2,fDistr);
  PostData(3,fNentries);
  PostData(5,fCounter);  
  return;
}

//________________________________________________________________________
void AliAnalysisTaskSED0Mass::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event:
  // heavy flavor candidates association to MC truth
  //cout<<"I'm in UserExec"<<endl;


  //cuts order
  //       printf("    |M-MD0| [GeV]    < %f\n",fD0toKpiCuts[0]);
  //     printf("    dca    [cm]  < %f\n",fD0toKpiCuts[1]);
  //     printf("    cosThetaStar     < %f\n",fD0toKpiCuts[2]);
  //     printf("    pTK     [GeV/c]    > %f\n",fD0toKpiCuts[3]);
  //     printf("    pTpi    [GeV/c]    > %f\n",fD0toKpiCuts[4]);
  //     printf("    |d0K|  [cm]  < %f\n",fD0toKpiCuts[5]);
  //     printf("    |d0pi| [cm]  < %f\n",fD0toKpiCuts[6]);
  //     printf("    d0d0  [cm^2] < %f\n",fD0toKpiCuts[7]);
  //     printf("    cosThetaPoint    > %f\n",fD0toKpiCuts[8]);
  

  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());

  TString bname;
  if(fArray==0){ //D0 candidates
    // load D0->Kpi candidates
    //cout<<"D0 candidates"<<endl;
    bname="D0toKpi";
  } else { //LikeSign candidates
    // load like sign candidates
    //cout<<"LS candidates"<<endl;
    bname="LikeSign2Prong";
  }

  TClonesArray *inputArray=0;
 
  if(!aod && AODEvent() && IsStandardAOD()) {
    // In case there is an AOD handler writing a standard AOD, use the AOD 
    // event in memory rather than the input (ESD) event.    
    aod = dynamic_cast<AliAODEvent*> (AODEvent());
    // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
    // have to taken from the AOD event hold by the AliAODExtension
    AliAODHandler* aodHandler = (AliAODHandler*) 
      ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());

    if(aodHandler->GetExtensions()) {
      AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
      AliAODEvent* aodFromExt = ext->GetAOD();
      inputArray=(TClonesArray*)aodFromExt->GetList()->FindObject(bname.Data());
    }
  } else if(aod) {
    inputArray=(TClonesArray*)aod->GetList()->FindObject(bname.Data());
  }


  if(!inputArray || !aod) {
    printf("AliAnalysisTaskSED0Mass::UserExec: input branch not found!\n");
    return;
  }
  
  // fix for temporary bug in ESDfilter
  // the AODs with null vertex pointer didn't pass the PhysSel
  if(!aod->GetPrimaryVertex() || TMath::Abs(aod->GetMagneticField())<0.001) return;


  TClonesArray *mcArray = 0;
  AliAODMCHeader *mcHeader = 0;

  if(fReadMC) {
    // load MC particles
    mcArray = (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if(!mcArray) {
      printf("AliAnalysisTaskSED0Mass::UserExec: MC particles branch not found!\n");
      return;
    }
    
    // load MC header
    mcHeader = (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!mcHeader) {
      printf("AliAnalysisTaskSED0Mass::UserExec: MC header branch not found!\n");
      return;
    }
  }
  
  //printf("VERTEX Z %f %f\n",vtx1->GetZ(),mcHeader->GetVtxZ());
  
  //histogram filled with 1 for every AOD
  fNentries->Fill(0);
  fCounter->StoreEvent(aod,fCuts,fReadMC); 
  //fCounter->StoreEvent(aod,fReadMC); 

  // trigger class for PbPb C0SMH-B-NOPF-ALLNOTRD, C0SMH-B-NOPF-ALL
  TString trigclass=aod->GetFiredTriggerClasses();
  if(trigclass.Contains("C0SMH-B-NOPF-ALLNOTRD") || trigclass.Contains("C0SMH-B-NOPF-ALL")) fNentries->Fill(14);

  if(!fCuts->IsEventSelected(aod)) {
    if(fCuts->GetWhyRejection()==1) // rejected for pileup
      fNentries->Fill(13);
    if(fSys==1 && (fCuts->GetWhyRejection()==2 || fCuts->GetWhyRejection()==3)) fNentries->Fill(15);
    return;
  }

  // Check the Nb of SDD clusters
  if (fIsRejectSDDClusters) { 
    Bool_t skipEvent = kFALSE;
    Int_t ntracks = 0;
    if (aod) ntracks = aod->GetNTracks();
    for(Int_t itrack=0; itrack<ntracks; itrack++) { // loop on tacks
      //    ... get the track
      AliAODTrack * track = aod->GetTrack(itrack);
      if(TESTBIT(track->GetITSClusterMap(),2) || TESTBIT(track->GetITSClusterMap(),3) ){
	skipEvent=kTRUE;
	fNentries->Fill(16);
	break;
      }
    }
    if (skipEvent) return;
  }
  
  // AOD primary vertex
  AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();

  Bool_t isGoodVtx=kFALSE;

  //vtx1->Print();
  TString primTitle = vtx1->GetTitle();
  if(primTitle.Contains("VertexerTracks") && vtx1->GetNContributors()>0) {
    isGoodVtx=kTRUE;
    fNentries->Fill(3);
  }

  // loop over candidates
  Int_t nInD0toKpi = inputArray->GetEntriesFast();
  if(fDebug>2) printf("Number of D0->Kpi: %d\n",nInD0toKpi);

  // FILE *f=fopen("4display.txt","a");
  // fprintf(f,"Number of D0->Kpi: %d\n",nInD0toKpi);
  Int_t nSelectedloose=0,nSelectedtight=0;  
  for (Int_t iD0toKpi = 0; iD0toKpi < nInD0toKpi; iD0toKpi++) {
    AliAODRecoDecayHF2Prong *d = (AliAODRecoDecayHF2Prong*)inputArray->UncheckedAt(iD0toKpi);
 
    if(d->GetSelectionMap()) if(!d->HasSelectionBit(AliRDHFCuts::kD0toKpiCuts)){
	fNentries->Fill(2);	
	continue; //skip the D0 from Dstar
      }

    // Bool_t unsetvtx=kFALSE;
    // if(!d->GetOwnPrimaryVtx()) {
    //   d->SetOwnPrimaryVtx(vtx1); // needed to compute all variables
    //   unsetvtx=kTRUE;
    // }
  
    
    if ( fCuts->IsInFiducialAcceptance(d->Pt(),d->Y(421)) ) {
      nSelectedloose++;
      nSelectedtight++;      
      if(fSys==0){
	if(fCuts->IsSelected(d,AliRDHFCuts::kTracks,aod))fNentries->Fill(6);       
      }
      Int_t ptbin=fCuts->PtBin(d->Pt());
      if(ptbin==-1) {fNentries->Fill(4); continue;} //out of bounds
      fIsSelectedCandidate=fCuts->IsSelected(d,AliRDHFCuts::kCandidate,aod); //selected
      if(fFillVarHists) {
	//if(!fCutOnDistr || (fCutOnDistr && fIsSelectedCandidate)) {
	fDaughterTracks.AddAt((AliAODTrack*)d->GetDaughter(0),0);
	fDaughterTracks.AddAt((AliAODTrack*)d->GetDaughter(1),1);
	//check daughters
	if(!fDaughterTracks.UncheckedAt(0) || !fDaughterTracks.UncheckedAt(1)) {
	  AliDebug(1,"at least one daughter not found!");
	  fNentries->Fill(5);
	  fDaughterTracks.Clear();
	  continue;
	}
	//}
	FillVarHists(aod,d,mcArray,fCuts,fDistr);
      }

      FillMassHists(d,mcArray,fCuts,fOutputMass);
    }
  
    fDaughterTracks.Clear();
    //if(unsetvtx) d->UnsetOwnPrimaryVtx();
  } //end for prongs
  fCounter->StoreCandidates(aod,nSelectedloose,kTRUE);  
  fCounter->StoreCandidates(aod,nSelectedtight,kFALSE);  
  // Post the data
  PostData(1,fOutputMass);
  if(fFillVarHists) PostData(2,fDistr);
  PostData(3,fNentries);
  PostData(5,fCounter);
  return;
}

//____________________________________________________________________________
void AliAnalysisTaskSED0Mass::FillVarHists(AliAODEvent* aod,AliAODRecoDecayHF2Prong *part, TClonesArray *arrMC, AliRDHFCutsD0toKpi *cuts, TList *listout){
  //
  // function used in UserExec to fill variable histograms:
  //


  Int_t pdgDgD0toKpi[2]={321,211};
  Int_t lab=-9999;
  if(fReadMC) lab=part->MatchToMC(421,arrMC,2,pdgDgD0toKpi); //return MC particle label if the array corresponds to a D0, -1 if not (cf. AliAODRecoDecay.cxx)
  //Double_t pt = d->Pt(); //mother pt
  Int_t isSelectedPID=3;
  if(!fReadMC || (fReadMC && fUsePid4Distr)) isSelectedPID=cuts->IsSelectedPID(part); //0 rejected,1 D0,2 Dobar, 3 both
  if (isSelectedPID==0)fNentries->Fill(7);
  if (isSelectedPID==1)fNentries->Fill(8);
  if (isSelectedPID==2)fNentries->Fill(9);
  if (isSelectedPID==3)fNentries->Fill(10);
  //fNentries->Fill(8+isSelectedPID);

  if(fCutOnDistr && !fIsSelectedCandidate) return; 
  //printf("\nif no cuts or cuts passed\n");


  //add distr here
  UInt_t pdgs[2];
    
  Double_t mPDG=TDatabasePDG::Instance()->GetParticle(421)->Mass();
  pdgs[0]=211;
  pdgs[1]=321;
  Double_t minvD0 = part->InvMassD0();
  pdgs[1]=211;
  pdgs[0]=321;
  Double_t minvD0bar = part->InvMassD0bar();
  //cout<<"inside mass cut"<<endl;

  Double_t invmasscut=0.03;

  TString fillthispi="",fillthisK="",fillthis="";

  Int_t ptbin=cuts->PtBin(part->Pt());

  Double_t dz1[2],dz2[2],covar1[3],covar2[3];//,d0xd0proper,errd0xd0proper;
  dz1[0]=-99; dz2[0]=-99;
  Double_t d0[2];
  Double_t decl[2]={-99,-99};
  Bool_t recalcvtx=kFALSE;



  if(fCuts->GetIsPrimaryWithoutDaughters()){
    recalcvtx=kTRUE;
    //recalculate vertex
    AliAODVertex *vtxProp=0x0;
    vtxProp=GetPrimaryVtxSkipped(aod);
    if(vtxProp) {
      part->SetOwnPrimaryVtx(vtxProp);
      //Bool_t unsetvtx=kTRUE;
      //Calculate d0 for daughter tracks with Vtx Skipped
      AliESDtrack *esdtrack1=new AliESDtrack((AliVTrack*)fDaughterTracks.UncheckedAt(0));
      AliESDtrack *esdtrack2=new AliESDtrack((AliVTrack*)fDaughterTracks.UncheckedAt(1));
      esdtrack1->PropagateToDCA(vtxProp,aod->GetMagneticField(),1.,dz1,covar1);
      esdtrack2->PropagateToDCA(vtxProp,aod->GetMagneticField(),1.,dz2,covar2);
      delete vtxProp; vtxProp=NULL;
      delete esdtrack1;
      delete esdtrack2;
    }

    d0[0]=dz1[0];
    d0[1]=dz2[0];

    decl[0]=part->DecayLength2();
    decl[1]=part->NormalizedDecayLength2();
    part->UnsetOwnPrimaryVtx();
  
  }

  Double_t cosThetaStarD0 = 99;
  Double_t cosThetaStarD0bar = 99;
  Double_t cosPointingAngle = 99;
  Double_t normalizedDecayLength2 = -1, normalizedDecayLengthxy=-1;
  Double_t decayLength2 = -1, decayLengthxy=-1;
  Double_t ptProng[2]={-99,-99};
  Double_t d0Prong[2]={-99,-99};
  

  //disable the PID
  if(!fUsePid4Distr) isSelectedPID=0;
  if((lab>=0 && fReadMC) || (!fReadMC && isSelectedPID)){ //signal (from MC or PID)
   
    //check pdg of the prongs
    AliAODTrack *prong0=(AliAODTrack*)fDaughterTracks.UncheckedAt(0);
    AliAODTrack *prong1=(AliAODTrack*)fDaughterTracks.UncheckedAt(1);

    if(!prong0 || !prong1) {
      return;
    }
    Int_t labprong[2];
    if(fReadMC){
      labprong[0]=prong0->GetLabel();
      labprong[1]=prong1->GetLabel();
    }
    AliAODMCParticle *mcprong=0;
    Int_t pdgProng[2]={0,0};

    for (Int_t iprong=0;iprong<2;iprong++){
      if(fReadMC && labprong[iprong]>=0) {
	mcprong= (AliAODMCParticle*)arrMC->At(labprong[iprong]);
	pdgProng[iprong]=mcprong->GetPdgCode();
      }
    }
    
    if(fSys==0){
      //no mass cut ditributions: ptbis
	
      fillthispi="hptpiSnoMcut_";
      fillthispi+=ptbin;

      fillthisK="hptKSnoMcut_";
      fillthisK+=ptbin;

      if ((TMath::Abs(pdgProng[0]) == 211 && TMath::Abs(pdgProng[1]) == 321)
	  || (isSelectedPID==1 || isSelectedPID==3)){
	((TH1F*)listout->FindObject(fillthispi))->Fill(prong0->Pt());
	((TH1F*)listout->FindObject(fillthisK))->Fill(prong1->Pt());
      }

      if ((TMath::Abs(pdgProng[0]) == 321 && TMath::Abs(pdgProng[1]) == 211)
	  || (isSelectedPID==2 || isSelectedPID==3)){
	((TH1F*)listout->FindObject(fillthisK))->Fill(prong0->Pt());
	((TH1F*)listout->FindObject(fillthispi))->Fill(prong1->Pt());
      }
    }

    //no mass cut ditributions: mass
    fillthis="hMassS_";
    fillthis+=ptbin;
      
    if ((fReadMC && ((AliAODMCParticle*)arrMC->At(lab))->GetPdgCode() == 421)
	|| (!fReadMC && (isSelectedPID==1 || isSelectedPID==3))){//D0
      ((TH1F*)listout->FindObject(fillthis))->Fill(minvD0);
    }
    else { //D0bar
      if(fReadMC || (!fReadMC && isSelectedPID > 1))
	((TH1F*)listout->FindObject(fillthis))->Fill(minvD0bar);
    }

    //apply cut on invariant mass on the pair
    if(TMath::Abs(minvD0-mPDG)<invmasscut || TMath::Abs(minvD0bar-mPDG)<invmasscut){

      cosThetaStarD0 = part->CosThetaStarD0();
      cosThetaStarD0bar = part->CosThetaStarD0bar();
      cosPointingAngle = part->CosPointingAngle();
      normalizedDecayLength2 = part->NormalizedDecayLength2();
      decayLength2 = part->DecayLength2();
      decayLengthxy = part->DecayLengthXY();
      normalizedDecayLengthxy=decayLengthxy/part->DecayLengthXYError();

      ptProng[0]=prong0->Pt(); ptProng[1]=prong1->Pt();
      d0Prong[0]=part->Getd0Prong(0); d0Prong[1]=part->Getd0Prong(1);

      if(fArray==1) cout<<"LS signal: ERROR"<<endl;
      for (Int_t iprong=0; iprong<2; iprong++){
	AliAODTrack *prong=(AliAODTrack*)fDaughterTracks.UncheckedAt(iprong);
	if (fReadMC) labprong[iprong]=prong->GetLabel();
	  
	//cout<<"prong name = "<<prong->GetName()<<" label = "<<prong->GetLabel()<<endl;
	Int_t pdgprong=0;
	if(fReadMC && labprong[iprong]>=0) {
	  mcprong= (AliAODMCParticle*)arrMC->At(labprong[iprong]);
	  pdgprong=mcprong->GetPdgCode();
	}

	Bool_t isPionHere[2]={(isSelectedPID==1 || isSelectedPID==3) ? kTRUE : kFALSE,(isSelectedPID==2 || isSelectedPID==3) ? kTRUE : kFALSE};

	if(TMath::Abs(pdgprong)==211 || isPionHere[iprong]) {
	  //cout<<"pi"<<endl;
	   
	  if(fSys==0){ 
	    fillthispi="hptpiS_";
	    fillthispi+=ptbin;
	    ((TH1F*)listout->FindObject(fillthispi))->Fill(ptProng[iprong]);
	  }
	  fillthispi="hd0piS_";
	  fillthispi+=ptbin;
	  ((TH1F*)listout->FindObject(fillthispi))->Fill(d0Prong[iprong]);
	  if(recalcvtx) {

	    fillthispi="hd0vpiS_";
	    fillthispi+=ptbin;
	    ((TH1F*)listout->FindObject(fillthispi))->Fill(d0[iprong]);
	  }

	}
	  
	if(TMath::Abs(pdgprong)==321 || !isPionHere[iprong]) {
	  //cout<<"kappa"<<endl;
	  if(fSys==0){
	    fillthisK="hptKS_";
	    fillthisK+=ptbin;
	    ((TH1F*)listout->FindObject(fillthisK))->Fill(ptProng[iprong]);
	  }
	  fillthisK="hd0KS_";
	  fillthisK+=ptbin;
	  ((TH1F*)listout->FindObject(fillthisK))->Fill(d0Prong[iprong]);
	  if (recalcvtx){
	    fillthisK="hd0vKS_";
	    fillthisK+=ptbin;
	    ((TH1F*)listout->FindObject(fillthisK))->Fill(d0[iprong]);
	  }
	}

	if(fSys==0){
	  fillthis="hcosthpointd0S_";
	  fillthis+=ptbin;	  
	  ((TH1F*)listout->FindObject(fillthis))->Fill(cosPointingAngle,d0Prong[iprong]);
	}
      } //end loop on prongs

      fillthis="hdcaS_";
      fillthis+=ptbin;	  
      ((TH1F*)listout->FindObject(fillthis))->Fill(part->GetDCA());

      fillthis="hcosthetapointS_";
      fillthis+=ptbin;	  
      ((TH1F*)listout->FindObject(fillthis))->Fill(cosPointingAngle);

      fillthis="hcosthetapointxyS_";
      fillthis+=ptbin;	  
      ((TH1F*)listout->FindObject(fillthis))->Fill(part->CosPointingAngleXY());


      fillthis="hd0d0S_";
      fillthis+=ptbin;
      ((TH1F*)listout->FindObject(fillthis))->Fill(part->Prodd0d0());

      fillthis="hdeclS_";
      fillthis+=ptbin;
      ((TH1F*)listout->FindObject(fillthis))->Fill(decayLength2);

      fillthis="hnormdeclS_";
      fillthis+=ptbin;
      ((TH1F*)listout->FindObject(fillthis))->Fill(normalizedDecayLength2);

      fillthis="hdeclxyS_";
      fillthis+=ptbin;
      ((TH1F*)listout->FindObject(fillthis))->Fill(decayLengthxy);

      fillthis="hnormdeclxyS_";
      fillthis+=ptbin;
      ((TH1F*)listout->FindObject(fillthis))->Fill(normalizedDecayLengthxy);

      fillthis="hdeclxyd0d0S_";
      fillthis+=ptbin;
      ((TH2F*)listout->FindObject(fillthis))->Fill(decayLengthxy,d0Prong[0]*d0Prong[1]);

      fillthis="hnormdeclxyd0d0S_";
      fillthis+=ptbin;
      ((TH2F*)listout->FindObject(fillthis))->Fill(normalizedDecayLengthxy,d0Prong[0]*d0Prong[1]);

      if(recalcvtx) {
	fillthis="hdeclvS_";
	fillthis+=ptbin;
	((TH1F*)listout->FindObject(fillthis))->Fill(decl[0]);

	fillthis="hnormdeclvS_";
	fillthis+=ptbin;
	((TH1F*)listout->FindObject(fillthis))->Fill(decl[1]);

	fillthis="hd0d0vS_";
	fillthis+=ptbin;
	((TH1F*)listout->FindObject(fillthis))->Fill(d0[0]*d0[1]);
      }

      if(fSys==0){
	fillthis="hcosthetastarS_";
	fillthis+=ptbin;
	if ((fReadMC && ((AliAODMCParticle*)arrMC->At(lab))->GetPdgCode() == 421)) ((TH1F*)listout->FindObject(fillthis))->Fill(cosThetaStarD0);
	else {
	  if (fReadMC || isSelectedPID>1)((TH1F*)listout->FindObject(fillthis))->Fill(cosThetaStarD0bar);
	  if(isSelectedPID==1 || isSelectedPID==3)((TH1F*)listout->FindObject(fillthis))->Fill(cosThetaStarD0);
	}

	fillthis="hcosthpointd0d0S_";
	fillthis+=ptbin;	  
	((TH2F*)listout->FindObject(fillthis))->Fill(cosPointingAngle,part->Prodd0d0());
      }

    } //end mass cut
    
  } else{ //Background or LS
    //if(!fReadMC){
    //cout<<"is background"<<endl;
     
    //no mass cut distributions: mass, ptbis
    fillthis="hMassB_";
    fillthis+=ptbin;

    if (!fCutOnDistr || (fCutOnDistr && (fIsSelectedCandidate==1 || fIsSelectedCandidate==3))) ((TH1F*)listout->FindObject(fillthis))->Fill(minvD0);
    if (!fCutOnDistr || (fCutOnDistr && fIsSelectedCandidate>1)) ((TH1F*)listout->FindObject(fillthis))->Fill(minvD0bar);
    if(fSys==0){
      fillthis="hptB1prongnoMcut_";
      fillthis+=ptbin;
      
      ((TH1F*)listout->FindObject(fillthis))->Fill(((AliAODTrack*)fDaughterTracks.UncheckedAt(0))->Pt());
      
      fillthis="hptB2prongsnoMcut_";
      fillthis+=ptbin;
      ((TH1F*)listout->FindObject(fillthis))->Fill(((AliAODTrack*)fDaughterTracks.UncheckedAt(0))->Pt());
      ((TH1F*)listout->FindObject(fillthis))->Fill(((AliAODTrack*)fDaughterTracks.UncheckedAt(0))->Pt());
    }
    //apply cut on invariant mass on the pair
    if(TMath::Abs(minvD0-mPDG)<invmasscut || TMath::Abs(minvD0bar-mPDG)<invmasscut){
      if(fSys==0){
	ptProng[0]=((AliAODTrack*)fDaughterTracks.UncheckedAt(0))->Pt(); ptProng[1]=((AliAODTrack*)fDaughterTracks.UncheckedAt(0))->Pt();
	cosThetaStarD0 = part->CosThetaStarD0();
	cosThetaStarD0bar = part->CosThetaStarD0bar();
      }

      cosPointingAngle = part->CosPointingAngle();
      normalizedDecayLength2 = part->NormalizedDecayLength2();
      decayLength2 = part->DecayLength2();
      decayLengthxy = part->DecayLengthXY();
      normalizedDecayLengthxy=decayLengthxy/part->DecayLengthXYError();
      d0Prong[0]=part->Getd0Prong(0); d0Prong[1]=part->Getd0Prong(1);


      AliAODTrack *prongg=(AliAODTrack*)fDaughterTracks.UncheckedAt(0);
      if(!prongg) {
	if(fDebug>2) cout<<"No daughter found";
	return;
      }
      else{
	if(fArray==1){
	  if(prongg->Charge()==1) {
	    //fTotPosPairs[ptbin]++;
	    ((TH1F*)fOutputMass->FindObject("hpospair"))->Fill(ptbin);
	  } else {
	    //fTotNegPairs[ptbin]++;
	    ((TH1F*)fOutputMass->FindObject("hnegpair"))->Fill(ptbin);
	  }
	}
      }
	
      fillthis="hd0B_";
      fillthis+=ptbin;
      ((TH1F*)listout->FindObject(fillthis))->Fill(d0Prong[0]);
      ((TH1F*)listout->FindObject(fillthis))->Fill(d0Prong[1]);

      if(fReadMC){
	Int_t pdgMother[2]={0,0};
	Double_t factor[2]={1,1};

	for(Int_t iprong=0;iprong<2;iprong++){
	  AliAODTrack *prong=(AliAODTrack*)fDaughterTracks.UncheckedAt(iprong);
	  lab=prong->GetLabel();
	  if(lab>=0){
	    AliAODMCParticle* mcprong=(AliAODMCParticle*)arrMC->At(lab);
	    if(mcprong){
	      Int_t labmom=mcprong->GetMother();
	      if(labmom>=0){
		AliAODMCParticle* mcmother=(AliAODMCParticle*)arrMC->At(labmom);
		if(mcmother) pdgMother[iprong]=mcmother->GetPdgCode();
	      }
	    }
	  }

	  if(fSys==0){

	    fillthis="hd0moresB_";
	    fillthis+=ptbin;
	  
	    if(TMath::Abs(pdgMother[iprong])==310 || TMath::Abs(pdgMother[iprong])==130 || TMath::Abs(pdgMother[iprong])==321){ //K^0_S, K^0_L, K^+-
	      if(ptProng[iprong]<=1)factor[iprong]=1./.7;
	      else factor[iprong]=1./.6;
	      fNentries->Fill(11);
	    }
	    
	    if(TMath::Abs(pdgMother[iprong])==3122) { //Lambda
	      factor[iprong]=1./0.25;
	      fNentries->Fill(12);
	    }
	    fillthis="hd0moresB_";
	    fillthis+=ptbin;

	    ((TH1F*)listout->FindObject(fillthis))->Fill(d0Prong[iprong],factor[iprong]);

	    if(recalcvtx){
	      fillthis="hd0vmoresB_";
	      fillthis+=ptbin;
	      ((TH1F*)listout->FindObject(fillthis))->Fill(d0[iprong],factor[iprong]);
	    }
	  }
	} //loop on prongs

	if(fSys==0){
	  fillthis="hd0d0moresB_";
	  fillthis+=ptbin;
	  ((TH1F*)listout->FindObject(fillthis))->Fill(part->Prodd0d0(),factor[0]*factor[1]);

	  fillthis="hcosthetapointmoresB_";
	  fillthis+=ptbin;
	  ((TH1F*)listout->FindObject(fillthis))->Fill(cosPointingAngle,factor[0]*factor[1]);

	  if(recalcvtx){
	    fillthis="hd0d0vmoresB_";
	    fillthis+=ptbin;
	    ((TH1F*)listout->FindObject(fillthis))->Fill(d0[0]*d0[1],factor[0]*factor[1]);
	  }
	}
      } //readMC

      if(fSys==0){	    
	//normalise pt distr to half afterwards
	fillthis="hptB_";
	fillthis+=ptbin;
	((TH1F*)listout->FindObject(fillthis))->Fill(ptProng[0]);
	((TH1F*)listout->FindObject(fillthis))->Fill(ptProng[1]);

	fillthis="hcosthetastarB_";
	fillthis+=ptbin;
	if (!fCutOnDistr || (fCutOnDistr && (fIsSelectedCandidate==1 || fIsSelectedCandidate==3)))((TH1F*)listout->FindObject(fillthis))->Fill(cosThetaStarD0);
	if (!fCutOnDistr || (fCutOnDistr && fIsSelectedCandidate>1))((TH1F*)listout->FindObject(fillthis))->Fill(cosThetaStarD0bar);	


	fillthis="hd0p0B_";
	fillthis+=ptbin;
	((TH1F*)listout->FindObject(fillthis))->Fill(d0Prong[0]);
	fillthis="hd0p1B_";
	fillthis+=ptbin;
	((TH1F*)listout->FindObject(fillthis))->Fill(d0Prong[1]);
	
	fillthis="hcosthpointd0d0B_";
	fillthis+=ptbin;
	((TH2F*)listout->FindObject(fillthis))->Fill(cosPointingAngle,part->Prodd0d0());
	
	fillthis="hcosthpointd0B_";
	fillthis+=ptbin;	  
	((TH1F*)listout->FindObject(fillthis))->Fill(cosPointingAngle,d0Prong[0]);
	((TH1F*)listout->FindObject(fillthis))->Fill(cosPointingAngle,d0Prong[1]);
	  

	if(recalcvtx){

	  fillthis="hd0vp0B_";
	  fillthis+=ptbin;
	  ((TH1F*)listout->FindObject(fillthis))->Fill(d0[0]);
	  fillthis="hd0vp1B_";
	  fillthis+=ptbin;
	  ((TH1F*)listout->FindObject(fillthis))->Fill(d0[1]);
	  
	  fillthis="hd0vB_";
	  fillthis+=ptbin;
	  ((TH1F*)listout->FindObject(fillthis))->Fill(d0[0]);
	  ((TH1F*)listout->FindObject(fillthis))->Fill(d0[1]);

	}

      }  

      fillthis="hdcaB_";
      fillthis+=ptbin;
      ((TH1F*)listout->FindObject(fillthis))->Fill(part->GetDCA());

      fillthis="hd0d0B_";
      fillthis+=ptbin;
      ((TH1F*)listout->FindObject(fillthis))->Fill(d0Prong[0]*d0Prong[1]);

      if(recalcvtx){
	fillthis="hd0d0vB_";
	fillthis+=ptbin;
	((TH1F*)listout->FindObject(fillthis))->Fill(d0[0]*d0[1]);
      }

      fillthis="hcosthetapointB_";
      fillthis+=ptbin;
      ((TH1F*)listout->FindObject(fillthis))->Fill(cosPointingAngle);

      fillthis="hcosthetapointxyB_";
      fillthis+=ptbin;	  
      ((TH1F*)listout->FindObject(fillthis))->Fill(part->CosPointingAngleXY());

      fillthis="hdeclB_";
      fillthis+=ptbin;
      ((TH1F*)listout->FindObject(fillthis))->Fill(decayLength2);

      fillthis="hnormdeclB_";
      fillthis+=ptbin;
      ((TH1F*)listout->FindObject(fillthis))->Fill(normalizedDecayLength2);

      fillthis="hdeclxyB_";
      fillthis+=ptbin;
      ((TH1F*)listout->FindObject(fillthis))->Fill(decayLengthxy);

      fillthis="hnormdeclxyB_";
      fillthis+=ptbin;
      ((TH1F*)listout->FindObject(fillthis))->Fill(normalizedDecayLengthxy);

      fillthis="hdeclxyd0d0B_";
      fillthis+=ptbin;
      ((TH2F*)listout->FindObject(fillthis))->Fill(decayLengthxy,d0Prong[0]*d0Prong[1]);

      fillthis="hnormdeclxyd0d0B_";
      fillthis+=ptbin;
      ((TH2F*)listout->FindObject(fillthis))->Fill(normalizedDecayLengthxy,d0Prong[0]*d0Prong[1]);	

      if(recalcvtx) {

	fillthis="hdeclvB_";
	fillthis+=ptbin;
	((TH1F*)listout->FindObject(fillthis))->Fill(decl[0]);

	fillthis="hnormdeclvB_";
	fillthis+=ptbin;
	((TH1F*)listout->FindObject(fillthis))->Fill(decl[1]);


      }
    }//mass cut	
  }//else (background)
  
  return;
}
//____________________________________________________________________________

void AliAnalysisTaskSED0Mass::FillMassHists(AliAODRecoDecayHF2Prong *part, TClonesArray *arrMC, AliRDHFCutsD0toKpi* cuts, TList *listout){
  //
  // function used in UserExec to fill mass histograms:
  //


  Double_t mPDG=TDatabasePDG::Instance()->GetParticle(421)->Mass();

  //cout<<"is selected = "<<fIsSelectedCandidate<<endl;

  //cout<<"check cuts = "<<endl;
  //cuts->PrintAll();
  if (!fIsSelectedCandidate){
    //cout<<"Not Selected"<<endl;
    //cout<<"Rejected because "<<cuts->GetWhy()<<endl;
    return;
  }


  if(fDebug>2)  cout<<"Candidate selected"<<endl;

  Double_t invmassD0 = part->InvMassD0(), invmassD0bar = part->InvMassD0bar();
  //printf("SELECTED\n");
  Int_t ptbin=cuts->PtBin(part->Pt());

  // AliAODTrack *prong=(AliAODTrack*)fDaughterTracks.UncheckedAt(0);
  // if(!prong) {
  //   AliDebug(2,"No daughter found");
  //   return;
  // }
  // else{
    // if(prong->Charge()==1) {
    //   ((TH1F*)fDistr->FindObject("hpospair"))->Fill(fCuts->GetNPtBins()+ptbin);
    //   //fTotPosPairs[ptbin]++;
    // } else {
    //   ((TH1F*)fDistr->FindObject("hnegpair"))->Fill(fCuts->GetNPtBins()+ptbin);
    //   //fTotNegPairs[ptbin]++;
    // }
  //  }
 
  // for(Int_t it=0;it<2;it++){
 
  //    //request on spd points to be addes
  //   if(/*nSPD==2 && */part->Pt() > 5. && (TMath::Abs(invmassD0-mPDG)<0.01 || TMath::Abs(invmassD0bar-mPDG)<0.01)){
  //     FILE *f=fopen("4display.txt","a");
  //     fprintf(f,"pt: %f \n Rapidity: %f \t Period Number: %x \t Run Number: %d \t BunchCrossNumb: %d \t OrbitNumber: %d\n",part->Pt(),part->Y(421),aod->GetPeriodNumber(),aod->GetRunNumber(),aod->GetBunchCrossNumber(),aod->GetOrbitNumber());
  //     fclose(f);
  //     //printf("PrimVtx NContributors: %d \n Prongs Rel Angle: %f \n \n",ncont,relangle);
  //   }
  // }
 
  TString fillthis="";
  Int_t pdgDgD0toKpi[2]={321,211};
  Int_t labD0=-1;
  if (fReadMC) labD0 = part->MatchToMC(421,arrMC,2,pdgDgD0toKpi); //return MC particle label if the array corresponds to a D0, -1 if not (cf. AliAODRecoDecay.cxx)

  //count candidates selected by cuts
  fNentries->Fill(1);
  //count true D0 selected by cuts
  if (fReadMC && labD0>=0) fNentries->Fill(2);

  if ((fIsSelectedCandidate==1 || fIsSelectedCandidate==3) && fFillOnlyD0D0bar<2) { //D0
    if(fReadMC){
      if(labD0>=0) {
	if(fArray==1) cout<<"LS signal ERROR"<<endl;

	AliAODMCParticle *partD0 = (AliAODMCParticle*)arrMC->At(labD0);
	Int_t pdgD0 = partD0->GetPdgCode();
	//	cout<<"pdg = "<<pdgD0<<endl;
	if (pdgD0==421){ //D0
	  //	  cout<<"Fill S with D0"<<endl;
	  fillthis="histSgn_";
	  fillthis+=ptbin;
	  ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0);
	  if(fSys==0){
	    if(TMath::Abs(invmassD0 - mPDG) < 0.027 && fFillVarHists){
	      fillthis="histSgn27_";
	      fillthis+=ptbin;
	      ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0);
	    }
	  }
	} else{ //it was a D0bar
	  fillthis="histRfl_";
	  fillthis+=ptbin;
	  ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0);
	}
      } else {//background
	fillthis="histBkg_";
	fillthis+=ptbin;
	((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0);
      }

    }else{
      fillthis="histMass_";
      fillthis+=ptbin;
      //      cout<<"Filling "<<fillthis<<endl;

      //      printf("Fill mass with D0");
      ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0);
    }
     
  }
  if (fIsSelectedCandidate>1 && (fFillOnlyD0D0bar==0 || fFillOnlyD0D0bar==2)) { //D0bar
    if(fReadMC){
      if(labD0>=0) {
	if(fArray==1) cout<<"LS signal ERROR"<<endl;
	AliAODMCParticle *partD0 = (AliAODMCParticle*)arrMC->At(labD0);
	Int_t pdgD0 = partD0->GetPdgCode();
	//	cout<<" pdg = "<<pdgD0<<endl;
	if (pdgD0==-421){ //D0bar
	  fillthis="histSgn_";
	  fillthis+=ptbin;
	  ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0bar);
	  // if (TMath::Abs(invmassD0bar - mPDG) < 0.027){
	  //   fillthis="histSgn27_";
	  //   fillthis+=ptbin;
	  //   ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0bar);
	  // }

	  
	} else{
	  fillthis="histRfl_";
	  fillthis+=ptbin;
	  ((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0bar);
	}
      } else {//background or LS
	fillthis="histBkg_";
	fillthis+=ptbin;
	((TH1F*)(listout->FindObject(fillthis)))->Fill(invmassD0bar);
      }
    }else{
      fillthis="histMass_";
      fillthis+=ptbin;
      //      printf("Fill mass with D0bar");
      ((TH1F*)listout->FindObject(fillthis))->Fill(invmassD0bar);

    }
  }

  return;
}

//__________________________________________________________________________
AliAODVertex* AliAnalysisTaskSED0Mass::GetPrimaryVtxSkipped(AliAODEvent *aodev){
  //Calculate the primary vertex w/o the daughter tracks of the candidate
  
  Int_t skipped[2];
  Int_t nTrksToSkip=2;
  AliAODTrack *dgTrack = (AliAODTrack*)fDaughterTracks.UncheckedAt(0);
  if(!dgTrack){
    AliDebug(2,"no daughter found!");
    return 0x0;
  }
  skipped[0]=dgTrack->GetID();
  dgTrack = (AliAODTrack*)fDaughterTracks.UncheckedAt(1);
  if(!dgTrack){
    AliDebug(2,"no daughter found!");
    return 0x0;
  }
  skipped[1]=dgTrack->GetID();

  AliESDVertex *vertexESD=0x0;
  AliAODVertex *vertexAOD=0x0;
  AliVertexerTracks *vertexer = new AliVertexerTracks(aodev->GetMagneticField());
  
  //
  vertexer->SetSkipTracks(nTrksToSkip,skipped);
  vertexer->SetMinClusters(4);  
  vertexESD = (AliESDVertex*)vertexer->FindPrimaryVertex(aodev); 
  if(!vertexESD) return vertexAOD;
  if(vertexESD->GetNContributors()<=0) { 
    AliDebug(2,"vertexing failed"); 
    delete vertexESD; vertexESD=NULL;
    return vertexAOD;
  }
  
  delete vertexer; vertexer=NULL;
  
  
  // convert to AliAODVertex
  Double_t pos[3],cov[6],chi2perNDF;
  vertexESD->GetXYZ(pos); // position
  vertexESD->GetCovMatrix(cov); //covariance matrix
  chi2perNDF = vertexESD->GetChi2toNDF();
  delete vertexESD; vertexESD=NULL;
  
  vertexAOD = new AliAODVertex(pos,cov,chi2perNDF);
  return vertexAOD;
  
}


//________________________________________________________________________
void AliAnalysisTaskSED0Mass::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  if(fDebug > 1) printf("AnalysisTaskSED0Mass: Terminate() \n");


  fOutputMass = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputMass) {     
    printf("ERROR: fOutputMass not available\n");
    return;
  }
  if(fFillVarHists){
    fDistr = dynamic_cast<TList*> (GetOutputData(2));
    if (!fDistr) {
      printf("ERROR: fDistr not available\n");
      return;
    }
  }

  fNentries = dynamic_cast<TH1F*>(GetOutputData(3));
  
  if(!fNentries){
    printf("ERROR: fNEntries not available\n");
    return;
  }
  fCuts = dynamic_cast<AliRDHFCutsD0toKpi*>(GetOutputData(4));
  if(!fCuts){
    printf("ERROR: fCuts not available\n");
    return;
  }
  fCounter = dynamic_cast<AliNormalizationCounter*>(GetOutputData(5));    
  if (!fCounter) {
    printf("ERROR: fCounter not available\n");
    return;
  }
  Int_t nptbins=fCuts->GetNPtBins();
  for(Int_t ipt=0;ipt<nptbins;ipt++){ 

    if(fArray==1 && fFillVarHists){ 
      fLsNormalization = 2.*TMath::Sqrt(((TH1F*)fOutputMass->FindObject("hpospair"))->Integral(nptbins+ipt+1,nptbins+ipt+2)*((TH1F*)fOutputMass->FindObject("hnegpair"))->Integral(nptbins+ipt+1,nptbins+ipt+2)); //after cuts


      if(fLsNormalization>1e-6) {
	
	TString massName="histMass_";
	massName+=ipt;
	((TH1F*)fOutputMass->FindObject(massName))->Scale((1/fLsNormalization)*((TH1F*)fOutputMass->FindObject(massName))->GetEntries());

      }
    

      fLsNormalization = 2.*TMath::Sqrt(((TH1F*)fOutputMass->FindObject("hpospair"))->Integral(ipt+1,ipt+2)*((TH1F*)fOutputMass->FindObject("hnegpair"))->Integral(ipt+1,ipt+2)); 
      //fLsNormalization = 2.*TMath::Sqrt(fTotPosPairs[4]*fTotNegPairs[4]);

      if(fLsNormalization>1e-6) {

	TString nameDistr="hptB_";
	nameDistr+=ipt;
	((TH1F*)fDistr->FindObject(nameDistr))->Scale((1/fLsNormalization)*((TH1F*)fDistr->FindObject(nameDistr))->GetEntries());
	nameDistr="hdcaB_";
	nameDistr+=ipt;
	((TH1F*)fDistr->FindObject(nameDistr))->Scale((1/fLsNormalization)*((TH1F*)fDistr->FindObject(nameDistr))->GetEntries());
	nameDistr="hcosthetastarB_";
	nameDistr+=ipt;
	((TH1F*)fDistr->FindObject(nameDistr))->Scale((1/fLsNormalization)*((TH1F*)fDistr->FindObject(nameDistr))->GetEntries());
	nameDistr="hd0B_";
	nameDistr+=ipt;
	((TH1F*)fDistr->FindObject(nameDistr))->Scale((1/fLsNormalization)*((TH1F*)fDistr->FindObject(nameDistr))->GetEntries());
	nameDistr="hd0d0B_";
	nameDistr+=ipt;
	((TH1F*)fDistr->FindObject(nameDistr))->Scale((1/fLsNormalization)*((TH1F*)fDistr->FindObject(nameDistr))->GetEntries());
	nameDistr="hcosthetapointB_";
	nameDistr+=ipt;
	((TH1F*)fDistr->FindObject(nameDistr))->Scale((1/fLsNormalization)*((TH1F*)fDistr->FindObject(nameDistr))->GetEntries());
	if(fSys==0){
	  nameDistr="hcosthpointd0d0B_";
	  nameDistr+=ipt;
	  ((TH2F*)fDistr->FindObject(nameDistr))->Scale((1/fLsNormalization)*((TH2F*)fDistr->FindObject(nameDistr))->GetEntries());
	}
      }
    }
  }
  TString cvname,cstname;

  if (fArray==0){
    cvname="D0invmass";
    cstname="cstat0";
  } else {
    cvname="LSinvmass";
    cstname="cstat1";
  }

  TCanvas *cMass=new TCanvas(cvname,cvname);
  cMass->cd();
  ((TH1F*)fOutputMass->FindObject("histMass_3"))->Draw();

  TCanvas* cStat=new TCanvas(cstname,Form("Stat%s",fArray ? "LS" : "D0"));
  cStat->cd();
  cStat->SetGridy();
  fNentries->Draw("htext0");

  // TCanvas *ccheck=new TCanvas(Form("cc%d",fArray),Form("cc%d",fArray));
  // ccheck->cd();

  return;
}

