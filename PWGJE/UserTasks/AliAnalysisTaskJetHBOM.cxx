// **************************************
// Task used for the correction of detector effects for background fluctuations in jet spectra by the HBOM method
// *******************************************


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

 
#include <TROOT.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include <TChain.h>
#include <TRefArray.h>
#include <TFile.h>
#include <TKey.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TProfile.h>
#include <TF1.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>
#include  "TDatabasePDG.h"

#include "AliAnalysisTaskJetHBOM.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAODExtension.h"
#include "AliAODTrack.h"
#include "AliAODJet.h"
#include "AliAODMCParticle.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliInputEventHandler.h"
#include "AliAODJetEventBackground.h"


using std::cout;
using std::endl;
using std::vector;

ClassImp(AliAnalysisTaskJetHBOM)

AliAnalysisTaskJetHBOM::~AliAnalysisTaskJetHBOM(){
  //
  // Destructor
  //

  delete fRef;
  fRef = 0;
  delete fRandom;
  fRandom = 0;

  if(fTCARandomConesOut)fTCARandomConesOut->Delete();
  delete fTCARandomConesOut;
  fTCARandomConesOut = 0;

}

AliAnalysisTaskJetHBOM::AliAnalysisTaskJetHBOM(): 
  AliAnalysisTaskSE(),
  fAOD(0x0),
  fAODExtension(0x0),
  fRef(0),
  fUseAODTrackInput(kFALSE),
  fUseAODMCInput(kFALSE),
  fEventSelection(kFALSE),     
  fFilterMask(0),
  fFilterMaskBestPt(0),
  fFilterType(0),
  fJetTypes(kJet),
  fTrackTypeRec(kTrackUndef),
  fTrackTypeGen(kTrackUndef),  
  fNSkipLeadingRan(0),
  fNSkipLeadingCone(0),
  fNRandomCones(0),
  randCone_pos(0),
  randCone_Eta(0),
  randCone_Phi(0),
  fNHBOM(0),
  fTrackEtaWindow(0.9),    
  fTrackPtCut(0.),							
  fJetOutputMinPt(0.150),
  fMaxTrackPtInJet(100.),
  fVtxZCut(8),
  fVtxR2Cut(1),
  fCentCutUp(0),
  fCentCutLo(0),
  fNonStdBranch(""),
  fBackgroundBranch(""),
  fNonStdFile(""),
  fMomResH1(0x0),
  fMomResH2(0x0),
  fMomResH3(0x0),
  fMomResH1Fit(0x0),
  fMomResH2Fit(0x0),
  fMomResH3Fit(0x0),
  fhEffH1(0x0),
  fhEffH2(0x0),
  fhEffH3(0x0),
  fUseTrMomentumSmearing(kFALSE),
  fUseDiceEfficiency(kFALSE),
  fRparam(1.0), 
  fAlgorithm(fastjet::kt_algorithm),
  fStrategy(fastjet::Best),
  fRecombScheme(fastjet::BIpt_scheme),
  fAreaType(fastjet::active_area), 
  fGhostArea(0.01),
  fActiveAreaRepeats(1),
  fGhostEtamax(1.5),
  background(0),
  fTCARandomConesOut(0x0),
  fRandom(0),
  fh1Xsec(0x0),
  fh1Trials(0x0),
  fh1PtHard(0x0),
  fh1PtHardNoW(0x0),  
  fh1PtHardTrials(0x0),
  fh1Nch(0x0),
  fh1CentralityPhySel(0x0), 
  fh1Centrality(0x0), 
  fh1DeltapT(0x0),
  fh1Rho(0x0),
  fh1RhoSigma(0x0),
  fh1PtRandCone(0x0),
  fh1efficiencyPt(0x0),
  fh2efficiencyPhi(0x0),
  fh1ZPhySel(0x0), 
  fh1Z(0x0), 
  fHistList(0x0)  
{
  
}

AliAnalysisTaskJetHBOM::AliAnalysisTaskJetHBOM(const char* name):
  AliAnalysisTaskSE(name),
  fAOD(0x0),
  fAODExtension(0x0),
  fRef(0),
  fUseAODTrackInput(kFALSE),
  fUseAODMCInput(kFALSE),
  fEventSelection(kFALSE),							  
  fFilterMask(0),
  fFilterMaskBestPt(0),
  fFilterType(0),
  fJetTypes(kJet),
  fTrackTypeRec(kTrackUndef),
  fTrackTypeGen(kTrackUndef),
  fNSkipLeadingRan(0),
  fNSkipLeadingCone(0),
  fNRandomCones(0),
  randCone_pos(0),
  randCone_Eta(0),
  randCone_Phi(0),
  fNHBOM(0),
  fTrackEtaWindow(0.9),    
  fTrackPtCut(0.),							
  fJetOutputMinPt(0.150),
  fMaxTrackPtInJet(100.),
  fVtxZCut(8),
  fVtxR2Cut(1),
  fCentCutUp(0),
  fCentCutLo(0),
  fNonStdBranch(""),
  fBackgroundBranch(""),
  fNonStdFile(""),
  fMomResH1(0x0),
  fMomResH2(0x0),
  fMomResH3(0x0),
  fMomResH1Fit(0x0),
  fMomResH2Fit(0x0),
  fMomResH3Fit(0x0),
  fhEffH1(0x0),
  fhEffH2(0x0),
  fhEffH3(0x0),
  fUseTrMomentumSmearing(kFALSE),
  fUseDiceEfficiency(kFALSE),
  fRparam(1.0), 
  fAlgorithm(fastjet::kt_algorithm),
  fStrategy(fastjet::Best),
  fRecombScheme(fastjet::BIpt_scheme),
  fAreaType(fastjet::active_area), 
  fGhostArea(0.01),
  fActiveAreaRepeats(1),
  fGhostEtamax(1.5),
  background(0),
  fTCARandomConesOut(0x0),
  fRandom(0),
  fh1Xsec(0x0),
  fh1Trials(0x0),
  fh1PtHard(0x0),
  fh1PtHardNoW(0x0),  
  fh1PtHardTrials(0x0),
  fh1Nch(0x0),
  fh1CentralityPhySel(0x0), 
  fh1Centrality(0x0), 
  fh1DeltapT(0x0),
  fh1Rho(0x0),
  fh1RhoSigma(0x0),
  fh1PtRandCone(0x0),
  fh1efficiencyPt(0x0),
  fh2efficiencyPhi(0x0),
  fh1ZPhySel(0x0), 
  fh1Z(0x0), 
  fHistList(0x0)
{
 DefineOutput(1, TList::Class());  
}



Bool_t AliAnalysisTaskJetHBOM::Notify()
{
  //
  // Implemented Notify() to read the cross sections
  // and number of trials from pyxsec.root
  // 
  return kTRUE;
}

void AliAnalysisTaskJetHBOM::UserCreateOutputObjects()
{

  //
  // Create the output container
  //

  fRandom = new TRandom3(0);


  // Connect the AOD


  if (fDebug > 1) printf("AnalysisTaskJetHBOM::UserCreateOutputObjects() \n");

  

  if(fNonStdBranch.Length()!=0)
    {
      // only create the output branch if we have a name
      // Create a new branch for jets...
      //  -> cleared in the UserExec....
      // here we can also have the case that the brnaches are written to a separate file
      
      // create the branch for the random cones with the same R 
      TString cName = Form("%sRandomConeSkip%02d",fNonStdBranch.Data(),fNSkipLeadingCone);
      if(fUseDiceEfficiency || fUseTrMomentumSmearing)
	cName = Form("%sDetector%d%d_RandomConeSkip%02d",fNonStdBranch.Data(),fUseTrMomentumSmearing,fUseDiceEfficiency,fNSkipLeadingCone);
    
      //create array for the random cones; Until now only one cone per event is used
      if(!AODEvent()->FindListObject(cName.Data())){
	fTCARandomConesOut = new TClonesArray("AliAODJet", 0);
	fTCARandomConesOut->SetName(cName.Data());
	AddAODBranch("TClonesArray",&fTCARandomConesOut,fNonStdFile.Data());
      }
      
      if(fNonStdFile.Length()!=0){
	// 
	// case that we have an AOD extension we need to fetch the jets from the extended output
	// we identify the extension aod event by looking for the branchname
	AliAODHandler *aodH = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());
	// case that we have an AOD extension we need can fetch the background maybe from the extended output                                                                  
	fAODExtension = (aodH?aodH->GetExtension(fNonStdFile.Data()):0);
      }
    }

  //  FitMomentumResolution();


  if(!fHistList)fHistList = new TList();
  fHistList->SetOwner();
  PostData(1, fHistList); // post data in any case once

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  //
  //  Histogram
    
  const Int_t nBinPt = 100;
  Double_t binLimitsPt[nBinPt+1];
  for(Int_t iPt = 0;iPt <= nBinPt;iPt++){
    if(iPt == 0){
      binLimitsPt[iPt] = 0.0;
    }
    else {// 1.0
      binLimitsPt[iPt] =  binLimitsPt[iPt-1] + 2.0;
    }
  }
  
  const Int_t nBinPhi = 90;
  Double_t binLimitsPhi[nBinPhi+1];
  for(Int_t iPhi = 0;iPhi<=nBinPhi;iPhi++){
    if(iPhi==0){
      binLimitsPhi[iPhi] = -1.*TMath::Pi();
    }
    else{
      binLimitsPhi[iPhi] = binLimitsPhi[iPhi-1] + 1/(Float_t)nBinPhi * TMath::Pi()*2;
    }
  }



  const Int_t nBinEta = 40;
  Double_t binLimitsEta[nBinEta+1];
  for(Int_t iEta = 0;iEta<=nBinEta;iEta++){
    if(iEta==0){
      binLimitsEta[iEta] = -2.0;
    }
    else{
      binLimitsEta[iEta] = binLimitsEta[iEta-1] + 0.1;
    }
  }

  const int nChMax = 5000;

  fh1Xsec = new TProfile("fh1Xsec","xsec from pyxsec.root",1,0,1);
  fh1Xsec->GetXaxis()->SetBinLabel(1,"<#sigma>");

  fh1Trials = new TH1F("fh1Trials","trials root file",1,0,1);
  fh1Trials->GetXaxis()->SetBinLabel(1,"#sum{ntrials}");


  fh1PtHard = new TH1F("fh1PtHard","PYTHIA Pt hard;p_{T,hard}",nBinPt,binLimitsPt);
  fh1PtHardNoW = new TH1F("fh1PtHardNoW","PYTHIA Pt hard no weight;p_{T,hard}",nBinPt,binLimitsPt);
  fh1PtHardTrials = new TH1F("fh1PtHardTrials","PYTHIA Pt hard weight with trials;p_{T,hard}",nBinPt,binLimitsPt);

  fh1Nch = new TH1F("fh1Nch","charged multiplicity; N_{ch}",nChMax,-0.5,nChMax-0.5);

  fh1Centrality = new TH1F("fh1Centrality",";cent (%)",111,-0.5,110.5);
  fh1CentralityPhySel = new TH1F("fh1CentralityPhySel",";cent (%)",111,-0.5,110.5);

  fh1Z = new TH1F("fh1Z",";zvtx",100,-25,25);
  fh1ZPhySel = new TH1F("fh1ZPhySel",";zvtx",100,-25,25);
  
  fh1DeltapT = new TH1F("fh1DeltapT","DeltapT",100,-50,50);
  fh1Rho = new TH1F("fh1Rho","Rho",100,0,200);
  fh1RhoSigma = new TH1F("fh1RhoSigma","SigmaRho",40,0,40);
  fh1PtRandCone = new TH1F("fh1PtRandCone","pt",100,0,200);

  const Int_t saveLevel = 3; // large save level more histos
  if(saveLevel>0){
    fHistList->Add(fh1Xsec);
    fHistList->Add(fh1Trials);

    fHistList->Add(fh1Nch);
    fHistList->Add(fh1Centrality);
    fHistList->Add(fh1CentralityPhySel);
    fHistList->Add(fh1Z);
    fHistList->Add(fh1ZPhySel);
    fHistList->Add(fh1DeltapT);
    fHistList->Add(fh1Rho);
    fHistList->Add(fh1RhoSigma);
    fHistList->Add(fh1PtRandCone);
  }

  // =========== Switch on Sumw2 for all histos ===========
  for (Int_t i=0; i<fHistList->GetEntries(); ++i) {
    TH1 *h1 = dynamic_cast<TH1*>(fHistList->At(i));
    if (h1){
      h1->Sumw2();
      continue;
    }
    THnSparse *hn = dynamic_cast<THnSparse*>(fHistList->At(i));
    if(hn)hn->Sumw2();
  }
  TH1::AddDirectory(oldStatus);
}

void AliAnalysisTaskJetHBOM::Init()
{
  //
  // Initialization
  //

  if (fDebug > 1) printf("AnalysisTaskJetHBOM::Init() \n");

  FitMomentumResolution();

}

void AliAnalysisTaskJetHBOM::UserExec(Option_t */*option*/)
{

  // handle and reset the output jet branch 

  if(fTCARandomConesOut)fTCARandomConesOut->Delete();

  //
  // Execute analysis for current event
  //
  AliESDEvent *fESD = 0;
  if(fUseAODTrackInput){    
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD){
      Printf("%s:%d AODEvent not found in Input Manager %d",(char*)__FILE__,__LINE__,fUseAODTrackInput);
      return;
    }
    // fetch the header
  }
  else{
    //  assume that the AOD is in the general output...
    fAOD  = AODEvent();
    if(!fAOD){
      Printf("%s:%d AODEvent not found in the Output",(char*)__FILE__,__LINE__);
      return;
    }
    if(fDebug>0){
      fESD = dynamic_cast<AliESDEvent*> (InputEvent());
    }
  }

  //Check if information is provided detector level effects
  if(!fMomResH1 || !fMomResH2 || !fMomResH3) fUseTrMomentumSmearing = kFALSE;
  if(!fhEffH1 || !fhEffH2 || !fhEffH3)       fUseDiceEfficiency = kFALSE;
  
  Bool_t selectEvent =  false;
  Bool_t physicsSelection = true;// handled by the framework(fInputHandler->IsEventSelected()&AliVEvent::kMB)==AliVEvent::kMB;

  Float_t cent = 0;
  Float_t zVtx  = 0;

  if(fAOD){
    const AliAODVertex *vtxAOD = fAOD->GetPrimaryVertex();
    TString vtxTitle(vtxAOD->GetTitle());
    zVtx = vtxAOD->GetZ();

    cent = fAOD->GetHeader()->GetCentrality();
    if(physicsSelection){
      fh1CentralityPhySel->Fill(cent);
      fh1ZPhySel->Fill(zVtx);
    }
    // zVertex and centrality selection
    if(fEventSelection){
      if(vtxAOD->GetNContributors()>2&&!vtxTitle.Contains("TPCVertex")){
	Float_t yvtx = vtxAOD->GetY();
	Float_t xvtx = vtxAOD->GetX();
	Float_t r2   = yvtx*yvtx+xvtx*xvtx;  
	if(TMath::Abs(zVtx)<fVtxZCut&&r2<fVtxR2Cut){//usual fVtxZCut=10 and fVtxR2Cut=1  // apply vertex cut later on 
	  if(physicsSelection){
	    selectEvent = true;
	  }
	}
      }
      //centrality selection
      if(fCentCutUp>0){
	if(cent<fCentCutLo||cent>fCentCutUp){
	  selectEvent = false;
	}
      }
    }else{
      selectEvent = true;
    }
  }
  

  if(!selectEvent){
    PostData(1, fHistList);
    return;
  }
  fh1Centrality->Fill(cent);  
  fh1Z->Fill(zVtx);
  fh1Trials->Fill("#sum{ntrials}",1);
  

  if (fDebug > 10)Printf("%s:%d",(char*)__FILE__,__LINE__);

  // ==== General variables needed



  // we simply fetch the tracks/mc particles as a list of AliVParticles

  //reconstructed particles
  TList recParticles;
  Int_t nT = GetListOfTracks(&recParticles,fTrackTypeRec);
  Float_t nCh = recParticles.GetEntries(); 
  fh1Nch->Fill(nCh);
  if(fDebug>2)Printf("%s:%d Selected Rec tracks: %d %d",(char*)__FILE__,__LINE__,nT,recParticles.GetEntries());
  //nT = GetListOfTracks(&genParticles,fTrackTypeGen);
  //if(fDebug>2)Printf("%s:%d Selected Gen tracks: %d %d",(char*)__FILE__,__LINE__,nT,genParticles.GetEntries());


  //apply efficiency fNHBOM times
  if(fNHBOM>0){
    for(int particle=0;particle<recParticles.GetEntries();particle++){
      // hier Effizienzen laden und überprüfen ob das Teilchen nachgewiesen wird.
      AliVParticle *vp = (AliVParticle*)recParticles.At(particle);
      
      Double_t pT = vp->Pt();
      Double_t phi = vp->Phi();
      
      //load efficiency
      Double_t efficiencyPt = fh1efficiencyPt->GetBinContent(fh1efficiencyPt->FindBin(pT));
      Double_t efficiencyPhi = fh2efficiencyPhi->GetBinContent(fh2efficiencyPhi->FindBin(phi,pT));
      Double_t eff = efficiencyPt*efficiencyPhi; //over all efficiency
      
      // if ran<eff -> particle is detected eff^fNHBOM = efficiency to detect particle fNHBOM times
      Double_t ran = fRandom->Rndm();
      if(ran>TMath::Power(eff,fNHBOM)){
	recParticles.Remove(vp);
      }
    }
  }


  // find the jets....

  vector<fastjet::PseudoJet> inputParticlesRec;
  vector<fastjet::PseudoJet> inputParticlesRecRan;
  
  //randomize particles
  AliAODJet vTmpRan(1,0,0,1);
  for(int i = 0; i < recParticles.GetEntries(); i++){
    AliVParticle *vp = (AliVParticle*)recParticles.At(i);

    // Carefull energy is not well determined in real data, should not matter for p_T scheme?
    // we take total momentum here

      //Add particles to fastjet in case we are not running toy model
      fastjet::PseudoJet jInp(vp->Px(),vp->Py(),vp->Pz(),vp->P());
      jInp.set_user_index(i);
      inputParticlesRec.push_back(jInp);

    // the randomized input changes eta and phi, but keeps the p_T
      Double_t pT = vp->Pt();
      Double_t eta = 2.*fTrackEtaWindow * fRandom->Rndm() - fTrackEtaWindow;
      Double_t phi = 2.* TMath::Pi() * fRandom->Rndm();
      
      Double_t theta = 2.*TMath::ATan(TMath::Exp(-eta));  
      Double_t pZ = pT/TMath::Tan(theta);

      Double_t pX = pT * TMath::Cos(phi);
      Double_t pY = pT * TMath::Sin(phi);
      Double_t p  = TMath::Sqrt(pT*pT+pZ*pZ); 
      fastjet::PseudoJet jInpRan(pX,pY,pZ,p);

      jInpRan.set_user_index(i);
      inputParticlesRecRan.push_back(jInpRan);
      vTmpRan.SetPxPyPzE(pX,pY,pZ,p);

  }// recparticles

  if(inputParticlesRec.size()==0){
    if(fDebug)Printf("%s:%d No input particles found, skipping event",(char*)__FILE__,__LINE__);
    PostData(1, fHistList);
    return;
  }
  
  // run fast jet
  // employ setters for these...

  fastjet::GhostedAreaSpec ghostSpec(fGhostEtamax, fActiveAreaRepeats, fGhostArea);
  fastjet::AreaType areaType =   fastjet::active_area;
  fastjet::AreaDefinition areaDef = fastjet::AreaDefinition(areaType,ghostSpec);
  fastjet::JetDefinition jetDef(fAlgorithm, fRparam, fRecombScheme, fStrategy);
  fastjet::ClusterSequenceArea clustSeq(inputParticlesRec, jetDef,areaDef);
  
  //range where to compute background
  Double_t phiMin = 0, phiMax = 0, rapMin = 0, rapMax = 0;
  phiMin = 0;
  phiMax = 2*TMath::Pi();
  rapMax = fGhostEtamax - fRparam;
  rapMin = - fGhostEtamax + fRparam;
  fastjet::RangeDefinition range(rapMin,rapMax, phiMin, phiMax);
 

  const vector <fastjet::PseudoJet> &inclusiveJets = clustSeq.inclusive_jets();
  const vector <fastjet::PseudoJet> &sortedJets = sorted_by_pt(inclusiveJets);

 
 // loop over all jets and fill information, first one is the leading jet

  if(inclusiveJets.size()>0){
    
    //background estimates:all bckg jets wo the 2 hardest
    vector<fastjet::PseudoJet> jets2=sortedJets;
    if(jets2.size()>2) jets2.erase(jets2.begin(),jets2.begin()+2); //removes the two jets with the highest pT; +2 is correct ro remove 2 jets
    Double_t bkg1=0;
    Double_t sigma1=0.;
    Double_t meanarea1=0.;
    clustSeq.get_median_rho_and_sigma(jets2, range, true, bkg1, sigma1, meanarea1, true);
    fh1RhoSigma->Fill(sigma1);// fluctuation of the background
    background = bkg1;//sets background variable of the task to the correct value


    // generate random cones
    if(fTCARandomConesOut){
      // create a random jet within the acceptance
      Double_t etaMax = fTrackEtaWindow - fRparam;//0.9 - 0.4
      Int_t nCone = 0;
      Double_t pTC = 1; // small number
      Double_t etaC = etaMax*2.*(fRandom->Rndm()-0.5); // +- etamax
      Double_t phiC = fRandom->Rndm()*2.*TMath::Pi(); // 0 - 2pi
      // use fixed position for random Cones
      if(randCone_pos){
	etaC = randCone_Eta;
	phiC = randCone_Phi;
      }
      // massless jet
      Double_t thetaC = 2.*TMath::ATan(TMath::Exp(-etaC));  
      Double_t pZC = pTC/TMath::Tan(thetaC);
      Double_t pXC = pTC * TMath::Cos(phiC);
      Double_t pYC = pTC * TMath::Sin(phiC);
      Double_t pC  = TMath::Sqrt(pTC*pTC+pZC*pZC); 
      AliAODJet tmpRecC (pXC,pYC,pZC, pC); 
      
      tmpRecC.SetBgEnergy(0,0); // this is use as temporary storage of the summed p_T below
      if(fTCARandomConesOut)new ((*fTCARandomConesOut)[nCone++]) AliAODJet(tmpRecC);
  
      // loop over the reconstructed particles and add up the pT in the random cones
      // maybe better to loop over randomized particles not in the real jets...
      // but this by definition brings dow average energy in the whole  event
      AliAODJet vTmpRanR(1,0,0,1);
      for(int i = 0; i < recParticles.GetEntries(); i++){
	AliVParticle *vp = (AliVParticle*)recParticles.At(i);
	//add up energy in cone
	if(fTCARandomConesOut){
	  AliAODJet *jC = (AliAODJet*)fTCARandomConesOut->At(0);  
	  if(jC&&jC->DeltaR(vp)<fRparam){
	    if(vp->Pt()>fMaxTrackPtInJet)jC->SetTrigger(AliAODJet::kHighTrackPtTriggered);
	    jC->SetBgEnergy(jC->ChargedBgEnergy()+vp->Pt(),0);
	  }
	}// add up energy in cone
	
      }// loop over recparticles
    }  //fTCARandomConesOut
    Float_t jetArea = fRparam*fRparam*TMath::Pi();
    if(fTCARandomConesOut){
      // rescale the momentum vectors for the random cones
      AliAODJet *rC = (AliAODJet*)fTCARandomConesOut->At(0);
      if(rC){
	Double_t etaC = rC->Eta();
	Double_t phiC = rC->Phi();
	// massless jet, unit vector
	Double_t pTC = rC->ChargedBgEnergy();
	if(pTC<=0)pTC = 0.001; // for almost empty events
	Double_t thetaC = 2.*TMath::ATan(TMath::Exp(-etaC));  
	Double_t pZC = pTC/TMath::Tan(thetaC);
	Double_t pXC = pTC * TMath::Cos(phiC);
	Double_t pYC = pTC * TMath::Sin(phiC);
	Double_t pC  = TMath::Sqrt(pTC*pTC+pZC*pZC); 
	rC->SetPxPyPzE(pXC,pYC,pZC, pC); 
	rC->SetBgEnergy(0,0);
	rC->SetEffArea(jetArea,0);
      }
    }
  }//inclusive Jets > 0

 //Calculate delta pT
 AliAODJet *randCone = (AliAODJet*)fTCARandomConesOut->At(0);
 if(randCone){
   //background is the backbround density per area and area=pi*0.4^2 -> backgroundCone is the background energy under the cone
   Float_t backgroundCone = background * randCone->EffectiveAreaCharged();
   //calculates difference between expected and measured energy density
   Float_t ptSub = randCone->Pt() - backgroundCone;
   fh1DeltapT->Fill(ptSub);// delta pT
   fh1Rho->Fill(background);// background rho
   fh1PtRandCone->Fill(randCone->Pt());// pT of random cone
 }else{
   if(fDebug)Printf("%s:%d No random Cone found",(char*)__FILE__,__LINE__);
 }
 

 if (fDebug > 2){
   if(fTCARandomConesOut)Printf("%s:%d RC %d",(char*)__FILE__,__LINE__,fTCARandomConesOut->GetEntriesFast());
 }
 PostData(1, fHistList);
}

void AliAnalysisTaskJetHBOM::Terminate(Option_t */*option*/)
{
  //
  // Terminate analysis
  //
    if (fDebug > 1) printf("AnalysisJetHBOM: Terminate() \n");

    if(fMomResH1Fit) delete fMomResH1Fit;
    if(fMomResH2Fit) delete fMomResH2Fit;
    if(fMomResH3Fit) delete fMomResH3Fit;
    
}


Int_t  AliAnalysisTaskJetHBOM::GetListOfTracks(TList *list,Int_t type){

  //
  // get list of tracks/particles for different types
  //

  if(fDebug>2)Printf("%s:%d Selecting tracks with %d",(char*)__FILE__,__LINE__,type);

  Int_t iCount = 0;
   if(type==kTrackAOD || type==kTrackAODextra || type==kTrackAODextraonly){
    if(type!=kTrackAODextraonly) {
      AliAODEvent *aod = 0;
      if(fUseAODTrackInput)aod = dynamic_cast<AliAODEvent*>(InputEvent());
      else aod = AODEvent();
      if(!aod){
	if(fDebug>2)Printf("%s:%d No AOD",(char*)__FILE__,__LINE__);
	return iCount;
      }

      for(int it = 0;it < aod->GetNumberOfTracks();++it){
	AliAODTrack *tr = aod->GetTrack(it);
	Bool_t bGood = false;
	if(fFilterType == 0)bGood = true;
	else if(fFilterType == 1)bGood = tr->IsHybridTPCConstrainedGlobal();
	else if(fFilterType == 2)bGood = tr->IsHybridGlobalConstrainedGlobal();
	if((fFilterMask>0)&&((!tr->TestFilterBit(fFilterMask)||(!bGood)))){
	  if(fDebug>10)Printf("%s:%d Not matching filter %d/%d %d/%d",(char*)__FILE__,__LINE__,it,aod->GetNumberOfTracks(),fFilterMask,tr->GetFilterMap());	
	  continue;
	}
	if(TMath::Abs(tr->Eta())>fTrackEtaWindow){
	  if(fDebug>10)Printf("%s:%d Not matching eta %d/%d",(char*)__FILE__,__LINE__,it,aod->GetNumberOfTracks());	
	  continue;
	}
	if(tr->Pt()<fTrackPtCut){
	  if(fDebug>10)Printf("%s:%d Not matching pt %d/%d",(char*)__FILE__,__LINE__,it,aod->GetNumberOfTracks());	
	  continue;
	}
	if(fDebug>10)Printf("%s:%d MATCHED %d/%d",(char*)__FILE__,__LINE__,it,aod->GetNumberOfTracks());	
	list->Add(tr);
	iCount++;
      }
    }
    if(type==kTrackAODextra || type==kTrackAODextraonly) {
      AliAODEvent *aod = 0;
      if(fUseAODTrackInput)aod = dynamic_cast<AliAODEvent*>(InputEvent());
      else aod = AODEvent();
      
      if(!aod){
	return iCount;
      }
      TClonesArray *aodExtraTracks = dynamic_cast<TClonesArray*>(aod->FindListObject("aodExtraTracks"));
      if(!aodExtraTracks)return iCount;
      for(int it =0; it<aodExtraTracks->GetEntries(); it++) {
	AliVParticle *track = dynamic_cast<AliVParticle*> ((*aodExtraTracks)[it]);
	if (!track) continue;

	AliAODTrack *trackAOD = dynamic_cast<AliAODTrack*> (track);
	if(!trackAOD)continue;
	Bool_t bGood = false;
	if(fFilterType == 0)bGood = true;
	else if(fFilterType == 1)bGood = trackAOD->IsHybridTPCConstrainedGlobal();
	else if(fFilterType == 2)bGood = trackAOD->IsHybridGlobalConstrainedGlobal();
	if((fFilterMask>0)&&((!trackAOD->TestFilterBit(fFilterMask)||(!bGood))))continue;
        if(TMath::Abs(trackAOD->Eta())>fTrackEtaWindow) continue;
	if(trackAOD->Pt()<fTrackPtCut) continue;
	list->Add(trackAOD);
	iCount++;
      }
    }
  }
  else if (type ==  kTrackKineAll||type == kTrackKineCharged){
    AliMCEvent* mcEvent = MCEvent();
    if(!mcEvent)return iCount;
    // we want to have alivpartilces so use get track
    for(int it = 0;it < mcEvent->GetNumberOfTracks();++it){
      if(!mcEvent->IsPhysicalPrimary(it))continue;
      AliMCParticle* part = (AliMCParticle*)mcEvent->GetTrack(it);
      if(type == kTrackKineAll){
	if(part->Pt()<fTrackPtCut)continue;
	list->Add(part);
	iCount++;
      }
      else if(type == kTrackKineCharged){
	if(part->Particle()->GetPDG()->Charge()==0)continue;
	if(part->Pt()<fTrackPtCut)continue;
	list->Add(part);
	iCount++;
      }
    }
  }
  else if (type == kTrackAODMCCharged || type == kTrackAODMCAll || type == kTrackAODMCChargedAcceptance) {
    AliAODEvent *aod = 0;
    if(fUseAODMCInput)aod = dynamic_cast<AliAODEvent*>(InputEvent());
    else aod = AODEvent();
    if(!aod)return iCount;
    TClonesArray *tca = dynamic_cast<TClonesArray*>(aod->FindListObject(AliAODMCParticle::StdBranchName()));
    if(!tca)return iCount;
    for(int it = 0;it < tca->GetEntriesFast();++it){
      AliAODMCParticle *part = (AliAODMCParticle*)(tca->At(it));
      if(!part->IsPhysicalPrimary())continue;
      if(type == kTrackAODMCAll){
	if(part->Pt()<fTrackPtCut)continue;
	list->Add(part);
	iCount++;
      }
      else if (type == kTrackAODMCCharged || type == kTrackAODMCChargedAcceptance ){
	if(part->Charge()==0)continue;
	if(part->Pt()<fTrackPtCut)continue;
	if(kTrackAODMCCharged){
	  list->Add(part);
	}
	else {
	  if(TMath::Abs(part->Eta())>fTrackEtaWindow)continue;
	  list->Add(part);
	}
	iCount++;
      }
    }
  }// AODMCparticle
  list->Sort();
  return iCount;
}

void AliAnalysisTaskJetHBOM::SetMomentumResolutionHybrid(TProfile *p1, TProfile *p2, TProfile *p3) {

  //
  // set mom res profiles
  //

  fMomResH1 = (TProfile*)p1->Clone("fMomResH1");
  fMomResH2 = (TProfile*)p2->Clone("fMomResH2");
  fMomResH3 = (TProfile*)p3->Clone("fMomResH3");
}

void AliAnalysisTaskJetHBOM:: SetEfficiencyHybrid(TH1 *h1, TH1 *h2, TH1 *h3) {
  //
  // set tracking efficiency histos
  //

  fhEffH1 = (TH1*)h1->Clone("fhEffH1");
  fhEffH2 = (TH1*)h2->Clone("fhEffH2");
  fhEffH3 = (TH1*)h3->Clone("fhEffH3");
}

Double_t AliAnalysisTaskJetHBOM::GetMomentumSmearing(Int_t cat, Double_t pt) {

  //
  // Get smearing on generated momentum
  //

  //printf("GetMomentumSmearing for cat %d and pt = %f \n",cat,pt);

  TProfile *fMomRes = 0x0;
  if(cat==1) fMomRes = (TProfile*)fMomResH1->Clone("fMomRes");
  if(cat==2) fMomRes = (TProfile*)fMomResH2->Clone("fMomRes");
  if(cat==3) fMomRes = (TProfile*)fMomResH3->Clone("fMomRes");

  if(!fMomRes) {
    return 0.;
  }


  Double_t smear = 0.;

  if(pt>20.) {
    if(cat==1 && fMomResH1Fit) smear = fMomResH1Fit->Eval(pt);
    if(cat==2 && fMomResH2Fit) smear = fMomResH2Fit->Eval(pt);
    if(cat==3 && fMomResH3Fit) smear = fMomResH3Fit->Eval(pt);
  }
  else {

    Int_t bin = fMomRes->FindBin(pt);

    smear = fRandom->Gaus(fMomRes->GetBinContent(bin),fMomRes->GetBinError(bin));

  }

  if(fMomRes) delete fMomRes;
  
  return smear;
}

void AliAnalysisTaskJetHBOM::FitMomentumResolution() {
  //
  // Fit linear function on momentum resolution at high pT
  //

  if(!fMomResH1Fit && fMomResH1) {
    fMomResH1Fit = new TF1("fMomResH1Fit","[0]+[1]*x",0.,200.);
    fMomResH1->Fit(fMomResH1Fit,"LL V0","",5.,30.);
    fMomResH1Fit ->SetRange(5.,100.);
  }

  if(!fMomResH2Fit && fMomResH2) {
    fMomResH2Fit = new TF1("fMomResH2Fit","[0]+[1]*x",0.,200.);
    fMomResH2->Fit(fMomResH2Fit,"LL V0","",5.,30.);
    fMomResH2Fit ->SetRange(5.,100.);
  }

  if(!fMomResH3Fit && fMomResH3) {
    fMomResH3Fit = new TF1("fMomResH3Fit","[0]+[1]*x",0.,200.);
    fMomResH3->Fit(fMomResH3Fit,"LL V0","",5.,30.);
    fMomResH3Fit ->SetRange(5.,100.);
  }

}

