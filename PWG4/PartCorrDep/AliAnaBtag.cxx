/**************************************************************************
 * Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
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
// Class for the electron identification and B-tagging.
// Clusters from EMCAL matched to tracks
// and kept in the AOD. Few histograms produced.
//
// -- Author: T.R.P.Aronsson (Yale) J.L. Klay (Cal Poly), M. Heinz (Yale)
//////////////////////////////////////////////////////////////////////////////
  
// --- ROOT system --- 
#include <TH2F.h>
#include <TH3F.h>
#include <TParticle.h>
#include <TClonesArray.h>
//#include <TObjString.h>
//#include <Riostream.h>

// --- Analysis system --- 
#include "AliAnaBtag.h" 
#include "AliCaloTrackReader.h"
#include "AliMCAnalysisUtils.h"
#include "AliAODCaloCluster.h"
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
//#include "iostream.h"
ClassImp(AliAnaBtag)
  
//____________________________________________________________________________
AliAnaBtag::AliAnaBtag() 
: AliAnaPartCorrBaseClass(),fCalorimeter(""),
  fpOverEmin(0.),fpOverEmax(0.),fResidualCut(0.),fMinClusEne(0.),
  fDrCut(0.),fPairDcaCut(0.),fDecayLenCut(0.),fImpactCut(0.),
  fAssocPtCut(0.),fMassCut(0.),fSdcaCut(0.),fITSCut(0),
  fNTagTrkCut(0),fIPSigCut(0.),fJetEtaCut(0.3),fJetPhiMin(1.8),fJetPhiMax(2.9),
  fhEmcalElectrons(0),fhTRDElectrons(0),fhTPCElectrons(0),fhDVM1(0),fhDVM2(0),fhJets(0),fhJetsAllEtaPhi(0),fhJetsLeadingBElectronEtaPhi(0),fhDVM1EtaPhi(0),fhBJetElectronDetaDphi(0),fhClusterEnergy(0),fhTestalle(0),fhResidual(0),fhElectrons(0),fhTracks(0)
{
  //default ctor
  
  //Initialize parameters
  InitParameters();

}

//____________________________________________________________________________
AliAnaBtag::~AliAnaBtag() 
{
  //dtor

}


//________________________________________________________________________
TList *  AliAnaBtag::GetCreateOutputObjects()
{  
  // Create histograms to be saved in output file and 
  // store them in outputContainer
  TList * outputContainer = new TList() ; 
  outputContainer->SetName("ElectronHistos") ; 


	fhEmcalElectrons = new TH1F("fhEmcalElectrons","",400,0,400);
    outputContainer->Add(fhEmcalElectrons);
    
    fhTRDElectrons = new TH1F("fhTRDElectrons","",400,0,400);
    outputContainer->Add(fhTRDElectrons);

    fhTPCElectrons = new TH1F("fhTPCElectrons","",400,0,400);
    outputContainer->Add(fhTPCElectrons);

    fhDVM1 = new TH1F("fhDVM1","",400,0,400);
    outputContainer->Add(fhDVM1);

    fhDVM2 = new TH1F("fhDVM2","",400,0,400);
    outputContainer->Add(fhDVM2);


    fhJets = new TH2F("fhJets","",400,0,400,20,0,20);
    outputContainer->Add(fhJets);

    fhJetsAllEtaPhi = new TH2F("fhJetsAllEtaPhi","",100,-2,2,100,-2,8);
    outputContainer->Add(fhJetsAllEtaPhi);

    fhJetsLeadingBElectronEtaPhi = new TH2F("fhJetsLeadingBElectronEtaPhi","",100,-5,5,200,-10,10);
    outputContainer->Add(fhJetsLeadingBElectronEtaPhi);

    fhDVM1EtaPhi = new TH2F("fhDVM1EtaPhi","",100,-2,2,100,-2,8);
    outputContainer->Add(fhDVM1EtaPhi);

    fhBJetElectronDetaDphi = new TH2F("fhBJetElectronDetaDphi","",100,-5,5,200,-10,10);
    outputContainer->Add(fhBJetElectronDetaDphi);

    fhClusterEnergy = new TH1F("fhClusterEnergy","",100,0,10);
    outputContainer->Add(fhClusterEnergy);

    fhTestalle = new TH1F("fhTestalle","",400,0,400);
    outputContainer->Add(fhTestalle);

    fhResidual = new TH1F("fhResidual","",500,0,5);
    outputContainer->Add(fhResidual);

    fhElectrons = new TH2F("fhElectrons","",200,0,100,20,0,20);
    outputContainer->Add(fhElectrons);

    fhTracks = new TH2F("fhTracks","",200,0,100,20,0,20);
    outputContainer->Add(fhTracks);


  
  return outputContainer ;
  
}

//____________________________________________________________________________
void AliAnaBtag::Init()
{

  //do some initialization
  if(fCalorimeter == "PHOS") {
    printf("AliAnaBtag::Init() - !!STOP: You want to use PHOS in analysis but this is not (yet) supported!!\n!!Check the configuration file!!\n");
    fCalorimeter = "EMCAL";
  }
  if(fCalorimeter == "EMCAL" && !GetReader()->IsEMCALSwitchedOn()){
    printf("AliAnaBtag::Init() - !!STOP: You want to use EMCAL in analysis but it is not read!!\n!!Check the configuration file!!\n");
    abort();
  }

}


//____________________________________________________________________________
void AliAnaBtag::InitParameters()
{
  
  //Initialize the parameters of the analysis.
  SetOutputAODClassName("AliAODPWG4Particle");
  SetOutputAODName("PWG4Particle");

  AddToHistogramsName("Btag_");

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
void  AliAnaBtag::MakeAnalysisFillAOD() 
{
  //This reads in tracks, extrapolates to EMCAL, does p/E selectrons, identifies electron candidates
  //After candidates are obtained, btagging and saving into AOD.

  Double_t bfield = 0.;
  if(GetReader()->GetDataType() != AliCaloTrackReader::kMC) bfield = GetReader()->GetBField();

  //Select the calorimeter of the electron
  if(fCalorimeter != "EMCAL") {
    printf("This class not yet implemented for PHOS\n");
    abort();
  }
  
  TObjArray *cl = GetAODEMCAL();
  
  if(!GetAODCTS() || GetAODCTS()->GetEntriesFast() == 0) return ;
  Int_t ntracks = GetAODCTS()->GetEntriesFast();
  if(GetDebug() > 0)
    printf("AliAnaBtag::MakeAnalysisFillAOD() - In CTS aod entries %d\n", ntracks);

  Int_t iCluster = -999;
  Int_t ntot = cl->GetEntriesFast();

  //CLUSTER STUFF
  if(1){
    for(Int_t iclus = 0; iclus < ntot; iclus++) {
      AliAODCaloCluster * clus = (AliAODCaloCluster*) (cl->At(iclus));
      if(!clus) continue;
      fhClusterEnergy->Fill(clus->E());
    }
  }


  for (Int_t itrk =  0; itrk <  ntracks; itrk++) {////////////// track loop
    iCluster = -999; //start with no match
    AliAODTrack * track = (AliAODTrack*) (GetAODCTS()->At(itrk)) ;

    Double_t imp[2] = {-999.,-999.}; Double_t cov[3] = {-999.,-999.,-999.};
    Bool_t dcaOkay = GetDCA(track,imp,cov);  //homegrown dca calculation until AOD is fixed
    if(!dcaOkay) printf("AliAnaBtag::Problem computing DCA to primary vertex for track %d.  Skipping it...\n",itrk);

    fhTracks->Fill(track->Pt(),1);

    AliAODPid* pid = (AliAODPid*) track->GetDetPid();
    if(pid == 0) {
      if(GetDebug() > 0) printf("AliAnaBtag::MakeAnalysisFillAOD() - No PID object - skipping track %d",itrk);
      continue;
    } 
    else {
      Double_t emcpos[3];
      pid->GetEMCALPosition(emcpos);
      Double_t emcmom[3];
      pid->GetEMCALMomentum(emcmom);
      
      TVector3 pos(emcpos[0],emcpos[1],emcpos[2]);
      TVector3 mom(emcmom[0],emcmom[1],emcmom[2]);
      Double_t tphi = pos.Phi();
      Double_t teta = pos.Eta();
      Double_t tmom = mom.Mag();

      Bool_t in = kFALSE;
      if(mom.Phi()*180./TMath::Pi() > 80. && mom.Phi()*180./TMath::Pi() < 190. &&
	 mom.Eta() > -0.7 && mom.Eta() < 0.7) in = kTRUE; //Checks momentum?
      //Also check the track
      if(track->Phi()*180./TMath::Pi() > 80. && track->Phi()*180./TMath::Pi() < 190. &&
	 track->Eta() > -0.7 && track->Eta() < 0.7) in = kTRUE;

      if(GetDebug() > 1) printf("AliAnaBtag::MakeAnalysisFillAOD() - Track(Extrap) pt %2.2f(%2.2f), phi %2.2f(%2.2f), eta %2.2f(%2.2f) in fiducial cut %d\n",track->Pt(), mom.Pt(), track->Phi(), mom.Phi(), track->Eta(),mom.Eta(), in);

      if(mom.Pt() > GetMinPt() && in) {
	fhTracks->Fill(track->Pt(),3);
	Double_t dEdx = pid->GetTPCsignal();
		
	//Int_t ntot = cl->GetEntriesFast();
	Double_t res = 999.;
	Double_t pOverE = -999.;
	
	Int_t pidProb = track->GetMostProbablePID();
	Bool_t tpcEle = kFALSE; if(dEdx > 70.) tpcEle = kTRUE;
	Bool_t trkEle = kFALSE; if(pidProb == AliAODTrack::kElectron) trkEle = kTRUE;
	Bool_t trkChgHad = kFALSE; if(pidProb == AliAODTrack::kPion || pidProb == AliAODTrack::kKaon || pidProb == AliAODTrack::kProton) trkChgHad = kTRUE;


	//Track Matching!
	Bool_t emcEle = kFALSE;      
	double minRes=100.;
	Double_t minR  = 99;
        Double_t minPe =-1;
        Double_t minEp =-1;
        Double_t minMult = -1;
        Double_t minPt   = -1;

	for(Int_t iclus = 0; iclus < ntot; iclus++) {
	  AliAODCaloCluster * clus = (AliAODCaloCluster*) (cl->At(iclus));
	  if(!clus) continue;

	  //As of 11-Oct-2009
	  //only select "good" clusters	  
// 	  if (clus->GetNCells()       < 2    ) continue;
//           if (clus->GetNCells()       > 30   ) continue;
//           if (clus->E()               < fMinClusEne ) continue;
//           if (clus->GetDispersion()   > 1    ) continue;
//           if (clus->GetM20()          > 0.4  ) continue;
//           if (clus->GetM02()          > 0.4  ) continue;
//           if (clus->GetM20()          < 0.03 ) continue;
//           if (clus->GetM02()          < 0.03 ) continue;
// new optimized from ben. 2010May
	  if (clus->GetNCells()       < 2    ) continue;
          if (clus->GetNCells()       > 35   ) continue;
          if (clus->E()               < 0 ) continue;
          if (clus->GetDispersion()   > 1.08    ) continue;
          if (clus->GetM20()          > 0.42  ) continue;
          if (clus->GetM02()          > 0.4  ) continue;
          if (clus->GetM20()          < 0 ) continue;
          if (clus->GetM02()          < 0.06 ) continue;
	  
	  
	  Double_t x[3];
	  clus->GetPosition(x);
	  TVector3 cluspos(x[0],x[1],x[2]);
	  Double_t deta = teta - cluspos.Eta();
	  Double_t dphi = tphi - cluspos.Phi();
	  if(dphi > TMath::Pi()) dphi -= 2*TMath::Pi();
	  if(dphi < -TMath::Pi()) dphi += 2*TMath::Pi();

	  res = sqrt(dphi*dphi + deta*deta);
	  
	  if(res<minRes){
	    minRes=res;
	    
	  }
	  
	  if(res < 0.0275) { //	  if(res < fResidualCut) { //Optimized from Ben
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


	  } else {
	      //unmatched
	  }//res cut
	}//calo cluster loop

	fhResidual->Fill(minRes);

	if(minPe > 0.9 && minPe < 1.08) emcEle = kTRUE;//	if(minPe > fpOverEmin && minPe < fpOverEmax) emcEle = kTRUE;
	

	
	if(emcEle)
	  fhEmcalElectrons->Fill(track->Pt());
	if(trkEle)
	  fhTRDElectrons->Fill(track->Pt());
	if(tpcEle)
	  fhTPCElectrons->Fill(track->Pt());
	
	
	//Take all emcal electrons, but the others only if pT < 10 GeV
	if(emcEle && (tpcEle || trkEle) ) {
	  fhTestalle->Fill(track->Pt());
	  //B-tagging
	  if(GetDebug() > 1) printf("Found Electron - do b-tagging\n");
	  Int_t dvmbtag = GetDVMBtag(track);
	  
	  
	  if(dvmbtag>0){
	    fhDVM1->Fill(track->Pt());
            fhDVM1EtaPhi->Fill(track->Eta(),track->Phi());	  
	    
          }
	  if(dvmbtag>1)
	    fhDVM2->Fill(track->Pt());
	  

	  //Create particle to save in AOD
	  Double_t eMass = 0.511/1000; //mass in GeV
	  Double_t eleE = sqrt(track->P()*track->P() + eMass*eMass);
	  AliAODPWG4Particle tr = AliAODPWG4Particle(track->Px(),track->Py(),track->Pz(),eleE);
	  tr.SetLabel(track->GetLabel());
	  tr.SetCaloLabel(iCluster,-1); //sets the indices of the original caloclusters
	  tr.SetTrackLabel(track->GetID(),-1); //sets the indices of the original tracks tr.SetTrackLabel(track->GetID(),-1) instead of itrk;

	  if(emcEle) {//PID determined by EMCAL
	    tr.SetDetector(fCalorimeter);
	  } else {
	    tr.SetDetector("CTS"); //PID determined by CTS
	  }

	  if(GetReader()->GetAODCTSNormalInputEntries() <= itrk) tr.SetInputFileIndex(1);
	  //Make this preserve sign of particle
	  if(track->Charge() < 0) tr.SetPdg(11); //electron is 11
	  else  tr.SetPdg(-11); //positron is -11

	  tr.SetBtag(dvmbtag);
	  

	  //Check origin of the candiates
	  if(IsDataMC()){
	    tr.SetTag(GetMCAnalysisUtils()->CheckOrigin(tr.GetLabel(),GetReader(),tr.GetInputFileIndex()));
	    if(GetDebug() > 0) printf("AliAnaBtag::MakeAnalysisFillAOD() - Origin of candidate %d\n",tr.GetTag());
	  }//Work with stack also   
	  
	  AddAODParticle(tr);
	  
	  if(GetDebug() > 1) printf("AliAnaBtag::MakeAnalysisFillAOD() - Electron selection cuts passed: pT %3.2f, pdg %d\n",tr.Pt(),tr.GetPdg());	
	}//electron
      }//pt, fiducial selection
    }//pid check
  }//track loop                         
  



  if(GetDebug() > 1) printf("AliAnaBtag::MakeAnalysisFillAOD()  End fill AODs \n");  
  
}

//__________________________________________________________________
void  AliAnaBtag::MakeAnalysisFillHistograms() 
{
  //Do analysis and fill histograms

  AliStack * stack = 0x0;
//   TParticle * primary = 0x0;




  if(IsDataMC()) {
    if(GetReader()->ReadStack()){
      stack =  GetMCStack() ;      
      if(!stack)
	printf("AliAnaBtag::MakeAnalysisFillHistograms() *** no stack ***: \n");
      
    }
  }// is data and MC

  ////////////////////////////////////
  //Loop over jets and check for b-tag
  ////////////////////////////////////
  double maxjetEta=-4.;
  double maxjetPhi=-4.;

  Int_t njets = (GetReader()->GetOutputEvent())->GetNJets();
  if(njets > 0) {
    if(GetDebug() > 0) printf("AliAnaBtag::MakeAnalysisFillHistograms() - Jet AOD branch has %d jets.  Performing b-jet tag analysis\n",njets);

    for(Int_t ijet = 0; ijet < njets ; ijet++) {
      AliAODJet * jet = (AliAODJet*)(GetReader()->GetOutputEvent())->GetJet(ijet) ;

      if(ijet==0){
        maxjetEta=jet->Eta();
        maxjetPhi=jet->Phi();
      }

      fhJets->Fill(jet->Pt(),1);
      fhJetsAllEtaPhi->Fill(jet->Eta(),jet->Phi());

      if(jet->Pt() < 0.) continue; //This has to be adjusted depending on pp or AA!
      fhJets->Fill(jet->Pt(),3); //All jets after pt cut

      //Geometric EMCAL cut
      if(TMath::Abs(jet->Eta()) > fJetEtaCut) continue;
      if(jet->Phi() < fJetPhiMin || jet->Phi() > fJetPhiMax) continue;
      fhJets->Fill(jet->Pt(),4); //All jets after geometric cut

      Bool_t leadJet  = kFALSE;
      if (ijet==0){ 
	fhJets->Fill(jet->Pt(),5); //Leading jets
	leadJet= kTRUE;
      }


      Bool_t dvmJet = kFALSE;  
      TRefArray* rt = jet->GetRefTracks();
      Int_t ntrk = rt->GetEntries();

      for(Int_t itrk = 0; itrk < ntrk; itrk++) {
      	AliAODTrack* jetTrack = (AliAODTrack*)jet->GetTrack(itrk);
	Bool_t isDVM = CheckIfBjet(jetTrack);
	if(isDVM) dvmJet = kTRUE;
      }

      if(dvmJet)
	fhJets->Fill(jet->Pt(),6);



      if(IsDataMC()) {
	//determine tagging efficiency & mis-tagging rate
	//using b-quarks from stack
	Bool_t isTrueBjet = IsMcBJet(jet->Eta(), jet->Phi());
	Bool_t isTrueDjet = IsMcDJet(jet->Eta(), jet->Phi());
	if (isTrueBjet && GetDebug() > 0) printf("== True Bjet==\n");
	if (isTrueDjet && GetDebug() > 0) printf("== True Charm-jet==\n");
	if (dvmJet && GetDebug() > 0)     printf("== found DVM jet==\n");
      }

    } //jet loop
  } //jets exist
  



  //Electron loop, read back electrons, fill histos
  Int_t naod = GetOutputAODBranch()->GetEntriesFast();
  if(GetDebug() > 0) printf("AliAnaBtag::MakeAnalysisFillHistograms() - aod branch entries %d\n", naod);
  
  for(Int_t iaod = 0; iaod < naod ; iaod++){
    AliAODPWG4Particle* ele =  (AliAODPWG4Particle*) (GetOutputAODBranch()->At(iaod));
    Int_t pdg = ele->GetPdg();


    if(TMath::Abs(pdg) != AliCaloPID::kElectron) continue; //not necessary..

    //MC tag of this electron
    //    Int_t mctag = ele->GetTag();


    fhElectrons->Fill(ele->Pt(),1); //All electrons
    Bool_t photonic = kFALSE;
    photonic = PhotonicV0(ele->GetTrackLabel(0)); //check against V0s
    if(!photonic) fhElectrons->Fill(ele->Pt(),3); //nonphotonic electrons
    if(photonic) fhElectrons->Fill(ele->Pt(),4);  //photonic electrons

    //Fill electron histograms 
    Float_t phiele = ele->Phi();
    Float_t etaele = ele->Eta();


    if(ele->GetBtag()>0){ // removed bit tag shit
      fhElectrons->Fill(ele->Pt(),5);
      if(!photonic) fhElectrons->Fill(ele->Pt(),6);
      if(photonic) fhElectrons->Fill(ele->Pt(),7);
      fhJetsLeadingBElectronEtaPhi->Fill(maxjetEta,maxjetPhi); 
      double deta=etaele-maxjetEta;
      double dphi=phiele-maxjetPhi;
      
      //double r = sqrt((dphi*dphi)+(deta*deta));
      fhBJetElectronDetaDphi->Fill(deta,dphi);
      
    }
    
  }//electron aod loop

}

//__________________________________________________________________
Int_t AliAnaBtag::GetDVMBtag(AliAODTrack * tr )
{
  //This method uses the Displaced Vertex between electron-hadron
  //pairs and the primary vertex to determine whether an electron is
  //likely from a B hadron.

  Int_t ncls1 = 0;
  for(Int_t l = 0; l < 6; l++) if(TESTBIT(tr->GetITSClusterMap(),l)) ncls1++;


  if (ncls1 < fITSCut) return 0;

  Double_t imp[2] = {-999.,-999.}; Double_t cov[3] = {-999.,-999.,-999.};
  Bool_t dcaOkay = GetDCA(tr,imp,cov);  //homegrown dca calculation until AOD is fixed                  
  if(!dcaOkay) {
    printf("AliAnaBtag::Problem computing DCA to primary vertex for track %d",tr->GetID());
    return 0;
  }

  if (TMath::Abs(imp[0])   > fImpactCut ) return 0;
  if (TMath::Abs(imp[1])   > fImpactCut ) return 0;

  Int_t nvtx1 = 0;
  Int_t nvtx2 = 0;
  Int_t nvtx3 = 0;

  for (Int_t k2 =0; k2 < GetAODCTS()->GetEntriesFast() ; k2++) {
    //loop over assoc
    AliAODTrack* track2 = (AliAODTrack*)GetAODCTS()->At(k2);
    Int_t id1 = tr->GetID();
    Int_t id2 = track2->GetID();
    if(id1 == id2) continue;

    Int_t ncls2 = 0;
    for(Int_t l = 0; l < 6; l++) if(TESTBIT(track2->GetITSClusterMap(),l)) ncls2++;
    if (ncls2 < fITSCut) continue;

    if(tr->Pt()<6.&&track2->Pt() < 0.4) continue;
    if(tr->Pt()>6.&&tr->Pt()<10.&&track2->Pt() < 0.2) continue;
    if(tr->Pt()>10.&&track2->Pt() < 0.6) continue;


    Double_t dphi = tr->Phi() - track2->Phi();
    if(dphi > TMath::Pi()) dphi -= 2*TMath::Pi();
    if(dphi < -TMath::Pi()) dphi += 2*TMath::Pi();
    Double_t deta = tr->Eta() - track2->Eta();
    Double_t dr = sqrt(deta*deta + dphi*dphi);

    if(dr > fDrCut) continue;
    
    if(tr->Pt()<6.){
      Double_t sDca = ComputeSignDca(tr, track2, 1.4,0.025);
	if(sDca > 0.06) nvtx2++;
    } 
    if(tr->Pt()>6.&&tr->Pt()<10.){
      Double_t sDca = ComputeSignDca(tr, track2, 1.7,0.012);
        if(sDca > 0.03) nvtx2++;
    } 
    if(tr->Pt()>10.){
      Double_t sDca = ComputeSignDca(tr, track2, 1.5,0.14);
        if(sDca > 0.04) nvtx2++;
    } 


  } //loop over hadrons

  if(GetDebug() > 0) {
    if (nvtx1>0) printf("result1 of btagging: %d \n",nvtx1);
    if (nvtx2>0) printf("result2 of btagging: %d \n",nvtx2);
    if (nvtx3>0) printf("result3 of btagging: %d \n",nvtx3);
  }




  return nvtx2;

}

//__________________________________________________________________
Double_t AliAnaBtag::ComputeSignDca(AliAODTrack *tr, AliAODTrack *tr2 , float masscut, double pdcacut)
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
    GetReader()->GetVertex(vertex); //If only one file, get the vertex from there
    //FIXME:  Add a check for whether file 2 is PYTHIA or HIJING
    //If PYTHIA, then set the vertex from file 2, if not, use the
    //vertex from file 1
    if(GetReader()->GetSecondInputAODTree()) GetReader()->GetSecondInputAODVertex(vertex);
  }
  
  TVector3 primV(vertex[0],vertex[1],vertex[2]) ;

  if(GetDebug() > 5) printf(">>ComputeSdca:: primary vertex = %2.2f,%2.2f,%2.2f \n",vertex[0],vertex[1],vertex[2]) ;

  AliExternalTrackParam *param1 = new AliExternalTrackParam(tr);
  AliExternalTrackParam *param2 = new AliExternalTrackParam(tr2);

  Double_t bfield[3];
  param1->GetBxByBz(bfield);
  bfield[0]=0;
  bfield[1]=0;
  bfield[2]=5.;
  Double_t bz = 5.; // = param1->GetBz();
  //cout<<"In ComputeSignDCA, bfield[3] and bz: "<<bfield[0]<<" "<<bfield[1]<<" "<<bfield[2]<<" "<<bz<<endl;
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



  if (emomAtB.Mag()>0 && pairdca < pdcacut && decaylength < fDecayLenCut ) {
    TVector3 sumMom = emomAtB+hmomAtB;
    Double_t ener1 = sqrt(pow(emomAtB.Mag(),2) + massE*massE);
    Double_t ener2 = sqrt(pow(hmomAtB.Mag(),2) + massK*massK);
    Double_t ener3 = sqrt(pow(hmomAtB.Mag(),2) + massE*massE);
    Double_t mass = sqrt(pow((ener1+ener2),2) - pow(sumMom.Mag(),2));
    Double_t massPhot = sqrt(pow((ener1+ener3),2) - pow(sumMom.Mag(),2));
    Double_t sDca = decayvector.Dot(emomAtB)/emomAtB.Mag();

//     if (masscut<1.1) fhDVMBtagQA2->Fill(sDca, mass);

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
Bool_t AliAnaBtag::PhotonicV0(Int_t id) 
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
Bool_t AliAnaBtag::GetDCA(const AliAODTrack* track,Double_t impPar[2], Double_t cov[3]) 
{
  //Use the Event vertex and AOD track information to get
  //a real impact parameter for the track


  Double_t maxD = 100000.; //max transverse IP
  if(GetReader()->GetDataType() != AliCaloTrackReader::kMC) {
    AliVEvent* ve = (AliVEvent*)GetReader()->GetInputEvent();
    AliVVertex *vv = (AliVVertex*)ve->GetPrimaryVertex();
    AliESDtrack esdTrack(track);
    Double_t bfield[3];
    esdTrack.GetBxByBz(bfield);
    bfield[0]=0;
    bfield[1]=0;
    bfield[2]=5.;
    
    Bool_t gotit = esdTrack.PropagateToDCABxByBz(vv,bfield,maxD,impPar,cov);
    return gotit;
  }

  return kFALSE;

}
//__________________________________________________________________
Bool_t AliAnaBtag::CheckIfBjet(const AliAODTrack* track)
{
  Bool_t bjet = kFALSE;
  Int_t trackId = track->GetID(); //get the index in the reader
  Int_t naod = GetOutputAODBranch()->GetEntriesFast();
  if(GetDebug() > 3) printf("AliAnaBtag::CheckIfBjet() - aod branch entries %d\n", naod);

  for(Int_t iaod = 0; iaod < naod ; iaod++){
    AliAODPWG4Particle* ele =  (AliAODPWG4Particle*) (GetOutputAODBranch()->At(iaod));
    Int_t electronlabel = ele->GetTrackLabel(0);
    if(electronlabel != trackId) continue;  //skip to the next one if they don't match
    if(ele->GetBtag()>0)
      bjet = kTRUE;
  } 

  return bjet;
} 



//__________________________________________________________________
AliAODMCParticle* AliAnaBtag::GetMCParticle(Int_t ipart) 
{
  //Get the MC particle at position ipart

  AliAODMCParticle* aodprimary = 0x0;
  TClonesArray * mcparticles0 = 0x0;
  TClonesArray * mcparticles1 = 0x0;

  if(GetReader()->ReadAODMCParticles()){
    //Get the list of MC particles                                                                                                                           
    mcparticles0 = GetReader()->GetAODMCParticles(0);
    if(!mcparticles0 && GetDebug() > 0) {
      printf("AliAnaBtag::MakeAnalysisFillHistograms() -  Standard MCParticles not available!\n");
    }
    if(GetReader()->GetSecondInputAODTree()){
      mcparticles1 = GetReader()->GetAODMCParticles(1);
      if(!mcparticles1 && GetDebug() > 0) {
	printf("AliAnaBtag::MakeAnalysisFillHistograms() -  Second input MCParticles not available!\n");
      }
    }

    Int_t npart0 = mcparticles0->GetEntriesFast();
    Int_t npart1 = 0;
    if(mcparticles1) npart1 = mcparticles1->GetEntriesFast();
    if(ipart < npart0) aodprimary = (AliAODMCParticle*)mcparticles0->At(ipart);
    else aodprimary = (AliAODMCParticle*)mcparticles1->At(ipart-npart0);
    if(!aodprimary) {
      printf("AliAnaBtag::GetMCParticle() *** no primary ***:  label %d \n", ipart);
      return 0x0;
    }

  } else {
    printf("AliAnaBtag::GetMCParticle() - Asked for AliAODMCParticle but we have a stack reader.\n");
  }
  return aodprimary;

}

//__________________________________________________________________
Bool_t  AliAnaBtag::IsMcBJet(Double_t jeta, Double_t jphi)
{
  //Check the jet eta,phi against that of the b-quark
  //to decide whether it is an MC B-jet
  Bool_t bjet=kFALSE;

  AliStack* stack = 0x0;
  
  for(Int_t ipart = 0; ipart < 100; ipart++) {

    Double_t pphi = -999.;
    Double_t peta = -999.;
    Int_t pdg = 0;
    if(GetReader()->ReadStack()) {
      stack = GetMCStack();
      if(!stack) {
	printf("AliAnaBtag::IsMCBJet() *** no stack ***: \n");
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
Bool_t  AliAnaBtag::IsMcDJet(Double_t jeta, Double_t jphi)
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
	printf("AliAnaBtag::IsMCDJet() *** no stack ***: \n");
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
void AliAnaBtag::Print(const Option_t * opt) const
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


//__________________________________________________________________
void  AliAnaBtag::Terminate(TList* outputList)
{
 
  //Do some plots to end
  //Recover histograms from output histograms list, needed for
  //distributed analysis.                
  //ReadHistograms(outputList);

  printf(" AliAnaBtag::Terminate()  *** %s Report: %d outputs\n", GetName(), outputList->GetEntries()) ;

}

