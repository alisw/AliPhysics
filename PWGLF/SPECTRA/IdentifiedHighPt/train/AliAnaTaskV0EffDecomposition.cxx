/*CHANGES:

 * 18 Jul 2016: changed output filename in AddTask.C
 * 05 Aug 2016: implement standard cuts, added also in .h
 * 17 Aug 2016: added #includes that were missing
 * 17 Aug 2016: removed fminpt
 * 19 Aug 2016: it compiles locally but gives no output
 * 20 Aug 2016: remove mass cut and remove vtxSPD
 * 07 Sep 2016: histogram cuts
 * 10 Sep 2016: 
	       a) mass cut
	       b) aodV0->DcaPosToPrimVertex(); and aodV0->DcaNegToPrimVertex();
	       c) and dcaV0
 * 11 Sep 2016: testing filter bits 1, 128, 768, 1024, or hard coded -> run without filter bit on tracks for now 
 * 14 Sep 2016: 
               a) adding a weight to the ghost V0 and tracks -> solving the track pair efficiency issue, but need to add TPCrefit later to get rid of ghost tracks (for now just subtract them)
	       b) clean-up and change hard coded values (generalize)
 * 23 Sep 2016: 
               a) finer pt binning
               b) pt window in cut distributions
 * 31 jul 2017: test push this comment	  

//to run locally:
aliroot
.L runAAFDev.C
runAAF(2, "local aod PbPb MC", "test_aod_mc.txt", 4)

*/


#include "AliAnaTaskV0EffDecomposition.h"
#include <TList.h>
#include <TH1.h>
#include <AliAnalysisManager.h>
#include <AliMCEventHandler.h>
#include <AliMCEvent.h>

//#include <bitset>
//#include <TTree.h>
//#include <TMath.h> 
//#include <TParticle.h>
//#include <TFile.h>
//#include <AliAnalysisFilter.h>
//#include <AliLog.h>
//#include <AliExternalTrackParam.h> 
//#include <AliAODVZERO.h>
//#include <AliStack.h>
//#include <AliHeader.h>
//#include <AliGenPythiaEventHeader.h>
//#include <AliGenDPMjetEventHeader.h>
//#include "AliCentrality.h" 
// #include <AliKFVertex.h>
// #include <AliAODVertex.h>
// #include <AliAODTrack.h> 
// #include <AliAODPid.h> 
// #include <AliAODMCHeader.h> 
 
#include <iostream>


// STL includes
#include <iostream>
using namespace std;


ClassImp(AliAnaTaskV0EffDecomposition)

//_____________________________________________________________________________
AliAnaTaskV0EffDecomposition::AliAnaTaskV0EffDecomposition():
fLowPtFraction(0.01),
  fMassCut(0.1),
  fPdca(0.1),
  fNdca(0.1),
  fMinPtV0(0.1), 
  fDecayRmax(100.),
  fDecayRmin(5.0),  
  fDcaDaugh(1.0),
  fDcaV0(-100),
  fV0pndca(0.1),
  fCospt(0.998),
  AliAnalysisTaskSE(),
  fAOD(0x0),
  fMCArray(0x0),
  fTrackFilterBit(1), // default is TPC filter bit
  fTrigBit(AliVEvent::kMB), // default is MB
  fVtxCut(10),        // default is 10 cm
  fEtaCut(0.8),       // default is 0.8
  fMinCent(0),        // default is 0%
  fMaxCent(10),       // default is 10%
  fPdgV0(3122),       // default is Lambda
  fPdgPos(2212),      // default is p
  fPdgNeg(-211),      // default is pi-
  fListOfObjects(0x0),
  fCt(3.0),
  fNcl(60.),
  fChi2perNDF(4.0),
  hV0Gen(0x0),
  hV0Rec(0x0),
  hDaughterRec(0x0),
  hV0Ghost(0x0),
  hTrackGhost(0x0),
  hV0ButNoTracks(0x0),
  hCentr(0x0),            
  hCt(0x0),            
  hDCAdaugh(0x0),
  hDecayR(0x0),
  hCosPA(0x0),
  hInvMass(0x0),
  hNcl(0x0),
  hChi2perNDF(0x0),
  hPdca(0x0),
  hNdca(0x0),
  hDcaV0(0x0),
  hCtHigh(0x0),            
  hDCAdaughHigh(0x0),
  hDecayRHigh(0x0),
  hCosPAHigh(0x0),
  hInvMassHigh(0x0),
  hNclHigh(0x0),
  hChi2perNDFHigh(0x0),
  hPdcaHigh(0x0),
  hNdcaHigh(0x0),
  hDcaV0High(0x0)
{
  // Default constructor (should not be used)
}

//______________________________________________________________________________
AliAnaTaskV0EffDecomposition::AliAnaTaskV0EffDecomposition(const char *name):
  fLowPtFraction(0.01),
  fMassCut(0.1),
  fPdca(0.1),
  fNdca(0.1),
  fMinPtV0(0.1),
  fDecayRmax(100.),
  fDecayRmin(5.0),  
  fDcaDaugh(1.0),
  fDcaV0(-100),
  fV0pndca(0.1),
  fCospt(0.998),
  AliAnalysisTaskSE(name),
  fAOD(0x0),
  fMCArray(0x0),
  fTrackFilterBit(1), // default is TPC filter bit
  fTrigBit(AliVEvent::kMB), // default is MB
  fVtxCut(10),        // default is 10 cm
  fEtaCut(0.8),       // default is 0.8
  fMinCent(0),        // default is 0%
  fMaxCent(10),       // default is 10%
  fPdgV0(3122),       // default is Lambda
  fPdgPos(2212),      // default is p
  fPdgNeg(-211),      // default is pi-
  fListOfObjects(0x0),
  fCt(3.0),
  fNcl(60.),
  fChi2perNDF(4.0),
  hV0Gen(0x0),
  hV0Rec(0x0),
  hDaughterRec(0x0),
  hV0Ghost(0x0),
  hTrackGhost(0x0),
  hV0ButNoTracks(0x0),
  hCt(0x0),            
  hCentr(0x0),            
  hDCAdaugh(0x0),
  hDecayR(0x0),
  hCosPA(0x0),
  hInvMass(0x0),
  hNcl(0x0),
  hChi2perNDF(0x0),
  hPdca(0x0),
  hNdca(0x0),
  hDcaV0(0x0),
 hCtHigh(0x0),            
  hDCAdaughHigh(0x0),
  hDecayRHigh(0x0),
  hCosPAHigh(0x0),
  hInvMassHigh(0x0),
  hNclHigh(0x0),
  hChi2perNDFHigh(0x0),
  hPdcaHigh(0x0),
  hNdcaHigh(0x0),
  hDcaV0High(0x0)
{
  // Output slot #1 writes into a TList
  DefineOutput(1, TList::Class());
}

//_____________________________________________________________________________
AliAnaTaskV0EffDecomposition::~AliAnaTaskV0EffDecomposition()
{
  // Destructor
  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor
}

//______________________________________________________________________________
void AliAnaTaskV0EffDecomposition::UserCreateOutputObjects()
{ 
  // This method is called once per worker node
  // Here we define the output: histograms and debug tree if requested 

  OpenFile(1);
  fListOfObjects = new TList();
  fListOfObjects->SetOwner();
  
  //
  // Histograms
  //  

  const Int_t nPtBins = 37;
  Double_t ptBins[nPtBins+1] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0, 5.5, 6.5, 8.0, 10.0, 14.0};
 
  hV0Gen = new TH1D("hV0Gen", "Number of generated V0s; p_{T} [GeV/c]; Counts", nPtBins, ptBins);
  fListOfObjects->Add(hV0Gen);

  hV0Rec = new TH1D("hV0Rec", "Number of reconstructed V0s; p_{T} [GeV/c]; Counts", nPtBins, ptBins);
  fListOfObjects->Add(hV0Rec);

  hDaughterRec = new TH1D("hDaughterRec", "Number of reconstructed track pairs; V0 p_{T} [GeV/c]; Counts", nPtBins, ptBins);
  fListOfObjects->Add(hDaughterRec);

  hV0Ghost = new TH1D("hV0Ghost", "Number of reconstructed V0 ghosts; V0 p_{T} [GeV/c]; Counts", nPtBins, ptBins);
  fListOfObjects->Add(hV0Ghost);

  hTrackGhost = new TH1D("hTrackGhost", "Number of reconstructed track pair ghosts; V0 p_{T} [GeV/c]; Counts", nPtBins, ptBins);
  fListOfObjects->Add(hTrackGhost);

  hV0ButNoTracks = new TH1D("hV0ButNoTracks", "Strange V0s: V0 rec but no daughter tracks; V0 p_{T} [GeV/c]; Counts", nPtBins, ptBins);
  fListOfObjects->Add(hV0ButNoTracks);
  
  hCentr = new TH1D("hCentr", "Centrality distribution", 110, 0, 11);
  fListOfObjects->Add(hCentr);


  hCt = new TH1D("hCt", "DCA to primary vertex; DCA [cm]; Counts", 300, 0, 50);
  fListOfObjects->Add(hCt);

  hDCAdaugh = new TH1D("hDCAdaugh", "DCA between daughters; DCA_{d-d} [cm]; Counts", 200, 0, 2);
  fListOfObjects->Add(hDCAdaugh);

  hDecayR = new TH1D("hDecayR", "V0 decay radius; Decay radius [cm]; Counts", 240, 0, 120);
  fListOfObjects->Add(hDecayR);

  hCosPA = new TH1D("hCosPA", "Cosine of pointing angle; cos(PA); Counts", 600, 0.998, 1);
  fListOfObjects->Add(hCosPA);

  hInvMass = new TH1D("hInvMass", "Inv Mass; inv mass (GeV/c2); Counts", 600, -0.3, 0.3);
  fListOfObjects->Add(hInvMass);

  hNcl = new TH1D("hNcl", "Number of TPC clusters; ncl; Counts", 170, 0, 170);
  fListOfObjects->Add(hNcl);

  hChi2perNDF = new TH1D("hChi2perNDF", "Chi2 per NDF; Chi2; Counts", 100, 0, 100);
  fListOfObjects->Add(hChi2perNDF);

  hPdca = new TH1D("hPdca", "DCA of positive daughter to primary vertex; DCA_{pd-PV}; Counts", 1000, 0., 100);
  fListOfObjects->Add(hPdca);

  hNdca = new TH1D("hNdca", "DCA of negative daughter to primary vertex; DCA_{nd-PV}; Counts", 1000, 0., 100);
  fListOfObjects->Add(hNdca);

  hDcaV0 = new TH1D("hDcaV0", "DCA of V0 to primary vertex; DCA_{V0-PV}; Counts", 1000, 0., 10);
  fListOfObjects->Add(hDcaV0);


   

  hCtHigh = new TH1D("hCtHigh", "DCA to primary vertex; DCA [cm]; Counts", 300, 0, 50);
  fListOfObjects->Add(hCtHigh);

  hDCAdaughHigh = new TH1D("hDCAdaughHigh", "DCA between daughters; DCA_{d-d} [cm]; Counts", 200, 0, 2);
  fListOfObjects->Add(hDCAdaughHigh);

  hDecayRHigh = new TH1D("hDecayRHigh", "V0 decay radius; Decay radius [cm]; Counts", 240, 0, 120);
  fListOfObjects->Add(hDecayRHigh);

  hCosPAHigh = new TH1D("hCosPAHigh", "Cosine of pointing angle; cos(PA); Counts", 600, 0.998, 1);
  fListOfObjects->Add(hCosPAHigh);

  hInvMassHigh = new TH1D("hInvMassHigh", "Inv Mass; inv mass (GeV/c2); Counts", 600, -0.3, 0.3);
  fListOfObjects->Add(hInvMassHigh);

  hNclHigh = new TH1D("hNclHigh", "Number of TPC clusters; ncl; Counts", 170, 0, 170);
  fListOfObjects->Add(hNclHigh);

  hChi2perNDFHigh = new TH1D("hChi2perNDFHigh", "Chi2 per NDF; Chi2; Counts", 100, 0, 100);
  fListOfObjects->Add(hChi2perNDFHigh);

  hPdcaHigh = new TH1D("hPdcaHigh", "DCA of positive daughter to primary vertex; DCA_{pd-PV}; Counts", 1000, 0., 100);
  fListOfObjects->Add(hPdcaHigh);

  hNdcaHigh = new TH1D("hNdcaHigh", "DCA of negative daughter to primary vertex; DCA_{nd-PV}; Counts", 1000, 0., 100);
  fListOfObjects->Add(hNdcaHigh);

  hDcaV0High = new TH1D("hDcaV0High", "DCA of V0 to primary vertex; DCA_{V0-PV}; Counts", 1000, 0., 10);
  fListOfObjects->Add(hDcaV0High);



  // Post output data.
  PostData(1, fListOfObjects);
}

//______________________________________________________________________________
void AliAnaTaskV0EffDecomposition::UserExec(Option_t *) 
{
  // Main loop
  
  //
  // Comment: This method matches completely the same method for the high pT tracks
  //


  //
  // First we make sure that we have valid input(s)!
  //
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!fAOD){
    Printf("%s:%d AODEvent not found in Input Manager",(char*)__FILE__,__LINE__);
    this->Dump();
    return;
  }

  fMCArray = (TClonesArray*)fAOD->FindListObject("mcparticles");
  if(!fMCArray){
    Printf("%s:%d AOD MC array not found in Input Manager",(char*)__FILE__,__LINE__);
    this->Dump();
    return;
  }    

  // Here we check: 
  // 1) that the event was triggered
  if(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))
     ->IsEventSelected() & fTrigBit ) {

    // 2) if the reconstrcuted vertex was in the correct range
    Float_t zvtx = -999;
    const AliVVertex* primaryVertex = fAOD->GetPrimaryVertex(); 
    if(primaryVertex->GetNContributors()>0)
      zvtx = primaryVertex->GetZ();
    if (TMath::Abs(zvtx) < fVtxCut) {	
      
	  // 3) if the centrality is in the correct range
	  Float_t centrality = -10;
	  AliCentrality *centObject = fAOD->GetCentrality();
	  if(centObject)
	    centrality = centObject->GetCentralityPercentile("V0M"); 
	  if((centrality < fMaxCent) && (centrality>=fMinCent)) {
	    
	    hCentr->Fill(centrality);

	    // as we focus here on Pb-Pb data I decide that the simplest would be
	    // to demand that all 3 are ok before going on
	    
	    // Fill MC gen histograms
	    ProcessMCTruthAOD();
	    
	    // Fill V0 histograms (loop over V0s)
	    AnalyzeV0AOD();
	    
	    // Fill track pair histograms (loop over tracks)
	    AnalyzeDaughtersAOD();	
	    
	    // In this final step we try to see for each generated V0 how many
	    // times it was reconstructed as a V0 and/or a track pair
	    AnalyzeRecMothersAOD();
      }
    }
  }
  
  // Post output data.
  PostData(1, fListOfObjects);
}

//_____________________________________________________________________________
void AliAnaTaskV0EffDecomposition::ProcessMCTruthAOD() 
{
  // Fill the special MC histogram with the MC truth info

  const Int_t nTracksMC = fMCArray->GetEntriesFast();
  
  for (Int_t iTracks = 0; iTracks < nTracksMC; iTracks++) {
    
    AliAODMCParticle* trackMC = dynamic_cast<AliAODMCParticle*>(fMCArray->At(iTracks));
    // We will use the Generator Index as a counter to see what we reconstruct
    trackMC->SetGeneratorIndex(0);

 
    // TString genname;
    // Bool_t yesno=fMCArray->GetCocktailGenerator(iTracks, genname);
    // if(yesno) {
    //   if(!genname.Contains("Hijing")) continue;
    // }
    

    Int_t pdgCode = trackMC->PdgCode();
    if(pdgCode != fPdgV0)
      continue;
    
    // Select only primaries
    if(!(trackMC->IsPhysicalPrimary()))
      continue;
    
    if (TMath::Abs(trackMC->Eta()) > fEtaCut )
      continue;

    Float_t ptMC      = trackMC->Pt();
    hV0Gen->Fill(ptMC);
  }
}




//____________________________________________________________________
AliAODMCParticle* AliAnaTaskV0EffDecomposition::FindPrimaryMotherAOD(AliAODMCParticle* startParticle, Int_t& nSteps)
{
  //
  // Finds the first mother among the primary particles of the particle
  // identified by <label>, i.e. the primary that "caused" this particle
  //
  // Taken from AliPWG0Helper class
  //

  nSteps = 0;

  AliAODMCParticle* mcPart = startParticle;

  while (mcPart) {
    
    if(mcPart->IsPrimary())
      return mcPart;
    
    Int_t mother = mcPart->GetMother();
    
    mcPart = dynamic_cast<AliAODMCParticle*>(fMCArray->At(mother));
    nSteps++; // 1 level down
  }
  
  return 0;
}

//_____________________________________________________________________
void AliAnaTaskV0EffDecomposition::AnalyzeV0AOD() {
  
  Int_t nV0s = fAOD->GetNumberOfV0s();
  if(nV0s < 1)
    return;
  
  // Check that a primary vertex was reconstructed
  AliAODVertex *myBestPrimaryVertex = fAOD->GetPrimaryVertex();
  if (!myBestPrimaryVertex) return;
  
  
  for (Int_t iV0 = 0; iV0 < nV0s; iV0++) {
    
    AliAODv0 *aodV0 = fAOD->GetV0(iV0);
    if (!aodV0) continue;
 

    Double_t ct = ((aodV0->DecayLength(myBestPrimaryVertex))*(1.115683))/aodV0->P();


    if(aodV0->Pt()>1.0 && aodV0->Pt()<1.3){
      hDecayR->Fill(aodV0->RadiusV0());           
      hCosPA->Fill(aodV0->CosPointingAngle(myBestPrimaryVertex));
      hCt->Fill(ct);     
      hDCAdaugh->Fill(aodV0->DcaV0Daughters());   
      hInvMass->Fill(aodV0->MassLambda()-1.116);
      hPdca->Fill(aodV0->DcaPosToPrimVertex());
      hNdca->Fill(aodV0->DcaNegToPrimVertex());
      hDcaV0->Fill(aodV0->DcaV0ToPrimVertex());
    }

    if(aodV0->Pt()>4.5 && aodV0->Pt()<6.5){
      hDecayRHigh->Fill(aodV0->RadiusV0());           
      hCosPAHigh->Fill(aodV0->CosPointingAngle(myBestPrimaryVertex));
      hCtHigh->Fill(ct);     
      hDCAdaughHigh->Fill(aodV0->DcaV0Daughters());   
      hInvMassHigh->Fill(aodV0->MassLambda()-1.116);
      hPdcaHigh->Fill(aodV0->DcaPosToPrimVertex());
      hNdcaHigh->Fill(aodV0->DcaNegToPrimVertex());
      hDcaV0High->Fill(aodV0->DcaV0ToPrimVertex());
    }


 
  if(TMath::Abs(aodV0->Eta())>fEtaCut)continue;
    if(aodV0->RadiusV0()<fDecayRmin || aodV0->RadiusV0()>fDecayRmax)continue;
    // if(TMath::Abs(aodV0->MassLambda()-1.115683)>fMassCut)continue;
    if(aodV0->CosPointingAngle(myBestPrimaryVertex)<fCospt)continue;
    if(aodV0->DcaV0Daughters()>fDcaDaugh)continue;
    if(ct>fCt*7.89)continue;
    if(aodV0->DcaPosToPrimVertex()<fPdca || aodV0->DcaNegToPrimVertex()<fNdca)continue;
    //    if(aodV0->DcaV0ToPrimVertex()<fDcaV0)continue;
    
    //Reject on-the-fly tracks too
    if(aodV0->GetOnFlyStatus() != 0 ) continue;

  
    //Check if V0 is primary by analyzing its daughters:

    // AliAODTrack (V0 Daughters)
    AliAODVertex* vertex = aodV0->GetSecondaryVtx();
    if (!vertex) {
      Printf("ERROR: Could not retrieve vertex");
      continue;
    }
    
    AliAODTrack *p_track = (AliAODTrack*)vertex->GetDaughter(0);
    AliAODTrack *n_track = (AliAODTrack*)vertex->GetDaughter(1);
    if (!p_track || !n_track) {
      Printf("ERROR: Could not retrieve one of the daughter track");
      continue;
    }
    
    // Remove like-sign
    if (p_track->Charge() == n_track->Charge()) {
      continue;
    } 
    
    // Make sure charge ordering is ok
    if (p_track->Charge() < 0) {
      AliAODTrack* helpTrack = p_track;
      p_track = n_track;
      n_track = helpTrack;
    } 
    
      if(TMath::Abs(p_track->Eta()) > fEtaCut || TMath::Abs(n_track->Eta()) > fEtaCut)continue;
    
    //not for mc, right?
    // if(!aodTrack->IsOn(AliAODTrack::kTPCrefit))continue;
    //but maybe;???:
    // if(!p_track->IsOn(AliAODTrack::kTPCrefit))continue;



    AliAODMCParticle* p_mother = ValidateTrack(p_track, fPdgPos);
    if(!p_mother)
      continue;
    AliAODMCParticle* n_mother = ValidateTrack(n_track, fPdgNeg);
    if(!n_mother)
      continue;
    
    // check that mother is the same
    if(p_mother != n_mother)
      continue;

    // check that mother has good eta
    if (TMath::Abs(p_mother->Eta()) > fEtaCut )
      continue;
    
    // One could also consider to fill the pT of the reconstructed mother
    hV0Rec->Fill(p_mother->Pt());

    p_mother->SetGeneratorIndex(p_mother->GetGeneratorIndex() + 1);
  }//end loop over v0's
}


//________________________________________________________________________
AliAODMCParticle* AliAnaTaskV0EffDecomposition::ValidateTrack(AliAODTrack* track, Int_t pdgDaughter)
{
  // Validate V0 daughter track
  
  //Apply track quality cuts
  // if(!track->TestFilterBit(fTrackFilterBit)) {
  //   return 0;
  // }


 if(track->Pt()>1.0 && track->Pt()<1.3){
    hNcl->Fill(track->GetTPCNcls());
    hChi2perNDF->Fill(track->Chi2perNDF());
  }

  if(track->Pt()>4.5 && track->Pt()<6.5){
    hNclHigh->Fill(track->GetTPCNcls());
    hChi2perNDFHigh->Fill(track->Chi2perNDF());
  }
  
  //standard primary track cuts:
  if(track->GetTPCNcls()<fNcl){return 0;}
  if(track->Chi2perNDF()>fChi2perNDF){return 0;}
  
  // track->SetMinNClustersTPC(50);
  // track->SetMaxChi2PerClusterTPC(4);
  // track->SetMaxDCAToVertexZ(3.2);
  // track->SetMaxDCAToVertexXY(2.4);
  // track->SetDCAToVertex2D(kTRUE);
  // track->SetAcceptKinkDaughters(kFALSE);  



  const Int_t label = TMath::Abs(track->GetLabel());
  AliAODMCParticle* mcTrack = dynamic_cast<AliAODMCParticle*>(fMCArray->At(label));
  if (!mcTrack)
    return 0;
  if(mcTrack->IsPhysicalPrimary())
    return 0;
  Int_t pdgCode = mcTrack->GetPdgCode();
  if(pdgCode != pdgDaughter)
    return 0;
  
  // TString genname;
  // Bool_t yesno=fMCArray->GetCocktailGenerator(label, genname);
  //   if(yesno) {
  //     if(!genname.Contains("Hijing")) continue;
  //   }
 
  // mother_steps is the number of steps we have to go backe to find the
  // primary mother. Here we only accept primary V0s so it has to be 1
  Int_t mother_steps = 0;
  AliAODMCParticle* mother = FindPrimaryMotherAOD(mcTrack, mother_steps);
  if(!mother) 
    return 0;
  if(mother_steps != 1)
    return 0;
  Int_t pdgMother = mother->GetPdgCode();
  if(pdgMother != fPdgV0)
    return 0;

  if(!(mother->IsPhysicalPrimary()))
    return 0;

  return mother;
}



//________________________________________________________________________
void AliAnaTaskV0EffDecomposition::AnalyzeDaughtersAOD()
{
  // check to see if two daughter tracks of the same mother was reconstructed
  // (here there is no V0 requirement)
  const Int_t nAODTracks = fAOD->GetNumberOfTracks();
  // Bool_t wasUsed[nAODTracks];
  // for(Int_t iT = 0; iT < nAODTracks; iT++)
  //   wasUsed[iT] = kFALSE;
  
  for(Int_t iT = 0; iT < nAODTracks; iT++) {
    
    // if(wasUsed[iT] == kTRUE)
    //   continue; // already part of a matched pair
    AliAODTrack* track1 = (AliAODTrack*)fAOD->GetTrack(iT);
    
    Int_t charge  = track1->Charge();

    Int_t pdgDaughter1 = fPdgPos;
    Int_t pdgDaughter2 = fPdgNeg;

    if(charge < 0) {

       pdgDaughter1 = fPdgNeg;
       pdgDaughter2 = fPdgPos;
    }
    
    AliAODMCParticle* mother1 = ValidateTrack(track1, pdgDaughter1);
    if(!mother1)
      continue;

    
    // loop over the remaining tracks and see of the other daughter track was
    // also reconstructed
    for(Int_t jT = iT+1; jT < nAODTracks; jT++) {
      
      // if(wasUsed[jT] == kTRUE)
      // 	continue; // already part of a matched pair
      AliAODTrack* track2 = (AliAODTrack*)fAOD->GetTrack(jT);
    
      AliAODMCParticle* mother2 = ValidateTrack(track2, pdgDaughter2);
      if(!mother2)
	continue;
      
      // check that mother is the same
      if(mother1 != mother2)
	continue;

      //      wasUsed[jT] = kTRUE;
      
      // check that mother has good eta
      if (TMath::Abs(mother1->Eta()) < fEtaCut ) {
	
	hDaughterRec->Fill(mother1->Pt());
	mother1->SetGeneratorIndex(mother1->GetGeneratorIndex() + 100);

      }
    }    
  }
}
  
//_____________________________________________________________________________
void AliAnaTaskV0EffDecomposition::AnalyzeRecMothersAOD() 
{
  // See how many mothers were reconstruted and how many times

  const Int_t nTracksMC = fMCArray->GetEntriesFast();
  
  for (Int_t iTracks = 0; iTracks < nTracksMC; iTracks++) {
    
    AliAODMCParticle* trackMC = dynamic_cast<AliAODMCParticle*>(fMCArray->At(iTracks));

    Int_t pdgCode = trackMC->GetPdgCode();
    if(pdgCode != fPdgV0)
      continue;
    
    // Select only primaries
    if(!(trackMC->IsPhysicalPrimary()))
      continue;
    
    if (TMath::Abs(trackMC->Eta()) > fEtaCut )
      continue;
    
    Float_t ptMC      = trackMC->Pt();
      
    Int_t nV0rec = trackMC->GetGeneratorIndex()%100;
    Int_t nTrackrec = Int_t(trackMC->GetGeneratorIndex()/100);

    if(nV0rec > 1)
      hV0Ghost->Fill(ptMC, nV0rec-1);
    if(nTrackrec > 1)
      hTrackGhost->Fill(ptMC, nTrackrec-1);
    if(nV0rec==1 && nTrackrec==0)
      hV0ButNoTracks->Fill(ptMC);      
  }
}

