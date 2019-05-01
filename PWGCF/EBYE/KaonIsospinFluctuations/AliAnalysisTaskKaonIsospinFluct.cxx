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



#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TList.h"
#include "TPDGCode.h"
#include "THnSparse.h"
#include "TRandom3.h"


#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliStack.h"
#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliMCVertex.h"
#include "AliVTrack.h"
#include "AliExternalTrackParam.h"
#include "AliESDfriendTrack.h"

#include "AliAODEvent.h"
#include "AliAODv0.h"
#include "AliAODcascade.h"

#include "AliPID.h"   
#include "AliPIDResponse.h"
#include "AliAODPid.h"


#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"

#include "AliInputEventHandler.h"
#include "AliAODpidUtil.h"
#include "AliAnalysisTaskKaonIsospinFluct.h"

ClassImp(AliAnalysisTaskKaonIsospinFluct)

//________________________________________________________________________
AliAnalysisTaskKaonIsospinFluct::AliAnalysisTaskKaonIsospinFluct() // All data members should be initialised here
:AliAnalysisTaskSE(), fcutCosPa(0.999),fcutcTauMin(-999), fcutNcTauMax(4.0), fcutRapidity(0.5), fcutArmenteros(0.22),fcutDCA(1.0), fcutNsigma(3.0),
fOutput(0),
fPIDResponse(0),
fAOD(0), 
fPrimaryVtx(0),fTrackBuffSize(18000), 
fHistGoodEvent(0),
fHistPt(0), 
fHistEta(0),
fHistNV0(0),
fHistZVertexCent(0),
fHistCosPaMK0(0),	
fHistcTauMK0(0),		
fHistDcaMK0(0),	
fHistRapMK0(0),	
fHistArmPodK0(0),
fHistoCorrelation(0),
fKshortSparse(0),
fdEdXPID(0),
fdEdXnoPID(0),  
fnsigmakaon(0),
fnsigmaproton(0),
fnsigmapion(0)



// The last in the above list should not have a comma after it
{
    // Dummy constructor ALWAYS needed for I/O.
}

//________________________________________________________________________
AliAnalysisTaskKaonIsospinFluct::AliAnalysisTaskKaonIsospinFluct(const char *name) // All data members should be initialised here
:AliAnalysisTaskSE(name), fcutCosPa(0.999),fcutcTauMin(-999), fcutNcTauMax(4.0), fcutRapidity(0.5),fcutArmenteros(0.22),
fOutput(0), fcutNsigma(3.0),
fcutDCA(1.0),
fPIDResponse(0),
fAOD(0), 
fPrimaryVtx(0),fTrackBuffSize(18000), 
fHistGoodEvent(0),
fHistPt(0), 
fHistEta(0),
fHistNV0(0),
fHistZVertexCent(0), 
fHistCosPaMK0(0),	
fHistcTauMK0(0),		
fHistDcaMK0(0),	
fHistRapMK0(0),	
fHistArmPodK0(0),
fHistoCorrelation(0),
fKshortSparse(0),
fdEdXPID(0),
fdEdXnoPID(0), 
fnsigmakaon(0),
fnsigmaproton(0),
fnsigmapion(0)



// The last in the above list should not have a comma after it
{

    DefineOutput(1, TList::Class());                                            // for output list
}

//________________________________________________________________________
AliAnalysisTaskKaonIsospinFluct::~AliAnalysisTaskKaonIsospinFluct()
{

    if (fOutput) delete fOutput;    

}

//________________________________________________________________________
void AliAnalysisTaskKaonIsospinFluct::UserCreateOutputObjects()
{
    // Create histograms
    // Called once (on the worker node)
	
  

  fOutput = new TList();
  fOutput->SetOwner();  // IMPORTANT!
	

  fHistGoodEvent = new TH1F("hGoodEvent","No of events passing the cuts.",10,   0.5,9.5);
    
	// Create histograms
  Int_t ptbins = 100;
  Float_t ptlow = 0.1, ptup = 40.0;
  fHistPt = new TH1F("fHistPt", "P_{T} distribution for reconstructed K0s", ptbins, ptlow, ptup);
  fHistPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  fHistPt->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  fHistPt->SetMarkerStyle(kFullCircle);
	
  Int_t etabins = 40;
  Float_t etalow = -1.0, etaup = 1.0;
  fHistEta = new TH1F("fHistEta","#eta distribution for reconstructed K0s",etabins, etalow, etaup);
  fHistEta->GetXaxis()->SetTitle("#eta");
  fHistEta->GetYaxis()->SetTitle("counts");
	

  fHistNV0 = new TH1F("fHistNV0","Number of V0s per event",100, 0, 5000);
   
    
  Int_t K0Bins[4]   = {100, 300,  140, 24 };
  Double_t K0Min[4] = {-0.5,  0., 0.414, -1.2};
  Double_t K0Max[4] = { 99.5 ,30.,0.582,1.2};
  
  const Char_t *K0Title[] = {"centrality","Pt","Mass","Eta"};
  
  fKshortSparse = new THnSparseD("fThnKshort", "", 4, K0Bins, K0Min, K0Max);
  for (Int_t iaxis = 0; iaxis < 4; iaxis++)
    fKshortSparse->GetAxis(iaxis)->SetTitle(K0Title[iaxis]);
  
  
  fHistZVertexCent = new TH2F("fHistZVertexCent"," Vz of primary vertex Vs Centrality; V_{Z} {cm} ; Centrality {%}", 60, -15, 15,100,0,100);
  
    
  fHistCosPaMK0		 = new	TH2D("fHistCosPaMK0","	Reconstructed Mass vs CosPa for K0Short Candidates;CosPA; K0 Mass (GeV/c^2)",200,0.99,1.001,140,0.414,0.582);
    
    
  fHistcTauMK0		 = new	TH2D("fHistcTauMK0","	Reconstructed Mass vs cTau for K0Short Candidates; cTau; K0 Mass (GeV/c^2)",200,0,50,140,0.414,0.582);
	
	
  fHistDcaMK0		 = new	TH2D("fHistDcaMK0","	Reconstructed Mass vs Dca for K0Short Candidates; DCA; K0 Mass (GeV/c^2)",200,0,1.5,140,0.414,0.582);
	
	
  fHistRapMK0		 = new	TH2D("fHistRapMK0","	Reconstructed Mass vs Rap for K0Short Candidates; Rapidity; K0 Mass (GeV/c^2)",200,-1,1,140,0.414,0.582);
	
  fHistArmPodK0 = new	TH2D("fHistArmPodK0","Armenteros plot for K0 candidates; Alpha; PtArm",100,-1,1,50,0,0.3);

  fdEdXPID = new TH2D("dEdXPID", "TPC dE/dX", 500, 0.2, 1.5, 500, 0., 1000);
  fdEdXPID->GetXaxis()->SetTitle("Momentum [GeV/c]");
  fdEdXPID->GetYaxis()->SetTitle("Energy loss [MeV]");

  fdEdXnoPID = new TH2D("dEdXnoPID", "TPC dE/dX", 500, 0.2, 1.5, 500, 0., 1000);
  fdEdXnoPID->GetXaxis()->SetTitle("Momentum [GeV/c]");
  fdEdXnoPID->GetYaxis()->SetTitle("Energy loss [MeV]");
	
  fnsigmakaon = new TH2D("K distribution","K distribution",400,0.0,4.0,100,-5.0,5.0);
  fnsigmakaon->GetXaxis()->SetTitle("Momentum[GeV/c]");
  fnsigmakaon->GetYaxis()->SetTitle("nsigmaK");
        
  fnsigmaproton = new TH2D("P distribution","P distribution",400,0.0,4.0,100,-5.0,5.0);
  fnsigmaproton->GetXaxis()->SetTitle("Momentum[GeV/c]");
  fnsigmaproton->GetYaxis()->SetTitle("nsigmaP");
	
  fnsigmapion = new TH2D("pi distribution","pi distribution",400,0.0,4.0,100,-5.0,5.0);
  fnsigmapion->GetXaxis()->SetTitle("Momentum[GeV/c]");
  fnsigmapion->GetYaxis()->SetTitle("nsigmapi");
	
  Int_t fgSparseDataBins[8]   = {100,  2000,   2000,   2000,   2000,   2000, 2000, 2000};
  Double_t fgSparseDataMin[8] = {-0.5, -0.5,   -0.5,   -0.5,   -0.5,   -0.5, -0.5, -0.5};
  Double_t fgSparseDataMax[8] = {99.5, 1999.5, 1999.5, 1999.5, 1999.5, 1999.5, 1999.5, 1999.5};
	
  const Char_t *fgkSparseDataTitle[] = {"centrality","K_{s}^{0}","K^{+}","K^{-}", "#pi^{+}", "#pi^{-}", "K=K^{+}+K^{-}", "#pi=#pi^{+}+#pi^{-}"};
	
  fHistoCorrelation = new THnSparseD("fThnCorr", "", 8, fgSparseDataBins, fgSparseDataMin, fgSparseDataMax);
  for (Int_t iaxis = 0; iaxis < 8; iaxis++)
  fHistoCorrelation->GetAxis(iaxis)->SetTitle(fgkSparseDataTitle[iaxis]);
	
  fOutput->Add(fHistGoodEvent);
  /* fOutput->Add(fHistPt);
     fOutput->Add(fHistEta);
     fOutput->Add(fHistNV0);
     //fOutput->Add(fHistCentrality);
     
     fOutput->Add(fHistBBK0Pos);
     fOutput->Add(fHistBBK0Neg);
     fOutput->Add(fHistBBPion);
     fOutput->Add(fHistZVertexCent);
     fOutput->Add(fHistCosPaMK0);	
     fOutput->Add(fHistcTauMK0);		
     fOutput->Add(fHistDcaMK0);	
     fOutput->Add(fHistRapMK0);	
     fOutput->Add(fHistArmPodK0);
     fOutput->Add(fHistoCorrelation);
     fOutput->Add(fKshortSparse);
     fOutput->Add(fLambdaSparse);
     fOutput->Add(fALambdaSparse);*/
   
  fOutput->Add(fHistoCorrelation);
  fOutput->Add(fKshortSparse);
  fOutput->Add(fHistArmPodK0);
  //fOutput->Add(fdEdXPID);
  //fOutput->Add(fdEdXnoPID);
  /*fOutput->Add(fnsigmakaon);
    fOutput->Add(fnsigmaproton);
    fOutput->Add(fnsigmapion);*/
   
    
	// NEW HISTO added to fOutput here
  PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
}



//________________________________________________________________________

static Bool_t AcceptTrack(const AliAODTrack *trk){

  if (!trk) return kFALSE;
  if(!trk->TestFilterBit(272)) return kFALSE;
  if(trk->Charge() == 0) return kFALSE;
  if (trk->Eta() < -0.5 || trk->Eta() > 0.5) return kFALSE;
  if (trk->Pt() < 0.2 || trk->Pt() > 1.5)  return kFALSE;
	
  return kTRUE;
}

//________________________________________________________________________
static Bool_t AcceptV0_general(const AliAODv0 *v1, const AliAODEvent *aod, double fcutCosPa, double fcutNImpact, double fcutDCA, double fcutEta) 
{
	
  if (v1->GetOnFlyStatus()) return kFALSE;
	
  int nnum = 1, pnum = 0;	
 
  const AliAODTrack *ptrack1=(AliAODTrack *)v1->GetDaughter(pnum);
  const AliAODTrack *ntrack1=(AliAODTrack *)v1->GetDaughter(nnum);

	
  Float_t impact=v1->DcaNegToPrimVertex();
  if (TMath::Abs(impact)<0.1) return kFALSE;
  if (TMath::Abs(impact)<fcutNImpact && fcutNImpact != -999) return kFALSE;
  impact=v1->DcaPosToPrimVertex();
  if (TMath::Abs(impact)<0.1) return kFALSE;
  if (TMath::Abs(impact)<fcutNImpact && fcutNImpact != -999) return kFALSE;
	
  Double_t dca=v1->DcaV0Daughters();
  if (TMath::Abs(dca)>fcutDCA && fcutDCA != -999) return kFALSE;
	
  Double_t cpa=v1->CosPointingAngle(aod->GetPrimaryVertex());
  if (cpa<fcutCosPa && fcutCosPa != -999) return kFALSE;
	
  Double_t etaN = v1->PseudoRapNeg();
  Double_t etaP = v1->PseudoRapPos();
  if ((TMath::Abs(etaN)>fcutEta || TMath::Abs(etaP)>fcutEta) && fcutEta != -999) return kFALSE;
	
  return kTRUE;
}

//________________________________________________________________________
static Bool_t AcceptV0_particle(const AliAODv0 *v1, double fcutcTauMin, double fcutRapidity, Double_t decayL, double fcutNcTauMax) 
{
  Double_t cTau = 0;
  cTau = decayL*(v1->MassK0Short())/(v1->P());
 
  if (cTau < fcutcTauMin && fcutcTauMin != -999 ) return kFALSE;
  if (cTau > (fcutNcTauMax*2.68) && fcutNcTauMax != -999) return kFALSE;
  
  Double_t rap = 0;
  rap = v1->RapK0Short();
	  
  if (TMath::Abs(rap)>fcutRapidity && fcutRapidity != -999) return kFALSE;
  return kTRUE;
}

//________________________________________________________________________

void AliAnalysisTaskKaonIsospinFluct::UserExec(Option_t *) 
{
  // Main loop                                                                
  // Called for each event                                                    


  // Fill a control histogram                                                  
  double fcutNImpact(-999), fcutDCA(-999);
  fHistGoodEvent->Fill(0.0);

  // Get the event                                                             
                                                                                
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fAOD) {
    printf("ERROR: fAOD not available\n");
    return;

  }
    // Fill a control histogram                                                
                                                                               
    fHistGoodEvent->Fill(1.0);
    AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();

    if (!(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kMB)) return;


    // centrality selection
    
    AliCentrality* centrality = fAOD->GetCentrality();
    if (!centrality) {
      printf ("ERROR: couldn't get the AliCentrality\n");
      return;
    }
    Double_t centPercentile = centrality->GetCentralityPercentile("V0M");
                                                                               
    fHistGoodEvent->Fill(2.0);

    // primary vertex
    
    fPrimaryVtx = fAOD->GetPrimaryVertex();
    if (!fPrimaryVtx){
      printf ("ERROR: no primary vertex\n");
      return;
    }
                                   
    fHistGoodEvent->Fill(3.0);

    Double_t xv=fPrimaryVtx->GetX();
    Double_t yv=fPrimaryVtx->GetY();
    Double_t zv=fPrimaryVtx->GetZ();

    // Vz cut
    
    if (TMath::Abs(zv) > 10.0)
      {
	return;
      }

    fHistGoodEvent->Fill(4.0);
    fHistZVertexCent->Fill(zv, centPercentile);
    
    Int_t nIDs[2] = {-9999999};

    TObjArray *posarray = new TObjArray();
     TObjArray *negarray = new TObjArray();
     
    Int_t shared_dau = 0;
    Int_t kshort = 0;
    Int_t kshort0 = 0;
    Int_t kshort1 = 0;
    Int_t kshort2 = 0;
    Int_t kshort3 = 0;
    Int_t nkshort0 = 0;
    Int_t nkshort = 0;
    Int_t nkshortB = 0;
    Int_t nkshortlB = 0;
    Int_t nkshortrB = 0;
    Int_t ptks1 = 0; Int_t ptks2 = 0; Int_t ptks3 = 0;
    Int_t ptks4 = 0; Int_t ptks5 = 0; Int_t ptks6 = 0;
    Int_t ptks7 = 0; Int_t ptks8 = 0; Int_t ptks9 = 0;
    Int_t ptks10 = 0; Int_t ptks11 = 0; Int_t ptks12 = 0;
    Int_t ptks13 = 0; Int_t ptks14 = 0;

    // count the number of K0's in the event

    
    Int_t nv0s = fAOD->GetNumberOfV0s();

    fHistNV0->Fill(nv0s);
    
    //cout<<"number of V0"<<     nv0s<<endl;

    for(Int_t i = 0; i < nv0s; i++) 
      {
	AliAODv0 *v0=fAOD->GetV0(i); // pointer to reconstructed v0          
	if(!v0) 
	  { 
	    cout<<"No V0 "<<endl;
	    continue; 
	  }
	  
     // V0s not consistent with Kshort

	Int_t ppid = -999; Int_t npid = -999;	  
	  
   
	Bool_t k0Candidate = kTRUE;
	Bool_t k0Cand = kFALSE;
	

	
	if (v0->MassK0Short() < 0.414 || v0->MassK0Short() > 0.582)
	  { k0Candidate = kFALSE;}
	 
	if(!k0Candidate) continue;

	Double_t cosPA=v0->CosPointingAngle(fAOD->GetPrimaryVertex());
	Double_t xyz[3]; 
	  
	v0->GetSecondaryVtx(xyz);
	  
	Double_t decayL = TMath::Sqrt((xyz[0]-xv)*(xyz[0]-xv)+(xyz[1]-yv)*(xyz[1]-yv)+(xyz[2]-zv)*(xyz[2]-zv));
	Double_t decayRad = TMath::Sqrt((xyz[0]-xv)*(xyz[0]-xv)+(xyz[1]-yv)*(xyz[1]-yv);

					/*	if ( decayRad < 0.5 || decayrad > 100.0) {
					return;
					}*/
		  
	if(!AcceptV0_general(v0,fAOD,fcutCosPa,fcutNImpact,fcutDCA,fcutEta))
	  { 
	    continue; 
	  }
	  
	int nnum = 1, pnum = 0;


	AliAODTrack *ptrack1=(AliAODTrack *)v0->GetDaughter(pnum);
	AliAODTrack *ntrack1=(AliAODTrack *)v0->GetDaughter(nnum);

	if(ntrack1->Charge()==0 || ptrack1->Charge()==0) continue;

	if(!AcceptV0_particle(v0,fcutcTauMin, fcutRapidity, decayL, fcutNcTauMax))
	  { k0Candidate = kFALSE; }


	if( !k0Candidate ) continue;
					
    	//Get the pid of the daughter tracks

	if (ptrack1->P() < 0.5){
	ppid = ProcessTPC(ptrack1);
	}
	else if (ptrack1->P() < 0.7){
	ppid = ProcessHybridPro(ptrack1);
	}
	else if (ptrack1->P() < 100.0) {
	ppid = ProcessHybrid(ptrack1);
    }
	  
	  
        if (ntrack1->P() < 0.5){
	npid = ProcessTPC(ntrack1);
	}
	else if (ntrack1->P() < 0.7){
	npid = ProcessHybridPro(ntrack1);
	}
	else if (ntrack1->P() < 100.0) {
	npid = ProcessHybrid(ntrack1);
	}
   	    
	if ((ppid == 1) && ( npid == 1)){

	  k0Cand = kTRUE;

	  nkshort++;            
      
    }

    Double_t dca=v0->DcaV0Daughters();
    double ArmenterosAlpha = v0->Alpha();
    double ArmenterosPt	   = v0->QtProng();
    
    if( ArmenterosPt <= TMath::Abs(fcutArmenteros*ArmenterosAlpha) && fcutArmenteros !=-999 )
    {k0Candidate = false;}
    if(!k0Candidate ) continue;
    
    	//pT cut for Kshort
	if (TMath::Sqrt(v0->Pt2V0()) < 0.4 || TMath::Sqrt(v0->Pt2V0()) > 1.5 )
	  {
	    return;
	  }
	if (v0->PseudoRapV0() < -0.5 || v0->PseudoRapV0() > 0.5)
	  {return;}
 					
       if(k0Cand){

       if (v0->MassK0Short() > 0.490 && v0->MassK0Short() < 0.506)
           {
            kshort3++;
	     posarray->Add(ptrack1);
	     negarray->Add(ntrack1);
	    
	      if ( TMath::Sqrt(v0->Pt2V0()) >= 0.19 && TMath::Sqrt(v0->Pt2V0()) < 0.21 ){ptks1++;}
	      if ( TMath::Sqrt(v0->Pt2V0()) >= 0.21 && TMath::Sqrt(v0->Pt2V0()) < 0.31 ){ptks2++;}
	      if ( TMath::Sqrt(v0->Pt2V0()) >= 0.31 && TMath::Sqrt(v0->Pt2V0()) < 0.41 ){ptks3++;}
	      if ( TMath::Sqrt(v0->Pt2V0()) >= 0.41 && TMath::Sqrt(v0->Pt2V0()) < 0.51 ){ptks4++;}
	     
	      if ( TMath::Sqrt(v0->Pt2V0()) >= 0.51 && TMath::Sqrt(v0->Pt2V0()) < 0.61 ){ptks5++;}
	      
	      if ( TMath::Sqrt(v0->Pt2V0()) >= 0.61 && TMath::Sqrt(v0->Pt2V0()) < 0.71 ){ptks6++;}
	      
	      if ( TMath::Sqrt(v0->Pt2V0()) >= 0.71 && TMath::Sqrt(v0->Pt2V0()) < 0.81 ){ptks7++;}
	     
	      if ( TMath::Sqrt(v0->Pt2V0()) >= 0.81 && TMath::Sqrt(v0->Pt2V0()) < 0.91 ){ptks8++;}
	      
	      if ( TMath::Sqrt(v0->Pt2V0()) >= 0.91 && TMath::Sqrt(v0->Pt2V0()) < 1.01 ){ptks9++;}
	     
	      if ( TMath::Sqrt(v0->Pt2V0()) >= 1.01 && TMath::Sqrt(v0->Pt2V0()) < 1.11 ){ptks10++;}
	      
	      if ( TMath::Sqrt(v0->Pt2V0()) >= 1.11 && TMath::Sqrt(v0->Pt2V0()) < 1.21 ){ptks11++;}
	      
	      if ( TMath::Sqrt(v0->Pt2V0()) >= 1.21 && TMath::Sqrt(v0->Pt2V0()) < 1.31 ){ptks12++;}
	     
	      if ( TMath::Sqrt(v0->Pt2V0()) >= 1.31 && TMath::Sqrt(v0->Pt2V0()) < 1.41 ){ptks13++;}
	      
	      if ( TMath::Sqrt(v0->Pt2V0()) >= 1.41 && TMath::Sqrt(v0->Pt2V0()) < 1.51 ){ptks14++;}
	      
	     
	   }
    
      if (v0->MassK0Short() > 0.420 && v0->MassK0Short() < 0.47)
	{nkshortlB++;}

      if (v0->MassK0Short() > 0.53 && v0->MassK0Short() < 0.575)
	{nkshortrB++;}
	  	  	       
      fHistCosPaMK0->Fill(cosPA,v0->MassK0Short()); 
      fHistcTauMK0->Fill(decayL*(v0->MassK0Short())/(v0->P()),v0->MassK0Short()); 	
      fHistDcaMK0->Fill(dca,v0->MassK0Short()); 	
      fHistRapMK0->Fill(v0->RapK0Short(),v0->MassK0Short()); 
      fHistArmPodK0->Fill(ArmenterosAlpha,ArmenterosPt);
	      	      
      Double_t k0sparse[4];
      k0sparse[0]   = centPercentile;
      k0sparse[1]   = TMath::Sqrt(v0->Pt2V0());
      k0sparse[2]   = v0->MassK0Short();
      k0sparse[3]   = v0->PseudoRapV0();
	      
      fKshortSparse->Fill(k0sparse);
    }//Kshort loop end
	      
      }  //v0 loop ends

    //calculating efficiency corrected kshort



     Double_t PtKshort_uncorr[11] = {ptks4 ,ptks5 ,ptks6 ,ptks7 ,ptks8 ,ptks9 ,ptks10 ,ptks11 ,ptks12 ,ptks13 ,ptks14};
     Double_t effkshort[11] = {0.0186395, 0.029414, 0.0414794, 0.0547837, 0.0673469, 0.0792357, 0.090228, 0.10008, 0.108047, 0.115824, 0.119694};

	 Int_t kshort2sqr  = 0;
	 Int_t PtKshort_corr[11] = {0};
	 Int_t PtKshort_corrsqr[11] = {0};
	 
	      

	 /*	 for (int i = 0; i < 11; i++){

	 PtKshort_corr[i] = PtKshort_uncorr[i]/effkshort[i];//calculate number of kshort corrected
	
	 }
	 for (int ii = 0; ii <11; ii++)
		  {
		    kshort2+=PtKshort_corr[ii];
		  }
       	 for (int i = 0; i < 11; i++)
         {
	 
	 for (int j = 0; j <11; j++)
	   {
	   if(i == j)
	     {
	  kshort2sqr  +=  (PtKshort_uncorr[j]*(PtKshort_uncorr[j] -1))/(effkshort[j]*effkshort[j]);
	       }
	     else
	       kshort2sqr  += PtKshort_corr[i]*(PtKshort_corr[j] -1);
	   }
	 if (kshort2sqr < 0)  {return;}
       }*/
     //cout<<"ptksouncorr = "<<PtKshort_uncorr[11]<< "ptksocorr = "<<PtKshort_corr[11]<<endl;	

     //Check for daughter sharing
     for(Int_t ip=0; ip<kshort3; ip++) {
      
	  AliAODTrack *pos = (AliAODTrack*)posarray->At(ip);
	  AliAODTrack *neg = (AliAODTrack*)negarray->At(ip);
	  
	  for(Int_t ipp=kshort3-1; ipp>ip;ipp--) {
	    AliAODTrack *pos1 = (AliAODTrack *)posarray->At(ipp);
	    AliAODTrack *neg1 = (AliAODTrack *)negarray->At(ipp);

	    Int_t idpos= pos->GetID();
	    Int_t idpos1= pos1->GetID();
	    //  cout<<"****************************"<<endl;
	    Int_t idneg= neg->GetID();
	    Int_t idneg1= neg1->GetID();
	    // cout<<idpos<<"*** idpos"<<idneg<<"**** idneg"<<idpos1<<"*** idpos1"<<idneg1<<endl;
	    if (idpos == idpos1 || idneg == idneg1) {shared_dau++;}
      }
      }

     kshort2 = kshort3 - shared_dau; //rejecting the shared daughters.	 
	    
    //track loop

    Int_t ntracks = fAOD->GetNumberOfTracks();
    
    //cout<<ntracks<<"*****"<<endl;
    
    Int_t nkaons = 0;
    Int_t npions = 0;
    Int_t nprotons = 0;
    
    Int_t nkaonsp = 0;
    Int_t npionsp = 0;
    Int_t nprotonsp = 0;
    
    Int_t nkaonsn = 0;
    Int_t npionsn = 0;
    Int_t nprotonsn = 0;
    Int_t ptkp1 = 0; Int_t ptkp2 = 0; Int_t ptkp3 = 0;
    Int_t ptkp4 = 0; Int_t ptkp5 = 0; Int_t ptkp6 = 0;
    Int_t ptkp7 = 0; Int_t ptkp8 = 0; Int_t ptkp9 = 0;
    Int_t ptkp10 = 0; Int_t ptkp11 = 0; Int_t ptkp12 = 0;
    Int_t ptkp13 = 0; Int_t ptkp14 = 0;

    Int_t ptkn1 = 0; Int_t ptkn2 = 0; Int_t ptkn3 = 0;
    Int_t ptkn4 = 0; Int_t ptkn5 = 0; Int_t ptkn6 = 0;
    Int_t ptkn7 = 0; Int_t ptkn8 = 0; Int_t ptkn9 = 0;
    Int_t ptkn10 = 0; Int_t ptkn11 = 0; Int_t ptkn12 = 0;
    Int_t ptkn13 = 0; Int_t ptkn14 = 0;
    

    Int_t a = -9999;
  
    for(Int_t itrk = 0; itrk < ntracks; itrk++){
      
    AliAODTrack* track = dynamic_cast<AliAODTrack *>(fAOD->GetTrack(itrk));
	 
    if (!track) continue;
    if (!AcceptTrack(track)) continue;
  
    fHistPt->Fill(track->Pt());
    fHistEta->Fill(track->Eta());

    if (track->P() < 0.5){
    a = ProcessTPC(track);
    }
    else if (track->P() < 0.7){
      a = ProcessHybridPro(track);
    }
    else if (track->P() < 100.0) {

    a =	ProcessHybrid(track);
    }
    
      
    if(a < 1 || a > 3) continue;

    if(a == 1) {
      npions++;
	
      if (track->Charge()>0) {npionsp++ ;}
      if (track->Charge()<0) {npionsn++ ;}
    }
      
    if(a == 2) {

              

      if (track->Charge()>0) {
	     if ( track->Pt() >= 0.19 && track->Pt() < 0.21  ){ptkp1++;}
	     if ( track->Pt() >= 0.21 && track->Pt() < 0.31  ){ptkp2++;}
	     if ( track->Pt() >= 0.31 && track->Pt() < 0.41  ){ptkp3++;}
	     if ( track->Pt() >= 0.41 && track->Pt() < 0.51  ){ptkp4++;}
	     if ( track->Pt() >= 0.51 && track->Pt() < 0.61  ){ptkp5++;}
	     if ( track->Pt() >= 0.61 && track->Pt() < 0.71  ){ptkp6++;}
	     if ( track->Pt() >= 0.71 && track->Pt() < 0.81  ){ptkp7++;}
	     if ( track->Pt() >= 0.81 && track->Pt() < 0.91  ){ptkp8++;}
	     if ( track->Pt() >= 0.91 && track->Pt() < 1.01  ){ptkp9++;}
	     if ( track->Pt() >= 1.01 && track->Pt() < 1.11  ){ptkp10++;}
	     if ( track->Pt() >= 1.11 && track->Pt() < 1.21  ){ptkp11++;}
	     if ( track->Pt() >= 1.21 && track->Pt() < 1.31  ){ptkp12++;}
	     if ( track->Pt() >= 1.31 && track->Pt() < 1.41  ){ptkp13++;}
	     if ( track->Pt() >= 1.41 && track->Pt() < 1.51  ){ptkp14++;}
	      
	      
	     
	nkaonsp++ ;
      }
      
      if (track->Charge()<0) {
	      if ( track->Pt() >= 0.19 && track->Pt() < 0.21  ){ptkn1++;}
	      if ( track->Pt() >= 0.21 && track->Pt() < 0.31  ){ptkn2++;}
	      if ( track->Pt() >= 0.31 && track->Pt() < 0.41  ){ptkn3++;}
	      if ( track->Pt() >= 0.41 && track->Pt() < 0.51  ){ptkn4++;}
	      if ( track->Pt() >= 0.51 && track->Pt() < 0.61  ){ptkn5++;}
	      if ( track->Pt() >= 0.61 && track->Pt() < 0.71  ){ptkn6++;}
	      if ( track->Pt() >= 0.71 && track->Pt() < 0.81  ){ptkn7++;}
	      if ( track->Pt() >= 0.81 && track->Pt() < 0.91  ){ptkn8++;}
	      if ( track->Pt() >= 0.91 && track->Pt() < 1.01  ){ptkn9++;}
	      if ( track->Pt() >= 1.01 && track->Pt() < 1.11  ){ptkn10++;}
	      if ( track->Pt() >= 1.11 && track->Pt() < 1.21  ){ptkn11++;}
	      if ( track->Pt() >= 1.21 && track->Pt() < 1.31  ){ptkn12++;}
	      if ( track->Pt() >= 1.31 && track->Pt() < 1.41  ){ptkn13++;}
	      if ( track->Pt() >= 1.41 && track->Pt() < 1.51  ){ptkn14++;}
	      
	      
	      
	nkaonsn++ ;
      }
    }

    if(a == 3) {
  
      nprotons++;
     
      if (track->Charge()>0) {nprotonsp++ ;}
      if (track->Charge()<0) {nprotonsn++ ;}
    }

    } // track loop ends here

    Double_t PtKplus_uncorr[13] = {ptkp2, ptkp3, ptkp4, ptkp5, ptkp6, ptkp7, ptkp8, ptkp9, ptkp10, ptkp11, ptkp12, ptkp13, ptkp14};
     Double_t effkplus[13] = {0.201674, 0.409359, 0.508769, 0.559053, 0.606346, 0.638419,0.665882, 0.686093,0.702675, 0.719366, 0.733534, 0.74120, 0.741186};

	      Int_t kplus2sqr  = 0;
	      Int_t PtKplus_corr[13] = {0};
	      Int_t PtKplus_corrsqr[13] = {0};

	      /*    for (int i = 0; i < 13; i++){

	    // if  PtKplus_corr[i] == 0 {continue;}

	     PtKplus_corr[i] = PtKplus_uncorr[i]/effkplus[i];//calculate number of kplus corrected
	      }

	     for (int ii = 0; ii <13; ii++)
		  {
		    nkaonsp+=PtKplus_corr[ii];//to be added
		    
		  }

	    for (int i = 0; i < 13; i++)
		{

	    for (int j = 0; j <13; j++)
		{
	        if(i == j)
	     	{
	    kplus2sqr  +=  (PtKplus_uncorr[j]*(PtKplus_uncorr[j] -1))/(effkplus[j]*effkplus[j]);
	      	}
	      	else
	    kplus2sqr  += PtKplus_corr[i]*(PtKplus_corr[j] - 1);
		    }
		}

		 Int_t PtKminus_uncorr[13] = {ptkn2, ptkn3, ptkn4, ptkn5, ptkn6, ptkn7, ptkn8, ptkn9, ptkn10, ptkn11, ptkn12, ptkn13, ptkn14};
	     Double_t effkminus[13] = {0.177988, 0.377285, 0.478507, 0.535829, 0.590221, 0.619686, 0.639414, 0.661616, 0.683636, 0.706696, 0.71968, 0.731811, 0.732079};

	      Int_t kminus2sqr  = 0;
	      Int_t PtKminus_corr[13] = {0};
	      Int_t PtKminus_corrsqr[13] = {0};

	       for (int i = 0; i < 13; i++){
		PtKminus_corr[i] = PtKminus_uncorr[i]/effkminus[i];//calculate number of kminus corrected
	
	      }
	      for (int ii = 0; ii <13; ii++)
		  {
		    nkaonsn+=PtKminus_corr[ii]; //to be added

		  }

	      for (int i = 0; i < 13; i++)
		{

		  for (int j = 0; j <13; j++)
		    {
		      if(i == j)
			{
			  kminus2sqr  +=  (PtKminus_uncorr[j]*(PtKminus_uncorr[j] -1))/(effkminus[j]*effkminus[j]);
			}
			else
			  kminus2sqr  += PtKminus_corr[i]*(PtKminus_corr[j] - 1);
		    }
		}

	      Int_t PtKaon_uncorr[13] = {0};
	      
	      for (int i = 0; i < 13; i++){
		PtKaon_uncorr[i] = ( PtKplus_uncorr[i] + PtKminus_uncorr[i]);
	      }

    Double_t effkaon[13] = {0.164324, 0.393735, 0.506215, 0.576249, 0.62681, 0.66166, 0.687436, 0.708101, 0.72790, 0.747362, 0.760914, 0.76933, 0.766969};   

              Int_t kaon2sqr  = 0;
	      Int_t PtKaon_corr[13] = {0};
	      Int_t PtKaon_corrsqr[13] = {0};
              Int_t fnkaons = 0;
	      Int_t fnkaons2 = 0;
	     
	        
	      for (int i = 0; i < 13; i++){
		PtKaon_corr[i] =  PtKaon_uncorr[i]/effkaon[i];//calculate number of kminus corrected
	      }
	      for (int ii = 0; ii <13; ii++){
	        fnkaons+=PtKaon_corr[ii];//to be added
	      }

	      for (int i = 0; i < 13; i++)
		{

		  for (int j = 0; j <13; j++)
		    {
		      if(i == j)
			{
			  kaon2sqr  +=  (PtKaon_uncorr[j]*(PtKaon_uncorr[j] -1))/(effkaon[j]*effkaon[j]);
			}
			else
			  kaon2sqr  += PtKaon_corr[i]*(PtKaon_corr[j] - 1);
		     
		    }
		  if(kaon2sqr < 0) {return;}
		  }*/
      
    Int_t fnkaons =  nkaonsp + nkaonsn;
    Int_t fnpions =  npionsp + npionsn;
    // Int_t fnkaonssqr = kplus2sqr + kminus2sqr;

   
 
   Double_t vsparse[8];

   vsparse[0]   = centPercentile;
   vsparse[1]   = kshort2;//kshort
   vsparse[2]   = nkaonsp;
   vsparse[3]   = nkaonsn;
   vsparse[4]   = npionsp;
   vsparse[5]   = npionsn;
   vsparse[6]   = fnkaons;//total kaons
   vsparse[7]   = fnpions;
   //vsparse[8]   = fnkaonssqr;//total kaon square
   //vsparse[9]   = kshort2sqr;//kshort sqr
   //vsparse[8]   = kshortlB;
   // vsparse[9]   = kshortrB;
   
   
    fHistoCorrelation->Fill(vsparse);
 
	PostData(1, fOutput);
}


//________________________________________________________________________
void AliAnalysisTaskKaonIsospinFluct::Terminate(Option_t *) 
{
	// Draw result to screen, or perform fitting, normalizations
	// Called once at the end of the query
	
	fOutput = dynamic_cast<TList*> (GetOutputData(1));
	if(!fOutput) { Printf("ERROR: could not retrieve TList fOutput"); return; }
	
	
	
	fdEdX1 = dynamic_cast<TH2D*> (fOutput->FindObject("fdEdX1"));
	if (!fdEdX1) { Printf("ERROR: could not retrieve fdEdX1"); return;}
	
	TCanvas *c = new TCanvas("AliAnalysisTaskKaonIsospinFluct","P_{T} & #eta",10,10,1020,510);
	c->Divide(2,1);
	c->cd(1)->SetLogy();
	fHistPt->DrawCopy("E");
      

}




Int_t AliAnalysisTaskKaonIsospinFluct::ProcessTPC(AliAODTrack *track){
 
   Double_t nsigmaelectron = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kElectron)) ;

   Double_t nsigmakaon = fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon) ;
 
   Double_t nsigmaproton =TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton)) ;

   Double_t nsigmapion = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion)) ;
 
  
 
   if( nsigmaelectron  < 2.0 ) { return 6;}  

   //   if(( nsigmakaon==nsigmapion ) && ( nsigmakaon==nsigmaproton )) return ;
   Double_t kaonsym = TMath::Abs(nsigmakaon);
   
   if(track->P() >= 0.39 && track->P() <= 0.47){
     
     if(nsigmakaon < -0.5 && nsigmakaon > 2.0 && nsigmapion > 2.0 && nsigmaproton > 2.0){
     return 2;

   }
   }

   //fill nsigma kaons also
   
   if (track->P() < 0.39 || track->P() >0.47){
     
   if( kaonsym   <= 2.0  &&  nsigmapion > 2.0 && nsigmaproton > 2.0){
      fnsigmakaon->Fill((track->P()),nsigmakaon); 
     return 2;
   }
   }
   
   if( nsigmapion   <= 2.0 ){
     fnsigmapion->Fill((track->P()),nsigmapion);
     return 1;
   }

  
   if( nsigmaproton   <= 2.0 ){
     fnsigmaproton->Fill((track->P()),nsigmaproton);	
   
     return 3;
   }

   
   //why not filling nsigmakaon
   return 10;

}


   Int_t AliAnalysisTaskKaonIsospinFluct::ProcessHybridPro(AliAODTrack *track){

   Double_t nsigmaelectron = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kElectron)) ;
   
   Double_t nsigmakaon = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon)) ;

   Double_t nsigmaproton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton)) ;

   Double_t nsigmapion = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion)) ;


   Double_t nsigmaprotonTOF=999.,nsigmakaonTOF=999.,nsigmapionTOF=999.;
 
   
   
   Bool_t fHasTOFPID;
  
  
    
  if(track->GetStatus() & AliVTrack::kTOFpid){ fHasTOFPID=kTRUE; }
  
  else{ fHasTOFPID=kFALSE; }
  
  
  if (fHasTOFPID){
  
    nsigmapionTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion) ;
    nsigmakaonTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon) ;
    nsigmaprotonTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton) ;
    
   
    if( ( TMath::Abs(nsigmapionTOF)   <= 2.0 ) && ( TMath::Abs(nsigmapion) <= 2.0) ) 
      {
	fnsigmapion->Fill((track->P()),nsigmapion);
	return 1;
      }    
    
    if( ( TMath::Abs(nsigmakaonTOF)   <= 2.0 ) && ( TMath::Abs(nsigmakaon) <= 2.0) ) 
      {
	fnsigmakaon->Fill((track->P()),nsigmakaon);
	  return 2;
	}    

    if( nsigmaproton   <= 2.0 ){
      fnsigmaproton->Fill((track->P()),nsigmaproton);
      return 3;
    }
    
    
  }

  else{

    if( nsigmaproton   <= 2.0 ){
      fnsigmaproton->Fill((track->P()),nsigmaproton);
      return 3;
    }
  }
   

  return 10;

   }

   Int_t AliAnalysisTaskKaonIsospinFluct::ProcessHybrid(AliAODTrack *track){
     
   Double_t nsigmaelectron = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kElectron)) ;

   Double_t nsigmakaon = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon)) ;

   Double_t nsigmaproton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton)) ;

   Double_t nsigmapion = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion)) ;


   Double_t nsigmaprotonTOF=999.,nsigmakaonTOF=999.,nsigmapionTOF=999.;
 
   Bool_t fHasTOFPID;
  
  
    
   if(track->GetStatus() & AliVTrack::kTOFpid){ fHasTOFPID=kTRUE; }
   
   else{ fHasTOFPID=kFALSE; }
  
  
   if (fHasTOFPID){
  
    nsigmapionTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kPion) ;
    nsigmakaonTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kKaon) ;
    nsigmaprotonTOF = fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton) ;
    
      
    if( ( TMath::Abs(nsigmapionTOF)   <= 2.0 ) && ( TMath::Abs(nsigmapion) <= 2.0) ) 
	{
          fnsigmapion->Fill((track->P()),nsigmapion);
	  return 1;
	}    

    if( ( TMath::Abs(nsigmakaonTOF)   <= 2.0 ) && ( TMath::Abs(nsigmakaon) <= 2.0) ) 
	{
          fnsigmakaon->Fill((track->P()),nsigmakaon);
	  return 2;
	}    

    if( ( TMath::Abs(nsigmaprotonTOF)   <= 2.0 ) && ( TMath::Abs(nsigmaproton) <= 2.0) ) 
      {
	fnsigmaproton->Fill((track->P()),nsigmaproton);
	return 3;
      }    
    
  }

   
   return 10;

   }

   


