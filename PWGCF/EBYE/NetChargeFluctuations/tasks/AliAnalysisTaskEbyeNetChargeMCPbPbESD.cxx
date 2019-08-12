#include "Riostream.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2.h"
#include "TH1D.h"
#include "TH1I.h"
#include "TROOT.h"
#include "TMath.h"
#include "TF1.h"
#include "THn.h"
#include "TTree.h"
#include "TList.h"
#include "TObjArray.h"
#include "TTreeStream.h"
#include "TDatabasePDG.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliESDUtils.h"
#include "AliVEvent.h"
#include "AliInputEventHandler.h"
#include "AliESDInputHandler.h"
#include "AliHeader.h"
#include "AliMultSelection.h"
#include "AliMultSelectionBase.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCVertex.h"
#include "AliCentrality.h"
#include "AliPhysicsSelection.h"
#include "AliStack.h"
#include "AliLog.h"
#include "AliGenEventHeader.h"
#include "AliAnalysisUtils.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliVEventHandler.h"

#include "AliGenHijingEventHeader.h"
#include <AliStack.h>

#include "AliAnalysisTaskEbyeNetChargeMCPbPbESD.h"

class AliAnalysisTaskEbyeNetChargeMCPbPbESD;

#include "TFile.h"
#include <iostream>
#include <fstream>

using namespace std;

ClassImp(AliAnalysisTaskEbyeNetChargeMCPbPbESD)

//================================================
//--------------Constructor-----------------------
//================================================

AliAnalysisTaskEbyeNetChargeMCPbPbESD::AliAnalysisTaskEbyeNetChargeMCPbPbESD()
: AliAnalysisTaskSE(),
fOutputList(0),
fTreeSRedirector(0x0),
fEtaDown(-0.8),
fEtaUp(0.8),
fHistCentralityMultSelection(0),
fHistCentralityImpPar(0),
fMCImpactParameter(0),
fCentrality(0),
fHistVertexStats(0),
fHistZVertexCent(0),
fEventStatistics(0),
fHistVx(0),
fHistVy(0),
fHistVz(0),
fRunNumber(-1),
fArrayMC(0),
hGenPt(0),
hGenPhi(0),
hGenEta(0),
fEventCuts(0),
EventNumber(0),
fEta(0),
fpT(0),
fPhi(0),
fhCent(0),
fEtaMCgen(0),
fpTMCgen(0),
fPhiMCgen(0),
fChargegen(0),
genPos(0),
genNeg(0),
genNch(0),
genMomentsCross(0),
genMomentsPos(0),
genMomentsNeg(0),
fTreeMCgen(0x0)
{
    // default constructor
    
}

//-------------------------------------------------

AliAnalysisTaskEbyeNetChargeMCPbPbESD::AliAnalysisTaskEbyeNetChargeMCPbPbESD(const char* name)
:AliAnalysisTaskSE(name),
fOutputList(0),
fTreeSRedirector(0x0),
fEtaDown(-0.8),
fEtaUp(0.8),
fHistCentralityMultSelection(0),
fHistCentralityImpPar(0),
fMCImpactParameter(0),
fCentrality(0),
fHistVertexStats(0),
fHistZVertexCent(0),
fEventStatistics(0),
fHistVx(0),
fHistVy(0),
fHistVz(0),
fRunNumber(-1),
fArrayMC(0),
hGenPt(0),
hGenPhi(0),
hGenEta(0),
fEventCuts(0),
EventNumber(0),
fEta(0),
fpT(0),
fPhi(0),
fhCent(0),
fEtaMCgen(0),
fpTMCgen(0),
fPhiMCgen(0),
fChargegen(0),
genPos(0),
genNeg(0),
genNch(0),
genMomentsCross(0),
genMomentsPos(0),
genMomentsNeg(0),
fTreeMCgen(0x0)
{
    
    cout << "============================================================" << endl;
    cout << "=================End of Constructor====================="  << endl;
    cout << "============================================================" << endl;
    
    //constructor
    DefineOutput(1, TList::Class());
    DefineOutput(2, TTree::Class());

    
}
//==========================================================
//--------------Destructor----------------------------------
//==========================================================
AliAnalysisTaskEbyeNetChargeMCPbPbESD::~AliAnalysisTaskEbyeNetChargeMCPbPbESD()
{
    //destructor
    if(fOutputList) 			        delete fOutputList;
    if (fhCent)               		    delete fhCent;
    if (fHistCentralityMultSelection)   delete fHistCentralityMultSelection;
    if (fEventStatistics)               delete fEventStatistics;
    if (fHistCentralityImpPar)delete fHistCentralityImpPar;

}
//============================================================
//--------------UserCreateOutputObjects-----------------------
//============================================================
void AliAnalysisTaskEbyeNetChargeMCPbPbESD::UserCreateOutputObjects()
{
    
    fOutputList = new TList();
    fOutputList->SetOwner(kTRUE);
    
    // Event statistics
    fEventStatistics = new TH1I("fEventStatistics","",10,0,10);
    fOutputList->Add(fEventStatistics);
    
    hGenPt = new TH1F("hGenPt","generated p_{T};p_{T} (GeV/c);",200,0.0,10.0);
    hGenPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
    hGenPt->GetYaxis()->SetTitle("Counts");
    fOutputList->Add(hGenPt);
    hGenPhi = new TH1F("hGenPhi","generated #varphi;#varphi;",160,0,TMath::TwoPi());
    hGenPhi->GetXaxis()->SetTitle("#Phi ");
    hGenPhi->GetYaxis()->SetTitle("Counts");
    fOutputList->Add(hGenPhi);
    hGenEta = new TH1F("hGenEta","generated #eta;#eta;",40, -2.0, 2.0);
    hGenEta->GetXaxis()->SetTitle("#eta ");
    hGenEta->GetYaxis()->SetTitle("Counts");
    fOutputList->Add(hGenEta);
    
    //Vertex distributions
    fHistVx = new TH1D("fHistVx","Primary vertex distribution - x coordinate;V_{x} (cm)",400,-20, 20);
    fOutputList->Add(fHistVx);
    fHistVy = new TH1D("fHistVy","Primary vertex distribution - y coordinate;V_{y} (cm)",400,-20, 20);
    fOutputList->Add(fHistVy);
    fHistVz = new TH1D("fHistVz","Primary vertex distribution - z coordinate;V_{z} (cm)",400,-20,20);
    fOutputList->Add(fHistVz);
    
    fHistCentralityMultSelection = new TH1D("fHistCentralityMultSelection","Centrality Percentile ;Centrality;Entries",100,0.,100.);  //new sk
    fHistCentralityMultSelection->GetXaxis()->SetTitle("Centrality (%)");
    fOutputList->Add(fHistCentralityMultSelection);
    
    fHistCentralityImpPar  = new TH1F("fCentralityImpPar",     "control histogram for centrality imppar"    , 100, 0., 100.);
    fHistCentralityMultSelection->GetXaxis()->SetTitle("CentralityImpPar (%)");
    fOutputList->Add(fHistCentralityImpPar);

    
    fHistZVertexCent = new TH2F("fHistZVertexCent"," Vz of primary vertex Vs Centrality; V_{Z} {cm} ; Centrality {%}", 60, -15, 15,100,0,100);
    fOutputList->Add(fHistZVertexCent);
  

    // ************************************************************************
    
    // Tree for pt eta and centrality checks
    fTreeSRedirector    = new TTreeSRedirector();
    fTreeMCgen          = ((*fTreeSRedirector)<<"MCgen").GetTree();
    
    cout << "============================================================" << endl;
    cout << "=================End of User Create Object=================="  << endl;
    cout << "============================================================" << endl;
    
    //Send data to the container
    PostData(1, fOutputList);
    PostData(2, fTreeMCgen);
    
}
//============================================================
//----------------------UserExec------------------------------
//============================================================
void AliAnalysisTaskEbyeNetChargeMCPbPbESD::UserExec(Option_t *){
    
    
  //  cout << "============================================================" << endl;
  //  cout << "================ User Exec Object=================="  << endl;
 //   cout << "============================================================" << endl;


    AliAnalysisManager* anMan = AliAnalysisManager::GetAnalysisManager();
    AliESDInputHandler *handler = (AliESDInputHandler*)anMan->GetInputEventHandler();

    AliMCEventHandler* eventHandler = 0x0;
                                                                                                                                                           
    // ========================== MC EVENT INF O====================================>>>                                                                      
    eventHandler = (AliMCEventHandler*)anMan->GetMCtruthEventHandler();
    if (!eventHandler) { printf("ERROR: Could not retrieve MC event handler\n"); return; }

    AliMCEvent *fMCEvent = 0x0;
    fMCEvent = eventHandler->MCEvent();
    if (!fMCEvent) { printf("ERROR: Could not retrieve MC event\n"); return; }
 
//========= Generator Info ==========================                                                                                                                                            
    if(fMCEvent) {
      AliGenEventHeader* genHeader = fMCEvent->GenEventHeader();
      if(!genHeader){
        printf("  Event generator header not available!!!\n");
        return;
      }
      fMCImpactParameter = ((AliGenHijingEventHeader*) genHeader)->ImpactParameter();
    }
                  
    AliStack *stack = fMCEvent->Stack();

    //===========================Centrality calculation =====================
 
   // cout<< " hello here  Centrality selection "<< endl;   

    Double_t impParArr[10] = {0.0, 3.72, 5.23, 7.31, 8.88, 10.20, 11.38, 12.47, 13.50, 14.5};
    
    if (fMCImpactParameter>=impParArr[0] && fMCImpactParameter<impParArr[1]) fCentrality=2.5;
    if (fMCImpactParameter>=impParArr[1] && fMCImpactParameter<impParArr[2]) fCentrality=7.5;
    if (fMCImpactParameter>=impParArr[2] && fMCImpactParameter<impParArr[3]) fCentrality=15.;
    if (fMCImpactParameter>=impParArr[3] && fMCImpactParameter<impParArr[4]) fCentrality=25.;
    if (fMCImpactParameter>=impParArr[4] && fMCImpactParameter<impParArr[5]) fCentrality=35.;
    if (fMCImpactParameter>=impParArr[5] && fMCImpactParameter<impParArr[6]) fCentrality=45.;
    if (fMCImpactParameter>=impParArr[6] && fMCImpactParameter<impParArr[7]) fCentrality=55.;
    if (fMCImpactParameter>=impParArr[7] && fMCImpactParameter<impParArr[8]) fCentrality=65.;
    if (fMCImpactParameter>=impParArr[8] && fMCImpactParameter<impParArr[9]) fCentrality=75.;
    if (fMCImpactParameter<impParArr[0]  || fMCImpactParameter>impParArr[9]) fCentrality=-10.;

  //  cout << " Centrality = " << fCentrality << " ------ Impact Parameter = " << fMCImpactParameter << endl;

    fHistCentralityImpPar->Fill(fCentrality);

    // =========Different Acceptance cuts for Eta momentum and Centrality=============
    Double_t momDownArray[1]    = {0.2};
    Double_t momUpArray[1]      = {5.0};
    Double_t etaDownArray[8]    = {-0.1, -0.2, -0.3, -0.4, -0.5, -0.6, -0.7, -0.8};
    Double_t etaUpArray[8]      = {0.1,  0.2,  0.3,  0.4,  0.5,  0.6,  0.7,  0.8};
    Double_t centDownArray[9]   = {0., 5.,  10., 20., 30., 40., 50., 60., 70.};
    Double_t centUpArray[9]     = {5., 10., 20., 30., 40., 50., 60., 70., 80.};
    
    for (Int_t imom=0; imom<1; imom++){
        for (Int_t ieta=0; ieta<8; ieta++){
            for (Int_t icent=0; icent<9; icent++){

             Double_t centBin = (centDownArray[icent]+centUpArray[icent])/2.;   
             Double_t etaBin  = (TMath::Abs(etaDownArray[ieta])+etaUpArray[ieta]);

              genPos=0.,genNeg=0.,genNch=0.;
 	      Int_t subsample1     = 0 ;
              Int_t trCountMCgen   = 0;
              subsample1 = EventNumber%30;

            //=========================== track loop ==============================
         
	        for (Int_t iTracks = 0; iTracks < fMCEvent->GetNumberOfTracks(); iTracks++) {
		
		pdgMomPhysicalPrimgen = 0, pdgMomPrimgen = 0,pdgPhysicalPrimgen = 0, pdgPrimgen = 0;

		AliMCParticle* mctrack = dynamic_cast<AliMCParticle*>(fMCEvent->GetTrack(iTracks));
		if(!mctrack)
		  continue;
		
		TParticle * part = mctrack->Particle();
		if( !part )
		  continue;                    
	
		if( !stack->IsPhysicalPrimary(iTracks) )
		  continue;

		fpTMCgen     = mctrack->Pt();
		fPhiMCgen    = mctrack->Phi();
		fChargegen   = mctrack->Charge();
		fEtaMCgen    = mctrack->Eta();

		//cout<<"pt  -->"<< fpTMCgen<< "   "<<"phi  -->"<<fPhiMCgen<< "   "<<"eta  -->"<< fChargegen<< endl;
		
		if (fChargegen == 0) continue;
		if (fpTMCgen < 0.2) continue;

		   hGenPt->Fill(fpT);
                    hGenPhi->Fill(fPhi);
                    hGenEta->Fill(fEta);

                 // MC eta cut
                    if ((fEtaMCgen < etaDownArray[ieta]) || (fEtaMCgen > etaUpArray[ieta])) continue;
             	
    //===================apply cuts on pt eta and centrality=====================
                if ((fpTMCgen>=momDownArray[imom])
                        &&(fpTMCgen<momUpArray[imom])
 			&&(fCentrality>=centDownArray[icent])
                        &&(fCentrality<centUpArray[icent])){

             // calculate first moments
		if (fChargegen < 0 || fChargegen > 0) genNch++;
		if(fChargegen > 0) genPos++;
		if(fChargegen < 0) genNeg++;
                 trCountMCgen++;

} // end loop of apply cuts
		
	      } //============end of track loop====================
	      
  		// calculate second moments
                genMomentsCross  = genPos * genNeg;
                genMomentsPos    = genPos * genPos;
                genMomentsNeg    = genNeg * genNeg;
	  
		  
                // Fill the pt eta and centrality check in the tree
                if ( trCountMCgen>0 ){
                    if(!fTreeSRedirector) return;
                    (*fTreeSRedirector)<<"MCgen"<<
                    "pos="                 << genPos                    <<    // MC gen positive particles
                    "neg="                 << genNeg                    <<    // MC gen negative particles
                    "isample1="      	   << subsample1                <<    // sample id for subsample method
                    "momentscross="        << genMomentsCross           <<    // MC gen second moments of pos and neg
                    "momentspos="          << genMomentsPos             <<    // MC gen second moments of pos
                    "momentsneg="          << genMomentsNeg             <<    // MC gen second moments of neg
                    "etaDown="             << etaDownArray[ieta]        <<    // lower edge of eta bin
                    "etaUp="               << etaUpArray[ieta]          <<    // upper edge of eta bin
                    "centDown="            << centDownArray[icent]      <<    // lower edge of cent bin
                    "centUp="              << centUpArray[icent]        <<    // upper edge of cent bin
                    "momDown="             << momDownArray[imom]        <<    // lower edge of mom bin
                    "momUp="               << momUpArray[imom]          <<    // upper edge of mom bin
                    "ieta="                << ieta                      <<    // eta index loop
                    "icent="               << icent                     <<    // cent index loop
                    "imom="                << imom                      <<    // mom index loop
                    "centBin="             << centBin                   <<    // centrality bin
                    "etabin="              << etaBin                    <<    // eta bin
                    "cent="                << fCentrality               <<
                    "impPar="              << fMCImpactParameter        <<
                    "\n";
                    
                } //===========tree filling========
   
       } //===========end of momentum loop
        } //===========end of eta loop
 } //===========end of centrality loop
    EventNumber++;
    
    PostData(1, fOutputList);
    PostData(2, fTreeMCgen);
    
}

//------------------------------------------------------------------------------------------
void AliAnalysisTaskEbyeNetChargeMCPbPbESD::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
    Info("AliAnalysisEbyeNetChargeFluctuations"," Task Successfully finished");
    AliInfo(Form("Found  %d MC events",EventNumber));
    
    
}
//____________________________________________________________________________________________
