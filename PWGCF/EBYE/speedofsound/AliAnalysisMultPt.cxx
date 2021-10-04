#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <math.h>
#include "TChain.h"
#include "TFile.h"
#include "TKey.h"
#include "TObject.h"
#include "TObjString.h"
#include "TList.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TBits.h"
#include "TRefArray.h"
#include "TStopwatch.h"
#include "TAxis.h"
#include "TRandom.h"
#include "TList.h"
#include "TComplex.h"
#include "TMath.h"
#include <cmath>
#include <cstdlib>
#include "AliAnalysisTask.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAODpidUtil.h"
#include "AliPID.h"
#include "AliAODPid.h"
#include "AliCentrality.h"
#include "AliPIDResponse.h"
#include "AliAODVZERO.h"
#include "AliCollisionGeometry.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliGenEventHeader.h"
#include "AliHeader.h"
#include "AliMultSelection.h"
#include "AliAnalysisUtils.h"
#include "AliAnalysisMultPt.h"

//#include "AliAODVertex.h"




class AliAnalysisMultPt;    // your analysis class
using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisMultPt) // classimp: necessary for root
//________________________________________________________________________

AliAnalysisMultPt::AliAnalysisMultPt() : AliAnalysisTaskSE(),
fAOD(0),fMCEvent(0),
fOutputList(0),
fHistMult(0),
fHistPt(0),
fHistMultPt(0),
fHistRatio(0),
fHistProj(0),
fHistMultMC(0),
fHistPtMC(0),
fHistMultPtMC(0),
fHistProjMC(0),
fIsMC(0),
mixer(0)
{
    mMode = 0;

     
}
//________________________________________________________________________
AliAnalysisMultPt::AliAnalysisMultPt(const char* name) :
AliAnalysisTaskSE(name),
fAOD(0),fMCEvent(0),
fOutputList(0),
fHistMult(0),
fHistPt(0),
fHistMultPt(0),
fHistProj(0),
fHistMultMC(0),
fHistPtMC(0),
fHistMultPtMC(0),
fHistProjMC(0),
fHistRatio(0),
fIsMC(0),
mixer(0)

{


    //DefineInput(0, TChain::Class());
    mMode = 0;

    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
    //DefineOutput(1, TTree::Class());
    
}
//________________________________________________________________________
AliAnalysisMultPt::~AliAnalysisMultPt()
{
    // destructor
     if(fOutputList) {
        delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
    }
}

//________________________________________________________________________

void AliAnalysisMultPt::ProcessMCParticles()
{
    // process MC particles
    TClonesArray* AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if (AODMCTrackArray == NULL) return;

    // Loop over all primary MC particle
    for(Long_t i = 0; i < AODMCTrackArray->GetEntriesFast(); i++) {

      AliAODMCParticle* particle = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(i));
      if (!particle) continue;
      //cout << "PDG CODE = " << particle->GetPdgCode() << endl;
    }
}

void AliAnalysisMultPt::UserCreateOutputObjects()
    {
        // create a new TList that OWNS its objects
        fOutputList = new TList();
        //fOutputList = new TTree();
        //fOutputList->SetOwner(true);
        fOutputList->SetOwner(kTRUE);
        

        //ostad
        // Create histograms
        // Called once
          //Init();// Initialize my settings
         // OpenFile(1,"RECREATE");
            //double PtBins[20] = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.2,1.4,1.6,1.8,2,2.5,3,3.5,4,4.5,5};

        
        // create our histo and add it to the list
        fHistPt = new TH1D("fHistPt", "fHistPt", 50, 0, 5);
        fOutputList->Add(fHistPt);
        
        fHistMult = new TH1D("fHistMult", "fHistMult", 350, 0, 3500);
        fOutputList->Add(fHistMult);
        
        fHistMultPt = new TH2D("fHistMultPt", "fHistMultPt",50, 0, 5, 350, 0, 3500);  //?????
        fOutputList->Add(fHistMultPt);

        fHistFit = new TH1D("fHistFit", "fHistFit", 350, 0, 3500);
        fOutputList->Add(fHistFit);
        
        fHistProj = new TH1D("fHistProj", "fHistProj", 50, 0, 5);
        fOutputList->Add(fHistProj);
        

        fHistPtMC = new TH1D("fHistPtMC", "fHistPtMC", 50, 0, 5);
        fOutputList->Add(fHistPtMC);
        
        fHistMultMC = new TH1D("fHistMultMC", "fHistMultMC", 350, 0, 3500);
        fOutputList->Add(fHistMultMC);
        
        fHistMultPtMC = new TH2D("fHistMultPtMC", "fHistMultPtMC",50, 0, 5, 350, 0, 3500);  //?????
        fOutputList->Add(fHistMultPtMC);

            
        //fHistRatio = new TH2D("fHistRatio", "fHistRatio", 50, 0, 5, 350, 0, 3500);  //?????
        //fOutputList->Add(fHistRatio);
        
        fHistRatio = new TH1D("fHistRatio", "fHistRatio", 50, 0,5);  //?????
        //TH1D *fHistRatio;
        fOutputList->Add(fHistRatio);
        
            
        fHistFitMC = new TH1D("fHistFitMC", "fHistFitMC", 350, 0, 3500);
        fOutputList->Add(fHistFitMC);
        
        fHistProjMC = new TH1D("fHistProjMC", "fHistProjMC", 50, 0, 5);
        fOutputList->Add(fHistProjMC);
 
        
        // add the list to our output file
        PostData(1,fOutputList);
        //PostData(1, mixer);
    }
//________________________________________________________________________

void AliAnalysisMultPt::UserExec(Option_t *)
{
    fMCEvent = MCEvent();
    // get an event from the analysis manager:
    fAOD = dynamic_cast<AliAODEvent *> (InputEvent());
    // check if there actually is an event:
    if(!fAOD) return;
    // ::Fatal("AliAnalysisMultPt::UserExec", "No AOD event found, check the event handler.");
    // let's loop over the tracks and fill our histogram,
    // first we get the number of tracks:
    Int_t iTracks(fAOD->GetNumberOfTracks());
    // and then loop over them:
    Double_t Mult=0.0; //For charged particle
    Double_t Mult1=0.0; //For positive
    Double_t Mult2=0.0; //For negative
    
    
    
    
//TCanvas *c1 = new TCanvas("c1","c1",800,600) ;
 //   c1->Divide(2,2);
    
    // loop over all these tracks:
    for(Int_t i(0); i < iTracks; i++) {
         // get a track (type AliAODTrack) from the event:
         AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
         if(!track) continue;                            // if we failed, skip this track
         //if(!track->TestFilterBit(768)) continue;
         //if(!track->TestFilterBit(96)) continue;
              Short_t Charge = track->Charge();
              if(TMath :: Abs(Charge)>0) {
                  if(track->Pt() < 0.2|| track->Pt() > 5) continue;
                  if(fabs(track->Eta()) > 0.8) continue;   //eta cut
                  if(track->TestFilterBit(96)){
                      Mult++;
                  //SUM OVER PT??????????
              }
              }
         }
    
    
    
    for(Int_t i(0); i < iTracks; i++) {
         // get a track (type AliAODTrack) from the event:
         AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
         // if we failed, skip this track:
         if(!track) continue;
              Short_t Charge = track->Charge();
        
              //Count the charge particles:
              if(TMath :: Abs(Charge)>0) {
                  if(track->Pt() < 0.2|| track->Pt() > 5) continue;
                  if(fabs(track->Eta()) > 0.8) continue;   //eta cut
                  if(track->TestFilterBit(96)){
                     fHistMultPt->Fill(track->Pt(), Mult);
              }
          }
         }
    
    TH1D *fHistProj = fHistMultPt->ProjectionX("fHistProj", 0, 3500); // where firstYbin = 0 and lastYbin = 5

//    c1->cd(1);
    //fHistProj->Draw();
    
    fHistMultPt->GetXaxis()->SetTitle("pT");
    fHistMultPt->GetYaxis()->SetTitle("Multiplicity");
    
/////////////////////////////////////////// Generated part:
    
     if(fMCEvent) ProcessMCParticles();

     TClonesArray *stack = 0;
     TList *lst = fAOD->GetList();
     stack = (TClonesArray*)lst->FindObject(AliAODMCParticle::StdBranchName());
     int nMCTracks;
     if (!stack) nMCTracks = 0;
     else nMCTracks = stack->GetEntries();
     
     
     
     Double_t mNMC  = 0.0;
     //Double_t  mSpTMC = 0.0;
     
     
     for (Int_t i(0); i < nMCTracks; i++) {
         
         AliAODMCParticle *particle=(AliAODMCParticle*)stack->UncheckedAt(i);
         if (!particle) continue;
         if(!particle->IsPrimary()) continue;
         if(!particle->IsPhysicalPrimary()) continue;
         Short_t Charge = particle->Charge();
         
             if(TMath :: Abs(Charge)>0) {
                 if(particle->Pt() < 0.2|| particle->Pt() > 5) continue;
                 if(fabs(particle->Eta()) > 0.8 ) continue;
                 if((fabs(particle->GetPdgCode())==211)||(fabs(particle->GetPdgCode())==2212)||(fabs(particle->GetPdgCode())==321)){
                 //mSpTMC+=particle->Pt();
                 mNMC++;
         }
       }
     }
    


 

     
     for (Int_t i(0); i < nMCTracks; i++) {
         
         AliAODMCParticle *particle=(AliAODMCParticle*)stack->UncheckedAt(i);
         if (!particle) continue;
         if(!particle->IsPrimary()) continue;
         if(!particle->IsPhysicalPrimary()) continue;
         Short_t Charge = particle->Charge();
         
             if(TMath :: Abs(Charge)>0) {
                 if(particle->Pt() < 0.2|| particle->Pt() > 5) continue;
                 if(fabs(particle->Eta()) > 0.8 ) continue;
                 if((fabs(particle->GetPdgCode())==211)||(fabs(particle->GetPdgCode())==2212)||(fabs(particle->GetPdgCode())==321)){
                 fHistMultPtMC->Fill(particle->Pt(), mNMC);
         }
        }
     }
            //fHistMultPtMC->Draw();
            TH1D *fHistProjMC = fHistMultPtMC->ProjectionX("fHistProjMC", 0, 3500); // where firstYbin = 0 and lastYbin = 5
            //fHistProjMC->Draw();
            //fHistProjMC-> SetLineColor(kRed);
            //fHistProj->Divide(fHistProjMC);
            //fHistRatio->Draw
    //c1->cd(2);
   // fHistProjMC->Draw();
            
            //fHistRatio->Fill(fHistProj/fHistProjMC);
            fHistMultPtMC->GetXaxis()->SetTitle("pT-MC");
            fHistMultPtMC->GetYaxis()->SetTitle("Multiplicity-MC");
       
            //TH2D *fHistRatio->(fHistMultPtMC->Divide(fHistMultPt));
           // fHistRatio->Draw();
           
           //TH1D *hratio = new TH1D(*fHistProj);
          //hratio->Divide(fHistProjMC);
          // hratio->Draw();

     //fHistProj->Draw();
  

         fHistRatio = (TH1D*)fHistProj->Clone();
         fHistRatio->Divide(fHistProjMC);
         fHistRatio->Draw();
     
    /*
            TH1D *hratio = new TH1D(*fHistProj);
            hratio->Divide(fHistProjMC);
            hratio->Draw();
     */
    
    //hratio->SetNameTitle("hratio", "h1 / h2");
    
        //cout<<"Charged particle = " << Mult<< endl;
       // fHistMult->Fill(Mult);
        //fHistMult->Draw();
        //fHistFit->Fill(); ??
        //fHistFit->Draw();
       //cout<<"Positive = " << Mult1 << endl;
       //cout<<"Negative = " << Mult2 << endl;
      // and save the data gathered in this iteration

       
   
   PostData(1, fOutputList);
}

//________________________________________________________________________
void AliAnalysisMultPt::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
    
