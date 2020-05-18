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

/* AliAnaysisTaskMyTask*/
// c++ headers
#include <iostream>
#include <string.h>

// root headers
#include "TH1I.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TObjString.h"
#include "TFile.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"

// aliroot headers
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliAODVZERO.h"
#include "AliAODZDC.h"
#include "AliESDVZERO.h"
#include "AliESDZDC.h"
#include "AliPIDResponse.h"
#include "AliAODTrack.h"
#include "AliAODPid.h"
#include "AliAODVertex.h"
#include "AliESDVertex.h"
#include "AliMultiplicity.h"
#include "AliESDtrack.h"
#include "AliESDMuonTrack.h"
#include "AliAODMCParticle.h"
#include "AliMCParticle.h"
#include "AliCentrality.h"
#include "AliKFVertex.h"
#include "AliExternalTrackParam.h"
//#include "AliTriggerAnalysis.h"
#include "AliAODMCHeader.h"
#include "TChain.h"
#include "TH1F.h"
#include "TList.h"
#include "AliAnalysisTask.h"
#include "AliAODInputHandler.h"
#include "Polarization.h"
#include "AliPIDResponse.h"
#include "TMath.h"

class Polarization;    
using namespace std;            

ClassImp(Polarization) // classimp: necessary for root


Polarization::Polarization() : AliAnalysisTaskSE(), 
    fAOD(0), fOutputList(0), fHistPt(0) , fHistP_TPC(0), fTree(0),fPt(0),fM(0),fPt0(0),fPt1(0),fPIDResponse(0),fKaonSigma1(0),fKaonSigma0(0),fPiSigma1(0),fPiSigma0(0), fTPCcluster1(0), fEta1(0),fTPCcluster2(0), fEta2(0),fHistEtaCut(0),fHistITS(0), fHistChargeCut(0),fHistPtCut(0),fHistPairCut(0),fHistPIDCut(0),
    fHistFilterbitCut(0),fHistTriggerCut(0),fHistCounter(0),fDCAxy2(0),fDCAz2(0),fDCAxy1(0),fDCAz1(0) , fTriggerClass(0),fdEdX(0),fPd(0),fPp(0),fPtd(0),fZDCdata(0),fZNAenergy(0), fZNCenergy(0),fZDCAtime(0),fZDCCtime(0),fMuSigma0(0),fMuSigma1(0),fRabs1(0),fRabs2(0),fTheta(0),fHelicityTheta(0),fCollinTheta(0),fPhi(0),fHelicityPhi(0),fCollinPhi(0), daughter1(0.,0.,0.,0.), daughter2(0.,0.,0.,0.) , parent(0.,0.,0.,0.)
  {
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
  }
//_____________________________________________________________________________
Polarization::Polarization(const char* name) : AliAnalysisTaskSE(name),
    fAOD(0), fOutputList(0), fHistPt(0),fHistP_TPC(0), fTree(0),fPt(0),fM(0),fPt0(0),fPt1(0) ,fPIDResponse(0),fKaonSigma0(0),fKaonSigma1(0),fPiSigma0(0),
    fPiSigma1(0), fTPCcluster1(0), fEta1(0),fTPCcluster2(0), fEta2(0),fHistEtaCut(0),fHistITS(0), fHistChargeCut(0),fHistPtCut(0),fHistPairCut(0),fHistPIDCut(0),
    fHistFilterbitCut(0),fHistTriggerCut(0),fHistCounter(0) ,fDCAxy1(0),fDCAz1(0),fDCAxy2(0),fDCAz2(0),fTriggerClass(0),fdEdX(0),fPd(0),fPp(0),fPtd(0),fZDCdata(0), fZNAenergy(0), fZNCenergy(0),fZDCAtime(0),fZDCCtime(0),fMuSigma0(0),fMuSigma1(0),fRabs1(0), fRabs2(0),fTheta(0),fHelicityTheta(0),fCollinTheta(0),fPhi(0),fHelicityPhi(0),fCollinPhi(0), daughter1(0.,0.,0.,0.), daughter2(0.,0.,0.,0.) , parent(0.,0.,0.,0.)
  {
    
    
    // constructor
    DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
                                        // this chain is created by the analysis manager, so no need to worry about it, 
                                        // it does its work automatically
    DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms 
                                        // you can add more output objects by calling DefineOutput(2, classname::Class())
                                        // if you add more output objects, make sure to call PostData for all of them, and to
                                        // make changes to your AddTask macro!
   DefineOutput(2, TTree::Class());
  }
//_____________________________________________________________________________
Polarization::~Polarization()
 {
  // destructor
  if(fOutputList) 
    {
  
    delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
    }
 }





//_____________________________________________________________________________
void Polarization::UserCreateOutputObjects()
  {
    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    if (man) 
      {
        AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
        if (inputHandler)   fPIDResponse = inputHandler->GetPIDResponse();
        if (!fPIDResponse) {
        cout<<"noPIDresopne"<<endl;
        return;
     }
   
  }
    
    
    fOutputList = new TList();          
    fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them
                                          // if requested (dont worry about this now)
  
      // example of a histogram
    fHistPt = new TH1F("fHistPt", "fHistPt", 100, 0, 10);       
    fOutputList->Add(fHistPt);  
    fHistCounter = new TH1I("fHistCounter","Histogram of Counter",100,0,100);  
    fOutputList->Add(fHistCounter);                              
    fHistP_TPC = new TH2F ("fHistP_TPC","fHistP_TPC",100,0,1,100,0,0.001);       // your histogram in the output file, add it to the list!
    fOutputList->Add(fHistP_TPC);
     
   //Defining the tree branches to fill therequired information.
   //
   
   
   
    fTree = new TTree("result","Momentum and TPC signal");
    
    
    fTree->Branch("daughter1","TLorentzVector",&daughter1);
    fTree->Branch("daughter2","TLorentzVector",&daughter2);
    fTree->Branch("parent","TLorentzVector",&parent);
    
    
    fTree->Branch("fPt", &fPt, "fPt/F");
    fTree->Branch("fTheta", &fTheta, "fTheta/F");
    fTree->Branch("fHelicityTheta", &fHelicityTheta, "fHelicityTheta/F");
    fTree->Branch("fCollinTheta", &fCollinTheta, "fCollinTheta/F");
    
    fTree->Branch("fPhi", &fPhi, "fPhi/F");
    fTree->Branch("fHelicityPhi", &fHelicityPhi, "fHelicityPhi/F");
    fTree->Branch("fCollinPhi", &fCollinPhi, "fCollinPhi/F");
    
    fTree->Branch("fM", &fM, "fM/F");
    fTree->Branch("fPt0", &fPt0, "fPt0/F");
    fTree->Branch("fPt1", &fPt1, "fPt1/F"); 
  //  fTree->Branch("fKaonSigma0", &fKaonSigma0, "fKaonSigma0/F");
  //  fTree->Branch("fKaonSigma1", &fKaonSigma1, "fKaonSigma1/F");
  //  fTree->Branch("fPiSigma0", &fPiSigma0, "fPiSigma0/F");
  //  fTree->Branch("fPiSigma1", &fPiSigma1, "fPiSigma1/F");
  //  fTree->Branch("fMuSigma0", &fMuSigma0, "fMuSigma0/F");
  //  fTree->Branch("fMuSigma1", &fMuSigma1, "fMuSigma1/F");
  //  fTree->Branch("fTPCcluster1", &fTPCcluster1, "fTPCcluster1/F");
  //  fTree->Branch("fEta1", &fEta1, "fEta1/F");
  //  fTree->Branch("fTPCcluster2", &fTPCcluster2, "fTPCcluster2/F");
  //  fTree->Branch("fEta2", &fEta2, "fEta2/F");
    fTree->Branch("fDCAxy1", &fDCAxy1, "fDCAxy1/F");
    fTree->Branch("fDCAxy2", &fDCAxy2, "fDCAxy2/F");
    fTree->Branch("fDCAz1", &fDCAz1, "fDCAz1/F");
    fTree->Branch("fDCAz2", &fDCAz2, "fDCAz2/F");
    fTree->Branch("fTriggerClass", &fTriggerClass,"fTriggerClass/I");
  //  fTree->Branch("fPp", &fPp, "fPp/F");
  //   fTree->Branch("fPt0", &fPd, "fPd/F");
  //  fTree->Branch("fPtd", &fPtd, "fPtd/F");
   // fTree->Branch("fdEdX", &fdEdX, "fdEdX/F");
    fTree->Branch("fZNAenergy", &fZNAenergy, "fZNAenergy/F"); 
    fTree->Branch("fZNCenergy", &fZNCenergy, "fZNCenergy/F");
    fTree->Branch("fZDCAtime", &fZDCAtime, "fZDCAtime/F");
    fTree->Branch("fZDCCtime", &fZDCCtime, "fZDCCtime/F");
  //  fTree->Branch("fHistCounter", &fHistCounter, "fHistCounter/I");
    
    PostData(1, fOutputList);           
   
  
    PostData(2,fTree);
  }
//_____________________________________________________________________________





void Polarization::UserExec(Option_t *)
  {



 
     fAOD = dynamic_cast<AliAODEvent*>(InputEvent());    
     if(!fAOD) return;                                  
       
     fHistCounter->Fill(1);
       
     TString trigger = fAOD-> GetFiredTriggerClasses();
    // cout<<trigger<<endl;
  //   if (trigger.Contains("CCUP8")) return; 
    // cout<< "trigger classes are " << trigger << endl; return;  
     if (!trigger.Contains("CMUP")) cout<<"no cmup"<<endl; return;
     if (trigger.Contains("CMUP10")){
         fTriggerClass =10;
         }
     if (trigger.Contains("CMUP11")){
         fTriggerClass =11;
         } 
     if  (trigger.Contains("CMUP13")){
         fTriggerClass =13;
         }  
     if  (trigger.Contains("CMUP26")){
         fTriggerClass =26;
         }
     if  (trigger.Contains("CMUP6")){
         fTriggerClass =6;
         }  
           
     cout<<trigger<<endl;     
    
      AliAODZDC *fZDCdata = fAOD->GetZDCData();
      fZNAenergy = fZDCdata->GetZNATowerEnergy()[0];
      fZNCenergy = fZDCdata->GetZNCTowerEnergy()[0];
      fZDCAtime = fZDCdata->GetZNATime();
      fZDCCtime = fZDCdata->GetZNCTime();
  //ZDC cuts but not applying here, planning to apply in the final result...
  //if(trigger.Contains("CCUP4-B"))fHistZDCCuts->Fill(1);
  /*if(fZNAenergy < 8200 && fZNCenergy < 8200) fHistZDCCuts->Fill(2);
  if(fZNAenergy < 683 && fZNCenergy < 683) fHistZDCCuts->Fill(3);
  if(fZDCAtime == 0 && fZDCCtime == 0) fHistZDCCuts->Fill(4);*/
  
  
     fHistCounter->Fill(2);            
     Int_t iTracks(fAOD->GetNumberOfTracks()); // see how many tracks there are in the event
     Int_t PairCounter = 0;
     Bool_t GoodTracks = kFALSE;
     Bool_t GoodTracks2 = kFALSE;
     AliAODTrack* savetrack1;
     AliAODTrack* savetrack2;
     for(Int_t i(0); i < iTracks; i++) 
      { 
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));                // loop ove rall these tracks
        
        fRabs1 = track->GetRAtAbsorberEnd();
        
      //  if(!(track->TestFilterBit(1<<0))) continue;
        fHistCounter->Fill(2); 
        if(track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1)) {GoodTracks = kTRUE;}
        
        else  {GoodTracks = kFALSE;}
        if (GoodTracks == kFALSE) fHistCounter->Fill(3);//continue; 
         
        if (track->Pt()>1) continue;           
      
        if (track->Eta()>-2.5||track->Eta()<-4) continue; //fHistCounter->Fill(4);//continue;
        
        
        
        if (17.5>fRabs1||fRabs1>89.5) continue;//fHistCounter->Fill(5);//continue;
        
        
        //debugger to check if the code is working
        if(17.5>fRabs1) cout<<"rabs cut not working on track 2"<<endl;
        if(89.5<fRabs1) cout<<"rabs cut not working on track 2"<<endl;
        if(-2.5<track->Eta()) cout<<"eta cut not working"<<endl;
        if(-4>track->Eta()) cout<<"eta cut not working on track 2"<<endl;  
          
          
        
        
        
                  
          for (Int_t j=i+1 ; j<iTracks ; j++) {
              AliAODTrack* track2 = static_cast<AliAODTrack*>(fAOD->GetTrack(j));
              fRabs2 = track2->GetRAtAbsorberEnd();
              
             
            //  fHistCounter->Fill(5);
              if(track2->HasPointOnITSLayer(0) || track2->HasPointOnITSLayer(1)) {GoodTracks2 = kTRUE;}
              else { GoodTracks2 = kFALSE ;}
              if (GoodTracks2 == kFALSE)   fHistCounter-> Fill(6);// continue;
              if (track2->Pt()>1) continue; 
              if (track2->Eta()>-2.5||track2->Eta()<-4) continue;
              if (17.5>fRabs2||fRabs2>89.5) continue;
             
             //some debuggers  to make sure the cut is working as I am little concerned they might not work
              if(17.5>fRabs2) cout<<"rabs cut not working on track 2"<<endl;
              if(89.5<fRabs2) cout<<"rabs cut not working on track 2"<<endl;
              if(-2.5<track2->Eta()) cout<<"eta cut not working"<<endl;
              if(-4>track2->Eta()) cout<<"eta cut not working on track 2"<<endl;
              
           
              fHistCounter->Fill(8);
              PairCounter = PairCounter+1;  
              savetrack1 = track;
              savetrack2 = track2;
             
             // i am braking it to make so that code do not run if pair is more than 1..
             // it makes code run faster 
              if (PairCounter<1) break;     
              
            
            
            
              }
       //end of for loop j tracks*/
       
      }//end of for loop i tracks
      //Selecting tracks with only 1 pair of tracks
     
     
      if (PairCounter!=1)  return;
      
      fHistCounter->Fill(9);
     
     
     
      Double_t dca[2] = {0.0,0.0}, cov[3] = {0.0,0.0,0.0};
      AliAODTrack* track1_clone=(AliAODTrack*)savetrack1->Clone("track1_clone");
      AliAODVertex *fAODVertex = fAOD->GetPrimaryVertex();
      if(!track1_clone->PropagateToDCA(fAODVertex,fAOD->GetMagneticField(),300.,dca,cov)) { 
       dca[0] = -999 ;
       dca[1] = -999 ;
       }
      fDCAxy1= dca[0];
      fDCAz1=  dca[1];
      delete track1_clone;
      AliAODTrack* track2_clone=(AliAODTrack*)savetrack2->Clone("track2_clone");
       if(!track2_clone->PropagateToDCA(fAODVertex,fAOD->GetMagneticField(),300.,dca,cov)){ 
       dca [0] = -999 ;
       dca [1] = -999 ;
       }
      fDCAxy2= dca[0];      
      fDCAz2= dca[1];
      delete track2_clone;
      if(TMath::Abs(dca[1]) > 2) return;
      //Charge Cut 
     
      if (savetrack1->Charge() * savetrack2->Charge()>= 0) return;
      
      fHistCounter->Fill(10);   
      fPiSigma0=fPIDResponse->NumberOfSigmasTPC(savetrack1, AliPID::kPion);
      fPiSigma1=fPIDResponse->NumberOfSigmasTPC(savetrack2, AliPID::kPion);
      fKaonSigma0=fPIDResponse->NumberOfSigmasTPC(savetrack1, AliPID::kKaon);
      fKaonSigma1=fPIDResponse->NumberOfSigmasTPC(savetrack2, AliPID::kKaon);
      fMuSigma0=fPIDResponse->NumberOfSigmasTPC(savetrack1, AliPID::kMuon);
      fMuSigma1=fPIDResponse->NumberOfSigmasTPC(savetrack2, AliPID::kMuon);
      
  
  
      fHistCounter->Fill(11);
      TLorentzVector d1;
      TLorentzVector d2; 
      d1.SetPtEtaPhiM(savetrack1->Pt(),savetrack1->Eta(),savetrack1->Phi(),0.105);
      d2.SetPtEtaPhiM(savetrack2->Pt(),savetrack2->Eta(),savetrack2->Phi(),0.105);
     // fPt0 = d1.Pt();
     // fPt1 = d2.Pt();
      d1= daughter1;
      d2= daughter2;
     
     
           
      // cout<<"mass of track"<< track->M()<<endl; 
      
      
      
      TLorentzVector p = d1+d2;
      
     
      p = parent;
      fPt = p.Pt();
      fM =  p.M(); 
      if (p.Rapidity()>-2.5||p.Rapidity()<-4.0) return;//fHistCounter->Fill(15);  
      if (p.Rapidity()> -2.5) cout<<"this cut is not working"<<endl;
      if (p.Rapidity()< -4) cout<<"this cut is not working"<<endl;
      if (fPt>0.25)  return;      
      
      
      fPhi =  p.Phi();
     
      fTheta = p.Theta();
      
      fHelicityTheta= CosThetaHelicityFrame(d1,d2,p);
      fCollinTheta= CosThetaCollinsSoper(d1,d2,p);
      
      fHelicityPhi= CosPhiHelicityFrame(d1,d2,p);
      fCollinPhi=  CosPhiCollinsSoper(d1,d2,p);
      
      
      
      
      
      
      
     // fPp  = TMath::Sqrt(p.Pt()*p.Pt()+p.Pz()*p.Pz());
      fTPCcluster1 = savetrack1->GetTPCNcls();
     // fEta1 = savetrack1->Eta();
      fTPCcluster2 = savetrack2->GetTPCNcls();
    //  fEta2 = savetrack2->Eta();
      
      
      fTree ->Fill();                                // continue until all the tracks are processed
      PostData(1, fOutputList);                           // stream the results the analysis of this event to
      PostData (2,fTree);                            // the output manager which will take care of writing
                                                          // it to a file
  
      }
      
      



Double_t Polarization::CosPhiHelicityFrame(  TLorentzVector muonPositive,
                                                          TLorentzVector muonNegative,
                                                          TLorentzVector possibleJPsi )
{
  /* - This function computes the helicity phi for the
     - helicity of the J/Psi.
     - The idea should be to get back to a reference frame where it
     - is easier to compute and to define the proper z-axis.
     -
   */

  /* - Half of the energy per pair of the colliding nucleons.
     -
  */
  Double_t HalfSqrtSnn   = 2510.;
  Double_t MassOfLead208 = 193.6823;
  Double_t MomentumBeam  = TMath::Sqrt( HalfSqrtSnn*HalfSqrtSnn*208*208 - MassOfLead208*MassOfLead208 );
  /* - Fill the Lorentz vector for projectile and target.
     - For the moment we do not consider the crossing angle.
     - Projectile runs towards the MUON arm.
     -
   */
  TLorentzVector pProjCM(0.,0., -MomentumBeam, HalfSqrtSnn*208); // projectile
  TLorentzVector pTargCM(0.,0.,  MomentumBeam, HalfSqrtSnn*208); // target
  /* - Translate the dimuon parameters in the dimuon rest frame
     -
   */
  TVector3       beta      = ( -1./possibleJPsi.E() ) * possibleJPsi.Vect();
  TLorentzVector pMu1Dimu  = muonPositive;
  TLorentzVector pMu2Dimu  = muonNegative;
  TLorentzVector pProjDimu = pProjCM;
  TLorentzVector pTargDimu = pTargCM;
  pMu1Dimu.Boost(beta);
  pMu2Dimu.Boost(beta);
  pProjDimu.Boost(beta);
  pTargDimu.Boost(beta);
  //
  // --- Determine the z axis for the calculation of the polarization angle
  // (i.e. the direction of the dimuon in the CM system)
  //
  TVector3 zaxis = (possibleJPsi.Vect()).Unit();
  TVector3 yaxis = ((pProjDimu.Vect()).Cross(pTargDimu.Vect())).Unit();
  TVector3 xaxis = (yaxis.Cross(zaxis)).Unit();
  //
  // --- Calculation of the azimuthal angle (Helicity)
  //
  Double_t phi = TMath::ATan2((pMu1Dimu.Vect()).Dot(yaxis),(pMu1Dimu.Vect()).Dot(xaxis));
  return   phi;
}      
      


Double_t Polarization::CosPhiCollinsSoper( TLorentzVector muonPositive,
                                                        TLorentzVector muonNegative,
                                                        TLorentzVector possibleJPsi )
{
  /* - This function computes the Collins-Soper PHI for the
     - helicity of the J/Psi.
     - The idea should be to get back to a reference frame where it
     - is easier to compute and to define the proper z-axis.
     -
   */

  /* - Half of the energy per pair of the colliding nucleons.
     -
   */
  Double_t HalfSqrtSnn   = 2510.;
  Double_t MassOfLead208 = 193.6823;
  Double_t MomentumBeam  = TMath::Sqrt( HalfSqrtSnn*HalfSqrtSnn*208*208 - MassOfLead208*MassOfLead208 );
  /* - Fill the Lorentz vector for projectile and target.
     - For the moment we do not consider the crossing angle.
     - Projectile runs towards the MUON arm.
     -
   */
  TLorentzVector pProjCM(0.,0., -MomentumBeam, HalfSqrtSnn*208); // projectile
  TLorentzVector pTargCM(0.,0.,  MomentumBeam, HalfSqrtSnn*208); // target
  /* - Translate the dimuon parameters in the dimuon rest frame
     -
   */
  TVector3       beta      = ( -1./possibleJPsi.E() ) * possibleJPsi.Vect();
  TLorentzVector pMu1Dimu  = muonPositive;
  TLorentzVector pMu2Dimu  = muonNegative;
  TLorentzVector pProjDimu = pProjCM;
  TLorentzVector pTargDimu = pTargCM;
  pMu1Dimu.Boost(beta);
  pMu2Dimu.Boost(beta);
  pProjDimu.Boost(beta);
  pTargDimu.Boost(beta);
  /* - Determine the z axis for the CS angle.
     -
   */
  TVector3 zaxisCS=(((pProjDimu.Vect()).Unit())-((pTargDimu.Vect()).Unit())).Unit();
  //
  // --- Determine the CS angle (angle between mu+ and the z axis defined above)
  //
  TVector3 yaxisCS=(((pProjDimu.Vect()).Unit()).Cross((pTargDimu.Vect()).Unit())).Unit();
  TVector3 xaxisCS=(yaxisCS.Cross(zaxisCS)).Unit();

  Double_t phi = TMath::ATan2((pMu1Dimu.Vect()).Dot(yaxisCS),((pMu1Dimu.Vect()).Dot(xaxisCS)));
  return   phi;
}



Double_t Polarization::CosThetaCollinsSoper( TLorentzVector muonPositive,
                                                          TLorentzVector muonNegative,
                                                          TLorentzVector possibleJPsi )
{
  /* - This function computes the Collins-Soper cos(theta) for the
     - helicity of the J/Psi.
     - The idea should be to get back to a reference frame where it
     - is easier to compute and to define the proper z-axis.
     -
   */

  /* - Half of the energy per pair of the colliding nucleons.
     -
   */
  Double_t HalfSqrtSnn   = 2510.;
  Double_t MassOfLead208 = 193.6823;
  Double_t MomentumBeam  = TMath::Sqrt( HalfSqrtSnn*HalfSqrtSnn*208*208 - MassOfLead208*MassOfLead208 );
  /* - Fill the Lorentz vector for projectile and target.
     - For the moment we do not consider the crossing angle.
     - Projectile runs towards the MUON arm.
     -
   */
  TLorentzVector pProjCM(0.,0., -MomentumBeam, HalfSqrtSnn*208); // projectile
  TLorentzVector pTargCM(0.,0.,  MomentumBeam, HalfSqrtSnn*208); // target
  /* - Translate the dimuon parameters in the dimuon rest frame
     -
   */
  TVector3       beta      = ( -1./possibleJPsi.E() ) * possibleJPsi.Vect();
  TLorentzVector pMu1Dimu  = muonPositive;
  TLorentzVector pMu2Dimu  = muonNegative;
  TLorentzVector pProjDimu = pProjCM;
  TLorentzVector pTargDimu = pTargCM;
  pMu1Dimu.Boost(beta);
  pMu2Dimu.Boost(beta);
  pProjDimu.Boost(beta);
  pTargDimu.Boost(beta);
  /* - Determine the z axis for the CS angle.
     -
   */
  TVector3 zaxisCS=(((pProjDimu.Vect()).Unit())-((pTargDimu.Vect()).Unit())).Unit();
  /* - Determine the CS angle (angle between mu+ and the z axis defined above)
     -
   */
  Double_t CosThetaCS = zaxisCS.Dot((pMu1Dimu.Vect()).Unit());
  return   CosThetaCS;
}


Double_t Polarization::CosThetaHelicityFrame( TLorentzVector muonPositive,
                                                           TLorentzVector muonNegative,
                                                           TLorentzVector possibleJPsi )
{
  /* - This function computes the Helicity cos(theta) for the
     - helicity of the J/Psi.
     - The idea should be to get back to a reference frame where it
     - is easier to compute and to define the proper z-axis.
     -
   */

  /* - Half of the energy per pair of the colliding nucleons.
     -
   */
  Double_t HalfSqrtSnn   = 2510.;
  Double_t MassOfLead208 = 193.6823;
  Double_t MomentumBeam  = TMath::Sqrt( HalfSqrtSnn*HalfSqrtSnn*208*208 - MassOfLead208*MassOfLead208 );
  /* - Fill the Lorentz vector for projectile and target.
     - For the moment we do not consider the crossing angle.
     - Projectile runs towards the MUON arm.
     -
   */
  TLorentzVector pProjCM(0.,0., -MomentumBeam, HalfSqrtSnn*208); // projectile
  TLorentzVector pTargCM(0.,0.,  MomentumBeam, HalfSqrtSnn*208); // target
  /* - Translate the dimuon parameters in the dimuon rest frame
     -
   */
  TVector3       beta      = ( -1./possibleJPsi.E() ) * possibleJPsi.Vect();
  TLorentzVector pMu1Dimu  = muonPositive;
  TLorentzVector pMu2Dimu  = muonNegative;
  TLorentzVector pProjDimu = pProjCM;
  TLorentzVector pTargDimu = pTargCM;
  pMu1Dimu.Boost(beta);
  pMu2Dimu.Boost(beta);
  pProjDimu.Boost(beta);
  pTargDimu.Boost(beta);
  //
  // --- Determine the z axis for the calculation of the polarization angle
  // (i.e. the direction of the dimuon in the CM system)
  //
  TVector3 zaxis = (possibleJPsi.Vect()).Unit();
  /* - Determine the He angle (angle between mu+ and the z axis defined above)
     -
   */
  Double_t CosThetaHE = zaxis.Dot((pMu1Dimu.Vect()).Unit());
  return   CosThetaHE;

}




      
//_____________________________________________________________________________
void Polarization::Terminate(Option_t *)
  {
    // terminate
   // called at the END of the analysis (when all events are processed)
  }
  //_____________________________________________________________________________
