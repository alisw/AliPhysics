/**************************************************************************
 * Authors: Martha Spyropoulou-Stassinaki and the  members 
 * of the Greek group at Physics Department of Athens University
 * Paraskevi Ganoti, Anastasia Belogianni and Filimon Roukoutakis 
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

//-----------------------------------------------------------------
//                 AliAnalysisKinkESDMC class
//       Example of an analysis task for kink topology study
//      Kaons from kink topology are 'identified' in this code
//-----------------------------------------------------------------

#include "TF1.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"

#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliESDkink.h"

#include "AliAnalysisKinkESDMC.h"

ClassImp(AliAnalysisKinkESDMC)
//________________________________________________________________________
AliAnalysisKinkESDMC::AliAnalysisKinkESDMC(const char *name) 
  : AliAnalysisTaskSE(name), fHistPtESD(0),fHistPt(0),fHistQtAll(0),fHistQt1(0),fHistQt2(0)
  , fHistPtKaon(0),fHistPtKPDG(0),fHistEta(0),fHistEtaK(0),fptKMC(0),fMultiplMC(0),fESDMult(0),fgenpt(0),frad(0),
  fKinkKaon(0), fKinkKaonBg(0), fM1kaon(0),  fgenPtEtR(0),fPtKink(0),  fptKink(0),
   fcodeH(0), fdcodeH(0), fAngMomK(0),fAngMomPi(0), fAngMomKC(0),  fMultESDK(0), fMultMCK(0),
  fRpr(0),fZpr(0),
   fZvXv(0),fZvYv(0), fXvYv(0), fPtPrKink(0),  f1(0), f2(0),
      fListOfHistos(0)

{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
 // DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliAnalysisKinkESDMC::UserCreateOutputObjects() 
{
  // Create histograms
  // Called once
  
   f1=new TF1("f1","((atan([0]*[1]*(1.0/(sqrt((x^2)*(1.0-([1]^2))-([0]^2)*([1]^2))))))*180.)/[2]",1.1,10.0);
   f1->SetParameter(0,0.493677);
   f1->SetParameter(1,0.9127037);
   f1->SetParameter(2,TMath::Pi());


   f2=new TF1("f2","((atan([0]*[1]*(1.0/(sqrt((x^2)*(1.0-([1]^2))-([0]^2)*([1]^2))))))*180.)/[2]",0.1,10.0);
   f2->SetParameter(0,0.13957018);
   f2->SetParameter(1,0.2731374);
   f2->SetParameter(2,TMath::Pi());
   //Open file  1= CAF 
    //OpenFile(1); 

  fHistPtESD = new TH1F("fHistPtESD", "P_{T} distribution",100, 0.0,10.0);
  fHistPtESD->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  fHistPtESD->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  fHistPtESD->SetMarkerStyle(kFullCircle);
  fHistPt = new TH1F("fHistPt", "P_{T} distribution",100, 0.0,10.0); 
  fHistQtAll = new TH1F("fHistQtAll", "Q_{T} distr All Kinks ",100, 0.0,.250); 
  fHistQt1= new TH1F("fHistQt1", "Q_{T} distribution",100, 0.0,.300); 
  fHistQt2= new TH1F("fHistQt2", "Q_{T} distribution",100, 0.0,.300); 
  fHistPtKaon = new TH1F("fHistPtKaon", "P_{T}Kaon distribution",100, 0.0,10.0); 
  fHistPtKPDG = new TH1F("fHistPtKPDG", "P_{T}Kaon distribution",100, 0.0,10.0); 
  fHistEta= new TH1F("fHistEta", "Eta distribution", 26,-1.3, 1.3); 
  fHistEtaK= new TH1F("fHistEtaK", "EtaK distribution", 26,-1.3, 1.3); 
  fptKMC= new TH1F("fptKMC", "P_{T}Kaon generated",100, 0.0,10.0); 
  fMultiplMC= new TH1F("fMultiplMC", "charge multiplicity MC",60, 0.5,120.5); 
  fESDMult= new TH1F("fESDMult", "charge multipliESD", 60, 0.5,120.5); 
  fgenpt= new TH1F("fgenpt", "genpt   K distribution",100, 0.0,10.0); 
  frad= new TH1F("frad", "radius  K generated",100,0.,400.); 
  fKinkKaon= new TH1F("fKinkKaon", "P_{T}Kaon kinks identi",100, 0.0,10.0); 
  fKinkKaonBg= new TH1F("fKinkKaonBg", "P_{T}Kaon kinks backgr",100, 0.0,10.0); 
  fM1kaon= new TH1F("fM1kaon","Invar m(kaon) from kink->mu+netrino decay",80,0.0, 0.8); 
  fgenPtEtR= new TH1F("fgenPtEtR", "P_{T}Kaon distribution",100, 0.0,10.0); 
  fPtKink= new TH1F("fPtKink", "P_{T}Kaon Kink  bution",100, 0.0,10.0); 
  fptKink= new TH1F("fptKink", "P_{T}Kaon Kink  bution",100, 0.0,10.0); 
  fcodeH   = new TH2F("fcodeH", "code vrs dcode dist. kinks,K",100,0.,2500.,100,0.,2500.);
  fdcodeH = new TH2F("fdcodeH", "code vrs dcode dist. kinks,K",100,0.,2500.,100,0.,2500.);
  fAngMomK= new TH2F("fAngMomK","Decay angle vrs Mother Mom,K",100,0.0,5.0,80,0.,80.);
  fAngMomPi= new TH2F("fAngMomPi","Decay angle vrs Mother Mom,Pi",100,0.0,5.0,80,0.,80.);
  fAngMomKC= new TH2F("fAngMomKC","Decay angle vrs Mother Mom,K",100,0.0,5.0,80,0.,80.);
  fMultESDK=new TH1F("fMultESDK", "charge multipliESD kaons", 60, 0.5,240.5); 
  fMultMCK=new TH1F("fMultMCK", "charge multipli MC kaons",100, 0.5,600.5); 
  fRpr = new TH1D("fRpr", "rad distribution  PID pr",50,0.0, 2.5);
  fZpr = new TH1D("fZpr", "z distribution PID pr  ",60,-15.,15.);
  fZvXv= new TH2F("fZvXv","Xv-Zv main vtx",60,-0.5,0.5,60, -15., 15.0);
  fZvYv= new TH2F("fZvYv","Yv-Zv main vtx",60,-0.5,0.5, 60, -15., 15.);
  fXvYv= new TH2F("fXvYv","Xv-Yv main vtx", 60,-1.5,1.5, 60, -1.5, 1.5);
  fPtPrKink=new TH1F("fPtPrKink","pt of ESD  kaonKink tracks",100, 0.0,10.0);

   fListOfHistos=new TList();

   fListOfHistos->Add(fHistPtESD);
   fListOfHistos->Add(fHistPt);
   fListOfHistos->Add(fHistQtAll);
   fListOfHistos->Add(fHistQt1);
   fListOfHistos->Add(fHistQt2);
   fListOfHistos->Add(fHistPtKaon);
   fListOfHistos->Add(fHistPtKPDG);
   fListOfHistos->Add(fHistEta);
   fListOfHistos->Add(fHistEtaK);
   fListOfHistos->Add(fptKMC);
   fListOfHistos->Add(fMultiplMC);
   fListOfHistos->Add(fESDMult);
   fListOfHistos->Add(fgenpt);
   fListOfHistos->Add(frad);
   fListOfHistos->Add(fKinkKaon);
   fListOfHistos->Add(fKinkKaonBg);
   fListOfHistos->Add(fM1kaon);
   fListOfHistos->Add(fgenPtEtR);
   fListOfHistos->Add(fPtKink);
   fListOfHistos->Add(fptKink);
   fListOfHistos->Add(fcodeH);
   fListOfHistos->Add(fdcodeH);
   fListOfHistos->Add(fAngMomK);
   fListOfHistos->Add(fAngMomPi);
   fListOfHistos->Add(fAngMomKC);
   fListOfHistos->Add(fMultESDK);
   fListOfHistos->Add(fMultMCK);
   fListOfHistos->Add(fRpr);
   fListOfHistos->Add(fZpr);
   fListOfHistos->Add(fZvXv);
   fListOfHistos->Add(fZvYv);
   fListOfHistos->Add(fXvYv);
   fListOfHistos->Add(fPtPrKink);
}

//________________________________________________________________________
void AliAnalysisKinkESDMC::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event

  // Process MC truth, therefore we receive the AliAnalysisManager and ask it for the AliMCEventHandler
  // This handler can return the current MC event
  
   AliVEvent *event = InputEvent();
  if (!event) {
     Printf("ERROR: Could not retrieve event");
     return;
  }

  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(event);
  if (!esd) {
     Printf("ERROR: Could not retrieve esd");
     return;
  }

  AliMCEvent* mcEvent = MCEvent();
  if (!mcEvent) {
     Printf("ERROR: Could not retrieve MC event");
     return;
  }

  Printf("MC particles: %d", mcEvent->GetNumberOfTracks());

  AliStack* stack=mcEvent->Stack();
  
//primary tracks  in MC
         Int_t  nPrim = stack->GetNprimary();
//
  const AliESDVertex *vertex=GetEventVertex(esd);
    if(!vertex) return;
//
  Double_t vpos[3];
  vertex->GetXYZ(vpos);
    fZpr->Fill(vpos[2]);         

  Double_t vtrack[3], ptrack[3];
  
     
 // Int_t nESDTracK = 0;

   Int_t nGoodTracks =  esd->GetNumberOfTracks();
    fESDMult->Fill(nGoodTracks);
     
//
// track loop
   for (Int_t iTracks = 0; iTracks < esd->GetNumberOfTracks(); iTracks++) {
// loop only on primary tracks

    AliESDtrack* track = esd->GetTrack(iTracks);
    if (!track) {
      Printf("ERROR: Could not receive track %d", iTracks);
      continue;
    }
    

    fHistPt->Fill(track->Pt());

    
    Int_t label = track->GetLabel();
    label = TMath::Abs(label);
    TParticle * part = stack->Particle(label);
    if (!part) continue;
// loop only on Primary tracks
      if (label > nPrim) continue; /// primary tracks only   , 26/8/09  EFF study

      //    pt cut at 300 MeV
        if ( (track->Pt())<.300)continue;

//    UInt_t status=track->GetStatus();

  //  if((status&AliESDtrack::kITSrefit)==0) continue;
  //  if((status&AliESDtrack::kTPCrefit)==0) continue;
  // ayksanei to ratio BG/real    if((track->GetTPCchi2()/track->GetTPCclusters(0))>4.5) continue;
     if((track->GetTPCchi2()/track->GetTPCclusters(0))>4.5) continue;

      Double_t extCovPos[15];
      track->GetExternalCovariance(extCovPos);    
    //  if(extCovPos[0]>2) continue;
    //  if(extCovPos[2]>2) continue;    
    //  if(extCovPos[5]>0.5) continue;  
    //  if(extCovPos[9]>0.5) continue;
    //  if(extCovPos[14]>2) continue;


    track->GetXYZ(vtrack);
 fXvYv->Fill(vtrack[0],vtrack[1]);  
 fZvYv->Fill(vtrack[0],vtrack[2]);  
 fZvXv->Fill(vtrack[1],vtrack[2]);  

// track momentum
     track->GetPxPyPz(ptrack);
    
    TVector3 trackMom(ptrack[0],ptrack[1],ptrack[2]);
    
    Double_t trackEta=trackMom.Eta();
    
    
    
    Float_t nSigmaToVertex = GetSigmaToVertex(track);  
      
    Float_t bpos[2];
    Float_t bCovpos[3];
    track->GetImpactParameters(bpos,bCovpos);
    
    if (bCovpos[0]<=0 || bCovpos[2]<=0) {
     Printf("Estimated b resolution lower or equal zero!");
     bCovpos[0]=0; bCovpos[2]=0;
    }

    Float_t dcaToVertexXYpos = bpos[0];
    Float_t dcaToVertexZpos = bpos[1];
    
    fRpr->Fill(dcaToVertexZpos);
//  test for secondary kinks   , 5/7/2009  
    if(nSigmaToVertex>=4) continue;
 //    if((dcaToVertexXYpos>3.0)||(dcaToVertexZpos>3.0)) continue;  //arxiko 
     if((dcaToVertexXYpos>4.0)||(dcaToVertexZpos>5.0)) continue;   // 27/8
    

 //  cut on eta 
        if(  (TMath::Abs(trackEta )) > 0.9 ) continue;
    fHistPtESD->Fill(track->Pt());

   // Add Kink analysis
   
   	    	Int_t indexKinkPos=track->GetKinkIndex(0);
//  loop on kinks
		if(indexKinkPos<0){
               fPtKink->Fill(track->Pt()); /// pt from track

	// select kink class	

	  AliESDkink *kink=esd->GetKink(TMath::Abs(indexKinkPos)-1);
//
	
	  Int_t eSDfLabel1=kink->GetLabel(0);
	  TParticle *particle1= stack->Particle(TMath::Abs(eSDfLabel1));
          Int_t code1= particle1->GetPdgCode();
	  
	  Int_t eSDfLabeld=kink->GetLabel(1);
	  TParticle *particled= stack->Particle(TMath::Abs(eSDfLabeld));
          Int_t dcode1= particled->GetPdgCode();

	   const TVector3 motherMfromKink(kink->GetMotherP());
	   const TVector3 daughterMKink(kink->GetDaughterP());
	   Float_t qT=kink->GetQt();

           fHistQtAll->Fill(qT) ;  //  Qt   distr
                  
           fptKink->Fill(motherMfromKink.Pt()); /// pt from kink

           Float_t kinkAngle=TMath::RadToDeg()*kink->GetAngle(2);

//       fake kinks, small Qt and small kink angle
    if(( kinkAngle<1.)&&(TMath::Abs(code1)==321)&&(TMath::Abs(dcode1)==13)) fHistQt1  ->Fill(qT) ;  //  Qt   distr

           if(qT<0.012) continue;  // remove fake kinks

//remove the double taracks 
           if( (kinkAngle<1.)  ) continue;

//

      if((kink->GetR()>120.)&&(kink->GetR()<220.)&&(TMath::Abs(trackEta)<0.9)&&(label<nPrim)) {
         if(       ( (TMath::Abs(code1)==321)&&(TMath::Abs(dcode1)==13))||    
           ( (TMath::Abs(code1)==321)&&(TMath::Abs(dcode1)==11))  ||    
           ( (TMath::Abs(code1)==321)&&(TMath::Abs(dcode1)==211))  ) { 
               fHistPtKPDG->Fill(track->Pt());  // ALL KAONS (pdg) inside ESD  kink sample
           fHistEta->Fill(trackEta) ;  //   Eta distr of PDG kink ESD  kaons
    fMultESDK->Fill(nGoodTracks);
      fHistQt2->Fill(qT);  // PDG ESD kaons            
     }
     }

//          maximum decay angle at a given mother momentum
	   Double_t maxDecAngKmu=f1->Eval(motherMfromKink.Mag(),0.,0.,0.);
	   Double_t maxDecAngpimu=f2->Eval(motherMfromKink.Mag(),0.,0.,0.);
//  two dimensional plot 
                if(TMath::Abs(code1)==321) fAngMomK->Fill(motherMfromKink.Mag(), kinkAngle); 
                if(TMath::Abs(code1)==211) fAngMomPi->Fill(motherMfromKink.Mag(), kinkAngle); 
//
// invariant mass of mother track decaying to mu
	 Float_t energyDaughterMu=TMath::Sqrt(daughterMKink.Mag()*daughterMKink.Mag()+0.105658*0.105658);
	 Float_t p1XM= motherMfromKink.Px();
         Float_t p1YM= motherMfromKink.Py();
         Float_t p1ZM= motherMfromKink.Pz();
         Float_t p2XM= daughterMKink.Px();
         Float_t p2YM= daughterMKink.Py();
         Float_t p2ZM= daughterMKink.Pz();
         Float_t p3Daughter=TMath::Sqrt(((p1XM-p2XM)*(p1XM-p2XM))+((p1YM-p2YM)*(p1YM-p2YM))+((p1ZM-p2ZM)*(p1ZM-p2ZM)));
         Double_t invariantMassKmu= TMath::Sqrt((energyDaughterMu+p3Daughter)*(energyDaughterMu+p3Daughter)-motherMfromKink.Mag()*motherMfromKink.Mag());
  
         fM1kaon->Fill(invariantMassKmu);

//  kaon selection from kinks
//if((kinkAngle>maxDecAngpimu)&&(qT>0.05)&&(qT<0.25)&&((kink->GetR()>120.)&&(kink->GetR()<220.))&&(TMath::Abs(trackEta)<0.9)&&(invariantMassKmu<0.6)) {
if((kinkAngle>maxDecAngpimu)&&(qT>0.04)&&(qT<0.30)&&((kink->GetR()>120.)&&(kink->GetR()<220.))&&(TMath::Abs(trackEta)<0.9)&&(invariantMassKmu<0.6)) {

        if( (kinkAngle>maxDecAngKmu*1.1) && ( motherMfromKink.Mag()> 1.2 ) ) fcodeH->Fill(TMath::Abs(code1), TMath::Abs(dcode1));
                if ( (kinkAngle>maxDecAngKmu*1.1) && ( motherMfromKink.Mag()> 1.2 ) )  fAngMomKC->Fill(motherMfromKink.Mag(), kinkAngle);


                 if ( (kinkAngle>maxDecAngKmu*1.1) && ( motherMfromKink.Mag()> 1.2 ) ) continue; // maximum angle selection revised 22/10/2009

                                    fHistPtKaon->Fill(track->Pt());   //all PID kink-kaon

//  kaons from k to mun and k to pipi and to e decay 
         if(       ( (TMath::Abs(code1)==321)&&(TMath::Abs(dcode1)==13))||    
           ( (TMath::Abs(code1)==321)&&(TMath::Abs(dcode1)==11))  ||    
           ( (TMath::Abs(code1)==321)&&(TMath::Abs(dcode1)==211))  ) { 

         if ( label<=nPrim ) fPtPrKink->Fill(track->Pt());
             fKinkKaon->Fill(track->Pt());        
                fHistEtaK->Fill(trackEta);
                             }
         else {
             fKinkKaonBg->Fill(track->Pt());     
 fdcodeH->Fill( TMath::Abs(code1), TMath::Abs(dcode1));   // put it here,  22/10/2009
          }   // primary and all +BG    

        }  //  kink selection 
                  

	}  //End Kink Information    
  

  } //track loop 

  // loop over mc particles

  // variable to count tracks
  Int_t nMCTracks = 0;
  Int_t nMCKinkKs = 0;

for (Int_t iMc = 0; iMc < nPrim; iMc++)
  {
   // AliDebug(AliLog::kDebug+1, Form("MC Loop: Processing particle %d.", iMc));

    TParticle* particle = stack->Particle(iMc);

    if (!particle)
    {
    //  AliDebug(AliLog::kError, Form("UNEXPECTED: particle with label %d not found in stack (mc loop).", iMc));
      continue;
    }

 //   
    Float_t eta = particle->Eta();

    if (TMath::Abs(eta) > 0.9)
      continue;

            Double_t ptK = particle->Pt();

   if( ptK <0.300) continue;

      Int_t code = particle->GetPdgCode();
            Int_t  mcProcess=-1011;
      if ((code==321)||(code==-321)){
	    

	      fptKMC->Fill(ptK);

    
  // fMultiplMC->Fill(nPrim);
       Int_t nMCKpi =0;

	    Int_t firstD=particle->GetFirstDaughter();
	    Int_t lastD=particle->GetLastDaughter();
//loop on secondaries
	     for (Int_t k=firstD;k<=lastD;k++) {
              if ( k > 0 )    {
	     TParticle*    daughter1=stack->Particle(k);   // 27/8   
	     Int_t dcode = daughter1->GetPdgCode();

        mcProcess=daughter1->GetUniqueID();
     if (mcProcess==4) {        
    

		 // frad->Fill(daughter1->R());

		 if (((daughter1->R())>120)&&((daughter1->R())<220) ){
        if (( ( code==321 )&& ( dcode ==-13  ))|| (( code == -321 )&& ( dcode == 13))) fgenPtEtR->Fill( ptK );//to muon
      //  if (( ( code==321 )&& ( dcode ==211  ))|| (( code == -321 )&& ( dcode ==-211))) fgenPtEtR->Fill( ptK );//to pion
        if (( ( code==321 )&& ( dcode ==-11  ))|| (( code == -321 )&& ( dcode ==11))) fgenPtEtR->Fill( ptK );// to electr
        if (( ( code==321 )&& ( dcode ==211  ))|| (( code == -321 )&& ( dcode ==-211)))    nMCKpi++ ; 
 
  //fMultMCK->Fill(nPrim);
		 frad->Fill(daughter1->R());
    nMCKinkKs++;
       }
//  next only  to mu decay
		 if (((code==321)&&(dcode==-13))||((code==-321)&&(dcode==13))){

		 
		    if (((daughter1->R())>120)&&((daughter1->R())<220)){
		 //   fgenpt->Fill(ptK);
		    }
		    }
		    }//    decay
         } // positive k
       }//  daughters
              if( nMCKpi == 1) fgenPtEtR->Fill(ptK);  //  k to pipi
              if( nMCKpi > 1) fgenpt->Fill(ptK);// k to pipipi

      }   /// kaons


      
    nMCTracks++;
  }// end of mc particle

//  printf("hello mc");
  fMultiplMC->Fill(nMCTracks);
  if( nMCKinkKs>0) fMultMCK->Fill(nPrim);

  PostData(1, fListOfHistos);

}      

//________________________________________________________________________
void AliAnalysisKinkESDMC::Terminate(Option_t *) 
{
  // Draw result to the screen 
  // Called once at the end of the query
   fHistPtKaon = dynamic_cast<TH1F*> (((TList*)GetOutputData(1))->FindObject("fHistPtKaon"));
   if (!fHistPtKaon) {
     Printf("ERROR: fHistPtKaon not available");
     return;
   }
   
   TCanvas *canCheckKinkESDMC = new TCanvas("AliAnalysisKinkESDMC","Check KinkESDMC",1);
   canCheckKinkESDMC->Draw();
   canCheckKinkESDMC->cd();
   fHistPtKaon->DrawCopy("E");

}
//____________________________________________________________________//

Float_t AliAnalysisKinkESDMC::GetSigmaToVertex(AliESDtrack* esdTrack) const {
  // Calculates the number of sigma to the vertex.
  
  Float_t b[2];
  Float_t bRes[2];
  Float_t bCov[3];

    esdTrack->GetImpactParameters(b,bCov);
  
  if (bCov[0]<=0 || bCov[2]<=0) {
    //AliDebug(1, "Estimated b resolution lower or equal zero!");
    bCov[0]=0; bCov[2]=0;
  }
  bRes[0] = TMath::Sqrt(bCov[0]);
  bRes[1] = TMath::Sqrt(bCov[2]);
  
  if (bRes[0] == 0 || bRes[1] ==0) return -1;
  
  Float_t d = TMath::Sqrt(TMath::Power(b[0]/bRes[0],2) + TMath::Power(b[1]/bRes[1],2));
  
  if (TMath::Exp(-d * d / 2) < 1e-10) return 1000;
  
  d = TMath::ErfInverse(1 - TMath::Exp(-d * d / 2)) * TMath::Sqrt(2);
  
  return d;
  }
//____________________________________________________________________//

const AliESDVertex* AliAnalysisKinkESDMC::GetEventVertex(const AliESDEvent* esd) const
  // Get the vertex from the ESD and returns it if the vertex is valid
  
{
  // Get the vertex 
  
  const AliESDVertex* vertex = esd->GetPrimaryVertex();

  if((vertex->GetStatus()==kTRUE)&&(vertex->GetNContributors()>2)) return vertex;
  else
  { 
     vertex = esd->GetPrimaryVertexSPD();
     if((vertex->GetStatus()==kTRUE)&&(vertex->GetNContributors()>2)) return vertex;
     else
     return 0;
  }
}
