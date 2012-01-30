/**************************************************************************
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
//                 AliAnalysisKinkESDat class
//       Example of an analysis task for kink topology study
//      Kaons from kink topology are 'identified' in this code
//-----------------------------------------------------------------

#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TParticle.h"
#include <TVector3.h>
#include "TF1.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliAnalysisKinkESDat.h"
#include "AliStack.h"
#include "AliESDpid.h"
#include "AliPID.h"
#include "AliESDkink.h"
#include "AliESDtrack.h"
#include "AliPhysicsSelectionTask.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliVEvent.h"
//#include "AddTaskTender.h"
//  #include "AliTPCpidESD.h"
#include "AliESDtrackCuts.h"
ClassImp(AliAnalysisKinkESDat)


//________________________________________________________________________
AliAnalysisKinkESDat::AliAnalysisKinkESDat(const char *name) 
  : AliAnalysisTaskSE(name), fHistPtESD(0),fHistPt(0),fHistQtAll(0),fHistQt1(0),fHistQt2(0)
  , fHistPtKaon(0),fHistPtKPDG(0),fHistEta(0),fHistEtaK(0),fptKMC(0),fMultiplMC(0),fESDMult(0),fgenpt(0),frad(0),
  fKinkKaon(0),fKinKRbn(0), fKinkKaonBg(0), fM1kaon(0),  fgenPtEtR(0),fPtKink(0),  fptKink(0),
   fcodeH(0), fdcodeH(0), fAngMomK(0),fAngMomPi(0), fAngMomKC(0),  fMultESDK(0), fMultMCK(0),
 fSignPtNcl(0), fSignPtEta(0), fEtaNcl(0), fSignPt(0), fChi2NclTPC(0), fRatChi2Ncl(0), fRadiusNcl(0), fTPCSgnlP(0),
   fTPCSgnlPa(0), fRpr(0),fZpr(0), fdcatoVxXY(0), fnSigmToVx(0), fKinkMothDau(0),
 fZvXv(0),fZvYv(0), fXvYv(0), fPtPrKink(0), fHistPtKaoP(0), fHistPtKaoN(0),frapiKESD(0), flifetime(), fradLK(0),
    fradPtRpDt(0), fInvMuNuAll(0), fQtInvM(0), 
         fDCAkink(0), fPosiKink(0),  fPosiKinkK(0),fPosiKinKXZ(0), fPosiKinKYZ(0), fPosiKinKBg(0), fQtMothP(0),
 f1(0), f2(0),
      fListOfHistos(0),fLowMulcut(-1),fUpMulcut(-1),fCutsMul(0)

{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
 // DefineInput(0, TChain::Class());
//----------------------Marek multiplicity bins 
 fCutsMul=new AliESDtrackCuts("Mul","Mul");
        fCutsMul->SetMinNClustersTPC(70);
        fCutsMul->SetMaxChi2PerClusterTPC(4);
        fCutsMul->SetAcceptKinkDaughters(kFALSE);
        fCutsMul->SetRequireTPCRefit(kTRUE);
        // ITS
        fCutsMul->SetRequireITSRefit(kTRUE);
        fCutsMul->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                                AliESDtrackCuts::kAny);
        fCutsMul->SetMaxDCAToVertexXYPtDep("0.0182+0.0350/pt^1.01");

        fCutsMul->SetMaxDCAToVertexZ(2);
        fCutsMul->SetDCAToVertex2D(kFALSE);
        fCutsMul->SetRequireSigmaToVertex(kFALSE);

        fCutsMul->SetEtaRange(-0.8,+0.8);
        fCutsMul->SetPtRange(0.15, 1e10);


  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliAnalysisKinkESDat::UserCreateOutputObjects() 
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
 //  Double_t gPt[31] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,0.9,1.0,
   //                     1.1, 1.2, 1.3, 1.4, 1.5, 1.6,1.7,1.8,1.9,  2.0,
     //                    2.2, 2.4, 2.6, 2.8, 3.0, 3.3, 3.6,3.9, 4.2, 4.5, 4.8};
   Double_t gPt7[43] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,0.9,1.0,
                        1.1, 1.2, 1.3, 1.4, 1.5, 1.6,1.7,1.8,1.9,  2.0,
                         2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7,2.8, 2.9, 3.0, 
                         3.2, 3.4, 3.6, 3.8, 4.0, 4.4, 4.8,5.2, 5.6, 6.0,  7.0, 8.0,10.0 };

  fHistPtESD = new TH1F("fHistPtESD", "P_{T} distribution",50, 0.0,5.0);
  fHistPtESD->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  fHistPtESD->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  fHistPtESD->SetMarkerStyle(kFullCircle);
  fHistPt = new TH1F("fHistPt", "P_{T} distribution",50, 0.0,5.0); 
  fHistQtAll = new TH1F("fHistQtAll", "Q_{T} distr All Kinks ",100, 0.0,.300); 
  fHistQt1= new TH1F("fHistQt1", "Q_{T} distribution",100, 0.0,.300); 
  fHistQt2= new TH1F("fHistQt2", "Q_{T} distribution",100, 0.0,.300); 
  //fHistPtKaon = new TH1F("fHistPtKaon", "P_{T}Kaon distribution",50, 0.0,5.0); 
  fHistPtKaon = new TH1F("fHistPtKaon", "P_{T}Kaon distribution",100, 0.0,10.0); 
  fHistPtKPDG = new TH1F("fHistPtKPDG", "P_{T}Kaon distribution",50, 0.0,5.0); 
  fHistEta= new TH1F("fHistEta", "Eta distribution", 26,-1.3, 1.3); 
  fHistEtaK= new TH1F("fHistEtaK", "EtaK distribution", 26,-1.3, 1.3); 
  fptKMC= new TH1F("fptKMC", "P_{T}Kaon generated",100, 0.0,10.0); 
  fMultiplMC= new TH1F("fMultiplMC", "charge multiplicity MC",100, 0.0,300.0); 
  fESDMult= new TH1F("fESDMult", "charge multipliESD",100, 0.0,300.0); 
  fgenpt= new TH1F("fgenpt", "genpt   K distribution",50, 0.0,5.0); 
   //frad= new TH1F("frad", "radius  K generated",100, 50., 250.0);
   frad= new TH1F("frad", "radius  K generated",100, 0.,1000.0);
  // fKinkKaon= new TH1F("fKinkKaon", "P_{T}Kaon kinks identi",50, 0.0,5.0); 
  fKinkKaon= new TH1F("fKinkKaon", "P_{T}Kaon kinks identi",100, 0.0,10.0); 
  fKinKRbn= new TH1F("fKinKRbn", "p_{t}Kaon kinks identi[GeV/c],Entries",42,gPt7); 
  fKinkKaonBg= new TH1F("fKinkKaonBg", "P_{T}Kaon kinks backgr",50, 0.0,5.0); 
  fM1kaon= new TH1F("fM1kaon","Invar m(kaon) from kink->mu+netrino decay",80,0.0, 0.8); 
  fgenPtEtR= new TH1F("fgenPtEtR", "P_{T}Kaon distribution",50, 0.0,5.0); 
  fPtKink= new TH1F("fPtKink", "P_{T}Kaon Kink  bution",50, 0.0,5.0); 
  fptKink= new TH1F("fptKink", "P_{T}Kaon Kink  bution",50, 0.0,5.0); 
  fcodeH   = new TH2F("fcodeH", "code vrs dcode dist. kinks,K",100,0.,2500.,100,0.,2500.);
  fdcodeH = new TH2F("fdcodeH", "code vrs dcode dist. kinks,K",100,0.,2500.,100,0.,2500.);
  fAngMomK= new TH2F("fAngMomK","Decay angle vrs Mother Mom,K",100,0.0,5.0,80,0.,80.);
  fAngMomPi= new TH2F("fAngMomPi","Decay angle vrs Mother Mom,Pi",100,0.0,5.0,80,0.,80.);
  fAngMomKC= new TH2F("fAngMomKC","Decay angle vrs Mother Mom,K",100,0.0,5.0,80,0.,80.);
  fMultESDK=new TH1F("fMultESDK", "charge multipliESD kaons",100, 0.0,100.0); 
  fMultMCK=new TH1F("fMultMCK", "charge multipli MC kaons",100, 0.0,100.0); 
  fSignPtNcl= new TH2F("fSignPtNcl","SignPt vrs Ncl,K",80,-4.,4.0,70,20.,160.);
  fSignPtEta= new TH2F("fSignPtEta","SignPt vrs Eta,K",80,-4.0,4.0,30,-1.5,1.5);
  fEtaNcl= new TH2F("fEtaNcl","Eta vrs nclust,K",30,-1.5,1.5, 70,20, 160);
  //fSignPt= new TH1F("fSignPt","SignPt ,K",40,-4.0,4.0);
  fSignPt= new TH1F("fSignPt","SignPt ,K",80,-4.0,4.0);
  fChi2NclTPC= new TH2F("fChi2NclTPC","Chi2vrs nclust,K",100,0.,500., 70,20, 160);
  fRatChi2Ncl= new TH1F("fRatChi2Ncl","Ratio chi2/nclusters in TPC,K",50,0.0,5.0);
  fRadiusNcl = new TH2F("fRadiusNcl","kink radius vrs Nclust,K",75,100.,250., 70,20, 160);
    fTPCSgnlP = new TH2F("fTPCSgnlP","TPC signal de/dx Mom,K",100,0.0,4.0,100,0.,250.);
  fTPCSgnlPa= new TH2F("fTPCSgnlPa","TPC signal de/dx Mom,K",100,0.0,4.,100, 0.,250.);
  fRpr = new TH1D("fRpr", "rad distribution  PID pr",100,-10.0, 10.0);
  fZpr = new TH1D("fZpr", "z distribution PID pr  ",80,-20.,20.);
  fdcatoVxXY = new TH1D("fdcatoVxXY", "dca  distribution PID  ",20,-1.,1.);
  fnSigmToVx = new TH1D("fnSigmToVx", "dca  distribution PID  ",80,0.,8.);
  fKinkMothDau= new TH2F("fKinkMothDau","TPC kink Moth Daugh ,K",50,0.0,2.5,50, 0., 2.5);
  fZvXv= new TH2F("fZvXv","Xv-Zv main vtx",60,-0.5,0.5,60, -15., 15.0);
  fZvYv= new TH2F("fZvYv","Yv-Zv main vtx",60,-0.5,0.5, 60, -15., 15.);
  fXvYv= new TH2F("fXvYv","Xv-Yv main vtx", 60,-1.5,1.5, 60, -1.5, 1.5);
  fPtPrKink=new TH1F("fPtPrKink","pt of ESD  kaonKink tracks",100, 0.0,10.0);
  // fHistPtKaoP = new TH1F("fHistPtKaoP", "P_{T}KaonP  distribution",50, 0.0,5.0); 
  fHistPtKaoP = new TH1F("fHistPtKaoP", "P_{T}KaonP  distribution",100, 0.0,10.0); 
  // fHistPtKaoN = new TH1F("fHistPtKaoN", "P_{T}KaonN  distribution",50, 0.0,5.0); 
  fHistPtKaoN = new TH1F("fHistPtKaoN", "P_{T}KaonN  distribution",100, 0.0,10.0); 

  frapiKESD=new TH1F("frapiKESD","rapid Kdistribution", 26,-1.3, 1.3); 
  flifetime= new TH1F("flifetime", "ct study of K-kinks",100,0.,1000.); 
  fradLK= new TH1F("fradLK", "Length of   K generated",100,0.,1000.); 
  fradPtRpDt=new TH3F("fradPtRpDt","rad pt rap dat",28,100.,240., 20, 0., 5., 20, -1., 1. );
  fInvMuNuAll= new TH1F("fInvMuNuAll", " Inv Mass MuNu all kink",80,0.,0.8); 
  fQtInvM= new TH2F("fQtInvM", "Q_{T} Versus Inv MuNu ",80, 0., 0.80 , 100 , 0., 0.300); 
    fDCAkink = new TH1F("fDCAkink ", "DCA kink vetrex ",50, 0.0,1.0);

  fPosiKink= new TH2F("fPosiKink", "Y vrx kink Vrex ",100, -300.0,300.0,100, -300, 300.);
  fPosiKinkK= new TH2F("fPosiKinkK", "Y vrx kink VrexK ",100, -300.0,300.0,100, -300, 300.);
  fPosiKinKXZ= new TH2F("fPosiKinKXZ", "Y vrx kink VrexK ",100, -300.0,300.0,100, -300, 300.);
  fPosiKinKYZ= new TH2F("fPosiKinKYZ", "Y vrx kink VrexK ",100, -300.0,300.0,100, -300, 300.);
  fPosiKinKBg= new TH2F("fPosiKinKBg", "z vrx kink rad    ",100, -300.0,300.0,100,100., 300.);
  fQtMothP = new TH2F("fQtMothP", " Qt vrs Mother P", 100, 0., 5.0,100, 0.,0.300);

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
   fListOfHistos->Add(fKinKRbn);
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
   fListOfHistos->Add(fSignPtNcl);
   fListOfHistos->Add(fSignPtEta);
   fListOfHistos->Add(fEtaNcl);
   fListOfHistos->Add(fSignPt);
   fListOfHistos->Add(fChi2NclTPC);
   fListOfHistos->Add(fRatChi2Ncl);
   fListOfHistos->Add(fRadiusNcl);
   fListOfHistos->Add(fTPCSgnlP);
   fListOfHistos->Add(fTPCSgnlPa);
   fListOfHistos->Add(fRpr);
   fListOfHistos->Add(fZpr);
   fListOfHistos->Add(fdcatoVxXY);
   fListOfHistos->Add(fnSigmToVx);
   fListOfHistos->Add(fKinkMothDau);
   fListOfHistos->Add(fZvXv);
   fListOfHistos->Add(fZvYv);
   fListOfHistos->Add(fXvYv);
   fListOfHistos->Add(fPtPrKink);
   fListOfHistos->Add(fHistPtKaoP);
   fListOfHistos->Add(fHistPtKaoN);
   fListOfHistos->Add(frapiKESD);
   fListOfHistos->Add(flifetime);
   fListOfHistos->Add(fradLK);
   fListOfHistos->Add(fradPtRpDt);
   fListOfHistos->Add(fInvMuNuAll);
   fListOfHistos->Add(fQtInvM);
   fListOfHistos->Add(fDCAkink);
   fListOfHistos->Add(fPosiKink);
   fListOfHistos->Add(fPosiKinkK);
   fListOfHistos->Add(fPosiKinKXZ);
   fListOfHistos->Add(fPosiKinKYZ);
   fListOfHistos->Add(fPosiKinKBg);
   fListOfHistos->Add(fQtMothP);


}
//=======================new thing
//     Float_t nCrossedRowsTPC = esdTrack->GetTPCClusterInfo(2,1);
//  Float_t  ratioCrossedRowsOverFindableClustersTPC = 1.0;
//  if (esdTrack->GetTPCNclsF()>0) {
//    ratioCrossedRowsOverFindableClustersTPC = nCrossedRowsTPC / esdTrack->GetTPCNclsF();
//   }
//________________________________________________________________________
void AliAnalysisKinkESDat::UserExec(Option_t *) 
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
//==================check of Physics selection?
       Bool_t isSelected =
((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()&AliVEvent::kMB;

         if ( isSelected ==kFALSE) return;   //  24/6/11 apo MF
//
//
   Int_t nESDTracks =  esd->GetNumberOfTracks();
      fMultMCK->Fill(nESDTracks);
//===============Marek multiplicity
   Float_t refmultiplicity=fCutsMul->CountAcceptedTracks(esd);
        if(fLowMulcut>-1)
        {
                if(refmultiplicity<fLowMulcut)
                        return;
        }
        if(fUpMulcut>-1)
        {
                if(refmultiplicity>fUpMulcut)
                        return;
        }



       fMultESDK->Fill(refmultiplicity);

//
//   Int_t nESDTracks =  esd->GetNumberOfTracks();
  //   if ( nESDTracks>0 ) fMultMCK->Fill(nESDTracks);



  const AliESDVertex *vertex=GetEventVertex(esd);    // 22/8
    if(!vertex) return;
      fMultiplMC->Fill(nESDTracks);
//
  Double_t vpos[3];
  vertex->GetXYZ(vpos);
    fZpr->Fill(vpos[2]);         
     if (TMath::Abs( vpos[2] ) > 10. ) return;   

    

  Double_t vtrack[3], ptrack[3];
  
     
 Int_t nESDTracK = 0;
 Int_t nESDTrKink = 0;

   Int_t nGoodTracks =  esd->GetNumberOfTracks();
    fESDMult->Fill(nGoodTracks);
      
  Double_t fAlephParameters[5] = {0.0283086,
                                  2.63394e+01,
                                  5.04114e-11,
                                  2.12543e+00,
                                  4.88663e+00};
       Double_t nsigma = 100.0;
       AliESDpid *fESDpid = new AliESDpid();                  
    fESDpid->GetTPCResponse().SetBetheBlochParameters(fAlephParameters[0],
                                                    fAlephParameters[1],
                                                   fAlephParameters[2],
                                                 fAlephParameters[3],
                                               fAlephParameters[4]);
//
// track loop
   for (Int_t iTracks = 0; iTracks < esd->GetNumberOfTracks(); iTracks++) {

    AliESDtrack* track = esd->GetTrack(iTracks);
    if (!track) {
      Printf("ERROR: Could not receive track %d", iTracks);
      continue;
    }
    
    fHistPt->Fill(track->Pt());


     //    sigmas
    nsigma = TMath::Abs(fESDpid->NumberOfSigmasTPC(track,AliPID::kKaon));



      Int_t tpcNCl = track->GetTPCclusters(0);  
      Double_t tpcSign = track->GetSign();  
    
    Int_t label = track->GetLabel();
    label = TMath::Abs(label);


    UInt_t status=track->GetStatus();

    if((status&AliESDtrack::kITSrefit)==0) continue;   
    if((status&AliESDtrack::kTPCrefit)==0) continue;
      if((track->GetTPCchi2()/track->GetTPCclusters(0))>3.8) continue;  

      Double_t extCovPos[15];
      track->GetExternalCovariance(extCovPos);    
      if(extCovPos[0]>2) continue;
     if(extCovPos[2]>2) continue;    
     if(extCovPos[5]>0.5) continue;  
     if(extCovPos[9]>0.5) continue;
     if(extCovPos[14]>2) continue;


    track->GetXYZ(vtrack);
 fXvYv->Fill(vtrack[0],vtrack[1]);  
 fZvYv->Fill(vtrack[0],vtrack[2]);  
 fZvXv->Fill(vtrack[1],vtrack[2]);  

// track momentum, rapidity calculation
     track->GetPxPyPz(ptrack);
    
    TVector3 trackMom(ptrack[0],ptrack[1],ptrack[2]);
    
          Double_t   etracK= TMath::Sqrt(trackMom.Mag()*trackMom.Mag() + 0.493677 *0.493677  );
         Double_t rapiditK = 0.5 * (TMath::Log(  (etracK + ptrack[2]  ) / ( etracK - ptrack[2])  ))  ;
    
    Double_t trackEta=trackMom.Eta();
     Double_t trMoment=trackMom.Mag();       
    Double_t trackPt = track->Pt();
    
    
      
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
 
   if((TMath::Abs(dcaToVertexXYpos)>0.3)||(TMath::Abs(dcaToVertexZpos)>3.0))  nESDTrKink++;  //  count of second  23Jul11    
//
     //   if((dcaToVertexXYpos>0.3)||(dcaToVertexZpos>0.3)) continue;   //    allagi-dokini    3/6                 
     // if((TMath::Abs(dcaToVertexXYpos)>0.4)||(dcaToVertexZpos>2.5)) continue;   //    allagi  23Jul11               
     if((TMath::Abs(dcaToVertexXYpos)>0.3)||(TMath::Abs(dcaToVertexZpos)>2.5)) continue;   //    allagi  23Jul11               
    
//  track Mult. after selection 
    nESDTracK++;        
  //    
//=========================================
    fHistPtESD->Fill(track->Pt());

   // Add Kink analysis           =============================
   
   	    	Int_t indexKinkPos=track->GetKinkIndex(0);
//  loop on kinks
		if(indexKinkPos<0){     ////mother kink
               fPtKink->Fill(track->Pt()); /// pt from track

	// select kink class	

	  AliESDkink *kink=esd->GetKink(TMath::Abs(indexKinkPos)-1);
//
         // DCA kink
          Double_t  Dist2 = kink->GetDistance();
          fDCAkink->Fill( Dist2   );
//
	
      
	   const TVector3 vposKink(kink->GetPosition());
 fPosiKink ->Fill( vposKink[0], vposKink[1]  );
//   Double_t  lengthK = TMath::Sqrt( vposKink[0]*vposKink[0] + vposKink[1]*vposKink[1] + vposKink[2]*vposKink[2] ) ;
   // Double_t dxKink = vpos[0]-vposKink[0], dyKink=vpos[1]-vposKink[1], dzKink=vpos[2]-vposKink[2]; 
   Double_t  dzKink=vpos[2]-vposKink[2]; 
   // Double_t lifeKink= TMath::Sqrt( dxKink*dxKink + dyKink*dyKink + dzKink*dzKink ) ;
//
            Double_t tanLamda = track->GetTgl();  // 25/6/2010

   Double_t lifeKink= (TMath::Abs( dzKink ))*( TMath::Sqrt(1.+ tanLamda*tanLamda) ) / (TMath::Abs( tanLamda)) ;

	   const TVector3 motherMfromKink(kink->GetMotherP());
	   const TVector3 daughterMKink(kink->GetDaughterP());

	   Float_t qT=kink->GetQt();
            Float_t motherPt=motherMfromKink.Pt();
       //     Float_t etaMother=motherMfromKink.Eta();

           fHistQtAll->Fill(qT) ;  //  Qt   distr
                  
           fptKink->Fill(motherMfromKink.Pt()); /// pt from kink

           Float_t kinkAngle=TMath::RadToDeg()*kink->GetAngle(2);

          if(  (TMath::Abs(rapiditK )) > 0.7 ) continue;
        if ( (track->Pt())<.250)continue;

                fQtMothP->Fill( track->P(), qT);

    if ( qT> 0.04)  fHistQt1  ->Fill(qT) ;  //  Qt   distr



            fHistEta->Fill(trackEta) ;  //   Eta distr of PDG kink ESD  kaons
      fHistQt2->Fill(qT);  // PDG ESD kaons            

//          maximum decay angle at a given mother momentum
	   //Double_t maxDecAngKmu=f1->Eval(motherMfromKink.Mag(),0.,0.,0.);
	   Double_t maxDecAngKmu=f1->Eval(track->P()          ,0.,0.,0.);
	   Double_t maxDecAngpimu=f2->Eval(    track->P(),       0.,0.,0.);
         if( (kinkAngle<2.)  ) continue;

           
      //  BG  ?????==============
              if ( TMath::Abs(vposKink[2]) >  225. ) continue ;
              if ( TMath::Abs(vposKink[2]) <  0.5 ) continue ;
//                if(( vposKink[2] >0. )&& (vposKink[2]< 5.)  ) continue; 
//
            fKinkKaonBg->Fill(motherPt);     
 fAngMomPi->Fill( track->P(),           kinkAngle); 
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
           fQtInvM -> Fill ( invariantMassKmu,  qT);
           fInvMuNuAll->Fill(invariantMassKmu);
   //     Float_t ptKink=TMath::Sqrt(p1XM*p1XM + p1YM*p1YM);
  
      if( ( kink->GetR()> 120 ) && ( kink->GetR() < 210 )  )  {
      if (qT>0.12)  fAngMomKC->Fill(track->P(), kinkAngle); 
          if ( qT>0.12) fM1kaon->Fill(invariantMassKmu);
             if ( qT > 0.12) 
         fRadiusNcl->Fill( (kink->GetR()) ,(track->GetTPCclusters(0)  ) ) ;
  }    
        //  fPosiKinKBg->Fill( vposKink[2], kink->GetR() );
          if(  ( tpcNCl<30) ) continue;
         if( ( ( track->GetTPCclusters(0) ) / (kink->GetR() ) ) > 0.63 ) continue;
            if( ( ( track->GetTPCclusters(0) ) / (kink->GetR() ) ) < 0.20 ) continue;

//
               fHistPtKPDG->Fill(track->Pt());  // ALL KAONS (pdg) inside ESD  kink sample
//if((kinkAngle>maxDecAngpimu)&&(qT>0.04)&&(qT<0.30)&&((kink->GetR()>120.)&&(kink->GetR()<200.))&&(TMath::Abs(trackEta)<0.9)&&(invariantMassKmu<0.6)) {
 // if((kinkAngle>maxDecAngpimu)&&(qT>0.04)&&(qT<0.30)&&((kink->GetR()>=120.)&&(kink->GetR()<=200.))&&(TMath::Abs(etaMother)<0.9)&&(invariantMassKmu<0.6)){
    if((kinkAngle>maxDecAngpimu)&&(qT>0.12)&&(qT<0.30)&&((kink->GetR()>=120.)&&(kink->GetR()<=210.))&&(TMath::Abs(rapiditK)<0.7)&&(invariantMassKmu<0.6)){
 // 16/10  if((kinkAngle>maxDecAngpimu)&&(qT>0.04)&&(qT<0.30)&&((kink->GetR()>=120.)&&(kink->GetR()<=210.))&&(TMath::Abs(rapiditK)<0.7)&&(invariantMassKmu<0.6)){
//  if((kinkAngle>maxDecAngpimu)&&(qT>0.04)&&(qT<0.30)&&((kink->GetR()>=133.)&&(kink->GetR()<=179.))&&(TMath::Abs(rapiditK)<0.5)&&(invariantMassKmu<0.6)){      // STAR 

        if( (kinkAngle<maxDecAngpimu*1.2) ) continue; 
                 if ( (kinkAngle>maxDecAngKmu*.98) && ( track->P() >1.2 )) continue;  ///5/5/2010

/*
*/

         fTPCSgnlPa->Fill( trMoment ,(track->GetTPCsignal()  ) ) ;
                          if ( nsigma               > 3.5) continue;
                                     fHistPtKaon->Fill(track->Pt());   //all PID kink-kaon
              if(tpcSign >0.)        fHistPtKaoP->Fill( track->Pt()         ) ;   //all PID kink-kaon
               if ( tpcSign <0.)    fHistPtKaoN->Fill(  track->Pt()        ) ;   //all PID kink-kaon
                frad->Fill(kink->GetR());  // kink 
               fradLK->Fill(lifeKink    );  // kink 
             fHistEtaK->Fill(trackEta);
            frapiKESD ->Fill(rapiditK);  //  rapidityof kaons 
          fPosiKinKBg->Fill( vposKink[2], kink->GetR() );

                     Float_t signPt= tpcSign*trackPt;
                  fSignPtNcl->Fill( signPt  ,   tpcNCl   );   ///  28/4/2010
                  fSignPtEta->Fill( signPt  , rapiditK  );
                  fEtaNcl->Fill( rapiditK, tpcNCl    );
                  fSignPt->Fill( signPt );
                  fChi2NclTPC->Fill( (track->GetTPCchi2() ) ,  tpcNCl );
         fRatChi2Ncl-> Fill (  (track->GetTPCchi2()/track->GetTPCclusters(0)  )) ;
    fdcatoVxXY->Fill(dcaToVertexXYpos);
                  fKinkMothDau->Fill(track->P(),daughterMKink.Mag());
         fTPCSgnlP->Fill(track->P(), (track->GetTPCsignal()  ) ) ;
		        flifetime->Fill(( lifeKink*.493667   )  /track->P()   ) ;
             fKinkKaon->Fill(track->Pt());        

             fKinKRbn->Fill(track->Pt());        
             fptKMC   ->Fill(  track->Pt()    );        
              fradPtRpDt->Fill( kink->GetR(), 1./track->Pt(), rapiditK); 
               fAngMomK->Fill(    track->P(),        kinkAngle); 
       fPosiKinkK->Fill( vposKink[0], vposKink[1]  );
          fPosiKinKXZ->Fill( vposKink[2], vposKink[0]  );
          fPosiKinKYZ->Fill( vposKink[2], vposKink[1]  );

        }  //  kink selection 
                  

	}  //End Kink Information    
  

  } //track loop 

 //     fMultiplMC->Fill(nESDTracK );

  PostData(1, fListOfHistos);

}      

//________________________________________________________________________
void AliAnalysisKinkESDat::Terminate(Option_t *) 
{
  // Draw result to the screen 
  // Called once at the end of the query

}

//____________________________________________________________________//


const AliESDVertex* AliAnalysisKinkESDat::GetEventVertex(const AliESDEvent* esd) const
  // Get the vertex from the ESD and returns it if the vertex is valid
  
{
  // Get the vertex 
  
// 24/3  const AliESDVertex* vertex = esd->GetPrimaryVertex();
   const AliESDVertex* vertex = esd->GetPrimaryVertexTracks();

  // if((vertex->GetStatus()==kTRUE)&&(vertex->GetNContributors()>2)) return vertex;
  if((vertex->GetStatus()==kTRUE)) return vertex;
  else
  { 
     vertex = esd->GetPrimaryVertexSPD();
      if((vertex->GetStatus()==kTRUE)&&(vertex->GetNContributors()>0)) return vertex;
//     if((vertex->GetStatus()==kTRUE)) return vertex;
     else
     return 0;
  }
}