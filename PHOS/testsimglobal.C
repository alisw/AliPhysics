#include "AliPHOSGetter.h"
#include "TH2.h"
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "AliPHOSHit.h"
#include "TF1.h"
#include "TFormula.h"
#include "TFolder.h"
#include "TStopwatch.h"
#include "TObjArray.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSDigit.h"
#include "AliPHOSSDigitizer.h"
#include "AliPHOSDigitizer.h"
#include "AliPHOSClusterizer.h"
#include "AliPHOSClusterizerv1.h"
#include "AliPHOSTrackSegmentMaker.h"
#include "AliPHOSTrackSegmentMakerv1.h"
#include "AliPHOSPID.h"
#include "AliPHOSPIDv1.h"


void testsimglobal(Int_t nevent = 1, const char *config="testconfig.C")
{
  //
  // Simple macro to run aliroot in a batch mode
  //
  cerr<<" __________________________________________________________________ "<<endl;
  cerr<< " " <<endl;
  cerr<<"             MESS ==> Beginning of the simulation examination."<<endl;
  cerr<<" __________________________________________________________________ "<<endl;
  // Definition of the alarm bounds
     const Float_t maxAlaHitsM = 12.79  ;  // total multiplicity
     const Float_t maxAlaTotEn = 19.34  ;  // total energy
     const Float_t maxAlaTotEnB = 18.35 ;  // per block multiplicity
     const Float_t maxAlaHitsMB = 11.1  ;  // per block energy
   
       // boolean which will test if there is an alarm
   Int_t boolerror = 0;
   
   // define the array in which the events that have not reach th EMCA will be put.
  
   

   TStopwatch timer;
   timer.Start();

   // Get the number of events generated in the simulation

  AliPHOSGetter * gime = AliPHOSGetter::GetInstance("testPHOS.root") ; 
  Int_t maxevent = gime->MaxEvent() ; 
  gime->Event(0,"Q");
 
  // Examine the alarms

   TObjArray * alahm = (TObjArray*)(gime->Alarms()->FindObject("HitsM"));
   Float_t ratiohm = 100.0*(Float_t)alahm->GetEntries()/(Float_t)maxevent ;  
 
  
   TObjArray * alaet = (TObjArray*)(gime->Alarms()->FindObject("TotEn"));
   Float_t ratioet = 100.0*(Float_t)alaet->GetEntries()/(Float_t)maxevent ; 
  

   // Define the alarms per block and examine them
   char namemul[80], namen[80];
   TObjArray* alahmb[5];
   TObjArray* alaenb[5];
   Float_t ratiohmb[5], ratioenb[5];
   for (Int_t i = 0 ; i < 5 ; i++)
     {
        sprintf(namemul,"HitsMB%d",i+1);
        sprintf(namen,"TotEnB%d",i+1);
        alahmb[i] = (TObjArray*) (gime->Alarms()->FindObject(namemul));
	ratiohmb[i] = 100.0*(Float_t)alahmb[i]->GetEntries()/(Float_t)maxevent;        
        alaenb[i] = (TObjArray*)(gime->Alarms()->FindObject(namen));
        ratioenb[i] = 100.0*(Float_t)alaenb[i]->GetEntries()/(Float_t)maxevent;      
        if (ratiohmb[i]>maxAlaHitsMB){
          boolerror = 1; 
          cerr<<" _____________________________________________________________ "<<endl;
          cerr<< " " <<endl;
          cerr << "             MESS ==> Examination detected an error in "<<namemul << endl ;    
          cerr<<" _____________________________________________________________ "<<endl;
	}
        if (ratioenb[i]>maxAlaTotEnB) {
          boolerror = 1;
          cerr<<" _____________________________________________________________ "<<endl;
          cerr<< " " <<endl;
          cerr << "             MESS ==> Examination detected an error in "<<namen << endl ;
          cerr<<" _____________________________________________________________ "<<endl;
	}
	    
     }
 

  timer.Stop();
  timer.Print();

      
  if (ratiohm>maxAlaHitsM){
     boolerror = 1;
     cerr<<" _____________________________________________________________ "<<endl;
     cerr<< " " <<endl;
     cerr << "             MESS ==> Examination detected an error in HitsM." << endl ;
     cerr<<" _____________________________________________________________ "<<endl;
  }

  if (ratioet>maxAlaTotEn){
     boolerror = 1;
     cerr<<" _____________________________________________________________ "<<endl;
     cerr<< " " <<endl;
     cerr << "             MESS ==> Examination detected an error in TotEn." << endl ;
     cerr<<" _____________________________________________________________ "<<endl;
  }

 
  // Condition that will launch the general loop that builds the histograms in order to allow a further analysis.
  
    if ( boolerror == 1 ) {
  cerr<<" _____________________________________________________________ "<<endl;
  cerr<< " " <<endl;
  cerr << "             MESS ==> Examination sets up the file that will be sent to PHOS director (30s). " << endl ;
  cerr<<" _____________________________________________________________ "<<endl;    
      Int_t index = 0 ; 
      Int_t nhits = 0 ; 
    
      TH1F * his = new TH1F("Total Multiplicity", "Total Multiplicity in PHOS", 200, 0., 200.) ;
      TH1F * hisnrg = new TH1F("Total Energy", "Total energy in PHOS", 200, 0., 12.) ;
    
     // name and define the different histograms per block involved in the analysis
      TClonesArray hisba("TH1F") ;
      TClonesArray hisbanrg("TH1F") ;  
      for (Int_t i = 0 ; i < 6 ; i++)
       {
         char name[80], title[80] ; 
         sprintf(name,"multiplicity for block %d",i) ; 
         sprintf(title,"multiplicity per blocks, block %d",i) ;
         new(hisba[i]) TH1F(name, title, 100, 0., 100.) ; 
      

         char namenrg[80], titlenrg[80] ; 
         sprintf(namenrg,"total energy for block %d",i) ; 
         sprintf(titlenrg,"total energy per block, block %d",i) ; 	    
         new(hisbanrg[i]) TH1F(namenrg, titlenrg, 200, 0., 12.) ; 
       }
      // define the global histograms, the hit pointer and give the means to get the actual block reached by the photon
      AliPHOSHit * hit;
      TH1F * hist ; 
      TH1F * histnrg;
      TH2F * hbiz = new TH2F ("hbiz","hbiz",200.,0.,200.,200.,0.,12.);
      AliPHOSGeometry * geom = AliPHOSGeometry::GetInstance() ; 


      // the very general loop
      for (index = 0 ; index < maxevent ; index++)
       {
	 //get the number of the event
         gime->Event(index) ; 
         // get the number of cells reached during this event and fill the total multiplicity histogram
         Int_t n = gime->Hits()->GetEntries() ;
         nhits += n ; 
         his->Fill(n) ;  
         // Get the data per block      
         TClonesArray * hita = (TClonesArray *) gime -> Hits();
         TIter next(hita);
         Float_t Et = 0.;
         Int_t id = 0, block = 0;
         Int_t nhit[6], rid[4];
         Float_t etblock[6]; 
         for ( Int_t i = 0; i < 6 ; i++) { 
	   nhit[i] = 0 ; 
	   etblock[i] = 0 ;
	 }
     
	 while ( (hit = (AliPHOSHit *) next()) ) {
	         id = hit->GetId();
		 if  (geom->IsInEMC(id) ) { 
		   Et += ( hit -> GetEnergy());
		   geom->AbsToRelNumbering(id,rid) ;
		   block = rid[0];
		   nhit[block]++  ;
		   etblock[block] +=  ( hit -> GetEnergy());
		 }
                }


	 //Fill all the histograms but total multiplicity, already done
	 hist = static_cast<TH1F*>(hisba.At(block)) ; 
	 hist->Fill(nhit[block]) ; 
	 histnrg = static_cast<TH1F*>(hisbanrg.At(block)) ;
	 histnrg->Fill(etblock[block]);
      	 hisnrg -> Fill(Et);
         hbiz->Fill(n,Et);	         
       }   
        
      nhits /= maxevent ; 
      cerr << "av = " << nhits << endl ;
      TFile * file = gROOT -> GetFile("testPHOS.root");
      file -> Write();
      his->Draw() ;
      hisnrg->Draw() ;
      hbiz->Draw();

      //Put the histograms in the root file
      for (Int_t i = 0 ; i < 6 ; i++)	{ 
         hist = static_cast<TH1F*>(hisba.At(i)) ;
         histnrg = static_cast<TH1F*>(hisbanrg.At(i));
     
        
	 // cout << hist << endl << i << endl ; 
	 hist->Draw();
	 histnrg->Draw(); 
	 
   
       }
      
      file->Write() ;
      hisba.Delete() ;
      hisbanrg.Delete() ; 
      file->Close();
        gSystem->Exec("uuencode $ALICE_ROOT/PHOS/testPHOS.root testPHOS.root | mail -s 'PHOS INSTALLATION ERROR  ' almajor@caramail.com"); 
   
    }
  cerr<<" _____________________________________________________________ "<<endl;
  cerr<< " " <<endl;
  cerr << "             MESS ==> Examination ended successfully. " << endl ;
  cerr<<" _____________________________________________________________ "<<endl;   



  const Float_t maxSDigits = 62.89 ;
  const Float_t widSDigits = TMath::Sqrt(maxSDigits) ;
  const Float_t maxDigits = 3489.41 ;
  const Float_t widDigits = TMath::Sqrt(maxDigits) ;
  const Float_t maxRecPoints = 222.83 ;
  const Float_t widRecPoints = TMath::Sqrt(maxRecPoints) ;
  const Float_t maxTrackSegments = 1 ;
  const Float_t maxRecParticles = 1 ;
  TString reconame = "test suite" ;
  if (boolerror == 1){
  cerr<<" ___________________________________________________________________ "<<endl;
  cerr<<" "<<endl;
  cerr<<"             MESS ==> Please find the error and try the whole simulation again. "<<endl;
  cerr<<" ___________________________________________________________________ "<<endl;
  }

  else {

  cerr<<" ___________________________________________________________________ "<<endl;
  cerr<<" "<<endl;
  cerr<<"             MESS ==> Beginning of the PHOS reconstruction. "<<endl;
  cerr<<" ___________________________________________________________________ "<<endl;

  //SDigits process
   AliPHOSSDigitizer *sd = new AliPHOSSDigitizer("testPHOS.root",reconame.Data());
   AliPHOSGetter * gime = AliPHOSGetter::GetInstance() ;
   sd->ExecuteTask("deb"); 
   Float_t nSDigits =  (Float_t) (gime->SDigitizer()->GetSDigitsInRun()) / gime->MaxEvent(); 
   if ( nSDigits < maxSDigits-widSDigits || nSDigits > maxSDigits+widSDigits ) {
    boolerror = 1 ;  
    cerr<<"__________________________________________________________________"<<endl;
    cerr<<" "<<endl;
    cerr<<"             MESS ==> Error detected in the SDigits process. Sending error file to PHOS director."<<endl;
    cerr<<"__________________________________________________________________"<<endl;
    gSystem->Exec("uuencode $ALICE_ROOT/PHOS/testPHOS.root testPHOS.root | mail -s 'PHOS INSTALLATION ERROR' almajor@caramail.com");
    }

  cerr<<"__________________________________________________________________"<<endl;
  cerr<<" "<<endl;
  cerr<<"             MESS ==> SDigits process ended successfully."<<endl;
  cerr<<"__________________________________________________________________"<<endl;

  //Digits process
     if (boolerror == 0){
        AliPHOSDigitizer *d = new AliPHOSDigitizer("testPHOS.root",reconame.Data());
        d->ExecuteTask("deb"); 
         Float_t nDigits = (Float_t)(gime->Digitizer()->GetDigitsInRun()) / gime->MaxEvent();
  
              if ( nDigits < maxDigits-widDigits || nDigits > maxDigits+widDigits ) {
                  boolerror =1 ;
                  cerr<<"__________________________________________________________________"<<endl;
                  cerr<<" "<<endl;
                  cerr<<"             MESS ==> Error detected in the Digits process. Sending error file to PHOS director."<<endl;
                  cerr<<"__________________________________________________________________"<<endl;
   
    gSystem->Exec("uuencode $ALICE_ROOT/PHOS/testPHOS.root testPHOS.root | mail -s 'PHOS INSTALLATION ERROR' almajor@caramail.com");
	      }
  
    cerr<<"__________________________________________________________________"<<endl;
    cerr<<" "<<endl;
    cerr<<"             MESS ==> Digits process ended successfully."<<endl;
    cerr<<"__________________________________________________________________"<<endl;
  }
  
//RecPoints process
 if (boolerror == 0){
 AliPHOSClusterizer * cluster = new  AliPHOSClusterizerv1("testPHOS.root", reconame.Data());
 
  cluster->ExecuteTask("deb");
  Float_t nRecPoints =  (Float_t) (gime->Clusterizer(reconame.Data())->GetRecPointsInRun()) / gime->MaxEvent();
 
   if ( nRecPoints < maxRecPoints-widRecPoints || nRecPoints > maxRecPoints+widRecPoints ) {
    boolerror = 1;
    cerr<<"__________________________________________________________________"<<endl;
    cerr<<" "<<endl;
    cerr<<"             MESS ==> Error detected in the Clusterizing process. Sending error file to PHOS director."<<endl;
    cerr<<"__________________________________________________________________"<<endl;
  
      gSystem->Exec("uuencode $ALICE_ROOT/PHOS/testPHOS.root testPHOS.root | mail -s 'PHOS INSTALLATION ERROR' almajor@caramail.com");
    }
   
  cerr<<"__________________________________________________________________"<<endl;
  cerr<<" "<<endl;
  cerr<<"             MESS ==> Cluster process ended successfully."<<endl;
  cerr<<"__________________________________________________________________"<<endl;
 }

//TrackSegments process
 if(boolerror == 0){ 
  AliPHOSTrackSegmentMaker * tracks = new AliPHOSTrackSegmentMakerv1("testPHOS.root",reconame.Data());
  tracks->ExecuteTask("deb");
  Float_t nTrackSegments =  (Float_t) (gime->TrackSegmentMaker(reconame.Data())->GetTrackSegmentsInRun()) / gime->MaxEvent();
 
   if ( nTrackSegments < maxTrackSegments-0.25 || nTrackSegments > maxTrackSegments+0.25 ) {
    boolerror = 1;
    cerr<<"__________________________________________________________________"<<endl;
    cerr<<" "<<endl;
    cerr<<"             MESS ==> Error detected in the TrackSegments process. Sending error file to PHOS director."<<endl;
    cerr<<"__________________________________________________________________"<<endl;
    
      gSystem->Exec("uuencode $ALICE_ROOT/PHOS/testPHOS.root testPHOS.root | mail -s 'PHOS INSTALLATION ERROR' almajor@caramail.com");
    }
   
  cerr<<"__________________________________________________________________"<<endl;
  cerr<<" "<<endl;
  cerr<<"             MESS ==> TrackSegments process ended successfully."<<endl;
  cerr<<"__________________________________________________________________"<<endl;
 }
//RecParticles process
 if(boolerror == 0){
  AliPHOSPID * pid = new AliPHOSPIDv1("testPHOS.root",reconame.Data()); pid->ExecuteTask("deb");
 
  Float_t nRecParticles =  (Float_t) (gime->PID(reconame.Data())->GetRecParticlesInRun())/gime->MaxEvent();
 
 
   if ( nRecParticles < maxRecParticles-0.25 || nRecParticles > maxRecParticles+0.25 ) {
    boolerror = 1;
    cerr<<"__________________________________________________________________"<<endl;
    cerr<<" "<<endl;
    cerr<<"             MESS ==> Error detected in the RecParticles process. Sending error file to PHOS director.Stop reconstruction."<<endl;
    cerr<<"__________________________________________________________________"<<endl;
    
    gSystem->Exec("uuencode $ALICE_ROOT/PHOS/testPHOS.root testPHOS.root | mail -s 'PHOS INSTALLATION ERROR' almajor@caramail.com");
    }
 
  cerr<<"__________________________________________________________________"<<endl;
  cerr<<" "<<endl;
  cerr<<"             MESS ==> RecParticles process ended successfully."<<endl;
  cerr<<"__________________________________________________________________"<<endl;
 }
 if(boolerror == 0){
  cerr<<"_____________________________________________________________________"<<endl;
  cerr<<" "<<endl;
  cerr<<"             MESS ==> reconstruction ended successfully."<<endl;
  cerr<<"_____________________________________________________________________"<<endl;
 }
 else {
 cerr<<"_____________________________________________________________________"<<endl;
  cerr<<" "<<endl;
  cerr<<"             MESS ==> Due to an error, reconstruction has not been completed.Please try again when the error has been found."<<endl;
  cerr<<"_____________________________________________________________________"<<endl;
  }
 }
}
