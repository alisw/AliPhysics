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




void testsimexam (Int_t nevent=1, const char *config="testconfig.C")
{
  //
  // Simple macro to run aliroot in a batch mode
  //
  cerr<<" __________________________________________________________________ "<<endl;
  cerr<< " " <<endl;
  cerr<<"           MESS ==> Beginning of the simulation examination."<<endl;
  cerr<<" __________________________________________________________________ "<<endl;
  // Definition of the alarm bounds
     const Float_t maxAlaHitsM = 12.79  ;  // total multiplicity
     const Float_t maxAlaTotEn = 19.34  ;  // total energy
     const Float_t maxAlaTotEnB = 18.35 ;  // per block multiplicity
     const Float_t maxAlaHitsMB = 11.1  ;  // per block energy
   
       // boolean which will test if there is an alarm
  
   Int_t boolala; 
   boolala = 0;
   // define the array in which the events that have not reach th EMCA will be put.
  
   

   TStopwatch timer;
   timer.Start();

   // Get the number of events generated in the simulation

  AliPHOSGetter * gime = AliPHOSGetter::GetInstance("testPHOS.root") ; 
  gime->Event(0,"Q");
  Int_t maxevent = gime->MaxEvent() ; 
 
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
          boolala = 1; 
          cerr<<" _____________________________________________________________ "<<endl;
          cerr<< " " <<endl;
          cerr << "             MESS ==> Examination detected an error in "<<namemul << endl ;    
          cerr<<" _____________________________________________________________ "<<endl;
	}
        if (ratioenb[i]>maxAlaTotEnB) {
          boolala = 1;
          cerr<<" _____________________________________________________________ "<<endl;
          cerr<< " " <<endl;
          cerr << "             MESS ==> Examination detected an error in "<<namen << endl ;
          cerr<<" _____________________________________________________________ "<<endl;
	}
	    
     }
 

  timer.Stop();
  timer.Print();

      
  if (ratiohm>maxAlaHitsM){
     boolala = 1;
     cerr<<" _____________________________________________________________ "<<endl;
     cerr<< " " <<endl;
     cerr << "             MESS ==> Examination detected an error in HitsM." << endl ;
     cerr<<" _____________________________________________________________ "<<endl;
  }

  if (ratioet>maxAlaTotEn){
     boolala = 1;
     cerr<<" _____________________________________________________________ "<<endl;
     cerr<< " " <<endl;
     cerr << "             MESS ==> Examination detected an error in TotEn." << endl ;
     cerr<<" _____________________________________________________________ "<<endl;
  }

 
  // Condition that will launch the general loop that builds the histograms in order to allow a further analysis.
  
    if ( boolala == 1 ) {
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
      //  gSystem->Exec("uuencode $ALICE_ROOT/PHOS/testPHOS.root testPHOS.root | mail -s 'PHOS INSTALLATION ERROR  ' schutz@in2p3.fr"); 
   
    }
  cerr<<" _____________________________________________________________ "<<endl;
  cerr<< " " <<endl;
  cerr << "             MESS ==> Examination ended successfully. " << endl ;
  cerr<<" _____________________________________________________________ "<<endl;   
}
