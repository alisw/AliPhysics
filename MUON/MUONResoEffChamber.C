#if !defined(__CINT__) || defined(__MAKECINT__)

// ROOT includes
#include "TBranch.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TParticle.h"
#include "TFile.h"

#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TMath.h"

#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraph2D.h"
#include "TStyle.h"
#include "TFitter.h"
#include "TRandom.h"

// STEER includes
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliHeader.h"
#include "AliLoader.h"
#include "AliTracker.h"
#include "AliStack.h"
#include "AliMagFMaps.h"


// MUON includes
#include "AliMUON.h"
#include "AliMUONData.h"
#include "AliMUONConstants.h"

#include "AliMUONHit.h"
#include "AliMUONHitForRec.h"
#include "AliMUONRawCluster.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackParam.h"
#include "AliMUONTrackExtrap.h"

#include "AliMpVSegmentation.h"
#include "AliMpIntPair.h"
#include "AliMpDEManager.h"
#endif


//Macro to calculate the resolution and the efficiency of chamber(s) chosen in function
//of Phi (angle on the chamber between the X axis and the straight line created by the
//center of the chamber and the impact of particles on this chamber, the azimuthal angle)
//and Theta (the polar angle), or in function of ThetaI (angle of incidence of particles
//on the chamber)




  const Double_t kInvPi = 1./3.14159265;


//Chamber number:
  Int_t chamberNbr;
//Number of events:
  Int_t nEvents, iEvent;
//Number of tracks:
  Int_t nTracks, iTrack;
//Number of hits:
  Int_t nHits,iHit;
//Chamber(s) chosen
  Int_t firstChamber, lastChamber;

 
  AliMUONTrack      * track  ;
  AliMUONHitForRec  * hit = 0;
  AliMUONTrackParam * trackParam = 0;
 
  TClonesArray      * tracks  ;
  TClonesArray      * trackParams;
  TClonesArray      * hits  ;






/*****************************************************************************************************************/
/*****************************************************************************************************************/
                                               /*EFFICIENCY*/

void efficiency( Int_t event2Check=0, char * filename="galice.root" )
{
  cout<<"\nChamber(s) chosen;\nFirst chamber:";
  cin>>firstChamber;
  cout<<"Last chamber:";
  cin>>lastChamber;
  cout<<"\n\n";

//Creating Run Loader and openning file containing Hits
//--------------------------------------------------------------------------------------------------------------
  AliRunLoader * RunLoader = AliRunLoader::Open(filename,"MUONFolder","READ");
  if (RunLoader == 0x0) {printf(">>> Error : Error Opening %s file \n",filename);return;}


//Getting MUONLoader
//--------------------------------------------------------------------------------------------------------------
  AliLoader * MUONLoader = RunLoader -> GetLoader("MUONLoader");
  MUONLoader -> LoadTracks("READ");
  MUONLoader -> LoadRecPoints("READ");


//Creating a MUON data container
//--------------------------------------------------------------------------------------------------------------
  AliMUONData muondata(MUONLoader,"MUON","MUON");


  nEvents = (Int_t) RunLoader -> GetNumberOfEvents();

    //Loop on events
    Int_t trackNb = 0;
    Int_t chamber[10] = {0};
    Int_t detEltNew, detElt;
    Int_t detEltOld = 0;

    for ( iEvent = 0; iEvent < nEvents; ++iEvent )
    {
      if ( event2Check!=0 ) 
	  iEvent=event2Check;
      printf("\r>>> Event %d ",iEvent);
      RunLoader -> GetEvent(iEvent);
      //Addressing
      muondata.SetTreeAddress("RT");
      muondata.GetRecTracks();
      tracks = muondata.RecTracks();

      //Loop on track
      nTracks = (Int_t) tracks -> GetEntriesFast();
      for ( iTrack = 0; iTrack < nTracks; ++iTrack )
      { 
	track     = (AliMUONTrack*) tracks -> At(iTrack);
	hits = track -> GetHitForRecAtHit(); 
	detEltOld = 0; 
	//Loop on hit
	nHits = (Int_t) hits -> GetEntriesFast();

	for ( iHit = 0; iHit < nHits; ++iHit )
	{ 
	  hit        = (AliMUONHitForRec*) hits -> At(iHit);
	  chamberNbr =                     hit  -> GetChamberNumber();
	  detElt     =                     hit  -> GetDetElemId();
	  detEltNew  = int(detElt/100);
	  if( chamberNbr >= firstChamber-1 ) {
	    if( chamberNbr < lastChamber ) {
	      if( detEltNew == detEltOld ) 
		  continue;
	      else {
		chamber[chamberNbr] = chamber[chamberNbr] + 1; 
		detEltOld = detEltNew;
	      }
	    }
	  }
	}	    
	//End loop on hit

      }
      //End loop on track 
      muondata.ResetRecTracks();
      if (event2Check != 0)
	  iEvent = nEvents;
      trackNb = trackNb + nTracks;
    }
    //End loop on event
//--------------------------------------------------------------------------------------------------------------

    cout<<"\n\n\n";
    for (Int_t i = firstChamber-1; i < lastChamber; ++i )
    { 
      printf ("\nChamber %d has responded %d times on %d tracks", i+1, chamber[i], trackNb);
      if (trackNb == 0) 
	  cout<<"\nEfficiency = ? (IS UNKNOWN)\n";
      else {
	Double_t eff = chamber[i]*100./trackNb; cout<<"\nEfficiency = "<<eff<<" %\n";
      }
    }
    cout<<"\n\n\n";
    MUONLoader->UnloadTracks();
}





/*****************************************************************************************************************/
/*****************************************************************************************************************/
                               /*EFFICIENCY IN FUNCTION OF THETA AND PHI*/

void efficiencyThetaPhi( Int_t event2Check=0, char * filename="galice.root" )
{
  cout<<"\nChamber(s) chosen;\nFirst chamber:";
  cin>>firstChamber;
  cout<<"Last chamber:";
  cin>>lastChamber;
  cout<<"\n\n";

  Int_t eff    [54] = {0};
  Int_t trackN [54] = {0};
  Int_t chamber;
  Int_t oldChamber;

//Creating Run Loader and openning file containing Hits
//--------------------------------------------------------------------------------------------------------------
  AliRunLoader * RunLoader = AliRunLoader::Open(filename,"MUONFolder","READ");
  if (RunLoader == 0x0) {printf(">>> Error : Error Opening %s file \n",filename);return;}


//Getting MUONLoader
//--------------------------------------------------------------------------------------------------------------
  AliLoader * MUONLoader = RunLoader -> GetLoader("MUONLoader");
  MUONLoader -> LoadTracks("READ");
  MUONLoader -> LoadRecPoints("READ");


//--------------------------------------------------------------------------------------------------------------
//Set mag field; waiting for mag field in CDB 
  printf("Loading field map...\n");
  if (!AliTracker::GetFieldMap()) {
  AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 1, 1., 10., AliMagFMaps::k5kG);
  AliTracker::SetFieldMap(field, kFALSE);}


//--------------------------------------------------------------------------------------------------------------
//Set Field Map for track extrapolation
  AliMUONTrackExtrap::SetField(AliTracker::GetFieldMap());


//Creating a MUON data container
//--------------------------------------------------------------------------------------------------------------
  AliMUONData muondata(MUONLoader,"MUON","MUON");


  nEvents = (Int_t) RunLoader -> GetNumberOfEvents();

    //Loop on events
    for ( iEvent = 0; iEvent < nEvents; ++iEvent )
    { 
      if ( event2Check!=0 )
	  iEvent=event2Check;
      printf("\r>>> Event %d ",iEvent);
      RunLoader->GetEvent(iEvent);
      
      //Addressing
      muondata.SetTreeAddress("RT");     
      muondata.GetRecTracks();
      tracks = muondata.RecTracks();

      //Loop on track
      nTracks = (Int_t) tracks -> GetEntriesFast();
      for ( iTrack = 0; iTrack < nTracks; ++iTrack )
      {
	track       = (AliMUONTrack*) tracks ->At(iTrack)          ;
	trackParams =                 track  ->GetTrackParamAtHit();
	hits        =                 track  ->GetHitForRecAtHit() ;
	chamber     = firstChamber-1;
	oldChamber  = -1;
	Int_t k     =  1;
  
	//Loop on hits
	nHits = (Int_t) hits->GetEntriesFast();
	for ( iHit = 0; iHit<nHits; ++iHit )
	{ 
	  trackParam = (AliMUONTrackParam*) trackParams -> At(iHit);
	  hit        = (AliMUONHitForRec* ) hits        -> At(iHit);
	  chamberNbr = hit -> GetChamberNumber();
	       
	  if ( chamberNbr == oldChamber ) 
	      continue;
	  else { 
	    oldChamber = chamberNbr;
	    if ( chamberNbr > chamber - k ) {
	      if ( chamber < lastChamber ) {
		if ( chamberNbr == chamber ) {
		  //Positions
		  Double_t traX, traY, traZ;
		  Double_t hitX, hitY, hitZ;
		  Double_t aveX, aveY      ;

		  //Angle (Phi)
		  Double_t phi   = 0.;
		  Double_t theta = 0.;
		  Int_t    iPhi   = 0 ;
		  Int_t    iTheta = 0 ;					     
						  
		  traX = trackParam -> GetNonBendingCoor();
		  traY = trackParam -> GetBendingCoor()   ;
		  traZ = trackParam -> GetZ()             ;
						  
		  hitX = hit        -> GetNonBendingCoor();
		  hitY = hit        -> GetBendingCoor()   ;
		  hitZ = hit        -> GetZ()             ;
						  
		  aveX = 1./2.*(traX + hitX);
		  aveY = 1./2.*(traY + hitY);
		
		  //The calculation of phi:
		  phi   = 180. * kInvPi * (TMath::ATan2( aveY, aveX ));

		  //The calculation of theta, theta is in fact 180 - Theta(The polar angle):
		  theta = 180. * kInvPi * (TMath::ATan2( sqrt(aveX*aveX+aveY*aveY), -hitZ ));

		  if ( phi < 0.) phi = 360 - abs(phi);
		  iPhi = int( phi/72. );
		  iTheta = int( theta );
		  if( theta > 10 ) iTheta = 9;
		  if( theta < 1  ) iTheta = 1;
					       
		  eff    [9*iPhi+iTheta-1] = eff    [9*iPhi+iTheta-1] + 1;
		  trackN [9*iPhi+iTheta-1] = trackN [9*iPhi+iTheta-1] + 1;  
		  chamber                  = chamber + 1;
		}   

		else {
		  //Positions
		  Double_t chamberZpos;
		  Double_t exXpos, exYpos;
						  
		  //Angles
		  Double_t phi   = 0.;
		  Double_t theta = 0.;
		  Int_t    iPhi   = 0 ;
		  Int_t    iTheta = 0 ;
						  
		  chamberZpos = AliMUONConstants::DefaultChamberZ(chamber);
		  AliMUONTrackExtrap::ExtrapToZ(trackParam,chamberZpos);
		  exXpos      = (Double_t) trackParam->GetNonBendingCoor();
		  exYpos      = (Double_t) trackParam->GetBendingCoor();
		  
		  //The calculation of phi:
		  phi   = 180. * kInvPi * (TMath::ATan2( exYpos, exXpos ));				  

		  //The calculation of theta, theta is in fact 180 - Theta(The polar angle):
		  theta = 180. * kInvPi * (TMath::ATan2( sqrt(exXpos*exXpos+exYpos*exYpos), - chamberZpos ));

		  if ( phi < 0.) phi = 360 - abs(phi);
		  iPhi = int( phi/72. );
		  iTheta = int( theta );
		  if( theta > 10 ) iTheta = 9;
		  if( theta < 1  ) iTheta = 1;

		  eff    [9*iPhi+iTheta-1] = eff    [9*iPhi+iTheta-1] + 0;
		  trackN [9*iPhi+iTheta-1] = trackN [9*iPhi+iTheta-1] + 1;
		  chamber                  = chamber + 1;
		  iHit                     = iHit - 1;
		}
	      }
	    }
	  }
						   
	  if ( iHit == nHits-1 ) {
	    if    ( chamber == lastChamber )
		continue;
	    else {
	      oldChamber = -1;
	      k          =  5;
	      iHit       = iHit-1;
	    }
	  }                                   
	} 
	//End Loop on hits

      } 
      //End Loop on tracks

      muondata.ResetRecTracks();
      if ( event2Check!=0 )
	  iEvent=nEvents;
    }
    //End Loop on events


    TCanvas * c1 = new TCanvas();
    TGraph2D* effPhiTheta = new TGraph2D();
    Double_t  efficiency = 0;

    if ( firstChamber == lastChamber )
    {
      for ( Int_t ph = 0; ph < 5; ++ph )
      {
	for ( Int_t th = 1; th < 10; ++th )
	{
	  Int_t i = 9*ph+th-1;
	  cout<<"\nFor Phi = ["<<ph*72<<","<<ph*72+72<<"] and Theta = ["<<179-th<<","<<180-th<<"]:";
	  cout<<"\nThe chamber "<<firstChamber<<" has responded "<<eff[i]<<" times on "<<trackN[i]<<" tracks\n";

	  Double_t e = eff   [i] ;
	  Double_t n = trackN[i] ;
	  Double_t p = ph*72.+36.;
	  Double_t t = th*1. +0.5;
	
	  if ( trackN[i] == 0 ) {
	    efficiency = 0.;
	    cout<<"Efficiency = ? % ( IS UNKNOWN )\n";
	  }
	  else {
	    efficiency = 100.*e/n;
	    cout<<"Efficiency = "<<efficiency<<" %\n";
	  }

	  effPhiTheta -> SetPoint( i, p, t, efficiency);
	}
      }
    }

    else{
      for ( Int_t ph = 0; ph < 5; ++ph )
      {
	for ( Int_t th = 1; th < 10; ++th )
	{
	  Int_t i = 9*ph+th-1;
	  cout<<"\nFor Phi = ["<<ph*72<<","<<ph*72+72<<"] and Theta = ["<<179-th<<","<<180-th<<"]:";
	  cout<<"\nThe chambers "<<firstChamber<<" to "<<lastChamber<<" have responded "<<eff[i]<<" times on "<<trackN[i]<<" tracks\n";

	  Double_t e = eff   [i] ;
	  Double_t n = trackN[i] ;
	  Double_t p = ph*72.+36.;
	  Double_t t = th*1. +0.5;
	
	  if ( trackN[i] == 0 ) {
	    efficiency = 0.;
	    cout<<"Efficiency = ? % ( IS UNKNOWN )\n";
	  }
	  else {
	    efficiency = 100.*e/n;
	    cout<<"Efficiency = "<<efficiency<<" %\n";
	  }
	  
	  effPhiTheta -> SetPoint( i, p, t, efficiency);
	}
      }
    }

    gStyle->SetPalette(1);
    effPhiTheta -> Draw("surf1");

    cout<<"\n\n\n";
    MUONLoader->UnloadTracks();
    c1->Update();

}







/*****************************************************************************************************************/
/*****************************************************************************************************************/
                          /*EFFICIENCY IN FUNCTION OF THETA OF INCIDENCE*/


void efficiencyThetaI( Int_t event2Check=0, char * filename="galice.root" )
{
  cout<<"\nChamber(s) chosen;\nFirst chamber:";
  cin>>firstChamber;
  cout<<"Last chamber:";
  cin>>lastChamber;
  cout<<"\n\n";

  Int_t eff    [12] = {0};
  Int_t trackN [12] = {0};
  Int_t chamber;
  Int_t oldChamber;


//Creating Run Loader and openning file containing Hits
//--------------------------------------------------------------------------------------------------------------
  AliRunLoader * RunLoader = AliRunLoader::Open(filename,"MUONFolder","READ");
  if (RunLoader == 0x0) {printf(">>> Error : Error Opening %s file \n",filename);return;}


//Getting MUONLoader
//--------------------------------------------------------------------------------------------------------------
  AliLoader * MUONLoader = RunLoader -> GetLoader("MUONLoader");
  MUONLoader -> LoadTracks("READ");
  MUONLoader -> LoadRecPoints("READ");


//--------------------------------------------------------------------------------------------------------------
//Set mag field; waiting for mag field in CDB 
  printf("Loading field map...\n");
  if (!AliTracker::GetFieldMap()) {
  AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 1, 1., 10., AliMagFMaps::k5kG);
  AliTracker::SetFieldMap(field, kFALSE);}


//--------------------------------------------------------------------------------------------------------------
//Set Field Map for track extrapolation
  AliMUONTrackExtrap::SetField(AliTracker::GetFieldMap());


//Creating a MUON data container
//--------------------------------------------------------------------------------------------------------------
  AliMUONData muondata(MUONLoader,"MUON","MUON");


  nEvents = (Int_t) RunLoader -> GetNumberOfEvents();

    //Loop on events
    for ( iEvent = 0; iEvent < nEvents; ++iEvent )
    { 
      if ( event2Check!=0 )
	  iEvent=event2Check;
      printf("\r>>> Event %d ",iEvent);
      RunLoader->GetEvent(iEvent);
      
      //Addressing
      muondata.SetTreeAddress("RT");     
      muondata.GetRecTracks();
      tracks = muondata.RecTracks();

      //Loop on track
      nTracks = (Int_t) tracks -> GetEntriesFast();
      for ( iTrack = 0; iTrack < nTracks; ++iTrack )
      { 
	track       = (AliMUONTrack*) tracks ->At(iTrack)          ;
	trackParams =                 track  ->GetTrackParamAtHit();
	hits        =                 track  ->GetHitForRecAtHit() ;
	chamber     = firstChamber - 1;
	oldChamber  = -1;
	Int_t k     =  1;

	//Loop on hits
	nHits = (Int_t) hits -> GetEntriesFast();
	for ( iHit = 0; iHit < nHits; ++iHit )
	{ 
	  trackParam = (AliMUONTrackParam*) trackParams -> At(iHit);
	  hit        = (AliMUONHitForRec*)  hits        -> At(iHit);
	  chamberNbr = hit -> GetChamberNumber();

	  if ( chamberNbr == oldChamber )
	      continue;
	  else {
	    oldChamber = chamberNbr;
	    if ( chamberNbr > chamber - k ) {
	      if ( chamber < lastChamber ) {
		if ( chamberNbr == chamber ) {
		  //Momentum
		  Double_t Px, Py, Pz, Pr;

		  //Angle
		  Double_t theta;
		  Int_t    iTheta;
						  
		  Px = trackParam -> Px()   ;
		  Py = trackParam -> Py()   ;
		  Pz = trackParam -> Pz()   ;
		  Pr = sqrt( Px*Px + Py*Py );

		  //The calculation of theta, the angle of incidence:						   
		  theta = 180. * kInvPi * (TMath::ATan( - Pz / Pr));

		  if      ( theta < 79 ) iTheta = 11;
		  else if ( theta < 90 ) iTheta = int( theta - 79.);
		  else                   iTheta = 11;

		  eff    [iTheta] = eff    [iTheta] + 1;
		  trackN [iTheta] = trackN [iTheta] + 1;
		  chamber         = chamber + 1;
		}

		else {
		  //Positions
		  Double_t chamberZpos;

		  //Momentum
		  Double_t Px, Py, Pz, Pr;

		  //Angles
		  Double_t theta = 0.;
		  Int_t    iTheta = 0 ;
						   
		  chamberZpos = AliMUONConstants::DefaultChamberZ(chamber);
		  AliMUONTrackExtrap::ExtrapToZ(trackParam,chamberZpos);
						   
		  Px = trackParam -> Px()   ;
		  Py = trackParam -> Py()   ;
		  Pz = trackParam -> Pz()   ;
		  Pr = sqrt( Px*Px + Py*Py );
						  
		  //The calculation of thetaI, the angle of incidence:
		  theta = 180. * kInvPi * (TMath::ATan( - Pz / Pr));

		  if      ( theta < 79 ) iTheta = 11;
		  else if ( theta < 90 ) iTheta = int( theta - 79.);
		  else                   iTheta = 11;
						   
		  eff    [iTheta] = eff    [iTheta] + 0;
		  trackN [iTheta] = trackN [iTheta] + 1;
		  chamber         = chamber + 1;
		  iHit            = iHit - 1;
		}
	      }
	    }
	  }
						   
	  if ( iHit == nHits-1 ) {
	    if ( chamber == lastChamber )
		continue;
	    else {
	      oldChamber = -1;
	      k          =  5;
	      iHit       = iHit-1;
	    }
	  }
	}
	//End loop on hits

      } 
      //End Loop on tracks

      muondata.ResetRecTracks();
      if ( event2Check!=0 )
	  iEvent=nEvents;
    }
    //End Loop on events

    Double_t t [11];
    Double_t efficiency [11];
    Int_t i = 11;

    Int_t th;
    TGraph * effTheta = new TGraph ();

    if ( firstChamber == lastChamber ) {
      for ( th = 0; th < 11; ++th )
      {
	cout<<"\nFor Theta (angle of incidence) = ["<<th+79<<","<<th+80<<"]:\n";
	cout<<"The chamber "<<firstChamber<<" has responded "<<eff [th]<<" times on "<<trackN [th]<<" tracks\n";
  
	t [th]  = th + 79.5  ;
	Double_t e = eff    [th];
	Double_t n = trackN [th];
 
	if ( n == 0. ) {
	  efficiency [th] = 0.;
	  cout<<"Efficiency = ? % (IS UNKNOWN)          \n";
	}
	else {
	  efficiency [th] = 100.*e/n;
	  cout<<"Efficiency = "<<efficiency [th]<<" %\n";
	}
      }
    }

    else{
      for ( th = 0; th < 11; ++th )
      {
	cout<<"\nFor Theta (angle of incidence) = ["<<th+79<<","<<th+80<<"]:\n";
	cout<<"The chambers "<<firstChamber<<" to "<<lastChamber<<" have responded "<<eff [th]<<" times on "<<trackN [th]<<" tracks\n";
  
	t [th]  = th + 79.5  ;
	Double_t e = eff    [th];
	Double_t n = trackN [th];
 
	if ( n == 0. ) {
	  efficiency [th] = 0.;
	  cout<<"Efficiency = ? % (IS UNKNOWN)          \n";
	}
	else {
	  efficiency [th] = 100.*e/n;
	  cout<<"Efficiency = "<<efficiency [th]<<" %\n";
	}
      }
    }

    TCanvas * c1 = new TCanvas ();
    effTheta = new TGraph( i, t, efficiency );

    effTheta -> Draw("ALP");
 
    MUONLoader->UnloadTracks();
    c1->Update();
}






/*****************************************************************************************************************/
/*****************************************************************************************************************/
                                                /*RESOLUTION*/

void resolution( Int_t event2Check=0, char * filename="galice.root" )
{
  cout<<"\nChamber(s) chosen;\nFirst chamber:";
  cin>>firstChamber;
  cout<<"Last chamber:";
  cin>>lastChamber;
  cout<<"\n\n";

  TH1F * hDeltax;
  TH1F * hDeltay;
  TH2  * hDelta3D;

  hDeltax    = new TH1F ("hDeltax " , "No Interest", 100, -0.3, 0.3);
  hDeltay    = new TH1F ("hDeltay " , "No Interest", 500, -0.1, 0.1);
  hDelta3D   = new TH2F ("hDelta3D" , "No Interest", 100, -0.3, 0.3, 500, -0.1, 0.1);


//Creating Run Loader and openning file containing Hits
//--------------------------------------------------------------------------------------------------------------
  AliRunLoader * RunLoader = AliRunLoader::Open(filename,"MUONFolder","READ");
  if (RunLoader == 0x0) {printf(">>> Error : Error Opening %s file \n",filename);return;}


//Getting MUONLoader
//--------------------------------------------------------------------------------------------------------------
  AliLoader * MUONLoader = RunLoader -> GetLoader("MUONLoader");
  MUONLoader -> LoadTracks("READ");
  MUONLoader -> LoadRecPoints("READ");


//Creating a MUON data container
//--------------------------------------------------------------------------------------------------------------
  AliMUONData muondata(MUONLoader,"MUON","MUON");


  nEvents = (Int_t) RunLoader -> GetNumberOfEvents();
    
    //Loop on events
     for ( iEvent = 0; iEvent < nEvents; ++iEvent )
    { 
      if (event2Check!=0)
	  iEvent=event2Check;
      printf("\r>>> Event %d ",iEvent);
      RunLoader->GetEvent(iEvent);
      
      //Addressing
      muondata.SetTreeAddress("RT");     
      muondata.GetRecTracks();
      tracks = muondata.RecTracks();

      //Loop on track
      nTracks = (Int_t) tracks -> GetEntriesFast();
      for ( iTrack = 0; iTrack < nTracks; ++iTrack )
      {
	track       = (AliMUONTrack*) tracks ->At(iTrack)          ;
	trackParams =                 track  ->GetTrackParamAtHit();
	hits        =                 track  ->GetHitForRecAtHit() ;
       
	//Loop on hits
	nHits = (Int_t) hits -> GetEntriesFast();
	for ( iHit = 0; iHit < nHits; ++iHit )
	{ 
	  trackParam = (AliMUONTrackParam*) trackParams -> At(iHit);
	  hit        = (AliMUONHitForRec*)  hits        -> At(iHit);
	  chamberNbr = hit -> GetChamberNumber();
	  if ( chamberNbr >= firstChamber-1 ) {
	    if ( chamberNbr < lastChamber ) {
	      //Positions
	      Double_t traX, traY;
	      Double_t hitX, hitY;
						     
	      //Resolution
	      Double_t deltaX, deltaY;
						     
	      traX = trackParam -> GetNonBendingCoor();
	      traY = trackParam -> GetBendingCoor();
	      hitX = hit        -> GetNonBendingCoor();
	      hitY = hit        -> GetBendingCoor();
						     
	      deltaX = traX - hitX;
	      deltaY = traY - hitY;
						     
	      hDeltax  -> Fill (deltaX);
	      hDeltay  -> Fill (deltaY);
	      hDelta3D -> Fill (deltaX, deltaY);
	    }
	  }
	}
	//End loop on hits
      }
      //End loop on tracks
      muondata.ResetRecTracks();
      if ( event2Check!=0 )
	  iEvent=nEvents;
    }
    //End loop on events

    char hXtitle[80];
    char hYtitle[80];
    char h3title[80];

    if ( firstChamber == lastChamber ) {
      sprintf(hXtitle,"DeltaX Chamber %d",firstChamber);
      sprintf(hYtitle,"DeltaY Chamber %d",lastChamber);
      sprintf(h3title,"DeltaX and Delta Y Chamber %d",lastChamber);
    }
    else{
      sprintf(hXtitle,"DeltaX Chamber %d to %d",firstChamber,lastChamber);
      sprintf(hYtitle,"DeltaY Chamber %d to %d",firstChamber,lastChamber);
      sprintf(h3title,"DeltaX and Delta Y Chamber %d to %d",firstChamber,lastChamber);
    }
  

    hDeltax  -> SetTitle (hXtitle);
    hDeltay  -> SetTitle (hYtitle);
    hDelta3D -> SetTitle (h3title);

    Double_t rmsX = hDeltax -> GetRMS     ();
    Double_t errX = hDeltax -> GetRMSError();

    TF1 *fY = new TF1("fY","gaus",-0.3, 0.3);
    hDeltay -> Fit("fY","R","E");
    Double_t sigY = fY -> GetParameter(2);
    Double_t errY = fY -> GetParError (2);

    if ( firstChamber == lastChamber ) {
      cout<<"\nThe resolution (Non Bending plane) of the chamber "<<firstChamber<<" = "<<rmsX<<" cm +-"<<errX;
      cout<<"\nThe resolution (Bending plane) of the chamber "<<firstChamber<<" = "<<sigY<<" cm +- "<<errY;
    }
    else {
      cout<<"\nThe resolution (Non Bending plane) of the chambers "<<firstChamber<<" to "<<lastChamber<<" = "<<rmsX<<" cm +- "<<errX;
      cout<<"\nThe resolution (Bending plane) of the chambers "<<firstChamber<<" to "<<lastChamber<<" = "<<sigY<<" cm +- "<<errY;
    }

    cout<<"\n\n\n";
    MUONLoader->UnloadTracks();

}







/*****************************************************************************************************************/
/*****************************************************************************************************************/
                                       /*RESOLUTION IN FUNCTION OF PHI*/

void resolutionPhi( Int_t event2Check=0, char * filename="galice.root" )
{
  cout<<"\nChamber(s) chosen;\nFirst chamber:";
  cin>>firstChamber;
  cout<<"Last chamber:";
  cin>>lastChamber;
  cout<<"\n\n";

  TH1F * hDeltaX[5];
  TH1F * hDeltaY[5];

  hDeltaX[0] = new TH1F ("hDeltaX0" , "No Interest", 100, -0.3, 0.3);
  hDeltaX[1] = new TH1F ("hDeltaX1" , "No Interest", 100, -0.3, 0.3);
  hDeltaX[2] = new TH1F ("hDeltaX2" , "No Interest", 100, -0.3, 0.3);
  hDeltaX[3] = new TH1F ("hDeltaX3" , "No Interest", 100, -0.3, 0.3);
  hDeltaX[4] = new TH1F ("hDeltaX4" , "No Interest", 100, -0.3, 0.3);
  
  hDeltaY[0] = new TH1F ("hDeltaY0" , "No Interest", 500, -0.1, 0.1);
  hDeltaY[1] = new TH1F ("hDeltaY1" , "No Interest", 500, -0.1, 0.1);
  hDeltaY[2] = new TH1F ("hDeltaY2" , "No Interest", 500, -0.1, 0.1);
  hDeltaY[3] = new TH1F ("hDeltaY3" , "No Interest", 500, -0.1, 0.1);
  hDeltaY[4] = new TH1F ("hDeltaY4" , "No Interest", 500, -0.1, 0.1);


//Creating Run Loader and openning file containing Hits
//--------------------------------------------------------------------------------------------------------------
  AliRunLoader * RunLoader = AliRunLoader::Open(filename,"MUONFolder","READ");
  if (RunLoader == 0x0) {printf(">>> Error : Error Opening %s file \n",filename);return;}


//Getting MUONLoader
//--------------------------------------------------------------------------------------------------------------
  AliLoader * MUONLoader = RunLoader -> GetLoader("MUONLoader");
  MUONLoader -> LoadTracks("READ");
  MUONLoader -> LoadRecPoints("READ");


//Creating a MUON data container
//--------------------------------------------------------------------------------------------------------------
  AliMUONData muondata(MUONLoader,"MUON","MUON");


  nEvents = (Int_t) RunLoader -> GetNumberOfEvents();

    //Loop on events
    for ( iEvent = 0; iEvent < nEvents; ++iEvent )
    {
      if ( event2Check!=0 )
	  iEvent=event2Check;
      printf("\r>>> Event %d ",iEvent);
      RunLoader->GetEvent(iEvent);
      
      //Addressing
      muondata.SetTreeAddress("RT");     
      muondata.GetRecTracks();
      tracks = muondata.RecTracks();

      //Loop on track
      nTracks = (Int_t) tracks -> GetEntriesFast();
      for ( iTrack = 0; iTrack < nTracks; ++iTrack )
      {
	track       = (AliMUONTrack*) tracks ->At(iTrack)          ;
	trackParams =                 track  ->GetTrackParamAtHit();
	hits        =                 track  ->GetHitForRecAtHit() ;
       
	//Loop on hits
	nHits = (Int_t) hits -> GetEntriesFast();
	for ( iHit = 0; iHit < nHits; ++iHit )
	{ 
	  trackParam = (AliMUONTrackParam*) trackParams -> At(iHit);
	  hit        = (AliMUONHitForRec* ) hits        -> At(iHit);
	  chamberNbr = hit -> GetChamberNumber();
	  if ( chamberNbr >= firstChamber -1 ) {
	    if ( chamberNbr < lastChamber     ) {
	      //Positions
	      Double_t traX, traY;
	      Double_t hitX, hitY;
	      Double_t aveY, aveX;

	      //Angles
	      Double_t phi;
	      Int_t    iPhi;

	      //Resolution					    				     
	      Double_t deltaX;
	      Double_t deltaY;

	      traX = trackParam -> GetNonBendingCoor();
	      traY = trackParam -> GetBendingCoor()   ;
						      						      
	      hitX = hit        -> GetNonBendingCoor();
	      hitY = hit        -> GetBendingCoor()   ;
						      						      
	      aveX = 1./2.*(traX + hitX);
	      aveY = 1./2.*(traY + hitY);

	      //The calculation of phi:	
	      phi = 180. * kInvPi * (TMath::ATan2( aveY, aveX ));					      

	      if ( phi < 0.) phi = 360 - abs(phi);
	      iPhi = int( phi/72. );
						      
	      deltaX = traX - hitX;
	      deltaY = traY - hitY;

	      hDeltaX [iPhi] -> Fill( deltaX );
	      hDeltaY [iPhi] -> Fill( deltaY );
	    }
	  }
	    
	}
	//End loop on hits

      }
      //End loop on tracks
      muondata.ResetRecTracks();
      if ( event2Check!=0 )
	  iEvent=nEvents;

    }
    //End loop on events


    Int_t iPhi;
    Int_t iPhiMax = 5;
    Int_t phiMin, phiMax;
 
    Float_t phi[5];
    Float_t sigmaY[5];
    Float_t errSY [5];
    Float_t rmsX  [5];
    Float_t errSX [5];
    Float_t errPh [5];

    for ( iPhi = 0; iPhi < iPhiMax; ++iPhi )
    {
      char hXtitle[80];
      char hYtitle[80];

      phiMin = 72*iPhi     ;
      phiMax = 72*iPhi + 72;
	 
      TF1 *fY = new TF1("fY","gaus",-0.3, 0.3);
       
      if ( firstChamber == lastChamber ) {
	sprintf(hXtitle,"DeltaX (phi=[%d, %d]) Chamber %d",phiMin,phiMax,firstChamber);
	sprintf(hYtitle,"DeltaY (phi=[%d, %d]) Chamber %d",phiMin,phiMax,lastChamber);
      }
      else{
	sprintf(hXtitle,"DeltaX (phi=[%d, %d]) Chamber %d to %d",phiMin,phiMax,firstChamber,lastChamber);
	sprintf(hYtitle,"DeltaY (phi=[%d, %d]) Chamber %d to %d",phiMin,phiMax,firstChamber,lastChamber);
      }

      hDeltaY [iPhi] -> SetTitle(hYtitle);
      hDeltaX [iPhi] -> SetTitle(hXtitle);

      hDeltaY [iPhi]      -> Fit("fY","R","E");
      sigmaY  [iPhi] = fY -> GetParameter(2)  ;
      errSY   [iPhi] = fY -> GetParError(2)   ;

      rmsX    [iPhi] = hDeltaX [iPhi] -> GetRMS();
      errSX   [iPhi] = hDeltaX [iPhi] -> GetRMSError(); 

      phi     [iPhi] = 72*iPhi + 36 ;
      errPh   [iPhi] = 36;
    }

    //--------------------------------------------------------------------------------------------------------------
    //For plotting resolution in function of the angle of incidence

    TCanvas * c1 = new TCanvas("c1","",200,10,700,500);
    c1-> Divide(1,2);
    c1->cd(1);

    TGraphErrors * GraphX = new TGraphErrors( iPhiMax, phi, rmsX, errPh, errSX);
    GraphX->SetTitle("Resolution en X (Phi)");
    GraphX->Draw("ALP");

    c1->cd(2);
    TGraphErrors * GraphY = new TGraphErrors( iPhiMax, phi, sigmaY, errPh, errSY);
    GraphY->SetTitle("Resolution en Y (Phi)");
    GraphY->Draw("ALP");

    cout<<"\n\n\n";

    MUONLoader->UnloadTracks();

}







/*****************************************************************************************************************/
/*****************************************************************************************************************/
                                       /*RESOLUTION IN FUNCTION OF THETA*/

void resolutionTheta( Int_t event2Check=0, char * filename="galice.root" )
{
  cout<<"\nChamber(s) chosen;\nFirst chamber:";
  cin>>firstChamber;
  cout<<"Last chamber:";
  cin>>lastChamber;
  cout<<"\n\n";

  TH1F * hDeltaX[9];
  TH1F * hDeltaY[9];

  hDeltaX[0] = new TH1F ("hDeltaX0" , "No Interest", 100, -0.3, 0.3);
  hDeltaX[1] = new TH1F ("hDeltaX1" , "No Interest", 100, -0.3, 0.3);
  hDeltaX[2] = new TH1F ("hDeltaX2" , "No Interest", 100, -0.3, 0.3);
  hDeltaX[3] = new TH1F ("hDeltaX3" , "No Interest", 100, -0.3, 0.3);
  hDeltaX[4] = new TH1F ("hDeltaX4" , "No Interest", 100, -0.3, 0.3);
  hDeltaX[5] = new TH1F ("hDeltaX5" , "No Interest", 100, -0.3, 0.3);
  hDeltaX[6] = new TH1F ("hDeltaX6" , "No Interest", 100, -0.3, 0.3);
  hDeltaX[7] = new TH1F ("hDeltaX7" , "No Interest", 100, -0.3, 0.3);
  hDeltaX[8] = new TH1F ("hDeltaX8" , "No Interest", 100, -0.3, 0.3);
  
  hDeltaY[0] = new TH1F ("hDeltaY0" , "No Interest", 500, -0.1, 0.1);
  hDeltaY[1] = new TH1F ("hDeltaY1" , "No Interest", 500, -0.1, 0.1);
  hDeltaY[2] = new TH1F ("hDeltaY2" , "No Interest", 500, -0.1, 0.1);
  hDeltaY[3] = new TH1F ("hDeltaY3" , "No Interest", 500, -0.1, 0.1);
  hDeltaY[4] = new TH1F ("hDeltaY4" , "No Interest", 500, -0.1, 0.1);
  hDeltaY[5] = new TH1F ("hDeltaY5" , "No Interest", 500, -0.1, 0.1);
  hDeltaY[6] = new TH1F ("hDeltaY6" , "No Interest", 500, -0.1, 0.1);
  hDeltaY[7] = new TH1F ("hDeltaY7" , "No Interest", 500, -0.1, 0.1);
  hDeltaY[8] = new TH1F ("hDeltaY8" , "No Interest", 500, -0.1, 0.1);


//Creating Run Loader and openning file containing Hits
//--------------------------------------------------------------------------------------------------------------
  AliRunLoader * RunLoader = AliRunLoader::Open(filename,"MUONFolder","READ");
  if (RunLoader == 0x0) {printf(">>> Error : Error Opening %s file \n",filename);return;}


//Getting MUONLoader
//--------------------------------------------------------------------------------------------------------------
  AliLoader * MUONLoader = RunLoader -> GetLoader("MUONLoader");
  MUONLoader -> LoadTracks("READ");
  MUONLoader -> LoadRecPoints("READ");


//Creating a MUON data container
//--------------------------------------------------------------------------------------------------------------
  AliMUONData muondata(MUONLoader,"MUON","MUON");


  nEvents = (Int_t) RunLoader -> GetNumberOfEvents();

    //Loop on events
    for ( iEvent = 0; iEvent < nEvents; ++iEvent )
    {
      if ( event2Check!=0 )
	  iEvent=event2Check;
      printf("\r>>> Event %d ",iEvent);
      RunLoader->GetEvent(iEvent);
      
      //Addressing
      muondata.SetTreeAddress("RT");     
      muondata.GetRecTracks();
      tracks = muondata.RecTracks();

      //Loop on track
      nTracks = (Int_t) tracks -> GetEntriesFast();
      for ( iTrack = 0; iTrack < nTracks; ++iTrack )
      { 
	track       = (AliMUONTrack*) tracks ->At(iTrack)          ;
	trackParams =                 track  ->GetTrackParamAtHit();
	hits        =                 track  ->GetHitForRecAtHit() ;
       
	//Loop on hits
	nHits = (Int_t) hits -> GetEntriesFast();
	for ( iHit = 0; iHit < nHits; ++iHit )
	{ 
	  trackParam = (AliMUONTrackParam*) trackParams -> At(iHit);
	  hit        = (AliMUONHitForRec* ) hits        -> At(iHit);
	  chamberNbr = hit -> GetChamberNumber();
	  if ( chamberNbr >= firstChamber -1 ) {
	    if ( chamberNbr < lastChamber ) {
	      //Positions
	      Double_t traX, traY;
	      Double_t hitX, hitY, hitZ;
	      Double_t aveY, aveX;
						      
	      //Angles
	      Double_t theta;
	      Int_t    iTheta;

	      //Resolution					    				     
	      Double_t deltaX;
	      Double_t deltaY;

	      traX = trackParam -> GetNonBendingCoor();
	      traY = trackParam -> GetBendingCoor()   ;
						      						      
	      hitX = hit        -> GetNonBendingCoor();
	      hitY = hit        -> GetBendingCoor()   ;
	      hitZ = hit        -> GetZ();
						      
	      aveX = 1./2.*(traX + hitX);
	      aveY = 1./2.*(traY + hitY);
						      
	      //The calculation of theta, theta is in fact 180 - Theta(The polar angle):
	      theta = 180. * kInvPi * (TMath::ATan2( sqrt(aveX*aveX+aveY*aveY), -hitZ ));

	      iTheta = int( theta );
	      if( theta > 10 ) iTheta = 9;
	      if( theta < 1  ) iTheta = 1; 
						      
	      deltaX = traX - hitX;
	      deltaY = traY - hitY;

	      hDeltaX [iTheta-1] -> Fill( deltaX );
	      hDeltaY [iTheta-1] -> Fill( deltaY );
	    }
	  }
						    
	}
	//End loop on hits

      }
      //End loop on tracks
      muondata.ResetRecTracks();
      if ( event2Check!=0 )
	  iEvent=nEvents;

    }
    //End loop on events


    Int_t iTheta;
    Int_t iThetaMax = 9;
    Int_t thetaMin, thetaMax;
 
    Float_t theta [9];
    Float_t sigmaY[9];
    Float_t errSY [9];
    Float_t rmsX  [9];
    Float_t errSX [9];
    Float_t ErrTh [9];

    for ( iTheta = 0; iTheta < iThetaMax; ++iTheta )
    {
      char hXtitle[80];
      char hYtitle[80];

      //To find the polar angle
      thetaMin = 178 - iTheta;
      thetaMax = 179 - iTheta;
	 
      TF1 *fY = new TF1("fY","gaus",-0.3, 0.3);
       
      if ( firstChamber == lastChamber ) {
	sprintf(hXtitle,"DeltaX (theta=[%d, %d]) Chamber %d",thetaMin,thetaMax,firstChamber);
	sprintf(hYtitle,"DeltaY (theta=[%d, %d]) Chamber %d",thetaMin,thetaMax,lastChamber);
      }
      else{
	sprintf(hXtitle,"DeltaX (theta=[%d, %d]) Chamber %d to %d",thetaMin,thetaMax,firstChamber,lastChamber);
	sprintf(hYtitle,"DeltaY (theta=[%d, %d]) Chamber %d to %d",thetaMin,thetaMax,firstChamber,lastChamber);
      }

      hDeltaY [iTheta] -> SetTitle(hYtitle);
      hDeltaX [iTheta] -> SetTitle(hXtitle);

      hDeltaY [iTheta]      -> Fit("fY","R","E");
      sigmaY  [iTheta] = fY -> GetParameter(2)  ;
      errSY   [iTheta] = fY -> GetParError(2)   ;

      rmsX    [iTheta] = hDeltaX [iTheta] -> GetRMS();
      errSX   [iTheta] = hDeltaX [iTheta] -> GetRMSError(); 

      theta   [iTheta] = 178.5 - iTheta;
      ErrTh   [iTheta] = 0.5;
    }

    //--------------------------------------------------------------------------------------------------------------
    //For plotting resolution in function of the angle of incidence

    TCanvas * c1 = new TCanvas("c1","",200,10,700,500);
    c1-> Divide(1,2);
    c1->cd(1);

    TGraphErrors * GraphX = new TGraphErrors( iThetaMax, theta, rmsX, ErrTh, errSX);
    GraphX->SetTitle("Resolution en X (Theta)");
    GraphX->Draw("ALP");

    c1->cd(2);
    TGraphErrors * GraphY = new TGraphErrors( iThetaMax, theta, sigmaY, ErrTh, errSY);
    GraphY->SetTitle("Resolution en Y (Theta)");
    GraphY->Draw("ALP");

    cout<<"\n\n\n";
    MUONLoader->UnloadTracks();

}





/*****************************************************************************************************************/
/*****************************************************************************************************************/
                              /*RESOLUTION IN FUNCTION OF THETA OF INCIDENCE*/

void resolutionThetaI( Int_t event2Check=0, char * filename="galice.root" )
{
  cout<<"\nChamber(s) chosen;\nFirst chamber:";
  cin>>firstChamber;
  cout<<"Last chamber:";
  cin>>lastChamber;
  cout<<"\n\n";

  TH1F * hDeltaX[11];
  TH1F * hDeltaY[11];

  hDeltaX[0] = new TH1F ("hDeltaX0" , "No Interest", 100, -0.3, 0.3);
  hDeltaX[1] = new TH1F ("hDeltaX1" , "No Interest", 100, -0.3, 0.3);
  hDeltaX[2] = new TH1F ("hDeltaX2" , "No Interest", 100, -0.3, 0.3);
  hDeltaX[3] = new TH1F ("hDeltaX3" , "No Interest", 100, -0.3, 0.3);
  hDeltaX[4] = new TH1F ("hDeltaX4" , "No Interest", 100, -0.3, 0.3);
  hDeltaX[5] = new TH1F ("hDeltaX5" , "No Interest", 100, -0.3, 0.3);
  hDeltaX[6] = new TH1F ("hDeltaX6" , "No Interest", 100, -0.3, 0.3);
  hDeltaX[7] = new TH1F ("hDeltaX7" , "No Interest", 100, -0.3, 0.3);
  hDeltaX[8] = new TH1F ("hDeltaX8" , "No Interest", 100, -0.3, 0.3);
  hDeltaX[9] = new TH1F ("hDeltaX9" , "No Interest", 100, -0.3, 0.3);
  hDeltaX[10]= new TH1F ("hDeltaX10", "No Interest", 100, -0.3, 0.3);
  
  hDeltaY[0] = new TH1F ("hDeltaY0" , "No Interest", 500, -0.1, 0.1);
  hDeltaY[1] = new TH1F ("hDeltaY1" , "No Interest", 500, -0.1, 0.1);
  hDeltaY[2] = new TH1F ("hDeltaY2" , "No Interest", 500, -0.1, 0.1);
  hDeltaY[3] = new TH1F ("hDeltaY3" , "No Interest", 500, -0.1, 0.1);
  hDeltaY[4] = new TH1F ("hDeltaY4" , "No Interest", 500, -0.1, 0.1);
  hDeltaY[5] = new TH1F ("hDeltaY5" , "No Interest", 500, -0.1, 0.1);
  hDeltaY[6] = new TH1F ("hDeltaY6" , "No Interest", 500, -0.1, 0.1);
  hDeltaY[7] = new TH1F ("hDeltaY7" , "No Interest", 500, -0.1, 0.1);
  hDeltaY[8] = new TH1F ("hDeltaY8" , "No Interest", 500, -0.1, 0.1);
  hDeltaY[9] = new TH1F ("hDeltaY9" , "No Interest", 500, -0.1, 0.1);
  hDeltaY[10]= new TH1F ("hDeltaY10", "No Interest", 500, -0.1, 0.1);


//Creating Run Loader and openning file containing Hits
//--------------------------------------------------------------------------------------------------------------
  AliRunLoader * RunLoader = AliRunLoader::Open(filename,"MUONFolder","READ");
  if (RunLoader == 0x0) {printf(">>> Error : Error Opening %s file \n",filename);return;}


//Getting MUONLoader
//--------------------------------------------------------------------------------------------------------------
  AliLoader * MUONLoader = RunLoader -> GetLoader("MUONLoader");
  MUONLoader -> LoadTracks("READ");
  MUONLoader -> LoadRecPoints("READ");


//Creating a MUON data container
//--------------------------------------------------------------------------------------------------------------
  AliMUONData muondata(MUONLoader,"MUON","MUON");


  nEvents = (Int_t) RunLoader -> GetNumberOfEvents();

    //Loop on events
    for ( iEvent = 0; iEvent < nEvents; ++iEvent )
    { 
      if ( event2Check!=0 )
	  iEvent=event2Check;
      printf("\r>>> Event %d ",iEvent);
      RunLoader->GetEvent(iEvent);
      
      //Addressing
      muondata.SetTreeAddress("RT");     
      muondata.GetRecTracks();
      tracks = muondata.RecTracks();

      //Loop on track
      nTracks = (Int_t) tracks -> GetEntriesFast();
      for ( iTrack = 0; iTrack < nTracks; ++iTrack )
      {
	track       = (AliMUONTrack*) tracks ->At(iTrack)          ;
	trackParams =                 track  ->GetTrackParamAtHit();
	hits        =                 track  ->GetHitForRecAtHit() ;
       
	//Loop on hits
	nHits = (Int_t) hits -> GetEntriesFast();
	for ( iHit = 0; iHit < nHits; ++iHit )
	{
	  trackParam = (AliMUONTrackParam*) trackParams -> At(iHit);
	  hit        = (AliMUONHitForRec* ) hits        -> At(iHit);
	  chamberNbr = hit -> GetChamberNumber();
	  if ( chamberNbr >= firstChamber - 1 ) {
	    if ( chamberNbr < lastChamber ) {
	      //Momentum
	      Double_t Px, Py, Pz, Pr;

	      //Positions
	      Double_t traX, traY;
	      Double_t hitX, hitY;

	      //Resolution
	      Double_t deltaX;
	      Double_t deltaY;

	      //Angle
	      Double_t theta;
	      Int_t    iTheta;

	      Px = trackParam -> Px()   ;
	      Py = trackParam -> Py()   ;
	      Pz = trackParam -> Pz()   ;
	      Pr = sqrt( Px*Px + Py*Py );

	      traX = trackParam -> GetNonBendingCoor();
	      traY = trackParam -> GetBendingCoor();

	      hitX = hit        -> GetNonBendingCoor();
	      hitY = hit     	-> GetBendingCoor();					    

	      //The calculation of theta, the angle of incidence:
	      theta = 180. * kInvPi * (TMath::ATan( - Pz / Pr))
;
	      if ( theta < 79 ) iTheta = 0;
	      else iTheta = int( theta - 79.);

	      deltaX = traX - hitX;
	      deltaY = traY - hitY;

	      hDeltaX [iTheta] -> Fill (deltaX);
	      hDeltaY [iTheta] -> Fill (deltaY);
	    }
	  }
						    
	}
	//End loop on hits

      }
      //End loop on tracks
      muondata.ResetRecTracks();
      if ( event2Check!=0 )
	  iEvent=nEvents;

    }
    //End loop on events


    Int_t iTheta;
    Int_t iThetaMax = 11;
    Int_t thetaMin, thetaMax;
 
    Float_t theta [11];
    Float_t sigmaY[11];
    Float_t errSY [11];
    Float_t rmsX  [11];
    Float_t errSX [11];
    Float_t errTh [11];

    for ( iTheta = 0; iTheta < iThetaMax; ++iTheta )
    {
      char hXtitle[80];
      char hYtitle[80];

      thetaMin = iTheta + 79;
      thetaMax = iTheta + 80;
	 
      TF1 *fY = new TF1("fY","gaus",-0.3, 0.3);
       
      if ( firstChamber == lastChamber ) {
	sprintf(hXtitle,"DeltaX (theta=[%d, %d]) Chamber %d",thetaMin,thetaMax,firstChamber);
	sprintf(hYtitle,"DeltaY (theta=[%d, %d]) Chamber %d",thetaMin,thetaMax,lastChamber);
      }
      else{
	sprintf(hXtitle,"DeltaX (theta=[%d, %d]) Chamber %d to %d",thetaMin,thetaMax,firstChamber,lastChamber);
	sprintf(hYtitle,"DeltaY (theta=[%d, %d]) Chamber %d to %d",thetaMin,thetaMax,firstChamber,lastChamber);
      }

      hDeltaY [iTheta] -> SetTitle(hYtitle);
      hDeltaX [iTheta] -> SetTitle(hXtitle);

      hDeltaY [iTheta]      -> Fit("fY","R","E");
      sigmaY  [iTheta] = fY -> GetParameter(2)  ;
      errSY   [iTheta] = fY -> GetParError(2)   ;

      rmsX    [iTheta] = hDeltaX [iTheta] -> GetRMS();
      errSX   [iTheta] = hDeltaX [iTheta] -> GetRMSError(); 

      theta   [iTheta] = iTheta + 79.5 ;
      errTh   [iTheta] = 0.5;
    }

    //--------------------------------------------------------------------------------------------------------------
    //For plotting resolution in function of the angle of incidence

    TCanvas * c1 = new TCanvas("c1","",200,10,700,500);
    c1-> Divide(1,2);
    c1->cd(1);

    TGraphErrors * GraphX = new TGraphErrors( iThetaMax, theta, rmsX, errTh, errSX);
    GraphX->SetTitle("Resolution en X (Theta)");
    GraphX->Draw("ALP");

    c1->cd(2);
    TGraphErrors * GraphY = new TGraphErrors( iThetaMax, theta, sigmaY, errTh, errSY);
    GraphY->SetTitle("Resolution en Y (Theta)");
    GraphY->Draw("ALP");

    cout<<"\n\n\n";
    MUONLoader->UnloadTracks();

}




