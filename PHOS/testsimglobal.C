// Root
#include "TH2.h"
#include "TH1.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TString.h"
#include "TF1.h"
#include "TFormula.h"
#include "TFolder.h"
#include "TStopwatch.h"
#include "TObjArray.h"

//AliRoot
#include "AliPHOSHit.h"
#include "AliPHOSLoader.h"
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

Bool_t sim_exam() ;
Bool_t sdigit() ;
Bool_t digit() ;
Bool_t recpoint() ;
Bool_t track() ;
Bool_t particle() ;
void write_info(TString) ;


//____________________________________________________________________________ 
void testsimglobal() 
{
  
  gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include -I.");
  Bool_t error = kFALSE ;
  TString mess("") ;
  
  if ( sim_exam() ){
    mess = "Examination ended successfully." ;
    write_info(mess) ;
  }
  else
    error = kTRUE ;
  
  if (!error){
    mess = "Beginning of the PHOS reconstruction."  ;
    write_info(mess) ;
    if(sdigit()){
      mess = "SDigits process ended successfully."  ;
      write_info(mess) ;
    }
    else
      error = kTRUE ;
  }
  
  if (!error){
    if (digit()){
      mess = "Digits process ended successfully."  ;
      write_info(mess) ;
    }
    else
      error = kTRUE ;
  }
  
  if (!error){
    if (recpoint()){
      mess = "Cluster process ended successfully."  ;
      write_info(mess) ;
    }
    else
      error = kTRUE ;
  }
  
  if (!error){
    if (track()){
      mess = "TrackSegments process ended successfully."  ;
      write_info(mess) ;
    }
    else
      error = kTRUE ;
  }
  
  if (!error){
    if (particle()){
      mess = "RecParticles process ended successfully."  ;
      write_info(mess) ;
    }
    else
      error = kTRUE ;
  }
  
  if(!error){
    mess = "reconstruction ended successfully."  ;
    write_info(mess) ;
  }
  else {
    gSystem->Exec("uuencode $ALICE_ROOT/PHOS/testPHOS.root testPHOS.root | mail -s 'PHOS INSTALLATION ERROR' schutz@in2p3.fr") ;
    mess = "Error message to PHOS director has been sent, please wait for his answer before launching the whole simulation again."  ;
    write_info(mess) ;
  }
}


//____________________________________________________________________________ 
Bool_t sim_exam()
{
  
  // Definition of the alarm bounds for the examination
  const Float_t maxAlaHitsM = 12.79  ;  // total multiplicity
  const Float_t maxAlaTotEn = 19.34  ;    // total energy
  const Float_t maxAlaTotEnB = 18.35 ;     // per block multiplicity
  const Float_t maxAlaHitsMB = 11.1  ;  // per block energy
  
  TString mess("") ;

  mess = "Beginning of the simulation examination." ;
  write_info(mess) ;
  
  // define the array in which the events that have not reach th EMCA will be put.
  Bool_t error = kFALSE ;
  
  TStopwatch timer ;
  timer.Start() ;
   
  // Get the number of events generated in the simulation
  
  AliRunLoader* rl = AliRunLoader::Open("testPHOS.root");
  
  AliPHOSLoader* gim = (AliPHOSLoader*)rl->GetLoader("PHOSLoader");
  
  Int_t maxevent = rl->GetNumberOfEvents() ;
  
  
  // Examine the alarms
  TObjArray * alahm = dynamic_cast<TObjArray*>(gime->Alarms()->FindObject("HitsM")) ;
  Float_t ratiohm = 100.0*static_cast<Float_t>(alahm->GetEntries())/static_cast<Float_t>(maxevent) ;
  
  TObjArray * alaet = dynamic_cast<TObjArray*>(gime->Alarms()->FindObject("TotEn")) ;
  Float_t ratioet = 100.0*static_cast<Float_t>(alaet->GetEntries())/static_cast<Float_t>(maxevent) ;
  
  // Define the alarms per block and examine them
  char namemul[80], namen[80] ;
  TObjArray* alahmb[5] ;
  TObjArray* alaenb[5] ;
  Float_t ratiohmb[5], ratioenb[5] ;
  
  Int_t i = 0 ;
  for (i = 0 ; i < 5 ; i++){
    sprintf(namemul,"HitsMB%d",i+1) ;
    sprintf(namen,"TotEnB%d",i+1) ;
    alahmb[i]   = (TObjArray*) (gime->Alarms()->FindObject(namemul)) ;
    ratiohmb[i] = 100.0*(Float_t)alahmb[i]->GetEntries()/(Float_t)maxevent ;
    alaenb[i]   = (TObjArray*)(gime->Alarms()->FindObject(namen)) ;
    ratioenb[i] = 100.0*(Float_t)alaenb[i]->GetEntries()/(Float_t)maxevent ;
    
    if (ratiohmb[i] > maxAlaHitsMB){
      error = kTRUE ;
      mess = "Examination detected an error in " + TString(namemul) ;
      write_info(mess) ;
    }
    
    if (ratioenb[i] > maxAlaTotEnB) {
      error = kTRUE ;
      mess = "Examination detected an error in " + TString(namen)  ;
      write_info(mess) ;
    }
  }
  
  
  timer.Stop() ;
  timer.Print() ;
  
  
  if (ratiohm > maxAlaHitsM){
    error = kTRUE ;
    mess = "Examination detected an error in HitsM." ;
    write_info(mess) ;
  }
  
  if (ratioet>maxAlaTotEn){
    error = kTRUE ;
    mess = "Examination detected an error in TotEn." ;
    write_info(mess) ;
  }
  
  // Condition that will launch the general loop that builds the histograms in order to allow a further analysis.
  
  if (!error)
    return kTRUE ;
  else {
    mess = "Examination sets up the file that will be sent to PHOS director (30s)." ;
    write_info(mess) ;
    
    Int_t index = 0 ;
    Int_t nhits = 0 ;
    
    TH1F * his    = new TH1F("Total Multiplicity", "Total Multiplicity in PHOS", 200, 0., 200.) ;
    TH1F * hisnrg = new TH1F("Total Energy", "Total energy in PHOS", 200, 0., 12.) ;
    
    // name and define the different histograms per block involved in the analysis
    
    TClonesArray hisba("TH1F") ;
    TClonesArray hisbanrg("TH1F") ;
    Int_t i = 0 ;
    
    char name[80], title[80] ;
    for (i = 0 ; i < 6 ; i++) {
      sprintf(name,"multiplicity for block %d",i) ;
      sprintf(title,"multiplicity per blocks, block %d",i) ;
      new(hisba[i]) TH1F(name, title, 100, 0., 100.) ;
      
      sprintf(name,"total energy for block %d",i) ;
      sprintf(title,"total energy per block, block %d",i) ;
      new(hisbanrg[i]) TH1F(name, title, 200, 0., 12.) ;
    }
    
    // define the global histograms, the hit pointer and give the means to get the actual block reached by the photon
    
    AliPHOSHit * hit ;
    TH1F * hist ;
    TH1F * histnrg ;
    TH2F * hbiz = new TH2F ("hbiz","hbiz", 200, 0., 200., 200, 0., 12.) ;
    const AliPHOSGeometry * geom = gime->PHOSGeometry() ;
    
    // the very general loop
    
    for (index = 0 ; index < maxevent ; index++) {
      //get the number of the event
      gime->Event(index) ;
      // get the number of cells reached during this event and fill the total multiplicity histogram
      Int_t n = gime->Hits()->GetEntries() ;
      nhits += n ;
      his->Fill(n) ;
      // Get the data per block
      const TClonesArray * hita = static_cast<const TClonesArray *>(gime -> Hits()) ;
      TIter next(hita) ;
      Float_t Et = 0. ;
      Int_t id = 0, block = 0 ;
      Int_t nhit[6], rid[4] ;
      Float_t etblock[6] ;
      Int_t i = 0 ;
      
      for ( i = 0; i < 6 ; i++) {
	nhit[i] = 0 ;
	etblock[i] = 0 ;
      }
      
      while ( (hit = static_cast<AliPHOSHit *>(next())) ) {
	id = hit->GetId() ;
	if  (geom->IsInEMC(id) ) {
	  Et += ( hit -> GetEnergy()) ;
	  geom->AbsToRelNumbering(id,rid) ;
	  block = rid[0] ;
	  nhit[block]++ ;
	  etblock[block] +=  ( hit -> GetEnergy()) ;
	}
      }
      
      //Fill all the histograms but total multiplicity, already done
      hist = static_cast<TH1F*>(hisba.At(block)) ;
      hist->Fill(nhit[block]) ;
      histnrg = static_cast<TH1F*>(hisbanrg.At(block)) ;
      histnrg->Fill(etblock[block]) ;
      hisnrg -> Fill(Et) ;
      hbiz->Fill(n,Et) ;
    }
         
    TFile * file = gROOT -> GetFile("testPHOS.root") ;
    file -> Write() ;
    file->Close() ;
    return kFALSE ;
  }
}


//____________________________________________________________________________ 
Bool_t sdigit()
{
  //SDigits process
  
  const Float_t maxSDigits = 62.89 ;
  const Float_t widSDigits = TMath::Sqrt(maxSDigits) ;

  TString mess("") ;
  TString reconame = "test suite" ;
  
  AliPHOSSDigitizer *sd = new AliPHOSSDigitizer("testPHOS.root",reconame.Data()) ;
  AliPHOSGetter * gime = AliPHOSGetter::GetInstance() ;
  
  sd->ExecuteTask("deb") ;
  
  Float_t nSDigits =  static_cast<Float_t>(gime->SDigitizer()->GetSDigitsInRun()) / static_cast<Float_t>(gime->MaxEvent()) ;
  if ( nSDigits < maxSDigits-widSDigits ||
       nSDigits > maxSDigits+widSDigits ) {
    mess = "Error detected in the SDigits process. Sending error file to PHOS director." ;
    cout <<  "sdigit() : nsDigits = " << nSDigits 
	 << " maxSDigits,widSDigits= " << maxSDigits << "," << widSDigits << endl ;    
    write_info(mess) ;
    return kFALSE ;
  }
  else
    return kTRUE ;
}


//____________________________________________________________________________ 
Bool_t digit()
{
  
  //Digits process
  AliPHOSGetter * gime = AliPHOSGetter::GetInstance("testPHOS.root") ;
  TString reconame = "test suite" ;
  const Float_t maxDigits = 2860. ;
  const Float_t widDigits = TMath::Sqrt(maxDigits) ;
  
  TString mess("") ;

  AliPHOSDigitizer *d = new AliPHOSDigitizer("testPHOS.root",reconame.Data()) ;
  
  d->ExecuteTask("deb") ;
  
  Float_t nDigits = static_cast<Float_t>(gime->Digitizer()->GetDigitsInRun()) / static_cast<Float_t>(gime->MaxEvent()) ;
  
  if ( nDigits < maxDigits-widDigits || nDigits > maxDigits+widDigits ) {
    cout <<  "digit() : nDigits = " << nDigits 
	 << " maxDigits,widDigits= " << maxDigits << "," << widDigits << endl ;    
    mess = "Error detected in the Digits process. Sending error file to PHOS director." ;
    write_info(mess) ;
    return kFALSE ;
  }
  else
    return kTRUE ;
}


//____________________________________________________________________________ 
Bool_t recpoint()
{
  
  //RecPoints process
  AliPHOSGetter * gime = AliPHOSGetter::GetInstance("testPHOS.root") ;
  TString reconame = "test suite" ;
  
  const Float_t maxRecPoints = 1.0 ;
  const Float_t widRecPoints = TMath::Sqrt(maxRecPoints) ;
  
  TString mess("") ;

  AliPHOSClusterizer * cluster = new  AliPHOSClusterizerv1("testPHOS.root", reconame.Data()) ;
  
  cluster->ExecuteTask("deb") ;
  
  Float_t nRecPoints =  static_cast<Float_t>(gime->Clusterizer(reconame.Data())->GetRecPointsInRun()) /
    static_cast<Float_t>(gime->MaxEvent()) ;
  
  if ( nRecPoints < maxRecPoints-widRecPoints
       || nRecPoints > maxRecPoints+widRecPoints ) {
    cout <<  "recpoint() : nRecPoints = " << nRecPoints 
	 << " maxRecPoints,widRecPoints= " << maxRecPoints << "," << widRecPoints << endl ;    
    mess = "Error detected in the Clusterizing process. Sending error file to PHOS director." ;
    write_info(mess) ;
    return kFALSE ;
  }
  else
    return kTRUE ;
}


//____________________________________________________________________________ 
Bool_t track()
{
  
  //TrackSegments process
  AliPHOSGetter * gime = AliPHOSGetter::GetInstance("testPHOS.root") ;
  const Float_t maxTrackSegments = 1 ;

  TString mess("") ;
  TString reconame = "test suite" ;
  
  AliPHOSTrackSegmentMaker * tracks = new AliPHOSTrackSegmentMakerv1("testPHOS.root",reconame.Data()) ;
  
  tracks->ExecuteTask("deb") ;
  
  Float_t nTrackSegments =  static_cast<Float_t> (gime->TrackSegmentMaker(reconame.Data())->GetTrackSegmentsInRun()) /
    static_cast<Float_t>(gime->MaxEvent()) ;
  
  if ( nTrackSegments < maxTrackSegments-0.25 ||
       nTrackSegments > maxTrackSegments+0.25 ) {
    cout <<  "track() : nTrackSegments = " << nTrackSegments
	 << " maxTrackSegments,widTrackSegments= " << maxTrackSegments << "," << "0.25" << endl ;    
    mess = "Error detected in the TrackSegments process. Sending error file to PHOS director." ;
    write_info(mess) ;
    return kFALSE ;
  }
  else
    return kTRUE ;
}


//____________________________________________________________________________ 
Bool_t particle()
{
  
  //RecParticles process
  AliPHOSGetter * gime = AliPHOSGetter::GetInstance("testPHOS.root") ;
  const Float_t maxRecParticles = 1 ;

  TString mess("") ;
  TString reconame = "test suite" ;
  
  AliPHOSPID * pid = new AliPHOSPIDv1("testPHOS.root",reconame.Data()) ;
  
  pid->ExecuteTask("deb") ;
  
  Float_t nRecParticles =  static_cast<Float_t> (gime->PID(reconame.Data())->GetRecParticlesInRun()) /
    static_cast<Float_t>(gime->MaxEvent()) ;
  
  
  if ( nRecParticles < maxRecParticles-0.25 ||
       nRecParticles > maxRecParticles+0.25 ) {
    cout <<  "particle() : nRecParticles = " << nRecParticles 
	 << " maxRecParticles,widRecParticles= " << maxRecParticles << "," << "0.25" << endl ;    
    mess = "Error detected in the RecParticles process. Sending error file to PHOS director.Stop reconstruction." ;
    write_info(mess) ;
    return kFALSE ;
  }
  else
    return kTRUE ;
}


//____________________________________________________________________________ 
void write_info(TString mess)
{
  cerr << " _____________________________________________________________ " << endl
       << " " << endl
       << "MESS ==> " << mess <<endl
       << " _____________________________________________________________ " <<endl ;
}

