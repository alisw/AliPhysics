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

/* $Id:  */

/* $Log:
   29.05.2001 Yuri Kharlov:
              Everywhere reading the treese TTree->GetEvent(i)
              is replaced by reading the branches TBranch->GetEntry(0)
*/

//_________________________________________________________________________
//  A singleton. This class should be used in the analysis stage to get 
//  reconstructed objects: Digits, RecPoints, TrackSegments and RecParticles,
//  instead of directly reading them from galice.root file. This container 
//  ensures, that one reads Digits, made of these particular digits, RecPoints, 
//  made of these particular RecPoints, TrackSegments and RecParticles. 
//  This becomes non trivial if there are several identical branches, produced with
//  different set of parameters. 
//
//  An example of how to use (see also class AliPHOSAnalyser):
//  AliPHOSGetter * gime = AliPHOSGetter::GetInstance("galice.root","test") ;
//  for(Int_t irecp = 0; irecp < gime->NRecParticles() ; irecp++)
//     AliPHOSRecParticle * part = gime->RecParticle(1) ;
//     ................
//  please->GetEvent(event) ;    // reads new event from galice.root
//                  
//*-- Author: Yves Schutz (SUBATECH) & Dmitri Peressounko (RRC KI & SUBATECH)
//*--         Completely redesigned by Dmitri Peressounko March 2001  
//
//*-- YS June 2001 : renamed the original AliPHOSIndexToObject and make
//*--         systematic usage of TFolders without changing the interface        
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---

#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TObjString.h"
#include "TFolder.h"

// --- Standard library ---
#include <iostream.h>

// --- AliRoot header files ---

#include "AliRun.h"
#include "AliPHOSGetter.h"
#include "AliPHOS.h"
#include "AliPHOSDigitizer.h"
#include "AliPHOSSDigitizer.h"
#include "AliPHOSClusterizer.h"
#include "AliPHOSClusterizerv1.h"
#include "AliPHOSTrackSegmentMaker.h"
#include "AliPHOSTrackSegmentMakerv1.h"
#include "AliPHOSTrackSegment.h"
#include "AliPHOSPID.h" 
#include "AliPHOSPIDv1.h" 
#include "AliPHOSGeometry.h"

ClassImp(AliPHOSGetter)
  
  AliPHOSGetter * AliPHOSGetter::fgObjGetter = 0 ; 

//____________________________________________________________________________ 
AliPHOSGetter::AliPHOSGetter(const char* headerFile, const char* branchTitle )
{
  //Initialize  all lists

  fHeaderFile         = headerFile ; 
  fSDigitsTitle       = branchTitle ; 
  fDigitsTitle        = branchTitle ; 
  fRecPointsTitle     = branchTitle ; 
  fRecParticlesTitle  = branchTitle ; 
  fTrackSegmentsTitle = branchTitle ; 

  fPrimaries = new TObjArray(1) ;

  //open headers file
  TFile * file = (TFile*) gROOT->GetFile(fHeaderFile.Data() ) ;

 if(file == 0){    //if file was not opened yet, read gAlice
    if(fHeaderFile.Contains("rfio")) // if we read file using HPSS
      file =	TFile::Open(fHeaderFile.Data(),"update") ;
    else
      file = new TFile(fHeaderFile.Data(),"update") ;
  
  if (!file->IsOpen()) {
    cerr << "ERROR : AliPHOSGetter::AliPHOSGetter -> Cannot open " << fHeaderFile.Data() << endl ; 
    abort() ; 
  }
    
     gAlice = (AliRun *) file->Get("gAlice") ;
    
    if (!gAlice) {
      cerr << "ERROR : AliPHOSGetter::AliPHOSGetter -> Cannot find gAlice in " << fHeaderFile.Data() << endl ; 
      abort() ; 
    }
  }
}

//____________________________________________________________________________ 
void AliPHOSGetter::CreateWhiteBoard() const
{
  // Posts a few item to the white board (folders)
  
  // -- the geometry
  Post(fHeaderFile, "G") ; 
  
  // -- the hits
  Post(fHeaderFile, "H") ; 
    
}

//____________________________________________________________________________ 
AliPHOSGetter * AliPHOSGetter::GetInstance()
{
  // Returns the pointer of the unique instance already defined
  
  AliPHOSGetter * rv = 0 ;
  if ( fgObjGetter )
    rv = fgObjGetter ;
  else
    cout << "AliPHOSGetter::GetInstance ERROR: not yet initialized" << endl ;

  return rv ;
}

//____________________________________________________________________________ 
AliPHOSGetter * AliPHOSGetter::GetInstance(const char* headerFile,
					   const char* branchTitle)
{
  // Creates and returns the pointer of the unique instance
  // Must be called only when the environment has changed 


  if ( fgObjGetter )      // delete it if already exists
    delete fgObjGetter ; 

  fgObjGetter = new AliPHOSGetter(headerFile,branchTitle) ; 
  
  // Posts a few item to the white board (folders)
  fgObjGetter->CreateWhiteBoard() ;

  // Get the first event into the arrays posted on the white board
  // branchTitle = 0, means simulation run and no event yet 
  if (branchTitle) 
    fgObjGetter->Event(0) ;

  return fgObjGetter ; 
  
}

//____________________________________________________________________________ 
 const  AliPHOS * AliPHOSGetter::PHOS() const 
{
  // returns the PHOS object 
  return ( (AliPHOS*)gAlice->GetDetector("PHOS") ); 
}  

//____________________________________________________________________________ 
  const AliPHOSGeometry *  AliPHOSGetter::PHOSGeometry() const 
{
  // retrieves the geometr from the folder

  TString path("YSAlice/WhiteBoard/Geometry/PHOS/") ;
  path += PHOS()->GetTitle() ;
  return (AliPHOSGeometry*)gROOT->FindObjectAny(path.Data()) ; 
}

//____________________________________________________________________________ 
void AliPHOSGetter::Post(const char * headerFile, const char * opt, const char * name, const Int_t event) const 
{
  // Adds a new folder for summable digits 
 
  TString foldertitle ; 
  if ( event >= 0 ) 
    foldertitle += event ; 
  else 
    foldertitle = "" ; 

  TFolder * aliceF  = (TFolder*)gROOT->FindObjectAny("YSAlice") ; 
 
  if ( strcmp(opt, "G") == 0 ) { // Geometry
    // the hierarchy is //YSALICE/WhiteBoard/Geometry/PHOS/
    AliPHOSGeometry * geom = AliPHOSGeometry::GetInstance(PHOS()->GetTitle(),"") ;  
    TFolder * geomF  = (TFolder*)aliceF->FindObject("WhiteBoard/Geometry/PHOS") ; 
    if ( !geomF ) {
      cerr << "ERROR: AliPHOSGetter::Post G -> Folder WhiteBoard/Geometry/PHOS/" << " not found!" << endl;
      abort() ; 
    }    
    else 
      geomF->Add((TObject*)geom) ; 
   
  } else if ( strcmp(opt, "H") == 0 ) { // Hits
    // the hierarchy is //YSALICE/WhiteBoard/SDigits/PHOS/Hits
    TClonesArray * hits = new TClonesArray("AliPHOSHit",1000) ;
    hits->SetName("Hits") ; 
    TFolder * hitsF  = (TFolder*)aliceF->FindObject("WhiteBoard/Hits/PHOS") ; 
    if ( !hitsF ) {
      cerr << "ERROR: AliPHOSGetter::Post H -> Folder WhiteBoard/Hits/PHOS/" << " not found!" << endl;
      abort() ; 
    }    
    else 
      hitsF->Add(hits) ; 

  } else if ( strcmp(opt, "S") == 0 ) { // summable digits
    // the hierarchy is //YSALICE/WhiteBoard/SDigits/PHOS/headerFile/sdigitsname
    // because you can have sdigits from several hit files for mixing
    TClonesArray * sdigits = new TClonesArray("AliPHOSDigit",1000) ;
    TString sdigitsName ; 
    if (name) 
      sdigitsName = name ;
    else 
      sdigitsName = "SDigits" ;	
    sdigits->SetName( sdigitsName.Data() ) ; 
    TFolder * sdigitsF  = (TFolder*)aliceF->FindObject("WhiteBoard/SDigits/PHOS") ;
    TFolder * sdigitsF2 = 0 ; 
    TString subdir(headerFile) ;
    if ( !(sdigitsF2=(TFolder*)sdigitsF->FindObject(subdir)) ) 
      sdigitsF2 = sdigitsF->AddFolder(subdir, foldertitle); 
    else {
      if ( sdigitsF2->FindObject( sdigitsName.Data() ) ) {
	cerr <<"INFO: AliPHOSGetter::Post S -> Folder " << subdir << ", " << foldertitle
	     << " already exists!" << endl ;  
	return ; 
      }
    }
    if ( !sdigitsF2 ) {
      cerr << "ERROR: AliPHOSGetter::Post S -> Folder WhiteBoard/SDigits/PHOS/" << subdir << " not created!" << endl;
      abort() ; 
    }
    else 
      sdigitsF2->Add(sdigits) ;

  } else if ( strcmp(opt, "Ser") == 0 ) { // sdigizer
    // the hierarchy is //YSALICE/tasks/Digitizer/PHOS/sdigitsname
    AliPHOSSDigitizer * sdigitizer = new AliPHOSSDigitizer() ; 
    TString sdigitsName ; 
    if (name) 
      sdigitsName = name ;
    else 
      sdigitsName = "SDigitizer" ;
    sdigitizer->SetName( sdigitsName.Data() ) ; 
    TTask * sdigitsF  = (TTask*)aliceF->FindObject("tasks/SDigitizer") ; 
    if ( !sdigitsF ) {
      cerr << "ERROR: AliPHOSGetter::Post Ser -> Task tasks/SDigitizer" << " not found!" << endl;
      abort() ; 
    }        
    TTask * phos = (TTask*)sdigitsF->GetListOfTasks()->FindObject("PHOS") ; 
    if ( !phos )  {
      cerr <<"ERROR: AliPHOSGetter::Post Ser -> tasks/SDigitizer/PHOS" << " not found!" << endl; 
      abort() ; 
    } else {
      AliPHOSSDigitizer * phossd = (AliPHOSSDigitizer*)phos->GetListOfTasks()->FindObject(sdigitsName.Data()) ; 
      if (phossd) { 
	cout << "INFO: AliPHOSGetter::Post Ser -> Task " << sdigitsName.Data() << " already exists" << endl ; 
	return ; 
      } else 
	phos->Add(sdigitizer) ;
    } 

  } else if ( strcmp(opt, "D") == 0 ) { // digits
    // the hierarchy is //YSALICE/WhiteBoard/Digits/PHOS/digitsname
    TClonesArray * digits  = new TClonesArray("AliPHOSDigit",20000) ;
    TString digitsName ; 
    if (name) 
      digitsName = name ;
    else 
      digitsName = "Digits" ;	
    digits->SetName( digitsName.Data() ) ; 
    TFolder * digitsF  = (TFolder*)aliceF->FindObject("WhiteBoard/Digits/PHOS") ;
    if ( !digitsF ) {
      cerr << "ERROR: AliPHOSGetter::Post D -> Folder WhiteBoard/Digits/PHOS/" << " not found!" << endl;
      abort() ; 
    }    
    digitsF->SetTitle(foldertitle) ; 
    if ( digitsF->FindObject( digitsName.Data() ) ) {
      cerr <<"INFO: AliPHOSGetter::Post D -> Object " << digitsName.Data() 
	   << " already exists!" << endl ;  
      return ; 
    } 
    else 
      digitsF->Add(digits) ;  

  } else if ( strcmp(opt, "Der") == 0 ) { // sdigizer
    // the hierarchy is //YSALICE/tasks/Digitizer/PHOS/digitsname
    AliPHOSDigitizer * digitizer = new AliPHOSDigitizer() ; 
    TString digitsName ; 
    if (name) 
      digitsName = name ;
    else 
      digitsName = "Digitizer" ;
    digitizer->SetName( digitsName.Data() ) ; 
    TTask * digitsF  = (TTask*)aliceF->FindObject("tasks/Digitizer") ; 
    if ( !digitsF ) {
      cerr << "ERROR: AliPHOSGetter::Post Der -> Task tasks/Digitizer" << " not found!" << endl;
      abort() ; 
    }        
    TTask * phos = (TTask*)digitsF->GetListOfTasks()->FindObject("PHOS") ; 
    if ( !phos )  {
      cerr <<"ERROR: AliPHOSGetter::Post Der -> tasks/Digitizer/PHOS" << " not found!" << endl; 
      abort() ; 
    } else {
      AliPHOSDigitizer * phosd = (AliPHOSDigitizer*)phos->GetListOfTasks()->FindObject(digitsName.Data()) ; 
      if (phosd) { 
	cout << "INFO: AliPHOSGetter::Post Der -> Task " << digitsName.Data() << " already exists" << endl ; 
	return ; 
      } else 
	phos->Add(digitizer) ;
    }
 
  } else if ( strcmp(opt, "R") == 0 ) { // RecPoints
    // the hierarchy is //YSALICE/WhiteBoard/RecPoints/PHOS/emc/recpointsname
    //                  //YSALICE/WhiteBoard/RecPoints/PHOS/cpv/recpointsname
    TObjArray * emcrp  = new TObjArray(100) ;
    TObjArray * cpvrp  = new TObjArray(100) ;
    TString recpointsName ; 
    if (name) 
      recpointsName = name ;
    else       
      recpointsName = "RecPoints" ;	
    emcrp->SetName( recpointsName.Data() ) ; 
    cpvrp->SetName( recpointsName.Data() ) ; 
    TFolder * emcrpF  = (TFolder*)aliceF->FindObject("WhiteBoard/RecPoints/PHOS/emc") ; 
    TFolder * cpvrpF  = (TFolder*)aliceF->FindObject("WhiteBoard/RecPoints/PHOS/cpv") ;
    emcrpF->SetTitle(foldertitle) ; 
    cpvrpF->SetTitle(foldertitle) ; 
    if ( !emcrpF || !cpvrpF ) {
      cerr << "ERROR: AliPHOSGetter::Post R -> Folder WhiteBoard/RecPoints/PHOS/emc(cpv)" << " not found!" << endl;
      abort() ; 
    }    
    // TString title("PHOS Digits") ; 
    if ( emcrpF->FindObject( recpointsName.Data() ) ||  cpvrpF->FindObject( recpointsName.Data() ) ) {
      cerr <<"INFO: AliPHOSGetter::Post R -> Object " << recpointsName.Data() 
	   << " already exists!" << endl ;  
      return ; 
    } 
    else {
      emcrpF->Add(emcrp) ;  
      cpvrpF->Add(cpvrp) ;
    }  

  } else if ( strcmp(opt, "Rer") == 0 ) { // clusterizer
    // the hierarchy is //YSALICE/tasks/Reconstructionner/PHOS/recpointsname
    AliPHOSClusterizer * clusterizer; 
    if ( strstr(name, "clu-v1") != 0 )   
	 clusterizer = new AliPHOSClusterizerv1() ;
    else {
      cerr << "ERROR: AliPHOSGetter::Post Rer -> " << name << " unknown clusterizer version" << endl ;  
      abort() ; 
    }
    TString recpointsName ; 
    if (name) 
      recpointsName = name ;
    else 
      recpointsName = "Clusterizer" ;
    clusterizer->SetName( recpointsName.Data() ) ; 
    TTask * reF  = (TTask*)aliceF->FindObject("tasks/Reconstructioner") ; 
    if ( !reF ) {
      cerr << "ERROR: AliPHOSGetter::Post Rer -> Task tasks/Reconstructioner" << " not found!" << endl;
      abort() ; 
    }        
    TTask * phos = (TTask*)reF->GetListOfTasks()->FindObject("PHOS") ; 
    if ( !phos )  {
      cerr <<"ERROR: AliPHOSGetter::Post Rer -> tasks/Reconstructioner/PHOS" << " not found!" << endl; 
      abort() ; 
    } else {
       AliPHOSClusterizer * phoscl = (AliPHOSClusterizer*)phos->GetListOfTasks()->FindObject(recpointsName.Data()) ; 
      if (phoscl) { 
	cout << "INFO: AliPHOSGetter::Post Rer -> Task " << recpointsName.Data() << " already exists" << endl ; 
	return ; 
      } else 
	phos->Add(clusterizer) ;
    }

  } else if ( strcmp(opt, "T") == 0 ) { //TrackSegments
    // the hierarchy is //YSALICE/WhiteBoard/TrackSegments/PHOS/tracksegmentsname

    TClonesArray * tracksegments  = new TClonesArray("AliPHOSTrackSegment", 200) ;
    TString tracksegmentsName ; 
    if (name) 
      tracksegmentsName = name ;
    else       
      tracksegmentsName = "TrackSegments" ;	
    tracksegments->SetName( tracksegmentsName.Data() ) ; 
    TFolder * tracksegmentsF  = (TFolder*)aliceF->FindObject("WhiteBoard/TrackSegments/PHOS") ; 
    tracksegmentsF->SetTitle(foldertitle) ;  
    if ( !tracksegmentsF) {
      cerr << "ERROR: AliPHOSGetter::Post T -> Folder WhiteBoard/TrackSegments/PHOS" << " not found!" << endl;
      abort() ; 
    }    
    if ( tracksegmentsF->FindObject( tracksegmentsName.Data() ) ) {
      cerr <<"INFO: AliPHOSGetter::Post T -> Object " << tracksegmentsName.Data() 
	   << " already exists!" << endl ;  
      return ; 
    } 
    else 
      tracksegmentsF->Add(tracksegments) ;  

  } else if ( strcmp(opt, "Ter") == 0 ) { // TrackSegmentsMaker
    // the hierarchy is //YSALICE/tasks/Reconstructionner/PHOS/tracksegmentsname
    AliPHOSTrackSegmentMaker * tracksegmentmaker ; 
    if ( strstr(name, "tsm-v1") != 0 )   
	 tracksegmentmaker = new AliPHOSTrackSegmentMakerv1() ;
    else {
      cerr << "ERROR: AliPHOSGetter::Post Ter -> " << name << " unknown track segment maker version" << endl ;  
      abort() ; 
    }
    TString tracksegmentsName ; 
    if (name) 
      tracksegmentsName = name ;
    else 
      tracksegmentsName = "TrackSegmentMaker" ;
    tracksegmentmaker->SetName( tracksegmentsName.Data() ) ; 
    TTask * reF  = (TTask*)aliceF->FindObject("tasks/Reconstructioner") ; 
    if ( !reF ) {
      cerr << "ERROR: AliPHOSGetter::Post Ter -> Task tasks/Reconstructioner" << " not found!" << endl;
      abort() ; 
    }        
    TTask * phos = (TTask*)reF->GetListOfTasks()->FindObject("PHOS") ; 
    if ( !phos )  {
      cerr <<"ERROR: AliPHOSGetter::Post Ter -> tasks/Reconstructioner/PHOS" << " not found!" << endl; 
      abort() ; 
    } else {
       AliPHOSTrackSegmentMaker * phosts = (AliPHOSTrackSegmentMaker*)phos->GetListOfTasks()->FindObject(tracksegmentsName.Data()) ; 
      if (phosts) { 
	cout << "INFO: AliPHOSGetter::Post Ter -> Task " << tracksegmentsName.Data() << " already exists" << endl ; 
	return ; 
      } else 
	phos->Add(tracksegmentmaker) ;
    }
  
  } else if ( strcmp(opt, "P") == 0 ) { // RecParticles
    // the hierarchy is //YSALICE/WhiteBoard/RecParticles/PHOS/recparticlesname

    TClonesArray * recparticles  = new TClonesArray("AliPHOSRecParticle", 200) ;
    TString recparticlesName ; 
    if (name) 
      recparticlesName = name ;
    else       
      recparticlesName = "RecParticles" ;	
    recparticles->SetName( recparticlesName.Data() ) ; 
    TFolder * recparticlesF  = (TFolder*)aliceF->FindObject("WhiteBoard/RecParticles/PHOS") ; 
    recparticlesF->SetTitle(foldertitle) ;  
    if ( !recparticlesF) {
      cerr << "ERROR: AliPHOSGetter::Post P -> Folder WhiteBoard/RecParticles/PHOS" << " not found!" << endl;
      abort() ; 
    }    
    if ( recparticlesF->FindObject( recparticlesName.Data() ) ) {
      cerr <<"INFO: AliPHOSGetter::Post P -> Object " << recparticlesName.Data() 
	   << " already exists!" << endl ;  
      return ; 
    } 
    else 
      recparticlesF->Add(recparticles) ;  

  } else if ( strcmp(opt, "Per") == 0 ) { // PID Maker
    // the hierarchy is //YSALICE/tasks/Reconstructionner/PHOS/recparticlesname
    AliPHOSPID * pid ; 
    if ( strstr(name, "pid-v1") != 0 )   
	 pid = new AliPHOSPIDv1() ;
    else {
      cerr << "ERROR: AliPHOSGetter::Post Per -> " << name << " unknown PID maker version" << endl ;  
      abort() ; 
    }
    TString recparticlesName ; 
    if (name) 
      recparticlesName = name ;
    else 
      recparticlesName = "PID" ;
    pid->SetName( recparticlesName.Data() ) ; 
    TTask * reF  = (TTask*)aliceF->FindObject("tasks/Reconstructioner") ; 
    if ( !reF ) {
      cerr << "ERROR: AliPHOSGetter::Post Per -> Task tasks/Reconstructioner" << " not found!" << endl;
      abort() ; 
    }        
    TTask * phos = (TTask*)reF->GetListOfTasks()->FindObject("PHOS") ; 
    if ( !phos )  {
      cerr <<"ERROR: AliPHOSGetter::Post Per -> tasks/Reconstructioner/PHOS" << " not found!" << endl; 
      abort() ; 
    } else {
       AliPHOSPID * phospid = (AliPHOSPID*)phos->GetListOfTasks()->FindObject(recparticlesName.Data()) ; 
      if (phospid) { 
	cout << "INFO: AliPHOSGetter::Post Per -> Task " << recparticlesName.Data() << " already exists" << endl ; 
	return ; 
      } else 
	phos->Add(pid) ;
    }
  } 
  else if ( strcmp(opt, "QA") == 0 ) { // Alarms
    // the hierarchy is //YSALICE/WhiteBoard/Alarms/PHOS/

    TFolder * alarmsF  = new TFolder() ;
    TString alarmsName ; 
    if (name) 
      alarmsName = name ;
    else       
      alarmsName = "Alarm with no name" ;	
    alarmsF->SetName( alarmsName.Data() ) ; 
    alarmsF->SetTitle(foldertitle) ;  
    TFolder * qaaF  = (TFolder*)aliceF->FindObject("WhiteBoard/QAAlarms") ; 
    if ( !qaaF) {
      cerr << "ERROR: AliPHOSGetter::Post QA -> Folder WhiteBoard/QAAlarms/" << " not found!" << endl;
      return ; 
    }    
    if ( qaaF->FindObject( alarmsName.Data() ) ) 
      qaaF->RecursiveRemove(  qaaF->FindObject( alarmsName.Data() ) ) ; 
    
    qaaF->Add(alarmsF) ;  

  }
}

//____________________________________________________________________________ 
const TParticle * AliPHOSGetter::Primary(Int_t index) const
{
  // Return primary particle numbered by <index>

  if(index < 0) 
    return 0 ;
  
  Int_t primaryIndex = index % 10000000 ; 
  Int_t primaryList = (Int_t ) ((index-primaryIndex)/10000000.)  ;
  
  if ( primaryList > 0  ) {
    cout << " Getter does not support currently Mixing of primary " << endl ;
    cout << "   can not return primary: " << index<< " (list "<< primaryList<< " primary # " << primaryIndex << " )"<<endl ;
    return 0;
  }
  
  return gAlice->Particle(primaryIndex) ;
  
}

//____________________________________________________________________________ 
void AliPHOSGetter::ReadTreeD()
{
  // Read the digit tree gAlice->TreeD()  
  if(gAlice->TreeD()== 0){
    cerr <<   "ERROR: AliPHOSGetter::ReadTreeD: can not read TreeD " << endl ;
  return ;
  }
  
  TObjArray * lob = (TObjArray*)gAlice->TreeD()->GetListOfBranches() ;
  TIter next(lob) ; 
  TBranch * branch = 0 ; 
  TBranch * digitsbranch = 0 ; 
  TBranch * digitizerbranch = 0 ; 
  Bool_t phosfound = kFALSE, digitizerfound = kFALSE ; 
  
  while ( (branch = (TBranch*)next()) && (!phosfound || !digitizerfound) ) {
    if ( (strcmp(branch->GetName(), "PHOS")==0) && (strcmp(branch->GetTitle(), fDigitsTitle)==0) ) {
      digitsbranch = branch ; 
      phosfound = kTRUE ;
    }
    else if ( (strcmp(branch->GetName(), "AliPHOSDigitizer")==0) && (strcmp(branch->GetTitle(), fDigitsTitle)==0) ) {
      digitizerbranch = branch ; 
      digitizerfound = kTRUE ; 
    }
  }

  if ( !phosfound || !digitizerfound ) {
    cerr << "WARNING: AliPHOSGetter::ReadTreeD -> Cannot find Digits and/or Digitizer with name " << fDigitsTitle << endl ;
    return ; 
  }   
 

 // Post the Digits
  Post(fHeaderFile, "D", fDigitsTitle) ; 
  
   // Post the Digitizer
  Post(fHeaderFile, "Der", fDigitsTitle) ; 

  TClonesArray * digits = Digits(fDigitsTitle) ; 
  digits->Clear() ; 
  digitsbranch   ->SetAddress(&digits) ;

  AliPHOSDigitizer * digitizer = Digitizer(fDigitsTitle) ;
  digitizerbranch->SetAddress(&digitizer) ;

  digitsbranch   ->GetEntry(0) ;
  digitizerbranch->GetEntry(0) ;
  
}

//____________________________________________________________________________ 
void AliPHOSGetter::ReadTreeH()
{
  // Read the first entry of PHOS branch in hit tree gAlice->TreeH()

  if(gAlice->TreeH()== 0){
    cerr <<   "ERROR: AliPHOSGetter::ReadTreeH: -> Cannot read TreeH " << endl ;
    return ;
  }
  
  TBranch * hitsbranch = (TBranch*)gAlice->TreeH()->GetListOfBranches()->FindObject("PHOS") ;
  if ( !hitsbranch ) {
    cerr << "WARNING:  AliPHOSGetter::ReadTreeH -> Cannot find branch PHOS" << endl ; 
  } else {
    TClonesArray * hits = Hits() ; 
    hits->Clear() ; 
    hitsbranch->SetAddress(&hits) ;
    hitsbranch->GetEntry(0) ;
  }
}

//____________________________________________________________________________ 
void AliPHOSGetter::ReadTreeQA()
{
  // Read the digit tree gAlice->TreeQA()
  // so far only PHOS knows about this Tree  

  if(PHOS()->TreeQA()== 0){
    cerr <<   "ERROR: AliPHOSGetter::ReadTreeQA: can not read TreeQA " << endl ;
    return ;
  }
  
  TBranch * qabranch = PHOS()->TreeQA()->GetBranch("PHOS") ; 
  if (!qabranch) { 
    cerr << "WARNING: AliPHOSGetter::ReadTreeQA -> Cannot find QA Alarms for PHOS" << endl ;
    return ; 
  }   

 // Post the QA Alarms
  Post(fHeaderFile, "QA", "PHOS") ; 
  TFolder * alarmsF = Alarms() ; 
  alarmsF->Clear() ; 
  qabranch->SetAddress(&alarmsF) ;
  qabranch->GetEntry(0) ;
  
}

//____________________________________________________________________________ 
void AliPHOSGetter::ReadTreeR()
{
  // Read the reconstrunction tree gAlice->TreeR()

  if(gAlice->TreeR()== 0){
    cout <<   "ERROR: AliPHOSGetter::ReadTreeR: can not read TreeR " << endl ;
   return ;
  }
  // RecPoints 
  TObjArray * lob = (TObjArray*)gAlice->TreeR()->GetListOfBranches() ;
  TIter next(lob) ; 
  TBranch * branch = 0 ; 
  TBranch * emcbranch = 0 ; 
  TBranch * cpvbranch = 0 ; 
  TBranch * clusterizerbranch = 0 ; 
  Bool_t phosemcrpfound = kFALSE, phoscpvrpfound = kFALSE, clusterizerfound = kFALSE ; 
  
  while ( (branch = (TBranch*)next()) && (!phosemcrpfound || !phoscpvrpfound || !clusterizerfound) ) {
    if ( (strcmp(branch->GetName(), "PHOSEmcRP")==0) && (strcmp(branch->GetTitle(), fRecPointsTitle)==0) ) {
      emcbranch = branch ; 
      phosemcrpfound = kTRUE ;
    }
    else if ( (strcmp(branch->GetName(), "PHOSCpvRP")==0) && (strcmp(branch->GetTitle(), fRecPointsTitle)==0) ) {
      cpvbranch = branch ; 
      phoscpvrpfound = kTRUE ;
    }
    else if ( (strcmp(branch->GetName(), "AliPHOSClusterizer")==0) && (strcmp(branch->GetTitle(), fRecPointsTitle)==0) ) {
      clusterizerbranch = branch ; 
      clusterizerfound = kTRUE ; 
    }
  }

  if ( !phosemcrpfound || !phoscpvrpfound || !clusterizerfound ) {
    cerr << "WARNING: AliPHOSGetter::ReadTreeR -> Cannot find RecPoints and/or Clusterizer with name " << fRecPointsTitle << endl ;
    return ; 
  }   
 
 // Post the RecPoints
  Post(fHeaderFile, "R", fRecPointsTitle) ; 

  // Post the Clusterizer
  //  Need the version first
  AliPHOSClusterizer * clusterizer = 0 ; 
  clusterizerbranch->SetAddress(&clusterizer) ;
  clusterizerbranch->GetEntry(0) ;
  TString clusterizerName(fRecPointsTitle) ; 
  clusterizerName.Append(clusterizer->Version()) ; 
  delete clusterizer ;
  Post(fHeaderFile, "Rer", clusterizerName) ; 

  TObjArray * emcRecPoints = EmcRecPoints(fRecPointsTitle) ;
  emcRecPoints->Clear() ; 
  emcbranch->SetAddress(&emcRecPoints) ;

  TObjArray * cpvRecPoints = CpvRecPoints(fRecPointsTitle) ;
  cpvRecPoints->Clear() ; 
  cpvbranch->SetAddress(&cpvRecPoints) ;

  clusterizer = Clusterizer(clusterizerName) ;
  clusterizerbranch->SetAddress(&clusterizer) ;

  emcbranch        ->GetEntry(0) ;
  cpvbranch        ->GetEntry(0) ;
  clusterizerbranch->GetEntry(0) ;
 
  // TrackSegments
  next.Reset() ; 
  TBranch * tsbranch = 0 ; 
  TBranch * tsmakerbranch = 0 ; 
  Bool_t phostsfound = kFALSE, tsmakerfound = kFALSE ; 
  
  while ( (branch = (TBranch*)next()) && (!phostsfound || !tsmakerfound) ) {
    if ( (strcmp(branch->GetName(), "PHOSTS")==0) && (strcmp(branch->GetTitle(), fTrackSegmentsTitle)==0) ) {
      tsbranch = branch ; 
      phostsfound = kTRUE ;
    }
    else if ( (strcmp(branch->GetName(), "AliPHOSTrackSegmentMaker")==0) && (strcmp(branch->GetTitle(), fTrackSegmentsTitle)==0) ) {
      tsmakerbranch = branch ; 
      tsmakerfound  = kTRUE ; 
    }
  }

  if ( !phostsfound || !tsmakerfound ) {
    cerr << "WARNING: AliPHOSGetter::ReadTreeR -> Cannot find TrackSegments and/or TrackSegmentMaker with name " << fTrackSegmentsTitle << endl ;
    return ; 
  } 

 // Post the TrackSegments
  Post(fHeaderFile, "T", fTrackSegmentsTitle) ; 

  // Post the TrackSegment Maker
  //  Need the version first
  AliPHOSTrackSegmentMaker * tsmaker = 0 ; 
  tsmakerbranch->SetAddress(&tsmaker) ;
  tsmakerbranch->GetEntry(0) ;
  TString tsmakerName(fTrackSegmentsTitle) ; 
  tsmakerName.Append(tsmaker->Version()) ; 
  delete tsmaker ;
  Post(fHeaderFile, "Ter", tsmakerName) ; 

  TClonesArray * tracksegments = TrackSegments(fTrackSegmentsTitle) ;
  tracksegments->Clear() ; 
  tsbranch->SetAddress(&tracksegments) ;
 
  tsmaker = TrackSegmentMaker(tsmakerName) ; 
  tsmakerbranch->SetAddress(&tsmaker) ;

  tsmakerbranch    ->GetEntry(0) ;
  tsbranch         ->GetEntry(0) ;
   
  // RecParticles
  next.Reset() ; 
  TBranch * rpabranch = 0 ; 
  TBranch * pidbranch = 0 ; 
  Bool_t phosrpafound = kFALSE, pidfound = kFALSE ; 
  
  while ( (branch = (TBranch*)next()) && (!phosrpafound || !pidfound) ) {
    if ( (strcmp(branch->GetName(), "PHOSRP")==0) && (strcmp(branch->GetTitle(), fRecParticlesTitle)==0) ) {
      rpabranch = branch ; 
      phosrpafound = kTRUE ;
    }
    else if ( (strcmp(branch->GetName(), "AliPHOSPID")==0) && (strcmp(branch->GetTitle(), fRecParticlesTitle)==0) ) {
      pidbranch = branch ; 
      pidfound  = kTRUE ; 
    }
  }

  if ( !phosrpafound || !pidfound ) {
    cerr << "WARNING: AliPHOSGetter::ReadTreeR -> Cannot find RecParticles and/or PID with name " << fRecParticlesTitle << endl ;
    return ; 
  } 

  // Post the RecParticles
  Post(fHeaderFile, "P", fRecParticlesTitle) ; 

  // Post the PID
  //  Need the version first
  AliPHOSPID * pid = 0 ; 
  pidbranch->SetAddress(&pid) ;
  pidbranch->GetEntry(0) ;
  TString pidName(fRecParticlesTitle) ; 
  pidName.Append(pid->Version()) ; 
  delete pid ;

  Post(fHeaderFile, "Per", pidName) ; 

  TClonesArray * recParticles = RecParticles(fRecParticlesTitle) ; 
  recParticles->Clear() ; 
  rpabranch->SetAddress(&recParticles) ;

  pid = PID(pidName) ; 
  pidbranch->SetAddress(&pid) ;
  
  pidbranch        ->GetEntry(0) ;
  rpabranch        ->GetEntry(0) ;
  
}

//____________________________________________________________________________ 
void AliPHOSGetter::ReadTreeS()
{
  // Read the summable digits tree gAlice->TreeS()  

  if(gAlice->TreeS()== 0){
    cerr <<   "ERROR: AliPHOSGetter::ReadTreeS -> Cannot find TreeS " << endl ;
    return ;
  }
  
  TObjArray * lob = (TObjArray*)gAlice->TreeS()->GetListOfBranches() ;
  TIter next(lob) ; 
  TBranch * branch = 0 ; 
  TBranch * sdigitsbranch = 0 ; 
  TBranch * sdigitizerbranch = 0 ; 
  Bool_t phosfound = kFALSE, sdigitizerfound = kFALSE ; 
  
  while ( (branch = (TBranch*)next()) && (!phosfound || !sdigitizerfound) ) {
    if ( (strcmp(branch->GetName(), "PHOS")==0) && (strcmp(branch->GetTitle(), fSDigitsTitle)==0) ) {
      sdigitsbranch = branch ; 
      phosfound = kTRUE ;
    }
    else if ( (strcmp(branch->GetName(), "AliPHOSSDigitizer")==0) && (strcmp(branch->GetTitle(), fSDigitsTitle)==0) ) {
      sdigitizerbranch = branch ; 
      sdigitizerfound = kTRUE ; 
    }
  }

  if ( !phosfound || !sdigitizerfound ) {
    cerr << "WARNING: AliPHOSGetter::ReadTreeS -> Cannot find SDigits and/or SDigitizer with name " << fSDigitsTitle << endl ;
    return ; 
  }   

  // -- the SDigits 
  Post(fHeaderFile, "S", fSDigitsTitle) ; 

  // Post the SDigitizer
  Post(fHeaderFile, "Ser", fSDigitsTitle) ; 
  
  TClonesArray * sdigits = SDigits(fSDigitsTitle) ; 
  sdigits->Clear() ; 
  sdigitsbranch->SetAddress(&sdigits) ;

  AliPHOSSDigitizer * sdigitizer = SDigitizer(fSDigitsTitle) ;
  sdigitizerbranch->SetAddress(&sdigitizer) ;

  sdigitsbranch->GetEvent(0) ;
  sdigitizerbranch->GetEvent(0) ;
    
}

//____________________________________________________________________________ 
void AliPHOSGetter::ReadPrimaries()
{
  // Reads specific branches of primaries
  
  fNPrimaries = gAlice->GetNtrack();
  
  //   //Check, is it necessary to open new files
  //   TArrayI* events = fDigitizer->GetCurrentEvents() ; 
  //   TClonesArray * filenames = fDigitizer->GetHeadersFiles() ;
//   Int_t input ;
//   for(input = 0; input < filenames->GetEntriesFast(); input++){

//     TObjString * filename = (TObjString *) filenames->At(input) ;

//     //Test, if this file already open
//     TFile *file = (TFile*) gROOT->GetFile( filename->GetString() ) ;
//     if(file == 0)
//       file = new TFile( filename->GetString()) ;
//     file->cd() ;
    
//     // Get Kine Tree from file
// //     char treeName[20];
// //     sprintf(treeName,"TreeK%d",events->At(input));
// //     TTree * treeK = (TTree*)gDirectory->Get(treeName);
// //     if (treeK) 
// //       treeK->SetBranchAddress("Particles", &fParticleBuffer);
// //     else    
// //       cout << "AliPHOSGetter: cannot find Kine Tree for event:" << events->At(input) << endl;

// //     // Create the particle stack
// //     if(!fParticles) fParticles = new TClonesArray("TParticle",1000);
// //     // Build the pointer list
// //     if(fParticleMap) {     <----
// //       fParticleMap->Clear();
// //       fParticleMap->Expand(treeK->GetEntries());
// //     } else
// //       fParticleMap = new TObjArray(treeK->GetEntries());
    
//     // From gAlice->Particle(i) 


// //   if(!(*fParticleMap)[i]) {
// //     Int_t nentries = fParticles->GetEntries();
    
// //     // algorithmic way of getting entry index
// //     // (primary particles are filled after secondaries)
// //     Int_t entry;
// //     if (i<fHeader.GetNprimary())
// //       entry = i+fHeader.GetNsecondary();
// //     else 
// //       entry = i-fHeader.GetNprimary();
      
// //     // only check the algorithmic way and give
// //     // the fatal error if it is wrong
// //     if (entry != fParticleFileMap[i]) {
// //       Fatal("Particle",
// //         "!!!! The algorithmic way is WRONG: !!!\n entry: %d map: %d",
// // 	entry, fParticleFileMap[i]); 
// //     }  
      
// //     fTreeK->GetEntry(fParticleFileMap[i]);
// //     new ((*fParticles)[nentries]) TParticle(*fParticleBuffer);
// //     fParticleMap->AddAt((*fParticles)[nentries],i);
// //   }
// //   return (TParticle *) (*fParticleMap)[i];

   
    
//   }


//   //scan over opened files and read corresponding TreeK##

  return ;
}
//____________________________________________________________________________ 
void AliPHOSGetter::Event(Int_t event)
{
  // Reads the content of all Tree's S, D and R
  
  if ( event > gAlice->TreeE()->GetEntries() ) {
    cerr << "ERROR: AliPHOSGetter::Event -> There are only " << gAlice->TreeE()->GetEntries() << " events in this file" << endl ; 
    return ;
  }
  
  gAlice->GetEvent(event) ;
  
  ReadTreeH() ;
  ReadTreeS() ;
  ReadTreeD() ;
  ReadTreeR() ;
  ReadTreeQA() ;
  ReadPrimaries() ;

}

//____________________________________________________________________________ 
const TObject * AliPHOSGetter::ReturnO(TString what, TString name, TString file) const 
{
  // get the object named "what" from the folder
  // folders are named like //YSAlice/WhiteBoard/what/PHOS/name

  if ( file.IsNull() ) 
    file = fHeaderFile ; 
  TString path("WhiteBoard/") ;
  if ( name.IsNull() ) {
    if ( what.CompareTo("Hits") == 0 ) {
      path += what ; 
      path += "/PHOS/"; 
      path += what ; 
    }
    else if ( what.CompareTo("SDigits") == 0 ) { 
      path += what ; 
      path += "/PHOS/"; 
      path += file ; 
      path += "/" ; 
      path += fSDigitsTitle ; 
    }
    else if ( what.CompareTo("Digits") == 0 ){
      path += what ; 
      path += "/PHOS/"; 
      path += fDigitsTitle ;
    } 
    else if ( what.CompareTo("EmcRecPoints") == 0 ) {
      path += "RecPoints/PHOS/";  
      path += "emc/" ; 
      path += fRecPointsTitle ; 
    }
    else if ( what.CompareTo("CpvRecPoints") == 0 ) {
      path += "RecPoints/PHOS/";  
      path += "cpv/" ; 
      path += fRecPointsTitle ; 
    }
    else if ( what.CompareTo("TrackSegments") == 0 ) {
      path += "TrackSegments/PHOS/";  
      path += fTrackSegmentsTitle ; 
    }  
    else if ( what.CompareTo("RecParticles") == 0 ) {
      path += "RecParticles/PHOS/";  
      path += fRecParticlesTitle ; 
    }  
     else if ( what.CompareTo("Alarms") == 0 ) {
      path += "QAAlarms/PHOS";   
    }  
  } 
  else {
    if ( what.CompareTo("SDigits") == 0 ) { 
      path += what ; 
      path += "/PHOS/"; 
      path += file ; 
      path += "/" ; 
      path += name ;
    } 
    else if ( what.CompareTo("Digits") == 0 ) {
      path += what ; 
      path += "/PHOS/"; 
      path += name ; 
    }
    else if ( what.CompareTo("EmcRecPoints") == 0 ) {
      path += "RecPoints/PHOS/";  
      path += "emc/" ; 
      path += name ; 
    }
    else if ( what.CompareTo("CpvRecPoints") == 0 ) {
      path += "RecPoints/PHOS/";  
      path += "cpv/" ; 
      path += name ; 
    }  
    else if ( what.CompareTo("TrackSegments") == 0 ) {
      path += "TrackSegments/PHOS/";  
      path += name ; 
    } 
    else if ( what.CompareTo("RecParticles") == 0 ) {
      path += "RecParticles/PHOS/";  
      path += name ; 
    } 
    else if ( what.CompareTo("Alarms") == 0 ) {
      path += "QAAlarms/PHOS/";  
      path += name ; 
    } 
  }
  path.Prepend("YSAlice/") ;
  TObject * phosO  = (TObject*)gROOT->FindObjectAny(path) ; 
  if (!phosO) {
    cerr << "ERROR : AliPHOSGetter::ReturnO -> Object " << path << " not found!" << endl ; 
    abort() ;
  }
  return phosO ;
}
  
//____________________________________________________________________________ 
const TTask * AliPHOSGetter::ReturnT(TString what, TString name) const 
{
  // get the TTask named "what" from the folder
  // folders are named like //YSAlice/Tasks/what/PHOS/name

  TString path("tasks") ;
 
  if ( what.CompareTo("SDigitizer") == 0 ) 
    path += "/SDigitizer" ;
  else if ( what.CompareTo("Digitizer") == 0 ) 
    path += "/Digitizer" ; 
  else if ( what.CompareTo("Clusterizer") == 0 ) 
    path += "/Reconstructioner" ; 
  else if ( what.CompareTo("TrackSegmentMaker") == 0 ) 
    path += "/Reconstructioner" ; 
  else if ( what.CompareTo("PID") == 0 ) 
    path += "/Reconstructioner" ; 

  TFolder * aliceF  = (TFolder*)gROOT->FindObjectAny("YSAlice") ; 
  TTask * aliceT  = (TTask*)aliceF->FindObject(path) ; 

  if (!aliceT) {
    cerr << "ERROR: AliPHOSGetter::ReturnT -> Task " << path << " not found!" << endl ;  
    abort() ;
  }

  TTask * phosT = (TTask*)aliceT->GetListOfTasks()->FindObject("PHOS") ; 
  if (!phosT) { 
    cerr << "ERROR: AliPHOSGetter::ReturnT -> Task " << path << "/PHOS not found!" << endl ;  
    abort() ;
  }
  TList * l = phosT->GetListOfTasks() ; 
 
  if (what.CompareTo("SDigitizer") == 0) {  
    if ( name.IsNull() )
      name =  fSDigitsTitle ; 
  } else  if (what.CompareTo("Digitizer") == 0){ 
    if ( name.IsNull() )
      name =  fDigitsTitle ;
  } else  if (what.CompareTo("Clusterizer") == 0){ 
    if ( name.IsNull() )
      name =  fRecPointsTitle ;
  }
  else  if (what.CompareTo("TrackSegmentMaker") == 0){ 
    if ( name.IsNull() )
      name =  fTrackSegmentsTitle ;
  }
  else  if (what.CompareTo("PID") == 0){ 
    if ( name.IsNull() )
      name =  fRecParticlesTitle ;
  }
  
  TTask * task = (TTask*)l->FindObject(name) ; 

  if (!task)
    cout << "WARNING: AliPHOSGetter::ReturnT -> Task " << path << "/" << name << " not found!" << endl ; 
  
  return task ;
}
