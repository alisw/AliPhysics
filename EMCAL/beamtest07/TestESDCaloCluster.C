#if !defined(__CINT__) || defined(__MAKECINT__)

//Root include files
#include <Riostream.h>
#include <TFile.h>
#include <TChain.h>
#include <TParticle.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TObjArray.h>
#include <TSystem.h>
#include <TString.h> 
#include <TH1F.h> 
#include <TVector.h> 
#include <TParticle.h> 

//AliRoot include files
#include "AliESD.h"
#include "AliESDEvent.h"
#include "AliEMCALLoader.h"
#include "AliESDCaloCluster.h"
#include "AliEMCALRecPoint.h"
#include "AliPID.h"
#include "AliLog.h" 

//macros

#endif

TChain * AliReadESDfromdisk(const UInt_t eventsToRead, 
				   const TString dirName, 
				   const TString esdTreeName, 
				   const char *  pattern) 
{
  // Reads ESDs from Disk
 TChain *  rv = 0; 
  
  // create a TChain of all the files 
  TChain * cESDTree = new TChain(esdTreeName) ; 

  // read from the directory file until the require number of events are collected
  void * from = gSystem->OpenDirectory(dirName) ;
   if (!from) 
     rv = 0 ;   
  else{ // reading file names from directory
    const char * subdir ; 
    // search all subdirectories witch matching pattern
    while( (subdir = gSystem->GetDirEntry(from))  && 
	   (cESDTree->GetEntries() < eventsToRead)) {
      if ( strstr(subdir, pattern) != 0 ) { 
	char file[200] ; 
        sprintf(file, "%s%s/AliESDs.root", dirName.Data(), subdir); 	
	cESDTree->Add(file) ;
      }
    } // while file
  
    rv = cESDTree ; 
    
  } // reading file names from directory
  return rv ; 
}

//======================================================================
TChain * AliReadESD(const UInt_t eventsToRead,
		  const TString dirName, 
		  const TString esdTreeName, 
		  const char *  pattern)  
{
  // Read AliESDs files and return a Chain of events
 
  if ( dirName == "" ) 
    return 0 ; 
  if ( esdTreeName == "" ) 
    return AliReadESDfromdisk(eventsToRead, dirName,"","") ;//Last 2 arguments are not necessary but pdsf compiler complains "","'
  else if ( strcmp(pattern, "") == 0 )
    return AliReadESDfromdisk(eventsToRead, dirName, esdTreeName,"") ;//Last argument is not necessary but pdsf compiler complains "","'
  else 
    return AliReadESDfromdisk(eventsToRead, dirName, esdTreeName, pattern) ;	    
}

//=====================================================================
//  Do:
//  .L TestESDCaloCluster.C++
//  TestESDCaloCluster(number of events to process)
//=====================================================================
void TestESDCaloCluster(const UInt_t eventsToProcess = 5, 
			TString dirName = "./", 
			const TString esdTreeName = "esdTree", 
			const char *  pattern = ".") 
{ 
  
  //Create chain of esd trees
  //AliReadESD(eventsToProcess, directoryName,esdTreeName,patternOfDirectory) ;
  //By default the root files are in the same directory 
  TChain * t = AliReadESD(eventsToProcess, dirName,esdTreeName,pattern) ; 
 
  
  // ESD
  AliESDEvent * esd = new AliESDEvent();
  esd->ReadFromTree(t);
  
  //Define few variables to be used in macro
  TString alirunName = "" ; 

  //Define example histograms
  TH1F * hEnergy = new TH1F("hEnergy","Energy Distribution",100,0.,100.);   
  TH1F * hEta    = new TH1F("hEta","Eta Distribution",100,-0.7,0.7);
  TH1F * hPhi    = new TH1F("hPhi","Phi Distribution",100,0,2*TMath::Pi());

  Int_t beg, end ;
  UInt_t event ;
  Float_t pos[3] ; 

  for (event = 0; event < eventsToProcess; event++) {//event loop
    //AliInfo( Form("Event %d \n",event) );  
    Int_t nbytes = t->GetEntry(event); // store event in esd
    //cout<<"nbytes "<<nbytes<<endl;
    if ( nbytes == 0 ) //If nothing in ESD 
      break ; 
    
    // Check that name of file is correct
    if (alirunName != t->GetFile()->GetName()) {        
      alirunName = t->GetFile()->GetName() ; 
      alirunName.ReplaceAll("galice.root", "AliESDs.root") ;
    }
    
    //get reconstructed vertex position 
    
    //Double_t vertex_position[3] ; 
    //    esd->GetVertex()->GetXYZ(vertex_position) ; 
    
    cout<<"Event >>>>>>>>>>> "<<event<<endl;
    

    //select EMCAL tracks only 
    end = esd->GetNumberOfCaloClusters() ;  
    beg = 0;
    Int_t nphP ;
    //cout<<"begin "<<beg<<" end "<<end<<endl;
    Int_t iclus = 0;
    for (nphP =  beg; nphP <  end; nphP++) {//////////////EMCAL cluster loop
      AliESDCaloCluster * clus = esd->GetCaloCluster(nphP) ; // retrieve cluster from esd
      
      //Get the cluster parameters
      Float_t energy   = clus->E() ;
      cout << "Cluster " << nphP << " energy = " << energy << endl;

	//	Int_t mult       = clus->GetNumberOfDigits() ;
	//	Float_t disp     = clus->GetClusterDisp() ;
	// 	UShort_t * amp   = clus->GetDigitAmplitude() ;
	// 	UShort_t * time  = clus->GetDigitTime() ;
	// 	UShort_t * index = clus->GetDigitIndex() ;
	//	Int_t iprim = clus->GetPrimaryIndex();

      /*
      clus->GetPosition(pos) ;
      TVector3 vpos(pos[0],pos[1],pos[2]) ;
      //Print values on screen
      cout<<"ESD cluster "<<iclus <<"; Energy "<<energy
	  <<"; Phi "<<vpos.Phi()*180/TMath::Pi()<<"; Eta "<<vpos.Eta()<<endl;
      //cout<<"Dispersion "<<disp<<"; multiplicity "<<mult<<endl;
      //Get the parent main values
      */
      //Fill histograms
      hEnergy->Fill(energy) ;
      //      hEta->Fill(vpos.Eta()) ;
      //hPhi->Fill(vpos.Phi()) ;

    }// track loop 
  }//////////////////////////////////////////event loop
  hEnergy->Draw();//Draw histogram on screen
  //Write histograms in Root file 
  TFile outf("histo.root","recreate") ;
  hEnergy->Write() ;
  hPhi->Write() ;
  hEta->Write() ;
  outf.Close() ;
}


