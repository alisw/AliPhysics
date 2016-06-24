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
#include <TRefArray.h>
#include <TArrayS.h>

//AliRoot include files
#include "AliRunLoader.h"
#include "AliStack.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliESDCaloCluster.h"
#include "AliESDCaloCells.h"
#include "AliPID.h"
#include "AliLog.h" 

#endif

//example macro to extract information from AliESDCaloClusters and AliESDCaloCells.
// Author: Gustavo.Conesa.Balbastre@cern.ch (INFN-LNF)

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
//  .L TestESDCaloClusterAndCell.C++
//  TestESDCaloClusterAndCell(calo, number of events to process)
//=====================================================================
void TestESDCaloClusterAndCell(TString det = "EMCAL", const UInt_t eventsToProcess = 5, 
			TString dirName = "./", 
			const TString esdTreeName = "esdTree", 
			const char *  pattern = ".") 
{ 

  if(det !="EMCAL" && det !="PHOS" )
 {
   cout << "Wrong detector name "<<det<<endl; 
   return;
 }

  //Create chain of esd trees
  //By default the root files are in the same directory 
  TChain * t = AliReadESD(eventsToProcess, dirName,esdTreeName,pattern) ; 
 
  
  // ESD
  AliESDEvent * esd = new AliESDEvent(); 
  esd->ReadFromTree(t); // This checks also if the branch address is already set

  
  //Define few variables to be used in macro
  TString alirunName = "" ;
  UInt_t event ;
  Float_t pos[3] ; 

  for (event = 0; event < eventsToProcess; event++) {//event loop
    //AliInfo( Form("Event %d \n",event) );  
    Int_t nbytes = t->GetEntry(event); // store event in esd
    if ( nbytes == 0 ) //If nothing in ESD 
      break ; 
    
    // Check that name of file is correct
    if (alirunName != t->GetFile()->GetName()) {        
      alirunName = t->GetFile()->GetName() ; 
      alirunName.ReplaceAll("galice.root", "AliESDs.root") ;
    }
    
    AliRunLoader *rl = AliRunLoader::Open("galice.root",AliConfig::GetDefaultEventFolderName(),  "read");
    rl->LoadKinematics();  
    rl->GetEvent(event);
    AliStack *sta=rl->Stack();

    //get reconstructed vertex position     
    Double_t vertex_position[3] ; 
    esd->GetVertex()->GetXYZ(vertex_position) ; 
    
    //cout<<"Event >>>>>>>>>>>  "<<event<<" vertex >>> "<<vertex_position[0]<<" "<<vertex_position[1]<<" "<<vertex_position[2]<<endl;
    
    //Get CaloCells

    AliESDCaloCells &cells= *(esd->GetEMCALCells());
    if(det == "PHOS")  AliESDCaloCells &cells = *(esd->GetPHOSCells());
    Int_t ncell = cells.GetNumberOfCells() ;  
    Int_t type = cells.GetType();
    cout<<">>> Event "<<event<<" ncells "<<ncell<< " type "<<type<<endl;
    // Uncomment to see the full list of digit amplitudes and times.
//     for (Int_t icell=  0; icell <  ncell; icell++) {
//       cout<<"icell "<<icell<<
// 	" ID "<<cells.GetCellNumber(icell)<<
// 	" amp "<<cells.GetAmplitude(icell)<<
// 	" time "<<cells.GetTime(icell)<<endl;

//     }// cell loop

    //Get the CaloClusters

    //select EMCAL clusters only 
    TRefArray * caloClusters  = new TRefArray();

    if(det == "EMCAL") esd->GetEMCALClusters(caloClusters);
    else if(det == "PHOS") esd->GetPHOSClusters(caloClusters);
    
    Int_t nclus = caloClusters->GetEntries() ;  
   
    cout<<"Event: "<<event<<"; Number of clusters "<<nclus<<endl;   
    for (Int_t iclus =  0; iclus <  nclus; iclus++) {//////////////EMCAL cluster loop
      AliESDCaloCluster * clus = (AliESDCaloCluster *) caloClusters->At(iclus) ; // retrieve cluster from esd
      //cout<<">>>> Cluster >>>"<<iclus<<endl;
	
      //Get the cluster info
      
      Float_t energy   = clus->E() ;
      //Float_t disp     = clus->GetClusterDisp() ;
      Int_t iprim = clus->GetLabel();
      
      clus->GetPosition(pos) ;
      TVector3 vpos(pos[0],pos[1],pos[2]) ;
      TLorentzVector p ;
      clus->GetMomentum(p,vertex_position);
      
      Int_t mult       = clus->GetNCells() ;
      UShort_t * index = clus->GetCellsAbsId() ;
      Double_t * fraction = clus->GetCellsAmplitudeFraction() ;
      //Print cluster info
      cout<<"Cluster "<<iclus <<"; digits mult "<<mult<<"; type "<<(Int_t )clus->GetType()
	  <<"; Energy "<<energy
	  <<"; Phi "<<(vpos.Phi()*180/TMath::Pi())<<"; Eta "<<vpos.Eta()
	  <<"; label "<<iprim<<endl;      
      //Print primary info
      if(iprim>=0 && iprim < sta->GetNtrack()){
	TParticle * particle = sta->Particle(clus->GetLabel());
	//Print primary values
	cout<<" Primary: "<<particle->GetName()<< "; Energy "<<particle->Energy()<<endl;    
      }
      else if( iprim >= sta->GetNtrack()) cout <<"PROBLEM, label is too large : "<<iprim<<" >= particles in stack "<< sta->GetNtrack() <<endl;
      else cout<<"Negative label!!!  : "<<iprim<<endl;
      
      //Get CaloCells of cluster and print	 
      for(Int_t i = 0; i < mult ; i++){
	Int_t absId =   index[i]; // or clus->GetCellNumber(i) ;
	Double_t ampFract =  fraction[i];
	Float_t amp = cells.GetCellAmplitude(absId) ;
	Float_t time = cells.GetCellTime(absId);
	cout<<"AbsID : "<< absId << "; Amplitude "<< amp << "; Fraction "<<ampFract<<"; Time " <<time<<endl;
      }
      
    }// clusters
    
  }// event loop
  
}


