/// \file ana.C
/// \ingroup CaloTrackCorrMacros
/// \brief Example of execution macro
///
/// Example macro to do analysis with the
/// analysis classes in CaloTrackCorrelations,
/// in local, grid or plugin modes.
///
/// Pay attention to the options and definitions
/// set in the lines below
///
/// \author Gustavo Conesa Balbastre <Gustavo.Conesa.Balbastre@cern.ch>, (LPSC-CNRS)
///

#if !defined(__CINT__) || defined(__MAKECINT__)

R__ADD_INCLUDE_PATH($ALICE_ROOT)
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)

// ROOT
#include <Riostream.h>
#include <TString.h>
#include <TFile.h>
#include <TChain.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TApplication.h>
#include <TKey.h>
#include <TProfile.h>
#include <TH1F.h>

#include "TGrid.h"
#include "TGridCollection.h"
#include "TGridResult.h"
//#include "TProof.h"
//#include "TFileCollection.h"
//#include "TFileInfo.h"

// ALIROOT
#include "AliLog.h"
#include "AliAnalysisGrid.h"
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliAnalysisTaskSE.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODInputHandler.h"
#include "AliAODHandler.h"
#include "AliMultiInputEventHandler.h"
#include "AliAnalysisDataContainer.h"

// Main AddTasks and associated classes
#include "AliAnalysisTaskCaloTrackCorrelation.h"
#include "AliAnaCaloTrackCorrMaker.h"
#include "PWGGA/CaloTrackCorrelations/macros/AddTaskGammaHadronCorrelationSelectAnalysis.C" 
//#include "PWGGA/CaloTrackCorrelations/macros/AddTaskMultipleTrackCutIsoConeAnalysis.C" 
// comment, does not compile with AddTaskGammaHadronCorrelationSelectAnalysis.C
//#include "PWGGA/CaloTrackCorrelations/macros/AddTaskPi0IMGammaCorrQA.C" 
// comment, does not compile with AddTaskGammaHadronCorrelationSelectAnalysis.C

#include "AliPhysicsSelection.h"
#include "AliPhysicsSelectionTask.h"
#include "OADB/macros/AddTaskPhysicsSelection.C"

//#include "AliCentralitySelectionTask.h"
//#include "AddTaskCentrality.C"

#include "AliMultSelectionTask.h" 
#include "OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C"

//#include "AliVZEROEPSelectionTask.h"
//#include "AliEPSelectionTask.h"
//#include "AddTaskVZEROEPSelection.C"
//#include "ANALYSIS/macros/AddTaskEventplane.C"

//#include "CreateAlienHandler.C"

//#include "AliAnalysisTaskCounter.h"
//#include "PWGGA/CaloTrackCorrelations/macros/AddTaskCounter.C"

//#include "AliTender.h"
//#include "AliEmcalTenderTask.h"
//#include "AliEMCALTenderSupply.h"
//#include "AliEMCALRecParam.h"
//#include "PWG/EMCAL/macros/AddTaskEMCALTender.C"

#include "AliTaskCDBconnect.h"
#include "PWGPP/PilotTrain/AddTaskCDBconnect.C"

#include "AliAnalysisTaskEMCALClusterize.h"
#include "PWGPP/EMCAL/macros/AddTaskEMCALClusterize.C"

#include "AliEmcalCorrectionTask.h"
#include "PWG/EMCAL/macros/AddTaskEmcalCorrectionTask.C"

#endif

//---------------------------------------------------------------------------
/// Different analysis modes
enum anaModes
{
  mLocal  = 0, /// Analyze locally files in your computer.
  mPROOF  = 1, /// Analyze files on GRID with Plugin
  mPlugin = 2, /// Analyze files on GRID with Plugin
  mGRID   = 3  /// Analyze files on GRID, jobs launched from aliensh
};

//---------------------------------------------------------------------------
/// Settings to read locally several files, only for "mLocal" mode
/// The different values are default, they can be set with environmental 
/// variables: INDIR, PATTERN, NFILES, respectivelly

char * kInDir   = (char*)"/user/data/files/"; /// Directory path to files
char * kPattern = (char*)"";                  /// Common pattern in directory name containing files 
Int_t  kFile    = 10;                         /// Maximum number of files to analyze

//---------------------------------------------------------------------------
// Old PROOF settings, not used, here for historical reference
// Dataset for proof analysis, mode=mPROOF
// char * kDataset = (char*)"/alice/vernet/PbPb_LHC10h_ESD";
//
//char *  kDatasetPROOF     = (char*)"/alice/vernet/LHC11b_149646";
//Int_t   kDatasetNMaxFiles = 20;
//TString ccin2p3UserName   = "arbor" ;
//TString alienUserName     = "narbor" ;

//---------------------------------------------------------------------------
/// Collection file for grid analysis
char * kXML = (char*)"collection.xml"; /// Global name for the xml collection file with data on grid

//---------------------------------------------------------------------------
/// Scale histograms from file. Change to kTRUE when xsection file exists
/// Put name of file containing xsection 
/// Put number of events per ESD file
/// This is an specific case for normalization of Pythia files.
const char * kXSFileName = (char*)"pyxsec.root"; /// Name of file with pT-hard cross sections

// Container of cross section if it is in file pyxsec_hist.root
TArrayF* xsArr;
TArrayI* trArr;

/// Check during configuration the cross section, does not do anything to analysis
Bool_t bCheckXS = kFALSE; 

//---------------------------------------------------------------------------
///  Set some default values, but used values are set in the code!
Bool_t  kMC        = kFALSE; /// With real data kMC = kFALSE
TString kInputData = "ESD";  /// ESD, AOD, MC, deltaAOD
Int_t   kYear      = 2011;   /// Year of data
TString kPeriod    ="LHC11c";/// Run period
TString kCollision = "pp";   /// Collision type: pp, pPb, PbPb
Bool_t  outAOD     = kFALSE; /// Create output AOD, needed by some.
TString kTreeName;           /// "esdTree" or "aodTree" or "TE" for pure MC kinematics analysis    
TString kPass      = "";     /// "passX"
Int_t   kRun       = 0;      /// Run number

//---------------------------------------------------------------------------
/// Activate here the desired analysis combinations, what correction/clusterization.
Bool_t  bEMCCluster = kFALSE; /// Use the EMCal clusterization task 
Bool_t  bEMCCorrFra = kFALSE; /// Use the EMCal correction framework 
Int_t   xTalkEmul   = 0;      /// Activate cross-talk emulation, 0 -no, 1 do not subtract induced energy from reference cell, 2 subtract (preferred)

Bool_t  bAnalysis   = kTRUE;  /// Do photon/pi0 isolation correlation analysis 
//Bool_t  bAnalysisQA = kTRUE; /// Execute analysis QA train wagon, comment, does not compile with bAnalysis 
Bool_t  bMultiplicity= kFALSE;/// Execute multiplicity task 

//_________________________________
/// Load par files, create analysis libraries
/// For testing, if par file already decompressed and modified
/// classes then do not decompress.
//_________________________________
void SetupPar(char* pararchivename)
{
  TString cdir(Form("%s", gSystem->WorkingDirectory() )) ; 
  TString parpar(Form("%s.par", pararchivename)) ; 
  
  if ( gSystem->AccessPathName(pararchivename) ) 
  {  
    TString processline = Form(".! tar xvzf %s",parpar.Data()) ;
    gROOT->ProcessLine(processline.Data()); // REVIEW FOR ROOT6
  }
  
  TString ocwd = gSystem->WorkingDirectory();
  gSystem->ChangeDirectory(pararchivename);
  
  // check for BUILD.sh and execute
  if ( !gSystem->AccessPathName("PROOF-INF/BUILD.sh") ) 
  {
    printf("*******************************\n");
    printf("*** Building PAR archive    ***\n");
    cout<<pararchivename<<endl;
    printf("*******************************\n");
    
    if ( gSystem->Exec("PROOF-INF/BUILD.sh") ) 
    {
      Error("runProcess","Cannot Build the PAR Archive! - Abort!");
      return ;
    }
  }
  // check for SETUP.C and execute
  if ( !gSystem->AccessPathName("PROOF-INF/SETUP.C") )
  {
    printf("*******************************\n");
    printf("*** Setup PAR archive       ***\n");
    cout<<pararchivename<<endl;
    printf("*******************************\n");
    gROOT->Macro("PROOF-INF/SETUP.C"); // REVIEW FOR ROOT6
  }
  
  gSystem->ChangeDirectory(ocwd.Data());
  printf("Current dir: %s\n", ocwd.Data());
}

//______________________________________
/// Sets input data and tree strings.
//______________________________________
void CheckInputData(const anaModes mode)
{  
  TString ocwd = gSystem->WorkingDirectory();

  //---------------------------------------
  // Local files analysis
  //---------------------------------------
  if ( mode == mLocal )
  {    
    // If you want to add several ESD files sitting in a common directory INDIR
    // Specify as environmental variables the directory (INDIR), the number of files 
    // to analyze (NFILES) and the pattern name of the directories with files (PATTERN)
    
    if ( gSystem->Getenv("INDIR") )  
      kInDir = (char*)gSystem->Getenv("INDIR") ; 
    else 
      cout<<"INDIR not set, use default: "<<kInDir<<endl;  
    
    TString sindir(kInDir);
    if      ( sindir.Contains("pass1") ) kPass = "pass1";
    else if ( sindir.Contains("pass2") ) kPass = "pass2";
    else if ( sindir.Contains("pass3") ) kPass = "pass3";
    else if ( sindir.Contains("pass4") ) kPass = "pass4";
    
    if ( gSystem->Getenv("PATTERN") )   
      kPattern = (char*) gSystem->Getenv("PATTERN") ; 
    else  
      cout<<"PATTERN not set, use default: "<<kPattern<<endl;
    
    cout<<"INDIR   : "<<kInDir<<endl;
    cout<<"NFILES  : "<<kFile<<endl;
    
    char fileE[120] ;   
    char fileA[120] ;   
    char fileG[120] ;
    char fileEm[120] ;   
    for (Int_t event = 0 ; event < kFile ; event++) 
    {
      sprintf(fileE,  "%s/%s%d/AliESDs.root",    kInDir,kPattern,event) ; 
      sprintf(fileA,  "%s/%s%d/AliAOD.root",     kInDir,kPattern,event) ; 
      sprintf(fileG,  "%s/%s%d/galice.root",     kInDir,kPattern,event) ; 
      sprintf(fileEm, "%s/%s%d/embededAOD.root", kInDir,kPattern,event) ; 
      
      TFile * fESD = TFile::Open(fileE) ; 
      TFile * fAOD = TFile::Open(fileA) ; 
      
      // Check if file exists and add it, if not skip it
      if ( fESD ) 
      {
        kTreeName  = "esdTree";
        kInputData = "ESD";
        TFile * fG = TFile::Open(fileG);
        if(fG) { kMC = kTRUE; fG->Close();}
        else     kMC = kFALSE;
        
        // Get run number
        TTree* esdTree = (TTree*)fESD->Get("esdTree");
        AliESDEvent* esd = new AliESDEvent();
        esd->ReadFromTree(esdTree);
        esdTree->GetEvent(0);
        kRun = esd->GetRunNumber();
        
        return;
      }
      else if ( fAOD )
      {
        kTreeName  = "aodTree";
        kInputData = "AOD";
        if(((TTree*) fAOD->Get("aodTree"))->GetBranch("mcparticles")) kMC=kTRUE;
        else kMC = kFALSE;
        
        // Get run number
        TTree* aodTree = (TTree*)fAOD->Get("aodTree");
        AliAODEvent* aod = new AliAODEvent();
        aod->ReadFromTree(aodTree);
        aodTree->GetEvent(0);
        kRun = aod->GetRunNumber();
        return;
      }
      else if ( TFile::Open(fileEm) )
      {
        kTreeName  = "aodTree";
        kInputData = "AOD";
        kMC        = kTRUE;
        
        return;
      }
      else if ( TFile::Open(fileG) )
      {
        kTreeName  = "TE";
        kInputData = "MC";
        kMC        = kTRUE;
        return;
      }
      
      if ( fESD ) fESD->Close();
      if ( fAOD ) fAOD->Close();
    }
    
  }// local files analysis
  
  //------------------------------
  // GRID xml files
  //-----------------------------
  else if ( mode == mGRID )
  {
    // Get colection file. It is specified by the environmental
    // variable XML, if non provided, collection.xml is expected.
    if ( gSystem->Getenv("XML") )
      kXML = (char*) gSystem->Getenv("XML");
    
    if ( !TFile::Open(kXML) ) 
    {
      printf("No collection file with name -- %s -- was found\n",kXML);
      return ;
    }
    else 
      cout<<"XML file "<<kXML<<endl;
    
    // Load necessary libraries and connect to the GRID
    gSystem->Load("libNetx.so") ; 
    gSystem->Load("libRAliEn.so"); 
    TGrid::Connect("alien://") ;
        
    // Feed Grid with collection file
    TGridCollection * collection = gGrid->OpenCollection(kXML);
    if ( !collection )
    {
      printf("%s not found\n", kXML) ; 
      return ; 
    }
    
    TGridResult* result = collection->GetGridResult("",0 ,0);
    
    for (Int_t index = 0; index < result->GetEntries(); index++) 
    {
      TString alienURL = result->GetKey(index, "turl") ; 
      cout << "================== " << alienURL << endl ; 
      
      if     (alienURL.Contains("pass1")) kPass = "pass1";
      else if(alienURL.Contains("pass2")) kPass = "pass2";
      else if(alienURL.Contains("pass3")) kPass = "pass3";
      else if(alienURL.Contains("pass4")) kPass = "pass4";
      
      kRun = AliAnalysisManager::GetRunFromAlienPath(alienURL.Data());
      printf("Run number from alien path = %d\n",kRun);
      
      TFile * fAOD = 0 ; 
      // Check if file exists and add it, if not skip it
      if ( alienURL.Contains("AliESDs.root") )  
      {
        kTreeName  = "esdTree";
        kInputData = "ESD";
        alienURL.ReplaceAll("AliESDs.root","galice.root");
        if(TFile::Open(alienURL)) kMC=kTRUE;
        else kMC = kFALSE;
        return;
      }
      else if ( alienURL.Contains("AliAOD.root") )
      {
        kTreeName  = "aodTree";
        kInputData = "AOD";
        fAOD = TFile::Open(alienURL);
        if(((TTree*) fAOD->Get("aodTree"))->GetBranch("mcparticles")) kMC=kTRUE;
        else kMC = kFALSE;
        return;
      }
      else if ( alienURL.Contains("embededAOD.root") )
      {
        kTreeName  = "aodTree";
        kInputData = "AOD";
        kMC=kTRUE;
        return;
      }
      else if ( alienURL.Contains("galice.root") )
      {
        kTreeName  = "TE";
        kInputData = "MC";
        kMC=kTRUE;
        return;
      } 
    }
  }// xml analysis
   //------------------------------
   //PROOF files
   //-----------------------------
//  else if(mode == mPROOF)
//  {
//    
//    TFileCollection* coll  = gProof->GetDataSet(kDatasetPROOF)->GetStagedSubset();
//    
//    TIter iter(coll->GetList());
//    
//    TFileInfo* fileInfo = 0;
//    while ((fileInfo = dynamic_cast<TFileInfo*> (iter())))
//    {
//      if (fileInfo->GetFirstUrl()) 
//  {
//        TString ProofURL = fileInfo->GetFirstUrl()->GetUrl();
//        cout << "================== " << ProofURL << endl ; 
//        
//        if     (ProofURL.Contains("pass1")) kPass = "pass1";
//        else if(ProofURL.Contains("pass2")) kPass = "pass2";
//        else if(ProofURL.Contains("pass3")) kPass = "pass3";
//        
//        kRun = AliAnalysisManager::GetRunFromAlienPath(ProofURL.Data());
//        printf("Run number from alien path = %d\n",kRun);
//        
//        TFile * fAOD = 0 ; 
//        //Check if file exists and add it, if not skip it
//        if (ProofURL.Contains("AliESDs.root"))  
//        {
//          kTreeName  = "esdTree";
//          kInputData = "ESD";
//          alienURL.ReplaceAll("AliESDs.root","galice.root");
//          if(TFile::Open(ProofURL)) kMC=kTRUE;
//          else kMC = kFALSE;
//          
//          return;
//        }
//        else if(ProofURL.Contains("AliAOD.root"))
//        {
//          kTreeName  = "aodTree";
//          kInputData = "AOD";
//          fAOD = TFile::Open(ProofURL);
//          if(((TTree*) fAOD->Get("aodTree"))->GetBranch("mcparticles")) kMC=kTRUE;
//          else kMC = kFALSE;
//          return;
//        }
//        else if(ProofURL.Contains("embededAOD.root"))
//        {
//          kTreeName  = "aodTree";
//          kInputData = "AOD";
//          kMC=kTRUE;
//          return;
//        }
//        else if(ProofURL.Contains("galice.root"))
//        {
//          kTreeName  = "TE";
//          kInputData = "MC";
//          kMC=kTRUE;
//          return;
//        } 
//      }
//    }
//  }// proof analysis
  
  gSystem->ChangeDirectory(ocwd.Data());
}

//_____________________________________________________________________
/// Fills chain with data files paths.
//_____________________________________________________________________
void CreateChain(const anaModes mode, TChain * chain, TChain * chainxs)
{
  TString ocwd = gSystem->WorkingDirectory();
  
  if ( kInputData == "AOD" )
  {
    xsArr = new TArrayF;
    trArr = new TArrayI;
  } 
  
  //---------------------------------------
  // Local files analysis
  //---------------------------------------
  if ( mode == mLocal )
  {    
    // If you want to add several ESD files sitting in a common directory INDIR
    // Specify as environmental variables the directory (INDIR), the number of files 
    // to analyze (NFILES) and the pattern name of the directories with files (PATTERN)
    
    if ( gSystem->Getenv("INDIR") )  
      kInDir = (char*)gSystem->Getenv("INDIR") ; 
    else 
      cout<<"INDIR not set, use default: "<<kInDir<<endl;  
    
    if ( gSystem->Getenv("PATTERN") )   
      kPattern = (char*) gSystem->Getenv("PATTERN") ; 
    else  
      cout<<"PATTERN not set, use default: "<<kPattern<<endl;
    
    if ( gSystem->Getenv("NFILES") )
      kFile = atoi(gSystem->Getenv("NFILES")) ;
    else 
      cout<<"NFILES not set, use default: "<<kFile<<endl;
    
    // Check if env variables are set and are correct
    if ( kInDir  && kFile) 
    {
      printf("Get %d files from directory %s\n",kFile,kInDir);
      // Check if ESDs directory exist
      if ( ! gSystem->cd(kInDir) ) 
      {
        printf("%s does not exist\n", kInDir) ;
        return ;
      }
      
      // if(gSystem->Getenv("XSFILE"))  
      // kXSFileName = gSystem->Getenv("XSFILE") ; 
      // else cout<<" XS file name not set, use default: "<<kXSFileName<<endl; 
      
      char * kGener = (char*) gSystem->Getenv("GENER");
      if ( kGener ) 
      {
        cout<<"GENER "<<kGener<<endl;
        if     (!strcmp(kGener,"PYTHIA")) kXSFileName = "pyxsec.root";
        else if(!strcmp(kGener,"HERWIG")) kXSFileName = "hexsec.root";
        else cout<<" UNKNOWN GENER, use default: "<<kXSFileName<<endl;
      }
      else 
        cout<<" GENER not set, use default xs file name: "<<kXSFileName<<endl;
      
      if ( kInputData == "AOD" )
      {
        kXSFileName = "pyxsec_hists.root";
      }      
      
      cout<<"INDIR   : "<<kInDir     <<endl;
      cout<<"NFILES  : "<<kFile      <<endl;
      cout<<"PATTERN : "<<kPattern   <<endl;
      cout<<"XSFILE  : "<<kXSFileName<<endl;
      
      TString datafile="";
      if     (kInputData == "ESD")        datafile = "AliESDs.root" ;
      else if(kInputData.Contains("AOD")) datafile = "AliAOD.root"  ;
      else if(kInputData == "MC")         datafile = "galice.root"  ;
      
      //Loop on ESD/AOD/MC files, add them to chain
      Int_t event =0;
      Int_t skipped=0 ; 
      char file[120] ;
      char filexs[120] ;
      
      if ( bCheckXS )
      {
        xsArr->Set(kFile);
        trArr->Set(kFile);
      }
      
      for (event = 0 ; event < kFile ; event++) 
      {
        sprintf(file,   "%s/%s%d/%s", kInDir,kPattern,event,datafile.Data()) ; 
        if ( bCheckXS )
          sprintf(filexs, "%s/%s%d/%s", kInDir,kPattern,event,kXSFileName) ;
        
        TFile * fData = TFile::Open(file) ; 
        // Check if file exists and add it, if not skip it
        if ( fData ) 
        {
          if ( fData->Get(kTreeName) ) 
          { 
            printf("++++ Adding %s\n", file) ;
            chain->AddFile(file);
            
            if ( bCheckXS )
            {
              if ( kInputData != "AOD" )
              {
                chainxs->Add(filexs) ;
              }
              else
              {
                TFile*  fxsec = TFile::Open(filexs);
                if(fxsec)
                {
                  TKey* key = (TKey*)fxsec->GetListOfKeys()->At(0);
                  if(!key)
                  {
                    fxsec->Close();
                    printf("No key!");
                    continue;
                  }
                  
                  TList *list = dynamic_cast<TList*>(key->ReadObj());
                  if(!list)
                  {
                    fxsec->Close();
                    printf("No list!");
                    continue;
                  }
                  
                  Float_t xsection = ((TProfile*)list->FindObject("h1Xsec"))  ->GetBinContent(1);
                  Int_t   ntrials  = ((TH1F*)    list->FindObject("h1Trials"))->GetBinContent(1);
                  fxsec->Close();
                  
                  xsArr->SetAt(xsection,event);
                  trArr->SetAt(ntrials,event);
                  
                  printf("recovered xs %f, ntrials %d, event %d\n",xsection,ntrials, event);
                  //chainxs->Add(tree);
                  //fileTMP->Close();
                } // fxsec exists
              } // xs in AODs
            } // check XS
            
          } // data tree chake
        }
        else 
        { 
          printf("---- Skipping %s\n", file) ;
          skipped++ ;
        }
      }
    }
    else 
    {
      TString input = "AliESDs.root" ;
      cout<<">>>>>> No list added, take a single file <<<<<<<<< "<<input<<endl;
      chain->AddFile(input);
    }
    
  }// local files analysis
  
  //------------------------------
  // GRID xml files
  //------------------------------
  else if ( mode == mGRID )
  {
    // Get colection file. It is specified by the environmental
    // variable XML
    
    // Feed Grid with collection file
    TGridCollection * collection = gGrid->OpenCollection(kXML);
    if (! collection ) 
    {
      printf("%s not found \n", kXML) ; 
      return ; 
    }
    
    TGridResult* result = collection->GetGridResult("",0 ,0);
    
    // Makes the chain 
    Int_t nFiles = result->GetEntries();
    Int_t nXSFilesFound = 0;
    if ( bCheckXS )
    {
      xsArr->Set(nFiles);
      trArr->Set(nFiles);
    }
    
    printf("*** Filling the Chain with %d files***\n",nFiles);
    for (Int_t index = 0; index < nFiles; index++) 
    {
      TString alienURL = result->GetKey(index, "turl") ; 
      cout << "================== " << alienURL << endl ; 
      chain->Add(alienURL) ; 
      
      if ( bCheckXS )
      {
        if ( kInputData != "AOD" )
        {
          alienURL.ReplaceAll("AliESDs.root",kXSFileName);
          alienURL.ReplaceAll("AliAOD.root" ,kXSFileName);
          chainxs->Add(alienURL) ;
        }
        else
        {
          alienURL.ReplaceAll("AliESDs.root","pyxsec_hists.root");
          alienURL.ReplaceAll("AliAOD.root", "pyxsec_hists.root");
          
          TFile*  fxsec = TFile::Open(alienURL);
          if ( fxsec )
          {
            TKey* key = (TKey*)fxsec->GetListOfKeys()->At(0);
            if(!key)
            {
              fxsec->Close();
              printf("No key!");
              continue;
            }
            
            TList *list = dynamic_cast<TList*>(key->ReadObj());
            if ( !list )
            {
              fxsec->Close();
              printf("No list!");
              continue;
            }
            
            Float_t xsection = ((TProfile*)list->FindObject("h1Xsec"))  ->GetBinContent(1);
            Int_t   ntrials  = ((TH1F*)    list->FindObject("h1Trials"))->GetBinContent(1);
            fxsec->Close();
            
            xsArr->SetAt(xsection,index);
            trArr->SetAt(ntrials,index);
            nXSFilesFound++;
            printf("recovered xs %f, ntrials %d, index %d\n",xsection,ntrials, index);
            
          } // fxsec exists
          
        } // xs in AODs
      }
      
    } // loop files
    
    if ( bCheckXS )
    {
      xsArr->Set(nXSFilesFound);
      trArr->Set(nXSFilesFound);
    }
    
  }// xml analysis
  
  //------------------------------
  // PROOF
  //------------------------------
//  else if (mode == mPROOF) 
//  {
//    
//    TFileCollection* ds= gProof->GetDataSet(kDatasetPROOF)->GetStagedSubset();
//    
//#if defined(__CINT__)
//    gROOT->LoadMacro("/afs/in2p3.fr/group/alice/laf/dataset_management/CreateChainFromDataSet.C");
//#endif
//    chain = CreateChainFromDataSet(ds, kTreeName , kDatasetNMaxFiles);
//    printf("chain has %d entries\n",chain->GetEntries());
//  }
  
  gSystem->ChangeDirectory(ocwd.Data());
}

//______________________________
/// Access one data file and set the year,
/// collision type and run number.
/// It is possible to set them via external parameters.
//______________________________
void CheckEnvironmentVariables()
{
  Bool_t bRecalibrate = kFALSE;
  Bool_t bBadChannel  = kFALSE;
  
  TString sRun = "";

  for (int i=0; i< gApplication->Argc();i++)
  {
#ifdef VERBOSEARGS
    printf("Arg  %d:  %s\n",i,gApplication->Argv(i));
#endif
    
    sRun = "";
    
    if (!(strcmp(gApplication->Argv(i),"--recalibrate")))
      bRecalibrate = atoi(gApplication->Argv(i+1));
    
    if (!(strcmp(gApplication->Argv(i),"--badchannel")))
      bBadChannel = atoi(gApplication->Argv(i+1));

    if (!(strcmp(gApplication->Argv(i),"--year")))
      kYear = atoi(gApplication->Argv(i+1));
    
    if (!(strcmp(gApplication->Argv(i),"--run")))
    {
      sRun = gApplication->Argv(i+1);
      if ( sRun.Contains("LHC10") ) 
      {
        kYear = 2010;
      }
      else
      {
        if ( kRun <=0 )
        {
          kRun = atoi(gApplication->Argv(i+1));
        }
        else 
          printf("** Run number already set  to %d, do not set to %d\n",kRun,atoi(gApplication->Argv(i+1)));
      } // numeric run
    } // --run available
    
  }// args loop
  
  // Check run number and decide kYear and kCollision, by default kCollision = "pp"
  if ( !sRun.Contains("LHC") )
  {
    if     ( kRun < 140000) 
    {
      kYear   = 2010;
      kPeriod = "LHC10";
      if( kRun >= 136851 ) 
      {
        kPeriod    = "LHC10h";
        kCollision = "PbPb";
      }
    }
    else if( kRun < 170600)
    {
      kYear   = 2011;
      kPeriod = "LHC11";
      if( kRun >= 166500 ) 
      {
        kPeriod    = "LHC11h";
        kCollision = "PbPb";
      }
    }
    else if( kRun < 200000 )
    {
      kYear   = 2012;
      kPeriod = "LHC12";
      if( kRun >= 194000 ) 
      {
        kPeriod    = "LHC13";
        kCollision = "pPb";
      }
    }
    else if( kRun < 247000 )
    {
      kYear   = 2015;
      kPeriod = "LHC15";
      if( kRun >= 244820 ) 
      {
        kPeriod    = "LHC15o";
        kCollision = "PbPb";
      }
    }
    else if( kRun < 268875 )
    {
      kYear   = 2016;
      kPeriod = "LHC16";
      if( kRun >= 265015 ) 
      {
        kPeriod    = "LHC16q";
        kCollision = "pPb";
      }
    }
    else if( kRun < 283616 )
    {
      kYear   = 2017;
      kPeriod = "LHC17";
      if( kRun == 280235 || kRun == 280234 ) 
      {
        kPeriod    = "LHC17n";
        kCollision = "PbPb"; // XeXe
      }
    }
    else
    {
      kYear   = 2018; 
      kPeriod = "LHC18";
      // To be defined
      //if( kRun >= XXXX )
      //{
      //  kPeriod    = "LHC18X";
      //  kCollision = "PbPb";
      //}
    }
  }
    
  printf("*********************************************\n");
  //printf("*** Settings recalibrate %d, remove bad channels %d, year %d, collision %s, run %d ***\n",
  //       bRecalibrate,bBadChannel, kYear,kCollision.Data(), kRun);
  printf("*** Settings year %d, collision %s, run %d ***\n",kYear,kCollision.Data(), kRun);
  printf("*********************************************\n");
  
}

//______________________________________________________________________________
/// Read the PYTHIA statistics from the file pyxsec.root created by
/// the function WriteXsection():
/// integrated cross section (xsection) and
/// the  number of Pyevent() calls (ntrials)
/// and calculate the weight per one event xsection/ntrials
/// The spectrum calculated by a user should be
/// multiplied by this weight, something like this:
/// TH1F *userSpectrum ... // book and fill the spectrum
/// userSpectrum->Scale(weight)
///
/// Yuri Kharlov 19 June 2007
/// Gustavo Conesa 15 April 2008
/// Add recovery of xs from pyxsec_hists.root file 15/jan/2015
//______________________________________________________________________________
Bool_t GetAverageXsection(TTree * tree, Double_t & xs, Float_t & ntr, Int_t & n)
{
  Double_t xsection = 0 ;
  UInt_t    ntrials = 0 ;
  Int_t      nfiles = 0 ;
    
  xs  = 0;
  ntr = 0;
  n   = 0;
  
  if( kInputData != "AOD" &&  tree )
  {
    nfiles =  tree->GetEntries()  ;

    tree->SetBranchAddress("xsection",&xsection);
    tree->SetBranchAddress("ntrials" ,&ntrials );
    for(Int_t i = 0; i < nfiles; i++)
    {
      tree->GetEntry(i);
      if(xsection > 0)
      {
        xs  += xsection ;
        ntr += ntrials ;
        n++;
      }
      printf("\t i %d xsection %e, trials %d\n",i, xsection,ntrials);
    } // loop
  }
  else if( kInputData == "AOD" && xsArr )
  {

    nfiles = xsArr->GetSize();
    
    for(Int_t i = 0; i < nfiles; i++)
    {
      if(xsArr->GetAt(i) > 0)
      {
        xs  += xsArr->GetAt(i) ;
        ntr += trArr->GetAt(i) ;
        n++;
      }
      printf("\t i %d xsection %e, trials %f\n",i, xsArr->GetAt(i),trArr->GetAt(i));
    } // loop
  }
  else return kFALSE;
  
  printf("\t total xs %e, ntr %e, n used files %d, n total files %d\n",xs,ntr,n,nfiles);
    
  if ( n <= 0 ) 
  {
    printf("CAREFUL: No files %d\n",n);
    return kFALSE;
  }
  
  xs =   xs /  n;
  ntr =  ntr / n;
  cout << "-----------------------------------------------------------------"<<endl;
  cout << "Average of "<< n <<" files: xsection " <<xs<<" ntrials "<<ntr<<endl;
  cout << "-----------------------------------------------------------------"<<endl;
  
  return kTRUE;
}

//_____________________________
/// Load analysis libraries.
/// Add here the call to the SetUp par files if need to work
/// on grid with modified code.
//_____________________________
void  LoadLibraries(Int_t /*mode*/)
{
  //  if (mode == mPROOF)
  //  {
  //    //TProof::Mgr("ccalpmaster")->SetROOTVersion("ALICE_v5-27-06b");
  //#if defined(__CINT__)
  //    gROOT->LoadMacro("/afs/in2p3.fr/group/alice/laf/EnableAliRootForLAF.C");
  //#endif
  //    TProof* proof = EnableAliRootForLAF("ccaplmaster",nPROOFWorkers.Data(),ccin2p3UserName.Data(),alienUserName.Data(),"",kFALSE,kTRUE,kTRUE,"OADB:ANALYSIS:ANALYSISalice:AOD:ESD:CORRFW:STEERBase:EMCALUtils:PHOSUtils:PWGCaloTrackCorrBase:PWGGACaloTrackCorrelations:PWGPPEMCAL");
  //    
  //    //  TProof* proof = TProof::Open("ccaplmaster",Form("workers=%s",nPROOFWorkers.Data()));
  //    
  //    //     //proof->ClearPackages();
  //    //     proof->UploadPackage("STEERBase");
  //    //     proof->UploadPackage("ESD");
  //    //     proof->UploadPackage("AOD");
  //    //     proof->UploadPackage("ANALYSIS");
  //    //     proof->UploadPackage("OADB");
  //    //     proof->UploadPackage("ANALYSISalice");
  //    //     proof->UploadPackage("CORRFW");
  //    //     //proof->UploadPackage("JETAN");
  //    //     proof->UploadPackage("PHOSUtils");
  //    //     proof->UploadPackage("EMCALUtils");
  //    //     proof->UploadPackage("PWGCaloTrackCorrBase");
  //    //     proof->UploadPackage("PWGGACaloTrackCorrelations");
  //    //     proof->UploadPackage("PWGPPEMCAL");
  //    
  //    //     proof->EnablePackage("STEERBase");
  //    //     proof->EnablePackage("ESD");
  //    //     proof->EnablePackage("AOD");
  //    //     proof->EnablePackage("ANALYSIS");
  //    //     proof->EnablePackage("OADB");
  //    //     proof->EnablePackage("ANALYSISalice");
  //    //     proof->EnablePackage("CORRFW");
  //    //     //proof->EnablePackage("JETAN");
  //    //     proof->EnablePackage("PHOSUtils");
  //    //     proof->EnablePackage("EMCALUtils");
  //    //     proof->EnablePackage("PWGCaloTrackCorrBase");
  //    //     proof->EnablePackage("PWGGACaloTrackCorrelations");
  //    //     proof->EnablePackage("PWGPPEMCAL");
  //    return;
  //  }  
  
  //--------------------------------------
  // Load the needed libraries via par files if modified
  //--------------------------------------
  //SetupPar("EMCALUtils");
  //SetupPar("EMCALraw");
  //SetupPar("EMCALbase");
  //SetupPar("EMCALsim");
  //SetupPar("EMCALrec");
  
  //SetupPar("PWGPPEMCAL");
  
  //SetupPar("PWGCaloTrackCorrBase");
  //SetupPar("PWGGACaloTrackCorrelations"); 
  
  //SetupPar("PWGGAGammaConv"); 
  
  // needed for plugin?
  gSystem->AddIncludePath("-I$ALICE_ROOT");
  gSystem->AddIncludePath("-I$ALICE_PHYSICS");
  gSystem->AddIncludePath("-I./");     
}

//___________________________
/// Main execution method. It:
/// * 1) loads the needed libraries in method LoadLibraries
/// * 2) depending on the files path, run etc, the variables year, collision type, data type, are obtained in methods CheckInputData and CheckEnvironmentVariables
/// * 3) put the data files in a list to be passed to the analysis frame in method CreateChain
/// * 4) In case of MC pt-Hard bin simulations, the file containing the cross sections is read and scaling parameter is obtained via the method GetAverageXsection
/// * 5) The analysis frame is initialized via de analysis manager
/// * 6) Different general analysis are initialized: Physics selection, centrality etc.
/// * 7) Specialized analysis are initialized: AliAnalysistaskCounter, AliAnalysisTaskEMCALClusterizer/EMCal correction framework, 
///      AliAnalysisTaskCaloTrackCorrelations and executed for different settings.
/// * 8) The output/input containers are passed to the analysis manager
/// * 9) The analysis is executed
///
/// \param mode: analysis mode defined in enum anaModes
//________________________
void ana ( anaModes mode = mGRID )
{  
  //--------------------------------------------------------------------
  // Load analysis libraries
  LoadLibraries(mode) ;
  //gSystem->ListLibraries();
  
  //-----------------------------------------------------------------------------
  // Create chain from ESD and from cross sections files, look below for options.
  
  // Set kInputData and kTreeName looking to the kINDIR
  CheckInputData(mode);
  
  // Check global analysis settings  
  CheckEnvironmentVariables();
  
  printf("*********************************************\n");
  printf("*** Input data < %s >, pass %s, tree < %s >, MC?  < %d > ***\n",
         kInputData.Data(),kPass.Data(),kTreeName.Data(),kMC);
  printf("*********************************************\n");
  
  TChain * chain   = new TChain(kTreeName) ;
  TChain * chainxs = new TChain("Xsection") ;
  CreateChain(mode, chain, chainxs); 
  
  if ( !chain )
  { 
    printf("STOP, no chain available\n"); 
    return;
  } 
  
  printf("*********************************************\n");
  printf("number of entries in chain # %lld \n", chain->GetEntries()) ;   
  printf("*********************************************\n");
  
  // Recover the cross section and print average value
  if ( kMC && bCheckXS )
  {
    //Get the cross section
    Double_t scale    = -1;
    Double_t xsection = 0;
    Float_t  ntrials  = 0;
    Int_t    nfiles   = 0;
    printf("===== kMC %d, chainxs %p\n",kMC,chainxs);

    printf("Average Cross section:\n");
    Bool_t ok = GetAverageXsection(chainxs, xsection, ntrials, nfiles);
    
    printf("\t ok %d n xs files %d ntrials %f \n",ok, nfiles, ntrials);
    if ( nfiles > 0 && ntrials > 0 )
    {
      if ( ok )
      {
        Int_t  nEventsPerFile = chain->GetEntries() / nfiles;
        
        Double_t trials = ntrials / nEventsPerFile ;
        
        scale = xsection / trials;
        
        printf("\t Get Cross section : nfiles  %d, nevents %lld, nevents per file %d \n",
               nfiles, chain->GetEntries(),nEventsPerFile);
        printf("\t \t ntrials %2.2f, trials %2.2f, xs %2.2e, scale factor %2.2e\n", 
               ntrials,trials,xsection,scale);
        
        if ( nfiles != chain->GetNtrees() ) 
          printf("\t CAREFUL: Number of files in data chain %d, in cross section chain %d \n",
                 chain->GetNtrees(),nfiles);
      } // ok
      
      // comment out this line in case the simulation did not have the 
      // cross section files produced in the directory
      if ( scale <= 0  || !ok )
      { 
        printf( "\t STOP, cross section not available! nfiles %lld \n", 
               chainxs->GetEntries() ) ; 
        return ; 
      }
    }
  }
    
  AliLog::SetGlobalLogLevel(AliLog::kError);//Minimum prints on screen
  
  //------------------------------------------
  //  Alien handler part
  //------------------------------------------
//  AliAnalysisGrid * alienHandler =0x0;
//  if ( mode==mPlugin )
//  {
//    // Create and configure the alien handler plugin
//#if defined(__CINT__)
//    gROOT->LoadMacro("CreateAlienHandler.C");
//#endif
//    alienHandler = CreateAlienHandler();
//    if ( !alienHandler ) return;
//  }  
  
  //--------------------------------------
  //--------------------------------------
  // If automatic check does not work, 
  // force here the collision type, period etc
//  kYear      = 2011;
//  kPeriod    = "LHC11c";
//  kCollision = "pp";
//  kMC        = kTRUE;
//  kInputData = "AOD";
  //--------------------------------------
  //--------------------------------------

  //--------------------------------------
  // Make the analysis manager
  //-------------------------------------
  AliAnalysisManager *mgr  = new AliAnalysisManager("Manager", "Manager");
  //AliAnalysisManager::SetUseProgressBar(kTRUE);
  //mgr->SetSkipTerminate(kTRUE);
  //mgr->SetNSysInfo(1);
  
  // Connect plugin to the analysis manager
//  if ( mode == mPlugin )
//  {
//    mgr->SetGridHandler(alienHandler);
//  }
  
  // MC handler
  if ( (kMC || kInputData == "MC") && !kInputData.Contains("AOD") )
  {
    AliMCEventHandler* mcHandler = new AliMCEventHandler();
    mcHandler->SetReadTR(kFALSE);//Do not search TrackRef file
    mgr->SetMCtruthEventHandler(mcHandler);
    if ( kInputData == "MC" ) 
    {
      cout<<"MC INPUT EVENT HANDLER"<<endl;
      mgr->SetInputEventHandler(NULL);
    }
  }
  
  // AOD output handler, very special analysis
  if ( kInputData != "deltaAOD" && outAOD)
  {
    cout<<"Init output handler"<<endl;
    AliAODHandler* aodoutHandler   = new AliAODHandler();
    aodoutHandler->SetOutputFileName("aod.root");
    ////aodoutHandler->SetCreateNonStandardAOD();
    mgr->SetOutputEventHandler(aodoutHandler);
  }
  
  //=========
  // Input
  
  if ( kInputData == "ESD" )
  {
    // ESD handler
    AliESDInputHandler *esdHandler = new AliESDInputHandler();
    esdHandler->SetReadFriends(kFALSE);
    mgr->SetInputEventHandler(esdHandler);
    cout<<"ESD handler "<<mgr->GetInputEventHandler()<<endl;
  }
  else if ( kInputData.Contains("AOD") )
  {
    // AOD handler
    AliAODInputHandler *aodHandler = new AliAODInputHandler();
    mgr->SetInputEventHandler(aodHandler);
    if(kInputData == "deltaAOD") aodHandler->AddFriend((char*)"deltaAODCaloTrackCorr.root");
    cout<<"AOD handler "<<mgr->GetInputEventHandler()<<endl;
  }
  
  //mgr->RegisterExternalFile("deltaAODCaloTrackCorr.root");
  //mgr->SetDebugLevel(1); // For debugging, do not uncomment if you want no messages.
  
  TString outputFile = AliAnalysisManager::GetCommonFileName(); 
  
  //-------------------------------------------------------------------------
  // Define task, put here any other task that you want to use.
  //-------------------------------------------------------------------------
  
  // Physics selection
  if ( !kMC )
  {
#if defined(__CINT__)
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C"); 
#endif
    AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(kMC,kTRUE); 
  }
  
  // Centrality, valid for Run1, but superseeded by new task below
//  if ( bMultiplicity && kCollision.Contains("Pb") )
//  {
//    if ( kYear < 200000 && kInputData=="ESD" )
//    {
//
//#if defined(__CINT__)
//     gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
//#endif
//     AliCentralitySelectionTask *taskCentrality = AddTaskCentrality();
//    }
//  }
  
  if ( bMultiplicity )
  {
    // New centrality/multiplicity selector
#if defined(__CINT__)
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
#endif
    AliMultSelectionTask * task = AddTaskMultSelection(kFALSE); // user mode:
    
    // use the default calibration for runs which have not yet been calibrated
    task->SetUseDefaultCalib  (kTRUE); // data
    task->SetUseDefaultMCCalib(kTRUE); // MC
  }
  
//  if ( kCollision=="PbPb" )
//  {
//#if defined(__CINT__)
//    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskVZEROEPSelection.C");
//    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskEventplane.C");
//#endif
//
//    AliVZEROEPSelectionTask  * EPV0 = AddTaskVZEROEPSelection();  
//    
//    AliEPSelectionTask * EP = AddTaskEventplane();
//  } 
  
  // OCDB connect
  //
  if ( bEMCCorrFra || bEMCCluster )
  {
#if defined(__CINT__)
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/PilotTrain/AddTaskCDBconnect.C");
    //gROOT->LoadMacro("AddTaskCDBconnect.C");
#endif

    AddTaskCDBconnect();
    ((AliTaskCDBconnect*)(AliAnalysisManager::GetAnalysisManager()->GetTask("CDBconnect")))->SetFallBackToRaw(kTRUE);
  }
  
  // EMCAL correction framework
  //
  if ( bEMCCorrFra && !bEMCCluster )
  {
    printf("INIT EMCal corrections\n");
#if defined(__CINT__)
    gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalCorrectionTask.C");
#endif

    AliEmcalCorrectionTask * emcorr = AddTaskEmcalCorrectionTask();
    
    //emcorr->SetUserConfigurationFilename("./EMCalCorrConfig_Data_ClV1_Run2TCalib_Test.yaml");
    // Data or MC specific configurations
    if ( !kMC )
    {
      emcorr->SetUserConfigurationFilename("$ALICE_PHYSICS_SRC/PWGGA/CaloTrackCorrelations/yaml/EMCalCorrConfig_Gamma_Data.yaml");
    }
    else
    {
      // Without cross-talk
      if ( xTalkEmul  == 0 )
        emcorr->SetUserConfigurationFilename("$ALICE_PHYSICS_SRC/PWGGA/CaloTrackCorrelations/yaml/EMCalCorrConfig_MC_ClV1.yaml");
      // With cross-talk
      else 
      {
        //emcorr->SetUserConfigurationFilename("$ALICE_PHYSICS_SRC/PWGGA/CaloTrackCorrelations/yaml/EMCalCorrConfig_MC_Run1_ClV1_xTalk.yaml");
        emcorr->SetUserConfigurationFilename("$ALICE_PHYSICS_SRC/PWGGA/CaloTrackCorrelations/yaml/EMCalCorrConfig_MC_Run1_ClV1_xTalk_ECellCut.yaml");
        //emcorr->SetUserConfigurationFilename("$ALICE_PHYSICS_SRC/PWGGA/CaloTrackCorrelations/yaml/EMCalCorrConfig_MC_Run1_ClV1_xTalk_ECellCut_Leak5MeV.yaml");
      }
    }
    
    //emcorr->SelectCollisionCandidates( AliVEvent::kAnyINT | AliVEvent::kEMC7 | AliVEvent::kEMC8 | AliVEvent::kEMC1 | AliVEvent::kEMCEGA | AliVEvent::kEMCEJE );
    
    emcorr->SetForceBeamType(AliEmcalCorrectionTask::kpp);  
    emcorr->Initialize();
  }
  
  // Clusterization task
  // For experts
  //
  TString clustersArray = "";
  TString cellsArray    = "";
  if ( !bEMCCorrFra && bEMCCluster )
  {
    printf("INIT EMCal Clusterizer\n");
#if defined(__CINT__)
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/EMCAL/macros/AddTaskEMCALClusterize.C"); 
#endif
 
    TString sClust    = "V1Unfold"; // Options: V1, V2, V1Unfold, NxN
    Int_t   clTM      = 2;       // Do track matching, 0 no, 1 TPC only, 2 hybrid
    Bool_t  exo       = kTRUE;   // Remove exotic cells
    Bool_t  clnonlin  = kTRUE;   // Apply non linearity (clusterizer), CAREFUL check that not done in analysis
    Int_t   minEcell  = 100;     // 50  MeV (10 MeV used in reconstruction)
    Int_t   minEseed  = 500;     // 100 MeV
    Int_t   dTime     = 10000;   // open
    Int_t   wTime     = 10000;   // open
    Int_t   unfMinE   = 15;      // Remove cells with less than 15 MeV from cluster after unfolding
    Int_t   unfFrac   = 1;       // Remove cells with less than 1% of cluster energy after unfolding
    Bool_t  updateCell= kFALSE;  // Calibrate cells and modify them on the fly
    Bool_t  filterEvents = kFALSE; // Filter events with activity in EMCal
    Int_t   cenBin[]  = {-1,-1}; // Centrality bin min-max of accepted events. {-1,-1} take all
    // Calibration, bad map ...
    
    Bool_t calibEE = kTRUE; // It is set automatically, but here we force to use ir or not in any case
    Bool_t calibTT = kTRUE; // It is set automatically, but here we force to use ir or not in any case
    Bool_t badMap  = kTRUE; // It is set automatically, but here we force to use it or not in any case  
       
    AliAnalysisTaskEMCALClusterize * cl = 
    AddTaskEMCALClusterize(clustersArray, outAOD, kMC, exo,sClust,"", clTM,
                           minEcell,minEseed,dTime,wTime,unfMinE,unfFrac,
                           calibEE,badMap,calibTT,clnonlin,
                           cenBin[0],cenBin[1],-1,1,1,filterEvents,xTalkEmul,updateCell);
    
    //cl->GetRecoUtils()->SwitchOffRunDepCorrection(); // Off for Run2 for the moment 

    //  Force option different than clTM 
    //  cl->GetRecoUtils()->SwitchOnAODHybridTracksMatch();
    //  cl->GetRecoUtils()->SwitchOnAODTPCOnlyTracksMatch();
    //  cl->GetRecoUtils()->SetAODTrackFilterMask(128);
    cl->GetRecoUtils()->SetRequireTrackDCA(kFALSE); // careful with this on old MC

    //  
    if ( kMC )
    {
      cl->SwitchOnUseClusterMCLabelForCell(0) ;  // For Old Run1 MC
      cl->GetRecoUtils()->SetNonLinearityFunction(AliEMCALRecoUtils::kPi0MCv3);
    }
      
    //    cl->GetRecoUtils()->SetWarmChannelAsGood();
    //    cl->GetRecoUtils()->SetDeadChannelAsGood();
    //    cl->GetRecoUtils()->SetHotChannelAsGood();
    
    clustersArray = Form("%s_Ecell%d_Eseed%d",sClust.Data(),minEcell,minEseed);
    cl->SetAODBranchName(clustersArray);

    if ( !updateCell )
    {
      cellsArray = "Cells_Updated";
      if ( xTalkEmul > 0 )
        cellsArray = "Cells_xTalkEmulation";
    }
    cl->SetAODCellsName (cellsArray);
    
    //      cl->SetMaxEvent(20);
    //      cl->SetDebugLevel(100);
  }

  /*  
   // -----------------
   // Photon conversion
   // ----------------- 
   
   if(kInputData=="ESD"){
   printf("* Configure photon conversion analysis in macro \n");
   TString arguments = "-run-on-train -use-own-xyz  -force-aod -mc-off ";
#if defined(__CINT__)
   gROOT->LoadMacro("$ALICE_PHYSICS/PWGGA/GammaConversion/macros/ConfigGammaConversion.C");
#endif
   AliAnalysisTaskGammaConversion * taskGammaConversion = 
   ConfigGammaConversion(arguments,mgr->GetCommonInputContainer());
   taskGammaConversion->SelectCollisionCandidates();
   
   // Gamma Conversion AOD to AODPWG4Particle
   AliAnalysisTaskGCPartToPWG4Part * taskGCToPC = new AliAnalysisTaskGCPartToPWG4Part("GCPartToPWG4Part");
   taskGCToPC->SetGammaCutId("90035620401003321022000000090");
   mgr->AddTask(taskGCToPC);
   mgr->ConnectInput  (taskGCToPC, 0, mgr->GetCommonInputContainer() );
   mgr->ConnectOutput (taskGCToPC, 0, mgr->GetCommonOutputContainer()); 
   }
   */  

  // -----------------
  // CaloTrack Correlations Task
  // -----------------

  // Common settings for Correlation and QA tasks
  Bool_t   calibrate     = kFALSE;
  Int_t    minCen        = -1;
  Int_t    maxCen        = -1;
  Int_t    debug         = -1;
  
  // Possible triggered events
  TString  lTrig[]       = {"default","EMCAL_L0","EMCAL_L1","EMCAL_L2"}; 
  Int_t    nTrig         = 4;
  Int_t    trig0         = 0;
  Int_t    fixTrig       = -1;
  if ( fixTrig >= 0 ) 
  {
    trig0 = fixTrig;
    nTrig = fixTrig+1;
  }
  
  if ( kYear == 2011 )
    nTrig = 2;

  if ( kYear == 2010 )
    nTrig = 1;  
  
  // -----------------
  // Photon/Pi0/Isolation/Correlation etc
  // -----------------
  
  if ( bAnalysis )
  {
    Int_t    rejectEMCTrig = 0;
    Bool_t   nonLinOn      = kFALSE;
    Float_t  shshMax       = 0.27;
    Float_t  isoCone       = 0.4;
    Float_t  isoConeMin    = -1;
    Float_t  isoPtTh       = 1;
    Int_t    isoMethod     = AliIsolationCut::kSumPtIC;
    Int_t    isoContent    = AliIsolationCut::kNeutralAndCharged;
    Int_t    leading       = 0;
    Int_t    tm            = 2;
    Bool_t   mixOn         = kFALSE;
    TString  outputfile    = "";
    Bool_t   printSettings = kFALSE;
    
    TString  cutSelected      = "SPDPileUp";//"_MCEnScale_ITSonly";
     // Activate photon selection and invariant mass analysis
    TString  analysisSelected = "Photon_InvMass";//_MergedPi0_Isolation_Correlation_ClusterShape_PerSM_PerTCard";
    // More options:
    // "Photon_InvMass_MergedPi0_Isolation_Correlation_ClusterShape_PerSM_PerTCard_QA_Charged_Bkg"; 
    
#if defined(__CINT__)
    gROOT->LoadMacro("$ALICE_PHYSICS/PWGGA/CaloTrackCorrelations/macros/AddTaskGammaHadronCorrelationSelectAnalysis.C");
#endif
    
    for(Int_t itrig = trig0; itrig < nTrig; itrig++)
    {
      if ( itrig > 0 && kMC ) continue; // Any MC has only one kind of trigger
      
      AliAnalysisTaskCaloTrackCorrelation * emc = AddTaskGammaHadronCorrelationSelectAnalysis
      ("EMCAL",kMC,kYear,kCollision,kPeriod,rejectEMCTrig,clustersArray,cutSelected,calibrate,nonLinOn, analysisSelected,
       shshMax,isoCone,isoConeMin,isoPtTh,isoMethod ,isoContent,leading,
       tm,minCen,maxCen,mixOn,outputfile,printSettings,debug,lTrig[itrig]);
      
      emc->GetAnalysisMaker()->GetReader()->SetEMCALCellsListName(cellsArray);
      emc->GetAnalysisMaker()->GetReader()->SwitchOffRejectNoTrackEvents();
      
      // Careful, need time calibration to use time cuts defined in macro
      if ( !bEMCCluster && !bEMCCorrFra && !calibrate )
      {
        emc->GetAnalysisMaker()->GetReader()->SwitchOffUseEMCALTimeCut();
        emc->GetAnalysisMaker()->GetReader()->SetEMCALTimeCut(-1e10,1e10); // Open time cut
      }
      
      //  emc ->SelectCollisionCandidates( AliVEvent::kINT7 | AliVEvent::kCentral  | AliVEvent::kSemiCentral | AliVEvent::kMB ); // Done internally, here as example
      //  emc->GetAnalysisMaker()->GetReader()->SetNameOfMCEventHederGeneratorToAccept("Pythia");
      //  emc->GetAnalysisMaker()->GetReader()->SwitchOffShowerShapeSmearing();
      //  emc->GetAnalysisMaker()->GetReader()->SetSmearingFunction(AliCaloTrackReader::kNoSmearing);
      
      //  emc ->GetAnalysisMaker()->GetReader()->SwitchOnAliCentrality () ;
      //  emc->SetLastEvent(maxEvent);
      
      // // Example on how to modify settings of a sub-wagon if not in corresponding macro
      //  TList * anaList = emc->GetAnalysisMaker()->GetListOfAnalysisContainers();
      //  AliAnaClusterShapeCorrelStudies * shapeAna = (AliAnaClusterShapeCorrelStudies*) anaList->At(9);
      //  shapeAna->SetNCellBinLimits(3); // no analysis on predefined bins in nCell
      //  shapeAna->SetDistToBadMin(2);
      //  shapeAna->SwitchOnStudyColRowFromCellMax() ;
      //  shapeAna->Print("");
      
      if ( kYear < 2014 ) continue;
      
      TString dcalTrig = lTrig[itrig];
      dcalTrig.ReplaceAll("EM","D");
      
      AliAnalysisTaskCaloTrackCorrelation * dmc = AddTaskGammaHadronCorrelationSelectAnalysis
      ("DCAL",kMC,kYear,kCollision,kPeriod,rejectEMCTrig,clustersArray,cutSelected,calibrate,nonLinOn, analysisSelected,
       shshMax,isoCone,isoConeMin,isoPtTh,isoMethod ,isoContent,leading,
       tm,minCen,maxCen,mixOn,outputfile,printSettings,-1,lTrig[itrig]);
      dmc->GetAnalysisMaker()->GetReader()->SetEMCALCellsListName(cellsArray);
      dmc->GetAnalysisMaker()->GetReader()->SwitchOffRejectNoTrackEvents();
      
      // Careful, need time calibration to use time cuts defined in macro
      if ( !bEMCCluster && !bEMCCorrFra && !calibrate )
      {
        dmc->GetAnalysisMaker()->GetReader()->SwitchOffUseEMCALTimeCut();
        dmc->GetAnalysisMaker()->GetReader()->SetEMCALTimeCut(-1e10,1e10); // Open time cut
      }
    } // trigger loop
  } // bAnalysis
  
  // -----------------
  // QA train analysis, comment out since it cannot compile with bAnalysis uncommented
  // -----------------

//  if ( bAnalysisQA )
//  {
//    Int_t    minTime       = -1000;
//    Int_t    maxTime       =  1000;
//    Bool_t   qaan          = kTRUE;
//    Bool_t   hadronan      = kTRUE;
//    
//#if defined(__CINT__)
//    gROOT->LoadMacro("$ALICE_PHYSICS/PWGGA/CaloTrackCorrelations/macros/QA/AddTaskPi0IMGammaCorrQA.C");
//#endif
//    
//    // To test the train environment variables
//    //
//    {
//      char col [1024];
//      char tag [1024];
//      char mc  [1024];
//      
//      sprintf(col,"%s",kCollision.Data());
//      sprintf(tag,"%s",kPeriod   .Data());
//      if ( kMC ) sprintf(mc,"MC" );
//      else       sprintf(mc,"RAW");
//
//      gSystem->Setenv("ALIEN_JDL_LPMINTERACTIONTYPE",col);
//      gSystem->Setenv("ALIEN_JDL_LPMPRODUCTIONTAG"  ,tag);
//      gSystem->Setenv("ALIEN_JDL_LPMPRODUCTIONTYPE" ,mc );
//    }
//    
//    for(Int_t itrig = trig0; itrig < nTrig; itrig++)
//    {
//      if ( itrig > 0 && kMC ) continue; // Any MC has only one kind of trigger
//      
//      AliAnalysisTaskCaloTrackCorrelation * qaTrain = AddTaskPi0IMGammaCorrQA
//      ("EMCAL",kMC,"","",qaan,hadronan,calibrate,minTime,maxTime,
//        minCen,maxCen,debug,lTrig[itrig]);
//    }
//    
//    // Other QA
//    // Detector QA
//#if defined(__CINT__)
//    //  gROOT->LoadMacro("$ALICE_PHYSICS/PWGGA/CaloTrackCorrelations/macros/QA/AddTaskCalorimeterQA.C");  
//#endif
//    //  AliAnalysisTaskCaloTrackCorrelation * qatask = AddTaskCalorimeterQA(kInputData,kYear,kPrint,kMC); 
//    // 
//    // Very old Trigger QA, not in use
//#if defined(__CINT__)
//    //  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/EMCAL/macros/AddTaskEMCALTriggerQA.C");  
//#endif
//    //  AliAnalysisTaskEMCALTriggerQA * qatrigtask = AddTaskEMCALTriggerQA(); 
//  } // bAnalysis QA
//  
      // Simple event counting tasks
//      
//    #if defined(__CINT__)
//      gROOT->LoadMacro("$ALICE_PHYSICS/PWGGA/CaloTrackCorrelations/macros/AddTaskCounter.C");   
//    #endif
//    
//      AliAnalysisTaskCounter* count    = AddTaskCounter("",kMC);   // All, fill histo with cross section and trials if kMC is true
//      AliAnalysisTaskCounter* countmb  = AddTaskCounter("MB"); // Min Bias
//      AliAnalysisTaskCounter* countany = AddTaskCounter("Any"); 
//      AliAnalysisTaskCounter* countint = AddTaskCounter("AnyINT");// Min Bias
//      
//      if ( !kMC )
//      {
//        AliAnalysisTaskCounter* countemg = AddTaskCounter("EMCEGA"); 
//        AliAnalysisTaskCounter* countemj = AddTaskCounter("EMCEJE"); 
//        if ( kCollision=="PbPb" )
//        {
//          AliAnalysisTaskCounter* countcen = AddTaskCounter("Central"); 
//          AliAnalysisTaskCounter* countsce = AddTaskCounter("SemiCentral"); 
//          AliAnalysisTaskCounter* countssce= AddTaskCounter("SemiOrCentral"); 
//          AliAnalysisTaskCounter* countphP = AddTaskCounter("PHOSPb"); 
//        }
//        else
//        {
//          AliAnalysisTaskCounter* countem1 = AddTaskCounter("EMC1"); // Trig Th > 1.5 GeV approx
//          AliAnalysisTaskCounter* countem7 = AddTaskCounter("EMC7"); // Trig Th > 4-5 GeV 
//          AliAnalysisTaskCounter* countphp = AddTaskCounter("PHOS"); 
//        }
//      }  
//     
      
  
  
  //-----------------------
  // Run the analysis
  //-----------------------    
  mgr->InitAnalysis();
  mgr->PrintStatus();
  
  if      (mode == mPlugin) mgr->StartAnalysis("grid");
  else if (mode == mPROOF ) mgr->StartAnalysis("proof",chain);
  else                      mgr->StartAnalysis("local",chain);
  
  cout <<" Analysis ended sucessfully "<< endl ;
}


