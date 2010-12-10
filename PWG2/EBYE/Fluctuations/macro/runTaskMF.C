/****************************************************************

Lunch the Analysis Task for Multiplicity Fluctuation 
Author: Satyajit Jena
Email:  sjena@cern.ch

Running in Local = analmode = local, 
input list of files, number of esds
useMC and format

and AFF = analmode = "proof" and
data set.

*******************************************************************/

void runTaskMF(TString analmode = "local")
{

 TString dataset = "/PMD/sjena/LHC10a10_2360Gev_pythia";
 Bool_t  useMC = kFALSE;
 TString format="esd";
 TString output = "test.root";


 if(analmode.CompareTo("local") == 0)
    {
      if(!LoadLibraries(analmode))
	{
	  printf("Library Not loaded\n");
	  return;
	}
      runLocal("file.txt",5, useMC, format, output);
      
    }
 else  if(analmode.CompareTo("proof") == 0)
   {
     if(!LoadLibraries(analmode)) return;
     runproof(dataset, useMC);
          
     
   }
 else printf("load error\n");

}

//___________________________________________________________________________________
Bool_t LoadLibraries(TString mode = "local")
{

  if(mode.CompareTo("local") == 0)
    {
      gSystem->Load("libSTEERBase.so");
      gSystem->Load("libVMC.so");
      gSystem->Load("libESD.so");
      gSystem->Load("libAOD.so");
      gSystem->Load("libANALYSIS.so");
      gSystem->Load("libANALYSISalice.so");
      gSystem->Load("libPWG2ebye.so");
      gSystem->AddIncludePath("-I$ALICE_ROOT/include");
      printf("Library is Loaded \n");

      return kTRUE;

    }
  else if(mode.CompareTo("proof") == 0)
    {
      TString proofCluster="alice-caf.cern.ch";
      TString alirootVer = "VO_ALICE@AliRoot::v4-20-11-AN";
      
      gEnv->SetValue("XSec.GSI.DelegProxy","2");
      TProof *p = TProof::Open(proofCluster.Data());
      
      p->ShowPackages();
      p->EnablePackage(alirootVer.Data());
     
      return kTRUE;
      
    }
  else 
    {
      printf(" ERROR: Give proper running mode \n");
      return kFALSE;
    }
  
}
//______________________________________________________________________
// Running local
void runLocal(TString input,int nfile,  Bool_t useMC,  TString format,TString OutPutFile)
{

  gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
  TChain* chain = CreateESDChain(input, nfile);
  
  AliAnalysisManager *mgr = new AliAnalysisManager("Analysis");
  AliVEventHandler* esdH = new AliESDInputHandler();
  mgr->SetInputEventHandler(esdH);

//  gROOT->LoadMacro("$ALICE_ROOT/PWG2/EBYE/Fluctuations/AliEbyEEventSelector.cxx++g");  
//  gROOT->LoadMacro("$ALICE_ROOT/PWG2/EBYE/Fluctuations/AliEbyEMultiplicityFluctuationAnalysis.cxx++g");  
//  gROOT->LoadMacro("$ALICE_ROOT/PWG2/EBYE/Fluctuations/AliEbyEMFAnalysisTask.cxx++g");  
  gROOT->LoadMacro("$ALICE_ROOT/PWG2/EBYE/Fluctuations/macro/AddTaskMF.C");
  AddTaskMF();

  mgr->SetDebugLevel(0);
  
  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();

  mgr->InitAnalysis();
  mgr->StartAnalysis("local", chain);


}


//________________________________________________________
//Running in proof

void runproof(TString dataset, Bool_t useMC)
{

     AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
     if (!mgr) mgr = new AliAnalysisManager("FirstCheck");
     
     AliVEventHandler* esdH = new AliESDInputHandler();
     mgr->SetInputEventHandler(esdH);

     gProof->LoadMacro("$ALICE_ROOT/PWG2/EBYE/Fluctuations/AliEbyEEventSelector.cxx++g");  
     gProof->LoadMacro("$ALICE_ROOT/PWG2/EBYE/Fluctuations/AliEbyEMFAnalysisTask.cxx++g");  
     gProof->LoadMacro("$ALICE_ROOT/PWG2/EBYE/Fluctuations/macro/AddTaskMF.C");
     AddTaskMF();
  
     mgr->SetDebugLevel(0);
     
     if (!mgr->InitAnalysis()) 
    return;
     
     mgr->PrintStatus();
     mgr->InitAnalysis();
     
     mgr->StartAnalysis("proof",dataset.Data());
  
}

//_________________
//Runningin grid
