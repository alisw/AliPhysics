// VARIABLES RANGES
const Double_t ymin          =  -4.0 	    ;		
const Double_t ymax          =  -2.5 	    ;
const Double_t ptmin         =   0.0 	    ;
const Double_t ptmax         =  20.   	    ;
const Double_t cCSmin        =  -1.	    ;
const Double_t cCSmax        =   1.	    ;
const Double_t cHEmin        =  -1.         ;
const Double_t cHEmax        =   1.         ;
const Double_t pCSmin        =   0.         ;           
const Double_t pCSmax        =   TMath::Pi();  
const Double_t pHEmin        =   0.         ;           
const Double_t pHEmax        =   TMath::Pi(); 
const Double_t massmin       =   0.         ; 
const Double_t massmax       =  12.         ; 		
const Double_t trigmin       =   0.	    ;
const Double_t trigmax       =   4.         ;
const Double_t ptmuminMIN    =   0.	    ;
const Double_t ptmuminMAX    =   100.	    ;
const Double_t ptmumaxMIN    =   0.         ; 
const Double_t ptmumaxMAX    =   100.       ; 
const Double_t thetamuminMIN =   0.	    ;
const Double_t thetamuminMAX =   180.	    ;
const Double_t thetamumaxMIN =   0.         ;
const Double_t thetamumaxMAX =   180.       ;
const Double_t pmuminMIN     =   0.	    ;
const Double_t pmuminMAX     =   100.	    ;
const Double_t pmumaxMIN     =   0.         ;
const Double_t pmumaxMAX     =   100.       ;
const Double_t trigsideMIN   =   0	    ;
const Double_t trigsideMAX   =   4	    ;
 


AliAnalysisTaskDimuonCFContainerBuilder *AddTaskDimuonCFContainerBuilder(Bool_t readAOD=kTRUE, Bool_t readMC=kTRUE, 
						Bool_t isaccept = kTRUE, Double_t beamEn=3500)
{

   // Check and Info printings
   //==============================================================================
    if(!readMC && isaccept) {
   	printf("ERROR: incompatible choice readMC-isaccept. If isaccept you must readMC!\n"); 
	return NULL;
    } else if (readMC && isaccept) printf("Creating task for filling a CFcontainer with acceptance data.\n");
    else if (readMC && !isaccept) printf("Creating task for filling a CFcontainer with simulated data.\n");
    else printf("Creating task for filling a CFcontainer with real data.\n");
       

   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
      ::Error("AddDataTaskMuonPolarCF", "No analysis manager to connect to.");
      return NULL;
    }   
   

   // MC handler if needed
   //==============================================================================
    if(!readAOD && readMC){
     AliMCEventHandler *mcH = (AliMCEventHandler*)mgr->GetMCtruthEventHandler();
     if (!mcH) {
       ::Error("AddDataTaskMuonPolarCF", "No MC handler connected");
       return NULL;
     }	
    }


   // DEFINING CONTAINER
   //==============================================================================

     UInt_t y  	      = 0;		// Association of variables with int numbers
     UInt_t pt        = 1;
     UInt_t costHE    = 2;
     UInt_t phiHE     = 3;
     UInt_t costCS    = 4;
     UInt_t phiCS     = 5;
     UInt_t mass      = 6;
     UInt_t trig      = 7;
     UInt_t ptmumin   = 8;
     UInt_t ptmumax   = 9;
     UInt_t thetamumin=10;
     UInt_t thetamumax=11;
     UInt_t pmumin    =12;
     UInt_t pmumax    =13;
     UInt_t trigside  =14;

     UInt_t nstep = 8;			// Number of layers (always 8 - not always filled). 4 with CINT and 4 with CMU

     const Int_t nvar   = 15 ;		// Number of variables of the grid

     const Int_t nbin1  = 10 ;		// Number of bins for each variable
     const Int_t nbin2  = 10 ;  
     const Int_t nbin3  = 20 ;  
     const Int_t nbin4  = 20 ;  
     const Int_t nbin5  = 20 ;  
     const Int_t nbin6  = 20 ;  
     const Int_t nbin7  = 240;  
     const Int_t nbin8  = 40 ;  
     const Int_t nbin9  = 100;  
     const Int_t nbin10 = 100;  
     const Int_t nbin11 = 180;  
     const Int_t nbin12 = 180;  
     const Int_t nbin13 = 100;  
     const Int_t nbin14 = 100;  
     const Int_t nbin15 = 4  ;  

     Int_t iBin[nvar];			// Array containing the number of bins for each variable
     iBin[0] =nbin1;
     iBin[1] =nbin2;
     iBin[2] =nbin3;
     iBin[3] =nbin4;
     iBin[4] =nbin5;
     iBin[5] =nbin6;
     iBin[6] =nbin7;
     iBin[7] =nbin8;
     iBin[8] =nbin9;
     iBin[9] =nbin10;
     iBin[10]=nbin11;
     iBin[11]=nbin12;
     iBin[12]=nbin13;
     iBin[13]=nbin14;
     iBin[14]=nbin15;
  
     Double_t *binLim1 = new Double_t[nbin1+1];		// Arrays for lower bounds 
     Double_t *binLim2 = new Double_t[nbin2+1];
     Double_t *binLim3 = new Double_t[nbin3+1];
     Double_t *binLim4 = new Double_t[nbin4+1];
     Double_t *binLim5 = new Double_t[nbin5+1];
     Double_t *binLim6 = new Double_t[nbin6+1];
     Double_t *binLim7 = new Double_t[nbin7+1];
     Double_t *binLim8 = new Double_t[nbin8+1];
     Double_t *binLim9 = new Double_t[nbin9+1];
     Double_t *binLim10= new Double_t[nbin10+1];
     Double_t *binLim11= new Double_t[nbin11+1];
     Double_t *binLim12= new Double_t[nbin12+1];
     Double_t *binLim13= new Double_t[nbin13+1];
     Double_t *binLim14= new Double_t[nbin14+1];
     Double_t *binLim15= new Double_t[nbin15+1];
     for(Int_t i=0; i<=nbin1; i++) binLim1[i] =(Double_t)ymin+(ymax-ymin)/nbin1*(Double_t)i;
     for(Int_t i=0; i<=nbin2; i++) binLim2[i] =(Double_t)ptmin+(ptmax-ptmin)/nbin2*(Double_t)i; 
     for(Int_t i=0; i<=nbin3; i++) binLim3[i] =(Double_t)cHEmin+(cHEmax-cHEmin)/nbin3*(Double_t)i ; 
     for(Int_t i=0; i<=nbin4; i++) binLim4[i] =(Double_t)pHEmin+(pHEmax-pHEmin)/nbin4*(Double_t)i ; 
     for(Int_t i=0; i<=nbin5; i++) binLim5[i] =(Double_t)cCSmin+(cCSmax-cCSmin)/nbin5*(Double_t)i ; 
     for(Int_t i=0; i<=nbin6; i++) binLim6[i] =(Double_t)pCSmin+(pCSmax-pCSmin)/nbin6*(Double_t)i ; 
     for(Int_t i=0; i<=nbin7; i++) binLim7[i] =(Double_t)massmin+(massmax-massmin)/nbin7*(Double_t)i ; 
     for(Int_t i=0; i<=nbin8; i++) binLim8[i] =(Double_t)trigmin+(trigmax-trigmin)/nbin8*(Double_t)i ; 
     for(Int_t i=0; i<=nbin9; i++) binLim9[i] =(Double_t)ptmuminMIN+(ptmuminMAX-ptmuminMIN)/nbin9*(Double_t)i ; 
     for(Int_t i=0; i<=nbin10;i++) binLim10[i]=(Double_t)ptmumaxMIN+(ptmumaxMAX-ptmumaxMIN)/nbin10*(Double_t)i ; 
     for(Int_t i=0; i<=nbin11;i++) binLim11[i]=(Double_t)thetamuminMIN+(thetamuminMAX-thetamuminMIN)/nbin11*(Double_t)i ; 
     for(Int_t i=0; i<=nbin12;i++) binLim12[i]=(Double_t)thetamumaxMIN+(thetamumaxMAX-thetamumaxMIN)/nbin12*(Double_t)i ; 
     for(Int_t i=0; i<=nbin13;i++) binLim13[i]=(Double_t)pmuminMIN+(pmuminMAX-pmuminMIN)/nbin13*(Double_t)i ; 
     for(Int_t i=0; i<=nbin14;i++) binLim14[i]=(Double_t)pmumaxMIN+(pmumaxMAX-pmumaxMIN)/nbin14*(Double_t)i ; 
     for(Int_t i=0; i<=nbin15;i++) binLim15[i]=(Double_t)trigsideMIN+(trigsideMAX-trigsideMIN)/nbin15*(Double_t)i ; 

     AliCFContainer* container = new AliCFContainer("container","Container for Dimuons",nstep,nvar,iBin);
  
     container -> SetBinLimits(y,binLim1);		// setting the bin limits
     container -> SetBinLimits(pt,binLim2);
     container -> SetBinLimits(costHE,binLim3);
     container -> SetBinLimits(phiHE,binLim4);
     container -> SetBinLimits(costCS,binLim5);
     container -> SetBinLimits(phiCS,binLim6);
     container -> SetBinLimits(mass,binLim7);
     container -> SetBinLimits(trig,binLim8);
     container -> SetBinLimits(ptmumin,binLim9);
     container -> SetBinLimits(ptmumax,binLim10);
     container -> SetBinLimits(thetamumin,binLim11);
     container -> SetBinLimits(thetamumax,binLim12);
     container -> SetBinLimits(pmumin,binLim13);
     container -> SetBinLimits(pmumax,binLim14);
     container -> SetBinLimits(trigside,binLim15);


   // CF Manager
   //==============================================================================
     AliCFManager* man = new AliCFManager() ;
     man->SetParticleContainer(container);

  
   // Outputs: list of histograms + CFContainer
   //==============================================================================
     TString outputfile = AliAnalysisManager::GetCommonFileName();
     outputfile += ":PWG3Muon_DimuonCFContainer";

     AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("Histos",TList::Class(),AliAnalysisManager::kOutputContainer,outputfile);
     AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("DimuonCFContainer",AliCFContainer::Class(),AliAnalysisManager::kOutputContainer,outputfile);

  // The task with the associtated CF manager
   //==============================================================================
     AliAnalysisTaskDimuonCFContainerBuilder *task = new AliAnalysisTaskDimuonCFContainerBuilder("AliAnalysisTaskDimuonCFContainerBuilder",readAOD,readMC,isaccept,beamEn);
     task->SetCFManager(man);


   // Additional settings for the task (including some cuts)
   //==============================================================================
    //task->SetDistinguishTrigClass(kTRUE);
    //task->SetReadMCinfo(kTRUE);
    //Double_t ptlimits[2]={1.,1000.};
    //task->SetPtSingMuLimits(ptlimits);
    //task->SetCutonZvtxSPD(kTRUE);
    //Double_t vtxlimits[2]={-10.,10.};
    //task->SetZprimVertLimits(vtxlimits);


   // Adding the task to the analysis manager and connecting inputs and outputs
   //==============================================================================
    mgr->AddTask(task);
    mgr->ConnectInput(task,  0, mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task,1,coutput1);
    mgr->ConnectOutput(task,2,coutput2);

   return task;
}
