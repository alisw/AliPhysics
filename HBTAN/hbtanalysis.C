#ifndef __MAKECINT__ 
 #ifndef __CINT__
  #include "alles.h"
  #include <TPDGCode.h>
  #include "AliRunAnalysis.h"  
  #include "AliHBTAnalysis.h"  
  #include "AliReader.h"
  #include "AliReaderESD.h"
  #include "AliReaderKineTree.h"
  #include "AliParticleCut.h"
  #include "AliAOD.h"
  #include "AliAODPairCut.h"
  #include "AliAODParticleCut.h"
  #include "AliHBTQResolutionFctns.h"
  #include "AliHBTQDistributionFctns.h"
  #include "AliHBTMonDistributionFctns.h"
  #include "AliHBTTwoTrackEffFctn.h"
  #include "AliHBTCorrelFctn.h"
  #include "AliHBTWeightFctn.h"
  #include "AliHBTWeightTheorFctn.h"
  #include "AliHBTWeights.h"
  #include "AliHBTLLWeights.h"
  #include "TSystem.h"
  #include "TObjString.h"
  #include "TString.h"
  #include "AliPDG.h"
  #include "AliHBTMonDistributionFctns.h"
  #include "AliHBTMonResolutionFctns.h"
 #endif 
#endif

//This macro can not be executed from 
AliHBTAnalysis* buildanalysis(Int_t pid1, Int_t pid2, Bool_t threeDcuts, Float_t qinvmax, const char* anasuffix = "");


void hbtanalysis(Option_t* datatype="Kine", Option_t* processopt="TracksAndParticles",
                Int_t first = -1,Int_t last = -1, 
                char *prefname = "")
 {
  const char* basedir=".";
  const char* serie=".";
  const char* field = ".";
 
  Int_t PID[10];
  PID[0]=kPiPlus;
  PID[1]=kPiMinus;
  PID[1]=0;

  Bool_t ident = kTRUE;//identical particles 
  Bool_t nonident = kFALSE;//non identical particles analysis
// Dynamically link some shared libs                    
//  gROOT->LoadMacro("loadlibs.C");
//  loadlibs();
  gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libANALYSIS");
  gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libHBTAN");
  
  /***********************************************************/
  //Create analysis and reader
  
  TString dataname;
  AliReader* reader;
  Int_t kine = strcmp(datatype,"Kine");
  Int_t ESD = strcmp(datatype,"ESD");
  Int_t intern = strcmp(datatype,"Intern");
  
  if(!kine)
   {
    reader = new AliReaderKineTree();
    processopt="Particles"; //this reader by definition reads only simulated particles
    dataname = "Kine";
   }
  else if(!ESD)
   {
    reader = new AliReaderESD();
    dataname = "TPC";
   }
  else if(!intern)
   {
//    dataname = "TPC";
//    reader = new AliReaderAOD(dataname+".data.root");
   }
  else
   {
    cerr<<"Option "<<datatype<<"  not recognized. Exiting"<<endl;
    return;
   }

  cout<<"hbtanalysis.C: Reader Set, "<<reader->Class()->GetName()<<" dataname: "<<dataname<<endl;
  cout<<"hbtanalysis.C: Setting dirs\n";
  TObjArray* dirs=0;
  if ( (first >= 0) && (last>=0) && ( (last-first)>=0 ) )
   {//read from many dirs dirs
     char buff[50];
     dirs = new TObjArray(last-first+1);
     for (Int_t i = first; i<=last; i++)
      {
//        sprintf(buff,"%s/%s/%s/%05.5d",basedir,field,serie,i);
        sprintf(buff,"%s/%s/%s/%d",basedir,field,serie,i);
        TObjString *odir= new TObjString(buff);
        dirs->Add(odir);
      }
    }
   reader->SetDirs(dirs);
   reader->SetName("READER");
   reader->SetTitle("Reader Used in HBT Analysis");
   cout<<"hbtanalysis.C: Dirs set\n";
  /***********************************************************/    
  /*****   R E A D E R ***************************************/    
  /***********************************************************/    
  
  //we want have only low pt pi+ so set a cut to reader
  AliAODParticleCut* readerpartcut= new AliAODParticleCut();
  readerpartcut->SetPtRange(0.0,10000.0);
//  readerpartcut->SetThetaRange(0.25*TMath::Pi(), 0.75*TMath::Pi());//45->135 deg.
  readerpartcut->SetYRange(-1.0,1.0);
  
  Int_t l = 0;
  while (PID[l] != 0)
   {
     readerpartcut->SetPID(PID[l]);
     reader->AddParticleCut(readerpartcut);//read this particle type with this cut
     l++;
   }
  
  TString filename = prefname + dataname + ".anal.";
  TString fn;
  TParticlePDG* pdg1;  
  TParticlePDG* pdg2;  
  TString p2name;  
  TString p1name ;

  /************************************************************/
  /****   P R O C E S S                 ***********************/
  /************************************************************/
  AliHBTAnalysis* hbtanalysis = 0x0;
  
  l = 0;
  Int_t k;
  while (PID[l] != 0)
   {
     pdg1 = TDatabasePDG::Instance()->GetParticle(PID[l]);
     if (pdg1 == 0x0)
       {
         cerr<<"ERROR: hbtanalysis.C: Can not find particle type "<<PID[l]<<" in PDG DataBase\n";
         l++;
         continue;
       }
     p1name= pdg1->GetName();
     k=l;
     while (PID[k] != 0)
      {
        if ( (k==l)&&(ident==kFALSE) )  { k++;continue;}
        if ( (k!=l)&&(nonident==kFALSE) ) { k++;continue; }
        
        pdg2 = TDatabasePDG::Instance()->GetParticle(PID[k]);
        if (pdg1 == 0x0)
         {
           cerr<<"ERROR: hbtanalysis.C: Can not find particle type "<<PID[k]<<" in PDG DataBase\n";
           k++;
           continue;
         }
        p2name = pdg2->GetName();
        cout<<"hbtanalysis.C: buildanalysis for "<<PID[l]<<"("<<p1name<<")  "
                                  <<PID[k]<<"("<<p2name<<")\n";
        fn = filename + p1name + p2name + ".root";
        cout<<"hbtanalysis.C: File Name is: "<<fn<<endl;
        
        AliRunAnalysis* analysis = new AliRunAnalysis();
        analysis->SetReader(reader);
        
        cout<<"hbtanalysis.C: buildanalysis "<<endl;
        hbtanalysis = buildanalysis(PID[l],PID[k],kFALSE,1.0,"-a1");
        cout<<"hbtanalysis.C: buildanalysis Done"<<endl;
        if (hbtanalysis == 0x0)
         {
           cout<<"ERROR: hbtanalysis.C: buildanalysis function returned NULL pointer\n";
           cout<<"ERROR: Deleting Run Analysis ..."<<endl;
           delete analysis;
           cout<<"ERROR: Deleting Run Analysis Done"<<endl;
           k++;
           continue; 
         }
        TFile* histoOutput = TFile::Open(fn,"recreate");
//        hbtanalysis->SetOutputFileName(fn);
        
        AliHBTLLWeights::Instance()->SetParticlesTypes(PID[l],PID[k]);
        AliHBTLLWeights::Instance()->SetColoumb(kFALSE);
        AliHBTLLWeights::Instance()->SetStrongInterSwitch(kFALSE);
        AliHBTLLWeights::Instance()->Init();
        AliHBTLLWeights::Instance()->Set();
        
//        AliHBTCrab::Instance()->Init(PID[l],PID[k]);
//        AliHBTCrab::Instance()->Set();
        
        hbtanalysis->SetProcessOption(AliHBTAnalysis::kSimulated);
        analysis->Add(hbtanalysis);

        histoOutput->cd();
        cout<<"hbtanalysis.C: Run\n";
        analysis->Run();
        
        cout<<"hbtanalysis.C: Write\n";
        reader->Write();
        cout<<"hbtanalysis.C: Delete analysis\n";
        hbtanalysis->SetOwner(kTRUE);
        delete analysis;
        histoOutput->Close();
        k++;
      }
     l++;
   }
  
  delete readerpartcut;
  if (dirs) 
   {
    dirs->SetOwner();
    delete dirs;
   }

 }

AliHBTAnalysis* buildanalysis(Int_t pid1, Int_t pid2, Bool_t threeDcuts, Float_t qinvmax, const char* anasuffix)
 {
  /************************************************************/
  /****   Q INV Correlation Function   ************************/
  /************************************************************/
  
  
  AliHBTAnalysis * analysis = new AliHBTAnalysis();
  analysis->SetBufferSize(1);
  
  
  AliAODPairCut *paircut = new AliAODPairCut();
  
  TString anasfx(anasuffix);
  
  Float_t qinvmin = 0.0;
  paircut->SetQInvRange(qinvmin,1.5*qinvmax);  

  TParticlePDG* pdg1 = TDatabasePDG::Instance()->GetParticle(pid1);
  TString system = pdg1->GetName();
  system+=" ";
  if (pid2)
   {
    cout<<"Setting cuts on particles: "<<pid1<<"  "<<pid2<<endl;
    AliAODParticleCut* partcut = new AliAODParticleCut();
    partcut->SetPID(pid1);
    paircut->SetFirstPartCut(partcut);
    partcut->SetPID(pid2);
    paircut->SetSecondPartCut(partcut);
    delete partcut;
    TParticlePDG* pdg2 = TDatabasePDG::Instance()->GetParticle(pid2);
    system+=pdg2->GetName();
   }
  else
   {
    system+=pdg1->GetName();
   }
  analysis->SetGlobalPairCut(paircut);
  paircut->SetName("AnalysisGlobalPairCut");
  paircut->Write();

  cout<<"Building functions"<<endl;
  
  /************************************************************/
  /****   Q INV Correlation Function   ************************/
  /************************************************************/
  //general cross check
  AliHBTQInvCorrelFctn * qinvcfP = new AliHBTQInvCorrelFctn(1000,qinvmax);
  qinvcfP->SetNumberOfBinsToScale(200);
  qinvcfP->Rename("qinvcfP"+anasfx,"Particle (simulated) Q_{inv} Correlation Function "+ system + anasfx);
  analysis->AddParticleFunction(qinvcfP);

//===========================================================================================  

  AliHBTWeightTheorQInvFctn  *wqinvcf = new AliHBTWeightTheorQInvFctn(1000,qinvmax);
  wqinvcf->SetNumberOfBinsToScale(200);
  wqinvcf->Rename("wqinvcf","Particle (simulated) Lednicky Q_{inv} Correlation Function "+system);
  analysis->AddParticleFunction(wqinvcf);
   
  AliHBTWeightTheorOSLFctn *wqoslcf = new AliHBTWeightTheorOSLFctn(60,qinvmax,0.0,60,qinvmax,0.0,60,qinvmax,0.0);
  wqoslcf->Rename("wqoslcf","Particle (simulated) Lednicky OSL Theoret Correlation Function "+system);
  wqoslcf->SetNumberOfBinsToScale(20,20,20);
  analysis->AddParticleFunction(wqoslcf);

//===========================================================================================  
  AliAODPairCut *cfPairCutKt = new AliAODPairCut();

//===========================================================================================  
  cfPairCutKt->SetKtRange(0.0,0.2);
  
  AliHBTWeightTheorQInvFctn  *wqinvcfK11 = new AliHBTWeightTheorQInvFctn(1000,qinvmax);
  wqinvcfK11->SetNumberOfBinsToScale(200);
  wqinvcfK11->SetPairCut(cfPairCutKt);
  wqinvcfK11->Rename("wqinvcfK11","Particle (simulated) Lednicky Q_{inv} Correlation Function 0.0<p_{t}<0.2 "+system);
  analysis->AddParticleFunction(wqinvcfK11);
   
  AliHBTWeightTheorOSLFctn *wqoslcfK11 = new AliHBTWeightTheorOSLFctn(60,qinvmax,0.0,60,qinvmax,0.0,60,qinvmax,0.0);
  wqoslcfK11->Rename("wqoslcfK11","Particle (simulated) Lednicky OSL Theoret Correlation Function 0.0<p_{t}<0.2 "+system);
  wqoslcfK11->SetNumberOfBinsToScale(20,20,20);
  wqoslcfK11->SetPairCut(cfPairCutKt);
  analysis->AddParticleFunction(wqoslcfK11);
//===========================================================================================  

  cfPairCutKt->SetKtRange(0.2,0.4);
  
  AliHBTWeightTheorQInvFctn  *wqinvcfK12 = new AliHBTWeightTheorQInvFctn(1000,qinvmax);
  wqinvcfK12->SetNumberOfBinsToScale(200);
  wqinvcfK12->SetPairCut(cfPairCutKt);
  wqinvcfK12->Rename("wqinvcfK12","Particle (simulated) Lednicky Q_{inv} Correlation Function 0.2<p_{t}<0.4 "+system);
  analysis->AddParticleFunction(wqinvcfK12);
   
  AliHBTWeightTheorOSLFctn *wqoslcfK12 = new AliHBTWeightTheorOSLFctn(60,qinvmax,0.0,60,qinvmax,0.0,60,qinvmax,0.0);
  wqoslcfK12->Rename("wqoslcfK12","Particle (simulated) Lednicky OSL Theoret Correlation Function 0.2<p_{t}<0.4 "+system);
  wqoslcfK12->SetNumberOfBinsToScale(20,20,20);
  wqoslcfK12->SetPairCut(cfPairCutKt);
  analysis->AddParticleFunction(wqoslcfK12);
//=========================================================================================== 

  cfPairCutKt->SetKtRange(0.4,0.6);

  AliHBTWeightTheorQInvFctn  *wqinvcfK21 = new AliHBTWeightTheorQInvFctn(1000,qinvmax);
  wqinvcfK21->SetNumberOfBinsToScale(200);
  wqinvcfK21->SetPairCut(cfPairCutKt);
  wqinvcfK21->Rename("wqinvcfK21","Particle (simulated) Lednicky Q_{inv} Correlation Function 0.4<p_{t}<0.6 "+system);
  analysis->AddParticleFunction(wqinvcfK21);

  AliHBTWeightTheorOSLFctn *wqoslcfK21 = new AliHBTWeightTheorOSLFctn(60,qinvmax,0.0, 60,qinvmax,0.0, 60,qinvmax,0.0);
  wqoslcfK21->Rename("wqoslcfK21","Particle (simulated) Lednicky OSL Theoret Correlation Function 0.4<p_{t}<0.6 "+system);
  wqoslcfK21->SetNumberOfBinsToScale(20,20,20);
  wqoslcfK21->SetPairCut(cfPairCutKt);
  analysis->AddParticleFunction(wqoslcfK21);
//===========================================================================================  

  cfPairCutKt->SetKtRange(0.6,0.8);

  AliHBTWeightTheorQInvFctn  *wqinvcfK22 = new AliHBTWeightTheorQInvFctn(1000,qinvmax);
  wqinvcfK22->SetNumberOfBinsToScale(200);
  wqinvcfK22->SetPairCut(cfPairCutKt);
  wqinvcfK22->Rename("wqinvcfK22","Particle (simulated) Lednicky Q_{inv} Correlation Function 0.6<p_{t}<0.8 "+system);
  analysis->AddParticleFunction(wqinvcfK22);

  AliHBTWeightTheorOSLFctn *wqoslcfK22 = new AliHBTWeightTheorOSLFctn(60,qinvmax,0.0, 60,qinvmax,0.0, 60,qinvmax,0.0);
  wqoslcfK22->Rename("wqoslcfK22","Particle (simulated) Lednicky OSL Theoret Correlation Function 0.6<p_{t}<0.8 "+system);
  wqoslcfK22->SetNumberOfBinsToScale(20,20,20);
  wqoslcfK22->SetPairCut(cfPairCutKt);
  analysis->AddParticleFunction(wqoslcfK22);
//===========================================================================================  

  cfPairCutKt->SetKtRange(0.8,1.2);

  AliHBTWeightTheorQInvFctn  *wqinvcfK3 = new AliHBTWeightTheorQInvFctn(1000,qinvmax);
  wqinvcfK3->SetNumberOfBinsToScale(200);
  wqinvcfK3->SetPairCut(cfPairCutKt);
  wqinvcfK3->Rename("wqinvcfK3","Particle (simulated) Lednicky Q_{inv} Correlation Function 0.8<p_{t}<1.2 "+system);
  analysis->AddParticleFunction(wqinvcfK3);

  AliHBTWeightTheorOSLFctn *wqoslcfK3 = new AliHBTWeightTheorOSLFctn(60,qinvmax,0.0, 60,qinvmax,0.0, 60,qinvmax,0.0);
  wqoslcfK3->Rename("wqoslcfK3","Particle (simulated) Lednicky OSL Theoret Correlation Function 0.8<p_{t}<1.2 "+system);
  wqoslcfK3->SetNumberOfBinsToScale(20,20,20);
  wqoslcfK3->SetPairCut(cfPairCutKt);
  analysis->AddParticleFunction(wqoslcfK3);
//===========================================================================================  

  cfPairCutKt->SetKtRange(1.2,100000.);

  AliHBTWeightTheorQInvFctn  *wqinvcfK4 = new AliHBTWeightTheorQInvFctn(1000,qinvmax);
  wqinvcfK4->SetNumberOfBinsToScale(200);
  wqinvcfK4->SetPairCut(cfPairCutKt);
  wqinvcfK4->Rename("wqinvcfK4","Particle (simulated) Lednicky Q_{inv} Correlation Function 1.2<p_{t} "+system);
  analysis->AddParticleFunction(wqinvcfK4);

  AliHBTWeightTheorOSLFctn *wqoslcfK4 = new AliHBTWeightTheorOSLFctn(60,qinvmax,0.0, 60,qinvmax,0.0, 60,qinvmax,0.0);
  wqoslcfK4->Rename("wqoslcfK4","Particle (simulated) Lednicky OSL Theoret Correlation Function 1.2<p_{t} "+system);
  wqoslcfK4->SetNumberOfBinsToScale(20,20,20);
  wqoslcfK4->SetPairCut(cfPairCutKt);
  analysis->AddParticleFunction(wqoslcfK4);
//===========================================================================================  
  

  AliHBTMonVyDistributionVsVxFctn* vydistrvsvx = new  AliHBTMonVyDistributionVsVxFctn(800,1.0e-11,-1.0e-11,800,1.0e-11,-1.0e-11);
  AliHBTMonRtDistributionVsVzFctn* rtdistrvsvz = new  AliHBTMonRtDistributionVsVzFctn(800,1.0e-11,-1.0e-11,800,1.0e-11,0.0);
  analysis->AddParticleMonitorFunction(vydistrvsvx);
  analysis->AddParticleMonitorFunction(rtdistrvsvz);

  AliHBTMonPhiDistributionVsPtFctn*  phidistr  = new AliHBTMonPhiDistributionVsPtFctn(800,0,4,200,6.3,-3.2);
  AliHBTMonThetaDistributionVsPtFctn* thetadistr = new AliHBTMonThetaDistributionVsPtFctn();
  analysis->AddParticleMonitorFunction(phidistr);
  analysis->AddParticleMonitorFunction(thetadistr);

  cout<<"Building functions Done "<<endl;
  delete paircut;
  
  return analysis;

//===========================================================================================  
//  EEEE   N   N   DDD   =====================================================================  
//  E      NN  N   D  D  =====================================================================  
//  EEE    N N N   D  D  =====================================================================  
//  E      N  NN   D  D  =====================================================================  
//  EEEE   N   N   DDD   =====================================================================  
//===========================================================================================  
  /************************************************************/
  /****   Q OUT Correlation Function   ************************/
  /************************************************************/
  AliAODPairCut *outPairCut = new AliAODPairCut();
  outPairCut->SetQSideCMSLRange(-0.01,0.01);
  outPairCut->SetQLongCMSLRange(-0.01,0.01);
  outPairCut->SetKtRange(0.0,0.8);

  AliHBTWeightTheorQOutFctn *wqcfoutP = new AliHBTWeightTheorQOutFctn(1000,qinvmax,-qinvmax);
  wqcfoutP->SetNumberOfBinsToScale(200);
  wqcfoutP->SetPairCut(outPairCut);
  analysis->AddParticleFunction(wqcfoutP);

  /************************************************************/
  /****   Q SIDE Correlation Function   ************************/
  /************************************************************/
  AliAODPairCut *sidePairCut = new AliAODPairCut();
  sidePairCut->SetQOutCMSLRange(-0.01,0.01);
  sidePairCut->SetQLongCMSLRange(-0.01,0.01);
  sidePairCut->SetKtRange(0.0,0.8);
  
  AliHBTWeightTheorQSideFctn *wqcfsideP = new AliHBTWeightTheorQSideFctn(1000,qinvmax,-qinvmax);
  wqcfsideP->SetNumberOfBinsToScale(200);
  wqcfsideP->SetPairCut(sidePairCut);
  analysis->AddParticleFunction(wqcfsideP);

  /************************************************************/
  /****   Q LONG Correlation Function   ************************/
  /************************************************************/
  AliAODPairCut *longPairCut = new AliAODPairCut();
  longPairCut->SetQOutCMSLRange(-0.01,0.01);
  longPairCut->SetQSideCMSLRange(-0.01,0.01);
  longPairCut->SetKtRange(0.0,0.8);
  
  AliHBTWeightTheorQLongFctn *wqcflongP = new AliHBTWeightTheorQLongFctn(1000,qinvmax,-qinvmax);
  wqcflongP->SetNumberOfBinsToScale(200);
  wqcfsideP->SetPairCut(longPairCut);
  analysis->AddParticleFunction(wqcflongP);

  
 
  AliHBTMonPxDistributionFctn* pxdistr = new AliHBTMonPxDistributionFctn();
  AliHBTMonPyDistributionFctn* pydistr = new AliHBTMonPyDistributionFctn();
  AliHBTMonPzDistributionFctn* pzdistr = new AliHBTMonPzDistributionFctn();
  AliHBTMonPDistributionFctn*  pdistr  = new AliHBTMonPDistributionFctn();
  analysis->AddParticleMonitorFunction(pxdistr);
  analysis->AddParticleMonitorFunction(pydistr);
  analysis->AddParticleMonitorFunction(pzdistr);
  analysis->AddParticleMonitorFunction(pdistr);
  
  AliHBTMonVxDistributionFctn* vxdistr = new AliHBTMonVxDistributionFctn(800,1.0e-11,-1.0e-11);
  AliHBTMonVyDistributionFctn* vydistr = new AliHBTMonVyDistributionFctn(800,1.0e-11,-1.0e-11);
  AliHBTMonVzDistributionFctn* vzdistr = new AliHBTMonVzDistributionFctn(800,1.0e-11,-1.0e-11);
  AliHBTMonRDistributionFctn*  vrdistr = new AliHBTMonRDistributionFctn (800,1.0e-11,-1.0e-11);
  analysis->AddParticleMonitorFunction(vxdistr);
  analysis->AddParticleMonitorFunction(vydistr);
  analysis->AddParticleMonitorFunction(vzdistr);
  analysis->AddParticleMonitorFunction(vrdistr);
  
  
  AliAODParticleCut* vxycut = new AliAODParticleCut();
  AliHBTMonVyDistributionVsVxFctn* vydistrvsvx1 = new  AliHBTMonVyDistributionVsVxFctn(800,1.0e-11,-1.0e-11,800,1.0e-11,-1.0e-11);
  vydistrvsvx1->Rename("vydistrvsvx1","Vx:Vy 0.2 < Pt < 1.8");
  vxycut->SetPtRange(0.2,0.8);
  vydistrvsvx1->SetParticleCut(vxycut);
  analysis->AddParticleMonitorFunction(vydistrvsvx1);

  AliHBTMonVyDistributionVsVxFctn* vydistrvsvx2 = new  AliHBTMonVyDistributionVsVxFctn(800,1.0e-11,-1.0e-11,800,1.0e-11,-1.0e-11);
  vydistrvsvx2->Rename("vydistrvsvx2","Vx:Vy 0.8 < Pt < 1.2");
  vxycut->SetPtRange(0.8,1.2);
  vydistrvsvx2->SetParticleCut(vxycut);
  analysis->AddParticleMonitorFunction(vydistrvsvx2);

  AliHBTMonVyDistributionVsVxFctn* vydistrvsvx3 = new  AliHBTMonVyDistributionVsVxFctn(800,1.0e-11,-1.0e-11,800,1.0e-11,-1.0e-11);
  vydistrvsvx3->Rename("vydistrvsvx3","Vx:Vy 1.2 < Pt ");
  vxycut->SetPtRange(1.2,10000.0);
  vydistrvsvx3->SetParticleCut(vxycut);
  analysis->AddParticleMonitorFunction(vydistrvsvx3);
  
//  AliHBTTwoKStarCorrelFctn* kstarP = new AliHBTTwoKStarCorrelFctn(600,qinvmax);
//  kstarP->Rename("kstarP","Particle (simulated) 2K^{*} Correlation Function "+system);
//  analysis->AddParticleFunction(kstarP);
  
  
  delete paircut;
//  delete outPairCut;
//  delete sidePairCut;
//  delete longPairCut;
  return analysis;
 }
