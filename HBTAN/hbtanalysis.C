  #include "AliHBTAnalysis.h"
//  #include "alles.h"
  #include "AliHBTReader.h"
  #include "AliHBTReaderKineTree.h"
  #include "AliHBTReaderITSv2.h"
  #include "AliHBTReaderITSv1.h"
  #include "AliHBTReaderTPC.h"
  #include "AliHBTReaderInternal.h"
  #include "AliHBTParticleCut.h"
  #include "AliHBTEvent.h"
  #include "AliHBTPairCut.h"
  #include "AliHBTQResolutionFctns.h"
  #include "AliHBTQDistributionFctns.h"
  #include "AliHBTTwoTrackEffFctn.h"
  #include "AliHBTCorrelFctn.h"
  #include "TSystem.h"
  #include "TObjString.h"
  #include "TString.h"
  #include "TPDGCode.h"
  #include "AliHBTMonDistributionFctns.h"
  #include "AliHBTMonResolutionFctns.h"
  
void hbtanalysis(Option_t* datatype, Option_t* processopt="TracksAndParticles", 
                Int_t first = -1,Int_t last = -1, 
                char *outfile = "hbtanalysis.root")
 {
   
  AliHBTTrackPoints::SetDebug(0);
  AliHBTParticle::SetDebug(0);
  AliLoader::SetDebug(0);
  
  AliHBTParticleCut c1;
  AliHBTParticleCut c2 = c1;

//HBT Anlysis Macro
//Anlyzes TPC recontructed tracks and simulated particles that corresponds to them

//datatype defines type of data to be read
//  Kine  - analyzes Kine Tree: simulated particles
//  TPC   - analyzes TPC   tracking + particles corresponding to these tracks
//  ITSv1 - analyzes ITSv1 ----------------------//--------------------------
//  ITSv2 - analyzes ITSv2 ----------------------//--------------------------

//processopt defines option passed to AliHBTAnlysis::Process method
// default: TracksAndParticles - process both recontructed tracks and sim. particles corresponding to them 
//          Tracks - process only recontructed tracks
//          Particles - process only simulated particles

//Reads data from diroctories from first to last(including)
// For examples if first=3 and last=5 it reads from
//  ./3/
//  ./4/
//  ./5/
//if first or last is negative (or both), it reads from current directory
//
//these names I use when analysis is done directly from CASTOR files via RFIO
//  const char* basedir="rfio:/castor/cern.ch/user/s/skowron/";
//  const char* serie="standard";
//  const char* field = "0.4";

  const char* basedir=".";
  const char* serie="";
  const char* field = "";
 
  Bool_t threeDcuts = kFALSE;

  // Dynamically link some shared libs                    
//  gROOT->LoadMacro("loadlibs.C");
//  loadlibs();
  
  gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libHBTAN");
  

  /***********************************************************/
  //Create analysis and reader
  AliHBTAnalysis * analysis = new AliHBTAnalysis();
  analysis->SetCutsOnTracks();
  
  AliHBTReader* reader;
  Int_t kine = strcmp(datatype,"Kine");
  Int_t ESD = strcmp(datatype,"ESD");
  Int_t TPC = strcmp(datatype,"TPC");
  Int_t ITSv1 = strcmp(datatype,"ITSv1");
  Int_t ITSv2 = strcmp(datatype,"ITSv2");
  Int_t intern = strcmp(datatype,"Intern");
  
  if(!kine)
   {
    reader = new AliHBTReaderKineTree();
    processopt="Particles"; //this reader by definition reads only simulated particles
   }
  else if(!ESD)
   {
    AliHBTReaderESD* esdreader = new AliHBTReaderESD();
    esdreader->ReadParticles(kTRUE);
    esdreader->SetNumberOfTrackPoints(5,30.);
    esdreader->SetClusterMap();
    reader = esdreader;
   }
  else if(!TPC)
   {
    reader = new AliHBTReaderTPC();
    //((AliHBTReaderTPC*)reader)->SetNumberOfTrackPoints(5,30.);
    ((AliHBTReaderTPC*)reader)->SetClusterMap();
   }
  else if(!ITSv1)
   {
    reader = new AliHBTReaderITSv1();
   }
  else if(!ITSv2)
   {
    reader = new AliHBTReaderITSv2();
   }
  else if(!intern)
   {
    reader = new AliHBTReaderInternal("data.root");
   }
  else
   {
    cerr<<"Option "<<datatype<<"  not recognized. Exiting"<<endl;
    return;
   }

  TObjArray* dirs=0;
  if ( (first >= 0) && (last>=0) && ( (last-first)>=0 ) )
   {//read from many dirs dirs
     char buff[50];
     dirs = new TObjArray(last-first+1);
     for (Int_t i = first; i<=last; i++)
      {
        sprintf(buff,"%s/%s/%s/%d",basedir,field,serie,i);
        TObjString *odir= new TObjString(buff);
        dirs->Add(odir);
      }
    }
   reader->SetDirs(dirs);
   
  /***********************************************************/    
  /*****   R E A D E R ***************************************/    
  /***********************************************************/    
  
  //we want have only low pt pi+ so set a cut to reader
  AliHBTParticleCut* readerpartcut= new AliHBTParticleCut();
  readerpartcut->SetPtRange(0.0,10.0);
  readerpartcut->SetPID(kPiPlus);
  reader->AddParticleCut(readerpartcut);//read this particle type with this cut
  
//  readerpartcut->SetPtRange(0.0,10.0);
//  readerpartcut->SetPID(kPiMinus);
//  reader->AddParticleCut(readerpartcut);//read this particle type with this cut
  
  analysis->SetReader(reader);

  AliHBTPairCut *paircut = new AliHBTPairCut();
  Float_t qinvmin = 0.0;
  Float_t qinvmax = 0.05;//50MeV
  paircut->SetQInvRange(qinvmin,qinvmax);  
  
 // paircut->SetAvSeparationRange(10.); //AntiMerging
  paircut->SetClusterOverlapRange(-.5,1.0);//AntiSplitting
  
//  AliHBTParticleCut* partcut= new AliHBTParticleCut();
//  partcut->SetPID(kPiPlus);
//  paircut->SetFirstPartCut(partcut);
//  partcut->SetPID(kPiMinus);
//  paircut->SetSecondPartCut(partcut);
  
  analysis->SetGlobalPairCut(paircut);
  
  /***********************************************************/    
  /*****   W E I G H T S        ******************************/    
  /***********************************************************/    

  AliHBTLLWeights::Instance()->SetParticlesTypes(kPiPlus,kPiPlus);
  AliHBTLLWeights::Instance()->SetColoumb(kFALSE);
  AliHBTLLWeights::Instance()->SetStrongInterSwitch(kFALSE);

  AliHBTLLWeights::Instance()->SetRandomPosition(kFALSE);
//AliHBTLLWeights::Instance()->SetR1dw(8.0);

  AliHBTLLWeights::Instance()->Init();
  AliHBTLLWeights::Instance()->Set();

  //example function
  AliHBTWeightTheorQInvFctn  *wqinvcfP = new AliHBTWeightTheorQInvFctn(100,qinvmax);
  wqinvcfP->SetNumberOfBinsToScale(30);
  wqinvcfP->Rename("wqinvcfP","Lednicky Q_{inv} Theoretical Correlation Function");
  analysis->AddParticleFunction(wqinvcfP);

  /************************************************************/
  /****   Q INV Correlation Function   ************************/
  /************************************************************/
   

  AliHBTQInvCorrelFctn * qinvcfT = new AliHBTQInvCorrelFctn(100,qinvmax);
  AliHBTQInvCorrelFctn * qinvcfP = new AliHBTQInvCorrelFctn(100,qinvmax);
  
  analysis->AddTrackFunction(qinvcfT);
  analysis->AddParticleFunction(qinvcfP);
  
  
  qinvcfP->Rename("qinvcfP","Particle (simulated) Qinv CF \\pi^{+} \\pi^{+}");
  qinvcfT->Rename("qinvcfT","Track (recontructed) Qinv CF \\pi^{+} \\pi^{+}");
  
  /************************************************************/
  /****   Q OUT Correlation Function   ************************/
  /************************************************************/
  
  AliHBTQOutLCMSCorrelFctn* qoutP = new AliHBTQOutLCMSCorrelFctn();
  qoutP->Rename("qoutP","Particle (simulated) Q_{out} CF \\pi^{+} \\pi^{+}");
  AliHBTQOutLCMSCorrelFctn* qoutT = new AliHBTQOutLCMSCorrelFctn(); 
  qoutT->Rename("qoutT","Track (recontructed) Q_{out} CF \\pi^{+} \\pi^{+}");

  AliHBTPairCut *outPairCut = new AliHBTPairCut();
  outPairCut->SetQOutCMSLRange(0.0,0.15);
  outPairCut->SetQSideCMSLRange(0.0,0.02);
  outPairCut->SetQLongCMSLRange(0.0,0.02);
  qoutP->SetPairCut(outPairCut);
  qoutT->SetPairCut(outPairCut);

  //analysis->AddTrackFunction(qoutP);
  //analysis->AddParticleFunction(qoutT);

  /************************************************************/
  /****   Q SIDE Correlation Function   ***********************/
  /************************************************************/
  
  AliHBTQSideLCMSCorrelFctn* qsideP = new AliHBTQSideLCMSCorrelFctn(100,qinvmax); 
  qsideP->Rename("qsideP","Particle (simulated) Q_{side} CF \\pi^{+} \\pi^{+}");
  AliHBTQSideLCMSCorrelFctn* qsideT = new AliHBTQSideLCMSCorrelFctn(100,qinvmax); 
  qsideT->Rename("qsideT","Track (recontructed) Q_{side} CF \\pi^{+} \\pi^{+}");

  AliHBTPairCut *sidePairCut = new AliHBTPairCut();
  sidePairCut->SetQOutCMSLRange(0.0,0.02);
  sidePairCut->SetQSideCMSLRange(0.0,0.15);
  sidePairCut->SetQLongCMSLRange(0.0,0.02);
  qsideP->SetPairCut(sidePairCut);
  qsideT->SetPairCut(sidePairCut);

  //analysis->AddTrackFunction(qsideP);
  //analysis->AddParticleFunction(qsideT);

  /************************************************************/
  /****   Q LONG Correlation Function   ***********************/
  /************************************************************/
    
  AliHBTQLongLCMSCorrelFctn* qlongP = new AliHBTQLongLCMSCorrelFctn(100,qinvmax);
  qlongP->Rename("qlongP","Particle (simulated) Q_{long} CF \\pi^{+} \\pi^{+}");
  AliHBTQLongLCMSCorrelFctn* qlongT = new AliHBTQLongLCMSCorrelFctn(100,qinvmax); 
  qlongT->Rename("qlongT","Track (recontructed) Q_{long} CF \\pi^{+} \\pi^{+}");

  AliHBTPairCut *longPairCut = new AliHBTPairCut();
  longPairCut->SetQOutCMSLRange(0.0,0.02);
  longPairCut->SetQSideCMSLRange(0.0,0.02);
  longPairCut->SetQLongCMSLRange(0.0,0.15);
  qlongP->SetPairCut(longPairCut);
  qlongT->SetPairCut(longPairCut);

  //analysis->AddTrackFunction(qlongT);
  //analysis->AddParticleFunction(qlongT);


  /************************************************************/
  /************************************************************/
  /************************************************************/
  /****   R E S O L U T I O N S         ***********************/
  /************************************************************/
  /************************************************************/
  /************************************************************/



  AliHBTQInvResolVsKtFctn qinvVsktF(200,1.0,0.0,300,0.015,-0.015);    //QInvLCMS  Res   Vs   Kt
  if (threeDcuts)  qinvVsktF.SetPairCut(paircut);
  //analysis->AddResolutionFunction(&qinvVsktF);

  AliHBTQOutResolVsKtFctn qoutVsktF(200,1.0,0.0,300,0.05,-0.05);    //QOutLCMS  Res   Vs   Kt
  if (threeDcuts)  qoutVsktF.SetPairCut(outPairCut);
  //analysis->AddResolutionFunction(&qoutVsktF);

  AliHBTQSideResolVsKtFctn qsideVsktF(200,1.0,0.0,300,0.015,-0.015);   //QSideLCMS Res   Vs   Kt
  if (threeDcuts)  qsideVsktF.SetPairCut(sidePairCut);
  //analysis->AddResolutionFunction(&qsideVsktF);

  AliHBTQLongResolVsKtFctn qlongVsktF(200,1.0,0.0,300,0.015,-0.015);   //QLongLCMS Res   Vs   Kt
  if (threeDcuts)  qlongVsktF.SetPairCut(longPairCut);
  //analysis->AddResolutionFunction(&qlongVsktF);

  AliHBTQInvResolVsQInvFctn qInvVsqinvF(200,qinvmax,qinvmin,300,0.015,-0.015); //QInvLCMS Res   Vs   QInvLCMS
  //analysis->AddResolutionFunction(&qInvVsqinvF);

  AliHBTQOutResolVsQInvFctn qoutVsqinvF(200,qinvmax,qinvmin,300,0.05,-0.05);  //QOutLCMS  Res   Vs   QInvLCMS
  if (threeDcuts)  qoutVsqinvF.SetPairCut(outPairCut);
  //analysis->AddResolutionFunction(&qoutVsqinvF);

  AliHBTQSideResolVsQInvFctn qsideVsqinvF(200,qinvmax,qinvmin,300,0.015,-0.015); //QSideLCMS Res   Vs   QInvLCMS
  if (threeDcuts)  qsideVsqinvF.SetPairCut(sidePairCut);
  //analysis->AddResolutionFunction(&qsideVsqinvF);

  AliHBTQLongResolVsQInvFctn qlongVsqinvF(200,qinvmax,qinvmin,300,0.015,-0.015); //QLongLCMS Res   Vs   QInvLCMS
  if (threeDcuts)  qlongVsqinvF.SetPairCut(longPairCut);
  //analysis->AddResolutionFunction(&qlongVsqinvF);


  AliHBTPairThetaResolVsQInvFctn pairThetaVsqinv(200,qinvmax,qinvmin,300,0.015,-0.015);
  if (threeDcuts) pairThetaVsqinv.SetPairCut(longPairCut);
  //analysis->AddResolutionFunction(&pairThetaVsqinv);

  AliHBTPairThetaResolVsKtFctn pairThetaVsKt(200,1.0,0.0,300,0.015,-0.015);
  if (threeDcuts) pairThetaVsKt.SetPairCut(longPairCut);
  //analysis->AddResolutionFunction(&pairThetaVsKt);

  AliHBTPairCut phipc;
  phipc.SetQLongCMSLRange(0.0,0.02);

  AliHBTPairPhiResolVsQInvFctn pairPhiVsqinv(200,qinvmax,qinvmin,300,0.015,-0.015);
  if (threeDcuts) pairPhiVsqinv.SetPairCut(&phipc);
  //analysis->AddResolutionFunction(&pairPhiVsqinv);

  AliHBTPairPhiResolVsKtFctn pairPhiVsKt(200,1.0,0.0,300,0.015,-0.015);
  if (threeDcuts) pairPhiVsKt.SetPairCut(&phipc);
  //analysis->AddResolutionFunction(&pairPhiVsKt);


  AliHBTQOutResolVsQOutFctn qoutVsqoutF;  //QOutLCMS  Res   Vs   QOut
  if (threeDcuts) qoutVsqoutF.SetPairCut(outPairCut);
  //analysis->AddResolutionFunction(&qoutVsqoutF);

  AliHBTQSideResolVsQSideFctn qsideVsqsideF;//QSideLCMS Res   Vs   QSide
  if (threeDcuts) qsideVsqsideF.SetPairCut(sidePairCut);
  //analysis->AddResolutionFunction(&qsideVsqsideF);

  
  AliHBTQLongResolVsQLongFctn qlongVsqlongF;//QLongLCMS Res   Vs   QLong
  if (threeDcuts) qlongVsqlongF.SetPairCut(longPairCut);
  //analysis->AddResolutionFunction(&qlongVsqlongF);


  /************************************************************/
  /************************************************************/
  /************************************************************/
  /****   D I S T R B U T I O N S         ***********************/
  /************************************************************/
  /************************************************************/
  /************************************************************/

  /************************************************************/
  /****   OUT VS KT                 ***********************/
  /************************************************************/
  AliHBTQOutDistributionVsKtFctn qoutVsKtDistP;
  qoutVsKtDistP.Rename("qoutVsKtDistP",qoutVsKtDistP.GetTitle());
  AliHBTQOutDistributionVsKtFctn qoutVsKtDistT;
  qoutVsKtDistT.Rename("qoutVsKtDistT",qoutVsKtDistT.GetTitle());
  if (threeDcuts)  qoutVsKtDistP.SetPairCut(outPairCut);
  if (threeDcuts)  qoutVsKtDistT.SetPairCut(outPairCut);
  //analysis->AddParticleFunction(&qoutVsKtDistP);
  //analysis->AddTrackFunction(&qoutVsKtDistT);
  /************************************************************/
  /****   Side VS KT                 ***********************/
  /************************************************************/
  AliHBTQSideDistributionVsKtFctn qSideVsKtDistP;
  qSideVsKtDistP.Rename("qSideVsKtDistP",qSideVsKtDistP.GetTitle());
  AliHBTQSideDistributionVsKtFctn qSideVsKtDistT;
  qSideVsKtDistT.Rename("qSideVsKtDistT",qSideVsKtDistT.GetTitle());
  if (threeDcuts)  qSideVsKtDistP.SetPairCut(sidePairCut);
  if (threeDcuts)  qSideVsKtDistT.SetPairCut(sidePairCut);
  //analysis->AddParticleFunction(&qSideVsKtDistP);
  //analysis->AddTrackFunction(&qSideVsKtDistT);
  /************************************************************/
  /****   Long VS KT                 ***********************/
  /************************************************************/
  AliHBTQLongDistributionVsKtFctn qLongVsKtDistP;
  qLongVsKtDistP.Rename("qLongVsKtDistP",qLongVsKtDistP.GetTitle());
  AliHBTQLongDistributionVsKtFctn qLongVsKtDistT;
  qLongVsKtDistT.Rename("qLongVsKtDistT",qLongVsKtDistT.GetTitle());
  if (threeDcuts)  qLongVsKtDistP.SetPairCut(longPairCut);
  if (threeDcuts)  qLongVsKtDistT.SetPairCut(longPairCut);
  //analysis->AddParticleFunction(&qLongVsKtDistP);
  //analysis->AddTrackFunction(&qLongVsKtDistT);

  /************************************************************/
  /************************************************************/
  /************************************************************/
  /****    M O N I T O R  R E S O L U T I O N S   *************/
  /************************************************************/
  /************************************************************/
  /************************************************************/
  
  AliHBTMonPxResolutionVsPtFctn PxVsPtResF(200,1.0,0.0,300,0.05,-0.05);
  //analysis->AddParticleAndTrackMonitorFunction(&PxVsPtResF);
  
  AliHBTMonPyResolutionVsPtFctn PyVsPtResF(200,1.0,0.0,300,0.05,-0.05);
  //analysis->AddParticleAndTrackMonitorFunction(&PyVsPtResF);
  
  AliHBTMonPzResolutionVsPtFctn PzVsPtResF(200,1.0,0.0,300,0.05,-0.05);
  //analysis->AddParticleAndTrackMonitorFunction(&PzVsPtResF);

  AliHBTMonPResolutionVsPtFctn PVsPtResF(200,1.0,0.0,300,0.05,-0.05);
  //analysis->AddParticleAndTrackMonitorFunction(&PVsPtResF);

  AliHBTMonPtResolutionVsPtFctn PtVsPtResF(400,4.0,0.0,300,0.07,-0.07);
  analysis->AddParticleAndTrackMonitorFunction(&PtVsPtResF);

  AliHBTMonPhiResolutionVsPtFctn PhiVsPtResF(200,1.0,0.0,300,0.02,-0.02);
  analysis->AddParticleAndTrackMonitorFunction(&PhiVsPtResF);

  AliHBTMonThetaResolutionVsPtFctn ThetaVsPtResF(200,1.0,0.0,300,0.02,-0.02);
  analysis->AddParticleAndTrackMonitorFunction(&ThetaVsPtResF);
  
  
  
  
  /************************************************************/
  /****   P R O C E S S                 ***********************/
  /************************************************************/
  analysis->SetBufferSize(2);
  analysis->Process(processopt);
  
  TFile histoOutput(outfile,"recreate"); 
  analysis->WriteFunctions();
  histoOutput.Close();
  
  delete qinvcfP;
  delete qinvcfT;
  delete paircut;
  delete readerpartcut;
  if (dirs) 
   {
    dirs->SetOwner();
    delete dirs;
   }
  delete reader;
  delete analysis;
  
 }

