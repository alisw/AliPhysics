void hbtanalysis(Int_t first = -1,Int_t last = -1, char *outfile = "hbtanalysis.root")
 {
//HBT Anlysis Macro
//Anlyzes TPC recontructed tracks and simulated particles that corresponds to them
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
 
  // Dynamically link some shared libs                    
  gROOT->LoadMacro("loadlibs.C");
  loadlibs();
  
  gSystem->Load("$(ALICE_ROOT)/lib/tgt_$(ALICE_TARGET)/libHBTAnalysis");
  

  /***********************************************************/
  //Create analysis and reader
  AliHBTAnalysis * analysis = new AliHBTAnalysis();
  
  AliHBTReader* reader;
  TObjArray* dirs=0;
  if ( (first < 0) || (last<0) || ( (last-first)<0 ) )
    {//read from current directory
     reader = new AliHBTReaderTPC();
    }
   else
    {//read from many dirs dirs
      char buff[50];
      dirs = new TObjArray(last-first+1);
      for (Int_t i = first; i<=last; i++)
       {
         sprintf(buff,"%s/%s/%s/%d",basedir,field,serie,i);
         TObjString *odir= new TObjString(buff);
         dirs->Add(odir);
       }
      reader = new AliHBTReaderTPC(dirs);
    }
   
  /***********************************************************/    
  
  //we want have only low pt pi+ so set a cut to reader
  AliHBTParticleCut* readerpartcut= new AliHBTParticleCut();
  readerpartcut->SetPtRange(0.0,3.0);
  readerpartcut->SetPID(kPiPlus);
  reader->AddParticleCut(readerpartcut);//read this particle type with this cut
  
  analysis->SetReader(reader);
  /************************************************************/
  
  AliHBTPairCut *paircut = new AliHBTPairCut();
  paircut->SetQInvRange(0.0,0.15);  
  analysis->SetGlobalPairCut(paircut);

  AliHBTQInvCorrelFctn * qinvcfT= new AliHBTQInvCorrelFctn();
  AliHBTQInvCorrelFctn * qinvcfP= new AliHBTQInvCorrelFctn();
  
  analysis->AddTrackFunction(qinvcfT);
  analysis->AddParticleFunction(qinvcfP);
  
  analysis->Process();
  
  qinvcfP->Rename("qinvcfP","Particle (simulated) Qinv CF \\pi^{+} \\pi^{+}");
  qinvcfT->Rename("qinvcfT","Track (recontructed) Qinv CF \\pi^{+} \\pi^{+}");

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

