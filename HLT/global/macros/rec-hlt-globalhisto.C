// $Id$
/**
 * @file rec-hlt-globalhisto.C
 * @brief Run reconstruction and fill the global histograms
 *
 * <pre>
 * Usage: aliroot -b -q -l \
 *     rec-hlt-globalhisto.C'("file", "cdb", minEvent, maxEvent)'
 * Example:
 *   aliroot -b -q -l \
 *     rec-hlt-globalhisto.C'("raw.root", "local://$ALICE_ROOT/OCDB")'
 * </pre>
 * 
 * Macro runs HLT reconstruction on raw data, all other AliRoot
 * modules switched off. The output of the GlobalEsdConverter is filled
 * into the GlobalHisto component, it's output is attached to a Root file
 * writer producing the file histo.root.
 *
 * A raw file, either simulated or real, is needed and the corresponding
 * OCDB has to be specified.
 *
 * Can also run on an AliESDs.root file instead of running the reconstruction
 *
 * aliroot -b -q -l \
 *   rec-hlt-globalhisto.C'("./","HLTesdTree","local://$ALICE_ROOT/OCDB/",5)' | tee log
 *
 * @author Matthias.Richter@ift.uib.no
 * @ingroup alihlt_qa
 */
void rec_hlt_globalhisto(const char *filename,
		         const char *cdbURI,
		         int minEvent=-1,
		         int maxEvent=-1)
{
  // connect to the GRID if we use a file or OCDB from the GRID
  TString struri=cdbURI;
  TString strfile=filename;
  if (struri.BeginsWith("raw://") ||
      strfile.Contains("://") && !strfile.Contains("local://")) {
    TGrid::Connect("alien");
  }

  // Set the CDB storage location
  AliCDBManager * man = AliCDBManager::Instance();
  man->SetDefaultStorage(cdbURI);
  if (struri.BeginsWith("local://")) {
    // set specific storage for GRP entry
    // search in the working directory and one level above, the latter
    // follows the standard simulation setup like e.g. in test/ppbench
    if (!gSystem->AccessPathName("GRP/GRP/Data")) {
      man->SetSpecificStorage("GRP/GRP/Data", "local://$PWD");
    } else if (!gSystem->AccessPathName("../GRP/GRP/Data")) {
      man->SetSpecificStorage("GRP/GRP/Data", "local://$PWD/..");      
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////
  //
  // setup the HLT system
  AliHLTSystem* pHLT=AliHLTPluginBase::GetInstance();

  AliHLTConfiguration globalhisto("globalhisto", "GlobalHisto", "GLOBAL-esd-converter",
				  "-histogram TrackPt    -size 1000 -expression Track_pt   "
				//"-histogram TrackPt    -size 1000 -expression Track_pt  -cut Track_pt<0.8&&Track_pt>0.6  "
				  "-histogram TrackPhi   -size 1000 -expression Track_phi   "
				  "-histogram TrackEta   -size 1000 -expression Track_eta   "
				  "-histogram TrackCount -size 1000 -expression trackcount  "
				  "-histogram TrackTheta -size 1000 -expression Track_theta "
				  "-histogram VertexX    -size 1000 -expression vertexX     "
				  "-histogram VertexY    -size 1000 -expression vertexY     "
				  "-histogram VertexZ    -size 1000 -expression vertexZ     "
				  "-histogram DCAr       -size 1000 -expression Track_DCAr  "
				  "-histogram DCAz       -size 1000 -expression Track_DCAz  "
				  "-histogram dEdx_vs_p  -size 1000 -expression Track_dEdx:Track_p  "
				  );
  AliHLTConfiguration writer("writer", "ROOTFileWriter", "globalhisto", "-datafile histo.root -overwrite -concatenate-events");

  TString hltOptions="loglevel=0x7c ignore-hltout chains=writer";

  // Reconstruction settings
  AliReconstruction rec;

  if (minEvent>=0 || maxEvent>minEvent) {
    if (minEvent<0) minEvent=0;
    if (maxEvent<minEvent) maxEvent=minEvent;
    rec.SetEventRange(minEvent,maxEvent);
  }

  rec.SetRunReconstruction("HLT");

  // QA options
  TString qaOptions="HLT";
  qaOptions+=":ALL";
  rec.SetRunQA(qaOptions) ;
  //rec.SetQARefDefaultStorage("local://$ALICE_ROOT/QAref") ;

  // AliReconstruction settings
  rec.SetWriteESDfriend(kTRUE);
  rec.SetRunVertexFinder(kFALSE);
  rec.SetRunMultFinder(kFALSE);
  rec.SetInput(filename);
  rec.SetOption("HLT", hltOptions);

  rec.SetRunPlaneEff(kFALSE);

  // switch off cleanESD
  rec.SetCleanESD(kFALSE);

  AliLog::Flush();
  rec.Run();

 }

void rec_hlt_globalhisto(const char *esdfilename,
		         const char *treename,
			 const char *cdbURI,
		         int nofEvents=-1
			 )

{
  // check the name of the tree
  TString strtree=treename;
  if (strtree.CompareTo("esdTree")==0) strtree="ESD";
  else if (strtree.CompareTo("HLTesdTree")==0) strtree="HLTESD";
  else {
    cerr << "invalid treename '" << treename << "', supported 'esdTree' and 'HLTesdTree'" << endl;
    return;
  }

  // connect to the GRID if we use a file or OCDB from the GRID
  TString struri=cdbURI;
  TString strfile=esdfilename;
  if (struri.BeginsWith("raw://") ||
      strfile.Contains("://") && !strfile.Contains("local://")) {
    TGrid::Connect("alien");
  }

  // open the ESD file and get the event count
  if (!strfile.EndsWith("/")) strfile+="/";
  strfile+="AliESDs.root";
  TFile* esdfile=TFile::Open(strfile);
  if (!esdfile || esdfile->IsZombie()) {
    cerr << "can not open file " << strfile << endl;
    return;
  }

  // get number of events
  TTree* pTree=NULL;
  esdfile->GetObject(treename, pTree);
  if (!pTree) {
    cerr << "can not find " << treename << " in file " << strfile << endl;
    return;
  }
  if (pTree->GetEntries()<=0) {
    cerr << "empty tree " << treename << " in file " << strfile << endl;
    return;
  }
  
  AliESDEvent* esd=new AliESDEvent;
  esd->ReadFromTree(pTree);
  pTree->GetEntry(0);

  if (nofEvents<0 || nofEvents>pTree->GetEntries())
    nofEvents=pTree->GetEntries();

  // Set the CDB storage location
  AliCDBManager * man = AliCDBManager::Instance();
  man->SetDefaultStorage(cdbURI);
  man->SetRun(esd->GetRunNumber());
  if (struri.BeginsWith("local://")) {
    // set specific storage for GRP entry
    // search in the working directory and one level above, the latter
    // follows the standard simulation setup like e.g. in test/ppbench
    if (!gSystem->AccessPathName("GRP/GRP/Data")) {
      man->SetSpecificStorage("GRP/GRP/Data", "local://$PWD");
    } else if (!gSystem->AccessPathName("../GRP/GRP/Data")) {
      man->SetSpecificStorage("GRP/GRP/Data", "local://$PWD/..");      
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////
  //
  // setup the HLT system
  AliHLTSystem *pHLT = AliHLTPluginBase::GetInstance();
  
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libPWG0base");
  gSystem->Load("libPWG0dep");

  TString arguments, histoinput;

  // ESD publisher component
  arguments=" -datapath "; arguments+=esdfilename;
  arguments+=" -entrytype "; arguments+=strtree;
  histoinput="ESD-publisher";

  AliHLTConfiguration esdpublisher(histoinput.Data(), "ESDMCEventPublisher", "", arguments.Data());

  //AliHLTConfiguration v0histo("v0","V0Histo",histoinput.Data(),"");
  
  AliHLTConfiguration globalhisto("globalhisto", "GlobalHisto", histoinput.Data(),
        			"-histogram TrackPt    -size 1000 -expression Track_pt "
        		        //"-histogram TrackPt    -size 1000 -expression Track_pt  -cut Track_pt<0.8&&Track_pt>0.6  "
         			"-histogram TrackPhi   -size 1000 -expression Track_phi "
         			"-histogram TrackEta   -size 1000 -expression Track_eta "
         			"-histogram TrackCount -size 1000 -expression trackcount "
				"-histogram Track_Nclusters -size 1000 -expression Track_Nclusters -cut Track_Nclusters>0 "
         			"-histogram TrackTheta -size 1000 -expression Track_theta "
         			"-histogram VertexX    -size 1000 -expression vertexX -cut nContributors>3 "
         			"-histogram VertexY    -size 1000 -expression vertexY -cut nContributors>3 "
         			"-histogram VertexZ    -size 1000 -expression vertexZ -cut nContributors>3 "
         			"-histogram DCAr       -size 1000 -expression Track_DCAr "
         			"-histogram DCAz       -size 1000 -expression Track_DCAz "
         			"-histogram dEdx_vs_p  -size 1000 -expression Track_dEdx:Track_p "
				"-histogram trendVertexX -size 1000 -expression vertexX:event -cut nContributors>3 "
				"-histogram trendVertexY -size 1000 -expression vertexY:event -cut nContributors>3 "
				"-histogram trendVertexZ -size 1000 -expression vertexZ:event -cut nContributors>3 "
				"-histogram XY -size 1000 -expression vertexX:vertexY -cut nContributors>3 "
				"-histogram DCA -size 1000 -expression TMath::Sqrt(Track_DCAr*Track_DCAr+Track_DCAz*Track_DCAz) "
				"-histogram V0_AP  -size 1000 -expression V0_pt:V0_AP -cut clust1>=60&&clust2>=60&&dev1>=2.5&&dev2>=2.5&&devPrim<=3.5&&length>=3.0*sigmaLength&&length>=0&&r<=50"				
				);
  
  AliHLTConfiguration writer("writer1", "ROOTFileWriter", "globalhisto", "-datafile   histo.root -overwrite -concatenate-events");
  //AliHLTConfiguration writer("writer2", "ROOTFileWriter", "v0",          "-datafile v0histo.root -overwrite -concatenate-events");

  // set option for the HLT system
  // arguments
  //  - libraries to be used as plugins
  //  - loglevel=0x79 : Important, Warning, Error, Fatal
  pHLT->ScanOptions("libAliHLTUtil.so libAliHLTGlobal.so loglevel=0x79");

  pHLT->BuildTaskList("writer1");
  //pHLT->BuildTaskList("writer2");
  pHLT->Run(nofEvents);
}
	
