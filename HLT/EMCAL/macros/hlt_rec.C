void hlt_rec(const char* input="./") {

  // For real data:
 // AliCDBManager::Instance()->SetDefaultStorage("raw://");
 
  gStyle->SetPalette(1);

  if(!gSystem->AccessPathName("galice.root")){
    cerr << "please delete the galice.root or run at different place." << endl;
    return;
  }
 
  if (!input) {
    cerr << "please specify input or run without arguments" << endl;
    return;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // init the HLT system in order to define the analysis chain below
  //
  
  AliHLTSystem* gHLT=AliHLTPluginBase::GetInstance();

  gHLT->LoadComponentLibraries("libAliHLTUtil.so");
  gHLT->LoadComponentLibraries("libAliHLTRCU.so");
  gHLT->LoadComponentLibraries("libAliHLTCalo.so");
  gHLT->LoadComponentLibraries("libAliHLTEMCAL.so");
  gHLT->LoadComponentLibraries("libAliHLTPHOS.so");
  gHLT->LoadComponentLibraries("libAliHLTGlobal.so");

  //////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // define the analysis chain to be run
  //

  // Define which modules and RCUs to use
  int moduleStart = 0;
  int moduleEnd = 3;
  int rcuStart = 0;
  int rcuEnd = 1;
  
  // The number of RCUs per module in EMCAL (used for calculating data specification)

  Int_t rcusPerModule = 2;
  Int_t ddlOffset = 4608; // The DDL offset for EMCAL (for PHOS the number is 1792)
  
  TString arg, fwInput, ecInput, dm;

  // Loop over the modules/RCUs
  for (int module = moduleStart; module <= moduleEnd; module++) 
    {
      TString clInput;
      
      for(int rcu = rcuStart; rcu <= rcuEnd; rcu++) 
	{
	  TString publisher, rawanalyzer, dm;

	  // Raw data publisher
	  publisher.Form("EMC-RP_%02d_%d", module, rcu);
	  
	 arg.Form("-verbose -minid %d -datatype 'DDL_RAW ' 'EMCA'  -dataspec 0x%x ", ddlOffset + module*(rcusPerModule) + rcu, 0x1 << (module*rcusPerModule + rcu));
	  AliHLTConfiguration pubConf(publisher.Data(), "AliRawReaderPublisher", NULL , arg.Data());
	  cout << arg.Data() << endl;

	  // Raw analyzer
	  arg = "";
	  rawanalyzer.Form("EMC-RA_%02d_%d", module, rcu);
	  AliHLTConfiguration rawConf(rawanalyzer.Data(), "EmcalRawCrude", publisher.Data(), arg.Data());
	   
	  // If you want to have a component subscribing to all the raw analyzers this will be it's input
	  if(fwInput.Length() > 0) fwInput += " ";
	  fwInput+=rawanalyzer;
	  cout << ">>>>>>>>>> fwInput:" << fwInput << endl;
	  
          // digit maker components
          dm.Form("EMC-DM_%02d_%d", module, rcu);
          arg="";
          arg.Form("-highgainfactor 0.005 -lowgainfactor 0.08");
          AliHLTConfiguration dmConf(dm.Data(), "EmcalDigitMaker", rawanalyzer.Data(), arg.Data());
          if(clInput.Length() > 0) clInput += " ";
          clInput+=dm;
	  cout << ">>>>>>>>>> clInput:" << clInput << endl;	
	
	} // end of RCU loop

   	TString cl, ca;

      	cl.Form("EMC-CL_%02d", module);
      	arg = "";
      	arg.Form("-digitthreshold 0.005 -recpointthreshold 0.1");
      	AliHLTConfiguration clConf(cl.Data(), "EmcalClusterizer", clInput.Data(), arg.Data());
      	if(ecInput.Length() > 0) ecInput += " ";
      	ecInput += cl;

    } // end of module loop

	TString ec;

     // The call the histo maker 
  	arg.Form("");
  	arg.Form("-pushfraction 5 -beverbose 1");
  	AliHLTConfiguration hfConf("emcalHisto","EmcalRawHistoMaker",fwInput.Data(), arg.Data());
	ecInput += " emcalHisto";

    	ec.Form("ESD-CONVERTER");
    	arg = "";
 	AliHLTConfiguration esdcconf(ec.Data(), "GlobalEsdConverter"   , ecInput.Data(), "");


  // The filewriter 
//  arg.Form("-datatype 'CHANNELT' 'EMCA' -datafile RAout.root -concatenate-blocks -concatenate-events");
//  AliHLTConfiguration fwConf("filewriter", "FileWriter", fwInput.Data(), arg.Data());
  
  // The call the histo maker 
//  arg.Form("");
//  arg.Form("-pushfraction 5 -beverbose 1");
//  AliHLTConfiguration hfConf("emcalHisto","EmcalRawHistoMaker",fwInput.Data(), arg.Data());
  
  // Write the root file 
  //AliHLTConfiguration rwConf("rootFileHisto","ROOTFileWriter", "emcalHisto ESD-CONVERTER","-datafile roothisto.root -concatenate-events -overwrite");



  
  //////////////////////////////////////////////////////////////////////
  //
  // Init and run the reconstruction
  // All but HLT reconstruction is switched off 
  //
  /////////////////////////////////////////////////////////////////////

  TString option="libAliHLTUtil.so libAliHLTRCU.so libAliHLTCalo.so libAliHLTEMCAL.so libAliHLTGlobal.so chains=";
  option+="ESD-CONVERTER";
  //option+="rootFileHisto loglevel=0x5f";

  AliReconstruction rec;
  Char_t outname[256];

  // uncomment for simulation
  rec.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  rec.SetSpecificStorage("GRP/GRP/Data", Form("local://%s",gSystem->pwd()));

  rec.SetRunReconstruction(":");
  rec.SetEventRange(0,10);
  rec.SetInput(input);
  
  //rec.SetInput("raw.root");
  rec.SetRunVertexFinder(kFALSE);

  rec.SetRunReconstruction("HLT");
  //rec.SetRunTracking(":");
  rec.SetLoadAlignFromCDB(0);
  rec.SetRunQA(":");

  rec.SetOption("HLT", option);
  
  rec.Run();

  TCanvas *c1 = new TCanvas("c1","",0,0,600,400);
  c1->Divide(2,2);

 // for (int i=0; i<3; i++) {
//  c1->cd(i+1);
//  fChannelETMap[i]->Draw("colz");

//  }

 c1->cd(1);
 fChannelETMap0->Draw("colz");
 c1->cd(2);
 fChannelETMap1->Draw("colz");
 c1->cd(3);
 fChannelETMap2->Draw("colz");
 c1->cd(4);
 fChannelETMap3->Draw("colz");

 sprintf(outname, "ChannleETMap.gif");
 c1->SaveAs(outname);


 TCanvas *c2 = new TCanvas("c2","",0,0,600,400);
 c2->Divide(2,2);


 c2->cd(1);
 fChannelEMap0->Draw("colz");
 c2->cd(2);
 fChannelEMap1->Draw("colz");
 c2->cd(3);
 fChannelEMap2->Draw("colz");
 c2->cd(4);
 fChannelEMap3->Draw("colz");

 sprintf(outname, "ChannleEMap.gif");
 c2->SaveAs(outname);

 TCanvas *c3 = new TCanvas("c3","",0,0,600,400);
 c3->Divide(2,2);

 c3->cd(1);
 fChannelTMap0->Draw("colz");
 c3->cd(2);
 fChannelTMap1->Draw("colz");
 c3->cd(3);
 fChannelTMap2->Draw("colz");
 c3->cd(4);
 fChannelTMap3->Draw("colz");

 sprintf(outname, "ChannleTMap.gif");
 c3->SaveAs(outname);

}
