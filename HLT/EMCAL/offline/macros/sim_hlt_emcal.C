
void sim_hlt_calo(const char *input = "./", const char *grp = "./", bool doPhos = true, bool doEmcal = true, bool doTM = true)
{

    //AliCDBManager::Instance()->SetRun(0);

 //   if (!gSystem->AccessPathName("galice.root")) {
   if(0){
        cerr << "please delete the galice.root or run at different place." << endl;
        return;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // init the HLT system in order to define the analysis chain below
    //
    AliHLTSystem* gHLT=AliHLTPluginBase::GetInstance();

    // Input to the ESD converter
    TString ecInput;
    TString option="libAliHLTUtil.so libAliHLTRCU.so libAliHLTTPC.so libAliHLTCalo.so libAliHLTPHOS.so libAliHLTEMCAL.so libAliHLTGlobal.so loglevel=0x7f chains=ESD-CONVERTER";

 const char* cdbEntry="GRP/Geometry/Data";
  AliCDBManager* pMan=AliCDBManager::Instance();
  if (pMan) {
    if (!pMan->IsDefaultStorageSet()) {
      pMan->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
      pMan->SetRun(0);
    }
    AliCDBEntry *pEntry = pMan->Get(cdbEntry);
    if (pEntry && 
	pEntry->GetObject()) {
    } else {
      HLTWarning("can not load CDB entry %s", cdbEntry);
    }
  }

//    if (doPhos)
    if(0)
    {
        ///////////////////////////////////////
        // The PHOS part of the chain
        //////////////////////////////////////

        int moduleStart = 2;
        int moduleEnd = 4;
        int rcuStart = 0;
        int rcuEnd = 3;
        int rcusPerModule = 4;
        int ddlOffset = 1792;

        for (int module = moduleStart; module <= moduleEnd; module++)
        {
            TString clInput;

             TString arg, cl;

            // Clusterizer
            cl.Form("PHS-CL_%02d", module);
            arg = "";
            AliHLTConfiguration clConf(cl.Data(), "PhosClusterizer", clInput.Data(), arg.Data());

            if (ecInput.Length() > 0) ecInput += " ";
            ecInput += cl;
        }

        // END OF PHOS
    }

    if (doEmcal)
    {
        int moduleStart = 0;
        int moduleEnd = 0;
        int rcuStart = 0;
        int rcuEnd = 1;
        int rcusPerModule = 2;
        int ddlOffset = 4608;

        histoInput = "";

        for (int module = moduleStart; module <= moduleEnd; module++)
        {
            TString clInput;

            for (int rcu = rcuStart; rcu <= rcuEnd; rcu++)
            {
 /*               TString arg, publisher, ra, dm;
                // raw data publisher components
                publisher.Form("EMC-RP_%02d_%d", module, rcu);
                arg.Form("-minid %d -datatype 'DDL_RAW ' 'EMCA'  -dataspec 0x%x ", ddlOffset + module*(rcusPerModule) + rcu, 0x1 << (module*rcusPerModule + rcu));
                AliHLTConfiguration pubConf(publisher.Data(), "AliRawReaderPublisher", NULL , arg.Data());
		
		
                // Raw analyzer
                arg = "";
                ra.Form("EMC-RA_%02d_%d", module, rcu);
                AliHLTConfiguration rawConf(ra.Data(), "EmcalRawCrude", publisher.Data(), arg.Data());

                // digit maker components
                dm.Form("EMC-DM_%02d_%d", module, rcu);
                arg="";
                AliHLTConfiguration dmConf(dm.Data(), "EmcalDigitMaker", ra.Data(), arg.Data());

                if (clInput.Length() > 0) clInput += " ";
                clInput+=dm;
  */          }
            
            TString arg, publisher, cl, ca;
	    publisher.Form("EMC-DP_%02d", module);
	    arg.Form("-detector EMCAL -module %d", module);
	    AliHLTConfiguration digPubConf(publisher.Data(), "CaloDigitPublisher", NULL, arg.Data());
	    
            // Clusterizer
            cl.Form("EMC-CL_%02d", module);
            arg = "";
            AliHLTConfiguration clConf(cl.Data(), "EmcalClusterizer", publisher.Data(), arg.Data());

            if (ecInput.Length() > 0) ecInput += " ";
            ecInput += cl;
        }


        // END OF EMCAL
    }

    // If there are no tracks it shouldn't do anything...
    AliHLTConfiguration tmconf("track-matcher", "TrackMatcher", ecInput.Data(), "");

    //if (doEmcal)
    if(0)
    {
        // EMCAL Histograms
        AliHLTConfiguration hconf("emcalHistocomp", "CaloPhysicsHistos", "track-matcher", "-emcal -invariantmass -clusterenergy -matchedtracks");
        AliHLTConfiguration fwconf("emcalHist", "ROOTFileWriter"   , "emcalHistocomp", "-datafile emcalHistograms -concatenate-events -overwrite");
    }

    //if (doPhos)
    if(0)
    {
        // PHOS Histograms
        AliHLTConfiguration hconf("phosHistocomp", "CaloPhysicsHistos", "track-matcher", "-phos -invariantmass -clusterenergy -matchedtracks");
        AliHLTConfiguration fwconf("phosHist", "ROOTFileWriter"   , "phosHistocomp", "-datafile phosHistograms -concatenate-events -overwrite");
    }
    
    TString arg, ec;


    ec.Form("ESD-CONVERTER");
    arg = "";

    AliHLTConfiguration esdcconf(ec.Data(), "GlobalEsdConverter"   , "track-matcher", "");
    
    AliSimulation sim;
        sim.SetRunHLT("libAliHLTUtil.so libAliHLTRCU.so libAliHLTTPC.so libAliHLTCalo.so libAliHLTPHOS.so libAliHLTEMCAL.so libAliHLTGlobal.so loglevel=0x7f chains=ESD-CONVERTER");
        sim.SetRunGeneration(kFALSE);
        sim.SetMakeDigits("");
        sim.SetMakeSDigits("");
        //sim.SetMakeDigits("EMCAL");
        //sim.SetMakeSDigits("EMCAL");
        sim.SetMakeDigitsFromHits("");
	sim.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
	sim.SetSpecificStorage("GRP/GRP/Data", Form("local://%s",gSystem->pwd()));

        sim.Run(5);

    
}
