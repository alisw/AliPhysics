
void rec_hlt_calo(const char *input = "./", const char *grp = "./", bool doPhos = true, bool doEmcal = true, bool doTM = true)
{

    AliCDBManager::Instance()->SetRun(0);

    if (!gSystem->AccessPathName("galice.root")) {
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


    if (doTM)
    {

        ///////////////////////////////////////
        // The TPC part of the chain
        //////////////////////////////////////

        int clusterFinderType=0; // 0 = v3; 1 = decoder; 2 = packed (offline v1)
        bool bUseCA=true;   // use the CA tracker and merger

        int iMinSlice=0;
        int iMaxSlice=35;
        int iMinPart=0;
        int iMaxPart=5;
        TString writerInput;
        TString mergerInput;
        TString histoInput;
        TString histogramHandlerInputClusterFinder;
        TString cdumpInput;
        for (int slice=iMinSlice; slice<=iMaxSlice; slice++) {
            TString trackerInput;
            for (int part=iMinPart; part<=iMaxPart; part++) {
                TString arg, publisher, cf;
                TString clusterHistoOutput;
                // raw data publisher components
                int ddlno=768;
                if (part>1) ddlno+=72+4*slice+(part-2);
                else ddlno+=2*slice+part;
                arg.Form("-minid %d -datatype 'DDL_RAW ' 'TPC '  -dataspec 0x%02x%02x%02x%02x -verbose", ddlno, slice, slice, part, part);
                publisher.Form("DP_%02d_%d", slice, part);
                AliHLTConfiguration pubconf(publisher.Data(), "AliRawReaderPublisher", NULL , arg.Data());

                // cluster finder components
                cf.Form("CF_%02d_%d", slice, part);
                if (clusterFinderType==1) {
                    AliHLTConfiguration cfconf(cf.Data(), "TPCClusterFinderDecoder", publisher.Data(), "-timebins 1001");
                } else if (clusterFinderType==2) {
                    AliHLTConfiguration cfconf(cf.Data(), "TPCClusterFinderPacked", publisher.Data(), "-timebins 1001 -sorted");
                } else {
                    AliHLTConfiguration cfconf(cf.Data(), "TPCClusterFinder32Bit", publisher.Data(), "");
                }
                if (trackerInput.Length()>0) trackerInput+=" ";
                trackerInput+=cf;
                if (writerInput.Length()>0) writerInput+=" ";
                writerInput+=cf;
                if (histoInput.Length()>0) histoInput+=" ";
                histoInput+=cf;
                if (cdumpInput.Length()>0) cdumpInput+=" ";
                cdumpInput+=cf;

                if (0) {
                    clusterHistoOutput.Form("CH_%02d_%d", slice, part);
                    AliHLTConfiguration cfconf(clusterHistoOutput.Data(), "TPCClusterHisto", cf.Data(), "");
                    if (histogramHandlerInputClusterFinder.Length()>0) histogramHandlerInputClusterFinder+=" ";
                    histogramHandlerInputClusterFinder+=clusterHistoOutput;
                }
            }
            TString tracker;
            // tracker components
            tracker.Form("TR_%02d", slice);
            if (bUseCA) {
                AliHLTConfiguration trackerconf(tracker.Data(), "TPCCATracker", trackerInput.Data(), "");
            } else {
                AliHLTConfiguration trackerconf(tracker.Data(), "TPCSliceTracker", trackerInput.Data(), "-pp-run");
            }

            if (writerInput.Length()>0) writerInput+=" ";
            writerInput+=tracker;
            if (mergerInput.Length()>0) mergerInput+=" ";
            mergerInput+=tracker;
            //add all slice tracks to histo input
            //if (histoInput.Length()>0) histoInput+=" ";
            //histoInput+=tracker;
        }

        // GlobalMerger component
        if (bUseCA) {
            AliHLTConfiguration mergerconf("globalmerger","TPCCAGlobalMerger",mergerInput.Data(),"");
        } else {
            AliHLTConfiguration mergerconf("globalmerger","TPCGlobalMerger",mergerInput.Data(),"");
        }



        //add all global tracks to histo input
        if (histoInput.Length()>0) histoInput+=" ";
        histoInput+="globalmerger";

        ecInput += " globalmerger";

        // END OF TPC
    }

    if (doPhos)
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

            for (int rcu = rcuStart; rcu <= rcuEnd; rcu++)
            {
                TString arg, publisher, ra, dm;
                // raw data publisher components
                publisher.Form("PHS-RP_%02d_%d", module, rcu);
                arg.Form("-minid %d -datatype 'DDL_RAW ' 'PHOS'  -dataspec 0x%x ", ddlOffset + module*(rcusPerModule) + rcu, 0x1 << (module*rcusPerModule + rcu));
                AliHLTConfiguration pubConf(publisher.Data(), "AliRawReaderPublisher", NULL , arg.Data());

                // Raw analyzer
                arg = "";
                ra.Form("PHS-RA_%02d_%d", module, rcu);
                AliHLTConfiguration rawConf(ra.Data(), "PhosRawCrude", publisher.Data(), arg.Data());

                // digit maker components
                dm.Form("PHS-DM_%02d_%d", module, rcu);
                arg="";
                AliHLTConfiguration dmConf(dm.Data(), "PhosDigitMaker", ra.Data(), arg.Data());

                if (clInput.Length() > 0) clInput += " ";
                clInput+=dm;
            }
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
                TString arg, publisher, ra, dm;
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
            }
            TString arg, cl, ca;

            // Clusterizer
            cl.Form("EMC-CL_%02d", module);
            arg = "";
            AliHLTConfiguration clConf(cl.Data(), "EmcalClusterizer", clInput.Data(), arg.Data());

            if (ecInput.Length() > 0) ecInput += " ";
            ecInput += cl;
        }


        // END OF EMCAL
    }

    // If there are no tracks it shouldn't do anything...
    AliHLTConfiguration tmconf("track-matcher", "TrackMatcher", ecInput.Data(), "");

    if (doEmcal)
    {
        // EMCAL Histograms
        AliHLTConfiguration hconf("emcalHistocomp", "CaloPhysicsHistos", "track-matcher", "-emcal -invariantmass -clusterenergy -matchedtracks");
        AliHLTConfiguration fwconf("emcalHist", "ROOTFileWriter"   , "emcalHistocomp", "-datafile emcalHistograms -concatenate-events -overwrite");
    }

    if (doPhos)
    {
        // PHOS Histograms
        AliHLTConfiguration hconf("phosHistocomp", "CaloPhysicsHistos", "track-matcher", "-phos -invariantmass -clusterenergy -matchedtracks");
        AliHLTConfiguration fwconf("phosHist", "ROOTFileWriter"   , "phosHistocomp", "-datafile phosHistograms -concatenate-events -overwrite");
    }
    
    TString arg, ec;


    ec.Form("ESD-CONVERTER");
    arg = "";

    AliHLTConfiguration esdcconf(ec.Data(), "GlobalEsdConverter"   , "track-matcher", "");

    AliReconstruction rec;

    rec.SetInput(input);
    rec.SetRunVertexFinder(kFALSE);
    rec.SetRunReconstruction("HLT");
    rec.SetRunTracking(":");
    rec.SetLoadAlignFromCDB(0);
    rec.SetRunQA(":");
    //  rec.SetRunLocalReconstruction("PHOS") ;
    rec.SetOption("HLT", option);
    rec.SetSpecificStorage("GRP/GRP/Data", Form("local://%s", grp));
    rec.Run();
}
