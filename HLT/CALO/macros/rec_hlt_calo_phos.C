
void rec_hlt_calo_phos(const char *input = "./", const char *grp = "./")
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

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // define the analysis chain to be run
    //
    int moduleStart = 2;
    int moduleEnd = 4;
    int rcuStart = 0;
    int rcuEnd = 3;
    int rcusPerModule = 4;
    int ddlOffset = 1792;

    TString option="libAliHLTUtil.so libAliHLTRCU.so libAliHLTCalo.so libAliHLTPHOS.so libAliHLTGlobal.so loglevel=0x7f chains=ESD-CONVERTER";
    TString ecInput;
    TString emInput;

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

    emInput = ecInput;

    TString arg, ec, em, hp, ef;


    ec.Form("ESD-CONVERTER");
    arg = "";

    AliHLTConfiguration esdcconf(ec.Data(), "GlobalEsdConverter"   , ecInput.Data(), "");

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
