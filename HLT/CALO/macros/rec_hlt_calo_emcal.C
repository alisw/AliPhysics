
void rec_hlt_calo_emcal(const char *input = "./", const char *grp = "./")
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
    int moduleEnd = 2;
    int rcuStart = 0;
    int rcuEnd = 1;
    int rcusPerModule = 2;
    int ddlOffset = 4608;

    TString option="libAliHLTUtil.so libAliHLTRCU.so libAliHLTCalo.so libAliHLTEMCAL.so libAliHLTGlobal.so loglevel=0x7f chains=ESD-CONVERTER";
    TString ecInput;
    TString emInput;

    for (int module = moduleStart; module <= moduleEnd; module++)
    {
        TString clInput;

        for (int rcu = rcuStart; rcu <= rcuEnd; rcu++)
        {
            TString arg, publisher, ra, dm;
            // raw data publisher components
            publisher.Form("EMC-RP_%02d_%d", module, rcu);
            arg.Form("-minid %d -datatype 'DDL_RAW ' 'EMCA'  -dataspec 0x%x ", 1792 + module*(4) + rcu, 0x1 << (module*4 + rcu));
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

    emInput = ecInput;

    TString arg, ec, em, hp, ef;


    ec.Form("ESD-CONVERTER");
    arg = "";

    AliHLTConfiguration esdcconf(ec.Data(), "GlobalEsdConverter"   , ecInput.Data(), "");

    AliReconstruction rec;

    rec.SetInput(input);

    rec.SetRunVertexFinder(kFALSE);
    rec.SetRunReconstruction("HLT EMCAL");
    rec.SetRunTracking(":");
    rec.SetLoadAlignFromCDB(0);
    rec.SetRunQA(":");

    rec.SetOption("HLT", option);
    rec.SetSpecificStorage("GRP/GRP/Data", Form("local://%s", grp));
    rec.Run();
}
