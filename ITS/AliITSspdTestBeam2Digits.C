void AliITSspdTestBeam2Digits(const Char_t *filename="run001dat"){
    // Macro to convert spd Test Beam data to ITS SPD Digits.
    // fCoord1 = Row
    // fCoord2 = Columb
    Int_t i;

    if(gClassTable->GetID("AliRun") <0){
        gRoot->ProcessLine(".x $(ALICE_ROOT)/macros/loadlibs.C");
    } else if(gAlice) {
        delete gAlice->GetRunLoader();
        delete gAlice;
        gAlice = 0;
    } // end if
    if(!gAlice) { // most create it separatly?
        gAlice = new AliRun("gAlice",
                            "The ALICE Off-line Reconstruction Framework");
    } // end if

    AliRunLoader *rl = AliRunLoader::Open("galice.root",
                           AliConfig::fgkDefaultEventFolderName,"new");
    gAlice->SetRunLoader(rl);
    rl->SetNumberOfEventsPerFile(1000);

    rl->SetEventFolderName();
    rl->MakeTree("E");
    AliITSvSPD02 *its = (AliITS*)gAlice->GetDetector("ITS");
    if(!its){
        its = new AliITSvSPD02("SPD testbeam Run001");
        gAlice->AddModule(its);
        its->Init();
        its->SetDefaults();
    } // end if
    rl->AddLoader(its);
    AliITSLoader *ldr = (AliITSLoader*) rl->GetLoader("ITSLoader");
    ldr->SetEventFolder(rl->GetEventFolder());

    //rl->WriteGeometry();

    AliITSspdTestBeam *spd = new AliITSspdTestBeam(filename);
    spd->Read();
    spd->Decode();
    spd->SetLoader(ldr);
    spd->SetITS(its);
    AliHeader *hdr;
    for(i=0;i<spd->GetNumberOfEvents();i++){
        spd->Digitize(i);
        rl->MakeHeader();
        hdr = rl->GetHeader();
        hdr->Reset(001,i,i);
    } // end for i
    rl->WriteHeader("OVERWRITE");
    rl->WriteRunLoader("OVERWRITE");
    rl->WriteAliRun("OVERWRITE");
    delete rl;
}
