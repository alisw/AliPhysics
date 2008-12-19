Bool_t AliITSCompareHitsRecPoints(Char_t *rfn="galice.root"){
    // This macro compares the average location of the hits associated
    // with a specific track and the best match RecPoint location.
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return:
    //   kTRUE if no errors were encountered, otherwise kFALSE.
    TProfile *pxg,*pyg,*pzg,*pxl,*pyl,*pzl;
    Double_t hg[3],hl[3],tof,rg[3],rl[3];
    Float_t rgf[3];

    // Dynamically link some shared libs
    if (gClassTable->GetID("AliRun") < 0) {
        gROOT->LoadMacro("loadlibs.C");
        loadlibs();
    } else {
        if(gAlice){
            delete AliRunLoader::GetRunLoader();
            delete gAlice;
            gAlice=0;
        } // end if gAlice
    }// end if gClassTable->GetID()

    gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSstandard.C");

    AliRunLoader *runl = AccessFile(rfn); // Set up to read in Data
    Int_t retval = runl->LoadHeader();
    if (retval){
        cerr<<"AliITSPrintRecPoints.C : LoadHeader returned error"<<endl;
        return kFALSE;
    } // end if retval

    AliITSLoader* ITSloader =  (AliITSLoader*) runl->GetLoader("ITSLoader");
    if(!ITSloader){
        cerr<<"AliITSPrintRecPoints.C :  ITS loader not found"<<endl;
        return kFALSE;
    } // end if !ITSloader

    if(!gGeoManager){
        gGeoManger = new TGeoManager();
        gGeoManager->Import("geometry.root","");
    } // end if !gGeoManager
    if(ITSloader->LoadHits("read")!=0){
        cerr<<"Error Loading Hits"<<endl;
        return kFALSE;
    }// end if ITSloader->LoadHits
    if(ITSloader->LoadRecPoints("read")!=0){
        cerr<<"Error Loading RecPoints"<<endl;
        return kFALSE;
    }// end if ITSloader->LoadRecPoints
    AliITS *ITS = 0;
    ITS = (AliITS*)(gAlice->GetDetector("ITS"));
    if(!ITS){
        cout << "Error: no ITS found. Aborting"<<endl;
        return kFALSE;
    } // end if !ITS
    if(!(ITS->GetDetTypeSim())){
        cout <<"No AliITSDetTypeSim object found in ITS"<<endl;
        return kFALSE;
    } // end if
    AliITSmodule     *m   = 0;
    AliITShit        *h   = 0;
    AliITSDetTypeSim *sim = ITS->GetDetTypeSim();
    AliITSgeom *gm=0;
    gm = ITS->GetITSgeom();
    if(!gm){
        cout <<"No AliITSgeom object found in ITS"<<endl;
        if(!gGeoManager){
            cout <<"No gGeoManger. Aborting"<<endl;
            return kFALSE;
        }else{
            ITS->UpdateInternalGeometry(gGeoManager);
        } // end if
    } // end if !AliITSgeom
    //
    Int_t nMods= gm->GetIndexMax(),nEvents=AliRunLoader::GetNumberOfEvents();
    Int_t mod=0,evnt=0,size=-1,irp=0,ih=0,trkindexOld=-1;
    Double_t xmod,nHitPerTrack;
    TTree            *rpt = 0;
    TClonesArray     *rpa = 0;
    AliITSRecPoint   *rp  = 0;
    AliITSDetTypeRec *rec = new AliITSDetTypeRec(ITSloader);
    rec->SetDefaults();
    // We are going to use the following
    // <hit_x[i]> == { Sum_j hit_x[i][j] }/N_hits per track
    // d_x[i] == <hit_x[i]>-recpoint_x[i]
    // <d_x>  == {sum_i d_x[i]}/N_RecPoints   (excluding noise and no merging)
    //         = {sum_i (<hit_x[i]>-recpoint_x[i])}/N_RecPoints
    //         = {sum_i ((sum_j hit_x[i][j])/N_hits_per_track -recpoinnt_x[i])}/N_Recpoints
    //         = sum_i sum_j {hit_x[i][j]/N_hits_per_track}/N_Recpoints -
    //           sum_i {recpoint_x[i]}/N_Recpoints
    pxg = new TProfile("XdiffGlobal","Mean displacement in gX direction",nMods,
                       -0.5,(Double_t)(nMods)+0.5," ");
    pxg->SetXTitle("module");
    pxg->SetYTitle("Global X [cm]");
    pyg = new TProfile("YdiffGlobal","Mean displacement in gY direction",nMods,
                       -0.5,(Double_t)(nMods)+0.5," ");
    pyg->SetXTitle("module");
    pyg->SetYTitle("Global Y [cm]");
    pzg = new TProfile("ZdiffGlobal","Mean displacement in gZ direction",nMods,
                       -0.5,(Double_t)(nMods)+0.5," ");
    pzg->SetXTitle("module");
    pzg->SetYTitle("Global Z [cm]");
    pxl = new TProfile("XdiffLocal","Mean displacement in lX direction",nMods,
                       -0.5,(Double_t)(nMods)+0.5," ");
    pxl->SetXTitle("module");
    pxl->SetYTitle("local X [cm]");
    pyl = new TProfile("YdiffLocal","Mean displacement in lY direction",nMods,
                       -0.5,(Double_t)(nMods)+0.5," ");
    pyl->SetXTitle("module");
    pyl->SetYTitle("local Y [cm]");
    pzl = new TProfile("ZdiffLocal","Mean displacement in lZ direction",nMods,
                       -0.5,(Double_t)(nMods)+0.5," ");
    pzl->SetXTitle("module");
    pzl->SetYTitle("local Z [cm]");
    //
    for(evnt=0;evnt<nEvents;evnt++){
        runl->GetEvent(evnt);
        ITS->InitModules(size,nMods);
        ITS->FillModules(evnt,0,-1," "," ");
        rec->SetTreeAddress();
        for(mod=0;mod<nMods;mod++){
            xmod = (Double_t) mod;
            m = ITS->GetModule(mod);
            rec->ResetRecPoints();
            rpt = ITSloader->TreeR();
            rpt->GetEvent(mod);
            rpa = rec->RecPoints();
            for(irp=0;irp<rpa->GetEntriesFast();irp++){
                rp = (AliITSRecPoint*)(rpa->At(irp));
                rp->GetGlobalXYZ(rgf);rg[0]=rgf[0];rg[1]=rgf[1];rg[2]=rgf[2];
                rl[0] = rp->GetDetLocalX();
                rl[1] = 0.0;
                rl[2] = rp->GetDetLocalZ();
                pxg->Fill(xmod,-rg[0],0.5);
                pyg->Fill(xmod,-rg[1],0.5);
                pzg->Fill(xmod,-rg[2],0.5);
                pxl->Fill(xmod,-rl[0],0.5);
                //pyl->Fill(xmod,-rl[1],0.5); // assumed to be zero always.
                pzl->Fill(xmod,-rl[2],0.5);
            }// rnf got itp
            trkindexOld=-1;
            nHitPerTrack = 0.0;
            for(ih=0;ih<m->GetNhits();ih++){ // We want the median hit location
                // for each track.
                h = m->GetHit(ih);
                if(m->GetHitTrackIndex(ih)!=trkindexOld){// Enterence location
                    trkindexOld = m->GetHitTrackIndex(ih);
                    nHitPerTrack = 1.0;
                    do{
                        nHitPerTrack += 1.0;
                    }while (m->GetHitTrackIndex(ih+nHitPerTrack-1)==trkindexOld);
                    h->GetPositionG0(hg[0],hg[1],hg[2],tof);
                    h->GetPositionL0(hl[0],hl[1],hl[2],tof);
                    pxg->Fill(xmod,hg[0],0.5/nHitPerTrack);
                    pyg->Fill(xmod,hg[1],0.5/nHitPerTrack);
                    pzg->Fill(xmod,hg[2],0.5/nHitPerTrack);
                    pxl->Fill(xmod,hl[0],0.5/nHitPerTrack);
                    pyl->Fill(xmod,hl[1],1.0/nHitPerTrack); // rl[1]=0 always.
                    pzl->Fill(xmod,hl[2],0.5/nHitPerTrack);
                } // end if
                if(nHitPerTrack==0.0) continue;
                h->GetPositionG(hg[0],hg[1],hg[2],tof);
                h->GetPositionL(hl[0],hl[1],hl[2],tof);
                pxg->Fill(xmod,hg[0],0.5/nHitPerTrack);
                pyg->Fill(xmod,hg[1],0.5/nHitPerTrack);
                pzg->Fill(xmod,hg[2],0.5/nHitPerTrack);
                pxl->Fill(xmod,hl[0],0.5/nHitPerTrack);
                pyl->Fill(xmod,hl[1],1.0/nHitPerTrack); // rl[1]=0 always.
                pzl->Fill(xmod,hl[2],0.5/nHitPerTrack);
            } // end for ih
        } // end for mod
        ITS->ClearModules();
    } // end for evnt
    //
    Int_t wh=800;
    TCanvas *c0 = new TCanvas("c0","Average displacements between hits and RecPoints"
                              ,1.5*wh,wh);
    c0->Divide(3,2);
    c0->cd(1);
    pxg->Draw();
    c0->cd(2);
    pyg->Draw();
    c0->cd(3);
    pzg->Draw();
    c0->cd(4);
    pxl->Draw();
    c0->cd(5);
    pyl->Draw();
    c0->cd(6);
    pzl->Draw();
    c0->Update();
    return kTRUE;
}
