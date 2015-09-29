void mkOCDBEmean(const char* fin="calib.PHOSCailb.root", Double_t eLow=0.06,
            Int_t mod=-1, Bool_t writeOCDB=kFALSE, Bool_t multiplyCC=kTRUE,
            const char* list="clist")
{
    //Calculates corrections ~1 to calibration coefficients (CC).
    //If mod<0, corrections calculated for all modules.
    //If writeOCDB=kTRUE, writes corrections to the local OCDB.
    //If multiplyCC=kTRUE, multiplies CCs by the CCs from the previous iteration.
    
    //Module numeration: 3,2,1.
    //Online:            2,3,4.
    
    TFile f(fin);
    
    TH1F* meanAmp[5] = {};
    TH1F* h1 = 0;
    Double_t meanPeak[5] = {};
    Double_t cc_prev;
    Double_t cc_i;
    
    char hname[128];
    char mname[128], mtitle[128];
    
    TH1F* cHist[5];
    char cname[128], ctitle[128];
    
    AliCDBManager::Instance()->SetDefaultStorage("local://OCDB");
    // const char* path = "alien://Folder=/alice/cern.ch/user/p/polishch/OCDB_mean_amps";
    // AliCDBManager::Instance()->SetDefaultStorage(path);
    
    AliPHOSCalibData cdb(0);
    // cdb.CreateNew();
    
    for(Int_t i=0; i<5; i++) {
        sprintf(cname,"ccHist%d",i+1);
        sprintf(ctitle,"Multiplicative corrections to the CC of Module %d", i+1);
        cHist[i] = new TH1F(cname,ctitle,1000,0.,2.);
    }
    
    TList *CellAmplitudes = (TList*)f.Get(list);
    if(!CellAmplitudes) { printf("TList %s not found!\n",list); }
    
    for(Int_t imod=1;imod<=5;imod++) {
        for(Int_t iX=1;iX<=64;iX++) {
            for(Int_t iZ=1;iZ<=56;iZ++) {
                sprintf(hname,"cell_m%d_x%d_z%d",imod,iX,iZ);
                if(CellAmplitudes)
                    h1 = (TH1F*)CellAmplitudes->FindObject(hname);
                else
                    h1 = (TH1F*)f.Get(hname);
                if(h1) {
                    if(eLow>0) h1->GetXaxis()->SetRangeUser(eLow,4.0);
                    if(!meanAmp[imod-1]) {
                        printf("Creating histogram for Module %d..\n",imod);
                        sprintf(mname,"meanAmpMod%d",imod);
                        sprintf(mtitle,"Mean Amplitudes in the cells of Module %d",imod);
                        meanAmp[imod-1] = new TH1F(mname,mtitle,500,0.,1.);
                    }
                    if(h1->GetMean()>0.) {
                        meanAmp[imod-1]->Fill(h1->GetMean());
                        meanPeak[imod-1] = meanAmp[imod-1]->GetBinCenter( meanAmp[imod-1]->GetMaximumBin());
                    }
                    
                }
                
            }
        }
    }
    
    
    for(Int_t imod=1;imod<=5;imod++) {
        printf("MeanPeak for Module [%d] =  %.3f\n",imod,meanPeak[imod-1]);
    }
    
    //Calculate CC
    
    for(Int_t imod=1;imod<=5;imod++) {
        for(Int_t iX=1;iX<=64;iX++) {
            for(Int_t iZ=1;iZ<=56;iZ++) {
                sprintf(hname,"cell_m%d_x%d_z%d",imod,iX,iZ);
                if(CellAmplitudes)
                    h1 = (TH1F*)CellAmplitudes->FindObject(hname);
                else
                    h1 = (TH1F*)f.Get(hname);
                if(h1) {
                    if(eLow>0) h1->GetXaxis()->SetRangeUser(eLow,4.0);
                    
                    if(h1->GetMean()>0.) {
                        if(multiplyCC)
                            cc_prev = cdb.GetADCchannelEmc(imod,iZ,iX);
                        else
                            cc_prev = 1.;
                        // cc_i = cc_prev*0.135/0.115*0.151/h1->GetMean();
                        cc_i = cc_prev*meanPeak[imod-1]/h1->GetMean();
                    }
                    else
                        cc_i = 1.;
                }
                else
                    cc_i = 1.;
                
                cdb.SetADCchannelEmc(imod,iZ,iX,cc_i);
                if(mod>0 && imod != mod) cdb.SetADCchannelEmc(imod,iZ,iX,1.);
                
                if(cc_i != 1.)
                    cHist[imod-1]->Fill(cc_i);
            }
        }
    }
    
    char outf[128];
    
    if(eLow>0)
        sprintf(outf,"mkOCDB_eLow=%.2f.root",eLow);
    else
        sprintf(outf,"mkOCDB.root");
    
    TFile fout(outf,"recreate");
    
    for(Int_t i=0; i<5; i++) {
        if(meanAmp[i]) meanAmp[i]->Write();
        cHist[i]->Write();
    }
    
    //Writing OCDB
    AliCDBMetaData *md= new AliCDBMetaData();
    md->SetResponsible("Boris Polishchuk");
    
    if(writeOCDB) {
        printf("Writing OCDB..");
        cdb.WriteEmc(0,999999,md);
        // cdb.WriteCpv(0,999999,md);
        // cdb.WriteEmcBadChannelsMap(0,999999,md);
        printf("done.\n");
    }
    
}
