
/* 
 * extractJetFlow.C macro
 * 
 * author: Redmer Alexander Bertens, 2013
 * rbertens@cern.ch, r.a.bertens@uu.nl, rbertens@nikhef.nl
 *
 * this macro performs a number of calculations which are necessary to understand the output
 * of the PWGJE PbPb Jet Flow wagons, most notably
 *
 * [1] get delta pt widths for different background subtraction schemes
 * [2] extract hybrid track (jet constituent track) and jet flow
 *     from different flow methods
 */

    // global variables
    TFile f("AnalysisResults.root");
    TFile w("extractedJetFlow.root", "RECREATE");
    TList* keys = f.GetListOfKeys();
    TList* lNoFit(0x0);         // placeholder pointer for output list no fit
    TList* lComb(0x0);          // placeholder pointer for output list combined fit
    TList* lInt(0x0);           // placeholder pointer for output list integrated v2
    TDirectoryFile* dirNoFit(0x0);      // placeholder dir for JetFlow
    TDirectoryFile* dirComb(0x0);        
    TDirectoryFile* dirInt(0x0);        
    // centralities
    Double_t _c[] = {0, 10, 30, 50, 90};
    TArrayD* centralities = new TArrayD((int)(sizeof(_c)/sizeof(_c[0])), _c);
    static const int maxCen((int)(sizeof(_c)/sizeof(_c[0])));// max number of centrality bins
    // pt array for jet flow analysis
    Double_t ptJ[] = {1, 10, 20, 30, 50, 100, 150, 200};
    TArrayD* _ptJ = new TArrayD(sizeof(ptJ)/sizeof(ptJ[0]), ptJ);
    // pt array for hybrid flow analysis
    Double_t ptH[] = {0., 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5};
    TArrayD* _ptH = new TArrayD(sizeof(ptH)/sizeof(ptH[0]), ptH);
    Double_t jetRadius(0.3);
    AliAnalysisTaskRhoVnModulation* rho = new AliAnalysisTaskRhoVnModulation();
    TH1F* r2V0A(0x0);      // container for second order ep resolution
    TH1F* r2V0C(0x0);   
    TH1F* r3V0A(0x0);      // container for third order ep resolution
    TH1F* r3V0C(0x0);
    // placeholders for delta pt histos from LHS fit
    TH1F* dPtTheoryNoVn(0x0);
    TH1F* dPtTheoryVn(0x0);
    TH1F* dPtkNoFit(0x0);
    TH1F* dPtkComb(0x0);
    TH1F* dPtkInt(0x0);
    // placeholders for delta pt histos from RMS
    TH1F* RMSdPtTheoryNoVn(0x0);
    TH1F* RMSdPtTheoryVn(0x0);
    TH1F* RMSdPtkNoFit(0x0);
    TH1F* RMSdPtkComb(0x0);
    TH1F* RMSdPtkInt(0x0);

//_____________________________________________________________________________
void extractJetFlow() {
    // macro to read output of v1.0 of jet flow tasks
    ReadOutputFile();
    w.cd();
    // Get the delta pt info from reading RMS and Mean from the histograms 
    if(lNoFit)  {
        w.mkdir(Form("DeltaPt_HISTO_%s", lNoFit->GetName()));
        w.cd(Form("DeltaPt_HISTO_%s", lNoFit->GetName()));
        RMSdPtkNoFit = GetDeltaPtRMS(lNoFit);
        RMSdPtkNoFit->Write();
    }
    if(lComb)   {
        w.mkdir(Form("DeltaPt_HISTO_%s", lComb->GetName()));
        w.cd(Form("DeltaPt_HISTO_%s", lComb->GetName()));
        RMSdPtkComb = GetDeltaPtRMS(lComb);
        RMSdPtkComb->Write();
    }
    if(lInt)    {
        w.mkdir(Form("DeltaPt_HISTO_%s", lInt->GetName()));
        w.cd(Form("DeltaPt_HISTO_%s", lInt->GetName()));
        RMSdPtkInt = GetDeltaPtRMS(lInt);
        RMSdPtkInt->Write();
    }
    // Get the delta pt info from doing a iterative LHS gaus fit
    if(lNoFit) {
        w.mkdir(Form("DeltaPt_LHSFIT_%s", lNoFit->GetName()));
        w.cd(Form("DeltaPt_LHSFIT_%s", lNoFit->GetName()));
        dPtkNoFit = GetDeltaPtSigma(lNoFit);
        dPtkNoFit->Write();
        GetDeltaPtMean(lNoFit)->Write();
    }
    if(lComb) {
        w.mkdir(Form("DeltaPt_LHSFIT_%s", lComb->GetName()));
        w.cd(Form("DeltaPt_LHSFIT_%s", lComb->GetName()));
        dPtkComb = GetDeltaPtSigma(lComb);
        dPtkComb->Write();
        GetDeltaPtMean(lComb)->Write();
    }
    if(lInt) {
        w.mkdir(Form("DeltaPt_LHSFIT_%s", lInt->GetName()));
        w.cd(Form("DeltaPt_LHSFIT_%s", lInt->GetName()));
        dPtkInt = GetDeltaPtSigma(lInt);
        dPtkInt->Write();
        GetDeltaPtMean(lInt)->Write();
    }
    // Get the delta pt predictions
    w.mkdir("DeltaPt_PREDICTION");
    GetPredictedDeltaPtSigma(lComb);
    // extract the flow
    TList* listNoFit = dirNoFit->GetListOfKeys();
    for(Int_t i(0); i < listNoFit->GetEntries(); i++) {
        TString string = listNoFit->At(i)->GetName();
        if(string.Contains("JetFlow")) {
            for(Int_t j(0); j < 5; j++) {
                if(string.EndsWith(Form("_%i_histograms", j*10))) {
                    TList* op = (TList*)dirNoFit->Get(string);
                    GetJetTrackFlow(op, j*10-5);
                }
            }
        }
    }
    // get the relative improvements
    GetRelativeImprovements();
    GetRelativeImprovementsFromRMS();
    // get flow of hyrbid tracks (jet constituents) and flow of jets
    TList* listNoFit = dirNoFit->GetListOfKeys();
    for(Int_t i(0); i < listNoFit->GetEntries(); i++) {
        TString string = listNoFit->At(i)->GetName();
        if(string.Contains("HybridFlow")) {
            for(Int_t j(0); j < 5; j++) {
                if(string.EndsWith(Form("_%i_histograms", j*10))) {
                    TList* op = (TList*)dirNoFit->Get(string);
                    GetHybridTrackFlow(op, j*10-5);
                }
            }
        }
    }
    TList* listComb = dirComb->GetListOfKeys();
    for(Int_t i(0); i < listComb->GetEntries(); i++) {
        TString string = listComb->At(i)->GetName();
        if(string.Contains("JetFlow")) {
            for(Int_t j(0); j < 5; j++) {
                if(string.EndsWith(Form("_%i_histograms", j*10))) {
                    TList* op = (TList*)dirComb->Get(string);
                    GetJetTrackFlow(op, j*10-5);
                }
            }
        }
    }
    TList* listInt = dirInt->GetListOfKeys();
    for(Int_t i(0); i < listInt->GetEntries(); i++) {
        TString string = listInt->At(i)->GetName();
        if(string.Contains("JetFlow")) {
            for(Int_t j(0); j < 5; j++) {
                if(string.EndsWith(Form("_%i_histograms", j*10))) {
                    TList* op = (TList*)dirInt->Get(string);
                    GetJetTrackFlow(op, j*10-5);
                }
            }
        }
    }
    // lock and write the output file
    w.Close();
}
//_____________________________________________________________________________
void LoadLibraries() {
    // Load common libraries - who knows what the future will bring ...
    // note that these need to be loaded prior to executing the macro (they're
    // necessary for the types specified as global variables)
    gSystem->Load("libTree");
    gSystem->Load("libVMC");
    gSystem->Load("libGeom");
    gSystem->Load("libGui");
    gSystem->Load("libXMLParser");
    gSystem->Load("libMinuit");
    gSystem->Load("libMinuit2");
    gSystem->Load("libProof");
    gSystem->Load("libPhysics");
    gSystem->Load("libSTEERBase");
    gSystem->Load("libESD");
    gSystem->Load("libAOD");
    gSystem->Load("libOADB");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libCDB");
    gSystem->Load("libRAWDatabase");
  //  gSystem->Load("libSTEER");
    gSystem->Load("libEVGEN");
    gSystem->Load("libANALYSISalice");
    gSystem->Load("libCORRFW");
    gSystem->Load("libTOFbase");
    //gSystem->Load("libTOFrec");
    gSystem->Load("libRAWDatabase.so");
    gSystem->Load("libRAWDatarec.so");
    gSystem->Load("libTPCbase.so");
    gSystem->Load("libTPCrec.so");
    gSystem->Load("libITSbase.so");
    gSystem->Load("libITSrec.so");
    gSystem->Load("libTRDbase.so");
    gSystem->Load("libTENDER.so");
    gSystem->Load("libSTAT.so");
    gSystem->Load("libTRDrec.so");
    gSystem->Load("libHMPIDbase.so");
    gSystem->Load("libPWGPP.so");
    gSystem->Load("libPWGHFbase");
    gSystem->Load("libPWGDQdielectron");
    gSystem->Load("libPWGHFhfe");
    gSystem->Load("libEMCALUtils");
    gSystem->Load("libPHOSUtils");
    gSystem->Load("libPWGCaloTrackCorrBase");
    gSystem->Load("libEMCALraw");
    gSystem->Load("libEMCALbase");
    gSystem->Load("libEMCALrec");
    gSystem->Load("libTRDbase");
    gSystem->Load("libVZERObase");
    gSystem->Load("libVZEROrec");
    gSystem->Load("libTENDER");
    gSystem->Load("libTENDERSupplies");
    gSystem->Load("libPWGEMCAL");
    gSystem->Load("libPWGGAEMCALTasks");
    gSystem->Load("libPWGTools");
    gSystem->Load("libPWGCFCorrelationsBase");
    gSystem->Load("libPWGCFCorrelationsDPhi");
    gSystem->Load("libCGAL");
    gSystem->Load("libfastjet");
    gSystem->Load("libSISConePlugin");
    gSystem->Load("libJETAN");
    gSystem->Load("libFASTJETAN");
    gSystem->Load("libPWGJEEMCALJetTasks");
    gSystem->Load("libPWGflowBase");
}
//_____________________________________________________________________________
void ReadOutputFile() {
    // loop over the keys
    for(Int_t i(0); i < keys->GetEntries(); i++) {
        if(TString(keys->At(i)->GetName()).EndsWith("PWGJE")) {
            if (TString(keys->At(i)->GetName()).Contains("kComb")) {
                printf(" > Found PWGJE output \n \t %s \n ", keys->At(i)->GetName());
                lComb = (TList*)f.Get(keys->At(i)->GetName());
            }
            if (TString(keys->At(i)->GetName()).Contains("kInt")) {
                printf(" > Found PWGJE output \n \t %s \n ", keys->At(i)->GetName());
                lInt = (TList*)f.Get(keys->At(i)->GetName());
            }
            if (TString(keys->At(i)->GetName()).Contains("kNoFit")) {
                printf(" > Found PWGJE output \n \t %s \n ", keys->At(i)->GetName());
                lNoFit = (TList*)f.Get(keys->At(i)->GetName());
            }
        }
        if(TString(keys->At(i)->GetName()).EndsWith("PWGCF")) {
            if(TString(keys->At(i)->GetName()).Contains("kComb")) {
                printf(" > Found PWGCF output \n \t %s \n ", keys->At(i)->GetName());
                dirComb = (TDirectoryFile*)f.Get(keys->At(i)->GetName());
            }
            if(TString(keys->At(i)->GetName()).Contains("kInt")) {
                printf(" > Found PWGCF output \n \t %s \n ", keys->At(i)->GetName());
                dirInt = (TDirectoryFile*)f.Get(keys->At(i)->GetName());
            }
            if(TString(keys->At(i)->GetName()).Contains("kNoFit")) {
                printf(" > Found PWGCF output \n \t %s \n ", keys->At(i)->GetName());
                dirNoFit = (TDirectoryFile*)f.Get(keys->At(i)->GetName());
            }
        }
    }   // end of loop over keys
}
//_____________________________________________________________________________
GetRelativeImprovements() {
    // get the relative improvements in delta pt widths
    // note that the error propagation towards the relative
    // improvement is NOT CORRECT !
    if(dPtTheoryVn && dPtTheoryNoVn ) {
        TH1F* impTheory = new TH1F("relative improvement #delta p_{T} #sigma, theory ", "relative improvement #delta p_{T} #sigma, theory", centralities->GetSize()-1, centralities->GetArray());
        impTheory->GetYaxis()->SetTitle("relative improvement [ #frac{#delta p_{T} #sigma no v_{n} - #delta p_{T} #sigma v_{n}}{#delta p_{T} #sigma no v_{n}} ]");
        impTheory->GetXaxis()->SetTitle("centrality percentile");
        for(Int_t i (0); i < centralities->GetSize(); i++) {
            Double_t a = dPtTheoryNoVn->GetBinContent(i+1);
            Double_t b = dPtTheoryVn->GetBinContent(i+1);
            if(b!=0) {
                impTheory->SetBinContent(i+1, (b-a)/b);
                impTheory->SetBinError(1+i, impTheory->GetBinError(i+1));
            }
        }
    }
    if(dPtkNoFit && dPtkComb) {
    TH1F* impComb = new TH1F("relative improvement #delta p_{T} #sigma, kCombined ", "relative improvement #delta p_{T} #sigma, kCombined", centralities->GetSize()-1, centralities->GetArray());
        impComb->GetYaxis()->SetTitle("relative improvement [ #frac{#delta p_{T} #sigma no v_{n} - #delta p_{T} #sigma v_{n}}{#delta p_{T} #sigma no v_{n}} ]");
        impComb->GetXaxis()->SetTitle("centrality percentile");
        for(Int_t i (0); i < centralities->GetSize(); i++) {
            Double_t a = dPtkComb->GetBinContent(i+1);
            Double_t b = dPtkNoFit->GetBinContent(i+1);
            if(b!=0) {
                impComb->SetBinContent(i+1, (b-a)/b);
                impComb->SetBinError(i+1, impComb->GetBinError(i+1));
            }
        }
    }
    if(dPtkNoFit && dPtkInt) {
    TH1F* impInt = new TH1F("relative improvement #delta p_{T} #sigma, kInt ", "relative improvement #delta p_{T} #sigma, kInt", centralities->GetSize()-1, centralities->GetArray());
        impInt->GetYaxis()->SetTitle("relative improvement [ #frac{#delta p_{T} #sigma no v_{n} - #delta p_{T} #sigma v_{n}}{#delta p_{T} #sigma no v_{n}} ]");
        impInt->GetXaxis()->SetTitle("centrality percentile");
        for(Int_t i (0); i < centralities->GetSize(); i++) {
            Double_t a = dPtkInt->GetBinContent(i+1);
            Double_t b = dPtkNoFit->GetBinContent(i+1);
            if(b!=0) {
                impInt->SetBinContent(i+1, (b-a)/b);
                impInt->SetBinError(1+i, impInt->GetBinError(i+1));
            }
        }
    }
    // write the output to file
    w.cd();
    w.mkdir("Relative improvement delta pt distributions");
    w.cd("Relative improvement delta pt distributions");
    impTheory->Write();
    impComb->Write();
    impInt->Write();
    
}
//_____________________________________________________________________________
GetRelativeImprovementsFromRMS() {
    // get the relative improvements in delta pt widths
    // note that the error propagation towards the relative
    // improvement is NOT CORRECT !
    if(dPtTheoryVn && dPtTheoryNoVn ) {
        TH1F* impTheory = new TH1F("relative improvement #delta p_{T} #sigma, theory", "relative improvement #delta p_{T} #sigma, theory ", centralities->GetSize()-1, centralities->GetArray());
        impTheory->GetYaxis()->SetTitle("relative improvement [ #frac{#delta p_{T} #sigma no v_{n} - #delta p_{T} #sigma v_{n}}{#delta p_{T} #sigma no v_{n}} ]");
        impTheory->GetXaxis()->SetTitle("centrality percentile");
        for(Int_t i (0); i < centralities->GetSize(); i++) {
            Double_t a = RMSdPtTheoryNoVn->GetBinContent(i+1);
            Double_t b = RMSdPtTheoryVn->GetBinContent(i+1);
            if(b!=0) {
                impTheory->SetBinContent(i+1, (b-a)/b);
                impTheory->SetBinError(1+i, impTheory->GetBinError(i+1));
            }
        }
    }
    if(RMSdPtkNoFit && RMSdPtkComb) {
    TH1F* impComb = new TH1F("relative improvement #delta p_{T} #sigma, kCombined RMS", "relative improvement #delta p_{T} #sigma, kCombined RMS", centralities->GetSize()-1, centralities->GetArray());
        impComb->GetYaxis()->SetTitle("relative improvement [ #frac{#delta p_{T} #sigma no v_{n} - #delta p_{T} #sigma v_{n}}{#delta p_{T} #sigma no v_{n}} ]");
        impComb->GetXaxis()->SetTitle("centrality percentile");
        for(Int_t i (0); i < centralities->GetSize(); i++) {
            Double_t a = RMSdPtkComb->GetBinContent(i+1);
            Double_t b = RMSdPtkNoFit->GetBinContent(i+1);
            if(b!=0) {
                impComb->SetBinContent(i+1, (b-a)/b);
                impComb->SetBinError(i+1, impComb->GetBinError(i+1));
            }
        }
    }
    if(RMSdPtkNoFit && RMSdPtkInt) {
    TH1F* impInt = new TH1F("relative improvement #delta p_{T} #sigma, kInt RMS", "relative improvement #delta p_{T} #sigma, kInt RMS", centralities->GetSize()-1, centralities->GetArray());
        impInt->GetYaxis()->SetTitle("relative improvement [ #frac{#delta p_{T} #sigma no v_{n} - #delta p_{T} #sigma v_{n}}{#delta p_{T} #sigma no v_{n}} ]");
        impInt->GetXaxis()->SetTitle("centrality percentile");
        for(Int_t i (0); i < centralities->GetSize(); i++) {
            Double_t a = RMSdPtkInt->GetBinContent(i+1);
            Double_t b = RMSdPtkNoFit->GetBinContent(i+1);
            if(b!=0) {
                impInt->SetBinContent(i+1, (b-a)/b);
                impInt->SetBinError(1+i, impInt->GetBinError(i+1));
            }
        }
    }
    // write the output to file
    w.cd();
    w.mkdir("Relative improvement delta pt distributions from RMS");
    w.cd("Relative improvement delta pt distributions from RMS");
    impTheory->Write();
    impComb->Write();
    impInt->Write();
    
}
//_____________________________________________________________________________
void GetHybridTrackFlow(TList* jf, Int_t c) {
    // get hybrid track flow from the vzero ep and qc anlaysis
    if(!jf) {
        printf(" > couldn't find output list with name %s \n", name.Data());
        return;
    }
    w.mkdir(Form("GetHybridTrackFlow_%s", jf->GetName()));
    w.cd(Form("GetHybridTrackFlow_%s", jf->GetName()));
    TProfile* v0a = (TProfile*)jf->FindObject("Differential v_{2}^{obs} VZEROA");
    TProfile* v0c = (TProfile*)jf->FindObject("Differential v_{2}^{obs} VZEROC");
    TProfile* qc2 = (TProfile*)jf->FindObject("Differential cumulants v2");
    TProfile* rc  = (TProfile*)jf->FindObject("Reference cumulants");
    if(qc2 && rc && v0a ) {
        TH1F* result = rho->GetDifferentialQC(rc, qc2, _ptH, 2);
        TString t = "qc2_";
        t+=jf->GetName();
        result->SetNameTitle(t.Data(), t.Data());
        for(Int_t i(0); i < result->GetXaxis()->GetNbins(); i++) {
//            Double_t M(rc->GetBinEntries(1)/2./700.);
//            Double_t Mp(qc2->GetBinEntries(1+i)/qc2->GetEntries());
//            Double_t errinv(rc->GetBinContent(1)*TMath::Sqrt(M*Mp));
//            cout << " errinv " << errinv << endl;
//            if(errinv > 0) errinv = TMath::Sqrt(errinv);
//            cout << " err " << errinv << endl;
            result->SetBinError(1+i, v0a->GetBinError(1+i));
        }
        result->Write();
    }
    if(v0a) {
        TH1F* result = rho->CorrectForResolutionDiff((TH1F*)v0a, AliAnalysisTaskRhoVnModulation::detectorType::kVZEROA, centralities, c, 2);
        result->Write();
    }
    if(v0c) {
        TH1F* result = rho->CorrectForResolutionDiff((TH1F*)v0c, AliAnalysisTaskRhoVnModulation::detectorType::kVZEROC, centralities, c, 2);
        result->Write();
    }
    // attempt to get the flow from the qc analysis
    TDirectoryFile* qc = (TDirectoryFile*)f.Get("QC");
    if(qc) {
        TList* qcl = (TList*)qc->Get(Form("QC_hybrid_flow_%s", app.Data()));
        if(qcl) {
            AliFlowCommonHistResults* flow = (AliFlowCommonHistResults*)qcl->FindObject("AliFlowCommonHistResults2ndOrderQC");
            if(flow) flow->GetHistDiffFlowPtPOI()->Write();
        }
    }
}
//_____________________________________________________________________________
void GetJetTrackFlow(TList* jf, Int_t c) {
    // get hybrid track flow from the vzero ep and qc anlaysis
    if(!jf) {
        printf(" > couldn't find output list with name %s \n", name.Data());
        return;
    }
    w.mkdir(Form("GetJetTrackFlow_%s", jf->GetName()));
    w.cd(Form("GetJetTrackFlow_%s", jf->GetName()));
    TProfile* v0a = (TProfile*)jf->FindObject("Differential v_{2}^{obs} VZEROA");
    TProfile* v0c = (TProfile*)jf->FindObject("Differential v_{2}^{obs} VZEROC");
    TProfile* qc2 = (TProfile*)jf->FindObject("Differential cumulants v2");
    TProfile* rc  = (TProfile*)jf->FindObject("Reference cumulants");
        if(qc2 && rc ) {
        TH1F* result = rho->GetDifferentialQC(rc, qc2, _ptJ, 2);
        TString t = "qc2_";
        t+=jf->GetName();
        result->SetNameTitle(t.Data(), t.Data());
        result->Write();
    }
    if(v0a) {
        TH1F* result = rho->CorrectForResolutionDiff((TH1F*)v0a, AliAnalysisTaskRhoVnModulation::detectorType::kVZEROA, centralities, c, 2);
        result->Write();
    }
    if(v0c) {
        TH1F* result = rho->CorrectForResolutionDiff((TH1F*)v0c, AliAnalysisTaskRhoVnModulation::detectorType::kVZEROC, centralities, c, 2);
        result->Write();
    }
}
//_____________________________________________________________________________
TH1F* GetDeltaPtRMS(TList* l) {
    // get the RMS value of delta pt
    TH1F* deltaPtRMS = new TH1F("#delta p_{T} RMS", "#delta p_{T} RMS", centralities->GetSize()-1, centralities->GetArray());
    deltaPtRMS->GetXaxis()->SetTitle("centrality percentile");
    deltaPtRMS->GetYaxis()->SetTitle("RMS [GeV/c]");
    for(Int_t i(0); i < maxCen; i++) {
        TString name = Form("fHistDeltaPtDeltaPhi2_%i", i);
        printf(" > searching for %s in list \n %s < \n ", name.Data(), l->GetName());
        TH2D* dpt((TH2D*)l->FindObject(name.Data()));
        if(!dpt) continue;
        deltaPtRMS->SetBinContent(i+1, dpt->GetRMS(2));
        deltaPtRMS->SetBinError(i+1, dpt->GetRMSError(2));
    }
    return deltaPtRMS;
}
//_____________________________________________________________________________
TH1F* GetDeltaPtSigma(TList* l) {
    // get the sigma of the delta pt distribution from a recursive LHS gauss fit
    TH1F* deltaPtMean = new TH1F("#delta p_{T} #sigma", "#delta p_{T} #sigma", centralities->GetSize()-1, centralities->GetArray());
    deltaPtMean->GetYaxis()->SetTitle("#sigma [GeV/c]");
    deltaPtMean->GetXaxis()->SetTitle("centrality percentile");
    for(Int_t i(0); i < maxCen; i++) {
        TString name = Form("fHistDeltaPtDeltaPhi2_%i", i);
        printf(" > searching for %s in list \n \t  %s < \n ", name.Data(), l->GetName());
        TH2D* dpt((TH2D*)l->FindObject(name.Data()));
        if(!dpt) continue;
        TH1D* temp = dpt->ProjectionY(/*"_py", 0, -1, "e"*/);    // do error propagation for projection
        Double_t s = temp->GetRMS();
        Double_t m = temp->GetMean();
        TF1* fit = new TF1(Form("sigma_%s", temp->GetName()), "gaus", m-3*s, m+0.5*s);
        for(Int_t j(0); j < 10; j++) {
            Double_t _m(m), _s(s);
            temp->Fit(fit, "QILR");       
            fit->SetRange(m-3*s, m+0.5*s);
            m = fit->GetParameter(1);
            s = fit->GetParameter(2);

        }
        deltaPtMean->SetBinContent(1+i, fit->GetParameter(2));
        deltaPtMean->SetBinError(1+i, fit->GetParError(2));
    }
    return deltaPtMean;
}
//_____________________________________________________________________________
TH1F* GetDeltaPtMean(TList* l) {
    // get the mean of the delta pt distribution from a recursive LHS gauss fit
    TH1F* deltaPtMean = new TH1F("#delta p_{T} mean", "#delta p_{T} meam", centralities->GetSize()-1, centralities->GetArray());
    deltaPtMean->GetYaxis()->SetTitle("mean [GeV/c]");
    deltaPtMean->GetXaxis()->SetTitle("centrality percentile");
    for(Int_t i(0); i < maxCen; i++) {
        TString name = Form("fHistDeltaPtDeltaPhi2_%i", i);
        printf(" > searching for %s in list \n \t %s < \n ", name.Data(), l->GetName());
        TH2D* dpt((TH2D*)l->FindObject(name.Data()));
        if(!dpt) continue;
        TH1D* temp = dpt->ProjectionY(/*"_py", 0, -1, "e"*/);    // do error propagation for projection
        Double_t s = temp->GetRMS();
        Double_t m = temp->GetMean();
        TF1* fit = new TF1(Form("mean_%s", temp->GetName()), "gaus", m-3*s, m+0.5*s);
        TH1F* qam = new TH1F(Form("QA_M_mean_%s", temp->GetName()), Form("QA_M_mean_%s", temp->GetName()), 10, 0, 10);
        TH1F* qas = new TH1F(Form("QA_S_mean_%s", temp->GetName()), Form("QA_S_mean_%s", temp->GetName()), 10, 0, 10);
        fit->SetParLimits(2, s/2., s*2.);
        for(Int_t j(0); j < 10; j++) {
            Double_t _m(m), _s(s);
            temp->Fit(fit, "QILR");       
            fit->SetRange(m-3*s, m+0.5*s);
            m = fit->GetParameter(1);
            s = fit->GetParameter(2);
            if(!m == 0) qam->SetBinContent(j+1, _m/m);
            if(!s == 0) qas->SetBinContent(j+1, _s/s);
        }
        deltaPtMean->SetBinContent(1+i, fit->GetParameter(1));
        deltaPtMean->SetBinError(1+i, fit->GetParError(1));
        temp->Write();
        qas->Write();
        qam->Write();
    }
    return deltaPtMean;
}
//_____________________________________________________________________________
void GetPredictedDeltaPtSigma(TList* l) {
    // get predicted delta pt sigma
    TH1F* deltaPtSigma = new TH1F("predicted #delta p_{T} #sigma ", "predicted #delta p_{T} #sigma", centralities->GetSize()-1, centralities->GetArray());
    deltaPtSigma->GetYaxis()->SetTitle("predicted #sigma [GeV/c]");
    deltaPtSigma->GetXaxis()->SetTitle("centrality percentile");
    TH1F* deltaPtSigmaNoV = new TH1F("predicted #delta p_{T} #sigma no vn", "predicted #delta p_{T} #sigma no vn", centralities->GetSize()-1, centralities->GetArray());
    deltaPtSigmaNoV->GetYaxis()->SetTitle("predicted #sigma [GeV/c]");
    deltaPtSigmaNoV->GetXaxis()->SetTitle("centrality percentile");
    // get the info from the methods of the class
    rho->SetOutputList((TList*)l->Clone());
    // get the resolution for the desired detector
    r2V0A = rho->GetResolutionFromOuptutFile(AliAnalysisTaskRhoVnModulation::detectorType::kVZEROA, 2, centralities);
    r3V0A = rho->GetResolutionFromOuptutFile(AliAnalysisTaskRhoVnModulation::detectorType::kVZEROA, 3, centralities);
    r2V0C = rho->GetResolutionFromOuptutFile(AliAnalysisTaskRhoVnModulation::detectorType::kVZEROC, 2, centralities);
    r3V0C = rho->GetResolutionFromOuptutFile(AliAnalysisTaskRhoVnModulation::detectorType::kVZEROC, 3, centralities);
    // grab the v2 and v3 values and do a resolution correction
    TH1F* v2 = new TH1F("v2", "v2", centralities->GetSize()-1, centralities->GetArray());
    TH1F* v3 = new TH1F("v3", "v3", centralities->GetSize()-1, centralities->GetArray());
    TProfile* pv2 = (TProfile*)l->FindObject("fProfV2");
    TProfile* pv3 = (TProfile*)l->FindObject("fProfV3");
    for(Int_t i(0); i < maxCen-1; i++) {
        v2->SetBinContent(1+i, pv2->GetBinContent(1+i));
        v2->SetBinError(1+i, pv2->GetBinError(1+i));
        v3->SetBinContent(1+i, pv3->GetBinContent(1+i));
        v3->SetBinError(1+i, pv3->GetBinError(1+i));
    }
    TH1F* cv2 = new TH1F("v2int","v2int",10,0,100);
    Double_t c_v2[] = {0, 0.036416,0.064765,0.084340,0.096771,0.104257,0.105902,0.104897,0.104811,0.104811,0.104811};
    cv2->SetContent(c_v2);
    TH1F* cv3 = new TH1F("v3int","v3int", 10, 0, 100);
    Double_t c_v3[] = {0, 0.0236149685, 0.02875255, 0.03241459, 0.03507416, 0.03730817, 0.03889757,0.04285879,0.05030896, 0, 0};
    cv3->SetContent(c_v3);
    for(Int_t i(0); i < centralities->GetSize()-1; i++) {
        TH1F* h = (TH1F*)l->FindObject(Form("fHistPicoTrackMult_%i", i));
        TH1F* m = (TH1F*)l->FindObject(Form("fHistPicoTrackPt_%i", i));
        if(!h||!m) {
            printf(" > Panic, one of the pt histos not found ! < \n");
            continue;
        }
        Double_t na(h->GetMean()*jetRadius*jetRadius*TMath::Pi()/(2.*.9*TMath::TwoPi()));
        Double_t naErr(h->GetMeanError()*jetRadius*jetRadius*TMath::Pi()/(2.*.9*TMath::TwoPi()));
        Double_t totErr(TMath::Sqrt((m->GetRMS()*m->GetRMS()+m->GetMean()*m->GetMean())*(m->GetRMS()*m->GetRMS()+m->GetMean()*m->GetMean())*naErr*naErr+4.*na*na*m->GetRMS()*m->GetRMS()*m->GetRMSError()*m->GetRMSError()+4.*na*na*m->GetMean()*m->GetMean()*m->GetMeanError()*m->GetMeanError()));
        Double_t err(m->GetRMS());
        Double_t mean(m->GetMean());
        Double_t _v2(cv2->GetBinContent(1+i));
        Double_t _v3(cv3->GetBinContent(1+i));
        // no v2
        Double_t dptnovn = TMath::Sqrt(na*(err*err+mean*mean));
        Double_t dpt = TMath::Sqrt(na*err*err+(na+2.*na*na*(_v2*_v2+_v3*_v3))*mean*mean);
        deltaPtSigma->SetBinContent(1+i, dpt);
        deltaPtSigma->SetBinError(1+i, totErr);
        deltaPtSigmaNoV->SetBinContent(1+i, dptnovn);
        deltaPtSigmaNoV->SetBinError(1+i, totErr);       // estimate without vn
    }
    w.mkdir("RhoTaskVnEstimates");
    w.cd("RhoTaskVnEstimates");
    r2V0A->Write();
    r3V0A->Write();
    cv2->Write();
    cv3->Write();
    w.cd("DeltaPt_PREDICTION");
    dPtTheoryVn = deltaPtSigma;
    RMSdPtTheoryVn = deltaPtSigma;
    dPtTheoryVn->Write();
    dPtTheoryNoVn = deltaPtSigmaNoV;
    RMSdPtTheoryNoVn = deltaPtSigmaNoV;
    dPtTheoryNoVn->Write();
}
//_____________________________________________________________________________
