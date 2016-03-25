

// Macro to calcuate efficiency QA
// Jitendra Kumar, Zaida

//Range of container variable
const Double_t etamin = -0.9;
const Double_t etamax =  0.9;
const Double_t ptmin = 0.5;
const Double_t ptmax = 24.0;
const Double_t phimin = -2*TMath::Pi();
const Double_t phimax = 2*TMath::Pi();
const Double_t thetamin = 0;
const Double_t thetamax = TMath::Pi();
const Double_t zvtxmin = -10.0;
const Double_t zvtxmax =  10.0;
const Double_t multmin = -1.0;
const Double_t multmax = 102;
const Double_t centmin = -1.0;
const Double_t centmax = 100;


const UInt_t ipt = 0;
const UInt_t ieta  = 1;
const UInt_t iphi = 2;
const UInt_t itheta = 3;
const UInt_t izvtx = 4;
const UInt_t imult = 5;
const UInt_t icent = 6;

Int_t nvars = 7;
Int_t* ivarSlice = new Int_t[nvars];
ivarSlice[0]= ipt;
ivarSlice[1] = ieta;
ivarSlice[2] = iphi;
ivarSlice[3] = itheta;
ivarSlice[4]= izvtx;
ivarSlice[5] = imult;
ivarSlice[6] = icent;

Double_t *mins = new Double_t[nvars];
Double_t *maxs = new Double_t[nvars];
mins[ipt] = ptmin;        maxs[ipt] = ptmax;
mins[ieta] = etamin;      maxs[ieta] = etamax;
mins[iphi] = phimin;      maxs[iphi] = phimax;
mins[itheta] = thetamin;  maxs[itheta] = thetamax;
mins[izvtx] = zvtxmin;    maxs[izvtx] = zvtxmax;
mins[imult] = multmin;    maxs[imult] = multmax;
mins[icent] = centmin;    maxs[icent] = centmax;



//Slicing container
AliCFContainer *CalcSingleTrackEffQAC(const char *inputfile="Run_117222/AnalysisResults_trainMarch10.root", const char *fParticleType = "NchFbit0"){
    
    
    TFile* infile = TFile::Open(inputfile,"read");
    TDirectoryFile *dir = (TDirectoryFile*)infile->Get(Form("PWGPP_CFSingleTrack"));
    if(!dir){Printf("ERROR: Directory not available  --Exiting.. !"); return;}
    
    AliCFContainer *effContainer = (AliCFContainer*) (dir->Get(Form("container%s",fParticleType)));
    if(!effContainer){Printf("ERROR: CF Container not available  --Exiting.. !"); return;}
    
    Printf("\n_________________________________________");
    AliCFContainer *ConSlice = (AliCFContainer*)effContainer->MakeSlice(nvars,ivarSlice,mins,maxs);
    cout<< "  ... slice done"<<endl;
    cout<< "  the new container has "<< ConSlice->GetNStep() << " steps"<<endl;

    infile->Close();
    return ConSlice;
    
}


//Main functions
Float_t CalcSingleTrackEffQAV(AliCFContainer *ConSliceC,
                              const char * fWhichVar="pt",
                              Float_t fVarMin=0,
                              Float_t fVarMax=24,
                              Int_t fWhicheff=2, //StepA/StepB
                              Int_t fWhichfbit=0,
                              const char *fParticleType = "NchFbit4",
                              Bool_t isAvg = kFALSE,
                              Int_t runN = 999999)
{
    AliCFContainer *ConSliceClone = (AliCFContainer *)ConSliceC->Clone(Form("Slice%s%i",fParticleType,isAvg));
    
    //1. EFF NUMERATOR
    AliCFGridSparse* gridSparseNum = 0x0;
    if( fWhicheff == 1 )  gridSparseNum = (AliCFGridSparse*)ConSliceClone->GetGrid(1); // GenAcc
    else if( fWhicheff == 2 )  gridSparseNum = (AliCFGridSparse*)ConSliceClone->GetGrid(6); // Rec (acc + cuts) draw reco properties
    else if( fWhicheff == 3 )  gridSparseNum = (AliCFGridSparse*)ConSliceClone->GetGrid(7); // RecPID
    else if( fWhicheff == 4 )  gridSparseNum = (AliCFGridSparse*)ConSliceClone->GetGrid(3); // Rec (no cuts)
    else if( fWhicheff == 5 )  gridSparseNum = (AliCFGridSparse*)ConSliceClone->GetGrid(4); // RecAcc
    else if( fWhicheff == 6 )  gridSparseNum = (AliCFGridSparse*)ConSliceClone->GetGrid(5); // Rec (acc + cuts) draw kine properties
    
    //2 EFF DENOMINATOR
    AliCFGridSparse* gridSparseDen = 0x0;
    if( fWhicheff == 1 )gridSparseDen = (AliCFGridSparse*)ConSliceClone->GetGrid(0); // LimAcc
    else gridSparseDen = (AliCFGridSparse*)ConSliceClone->GetGrid(1);              // GenAcc
    
    if(isAvg)return CalcAvgSingleTrackEff(gridSparseNum, gridSparseDen, fWhichVar, fWhicheff, fParticleType, fWhichfbit, fVarMin, fVarMax);
    else if (!isAvg) return CalcTotalSingleTrackEff(gridSparseNum, gridSparseDen, fWhichVar, fWhicheff,fWhichfbit,fParticleType, runN);
    else Printf ("\n \n \n No efficiency calulcated .. \n --> exit");

}



//Function for full range efficiecy
Float_t CalcTotalSingleTrackEff(AliCFGridSparse* gridSprsNum,
                                AliCFGridSparse* gridSprsDen,
                                const char *fWhichVar="",
                                Int_t fWcheff=1,
                                Int_t fWchfbit=0,
                                const char *fParticleType="",
                                Int_t runNumber = 999999)
{
    
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(1);
    gStyle->SetPalette(1);
    gStyle->SetCanvasColor(0);
    gStyle->SetFrameFillColor(0);
    
    Int_t iVarContainer = -99;
    Double_t xvalue = 0.0;
    if(fWhichVar=="pt" || fWhichVar=="PT" || fWhichVar=="Pt"){iVarContainer=0; xvalue = 12.6;}
    else if(fWhichVar=="eta" || fWhichVar=="ETA" || fWhichVar=="Eta"){iVarContainer=1; xvalue = -0.1;}
    else if(fWhichVar=="phi" || fWhichVar=="PHI" || fWhichVar=="Phi"){iVarContainer=2;}
    else if(fWhichVar=="theta" ||fWhichVar=="THETA" || fWhichVar=="Theta"){iVarContainer=3;}
    else if(fWhichVar=="zvtx" || fWhichVar=="ZVTX" || fWhichVar=="Zvtx"){iVarContainer=4; xvalue = 1.0;}
    else if(fWhichVar=="mult" || fWhichVar=="MULT" || fWhichVar=="Mult"){iVarContainer=5;}
    else if(fWhichVar=="cent" || fWhichVar=="CENT" || fWhichVar=="Cent"){iVarContainer=6;}
    else {Printf("ERROR: Input Variable is not available  --Exiting..");return;}
    
    const char *cLevelNum = "--";
    const char *cLevelDen = "GenAcc";
    if( fWcheff == 1 ) { cLevelNum = "GenAcc"; cLevelDen = "LimAcc";}
    else if( fWcheff == 2 )cLevelNum = "RecTrkQuaCut"; // Rec (acc + cuts) draw reco properties
    else if( fWcheff == 3 )cLevelNum = "RecPID"; // RecPID
    else if( fWcheff == 4 )cLevelNum = "Reco"; // Rec (no cuts)
    else if( fWcheff == 5 )cLevelNum = "RecKineAcc"; // RecAcc
    else if( fWcheff == 6 )cLevelNum = "RecTrkMCPart"; // Rec (acc + cuts) draw kine properties
    
    TString hEffName = "";
    hEffName.Form("%s_%s_Efficiency_%sBy%s_fbit%d", fParticleType, fWhichVar, cLevelNum, cLevelDen, fWchfbit);

    AliCFGridSparse* gridSprsNumClone = (AliCFGridSparse*)gridSprsNum->Clone(Form("Num%s_%s_Efficiency_%sBy%s_fbit%d", fParticleType, fWhichVar, cLevelNum, cLevelDen, fWchfbit));
    if(!gridSprsNumClone) return;
    
    AliCFGridSparse* gridSprsDenClone = (AliCFGridSparse*)gridSprsDen->Clone(Form("Den%s_%s_Efficiency_%sBy%s_fbit%d", fParticleType, fWhichVar, cLevelNum, cLevelDen, fWchfbit));
    if(!gridSprsDenClone) return;
    
    cout << "Eff Numerator   ---> " << Form("Num%s_%s_Efficiency_%sBy%s_fbit%d", fParticleType, fWhichVar, cLevelNum, cLevelDen, fWchfbit) << endl;
    cout << "Eff Denominator ---> " << Form("Den%s_%s_Efficiency_%sBy%s_fbit%d", fParticleType, fWhichVar, cLevelNum, cLevelDen, fWchfbit) << endl;
    
    TCanvas *c = new TCanvas(hEffName.Data(), hEffName.Data(), 2250,750);
    c->Divide(3,1);
   
    THnSparse* numData = (THnSparse*)gridSprsNumClone->GetGrid();
    numData->SetTitle(Form("new1gridSprsDenClone%i%s",runNumber,fParticleType));
    TH1D* heffnum = (TH1D*)numData->Projection(iVarContainer);
    c->cd(2);
    if(fWhichVar=="pt" || fWhichVar=="PT" || fWhichVar=="Pt")gPad->SetLogy();
    Float_t fVarMin = heffnum->GetXaxis()->GetXmin();
    Float_t fVarMax = heffnum->GetXaxis()->GetXmax();
    hEffName.Resize(0);
    hEffName.Form("%s: %s_%s_distribution (%.1f-%.1f)",cLevelDen, fParticleType, fWhichVar, fVarMin, fVarMax);
    SetHisto(heffnum,hEffName,fWhichVar, 9);
    heffnum->Draw();
    

    THnSparse* denData = (THnSparse*)gridSprsDenClone->GetGrid();
    TH1D* heffden = (TH1D*)denData->Projection(iVarContainer);
    c->cd(1);
    if(fWhichVar=="pt" || fWhichVar=="PT" || fWhichVar=="Pt")gPad->SetLogy();
    hEffName.Resize(0);
    hEffName.Form("%s: %s_%s_distribution (%.1f-%.1f)",cLevelNum, fParticleType, fWhichVar, fVarMin, fVarMax);
    SetHisto(heffden,hEffName,fWhichVar, 6);
    heffden->Draw();
    
    
    //Efficiencys
    TH1D* heff = (TH1D*)heffnum->Clone("heff");
    hEffName.Resize(0);
    hEffName.Form("Efficiency: %s_%s (%.1f-%.1f)",fParticleType, fWhichVar, fVarMin, fVarMax);
    heff->Divide(heffnum,heffden,1,1,"B");
    heff->SetTitle(hEffName.Data());
    heff->GetYaxis()->SetTitle("Efficiency");
    heff->SetMaximum(1.00);
    heff->SetMinimum(-0.20);
    heff->SetMarkerColor(4); //blue
    c->cd(3);
 
    
    gPad->SetFrameFillColor(5);
    heff->Draw();
    LatexName(xvalue, fWhichVar, runNumber, fWchfbit, fParticleType);
    c->Update();
    hEffName.Resize(0);
    hEffName.Form("%s_%s_Efficiency_%sBy%s_fbit%d_%d.pdf", fParticleType, fWhichVar, cLevelNum, cLevelDen, fWchfbit, runNumber);
    c->SaveAs(hEffName.Data());
    return 0;
}


//Function for avg efficiecy
Float_t CalcAvgSingleTrackEff(AliCFGridSparse* gridSprsNum,
                              AliCFGridSparse* gridSprsDen,
                              const char *fWhichVar="",
                              Int_t fWcheff=2,
                              const char *fParticleType="",
                              Int_t fWchfbit=0,
                              Float_t fVarMin=4.0,
                              Float_t fVarMax=8.0)

{
    
    
    AliCFGridSparse* gridSprsNumClone = (AliCFGridSparse*)gridSprsNum->Clone(Form("nNumClone%0.1f%0.1f%s",fVarMin, fVarMax,fParticleType));
    if(!gridSprsNumClone) return;
    
    AliCFGridSparse* gridSprsDenClone = (AliCFGridSparse*)gridSprsDen->Clone(Form("nDenClone%0.1f%0.1f%s",fVarMin, fVarMax,fParticleType));
    if(!gridSprsDenClone) return;
    
    Int_t iVarContainer = -99;
    if(fWhichVar=="pt" || fWhichVar=="PT" || fWhichVar=="Pt")iVarContainer=0;
    else if(fWhichVar=="eta" || fWhichVar=="ETA" || fWhichVar=="Eta")iVarContainer=1;
    else if(fWhichVar=="phi" || fWhichVar=="PHI" || fWhichVar=="Phi")iVarContainer=2;
    else if(fWhichVar=="theta" ||fWhichVar=="THETA" || fWhichVar=="Theta")iVarContainer=3;
    else if(fWhichVar=="zvtx" || fWhichVar=="ZVTX" || fWhichVar=="Zvtx")iVarContainer=4;
    else if(fWhichVar=="mult" || fWhichVar=="MULT" || fWhichVar=="Mult")iVarContainer=5;
    else if(fWhichVar=="cent" || fWhichVar=="CENT" || fWhichVar=="Cent")iVarContainer=6;
    else {Printf("ERROR: Input Variable is not available  --Exiting..");return;}
    
    const char *cLevelNum = "--";
    const char *cLevelDen = "GenAcc";
    if( fWcheff == 1 ) { cLevelNum = "GenAcc"; cLevelDen = "LimAcc";}
    else if( fWcheff == 2 )cLevelNum = "RecTrkQuaCut"; // Rec (acc + cuts) draw reco properties
    else if( fWcheff == 3 )cLevelNum = "RecPID"; // RecPID
    else if( fWcheff == 4 )cLevelNum = "Reco"; // Rec (no cuts)
    else if( fWcheff == 5 )cLevelNum = "RecKineAcc"; // RecAcc
    else if( fWcheff == 6 )cLevelNum = "RecTrkMCPart"; // Rec (acc + cuts) draw kine properties
    
    
    Int_t nLimits=1;
    Double_t* newLimits =new Double_t[nLimits+1];
    newLimits[0]=fVarMin;
    newLimits[1]=fVarMax;
    const Int_t nnewLimits = nLimits;
    
    // Numerators Rebinning
    THnSparse* numData = (THnSparse*)gridSprsNumClone->GetGrid();
    THnSparse* newnumData = (THnSparse*)numData->Clone(Form("NumNew%0.1f%0.1f%s",fVarMin, fVarMax,fParticleType));
    newnumData->Reset();
    TAxis* axis = (TAxis*)newnumData->GetAxis(iVarContainer);
    axis->Set(nnewLimits,newLimits);
    newnumData->SetBinEdges(iVarContainer,newLimits);
    newnumData->RebinnedAdd(numData, 1);
    TH1D* Rheffnum = (TH1D*)newnumData->Projection(iVarContainer);
    
    // Denominators Rebinning
    THnSparse* denData = (THnSparse*)gridSprsDenClone->GetGrid();
    THnSparse* newdenData = (THnSparse*)denData->Clone(Form("DenNew%0.1f%0.1f%s",fVarMin, fVarMax,fParticleType));
    newdenData->Reset();
    TAxis* axis2 = (TAxis*)newdenData->GetAxis(iVarContainer);
    axis2->Set(nnewLimits,newLimits);
    newdenData->SetBinEdges(iVarContainer,newLimits);
    newdenData->RebinnedAdd(denData, 1);
    TH1D* Rheffden = (TH1D*)newdenData->Projection(iVarContainer);
    
    //Avg Efficiency rebinned
    TH1D* Rbheff = (TH1D*)Rheffden->Clone(Form("heffAvgw%0.1f%0.1f%s",fVarMin, fVarMax,fParticleType));
    Rbheff->Divide(Rheffnum,Rheffden,1,1,"B");
    
    Printf("--- Avarage %s eff: (%s: %0.1f-%0.1f) = %f",fWhichVar, fParticleType, fVarMin, fVarMax, Rbheff->GetBinContent(1));
    
    return Rbheff->GetBinContent(1);
    
}


// Setting of histogram name etc.
void SetHisto(TH1D* h, TString title= "", const char *Var="pt", const Int_t mrkColor=5){
    
    //General
    h->SetMarkerStyle(20);
    h->SetMarkerSize(0.82);
    h->SetMarkerColor(mrkColor); //blue
    
    //Title
    h->SetTitle(title.Data());
    h->SetName(title.Data());
    if(Var=="eta" || Var=="ETA" || Var=="Eta")
    {
        h->SetMaximum(h->GetMaximum()*1.25);
        h->SetMinimum(h->GetMinimum()*0.10);
    }
    //h->SetTitleFont(43);
    //h->SetTitleSize(25);
    
    //X axis
    h->GetXaxis()->SetTitle(Form("%s", Var));
    //h->GetXaxis()->SetTitleSize(1);
    //h->GetXaxis()->SetTitleOffset(0.50);
    //h->GetXaxis()->SetLabelFont(43);
    //h->GetXaxis()->SetLabelSize(10);
    //h->GetXaxis()->SetLabelOffset(0.007);
    
    //Y-axis
    h->GetYaxis()->SetTitle("counts");
    //h->GetYaxis()->SetTitleFont(43);
    //h->GetYaxis()->SetLabelFont(43);
    //h->GetYaxis()->SetLabelSize(15);
    //h->GetYaxis()->SetLabelOffset(0.007);
    
}


//Latex Name on histograms
void LatexName(const Double_t xvalue, TString fWhichVar="", Int_t runNumber, Int_t fbit=0, const char *fPartype=""){
    
    TLatex Tl;
    Tl.SetTextAlign(12);
    Tl.SetTextSize(0.035);
    Tl.DrawLatex(xvalue,0.49, Form("Configurations"));
    Tl.DrawLatex(xvalue,0.48, ""); //new commit 21March 2016
    Tl.DrawLatex(xvalue,0.42, Form("1. runNumber: %d", runNumber));
    Tl.DrawLatex(xvalue,0.36, Form("2. Filterbit: %d", fbit));
    Tl.DrawLatex(xvalue,0.30, Form("3. Particle: %s ", fPartype));
    
    if(fWhichVar=="pt" || fWhichVar=="PT" || fWhichVar=="Pt"){
        Tl.DrawLatex(xvalue,0.24, Form("4. %.1f #leq #eta #leq %.1f ", etamin, etamax));
        Tl.DrawLatex(xvalue,0.18, Form("5. %.1f #leq Zvtx #leq %.1f ", zvtxmin, zvtxmax));
    }else if(fWhichVar=="eta" || fWhichVar=="ETA" || fWhichVar=="Eta"){
        Tl.DrawLatex(xvalue,0.24, Form("4. %.1f #leq p_{T} #leq %.1f ", ptmin, ptmax));
        Tl.DrawLatex(xvalue,0.18, Form("5. %.1f #leq Zvtx #leq %.1f ", zvtxmin, zvtxmax));
    }else if(fWhichVar=="zvtx" || fWhichVar=="ZVTX" || fWhichVar=="Zvtx"){
        Tl.DrawLatex(xvalue,0.24, Form("4. %.1f #leq #eta #leq %.1f ", etamin, etamax));
        Tl.DrawLatex(xvalue,0.18, Form("5. %.1f #leq p_{T} #leq %.1f ", ptmin, ptmax));
    }else {
        Tl.DrawLatex(xvalue,0.18, Form("4. %.1f #leq Zvtx #leq %.1f ", zvtxmin, zvtxmax));
        Tl.DrawLatex(xvalue,0.24, Form("5. %.1f #leq #eta #leq %.1f ", etamin,  etamax));
        Tl.DrawLatex(xvalue,0.30, Form("6. %.1f #leq p_{T} #leq %.1f", ptmin,   ptmax));
    }
    
}


// EOF


