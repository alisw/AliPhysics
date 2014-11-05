/* pt generation limits */
const Double_t ptMin = 0.;
const Double_t ptMax = 10.;
/* y generation limits */
const Double_t yMin = -0.1;
const Double_t yMax = 0.1;
/* phi generation limits */
const Double_t phiMin = -TMath::Pi();
const Double_t phiMax = TMath::Pi();

const Double_t pionMass = 139.57018e-3;
const Double_t kaonMass = 493.677e-3;
const Double_t protonMass = 938.272046e-3;

PhiModel(Int_t ngen = 1.e6)
{
    
    Double_t m = TDatabasePDG::Instance()->GetParticle("phi")->Mass();
    Double_t w = TDatabasePDG::Instance()->GetParticle("phi")->Width();
    Double_t dm[2] = {kaonMass, kaonMass};
    Int_t di[2] = {3, -3};
    
    ResonanceModel(m, w, 2, dm, di, "PhiModel.root", ngen);
}

KStarModel(Int_t ngen = 1.e6, Bool_t antiparticle = kFALSE)
{
    
    Double_t m = TDatabasePDG::Instance()->GetParticle("K*0")->Mass();
    Double_t w = TDatabasePDG::Instance()->GetParticle("K*0")->Width();
    Double_t dm[2] = {kaonMass, pionMass};
    Int_t di[2] = {3, -2};
    if (antiparticle) {di[0] = -di[0]; di[1] = -di[1];}
    
    if (antiparticle) ResonanceModel(m, w, 2, dm, di, "KStarBarModel.root", ngen);
    else ResonanceModel(m, w, 2, dm, di, "KStarModel.root", ngen);
}

Lambda1520Model(Int_t ngen = 1.e6, Bool_t antiparticle = kFALSE)
{
    
    Double_t m = 1519.54e-3; /* GeV */
    Double_t w = 15.73e-3; /* GeV */
    Double_t dm[2] = {protonMass, kaonMass};
    Int_t di[2] = {4, -3};
    if (antiparticle) {di[0] = -di[0]; di[1] = -di[1];}
    
    if (antiparticle) ResonanceModel(m, w, 2, dm, di, "Lambda1520BarModel.root", ngen);
    else ResonanceModel(m, w, 2, dm, di, "Lambda1520Model.root", ngen);
}

ResonanceModel(Double_t resonanceMass, Double_t resonanceWidth, Int_t nDaughters, Double_t *daughterMass, Int_t *daughterID, const Char_t *outname, Int_t ngen = 1.e6)
{
    
    /* functions */
    TF1 *fBW = new TF1("fBW", "[0] * TMath::BreitWigner(x, [1], [2])");
    fBW->SetParameter(0, ngen / 100.);
    fBW->SetParameter(1, resonanceMass);
    fBW->SetParameter(2, resonanceWidth);
    
    /* histograms */
    Int_t massBins = 50;
    Int_t ptBins = 50;
    Double_t massMin = resonanceMass - 10. * resonanceWidth;
    Double_t massMax = resonanceMass + 10. * resonanceWidth;
    
    TH1F *hGenMass = new TH1F("hGenMass", ";m (GeV/#it{c}^{2});", massBins, massMin, massMax);
    hGenMass->SetMarkerColor(2);
    hGenMass->SetMarkerStyle(25);
    hGenMass->Sumw2();
    
    TH1F *hRecMass = new TH1F("hRecMass", ";m (GeV/#it{c}^{2});", massBins, massMin, massMax);
    hRecMass->SetMarkerColor(4);
    hRecMass->SetMarkerStyle(20);
    hRecMass->Sumw2();
    
    TH1F *hGenPt = new TH1F("hGenPt", ";#it{p}_{T} (GeV/#it{c});", ptBins, ptMin, ptMax);
    hGenPt->SetMarkerColor(2);
    hGenPt->SetMarkerStyle(25);
    hGenPt->Sumw2();
    
    TH1F *hRecPt = new TH1F("hRecPt", ";#it{p}_{T} (GeV/#it{c});", ptBins, ptMin, ptMax);
    hRecPt->SetMarkerColor(4);
    hRecPt->SetMarkerStyle(20);
    hRecPt->Sumw2();
    
    /* init random generator */
    TRandom *gRandom = new TRandom();
    gRandom->SetSeed(123456789);
    
    /* resonance four-momentum */
    TLorentzVector resonance;
    TLorentzVector recoP;
    TGenPhaseSpace decay;
    
    /* benchmark */
    TBenchmark bm;
    bm.Start("time");
    
    /* loop over generation */
    Int_t igen = 0;
    while (igen < ngen) {
        
        /* generate mass */
        Double_t mass = gRandom->BreitWigner(resonanceMass, resonanceWidth);
        /* generate pt */
        Double_t pt = gRandom->Uniform(ptMin, ptMax);
        /* generate y */
        Double_t y = gRandom->Uniform(yMin, yMax);
        /* get eta */
        Double_t eta = y2eta(pt, resonanceMass, y);
        /* generate phi */
        Double_t phi = gRandom->Uniform(phiMin, phiMax);
        
        /* init resonance four-momentum */
        resonance.SetPtEtaPhiM(pt, eta, phi, mass);
        /* decay particle, if allowed */
        if (!decay.SetDecay(resonance, nDaughters, daughterMass)) continue;
        
        /* decay allowed */
        igen++;
        decay.Generate();
        
        /* fill gen plots */
        hGenMass->Fill(mass);
        hGenPt->Fill(pt);
        
        /* reconstructed invariant mass */
        recoP.SetPtEtaPhiM(0., 0., 0., 0.);
        TLorentzVector *daughterP;
        Bool_t daughterLost = kFALSE;
        for (Int_t idaugh = 0; idaugh < nDaughters; idaugh++) {
            daughterP = decay.GetDecay(idaugh);
            recoP += *daughterP;
            Double_t eff = GetEfficiency(daughterP, daughterID[idaugh]);
            if (gRandom->Uniform(0., 1.) > eff) daughterLost = kTRUE;
        }
        
        /* fill reco plots */
        if (!daughterLost) {
            hRecMass->Fill(recoP.M());
            hRecPt->Fill(pt);
        }
    }
    
    /* benchmark */
    bm.Stop("time");
    bm.Print("time");
    
    gStyle->SetHistMinimumZero();
    
    TCanvas *cGenRec = new TCanvas("cGenRec", "cRecGen", 600, 600);
    cGenRec->Divide(1, 2);
    cGenRec->cd(1);
    hGenMass->DrawCopy();
    hRecMass->DrawCopy("same");
    cGenRec->cd(2);
    hGenPt->DrawCopy();
    hRecPt->DrawCopy("same");
    
    TCanvas *cEff = new TCanvas("cEff", "cEff", 600, 600);
    cEff->Divide(1, 2);
    cEff->cd(1);
    TH1 *hEffMass = (TH1 *)hRecMass->Clone("hEffMass");
    hEffMass->Divide(hEffMass, hGenMass, 1., 1., "B");
    hEffMass->DrawCopy();
    cEff->cd(2);
    TH1 *hEffPt = (TH1 *)hRecPt->Clone("hEffPt");
    hEffPt->Divide(hEffPt, hGenPt, 1., 1., "B");
    hEffPt->DrawCopy();
    
    TFile *fout = TFile::Open(outname, "RECREATE");
    hGenMass->Write();
    hRecMass->Write();
    hGenPt->Write();
    hRecPt->Write();
    hEffPt->Write();
    fout->Close();
    
}

/*****************************************************************/

TFile *fTrackEff = NULL;
TFile *fMatchEff = NULL;

Double_t
GetEfficiency(TLorentzVector *lv, Int_t ipart)
{
    
    Char_t *pname[5] = {"", "", "pion", "kaon", "proton"};
    Char_t *cname[2] = {"positive", "negative"};
    
    Int_t icharge = ipart > 0 ? 0 : 1;
    ipart = TMath::Abs(ipart);

    Double_t pt = lv->Pt();
    Double_t eta = lv->Eta();
    Double_t phi = lv->Phi();

    if (!fTrackEff) fTrackEff = TFile::Open("TOF_trackingEfficiency.root");
    if (!fMatchEff) fMatchEff = TFile::Open("TOF_matchingEfficiency.root");
    TH1 *hTrackEff = fTrackEff->Get(Form("hTrackingEff_MB_%s_%s", pname[ipart], cname[icharge]));
    TH1 *hMatchEff = fMatchEff->Get(Form("hMatchEff_MB_%s_%s", pname[ipart], cname[icharge]));
    
    Double_t trackEff = hTrackEff->Interpolate(pt);
    Double_t matchEff = hMatchEff->Interpolate(pt);

    Double_t trackEff_GF = TrackingEff_geantflukaCorrection(ipart, icharge)->Eval(pt);
    Double_t matchEff_GF = TOFmatchMC_geantflukaCorrection(ipart, icharge)->Eval(pt);
    
    Double_t eff = trackEff;
    if (ipart == 3 && pt > 0.6) eff = trackEff * matchEff;
    if (ipart == 4 && pt > 1.1) eff = trackEff * matchEff;
    
    Double_t pidEff = 0.954;
    eff *= pidEff;

//    return trackEff;
    return eff;
}

/*****************************************************************/

Double_t 
y2eta(Double_t pt, Double_t mass, Double_t y){
    Double_t mt = TMath::Sqrt(mass * mass + pt * pt);
    return TMath::ASinH(mt / pt * TMath::SinH(y));
}
Double_t 
eta2y(Double_t pt, Double_t mass, Double_t eta){
    Double_t mt = TMath::Sqrt(mass * mass + pt * pt);
    return TMath::ASinH(pt / mt * TMath::SinH(eta));
}

/*****************************************************************/

TF1 *fTrackEff_GF[5][2] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};
TF1 *fMatchEff_GF[5][2] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};

TF1 *
TrackingEff_geantflukaCorrection(Int_t ipart, Int_t icharge)
{
    
    if (fTrackEff_GF[ipart][icharge]) return fTrackEff_GF[ipart][icharge];
    
    Char_t *pname[3] = {"pion", "kaon", "proton"};
    Char_t *cname[2] = {"positive", "negative"};
    
    TF1 *f;
    if (ipart == 3 && icharge == 1)
        f = new TF1(Form("fGeantFluka_%s_%s", pname[ipart - 2], cname[icharge]), "TrackingPtGeantFlukaCorrectionKaMinus(x)", 0., 5.);
    else if (ipart == 4 && icharge == 1)
        f = new TF1(Form("fGeantFluka_%s_%s", pname[ipart - 2], cname[icharge]), "TrackingPtGeantFlukaCorrectionPrMinus(x)", 0., 5.);
    else
        f = new TF1(Form("fGeantFluka_%s_%s", pname[ipart - 2], cname[icharge]), "TrackingPtGeantFlukaCorrectionNull(x)", 0., 5.);
    
    fTrackEff_GF[ipart][icharge] = f;
    return f;
}

Double_t
TrackingPtGeantFlukaCorrectionNull(Double_t pTmc)
{
    return 1.;
}

Double_t
TrackingPtGeantFlukaCorrectionPrMinus(Double_t pTmc)
{
    return (1 - 0.129758 *TMath::Exp(-pTmc*0.679612));
}

Double_t
TrackingPtGeantFlukaCorrectionKaMinus(Double_t pTmc)
{
    return TMath::Min((0.972865 + 0.0117093*pTmc), 1.);
}

/*****************************************************************/

TF1 *
TOFmatchMC_geantflukaCorrection(Int_t ipart, Int_t icharge)
{

    if (fMatchEff_GF[ipart][icharge]) return fMatchEff_GF[ipart][icharge];
    
    Char_t *pname[3] = {"pion", "kaon", "proton"};
    Char_t *cname[2] = {"positive", "negative"};

    TF1 *f;
    if (ipart == 3 && icharge == 1)
        f = new TF1(Form("fGeantFluka_%s_%s", pname[ipart - 2], cname[icharge]), "MatchingPtGeantFlukaCorrectionKaMinus(x)", 0., 5.);
    else if (ipart == 4 && icharge == 1)
        f = new TF1(Form("fGeantFluka_%s_%s", pname[ipart - 2], cname[icharge]), "MatchingPtGeantFlukaCorrectionPrMinus(x)", 0., 5.);
    else
        f = new TF1(Form("fGeantFluka_%s_%s", pname[ipart - 2], cname[icharge]), "MatchingPtGeantFlukaCorrectionNull(x)", 0., 5.);
    
    fMatchEff_GF[ipart][icharge] = f;
    return f;
}

Double_t
MatchingPtGeantFlukaCorrectionNull(Double_t pTmc)
{
    return 1.;
}

Double_t
MatchingPtGeantFlukaCorrectionPrMinus(Double_t pTmc)
{
    Float_t ptTPCoutP =pTmc*(1-6.81059e-01*TMath::Exp(-pTmc*4.20094));
    return (TMath::Power(1 - 0.129758*TMath::Exp(-ptTPCoutP*0.679612),0.07162/0.03471));
}

Double_t
MatchingPtGeantFlukaCorrectionKaMinus(Double_t pTmc)
{
    Float_t ptTPCoutK=pTmc*(1- 3.37297e-03/pTmc/pTmc - 3.26544e-03/pTmc);
    return TMath::Min((TMath::Power(0.972865 + 0.0117093*ptTPCoutK,0.07162/0.03471)), 1.);
}