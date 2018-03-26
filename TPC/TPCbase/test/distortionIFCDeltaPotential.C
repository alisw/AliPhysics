/*
  .L $AliRoot_SRC/TPC/TPCbase/test/distortionIFCDeltaPotential.C
  distortionIFCDeltaPotential(40,1000,17,17,90) // zShort, amplitude, rRow, zCol, phiSlice

*/

Double_t dFunctionVZ(Double_t *x, Double_t *par);
const Double_t gkVCE  = -50000;
const Double_t gkVROC = -70;
const Int_t gkMaxIter = 200;
const Float_t gkConvError = 1e-8;
const Double_t gkMaxEpsilon = 1e-2;


/// Input
/// \param zShort     -  position with shortcut at IFC
/// \param amplitude  -  deltaV at zShort
/// deltaV= ....
/// input density is equal to 0
/// \output - root file with distortion object and and residual trees
void distortionIFCDeltaPotential(Double_t zShort, Double_t amplitude, const Int_t rRow, const Int_t zColumn, const Int_t phiSlice){
    // set potential boundary in V
    TF1 * potentialBoundaryFunctionInZ = new TF1("dFunctionVZ", dFunctionVZ, -250.0, 250.0,2);
    potentialBoundaryFunctionInZ->SetParName(0,"zShort");
    potentialBoundaryFunctionInZ->SetParName(1,"amplitude");
    potentialBoundaryFunctionInZ->SetParameter(0,zShort);
    potentialBoundaryFunctionInZ->SetParameter(1,amplitude);
    potentialBoundaryFunctionInZ->Draw();


    // charge in histogram
    TH3 * chargeA = new TH3F ("charge A","charge A",250, 0.0, TMath::TwoPi(), 250, 85.0, 250.0, 250,0.0,250.0);
    TH3 * chargeC = new TH3F ("charge C","charge C",250, 0.0, TMath::TwoPi(), 250, 85.0, 250.0, 250,0.0,250.0);

    ::Info( "distortionIFCDeltaPotential","Begin");
    AliTPCSpaceCharge3DDriftLine * sc = new AliTPCSpaceCharge3DDriftLine(TString::Format("distortionIFCDeltaPotentialRRow%dZCol%dPhiSlice%d",rRow,zColumn,phiSlice).Data(),TString::Format("distortionIFCDeltaPotentialRRow%dZCol%dPhiSlice%d",rRow,zColumn,phiSlice).Data(), rRow,zColumn,phiSlice);

    // set the input boundary potential and charge densities
    sc->SetInputSpaceCharge(chargeA,0);
    sc->SetInputSpaceCharge(chargeC,1);
    sc->SetBoundaryIFCA(potentialBoundaryFunctionInZ);


    // set type of correction map
    //sc->SetCorrectionType(1);

    // create Poisson Solver
    AliTPCPoissonSolver *poissonSolver = new AliTPCPoissonSolver();
    sc->SetPoissonSolver(poissonSolver);
    sc->SetOmegaTauT1T2(0.35,1.,1.);
    // space charge parameters
    //Double_t c0 = 0.904466;
    //Double_t c1 = 0.904466;// set constant parameters
    //sc->Init();

    sc->InitSpaceCharge3DPoissonIntegralDz(rRow,zColumn,phiSlice,gkMaxIter,gkConvError);

    // create tree or distortion map
    ::Info( "distortionIFCDeltaPotential","Creating distortion tree");
    sc->CreateDistortionTree(5.0);

    // create histogram
    //TCanvas *can1 = new TCanvas("histogramDR");
    //sc->CreateHistoDistDRinXY(10,200,200)->Draw();

}


// function for dV(z)
Double_t dFunctionVZ(Double_t *x, Double_t *par)
{
    /// dU=0 at (0 and at 250)
    /// -k0*z     at z(0,z0)
    /// -k0*z0+k1*z   at z(z0,250)
    /// -k0*z0+k1*250=0
    ///  k1=k0*z0/250.

    Float_t xx =x[0];
    const Float_t z0 =par[0]; // U at central electrode
    const Float_t k0 =par[1]; // z0
    //const Float_t k1 = k0 * z0 / AliTPCPoissonSolver::fgkTPCZ0; // k1

    const Float_t k1 = k0  / (AliTPCPoissonSolver::fgkTPCZ0 - z0);
    Double_t f = 0.0;

    // z location near
    if (TMath::Abs(xx) < TMath::Abs(z0)) {
        // ideal case
        f = -k0 * (TMath::Abs(xx) / TMath::Abs(z0)) ;
    } else {
        f = -k0 + k1 * (TMath::Abs(xx) - z0);
    }

    return f;
}


