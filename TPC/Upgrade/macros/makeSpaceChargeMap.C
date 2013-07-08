//
// Origin: Christian Lippman, CERN, Christian.Lippmann@cern.ch
//

int makeSpaceChargeMap(Double_t multiplicity = 950., Double_t intRate = 5e4, Double_t eps = 10.,
		       Double_t gasfactor = 1., string filename = "SpaceChargeMap.root",
		       Double_t radialScaling = 2., Double_t epsilonScaling = 2./3.) {
  //
  // Charge distribution is splitted into two (RZ and RPHI) in order to speed up
  // the needed calculation time. It is dumped to 
  //
  // Explanation of variables:
  // 1) multiplicity: charghed particle dn/deta for top 80% centrality (660 for 2011,
  //    expect 950 for full energy)
  // 2) intRate: Total interaction rate (e.g. 50kHz for the upgrade)
  // 3) eps: Number of backdrifting ions per primary electron (0 for MWPC, e.g.10 for GEM)
  // 4) gasfactor: Use different gas. E.g. Ar/CO2 has twice the primary ionization, ion drift
  //    velocity factor 2.5 slower, so  gasfactor = 5.
  //

  TFile *f = new TFile(filename.c_str(), "RECREATE");
  
  // some grid, not too coarse
  Int_t nr   = 350;
  Int_t nphi = 180;
  Int_t nz   = 500;

  const Double_t fgkIFCRadius=  83.5;     // radius which renders the "18 rod manifold"
  const Double_t fgkOFCRadius= 254.5;     // Mean Radius of the Outer Field Cage
  const Double_t fgke0 = 8.854187817e-12; // vacuum permittivity [A·s/(V·m)]

  Double_t dr = (fgkOFCRadius-fgkIFCRadius)/(nr+1);
  Double_t dphi = TMath::TwoPi()/(nphi+1);
  Double_t dz = 500./(nz+1);
  Double_t safty = 0.; // due to a root bug which does not interpolate the boundary ..
  // .. (first and last bin) correctly

  // Charge distribution in ZR (rotational symmetric) ------------------

  TH2F *histoZR = new TH2F("chargeZR", "chargeZR",
                           nr, fgkIFCRadius-dr-safty, fgkOFCRadius+dr+safty,
                           nz, -250-dz-safty, 250+dz+safty);

  // For the normalization to same integral as radial exponent = 2
  Double_t radialExponent             = -2.; // reference = 2
  Double_t radiusInner                = histoZR->GetXaxis()->GetBinCenter(1) / 100.;//in [m]
  Double_t radiusOuter                = histoZR->GetXaxis()->GetBinCenter(nr) / 100.;//in [m]
  Double_t integralRadialExponent2    = TMath::Power(radiusOuter,radialExponent+1) * 1./(radialExponent+1) 
    - TMath::Power(radiusInner,radialExponent+1) * 1./(radialExponent+1);
  
  radialExponent                      = -radialScaling; // user set   
  Double_t integralRadialExponentUser = 0.;
  if(radialScaling > 1 + 0.000001 || radialScaling < 1 - 0.000001 ) // to avoid n = -1
    integralRadialExponentUser = TMath::Power(radiusOuter,radialExponent+1) * 1./(radialExponent+1) 
      - TMath::Power(radiusInner,radialExponent+1) * 1./(radialExponent+1);
  else
    integralRadialExponentUser = TMath::Log(radiusOuter) - TMath::Log(radiusInner);
    
  Double_t normRadialExponent         = integralRadialExponent2 / integralRadialExponentUser;
 
  for (Int_t ir=1;ir<=nr;++ir) {
    Double_t rp = histoZR->GetXaxis()->GetBinCenter(ir);
    for (Int_t iz=1;iz<=nz;++iz) {
      Double_t zp = histoZR->GetYaxis()->GetBinCenter(iz);
      
      // recalculation to meter
      Double_t lZ = 2.5; // approx. TPC drift length
      Double_t rpM = rp/100.; // in [m]
      Double_t zpM = TMath::Abs(zp/100.); // in [m]
 
      // calculation of "scaled" parameters
      Double_t a = multiplicity*intRate/76628;
      //Double_t charge = gasfactor * ( a / (rpM*rpM) * (1 - zpM/lZ) ); // charge in [C/m^3/e0], no IBF
      Double_t charge = normRadialExponent * gasfactor * ( a / (TMath::Power(rpM,radialScaling)) * (1 - zpM/lZ + epsilonScaling*eps) ); // charge in [C/m^3/e0], with IBF
      
      charge = charge*fgke0;          // [C/m^3]

      // from MC simulation (Stefan)
      // for 50kHz
      Double_t kon = (2.62243e-09); // charge in [C/m^3]
      // Add to normal charge: gain 2000 with {0.25,0.5%) ion feedback
      //charge += eps*(kon/(rpM*rpM));

      if (zp<0) charge *= 0.9; // Slightly less on C side due to front absorber

      histoZR->SetBinContent(ir, iz, charge); 
    }
  }
    
  histoZR->Write("SpaceChargeInRZ");

  // Charge distribution in RPhi (e.g. Floating GG wire) ------------
  
  TH3F *histoRPhi = new TH3F("chargeRPhi", "chargeRPhi",
                             nr, fgkIFCRadius-dr-safty, fgkOFCRadius+dr+safty,
                             nphi, 0-dphi-safty, TMath::TwoPi()+dphi+safty,
                             2, -1, 1); // z part - to allow A and C side differences
  
  // some 'arbitrary' GG leaks
  Int_t   nGGleaks = 5;
  Double_t secPosA[5]    = {3,6,6,11,13};         // sector
  Double_t radialPosA[5] = {125,100,160,200,230}; // radius in cm
  Double_t secPosC[5]    = {1,8,12,15,15};        // sector
  Double_t radialPosC[5] = {245,120,140,120,190}; // radius in cm

  for (Int_t ir=1;ir<=nr;++ir) {
    Double_t rp = histoRPhi->GetXaxis()->GetBinCenter(ir);
    for (Int_t iphi=1;iphi<=nphi;++iphi) {
      Double_t phip = histoRPhi->GetYaxis()->GetBinCenter(iphi);
      for (Int_t iz=1;iz<=2;++iz) {
        Double_t zp = histoRPhi->GetZaxis()->GetBinCenter(iz);
        
        Double_t charge = 0;
        
        for (Int_t igg = 0; igg<nGGleaks; igg++) { // loop over GG leaks
          
          // A side
          Double_t secPos = secPosA[igg]; 
          Double_t radialPos = radialPosA[igg];

          if (zp<0) { // C side
            secPos = secPosC[igg]; 
            radialPos = radialPosC[igg];
          }

          // some 'arbitrary' GG leaks
          if (  (phip<(TMath::Pi()/9*(secPos+1)) && phip>(TMath::Pi()/9*secPos) ) ) { // sector slice
            if ( rp>(radialPos-2.5) && rp<(radialPos+2.5))  // 5 cm slice
              //charge = 300;
	      charge = 0.;
	  }
          
        }               
       
        charge = charge*fgke0; // [C/m^3]
        histoRPhi->SetBinContent(ir,iphi,iz,charge); 
      }
    }
  }

  histoRPhi->Write("SpaceChargeInRPhi");

  f->Close();
  
}
