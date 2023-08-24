#ifndef __CLING__
#include <iostream>
#include <cmath>

#include <AliExternalTrackParam.h>
#include <AliGeomManager.h>
#include <TGeoManager.h>
#include <TString.h>
#include <TMath.h>
#include "TH1.h"
#include "TFile.h"
#include "TRandom3.h"
#endif

constexpr int kNtracks{50000};
constexpr double kPt{10000.};
constexpr double kB{5.};
constexpr bool kFullWeights{true};
constexpr double kMeanVertex{1.38};
constexpr double kSigmaVertex{6.76};
constexpr bool kUseEtaSmearing{false};

double MeanMaterialBudget(const double *start, const double *end, double *mparam)
{
    //
    // Calculate mean material budget and material properties between
    //    the points "start" and "end".
    //
    // "mparam" - parameters used for the energy and multiple scattering
    //  corrections:
    //
    // mparam[0] - mean density: sum(x_i*rho_i)/sum(x_i) [g/cm3]
    // mparam[1] - equivalent rad length fraction: sum(x_i/X0_i) [adimensional]
    // mparam[2] - mean A: sum(x_i*A_i)/sum(x_i) [adimensional]
    // mparam[3] - mean Z: sum(x_i*Z_i)/sum(x_i) [adimensional]
    // mparam[4] - length: sum(x_i) [cm]
    // mparam[5] - Z/A mean: sum(x_i*Z_i/A_i)/sum(x_i) [adimensional]
    // mparam[6] - number of boundary crosses
    //
    
    mparam[0]=0; mparam[1]=1; mparam[2] =0; mparam[3] =0;
    mparam[4]=0; mparam[5]=0; mparam[6] =0;
    //
    double bparam[6]; // total parameters
    double lparam[6]; // local parameters
    
    for (int i=0;i<6;i++) bparam[i]=0;
    
    if (!gGeoManager) {
        std::cerr << "No TGeo\n" << std::endl;
        return 0.;
    }
    //
    double length;
    double dir[3];
    length = std::sqrt((end[0]-start[0])*(end[0]-start[0])+
                       (end[1]-start[1])*(end[1]-start[1])+
                       (end[2]-start[2])*(end[2]-start[2]));
    mparam[4]=length;
    if (length<TGeoShape::Tolerance()) return 0.0;
    double invlen = 1./length;
    dir[0] = (end[0]-start[0])*invlen;
    dir[1] = (end[1]-start[1])*invlen;
    dir[2] = (end[2]-start[2])*invlen;
    
    // Initialize start point and direction
    TGeoNode *currentnode = 0;
    TGeoNode *startnode = gGeoManager->InitTrack(start, dir);
    if (!startnode) {
        std::cout << Form("start point out of geometry: x %f, y %f, z %f",
                          start[0],start[1],start[2]) << std::endl;
        return 0.0;
    }
    TGeoMaterial *material = startnode->GetVolume()->GetMedium()->GetMaterial();
    lparam[0]   = material->GetDensity();
    lparam[1]   = material->GetRadLen();
    lparam[2]   = material->GetA();
    lparam[3]   = material->GetZ();
    lparam[4]   = length;
    lparam[5]   = lparam[3]/lparam[2];
    if (material->IsMixture()) { // first material before SPD: BODY_Air (mixture)
        TGeoMixture * mixture = (TGeoMixture*)material;
        lparam[5] =0;
        double sum =0;
        for (int iel=0;iel<mixture->GetNelements();iel++){
            sum  += mixture->GetWmixt()[iel];
            lparam[5]+= mixture->GetZmixt()[iel]*mixture->GetWmixt()[iel]/mixture->GetAmixt()[iel];
        }
        lparam[5]/=sum;
        
    }
    
    // Locate next boundary within length without computing safety.
    // Propagate either with length (if no boundary found) or just cross boundary
    gGeoManager->FindNextBoundaryAndStep(length, kFALSE);
    double step = 0.0; // Step made
    double snext = gGeoManager->GetStep();
    // If no boundary within proposed length, return current density
    if (!gGeoManager->IsOnBoundary()) {
        mparam[0] = lparam[0];
        mparam[1] = lparam[4]/lparam[1];
        mparam[2] = lparam[2];
        mparam[3] = lparam[3];
        mparam[4] = lparam[4];
        return lparam[0];
    }
    
    // Try to cross the boundary and see what is next
    int nzero = 0;
    while (length>TGeoShape::Tolerance()) {
        currentnode = gGeoManager->GetCurrentNode();
        if (snext<2.*TGeoShape::Tolerance()) nzero++;
        else nzero = 0;
        // if (nzero>3) {
//            // This means navigation has problems on one boundary
//            // Try to cross by making a small step
//            
//            std:: cout << "Problem on a boundary!!" << std::endl;
//            
//            static int show_error = !(getenv("HLT_ONLINE_MODE") && strcmp(getenv("HLT_ONLINE_MODE"), "on") == 0);
//            if (show_error) std::cerr << "Cannot cross boundary" << std::endl;
//            mparam[0] = bparam[0]/step;
//            mparam[1] = bparam[1];
//            mparam[2] = bparam[2]/step;
//            mparam[3] = bparam[3]/step;
//            mparam[5] = bparam[5]/step;
//            mparam[4] = step;
//            mparam[0] = 0.;             // if crash of navigation take mean density 0
//            mparam[1] = 1000000;        // and infinite rad length
//            return bparam[0]/step;
        // }
        mparam[6]+=1.;
        step += snext;
        bparam[1]    += snext/lparam[1];
        bparam[2]    += snext*lparam[2];
        bparam[3]    += snext*lparam[3];
        bparam[5]    += snext*lparam[5];
        bparam[0]    += snext*lparam[0];
        
//        if(material->IsMixture()) {
//            std::cout << "material: " << material->GetName() << "  (Is Mixture!)" << std::endl;
//        } else {
//            std::cout << std::endl << std::endl << "material: " << material->GetName() << std::endl;
//        }
//        std::cout << "density: " << lparam[0] << "  rad. length: " << lparam[1] << "  A: " << lparam[2] << "  Z: " << lparam[3] << "  length: " << snext << "  X/X0: " << snext/lparam[1] << std::endl << std::endl;
        
        if (snext>=length) break;
        if (!currentnode) break;
        length -= snext;
        material = currentnode->GetVolume()->GetMedium()->GetMaterial();
        lparam[0] = material->GetDensity();
        lparam[1]  = material->GetRadLen();
        lparam[2]  = material->GetA();
        lparam[3]  = material->GetZ();
        lparam[5]   = lparam[3]/lparam[2];
        if (material->IsMixture()) {
            TGeoMixture * mixture = (TGeoMixture*)material;
            lparam[5]=0;
            double sum =0;
            for (int iel=0;iel<mixture->GetNelements();iel++){
                sum+= mixture->GetWmixt()[iel];
                lparam[5]+= mixture->GetZmixt()[iel]*mixture->GetWmixt()[iel]/mixture->GetAmixt()[iel];
            }
            lparam[5]/=sum;
        }
        
        gGeoManager->FindNextBoundaryAndStep(length, kFALSE);
        snext = gGeoManager->GetStep();
    }
    mparam[0] = bparam[0]/step;
    mparam[1] = bparam[1];
    mparam[2] = bparam[2]/step;
    mparam[3] = bparam[3]/step;
    mparam[5] = bparam[5]/step;
    return bparam[0]/step;
}

void InspectGeo() {
    AliGeomManager::LoadGeometry("geometry.root");
    
    double Rin   = 0.0;
    double Rout  = 1.0;
    
    double sum   = 0.0;
    double niter = 0.0;
    
    TH1F* hist_material_budget = new TH1F("hist_material_budget","Material budget vs r",400,0,400);
    hist_material_budget->GetXaxis()->SetTitle("r [cm]");
    hist_material_budget->GetYaxis()->SetTitle("Mat. budget [X_{0}]");
    
    TH1F* hist_a = new TH1F("hist_a","; r (cm); A",400,0,400);
    TH1F* hist_z = new TH1F("hist_z","; r (cm); Z",400,0,400);
    TH1F* hist_rho = new TH1F("density","density; r (cm); #rho",400,0,400);

    std::vector<std::array<double,3>> directions;
    directions.reserve(kNtracks);

    std::vector<double> material(kNtracks,0);
    std::vector<std::array<double,3>> currentPoint;
    currentPoint.reserve(kNtracks);
    double cov[21]{1.};
    for (int i{0}; i < kNtracks; ++i) {
        double eta = kUseEtaSmearing ? gRandom->Uniform(-0.9,0.9) : 0.;
        double theta= kUseEtaSmearing ? 2*TMath::ATan(TMath::Exp(-eta)) : TMath::PiOver2();
        double phi = gRandom->Uniform(0,TMath::TwoPi());
        double shift{gRandom->Gaus(kMeanVertex, kSigmaVertex)};

        while (std::abs(shift) > 10.)
            shift = gRandom->Gaus(kMeanVertex, kSigmaVertex);
        currentPoint.push_back({0.,0.,shift});

        directions.push_back({std::cos(phi) * std::sin(theta), std::sin(phi) * std::sin(theta), std::cos(theta)});
    }

    constexpr int step = 1;
    for (int i = 1; i<400; i++){
        
        sum = 0.0;
        double rho{0.}, a{0.}, z{0.};

        for (int it = 0; it < kNtracks; it+=step){
            const auto& dir{directions[it]};
            auto& pnt{currentPoint[it]};
            
            double endPoint[3]{pnt[0] + dir[0] * step, pnt[1] + dir[1] * step, pnt[2] + dir[2] * step};
            
            const double radius = std::hypot(endPoint[0], endPoint[1]);
            if (std::abs(radius - i) > 1.e-6) {
                std::cout << "Radius error " << radius << "\t" << i << ". Aborting." << std::endl;
                return;
            } 

            double param[7];
            
            MeanMaterialBudget(pnt.data(), endPoint, param);

            material[it] += param[1];
            sum   += material[it];

            rho += param[0];
            a += param[2] * (kFullWeights ? param[0] : 1.);
            z += param[3] * (kFullWeights ? param[0] : 1.);
    
            currentPoint[it][0] = endPoint[0];
            currentPoint[it][1] = endPoint[1];
            currentPoint[it][2] = endPoint[2];
        }

        double wA = a / (kFullWeights ? rho : kNtracks);
        double wZ = z / (kFullWeights ? rho : kNtracks);
        hist_a->SetBinContent(i, wA);
        hist_a->SetBinError(i, kNtracks / rho);
        
        hist_z->SetBinContent(i, wZ);
        hist_z->SetBinError(i, kNtracks / rho);

        hist_rho->SetBinContent(i, rho / kNtracks);

        hist_material_budget->Fill(i, sum / kNtracks);
        
        std::cout << "Mean material budget: " << sum/kNtracks << "\t A = " << wA << "\t Z = " << wZ << std::endl;
        std::cout << "R [cm] = : " << i << std::endl;
        Rin+=1.0;
        Rout+=1.0;
        
    }

    int tpctof[2]{168,370};
    for (int det{0}; det < 2; ++det) {
        double totA{0}, totZ{0}, rho{0};
        for (int iBin{1}; iBin <= tpctof[det]; ++iBin) {
            double density = 1 / hist_a->GetBinError(iBin);
            totA += hist_a->GetBinContent(iBin) * density;
            totZ += hist_z->GetBinContent(iBin) * density;
            rho += density;
        }
        totA /= rho; 
        totZ /= rho;
        std::cout << "A " << totA  << "\t Z " << totZ << std::endl;
    }
    
    TFile* fileout = new TFile("FastMeanX0R.root","RECREATE");
    hist_material_budget->Write("hist_material_budget");
    hist_a->Write();
    hist_z->Write();
    hist_rho->Write();
    fileout->Close();
    
    hist_material_budget->Draw("hist");
    
}
