static TF1* f_wdsx_z;
static TF1* f_ta;
static TF2* f_ta_rfi;

void glauber()
{
//
//  Calculates some geometrical properties of PbPb collisions 
//  in the Glauber Model
//
//  Wood-Saxon nuclear density function
//
    TCanvas *c1 = new TCanvas("c1","Wood Saxon",400,10,600,700);
    TF1* f_wdsx = new TF1("f_wdsx", wdsx, 0, 15., 4);
    f_wdsx->SetParameter(0,6.624);
    f_wdsx->SetParameter(1,0.549);
    f_wdsx->SetParameter(2,0.000);
    f_wdsx->SetParameter(3,7.69e-4);
    f_wdsx->Draw();
//
//  Wood Saxon-nuclear density (b-z)
//
    TCanvas *c2 = new TCanvas("c2","Wood Saxon",400,10,600,700);
    TF2* f_wdsx_bz = new TF2("f_wdsx_bz", wdsx_bz, 0, 15., 0., 15., 4);
    f_wdsx_bz->SetParameter(0,6.624);
    f_wdsx_bz->SetParameter(1,0.549);
    f_wdsx_bz->SetParameter(2,0.000);
    f_wdsx_bz->SetParameter(3,7.69e-4);
    f_wdsx_bz->Draw();
//
//  Wood Saxon-nuclear density (z, for fixed b)
//
    TCanvas *c3 = new TCanvas("c3","Wood Saxon",400,10,600,700);
    f_wdsx_z = new TF1("f_wdsx_z", wdsx_z, 0, 15., 5);
    f_wdsx_z->SetParameter(0,6.624);
    f_wdsx_z->SetParameter(1,0.549);
    f_wdsx_z->SetParameter(2,0.000);
    f_wdsx_z->SetParameter(3,7.69e-4);
    f_wdsx_z->SetParameter(4,0.);
    f_wdsx_z->Draw();
//
//  Thickness function 
//
    TCanvas *c4 = new TCanvas("c4","T_A",400,10,600,700);
    f_ta = new TF1("f_ta", ta, 0, 15., 0);
    f_ta->Draw();
//
//  Kernel of overlap function
//    

    TCanvas *c5 = new TCanvas("c5","T_A",400,10,600,700);
    f_ta_rfi = new TF2("f_ta_rfi", ta_rfi, 0, 15., 0., TMath::Pi(), 1);
    f_ta_rfi->SetParameter(0,0.);
    f_ta_rfi->Draw();
//
//  Overlap Function
//
    TCanvas *c6 = new TCanvas("c6","T_AA",400,10,600,700);
    TF1* f_taa = new TF1("f_taa", taa,0.,15., 0);
    f_taa->Draw();

}



Double_t taa(Double_t* x, Double_t* dum)
{
    printf("taa %f\n", x[0]);
    
    Double_t b    = x[0];
    Double_t y    = 10./(208*208)*hijing->Profile((Float_t)b);
    return y;
}


Double_t wdsx(Double_t* x, Double_t* par)
{
//
//  Wood Saxon Parameterisation
//  as a function of radius
//
    Double_t xx  = x[0];
    Double_t r0  = par[0];
    Double_t d   = par[1];
    Double_t w   = par[2];
    Double_t n   = par[3];
    
    Double_t y  = n * (1.+w*(xx/r0)*(xx/r0))/(1.+TMath::Exp((xx-r0)/d));
    return y;
}

Double_t wdsx_bz(Double_t* x, Double_t* par)
{
//
//  Wood Saxon Parameterisation
//  as a function of z and  b
//
    Double_t bb  = x[0];
    Double_t zz  = x[1];
    Double_t r0  = par[0];
    Double_t d   = par[1];
    Double_t w   = par[2];
    Double_t n   = par[3];
    Double_t xx  = TMath::Sqrt(bb*bb+zz*zz);
    Double_t y  = n * (1.+w*(xx/r0)*(xx/r0))/(1.+TMath::Exp((xx-r0)/d));
    return y;
}

Double_t wdsx_z(Double_t* x, Double_t* par)
{
//
//  Wood Saxon Parameterisation
//  as a function of z for fixed b
//
    Double_t bb  = par[4];
    Double_t zz  = x[0];
    Double_t r0  = par[0];
    Double_t d   = par[1];
    Double_t w   = par[2];
    Double_t n   = par[3];
    Double_t xx  = TMath::Sqrt(bb*bb+zz*zz);
    Double_t y  = n * (1.+w*(xx/r0)*(xx/r0))/(1.+TMath::Exp((xx-r0)/d));
    return y;
}

Double_t ta(Double_t* x, Double_t* par)
{
//
//  Thickness function 
//
    Double_t b  = x[0];
    f_wdsx_z->SetParameter(4,b);
    Double_t y  = 2.*f_wdsx_z->Integral(0.,15.);
    return y;
}

Double_t ta_rfi(Double_t* x, Double_t* par)
{
//
//  Kernel for overlap function
//
    Double_t b    = par[0];
    Double_t r1   = x[0];
    Double_t phi  = x[1];
    Double_t r2   = TMath::Sqrt(r1*r1+b*b-2.*r1*b*TMath::Cos(phi));
    Double_t y    = r1*f_ta->Eval(r1)*f_ta->Eval(r2);
//    Double_t y    = r1*f_ta->Eval(r1);
    return y;
}


Double_t taa(Double_t* x, Double_t* par)
{
//
//  Overlap function
//

    
    Double_t b    = x[0];
    f_ta_rfi->SetParameter(0,b);
    f_ta_rfi->Update();
    
    Double_t y    = 2.*
	f_ta_rfi->Integral(0.,15., 0., TMath::Pi(), 0.001);
    printf("taa %f %f\n", x[0], y);
    return y;
}










