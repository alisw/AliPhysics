 void testPDF(Double_t a = 208.)
{

    gSystem->Load("liblhapdf.so");
    StrucFunc_t pdf = kCTEQ6ll;
    

// Dynamically link some shared libs                    

    
//    
    char    parm[20][20];
    Double_t val[20];
//
//  Initialization
//
    printf("PDF initialized with: \n \n");
    strncpy(parm[0], "DEFAULT             ",  20);
    val[0]  = AliStructFuncType::PDFsetIndex(pdf);
    AliStructFuncType::PdfSet(parm, val);

//
//  Plots
//
// ================================================================
// Pdf 
//
    TCanvas *c1 = new TCanvas("c1","Gluon PDF",400,10,600,700);
    TF1* f_gl_2 = new TF1("f_gl_2", pdf_gl, -6., 0., 2);
    f_gl_2->SetParameter(0,2);
    f_gl_2->SetParameter(1,a);
    f_gl_2->SetMaximum(1000.);
    f_gl_2->SetLineColor(1);
    f_gl_2->Draw();

    TF1* f_gl_5 = new TF1("f_gl_5", pdf_gl, -6., 0., 2);
    f_gl_5->SetParameter(0,5);
    f_gl_5->SetParameter(1,a);
    f_gl_5->SetLineColor(2);
    f_gl_5->Draw("Same");

    TF1* f_gl_10 = new TF1("f_gl_10", pdf_gl, -6., 0., 2);
    f_gl_10->SetParameter(0,10);
    f_gl_10->SetParameter(1,a);
    f_gl_10->SetLineColor(3);
    f_gl_10->Draw("Same");

    TF1* f_gl_50 = new TF1("f_gl_50", pdf_gl, -6., 0., 2);
    f_gl_50->SetParameter(0,50);
    f_gl_50->SetParameter(1,a);
    f_gl_50->SetLineColor(4);
    f_gl_50->Draw("Same");

// ================================================================
// Gluon Modification
//
    TCanvas *c2 = new TCanvas("c2","Gluon PDF Nucl. Mod.",400,10,600,700);

    TF1* fmod_gl_2 = new TF1("fmod_gl_2", mod_gl, -6., 0., 2);
    fmod_gl_2->SetParameter(0,2.);
    fmod_gl_2->SetParameter(1,a);
    fmod_gl_2->SetMaximum(1.5);
    fmod_gl_2->SetLineColor(1);
    fmod_gl_2->Draw();
    fmod_gl_2->GetHistogram()->SetXTitle("log(x)");
    fmod_gl_2->GetHistogram()->SetYTitle("F^{g}_{A}/F^{g}_{p}");

    TF1* fmod_gl_5 = new TF1("fmod_gl_5", mod_gl, -6., 0., 2);
    fmod_gl_5->SetParameter(0,5.);
    fmod_gl_5->SetParameter(1,a);
    fmod_gl_5->SetLineColor(2);
    fmod_gl_5->Draw("Same");

    TF1* fmod_gl_10 = new TF1("fmod_gl_10", mod_gl, -6., 0., 2);
    fmod_gl_10->SetParameter(0,10.);
    fmod_gl_10->SetParameter(1,a);
    fmod_gl_10->SetLineColor(3);
    fmod_gl_10->Draw("Same");

    TF1* fmod_gl_50 = new TF1("fmod_gl_50", mod_gl, -6., 0., 2);
    fmod_gl_50->SetParameter(0,50.);
    fmod_gl_50->SetParameter(1,a);
    fmod_gl_50->SetLineColor(4);
    fmod_gl_50->Draw("Same");
    label();
    

//  =================================================================
//  Sea Quark Modification
//
    TCanvas *c3 = new TCanvas("c3","Sea Quark PDF Nucl. Mod.",400,10,600,700);
    TF1* fmod_sq_2 = new TF1("fmod_sq_2", mod_sq, -6., 0., 2);
    fmod_sq_2->SetParameter(0,2.);
    fmod_sq_2->SetParameter(1,a);
    fmod_sq_2->SetMaximum(1.5);
    fmod_sq_2->SetLineColor(1);
    fmod_sq_2->Draw();
    fmod_sq_2->GetHistogram()->SetXTitle("log(x)");
    fmod_sq_2->GetHistogram()->SetYTitle("F^{sea q}_{A}/F^{sea q}_{p}");

    TF1* fmod_sq_5 = new TF1("fmod_sq_5", mod_sq, -6., 0., 2);
    fmod_sq_5->SetParameter(0,5.);
    fmod_sq_5->SetParameter(1,a);
    fmod_sq_5->SetLineColor(2);
    fmod_sq_5->Draw("Same");

    TF1* fmod_sq_10 = new TF1("fmod_sq_10", mod_sq, -6., 0., 2);
    fmod_sq_10->SetParameter(0,10.);
    fmod_sq_10->SetParameter(1,a);
    fmod_sq_10->SetLineColor(3);
    fmod_sq_10->Draw("Same");

    TF1* fmod_sq_50 = new TF1("fmod_sq_50", mod_sq, -6., 0., 2);
    fmod_sq_50->SetParameter(0,50.);
    fmod_sq_50->SetParameter(1,a);
    fmod_sq_50->SetLineColor(4);
    fmod_sq_50->Draw("Same");
    label();
    
//  =================================================================
//  Valence Quark Modification
//
    TCanvas *c4 = new TCanvas("c4","Valence Quark PDF Nucl. Mod.",
			      400,10,600,700);
    TF1* fmod_vq_2 = new TF1("fmod_vq_2", mod_vq, -6., 0., 2);
    fmod_vq_2->SetParameter(0,2.);
    fmod_vq_2->SetParameter(1,a);
    fmod_vq_2->SetMaximum(1.5);
    fmod_vq_2->SetLineColor(1);
    fmod_vq_2->Draw();
    fmod_vq_2->GetHistogram()->SetXTitle("log(x)");
    fmod_vq_2->GetHistogram()->SetYTitle("F^{val q}_{A}/F^{val q}_{p}");

    TF1* fmod_vq_5 = new TF1("fmod_vq_5", mod_vq, -6., 0., 2);
    fmod_vq_5->SetParameter(0,5.);
    fmod_vq_5->SetParameter(1,a);
    fmod_vq_5->SetLineColor(2);
    fmod_vq_5->SetMaximum(1.5);
    fmod_vq_5->Draw("Same");

    TF1* fmod_vq_10 = new TF1("fmod_vq_10", mod_vq, -6., 0., 2);
    fmod_vq_10->SetParameter(0,10.);
    fmod_vq_10->SetParameter(1,a);
    fmod_vq_10->SetLineColor(3);
    fmod_vq_10->SetMaximum(1.5);
    fmod_vq_10->Draw("Same");

    TF1* fmod_vq_50 = new TF1("fmod_vq_50", mod_vq, -6., 0., 2);
    fmod_vq_50->SetParameter(0,50.);
    fmod_vq_50->SetParameter(1,a);
    fmod_vq_50->SetLineColor(4);
    fmod_vq_50->SetMaximum(1.5);
    fmod_vq_50->Draw("Same");
    label();
    
}
 

Double_t pdf_gl(Double_t* x, Double_t* par)
{
    Double_t xx = TMath::Power(10.,x[0]);
    Double_t y;
    Double_t upv, dnv, usea, dsea, str, chm, bot, top, gl;
    Double_t q = par[0];
    Double_t a = par[1];

    AliStructFuncType::StructA(xx, q, a, upv, dnv, usea, 
				   dsea, str, chm, bot, top, gl);
    
    y = gl;
    
    return y;
}

Double_t mod_gl(Double_t* x, Double_t* par)
{
    Double_t xx = TMath::Power(10.,x[0]);
    Double_t y;
    Double_t upv, dnv, usea, dsea, str, chm, bot, top, gl;
    Double_t q = par[0];
    Double_t a = par[1];

    AliStructFuncType::StructA(xx, q, a, upv, dnv, usea, 
				   dsea, str, chm, bot, top, gl);
    
    Double_t y1 = gl;

    AliStructFuncType::StructA(xx, q, 1., upv, dnv, usea, 
				   dsea, str, chm, bot, top, gl);
    
    return y1/gl;
}

Double_t mod_sq(Double_t* x, Double_t* par)
{
    Double_t xx = TMath::Power(10.,x[0]);
    Double_t y;
    Double_t upv, dnv, usea, dsea, str, chm, bot, top, gl;
    Double_t q = par[0];
    Double_t a = par[1];

    AliStructFuncType::StructA(xx, q, a, upv, dnv, usea, 
				   dsea, str, chm, bot, top, gl);
    
    Double_t y1 = usea+dsea;

    AliStructFuncType::StructA(xx, q, 1., upv, dnv, usea, 
				   dsea, str, chm, bot, top, gl);
    
    return y1/(usea+dsea);
}


Double_t mod_vq(Double_t* x, Double_t* par)
{
    Double_t xx = TMath::Power(10.,x[0]);
    Double_t y;
    Double_t upv, dnv, usea, dsea, str, chm, bot, top, gl;
    Double_t q = par[0];
    Double_t a = par[1];

    AliStructFuncType::StructA(xx, q, a, upv, dnv, usea, 
				   dsea, str, chm, bot, top, gl);
    
    Double_t y1 = upv+dnv;

    AliStructFuncType::StructA(xx, q, 1., upv, dnv, usea, 
				   dsea, str, chm, bot, top, gl);
    
    return y1/(upv+dnv);
}


void label()
{
    TLine *line = new TLine(-5.7,1.4,-5.2,1.4);
    line->SetLineColor(1);
    line->SetLineWidth(3);
    line->Draw();
    tex = new TLatex(-5.1, 1.38,"Q =   2 GeV");
    tex->SetTextSize(0.04);
    tex->SetLineWidth(2);
    tex->Draw();

    TLine *line = new TLine(-5.7,1.3,-5.2,1.3);
    line->SetLineColor(2);
    line->SetLineWidth(3);
    line->Draw();
    tex = new TLatex(-5.1, 1.28,"Q =   5 GeV");
    tex->SetTextSize(0.04);
    tex->SetLineWidth(2);
    tex->Draw();


    TLine *line = new TLine(-5.7,1.2,-5.2,1.2);
    line->SetLineColor(3);
    line->SetLineWidth(3);
    line->Draw();
    tex = new TLatex(-5.1, 1.18,"Q = 10 GeV");
    tex->SetTextSize(0.04);
    tex->SetLineWidth(2);
    tex->Draw();


    TLine *line = new TLine(-5.7,1.1,-5.2,1.1);
    line->SetLineColor(4);
    line->SetLineWidth(3);
    line->Draw();
    tex = new TLatex(-5.1, 1.08,"Q = 50 GeV");
    tex->SetTextSize(0.04);
    tex->SetLineWidth(2);
    tex->Draw();
}



