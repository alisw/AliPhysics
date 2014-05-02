#include <TF1.h>
#include <TH1.h>
#include <TMath.h>
#include <cmath>

struct PDF 
{
  
  static double cheb(double v, double* p, double* q, int n) {
    double num = 0;
    double den = 0;
    for (int i = n-1; i >= 0; i--) {
      num = (p[i] + num)*v;
      den = (q[i] + den)*v;
    }
    return num/den;
  }
  static double f1(double v) { 
    static double a1[3] = {0.04166666667,-0.01996527778, 0.02709538966};

    double u = std::exp(v+1);
    if (u < 1e-10) return 0;
    
    double ue = std::exp(-1/u);
    double us = std::sqrt(u);
    return 0.3989422803*(ue/us)*(1+(a1[0]+(a1[1]+a1[2]*u)*u)*u);
  }
  static double f2(double v) {
    static double p1[] = {0.4259894875,  -0.1249762550, 
			  0.03984243700, -0.006298287635,   
			  0.001511162253};
    static double q1[] = {1.0,           -0.3388260629, 
			  0.09594393323, -0.01608042283,    
			  0.003778942063};

    double u = std::exp(-v-1);
    return std::exp(-u)*std::sqrt(u)*cheb(v,p1,q1,5);
  }
  static double f3(double v) { 
    static double p2[] = {0.1788541609,   0.1173957403, 
			  0.01488850518, -0.001394989411,   
			  0.0001283617211};
    static double q2[] = {1.0,            0.7428795082, 
			  0.3153932961,   0.06694219548,    
			  0.008790609714};
    return cheb(v,p2,q2,5);
  }
  static double f4(double v) { 
    static double p3[] = {0.1788544503,   0.09359161662,
			  0.006325387654, 0.00006611667319,
			  -0.000002031049101};
    static double q3[] = {1.0,            0.6097809921, 
			  0.2560616665,   0.04746722384,    
			  0.006957301675};
    return cheb(v,p3,q3,5);
  }
  static double f5(double v) { 
    static double p4[] = {0.9874054407, 118.6723273,  
			  849.2794360, -743.7792444,      
			  427.0262186};
    static double q4[] = {1.0,           106.8615961,  
			  337.6496214,   2016.712389,      
			  1597.063511};
    double u   = 1/v;
    return u*u*cheb(u,p4,q4,5);
  }
  static double f6(double v) { 
    static double p5[] = {1.003675074,   167.5702434,  
			   4789.711289, 21217.86767,   
			   -22324.94910};
    static double q5[] = {1.0,          156.9424537,
			   3745.310488,  9834.698876,
			   66924.28357};
    double u   = 1/v;
    return u*u*cheb(u,p5,q5,5);
  }
  static double f7(double v) { 
    static double p6[] = {1.000827619,  664.9143136,  
			  62972.92665,  475554.6998,
			  -5743609.109};
    static double q6[] = {1.0,          651.4101098,  
			  56974.73333,  165917.4725,
			  -2815759.939};
    double u   = 1/v;
    return u*u*cheb(u,p6,q6,5);
  }
  static double f8(double v) { 
    static double a2[2] = {-1.845568670,-4.284640743};
    double u = 1/(v-v*std::log(v)/(v+1));
    return u*u*(1+(a2[0]+a2[1]*u)*u);
  }
    
  static Double_t wrap(Double_t* px, Double_t* pp) { 
    Int_t    i = pp[0];
    Double_t v = px[0];
    switch (i) { 
    case 1: return f1(v);
    case 2: return f2(v);
    case 3: return f3(v);
    case 4: return f4(v);
    case 5: return f5(v);
    case 6: return f6(v);
    case 7: return f7(v);
    case 8: return f8(v);
    }
    return 0;
  }

  static TF1* getComp(Int_t i) { 
    Double_t l = -10;
    Double_t h = -5.5;
    switch (i) { 
    case 2: l =   -5.5; h =   -1; break;
    case 3: l =   -1;   h =    1; break;
    case 4: l =    1;   h =    5; break;
    case 5: l =    5;   h =   12; break;
    case 6: l =   12;   h =   50; break;
    case 7: l =   50;   h =  300; break;
    case 8: l =  300;   h = 1000; break;
    default: i = 1; break;
    }
    
    l *= (l < 0 ? 1.5 : 0.4);
    h *= (h < 0 ? 0.5 : 1.5);

    TF1* f = new TF1(Form("f%d", i), &PDF::wrap, l, h, 1);
    f->SetParameter(0, i);
    f->SetLineColor(i);
    return f;
  }
  static void draw() { 
    TH1* h = new TH1F("h","h", 1000, -10, 1000);
    h->Draw();

    TF1* l = new TF1("l", "landau", -10, 1000);
    l->SetParameters(1,0,1);
    l->SetLineStyle(2);
    l->Draw("same");

    Double_t max = 0;
    for (Int_t i = 1; i <= 8; i++) { 
      TF1* f = getComp(i);
      f->Draw("same");
      max = TMath::Max(f->GetMaximum(),max);
    }
    h->SetMaximum(1.1*max);
  }
};


    
