//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Module: EvtBtoXsllUtil.cc
//
// Description: Routine to generate non-resonant B -> Xs l+ l- decays.
// It generates a dilepton mass spectrum according to
// F.Kruger and L.M.Sehgal, Phys. Lett. B380, 199 (1996)
// and then generates the two lepton momenta according to
// A.Ali, G.Hiller, L.T.Handoko and T.Morozumi, Phys. Rev. D55, 4105 (1997).
// Expressions for Wilson coefficients and power corrections are taken
// from A.Ali, E.Lunghi, C.Greub and G.Hiller, Phys. Rev. D66, 034002 (2002).
// Detailed formulae for shat dependence of these coefficients are taken
// from H.H.Asatryan, H.M.Asatrian, C.Greub and M.Walker, PRD65, 074004 (2002)
// and C.Bobeth, M.Misiak and J.Urban, Nucl. Phys. B574, 291 (2000).
// The resultant Xs particles may be decayed by JETSET.
//
// Modification history:
//
//    Stephane Willocq    Jan 19, 2001   Module created
//    Stephane Willocq    Nov  6, 2003   Update Wilson Coeffs & dG's
//    &Jeff Berryhill
//
//------------------------------------------------------------------------
//
#include "EvtGenBase/EvtPatches.hh"
//
#include <stdlib.h>
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenModels/EvtBtoXsllUtil.hh"
#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtDiLog.hh"

EvtComplex EvtBtoXsllUtil::GetC7Eff0(double sh, bool nnlo) 
{
  // This function returns the zeroth-order alpha_s part of C7

  if (!nnlo) return -0.313;

  double A7;

  // use energy scale of 2.5 GeV as a computational trick (G.Hiller)
  // at least for shat > 0.25
  A7 = -0.353 + 0.023;

  EvtComplex c7eff;
  if (sh > 0.25)
  { 
    c7eff = A7;
    return c7eff;
  }

  // change energy scale to 5.0 for full NNLO calculation below shat = 0.25
  A7 = -0.312 + 0.008;
  c7eff = A7;

  return c7eff;
}

EvtComplex EvtBtoXsllUtil::GetC7Eff1(double sh, double mbeff, bool nnlo) 
{
  // This function returns the first-order alpha_s part of C7

  if (!nnlo) return 0.0;
  double logsh;
  logsh = log(sh);

  EvtComplex uniti(0.0,1.0);

  EvtComplex c7eff = 0.0;
  if (sh > 0.25)
  { 
    return c7eff;
  }

  // change energy scale to 5.0 for full NNLO calculation below shat = 0.25
  double muscale = 5.0;
  double alphas = 0.215;
  //double A7 = -0.312 + 0.008;
  double A8 = -0.148;
  //double A9 = 4.174 + (-0.035);
  //double A10 = -4.592 + 0.379;
  double C1 = -0.487;
  double C2 = 1.024;
  //double T9 = 0.374 + 0.252;
  //double U9 = 0.033 + 0.015;
  //double W9 = 0.032 + 0.012;
  double Lmu = log(muscale/mbeff);

  EvtComplex F71;
  EvtComplex f71;
  EvtComplex k7100(-0.68192,-0.074998);
  EvtComplex k7101(0.0,0.0);
  EvtComplex k7110(-0.23935,-0.12289);
  EvtComplex k7111(0.0027424,0.019676);
  EvtComplex k7120(-0.0018555,-0.175);
  EvtComplex k7121(0.022864,0.011456);
  EvtComplex k7130(0.28248,-0.12783);
  EvtComplex k7131(0.029027,-0.0082265);
  f71 = k7100 + k7101*logsh + sh*(k7110 + k7111*logsh) +
        sh*sh*(k7120 + k7121*logsh) + 
        sh*sh*sh*(k7130 + k7131*logsh); 
  F71 = (-208.0/243.0)*Lmu + f71;

  EvtComplex F72;
  EvtComplex f72;
  EvtComplex k7200(4.0915,0.44999);
  EvtComplex k7201(0.0,0.0);
  EvtComplex k7210(1.4361,0.73732);
  EvtComplex k7211(-0.016454,-0.11806);
  EvtComplex k7220(0.011133,1.05);
  EvtComplex k7221(-0.13718,-0.068733);
  EvtComplex k7230(-1.6949,0.76698);
  EvtComplex k7231(-0.17416,0.049359);
  f72 = k7200 + k7201*logsh + sh*(k7210 + k7211*logsh) +
        sh*sh*(k7220 + k7221*logsh) + 
        sh*sh*sh*(k7230 + k7231*logsh); 
  F72 = (416.0/81.0)*Lmu + f72;
  
  EvtComplex F78;
  F78 = (-32.0/9.0)*Lmu + 8.0*EvtConst::pi*EvtConst::pi/27.0 + (-44.0/9.0) 
        + (-8.0*EvtConst::pi/9.0)*uniti +
        (4.0/3.0*EvtConst::pi*EvtConst::pi - 40.0/3.0)*sh +
        (32.0*EvtConst::pi*EvtConst::pi/9.0 - 316.0/9.0)*sh*sh +
        (200.0*EvtConst::pi*EvtConst::pi/27.0 - 658.0/9.0)*sh*sh*sh +
    (-8.0*logsh/9.0)*(sh + sh*sh + sh*sh*sh);
        
  c7eff = - alphas/(4.0*EvtConst::pi)*(C1*F71 + C2*F72 + A8*F78);

  return c7eff;
}


EvtComplex EvtBtoXsllUtil::GetC9Eff0(double sh, double /* mbeff */,
                                     bool nnlo, bool btod) 
{
  // This function returns the zeroth-order alpha_s part of C9

  if (!nnlo) return 4.344;
  double mch = 0.29;

  double A9;
  A9 = 4.287 + (-0.218);
  double C1;
  C1 = -0.697;
  double C2;
  C2 = 1.046;
  double T9;
  T9 = 0.114 + 0.280;
  double U9;
  U9 = 0.045 + 0.023;
  double W9;
  W9 = 0.044 + 0.016;

  EvtComplex uniti(0.0,1.0);

  EvtComplex hc;
  double xarg;
  xarg = 4.0*mch/sh;

  hc = -4.0/9.0*log(mch*mch) + 8.0/27.0 + 4.0*xarg/9.0;
  if (xarg < 1.0)
  { 
    hc = hc - 2.0/9.0*(2.0 + xarg)*sqrt(fabs(1.0 - xarg))*
      (log((sqrt(1.0 - xarg)+1.0)/(sqrt(1.0 - xarg) - 1.0)) - 
       uniti*EvtConst::pi);
  } 
  else
  {
    hc = hc - 2.0/9.0*(2.0 + xarg)*sqrt(fabs(1.0 - xarg))*
      2.0*atan(1.0/sqrt(xarg-1.0));
  }

  EvtComplex h1;
  xarg = 4.0/sh;
  h1 = 8.0/27.0 + 4.0*xarg/9.0;
  if (xarg < 1.0)
  { 
    h1 = h1 - 2.0/9.0*(2.0 + xarg)*sqrt(fabs(1.0 - xarg))*
      (log((sqrt(1.0 - xarg)+1.0)/(sqrt(1.0 - xarg) - 1.0)) - 
       uniti*EvtConst::pi);
  } 
  else
  {
    h1 = h1 - 2.0/9.0*(2.0 + xarg)*sqrt(fabs(1.0 - xarg))*
      2.0*atan(1.0/sqrt(xarg-1.0));
  }

  EvtComplex h0;
  h0 = 8.0/27.0 - 4.0*log(2.0)/9.0 + 4.0*uniti*EvtConst::pi/9.0;


  // X=V_{ud}^* V_ub / V_{td}^* V_tb * (4/3 C_1 +C_2) * (h(\hat m_c^2, hat s)-
  // h(\hat m_u^2, hat s))
  EvtComplex Vudstar(1.0 - 0.2279*0.2279/2.0, 0.0);
  EvtComplex Vub((0.118+0.273)/2.0, -1.0*(0.305+0.393)/2.0);
  EvtComplex Vtdstar(1.0 - (0.118+0.273)/2.0,(0.305+0.393)/2.0);
  EvtComplex Vtb(1.0,0.0);

  EvtComplex Xd;
  Xd = (Vudstar * Vub / Vtdstar * Vtb) * (4.0/3.0*C1 + C2) * (hc - h0);

  EvtComplex c9eff = 4.344;
  if (sh > 0.25)
  { 
    c9eff =  A9 + T9*hc + U9*h1 + W9*h0;
    if (btod)
    {
      c9eff += Xd; 
    }
    return c9eff;
  }

  // change energy scale to 5.0 for full NNLO calculation below shat = 0.25
  A9 = 4.174 + (-0.035);
  C1 = -0.487;
  C2 = 1.024;
  T9 = 0.374 + 0.252;
  U9 = 0.033 + 0.015;
  W9 = 0.032 + 0.012;

  Xd = (Vudstar * Vub / Vtdstar * Vtb) * (4.0/3.0*C1 + C2) * (hc - h0);

  c9eff = A9 + T9*hc + U9*h1 + W9*h0;

  if (btod)
  {
    c9eff += Xd; 
  }

  return c9eff;
}

EvtComplex EvtBtoXsllUtil::GetC9Eff1(double sh, double mbeff,
                                     bool nnlo, bool /*btod*/) 
{
  // This function returns the first-order alpha_s part of C9

  if (!nnlo) return 0.0;
  double logsh;
  logsh = log(sh);
  double mch = 0.29;

  EvtComplex uniti(0.0,1.0);

  EvtComplex c9eff = 0.0;
  if (sh > 0.25)
  { 
    return c9eff;
  }

  // change energy scale to 5.0 for full NNLO calculation below shat = 0.25
  double muscale = 5.0;
  double alphas = 0.215;
  double C1 = -0.487;
  double C2 = 1.024;
  double A8 = -0.148;
  double Lmu = log(muscale/mbeff);

  EvtComplex F91;
  EvtComplex f91;
  EvtComplex k9100(-11.973,0.16371);
  EvtComplex k9101(-0.081271,-0.059691);
  EvtComplex k9110(-28.432,-0.25044);
  EvtComplex k9111(-0.040243,0.016442);
  EvtComplex k9120(-57.114,-0.86486);
  EvtComplex k9121(-0.035191,0.027909);
  EvtComplex k9130(-128.8,-2.5243);
  EvtComplex k9131(-0.017587,0.050639);
  f91 = k9100 + k9101*logsh + sh*(k9110 + k9111*logsh) +
        sh*sh*(k9120 + k9121*logsh) + 
        sh*sh*sh*(k9130 + k9131*logsh); 
  F91 = (-1424.0/729.0 + 16.0*uniti*EvtConst::pi/243.0 
         + 64.0/27.0*log(mch))*Lmu - 16.0*Lmu*logsh/243.0 +
        (16.0/1215.0 - 32.0/135.0/mch/mch)*Lmu*sh +
        (4.0/2835.0 - 8.0/315.0/mch/mch/mch/mch)*Lmu*sh*sh +
    (16.0/76545.0 - 32.0/8505.0/mch/mch/mch/mch/mch/mch)*
    Lmu*sh*sh*sh -256.0*Lmu*Lmu/243.0 + f91;

  EvtComplex F92;
  EvtComplex f92;
  EvtComplex k9200(6.6338,-0.98225);
  EvtComplex k9201(0.48763,0.35815);
  EvtComplex k9210(3.3585,1.5026);
  EvtComplex k9211(0.24146,-0.098649);
  EvtComplex k9220(-1.1906,5.1892);
  EvtComplex k9221(0.21115,-0.16745);
  EvtComplex k9230(-17.12,15.146);
  EvtComplex k9231(0.10552,-0.30383);
  f92 = k9200 + k9201*logsh + sh*(k9210 + k9211*logsh) +
        sh*sh*(k9220 + k9221*logsh) + 
        sh*sh*sh*(k9230 + k9231*logsh); 
  F92 = (256.0/243.0 - 32.0*uniti*EvtConst::pi/81.0 
         - 128.0/9.0*log(mch))*Lmu + 32.0*Lmu*logsh/81.0 +
        (-32.0/405.0 + 64.0/45.0/mch/mch)*Lmu*sh +
        (-8.0/945.0 + 16.0/105.0/mch/mch/mch/mch)*Lmu*sh*sh +
    (-32.0/25515.0 + 64.0/2835.0/mch/mch/mch/mch/mch/mch)*
    Lmu*sh*sh*sh + 512.0*Lmu*Lmu/81.0 + f92;
  
  EvtComplex F98;
  F98 = 104.0/9.0 - 32.0*EvtConst::pi*EvtConst::pi/27.0 + 
        (1184.0/27.0 - 40.0*EvtConst::pi*EvtConst::pi/9.0)*sh +
        (14212.0/135.0 - 32.0*EvtConst::pi*EvtConst::pi/3.0)*sh*sh +
    (193444.0/945.0 - 560.0*EvtConst::pi*EvtConst::pi/27.0)*sh*sh*sh +
        16.0*logsh/9.0*(1.0 + sh + sh*sh + sh*sh*sh);

  c9eff = - alphas/(4.0*EvtConst::pi)*(C1*F91 + C2*F92 + A8*F98);

  return c9eff;
}

EvtComplex EvtBtoXsllUtil::GetC10Eff(double /*sh*/, bool nnlo) 
{

  if (!nnlo) return -4.669;
  double A10;
  A10 = -4.592 + 0.379;

  EvtComplex c10eff;
  c10eff = A10;

  return c10eff;
}

double EvtBtoXsllUtil::dGdsProb(double mb, double ms, double ml,
                                double s)
{
  // Compute the decay probability density function given a value of s
  // according to Ali-Lunghi-Greub-Hiller's 2002 paper
  // Note that the form given below is taken from
  // F.Kruger and L.M.Sehgal, Phys. Lett. B380, 199 (1996)
  // but the differential rate as a function of dilepton mass
  // in this latter paper reduces to Eq.(12) in ALGH's 2002 paper
  // for ml = 0 and ms = 0.

  bool btod = false;
  bool nnlo = true;

  double delta, lambda, prob;
  double f1, f2, f3, f4;
  double msh, mlh, sh;
  double mbeff = 4.8;

  mlh = ml / mb;
  msh = ms / mb;
  // set lepton and strange-quark masses to 0 if need to
  // be in strict agreement with ALGH 2002 paper
  //  mlh = 0.0; msh = 0.0;
  //  sh  = s  / (mb*mb);
  sh  = s  / (mbeff*mbeff);

  // if sh >1.0 code will return a nan. so just skip it
  if ( sh > 1.0 ) return 0.0;


  EvtComplex c7eff0 = EvtBtoXsllUtil::GetC7Eff0(sh,nnlo);
  EvtComplex c7eff1 = EvtBtoXsllUtil::GetC7Eff1(sh,mbeff,nnlo);
  EvtComplex c9eff0 = EvtBtoXsllUtil::GetC9Eff0(sh,mbeff,nnlo,btod);
  EvtComplex c9eff1 = EvtBtoXsllUtil::GetC9Eff1(sh,mbeff,nnlo,btod);
  EvtComplex c10eff = EvtBtoXsllUtil::GetC10Eff(sh,nnlo);

  double alphas = 0.119/
     (1 + 0.119*log(pow(4.8,2)/pow(91.1867,2))*23.0/12.0/EvtConst::pi);

  double omega7 = -8.0/3.0*log(4.8/mb)
                  -4.0/3.0*EvtDiLog::DiLog(sh) 
                  -2.0/9.0*EvtConst::pi*EvtConst::pi
                  -2.0/3.0*log(sh)*log(1.0-sh)
                  -log(1-sh)*(8.0+sh)/(2.0+sh)/3.0 
    -2.0/3.0*sh*(2.0 - 2.0*sh - sh*sh)*log(sh)/pow((1.0 - sh),2)/(2.0 + sh)
    -(16.0 - 11.0*sh - 17.0*sh*sh)/18.0/(2.0 + sh)/(1.0 - sh);
  double eta7 = 1.0 + alphas*omega7/EvtConst::pi;

  double omega79 = -4.0/3.0*log(4.8/mb)
                   -4.0/3.0*EvtDiLog::DiLog(sh) 
                   -2.0/9.0*EvtConst::pi*EvtConst::pi
                   -2.0/3.0*log(sh)*log(1.0-sh)
                   -1.0/9.0*(2.0+7.0*sh)*log(1.0 - sh)/sh
                   -2.0/9.0*sh*(3.0 - 2.0*sh)*log(sh)/pow((1.0 - sh),2) 
                   +1.0/18.0*(5.0 - 9.0*sh)/(1.0 - sh);
  double eta79 = 1.0 + alphas*omega79/EvtConst::pi;

  double omega9 = -2.0/9.0*EvtConst::pi*EvtConst::pi - 4.0/3.0*EvtDiLog::DiLog(sh)
                 - 2.0/3.0*log(sh)*log(1.0-sh)
                 - (5.0+4.0*sh)/(3.0*(1.0+2.0*sh)) * log(1.0-sh)
                 - 2.0*sh*(1.0+sh)*(1.0-2.0*sh)
                 /(3.0*pow(1.0-sh,2)*(1.0+2.0*sh)) * log(sh)
                 + (5.0+9.0*sh-6.0*sh*sh)/(6.0*(1.0-sh)*(1.0+2.0*sh));
  double eta9 = 1.0 + alphas*omega9/EvtConst::pi;

  EvtComplex c7eff = eta7*c7eff0 + c7eff1;
  EvtComplex c9eff = eta9*c9eff0 + c9eff1;
  c10eff *= eta9;

  double c7c7 = abs2(c7eff);
  double c7c9 = real((eta79*c7eff0 + c7eff1)*conj(eta79*c9eff0 + c9eff1));
  double c9c9plusc10c10  = abs2(c9eff) + abs2(c10eff);
  double c9c9minusc10c10 = abs2(c9eff) - abs2(c10eff);

  // Power corrections according to ALGH 2002
  double lambda_1 = -0.2;
  double lambda_2 = 0.12;
  double C1 = -0.487;
  double C2 = 1.024;
  double mc = 0.29 * mb;

  EvtComplex F;
  double r = s / (4.0 * mc * mc);
  EvtComplex uniti(0.0,1.0);
  F = 3.0 / (2.0 * r);
  if (r < 1)
  {
    F *= 1.0/sqrt(r*(1.0-r))*atan(sqrt(r/(1.0-r)))-1.0;
  }
  else
  {
    F *= 0.5/sqrt(r*(r-1.0))*(log((1.0-sqrt(1.0-1.0/r))/(1.0+sqrt(1.0-1.0/r)))
                              +uniti*EvtConst::pi)-1.0;
  }

  double G1 = 1.0 + lambda_1 / (2.0 * mb * mb)
                  + 3.0 * (1.0 - 15.0*sh*sh + 10.0*sh*sh*sh)
                        / ((1.0 - sh)*(1.0 -sh)*(1.0 + 2.0*sh))
                        * lambda_2 / (2.0*mb*mb);
  double G2 = 1.0 + lambda_1 / (2.0 * mb * mb)
                  - 3.0 * (6.0 + 3.0*sh - 5.0*sh*sh*sh)
                        / ((1.0 - sh)*(1.0 -sh)*(2.0 + sh))
                        * lambda_2 / (2.0*mb*mb);
  double G3 = 1.0 + lambda_1 / (2.0 * mb * mb)
                  - (5.0 + 6.0*sh - 7.0*sh*sh)
                     / ((1.0 - sh)*(1.0 -sh))
                     * lambda_2 / (2.0*mb*mb);
  double Gc = -8.0/9.0 * (C2 - C1/6.0) * lambda_2/(mc*mc) 
    * real(F*(conj(c9eff)*(2.0+sh)+conj(c7eff)*(1.0 + 6.0*sh - sh*sh)/sh));

  // end of power corrections section
  // now back to Kruger & Sehgal expressions

  double msh2=msh*msh;
  lambda = 1.0 + sh*sh + msh2*msh2 - 2.0*(sh + sh*msh2 + msh2);
  // negative lambda screw up sqrt below!
  if ( lambda < 0.0 ) return 0.0;

  f1 = pow(1.0-msh2,2) - sh*(1.0 + msh2);
  f2 = 2.0*(1.0 + msh2) * pow(1.0-msh2,2)
       - sh*(1.0 + 14.0*msh2 + pow(msh,4)) - sh*sh*(1.0 + msh2);
  f3 = pow(1.0-msh2,2) + sh*(1.0 + msh2) - 2.0*sh*sh
       + lambda*2.0*mlh*mlh/sh;
  f4 = 1.0 - sh + msh2;

  delta = (  12.0*c7c9*f1*G3 + 4.0*c7c7*f2*G2/sh ) * (1.0 + 2.0*mlh*mlh/sh)
            + c9c9plusc10c10*f3*G1 
            + 6.0*mlh*mlh*c9c9minusc10c10*f4
            + Gc;

  // avoid negative probs
  if ( delta < 0.0 ) delta=0.;
  // negative when sh < 4*mlh*mlh
  //               s < 4*ml*ml
  ///  prob =  sqrt(lambda*(1.0 - 4.0*mlh*mlh/sh)) * delta;
  prob =  sqrt(lambda*(1.0 - 4.0*ml*ml/s)) * delta;

  //   if ( !(prob>=0.0) && !(prob<=0.0) ) {
    //nan
     //     std::cout << lambda << " " << mlh << " " << sh << " " << delta << " " << mb << " " << mbeff << std::endl;
     // std::cout << 4.0*mlh*mlh/sh << " " << 4.0*ml*ml/s <<  " " << s-4.0*ml*ml << " " << ml << std::endl;
     //	std::cout << sh << " " << sh*sh << " " << msh2*msh2 << " " << msh << std::endl;
     //std::cout <<  (  12.0*c7c9*f1*G3 + 4.0*c7c7*f2*G2/sh ) * (1.0 + 2.0*mlh*mlh/sh)
     //	      <<" " << c9c9plusc10c10*f3*G1 
     //	      << " "<< 6.0*mlh*mlh*c9c9minusc10c10*f4
     //	      << " "<< Gc << std::endl;
     //std::cout << C2 << " " << C1 << " "<< lambda_2 << " " << mc <<  " " << real(F*(conj(c9eff)*(2.0+sh)+conj(c7eff)*(1.0 + 6.0*sh - sh*sh)/sh)) << " " << sh << " " << r << std::endl;
     //std::cout << c9eff << " " << eta9 << " " <<c9eff0 << " " <<  c9eff1 << " " << alphas << " " << omega9 << " " << sh << std::endl;

     //}
//  else{
//    if ( sh > 1.0) std::cout << "not a nan \n";
//  }
  return prob;
}

double EvtBtoXsllUtil::dGdsdupProb(double mb, double ms, double ml,
                                   double s,  double u)
{
  // Compute the decay probability density function given a value of s and u
  // according to Ali-Hiller-Handoko-Morozumi's 1997 paper
  // see Appendix E

  bool btod = false;
  bool nnlo = true;

  double prob;
  double f1sp, f2sp, f3sp;
  double mbeff = 4.8;

  //  double sh = s / (mb*mb);
  double sh  = s  / (mbeff*mbeff);

  // if sh >1.0 code will return a nan. so just skip it
  if ( sh > 1.0 ) return 0.0;

  EvtComplex c7eff0 = EvtBtoXsllUtil::GetC7Eff0(sh,nnlo);
  EvtComplex c7eff1 = EvtBtoXsllUtil::GetC7Eff1(sh,mbeff,nnlo);
  EvtComplex c9eff0 = EvtBtoXsllUtil::GetC9Eff0(sh,mbeff,nnlo,btod);
  EvtComplex c9eff1 = EvtBtoXsllUtil::GetC9Eff1(sh,mbeff,nnlo,btod);
  EvtComplex c10eff = EvtBtoXsllUtil::GetC10Eff(sh,nnlo);

  double alphas = 0.119/
     (1 + 0.119*log(pow(4.8,2)/pow(91.1867,2))*23.0/12.0/EvtConst::pi);

  double omega7 = -8.0/3.0*log(4.8/mb)
                  -4.0/3.0*EvtDiLog::DiLog(sh) 
                  -2.0/9.0*EvtConst::pi*EvtConst::pi
                  -2.0/3.0*log(sh)*log(1.0-sh)
                  -log(1-sh)*(8.0+sh)/(2.0+sh)/3.0 
    -2.0/3.0*sh*(2.0 - 2.0*sh - sh*sh)*log(sh)/pow((1.0 - sh),2)/(2.0 + sh)
    -(16.0 - 11.0*sh - 17.0*sh*sh)/18.0/(2.0 + sh)/(1.0 - sh);
  double eta7 = 1.0 + alphas*omega7/EvtConst::pi;

  double omega79 = -4.0/3.0*log(4.8/mb)
                   -4.0/3.0*EvtDiLog::DiLog(sh)
                   -2.0/9.0*EvtConst::pi*EvtConst::pi
                   -2.0/3.0*log(sh)*log(1.0-sh)
                   -1.0/9.0*(2.0+7.0*sh)*log(1.0 - sh)/sh
                   -2.0/9.0*sh*(3.0 - 2.0*sh)*log(sh)/pow((1.0 - sh),2) 
                   +1.0/18.0*(5.0 - 9.0*sh)/(1.0 - sh);
  double eta79 = 1.0 + alphas*omega79/EvtConst::pi;

  double omega9 = - 2.0/9.0*EvtConst::pi*EvtConst::pi - 4.0/3.0*EvtDiLog::DiLog(sh)
                 - 2.0/3.0*log(sh)*log(1.0-sh)
                 - (5.0+4.0*sh)/(3.0*(1.0+2.0*sh)) * log(1.0-sh)
                 - 2.0*sh*(1.0+sh)*(1.0-2.0*sh)
                 /(3.0*pow(1.0-sh,2)*(1.0+2.0*sh)) * log(sh)
                 + (5.0+9.0*sh-6.0*sh*sh)/(6.0*(1.0-sh)*(1.0+2.0*sh));
  double eta9 = 1.0 + alphas*omega9/EvtConst::pi;

  EvtComplex c7eff = eta7*c7eff0 + c7eff1;
  EvtComplex c9eff = eta9*c9eff0 + c9eff1;
  c10eff *= eta9;

  double c7c7  = abs2(c7eff);
  double c7c9  = real((eta79*c7eff0 + c7eff1)*conj(eta79*c9eff0 + c9eff1));
  double c7c10 = real((eta79*c7eff0 + c7eff1)*conj(eta9*c10eff));
  double c9c10 = real((eta9*c9eff0  + c9eff1)*conj(eta9*c10eff));
  double c9c9plusc10c10  = abs2(c9eff) + abs2(c10eff);

  f1sp = ( pow(mb*mb-ms*ms,2) - s*s) * c9c9plusc10c10 
         + 4.0*( pow(mb,4) - ms*ms*mb*mb - pow(ms,4)*(1.0 - ms*ms/(mb*mb))
         - 8.0*s*ms*ms - s*s*(1.0 + ms*ms/(mb*mb) ))*mb*mb*c7c7/s
    // kludged mass term
         *(1.0 + 2.0*ml*ml/s)
         - 8.0*(s*(mb*mb + ms*ms) - pow(mb*mb-ms*ms,2)) * c7c9
    // kludged mass term
         *(1.0 + 2.0*ml*ml/s);

  f2sp = 4.0*s*c9c10 + 8.0*(mb*mb + ms*ms)*c7c10;
  f3sp = - (c9c9plusc10c10)
         + 4.0*(1.0 + pow(ms/mb,4)) * mb*mb*c7c7/s
    // kludged mass term
         *(1.0 + 2.0*ml*ml/s);

  prob = (f1sp + f2sp*u + f3sp*u*u)/ pow(mb,3);
  if ( prob < 0.0 ) prob=0.;

  return prob;
}

double EvtBtoXsllUtil::FermiMomentum(double pf)
{
  // Pick a value for the b-quark Fermi motion momentum
  // according to Ali's Gaussian model

  double pb, pbmax, xbox, ybox;
  pb    = 0.0;
  pbmax = 5.0 * pf;

  while (pb == 0.0)
  {
    xbox = EvtRandom::Flat(pbmax);
    ybox = EvtRandom::Flat();
    if (ybox < FermiMomentumProb(xbox, pf)) { pb = xbox;}
  }

  return pb;
}

double EvtBtoXsllUtil::FermiMomentumProb(double pb, double pf)
{
  // Compute probability according to Ali's Gaussian model
  // the function chosen has a convenient maximum value of 1 for pb = pf

  double prsq = (pb*pb)/(pf*pf);
  double prob = prsq * exp(1.0 - prsq);

  return prob;
}

