// -*- C++ -*-

#include <TMath.h>
#include <TMatrixD.h>

#include "AliLog.h"
#include "AliDoubleGaussianBeamProfile.h"

/*
  All tables and equations mentioned in the code refer to
  Webb's thesis:  https://cds.cern.ch/record/2020875?ln=en

  Authors: jesus.guillermo.contreras.nuno@cern.ch, lucina.gabriela.espinoza.beltran@cern.ch
*/

ClassImp(AliDoubleGaussianBeamProfile);

namespace Detail {
  // ------------------------------------------------------------------------
  // definition of global variables to store input values
  // ------------------------------------------------------------------------
  // rotation of beams in laboratory frame:
  Double_t alpha_xz; // half corssing angle in rads
  Double_t alpha_yz; // half corssing angle in rads
  Double_t theta_xz; // eq. A.5
  Double_t theta_yz; // eq. A.5
  // parameters for double gaussian of beam 1
  Double_t sigma_xa1; //[mm] sigma in x of gaussian a in beam 1
  Double_t sigma_xb1; //[mm] sigma in x of gaussian b in beam 1
  Double_t sigma_ya1; //[mm] sigma in y of gaussian a in beam 1
  Double_t sigma_yb1; //[mm] sigma in y of gaussian b in beam 1
  Double_t sigma_za1; //[mm] sigma in z of gaussian a in beam 1
  Double_t sigma_zb1; //[mm] sigma in z of gaussian b in beam 1
  Double_t kappa_a1; // x,y correlation of gaussian a in beam 1
  Double_t kappa_b1; // x,y correlation of gaussian b in beam 1
  Double_t w1; // weight of gaussian a in beam 1
  // parameters for double gaussian of beam 2
  Double_t sigma_xa2; //[mm] sigma in x of gaussian a in beam 2
  Double_t sigma_xb2; //[mm] sigma in x of gaussian b in beam 2
  Double_t sigma_ya2; //[mm] sigma in y of gaussian a in beam 2
  Double_t sigma_yb2; //[mm] sigma in y of gaussian b in beam 2
  Double_t sigma_za2; //[mm] sigma in z of gaussian a in beam 2
  Double_t sigma_zb2; //[mm] sigma in z of gaussian b in beam 2
  Double_t kappa_a2; // x,y correlation of gaussian a in beam 2
  Double_t kappa_b2; // x,y correlation of gaussian b in beam 2
  Double_t w2; // weight of gaussian a in beam 2

  Bool_t   doDebug = kFALSE;
  Double_t scaleZ  = 5e-3;

  // ------------------------------------------------------------------------
  // ------------------------------------------------------------------------
  void DefineRotationMatrix(Double_t xz, Double_t yz, TMatrixD& m)
  // arguments are
  // --> the two crossing angles in rad!
  // --> the matrix to be defined
  {
#if 1
    Double_t elements[9] = {TMath::Cos(xz),
			    TMath::Sin(xz)*TMath::Sin(yz),
			    TMath::Sin(xz)*TMath::Cos(yz),
			    0,
			    TMath::Cos(yz),
			    -TMath::Sin(yz),
			    -TMath::Sin(xz),
			    TMath::Sin(yz)*TMath::Cos(xz),
			    TMath::Cos(xz)*TMath::Cos(yz)};
#else
    Double_t elements[9] = {TMath::Cos(yz),
			    0,
			    TMath::Sin(yz),

			    -TMath::Sin(xz)*TMath::Sin(yz),
			    TMath::Cos(xz),
			    -TMath::Sin(xz),

			    TMath::Sin(yz),
			    TMath::Sin(xz)*TMath::Cos(yz),
			    TMath::Cos(xz)*TMath::Cos(yz)};
#endif
    m.SetMatrixArray(elements);
  }

  // ------------------------------------------------------------------------
  // ------------------------------------------------------------------------
  void DefineCovMatrix(Double_t sigma_x, Double_t sigma_y, Double_t sigma_z,
		       Double_t rho_xy, TMatrixD& m)
  // arguments are
  // --> the three Gaussian widths, the correlation in x-y
  // --> the matrix to be defined
  {
    Double_t elements[9] = {sigma_x*sigma_x, sigma_x*sigma_y*rho_xy, 0,
			    sigma_x*sigma_y*rho_xy, sigma_y*sigma_y, 0,
			    0, 0, sigma_z*sigma_z};
    m.SetMatrixArray(elements);
  }

  // ------------------------------------------------------------------------
  // ------------------------------------------------------------------------
  void DefineVector(Double_t x, Double_t y, Double_t z, TMatrixD& m)
  // arguments are
  // --> the three elements of the vector
  // --> the vector (in Matrix notation!) to be defined
  {
    Double_t elements[3] = {x, y, z};
    m.SetMatrixArray(elements);
  }

  // ------------------------------------------------------------------------
  // ------------------------------------------------------------------------
  Double_t GetProduct2Gaussians(const TMatrixD& m1, const TMatrixD& a1, const TMatrixD& C1,
				const TMatrixD& m2, const TMatrixD& a2, const TMatrixD& C2,
				TMatrixD& alpha, TMatrixD& K)
  // compute the amplitude A, mean alpha and Covariance K of the gaussian obtained
  // by the product of two gaussians defined by m, a and C
  // The value of the amplitude is returned.
  {
    // ooooooooooooooooooooooooooooooooooooooooooooooooooo
    // variables directly derived from input matrices
    // ooooooooooooooooooooooooooooooooooooooooooooooooooo
    // --> transposed matrices
    const TMatrixD m1_T(TMatrixD::kTransposed,m1);
    const TMatrixD m2_T(TMatrixD::kTransposed,m2);
    const TMatrixD a1_T(TMatrixD::kTransposed,a1);
    const TMatrixD a2_T(TMatrixD::kTransposed,a2);
    // --> inverse matrices
    const TMatrixD C1_Inv(TMatrixD::kInverted,C1);
    const TMatrixD C2_Inv(TMatrixD::kInverted,C2);
    // --> determinants
    const Double_t C1_Det = C1.Determinant();
    const Double_t C2_Det = C2.Determinant();

    // ooooooooooooooooooooooooooooooooooooooooooooooooooo
    // products used in different equations
    // ooooooooooooooooooooooooooooooooooooooooooooooooooo
    const TMatrixD C1Ia1(C1_Inv,TMatrixD::kMult,a1);
    const TMatrixD C2Ia2(C2_Inv,TMatrixD::kMult,a2);
    const TMatrixD a1TC1I(a1_T,TMatrixD::kMult,C1_Inv);
    const TMatrixD a2TC2I(a2_T,TMatrixD::kMult,C2_Inv);
    const TMatrixD C1Ia1a1TC1I(C1Ia1,TMatrixD::kMult,a1TC1I);
    const TMatrixD C2Ia2a2TC2I(C2Ia2,TMatrixD::kMult,a2TC2I);
    const TMatrixD C1Ia1a2TC2I(C1Ia1,TMatrixD::kMult,a2TC2I);
    const TMatrixD C2Ia2a1TC1I(C2Ia2,TMatrixD::kMult,a1TC1I);
    const TMatrixD C1Ia1a1TC1Im1(C1Ia1a1TC1I,TMatrixD::kMult,m1);
    const TMatrixD C2Ia2a2TC2Im2(C2Ia2a2TC2I,TMatrixD::kMult,m2);
    const TMatrixD C1Ia1a2TC2Im2(C1Ia1a2TC2I,TMatrixD::kMult,m2);
    const TMatrixD C2Ia2a1TC1Im1(C2Ia2a1TC1I,TMatrixD::kMult,m1);
    const TMatrixD a1TC1Ia1(a1_T,TMatrixD::kMult,C1Ia1);
    const TMatrixD a2TC2Ia2(a2_T,TMatrixD::kMult,C2Ia2);
    const TMatrixD C1Im1(C1_Inv,TMatrixD::kMult,m1);
    const TMatrixD C2Im2(C2_Inv,TMatrixD::kMult,m2);
    const TMatrixD m1TC1Im1(m1_T,TMatrixD::kMult,C1Im1);
    const TMatrixD m2TC2Im2(m2_T,TMatrixD::kMult,C2Im2);
    const Double_t sigma_t_2_Inv = a1TC1Ia1(0,0)+a2TC2Ia2(0,0);
    const Double_t sigma_t_2 = 1.0/sigma_t_2_Inv;

    // ooooooooooooooooooooooooooooooooooooooooooooooooooo
    // compute covariance matrix
    // ooooooooooooooooooooooooooooooooooooooooooooooooooo
    TMatrixD K_Inv(C1_Inv,TMatrixD::kPlus,C2_Inv);
    Double_t K_Det = 0;
    K = K_Inv;
    K.InvertFast(&K_Det);

    // ooooooooooooooooooooooooooooooooooooooooooooooooooo
    // compute mean matrix
    // ooooooooooooooooooooooooooooooooooooooooooooooooooo
    TMatrixD m_Sum(C1Im1,TMatrixD::kPlus,C2Im2);
    alpha = K*m_Sum;
    TMatrixD alpha_T(TMatrixD::kTransposed,alpha);

    // ooooooooooooooooooooooooooooooooooooooooooooooooooo
    // compute amplitude
    // ooooooooooooooooooooooooooooooooooooooooooooooooooo
    const TMatrixD beta_Sum1(m1TC1Im1,TMatrixD::kPlus,m2TC2Im2);
    const TMatrixD beta_Sum2(1,1);
    TMatrixD beta_Sum3(1,1);
    beta_Sum3 = alpha_T*(K_Inv*alpha); // alpha_0?
    const Double_t beta = beta_Sum2(0,0)+beta_Sum3(0,0)-beta_Sum1(0,0);
    const Double_t A = TMath::Exp(0.5*beta)*TMath::Sqrt(sigma_t_2*K_Det/(C1_Det*C2_Det))/TMath::TwoPi();
    return A;
  }

  // ------------------------------------------------------------------------
  // ------------------------------------------------------------------------
  Double_t GetProduct2GaussiansWebb(TMatrixD& m1, TMatrixD& a1, TMatrixD& C1,
				    TMatrixD& m2, TMatrixD& a2, TMatrixD& C2,
				    TMatrixD& alpha, TMatrixD& K)
  // compute the amplitude A, mean alpha and Covariance K of the gaussian obtained
  // by the product of two gaussians defined by m, a and C
  // This is the implementation of formulas A.13, A.14 and A.15 of Webb's thesis
  // The value of the amplitude is returned.
  {
    // ooooooooooooooooooooooooooooooooooooooooooooooooooo
    // variables directly derived from input matrices
    // ooooooooooooooooooooooooooooooooooooooooooooooooooo
    // --> transposed matrices
    const TMatrixD m1_T(TMatrixD::kTransposed,m1);
    const TMatrixD m2_T(TMatrixD::kTransposed,m2);
    const TMatrixD a1_T(TMatrixD::kTransposed,a1);
    const TMatrixD a2_T(TMatrixD::kTransposed,a2);
    // --> inverse matrices
    const TMatrixD C1_Inv(TMatrixD::kInverted,C1);
    const TMatrixD C2_Inv(TMatrixD::kInverted,C2);
    // --> determinants
    const Double_t C1_Det = C1.Determinant();
    const Double_t C2_Det = C2.Determinant();

    // ooooooooooooooooooooooooooooooooooooooooooooooooooo
    // products used in different equations
    // ooooooooooooooooooooooooooooooooooooooooooooooooooo
    const TMatrixD C1Ia1(C1_Inv,TMatrixD::kMult,a1);
    const TMatrixD C2Ia2(C2_Inv,TMatrixD::kMult,a2);
    const TMatrixD a1TC1I(a1_T,TMatrixD::kMult,C1_Inv);
    const TMatrixD a2TC2I(a2_T,TMatrixD::kMult,C2_Inv);
    const TMatrixD C1Ia1a1TC1I(C1Ia1,TMatrixD::kMult,a1TC1I);
    const TMatrixD C2Ia2a2TC2I(C2Ia2,TMatrixD::kMult,a2TC2I);
    const TMatrixD C1Ia1a2TC2I(C1Ia1,TMatrixD::kMult,a2TC2I);
    const TMatrixD C2Ia2a1TC1I(C2Ia2,TMatrixD::kMult,a1TC1I);
    const TMatrixD C1Ia1a1TC1Im1(C1Ia1a1TC1I,TMatrixD::kMult,m1);
    const TMatrixD C2Ia2a2TC2Im2(C2Ia2a2TC2I,TMatrixD::kMult,m2);
    const TMatrixD C1Ia1a2TC2Im2(C1Ia1a2TC2I,TMatrixD::kMult,m2);
    const TMatrixD C2Ia2a1TC1Im1(C2Ia2a1TC1I,TMatrixD::kMult,m1);
    const TMatrixD a1TC1Ia1(a1_T,TMatrixD::kMult,C1Ia1);
    const TMatrixD a2TC2Ia2(a2_T,TMatrixD::kMult,C2Ia2);
    const TMatrixD C1Im1(C1_Inv,TMatrixD::kMult,m1);
    const TMatrixD C2Im2(C2_Inv,TMatrixD::kMult,m2);
    const TMatrixD m1TC1Im1(m1_T,TMatrixD::kMult,C1Im1);
    const TMatrixD m2TC2Im2(m2_T,TMatrixD::kMult,C2Im2);
    const Double_t sigma_t_2_Inv = a1TC1Ia1(0,0)+a2TC2Ia2(0,0);
    const Double_t sigma_t_2 = 1.0/sigma_t_2_Inv;

    // ooooooooooooooooooooooooooooooooooooooooooooooooooo
    // compute covariance matrix
    // ooooooooooooooooooooooooooooooooooooooooooooooooooo
    const TMatrixD C1I_plus_C2I(C1_Inv,TMatrixD::kPlus,C2_Inv);
    TMatrixD C_Sum(3,3);
    C_Sum = C1Ia1a1TC1I+C2Ia2a2TC2I+C1Ia1a2TC2I+C2Ia2a1TC1I;
    C_Sum = sigma_t_2*C_Sum;
    const TMatrixD K_Inv(C1I_plus_C2I,TMatrixD::kMinus,C_Sum);
    Double_t K_Det = 0;
    K = K_Inv;
    K.InvertFast(&K_Det);

    // ooooooooooooooooooooooooooooooooooooooooooooooooooo
    // compute mean matrix
    // ooooooooooooooooooooooooooooooooooooooooooooooooooo
    const TMatrixD m_Sum1(C1Im1,TMatrixD::kPlus,C2Im2);
    TMatrixD m_Sum2(3,1);
    m_Sum2 = C1Ia1a1TC1Im1+C2Ia2a2TC2Im2+C1Ia1a2TC2Im2+C2Ia2a1TC1Im1;
    m_Sum2 = sigma_t_2*m_Sum2;
    const TMatrixD m_minus(m_Sum1,TMatrixD::kMinus,m_Sum2);
    alpha = K*m_minus;
    const TMatrixD alpha_T(TMatrixD::kTransposed,alpha);

    // ooooooooooooooooooooooooooooooooooooooooooooooooooo
    // compute amplitude
    // ooooooooooooooooooooooooooooooooooooooooooooooooooo
    const TMatrixD beta_Sum1(m1TC1Im1,TMatrixD::kPlus,m2TC2Im2);
    const TMatrixD m1TC1Ia1a1TC1Im1(m1_T,TMatrixD::kMult,C1Ia1a1TC1Im1);
    const TMatrixD m2TC2Ia2a2TC2Im2(m2_T,TMatrixD::kMult,C2Ia2a2TC2Im2);
    const TMatrixD m1TC1Ia1a2TC2Im2(m1_T,TMatrixD::kMult,C1Ia1a2TC2Im2);
    const TMatrixD m2TC2Ia2a1TC1Im1(m2_T,TMatrixD::kMult,C2Ia2a1TC1Im1);
    TMatrixD beta_Sum2(1,1);
    beta_Sum2 = m1TC1Ia1a1TC1Im1+m2TC2Ia2a2TC2Im2+m1TC1Ia1a2TC2Im2+m2TC2Ia2a1TC1Im1; // eq A.16 corrected
    beta_Sum2 = sigma_t_2*beta_Sum2;
    TMatrixD beta_Sum3(1,1);
    beta_Sum3 = alpha_T*(K_Inv*alpha); // alpha_0?
    const Double_t beta = beta_Sum2(0,0)+beta_Sum3(0,0)-beta_Sum1(0,0);
    const Double_t A = TMath::Exp(0.5*beta)*TMath::Sqrt(sigma_t_2*K_Det/(C1_Det*C2_Det))/TMath::TwoPi();
    return A;
  }

  void ScaleCov(TMatrixD& c, Double_t f) {
    for (Int_t i=0; i<3; ++i) {
      c(2,i) *=f;
      c(i,2) *=f;
    }
  }

  // ------------------------------------------------------------------------
  // ------------------------------------------------------------------------
  void GetBeamSpotProfile(Double_t delta_x, Double_t delta_y, Double_t *profile)
  // for a given separation between beams get the beamspot profile
  // answer is stored in vector 'profile'
  {
    // ooooooooooooooooooooooooooooooooooooooooooooooooooo
    // prepare rotations
    // ooooooooooooooooooooooooooooooooooooooooooooooooooo
    // --> beam 1
    TMatrixD RotBeam1(3,3);
    DefineRotationMatrix(theta_xz,theta_yz,RotBeam1);
    TMatrixD RotBeam1_Inv(3,3);
    RotBeam1_Inv.Transpose(RotBeam1); // using that for rotations, the transpose = the inverse
    // --> beam 2
    TMatrixD RotBeam2(3,3);
    DefineRotationMatrix(-theta_xz,-theta_yz,RotBeam2);
    TMatrixD RotBeam2_Inv(3,3);
    RotBeam2_Inv.Transpose(RotBeam2); // using that for rotations, the transpose = the inverse

    // ooooooooooooooooooooooooooooooooooooooooooooooooooo
    // Define matrices of initial Gaussian distributions
    // Notation:
    //   beams labeled by numbers 1 and 2
    //   gaussians in a beam labeled by letters a and b
    //   gaussian a has weight w, gaussian b has (1-w)
    // ooooooooooooooooooooooooooooooooooooooooooooooooooo
    // --> beam 1: covariances for the two gaussians
    TMatrixD Cov_a1_B(3,3); // Covariance matrix of first gaussian in beam frame
    DefineCovMatrix(sigma_xa1,sigma_ya1,sigma_za1,kappa_a1,Cov_a1_B);
    TMatrixD Cov_a1(3,3); // Covariance matrix of first gaussian in lab frame
    Cov_a1 = RotBeam1_Inv*(Cov_a1_B*RotBeam1);
    TMatrixD Cov_b1_B(3,3); // Covariance matrix of second gaussian in beam frame
    DefineCovMatrix(sigma_xb1,sigma_yb1,sigma_zb1,kappa_b1,Cov_b1_B);
    TMatrixD Cov_b1(3,3); // Covariance matrix of second gaussian in lab frame
    Cov_b1 = RotBeam1_Inv*(Cov_b1_B*RotBeam1);
    // --> beam 1: only one mean
    // note that assumption of beam centered in lab frame is not valid for ALICE.
    TMatrixD mu0_1(3,1); // center of beam 1 at t = 0;
    DefineVector(0.5*delta_x,0.5*delta_y,0,mu0_1);
    TMatrixD a1(3,1); // eq. A4 (there named \vec{a})
#if 1
    DefineVector(-TMath::Sin(theta_xz),
		 TMath::Cos(theta_xz)*TMath::Sin(theta_yz),
		 scaleZ*TMath::Cos(theta_xz)*TMath::Cos(theta_yz),a1);
#else
    DefineVector(-TMath::Cos(theta_xz)*TMath::Sin(theta_yz),
		 -TMath::Sin(theta_xz),
		 TMath::Cos(theta_xz)*TMath::Cos(theta_yz),a1);
#endif
    // --> beam 2: covariances for the two gaussians
    TMatrixD Cov_a2_B(3,3); // Covariance matrix of first gaussian in beam frame
    DefineCovMatrix(sigma_xa2,sigma_ya2,sigma_za2,kappa_a2,Cov_a2_B);
    TMatrixD Cov_a2(3,3); // Covariance matrix of first gaussian in lab frame
    Cov_a2 = RotBeam2_Inv*(Cov_a2_B*RotBeam2);
    TMatrixD Cov_b2_B(3,3); // Covariance matrix of second gaussian in beam frame
    DefineCovMatrix(sigma_xb2,sigma_yb2,sigma_zb2,kappa_b2,Cov_b2_B);
    TMatrixD Cov_b2(3,3); // Covariance matrix of second gaussian in lab frame
    Cov_b2 = RotBeam2_Inv*(Cov_b2_B*RotBeam2);
    // --> beam 1: only one mean
    // note that assumption of beam centered in lab frame is not valid for ALICE.
    TMatrixD mu0_2(3,1); // center of beam 2 at t = 0;
    DefineVector(-0.5*delta_x,-0.5*delta_y,0,mu0_2);
    TMatrixD a2(3,1); // eq. A4 (there named \vec{b})
#if 1
    DefineVector(-TMath::Sin(theta_xz),
		 TMath::Cos(theta_xz)*TMath::Sin(theta_yz),
		 -scaleZ*TMath::Cos(theta_xz)*TMath::Cos(theta_yz),a2);
#else
    DefineVector(-TMath::Cos(-theta_xz)*TMath::Sin(-theta_yz),
		 -TMath::Sin(-theta_xz),
		 TMath::Cos(-theta_xz)*TMath::Cos(-theta_yz),a2);
#endif

    ScaleCov(Cov_a1, scaleZ);
    ScaleCov(Cov_b1, scaleZ);
    ScaleCov(Cov_a2, scaleZ);
    ScaleCov(Cov_b2, scaleZ);

    // ooooooooooooooooooooooooooooooooooooooooooooooooooo
    // Compute the time integrated product of the
    // two gaussians.
    // The product of two double gaussians produces
    // four products of two gaussians in total.
    // Notation: K is the covariance, alpha the mean and A the amplitude
    // ooooooooooooooooooooooooooooooooooooooooooooooooooo
    // first product: gaussian a1 * gaussian a2
    TMatrixD K_a1a2(3,3);
    TMatrixD alpha_a1a2(3,1);

    Double_t A_a1a2 = w1*w2*GetProduct2GaussiansWebb(mu0_1, a1, Cov_a1, mu0_2, a2, Cov_a2, alpha_a1a2, K_a1a2);
    // second product: gaussian a1 * gaussian b2
    TMatrixD K_a1b2(3,3);
    TMatrixD alpha_a1b2(3,1);
    Double_t A_a1b2 = w1*(1-w2)*GetProduct2GaussiansWebb(mu0_1, a1, Cov_a1, mu0_2, a2, Cov_b2, alpha_a1b2, K_a1b2);
    // third product: gaussian b1 * gaussian a2
    TMatrixD K_b1a2(3,3);
    TMatrixD alpha_b1a2(3,1);
    Double_t A_b1a2 = (1-w1)*w2*GetProduct2GaussiansWebb(mu0_1, a1, Cov_b1, mu0_2, a2, Cov_a2, alpha_b1a2, K_b1a2);
    // fourth product: gaussian b1 * gaussian b2
    TMatrixD K_b1b2(3,3);
    TMatrixD alpha_b1b2(3,1);
    Double_t A_b1b2 = (1-w1)*(1-w2)*GetProduct2GaussiansWebb(mu0_1, a1, Cov_b1, mu0_2, a2, Cov_b2, alpha_b1b2, K_b1b2);

    // ooooooooooooooooooooooooooooooooooooooooooooooooooo
    // Compute the luminosity
    // ooooooooooooooooooooooooooooooooooooooooooooooooooo
    const Double_t KF   = 2.0;
    const Double_t SumA = A_a1a2+A_a1b2+A_b1a2+A_b1b2;
    const Double_t Lumi = KF*SumA;

    // ooooooooooooooooooooooooooooooooooooooooooooooooooo
    // Compute the vector of mean positions
    // ooooooooooooooooooooooooooooooooooooooooooooooooooo
    TMatrixD mean_pos(3,1);
    mean_pos = (A_a1a2*alpha_a1a2)+(A_a1b2*alpha_a1b2)
      +(A_b1a2*alpha_b1a2)+(A_b1b2*alpha_b1b2);

    TMatrixD mean_pos_norm(3,1);
    mean_pos_norm = (1.0/SumA)*mean_pos;
    if (doDebug)
      Printf("A_mean_pos: %e | %e %e %e %e (%f %f)",
	     mean_pos_norm(0,0),
	     (A_a1a2*alpha_a1a2)(0,0),
	     (A_a1b2*alpha_a1b2)(0,0),
	     (A_b1a2*alpha_b1a2)(0,0),
	     (A_b1b2*alpha_b1b2)(0,0),
	     delta_x, delta_y);

    // ooooooooooooooooooooooooooooooooooooooooooooooooooo
    // Compute matrix of widths
    // ooooooooooooooooooooooooooooooooooooooooooooooooooo
    const TMatrixD p1(alpha_a1a2,TMatrixD::kMultTranspose,alpha_a1a2);
    const TMatrixD p2(alpha_a1b2,TMatrixD::kMultTranspose,alpha_a1b2);
    const TMatrixD p3(alpha_b1a2,TMatrixD::kMultTranspose,alpha_b1a2);
    const TMatrixD p4(alpha_b1b2,TMatrixD::kMultTranspose,alpha_b1b2);
    TMatrixD Cov1(K_a1a2,TMatrixD::kPlus,p1);
    TMatrixD Cov2(K_a1b2,TMatrixD::kPlus,p2);
    TMatrixD Cov3(K_b1a2,TMatrixD::kPlus,p3);
    TMatrixD Cov4(K_b1b2,TMatrixD::kPlus,p4);
    if (doDebug) {
      Printf("p1:"); p1.Print();
      Printf("p2:"); p2.Print();
      Printf("p3:"); p3.Print();
      Printf("p4:"); p4.Print();
      Printf("Cov1:"); Cov1.Print();
      Printf("Cov2:"); Cov2.Print();
      Printf("Cov3:"); Cov3.Print();
      Printf("Cov4:"); Cov4.Print();
      Printf("#################### A_a1a2,...=%e %e %e %e", A_a1a2, A_a1b2, A_b1a2, A_b1b2);
    }
    Cov1 = A_a1a2*Cov1;
    Cov2 = A_a1b2*Cov2;
    Cov3 = A_b1a2*Cov3;
    Cov4 = A_b1b2*Cov4;
    TMatrixD Cov_Sum1(3,3);
    Cov_Sum1 = Cov1+Cov2+Cov3+Cov4;
    Cov_Sum1 = (1.0/SumA)*Cov_Sum1;

    if (doDebug) {
      Printf("Cov1:"); Cov1.Print();
      Printf("Cov2:"); Cov2.Print();
      Printf("Cov3:"); Cov3.Print();
      Printf("Cov4:"); Cov4.Print();
    }

    TMatrixD Cov_Sum2(mean_pos,TMatrixD::kMultTranspose,mean_pos);
    Cov_Sum2 = (1.0/(SumA*SumA))*Cov_Sum2;
    const TMatrixD Cov(Cov_Sum1,TMatrixD::kMinus,Cov_Sum2);
    const Double_t Sx = TMath::Sqrt(Cov(0,0));
    const Double_t Sy = TMath::Sqrt(Cov(1,1));
    const Double_t Sz = TMath::Sqrt(Cov(2,2));
    const Double_t Sxy = Cov(0,1)/(Sx*Sy);

    // ooooooooooooooooooooooooooooooooooooooooooooooooooo
    // Fill values in output array
    // ooooooooooooooooooooooooooooooooooooooooooooooooooo
    profile[0] = Lumi * scaleZ*scaleZ*scaleZ;
    profile[1] = mean_pos_norm(0,0);
    profile[2] = mean_pos_norm(1,0);
    profile[3] = mean_pos_norm(2,0) / scaleZ;
    //  mean_pos_norm.Print();
    profile[4] = Sx;
    profile[5] = Sy;
    profile[6] = Sz / scaleZ;
    profile[7] = Sxy;
  }

}

Bool_t AliDoubleGaussianBeamProfile::Eval(Double_t sepX, Double_t sepY, const TVectorD &par, TVectorD &profile, Double_t scaleZ, Bool_t debug)
{
  if (par.GetNoElements()      < 20) AliFatalClass("par.GetNoElements() < 20");
  if (profile.GetNoElements() !=  8) AliFatalClass("profile.GetNoElements() !=  8");

  Detail::scaleZ  = scaleZ;
  Detail::doDebug = debug;

  Detail::sigma_xa1 = par[ 0];
  Detail::sigma_ya1 = par[ 1];
  Detail::sigma_za1 = par[ 2];
  Detail::kappa_a1  = par[ 3];

  Detail::sigma_xb1 = par[ 4]*par[ 0];
  Detail::sigma_yb1 = par[ 5]*par[ 1];
  Detail::sigma_zb1 = par[ 6]*par[ 2];
  Detail::kappa_b1  = par[ 7];

  Detail::w1        = par[ 8];

  Detail::sigma_xa2 = par[ 9];
  Detail::sigma_ya2 = par[10];
  Detail::sigma_za2 = par[11];
  Detail::kappa_a2  = par[12];

  Detail::sigma_xb2 = par[13]*par[ 9];
  Detail::sigma_yb2 = par[14]*par[10];
  Detail::sigma_zb2 = par[15]*par[11];
  Detail::kappa_b2  = par[16];

  Detail::w2        = par[17];

  Detail::alpha_xz  = par[18];
  Detail::alpha_yz  = par[19];

  Detail::theta_xz  = Detail::alpha_xz; // eq. A.5
  Detail::theta_yz  = TMath::ATan(-TMath::Cos(Detail::theta_xz)*TMath::Tan(Detail::alpha_yz));

  Detail::GetBeamSpotProfile(sepX, sepY, profile.GetMatrixArray());

  // conversion to the ALICE coordinate system
  profile[1] *= -1;
  profile[3] *= -1;
  profile[7] *= -1;

  return (!TMath::IsNaN(profile[1]) &&
	  !TMath::IsNaN(profile[2]) &&
	  !TMath::IsNaN(profile[3]) &&
	  !TMath::IsNaN(profile[4]) &&
	  !TMath::IsNaN(profile[5]) &&
	  !TMath::IsNaN(profile[6]) &&
	  !TMath::IsNaN(profile[7]));
}

Double_t AliDoubleGaussianBeamProfile::EvalProfile0(Double_t *x, Double_t *p) {
  const TVectorD par(20, p);
  TVectorD profile(8);
  Eval(x[0], x[1], par, profile, 5e-3);
  return (TMath::IsNaN(profile(0)) ? 0.0 : profile(0));
}

const char* AliDoubleGaussianBeamProfile::GetParName(Int_t i)
{
  static const char* parNames[] = {
    "#sigma^{X}_{1a}",
    "#sigma^{Y}_{1a}",
    "#sigma^{Z}_{1a}",
    "#rho^{XY}_{1a}",
    "S^{X}_{1b}",
    "S^{Y}_{1b}",
    "S^{Z}_{1b}",
    "#rho^{XY}_{1b}",
    "w_{1}",
    "#sigma^{X}_{2a}",
    "#sigma^{Y}_{2a}",
    "#sigma^{Z}_{2a}",
    "#rho^{XY}_{2a}",
    "S^{X}_{2b}",
    "S^{Y}_{2b}",
    "S^{Z}_{2b}",
    "#rho^{XY}_{2b}",
    "w_{2}",
    "#alpha^{XZ}",
    "#alpha^{YZ}"
  };
  const Int_t n = sizeof(parNames)/sizeof(const char*);
  if (i >= n) AliFatalClassF("%d > %d", i, n);
  return parNames[i];
}
