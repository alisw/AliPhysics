//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 2003      Caltech
//
// Module: EvtGen/EvtBBScalar
//
// Description:Implementation of the decay B- -> lambda p_bar pi according to
// hep-ph/0204185, hep-ph/0211240
// This model is intended to be applicable to all decays of the type B-> baryon baryon scalar
//
// Modification history:
//
//    Jan Strube     March 24th, 2006         Module created
//
//------------------------------------------------------------------------

#ifndef EVTBBSCALAR_HH
#define EVTBBSCALAR_HH

#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtScalarParticle.hh"
#include "EvtGenBase/EvtDiracParticle.hh"
#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtVector4C.hh"

#include <map>
#include <string>
#include <vector>
#include <bitset>

struct FormFactor {
    double value;
    double sigma1;
    double sigma2;
    double mV;
};

enum Baryons {
    Lambda, Proton, Neutron, Sigma0, Sigma_minus, Xi0, Xi_minus, nBaryons
};


class EvtBBScalar : public EvtDecayAmp {

public:
    EvtBBScalar();
    std::string getName();
    EvtDecayBase* clone();
    void decay(EvtParticle* p);
    void init();
    void initProbMax();

private:
    // used values of constants
    static const EvtComplex I;
    static const EvtComplex V_ub;
    static const EvtComplex V_us_star;
    static const EvtComplex a1;
    static const EvtComplex V_tb;
    static const EvtComplex V_ts_star;
    static const EvtComplex a4;
    static const EvtComplex a6;

    // used parameters in the calculation of the magnetic form factors
    static const double x[];
    static const double y[];
    // quark masses as used in the model
    static const double m_s;
    static const double m_u;
    static const double m_b;

    // used to choose the right value for the form factor depending on the type of scalar
    std::string _scalarType;
    mutable std::map<std::string, FormFactor> _f0Map;
    mutable std::map<std::string, FormFactor> _f1Map;

    // only consider F1+F2 here
    std::bitset<nBaryons> _baryonCombination;
    void setKnownBaryonTypes(const EvtId& baryon);
    
    double B_pi_f1(double t) const ;
    double B_pi_f0(double t) const ;
    double baryonF1F2(double t) const ;
    double G_p(double t) const ;
    double G_n(double t) const ;
    
    double baryon_gA(double t) const;
    double baryon_hA(double t) const;
    double baryon_gP(double t) const ;
    double baryon_fS(double t) const ;

    double D_A(double t) const ;
    double F_A(double t) const ;
    double D_P(double t) const ;
    double F_P(double t) const ;
    double D_S(double t) const ;
    double F_S(double t) const ;

    // (mB1 - mB2)/(mq1 - mq1)
    double _massRatio;
    double _baryonMassSum;
    double formFactorFit(double t, const std::vector<double>& params) const ;

    static const EvtComplex const_B;
    static const EvtComplex const_C;
    const EvtVector4C
    amp_A(const EvtVector4R& p4B, const EvtVector4R& p4Scalar);
    const EvtComplex
    amp_B(const EvtDiracParticle* baryon1, const EvtDiracSpinor& b1Pol
        , const EvtDiracParticle* baryon2, const EvtDiracSpinor& b2Pol
        , int index);
    const EvtComplex
    amp_B_vectorPart(const EvtDiracParticle* baryon1, const EvtDiracSpinor& b1Pol
                   , const EvtDiracParticle* baryon2, const EvtDiracSpinor& b2Pol
                   , int index);
    const EvtComplex
    amp_B_axialPart(const EvtDiracParticle* baryon1, const EvtDiracSpinor& b1Pol
                  , const EvtDiracParticle* baryon2, const EvtDiracSpinor& b2Pol
                  , int index);
    const EvtComplex
    amp_C(const EvtDiracParticle* baryon1, const EvtDiracSpinor& b1Pol
        , const EvtDiracParticle* baryon2, const EvtDiracSpinor& b2Pol
        , int index);
    const EvtComplex
    amp_C_scalarPart(const EvtDiracSpinor& b1Pol, const EvtDiracSpinor& b2Pol, double t);
    const EvtComplex
    amp_C_pseudoscalarPart(const EvtDiracSpinor& b1Pol, const EvtDiracSpinor& b2Pol, double t);

    // initialize phasespace and calculate the amplitude for one (i=0,1) state of the photon
    EvtComplex calcAmpliude(const EvtParticle* p, const unsigned int polState);
};

#endif
