// SusyCouplings.cc is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// Main authors of this file: N. Desai, P. Skands
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the
// supersymmetric couplings class.

#include "Pythia8/ParticleData.h"
#include "Pythia8/SusyCouplings.h"

namespace Pythia8 {

//==========================================================================

// The CoupSUSY class.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Allow verbose printout for debug purposes.
  const bool CoupSUSY::DBSUSY = false;

//--------------------------------------------------------------------------

// Initialize SM+SUSY couplings (only performed once).

void CoupSUSY::initSUSY (SusyLesHouches* slhaPtrIn, Info* infoPtrIn,
    ParticleData* particleDataPtrIn, Settings* settingsPtrIn) {

  // Save pointers.
  infoPtr         = infoPtrIn;
  slhaPtr         = slhaPtrIn;
  settingsPtr     = settingsPtrIn;
  particleDataPtr = particleDataPtrIn;

  // Only initialize SUSY parts if SUSY is actually switched on
  if (!slhaPtr->modsel.exists()) return;

  // Is NMSSM switched on?
  isNMSSM = (slhaPtr->modsel(3) != 1 ? false : true);
  settingsPtr->flag("SLHA:NMSSM",isNMSSM);
  int nNeut = (isNMSSM ? 5 : 4);
  int nChar = 2;

  // Initialize pole masses
  mZpole    = particleDataPtr->m0(23);
  wZpole    = particleDataPtr->mWidth(23);
  mWpole    = particleDataPtr->m0(24);
  wWpole    = particleDataPtr->mWidth(24);

  // Running masses and weak mixing angle
  // (default to pole values if no running available)
  mW        = mWpole;
  mZ        = mZpole;
  sin2W     = sin2thetaW();

  if (settingsPtr->mode("SUSY:sin2thetaWMode") == 3) {
    // Possibility to force on-shell sin2W definition, mostly intended for
    // cross-checks
    sin2W     = 1.0 - pow(mW/mZ,2);
  } else if (settingsPtr->mode("SUSY:sin2thetaWMode") == 2){
    // Possibility to use running sin2W definition, derived from gauge
    // couplings in running SLHA blocks (normally defined at SUSY scale)
    if(slhaPtr->gauge.exists(1) && slhaPtr->gauge.exists(2)
       && slhaPtr->hmix.exists(3)) {
      double gp=slhaPtr->gauge(1);
      double g =slhaPtr->gauge(2);
      double v =slhaPtr->hmix(3);
      mW      = g * v / 2.0;
      mZ      = sqrt(pow(gp,2)+pow(g,2)) * v / 2.0;
      //double tan2W   = pow2(gp)/pow2(g);
      //if (DBSUSY) cout << " tan2W = " << tan2W << endl;
      sin2W   = pow2(gp)/(pow2(g)+pow2(gp));
    } else {
      infoPtr->errorMsg("Warning in CoupSUSY::initSUSY: Block GAUGE"
        " not found or incomplete; using sin(thetaW) at mZ");
    }
  }

  // Overload the SM value of sin2thetaW
  s2tW = sin2W;

  sinW = sqrt(sin2W);
  cosW = sqrt(1.0-sin2W);

  // Tan(beta)
  // By default, use the running one in HMIX (if not found, use MINPAR)

  if(slhaPtr->hmix.exists(2))
    tanb = slhaPtr->hmix(2);
  else{
    tanb = slhaPtr->minpar(3);
    infoPtr->errorMsg("Warning in CoupSUSY::initSUSY: Block HMIX"
      " not found or incomplete; using MINPAR tan(beta)");
  }
  cosb = sqrt( 1.0 / (1.0 + tanb*tanb) );
  sinb = sqrt( max(0.0, 1.0 - cosb*cosb));

  double beta = acos(cosb);

  // Verbose output
  if (DBSUSY) {
    cout << " sin2W(Q) = " << sin2W << "  mW(Q) = " << mW
         << "  mZ(Q) = " << mZ << endl;
    cout << " vev(Q) = " << slhaPtr->hmix(3) << " tanb(Q) = " << tanb
         << endl;
    for (int i=1;i<=3;i++) {
      for (int j=1;j<=3;j++) {
        cout << " VCKM  [" << i << "][" << j << "] = "
             << scientific << setw(10) << VCKMgen(i,j) << endl;
      }
    }
  }

  // Higgs sector
  if(slhaPtr->alpha.exists()) {
    alphaHiggs = slhaPtr->alpha();
  }
  // If RPV, assume alpha = asin(RVHMIX(1,2)) (ignore Higgs-Sneutrino mixing)
  else if (slhaPtr->modsel(4) == 1) {
    alphaHiggs = asin(slhaPtr->rvhmix(1,2));
    infoPtr->errorMsg("Info from CoupSUSY::initSUSY:","Extracting angle"
      " alpha from RVHMIX", true);
  } else {
    infoPtr->errorMsg("Info from CoupSUSY::initSUSY:","Block ALPHA"
      " not found; using alpha = beta.", true);
    // Define approximate alpha by simple SM limit
    alphaHiggs = atan(tanb);
  }

  if(slhaPtr->hmix.exists(1) && slhaPtr->hmix.exists(4)){
    muHiggs = slhaPtr->hmix(1);
    mAHiggs = sqrt(slhaPtr->hmix(4));
  } else if (slhaPtr->rvamix.exists()){
    mAHiggs = particleDataPtr->m0(36);
    muHiggs = 0.0;
    infoPtr->errorMsg("Warning from CoupSUSY::initSUSY: Block HMIX"
      " not found or incomplete","; setting mu = 0.");
  } else{
    infoPtr->errorMsg("Warning from CoupSUSY::initSUSY: Block HMIX"
      " not found or incomplete","; setting mu = mA = 0.");
    muHiggs = 0.0;
    mAHiggs = 0.0;
  }

  // Pass SLHA input to 2HDM sector

  double sba = sin(beta-alphaHiggs);
  double cba = cos(beta-alphaHiggs);
  double cosa = cos(alphaHiggs);
  double sina = sin(alphaHiggs);

  // h0
  settingsPtr->parm("HiggsH1:coup2d", -sina/cosb);
  settingsPtr->parm("HiggsH1:coup2u",  cosa/sinb);
  settingsPtr->parm("HiggsH1:coup2l", cosa/sinb);
  settingsPtr->parm("HiggsH1:coup2Z", sba);
  settingsPtr->parm("HiggsH1:coup2W", sba);
  // H0
  settingsPtr->parm("HiggsH2:coup2d", cosa/cosb);
  settingsPtr->parm("HiggsH2:coup2u", sina/sinb);
  settingsPtr->parm("HiggsH2:coup2l", sina/sinb);
  settingsPtr->parm("HiggsH2:coup2Z", cba);
  settingsPtr->parm("HiggsH2:coup2W", cba);
  settingsPtr->parm("HiggsH2:coup2H1Z", 0.0);
  settingsPtr->parm("HiggsH2:coup2HchgW", sba);
  // A0
  settingsPtr->parm("HiggsA3:coup2d", tanb);
  settingsPtr->parm("HiggsA3:coup2u", cosb/sinb);
  settingsPtr->parm("HiggsA3:coup2l", cosb/sinb);
  settingsPtr->parm("HiggsA3:coup2Z", 0.0);
  settingsPtr->parm("HiggsA3:coup2W", 0.0);
  settingsPtr->parm("HiggsA3:coup2H1Z", cba);
  settingsPtr->parm("HiggsA3:coup2H2Z", sba);
  settingsPtr->parm("HiggsA3:coup2HchgW", 1.0);

  // H^+
  settingsPtr->parm("HiggsHchg:tanBeta", tanb);
  settingsPtr->parm("HiggsHchg:coup2H1W", cba);
  settingsPtr->parm("HiggsHchg:coup2H2W", sba);

  // Triple higgs couplings

  double cbpa = cos(beta+alphaHiggs);
  double sbpa = sin(beta+alphaHiggs);

  settingsPtr->parm("HiggsH1:coup2Hchg", cos(2*beta)*sbpa + 2*pow2(cosW)*sba);
  settingsPtr->parm("HiggsH2:coup2Hchg", -cos(2*beta)*cbpa + 2*pow2(cosW)*cba);
  settingsPtr->parm("HiggsH2:coup2H1H1", 2*sin(2*alphaHiggs)*sbpa
    - cos(2*alphaHiggs)*cbpa);
  settingsPtr->parm("HiggsH2:coup2A3A3", -cos(2*beta)*cbpa);
  settingsPtr->parm("HiggsH2:coup2A3H1", 0.0);
  settingsPtr->parm("HiggsA3:coup2H1H1", 0.0);
  settingsPtr->parm("HiggsA3:coup2Hchg", 0.0);

  // Shorthand for squark mixing matrices
  LHmatrixBlock<6> Ru(slhaPtr->usqmix);
  LHmatrixBlock<6> Rd(slhaPtr->dsqmix);
  LHmatrixBlock<6> imRu(slhaPtr->imusqmix);
  LHmatrixBlock<6> imRd(slhaPtr->imdsqmix);

  // Construct ~g couplings
  for (int i=1; i<=6; i++) {
    for (int j=1; j<=3; j++) {
      LsddG[i][j] = complex( Rd(i,j)  ,  imRd(i,j));
      RsddG[i][j] = complex(-Rd(i,j+3), -imRd(i,j+3));
      LsuuG[i][j] = complex( Ru(i,j)  ,  imRu(i,j));
      RsuuG[i][j] = complex(-Ru(i,j+3), -imRu(i,j+3));

      if (DBSUSY) {
        cout << " Lsddg  [" << i << "][" << j << "] = "
             << scientific << setw(10) << LsddG[i][j]
             << " RsddG  [" << i << "][" << j  << "] = "
             << scientific << setw(10) << RsddG[i][j] << endl;
        cout << " Lsuug  [" << i << "][" << j << "] = "
             << scientific << setw(10) << LsuuG[i][j]
             << " RsuuG  [" << i << "][" << j  << "] = "
             << scientific << setw(10) << RsuuG[i][j] << endl;
      }
    }
  }

  // Construct qqZ couplings
  for (int i=1 ; i<=6 ; i++) {

    // q[i] q[i] Z (def with extra factor 2 compared to [Okun])
    LqqZ[i] = af(i) - 2.0*ef(i)*sin2W ;
    RqqZ[i] =       - 2.0*ef(i)*sin2W ;

    // tmp: verbose output
    if (DBSUSY) {
      cout << " LqqZ  [" << i << "][" << i << "] = "
           << scientific << setw(10) << LqqZ[i]
           << " RqqZ  [" << i << "][" << i  << "] = "
           << scientific << setw(10) << RqqZ[i] << endl;
    }
  }

  // Construct ~q~qZ couplings
  for (int i=1 ; i<=6 ; i++) {

    // Squarks can have off-diagonal couplings as well
    for (int j=1 ; j<=6 ; j++) {

      // ~d[i] ~d[j] Z
      LsdsdZ[i][j] = 0.0;
      RsdsdZ[i][j] = 0.0;
      for (int k=1;k<=3;k++) {
        complex Rdik  = complex(Rd(i,k),  imRd(i,k)  );
        complex Rdjk  = complex(Rd(j,k),  imRd(j,k)  );
        complex Rdik3 = complex(Rd(i,k+3),imRd(i,k+3));
        complex Rdjk3 = complex(Rd(j,k+3),imRd(j,k+3));
        LsdsdZ[i][j] += LqqZ[1] * (Rdik*conj(Rdjk));
        RsdsdZ[i][j] += RqqZ[1] * (Rdik3*conj(Rdjk3));
      }

      // ~u[i] ~u[j] Z
      LsusuZ[i][j] = 0.0;
      RsusuZ[i][j] = 0.0;
      for (int k=1;k<=3;k++) {
        complex Ruik  = complex(Ru(i,k)  ,imRu(i,k)  );
        complex Rujk  = complex(Ru(j,k)  ,imRu(j,k)  );
        complex Ruik3 = complex(Ru(i,k+3),imRu(i,k+3));
        complex Rujk3 = complex(Ru(j,k+3),imRu(j,k+3));
        LsusuZ[i][j] += LqqZ[2] * (Ruik*conj(Rujk));
        RsusuZ[i][j] += RqqZ[2] * (Ruik3*conj(Rujk3));
      }

      // tmp: verbose output
      if (DBSUSY) {
        if (max(abs(LsdsdZ[i][j]),abs(RsdsdZ[i][j])) > 1e-6) {
          cout << " LsdsdZ[" << i << "][" << j << "] = "
               << scientific << setw(10) << LsdsdZ[i][j]
               << " RsdsdZ[" << i << "][" << j << "] = "
               << scientific << setw(10) << RsdsdZ[i][j] << endl;
        }
        if (max(abs(LsusuZ[i][j]),abs(RsusuZ[i][j]))> 1e-6) {
          cout << " LsusuZ[" << i << "][" << j << "] = "
               << scientific << setw(10) << LsusuZ[i][j]
               << " RsusuZ[" << i << "][" << j << "] = "
               << scientific << setw(10) << RsusuZ[i][j] << endl;
        }
      }
    }
  }

  for(int i = 1; i < 7; i++)
    for(int j = 1; j < 7; j++){
      Rsl[i][j] = slhaPtr->selmix(i,j);
    }


  for(int i = 1; i < 7; i++)
    for(int j = 1; j < 7; j++){
      if (i < 4 && j < 4) Rsv[i][j] = slhaPtr->snumix(i,j);
      else Rsv[i][j] = 0.0;
    }

  // In RPV, the slepton mixing matrices include Higgs bosons
  // Here we just extract the entries corresponding to the slepton PDG codes
  // I.e., slepton-Higgs mixing is not supported.

  // Charged sleptons
  if (slhaPtr->modsel(4) >= 1 && slhaPtr->rvlmix.exists()) {
    for (int i=1; i<=6; ++i)
      for (int j=1; j<=6; ++j)
        Rsl[i][j] = slhaPtr->rvlmix(i+1,j+2);
    // Check for Higgs-slepton mixing in RVLMIX
    bool hasCrossTerms = false;
    for (int i=2; i<=7; ++i)
      for (int j=1; j<=2; ++j)
        if (abs(slhaPtr->rvlmix(i,j)) > 1e-5) {
          hasCrossTerms = true;
          break;
        }
    if (hasCrossTerms)
      infoPtr->errorMsg("Warning from CoupSUSY::initSUSY: "
        "slepton-Higgs mixing not supported internally in PYTHIA");
  }

  // Neutral sleptons
  if (slhaPtr->modsel(4) >= 1 && slhaPtr->rvhmix.exists()) {
    for (int i=1; i<=3; ++i)
      for (int j=1; j<=3; ++j)
        Rsv[i][j] = slhaPtr->rvhmix(i+2,j+2);
    // Check for Higgs-sneutrino mixing in RVHMIX
    bool hasCrossTerms = false;
    for (int i=3; i<=7; ++i)
      for (int j=1; j<=2; ++j)
        if (abs(slhaPtr->rvhmix(i,j)) > 1e-5) {
          hasCrossTerms = true;
          break;
        }
    if (hasCrossTerms)
      infoPtr->errorMsg("Warning from CoupSUSY::initSUSY: "
        "sneutrino-Higgs mixing not supported internally in PYTHIA");
  }

  if(DBSUSY){
    cout<<"Rsl"<<endl;
    for(int i=1;i<=6;i++){
      for(int j=1;j<=6;j++){
        cout << scientific << setw(10) << Rsl[i][j]<<"  ";
      }
      cout<<endl;
    }
    cout<<"Rsv"<<endl;
    for(int i=1;i<=6;i++){
      for(int j=1;j<=6;j++){
        cout << scientific << setw(10) << Rsv[i][j]<<"  ";
      }
      cout<<endl;
    }
  }


  // Construct llZ couplings;
  for (int i=11 ; i<=16 ; i++) {

    LllZ[i-10] = af(i) - 2.0*ef(i)*sin2W ;
    RllZ[i-10] =       - 2.0*ef(i)*sin2W ;

    // tmp: verbose output
    if (DBSUSY) {
      cout << " LllZ  [" << i-10 << "][" << i-10 << "] = "
           << scientific << setw(10) << LllZ[i-10]
           << " RllZ  [" << i-10 << "][" << i-10  << "] = "
           << scientific << setw(10) << RllZ[i-10] << endl;
    }
  }

  // Construct ~l~lZ couplings
  // Initialize
  for(int i=1;i<=6;i++){
    for(int j=1;j<=6;j++){
      LslslZ[i][j] = 0.0;
      RslslZ[i][j] = 0.0;

      for (int k=1;k<=3;k++) {
        LslslZ[i][j] += LllZ[1] * Rsl[i][k] * Rsl[j][k];
        RslslZ[i][j] += RllZ[1] * Rsl[i][k+3] * Rsl[j][k+3];
      }

      // ~v[i] ~v[j] Z
      LsvsvZ[i][j] = 0.0;
      RsvsvZ[i][j] = 0.0;
      for (int k=1;k<=3;k++)
        LsvsvZ[i][j] += LllZ[2] * Rsv[i][k]*Rsv[j][k];
    }
  }

  for(int i=1;i<=6;i++){
    for(int j=1;j<=6;j++){
      if (DBSUSY) {
        if (max(abs(LsvsvZ[i][j]),abs(RsvsvZ[i][j])) > 1e-6) {
          cout << " LsvsvZ[" << i << "][" << j << "] = "
               << scientific << setw(10) << LsvsvZ[i][j]
               << " RsvsvZ[" << i << "][" << j << "] = "
               << scientific << setw(10) << RsvsvZ[i][j] << endl;
        }
        if (max(abs(LslslZ[i][j]),abs(RslslZ[i][j]))> 1e-6) {
          cout << " LslslZ[" << i << "][" << j << "] = "
               << scientific << setw(10) << LslslZ[i][j]
               << " RslslZ[" << i << "][" << j << "] = "
               << scientific << setw(10) << RslslZ[i][j] << endl;
        }
      }
    }
  }

  // Construct udW couplings
  // Loop over up [i] and down [j] quark generation
  for (int i=1;i<=3;i++) {
    for (int j=1;j<=3;j++) {

      // CKM matrix (use Pythia one if no SLHA)
      // (NB: could also try input one if no running one found, but
      // would then need to compute from Wolfenstein)
      complex Vij=VCKMgen(i,j);
      if (slhaPtr->vckm.exists()) {
        Vij=complex(slhaPtr->vckm(i,j),slhaPtr->imvckm(i,j));
      }

      // u[i] d[j] W
      LudW[i][j] = sqrt(2.0) * cosW * Vij;
      RudW[i][j] = 0.0;

      // tmp: verbose output
      if (DBSUSY) {
        cout << " LudW  [" << i << "][" << j << "] = "
             << scientific << setw(10) << LudW[i][j]
             << " RudW  [" << i << "][" << j << "] = "
             << scientific << setw(10) << RudW[i][j] << endl;
      }
    }
  }


  // Construct ~u~dW couplings
  // Loop over ~u[k] and ~d[l] flavours
  for (int k=1;k<=6;k++) {
    for (int l=1;l<=6;l++) {

      LsusdW[k][l]=0.0;
      RsusdW[k][l]=0.0;

      // Loop over u[i] and d[j] flavours
      for (int i=1;i<=3;i++) {
        for (int j=1;j<=3;j++) {

          // CKM matrix (use Pythia one if no SLHA)
          // (NB: could also try input one if no running one found, but
          // would then need to compute from Wolfenstein)
          complex Vij=VCKMgen(i,j);
          if (slhaPtr->vckm.exists()) {
            Vij=complex(slhaPtr->vckm(i,j),slhaPtr->imvckm(i,j));
          }

          // ~u[k] ~d[l] W (add one term for each quark flavour i,j)
          complex Ruki = complex(Ru(k,i),imRu(k,i));
          complex Rdlj = complex(Rd(l,j),imRd(l,j));
          LsusdW[k][l] += sqrt(2.0) * cosW * Vij * Ruki * conj(Rdlj);
          RsusdW[k][l] += 0.0;

        }
      }

      // tmp: verbose output
      if (DBSUSY) {
        if (max(abs(LsusdW[k][l]),abs(RsusdW[k][l]))> 1e-6) {
          cout << " LsusdW[" << k << "][" << l << "] = "
               << scientific << setw(10) << LsusdW[k][l]
               << " RsusdW[" << k << "][" << l << "] = "
               << scientific << setw(10) << RsusdW[k][l] << endl;
        }
      }

    }
  }



  // Construct lvW couplings
  for (int i=1;i<=3;i++){
    for (int j=1;j<=3;++j){
       LlvW[i][j] = (i==j) ? sqrt(2.0) * cosW : 0.0 ;
       RlvW[i][j] = 0.0;

      // tmp: verbose output
      if (DBSUSY) {
        cout << " LlvW  [" << i << "][" << j << "] = "
             << scientific << setw(10) << LlvW[i][j]
             << " RlvW  [" << i << "][" << j << "] = "
             << scientific << setw(10) << RlvW[i][j] << endl;
      }
    }
  }

  // Construct ~l~vW couplings
  for (int k=1;k<=6;k++) {
    for (int l=1;l<=6;l++) {
      LslsvW[k][l]=0.0;
      RslsvW[k][l]=0.0;

      if(l<=3) // Only left-handed sneutrinos
        for(int i=1;i<=3;i++)
          LslsvW[k][l] += sqrt(2.0) * cosW * Rsl[k][i] * Rsv[l][i];


      // tmp: verbose output
      if (DBSUSY) {
        cout << " LslsvW  [" << k << "][" << l << "] = "
             << scientific << setw(10) << LslsvW[k][l]
             << " RslsvW  [" << k << "][" << l << "] = "
             << scientific << setw(10) << RslsvW[k][l] << endl;
      }
    }
  }

  // Now we come to the ones with really many indices

  // RPV: check for lepton-neutralino mixing (not supported)
  if (slhaPtr->modsel(4) >= 1 && slhaPtr->rvnmix.exists()) {
    bool hasCrossTerms = false;
    for (int i=1; i<=3; ++i)
      for (int j=4; j<=7; ++j)
        if (abs(slhaPtr->rvnmix(i,j)) > 1e-5) {
          hasCrossTerms = true;
          break;
        }
    if (hasCrossTerms)
      infoPtr->errorMsg("Warning from CoupSUSY::initSUSY: "
        "Neutrino-Neutralino mixing not supported internally in PYTHIA");
  }

  // Construct ~chi0 couplings (allow for 5 neutralinos in NMSSM)
  for (int i=1;i<=nNeut;i++) {

    // Ni1, Ni2, Ni3, Ni4, Ni5
    complex ni1,ni2,ni3,ni4,ni5;

    // In RPV, ignore neutrino-neutralino mixing
    if (slhaPtr->modsel(4) >= 1 && slhaPtr->rvnmix.exists()) {
      ni1=complex( slhaPtr->rvnmix(i+3,4), 0.0 );
      ni2=complex( slhaPtr->rvnmix(i+3,5), 0.0 );
      ni3=complex( slhaPtr->rvnmix(i+3,6), 0.0 );
      ni4=complex( slhaPtr->rvnmix(i+3,7), 0.0 );
      ni5=complex( 0.0, 0.0);
    }
   else if (!isNMSSM) {
      ni1=complex( slhaPtr->nmix(i,1), slhaPtr->imnmix(i,1) );
      ni2=complex( slhaPtr->nmix(i,2), slhaPtr->imnmix(i,2) );
      ni3=complex( slhaPtr->nmix(i,3), slhaPtr->imnmix(i,3) );
      ni4=complex( slhaPtr->nmix(i,4), slhaPtr->imnmix(i,4) );
      ni5=complex( 0.0, 0.0);
    } else {
      ni1=complex( slhaPtr->nmnmix(i,1), slhaPtr->imnmnmix(i,1) );
      ni2=complex( slhaPtr->nmnmix(i,2), slhaPtr->imnmnmix(i,2) );
      ni3=complex( slhaPtr->nmnmix(i,3), slhaPtr->imnmnmix(i,3) );
      ni4=complex( slhaPtr->nmnmix(i,4), slhaPtr->imnmnmix(i,4) );
      ni5=complex( slhaPtr->nmnmix(i,5), slhaPtr->imnmnmix(i,5) );
    }

    // Change to positive mass convention
    complex iRot( 0., 1.);
    if (slhaPtr->mass(idNeut(i)) < 0.) {
      ni1 *= iRot;
      ni2 *= iRot;
      ni3 *= iRot;
      ni4 *= iRot;
      ni5 *= iRot;
    }

    // ~chi0 [i] ~chi0 [j] Z : loop over [j]
    for (int j=1; j<=nNeut; j++) {

      // neutralino [j] higgsino components
      complex nj3, nj4;
      if (slhaPtr->modsel(4) >= 1 && slhaPtr->rvnmix.exists()) {
        nj3=complex( slhaPtr->rvnmix(i+3,6), 0.0 );
        nj4=complex( slhaPtr->rvnmix(i+3,7), 0.0 );
      }
      else if (!isNMSSM) {
        nj3=complex( slhaPtr->nmix(j,3), slhaPtr->imnmix(j,3) );
        nj4=complex( slhaPtr->nmix(j,4), slhaPtr->imnmix(j,4) );
      } else {
        nj3=complex( slhaPtr->nmnmix(j,3), slhaPtr->imnmnmix(j,3) );
        nj4=complex( slhaPtr->nmnmix(j,4), slhaPtr->imnmnmix(j,4) );
      }
      // Change to positive mass convention
      if (slhaPtr->mass(idNeut(j)) < 0.) {
        nj3 *= iRot;
        nj4 *= iRot;
      }

      // ~chi0 [i] ~chi0 [j] Z : couplings
      OLpp[i][j] = -0.5 * ni3 * conj(nj3) + 0.5 * ni4 * conj(nj4);
      ORpp[i][j] = 0.5 * conj(ni3) * nj3 - 0.5 * conj(ni4) * nj4;

      // tmp: verbose output
      if (DBSUSY) {
        cout << " OL''  [" << i << "][" << j << "] = "
             << scientific << setw(10) << OLpp[i][j]
             << " OR''  [" << i << "][" << j << "] = "
             << scientific << setw(10) << ORpp[i][j] << endl;
      }

    }

    // ~chi0 [i] ~chi+ [j] W : loop over [j]
    for (int j=1; j<=nChar; j++) {

      // Chargino mixing
      complex uj1, uj2, vj1, vj2;
      if (slhaPtr->modsel(4)<1 || !slhaPtr->rvumix.exists()) {
        uj1=complex( slhaPtr->umix(j,1), slhaPtr->imumix(j,1) );
        uj2=complex( slhaPtr->umix(j,2), slhaPtr->imumix(j,2) );
        vj1=complex( slhaPtr->vmix(j,1), slhaPtr->imvmix(j,1) );
        vj2=complex( slhaPtr->vmix(j,2), slhaPtr->imvmix(j,2) );
      }
      // RPV: ignore lepton-chargino mixing
      else {
        uj1=complex( slhaPtr->rvumix(j+3,4), 0.0 );
        uj2=complex( slhaPtr->rvumix(j+3,5), 0.0 );
        vj1=complex( slhaPtr->rvvmix(j+3,4), 0.0 );
        vj2=complex( slhaPtr->rvvmix(j+3,5), 0.0 );
      }

      // ~chi0 [i] ~chi+ [j] W : couplings
      OL[i][j] = -1.0/sqrt(2.0)*ni4*conj(vj2)+ni2*conj(vj1);
      OR[i][j] = 1.0/sqrt(2.0)*conj(ni3)*uj2+conj(ni2)*uj1;

      // tmp: verbose output
      if (DBSUSY) {
        cout << " OL    [" << i << "][" << j << "] = "
             << scientific << setw(10) << OL[i][j]
             << " OR    [" << i << "][" << j << "] = "
             << scientific << setw(10) << OR[i][j] << endl;
      }
    }


    // ~qqX couplings
    // Quark Charges
    double ed  = -1.0/3.0;
    double T3d = -0.5;
    double eu  =  2.0/3.0;
    double T3u =  0.5;

    // Loop over quark [k] generation
    for (int k=1;k<=3;k++) {

      // Set quark masses
      // Initial guess 0,0,0,mc,mb,mt with the latter from the PDT
      double mu = particleDataPtr->m0(2*k);
      double md = particleDataPtr->m0(2*k-1);
      if (k == 1) { mu=0.0 ; md=0.0; }
      if (k == 2) { md=0.0 ; mu=0.0; }

      // Compute running mass from Yukawas and vevs if possible.
      if (slhaPtr->yd.exists() && slhaPtr->hmix.exists(3)) {
        double ykk=slhaPtr->yd(k,k);
        double v1=slhaPtr->hmix(3)/sqrt(1+pow(tanb,2));
        if (ykk > 0.0) md = ykk * v1 / sqrt(2.0) ;
      }
      if (slhaPtr->yu.exists() && slhaPtr->hmix.exists(3)) {
        double ykk=slhaPtr->yu(k,k);
        double v2=slhaPtr->hmix(3)/sqrt(1.0+1.0/pow(tanb,2));
        if (ykk > 0.0) mu = ykk * v2 / sqrt(2.0) ;
      }

      // tmp: verbose output
      if (DBSUSY) {
        cout  <<  " Gen = " << k << " mu = " << mu << " md = " << md
              << " yUU,DD = " << slhaPtr->yu(k,k) << ","
              << slhaPtr->yd(k,k) << endl;
      }

      // Loop over squark [j] flavour
      for (int j=1;j<=6;j++) {

        // Squark mixing
        complex Rdjk  = complex(Rd(j,k),  imRd(j,k)  );
        complex Rdjk3 = complex(Rd(j,k+3),imRd(j,k+3));
        complex Rujk  = complex(Ru(j,k),  imRu(j,k)  );
        complex Rujk3 = complex(Ru(j,k+3),imRu(j,k+3));

        // ~d[j] d[k] ~chi0[i]
        // Changed according to new notation
        LsddX[j][k][i] = ((ed-T3d)*sinW*ni1 + T3d*cosW*ni2)*conj(Rdjk)
          + md*cosW*ni3*conj(Rdjk3)/2.0/mW/cosb;
        RsddX[j][k][i] = -ed*sinW*conj(ni1)*conj(Rdjk3)
          + md*cosW*conj(ni3)*conj(Rdjk)/2.0/mW/cosb;

        // ~u[j] u[k] ~chi0[i]
        LsuuX[j][k][i] = ((eu-T3u)*sinW*ni1 + T3u*cosW*ni2)*conj(Rujk)
          + mu*cosW*ni4*conj(Rujk3)/2.0/mW/sinb;
        RsuuX[j][k][i] = -eu*sinW*conj(ni1)*conj(Rujk3)
          + mu*cosW*conj(ni4)*conj(Rujk)/2.0/mW/sinb;

        if (DBSUSY) {
          if (abs(LsddX[j][k][i]) > 1e-6) {
            // tmp: verbose output
            cout << " LsddX[" << j << "][" << k << "][" << i << "] = "
                 << scientific << setw(10) << LsddX[j][k][i] << endl;
          }
          if (abs(RsddX[j][k][i]) > 1e-6) {
            // tmp: verbose output
            cout << " RsddX[" << j << "][" << k << "][" << i << "] = "
                 << scientific << setw(10) << RsddX[j][k][i] << endl;
          }
          if (abs(LsuuX[j][k][i]) > 1e-6) {
            // tmp: verbose output
            cout << " LsuuX[" << j << "][" << k << "][" << i << "] = "
                 << scientific << setw(10) << LsuuX[j][k][i] << endl;
          }
          if (abs(RsuuX[j][k][i]) > 1e-6) {
            // tmp: verbose output
            cout << " RsuuX[" << j << "][" << k << "][" << i << "] = "
                 << scientific << setw(10) << RsuuX[j][k][i] << endl;
          }
        }
      }
    }

    // Start slepton couplings
    // Lepton Charges
    double el  = -1.0;
    double T3l = -0.5;
    double ev  =  0.0;
    double T3v =  0.5;

    // Need to define lepton mass
    for (int k=1;k<=3;k++) {
      // Set lepton masses
      double ml(0.0);
      if(k==3) ml = particleDataPtr->m0(15);

      for(int j=1;j<=6;j++){
        LsllX[j][k][i] = 0;
        RsllX[j][k][i] = 0;
        LsvvX[j][k][i] = 0;
        RsvvX[j][k][i] = 0;

        // ~l[j] l[k] ~chi0[i]
        // Changed according to new notation
        LsllX[j][k][i] = ((el-T3l)*sinW*ni1 + T3l*cosW*ni2)*Rsl[j][k]
          + ml*cosW*ni3/2.0/mW/cosb*Rsl[j][k+3];
        RsllX[j][k][i] = -el*sinW*conj(ni1)*Rsl[j][k+3]
          + ml*cosW*conj(ni3)/2.0/mW/cosb*Rsl[j][k];

        if(j<3 && j==k){
          // No sneutrino mixing; only left handed
          // ~v[j] v[k] ~chi0[i]
          LsvvX[j][k][i] = ((ev-T3v)*sinW*ni1 + T3v*cosW*ni2);
        }

        if (DBSUSY) {
          if (abs(LsllX[j][k][i]) > 1e-6) {
            // tmp: verbose output
            cout << " LsllX[" << j << "][" << k << "][" << i << "] = "
                 << scientific << setw(10) << LsllX[j][k][i] << endl;
          }
          if (abs(RsllX[j][k][i]) > 1e-6) {
            // tmp: verbose output
            cout << " RsllX[" << j << "][" << k << "][" << i << "] = "
                 << scientific << setw(10) << RsllX[j][k][i] << endl;
          }
          if (abs(LsvvX[j][k][i]) > 1e-6) {
            // tmp: verbose output
            cout << " LsvvX[" << j << "][" << k << "][" << i << "] = "
                 << scientific << setw(10) << LsvvX[j][k][i] << endl;
          }
        }
      }
    }
  }

  // RPV: check for lepton-chargino mixing (not supported)
  if (slhaPtr->modsel(4) >= 1 && slhaPtr->rvumix.exists()) {
    bool hasCrossTerms = false;
    for (int i=1; i<=3; ++i)
      for (int j=4; j<=5; ++j)
        if (abs(slhaPtr->rvumix(i,j)) > 1e-5
            || abs(slhaPtr->rvvmix(i,j)) > 1e-5) {
          hasCrossTerms = true;
          break;
        }
    if (hasCrossTerms)
      infoPtr->errorMsg("Warning from CoupSUSY::initSUSY: "
        "Lepton-Chargino mixing not supported internally in PYTHIA");
  }

  // Construct ~chi+ couplings
  // sqrt(2)
  double rt2 = sqrt(2.0);

  for (int i=1;i<=nChar;i++) {

    // Ui1, Ui2, Vi1, Vi2
    complex ui1,ui2,vi1,vi2;
    ui1=complex( slhaPtr->umix(i,1), slhaPtr->imumix(i,1) );
    ui2=complex( slhaPtr->umix(i,2), slhaPtr->imumix(i,2) );
    vi1=complex( slhaPtr->vmix(i,1), slhaPtr->imvmix(i,1) );
    vi2=complex( slhaPtr->vmix(i,2), slhaPtr->imvmix(i,2) );

    // ~chi+ [i] ~chi- [j] Z : loop over [j]
    for (int j=1; j<=nChar; j++) {

      // Chargino mixing
      complex uj1, uj2, vj1, vj2;
      uj1=complex( slhaPtr->umix(j,1), slhaPtr->imumix(j,1) );
      uj2=complex( slhaPtr->umix(j,2), slhaPtr->imumix(j,2) );
      vj1=complex( slhaPtr->vmix(j,1), slhaPtr->imvmix(j,1) );
      vj2=complex( slhaPtr->vmix(j,2), slhaPtr->imvmix(j,2) );

      // ~chi+ [i] ~chi- [j] Z : couplings
      OLp[i][j] = -vi1*conj(vj1) - 0.5*vi2*conj(vj2)
        + ( (i == j) ? sin2W : 0.0);
      ORp[i][j] = -conj(ui1)*uj1 - 0.5*conj(ui2)*uj2
        + ( (i == j) ? sin2W : 0.0);

      if (DBSUSY) {
        // tmp: verbose output
        cout << " OL'   [" << i << "][" << j << "] = "
             << scientific << setw(10) << OLp[i][j]
             << " OR'   [" << i << "][" << j << "] = "
             << scientific << setw(10) << ORp[i][j] << endl;
      }
    }

    // Loop over quark [l] flavour
    for (int l=1;l<=3;l++) {

      // Set quark [l] masses
      // Initial guess 0,0,0,mc,mb,mt with the latter from the PDT
      double mul = particleDataPtr->m0(2*l);
      double mdl = particleDataPtr->m0(2*l-1);
      if (l == 1 || l == 2) { mul=0.0 ; mdl=0.0; }

      // Compute running mass from Yukawas and vevs if possible.
      if (slhaPtr->yd.exists() && slhaPtr->hmix.exists(3)) {
        double yll=slhaPtr->yd(l,l);
        double v1=slhaPtr->hmix(3)/sqrt(1+pow(tanb,2));
        if (yll > 0.0) mdl = yll * v1 / sqrt(2.0) ;
      }
      if (slhaPtr->yu.exists() && slhaPtr->hmix.exists(3)) {
        double yll=slhaPtr->yu(l,l);
        double v2=slhaPtr->hmix(3)/sqrt(1.0+1.0/pow(tanb,2));
        if (yll > 0.0) mul = yll * v2 / sqrt(2.0) ;
      }

      // Loop over squark [j] flavour
      for (int j=1;j<=6;j++) {

        //Initialise to zero
        LsduX[j][l][i] = 0.0;
        RsduX[j][l][i] = 0.0;
        LsudX[j][l][i] = 0.0;
        RsudX[j][l][i] = 0.0;

        // Loop over off-diagonal quark [k] generation
        for (int k=1;k<=3;k++) {

          // Set quark [k] masses
          // Initial guess 0,0,0,0,mb,mt with the latter from the PDT
          double muk = particleDataPtr->m0(2*k);
          double mdk = particleDataPtr->m0(2*k-1);
          if (k == 1) { muk=0.0 ; mdk=0.0; }
          if (k == 2) { mdk=0.0 ; muk=0.0; }

          // Compute running mass from Yukawas and vevs if possible.
          if (slhaPtr->yd.exists() && slhaPtr->hmix.exists(3)) {
            double ykk=slhaPtr->yd(k,k);
            double v1=slhaPtr->hmix(3)/sqrt(1+pow(tanb,2));
            if (ykk > 0.0) mdk = ykk * v1 / sqrt(2.0) ;
          }
          if (slhaPtr->yu.exists() && slhaPtr->hmix.exists(3)) {
            double ykk=slhaPtr->yu(k,k);
            double v2=slhaPtr->hmix(3)/sqrt(1.0+1.0/pow(tanb,2));
            if (ykk > 0.0) muk = ykk * v2 / sqrt(2.0) ;
          }

          // CKM matrix (use Pythia one if no SLHA)
          // (NB: could also try input one if no running one found, but
          // would then need to compute from Wolfenstein)
          complex Vlk=VCKMgen(l,k);
          complex Vkl=VCKMgen(k,l);
          if (slhaPtr->vckm.exists()) {
            Vlk=complex(slhaPtr->vckm(l,k),slhaPtr->imvckm(l,k));
            Vkl=complex(slhaPtr->vckm(k,l),slhaPtr->imvckm(k,l));
          }

          // Squark mixing
          complex Rdjk  = complex(Rd(j,k),  imRd(j,k)  );
          complex Rdjk3 = complex(Rd(j,k+3),imRd(j,k+3));
          complex Rujk  = complex(Ru(j,k),  imRu(j,k)  );
          complex Rujk3 = complex(Ru(j,k+3),imRu(j,k+3));


          // ~d[j] u[l] ~chi+[i]
          LsduX[j][l][i] += (ui1*conj(Rdjk)
                             - mdk*ui2*conj(Rdjk3)/rt2/mW/cosb)*Vlk;
          RsduX[j][l][i] -= mul*conj(vi2)*Vlk*conj(Rdjk)/rt2/mW/sinb;

          // ~u[j] d[l] ~chi+[i]
          LsudX[j][l][i] += (conj(vi1)*Rujk
                             - muk*conj(vi2)*Rujk3/rt2/mW/sinb)*Vkl;
          RsudX[j][l][i] -= mdl*ui2*Vkl*Rujk/rt2/mW/cosb;

        }

        if (DBSUSY) {
          if (max(abs(LsduX[j][l][i]),abs(RsduX[j][l][i])) > 1e-6) {
            // tmp: verbose output
            cout << " LsduX[" << j << "][" << l << "][" << i << "] = "
                 << scientific << setw(10) << LsduX[j][l][i];
            cout << " RsduX[" << j << "][" << l << "][" << i << "] = "
                 << scientific << setw(10) << RsduX[j][l][i] << endl;
          }
          if (max(abs(LsudX[j][l][i]),abs(RsudX[j][l][i])) > 1e-6) {
            // tmp: verbose output
            cout << " LsudX[" << j << "][" << l << "][" << i << "] = "
                 << scientific << setw(10) << LsudX[j][l][i];
            cout << " RsudX[" << j << "][" << l << "][" << i << "] = "
                 << scientific << setw(10) << RsudX[j][l][i] << endl;
          }
        }
      }
    }

    // Loop over slepton [j] flavour
    for (int j=1;j<=6;j++) {
      for (int k=1;k<=3;k++) {

        LslvX[j][k][i] = 0.0;
        RslvX[j][k][i] = 0.0;
        LsvlX[j][k][i] = 0.0;
        RsvlX[j][k][i] = 0.0;

        // Set lepton [k] masses
        double ml = 0.0;
        if (k == 3) ml = particleDataPtr->m0(15);

        if(j==k || j==k+3){ // No lepton mixing

          // ~l[j] v[l] ~chi+[i]
          LslvX[j][k][i] += ui1- ml*ui2/rt2/mW/cosb;
          RslvX[j][k][i] -= ml*conj(vi2)/rt2/mW/sinb;

          // ~v[j] l[l] ~chi+[i]
          if(j<=3){ // No right handed sneutrinos
            LsvlX[j][k][i] += conj(vi1) - ml*conj(vi2)/rt2/mW/sinb;
          }
        }

        if (DBSUSY) {
          if (max(abs(LslvX[j][k][i]),abs(RslvX[j][k][i])) > 1e-6) {
            // tmp: verbose output
            cout << " LslvX[" << j << "][" << k << "][" << i << "] = "
                 << scientific << setw(10) << LslvX[j][k][i];
            cout << " RslvX[" << j << "][" << k << "][" << i << "] = "
                 << scientific << setw(10) << RslvX[j][k][i] << endl;
          }
          if (max(abs(LsvlX[j][k][i]),abs(RsvlX[j][k][i])) > 1e-6) {
            // tmp: verbose output
            cout << " LsvlX[" << j << "][" << k << "][" << i << "] = "
                 << scientific << setw(10) << LsvlX[j][k][i];
            cout << " RsvlX[" << j << "][" << k << "][" << i << "] = "
                 << scientific << setw(10) << RsvlX[j][k][i] << endl;
          }
        }
      }
    }
  }

  // Shorthand for RPV couplings
  // The input LNV lambda couplings
  LHtensor3Block<3> rvlle(slhaPtr->rvlamlle);
  // The input LNV lambda' couplings
  LHtensor3Block<3> rvlqd(slhaPtr->rvlamlqd);
  // The input BNV lambda'' couplings
  LHtensor3Block<3> rvudd(slhaPtr->rvlamudd);

  isLLE = false;
  isLQD = false;
  isUDD = false;

  // Symmetry properties
  if(rvlle.exists()){
    isLLE = true;
    for(int i=1;i<=3;i++){
      for(int j=i;j<=3;j++){
        //lambda(i,j,k)=-lambda(j,i,k)
        for(int k=1;k<=3;k++){
          if(i==j){
            rvLLE[i][j][k] = 0.0;
          }else{
            rvLLE[i][j][k] = rvlle(i,j,k);
            rvLLE[j][i][k] = -rvlle(i,j,k);
          }
        }
      }
    }
  }

  if(rvlqd.exists()){
    isLQD=true;
    for(int i=1;i<=3;i++){
      for(int j=1;j<=3;j++){
        for(int k=1;k<=3;k++){
            rvLQD[i][j][k] = rvlqd(i,j,k);
        }
      }
    }
  }

  //lambda''(k,j,i)=-lambda''(k,i,j)

  if(rvudd.exists()){
    isUDD = true;
    for(int i=1;i<=3;i++){
      for(int j=i;j<=3;j++){
        for(int k=1;k<=3;k++){
          if(i==j){
            rvUDD[k][i][j] = 0.0;
          }else{
            rvUDD[k][i][j] = rvudd(k,i,j);
            rvUDD[k][j][i] = -rvudd(k,i,j);
          }
        }
      }
    }
  }

  if(DBSUSY){
    for(int i=1;i<=3;i++){
      for(int j=1;j<=3;j++){
        for(int k=1;k<=3;k++){
          if(rvlle.exists())
            cout<<"LambdaLLE["<<i<<"]["<<j<<"]["<<k<<"]="<<rvLLE[i][j][k]<<" ";
          if(rvlqd.exists())
            cout<<"LambdaLQD["<<i<<"]["<<j<<"]["<<k<<"]="<<rvLQD[i][j][k]<<" ";
          if(rvudd.exists())
            cout<<"LambdaUDD["<<i<<"]["<<j<<"]["<<k<<"]="<<rvUDD[i][j][k]
                <<"\n";
        }
      }
    }
  }

  // Store the squark mixing matrix
  for(int i=1;i<=6;i++){
    for(int j=1;j<=3;j++){
      Rusq[i][j]   = complex(Ru(i,j),  imRu(i,j)  );
      Rusq[i][j+3] = complex(Ru(i,j+3),imRu(i,j+3));
      Rdsq[i][j]   = complex(Rd(i,j),  imRd(i,j)  );
      Rdsq[i][j+3] = complex(Rd(i,j+3),imRd(i,j+3));
    }
  }

  if(DBSUSY){
    cout<<"Ru"<<endl;
    for(int i=1;i<=6;i++){
      for(int j=1;j<=6;j++){
        cout << real(Rusq[i][j])<<"  ";
      }
      cout<<endl;
    }
    cout<<"Rd"<<endl;
    for(int i=1;i<=6;i++){
      for(int j=1;j<=6;j++){
        cout << real(Rdsq[i][j])<<"  ";
      }
      cout<<endl;
    }
  }

  // Let everyone know we are ready
  isInit = true;
}

//--------------------------------------------------------------------------

// Return neutralino flavour codes.

int CoupSUSY::idNeut(int idChi) {

  int id = 0;
  if      (idChi == 1) id = 1000022;
  else if (idChi == 2) id = 1000023;
  else if (idChi == 3) id = 1000025;
  else if (idChi == 4) id = 1000035;
  else if (idChi == 5) id = 1000045;
  return id;
}

//--------------------------------------------------------------------------

// Return chargino flavour codes.

int CoupSUSY::idChar(int idChi) {

  int id = 0;
  if      (idChi ==  1) id =  1000024;
  else if (idChi == -1) id = -1000024;
  else if (idChi ==  2) id =  1000037;
  else if (idChi == -2) id = -1000037;
  return id;
}

//--------------------------------------------------------------------------

// Return sup flavour codes.

int CoupSUSY::idSup(int iSup) {

  int id = 0;
  int sgn = ( iSup > 0 ) ? 1 : -1;
  iSup = abs(iSup);
  if      (iSup ==  1) id =  1000002;
  else if (iSup ==  2) id =  1000004;
  else if (iSup ==  3) id =  1000006;
  else if (iSup ==  4) id =  2000002;
  else if (iSup ==  5) id =  2000004;
  else if (iSup ==  6) id =  2000006;
  return sgn*id;
}

//--------------------------------------------------------------------------

// Return sdown flavour codes

int CoupSUSY::idSdown(int iSdown) {

  int id = 0;
  int sgn = ( iSdown > 0 ) ? 1 : -1;
  iSdown = abs(iSdown);
  if      (iSdown ==  1) id =  1000001;
  else if (iSdown ==  2) id =  1000003;
  else if (iSdown ==  3) id =  1000005;
  else if (iSdown ==  4) id =  2000001;
  else if (iSdown ==  5) id =  2000003;
  else if (iSdown ==  6) id =  2000005;
  return sgn*id;
}

//--------------------------------------------------------------------------

// Function to return slepton flavour codes

int CoupSUSY::idSlep(int iSlep) {

  int id = 0;
  int sgn = ( iSlep > 0 ) ? 1 : -1;
  iSlep = abs(iSlep);
  if      (iSlep ==  1) id =  1000011;
  else if (iSlep ==  2) id =  1000013;
  else if (iSlep ==  3) id =  1000015;
  else if (iSlep ==  4) id =  2000011;
  else if (iSlep ==  5) id =  2000013;
  else if (iSlep ==  6) id =  2000015;
  return sgn*id;
}

//--------------------------------------------------------------------------

// Return neutralino code; zero if not a (recognized) neutralino

int CoupSUSY::typeNeut(int idPDG) {
  int type = 0;
  int idAbs = abs(idPDG);
  if(idAbs == 1000022) type = 1;
  else if(idAbs == 1000023) type = 2;
  else if(idAbs == 1000025) type = 3;
  else if(idAbs == 1000035) type = 4;
  else if(isNMSSM && idAbs == 1000045) type = 5;
  return type;

}

//--------------------------------------------------------------------------

// Check whether particle is a Chargino

int CoupSUSY::typeChar(int idPDG) {
  int type = 0;
  if(abs(idPDG) == 1000024) type = 1;
  else if (abs(idPDG) == 1000037)type = 2;
  return type;
}

//==========================================================================

} // end namespace Pythia8
