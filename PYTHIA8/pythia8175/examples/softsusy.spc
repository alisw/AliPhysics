# SOFTSUSY1.9
# B.C. Allanach, Comput. Phys. Commun. 143 (2002) 305-331, hep-ph/0104145
Block SPINFO         # Program information
     1   SOFTSUSY    # spectrum calculator
     2   1.9.1       # version number
Block MODSEL  # Select model
     1    1   # sugra
Block SMINPUTS   # Standard Model inputs
     1    1.27934000e+02   # alpha_em^(-1)(MZ) SM MSbar
     2    1.16637000e-05   # G_Fermi
     3    1.17200000e-01   # alpha_s(MZ)MSbar
     4    9.11876000e+01   # MZ(pole)
     5    4.25000000e+00   # Mb(mb)
     6    1.74300000e+02   # Mtop(pole)
     7    1.77700000e+00   # Mtau(pole)
Block MINPAR  # SUSY breaking input parameters
     3    1.00000000e+01   # tanb
     4    1.00000000e+00   # sign(mu)
     1    1.00000000e+02   # m0
     2    2.50000000e+02   # m12
     5   -1.00000000e+02   # A0
# Low energy data in SOFTSUSY: MIXING=-1 TOLERANCE=1.00000000e-03
# mgut=2.45916471e+16 GeV
Block MASS   # Mass spectrum
#PDG code      mass              particle
        24     8.04191121e+01   # MW
        25     1.10762378e+02   # h0
        35     4.00599584e+02   # H0
        36     4.00231463e+02   # A0
        37     4.08513284e+02   # H+
   1000001     5.72700955e+02   # ~d_L
   1000002     5.67251814e+02   # ~u_L
   1000003     5.72700955e+02   # ~s_L
   1000004     5.67251814e+02   # ~c_L
   1000005     5.15211952e+02   # ~b_1
   1000006     3.95920984e+02   # ~t_1
   1000011     2.04276615e+02   # ~e_L
   1000012     1.88657729e+02   # ~nue_L
   1000013     2.04276615e+02   # ~mu_L
   1000014     1.88657729e+02   # ~numu_L
   1000015     1.36227147e+02   # ~stau_1
   1000016     1.87773326e+02   # ~nu_tau_L
   1000021     6.07604198e+02   # ~g
   1000022     9.72852615e+01   # ~neutralino(1)
   1000023     1.80961862e+02   # ~neutralino(2)
   1000024     1.80378828e+02   # ~chargino(1)
   1000025    -3.64435115e+02   # ~neutralino(3)
   1000035     3.83135773e+02   # ~neutralino(4)
   1000037     3.83371870e+02   # ~chargino(2)
   2000001     5.46070490e+02   # ~d_R
   2000002     5.46999685e+02   # ~u_R
   2000003     5.46070490e+02   # ~s_R
   2000004     5.46999685e+02   # ~c_R
   2000005     5.43966766e+02   # ~b_2
   2000006     5.85698733e+02   # ~t_2
   2000011     1.45526717e+02   # ~e_R
   2000013     1.45526717e+02   # ~mu_R
   2000015     2.08222793e+02   # ~stau_2
# Higgs mixing
Block alpha   # Effective Higgs mixing parameter
          -1.13732831e-01   # alpha
Block stopmix  # stop mixing matrix
  1  1     5.38083886e-01   # O_{11}
  1  2     8.42891293e-01   # O_{12}
  2  1     8.42891293e-01   # O_{21}
  2  2    -5.38083886e-01   # O_{22}
Block sbotmix  # sbottom mixing matrix
  1  1     9.47744273e-01   # O_{11}
  1  2     3.19031021e-01   # O_{12}
  2  1    -3.19031021e-01   # O_{21}
  2  2     9.47744273e-01   # O_{22}
Block staumix  # stau mixing matrix
  1  1     2.80956141e-01   # O_{11}
  1  2     9.59720609e-01   # O_{12}
  2  1     9.59720609e-01   # O_{21}
  2  2    -2.80956141e-01   # O_{22}
Block nmix  # neutralino mixing matrix
  1  1     9.86066377e-01   # N_{1,1}
  1  2    -5.46292061e-02   # N_{1,2}
  1  3     1.47649927e-01   # N_{1,3}
  1  4    -5.37424305e-02   # N_{1,4}
  2  1     1.02062420e-01   # N_{2,1}
  2  2     9.42721210e-01   # N_{2,2}
  2  3    -2.74985600e-01   # N_{2,3}
  2  4     1.58880154e-01   # N_{2,4}
  3  1    -6.04575099e-02   # N_{3,1}
  3  2     8.97030908e-02   # N_{3,2}
  3  3     6.95501068e-01   # N_{3,3}
  3  4     7.10335491e-01   # N_{3,4}
  4  1    -1.16624405e-01   # N_{4,1}
  4  2     3.16616055e-01   # N_{4,2}
  4  3     6.47194471e-01   # N_{4,3}
  4  4    -6.83587843e-01   # N_{4,4}
Block Umix  # chargino U mixing matrix 
  1  1     9.15531658e-01   # U_{1,1}
  1  2    -4.02245924e-01   # U_{1,2}
  2  1     4.02245924e-01   # U_{2,1}
  2  2     9.15531658e-01   # U_{2,2}
Block Vmix  # chargino V mixing matrix 
  1  1     9.72345994e-01   # V_{1,1}
  1  2    -2.33545003e-01   # V_{1,2}
  2  1     2.33545003e-01   # V_{2,1}
  2  2     9.72345994e-01   # V_{2,2}
Block gauge Q= 4.64231969e+02  
     1     3.60968173e-01   # g'(Q)MSSM DRbar
     2     6.46474399e-01   # g(Q)MSSM DRbar
     3     1.09626470e+00   # g3(Q)MSSM DRbar
Block yu Q= 4.64231969e+02  
  3  3     8.89731484e-01   # Yt(Q)MSSM DRbar
Block yd Q= 4.64231969e+02  
  3  3     1.39732269e-01   # Yb(Q)MSSM DRbar
Block ye Q= 4.64231969e+02  
  3  3     1.00914051e-01   # Ytau(Q)MSSM DRbar
Block hmix Q= 4.64231969e+02  # Higgs mixing parameters
     1     3.58339654e+02   # mu(Q)MSSM DRbar
     2     9.75145219e+00   # tan beta(Q)MSSM DRbar
     3     2.44923803e+02   # higgs vev(Q)MSSM DRbar
     4     1.67100152e+05   # mA^2(Q)MSSM DRbar
Block msoft Q=4.64231969e+02 # MSSM DRbar SUSY breaking parameters
     1     1.01439997e+02   # M_1(Q)
     2     1.91579315e+02   # M_2(Q)
     3     5.86586195e+02   # M_3(Q)
    21     3.23914077e+04   # mH1^2(Q)
    22    -1.29413007e+05   # mH2^2(Q)
    31     1.99042560e+02   # meL(Q)
    32     1.99042560e+02   # mmuL(Q)
    33     1.98204510e+02   # mtauL(Q)
    34     1.38811933e+02   # meR(Q)
    35     1.38811933e+02   # mmuR(Q)
    36     1.36392545e+02   # mtauR(Q)
    41     5.50815976e+02   # mqL1(Q)
    42     5.50815976e+02   # mqL2(Q)
    43     4.99361608e+02   # mqL3(Q)
    44     5.28861326e+02   # muR(Q)
    45     5.28861326e+02   # mcR(Q)
    46     4.18454191e+02   # mtR(Q)
    47     5.26100270e+02   # mdR(Q)
    48     5.26100270e+02   # msR(Q)
    49     5.22780488e+02   # mbR(Q)
Block au Q= 4.64231969e+02  
  1  1     0.00000000e+00   # Au(Q)MSSM DRbar
  2  2     0.00000000e+00   # Ac(Q)MSSM DRbar
  3  3    -5.04520155e+02   # At(Q)MSSM DRbar
Block ad Q= 4.64231969e+02  
  1  1     0.00000000e+00   # Ad(Q)MSSM DRbar
  2  2     0.00000000e+00   # As(Q)MSSM DRbar
  3  3    -7.97104366e+02   # Ab(Q)MSSM DRbar
Block ae Q= 4.64231969e+02  
  1  1     0.00000000e+00   # Ae(Q)MSSM DRbar
  2  2     0.00000000e+00   # Amu(Q)MSSM DRbar
  3  3    -2.56146632e+02   # Atau(Q)MSSM DRbar
