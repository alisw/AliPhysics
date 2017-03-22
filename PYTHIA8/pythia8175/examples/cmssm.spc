# SOFTSUSY3.3.0 SLHA compliant output
# B.C. Allanach, Comput. Phys. Commun. 143 (2002) 305-331, hep-ph/0104145
# CMSSM 10.1.1 input - see arXiv:1109.3859
Block SPINFO          # Program information
     1    SOFTSUSY    # spectrum calculator
     2    3.3.0       # version number
Block MODSEL  # Select model
     1    1   # sugra
Block SMINPUTS             # Standard Model inputs
     1    1.27934000e+02   # alpha_em^(-1)(MZ) SM MSbar
     2    1.16637000e-05   # G_Fermi
     3    1.17200000e-01   # alpha_s(MZ)MSbar
     4    9.11876000e+01   # MZ(pole)
     5    4.25000000e+00   # mb(mb)
     6    1.74300000e+02   # Mtop(pole)
     7    1.77700000e+00   # Mtau(pole)
Block MINPAR               # SUSY breaking input parameters
     1    1.25000000e+02   # m0
     2    5.00000000e+02   # m12
     3    1.00000000e+01   # tanb
     4    1.00000000e+00   # sign(mu)
     5    0.00000000e+00   # A0
# Low energy data in SOFTSUSY: MIXING=0 TOLERANCE=1.00000000e-03
# mgut=1.87825797e+16 GeV
Block MASS                      # Mass spectrum
# PDG code     mass             particle
        24     8.04613815e+01   # MW
        25     1.15412339e+02   # h0
        35     7.21892898e+02   # H0
        36     7.21606682e+02   # A0
        37     7.26310615e+02   # H+
   1000021     1.14687814e+03   # ~g
   1000022     2.04168238e+02   # ~neutralino(1)
   1000023     3.86197315e+02   # ~neutralino(2)
   1000024     3.86332261e+02   # ~chargino(1)
   1000025    -6.39638883e+02   # ~neutralino(3)
   1000035     6.53032956e+02   # ~neutralino(4)
   1000037     6.52679913e+02   # ~chargino(2)
   1000001     1.05561088e+03   # ~d_L
   1000002     1.05278695e+03   # ~u_L
   1000003     1.05560831e+03   # ~s_L
   1000004     1.05278437e+03   # ~c_L
   1000005     9.65462335e+02   # ~b_1
   1000006     8.03540912e+02   # ~t_1
   1000011     3.61778036e+02   # ~e_L
   1000012     3.52902884e+02   # ~nue_L
   1000013     3.61861990e+02   # ~mu_L
   1000014     3.52899536e+02   # ~numu_L
   1000015     2.22692429e+02   # ~stau_1
   1000016     3.51730554e+02   # ~nu_tau_L
   2000001     1.01102319e+03   # ~d_R
   2000002     1.01434248e+03   # ~u_R
   2000003     1.01102050e+03   # ~s_R
   2000004     1.01433974e+03   # ~c_R
   2000005     1.00775879e+03   # ~b_2
   2000006     1.01102814e+03   # ~t_2
   2000011     2.29874717e+02   # ~e_R
   2000013     2.29864278e+02   # ~mu_R
   2000015     3.62940076e+02   # ~stau_2
Block alpha                     # Effective Higgs mixing parameter
          -1.06797662e-01       # alpha
Block nmix                  # neutralino mixing matrix
  1  1     9.95886986e-01   # N_{1,1}
  1  2    -1.72694202e-02   # N_{1,2}
  1  3     8.22024839e-02   # N_{1,3}
  1  4    -3.39651202e-02   # N_{1,4}
  2  1     3.74375188e-02   # N_{2,1}
  2  2     9.71264892e-01   # N_{2,2}
  2  3    -1.95766260e-01   # N_{2,3}
  2  4     1.30071184e-01   # N_{2,4}
  3  1    -3.30727556e-02   # N_{3,1}
  3  2     4.82085025e-02   # N_{3,2}
  3  3     7.03472841e-01   # N_{3,3}
  3  4     7.08313557e-01   # N_{3,4}
  4  1    -7.55892561e-02   # N_{4,1}
  4  2     2.32426798e-01   # N_{4,2}
  4  3     6.78265645e-01   # N_{4,3}
  4  4    -6.92978905e-01   # N_{4,4}
Block Umix                  # chargino U mixing matrix 
  1  1     9.60950435e-01   # U_{1,1}
  1  2    -2.76720547e-01   # U_{1,2}
  2  1     2.76720547e-01   # U_{2,1}
  2  2     9.60950435e-01   # U_{2,2}
Block Vmix                  # chargino V mixing matrix 
  1  1     9.82582136e-01   # V_{1,1}
  1  2    -1.85828808e-01   # V_{1,2}
  2  1     1.85828808e-01   # V_{2,1}
  2  2     9.82582136e-01   # V_{2,2}
Block stopmix               # stop mixing matrix
  1  1     4.24985787e-01   # F_{11}
  1  2     9.05200023e-01   # F_{12}
  2  1     9.05200023e-01   # F_{21}
  2  2    -4.24985787e-01   # F_{22}
Block sbotmix               # sbottom mixing matrix
  1  1     9.75941340e-01   # F_{11}
  1  2     2.18033254e-01   # F_{12}
  2  1    -2.18033254e-01   # F_{21}
  2  2     9.75941340e-01   # F_{22}
Block staumix               # stau mixing matrix
  1  1     1.42739201e-01   # F_{11}
  1  2     9.89760335e-01   # F_{12}
  2  1     9.89760335e-01   # F_{21}
  2  2    -1.42739201e-01   # F_{22}
Block gauge Q= 8.73811408e+02  # SM gauge couplings
     1     3.62408332e-01   # g'(Q)MSSM DRbar
     2     6.43000889e-01   # g(Q)MSSM DRbar
     3     1.06090510e+00   # g3(Q)MSSM DRbar
Block yu Q= 8.73811408e+02  
  3  3     8.67028746e-01   # Yt(Q)MSSM DRbar
Block yd Q= 8.73811408e+02  
  3  3     1.35275115e-01   # Yb(Q)MSSM DRbar
Block ye Q= 8.73811408e+02  
  3  3     1.00502013e-01   # Ytau(Q)MSSM DRbar
Block hmix Q= 8.73811408e+02 # Higgs mixing parameters
     1     6.34225757e+02    # mu(Q)MSSM DRbar
     2     9.66949856e+00    # tan beta(Q)MSSM DRbar
     3     2.44043659e+02    # higgs vev(Q)MSSM DRbar
     4     5.40020729e+05    # mA^2(Q)MSSM DRbar
Block msoft Q= 8.73811408e+02  # MSSM DRbar SUSY breaking parameters
     1     2.09217931e+02      # M_1(Q)
     2     3.88513213e+02      # M_2(Q)
     3     1.11375823e+03      # M_3(Q)
    21     1.09175993e+05      # mH1^2(Q)
    22    -3.89577899e+05      # mH2^2(Q)
    31     3.54297503e+02      # meL(Q)
    32     3.54294169e+02      # mmuL(Q)
    33     3.53288859e+02      # mtauL(Q)
    34     2.22024882e+02      # meR(Q)
    35     2.22014061e+02      # mmuR(Q)
    36     2.18730815e+02      # mtauR(Q)
    41     1.01868349e+03      # mqL1(Q)
    42     1.01868086e+03      # mqL2(Q)
    43     9.38174894e+02      # mqL3(Q)
    44     9.80660335e+02      # muR(Q)
    45     9.80657548e+02      # mcR(Q)
    46     8.05547302e+02      # mtR(Q)
    47     9.76062679e+02      # mdR(Q)
    48     9.76059936e+02      # msR(Q)
    49     9.71058708e+02      # mbR(Q)
Block au Q= 8.73811408e+02  
  1  1    -1.13744947e+03      # Au(Q)MSSM DRbar
  2  2    -1.13744435e+03      # Ac(Q)MSSM DRbar
  3  3    -8.73322525e+02      # At(Q)MSSM DRbar
Block ad Q= 8.73811408e+02  
  1  1    -1.39629474e+03      # Ad(Q)MSSM DRbar
  2  2    -1.39628999e+03      # As(Q)MSSM DRbar
  3  3    -1.30365153e+03      # Ab(Q)MSSM DRbar
Block ae Q= 8.73811408e+02  
  1  1    -2.99482647e+02      # Ae(Q)MSSM DRbar
  2  2    -2.99477266e+02      # Amu(Q)MSSM DRbar
  3  3    -2.97855048e+02      # Atau(Q)MSSM DRbar
