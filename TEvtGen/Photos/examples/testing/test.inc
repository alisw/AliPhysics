EVENTS = 1000000

VIEWER=gv
COMPARE=../compare.sh

# Programs
PHOTOS_ONLY   =../photos_test.exe
PHOTOS_TAUOLA =../photos_tauola_test.exe

# Pythia modes
DEFAULT_MODE =0 # ee@200GeV    (for Wenu, Wmunu, WmunuNLO)
Z_MU_MODE    =1 # ee@91.187GeV (for Zmumu, ZmumuNLO, Ztautau, ScalarNLO)
PP_MODE      =2 # pp@14TeV     (for Zee, ttbar)
HIGGS_MODE   =3 # ee@500GeV    (for Htautau)

# Special modes
TOP_MODE     =1 # Turn on MC-Tester fix for Top decays
NLO_MODE     =2 # Turn on Photos NLO corrections

# Tauola generation modes
PIPI_MODE    =3 # ee@91.187GeV (for Ztautau)
RHORHO_MODE  =4 # ee@500GeV    (for Htautau)
K-STAR_MODE  =7 # ee@91.187GeV (for ScalarNLO)
