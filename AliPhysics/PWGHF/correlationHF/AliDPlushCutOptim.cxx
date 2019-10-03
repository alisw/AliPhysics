//
//  Class to crate a TTree for variable optimization
//
//-----------------------------------------------------------------------
//  Author Shyam Kumar : copied from D0 class
//  IIT Bombay
//  shyam.kumar@cern.ch
//-----------------------------------------------------------------------

#include <iostream>
#include "TObject.h"
#include "TMath.h"
#include "TFile.h"
#include "TDirectoryFile.h"
#include "TList.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TDatabasePDG.h"

#include "AliDPlushCutOptim.h"

//___________________________________________________________________________________________
AliDPlushCutOptim::AliDPlushCutOptim():
// default constructor
invMass(0),
cent(0),
pt(0),
pTk(0),
pTpi(0),
d0k(0),
d0pi(0),
cosPA(0),
sumd0square(0),
dca(0),
dlXY(0),
cosPAXY(0),
topom(0)
{

}
