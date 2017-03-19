//
//  Class to crate a TTree for variable optimization
//
//-----------------------------------------------------------------------
//  Author F.Colamaria
//  INFN Bari
//  fabio.colamaria@cern.ch
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

#include "AliD0hCutOptim.h"

//___________________________________________________________________________________________
AliD0hCutOptim::AliD0hCutOptim():
// default constructor
invMass(0),
cent(0),
pt(0),
dca(0),
cosThSt(0),
pTk(0),
pTpi(0),
d0k(0),
d0pi(0),
d0xd0(0),
cosThPt(0),
normLxy(0),
topom(0)
{

}
