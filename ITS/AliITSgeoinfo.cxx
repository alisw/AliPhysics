/*
#include <Riostream.h>
#include <stdlib.h>
#include <TMath.h>
#include <TRandom.h>
#include <TBranch.h>
#include <TVector.h>
#include <TClonesArray.h>
#include <TROOT.h>
#include <TObjectTable.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>

#include "TParticle.h"
#include "AliRun.h"
#include "AliITS.h"
#include "AliITSMap.h"
#include "AliITSDetType.h"
#include "AliITSClusterFinder.h"
#include "AliITSsimulationSPD.h"
#include "AliITSsimulationSDD.h"
#include "AliITSsimulationSSD.h"
#include "AliITSresponse.h"
#include "AliITSsegmentationSPD.h"
#include "AliITSresponseSPD.h"
#include "AliITSresponseSPDbari.h"
#include "AliITSsegmentationSDD.h"
#include "AliITSresponseSDD.h"
#include "AliITSsegmentationSSD.h"
#include "AliITSresponseSSD.h"
#include "AliITShit.h"
#include "AliITSgeomSPD.h"
#include "AliITSgeomSDD.h"
#include "AliITSgeomSSD.h"
#include "AliITSgeom.h"
#include "AliITSdigit.h"
#include "AliITSmodule.h"
#include "AliITSRecPoint.h"
#include "AliITSRawCluster.h"
#include "stdlib.h"
#include "AliKalmanTrack.h" 
#include "AliMagF.h"



#include "AliITStrack.h"
#include "AliITSiotrack.h"
#include "AliITStracking.h"
#include "AliITSRad.h"   
#include "../TPC/AliTPC.h"
#include "../TPC/AliTPCParam.h"
#include "../TPC/AliTPCtracker.h"

#include "AliITSTrackerV1.h"
*/
#include "AliITSgeoinfo.h"


ClassImp(AliITSgeoinfo)


//________________________________________________________________

AliITSgeoinfo::AliITSgeoinfo() {

  for(Int_t i=0; i<6; i++) {Nlad[i]=0; Ndet[i]=0; Avrad[i]=0.; Detx[i]=0.; Detz[i]=0.;}
  
}

