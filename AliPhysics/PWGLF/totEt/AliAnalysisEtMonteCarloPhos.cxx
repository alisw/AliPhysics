//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Base class for MC analysis, for PHOS
//  - MC output
//  implementation file 
//
//*-- Authors: Oystein Djuvsland (Bergen), David Silvermyr (ORNL)
//_________________________________________________________________________
#include "AliAnalysisEtMonteCarloPhos.h"
#include "AliAnalysisEtSelectorPhos.h"
#include "AliAnalysisEtCuts.h"
#include "AliESDtrack.h"
#include <iostream>
#include "AliPHOSGeoUtils.h"
#include "TFile.h"
#include "TH2I.h"
#include <AliPHOSGeometry.h>

using namespace std;

ClassImp(AliAnalysisEtMonteCarloPhos);


AliAnalysisEtMonteCarloPhos::AliAnalysisEtMonteCarloPhos():AliAnalysisEtMonteCarlo()
,fBadMapM2(0)
,fBadMapM3(0)
,fBadMapM4(0)
,fGeoUtils(0)
{
   fHistogramNameSuffix = TString("PhosMC");
}

AliAnalysisEtMonteCarloPhos::~AliAnalysisEtMonteCarloPhos()
{ // dtor
  delete fBadMapM2;
  delete fBadMapM3;
  delete fBadMapM4;
  delete fGeoUtils;
}


void AliAnalysisEtMonteCarloPhos::Init()
{ // Init
  AliAnalysisEtMonteCarlo::Init();
  fSelector = new AliAnalysisEtSelectorPhos(fCuts);
    
  fDetectorRadius = fCuts->GetGeometryPhosDetectorRadius();
  fSingleCellEnergyCut = fCuts->GetReconstructedPhosSingleCellEnergyCut();

 // ifstream f("badchannels.txt", ios::in);
  TFile *f = TFile::Open("badchannels.root", "READ");
  
  fBadMapM2 = (TH2I*)f->Get("bad_channels_m2");
   fBadMapM3 = (TH2I*)f->Get("bad_channels_m3");
   fBadMapM4 = (TH2I*)f->Get("bad_channels_m4");
// 
   fGeoUtils = AliPHOSGeometry::GetInstance("IHEP");
  
}


Bool_t AliAnalysisEtMonteCarloPhos::TooCloseToBadChannel(const AliESDCaloCluster &cluster) const
{ // too close to bad channel?

    Float_t gPos[3];
    cluster.GetPosition(gPos);
    Int_t relId[4];
    TVector3 glVec(gPos);
    fGeoUtils->GlobalPos2RelId(glVec, relId);

    TVector3 locVec;
    fGeoUtils->Global2Local(locVec, glVec, relId[0]);
//    std::cout << fGeoUtils << std::endl;
    //std::cout << relId[0] << " " << cluster.IsPHOS() <<  std::endl;
    //std::cout << locVec[0] << " " << " " << locVec[1] << " " << locVec[2] << std::endl;
    for (Int_t x = 0; x < fBadMapM2->GetNbinsX(); x++)
    {
        for (Int_t z = 0; z < fBadMapM2->GetNbinsY(); z++)
        {
            if (relId[0] == 3)
            {
                if (fBadMapM2->GetBinContent(x+1, z+1) != 0)
                {
		    Int_t tmpRel[4];
		    tmpRel[0] = 3;
		    tmpRel[1] = 0;
		    tmpRel[2] = x+1;
		    tmpRel[3] = z+1;
		    
		    Float_t tmpX;
		    Float_t tmpZ;
		    fGeoUtils->RelPosInModule(tmpRel, tmpX, tmpZ);

                    Float_t distance = TMath::Sqrt((tmpX-locVec[0])*(tmpX-locVec[0]) + (tmpZ - locVec[2])*(tmpZ-locVec[2]));
                    //Float_t distance = TMath::Sqrt((x-relId[3])*(x-relId[3]) + (z - relId[2])*(z-relId[2]));
		    if (distance < fCuts->GetPhosBadDistanceCut())
                    {
//		      std::cout << "d: " << distance << ", cut: " << fCuts->GetPhosBadDistanceCut() << std::endl;
                    
                        return kTRUE;
                    }
                }
            }
            if (relId[0] == 2)
            {
                if (fBadMapM3->GetBinContent(x+1, z+1) != 0)
                {
		    Int_t tmpRel[4];
		    tmpRel[0] = 2;
		    tmpRel[1] = 0;
		    tmpRel[2] = x+1;
		    tmpRel[3] = z+1;
		    
		    Float_t tmpX;
		    Float_t tmpZ;
		    fGeoUtils->RelPosInModule(tmpRel, tmpX, tmpZ);

                    Float_t distance = TMath::Sqrt((tmpX-locVec[0])*(tmpX-locVec[0]) + (tmpZ - locVec[2])*(tmpZ-locVec[2]));

//                    Float_t distance = TMath::Sqrt((x-locVec[0])*(x-locVec[0]) + (z - locVec[2])*(z-locVec[2]));
		    //Float_t distance = TMath::Sqrt((x-relId[3])*(x-relId[3]) + (z - relId[2])*(z-relId[2]));
                    if (distance < fCuts->GetPhosBadDistanceCut())
                    {
//		        std::cout << "d: " << distance << ", cut: " << fCuts->GetPhosBadDistanceCut() << std::endl;
                    
                        return kTRUE;
                    }
                }
            }
            if (relId[0] == 1)
            {
                if (fBadMapM4->GetBinContent(x+1, z+1) != 0)
                {
		    Int_t tmpRel[4];
		    tmpRel[0] = 1;
		    tmpRel[1] = 0;
		    tmpRel[2] = x+1;
		    tmpRel[3] = z+1;
		    
		    Float_t tmpX;
		    Float_t tmpZ;
		    fGeoUtils->RelPosInModule(tmpRel, tmpX, tmpZ);

                    Float_t distance = TMath::Sqrt((tmpX-locVec[0])*(tmpX-locVec[0]) + (tmpZ - locVec[2])*(tmpZ-locVec[2]));

//                    Float_t distance = TMath::Sqrt((x-locVec[0])*(x-locVec[0]) + (z - locVec[2])*(z-locVec[2]));
		    //Float_t distance = TMath::Sqrt((x-relId[3])*(x-relId[3]) + (z - relId[2])*(z-relId[2]));
                    if (distance < fCuts->GetPhosBadDistanceCut())
                    {
//		        std::cout << "d: " << distance << ", cut: " << fCuts->GetPhosBadDistanceCut() << std::endl;
                    
                        return kTRUE;
                    }
                }
            }

        }
    }

    return kFALSE;
}



void AliAnalysisEtMonteCarloPhos::CreateHistograms()
{ // add some extra histograms & objects to the ones from base class
  if(!fSelector){
    cout<<__FILE__<<" "<<"Creating new fSelector"<<endl;
    fSelector = new AliAnalysisEtSelectorPhos(fCuts);
  }
  AliAnalysisEtMonteCarlo::CreateHistograms();
}
