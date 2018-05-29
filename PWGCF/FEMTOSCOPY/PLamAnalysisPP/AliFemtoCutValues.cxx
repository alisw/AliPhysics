/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//Class that contains all cut values relevant for baryon-baryon femtoscopy

#include "AliFemtoCutValues.h"

ClassImp(AliFemtoCutValues)

AliFemtoCutValues::AliFemtoCutValues() :
  fEvtZvtxLow(-10.),//cm
  fEvtZvtxUp(10.),//cm
  fV0PIDThresholdPtTPCLow(0.3),//GeV/c
  fV0PIDThresholdPtTPCUp(4.3),//GeV/c
  fV0SelWidth(0.004),//GeV/c*2
  fV0EtaRange(0.8),
  fAntiCutK0sLow(0.48),//GeV/c*2
  fAntiCutK0sUp(0.515),//GeV/c*2
  fV0Nsigma(5.),//5
  fV0Pointing(0.99),
  fV0RxyLow(0.2),//cm
  fV0RxyUp(100.),//cm
  fV0Decayvtx(100.),//cm
  fV0DCAtrackPV(0.05),//cm
  fV0DCAtracksV0decay(1.5),//cm
  fV0PtBinsCPA(8),
  fV0PtBinsInvMass(8),
  fV0TPCCluster(70),
  fProtonPIDThresholdPtTPCLow(0.5),//GeV/c
  fProtonPIDThresholdPtTOFUp(4.05),//GeV/c
  fProtonPIDTPCTOFSwitch(0.75),//GeV/c
  fProtonPtBinsDCA(20),
  fProtonPtBinsPurity(33),
  fProtonTPCCluster(80),
  fProtonDCAxyCut(0.1),//cm
  fProtonDCAzCut(0.2),//cm
  fProtonEtaRange(0.8),
  fProtonNsigma(3.),
  fTrackFilterBit(128)
//fPDGDatabase(new TDatabasePDG())
{
  //Default constructor
  Info("AliFemtoCutValues","Calling default Constructor");
}
//_____________________________________________________________________________
AliFemtoCutValues::AliFemtoCutValues(systematics cutType) :
  fEvtZvtxLow(-10.),//cm
  fEvtZvtxUp(10.),//cm
  fV0PIDThresholdPtTPCLow(0.3),//GeV/c
  fV0PIDThresholdPtTPCUp(4.3),//GeV/c
  fV0SelWidth(0.004),//GeV/c*2
  fV0EtaRange(0.8),
  fAntiCutK0sLow(0.48),//GeV/c*2
  fAntiCutK0sUp(0.515),//GeV/c*2
  fV0Nsigma(5.),//5
  fV0Pointing(0.99),
  fV0RxyLow(0.2),//cm
  fV0RxyUp(100.),//cm
  fV0Decayvtx(100.),//cm
  fV0DCAtrackPV(0.05),//cm
  fV0DCAtracksV0decay(1.5),//cm
  fV0PtBinsCPA(8),
  fV0PtBinsInvMass(8),
  fV0TPCCluster(70),
  fProtonPIDThresholdPtTPCLow(0.5),//GeV/c
  fProtonPIDThresholdPtTOFUp(4.05),//GeV/c
  fProtonPIDTPCTOFSwitch(0.75),//GeV/c
  fProtonPtBinsDCA(20),
  fProtonPtBinsPurity(33),
  fProtonTPCCluster(80),
  fProtonDCAxyCut(0.1),//cm
  fProtonDCAzCut(0.2),//cm
  fProtonEtaRange(0.8),
  fProtonNsigma(3.),
  fTrackFilterBit(128)
//fPDGDatabase(new TDatabasePDG())
{
  if(cutType != kDefault) SetCutVariation(cutType);

  Info("AliFemtoCutValues","Calling extended Constructor");
}
//_____________________________________________________________________________
AliFemtoCutValues::AliFemtoCutValues(const AliFemtoCutValues &obj) :
  fEvtZvtxLow(obj.fEvtZvtxLow),
  fEvtZvtxUp(obj.fEvtZvtxUp),
  fV0PIDThresholdPtTPCLow(obj.fV0PIDThresholdPtTPCLow),
  fV0PIDThresholdPtTPCUp(obj.fV0PIDThresholdPtTPCUp),
  fV0SelWidth(obj.fV0SelWidth),
  fV0EtaRange(obj.fV0EtaRange),
  fAntiCutK0sLow(obj.fAntiCutK0sLow),
  fAntiCutK0sUp(obj.fAntiCutK0sUp),
  fV0Nsigma(obj.fV0Nsigma),
  fV0Pointing(obj.fV0Pointing),
  fV0RxyLow(obj.fV0RxyLow),
  fV0RxyUp(obj.fV0RxyUp),
  fV0Decayvtx(obj.fV0Decayvtx),
  fV0DCAtrackPV(obj.fV0DCAtrackPV),
  fV0DCAtracksV0decay(obj.fV0DCAtracksV0decay),
  fV0PtBinsCPA(obj.fV0PtBinsCPA),
  fV0PtBinsInvMass(obj.fV0PtBinsInvMass),
  fV0TPCCluster(obj.fV0TPCCluster),
  fProtonPIDThresholdPtTPCLow(obj.fProtonPIDThresholdPtTPCLow),
  fProtonPIDThresholdPtTOFUp(obj.fProtonPIDThresholdPtTOFUp),
  fProtonPIDTPCTOFSwitch(obj.fProtonPIDTPCTOFSwitch),
  fProtonPtBinsDCA(obj.fProtonPtBinsDCA),
  fProtonPtBinsPurity(obj.fProtonPtBinsPurity),
  fProtonTPCCluster(obj.fProtonTPCCluster),
  fProtonDCAxyCut(obj.fProtonDCAxyCut),
  fProtonDCAzCut(obj.fProtonDCAzCut),
  fProtonEtaRange(obj.fProtonEtaRange),
  fProtonNsigma(obj.fProtonNsigma)
{
  // copy constructor
}

//_____________________________________________________________________________
AliFemtoCutValues &AliFemtoCutValues::operator=(const AliFemtoCutValues &obj)
{
  //Assignment operator
  if(this == &obj) return *this;

  fEvtZvtxLow = obj.fEvtZvtxLow;
  fEvtZvtxUp = obj.fEvtZvtxUp;
  fV0PIDThresholdPtTPCLow = obj.fV0PIDThresholdPtTPCLow;
  fV0PIDThresholdPtTPCUp = obj.fV0PIDThresholdPtTPCUp;
  fV0SelWidth = obj.fV0SelWidth;
  fV0EtaRange = obj.fV0EtaRange;
  fAntiCutK0sLow = obj.fAntiCutK0sLow;
  fAntiCutK0sUp = obj.fAntiCutK0sUp;
  fV0Nsigma = obj.fV0Nsigma;
  fV0Pointing = obj.fV0Pointing;
  fV0RxyLow = obj.fV0RxyLow;
  fV0RxyUp = obj.fV0RxyUp;
  fV0Decayvtx = obj.fV0Decayvtx;
  fV0DCAtrackPV = obj.fV0DCAtrackPV;
  fV0DCAtracksV0decay = obj.fV0DCAtracksV0decay;
  fV0PtBinsCPA = obj.fV0PtBinsCPA;
  fV0PtBinsInvMass = obj.fV0PtBinsInvMass;
  fV0TPCCluster = obj.fV0TPCCluster;
  fProtonPIDThresholdPtTPCLow = obj.fProtonPIDThresholdPtTPCLow;
  fProtonPIDThresholdPtTOFUp = obj.fProtonPIDThresholdPtTOFUp;
  fProtonPIDTPCTOFSwitch = obj.fProtonPIDTPCTOFSwitch;
  fProtonPtBinsDCA = obj.fProtonPtBinsDCA;
  fProtonPtBinsPurity = obj.fProtonPtBinsPurity;
  fProtonTPCCluster = obj.fProtonTPCCluster;
  fProtonDCAxyCut = obj.fProtonDCAxyCut;
  fProtonDCAzCut = obj.fProtonDCAzCut;
  fProtonEtaRange = obj.fProtonEtaRange;
  fProtonNsigma = obj.fProtonNsigma;

  return (*this);
}
//_____________________________________________________________________________
Int_t AliFemtoCutValues::fFindPtBin(Double_t ptVal,Double_t ptLow,Double_t ptUp,Int_t nBins)
{
  Int_t ptBin = -9999;

  Float_t deltaPt = (ptUp - ptLow)/(Float_t)nBins;

  for(int i=0; i<nBins; i++)
    {
      if((ptLow + i*deltaPt) < ptVal && ptVal < (ptLow + (i+1)*deltaPt))
	{
	  ptBin = i;
	  break;
	}
    }

  return ptBin;
}
//_____________________________________________________________________________
void AliFemtoCutValues::SetCutVariation(systematics cutType)
{
  //This will update the cut values initialized for default values in the constructor

  if(cutType == kDefault) return;//do nothing
  //Proton Variations:
  else if(cutType == kProtonVariationLowerPtThresholdDown) fProtonPIDThresholdPtTPCLow = 0.4;
  else if(cutType == kProtonVariationLowerPtThresholdUp) fProtonPIDThresholdPtTPCLow = 0.6;
  else if(cutType == kProtonVariationEtaRangeUp) fProtonEtaRange = 0.9;
  else if(cutType == kProtonVariationEtaRangeDown) fProtonEtaRange = 0.7;
  else if(cutType == kProtonVariationNsigmaUp) fProtonNsigma = 4.;
  else if(cutType == kProtonVariationNsigmaDown) fProtonNsigma = 2.;
  else if(cutType == kProtonTPCClusterUp) fProtonTPCCluster = 90;
  else if(cutType == kProtonVariationFilterBitGlobal) fTrackFilterBit = 96;
  //V0 variations:
  else if(cutType == kV0VariationLowerPtThresholdDown) fV0PIDThresholdPtTPCLow = 0.24;
  else if(cutType == kV0VariationLowerPtThresholdUp) fV0PIDThresholdPtTPCLow = 0.36;
  else if(cutType == kV0VariationCosinePointingUp) fV0Pointing = 0.998;
  else if(cutType == kV0VariationNsigmaDown) fV0Nsigma = 4.;
  else if(cutType == kV0VariationTPCClusterUp) fV0TPCCluster = 80;
  else if(cutType == kV0VariationEtaRangeUp) fV0EtaRange = 0.9;
  else if(cutType == kV0VariationEtaRangeDown) fV0EtaRange = 0.7;
  else if(cutType == kV0VariationDCAatV0DecayDown) fV0DCAtracksV0decay = 1.2;
  else if(cutType == kV0VariationDCADaughtersToPVUp) fV0DCAtrackPV = 0.06;
  else std::cout << "Unknown variation, default values are used for analysis" << std::endl;
}
//_____________________________________________________________________________
AliFemtoCutValues::~AliFemtoCutValues()
{
  /*
  if(fPDGDatabase)
    {
      delete fPDGDatabase;
      fPDGDatabase = 0;
    }
  */
}
//_____________________________________________________________________________
