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

AliFemtoCutValues::AliFemtoCutValues() :
  fEvtZvtxLow(-10.),//cm
  fEvtZvtxUp(10.),//cm
  fV0PIDThresholdPtTPCLow(0.3),//GeV/c
  fV0PIDThresholdPtTPCUp(4.3),//GeV/c
  fV0SelWidth(0.004),//GeV/c*2
  fV0EtaRange(0.8),
  fAntiCutK0sLow(0.48),//GeV/c*2
  fAntiCutK0sUp(0.515),//GeV/c*2
  fV0Nsigma(5.),
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
  fProtonDCAxyCut(0.1),//cm
  fProtonDCAzCut(0.2),//cm
  fProtonEtaRange(0.8),
  fProtonNsigma(3.),
  fProtonPtBinsDCA(20),
  fProtonPtBinsPurity(33),
  fProtonTPCCluster(80),
  fPDGDatabase(new TDatabasePDG())
{
  //Default constructor
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
  fProtonDCAxyCut(obj.fProtonDCAxyCut),
  fProtonDCAzCut(obj.fProtonDCAzCut),
  fProtonEtaRange(obj.fProtonEtaRange),
  fProtonNsigma(obj.fProtonNsigma),
  fProtonPtBinsDCA(obj.fProtonPtBinsDCA),
  fProtonPtBinsPurity(obj.fProtonPtBinsPurity),
  fProtonTPCCluster(obj.fProtonTPCCluster),
  fPDGDatabase(new TDatabasePDG())
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
AliFemtoCutValues::~AliFemtoCutValues()
{
  if(fPDGDatabase)
    {
      delete fPDGDatabase;
      fPDGDatabase = 0;
    }
}
//_____________________________________________________________________________
