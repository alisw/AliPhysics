/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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

// Interface for all the histogram container classes in JCORRAN analysis

#include "AliJHistogramInterface.h"
#include "AliJCard.h"

/*
 * Default constructor
 */
AliJHistogramInterface::AliJHistogramInterface() :
  fCard(0x0),
  fHMG(NULL),
  fCentBin(),
  fVtxBin(),
  fPTtBin(),
  fPTaBin(),
  fXEBin(),
  fKLongBin(),
  fRGapBin(),
  fEtaGapBin(),
  fPhiGapBin(),
  fMassBin(),
  fTypBin(),
  fTypBin3(),
  fPairPtBin(),
  fhTrackSelection(),
  fhEvents(),
  fhEventTrigger(),
  fhZVertRaw(),
  fhVertexZTriggVtx()
{   // default constructor

}

/*
 * Constructor
 *
 *  AliJCard *cardP = JCard to be initialized for the histogram interface
 */
AliJHistogramInterface::AliJHistogramInterface(AliJCard *cardP) :
  fCard(cardP),
  fHMG(NULL),
  fCentBin(),
  fVtxBin(),
  fPTtBin(),
  fPTaBin(),
  fXEBin(),
  fKLongBin(),
  fRGapBin(),
  fEtaGapBin(),
  fPhiGapBin(),
  fMassBin(),
  fTypBin(),
  fTypBin3(),
  fPairPtBin(),
  fhTrackSelection(),
  fhEvents(),
  fhEventTrigger(),
  fhZVertRaw(),
  fhVertexZTriggVtx()
{   // constructor

    fHMG = new AliJHistManager( "HistManager","AliJHistos");
  
    fCentBin   .Set("Cent",   "C", "Cend:%2.0f-%2.0f%%" ).SetBin( fCard->GetVector("CentBinBorders"));
    fVtxBin    .Set("Vtx",    "V", "Vtx:%2.0f-%2.0f" ).SetBin(fCard->GetVector("zVertBins"));
    fPTtBin    .Set("PTt",    "T", "%.2f<p_{Tt}<%.2f").SetBin(fCard->GetVector("TriggPtBorders"));
    fPTaBin    .Set("PTa",    "A", "%.2f<p_{Ta}<%.2f").SetBin(fCard->GetVector("AssocPtBorders"));
    fXEBin     .Set("XE",     "X", "%.1f<x_{E}<%.1f" ).SetBin(fCard->GetVector("xEBorders"));
    fKLongBin  .Set("KLong", "L",  "%.1f<k_{#parallel}<%.1f").SetBin(fCard->GetVector("KlongBorders"));
    fRGapBin   .Set("RGap",  "R",  "%.1f<R_{gap}<%.1f").SetBin(fCard->GetVector("RGapThresholds"));
    fEtaGapBin .Set("EtaGap", "E", "%.1f<#eta_{gap}<%.1f").SetBin(fCard->GetVector("EtaGapThresholds"));
    fPhiGapBin .Set("PhiGap", "P", "%.1f<#phi_{gap}<%.1f" ).SetBin(fCard->GetVector("EtaGapThresholds"));
    fMassBin   .Set("Mass",   "M", "%.1f<M_{jj}<%.1f").SetBin(fCard->GetVector("PairInvariantMassBins"));
    fTypBin    .Set("Type",   "D", "", AliJBin::kSingle ).SetBin( "0 1" );
    fTypBin3   .Set("Type3",   "D", "", AliJBin::kSingle ).SetBin( "0 1 2 3" );
    fPairPtBin .Set("PairPt", "", AliJBin::kSingle ).SetBin( fCard->GetN("UpperPairPtCut") );
  
}

/*
 * Copy contructor
 */
AliJHistogramInterface::AliJHistogramInterface(const AliJHistogramInterface& obj) :
  fCard(obj.fCard),
  fHMG(obj.fHMG),
  fCentBin(obj.fCentBin),
  fVtxBin(obj.fVtxBin),
  fPTtBin(obj.fPTtBin),
  fPTaBin(obj.fPTaBin),
  fXEBin(obj.fXEBin),
  fKLongBin(obj.fKLongBin),
  fRGapBin(obj.fRGapBin),
  fEtaGapBin(obj.fEtaGapBin),
  fPhiGapBin(obj.fPhiGapBin),
  fMassBin(obj.fMassBin),
  fTypBin(obj.fTypBin),
  fTypBin3(obj.fTypBin3),
  fPairPtBin(obj.fPairPtBin),
  fhTrackSelection(obj.fhTrackSelection),
  fhEvents(obj.fhEvents),
  fhEventTrigger(obj.fhEventTrigger),
  fhZVertRaw(obj.fhZVertRaw),
  fhVertexZTriggVtx(obj.fhVertexZTriggVtx)
{   // copy constructor
  
}

/*
 * Destructor
 */
AliJHistogramInterface::~AliJHistogramInterface() {
	delete fHMG;
}

/*
 * Histogram creator for the histograms needed in AliJDataManager.cxx
 */
void AliJHistogramInterface::CreateDataManagerHistograms(){
  
  fHMG->cd();
  
  // z-vertex
  fhZVertRaw << TH1D("hZVertRaw","vertex 0", 120, -30., 30.) << "END";
  
  double   binsVertexMult[] = {0,1,2,3,4,5,10000};
  int   nbinsVertexMult  = sizeof(binsVertexMult)/sizeof(double)-1;
  double binsVertexZ[]    = {-10,-6,-3,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,3,6,10};
  int   nbinsVertexZ   = sizeof(binsVertexZ)/sizeof(double)-1;
  fhVertexZTriggVtx << TH2D("hVertexZTriggVtx","Vertex counts", nbinsVertexMult, binsVertexMult, nbinsVertexZ, binsVertexZ) << "END";
  
  // event counter
  fhEvents
    << TH1D("hEvents","Events passing cuts", 100, -0.5, 100-0.5 ) << "END";
  fhEventTrigger
    << TH1D("hEventTrigger","Trigger count", 50, -0.5, 50.-.5 ) << "END";
  fhTrackSelection
    << TH1D("hTrackSelection","Bit convention", 100, -0.5, 100-0.5) << "END";
  
}
