#include "AliEventContainer.h"

class AliEventContainer;    // your analysis class

ClassImp(AliEventContainer) // classimp: necessary for root

//_____________________________________________________________________________
AliEventContainer::AliEventContainer() :
TObject                     (), 
fV0MPercCorr                (-1),
fV0MPercentile              (-1),
fV0APercentile              (-1),
fV0CPercentile              (-1),

fSPDTrackletsPerc           (-1),
fSPDTrackletsInEta08Perc    (-1),
fSPDTrackletsInEta05Perc    (-1),

fSPDTracklets               (-1),
fSPDTrackletsInEta08        (-1),
fSPDTrackletsInEta05        (-1),

fV0ASignal                  (0.),
fV0CSignal                  (0.),

fIsCollisionCandidate       (kFALSE),
fTriggerMask                (0), 

//Event info
fMagField                   (0.),
fRunNumber                  (0.),
fEventNumber                (0.),

fHasVertex                  (kFALSE),
fProximityCut               (kFALSE),
fHasSPDANDTrkVtx            (kFALSE),
fVertexSelected2015pp       (kFALSE),

//---------------PILE-UP-INFO---------------------
fIsSPDClusterVsTrackletBG   (kFALSE),
fIsIncompleteDAQ            (kFALSE),
fIsINELgtZERO               (kFALSE),
fHasNoInconsistentVertices  (kFALSE),
fIsNotAsymmetricInVZERO     (kFALSE),
fHasGoodVertex2016			(kFALSE),
//---------------PILE-UP-INFO---------------------
fIsPileupMV                 (kFALSE),
fIsPileupOOB                (kFALSE),
fIsPileupFromSPD            (kFALSE),
fIsPileupFromSPDInMultBins  (kFALSE),
fIsPileupInMultBins         (kFALSE) 
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
    fPrimVertex     [0] = -999.;
    fPrimVertex     [1] = -999.;
    fPrimVertex     [2] = -999.;
    
    fPrimVertexCov  [0] = -1;
    fPrimVertexCov  [1] = -1;
    fPrimVertexCov  [2] = -1;
    fPrimVertexCov  [3] = -1;
    fPrimVertexCov  [4] = -1;
    fPrimVertexCov  [5] = -1;
}

//_____________________________________________________________________________
AliEventContainer::AliEventContainer(const AliEventContainer &source): // copy constructor 
TObject(source)
{
    // copy constructor
    fV0MPercCorr                = source.fV0MPercCorr;              
    fV0MPercentile              = source.fV0MPercentile;            
    fV0APercentile              = source.fV0APercentile;            
    fV0CPercentile              = source.fV0CPercentile;            

    fSPDTrackletsPerc           = source.fSPDTrackletsPerc;            
    fSPDTrackletsInEta08Perc    = source.fSPDTrackletsInEta08Perc;     
    fSPDTrackletsInEta05Perc    = source.fSPDTrackletsInEta05Perc;     
    
    fSPDTracklets               = source.fSPDTracklets;            
    fSPDTrackletsInEta08        = source.fSPDTrackletsInEta08;     
    fSPDTrackletsInEta05        = source.fSPDTrackletsInEta05;     

    fV0ASignal                  = source.fV0ASignal;               
    fV0CSignal                  = source.fV0CSignal;         

    fIsCollisionCandidate       = source.fIsCollisionCandidate;
    fTriggerMask                = source.fTriggerMask;             

    //Event info
    fMagField                   = source.fMagField;                
    fRunNumber                  = source.fRunNumber;
    fEventNumber                = source.fEventNumber;
    
    fHasVertex                  = source.fHasVertex;
    fProximityCut               = source.fProximityCut;
    fHasSPDANDTrkVtx            = source.fHasSPDANDTrkVtx;
    fVertexSelected2015pp       = source.fVertexSelected2015pp;

    //------------------------------------------------
    fPrimVertex     [0]         = source.fPrimVertex    [0];
    fPrimVertex     [1]         = source.fPrimVertex    [1];             
    fPrimVertex     [2]         = source.fPrimVertex    [2]; 
    
    fPrimVertexCov  [0]         = source.fPrimVertexCov [0];
    fPrimVertexCov  [1]         = source.fPrimVertexCov [1];             
    fPrimVertexCov  [2]         = source.fPrimVertexCov [2]; 
    fPrimVertexCov  [3]         = source.fPrimVertexCov [3];
    fPrimVertexCov  [4]         = source.fPrimVertexCov [4];             
    fPrimVertexCov  [5]         = source.fPrimVertexCov [5]; 
    
    //---------------PILE-UP-INFO---------------------
    fIsSPDClusterVsTrackletBG   = source.fIsSPDClusterVsTrackletBG;
    fIsIncompleteDAQ            = source.fIsIncompleteDAQ;         
    fIsINELgtZERO               = source.fIsINELgtZERO;            
    fHasNoInconsistentVertices  = source.fHasNoInconsistentVertices;
    fIsNotAsymmetricInVZERO     = source.fIsNotAsymmetricInVZERO;   
	fHasGoodVertex2016			= source.fHasGoodVertex2016;
    //---------------PILE-UP-INFO---------------------
    fIsPileupMV                 = source.fIsPileupMV;               
    fIsPileupOOB                = source.fIsPileupOOB;              
    fIsPileupFromSPD            = source.fIsPileupFromSPD;          
    fIsPileupFromSPDInMultBins  = source.fIsPileupFromSPDInMultBins;
    fIsPileupInMultBins         = source.fIsPileupInMultBins;       
}
//_____________________________________________________________________________
AliEventContainer::~AliEventContainer() // destructor
{
    // destructor
}
//_____________________________________________________________________________
void AliEventContainer::Reset()
{
    fV0MPercCorr                = -1;
    fV0MPercentile              = -1;
    fV0APercentile              = -1;
    fV0CPercentile              = -1;

    fSPDTrackletsPerc           = -1;
    fSPDTrackletsInEta08Perc    = -1;
    fSPDTrackletsInEta05Perc    = -1;
    
    fSPDTracklets               = -1;
    fSPDTrackletsInEta08        = -1;
    fSPDTrackletsInEta05        = -1;

    fV0ASignal                  = 0.;
    fV0CSignal                  = 0.;

    fIsCollisionCandidate       = kFALSE;
    fTriggerMask                = 0; 

    //Event info
    fMagField                   = 0.;
    fRunNumber                  = 0.;
    fEventNumber                = 0.;
    
    fHasVertex                  = kFALSE;
    fProximityCut               = kFALSE;
    fHasSPDANDTrkVtx            = kFALSE;
    fVertexSelected2015pp       = kFALSE;

    //------------------------------------------------
    fPrimVertex     [0]         = -999;
    fPrimVertex     [1]         = -999;
    fPrimVertex     [2]         = -999; 
    
    fPrimVertexCov  [0]         = -1;
    fPrimVertexCov  [1]         = -1;
    fPrimVertexCov  [2]         = -1;
    fPrimVertexCov  [3]         = -1;
    fPrimVertexCov  [4]         = -1;
    fPrimVertexCov  [5]         = -1;
    
    //---------------PILE-UP-INFO---------------------
    fIsSPDClusterVsTrackletBG   = kFALSE;
    fIsIncompleteDAQ            = kFALSE;
    fIsINELgtZERO               = kFALSE;
    fHasNoInconsistentVertices  = kFALSE;
    fIsNotAsymmetricInVZERO     = kFALSE;
	fHasGoodVertex2016			= kFALSE;
    //---------------PILE-UP-INFO---------------------
    fIsPileupMV                 = kFALSE;
    fIsPileupOOB                = kFALSE; 
    fIsPileupFromSPD            = kFALSE; 
    fIsPileupFromSPDInMultBins  = kFALSE; 
    fIsPileupInMultBins         = kFALSE;
}
//_____________________________________________________________________________
void AliEventContainer::Fill(AliMultSelection* MultSelection, AliAnalysisUtils* Utils, AliAODEvent *AODEvent) 
{
    //--------------SPD-MULTIPLICITY-----------------
    fSPDTrackletsPerc           = MultSelection ? MultSelection->GetMultiplicityPercentile("SPDTracklets")  : -1.;
    fSPDTrackletsInEta08Perc    = MultSelection ? MultSelection->GetMultiplicityPercentile("RefMult08")     : -1.;
    fSPDTrackletsInEta05Perc    = MultSelection ? MultSelection->GetMultiplicityPercentile("RefMult05")     : -1.;
    fSPDTracklets               = MultSelection ? MultSelection->GetEstimator("SPDTracklets")->GetValue()   : -1.;
    fSPDTrackletsInEta08        = MultSelection ? MultSelection->GetEstimator("RefMult08")->GetValue()      : -1.;
    fSPDTrackletsInEta05        = MultSelection ? MultSelection->GetEstimator("RefMult05")->GetValue()      : -1.;
    
    //--------------VZERO-MULTIPLICITY----------------
    fV0MPercentile          = MultSelection ? MultSelection->GetMultiplicityPercentile("V0M")           : -1.;
    fV0APercentile          = MultSelection ? MultSelection->GetMultiplicityPercentile("V0A")           : -1.;
    fV0CPercentile          = MultSelection ? MultSelection->GetMultiplicityPercentile("V0C")           : -1.;
    fV0ASignal              = MultSelection ? MultSelection->GetEstimator("V0A")->GetValue()            :  0.;
    fV0CSignal              = MultSelection ? MultSelection->GetEstimator("V0C")->GetValue()            :  0.;

    //-----------------EVENT-INFO---------------------
    fMagField               = AODEvent->GetMagneticField(); 
    fRunNumber              = AODEvent->GetRunNumber();     
    fEventNumber            = ( ( ((ULong64_t)AODEvent->GetPeriodNumber() ) << 36 ) |
                                ( ((ULong64_t)AODEvent->GetOrbitNumber () ) << 12 ) |
                                  ((ULong64_t)AODEvent->GetBunchCrossNumber() )  );
    
    
    const AliAODVertex *vertex = AODEvent->GetPrimaryVertexTracks();
    if (vertex->GetNContributors() < 1) {
        vertex = AODEvent->GetPrimaryVertexSPD();
        if (vertex->GetNContributors() < 1) fHasVertex = kFALSE;
        else fHasVertex = kTRUE;
        TString vtxTyp = vertex->GetTitle();
        Double_t cov[6]={0};
        vertex->GetCovarianceMatrix(cov);
        Double_t zRes = TMath::Sqrt(cov[5]);
        if (vtxTyp.Contains("vertexer:Z") && (zRes>0.25)) fHasVertex = kFALSE;
    }
    else fHasVertex = kTRUE;
    
    fHasSPDANDTrkVtx = kTRUE;
    fProximityCut    = kTRUE;
    fVertexSelected2015pp = SelectVertex2015pp(AODEvent,kTRUE,&fHasSPDANDTrkVtx,&fProximityCut);
    
    //fIsSPDClusterVsTrackletBG   = Utils->IsSPDClusterVsTrackletBG(AODEvent);
	//fIsIncompleteDAQ            = AODEvent->IsIncompleteDAQ();         
	fIsSPDClusterVsTrackletBG   = !MultSelection->GetThisEventPassesTrackletVsCluster();          
    fIsIncompleteDAQ            = !MultSelection->GetThisEventIsNotIncompleteDAQ();  
    
    fIsINELgtZERO               = MultSelection->GetThisEventINELgtZERO();   
    fHasNoInconsistentVertices  = MultSelection->GetThisEventHasNoInconsistentVertices();
    fIsNotAsymmetricInVZERO     = MultSelection->GetThisEventIsNotAsymmetricInVZERO();   
	fHasGoodVertex2016			= MultSelection->GetThisEventHasGoodVertex2016();
    
    fIsPileupMV                 = !MultSelection->GetThisEventIsNotPileupMV();      
    fIsPileupOOB                = Utils->IsOutOfBunchPileUp(AODEvent);     
	//fIsPileupFromSPD			= AODevent->IsPileupFromSPD();       
    fIsPileupFromSPD            = !MultSelection->GetThisEventIsNotPileup(); 
    fIsPileupFromSPDInMultBins  = AODEvent->IsPileupFromSPDInMultBins(); 
    fIsPileupInMultBins         = !MultSelection->GetThisEventIsNotPileupInMultBins();    
}
//_____________________________________________________________________________
Bool_t AliEventContainer::SelectVertex2015pp(AliAODEvent *aod,  Bool_t checkSPDres, Bool_t *SPDandTrkExists, Bool_t *PassProximityCut) 
{
  if (!aod) return kFALSE;
  const AliAODVertex * trkVertex = aod->GetPrimaryVertexTracks();
  const AliAODVertex * spdVertex = aod->GetPrimaryVertexSPD();
//   Bool_t hasSPD = spdVertex->GetStatus();
//   Bool_t hasTrk = trkVertex->GetStatus();
  Bool_t hasSPD = GetVertexStatus(spdVertex);
  Bool_t hasTrk = GetVertexStatus(trkVertex);
  //Note that AliVertex::GetStatus checks that N_contributors is > 0
  //reject events if both are explicitly requested and none is available
  //MOD: do not reject if SPD&Trk vtx. not there, but store it to the variable, if requested:
  if(SPDandTrkExists)
    (*SPDandTrkExists) = hasSPD && hasTrk;
  //Set initial value for proximity check. Only checking the proximity if SPDandTrkExists,
  //the default value should be 1, tracking vtx is not required
  if(PassProximityCut)
    (*PassProximityCut) = 1;
  
  //reject events if none between the SPD or track verteces are available
  //if no trk vertex, try to fall back to SPD vertex;
  if (!hasTrk) {
    if (!hasSPD) return kFALSE;
    //on demand check the spd vertex resolution and reject if not satisfied
    if (checkSPDres && !IsGoodSPDvertexRes(spdVertex)) return kFALSE;
  } else {
    if (hasSPD) {
      //if enabled check the spd vertex resolution and reject if not satisfied
      //if enabled, check the proximity between the spd vertex and trak vertex, and reject if not satisfied
      if (checkSPDres && !IsGoodSPDvertexRes(spdVertex)) return kFALSE;
      if(PassProximityCut)
	(*PassProximityCut) = TMath::Abs(spdVertex->GetZ() - trkVertex->GetZ())<=0.5;
  //if ((checkProximity && TMath::Abs(spdVertex->GetZ() - trkVertex->GetZ())>0.5)) return kFALSE; 
    }
  }

  //Not needed here, done separately
  //Cut on the vertex z position
  //const AliESDVertex * vertex = esd->GetPrimaryVertex();
  //if (TMath::Abs(vertex->GetZ())>10) return kFALSE;
  return kTRUE;
};
//_____________________________________________________________________________
Bool_t AliEventContainer::IsGoodSPDvertexRes(const AliAODVertex * spdVertex)
{
  if (!spdVertex) return kFALSE;
  Double_t covmatrix[6];
  spdVertex->GetCovarianceMatrix(covmatrix);
  if (spdVertex->IsFromVertexerZ() && !(/*spdVertex->GetDispersion()<0.04 &&*/ TMath::Sqrt(covmatrix[5]) <0.25)) return kFALSE;
  return kTRUE;
};
//_____________________________________________________________________________
Bool_t AliEventContainer::GetVertexStatus(const AliAODVertex *vertex)
{
    TString title = vertex->GetTitle();
    if(vertex->GetNContributors()>0 || (title.Contains("cosmics") && !title.Contains("failed"))) return 1;
    if(title.Contains("smearMC")) return 1;
    return 0;
};
