/*************************************************************************
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
//
// Class for electrons from beauty study
// Counting electrons from beauty
// by DCA cuts, background subtraction 
//
// Authors:
//   Hongyan Yang <hongyan@physi.uni-heidelberg.de>
//   Carlo Bombonati <Carlo.Bombonati@cern.ch>
//

#include "TMath.h"
#include "TList.h"
#include "AliLog.h"

#include <TParticle.h>
#include <TDatabasePDG.h>
#include "THnSparse.h"

#include "AliMCEvent.h"
#include "AliMCVertex.h"
#include "AliMCParticle.h"
#include "AliStack.h" 

#include "AliESDEvent.h"
#include "AliESDtrack.h"

#include "AliKFParticle.h"
#include "AliKFVertex.h"

#include "AliVertexerTracks.h"


#include "AliHFEdisplacedElectrons.h"

ClassImp(AliHFEdisplacedElectrons)

//__________________________________________________________
const Float_t AliHFEdisplacedElectrons::fgkDcaMinPtIntv[13] = {
  // 
  // define DCA low limits for single electrons cut in each Pt bin 
  // preliminary numbers 
  // these numbers should be replaced by the best numbers determined by 
  // cut efficiency study: signal/background  (not used right now)
  //
  450, // 0.0 - 0.5
  450, // 0.5 - 1.0
  450, // 1.0 - 1.5
  500, // 1.5 - 2.0
  400, // 2.0 - 2.5
  300, // 2.5 - 3.0
  300, // 3.0 - 4.0
  300, // 4.0 - 5.0
  200, // 5.0 - 7.0
  150, // 7.0 - 9.0
  150, // 9.0 - 12.0
  100, //12.0 - 16.0
  50}; //16.0 - 20.0
//__________________________________________________________
const Float_t AliHFEdisplacedElectrons::fgkPtIntv[14] = {
  //
  // define pT bins for spectra of single electrons 
  //
  0.0,0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0,7.0,9.0,12.0,16.0,20.0};

//__________________________________________________________
const Char_t *AliHFEdisplacedElectrons::fgkKineVar[3] = {
  "y", "phi", "pt"};

//__________________________________________________________
const Char_t *AliHFEdisplacedElectrons::fgkKineVarTitle[3] ={
  "rapidity;y;dN/dy",   "azimuthal;#phi;dN/d#phi",   "transverse momentum;p_{T};dN/dp_{T}", 
};

//__________________________________________________________
AliHFEdisplacedElectrons::AliHFEdisplacedElectrons():
  fDeDebugLevel(0)
  , fNclustersITS(0)
  , fMinNprimVtxContributor(1)
  , fTHnSparseDcaMcEleInfo(NULL)
  , fTHnSparseDcaEsdEleInfo(NULL)
  , fTHnSparseDcaDataEleInfo(NULL)
  , fDeOutputList(0x0)
{
  //
  // default constructor
  //

}

//__________________________________________________________
AliHFEdisplacedElectrons::AliHFEdisplacedElectrons(const AliHFEdisplacedElectrons &ref):
  TObject(ref)
  , fDeDebugLevel(ref.fDeDebugLevel)
  , fNclustersITS(ref.fNclustersITS)
  , fMinNprimVtxContributor(ref.fMinNprimVtxContributor)
  , fTHnSparseDcaMcEleInfo(ref.fTHnSparseDcaMcEleInfo)
  , fTHnSparseDcaEsdEleInfo(ref.fTHnSparseDcaEsdEleInfo)
  , fTHnSparseDcaDataEleInfo(ref.fTHnSparseDcaDataEleInfo)
  , fDeOutputList(ref.fDeOutputList)
  
{
  //
  // copy constructor
  //
}


//__________________________________________________________
AliHFEdisplacedElectrons&AliHFEdisplacedElectrons::operator=(const AliHFEdisplacedElectrons &ref)
{
  //
  // Assignment operator
  //

  if(this == &ref) return *this;
  AliHFEdisplacedElectrons::operator=(ref);

  fDeDebugLevel = ref.fDeDebugLevel;
  fNclustersITS = ref.fNclustersITS;
  fMinNprimVtxContributor = ref.fMinNprimVtxContributor;
  fTHnSparseDcaMcEleInfo = ref.fTHnSparseDcaMcEleInfo;
  fTHnSparseDcaEsdEleInfo = ref.fTHnSparseDcaEsdEleInfo;
  fTHnSparseDcaDataEleInfo = ref.fTHnSparseDcaDataEleInfo;
  fDeOutputList = ref.fDeOutputList;

  return *this;
}

//__________________________________________________________
AliHFEdisplacedElectrons::~AliHFEdisplacedElectrons()
{
  //
  // default constructor
  //

  if(fTHnSparseDcaMcEleInfo) 
    delete fTHnSparseDcaMcEleInfo;
  if(fTHnSparseDcaEsdEleInfo) 
    delete fTHnSparseDcaEsdEleInfo;
  if(fTHnSparseDcaDataEleInfo) 
    delete fTHnSparseDcaDataEleInfo;
    
  if(fDeOutputList){
    fDeOutputList->Clear();
    delete fDeOutputList;
    
  }

  Printf("analysis done\n");
}

//__________________________________________________________
void AliHFEdisplacedElectrons::InitAnalysis(){
  //
  // init analysis (no intialization yet)
  //
  

  Printf("initialize analysis\n");

}


//__________________________________________________________
void AliHFEdisplacedElectrons::CreateOutputs(TList* const displacedList){

  //
  //  create output fDeOutputList
  //

  // THnSparseF
  // 8+? interested electron sources: others, photon conv, direct photon,  pi0, eta, b, b->c, c + pion or missid, and missid pion
  // 41 possible DCA cuts XY: 
  // 13 pT bins
  // 10 rapidity bins
  // 10 azimuthal angle phi bins


  if(!displacedList) return;

  fDeOutputList = displacedList;
  fDeOutputList -> SetName("displacedElectrons");

  // electron source 
  Int_t nBinsEleSource = 10;
  Double_t minEleSource = -1.5; 
  Double_t maxEleSource = 8.5;  
  Double_t *binLimEleSource = new Double_t[nBinsEleSource+1];
  for(Int_t i=0; i<=nBinsEleSource; i++)
    binLimEleSource[i] = minEleSource + i*(maxEleSource-minEleSource)/nBinsEleSource;

  // dca bins:  XY 
  Int_t nBinsDcaXY = kNDcaMin-1;  // 41 bins
  Double_t minDcaXY = -20.5; 
  Double_t maxDcaXY = 20.5; 
  Double_t dcaXYBinWidth = (maxDcaXY-minDcaXY)/nBinsDcaXY;
  Double_t *binLimDcaXY = new Double_t[nBinsDcaXY+1];
  for(Int_t i=0; i<=nBinsDcaXY; i++)
    binLimDcaXY[i] = minDcaXY + i*dcaXYBinWidth;


  // dca bins:  Z
  Int_t nBinsDcaZ = kNDcaMin/2;  // 21 bins
  Double_t minDcaZ = -10.5; 
  Double_t maxDcaZ = 10.5; 
  Double_t dcaZBinWidth = (maxDcaZ-minDcaZ)/nBinsDcaZ;
  Double_t *binLimDcaZ = new Double_t[nBinsDcaZ+1];
  for(Int_t i=0; i<=nBinsDcaZ; i++)
    binLimDcaZ[i] = minDcaZ + i*dcaZBinWidth;
  
 // pt bins
  Int_t nBinsPt = kNPtIntv-1;
  Double_t *binLimPt = new Double_t[nBinsPt+1];
  for(Int_t i=0; i<=nBinsPt; i++)
    binLimPt[i] = fgkPtIntv[i];  // variable bins

  // rapidity bins
  Int_t nBinsRap = 10;
  Double_t minRap = -1.0; 
  Double_t maxRap = 1.0;
  Double_t *binLimRap = new Double_t[nBinsRap+1];
  for(Int_t i=0; i<=nBinsRap; i++)
    binLimRap[i] = minRap + i*(maxRap-minRap)/nBinsRap;

  // azumuthal phi angle
  Int_t nBinsPhi = 10;
  Double_t minPhi = 0; 
  Double_t maxPhi = 2*TMath::Pi();
  Double_t *binLimPhi = new Double_t[nBinsPhi+1];
  for(Int_t i=0; i<=nBinsPhi; i++)
    binLimPhi[i] = minPhi + i*(maxPhi-minPhi)/nBinsPhi;

  // production radii
  Int_t nBinsR = 30;
  Double_t minR = 0; 
  Double_t maxR = 15;
  Double_t *binLimR = new Double_t[nBinsR+1];
  for(Int_t i=0; i<=nBinsR; i++)
    binLimR[i] = minR + i*(maxR-minR)/nBinsR;

  // status
  Int_t nBinsStat = 11;
  Double_t minStat = 0; 
  Double_t maxStat = 11;
  Double_t *binLimStat = new Double_t[nBinsStat+1];
  for(Int_t i=0; i<=nBinsStat; i++)
    binLimStat[i] = minStat + i*(maxStat-minStat)/nBinsStat;

  fTHnSparseDcaMcEleInfo = 0x0;
  fTHnSparseDcaEsdEleInfo = 0x0;
  fTHnSparseDcaDataEleInfo = 0x0;

  // for MC only: MC electron ID
  const Int_t nVarMc = 8;
  Int_t iBinMc[nVarMc] = {nBinsEleSource, nBinsDcaXY, nBinsDcaZ, nBinsPt, nBinsRap, nBinsPhi, nBinsR, nBinsStat};
  
  fTHnSparseDcaMcEleInfo = new THnSparseF("dcaMcElectronInfo", 
					   "MC electrons;electron source ID;mc dcaXY [50 #mum];mc dcaZ [100 #mum];mc pT [GeV/c];mc y [rapidity];mc #phi [rad];mc R [cm];Status Code", 	
					   nVarMc, iBinMc);

  fTHnSparseDcaMcEleInfo->SetBinEdges(0, binLimEleSource); // electron source
  fTHnSparseDcaMcEleInfo->SetBinEdges(1, binLimDcaXY);  // dca xy cut
  fTHnSparseDcaMcEleInfo->SetBinEdges(2, binLimDcaZ);  // dca z cut
  fTHnSparseDcaMcEleInfo->SetBinEdges(3, binLimPt);   // pt
  fTHnSparseDcaMcEleInfo->SetBinEdges(4, binLimRap);  // rapidity
  fTHnSparseDcaMcEleInfo->SetBinEdges(5, binLimPhi);  // phi
  fTHnSparseDcaMcEleInfo->SetBinEdges(6, binLimR);
  fTHnSparseDcaMcEleInfo->SetBinEdges(7, binLimStat);
  fTHnSparseDcaMcEleInfo->Sumw2();  

  // for ESD with MC: HFE pid and MC pid
  const Int_t nVarEsd = 9;
  Int_t iBin[nVarEsd] = {nBinsEleSource,nBinsDcaXY, nBinsDcaZ, nBinsDcaXY, nBinsPt, nBinsRap, nBinsPhi, nBinsR, nBinsStat};
  
  fTHnSparseDcaEsdEleInfo 
    = new THnSparseF("dcaEsdElectronInfo", 
		     "ESD electrons;electron source ID;esd dcaXY [50 #mum];esd dcaZ [100 #mum];esd dcaXY KF [50 #mum];pT [GeV/c];y [rapidity];#phi [rad];R [cm];Status Code",		
		     nVarEsd, iBin);

  fTHnSparseDcaEsdEleInfo->SetBinEdges(0, binLimEleSource); // electron source
  fTHnSparseDcaEsdEleInfo->SetBinEdges(1, binLimDcaXY);  // dca xy without current track
  fTHnSparseDcaEsdEleInfo->SetBinEdges(2, binLimDcaZ);  // dca z without current track
  fTHnSparseDcaEsdEleInfo->SetBinEdges(3, binLimDcaXY);  // dca xy kf without current track
  fTHnSparseDcaEsdEleInfo->SetBinEdges(4, binLimPt);   // pt esd
  fTHnSparseDcaEsdEleInfo->SetBinEdges(5, binLimRap);  // rapidity esd
  fTHnSparseDcaEsdEleInfo->SetBinEdges(6, binLimPhi);  // phi esd
  fTHnSparseDcaEsdEleInfo->SetBinEdges(7, binLimR);
  fTHnSparseDcaEsdEleInfo->SetBinEdges(8, binLimStat);
  fTHnSparseDcaEsdEleInfo->Sumw2();
  

  // for ESD data: HFE pid
  const Int_t nVarData = 6;
  Int_t iBinData[nVarData] = {nBinsDcaXY, nBinsDcaZ, nBinsDcaXY, nBinsPt, nBinsRap, nBinsPhi};  
  
  fTHnSparseDcaDataEleInfo 
    = new THnSparseF("dcaDataElectronInfo", 
		     "Data electrons;dcaXY [50 #mum];dcaZ [100 #mum];dcaXYkf [50 #mum];pT [GeV/c];y [rapidity];#phi [rad];",
		     nVarData, iBinData);    
  fTHnSparseDcaDataEleInfo->SetBinEdges(0, binLimDcaXY);  // dca xy cut w/o
  fTHnSparseDcaDataEleInfo->SetBinEdges(1, binLimDcaZ);  // dca z cut w/o
  fTHnSparseDcaDataEleInfo->SetBinEdges(2, binLimDcaXY);  // dca xy kf cut 
  fTHnSparseDcaDataEleInfo->SetBinEdges(3, binLimPt);   // pt
  fTHnSparseDcaDataEleInfo->SetBinEdges(4, binLimRap);  // rapidity
  fTHnSparseDcaDataEleInfo->SetBinEdges(5, binLimPhi);  // phi
  fTHnSparseDcaDataEleInfo->Sumw2();

  fDeOutputList -> AddAt(fTHnSparseDcaMcEleInfo, kMcElectron);
  fDeOutputList -> AddAt(fTHnSparseDcaEsdEleInfo, kEsdElectron);
  fDeOutputList -> AddAt(fTHnSparseDcaDataEleInfo, kDataElectron);
  
  AliInfo("THnSparse histograms are created\n");
  fDeOutputList->Print();

}


//__________________________________________________________
void AliHFEdisplacedElectrons::FillMcOutput(const AliESDEvent *const fESD, AliMCEvent* const fMC, const AliMCParticle* const mctrack)
{

  // fill output
  //0. after mc event cut
  //1. after checking stack, mcpart etc are valid
  //2. after PID
  //3. after event and track selection 

  AliStack *stack = fMC->Stack();
  TParticle *part = mctrack->Particle();

  // obtain impact parameters in xy and z
  AliMCVertex *mcPrimVtx = (AliMCVertex *)fMC->GetPrimaryVertex();      
  Double_t mcPrimV[3];
  mcPrimV[0] = mcPrimVtx->GetX();
  mcPrimV[1] = mcPrimVtx->GetY();
  mcPrimV[2] = mcPrimVtx->GetZ();

  Double_t mcVtxXY = TMath::Abs(mcPrimV[0]*mcPrimV[0] + mcPrimV[1]*mcPrimV[1]);

  Float_t vx = part->Vx();  // in cm
  Float_t vy = part->Vy();  // in cm
  Float_t vz = part->Vz();   // in cm
  
  Float_t vxy = TMath::Sqrt(vx*vx+vy*vy);
  
  Float_t mcpx = part->Px();
  Float_t mcpy = part->Py();
  Float_t mcpt = TMath::Sqrt(mcpx*mcpx+mcpy*mcpy);
  Float_t mceta = part->Eta();
  Float_t mcphi = part->Phi();

  // begin new stuff
  Double_t mcR = part->R();
  Int_t mcStatus = part->GetStatusCode();
  // end new staff

  Int_t pdg = part->GetPdgCode();
  
  Int_t charge = 1;
  if(pdg==kPDGelectron || pdg==-kPDGpion) charge = -1;  
  
  // calculate mcDca ------------------------------------------------------------------ 
  const Float_t conv[2] = {1.783/1.6, 2.99792458};
  Float_t magneticField = 0;  // initialized as 5kG
  magneticField = fESD->GetMagneticField();  // in kG
  Float_t radiusMc = mcpt/(TMath::Abs(magneticField)/10.)*conv[0]*conv[1]; // pt in GeV/c, magnetic field in Tesla, radius in meter
  
  Float_t radius;
  radius = TMath::Abs(radiusMc);
  Float_t nx = mcpx/mcpt;
  Float_t ny = mcpy/mcpt;
  Double_t dxy = vxy - mcVtxXY;   // in cm
  Double_t dvx = vx - mcPrimV[0]; // in cm
  Double_t dvy = vy - mcPrimV[1]; // in cm
  Float_t mcDcaXYm = (radius - TMath::Sqrt(dxy*dxy/100./100. + radius*radius + 2*radius*charge*(dvx*ny-dvy*nx)/100.)) ;  // in meters
  Double_t mcDca[2] = {mcDcaXYm*100, vz};  // in cm
  Double_t mcDcaXY = mcDca[0]*1.0e4;  // conv dca in cm to dca in micron 
  Double_t mcDcaZ = mcDca[1]*1.0e4;
    
  const Int_t nvarMC=8;
  Double_t var[nvarMC];
  var[0] = -1;
  
  if(TMath::Abs(pdg)==kPDGelectron) {

    Int_t eleLabel = mctrack->GetLabel();
    Int_t sourcePdg = ElectronFromSource(stack, eleLabel);     
    
    if(sourcePdg==kPDGgamma){ // check direct photon or not 
      if(ElePhotonDirect(stack, eleLabel)!=-1)	
	var[0] = kEleDirectPhotonConv;    
      else 	
	var[0] = kElePhotonConv;     
      if(fDeDebugLevel>=10) 
	printf("photonconv=====this electron %d is from %d, source id=%.1f\n", eleLabel, sourcePdg, var[0]);   
    }   // photon or direct photon -> e
    
    if(sourcePdg==kPDGpi0){      
      var[0] = kElePi0;   
      if(fDeDebugLevel>=10) printf("pi0======this electron %d is from %d, source id=%.1f\n", eleLabel, sourcePdg, var[0]);    
    } 
    
    if(sourcePdg==kPDGeta){ 
      var[0] = kEleEta;     
      if(fDeDebugLevel>=10) 
	printf("eta=====this electron %d is from %d, source id=%.1f\n", eleLabel, sourcePdg, var[0]);    
    }   // from eta -> e
    
    if(TMath::Abs(sourcePdg%10000)/100==kPDGbeauty ||    // for special intermediate meson states: like 100553       
       TMath::Abs(sourcePdg)/1000==kPDGbeauty ||      
       TMath::Abs(sourcePdg)/100==kPDGbeauty ||     
       TMath::Abs(sourcePdg)==kPDGbeauty){     
      var[0]=kEleB;       
      if(fDeDebugLevel>=10) 
	printf("beauty======electron %d is from %d with id %.1f\n", eleLabel, ElectronFromSource(stack, eleLabel), var[0]);   
    } // direct beauty  -> e      
    if(TMath::Abs(sourcePdg%10000)/100==kPDGcharm ||    // for special intermediate meson states: like 1004**      
       TMath::Abs(sourcePdg)/1000==kPDGcharm ||       
       TMath::Abs(sourcePdg)/100==kPDGcharm ||        
       TMath::Abs(sourcePdg)==kPDGcharm){           
      // two cases: 
      //            electron from b->c->e     
      //            electron from c->e     
      if(ElectronFromCharm(stack, eleLabel)!=-1){
	var[0] = ElectronFromCharm(stack, eleLabel);
	if(fDeDebugLevel>=10)
	  printf("charm----->electron %d has mother %d is from %.1f\n", eleLabel, ElectronFromSource(stack, eleLabel), var[0]);
      } 
    }  // charm electrons: b->c->e or c->e
    
    if(fDeDebugLevel>=10) printf("others----->electron %d has mother %d is from %.1f\n", eleLabel, ElectronFromSource(stack, eleLabel), var[0]);
  } // electron source endif 
  
  else
    if(TMath::Abs(pdg)==kPDGpion)
      var[0] = kPion;
  
  if(TMath::Abs(mcDcaXY)<1000) var[1] = Int_t(mcDcaXY)/50;   // larger than 1mm should go to the last bin
  else
    var[1] = ((mcDcaXY>0)?1:-1)*20;
  
  if(TMath::Abs(mcDcaZ)<1000) var[2] = Int_t(mcDcaZ)/100;
  else
    var[2] = ((mcDcaZ>0)?1:-1)*10;
        
    var[3] = mcpt;      // pt 
    var[4] = mceta;     // eta
    var[5] = mcphi;     // phi    
    var[6] = mcR;       // production radius
    var[7] = mcStatus;  // internal status


    (static_cast<THnSparseF *>(fDeOutputList->At(kMcElectron)))->Fill(var);
}



//__________________________________________________________
void AliHFEdisplacedElectrons::FillEsdOutput(const AliESDEvent * const fESDEvent, AliESDtrack* const esdTrack, AliStack * const stack)
{
  // after esd event selection, esd track cuts, and hfe pid 
  // this is the case for ESD tracks, with MC information

  AliESDtrack *track = esdTrack;
  Double_t pt  = track->Pt();
  Double_t eta = track->Eta();
  Double_t phi = track->Phi();
  
  Int_t eleLabel = track->GetLabel();
  if(eleLabel<0 || eleLabel>stack->GetNtrack()) return;  

  TParticle *particle = stack->Particle(eleLabel);
  if(!particle) return;

  // begin new stuff
  Double_t mcR = particle->R();
  Int_t mcStatus = particle->GetStatusCode();
  // end new staff

  // obtain impact parameters in xy and z
  Float_t magneticField = 5;  // initialized as 5kG
  magneticField = fESDEvent->GetMagneticField();  // in kG 
  Double_t beampiperadius=3.;  
  const AliESDVertex *primVtx = fESDEvent->GetPrimaryVertex();      

  
  const Int_t nvarESD = 9;
  Double_t var[nvarESD];
  var[0] = -1;

  Int_t sourcePdg = -1;

  if(TMath::Abs(particle->GetPdgCode())==kPDGelectron){
    
    sourcePdg = ElectronFromSource(stack, eleLabel);     
    
    if(sourcePdg==kPDGgamma){ // check direct photon or not 
      if(ElePhotonDirect(stack, eleLabel)!=-1)	
	var[0] = kEleDirectPhotonConv;    
      else 	
	var[0] = kElePhotonConv;     
      if(fDeDebugLevel>=10) 
	printf("photonconv=====this electron %d is from %d, source id=%.1f\n", eleLabel, sourcePdg, var[0]);   
    } else if(sourcePdg==kPDGpi0){ // check dalitz     
      var[0] = kElePi0;   
      if(fDeDebugLevel>=10) printf("pi0======this electron %d is from %d, source id=%.1f\n", eleLabel, sourcePdg, var[0]);    
    } else if(sourcePdg==kPDGeta){ // check eta
      var[0] = kEleEta;     
      if(fDeDebugLevel>=10) 
	printf("eta=====this electron %d is from %d, source id=%.1f\n", eleLabel, sourcePdg, var[0]);    
    } else if(IsB(sourcePdg)) var[0] = kEleB; // check beauty
      else if(IsC(sourcePdg)) var[0] = CheckCharm(stack,eleLabel); // check charm 
    
  } else if(TMath::Abs(particle->GetPdgCode())==kPDGpion) var[0] = kEleMissIDpion;
  else var[0] = kEleMissID;
  // ---- PID END ---- 
  

  // excluding current track
  // ---- beginning --- method from Andrea D 28.05.2010
  AliVertexerTracks *vertexer = new AliVertexerTracks(magneticField);
  vertexer->SetITSMode();
  vertexer->SetMinClusters(fNclustersITS);
  Int_t skipped[2];
  skipped[0] = (Int_t)track->GetID();
  vertexer->SetSkipTracks(1,skipped);
  AliESDVertex *vtxESDSkip = (AliESDVertex*)vertexer->FindPrimaryVertex(fESDEvent);
  delete vertexer; vertexer = NULL;
  if(vtxESDSkip->GetNContributors()<fMinNprimVtxContributor) return;
  // -- ending --- method from Andrea D 28.05.2010 
  
  Double_t dz[2];   // error of dca in cm
  Double_t covardz[3];
  if(!track->PropagateToDCA(vtxESDSkip,magneticField, beampiperadius, dz, covardz)) return; // protection    
  
  Double_t dcaXY = dz[0]*1.0e4;  // conv dca in cm to dca in micron 
  Double_t dcaZ = dz[1]*1.0e4;  
  
  if(fDeDebugLevel>=10) printf("others----->electron %d has mother %d is from %.1f\n", eleLabel, ElectronFromSource(stack, eleLabel), var[0]);
  
  if(TMath::Abs(dcaXY)<1000) var[1] = Int_t(dcaXY)/50;   // larger than 1mm should go to the last bin
  else
    var[1] = ((dcaXY>0)?1:-1)*20;
  
  if(TMath::Abs(dcaZ)<1000) var[2] = Int_t(dcaZ)/100;
  else
    var[2] = ((dcaZ>0)?1:-1)*10;
  
  // calculate dca using AliKFParticle class------------------------------------------------------------------
  Float_t  kfDcaXY = 0;
  Int_t trkID = track->GetID();  
  AliKFParticle::SetField(magneticField);
  AliKFParticle kfParticle(*track, particle->GetPdgCode());  
  // prepare kfprimary vertex
  AliKFVertex kfESDprimary;
  // Reconstruct Primary Vertex (with ESD tracks)
  Int_t n=primVtx->GetNIndices();
  if (n>0 && primVtx->GetStatus()){
    kfESDprimary = AliKFVertex(*primVtx);        
    UShort_t *priIndex = primVtx->GetIndices();    
    for (Int_t i=0;i<n;i++){
      Int_t idx = Int_t(priIndex[i]);
      if (idx == trkID){
	kfESDprimary -= kfParticle;
	kfDcaXY = kfParticle.GetDistanceFromVertexXY(kfESDprimary)*1e4;
      }  // remove current track from this calculation
    }  // loop over all primary vertex contributors    
  }
  // end of KF dca
  
  if(TMath::Abs(kfDcaXY)<1000) 
    var[3] = Int_t(kfDcaXY)/50;
  else
    var[3] = (kfDcaXY>0?1:-1)*20;
  
  var[4] = pt;  // pt 
  var[5] = eta; // eta
  var[6] = phi; // phi    
  var[7] = mcR;
  var[8] = mcStatus;
  
  (static_cast<THnSparseF *>(fDeOutputList->At(kEsdElectron)))->Fill(var);
  
}

//__________________________________________________________
void AliHFEdisplacedElectrons::FillDataOutput(const AliESDEvent * const fESDEvent, AliESDtrack* const esdTrack)
{
  
  // this is pure data, without MC information at all
  // fill output: with HFE pid selection of electrons after all track quality cuts 

  if(!esdTrack) return;
  AliESDtrack *track = esdTrack;
  
  Double_t pt   = track->Pt();
  Double_t eta = track->Eta();
  Double_t phi = track->Phi();
  
  // obtain impact parameters in xy and y
  const AliESDVertex *primVtx = fESDEvent->GetPrimaryVertex();    
  
  Float_t magneticField = 5;  // initialized as 5kG
  magneticField = fESDEvent->GetMagneticField();  // in kG 
  Double_t beampiperadius=3.;     
  
 
  const Int_t nvarData = 6;
  Double_t varData[nvarData]; 
    
  //
  // excluding current track
  //
  
  //------ beginning --- method from Andrea D 28.05.2010
  AliVertexerTracks *vertexer = new AliVertexerTracks(fESDEvent->GetMagneticField());
  vertexer->SetITSMode();
  vertexer->SetMinClusters(fNclustersITS);
  Int_t skipped[2];
  skipped[0] = (Int_t)track->GetID();
  vertexer->SetSkipTracks(1,skipped);
  AliESDVertex *vtxESDSkip = (AliESDVertex*)vertexer->FindPrimaryVertex(fESDEvent);
  delete vertexer; vertexer = NULL;
  if(vtxESDSkip->GetNContributors()<fMinNprimVtxContributor) return;
  //------ ending --- method from Andrea D 28.05.2010 
 
  Double_t dz[2];   // error of dca in cm
  Double_t covardz[3];

  if(!track->PropagateToDCA(vtxESDSkip,magneticField, beampiperadius, dz, covardz)) return; // protection 
  
  Double_t dcaXY = dz[0]*1.0e4;  // conv dca in cm to dca in micron 
  Double_t dcaZ = dz[1]*1.0e4;
   
  if(TMath::Abs(dcaXY)<1000) varData[0] = Int_t(dcaXY)/50;   // larger than 1mm should go to the last bin
  else
    varData[0] = (dcaXY>0?1:-1)*20;
  if(TMath::Abs(dcaZ)<1000) varData[1] = Int_t(dcaZ)/100;
  else
    varData[1] = (dcaZ>0?1:-1)*10;


  // calculate dca using AliKFParticle class------------------------------------------------------------------
  Float_t  kfDcaXY = 0;
  Int_t trkID = track->GetID();  
  AliKFParticle::SetField(magneticField);
  AliKFParticle kfParticle(*track, -11*track->Charge());  
  // prepare kfprimary vertex
  AliKFVertex kfESDprimary;
  // Reconstruct Primary Vertex (with ESD tracks)
  Int_t n=primVtx->GetNIndices();
  if (n>0 && primVtx->GetStatus()){
    kfESDprimary = AliKFVertex(*primVtx);        
    UShort_t *priIndex = primVtx->GetIndices();    
    for (Int_t i=0;i<n;i++){
      Int_t idx = Int_t(priIndex[i]);
      if (idx == trkID){
	  kfESDprimary -= kfParticle;
	  kfDcaXY = kfParticle.GetDistanceFromVertexXY(kfESDprimary)*1e4;
      }  // remove current track from this calculation
    }  // loop over all primary vertex contributors    
  }
  if(TMath::Abs(kfDcaXY)<1000)
    varData[2] = Int_t(kfDcaXY)/50;
  else
    varData[2] = ((kfDcaXY)>0?1:-1)*20;
  // end of KF dca

  varData[3] = pt; // pt 
  varData[4] = eta; //eta
  varData[5] = phi; // phi
  
  (static_cast<THnSparseF *>(fDeOutputList->At(kDataElectron)))->Fill(varData);

}

//__________________________________________________________
Int_t AliHFEdisplacedElectrons::ElectronFromSource(AliStack * const stack, Int_t label) const
{

  //
  //  return electron source label via electron label
  //

  // kEleSource is supposed to be either photon conv, direct photon conv, pi0, eta, beauty --> 0, 1, 2, 3, 4

  Int_t eleLabel = label;

  if(eleLabel<0 || eleLabel>stack->GetNtrack()) return -1;

  TParticle *particle = stack->Particle(eleLabel);
  if(!particle) return -1;

  Int_t motherLabel = particle->GetFirstMother();
  if(motherLabel<0 || motherLabel>stack->GetNtrack()) return -1;
  TParticle *motherPart = stack->Particle(motherLabel);
  if(!motherPart) return -1;

  Int_t pdgCode = TMath::Abs(motherPart->GetPdgCode());

  if(pdgCode==kPDGelectron) {
    if(fDeDebugLevel>=10) printf("particle label: %d...(motherLabel=%d : motherPdg=%d)  grandmother's pdg code was returned...%d \n",
			       label, motherLabel, pdgCode, ElectronFromSource(stack, motherLabel));
    return ElectronFromSource(stack, motherLabel);
  }

  return pdgCode;    

}

//__________________________________________________________
Int_t AliHFEdisplacedElectrons::ElectronFromCharm(AliStack * const stack, Int_t label) const
{
  //
  //  separate electron: kEleC from c->eX, kEleBC from b->c->eX
  //

  Int_t motherLabel = label;
  TParticle *motherParticle = stack->Particle(motherLabel);  // mother part
  if(!motherParticle) return -1;
  Int_t gMotherLabel = motherParticle->GetFirstMother();  // grand mother
  if(gMotherLabel<0 || gMotherLabel>stack->GetNtrack()) return -1;
  
  TParticle *gMotherPart = stack->Particle(gMotherLabel);
  if(!gMotherPart) return -1;
  
  Int_t pdgCode = gMotherPart->GetPdgCode();
  if(TMath::Abs(pdgCode%10000)/100==kPDGbeauty ||    // for special intermediate meson states: like 100553
     TMath::Abs(pdgCode)/1000==kPDGbeauty || TMath::Abs(pdgCode)/100==kPDGbeauty || TMath::Abs(pdgCode)==kPDGbeauty) {
    if(fDeDebugLevel>=10)  printf("this electron label %d is from mother %d, and finally from %d\n", label, ElectronFromSource(stack, label), pdgCode);
    return kEleBC;
  }  // for sure it is from BC
  
  else 
    if(TMath::Abs(pdgCode%10000)/100==kPDGcharm ||  // for special intermediate meson states: like 100443, 10443, 204*3 (*=1, 2, 3, 4)
       TMath::Abs(pdgCode)/1000==kPDGcharm || 
       TMath::Abs(pdgCode)/100==kPDGcharm || 
       TMath::Abs(pdgCode)==kPDGcharm){
      
      if(CharmFromBeauty(stack, gMotherLabel)!=-1) return kEleBC;
      else return kEleC;
      
    }
  
    else
      return -1;
     
}

//__________________________________________________________
Int_t AliHFEdisplacedElectrons::CharmFromBeauty(AliStack * const stack, Int_t hfLabel) const
{

  //
  //  check if charm meson/hadron is from beauty decay
  //  -1, not from beauty  
  //  return the label of the mother of this charm hadron
  //

  Int_t charmLabel = hfLabel;
  
  //  AliStack *stack = fMC->Stack();
  if(charmLabel<0 || charmLabel>stack->GetNtrack()) return -1;
  TParticle *particle = stack->Particle(charmLabel);
  if(!particle) return -1;

  Int_t motherCharmLabel = particle->GetFirstMother();
  if(motherCharmLabel<0 || motherCharmLabel>stack->GetNtrack()) return -1;

  TParticle *motherCharmPart = stack->Particle(motherCharmLabel);
  if(!motherCharmPart) return -1;

  Int_t pdgCode = motherCharmPart->GetPdgCode();
  
  if(TMath::Abs(pdgCode%10000)/100==kPDGbeauty ||   // for special intermediate meson states: like 100553
     TMath::Abs(pdgCode)/1000==kPDGbeauty || 
     TMath::Abs(pdgCode)/100==kPDGbeauty || 
     TMath::Abs(pdgCode)==kPDGbeauty) 
    return motherCharmLabel;  
  else 
    if(TMath::Abs(pdgCode%10000)/100==kPDGcharm ||    // for special intermediate meson states: like 100443, 10443, 204*3 (*=1, 2, 3, 4)
       TMath::Abs(pdgCode)/1000==kPDGcharm || 
       TMath::Abs(pdgCode)/100==kPDGcharm) 
      return CharmFromBeauty(stack, motherCharmLabel);   // loop over to see if charm is from beauty -- if yes, return the label of mother of the charm
  
    else
      return -1;
}

//__________________________________________________________
Int_t AliHFEdisplacedElectrons::ElePhotonDirect(AliStack * const stack, Int_t label) const 
{
  //
  // electron is from photon, and check if this photon is direct
  //

  Int_t eleLabel = label;
  
  TParticle *particle = stack->Particle(eleLabel);
  if(particle->GetFirstMother()<0 || particle->GetFirstMother()>stack->GetNtrack())
    return -1;
  
  TParticle *photonPart = stack->Particle(particle->GetFirstMother());
  if(!photonPart) 
    return -1;
  Int_t motherPhotonLabel = photonPart->GetFirstMother();
  if(motherPhotonLabel<0 || motherPhotonLabel>stack->GetNtrack())
    return -1;

  TParticle *motherPhotonPart = stack->Particle(motherPhotonLabel);
  if(!motherPhotonPart) 
    return -1;
  
  Int_t pdgMotherPhoton = motherPhotonPart->GetPdgCode();
  if(TMath::Abs(pdgMotherPhoton)<=10 || TMath::Abs(pdgMotherPhoton)==21)
    return 1;

  else 
    return -1;


}

//__________________________________________________________
Int_t AliHFEdisplacedElectrons::GetMCpid(AliStack* const stack, Int_t partLabel) const
{

  //
  // Simply pdg code
  //
  
  Int_t label = partLabel;
//   AliStack* stack = fMC->Stack();
  if((label < 0) || (label >= stack->GetNtrack())) return -1;  

  // MC Information
  TParticle * particle = stack->Particle(label);
  if(!particle) return -1;
  Int_t pdg = particle->GetPdgCode();

  return pdg;
 
}

//__________________________________________________________
Int_t AliHFEdisplacedElectrons::GetMotherLabel(AliStack *const stack, Int_t eleLabel) const 
{
  //
  // Simply label of mother
  //


  Int_t label = eleLabel;
  //  AliStack* stack = fMC->Stack();
  if((label < 0) || (label >= stack->GetNtrack())) return -1;  

  // MC Information
  TParticle * particle = stack->Particle(label);
  if(!particle) return -1;

  return particle->GetFirstMother();
 
}


//__________________________________________________________
Float_t AliHFEdisplacedElectrons::GetRapidity(TParticle *part) const
{
  // return rapidity
  
  Float_t rapidity = -999;        
  if((part->Energy() - part->Pz())*(part->Energy() + part->Pz())>0)
    rapidity = 0.5*(TMath::Log((part->Energy()+part->Pz()) / (part->Energy()-part->Pz()))); 
  
  return rapidity;
}


//__________________________________________________________
Float_t AliHFEdisplacedElectrons::GetTrackRapidity(AliESDtrack * const track) const
{
  // return rapidity of electron
  
  Float_t px = track->Px();
  Float_t py = track->Py();
  Float_t pz = track->Pz();

  // electron mass 0.00051099906 GeV/c2
  TParticlePDG* electron = TDatabasePDG::Instance()->GetParticle(kPDGelectron);
  Double_t mass = electron->Mass();
  
  Float_t en = TMath::Sqrt(px*px + py*py + pz*pz + mass*mass);

  Float_t rapidity = -999;        
  if((en - pz)*(en + pz)>0)
    rapidity = 0.5*(TMath::Log((en - pz)*(en + pz))); 
  
  return rapidity;
}

//__________________________________________________________
Int_t AliHFEdisplacedElectrons::CheckCharm(AliStack * const stack, Int_t eleLabel)
{
  // checks the genealogy of the ele from charm for a beauty
  // this method needs the stack and the label of the electron
  // warning: it assumes that the ele comes from a charm

  TParticle *particle = stack->Particle(eleLabel);
  Int_t label = particle->GetFirstMother();

  while(label>0 && label<(stack->GetNtrack())){
    particle = stack->Particle(label);
    if(IsB(TMath::Abs(particle->GetPdgCode()))) return kEleBC;
    label = particle->GetFirstMother();
  }

  return kEleC; 
}

//__________________________________________________________
Bool_t AliHFEdisplacedElectrons::IsB(Int_t pdg) const
{
  // check if the pdg is that of a beauty particle
 
  if((TMath::Abs(pdg)%10000)/100==kPDGbeauty ||  
     TMath::Abs(pdg)/1000==kPDGbeauty ||      
     TMath::Abs(pdg)/100==kPDGbeauty ||     
     TMath::Abs(pdg)==kPDGbeauty) return kTRUE;
  else return kFALSE;
} 

//__________________________________________________________
Bool_t AliHFEdisplacedElectrons::IsC(Int_t pdg) const
{
  // check if the pdg is that of a charmed particle
   
  if((TMath::Abs(pdg)%10000)/100==kPDGcharm ||   
     TMath::Abs(pdg)/1000==kPDGcharm ||       
     TMath::Abs(pdg)/100==kPDGcharm ||        
     TMath::Abs(pdg)==kPDGcharm) return kTRUE;
  else return kFALSE; 
}
