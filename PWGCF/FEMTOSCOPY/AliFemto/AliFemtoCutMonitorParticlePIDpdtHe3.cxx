////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// AliFemtoCutMonitorParticlePIDpdtHe3 - the cut monitor for particles to study     //
// various aspects of the PID determination                                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
#include "AliFemtoCutMonitorParticlePIDpdtHe3.h"
#include "AliFemtoModelHiddenInfo.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TList.h>
#include <TMath.h>

AliFemtoCutMonitorParticlePIDpdtHe3::AliFemtoCutMonitorParticlePIDpdtHe3():
  AliFemtoCutMonitor()
  , fTOFParticle(0)
  , fIfUsePt(false)
  , fSaveDCA(0)
  , fTPCdEdx(nullptr)
  , fTOFNSigma(nullptr)
  , fTPCNSigma(nullptr)
  , fMass(nullptr)
{
  // Default constructor
  fTPCdEdx        = new TH2D("TPCdEdx", "TPC dEdx vs. momentum", 100, 0.0, 5.0, 250, 0.0, 500.0);
  fTOFNSigma      = new TH2D("TOFNSigma","TOF NSigma vs. momentum", 100, 0.0, 5.0, 100, -5.0, 5.0);
  fTPCNSigma      = new TH2D("TPCNSigma","TPC NSigma vs. momentum", 100, 0.0, 5.0, 100, -5.0, 5.0);
  fMass           = new TH2D("Mass", "m vs. p", 100, 0.0, 5.0, 250, -1.0, 10.0);
   
      fDCARPt     = new TH2D("DCARPt", "DCA in XY vs. Pt", 400, -3.0, 3.0, 100,0.0,5.0);
      fDCAZPt     = new TH2D("DCAZPt", "DCA in Z vs. Pt", 400, -3.0, 3.0, 100,0.0,5.0);
  

}

AliFemtoCutMonitorParticlePIDpdtHe3::AliFemtoCutMonitorParticlePIDpdtHe3(const char *aName, Int_t aTOFParticle, int MassBin, float LowMass, float UpMass):
  AliFemtoCutMonitor()
  , fTOFParticle(aTOFParticle)
  , fIfUsePt(false)
  , fSaveDCA(0)
  , fTPCdEdx(nullptr)
  , fTOFNSigma(nullptr)
  , fTPCNSigma(nullptr)
  , fMass(nullptr)
  , fDCARPt(nullptr)
  , fDCAZPt(nullptr) 
{
  // Normal constructor
  fTPCdEdx        = new TH2D(TString::Format("TPCdEdx%s", aName), "TPC dEdx vs. momentum", 200, 0.0, 4.0, 250, 0.0, 500.0);
  fTOFNSigma      = new TH2D(TString::Format("TOFNSigma%s", aName), "TOF NSigma vs. momentum", 100, 0.0, 5.0, 100, -5.0, 5.0);
  fTPCNSigma      = new TH2D(TString::Format("TPCNSigma%s", aName), "TPC NSigma vs. momentum", 100, 0.0, 5.0, 100, -5.0, 5.0);
  fMass           = new TH2D(TString::Format("Mass%s", aName), "m2 vs. p", 100, 0.0, 5.0,MassBin,LowMass,UpMass);
}

AliFemtoCutMonitorParticlePIDpdtHe3::AliFemtoCutMonitorParticlePIDpdtHe3(const AliFemtoCutMonitorParticlePIDpdtHe3 &aCut):
  AliFemtoCutMonitor(aCut)
  , fTOFParticle(aCut.fTOFParticle)
  , fIfUsePt(aCut.fIfUsePt)
  , fSaveDCA(aCut.fSaveDCA)
  , fTPCdEdx(new TH2D(*aCut.fTPCdEdx))
  , fTOFNSigma(new TH2D(*aCut.fTOFNSigma))
  , fTPCNSigma(new TH2D(*aCut.fTPCNSigma))
  , fMass(new TH2D(*aCut.fMass))
  , fDCARPt(new TH2D(*aCut.fDCARPt))
  , fDCAZPt(new TH2D(*aCut.fDCAZPt))

{
  // copy constructor
}

AliFemtoCutMonitorParticlePIDpdtHe3::~AliFemtoCutMonitorParticlePIDpdtHe3()
{ 
  // Destructor
  delete fTPCdEdx;
  delete fTOFNSigma;
  delete fTPCNSigma;
  delete fMass;
  delete fDCARPt;
  delete fDCAZPt;
}

AliFemtoCutMonitorParticlePIDpdtHe3& AliFemtoCutMonitorParticlePIDpdtHe3::operator=(const AliFemtoCutMonitorParticlePIDpdtHe3& aCut)
{
  // assignment operator
  if (this == &aCut)
    return *this;

  AliFemtoCutMonitor::operator=(aCut);

  fTOFParticle  = aCut.fTOFParticle;
  fIfUsePt      = aCut.fIfUsePt;
  fSaveDCA 	= aCut.fSaveDCA;
  *fTPCdEdx     = *aCut.fTPCdEdx;
  *fTOFNSigma   = *aCut.fTOFNSigma;
  *fTPCNSigma   = *aCut.fTPCNSigma;
  *fMass        = *aCut.fMass;

    if(fSaveDCA){
      *fDCARPt = *aCut.fDCARPt;
      *fDCAZPt = *aCut.fDCAZPt; 
  } 
  return *this;
}

AliFemtoString AliFemtoCutMonitorParticlePIDpdtHe3::Report()
{
  // Prepare report from the execution
  TString report = "*** AliFemtoCutMonitorParticlePIDpdtHe3 report";
  return AliFemtoString(report.Data());
}

void AliFemtoCutMonitorParticlePIDpdtHe3::Fill(const AliFemtoTrack* aTrack)
{
    // Fill in the monitor histograms with the values from the current track
    float tMom = (fIfUsePt) ? aTrack->Pt() : aTrack->P().Mag();
    float tdEdx = aTrack->TPCsignal();
    fTPCdEdx->Fill(tMom, tdEdx);

    if (fTOFParticle == 0) fTOFNSigma->Fill(tMom, aTrack->NSigmaTOFPi());
    if (fTOFParticle == 1) fTOFNSigma->Fill(tMom, aTrack->NSigmaTOFK());
    if (fTOFParticle == 2) fTOFNSigma->Fill(tMom, aTrack->NSigmaTOFP());
    if (fTOFParticle == 3) fTOFNSigma->Fill(tMom, aTrack->NSigmaTOFD());
    if (fTOFParticle == 4) fTOFNSigma->Fill(tMom, aTrack->NSigmaTOFT());
    if (fTOFParticle == 5) fTOFNSigma->Fill(tMom, aTrack->NSigmaTOFH());
    if (fTOFParticle == 6) fTOFNSigma->Fill(tMom, aTrack->NSigmaTOFE());

    if (fTOFParticle == 0) fTPCNSigma->Fill(tMom, aTrack->NSigmaTPCPi());
    if (fTOFParticle == 1) fTPCNSigma->Fill(tMom, aTrack->NSigmaTPCK());
    if (fTOFParticle == 2) fTPCNSigma->Fill(tMom, aTrack->NSigmaTPCP());
    if (fTOFParticle == 3) fTPCNSigma->Fill(tMom, aTrack->NSigmaTPCD());
    if (fTOFParticle == 4) fTPCNSigma->Fill(tMom, aTrack->NSigmaTPCT());
    if (fTOFParticle == 5) fTPCNSigma->Fill(tMom, aTrack->NSigmaTPCH());
    if (fTOFParticle == 6) fTPCNSigma->Fill(tMom, aTrack->NSigmaTPCE());

    double massTOF;//Mass square
    
    double massPDGPi    = 0.13957018;
    double massPDGK     = 0.493677;
    double massPDGP     = 0.938272013;
    double massPDGD     = 1.8756;
    double massPDGT     = 2.8089;
    double massPDGHe3   = 2.8084;

    float tMomTotal 	= aTrack->P().Mag();
    float c=1;
    double beta =aTrack->VTOF();
    if(beta!=0){
        massTOF= tMomTotal * tMomTotal/c/c*(1./(beta*beta)-1.);
        fMass->Fill(tMom,massTOF);
    }
    
    if(fSaveDCA){
      float tPt = aTrack->Pt();
      float dcar = aTrack->ImpactD();
      float dcaz = aTrack->ImpactZ();
      fDCARPt->Fill(dcar, tPt);
      fDCAZPt->Fill(dcaz, tPt);
    }    

}

void AliFemtoCutMonitorParticlePIDpdtHe3::Write()
{
    // Write out the relevant histograms
    fTPCdEdx->Write();
    fTOFNSigma->Write();
    fTPCNSigma->Write();
    fMass->Write();
    if(fSaveDCA){
        fDCARPt->Write();
        fDCAZPt->Write();
    }
}

TList *AliFemtoCutMonitorParticlePIDpdtHe3::GetOutputList()
{
    TList *tOutputList = new TList();
    tOutputList->Add(fTPCdEdx);
    tOutputList->Add(fTOFNSigma);
    tOutputList->Add(fTPCNSigma);
    tOutputList->Add(fMass);
   
    if(fSaveDCA){
        tOutputList->Add(fDCARPt);
        tOutputList->Add(fDCAZPt);
    }
       
        
    return tOutputList;
}
void AliFemtoCutMonitorParticlePIDpdtHe3::SetDCAInit(TString aName,int nbinsDCA,float lowDCA,float upDCA,
int nbinspT,float lowpT,float uppT){

	fDCARPt = new TH2D(aName+"DCAR","",nbinsDCA,lowDCA,upDCA,nbinspT,lowpT,uppT);
	fDCAZPt = new TH2D(aName+"DCAZ","",nbinsDCA,lowDCA,upDCA,nbinspT,lowpT,uppT);
}
