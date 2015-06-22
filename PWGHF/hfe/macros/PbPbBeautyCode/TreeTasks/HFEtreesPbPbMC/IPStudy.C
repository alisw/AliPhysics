#define IPStudy_cxx
// The class definition in DummySelectorReduced.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("DummySelectorReduced.C")
// Root > T->Process("DummySelectorReduced.C","some options")
// Root > T->Process("DummySelectorReduced.C+")
//

#include "IPStudy.h"
#include <TH2.h>
#include <TH1D.h>
#include <TStyle.h>
//#include <AliHFEreducedTrack.h>

void IPStudy::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

void IPStudy::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
   
   Double_t ptbinningX[19] = {0., 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 2., 2.5, 3., 4., 5., 6., 8., 12., 16., 20.};
   Double_t IPbinningY[4001];
   for(int i=0; i<=4000; i++){IPbinningY[i]=-0.2 + double(i)*0.0001;}
   Double_t IPbinningYSmall[401];
   for(int i=0; i<=400; i++){IPbinningYSmall[i]=-0.2 + double(i)*0.001;}
   Double_t IPerrbinningY[401];
   for(int i=0; i<=400; i++){IPerrbinningY[i]=double(i)*0.0001;}
   Double_t mctruthbinning[19] = {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5};
   Double_t largepTBinning[51];
   for(int i=0; i<=50; i++){largepTBinning[i]=double(i)*1.0;}
   
   fTrackDCApTnoCut = new TH3D("TrackDCApTnoCut","", 18, ptbinningX, 4000, IPbinningY, 6, mctruthbinning);
   fTrackDCApTnoCutnoTOF = new TH3D("TrackDCApTnoCutnoTOF","", 18, ptbinningX, 4000, IPbinningY, 6, mctruthbinning);
   fTrackDCApTCut = new TH3D("TrackDCApTCut","", 18, ptbinningX, 400, IPbinningYSmall, 6, mctruthbinning);
   fTrackDCApTnoCutV0 = new TH2D("TrackDCApTnoCutV0","", 18, ptbinningX, 400, IPbinningYSmall);
   fTrackDCApTnoCutV0weakerConditions = new TH2D("TrackDCApTnoCutV0weakerConditions","", 18, ptbinningX, 400, IPbinningYSmall);
   fTrackDCApTnoCutV0weakerConditionsNonFake = new TH2D("TrackDCApTnoCutV0weakerConditionsNonFake","", 18, ptbinningX, 400, IPbinningYSmall);
   fTrackDCApTnoCutV0weakerConditionsFake = new TH2D("TrackDCApTnoCutV0weakerConditionsFake","", 18, ptbinningX, 400, IPbinningYSmall);
   fTrackDCApTnoCutV0weakerConditionsDalitz = new TH2D("TrackDCApTnoCutV0weakerConditionsDalitz","", 18, ptbinningX, 400, IPbinningYSmall);
   fTrackDCApTnoCutStrange = new TH2D("TrackDCApTnoCutStrange","", 18, ptbinningX, 400, IPbinningYSmall);
   fTrackDCApTnoCutDalitzWOStrange = new TH2D("TrackDCApTnoCutDalitzWOStrange","", 18, ptbinningX, 400, IPbinningYSmall);
   fTrackDCApTnoCutConversionWOStrange = new TH2D("TrackDCApTnoCutConversionWOStrange","", 18, ptbinningX, 400, IPbinningYSmall);
   fTrackDCApTCutV0 = new TH2D("TrackDCApTCutV0","", 18, ptbinningX, 400, IPbinningYSmall);
   fTrackDCAerrpT = new TH3D("TrackDCAerrpT","", 18, ptbinningX, 400, IPerrbinningY, 6, mctruthbinning);
   
   fBeautyMotherCorrelation = new TH3D("BeautyMotherCorrelation","", 18, ptbinningX, 400, IPbinningYSmall, 50, largepTBinning);
   fBeautyMotherCorrelationRAA = new TH3D("BeautyMotherCorrelationRAA","", 18, ptbinningX, 400, IPbinningYSmall, 50, largepTBinning);
   fBeautyMotherCorrelationHalfRAA = new TH3D("BeautyMotherCorrelationHalfRAA","", 18, ptbinningX, 400, IPbinningYSmall, 50, largepTBinning);
   fCharmMotherCorrelation = new TH3D("CharmMotherCorrelation","", 18, ptbinningX, 400, IPbinningYSmall, 50, largepTBinning);
   fCharmMotherCorrelationRAA = new TH3D("CharmMotherCorrelationRAA","", 18, ptbinningX, 400, IPbinningYSmall, 50, largepTBinning);
   fCharmMotherCorrelationHalfRAA = new TH3D("CharmMotherCorrelationHalfRAA","", 18, ptbinningX, 400, IPbinningYSmall, 50, largepTBinning);
   
   fTrackDCAerr = new TH1D("TrackDCAerr","", 100,0.0,0.03);
   fhMultiplicity= new TH1D("hMultiplicity","", 1000,0.0,10000);
   fhMultiplicityPle= new TH1D("hMultiplicityPle","", 1000,0.0,10000);
   
   fTrackDCAerrpTConv = new TH3D("TrackDCAerrpTConv","", 20, 0,5, 100,0.0,0.03, 2, -0.5, 1.5);
   fTrackProdRnoCut = new TH2D("TrackProdRnoCut","", 20, 0,5, 500, 0, 20.0);
   fTrackProdRnoCutDalitz = new TH2D("TrackProdRnoCutDalitz","", 20, 0,5, 2000, 0, 1.0);
   fTrackProdRnoCutV0 = new TH2D("TrackProdRnoCutV0","", 20, 0,5, 500, 0, 20.0);
   fTrackProdRwithCut = new TH2D("TrackProdRwithCut","", 20, 0,5, 500, 0, 20.0);
   fTrackDCAerrpTFake = new TH2D("TrackDCAerrpTFake","", 20, 0,5, 100,0.0,0.03);
   fTrackDCAerrpTTrueConv = new TH2D("TrackDCAerrpTTrueConv","", 20, 0,5, 100,0.0,0.03);
   fTrackEta = new TH2D("TrackEta","", 100, -1.0, 1.0, 10, -0.5, 9.5);
   fConvPhi = new TH1D("ConvPhi","", 200, -2.*3.14159, 2.*3.14159);
   fTOFPhi = new TH1D("TOFPhi","", 200, -2.*3.14159, 2.*3.14159);
   fEventNr = new TH1D("EventNr", "", 5, -0.5, 4.5);
   fGaussTest = new TH1D("GaussTest", "", 100, -5.0, 5.0);
   
   fWidthOfGauss = new TF1("widthOfGauss", "0.00707/x^0.613 * 0.536", 0.01, 100.);
   
   fTPCTOFsigmas = new TH3D("TPCTOFsigmas", "", 20, 0.0, 10., 800, -15.,5., 200, -5., 5.);
   fTPCClusters   = new TH2D("TPCClusters", "", 40, 0.0, 10., 300, -0.5, 299.5);
   fITSClusters   = new TH2D("ITSClusters", "", 40, 0.0, 10., 7, -0.5, 6.5);
   
   fVertexXY = new TH2D("VertexXY","", 100, -0.5, 0.5, 100, -0.5, 0.5);
   
   fV0sConvpT = new TH1D("V0sConvpT", "", 20, 0., 5.);
   fV0sDalitzpT = new TH1D("V0sDalitzpT", "", 20, 0., 5.);
   fV0sHFpT = new TH1D("V0sHFpT", "", 20, 0., 5.);
   
   fConvEta  = new TH2D("ConvEta", "", 40, 0.0, 10., 100, -0.8, 0.8);
   fOtherEta = new TH2D("OtherEta", "", 40, 0.0, 10., 100, -0.8, 0.8);
   fK0pTDCA = new TH2D("K0pTDCA", "", 40, 0.0, 10., 200, -0.2, 0.2);
   fK0pTDCARCut = new TH2D("K0pTDCARCut", "", 40, 0.0, 10., 200, -0.2, 0.2);
   fRKaons = new TH1D("RKaons", "", 1000, 0., 100.);
   
   fClusterPIDCluster = new TH2D("ClusterPIDCluster", "", 160, -0.5, 159.5, 160, -0.5, 159.5);
   
   fTrackDCApTnoCut12 = new TH3D("TrackDCApTnoCut12","", 18, ptbinningX, 400, IPbinningYSmall, 6, mctruthbinning);
   fBeautyMotherCorrelation12 = new TH2D("BeautyMotherCorrelation12","", 18, ptbinningX, 400, IPbinningYSmall);
   fBeautyMotherCorrelationRAA12 = new TH2D("BeautyMotherCorrelationRAA12","", 18, ptbinningX, 400, IPbinningYSmall);
   fBeautyMotherCorrelationHalfRAA12 = new TH2D("BeautyMotherCorrelationHalfRAA12","", 18, ptbinningX, 400, IPbinningYSmall);
   fCharmMotherCorrelation12 = new TH2D("CharmMotherCorrelation12","", 18, ptbinningX, 400, IPbinningYSmall);
   fCharmMotherCorrelationRAA12 = new TH2D("CharmMotherCorrelationRAA12","", 18, ptbinningX, 400, IPbinningYSmall);
   fCharmMotherCorrelationHalfRAA12 = new TH2D("CharmMotherCorrelationHalfRAA12","", 18, ptbinningX, 400, IPbinningYSmall);
   
   fRD = new TRandom3(5);  // old : 4
   
   fOutput->Add(fTrackDCAerr);
   fOutput->Add(fhMultiplicity);
   fOutput->Add(fhMultiplicityPle);
   fOutput->Add(fTrackDCAerrpT);
   fOutput->Add(fTrackDCAerrpTConv);
   fOutput->Add(fTrackDCApTnoCut);
   fOutput->Add(fTrackDCApTnoCutnoTOF);
   fOutput->Add(fTrackDCApTCut);
   fOutput->Add(fTrackDCApTnoCutV0);
   fOutput->Add(fTrackDCApTnoCutV0weakerConditions);
   fOutput->Add(fTrackDCApTnoCutV0weakerConditionsNonFake);
   fOutput->Add(fTrackDCApTnoCutV0weakerConditionsFake);
   fOutput->Add(fTrackDCApTnoCutV0weakerConditionsDalitz);
   fOutput->Add(fTrackDCApTnoCutStrange);
   fOutput->Add(fTrackDCApTnoCutDalitzWOStrange);
   fOutput->Add(fTrackDCApTnoCutConversionWOStrange);
   fOutput->Add(fTrackDCApTCutV0);
   fOutput->Add(fTrackProdRnoCut);
   fOutput->Add(fTrackProdRnoCutDalitz);
   fOutput->Add(fTrackProdRnoCutV0);
   fOutput->Add(fTrackProdRwithCut);
   fOutput->Add(fTrackDCAerrpTFake);
   fOutput->Add(fTrackDCAerrpTTrueConv);
   fOutput->Add(fTrackEta);
   fOutput->Add(fConvPhi);
   fOutput->Add(fTOFPhi);
   fOutput->Add(fTPCTOFsigmas);
   fOutput->Add(fEventNr);
   fOutput->Add(fGaussTest);
   fOutput->Add(fTPCClusters);
   fOutput->Add(fITSClusters);
   fOutput->Add(fV0sConvpT);
   fOutput->Add(fV0sDalitzpT);
   fOutput->Add(fV0sHFpT);
   fOutput->Add(fVertexXY);
   fOutput->Add(fConvEta);
   fOutput->Add(fOtherEta);
   fOutput->Add(fK0pTDCA);
   fOutput->Add(fK0pTDCARCut);
   fOutput->Add(fRKaons);
   fOutput->Add(fClusterPIDCluster);
   fOutput->Add(fBeautyMotherCorrelation);
   fOutput->Add(fBeautyMotherCorrelationRAA);
   fOutput->Add(fBeautyMotherCorrelationHalfRAA);
   fOutput->Add(fCharmMotherCorrelation);
   fOutput->Add(fCharmMotherCorrelationRAA);
   fOutput->Add(fCharmMotherCorrelationHalfRAA);
   fOutput->Add(fTrackDCApTnoCut12);
   fOutput->Add(fBeautyMotherCorrelation12);
   fOutput->Add(fBeautyMotherCorrelationRAA12);
   fOutput->Add(fBeautyMotherCorrelationHalfRAA12);
   fOutput->Add(fCharmMotherCorrelation12);
   fOutput->Add(fCharmMotherCorrelationRAA12);
   fOutput->Add(fCharmMotherCorrelationHalfRAA12);
   
}

Bool_t IPStudy::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either DummySelectorReduced::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the data. When processing
   // keyed objects with PROOF, the object is already loaded and is available
   // via the fObject pointer.
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.
  
  GetEntry(entry);
  
  
  
  if(EventPassesCuts())
  {
    int nTracks = HFEevent->GetNumberOfTracks();
    fEventNr->Fill(0);
    AliHFEreducedTrack * hfeTrack; 
  
    fhMultiplicity->Fill(nTracks);
    fVertexXY->Fill(HFEevent->GetVX(), HFEevent->GetVY());
  
    double pT;
    double resolution;
    double RAARandom;
  
    bool CorrectResolution = true;
    
    for(int j=0;j<nTracks;j++)
    {
      fGaussTest->Fill(sampleGaus());
      
      hfeTrack=(AliHFEreducedTrack*)HFEevent->GetTrack(j);
      resolution = hfeTrack->HFEImpactParameter() / hfeTrack->HFEImpactParameterResolution();
      
      
      
      if(hfeTrack->IsV0pion())  // GetV0prodR() might be an additional condition
	fK0pTDCA->Fill(hfeTrack->Pt(), hfeTrack->HFEImpactParameter()*hfeTrack->Charge());
      if(hfeTrack->IsV0pion())  // GetV0prodR() might be an additional condition
	if(hfeTrack->GetV0prodR()<5.0)
	  fK0pTDCARCut->Fill(hfeTrack->Pt(), hfeTrack->HFEImpactParameter()*hfeTrack->Charge());
      if(hfeTrack->IsV0pion())
	  fRKaons->Fill(hfeTrack->GetV0prodR());
      //source->Fill(hfeTrack->MCSource());
      //TrackQuantity->Fill(hfeTrack->GetITSnclusters());
      if(TMath::Abs(hfeTrack->Eta())<0.8 && hfeTrack->HasITSrefit() && hfeTrack->HasTPCrefit() && hfeTrack->MCProdRadius()<3.1
	&& hfeTrack->HasITScluster(0) && hfeTrack->HasITScluster(1) && hfeTrack->GetITSnclusters()>3
        && hfeTrack->Pt() > 1.5 && hfeTrack->Pt() < 2.0)
      {
	fClusterPIDCluster->Fill(hfeTrack->GetTPCnclusters(), hfeTrack->GetTPCnclusterPID());
      }
	
	
      if(hfeTrack->MCSource()<2 && TMath::Abs(hfeTrack->Eta())<0.8)
	fTPCClusters->Fill(hfeTrack->Pt(), hfeTrack->GetTPCnclusters());
      if(hfeTrack->MCSource()<2 && TMath::Abs(hfeTrack->Eta())<0.8)
	fITSClusters->Fill(hfeTrack->Pt(), hfeTrack->GetITSnclusters());
      
      if(hfeTrack->IsV0electron() &&
	TMath::Abs(hfeTrack->Eta())<0.8 && hfeTrack->HasITSrefit() && hfeTrack->HasTPCrefit()
	&& hfeTrack->HasITScluster(0) && hfeTrack->HasITScluster(1) && hfeTrack->GetITSnclusters()>3   
        && hfeTrack->MCSource()<=3)  
      {
	  fTrackDCApTnoCutV0weakerConditions->Fill(hfeTrack->Pt(),hfeTrack->HFEImpactParameter()*hfeTrack->Charge()+0.458/*0.50*/*resolution*sampleGaus());
	  if(hfeTrack->MCProdRadius() < 4. && hfeTrack->MCProdRadius() > 2.)
	    fTrackDCApTnoCutV0weakerConditionsNonFake->Fill(hfeTrack->Pt(),hfeTrack->HFEImpactParameter()*hfeTrack->Charge()+0.458/*0.50*/*resolution*sampleGaus());
	  if(hfeTrack->MCProdRadius() > 4.)
	    fTrackDCApTnoCutV0weakerConditionsFake->Fill(hfeTrack->Pt(),hfeTrack->HFEImpactParameter()*hfeTrack->Charge()+0.458/*0.50*/*resolution*sampleGaus());
	  if(hfeTrack->MCProdRadius() < 2.)
	    fTrackDCApTnoCutV0weakerConditionsDalitz->Fill(hfeTrack->Pt(),hfeTrack->HFEImpactParameter()*hfeTrack->Charge()+0.458/*0.50*/*resolution*sampleGaus());
      }
      
      
      // Main part
      
      //if(TMath::Abs(hfeTrack->MCElectronSource()) >= 21 && TMath::Abs(hfeTrack->MCElectronSource()) <= 40 && passesCuts(hfeTrack) && hfeTrack->MCSource()<=6)
	//cout << "MCElectronSource: " << hfeTrack->MCElectronSource() << " MCSource: " << hfeTrack->MCSource() << endl;
      
      if(passesCuts(hfeTrack) && hfeTrack->MCSource()<=3){  // 0,1,2,3=charm,beauty, conversion, 
	
      double addRes = sampleGaus();
      double IP=hfeTrack->HFEImpactParameter()*hfeTrack->Charge()+0.458*resolution*addRes; // 10% resolution correction
      double IP12=hfeTrack->HFEImpactParameter()*hfeTrack->Charge()+0.50*resolution*addRes; // 12% resolution correction
      pT=hfeTrack->Pt();
      
      if(CorrectResolution)
        fTrackDCApTnoCutnoTOF->Fill(hfeTrack->Pt(),IP , hfeTrack->MCSource());
      else
	fTrackDCApTnoCutnoTOF->Fill(hfeTrack->Pt(),hfeTrack->HFEImpactParameter()*hfeTrack->Charge(), hfeTrack->MCSource());
      
      if(hfeTrack->HasTOFpid())
      {
      fTrackDCAerrpT->Fill(hfeTrack->Pt(), resolution, hfeTrack->MCSource());

      if(hfeTrack->GetTPCnclusters()>80 && hfeTrack->Pt()>0.5) fTrackEta->Fill(hfeTrack->MCEta(), hfeTrack->MCSource());
      
      if(hfeTrack->MCSource()==2) fConvEta->Fill(hfeTrack->Pt(), hfeTrack->Eta());
      else fOtherEta->Fill(hfeTrack->Pt(), hfeTrack->Eta());

      RAARandom = fRD->Rndm();
      
      if(hfeTrack->MCSource()==1){
        fBeautyMotherCorrelation->Fill(pT, IP, TMath::Abs(hfeTrack->MCElectronSourcePt()));
	fBeautyMotherCorrelation12->Fill(pT, IP12);}
      if(hfeTrack->MCSource()==1)
	if(RAARandom<(0.5/(1. + TMath::Exp((TMath::Abs(hfeTrack->MCElectronSourcePt())-7.)*0.7)) + 0.5 + (TMath::Abs(hfeTrack->MCElectronSourcePt())-15.)/300.)){
	  fBeautyMotherCorrelationRAA->Fill(pT, IP, TMath::Abs(hfeTrack->MCElectronSourcePt()));
	  fBeautyMotherCorrelationRAA12->Fill(pT, IP12);}
      if(hfeTrack->MCSource()==1)
	if(RAARandom<(1.+(0.5/(1. + TMath::Exp((TMath::Abs(hfeTrack->MCElectronSourcePt())-7.)*0.7)) + 0.5 + (TMath::Abs(hfeTrack->MCElectronSourcePt())-15.)/300.))/2.){
	  fBeautyMotherCorrelationHalfRAA->Fill(pT, IP, TMath::Abs(hfeTrack->MCElectronSourcePt()));
	  fBeautyMotherCorrelationHalfRAA12->Fill(pT, IP12);}
	
      if(hfeTrack->MCSource()==0){
        fCharmMotherCorrelation->Fill(pT, IP, TMath::Abs(hfeTrack->MCElectronSourcePt()));
	fCharmMotherCorrelation12->Fill(pT, IP12);}
      if(hfeTrack->MCSource()==0)
	if(RAARandom<(0.5*TMath::Gaus(hfeTrack->MCElectronSourcePt(), 2.5, 1.5) + 0.5 + (hfeTrack->MCElectronSourcePt()-15.)/300.)){
	  fCharmMotherCorrelationRAA->Fill(pT, IP, TMath::Abs(hfeTrack->MCElectronSourcePt()));
	  fCharmMotherCorrelationRAA12->Fill(pT, IP12);}
      if(hfeTrack->MCSource()==0)
	if(RAARandom<(1.+(0.5*TMath::Gaus(hfeTrack->MCElectronSourcePt(), 2.5, 1.5) + 0.5 + (hfeTrack->MCElectronSourcePt()-15.)/300.))/2.){
	  fCharmMotherCorrelationHalfRAA->Fill(pT, IP, TMath::Abs(hfeTrack->MCElectronSourcePt()));
	  fCharmMotherCorrelationHalfRAA12->Fill(pT, IP12);}
	
      if(TMath::Abs(hfeTrack->MCElectronSource()) >= 21 && TMath::Abs(hfeTrack->MCElectronSource()) <= 40)
	fTrackDCApTnoCutStrange->Fill(pT, IP);
      
	
      if(hfeTrack->MCSource()==2 && !(TMath::Abs(hfeTrack->MCElectronSource()) >= 21 && TMath::Abs(hfeTrack->MCElectronSource()) <= 40))
	fTrackDCApTnoCutConversionWOStrange->Fill(pT, IP);
      if(hfeTrack->MCSource()==3 && !(TMath::Abs(hfeTrack->MCElectronSource()) >= 21 && TMath::Abs(hfeTrack->MCElectronSource()) <= 40))
	fTrackDCApTnoCutDalitzWOStrange->Fill(pT, IP);
      
      
      
      if(CorrectResolution)
	//for(int m=0;m<1000;m++)
        fTrackDCApTnoCut->Fill(pT, IP, hfeTrack->MCSource());
      else
	fTrackDCApTnoCut->Fill(pT,hfeTrack->HFEImpactParameter()*hfeTrack->Charge(), hfeTrack->MCSource());
      fTrackDCApTnoCut12->Fill(pT, IP12, hfeTrack->MCSource());
      
      if(hfeTrack->IsV0electron())  // Already required that it is an electron (source<3)
	//for(int m=0;m<1000;m++)
	  fTrackDCApTnoCutV0->Fill(pT,IP);
	  
      if(hfeTrack->IsV0electron())
      {
	if(hfeTrack->MCSource()==2)fV0sConvpT->Fill(hfeTrack->Pt());
	if(hfeTrack->MCSource()==3)fV0sDalitzpT->Fill(hfeTrack->Pt());
	if(hfeTrack->MCSource()==1 || hfeTrack->MCSource()==0)fV0sHFpT->Fill(hfeTrack->Pt());
      }
	  
      if(hfeTrack->MCSource()==3)fTrackProdRnoCutDalitz->Fill(hfeTrack->MCPt(), hfeTrack->MCProdRadius());
      
	  if(hfeTrack->MCSource()==2)
	  {
	    if(hfeTrack->IsV0electron()) fTrackDCAerrpTConv->Fill(hfeTrack->Pt(), resolution, 0);
	    
	    fhMultiplicityPle->Fill(nTracks);
	    
	    fTrackDCAerrpTConv->Fill(hfeTrack->Pt(), resolution, 1);  // all
	    if(hfeTrack->MCProdRadius()>6.0)fConvPhi->Fill(hfeTrack->Phi());
	    if(hfeTrack->HasTOFpid())fTOFPhi->Fill(hfeTrack->Phi());
	    
	    fTrackProdRnoCut->Fill(hfeTrack->Pt(), hfeTrack->MCProdRadius());
	    if(hfeTrack->IsV0electron()) fTrackProdRnoCutV0->Fill(hfeTrack->Pt(), hfeTrack->MCProdRadius());
	    if(hfeTrack->MCProdRadius()>6.0) fTrackDCAerrpTFake->Fill(hfeTrack->Pt(), resolution);
	    else fTrackDCAerrpTTrueConv->Fill(hfeTrack->Pt(), resolution);
	  }
	  pT=hfeTrack->Pt();
          if(resolution<=(0.007313/pT+0.00173))
	  {
            fTrackDCApTCut->Fill(hfeTrack->Pt(),hfeTrack->HFEImpactParameter()*hfeTrack->Charge(), hfeTrack->MCSource());
	    if(hfeTrack->IsV0electron())
	       fTrackDCApTCutV0->Fill(hfeTrack->Pt(),hfeTrack->HFEImpactParameter()*hfeTrack->Charge());
	    if(hfeTrack->MCSource()==2)fTrackProdRwithCut->Fill(hfeTrack->MCPt(), hfeTrack->MCProdRadius());
  
	  }
	  
      } 
      }
      
    }
  }


   return kTRUE;
}

double IPStudy::sampleGaus(void)
{
  double r = fRD->Rndm();
  return TMath::ErfInverse(r*2.-1.)*TMath::Sqrt(2.);
}

Bool_t IPStudy::EventPassesCuts()
{
   //cout << "N contrib: " << HFEevent->GetNContribVertexSPD() << endl;
   //if(HFEevent->GetNContribVertexSPD()<20) return kFALSE;
   if(!HFEevent->GetNContribVertex()) return kFALSE;
   if(!HFEevent->GetNContribVertexSPD()) return false;
   if(TMath::Abs(HFEevent->GetVZ() - HFEevent->GetVZSPD()) > 0.5) return false;
   if(HFEevent->GetVertexZResolutionSPD() > 0.25) return false;
   if(TMath::Abs(HFEevent->GetVZ()) > 10.) return kFALSE;
   if(!(HFEevent->GetCentrality()>0.0 && HFEevent->GetCentrality()<60.)) return kFALSE;  // temporary
   return true;
}

Bool_t IPStudy::passesCuts(AliHFEreducedTrack * track)
{
  //     if(track.Pt()<0.5) return kFALSE;;
  //    if(abs(track.Eta())>0.6) return kFALSE;
      if(!track->HasITSrefit()) return kFALSE;
      if(!track->HasTPCrefit()) return kFALSE;
      if(track->IsKinkDaughter()) return kFALSE;
      if(track->IsKinkMother()) return kFALSE;   
      if(track->MCProdRadius()>100.0) return kFALSE;
      //if(track->MCProdRadius()>3.5) return kFALSE;
      //if((track->MCProdRadius()>8.0 || track->MCProdRadius()<6.0) && track->MCSource()==2) return kFALSE;
      if(track->GetITSnclusters()<4) return kFALSE;
      if(track->GetTPCnclusters()<110) return kFALSE;
      if(!(track->GetChi2PerTPCcluster()<4.)) return kFALSE;
      if(track->GetTPCnclusterPID()<80) return kFALSE;
      if(track->GetTPCclusterRatio()<=0.6) return kFALSE;
      if(!(abs(track->DCAr())<1.)) return kFALSE;
      if(!(abs(track->DCAz())<2.)) return kFALSE;
      if(!(track->HasITScluster(0))) return kFALSE;
      if(!(track->HasITScluster(1))) return kFALSE;
      if(!(TMath::Abs(track->Eta())<0.8)) return kFALSE;
      //if(!track->HasTOFpid()) return kFALSE; 
      //if(abs(track->GetTPCsigmaElCorrected())>3.) return kFALSE; 
      return kTRUE;
}



void IPStudy::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void IPStudy::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
  
  TString histname = Form("outputTreeMCenhNeg060.root");
  
   fTrackDCAerr = (TH1D *)fOutput->FindObject("TrackDCAerr");
   fhMultiplicity = (TH1D*)fOutput->FindObject("hMultiplicity");
   fhMultiplicityPle = (TH1D*)fOutput->FindObject("hMultiplicityPle");
   fTrackDCAerrpT = (TH3D*)fOutput->FindObject("TrackDCAerrpT"); 
   fTrackDCAerrpTConv = (TH3D*)fOutput->FindObject("TrackDCAerrpTConv"); 
   fTrackDCApTnoCut = (TH3D*)fOutput->FindObject("TrackDCApTnoCut"); 
   fTrackDCApTnoCutnoTOF = (TH3D*)fOutput->FindObject("TrackDCApTnoCutnoTOF");
   fTrackDCApTCut = (TH3D*)fOutput->FindObject("TrackDCApTCut"); 
   fTrackDCApTnoCutV0 = (TH2D*)fOutput->FindObject("TrackDCApTnoCutV0");
   fTrackDCApTnoCutV0weakerConditions = (TH2D*)fOutput->FindObject("TrackDCApTnoCutV0weakerConditions");
   fTrackDCApTnoCutV0weakerConditionsNonFake = (TH2D*)fOutput->FindObject("TrackDCApTnoCutV0weakerConditionsNonFake");
   fTrackDCApTnoCutV0weakerConditionsFake = (TH2D*)fOutput->FindObject("TrackDCApTnoCutV0weakerConditionsFake");
   fTrackDCApTnoCutV0weakerConditionsDalitz = (TH2D*)fOutput->FindObject("TrackDCApTnoCutV0weakerConditionsDalitz");
   fTrackDCApTnoCutStrange = (TH2D*)fOutput->FindObject("TrackDCApTnoCutStrange");
   fTrackDCApTnoCutDalitzWOStrange = (TH2D*)fOutput->FindObject("TrackDCApTnoCutDalitzWOStrange");
   fTrackDCApTnoCutConversionWOStrange = (TH2D*)fOutput->FindObject("TrackDCApTnoCutConversionWOStrange");
   fTrackDCApTCutV0 = (TH2D*)fOutput->FindObject("TrackDCApTCutV0"); 
   fTrackProdRnoCut = (TH2D*)fOutput->FindObject("TrackProdRnoCut"); 
   fTrackProdRnoCutDalitz = (TH2D*)fOutput->FindObject("TrackProdRnoCutDalitz"); 
   fTrackProdRnoCutV0 = (TH2D*)fOutput->FindObject("TrackProdRnoCutV0");
   fTrackProdRwithCut = (TH2D*)fOutput->FindObject("TrackProdRwithCut"); 
   fTrackDCAerrpTFake = (TH2D*)fOutput->FindObject("TrackDCAerrpTFake"); 
   fTrackDCAerrpTTrueConv = (TH2D*)fOutput->FindObject("TrackDCAerrpTTrueConv"); 
   fTrackEta = (TH2D*)fOutput->FindObject("TrackEta"); 
   fConvPhi = (TH1D*)fOutput->FindObject("ConvPhi"); 
   fTOFPhi = (TH1D*)fOutput->FindObject("TOFPhi");
   fEventNr = (TH1D*)fOutput->FindObject("EventNr");
   fGaussTest = (TH1D*)fOutput->FindObject("GaussTest");
   fTPCClusters = (TH2D*)fOutput->FindObject("TPCClusters");
   fITSClusters = (TH2D*)fOutput->FindObject("ITSClusters");
   fV0sConvpT = (TH1D*)fOutput->FindObject("V0sConvpT");
   fV0sDalitzpT = (TH1D*)fOutput->FindObject("V0sDalitzpT");
   fV0sHFpT = (TH1D*)fOutput->FindObject("V0sHFpT");
   fVertexXY = (TH2D*)fOutput->FindObject("VertexXY");
   fConvEta = (TH2D*)fOutput->FindObject("ConvEta");
   fOtherEta = (TH2D*)fOutput->FindObject("OtherEta");
   fK0pTDCA = (TH2D*)fOutput->FindObject("K0pTDCA");
   fK0pTDCARCut = (TH2D*)fOutput->FindObject("K0pTDCARCut");
   fRKaons = (TH1D*)fOutput->FindObject("RKaons");
   fClusterPIDCluster = (TH2D*)fOutput->FindObject("ClusterPIDCluster");
   fBeautyMotherCorrelation = (TH3D*)fOutput->FindObject("BeautyMotherCorrelation");
   fBeautyMotherCorrelationRAA = (TH3D*)fOutput->FindObject("BeautyMotherCorrelationRAA");
   fBeautyMotherCorrelationHalfRAA = (TH3D*)fOutput->FindObject("BeautyMotherCorrelationHalfRAA");
   fCharmMotherCorrelation = (TH3D*)fOutput->FindObject("CharmMotherCorrelation");
   fCharmMotherCorrelationRAA = (TH3D*)fOutput->FindObject("CharmMotherCorrelationRAA");
   fCharmMotherCorrelationHalfRAA = (TH3D*)fOutput->FindObject("CharmMotherCorrelationHalfRAA");
   fTrackDCApTnoCut12 = (TH3D*)fOutput->FindObject("TrackDCApTnoCut12");
   fBeautyMotherCorrelation12 = (TH2D*)fOutput->FindObject("BeautyMotherCorrelation12");
   fBeautyMotherCorrelationRAA12 = (TH2D*)fOutput->FindObject("BeautyMotherCorrelationRAA12");
   fBeautyMotherCorrelationHalfRAA12 = (TH2D*)fOutput->FindObject("BeautyMotherCorrelationHalfRAA12");
   fCharmMotherCorrelation12 = (TH2D*)fOutput->FindObject("CharmMotherCorrelation12");
   fCharmMotherCorrelationRAA12 = (TH2D*)fOutput->FindObject("CharmMotherCorrelationRAA12");
   fCharmMotherCorrelationHalfRAA12 = (TH2D*)fOutput->FindObject("CharmMotherCorrelationHalfRAA12");
   
  
   
   std::cout << "Writing output to " << histname.Data() << std::endl;
   TFile out(histname.Data(),"recreate");
   out.cd();
   fhMultiplicity->Write();
   fhMultiplicityPle->Write();
   fTrackDCApTnoCut->Write();
   fTrackDCApTnoCutnoTOF->Write();
   fTrackDCApTCut->Write();
   fTrackDCApTnoCutV0->Write();
   fTrackDCApTnoCutV0weakerConditions->Write();
   fTrackDCApTnoCutV0weakerConditionsNonFake->Write();
   fTrackDCApTnoCutV0weakerConditionsFake->Write();
   fTrackDCApTnoCutV0weakerConditionsDalitz->Write();
   fTrackDCApTnoCutStrange->Write();
   fTrackDCApTnoCutDalitzWOStrange->Write();
   fTrackDCApTnoCutConversionWOStrange->Write();
   fTrackDCApTCutV0->Write();
   fTrackDCAerrpT->Write();
   fTrackDCAerrpTConv->Write();
   fTrackProdRnoCut->Write();
   fTrackProdRnoCutDalitz->Write();
   fTrackProdRnoCutV0->Write();
   fTrackProdRwithCut->Write();
   fTrackDCAerrpTFake->Write();
   fTrackDCAerrpTTrueConv->Write();
   fTrackEta->Write();
   fConvPhi->Write();
   fTOFPhi->Write();
   fEventNr->Write();
   fGaussTest->Write();
   fTPCClusters->Write();
   fITSClusters->Write();
   fV0sConvpT->Write();
   fV0sDalitzpT->Write();
   fV0sHFpT->Write();
   fVertexXY->Write();
   fConvEta->Write();
   fOtherEta->Write();
   fK0pTDCA->Write();
   fK0pTDCARCut->Write();
   fRKaons->Write();
   fClusterPIDCluster->Write();
   fBeautyMotherCorrelation->Write();
   fBeautyMotherCorrelationRAA->Write();
   fBeautyMotherCorrelationHalfRAA->Write();
   fCharmMotherCorrelation->Write();
   fCharmMotherCorrelationRAA->Write();
   fCharmMotherCorrelationHalfRAA->Write();
   fTrackDCApTnoCut12->Write();
   fBeautyMotherCorrelation12->Write();
   fBeautyMotherCorrelationRAA12->Write();
   fBeautyMotherCorrelationHalfRAA12->Write();
   fCharmMotherCorrelation12->Write();
   fCharmMotherCorrelationRAA12->Write();
   fCharmMotherCorrelationHalfRAA12->Write();
   out.Close();
   std::cout << "output written " << std::endl;

}
