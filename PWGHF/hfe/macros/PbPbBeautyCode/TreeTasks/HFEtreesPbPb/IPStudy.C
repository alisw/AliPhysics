#define IPStudy_cxx
// Tree Analysis Task for Pb-Pb data reduced event trees

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
   Double_t ipBinningY[4001];
   for(int i=0; i<=4000; i++){ipBinningY[i]=-0.2 + double(i)*0.0001;}
   Double_t ipErrBinningY[401];
   for(int i=0; i<=400; i++){ipErrBinningY[i]=double(i)*0.0001;}
   Double_t mcTruthBinning[19] = {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5};
   fTrackDCApTnoCut = new TH3D("TrackDCApTnoCut","", 18, ptbinningX, 4000, ipBinningY, 6, mcTruthBinning);
   fTrackDCApTCut = new TH3D("TrackDCApTCut","", 18, ptbinningX, 4000, ipBinningY, 6, mcTruthBinning);
   fTrackDCAerrpT = new TH3D("TrackDCAerrpT","", 18, ptbinningX, 400, ipErrBinningY, 6, mcTruthBinning);
   fTrackDCApTnoCutHadron = new TH2D("TrackDCApTnoCutHadron","", 18, ptbinningX, 4000, ipBinningY);
   
   fTrackDCAerr = new TH1D("TrackDCAerr","", 100,0.0,0.03);
   fhMultiplicity= new TH1D("hMultiplicity","", 1000,0.0,10000);
   
   fTrackDCAerrpTConv = new TH3D("TrackDCAerrpTConv","", 20, 0,5, 100,0.0,0.03, 2, -0.5, 1.5);
   fTrackProdRnoCut = new TH2D("TrackProdRnoCut","", 20, 0,5, 500, 0, 20.0);
   fTrackProdRnoCutV0 = new TH2D("TrackProdRnoCutV0","", 20, 0,5, 500, 0, 20.0);
   fTrackProdRwithCut = new TH2D("TrackProdRwithCut","", 20, 0,5, 500, 0, 20.0);
   fTrackDCAerrpTFake = new TH2D("TrackDCAerrpTFake","", 20, 0,5, 100,0.0,0.03);
   fTrackDCAerrpTTrueConv = new TH2D("TrackDCAerrpTTrueConv","", 20, 0,5, 100,0.0,0.03);
   //TrackEta = new TH2D("TrackEta","", 100, -1.0, 1.0, 10, -0.5, 9.5);
   fEventNr = new TH1D("EventNr", "", 5, -0.5, 4.5);
   fNContribVertex = new TH2D("NContribVertex", "", 100, -0.5, 3999.5, 100, 0., 100.);
   fDeuteronTOF = new TH2D("DeuteronTOF","", 20, 0,5, 200,-20.0,20.0);
   
   //TPCTOFsigmas = new TH3D("TPCTOFsigmas", "", 50, 0.0, 10., 800, -15.,5., 200, -5., 5.);
   fBin0Numbers = new TH1D("Bin0Numbers", "", 5, -0.5, 4.5); // 0: All Events 1: pass Z-Cut, 2: pass all cuts, 3: also nVertex>=20
   fTPCPionV0 = new TH2D("TPCPionV0","", 20, 0,5, 500, 0, 200.0);
   fTPCElectronV0 = new TH2D("TPCElectronV0","", 20, 0,5, 500, 0, 200.0);
   fTPCSignal = new TH2D("TPCSignal","", 20, 0,5, 500, 0, 200.0);
   fITSSignal = new TH2D("ITSSignal","", 20, 0,5, 500, -20, 20.0);
   fTPCPionV0nsigma = new TH2D("TPCPionV0nsigma","", 20, 0,5, 500, -15, 5.0);
   fTPCElectronV0nsigma = new TH2D("TPCElectronV0nsigma","", 20, 0,5, 500, -15, 5.0);
   fTPCElectronV0nsigmaNoCuts = new TH2D("TPCElectronV0nsigmaNoCuts","", 20, 0,5, 500, -15, 5.0);
   fVertexXY = new TH2D("VertexXY","", 100, -0.5, 0.5, 100, -0.5, 0.5);
   fK0pTDCA = new TH2D("K0pTDCA", "", 40, 0.0, 10., 200, -0.2, 0.2);
   fK0pTDCARCut = new TH2D("K0pTDCARCut", "", 40, 0.0, 10., 200, -0.2, 0.2);
   
   fOutput->Add(fTrackDCAerr);
   fOutput->Add(fhMultiplicity);
   fOutput->Add(fTrackDCAerrpT);
   fOutput->Add(fTrackDCAerrpTConv);
   fOutput->Add(fTrackDCApTnoCut);
   fOutput->Add(fTrackDCApTCut);
   fOutput->Add(fTrackDCApTnoCutHadron);
   fOutput->Add(fTrackProdRnoCut);
   fOutput->Add(fTrackProdRnoCutV0);
   fOutput->Add(fTrackProdRwithCut);
   fOutput->Add(fTrackDCAerrpTFake);
   fOutput->Add(fTrackDCAerrpTTrueConv);
   //fOutput->Add(TrackEta);
   //fOutput->Add(TPCTOFsigmas);
   fOutput->Add(fEventNr);
   fOutput->Add(fDeuteronTOF);
   fOutput->Add(fNContribVertex);
   fOutput->Add(fBin0Numbers);
   fOutput->Add(fVertexXY);
   fOutput->Add(fK0pTDCA);
   fOutput->Add(fK0pTDCARCut);
   
   fOutput->Add(fTPCPionV0);
   fOutput->Add(fTPCElectronV0);
   fOutput->Add(fTPCPionV0nsigma);
   fOutput->Add(fTPCElectronV0nsigma);
   fOutput->Add(fTPCElectronV0nsigmaNoCuts);
   fOutput->Add(fTPCSignal);
   fOutput->Add(fITSSignal);
   
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
  
  if(HFEevent->GetCentrality()<20. && HFEevent->GetCentrality()>0.)
  {
  fBin0Numbers->Fill(0);
  if(TMath::Abs(HFEevent->GetVZ()) < 10.)
  {
    fBin0Numbers->Fill(1);
    if(HFEevent->GetNContribVertex() && HFEevent->GetNContribVertexSPD() &&  TMath::Abs(HFEevent->GetVZ() - HFEevent->GetVZSPD()) < 0.5 && HFEevent->GetVertexZResolutionSPD() < 0.25)
    {
      fBin0Numbers->Fill(2);
      if(HFEevent->GetNContribVertexSPD()>20)
	fBin0Numbers->Fill(3);
    }
  }
  }
  
  if(EventPassesCuts())
  {
  int nTracks = HFEevent->GetNumberOfTracks();
  fEventNr->Fill(0);
  fNContribVertex->Fill(HFEevent->GetNContribVertex(), HFEevent->GetCentrality());
  AliHFEreducedTrack * hfeTrack; 
  fVertexXY->Fill(HFEevent->GetVX(), HFEevent->GetVY());
  fhMultiplicity->Fill(nTracks);
  
      double pT;
      double resolution;
    
      for(int j=0;j<nTracks;j++)
    {
      
      
      hfeTrack=(AliHFEreducedTrack*)HFEevent->GetTrack(j);
      resolution = hfeTrack->HFEImpactParameter() / hfeTrack->HFEImpactParameterResolution();
      
      if(hfeTrack->IsV0pion())  // GetV0prodR() might be an additional condition
	fK0pTDCA->Fill(hfeTrack->Pt(), hfeTrack->HFEImpactParameter()*hfeTrack->Charge());
      if(hfeTrack->IsV0pion())  // GetV0prodR() might be an additional condition
	if(hfeTrack->GetV0prodR()<5.0)
	  fK0pTDCARCut->Fill(hfeTrack->Pt(), hfeTrack->HFEImpactParameter()*hfeTrack->Charge());
      //source->Fill(hfeTrack->MCSource());
      //TrackQuantity->Fill(hfeTrack->GetITSnclusters());
	if(hfeTrack->IsV0electron() && hfeTrack->HasITScluster(0) && hfeTrack->HasITScluster(1))
	  fTPCElectronV0nsigmaNoCuts->Fill(hfeTrack->Pt(), hfeTrack->GetTPCsigmaElCorrected());
	
      if(passesCuts(hfeTrack))
      {  // 0,1,2,3=charm,beauty, conversion, dalitz
	fDeuteronTOF->Fill(hfeTrack->Pt(), hfeTrack->GetTOFsigmaDeuteron());
	//TPCTOFsigmas->Fill(hfeTrack->Pt() ,hfeTrack->GetTPCsigmaElCorrected(), hfeTrack->GetTOFsigmaEl());
	
	if(abs(hfeTrack->GetTOFsigmaEl())<3. && hfeTrack->IsV0pion())
	  fTPCPionV0->Fill(hfeTrack->Pt(), hfeTrack->GetTPCdEdxCorrected());
	if(abs(hfeTrack->GetTOFsigmaEl())<3. && hfeTrack->IsV0electron() && abs(hfeTrack->GetITSsigmaEl())<5.)
	  fTPCElectronV0->Fill(hfeTrack->Pt(), hfeTrack->GetTPCdEdxCorrected());
	if(abs(hfeTrack->GetTOFsigmaEl())<3. && hfeTrack->IsV0pion())
	  fTPCPionV0nsigma->Fill(hfeTrack->Pt(), hfeTrack->GetTPCsigmaElCorrected());
	if(abs(hfeTrack->GetTOFsigmaEl())<3. && hfeTrack->IsV0electron() && abs(hfeTrack->GetITSsigmaEl())<5.)
	  fTPCElectronV0nsigma->Fill(hfeTrack->Pt(), hfeTrack->GetTPCsigmaElCorrected());
	if(abs(hfeTrack->GetTOFsigmaEl()))
	  fTPCSignal->Fill(hfeTrack->Pt(), hfeTrack->GetTPCdEdxCorrected());
	if(abs(hfeTrack->GetTOFsigmaEl()))
	  fITSSignal->Fill(hfeTrack->Pt(), hfeTrack->GetITSsigmaEl());
	
	//if(abs(hfeTrack->GetTOFsigmaEl())<3. && hfeTrack->GetTPCsigmaElCorrected()>0.163)
	//if(abs(hfeTrack->GetTOFsigmaEl())<3. && hfeTrack->GetTPCsigmaElCorrected()>0. && hfeTrack->GetTPCsigmaElCorrected()<3.)
	if(hfeTrack->GetTPCsigmaElCorrected()>-0.5 && hfeTrack->GetTPCsigmaElCorrected()<3.)
	{
	  fTrackDCAerrpT->Fill(hfeTrack->Pt(), resolution, 1);
	
	  //if(hfeTrack->GetTPCnclusters()>80 && hfeTrack->Pt()>0.5) TrackEta->Fill(hfeTrack->Eta(), 1);

	  if(hfeTrack->IsV0electron())
	    fTrackDCApTnoCut->Fill(hfeTrack->Pt(),hfeTrack->HFEImpactParameter()*hfeTrack->Charge(), 2);
	  else
	    fTrackDCApTnoCut->Fill(hfeTrack->Pt(),hfeTrack->HFEImpactParameter()*hfeTrack->Charge(), 1);

	  pT=hfeTrack->Pt();
          if(resolution<=(0.007313/pT+0.00173))
	  {
	    if(hfeTrack->IsV0electron())
	      fTrackDCApTCut->Fill(hfeTrack->Pt(),hfeTrack->HFEImpactParameter()*hfeTrack->Charge(), 2);
	    else
	      fTrackDCApTCut->Fill(hfeTrack->Pt(),hfeTrack->HFEImpactParameter()*hfeTrack->Charge(), 1);
  
	  }
	  
	}
	
	//if(abs(hfeTrack->GetTOFsigmaEl())<3. && hfeTrack->GetTPCsigmaElCorrected()>-5. && hfeTrack->GetTPCsigmaElCorrected()<-3.)
	if(hfeTrack->GetTPCsigmaElCorrected()>-5. && hfeTrack->GetTPCsigmaElCorrected()<-3.)
	{
	  fTrackDCApTnoCutHadron->Fill(hfeTrack->Pt(),hfeTrack->HFEImpactParameter()*hfeTrack->Charge());
	}
	
      }
      
    }
  }


   return kTRUE;
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
   if(!(HFEevent->GetCentrality()>0.0 && HFEevent->GetCentrality()<20.)) return kFALSE;  // should be okay
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
      //if(!(TMath::Abs(track->Eta())>0.6)) return kFALSE;
      if(!track->HasTOFpid()) return kFALSE; 
      if(abs(track->GetTOFsigmaEl())>3.) return kFALSE; 
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
  
   fTrackDCAerr = (TH1D *)fOutput->FindObject("TrackDCAerr");
   fhMultiplicity = (TH1D*)fOutput->FindObject("hMultiplicity");
   fTrackDCAerrpT = (TH3D*)fOutput->FindObject("TrackDCAerrpT"); 
   fTrackDCAerrpTConv = (TH3D*)fOutput->FindObject("TrackDCAerrpTConv"); 
   fTrackDCApTnoCut = (TH3D*)fOutput->FindObject("TrackDCApTnoCut");
   fTrackDCApTCut = (TH3D*)fOutput->FindObject("TrackDCApTCut"); 
   fTrackProdRnoCut = (TH2D*)fOutput->FindObject("TrackProdRnoCut"); 
   fTrackProdRnoCutV0 = (TH2D*)fOutput->FindObject("TrackProdRnoCutV0");
   fTrackProdRwithCut = (TH2D*)fOutput->FindObject("TrackProdRwithCut"); 
   fTrackDCAerrpTFake = (TH2D*)fOutput->FindObject("TrackDCAerrpTFake"); 
   fTrackDCAerrpTTrueConv = (TH2D*)fOutput->FindObject("TrackDCAerrpTTrueConv"); 
   //TrackEta = (TH2D*)fOutput->FindObject("TrackEta"); 
   //TPCTOFsigmas = (TH3D*)fOutput->FindObject("TPCTOFsigmas");
   fEventNr = (TH1D*)fOutput->FindObject("EventNr");
   fDeuteronTOF = (TH2D*)fOutput->FindObject("DeuteronTOF");
   fNContribVertex = (TH2D*)fOutput->FindObject("NContribVertex");
   fTrackDCApTnoCutHadron = (TH2D*)fOutput->FindObject("TrackDCApTnoCutHadron");
   fBin0Numbers = (TH1D*)fOutput->FindObject("Bin0Numbers");
   fVertexXY = (TH2D*)fOutput->FindObject("VertexXY");
   fK0pTDCA = (TH2D*)fOutput->FindObject("K0pTDCA");
   fK0pTDCARCut = (TH2D*)fOutput->FindObject("K0pTDCARCut");
   
   fTPCPionV0 = (TH2D*)fOutput->FindObject("TPCPionV0");
   fTPCElectronV0 = (TH2D*)fOutput->FindObject("TPCElectronV0");
   fTPCPionV0nsigma = (TH2D*)fOutput->FindObject("TPCPionV0nsigma");
   fTPCElectronV0nsigma = (TH2D*)fOutput->FindObject("TPCElectronV0nsigma");
   fTPCElectronV0nsigmaNoCuts = (TH2D*)fOutput->FindObject("TPCElectronV0nsigmaNoCuts");
   fTPCSignal = (TH2D*)fOutput->FindObject("TPCSignal");
   fITSSignal = (TH2D*)fOutput->FindObject("ITSSignal");
   
   
      
  
   TString histname = Form("outputTreePbPbDataEta-05pos020.root");
   std::cout << "Writing output to " << histname.Data() << std::endl;
   TFile out(histname.Data(),"recreate");
   out.cd();
   fhMultiplicity->Write();
   fTrackDCApTnoCut->Write();
   fTrackDCApTCut->Write();
   fTrackDCApTnoCutHadron->Write();
   fTrackDCAerrpT->Write();
   fTrackDCAerrpTConv->Write();
   fTrackProdRnoCut->Write();
   fTrackProdRnoCutV0->Write();
   fTrackProdRwithCut->Write();
   fTrackDCAerrpTFake->Write();
   fTrackDCAerrpTTrueConv->Write();
   //TrackEta->Write();
   //TPCTOFsigmas->Write();
   fEventNr->Write();
   fDeuteronTOF->Write();
   fNContribVertex->Write();
   fBin0Numbers->Write();
   fVertexXY->Write();
   
   fTPCPionV0->Write();
   fTPCElectronV0->Write();
   fTPCPionV0nsigma->Write();
   fTPCElectronV0nsigma->Write();
   fTPCElectronV0nsigmaNoCuts->Write();
   fTPCSignal->Write();
   fITSSignal->Write();
   fK0pTDCA->Write();
   fK0pTDCARCut->Write();
   out.Close();
   std::cout << "output written " << std::endl;

}
