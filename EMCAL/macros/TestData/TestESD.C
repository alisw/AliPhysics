/// \file TestESD.C
/// \brief Example to analyze calorimeter ESDs
///
/// Example macro to extract calorimeter related information from ESDs.
/// Very similar for AODs (see TestAOD.C).
/// It mostly prints clusters/cells information but it can also plot some example 
/// histograms with clusters energy, position, and track matching.
///
/// Different global bools can be set to inspect different parameters: cells, trigger, 
/// track-cluster matching, cells in clusters, PID or MC origin of the clusters.
///
/// \author Gustavo Conesa Balbastre, <Gustavo.Conesa.Balbastre@cern.ch>, LPSC-CNRS
///

#if !defined(__CINT__) || defined(__MAKECINT__)

//Root include files 
#include <Riostream.h>
#include <TFile.h>
//#include <TSystem.h>
#include <TH1F.h>
#include <TParticle.h>
#include <TRefArray.h>

//AliRoot include files 
#include "AliRunLoader.h"
#include "AliStack.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliESDCaloCluster.h"
#include "AliESDCaloCells.h"
#include "AliESDCaloTrigger.h"
#include "AliPID.h"
#include "AliEMCALGeometry.h"

#endif

// Change the bool depending on what information you want to print
// when all FALSE, prints minimum cluster information.
Bool_t kPrintKine         = kFALSE; /// Print MC related information. Do not use for raw data.
Bool_t kPrintCaloCells    = kFALSE; /// Print cells parameters
Bool_t kPrintCaloTrigger  = kFALSE; /// Print trigger patches information
Bool_t kPrintTrackMatches = kFALSE; /// Print cluster-track matching information
Bool_t kPrintClusterCells = kFALSE; /// Print cells in clusters information
Bool_t kPrintClusterPID   = kFALSE; /// Print clusters PID (bayesian) weights
Bool_t kPrintMisalMatrix  = kFALSE; /// Print the alignment matrices stored in ESDs
 
///
/// Main method to read information stored in AliESDCaloClusters and AliESDCaloCells
///
void TestESD() 
{	
  // Init some example histograms
  // ESD
  TH1F * hEta  = (TH1F*) new TH1F("hEta","reco #eta",1000, -0.71,0.71); 
  TH1F * hPhi  = (TH1F*) new TH1F("hPhi","reco #phi",360, 0,360); 
  TH1F * hE    = (TH1F*) new TH1F("hE"  ,"reco e",300, 0,30); 
  TH1F * hTime = (TH1F*) new TH1F("hTime"  ,"reco time",1000, 0,1000); 
  hEta ->SetXTitle("#eta");
  hPhi ->SetXTitle("#phi (deg)");
  hE   ->SetXTitle("E (GeV)");
  hTime->SetXTitle("time (ns)");
  
  // Monte Carlo
  TH1F * hMCEta = (TH1F*) new TH1F("hMCEta","MC #eta",1000, -0.71,0.71); 
  TH1F * hMCPhi = (TH1F*) new TH1F("hMCPhi","MC #phi",360, 0,360); 
  TH1F * hMCE   = (TH1F*) new TH1F("hMCE"  ,"MC e",300, 0,30); 
  hMCEta->SetXTitle("#eta");
  hMCPhi->SetXTitle("#phi (deg)");
  hMCE  ->SetXTitle("E (GeV)");
  
  // ESD - MonteCarlo
  TH1F * hMCDEta = (TH1F*) new TH1F("hMCDEta"," #eta cluster - #eta MC",500, -0.05,0.05); 
  TH1F * hMCDPhi = (TH1F*) new TH1F("hMCDPhi"," #phi cluster - #phi MC",500, -0.05,0.05); 
  TH1F * hMCDE   = (TH1F*) new TH1F("hMCDE"  ,"e cluster - e MC",200, -10,10); 
  hMCDEta->SetXTitle("#eta_{reco}-#eta_{MC}");
  hMCDPhi->SetXTitle("#phi_{reco}-#phi_{MC} (rad)");
  hMCDE  ->SetXTitle("E_{reco}-E_{MC} (GeV)");

  // ESD - Track-matching
  TH1F * hTMDEta   = (TH1F*) new TH1F("hTMDEta"," #eta cluster - eta track",500, -0.05,0.05); 
  TH1F * hTMDPhi   = (TH1F*) new TH1F("hTMDPhi"," #phi cluster - phi track",500, -0.05,0.05);     
  TH1F * hTMDEtaOut= (TH1F*) new TH1F("hTMDEtaOut"," #eta cluster - eta track-outer",500, -0.05,0.05); 
  TH1F * hTMDPhiOut= (TH1F*) new TH1F("hTMDPhiOut"," #phi cluster - phi track-outer",500, -0.05,0.05);   
  TH1F * hTMResEta = (TH1F*) new TH1F("hTMResEta"," Residual #eta cluster - eta track",500, -0.05,0.05); 
  TH1F * hTMResPhi = (TH1F*) new TH1F("hTMResPhi"," Residual #phi cluster - phi track",500, -0.05,0.05); 
  hTMDEta   ->SetXTitle("#eta_{cluster}-#eta_{track}");
  hTMDPhi   ->SetXTitle("#phi_{cluster}-#phi_{track} (rad)");  
  hTMDEtaOut->SetXTitle("#eta_{cluster}-#eta_{track-out}");
  hTMDPhiOut->SetXTitle("#phi_{cluster}-#phi_{track-out} (rad)");
  
  TH1F * hTMEOverP    = (TH1F*) new TH1F("hTMEOverP"," E_{cluster}/ p_{Track}",100,0,2); 
  TH1F * hTMEOverPOut = (TH1F*) new TH1F("hTMEOverPOut"," E_{cluster}/ p_{Track-out}",100,0,2); 
  hTMEOverP    ->SetXTitle("E_{cluster}/ p_{Track}");
  hTMEOverPOut ->SetXTitle("E_{cluster}/ p_{Track-out}");
  
  // L1 trigger 
  TH2I * hEGAPatch = new TH2I("hEGAPatch","EGA trigger",51,-0.5,50.5,65,-0.5,64.5);
  TH2I * hEJEPatch = new TH2I("hEJEPatch","EJE trigger",51,-0.5,50.5,65,-0.5,64.5);
  hEGAPatch->SetXTitle("column (#eta direction)");
  hEGAPatch->SetYTitle("row (#phi direction)");
  hEJEPatch->SetXTitle("column (#eta direction)");
  hEJEPatch->SetYTitle("row (#phi direction)");
  
  // Open the ESD file, get the tree with events
  TFile* f = new TFile("AliESDs.root");
  TTree* esdTree = (TTree*)f->Get("esdTree");
  AliESDEvent* esd = new AliESDEvent();
  esd->ReadFromTree(esdTree);
  
  // Init geometry and array that will contain the clusters, Run2
  AliEMCALGeometry *geom =  AliEMCALGeometry::GetInstance("EMCAL_COMPLETE12SMV1_DCAL_8SM") ;  	
  
  // Run1:
  //AliEMCALGeometry *geom =  AliEMCALGeometry::GetInstance("EMCAL_COMPLETE12SMV1") ;  	
  //AliEMCALGeometry *geom =  AliEMCALGeometry::GetInstance("EMCAL_COMPLETEV1") ;  	
  //AliEMCALGeometry *geom =  AliEMCALGeometry::GetInstance("EMCAL_FIRSTYEARV1") ;  	
  
  // Pass the geometry transformation matrix from ESDs to geometry
  for(Int_t mod = 0; mod < (geom->GetEMCGeometry())->GetNumberOfSuperModules(); mod++)
  {
    //printf("matrix %d\n",mod);
    if(esd->GetEMCALMatrix(mod)) 
    {
      if(kPrintMisalMatrix) 
      {
        printf("Misalign matrid for: mod %d, matrix %p\n",mod, esd->GetEMCALMatrix(mod));

        (esd->GetEMCALMatrix(mod))->Print("");
      }
      
      // Fill it in case cell index to global position is requested
      geom->SetMisalMatrix(esd->GetEMCALMatrix(mod),mod) ;
    }//matrix
  }//module

  // Loop of events

  Float_t pos[3];

  TRefArray* caloClusters = new TRefArray();

  Int_t nEvt = esdTree->GetEntries();
  
  for(Int_t iev = 0; iev < nEvt; iev++) 
  {
    cout << "<<<< Event: " << iev+1 << "/" << nEvt << " >>>>"<<endl;
    esdTree->GetEvent(iev);
        
    // In case you want to play with MC data, get stack
    AliStack *stack = 0;
    if(kPrintKine)
    {  
      AliRunLoader *rl = AliRunLoader::Open("galice.root",AliConfig::GetDefaultEventFolderName(),  "read");
      rl->LoadKinematics();  
      rl->GetEvent(iev);
      stack = rl->Stack();
    }  
    
    // Get reconstructed vertex position, only used if 
    // the momentum of the cluster is calculated, see comment below
    //
    //Double_t vertex_position[3];
    //esd->GetVertex()->GetXYZ(vertex_position);
    
    //------------------------------------------------------
    // Get Cells Array, all cells in event, print cells info 
    //------------------------------------------------------ 
    AliVCaloCells &cells = *(esd->GetEMCALCells());
    
    if(kPrintCaloCells)
    {  
      Int_t nTotalCells = cells.GetNumberOfCells() ;  
      //Int_t type        = cells.GetType();
      for (Int_t icell=  0; icell <  nTotalCells; icell++) 
      {
        cout<<"Cell   : "<<icell<<"/"<<nTotalCells<<" - ID: "<<cells.GetCellNumber(icell)<<"; Amplitude: "<<cells.GetAmplitude(icell)<<"; Time: "<<cells.GetTime(icell)*1e9;
        cout << "; MC label "<<cells.GetMCLabel(icell)<<"; Embeded E fraction "<<cells.GetEFraction(icell);
        cout<<endl;	       
      }// cell loop
    }
    
    //------------------------------------------------------
    // Calo Trigger, L1 
    //------------------------------------------------------
    
    if(kPrintCaloTrigger)
    { 
      // Bits to know what L1 trigger, 
      Int_t bitEGA = 6; // 4 for old data
      Int_t bitEJE = 8; // 5 for old data
      
      AliESDCaloTrigger& trg = *(esd->GetCaloTrigger("EMCAL"));
		  
      trg.Reset();
      while (trg.Next())
      {
        Int_t col, row;
        trg.GetPosition(col, row);
  
        if (col > -1 && row > -1) 
        {
          Int_t bit = 0;
          trg.GetTriggerBits(bit);
          
          Int_t ts = 0;
          trg.GetL1TimeSum(ts);

          // Check if it is EGA or EJE (it does not work??)
          Bool_t isEGA1 = ((bit>>  bitEGA   ) & 0x1);
          Bool_t isEGA2 = ((bit>> (bitEGA+1)) & 0x1);
          Bool_t isEJE1 = ((bit>>  bitEJE   ) & 0x1);
          Bool_t isEJE2 = ((bit>> (bitEJE+1)) & 0x1);
          
          if ( ts >= 0 )
          {
            printf("Bit %d (G1 %d,G2 %d,J1 %d,J2 %d), cell in patch position (ieta,iphi)=(%d,%d), L1 amplitude %d\n ",
                            bit, isEGA1, isEGA2, isEJE1, isEJE2, col, row, ts);
            if(bit>10) hEJEPatch->Fill(col,row);
            else       hEGAPatch->Fill(col,row);
          }
        }
      }
    }
    
    //------------------------------------------------------
    // Calo Clusters 
    //------------------------------------------------------
    
    // Get CaloClusters Array
    esd->GetEMCALClusters(caloClusters);
    
    // Loop over clusters
    Int_t nclus = caloClusters->GetEntries();
    for (Int_t icl = 0; icl < nclus; icl++) 
    {
      AliVCluster* clus = (AliVCluster*)caloClusters->At(icl);
      Float_t energy = clus->E();
      clus->GetPosition(pos);
      TVector3 vpos(pos[0],pos[1],pos[2]);
      
      // We can get a momentum TLorentzVector per cluster, 
      // corrected by the vertex position, see above 
      //
      //TLorentzVector p;
      //clus->GetMomentum(p,vertex_position);
      
      Double_t cphi = vpos.Phi();
      if(cphi < 0) cphi +=TMath::TwoPi();
      Double_t ceta = vpos.Eta();
      
      Int_t nMatched   = clus->GetNTracksMatched();
      Int_t trackIndex = clus->GetTrackMatchedIndex();
      Int_t nLabels    = clus->GetNLabels();
      Int_t labelIndex = clus->GetLabel();
      Int_t nCells     = clus->GetNCells();
      Int_t nlm        = clus->GetNExMax() ;
      
      // Fill some histograms
      hEta ->Fill(ceta);
      hPhi ->Fill(cphi*TMath::RadToDeg());
      hE   ->Fill(energy);
      hTime->Fill(clus->GetTOF()*1e9);
      
      // Print basic cluster information
      cout << "Cluster: " << icl+1 << "/" << nclus << " Energy: " << energy << "; Phi: " 
      << cphi*TMath::RadToDeg() << "; Eta: " << ceta 
      << "; NCells: " << nCells  << "; NLM: " << nlm
      << "; #Labels: " << nLabels << " Index: " 
      << labelIndex << "; Time "<<clus->GetTOF()*1e9<<" ns "<<endl;
      
      if(nMatched > 0)
      {
        printf("\t N matches %d, match index %d, Residual phi %2.4f, eta %2.4f\n",nMatched,trackIndex,clus->GetTrackDx(),clus->GetTrackDz());
        hTMResEta->Fill(clus->GetTrackDz());
        hTMResPhi->Fill(clus->GetTrackDx());
      }
      
      //Print primary info
      if(stack && kPrintKine) 
      {
        if(labelIndex >= 0 && labelIndex < stack->GetNtrack())
        {
          TParticle * particle = stack->Particle(labelIndex);
          // Fill histograms with primary info
          hMCEta  ->Fill(particle->Eta());
          hMCPhi  ->Fill(particle->Phi()*TMath::RadToDeg());
          hMCE    ->Fill(particle->Energy());
          hMCDEta ->Fill(ceta-particle->Eta());
          hMCDPhi ->Fill(cphi-particle->Phi());
          hMCDE   ->Fill(energy-particle->Energy());
          
          // Print primary values
          cout<<"         MC More  contributing primary: "<<particle->GetName()<<"; with kinematics: "<<endl;
          cout<<" \t     Energy: "<<particle->Energy()<<"; Phi: "<<particle->Phi()*TMath::RadToDeg()<<"; Eta: "<<particle->Eta()<<endl;   
          cout<<" \t     dE: "<<energy-particle->Energy()<<"; dPhi: "<<(cphi-particle->Phi())*TMath::RadToDeg()<<"; dEta: "<<ceta-particle->Eta()<<endl;         
          cout<<" \t     deposited energy fraction = "<<clus->GetClusterMCEdepFraction(0)<<endl;

          //printf("List of e fractions %p; cell edep fraction map %p\n",clus->GetClusterMCEdepFraction(),clus->GetCellsMCEdepFractionMap());

          for(Int_t i = 1; i < nLabels; i++)
          {
            //particle = stack->Particle((((AliESDCaloCluster*)clus)->GetLabelsArray())->At(i));
            particle = stack->Particle((clus->GetLabels())[i]);
            //or Int_t *labels = clus->GetLabels();
            //particle = stack->Particle(labels[i]);
            cout<<"         Other contributing primary: "<<particle->GetName()<< "; Energy "<<particle->Energy()<<
            "; deposited energy fraction = "<<clus->GetClusterMCEdepFraction(i)<<endl;
          }
        }
        else if( labelIndex >= stack->GetNtrack()) 
          cout <<"PROBLEM, label is too large : "<<labelIndex<<" >= particles in stack "<< stack->GetNtrack() <<endl;
        else 
          cout<<"Negative label!!!  : "<<labelIndex<<endl;
      } // play with stack
      
      // Matching results
      if(kPrintTrackMatches && trackIndex >= 0) 
      {
        AliESDtrack* track = esd->GetTrack(trackIndex);
        if(track)
        {
          Double_t tphi = track->Phi();
          Double_t teta = track->Eta();
          Double_t tmom = track->P();

          Double_t deta = teta - ceta;
          Double_t dphi = tphi - cphi;

          if(dphi >  TMath::Pi()) dphi -= 2*TMath::Pi();
          if(dphi < -TMath::Pi()) dphi += 2*TMath::Pi();
          
          Double_t dR = sqrt(dphi*dphi + deta*deta);

          Double_t pOverE = tmom/energy;

          cout << "\t Track Momentum        : " << tmom << " phi: " << tphi << " eta: " << teta << endl;

          hTMEOverP->Fill(pOverE);
          hTMDEta  ->Fill(deta);        
          hTMDPhi  ->Fill(dphi);
          
          if(track->GetOuterParam())
          {
            tphi = track->GetOuterParam()->Phi();
            teta = track->GetOuterParam()->Eta();
            tmom = track->GetOuterParam()->P();
            
            cout << "\t Track Momentum - Outer: " << tmom << " phi: " << tphi << " eta: " << teta << endl;
            
            deta = teta - ceta;
            dphi = tphi - cphi;
            
            if(dphi >  TMath::Pi()) dphi -= 2*TMath::Pi();
            if(dphi < -TMath::Pi()) dphi += 2*TMath::Pi();
            
            dR = sqrt(dphi*dphi + deta*deta);
            
            pOverE = tmom/energy;
            
            if(dR < 0.02 && pOverE < 1.8 && nCells > 1)
            {
              cout << "\t Excellent MATCH! dR = " << dR << " p/E = " << pOverE << "Residuals phi "<< dphi*TMath::RadToDeg()<< " eta "<< deta  <<  endl;
            }
            
            hTMEOverPOut->Fill(pOverE);
            hTMDEtaOut->Fill(deta);        
            hTMDPhiOut->Fill(dphi);  
          }
          else printf("!!! NO OUTER PARAM !!!\n");
        }
        else printf("!!! NO TRACK !!!\n");
      }// matching
      
      // Get PID weights and print them
      if(kPrintClusterPID)
      {
        const Double_t *pid = clus->GetPID();
        printf("PID weights: ph %0.2f, pi0 %0.2f, el %0.2f, conv el %0.2f, hadrons: pion %0.2f, kaon %0.2f, proton %0.2f , neutron %0.2f, kaon %0.2f \n",
               pid[AliVCluster::kPhoton],   pid[AliVCluster::kPi0],
               pid[AliVCluster::kElectron], pid[AliVCluster::kEleCon],
               pid[AliVCluster::kPion],     pid[AliVCluster::kKaon],   pid[AliVCluster::kProton],
               pid[AliVCluster::kNeutron],  pid[AliVCluster::kKaon0]);
      } // PID
      
      // Get CaloCells of cluster and print their info, position.
      if(kPrintClusterCells)
      {	
        UShort_t * index    = clus->GetCellsAbsId() ;
        Double_t * fraction = clus->GetCellsAmplitudeFraction() ;
        Int_t sm = -1;
        
        for(Int_t i = 0; i < nCells ; i++)
        {
          Int_t    absId    = index[i]; // or clus->GetCellNumber(i) ;
          Double_t ampFract = fraction[i];
          Float_t  amp      = cells.GetCellAmplitude(absId) ;
          Double_t time     = cells.GetCellTime(absId);
          
          cout<<"\t Cluster Cell: AbsID : "<< absId << " == "<<clus->GetCellAbsId(i) <<"; Amplitude "<< amp << "; Fraction "<<ampFract<<"; Time " <<time*1e9<<endl;

          if(kPrintKine)
          {
            Float_t eDepFrac[4];
            clus->GetCellMCEdepFractionArray(i,eDepFrac);
            Int_t map = 0;
            if(clus->GetCellsMCEdepFractionMap()) map = clus->GetCellsMCEdepFractionMap()[i];
            printf("\t \t cell MC label index %d, cell MC deposited energy map %d: p1 %2.2f, p2 %2.2f, p3 %2.2f p4 %2.2f\n",
                   cells.GetCellMCLabel(absId),map,eDepFrac[0],eDepFrac[1],eDepFrac[2],eDepFrac[3]);
          }
          
          // Geometry methods  
          Int_t iSupMod =  0 ;
          Int_t iTower  =  0 ;
          Int_t iIphi   =  0 ;
          Int_t iIeta   =  0 ;
          Int_t iphi    =  0 ;
          Int_t ieta    =  0 ;
          
          if(geom)
          {
            geom->GetCellIndex(absId,iSupMod,iTower,iIphi,iIeta); 
            //Gives SuperModule and Tower numbers
            geom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,
                                              iIphi, iIeta,iphi,ieta);
            //Gives label of cell in eta-phi position per each supermodule
            Float_t cellPhi = 0;
            Float_t cellEta = 0;
            geom->EtaPhiFromIndex(absId,cellEta,cellPhi);
            cout<< "                SModule "<<iSupMod<<"; Tower "<<iTower <<"; Eta "<<iIeta
            <<"; Phi "<<iIphi<<"; Index: Cell Eta "<<ieta<<"; Cell Phi "<<iphi
            <<"; Global: Cell Eta "<<cellEta<<"; Cell Phi "<<cellPhi*TMath::RadToDeg()
            <<endl;
            
            if      ( i == 0       ) sm = iSupMod;
            else if ( sm != iSupMod) printf("                ******CLUSTER SHARED BY 2 SuperModules!!!!\n");
            	
          }// geometry on
        }// cluster cell loop
      }// print cell clusters
    } //cluster loop
  } // event loop
  
  
  //Write histograms in a file
  TFile * fhisto = (TFile*) new TFile("histos.root","recreate");
  hEta->Write();
  hPhi->Write();
  hE  ->Write();
  hTime->Write();
  
  if(kPrintKine)
  {
    hMCEta->Write();
    hMCPhi->Write();
    hMCE  ->Write();
    hMCDEta->Write();
    hMCDPhi->Write();
    hMCDE  ->Write();
  }
  
  if(kPrintTrackMatches)
  {
    hTMEOverP->Write();
    hTMDEta  ->Write();
    hTMDPhi  ->Write();
 
    hTMEOverPOut->Write();
    hTMDEtaOut  ->Write();
    hTMDPhiOut  ->Write();
    
    hTMResEta->Write();
    hTMResPhi->Write();  
  }
  
  if(kPrintCaloTrigger)
  {
    hEJEPatch->Write();
    hEGAPatch->Write();
  }
  
  fhisto->Close();
}
