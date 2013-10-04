void printProfile(char* g3 = "geant3/timing_200_woZDC.root", char* g4 ="geant4/timing_200_woZDC.root")
{
  // Create monitors 
  AliTransportMonitor* mon3 = AliTransportMonitor::Import(g3);
  AliTransportMonitor* mon4 = AliTransportMonitor::Import(g4);
  //
  // Create timing histos and volume lists
  mon3->Print();
  tH3 = (TH1F*) volume_timing->Clone();
  TObjArray* volumes3 = mon3->GetVolumes(); 
  mon4->Print();
  tH4 = (TH1F*) volume_timing->Clone();
  TObjArray* volumes4 = mon4->GetVolumes(); 
  //
  // Normalise the geant 3 histo
  // Fraction of time in %
  Float_t ri = tH3->Integral();
  tH3->Scale(100./ri);
  //
  Float_t cumulant = 0.;  // cumulant
  Int_t   i = 1;          // position in list
  //
  // print up to 90% of total time
  //
  printf(" Pos.               Volume  t (%) ct (%)  Time_G3  Time_G4 perStepG3 perStepG4   StepsG3   StepsG4    te+e-     tgam     thad     te+e-    tgam     thad\n"); 
  while (cumulant < 90.) {
    Float_t rt = tH3->GetBinContent(i);
    cumulant += rt;
    // extract the corresponding geant3 and geant4 volumes
    char* volN = (tH3->GetXaxis())->GetBinLabel(i);
  AliTransportMonitor::AliTransportMonitorVol* vol4 = (AliTransportMonitor::AliTransportMonitorVol*) 
    volumes4->FindObject(volN);
  AliTransportMonitor::AliTransportMonitorVol* vol3 = 0;
    if (i > 1) {
      vol3 = (AliTransportMonitor::AliTransportMonitorVol*) 
	volumes3->FindObject(volN);
    } else {
      vol3 = (AliTransportMonitor::AliTransportMonitorVol*) 
	volumes3->At(1);
     }
    TString s(volN);
    TString s20(s(0,20));
    // Total time
    
    Float_t t4[5]; 
    t4[0] = vol4->GetTotalTime();
    Float_t t3[5]; 
    t3[0] = vol3->GetTotalTime();

    // For e+/e-
    t4[1] = ElectronTime(vol4);
    t3[1] = ElectronTime(vol3);
 
    // For photons
    t4[2] = PhotonTime(vol4);
    t3[2] = PhotonTime(vol3);
   
    // Hadrons
    t4[3] = HadronTime(vol4);
    t3[3] = HadronTime(vol3);
   
    printf("%5d %20s %6.2f %6.2f %8.3f %8.3f %8.3e %8.3e %8.3e %8.3e %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n", 
	   i,
	   s20.Data(),
	   rt, cumulant, t3[0], t4[0], t3[0]/(1.+vol3->GetNSteps()), 
	   t4[0]/(1.+vol4->GetNSteps()), vol3->GetNSteps(), vol4->GetNSteps(),
	   t3[1]/t3[0] * 100., t3[2]/t3[0] * 100., t3[3]/t3[0] * 100.,
	   t4[1]/t4[0] * 100., t4[2]/t4[0] * 100., t4[3]/t4[0] * 100.

	   );
    i++;
  }
}

Float_t ElectronTime(AliTransportMonitor::AliTransportMonitorVol* vol)
{
   Int_t nt = vol->GetNtypes();
   Float_t t = 0.;
   for (Int_t i = 0; i < nt; i++) {
     Int_t pdg = TMath::Abs(vol->GetPDG(i));
     if (pdg == 11) t += vol->GetTime(i);
   }

   return t;
}

Float_t PhotonTime(AliTransportMonitor::AliTransportMonitorVol* vol)
{
   Int_t nt = vol->GetNtypes();
   Float_t t = 0.;
   for (Int_t i = 0; i < nt; i++) {
     Int_t pdg = TMath::Abs(vol->GetPDG(i));
     if (pdg == 22) t += vol->GetTime(i);
   }
   return t;
}

Float_t HadronTime(AliTransportMonitor::AliTransportMonitorVol* vol)
{
   Int_t nt = vol->GetNtypes();
   Float_t t = 0.;
   for (Int_t i = 0; i < nt; i++) {
     Int_t pdg = TMath::Abs(vol->GetPDG(i));
     if (pdg == 211 || pdg == 2212 || pdg == 2112 || pdg == 130 || pdg == 321) t += vol->GetTime(i);
   }
   return t;
}


