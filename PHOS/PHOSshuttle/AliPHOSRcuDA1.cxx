#include "AliPHOSRcuDA1.h"
#include "TString.h"

ClassImp(AliPHOSRcuDA1)

//----------------------------------------------------------------
AliPHOSRcuDA1::AliPHOSRcuDA1(Int_t module, Int_t rcu) : TNamed(),
	   fHistoFile(0),fMod(module),fRCU(rcu),fWriteToFile(kTRUE),fHistoArray()

{
  // Create DA1 ("Calibration DA") object.
  // module is the PHOS module number (0..4).
  // Checks existence of histograms which might have been left
  // from the previous runs to continue their filling.
  // Histogram names: module_iX_iZ_gain for TH2 and  module_iX_iZ for TH1.
  // Root file name: PHOS_ModuleX_RCUY_Calib.root, where X - module number,
  // Y - RCU number. If no RCU specified (rcu<0), file name is simply
  // PHOS_ModuleX_Calib.root.
  
  char name[128]; TString sname;

  if(rcu<0) { 
    sname="PHOS_Module%d_Calib";
    snprintf(name,sname.Length(),sname.Data(),fMod);
  }
  else {
    sname="PHOS_Module%d_RCU%d_Calib";
    snprintf(name,sname.Length(),sname.Data(),fMod,fRCU);
  }

  SetName(name);
  SetTitle("Calibration Detector Algorithm");

  char rootname[128]; 
  TString srootname="%s.root";
  snprintf(rootname,srootname.Length(),srootname.Data(),GetName());

  fHistoFile =  new TFile(rootname,"update");
  fHistoArray.SetName("histo_container");

  char hname[128];
  TH1F* hist1=0;
  TH2F* hist2=0;
  TString shname;

  for(Int_t iX=0; iX<64; iX++) {
    for(Int_t iZ=0; iZ<56; iZ++) {
      
      shname = "%d_%d_%d";
      snprintf(hname,shname.Length(),shname.Data(),fMod,iX,iZ);

      hist1 = (TH1F*)fHistoFile->Get(hname);
      if(hist1) { 
	fHgLgRatio[iX][iZ] = hist1;
	fHistoArray.Add(hist1);
      }
      else
	fHgLgRatio[iX][iZ] = 0;

      for(Int_t iGain=0; iGain<2; iGain++) {
	
	shname = "%d_%d_%d_%d";
	snprintf(hname,shname.Length(),shname.Data(),fMod,iX,iZ,iGain);
	
	hist2 = (TH2F*)fHistoFile->Get(hname);
	if(hist2) {
	  fTimeEnergy[iX][iZ][iGain] = hist2;
	  fHistoArray.Add(hist2);
	}
	else
	  fTimeEnergy[iX][iZ][iGain] = 0;
      }

    }
  }
 
}

//-------------------------------------------------------------------
AliPHOSRcuDA1::AliPHOSRcuDA1(Int_t module, Int_t rcu, TObjArray* oldTimeEnergy) :
  TNamed(),fHistoFile(0),fMod(module),fRCU(rcu),fWriteToFile(kFALSE),
  fHistoArray()
{
  // Constructor. 
  // oldTimeEnergy is an array of histograms kept from the previous run.
  // By default the final histograms will not be saved to the root file.

  char name[128]; TString sname;

  if(rcu<0) {
    sname="PHOS_Module%d_Calib";
    snprintf(name,sname.Length(),sname.Data(),fMod);
  }
  else {
    sname="PHOS_Module%d_RCU%d_Calib";
    snprintf(name,sname.Length(),sname.Data(),fMod,fRCU);
  }
  
  SetName(name);
  SetTitle("Calibration Detector Algorithm");
  
  fHistoArray.SetName("histo_container");
  
  char hname[128];
  TH1F* hist1=0;
  TH2F* hist2=0;
  TString shname;

  for(Int_t iX=0; iX<64; iX++) {
    for(Int_t iZ=0; iZ<56; iZ++) {
      
      shname = "%d_%d_%d";
      snprintf(hname,shname.Length(),shname.Data(),fMod,iX,iZ);
      
      if(oldTimeEnergy)
	hist1 = (TH1F*)oldTimeEnergy->FindObject(hname);
      if(hist1) {
        fHgLgRatio[iX][iZ] = hist1;
	fHistoArray.Add(hist1);
      }
      else
        fHgLgRatio[iX][iZ] = 0;
      
      for(Int_t iGain=0; iGain<2; iGain++) {
	
        shname = "%d_%d_%d_%d";
        snprintf(hname,shname.Length(),shname.Data(),fMod,iX,iZ,iGain);
	
	if(oldTimeEnergy)
	  hist2 = (TH2F*)oldTimeEnergy->FindObject(hname);
        if(hist2) {
          fTimeEnergy[iX][iZ][iGain] = hist2;
	  fHistoArray.Add(hist2);
        }
        else
          fTimeEnergy[iX][iZ][iGain] = 0;
      }

    }
  }

}

//-------------------------------------------------------------------
AliPHOSRcuDA1::~AliPHOSRcuDA1()
{
  // Destructor
  
  UpdateHistoFile();
  if(fHistoFile) delete fHistoFile;
  
}

//-------------------------------------------------------------------
void AliPHOSRcuDA1::FillHistograms(Float_t e[64][56][2], Float_t t[64][56][2]) 
{
  // Fill energy vs time-of-flight and HG/LG ratio histograms in one event.
  // Energy and TOF data are encoded as e[X][Z][gain] and t[X][Z][gain], 
  // where X(0..63) and Z(0..55) - crystal position in the module, 
  // gain=0 - low gain, gain=1 - high gain.
  // Energy in ADC counts, time in samples.
  // If no energy or time read for particular channel, 
  // the correspondent array entry should be filled by zero.
  // WARNING: this function should be called once per event!

  TString hname;
  TString htitle;

  for(Int_t iX=0; iX<64; iX++) {
    for (Int_t iZ=0; iZ<56; iZ++) {
      
      // HG/LG
      if(e[iX][iZ][0]>10. && e[iX][iZ][1]>10. && e[iX][iZ][1]<900.) {
	
	if(fHgLgRatio[iX][iZ]) {
	  //printf("iX=%d iZ=%d,e[iX][iZ][1]=%.3f,e[iX][iZ][0]=%.3f, t1=%.3f,t0=%.3f\n",
	  // 		 iX,iZ,e[iX][iZ][1],e[iX][iZ][0], t[iX][iZ][1],t[iX][iZ][0]);
	  fHgLgRatio[iX][iZ]->Fill(e[iX][iZ][1]/e[iX][iZ][0]); 
	}
	else
	  {
	    hname.Clear(); htitle.Clear();
	    hname += fMod; hname += "_"; hname += iX; hname += "_"; hname += iZ;
	    htitle += "HG/LG ratio for crystal "; htitle += hname;
	    
	    fHgLgRatio[iX][iZ] = new TH1F(hname,htitle,400,14.,18.);
// 	    printf("iX=%d iZ=%d,e[iX][iZ][1]=%.3f,e[iX][iZ][0]=%.3f\n",
// 		   iX,iZ,e[iX][iZ][1],e[iX][iZ][0]);
	    fHgLgRatio[iX][iZ]->Fill(e[iX][iZ][1]/e[iX][iZ][0]);
	    fHistoArray.Add(fHgLgRatio[iX][iZ]);
	  }
      }
	  
      // Energy vs TOF
      for(Int_t iGain=0; iGain<2; iGain++) {

	if(e[iX][iZ][iGain]<10) continue;
	if(!(t[iX][iZ][iGain]>0)) continue;

	if(fTimeEnergy[iX][iZ][iGain]) 
	  fTimeEnergy[iX][iZ][iGain]->Fill(e[iX][iZ][iGain],t[iX][iZ][iGain]);
	else 
	  {
	    hname.Clear(); htitle.Clear();
	    hname += fMod; hname += "_"; hname += iX; hname += "_"; hname += iZ;
	    htitle = "Energy vs TOF for crystal "; htitle += hname; htitle += " and gain "; htitle += iGain;
	    hname += "_"; hname += iGain;
	    
	    fTimeEnergy[iX][iZ][iGain] = new TH2F(hname,htitle,100,0.,1024.,50,0.,50.);
	    fTimeEnergy[iX][iZ][iGain]->Fill(e[iX][iZ][iGain],t[iX][iZ][iGain]);
	    fHistoArray.Add(fTimeEnergy[iX][iZ][iGain]);
	  }
      }
      
    }
  }
  
}

//-------------------------------------------------------------------
void AliPHOSRcuDA1::UpdateHistoFile()
{
  // Write histograms to file

  if(!fWriteToFile) return;
  if(!fHistoFile) return;
  if(!fHistoFile->IsOpen()) return;

  fHistoArray.Write();
  fHistoFile->Purge();

}

//-----------------------------------------------------------------------
void AliPHOSRcuDA1::SetWriteToFile(Bool_t write)
{
  if(!write) { 
    fWriteToFile = write;
    return;
  }

  if(!fHistoFile) {
    TString rootname(GetName());
    rootname += ".root";
    fHistoFile =  new TFile(rootname.Data(),"update");
  } 
  
  fWriteToFile = write;
}
