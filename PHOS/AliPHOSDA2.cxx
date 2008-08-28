#include "AliPHOSDA2.h"
#include "TString.h"

ClassImp(AliPHOSDA2)

//----------------------------------------------------------------
AliPHOSDA2::AliPHOSDA2(int module) : TNamed(),
 fHistoFile(0),fFiredCells(0),fMod(module)

{
  // Create AliPHOSDA2 ("Bad channels finder") object.
  // module is the PHOS module number (0..4).
  // Quality histogram names: module_iX_iZ_gain.
  // Root file name: PHOS_ModuleX_BCM.root, where X - module number.
  
  char name[128];
  sprintf(name,"PHOS_Module%d_BCM",fMod);
  SetName(name);

  SetTitle("Detector Algorithm to check for PHOS channels quality");

  char rootname[128];
  sprintf(rootname,"%s.root",GetName());

  fHistoFile =  new TFile(rootname,"recreate"); // new file!
  
  for(Int_t iX=0; iX<64; iX++) {
    for(Int_t iZ=0; iZ<56; iZ++) {
      for(Int_t iGain=0; iGain<2; iGain++) {
	fHQuality[iX][iZ][iGain] = 0;
      }
    }
  }

  fMaps[0]=0;
  fMaps[1]=0;

  fFiredCells = new TH1I("fFiredCells","Number of fired cells per event",100,0,1000);

}

//-------------------------------------------------------------------
AliPHOSDA2::AliPHOSDA2(const AliPHOSDA2& da) : TNamed(da),
  fHistoFile(0),fFiredCells(0),fMod(da.fMod)
{
  // Copy constructor.

  char hname[128];
  TH1F* hist1=0;

  for(Int_t iX=0; iX<64; iX++) {
    for(Int_t iZ=0; iZ<56; iZ++) {
      for(Int_t iGain=0; iGain<2; iGain++) {

	sprintf(hname,"%d_%d_%d_%d",fMod,iX,iZ,iGain);
	hist1 = (TH1F*)da.fHistoFile->Get(hname);
	if(hist1) fHQuality[iX][iZ][iGain] = new TH1F(*hist1);
	else
	  fHQuality[iX][iZ][iGain] = 0;
      }
    }
  }
  
  if(da.fMaps[0]) 
    fMaps[0] = new TH2F(*da.fMaps[0]);
  else
    fMaps[0] = 0;

  if(da.fMaps[1]) 
    fMaps[1] = new TH2F(*da.fMaps[1]);
  else
    fMaps[1] = 0;
  
  fHistoFile = new TFile(da.GetName(),"recreate");
  fFiredCells = new TH1I(*da.fFiredCells);
  
}

//-------------------------------------------------------------------
AliPHOSDA2& AliPHOSDA2::operator= (const AliPHOSDA2& da)
{
  //Assignment operator.

  if(this != &da) {

    TString oldname(fHistoFile->GetName());
    TString newname(da.fHistoFile->GetName());

    if(oldname != newname) {
      delete fHistoFile;
      fHistoFile = new TFile(da.fHistoFile->GetName(),"update");
    }

    fMod = da.fMod;

    SetName(da.GetName());
    SetTitle(da.GetTitle());

    for(Int_t iX=0; iX<64; iX++) {
      for(Int_t iZ=0; iZ<56; iZ++) {
	for(Int_t iGain=0; iGain<2; iGain++) {
	  if (fHQuality[iX][iZ][iGain]) delete fHQuality[iX][iZ][iGain];
	  fHQuality[iX][iZ][iGain] = da.fHQuality[iX][iZ][iGain];
	}
      }
    }

    if(fMaps[0]) { 
      delete fMaps[0];
      fMaps[0] = da.fMaps[0];
    } 

    if(fMaps[1]) { 
      delete fMaps[1];
      fMaps[1] = da.fMaps[1];
    } 
    
    if(fFiredCells) {
      delete fFiredCells;
      fFiredCells = da.fFiredCells;
    }
    
  }
  
  return *this;
}


//-------------------------------------------------------------------
AliPHOSDA2::~AliPHOSDA2()
{
  // Destructor
  
  UpdateHistoFile();
  if(fHistoFile) delete fHistoFile;
  
}

//-------------------------------------------------------------------
void AliPHOSDA2::FillQualityHistograms(Float_t quality[64][56][2]) 
{
  // Fills quality histograms.
  // _By definition_, qood quality is 0<quality<1, 
  // all outside that is a bad quality.
  // If no quality value read for particular channel, 
  // the correspondent array entry should be filled by zero.
  // WARNING: this function should be called once per event!

  char hname[128];
  char htitl[128];

  for(Int_t iX=0; iX<64; iX++) {
    for (Int_t iZ=0; iZ<56; iZ++) {

      for(Int_t iGain=0; iGain<2; iGain++) {
	if(!quality[iX][iZ][iGain]) continue;
	
	if(fHQuality[iX][iZ][iGain]) 
	  fHQuality[iX][iZ][iGain]->Fill(quality[iX][iZ][iGain]);
	else {
	  sprintf(hname,"%d_%d_%d_%d",fMod,iX,iZ,iGain);
	  sprintf(htitl,"Quality for crystal %d_%d_%d and gain %d",fMod,iX,iZ,iGain);
	  fHQuality[iX][iZ][iGain] = new TH1F(hname,htitl,100,1.e-6,10.);
	  fHQuality[iX][iZ][iGain]->Fill(quality[iX][iZ][iGain]);
	}
      }

    }
  }

}

//-------------------------------------------------------------------
void  AliPHOSDA2::FillFiredCellsHistogram(Int_t nCells)
{
  fFiredCells->Fill(nCells);
}

//-------------------------------------------------------------------
void AliPHOSDA2::UpdateHistoFile()
{
  // Write histograms to file

  if(!fHistoFile) return;
  if(!fHistoFile->IsOpen()) return;

  char titl[128];

  if(fMaps[0]) 
    fMaps[0]->Reset();
  else {
    sprintf(titl,"Quality map for Low gain");
    fMaps[0] = new TH2F("gmaplow",  titl, 64,0.,64.,56,0.,56.);
  }

  if(fMaps[1]) 
    fMaps[1]->Reset();
  else {
    sprintf(titl,"Quality map for High gain");
    fMaps[1] = new TH2F("gmaphigh", titl, 64,0.,64.,56,0.,56.);
  }
    
  TH1F* hist1=0;

  for(Int_t iX=0; iX<64; iX++) {
    for(Int_t iZ=0; iZ<56; iZ++) {

      for(Int_t iGain=0; iGain<2; iGain++) {
	hist1 = fHQuality[iX][iZ][iGain];
	if(hist1) { 
	  hist1->Write(hist1->GetName(),TObject::kWriteDelete);
	  Double_t mean = hist1->GetMean();
	  fMaps[iGain]->SetBinContent(iX+1,iZ+1,mean);
	}
      } 

    }
  }

  fMaps[0]->Write(fMaps[0]->GetName(),TObject::kWriteDelete);
  fMaps[1]->Write(fMaps[1]->GetName(),TObject::kWriteDelete);

  fFiredCells->Write();

}

