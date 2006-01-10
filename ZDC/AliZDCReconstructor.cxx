/**************************************************************************
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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for ZDC reconstruction                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


#include <TF1.h>

#include "AliRunLoader.h"
#include "AliRawReader.h"
#include "AliESD.h"
#include "AliZDCDigit.h"
#include "AliZDCRawStream.h"
#include "AliZDCReco.h"
#include "AliZDCReconstructor.h"
#include "AliZDCCalibData.h"


ClassImp(AliZDCReconstructor)


//_____________________________________________________________________________
AliZDCReconstructor:: AliZDCReconstructor()
{
  // **** Default constructor
  // if(!fStorage) fStorage =  AliCDBManager::Instance()->GetStorage("local://DBlocal");

  //  ---      Number of generated spectator nucleons and impact parameter
  // --------------------------------------------------------------------------------------------------
  // [1] ### Results from a new production  -> 0<b<18 fm (Apr 2002)
  // Fit results for neutrons (Nspectator n true vs. EZN)
  fZNCen = new TF1("fZNCen",
      "(-2.287920+sqrt(2.287920*2.287920-4*(-0.007629)*(11.921710-x)))/(2*(-0.007629))",0.,164.);
  fZNPer = new TF1("fZNPer",
      "(-37.812280-sqrt(37.812280*37.812280-4*(-0.190932)*(-1709.249672-x)))/(2*(-0.190932))",0.,164.);
  // Fit results for protons (Nspectator p true vs. EZP)
  fZPCen = new TF1("fZPCen",
       "(-1.321353+sqrt(1.321353*1.321353-4*(-0.007283)*(3.550697-x)))/(2*(-0.007283))",0.,60.);
  fZPPer = new TF1("fZPPer",
      "(-42.643308-sqrt(42.643308*42.643308-4*(-0.310786)*(-1402.945615-x)))/(2*(-0.310786))",0.,60.);
  // Fit results for total number of spectators (Nspectators true vs. EZDC)
  fZDCCen = new TF1("fZDCCen",
      "(-1.934991+sqrt(1.934991*1.934991-4*(-0.004080)*(15.111124-x)))/(2*(-0.004080))",0.,225.);
  fZDCPer = new TF1("fZDCPer",
      "(-34.380639-sqrt(34.380639*34.380639-4*(-0.104251)*(-2612.189017-x)))/(2*(-0.104251))",0.,225.);
  // --------------------------------------------------------------------------------------------------
  // Fit results for b (b vs. EZDC)
  // [2] ### Results from a new production  -> 0<b<18 fm (Apr 2002)
  fbCen = new TF1("fbCen","-0.056923+0.079703*x-0.0004301*x*x+0.000001366*x*x*x",0.,220.);
  fbPer = new TF1("fbPer","17.943998-0.046846*x+0.000074*x*x",0.,220.);
  // --------------------------------------------------------------------------------------------------
  // Evaluating Nspectators and b from ZEM energy
  // [2] ### Results from a new production  -> 0<b<18 fm (Apr 2002)
  fZEMn  = new TF1("fZEMn","126.2-0.05399*x+0.000005679*x*x",0.,4000.);
  fZEMp  = new TF1("fZEMp","82.49-0.03611*x+0.00000385*x*x",0.,4000.);
  fZEMsp = new TF1("fZEMsp","208.7-0.09006*x+0.000009526*x*x",0.,4000.);
  fZEMb  = new TF1("fZEMb","16.06-0.01633*x+1.44e-5*x*x-6.778e-9*x*x*x+1.438e-12*x*x*x*x-1.112e-16*x*x*x*x*x",0.,4000.);
}

//_____________________________________________________________________________
AliZDCReconstructor::AliZDCReconstructor(const AliZDCReconstructor& 
                                         reconstructor):
  AliReconstructor(reconstructor)
{
// copy constructor

  Fatal("AliZDCReconstructor", "copy constructor not implemented");
}

//_____________________________________________________________________________
AliZDCReconstructor& AliZDCReconstructor::operator = 
  (const AliZDCReconstructor& /*reconstructor*/)
{
// assignment operator

  Fatal("operator =", "assignment operator not implemented");
  return *this;
}

//_____________________________________________________________________________
AliZDCReconstructor::~AliZDCReconstructor()
{
// destructor

  delete fZNCen;
  delete fZNPer;
  delete fZPCen;
  delete fZPPer;
  delete fZDCCen;
  delete fZDCPer;
  delete fbCen;
  delete fbPer;
  delete fZEMn;
  delete fZEMp;
  delete fZEMsp;
  delete fZEMb;
}


//_____________________________________________________________________________
void AliZDCReconstructor::Reconstruct(AliRunLoader* runLoader) const
{
  // *** Local ZDC reconstruction for digits
  
  // Get calibration data
  int runNumber = 0;
  AliZDCCalibData *calibda = GetCalibData(runNumber);  
  
  Float_t meanPed[47];
  for(Int_t jj=0; jj<47; jj++) meanPed[jj] = calibda->GetMeanPed(jj);

  AliLoader* loader = runLoader->GetLoader("ZDCLoader");
  if (!loader) return;
  loader->LoadDigits("read");
  loader->LoadRecPoints("recreate");
  AliZDCDigit digit;
  AliZDCDigit* pdigit = &digit;

  // Event loop
  for (Int_t iEvent = 0; iEvent < runLoader->GetNumberOfEvents(); iEvent++) {
    runLoader->GetEvent(iEvent);

    // load digits
    loader->LoadDigits();
    TTree* treeD = loader->TreeD();
    if (!treeD) continue;
    treeD->SetBranchAddress("ZDC", &pdigit);

    // loop over digits
    Float_t zncorr=0, zpcorr=0, zemcorr=0;
    for (Int_t iDigit = 0; iDigit < treeD->GetEntries(); iDigit++) {
      treeD->GetEntry(iDigit);
      if (!pdigit) continue;

      if(digit.GetSector(0) == 1)
       	 zncorr  += (Float_t) (digit.GetADCValue(0)-meanPed[digit.GetSector(1)]);	   // ped 4 high gain ZN ADCs
      else if(digit.GetSector(0) == 2)
	 zpcorr  += (Float_t) (digit.GetADCValue(0)-meanPed[digit.GetSector(1)+10]); // ped 4 high gain ZP ADCs
      else if(digit.GetSector(0) == 3){
	 if(digit.GetSector(1)==1)      zemcorr += (Float_t) (digit.GetADCValue(0)-meanPed[digit.GetSector(1)+20]); // ped 4 high gain ZEM1 ADCs
	 else if(digit.GetSector(1)==2) zemcorr += (Float_t) (digit.GetADCValue(0)-meanPed[digit.GetSector(1)+22]); // ped 4 high gain ZEM2 ADCs
      }
    }
    if(zncorr<0)  zncorr=0;
    if(zpcorr<0)  zpcorr=0;
    if(zemcorr<0) zemcorr=0;

    // reconstruct the event
    printf("\n \t ZDCReco from digits-> Ev.#%d ZN = %.0f, ZP = %.0f, ZEM = %.0f\n",iEvent,zncorr,zpcorr,zemcorr);
    ReconstructEvent(loader, zncorr, zpcorr, zemcorr);
  }

  loader->UnloadDigits();
  loader->UnloadRecPoints();
}

//_____________________________________________________________________________
void AliZDCReconstructor::Reconstruct(AliRunLoader* runLoader, 
                                      AliRawReader* rawReader) const
{
  // *** Local ZDC reconstruction for raw data
  
  // Calibration data
  int runNumber = 0;
  AliZDCCalibData *calibda = GetCalibData(runNumber);  
  
  Float_t meanPed[47];
  for(Int_t jj=0; jj<47; jj++) meanPed[jj] = calibda->GetMeanPed(jj);

  AliLoader* loader = runLoader->GetLoader("ZDCLoader");
  if (!loader) return;
  loader->LoadRecPoints("recreate");
  // Event loop
  Int_t iEvent = 0;
  while (rawReader->NextEvent()) {
    runLoader->GetEvent(iEvent++);

    // loop over raw data digits
    Float_t zncorr=0, zpcorr=0, zemcorr=0;
    AliZDCRawStream digit(rawReader);
    while (digit.Next()) {
      if(digit.IsADCDataWord()){
        if(digit.GetADCGain() == 0){
          if(digit.GetSector(0) == 1)	   zncorr  += (Float_t) (digit.GetADCValue()-meanPed[digit.GetSector(1)]); // pedestals for high gain ZN ADCs;
          else if(digit.GetSector(0) == 2) zpcorr  += (Float_t) (digit.GetADCValue()-meanPed[digit.GetSector(1)+10]); // pedestals for high gain ZP ADCs;
          else if(digit.GetSector(0) == 3) zemcorr += (Float_t) (digit.GetADCValue()-meanPed[digit.GetSector(1)+20]); // pedestals for high gain ZEM ADCs;
	}
      }
    }
    if(zncorr<0)  zncorr=0;
    if(zpcorr<0)  zpcorr=0;
    if(zemcorr<0) zemcorr=0;
    
    // reconstruct the event
    printf("\n\t ZDCReco from raw-> Ev.#%d ZN = %.0f, ZP = %.0f, ZEM = %.0f\n",iEvent,zncorr,zpcorr,zemcorr);
    ReconstructEvent(loader, zncorr, zpcorr, zemcorr);
  }

  loader->UnloadRecPoints();
}

//_____________________________________________________________________________
void AliZDCReconstructor::ReconstructEvent(AliLoader* loader, Float_t zncorr,
	Float_t zpcorr, Float_t zemcorr) const
{
  // ***** Reconstruct one event
  
  //  ---      ADCchannel -> photoelectrons
  // NB-> PM gain = 10^(5), ADC resolution = 6.4*10^(-7)
  // Move to V965 (E.S.,15/09/04) NB-> PM gain = 10^(5), ADC resolution = 8*10^(-7)
  Float_t znphe, zpphe, zemphe, convFactor = 0.08;
  znphe  = zncorr/convFactor;
  zpphe  = zpcorr/convFactor;
  zemphe = zemcorr/convFactor;
  //if AliDebug(1,Form("\n    znphe = %f, zpphe = %f, zemphe = %f\n",znphe, zpphe, zemphe);
  
  //  ---      Energy calibration
  // Conversion factors for hadronic ZDCs goes from phe yield to TRUE 
  // incident energy (conversion from GeV to TeV is included); while for EM 
  // calos conversion is from light yield to detected energy calculated by
  // GEANT NB -> ZN and ZP conversion factors are constant since incident
  // spectators have all the same energy, ZEM energy is obtained through a
  // fit over the whole range of incident particle energies 
  // (obtained with full HIJING simulations) 
  Float_t znenergy, zpenergy, zemenergy, zdcenergy;
  Float_t znphexTeV=329., zpphexTeV=369.;
  znenergy  = znphe/znphexTeV;
  zpenergy  = zpphe/zpphexTeV;
  zdcenergy = znenergy+zpenergy;
  zemenergy = -4.81+0.3238*zemphe;
  if(zemenergy<0) zemenergy=0;
  //  if AliDebug(1,Form("    znenergy = %f TeV, zpenergy = %f TeV, zdcenergy = %f GeV, "
  //			   "\n		zemenergy = %f TeV\n", znenergy, zpenergy, 
  //			   zdcenergy, zemenergy);
  
  //  if(zdcenergy==0)
  //    if AliDebug(1,Form("\n\n	###	ATTENZIONE!!! -> ev# %d: znenergy = %f TeV, zpenergy = %f TeV, zdcenergy = %f GeV, "
  //			     " zemenergy = %f TeV\n\n", fMerger->EvNum(), znenergy, zpenergy, zdcenergy, zemenergy); 

  //  ---      Number of incident spectator nucleons
  Int_t nDetSpecN, nDetSpecP;
  nDetSpecN = (Int_t) (znenergy/2.760);
  nDetSpecP = (Int_t) (zpenergy/2.760);
//  if AliDebug(1,Form("\n    nDetSpecN = %d, nDetSpecP = %d\n",nDetSpecN, nDetSpecP);

  Int_t nGenSpecN=0, nGenSpecP=0, nGenSpec=0;
  Double_t impPar=0;
  // Cut value for Ezem (GeV)
  // [2] ### Results from a new production  -> 0<b<18 fm (Apr 2002)
  Float_t eZEMCut = 420.;
  Float_t deltaEZEMSup = 690.; 
  Float_t deltaEZEMInf = 270.; 
  if(zemenergy > (eZEMCut+deltaEZEMSup)){
    nGenSpecN = (Int_t) (fZNCen->Eval(znenergy));
    nGenSpecP = (Int_t) (fZPCen->Eval(zpenergy));
    nGenSpec  = (Int_t) (fZDCCen->Eval(zdcenergy));
    impPar    = fbCen->Eval(zdcenergy);
  }
  else if(zemenergy < (eZEMCut-deltaEZEMInf)){
    nGenSpecN = (Int_t) (fZNPer->Eval(znenergy)); 
    nGenSpecP = (Int_t) (fZPPer->Eval(zpenergy));
    nGenSpec  = (Int_t) (fZDCPer->Eval(zdcenergy));
    impPar    = fbPer->Eval(zdcenergy);
  }
  else if(zemenergy >= (eZEMCut-deltaEZEMInf) && zemenergy <= (eZEMCut+deltaEZEMSup)){
    nGenSpecN = (Int_t) (fZEMn->Eval(zemenergy));
    nGenSpecP = (Int_t) (fZEMp->Eval(zemenergy));
    nGenSpec  = (Int_t)(fZEMsp->Eval(zemenergy));
    impPar    =  fZEMb->Eval(zemenergy);
  }
  // [2] ### Results from a new production  -> 0<b<18 fm (Apr 2002)
  if(znenergy>162.)  nGenSpecN = (Int_t) (fZEMn->Eval(zemenergy));
  if(zpenergy>59.75)  nGenSpecP = (Int_t) (fZEMp->Eval(zemenergy));
  if(zdcenergy>221.5) nGenSpec  = (Int_t)(fZEMsp->Eval(zemenergy));
  if(zdcenergy>220.)  impPar    =  fZEMb->Eval(zemenergy);
  
  if(nGenSpecN>125)    nGenSpecN=125;
  else if(nGenSpecN<0) nGenSpecN=0;
  if(nGenSpecP>82)     nGenSpecP=82;
  else if(nGenSpecP<0) nGenSpecP=0;
  if(nGenSpec>207)     nGenSpec=207;
  else if(nGenSpec<0)  nGenSpec=0;
  
  //  ---      Number of participants
  Int_t nPart, nPartTot;
  nPart = 207-nGenSpecN-nGenSpecP;
  nPartTot = 207-nGenSpec;
  printf("\t  ZDCeventReco-> ZNEn = %.0f GeV, ZPEn = %.0f GeV, ZEMEn = %.0f GeV\n",
  	znenergy, zpenergy, zemenergy);

  // create the output tree
  loader->MakeTree("R");
  TTree* treeR = loader->TreeR();
  AliZDCReco reco(znenergy, zpenergy, zdcenergy, zemenergy,
                  nDetSpecN, nDetSpecP, nGenSpecN, nGenSpecP, nGenSpec,
                  nPartTot, impPar);
  AliZDCReco* preco = &reco;
  const Int_t kBufferSize = 4000;
  treeR->Branch("ZDC", "AliZDCReco", &preco, kBufferSize);

  // write the output tree
  treeR->Fill();
  loader->WriteRecPoints("OVERWRITE");
}

//_____________________________________________________________________________
void AliZDCReconstructor::FillESD(AliRunLoader* runLoader, 
				  AliESD* esd) const
{
// fill energies and number of participants to the ESD

  AliLoader* loader = runLoader->GetLoader("ZDCLoader");
  if (!loader) return;
  loader->LoadRecPoints();

  TTree* treeR = loader->TreeR();
  if (!treeR) return;
  AliZDCReco reco;
  AliZDCReco* preco = &reco;
  treeR->SetBranchAddress("ZDC", &preco);

  treeR->GetEntry(0);
  esd->SetZDC(reco.GetZNenergy(), reco.GetZPenergy(), reco.GetZEMenergy(),
	      reco.GetNPart());

  loader->UnloadRecPoints();
}

//_____________________________________________________________________________
AliZDCCalibData* AliZDCReconstructor::GetCalibData(int runNumber) const
{

  //printf("\n\t AliZDCReconstructor::GetCalibData \n");
  //fStorage->PrintSelectionList();
  //AliCDBEntry *entry = fStorage->Get("DBlocal/ZDC/Calib/Data",runNumber);

  
  AliCDBStorage *fStorage = AliCDBManager::Instance()->GetStorage("local://$ALICE_ROOT");
  AliCDBEntry  *entry = fStorage->Get("ZDC/Calib/Data",0);
  
  AliZDCCalibData *calibda = (AliZDCCalibData*) entry->GetObject();

  //AliCDBManager::Instance()->Destroy();

  return calibda;

}
