#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliCDBManager.h"
#include "AliITSCalibrationSDD.h"
#include "AliITSgeomTGeo.h"
#include "AliCDBMetaData.h"
#include "AliCDBStorage.h"
#include "AliCDBId.h"
#include "AliCDBPath.h"
#include "AliCDBEntry.h"
#include <TObjArray.h>
#include <TRandom3.h>
#endif

void StoreCalibSDD(Int_t firstRun=0,Int_t lastRun=AliCDBRunRange::Infinity(), Bool_t isIdeal=kFALSE){
  ///////////////////////////////////////////////////////////////////////
  // Macro to generate and store the calibration files for SDD         //
  // Generates:                                                        //
  //  1 file with 260 AliITSCalibrationSDD objects with                //
  //    baselines, noise, gain, drift speed for each module (CalibSDD) //
  ///////////////////////////////////////////////////////////////////////
  
  if(!AliCDBManager::Instance()->IsDefaultStorageSet()) {
    AliCDBManager::Instance()->SetDefaultStorage("local://OCDB");
  }
  

  AliCDBMetaData *md= new AliCDBMetaData(); // metaData describing the object
  md->SetObjectClassName("TObjArray");
  md->SetResponsible("Francesco Prino");
  md->SetBeamPeriod(0);
  md->SetComment("Simulated data");

  AliCDBId idCalSDD("ITS/Calib/CalibSDD",firstRun, lastRun);
  TObjArray calSDD(260);
  calSDD.SetOwner(kFALSE);

  // BAD modules data 
  const Int_t nbadmod=6;
  Int_t idBadMod[nbadmod];
  idBadMod[0]=AliITSgeomTGeo::GetModuleIndex(3,4,1);
  idBadMod[1]=AliITSgeomTGeo::GetModuleIndex(3,4,2);
  idBadMod[2]=AliITSgeomTGeo::GetModuleIndex(3,4,3);
  idBadMod[3]=AliITSgeomTGeo::GetModuleIndex(4,1,1);
  idBadMod[4]=AliITSgeomTGeo::GetModuleIndex(4,12,1);
  idBadMod[5]=AliITSgeomTGeo::GetModuleIndex(4,21,3);
  
  // Modules with bad left side 
  const Int_t nbadleft=5;
  Int_t idBadLeft[nbadleft];
  idBadLeft[0]=AliITSgeomTGeo::GetModuleIndex(3,13,6);
  idBadLeft[1]=AliITSgeomTGeo::GetModuleIndex(4,4,8);
  idBadLeft[2]=AliITSgeomTGeo::GetModuleIndex(4,5,3);
  idBadLeft[3]=AliITSgeomTGeo::GetModuleIndex(4,18,4);
  idBadLeft[4]=AliITSgeomTGeo::GetModuleIndex(4,21,5);


 // BAD anodes data
  const Int_t nData = 209;
  Int_t anodeUp[209] = {0,36,0,12,20,32,0,0,12,76,28,8,16,0,0,0,8,0,0,0,20,4,0,0,0,0,0,0
			,0,0,8,0,0,0,0,0,0,0,0,0,0,0,12,0,8,0,4,4,0,160,0,0,0,252,16,0,8,8
			,0,0,12,0,0,0,12,0,15,20,0,0,0,0,0,0,0,12,12,0,0,0,0,0,8,8,3,240,212
			,12,0,8,0,12,0,0,12,0,8,24,0,0,12,0,0,0,0,40,0,0,40,12,28,0,0,12,12
			,0,0,0,0,20,0,0,0,12,0,24,0,0,0,0,0,0,0,8,16,0,0,0,0,0,256,0,0,0,0,0,20
			,0,12,0,0,0,0,24,0,0,0,0,0,0,0,20,0,0,16,0,0,0,0,24,0,0,0,8,0,16,40,0
			,0,0,0,0,0,0,0,0,4,0,32,12,8,28,0,76,0,0,0,12,60,144,0,0,0,0,16,0,16,0,3 };

  Int_t anodeDown[209] = {0,8,0,48,0,12,0,8,0,80,12,0,4,4,0,0,24,0,0,20,0,0,0,0,20,0,0,0,0,0,0,0
			  ,0,0,0,0,4,0,24,4,0,0,8,0,0,36,20,0,12,236,0,0,12,256,16,8,32,4,0,0,24
			  ,24,10,0,16,8,0,2,40,0,0,0,24,0,0,0,8,0,0,0,0,0,32,240,0,92,60,28,8,0,0
			  ,2,0,0,0,0,12,48,0,0,0,0,0,36,11,12,0,0,0,12,12,11,0,20,0,0,12,20,0,0,4
			  ,0,8,12,0,0,0,16,16,0,32,72,12,0,88,20,16,112,8,0,244,28,256,28,0,24,236
			  ,56,0,68,0,4,20,208,20,12,4,28,12,0,0,20,12,0,100,0,16,8,8,0,24,16,0,12,12
			  ,16,0,16,20,0,28,0,8,24,0,12,8,4,40,0,104,96,32,140,20,12,8,20,24,16,16,20
			  ,8,140,96,0,32,20,44};

  TRandom3 *gran = new TRandom3();
  
  for(Int_t mod=0;mod<260;mod++){
    AliITSCalibrationSDD* resd = new AliITSCalibrationSDD("simulated");
    resd->SetAMAt20MHz();
    Int_t nBadUp = 0;
    Int_t nBadDown = 0;
      

    // gain
    for(Int_t iWing=0; iWing<2;iWing++){
      for(Int_t iChip=0; iChip<4;iChip++){
	Float_t chipgain=gran->Gaus(1.,0.1);
	if(chipgain<0.1) chipgain=0.1;
	for(Int_t iChan=0; iChan<64;iChan++){
	  Float_t gain=gran->Gaus(chipgain,0.01);
	  if(gain<0.1) gain=0.1;
	  Int_t ian=resd->GetAnodeNumber(iWing,iChip,iChan);
	  resd->SetGain(ian,gain);
	}
      }
    }

    // bad channels
    if(!isIdeal){
      Int_t count=0;
      do {
	Int_t ranMod = Int_t(nData*gran->Uniform());
	nBadUp   = anodeUp[ranMod];
	nBadDown = anodeDown[ranMod];
      }while (nBadUp+nBadDown>25);

      resd->SetDeadChannels(nBadDown+nBadUp);
      Int_t remainingBad = nBadDown;
      while (remainingBad>0) {
	Int_t nBadChannel;
	if (remainingBad<5) {
	  nBadChannel = remainingBad;
	} else {
	  Int_t width = remainingBad-(5-1);
	  if (width>4) width = 4;
	  nBadChannel = 5 + Int_t(width*gran->Uniform());
	}
	

	Int_t iChannelPos = Int_t( (4*64-nBadChannel)*gran->Uniform() );
	//      Int_t iChip = iChannelPos/64;
	//      Int_t iChan = iChannelPos - iChip*64;
	if(resd->IsBadChannel(iChannelPos)) continue;	
	Int_t *clus = new Int_t[nBadChannel];
	Int_t ich = iChannelPos;
	for(Int_t i=0;i<nBadChannel;i++){
	  clus[i]=ich;
	  ich++;
	}
	
	for(Int_t i=0;i<nBadChannel;i++){
	  if(resd->IsBadChannel(clus[i])) break;
	  resd->SetBadChannel(count,clus[i]);
	  count++;	  
	}
	remainingBad -= nBadChannel;
	delete [] clus;
      }

      // May happen that, due to overlapping clusters, we
      // have less bad channels than requested
      // Let's put the remaining one per one ...
      Int_t nSeToBad = 0;
      for (Int_t i=0; i<4; i++){
	for(Int_t j=0;j<64;j++){
	  Int_t ian=resd->GetAnodeNumber(0,i,j);
	  if (resd->GetChannelGain(ian)<0.0001) nSeToBad++;
	}
      }
      while (nSeToBad<nBadDown) {
	Int_t i = Int_t(4*64*gran->Uniform());
	if(resd->IsBadChannel(i)==kFALSE){
	  resd->SetBadChannel(count,i);
	  count++;
	  nSeToBad++;
	}
      }
      
      remainingBad = nBadUp;
      while (remainingBad>0) {
	Int_t nBadChannel;
	if (remainingBad<5) {
	  nBadChannel = remainingBad;
	} else {
	  Int_t width = remainingBad-(5-1);
	  if (width>4) width = 4;
	  nBadChannel = 5 + Int_t(width*gran->Uniform());
	}
	
	Int_t iChannelPos = Int_t( (4*64-nBadChannel)*gran->Uniform() );
	//      Int_t iChip = iChannelPos/64;
	//      Int_t iChan = iChannelPos - iChip*64;
	if(resd->IsBadChannel(iChannelPos)) continue;
	
	Int_t *clus = new Int_t[nBadChannel];
	Int_t ich = iChannelPos;
	for(Int_t i=0;i<nBadChannel;i++){
	  clus[i]=ich+256;
	  ich++;
	}
	for(Int_t i=0;i<nBadChannel;i++){
	
	  if(resd->IsBadChannel(clus[i])) break;  
	  resd->SetBadChannel(count,clus[i]);
	  count++;
	  
	}
	remainingBad -= nBadChannel;
	delete [] clus;
      }
      
      nSeToBad = 0;
      for (Int_t i=0; i<4; i++){
	for(Int_t j=0;j<64;j++){
	  Int_t ian=resd->GetAnodeNumber(1,i,j);
	  if (resd->GetChannelGain(ian)<0.0001) nSeToBad++;
	}
      }

      while (nSeToBad<nBadUp) {
	Int_t i = Int_t(4*64*gran->Uniform());
	if(resd->IsBadChannel(i+256)==kFALSE){
	  resd->SetBadChannel(count,i+256);
	  count++;
	  nSeToBad++;
	}
      }
    }
    //baselines
    /*
      for(Int_t nan=0;nan<512;nan++){
	
      Int_t baseline = resd->GetBaseline(0);
      Double_t noise = resd->GetNoiseAfterElectronics(0);
      resd->SetBaseline(nan,baseline+(Int_t)gran->Gaus(0,4)); //baseline di def. e' 20, ho messo variazione di 4
      resd->SetNoiseAfterElectronics(nan,noise+(Double_t)gran->Gaus(0,0.6));		  
			  
      }
    */

    Int_t modId=mod+240;
    // Bad modules
    if(!isIdeal){
      for(Int_t ibadm=0; ibadm<nbadmod; ibadm++){
	if(modId==idBadMod[ibadm]) resd->SetBad();
      }
      // Modules with bad left side
      for(Int_t ibadl=0; ibadl<nbadleft; ibadl++){
	if(modId==idBadLeft[ibadl]) for(Int_t ichip=0;ichip<4;ichip++) resd->SetChipBad(ichip);
      }

      // Modules with bad pascal chips
      if( modId==AliITSgeomTGeo::GetModuleIndex(4,4,2) ){
	resd->SetChipBad(0);
	resd->SetChipBad(3);
      }
    }
    calSDD.Add(resd);
    printf("Added module %d\n",mod);
  }
    
  FILE* out = fopen("deadchannels.dat","w");
  for(Int_t i=0;i<260;i++){
    fprintf(out,"MODULE=%d\n",i);
    AliITSCalibrationSDD* cl = (AliITSCalibrationSDD*)calSDD.At(i);
    Int_t ndead=cl->GetDeadChannels();
    fprintf(out,"n %d\n",ndead);
    for(Int_t n=0;n<ndead;n++){
      fprintf(out,"%d\n",cl->GetBadChannel(n));
    }   
  }
  fclose(out);
  AliCDBManager::Instance()->GetDefaultStorage()->Put(&calSDD, idCalSDD, md);   
}
