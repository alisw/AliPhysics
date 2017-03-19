#include "TNamed.h"
#include "TRandom.h"

#include "AliPHOSTriggerUtils.h"
#include "AliVEvent.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliVCluster.h"
#include "AliVCaloCells.h"
#include "AliOADBContainer.h"
#include "AliPHOSGeoUtils.h"
#include "AliLog.h"

ClassImp(AliPHOSTriggerUtils)
//-----------------------------------------------------------------------    
AliPHOSTriggerUtils::AliPHOSTriggerUtils():TNamed(),
fbdrL(-1),
fbdrR( 3),
fRun(-1), 
fEvent(0x0),
fGeom(0x0),
fNtrg4x4(0)
{
    
   for(Int_t mod=0;mod<5;mod++)fPHOSBadMap[mod]=0x0 ;
    
}
//-----------------------------------------------------------------------    
AliPHOSTriggerUtils::AliPHOSTriggerUtils(const Text_t* name, const Text_t* title):
TNamed(name, title),
fbdrL(-1),
fbdrR( 3),
fRun(-1), 
fEvent(0x0),
fGeom(new AliPHOSGeoUtils("IHEP","")),
fNtrg4x4(0)
{

  fGeom = new AliPHOSGeoUtils("IHEP","");

  for(Int_t mod=0;mod<5;mod++)fPHOSBadMap[mod]=0x0 ;


}
//-----------------------------------------------------------------------    
void AliPHOSTriggerUtils::SetEvent(AliVEvent * event){
   //read bad map corresponding this run, 
   //prepre list of fired PHOS triggers
  
  fEvent=event ;

   Int_t runNumber=-999 ;
   AliAODEvent *aod = dynamic_cast<AliAODEvent*>(event) ;
   if(aod)
      runNumber = aod->GetRunNumber() ;
   else{
      AliESDEvent *esd = dynamic_cast<AliESDEvent*>(event) ;
      if(esd)
         runNumber = esd->GetRunNumber() ;
   }
  //Check if run number changed and one should (re)read bad map
  if(fRun!=runNumber){
   fRun = runNumber ;
     // 
    AliOADBContainer badmapContainer(Form("phosTriggerBadMap"));
    badmapContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSTrigBadMaps.root","phosTriggerBadMap");
//    badmapContainer.InitFromFile("$PHOSTrigBadMaps.root","phosTrigBadMap");
    TObjArray *maps = (TObjArray*)badmapContainer.GetObject(runNumber,"phosTriggerBadMap");
    if(!maps){
      AliError(Form("Can not read Trigger Bad map for run %d. \n",runNumber)) ;    
    }
    else{
      AliInfo(Form("Setting PHOS Trigger bad map with name %s \n",maps->GetName())) ;
      for(Int_t mod=0; mod<5;mod++){
        if(fPHOSBadMap[mod]) 
          delete fPHOSBadMap[mod] ;
        TH2I * h = (TH2I*)maps->At(mod) ;      
	if(h)
          fPHOSBadMap[mod]=new TH2I(*h) ;
      }
    }    
  }
  
  
 AliVCaloTrigger* trgESD = event->GetCaloTrigger("PHOS");
 trgESD->Reset();

 // *************************************************************************************************
 // Trigger
 // ************************************************************************
  fNtrg4x4 = 0;
  //Loop over 4x4 fired regions
  while(trgESD->Next()) {
    Int_t tmod,tabsId; // "Online" module number, bottom-left 4x4 edge cell absId
    trgESD->GetPosition(tmod,tabsId);
    Int_t trelid[4] ;
    fGeom->AbsToRelNumbering(tabsId,trelid);     
     
    Int_t l1Threshold=trgESD->GetL1TimeSum() ; // tdig->GetL1Threshold()  // L1 threshold: 0-High, 1-Medium, 2-Low
    Int_t trType=0; // 0-L0, 1-L1
    trgESD->GetTriggerBits(trType) ;
    Int_t trBits=0;
    if(trType==0) trBits=1 ; 
    else trBits=1<<(3-l1Threshold) ;  //Bits: 0: L0,  1:L1low, 2:L1med, 3:L1high 
  
    
    //Check if already exist
    Bool_t isNew=kTRUE ;
    for(Int_t i=0;i<fNtrg4x4; i++){
      if(fTrMod4x4[i] != trelid[0]) continue;
      if(fTrX4x4[i] != trelid[2]) continue;
      if(fTrZ4x4[i] != trelid[3])continue;
      //exist  
      fTrK4x4[i]|= trBits ;
      isNew=kFALSE ;
      break ;
    }
    if(isNew){
      fTrMod4x4[fNtrg4x4] = trelid[0];
      fTrX4x4[fNtrg4x4]  = trelid[2];
      fTrZ4x4[fNtrg4x4]  = trelid[3];
      fTrK4x4[fNtrg4x4] = trBits; // L0, L1low, L1med, L1high
      if ( fNtrg4x4<999 ) fNtrg4x4++;
      else AliError("PHOS trigger capacity exceeded, increase size of fTrMod4x4[],fTrX4x4[],fTrZ4x4[] arrays") ;
    }
    
  }  
}

//-----------------------------------------------------------------------    
Int_t AliPHOSTriggerUtils::IsFiredTrigger(AliVCluster * clu){ 
//Returns which trigger fired cluster:
//Bits: 0: L0,  1:L1low, 2:L1med, 3:L1high     
  
   if(!fEvent)
     return 0 ;
  
   //Maximum energy tower
   Int_t maxId, relid[4];
   Double_t eMax = -111;
     
   AliVCaloCells * phsCells=fEvent->GetPHOSCells() ;
   
   for (Int_t iDig=0; iDig< clu->GetNCells(); iDig++){
      Int_t cellAbsId1 = clu->GetCellAbsId(iDig);
      Double_t eCell = phsCells->GetCellAmplitude(cellAbsId1)*clu->GetCellAmplitudeFraction(iDig);
      if(eCell>eMax) 
        {
          eMax = eCell;
          maxId = cellAbsId1;
        }
    }
    
    fGeom->AbsToRelNumbering(maxId, relid);
    Int_t mod= relid[0] ;
    Int_t ix = relid[2];
    Int_t iz = relid[3];
        
    if(!TestBadMap(mod,ix,iz))
      return 0 ;
        
    for (Int_t itr=0; itr < fNtrg4x4; itr++ ){
      if(fTrMod4x4[itr] != mod ) continue; //trigger fired in other module
      if( (ix - fTrX4x4[itr] ) > fbdrL && (ix - fTrX4x4[itr] ) < fbdrR &&
          (iz - fTrZ4x4[itr] ) > fbdrL && (iz - fTrZ4x4[itr] ) < fbdrR  ){ //this is trigger
	return fTrK4x4[itr] ; 
//	    Int_t branch = FindBranch(tnX4x4[n4x4], tnZ4x4[n4x4]); //tru, gde srabotal trigger
      }
    }//tr44
    
    return 0 ;
}
//-----------------------------------------------------------------------  
Int_t AliPHOSTriggerUtils::IsFiredTriggerMC(AliVCluster * clu){
   //Estimate probability to fire PHOS trigger 
   //Bits: 0: L0,  1:L1low, 2:L1med, 3:L1high     
  
   if(!fEvent)
     return 0 ;
  
   //Maximum energy tower
   Int_t maxId, relid[4];
   Double_t eMax = -111;
     
   AliVCaloCells * phsCells=fEvent->GetPHOSCells() ;
   
   for (Int_t iDig=0; iDig< clu->GetNCells(); iDig++){
      Int_t cellAbsId1 = clu->GetCellAbsId(iDig);
      Double_t eCell = phsCells->GetCellAmplitude(cellAbsId1)*clu->GetCellAmplitudeFraction(iDig);
      if(eCell>eMax) 
        {
          eMax = eCell;
          maxId = cellAbsId1;
        }
    }
    
    fGeom->AbsToRelNumbering(maxId, relid);
    Int_t mod= relid[0] ;
    Int_t ix = relid[2];
    Int_t iz = relid[3];
        
    if(!TestBadMap(mod,ix,iz)){
      return 0 ;
    }
  
    //Now etimate probability from the parameterization
    
    Int_t result = 0;
    if(fRun<209122){ //Run1 Only L0
      if(gRandom->Uniform()<TriggerProbabilityLHC13bcdef(clu->E(),mod)) result=1 ;
    }
    else{
      for(Int_t bit=0; bit<4; bit++){
        if(gRandom->Uniform()<TriggerProbability(clu->E(),mod,bit))
            result|=1<<bit ;
      }
    }
    return result ;

}
//-----------------------------------------------------------------------  
Int_t AliPHOSTriggerUtils::FindBranch(Int_t nX, Int_t nZ) {
  //Calculate branch number from column and raw of a cell 
  
  Int_t sm=-1;

  if ( nX>=0  && nX <=15 && nZ <= 27 ) sm = 0;
  if ( nX>=16 && nX <=31 && nZ <= 27 ) sm = 1;        
  if ( nX>=32 && nX <=47 && nZ <= 27 ) sm = 2;
  if ( nX>=48 && nX <=64 && nZ <= 27 ) sm = 3;             

  if ( nX>=0  && nX <=15 && nZ >= 28 && nZ <= 56 ) sm = 4;
  if ( nX>=16 && nX <=31 && nZ >= 28 && nZ <= 56 ) sm = 5;
  if ( nX>=32 && nX <=47 && nZ >= 28 && nZ <= 56 ) sm = 6;
  if ( nX>=48 && nX <=64 && nZ >= 28 && nZ <= 56 ) sm = 7;

  return sm;
}
//-----------------------------------------------------------------------  
Bool_t AliPHOSTriggerUtils::TestBadMap(Int_t mod, Int_t ix,Int_t iz){
  //Test if this is good cell
  
  if(!fPHOSBadMap[mod])
    return kTRUE ;
  if( fPHOSBadMap[mod]->GetBinContent(ix,iz)>0)
    return kFALSE ;
  else 
    return kTRUE ;
  
  
}
//-----------------------------------------------------------------------  
Double_t AliPHOSTriggerUtils::TriggerProbability(Double_t eClu, Int_t module,Int_t bit){
 //Placeholder for turn-on curves parametrization
    
 return 1;    
}
//-----------------------------------------------------------------------  
Double_t AliPHOSTriggerUtils::TriggerProbabilityLHC13bcdef(Double_t eClu, Int_t module){
  //Parameterization provided by Yu. and V. Riabov
  
  if (module!=1 && module!=3) return 0;

  //Parameterization of efficiency function
  const Double_t ar[4][7]={0.,0.,0.,0.,0.,0.,0.,
   7.10000e-001,3.75156e+007,1.76476e+007,3.77306e+007,1.61226e+006, 3.76983e+007,9.02358e+000,
   0.,0.,0.,0.,0.,0.,0.,
   1.00000e+000, 1.15962e+008,1.04148e+008,9.06092e+000,1.37488e+008,2.25460e+011,1.28279e+001} ;
   
  //tabulated depence for low E part 
  const Int_t num[4]={0,22,0,24} ;
  const Double_t ent[4]={0.,6.83333,0.,7.5} ;
  const Double_t en[4][30]={
  0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 
  0.0000, 0.1667, 0.5000, 0.8333, 1.1667, 1.5000, 1.8333, 2.1667, 2.5000, 2.8333, 3.1667, 3.5000, 3.8333, 4.1667, 4.5000, 4.8333, 5.1667, 5.5000, 5.8333, 6.1667, 6.5000, 6.8333, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 
  0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 
  0.0000, 0.1667, 0.5000, 0.8333, 1.1667, 1.5000, 1.8333, 2.1667, 2.5000, 2.8333, 3.1667, 3.5000, 3.8333, 4.1667, 4.5000, 4.8333, 5.1667, 5.5000, 5.8333, 6.1667, 6.5000, 6.8333, 7.1667, 7.5000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000 } ;

  const Double_t eff[4][30]={
  0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 
  0.0000, 0.0001, 0.0002, 0.0006, 0.0020, 0.0056, 0.0112, 0.0201, 0.0314, 0.0442, 0.0580, 0.0748, 0.0959, 0.1154, 0.1371, 0.1498, 0.1705, 0.1844, 0.2194, 0.2465, 0.2798, 0.3366, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 
  0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 
  0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0001, 0.0001, 0.0003, 0.0007, 0.0023, 0.0086, 0.0229, 0.0404, 0.0561, 0.0701, 0.0734, 0.0831, 0.1090, 0.1279, 0.1529, 0.1950, 0.2313, 0.3302, 0.4269, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000 } ; 


  if (eClu >= ent[module]){
    Double_t denom= 1.+ ar[module][1]*TMath::Power(eClu,-ar[module][2]) + ar[module][3]*TMath::Power(eClu,-ar[module][4]) + ar[module][5]*TMath::Power(eClu,-ar[module][6]) ;
    if(denom>0)
      return ar[module][0]/denom ;
    else
      return 0. ;
  }
  else {
    for (Int_t i=1; i<num[module]; i++) {
      if (eClu <= en[module][i]) {
	return eff[module][i-1]+(eClu - en[module][i-1])*(eff[module][i]-eff[module][i-1])/(en[module][i]-en[module][i-1]);
      }
    }
  }
  return 0 ;  
  
}