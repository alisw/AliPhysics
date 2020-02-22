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
fFixedRun(kFALSE),
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
fFixedRun(kFALSE),
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
  if(fFixedRun){
    Bool_t toRead= kTRUE ;  
    for(Int_t mod=0; mod<5;mod++){
      if(fPHOSBadMap[mod]) 
        toRead = kFALSE ; //already read
    }
    if(toRead){
      AliOADBContainer badmapContainer(Form("phosTriggerBadMap"));
      badmapContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSTrigBadMaps.root","phosTriggerBadMap");
      TObjArray *maps = (TObjArray*)badmapContainer.GetObject(fRun,"phosTriggerBadMap");
      if(!maps){
        AliError(Form("Can not read Trigger Bad map for run %d. \n",fRun)) ;    
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
  }
  else{ //use runnumber from header
    //Check if run number changed and one should (re)read bad map
    if(fRun!=runNumber){
      fRun = runNumber ;
      // 
      AliOADBContainer badmapContainer(Form("phosTriggerBadMap"));
      badmapContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSTrigBadMaps.root","phosTriggerBadMap");
      TObjArray *maps = (TObjArray*)badmapContainer.GetObject(fRun,"phosTriggerBadMap");
      if(!maps){
        AliError(Form("Can not read Trigger Bad map for run %d. \n",fRun)) ;    
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
   Int_t maxId=-1, relid[4];
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
   Int_t maxId=-1, relid[4];
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
      if(fRun>=265015 && fRun<=267166){ //LHC16qrst
        Int_t ddl = WhichDDL(mod, ix) ;  
        Double_t rdm =gRandom->Uniform() ; 
printf("E=%f, ddl=%d, Prob=%e \n",clu->E(), ddl, TriggerProbabilityLHC16qrst(clu->E(),ddl)) ;       
        if(rdm<TriggerProbabilityLHC16qrst(clu->E(),ddl)) result=1 ;         
      }
      else{  
        if(fRun>=282008 && fRun<=282441){ //LHC17pq
          Int_t ddl = WhichDDL(mod, ix) ;  
          Double_t rdm =gRandom->Uniform() ; 
          if(rdm<TriggerProbabilityLHC17pq(clu->E(),ddl)) result=1 ;         
        }
        else{
          for(Int_t bit=0; bit<4; bit++){
            if(gRandom->Uniform()<TriggerProbability(clu->E(),mod,bit))
               result|=1<<bit ;
          }
        }
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
//______________________________________________________
Int_t AliPHOSTriggerUtils::WhichDDL(Int_t module, Int_t cellx)
{
  const Int_t Nmod=5;//totally, 5 PHOS modules are designed.
  Int_t ddl = -1;

  if(cellx<1 || 64<cellx) return -1;

  if(module<1 || 4<module){
    return -1;
  }
  else{
    ddl = (Nmod-module) * 4 + (cellx-1)/16;//convert offline module numbering to online.
    return ddl;
  }
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
//-----------------------------------------------------------------------  
Double_t AliPHOSTriggerUtils::TriggerProbabilityLHC16qrst(Double_t x, Int_t ddl){
  //Parameterization of turn-on curve for period LHC16qrst
  // Note that it is done for DDLs, not modules

  if(ddl<8 || ddl>19) return 0.;
  
  if(fRun>=265015 && fRun<=267166){ //LHC16qrst
    switch(ddl){  
    case 8 : return 9.573766e-01/(TMath::Exp((3.871535e+00-x)/2.876473e-01)+1.)+(1.-9.573766e-01)/(TMath::Exp((7.153766e+00-x)/2.876473e-01)+1.) ;
    case 9 : return 9.228485e-01/(TMath::Exp((3.798753e+00-x)/2.889749e-01)+1.)+(1.-9.228485e-01)/(TMath::Exp((8.212820e+00-x)/2.889749e-01)+1.) ;
    case 10 : return 9.737340e-01/(TMath::Exp((3.684112e+00-x)/2.749400e-01)+1.)+(1.-9.737340e-01)/(TMath::Exp((6.193852e+00-x)/2.749400e-01)+1.) ;
    case 11 : return 9.141516e-01/(TMath::Exp((3.775433e+00-x)/2.522816e-01)+1.)+(1.-9.141516e-01)/(TMath::Exp((5.245370e+00-x)/2.522816e-01)+1.) ;
    case 12 : return 8.921613e-01/(TMath::Exp((3.881223e+00-x)/2.823943e-01)+1.)+(1.-8.921613e-01)/(TMath::Exp((5.141498e+00-x)/2.823943e-01)+1.) ;
    case 13 : return 8.513185e-01/(TMath::Exp((3.970404e+00-x)/2.723763e-01)+1.)+(1.-8.513185e-01)/(TMath::Exp((5.808491e+00-x)/2.723763e-01)+1.) ;
    case 14 : return 9.620960e-01/(TMath::Exp((4.163946e+00-x)/3.106500e-01)+1.)+(1.-9.620960e-01)/(TMath::Exp((8.264108e+00-x)/3.106500e-01)+1.) ;
    case 15 : return 8.836921e-01/(TMath::Exp((3.968618e+00-x)/2.916503e-01)+1.)+(1.-8.836921e-01)/(TMath::Exp((9.500000e+00-x)/2.916503e-01)+1.) ;
    case 16 : return 7.232389e-01/(TMath::Exp((3.867679e+00-x)/2.633759e-01)+1.)+(1.-7.232389e-01)/(TMath::Exp((5.822319e+00-x)/2.633759e-01)+1.) ;
    case 17 : return 8.811230e-01/(TMath::Exp((3.934961e+00-x)/2.731523e-01)+1.)+(1.-8.811230e-01)/(TMath::Exp((7.521001e+00-x)/2.731523e-01)+1.) ;
    case 18 : return 8.920217e-01/(TMath::Exp((3.741835e+00-x)/2.869366e-01)+1.)+(1.-8.920217e-01)/(TMath::Exp((9.500000e+00-x)/2.869366e-01)+1.) ;
    case 19 : return 7.797249e-01/(TMath::Exp((3.878770e+00-x)/2.757002e-01)+1.)+(1.-7.797249e-01)/(TMath::Exp((7.344906e+00-x)/2.757002e-01)+1.) ;
    default : return 0;
    }
  }
  return 0.;
}

//-----------------------------------------------------------------------  
Double_t AliPHOSTriggerUtils::TriggerProbabilityLHC17pq(Double_t x, Int_t ddl){
  //Parameterization of turn-on curve for period LHC17pq
  // Note that it is done for DDLs, not modules

  if(ddl<8 || ddl>19) return 0.;
    
  if(fRun>=282008 && fRun<=282343){ //LHC17p
    switch(ddl){  
    case 8 : return 7.944762e-01/(TMath::Exp((3.641718e+00-x)/2.650322e-01)+1.)+(1.-7.944762e-01)/(TMath::Exp((4.694861e+00-x)/2.650322e-01)+1.) ;
    case 9 : return 9.166827e-01/(TMath::Exp((3.517378e+00-x)/2.592677e-01)+1.)+(1.-9.166827e-01)/(TMath::Exp((6.836709e+00-x)/2.592677e-01)+1.) ;
    case 10 : return 9.530612e-01/(TMath::Exp((3.496633e+00-x)/2.656076e-01)+1.)+(1.-9.530612e-01)/(TMath::Exp((4.780714e+00-x)/2.656076e-01)+1.) ;
    case 11 : return 8.475332e-01/(TMath::Exp((3.672446e+00-x)/2.545452e-01)+1.)+(1.-8.475332e-01)/(TMath::Exp((4.681494e+00-x)/2.545452e-01)+1.) ;
    case 12 : return 1.009702e+00/(TMath::Exp((3.789344e+00-x)/3.089832e-01)+1.)+(1.-1.009702e+00)/(TMath::Exp((5.224435e+00-x)/3.089832e-01)+1.) ;
    case 13 : return 9.043261e-01/(TMath::Exp((3.691091e+00-x)/2.750534e-01)+1.)+(1.-9.043261e-01)/(TMath::Exp((6.038003e+00-x)/2.750534e-01)+1.) ;
    case 14 : return 9.368423e-01/(TMath::Exp((3.863201e+00-x)/3.307754e-01)+1.)+(1.-9.368423e-01)/(TMath::Exp((6.490584e+00-x)/3.307754e-01)+1.) ;
    case 15 : return 8.929566e-01/(TMath::Exp((3.650563e+00-x)/2.849078e-01)+1.)+(1.-8.929566e-01)/(TMath::Exp((6.933385e+00-x)/2.849078e-01)+1.) ;
    case 16 : return 8.097278e-01/(TMath::Exp((3.664689e+00-x)/2.879874e-01)+1.)+(1.-8.097278e-01)/(TMath::Exp((6.285455e+00-x)/2.879874e-01)+1.) ;
    case 17 : return 7.103594e-01/(TMath::Exp((3.528692e+00-x)/2.930813e-01)+1.)+(1.-7.103594e-01)/(TMath::Exp((5.220063e+00-x)/2.930813e-01)+1.) ;
    case 18 : return 9.025508e-01/(TMath::Exp((3.406599e+00-x)/3.195572e-01)+1.)+(1.-9.025508e-01)/(TMath::Exp((7.277194e+00-x)/3.195572e-01)+1.) ;
    case 19 : return 7.898787e-01/(TMath::Exp((3.521535e+00-x)/2.833550e-01)+1.)+(1.-7.898787e-01)/(TMath::Exp((6.468302e+00-x)/2.833550e-01)+1.) ;
    default: return 0.;
    }
  }
  if(fRun>=282365 && fRun<=282441){ //LHC17q
    switch(ddl){  
    case 8 : return 7.670229e-01/(TMath::Exp((3.488723e+00-x)/2.742928e-01)+1.)+(1.-7.670229e-01)/(TMath::Exp((5.059602e+00-x)/2.742928e-01)+1.) ;
    case 9 : return 8.151913e-01/(TMath::Exp((3.395833e+00-x)/2.143634e-01)+1.)+(1.-8.151913e-01)/(TMath::Exp((4.792992e+00-x)/2.143634e-01)+1.) ;
    case 10 : return 8.780434e-01/(TMath::Exp((3.448331e+00-x)/2.344344e-01)+1.)+(1.-8.780434e-01)/(TMath::Exp((4.636605e+00-x)/2.344344e-01)+1.) ;
    case 11 : return 7.220132e+00/(TMath::Exp((3.872164e+00-x)/3.007561e-01)+1.)+(1.-7.220132e+00)/(TMath::Exp((3.893820e+00-x)/3.007561e-01)+1.) ;
    case 12 : return 5.536853e-01/(TMath::Exp((3.382345e+00-x)/2.497804e-01)+1.)+(1.-5.536853e-01)/(TMath::Exp((4.190075e+00-x)/2.497804e-01)+1.) ;
    case 13 : return 5.150823e-01/(TMath::Exp((3.342826e+00-x)/1.628658e-01)+1.)+(1.-5.150823e-01)/(TMath::Exp((4.432898e+00-x)/1.628658e-01)+1.) ;
    case 14 : return 5.384429e-01/(TMath::Exp((3.534127e+00-x)/3.184165e-01)+1.)+(1.-5.384429e-01)/(TMath::Exp((4.565801e+00-x)/3.184165e-01)+1.) ;
    case 15 : return 7.354045e-01/(TMath::Exp((3.564357e+00-x)/3.220581e-01)+1.)+(1.-7.354045e-01)/(TMath::Exp((5.170461e+00-x)/3.220581e-01)+1.) ;
    case 16 : return 6.502367e-01/(TMath::Exp((3.582700e+00-x)/3.271896e-01)+1.)+(1.-6.502367e-01)/(TMath::Exp((5.037315e+00-x)/3.271896e-01)+1.) ;
    case 17 : return 1.164583e+00/(TMath::Exp((4.573204e+00-x)/6.050801e-01)+1.)+(1.-1.164583e+00)/(TMath::Exp((5.757209e+00-x)/6.050801e-01)+1.) ;
    case 18 : return 2.307368e-01/(TMath::Exp((5.113748e+00-x)/2.183349e-01)+1.)+(1.-2.307368e-01)/(TMath::Exp((3.180296e+00-x)/2.183349e-01)+1.) ;
    case 19 : return 4.047165e-01/(TMath::Exp((3.188518e+00-x)/1.792187e-01)+1.)+(1.-4.047165e-01)/(TMath::Exp((4.656638e+00-x)/1.792187e-01)+1.) ;
    default: return 0.;
    }
  }
  return 0.;
}
