// macro to create and store in the OADB the root file containing the refractive index
// values to be used for nSigmas and probabiltiy calculation in the PID framework
// HmpParams: the refractive index is "manually" assigned . One for each chamber, constant vs y.
// HmpParamsFromOCDB: the refracitve index is calulated starting from the parameters stored in OCDB. One for each radiator, not constant vs y. 
// Contact: giacomo.volpe@cern.ch

HmpParams(Int_t runFirst, Int_t runLast)
{
  Double_t sizePcY = 48.*0.84;
    
  TObjArray *arrayDef = new TObjArray(7);
  TObjArray *array    = new TObjArray(7);
  
  TF1 *f1[7], *fIdx[7];
  
  for(Int_t iCh=0; iCh<7; iCh++){
    
    f1[iCh] = new TF1(Form("f1_%i",iCh),"1.284",0,3.*sizePcY);
    arrayDef->AddAt(f1[iCh],iCh);
    
    fIdx[iCh] = new TF1(Form("fIdx_ch%i",iCh),"1.290",0,3.*sizePcY);
       
    array->AddAt(fIdx[iCh],iCh);
  }  
    
  AliOADBContainer* con = new AliOADBContainer("HMPoadb");

  AliHMPIDPIDParams *pPar = new AliHMPIDPIDParams("HMPparams");
  pPar->SetHMPIDrefIndex(array);
  Printf("Populating HMPID OADB for period **");
  con->AppendObject(pPar,runFirst,runLast);
  
  AliHMPIDPIDParams *pParDefault = new AliHMPIDPIDParams("HMPparams");
  pParDefault->SetHMPIDrefIndex(arrayDef);
//  Printf("Populating HMPID OADB with default entry",period.Data());
  Printf("Populating HMPID OADB with default entry");
  con->AddDefaultObject(pParDefault);
      
  con->WriteToFile("HMPIDPIDParams.root");
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
HmpParamsFromOCDB( Int_t year, Int_t runNumber, Int_t runFirst, Int_t runLast)
{
  Double_t sizePcY = 48.*0.84;
    
  TGrid::Connect("alien://");
  
  AliCDBManager *man = AliCDBManager::Instance();
  
  man->SetDefaultStorage(Form("alien://folder=/alice/data/%i/OCDB",year));
  
  man->SetRun(runNumber);
      
  AliCDBEntry *entry = man->Get("HMPID/Calib/Nmean");
    
  TObjArray *arr = (TObjArray*)entry->GetObject();
     
  TF1 *fEmean = (TF1*)arr->At(42);
  
  TObjArray *arrayDef = new TObjArray(21);
  TObjArray *array    = new TObjArray(21);
  
  TF1 *f1[21], *fIdx[21];
  
  for(Int_t iCh=0; iCh<7; iCh++){
    
    for(Int_t iRad=0; iRad<3; iRad++){
    
       f1[3*iCh+iRad] = new TF1(Form("f1_%i_%i",iCh,iRad),"1.290",iRad*sizePcY,(iRad+1)*sizePcY);
       arrayDef->AddAt(f1[3*iCh+iRad],3*iCh+iRad);
    
       TF1 *fIn    = (TF1*)arr->At(6*iCh+2*iRad);
       TF1 *fOut   = (TF1*)arr->At(6*iCh+2*iRad+1);
     
       Double_t xmin, xmax;
     
       fIn->GetRange(xmin,xmax);
     
       Double_t tempIn  = fIn->Eval(0.5*(xmax+xmin));
       Double_t tempOut = fOut->Eval(0.5*(xmax+xmin));
    
       Double_t emean = fEmean->Eval(0.5*(xmax+xmin));  
       
       Double_t idxIn  = AliHMPIDParam::NIdxRad(emean,tempIn);
       Double_t idxOut = AliHMPIDParam::NIdxRad(emean,tempOut);
        
       Printf("ch = %i, rad = %i,Tin = %f, Tout = %f, emean = %f, ----- Tin(xmin) = %f, Tin(xmax) = %f, Tout(xmin) = %f, Tout(xmax) = %f",iCh,iRad,tempIn,tempOut,emean,fIn->Eval(xmin),fIn->Eval(xmax),fOut->Eval(xmin),fOut->Eval(xmax));
     
       if(tempOut<tempIn) tempOut = tempIn;
       
       Double_t gradT = (tempOut - tempIn)/sizePcY;
       
       TF1 *fTemp = new TF1(Form("Temp_%i_%i",iCh,iRad),Form("%f*x + %f",gradT,tempIn),iRad*sizePcY,(iRad+1)*sizePcY);
       
      // fIdx[3*iCh+iRad] = new TF1(Form("fIdx_ch%i_rad%i",iCh,iRad),Form("sqrt(1+0.554*(1239.84/6.675)*(1239.84/6.675)/((1239.84/6.675)*(1239.84/6.675)-5769))-0.0005*(Temp_%i_%i-20)",iCh,iRad),iRad*sizePcY,(iRad+1)*sizePcY);
       fIdx[3*iCh+iRad] = new TF1(Form("fIdx_ch%i_rad%i",iCh,iRad),"1.290",iRad*sizePcY,(iRad+1)*sizePcY);
       
       array->AddAt(fIdx[3*iCh+iRad],3*iCh+iRad);
    }
    
  }  
    
  //TString period = "LHC10b";
  //Int_t startTime;

  AliOADBContainer* con = new AliOADBContainer("HMPoadb");

  AliHMPIDPIDParams *pPar = new AliHMPIDPIDParams("HMPparams");
  pPar->SetHMPIDrefIndex(array);
  Printf("Populating HMPID OADB for period **");
  con->AppendObject(pPar,runFirst,runLast);
  
  AliHMPIDPIDParams *pParDefault = new AliHMPIDPIDParams("HMPparams");
  pParDefault->SetHMPIDrefIndex(arrayDef);
//  Printf("Populating HMPID OADB with default entry",period.Data());
  Printf("Populating HMPID OADB with default entry");
  con->AddDefaultObject(pParDefault);
      
  con->WriteToFile("HMPIDPIDParams.root");
}
