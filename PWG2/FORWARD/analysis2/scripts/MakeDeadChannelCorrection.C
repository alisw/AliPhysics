//_____________________________________________________________________
Float_t GetMaxR(Char_t ring) const {
  //Get max R of ring
  Float_t radius = 0;
  if(ring == 'I')
    radius = 17.2;
  else if(ring == 'O')
    radius = 28.0;
  else
    AliWarning("Unknown ring - must be I or O!");
  
  return radius;
}
//_____________________________________________________________________
Float_t GetMinR(Char_t ring) const{
  //Get min R of ring
  Float_t radius = 0;
  if(ring == 'I')
    radius = 4.5213;
  else if(ring == 'O')
    radius = 15.4;
  else
    AliWarning("Unknown ring - must be I or O!");
  
  return radius;

}
//_____________________________________________________________________
Float_t GetEtaFromStrip(UShort_t det, Char_t ring, UShort_t sec, UShort_t strip, Float_t zvtx)
{
  //Calculate eta from strip with vertex (redundant with AliESDFMD::Eta)
  Float_t   rad       = GetMaxR(ring)-GetMinR(ring);
  Float_t   nStrips   = (ring == 'I' ? 512 : 256);
  Float_t   segment   = rad / nStrips;
  Float_t   r         = GetMinR(ring) + segment*strip;
  Float_t   z         = 0;
  Int_t hybrid = sec / 2;
  
  if(det == 1) {
    if(!(hybrid%2)) z = 320.266; else z = 319.766;
  }
  if(det == 2 && ring == 'I' ) {
    if(!(hybrid%2)) z = 83.666; else z = 83.166;
  }
  if(det == 2 && ring == 'O' ) {
    if(!(hybrid%2)) z = 74.966; else z = 75.466;
  }
  if(det == 3 && ring == 'I' ) {
    if(!(hybrid%2)) z = -63.066; else z = -62.566;
  }
  if(det == 3 && ring == 'O' ) {
    if(!(hybrid%2)) z = -74.966; else z = -75.466;
  }
  
  Float_t   theta = TMath::ATan2(r,z-zvtx);
  Float_t   eta   = -1*TMath::Log(TMath::Tan(0.5*theta));
  
  return eta;
}
//_____________________________________________________________________
Float_t GetPhiFromSector(UShort_t det, Char_t ring, UShort_t sec) 
{
  //Get phi from sector
  Int_t nsec = (ring == 'I' ? 20 : 40);
  Float_t basephi = 0;
  if(det == 1) 
    basephi = 1.72787594; 
  if(det == 2 && ring == 'I')
    basephi = 0.15707963;
  if(det == 2 && ring == 'O')
    basephi = 0.078539818;
  if(det == 3 && ring == 'I')
    basephi = 2.984513044;
  if(det == 3 && ring == 'O')
    basephi = 3.06305289;
  
  Float_t step = 2*TMath::Pi() / nsec;
  Float_t phi = 0;
  if(det == 3)
    phi = basephi - sec*step;
  else
    phi = basephi + sec*step;
  
  if(phi < 0) 
    phi = phi +2*TMath::Pi();
  if(phi > 2*TMath::Pi() )
    phi = phi - 2*TMath::Pi();
  
  return phi;
}
//_____________________________________________________________________
void MakeDeadChannelCorrection(Int_t runnumber, Int_t nVtxBins=10, Float_t vtxLow=-10, Float_t vtxHigh=10){
  
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libPWG2forward2");
  
  Float_t delta = (vtxHigh - vtxLow) / (Float_t)nVtxBins;
  
  //TGrid::Connect("alien://",0,0,"t");
  AliCDBManager* cdb = AliCDBManager::Instance();
  //cdb->SetDefaultStorage("alien://Folder=/alice/data/2010/OCDB");
  cdb->SetDefaultStorage("local:///home/canute/ALICE/AliRoot/OCDB");
  cdb->SetRun(runnumber);
  
  TObjArray* fReadArray = new TObjArray();
  fReadArray->SetName("DeadChannels");
  TObjArray* fAllArray = new TObjArray();
  fAllArray->SetName("AllChannels");
  
  TH2D* hRead = 0; 
  TH2D* hAll = 0; 
  for(Int_t det =1; det<=3;det++)
    {
      TObjArray* detReadArray = new TObjArray();
      detReadArray->SetName(Form("FMD%d_Read",det));
      fReadArray->AddAtAndExpand(detReadArray,det);
      
      TObjArray* detArray = new TObjArray();
      detArray->SetName(Form("FMD%d_All",det));
      fAllArray->AddAtAndExpand(detArray,det);
      
      UShort_t nRings = (det==1 ? 1 : 2);
      for(Int_t ring = 0;ring<nRings;ring++)
	{
	  TObjArray* vtxArray = new TObjArray();
	  detArray->AddAtAndExpand(vtxArray,ring);
	  
	  TObjArray* vtxReadArray = new TObjArray();
	  detReadArray->AddAtAndExpand(vtxReadArray,ring);
	  
	  
	  Char_t ringChar = (ring == 0 ? 'I' : 'O');
	  Int_t nSec = (ringChar == 'I' ? 20 : 40);
	  
	  for(Int_t v=0; v<nVtxBins; v++) {
	    hRead = new TH2D(Form("Read_FMD%d%c_vtxbin%d",det,ringChar,v),Form("Read_FMD%d%c_vtxbin%d",det,ringChar,v),200,-4,6,nSec,0,2*TMath::Pi());
	    hAll  = new TH2D(Form("All_FMD%d%c_vtxbin%d",det,ringChar,v),Form("All_FMD%d%c_vtxbin%d",det,ringChar,v),200,-4,6,nSec,0,2*TMath::Pi());
	    hRead->Sumw2();
	    hAll->Sumw2();
	    vtxArray->AddAtAndExpand(hAll,v);
	    vtxReadArray->AddAtAndExpand(hRead,v);
	    
	  }
	  
	  
	}
    }
  
  
  
  
  AliFMDParameters* pars = AliFMDParameters::Instance();
  pars->Init();
  Int_t nDead = 0;
  
  for(UShort_t d=1;d<=1;d++) {
    UShort_t nRings = (d==1 ? 1 : 2);
    for (UShort_t ir = 0; ir < nRings; ir++) {
      
      Char_t   ringChar = (ir == 0 ? 'I' : 'O');
      UShort_t nsec     = (ir == 0 ? 20  : 40);
      UShort_t nstr     = (ir == 0 ? 512 : 256);
      
      std::cout<<Form("Now in FMD%d%c", d, ringChar)<<std::endl;
      for(Int_t v = 0; v<nVtxBins ; v++) { 
	    
	    TObjArray* detReadArray = (TObjArray*)fReadArray->At(d);
	    TObjArray* vtxReadArray = (TObjArray*)detReadArray->At(ir);
	    TH2D*      hRead        = (TH2D*)vtxReadArray->At(v);
	    
	    TObjArray* detArray = (TObjArray*)fAllArray->At(d);
	    TObjArray* vtxArray = (TObjArray*)detArray->At(ir);
	    TH2D*      hAll        = (TH2D*)vtxArray->At(v);
	    
	    Float_t vtxMean = (v+0.5)*delta-vtxHigh;
	    for(UShort_t sec =0; sec < nsec;  sec++) {
	
	      for(UShort_t strip = 0; strip < nstr; strip++) {
	  
	  
		if(pars->IsDead(d,ringChar,sec,strip)) {
		  nDead++;
		}
		else hRead->Fill(GetEtaFromStrip(d,ringChar,sec,strip,vtxMean),GetPhiFromSector(d,ringChar,sec));
	    
		hAll->Fill(GetEtaFromStrip(d,ringChar,sec,strip,vtxMean),GetPhiFromSector(d,ringChar,sec));
	    
	    
	  }
	}
      }
    }
  }
  
  Float_t reldead = (Float_t)nDead / 51200.;
  
  std::cout<<Form("Found %d dead channels or %f percent",nDead,100*reldead)<<std::endl;
  
  AliFMDCorrDeadChannels* deadObject = new AliFMDCorrDeadChannels();
  deadObject->SetVertexAxis(nVtxBins, vtxLow, vtxHigh);
  
  for(UShort_t d=1;d<=3;d++) {
    UShort_t nRings = (d==1 ? 1 : 2);
    for (UShort_t ir = 0; ir < nRings; ir++) {
      
      Char_t   ringChar = (ir == 0 ? 'I' : 'O');
      UShort_t nsec     = (ir == 0 ? 20  : 40);
      UShort_t nstr     = (ir == 0 ? 512 : 256);
      
      std::cout<<Form("Now saving in FMD%d%c", d, ringChar)<<std::endl;
      
      for(UShort_t vv = 0; vv<nVtxBins ; vv++) { 
	
	TObjArray* detReadArray = (TObjArray*)fReadArray->At(d);
	TObjArray* vtxReadArray = (TObjArray*)detReadArray->At(ir);
	TH2D*      hRead        = (TH2D*)vtxReadArray->At(vv);
	
	TObjArray* detArray = (TObjArray*)fAllArray->At(d);
	TObjArray* vtxArray = (TObjArray*)detArray->At(ir);
	TH2D*      hAll        = (TH2D*)vtxArray->At(vv);
	
	hRead->Divide(hRead,hAll,1,1,"B");
	UShort_t vtxbin = vv+1;
	deadObject->SetCorrection(d,ringChar,vtxbin,hRead);
	
      }
    }
  }
  
}
