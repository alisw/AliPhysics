{
Int_t evt = 0 ; 
RecAna * t = new RecAna("junk.root"); 
t->GetEvent(evt);   
TObjArray * lp = t->PHOSPpsdRP ; 
 cout << "Tree macro = " << lp << endl ; 
for (int i = 0 ; i < lp->GetEntries() ; i++ ) {
  AliPHOSPpsdRecPoint * rpp = (AliPHOSPpsdRecPoint *)lp->At(i)   ; 
  rpp.Print(); 
}
TObjArray * le = t->PHOSEmcRP ; 
for (int i = 0 ; i < le->GetEntries() ; i++ ) {
  AliPHOSEmcRecPoint * rp = (AliPHOSEmcRecPoint *)le->At(i)   ; 
  rp->Print(); 
}

AliPHOSIndexToObject * please = AliPHOSIndexToObject::GetInstance() ;

for (int i = 0 ; i < (t->PHOSTS_-1) ; i++) {
  cout << "TrackSegment # " << i << endl 
       << "====================" << endl ; 
  int index = t->PHOSTS_fEmcRecPoint[i] ;
  AliPHOSEmcRecPoint * emrp = (AliPHOSEmcRecPoint *) ( please->GimeRecPoint(index, TString("emc") ) ) ; 
  emrp->Print() ; 
  index = t->PHOSTS_fPpsdLowRecPoint[i] ;
  AliPHOSPpsdRecPoint * ppsdl = (AliPHOSPpsdRecPoint *) ( please->GimeRecPoint(index, TString("ppsd") ) ) ; 
  if (ppsdl)
    ppsdl->Print() ; 
  index = t->PHOSTS_fPpsdUpRecPoint[i] ;
  AliPHOSPpsdRecPoint * ppsdu = (AliPHOSPpsdRecPoint *) ( please->GimeRecPoint(index, TString("ppsd") ) ) ; 
  if (ppsdu)
    ppsdu->Print() ; 
}
for (int i = 0 ; i < (t->PHOSRP_-1) ; i++) {
  cout << "RecParticles # " << i << endl 
       << "====================" << endl ; 
  cout << "type = " << t->PHOSRP_fType[i] << " energy = " << t->PHOSRP_fE[i] << endl ;
}  
delete t ; 
}
