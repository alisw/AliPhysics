void ViewEMCAL(){
    gMC->Gsatt("XEN1","seen",0);  // EMCAL envelope/mother volume
    gMC->Gsatt("XALU","seen",1);  //
    gMC->Gsatt("XPST","seen",0);
    gMC->Gsatt("XPBX","seen",0);
    
}
