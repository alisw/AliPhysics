void DumpspdTestBeamData(const Char_t *filename,Int_t istart=0,Int_t iend=-1){
    AliITSspdTestBeam *spd = new AliITSspdTestBeam(filename);
    Int_t i;
         //fstream *fp = new fstream("dumpedata.txt",io::out);

    spd->Read();
    spd->Decode();
    if(iend<=0) iend = spd->GetNumberOfEvents();
    for(i=istart;i<iend;i++){
        cout << "********Event=" << i <<" *********"<< endl;
        spd->PrintEventData(i);
    } // end for
         //fp->Close();
}
