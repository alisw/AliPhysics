void AliITSDigits2Clusters(const char *filename="galice.root"){
    AliITSreconstruction *itsr = new AliITSreconstruction(filename);
    itsr->Exec(); 
    delete itsr;
}
