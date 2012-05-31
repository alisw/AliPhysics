
void cambiaNome(){
TSeqCollection *sec = gSystem->GetListOfFileHandlers();
sec->Print();
cout<< "Entries" << sec->GetEntries() << endl;
sec->Draw();
}
