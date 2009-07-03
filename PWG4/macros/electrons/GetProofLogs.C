{
 gEnv->SetValue("XSec.GSI.DelegProxy","2");
 logs=TProof::Mgr("kread@localhost")->GetSessionLogs(0);
 logs->Display();
 //logs->Grep("segmentation violation");
 logs->Save("*","logs.txt");
}
