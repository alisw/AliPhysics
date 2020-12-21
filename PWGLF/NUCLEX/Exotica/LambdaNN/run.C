void CreateHistos();
void DrawPlots();

void run(){

 gROOT->LoadMacro("Getoutput.cxx++");
 //CreateHistos();
 DrawPlots();

}

void CreateHistos(){
 Getoutput g;
 if( !g.LoadParams("params.txt")){
  printf("no param file, exiting... \n");
  return;
 }
 if(!g.LoadFile("LNNntupleTrd.root")){
  printf("proble with input file, exiting... \n");
  return;
 }
 g.BookOutputData();
 g.LoopOverV0();
 g.LoopOverV0(1);
 g.StoreOutputData("results.root");


}
void DrawPlots(){

 Getoutput d;
 if(!d.LoadOutputData("results3H1.root")){
  printf("issue wen loading analized data file, exiting...\n");
  return;

 }
 d.DrawResults();


}
