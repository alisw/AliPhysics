//void ViewTRD()
{
   geant3->Gsatt("TRD ","seen",0);
   geant3->Gsatt("UTRS","seen",0);
   geant3->Gsatt("UTRI","seen",0);

   geant3->Gsatt("UTCI","seen",0);
   geant3->Gsatt("UTCN","seen",0);
   geant3->Gsatt("UTCO","seen",0);

   geant3->Gsatt("UTII","seen",0);
   geant3->Gsatt("UTIN","seen",0);
   geant3->Gsatt("UTIO","seen",0);

   geant3->Gsatt("UTMI","seen",0);
   geant3->Gsatt("UTMN","seen",0);
   geant3->Gsatt("UTMO","seen",0);

   geant3->Gsatt("UT1I","seen",1);
   geant3->Gsatt("UT1N","seen",1);
   geant3->Gsatt("UT1O","seen",1);

   geant3->Gsatt("UT4I","seen",1);
   geant3->Gsatt("UT4N","seen",1);
   geant3->Gsatt("UT4O","seen",1);

}
