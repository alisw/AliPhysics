void generateFlipSlatFile(int line)
{
  char filename[80];

  const char* planes[] = { "Bending","NonBending" };

  for ( int iplane = 0; iplane  < 2; ++iplane)
    {
      for ( int i = 1; i <=4; ++i )
	{
	  if ( line > 4 )
	    {
	      // make a file consisting into flip-x right-slat.
	      sprintf(filename,"%dLL%d.%s.slat",i,line,planes[iplane]);
	      ofstream out(filename);
	      out << "FLIP_X " << i << "RL" << line << endl;
	      out.close();
	    }
	  else if ( line <= 4 )
	    {
	      // make a file consisting into flip-y top-slat.
	      const char* lr[] = { "L","R" };
	      for ( int j = 0; j < 2; ++j )
		{
		  sprintf(filename,"%d%sL%d.%s.slat",
			  i,lr[j],line,planes[iplane]);
		  ofstream out(filename);
		  out << "FLIP_Y " << i << lr[j] << "L" << (10-line) << endl;
		  out.close();
		}
	    }
	}
    }
}
