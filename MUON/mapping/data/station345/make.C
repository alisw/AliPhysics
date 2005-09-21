void make()
{
  std::ofstream out("padPosG.dat");

  int n = 1;

  for ( int ix = 0; ix < 16; ++ix )
    {
      for ( int iy = 0; iy < 3; ++ iy )
	{
	  out << n << "\t" << ix << "\t" << iy << std::endl;
	  ++n;
	}
    }
  out.close();

}
