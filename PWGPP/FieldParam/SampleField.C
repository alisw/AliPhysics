/**
 * Sample magnet field for specified geometry.
 */

TRandom3 *rnd;

void writeSamples(AliMagF *mag, int count, Double_t *dim)
{
	Double_t rpz[] = {0, 0, 0};
	Double_t xyz[] = {0, 0, 0};
	Double_t fld[] = {0, 0, 0};

	for (int i = 0; i < count; ++i)
	{
		rnd->RndmArray(3, rpz); // (0,1] uniform
		rpz[0] = rpz[0] * (dim[3]*dim[3] - dim[0]*dim[0]) + dim[0]*dim[0];
		rpz[1] = rpz[1] * (dim[4] - dim[1]) + dim[1];
		rpz[2] = rpz[2] * (dim[5] - dim[2]) + dim[2];
		xyz[0] = sqrt(rpz[0]) * cos(rpz[1]);
		xyz[1] = sqrt(rpz[0]) * sin(rpz[1]);
		xyz[2] = rpz[2];

		mag->Field(xyz, fld);
		printf("%f %f %f %f %f %f\n", xyz[0], xyz[1], xyz[2], fld[0], fld[1], fld[2]);
	}
}

void SampleField(int rows = 100000, int nkGauss = 5,
	             double r_min, double p_min, double z_min,
	             double r_max, double p_max, double z_max)
{
	rnd = new TRandom3(0);
	AliMagF *k2 = new AliMagF("k2", "k2", 1, 1, AliMagF::k2kG);
	AliMagF *k5 = new AliMagF("k5", "k5", 1, 1, AliMagF::k5kG);
	AliMagF *gauss;

	if (nkGauss == 5)
	{
		gauss = k5;
	}
	else if (nkGauss == 2)
	{
		gauss = k2;
	}
	else
	{
		fprintf(stderr, "%d kilo gauss field is not supported.\n", nkGauss);
	}

	Double_t dim[] = {r_min, p_min, z_min, r_max, p_max, z_max};
	writeSamples(gauss, rows, dim);
}
