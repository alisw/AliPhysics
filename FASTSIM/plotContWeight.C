void plotContWeight()
{
  AliQuenchingWeights afq;
  afq.InitMult();

  afq.PlotContWeights(1,4);
  afq.PlotContWeights(2,1.);
}

