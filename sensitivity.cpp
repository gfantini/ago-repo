#include "./coverage.cpp" // include many nice functions we can re-use right now

TGraph ScanLogL(double ExpSignal,double ExpBkg,Double_t& Nsignal_best_fit,bool produceplot = false,bool debug = false,double xmin = -10.,double xmax = 10.,int nbin = 200);

// Given the TRUE values ExpSignal / ExpBkg
// Given the TEST value for TestSignal
// -> Calls ScanLogL from coverage.cpp to build the profile likelihood graph t(s)
// -> Builds histo of profileL(TestSignal)-profileL(BestFit) repeating the procedure NumToyMC times
TH1D BuildTestStatisticsHisto(double ExpSignal,double ExpBkg,double TestSignal,const int NumToyMC,bool debug = false)
{
	TH1D hist_teststat("","",1000,0.,10.); // no-name histo prevents memory leak (?)
	TGraph profileL;
	int i,j;
	Double_t best_fit_signal;
	for(i=0;i<NumToyMC;i++)
	{
		if(debug && i%100 == 0)cout << "[BuildTestStatisticsHisto] DEBUG: cycle "<< i << "/" << NumToyMC << endl;
		profileL = ScanLogL(ExpSignal,ExpBkg,best_fit_signal,false,false);
		hist_teststat.Fill(profileL.Eval(TestSignal)-profileL.Eval(best_fit_signal));
	}
	return hist_teststat;
}


TGraph PrintThresholdAndSensitivity(const int Ntoys,const bool debug = false)
{
	
	// quantiles stuff
	const Int_t nq = 1; // # of quantiles to compute
	Double_t xq[nq];  // position where to compute the quantiles in [0,1]
	xq[0] = .9;
	Double_t yq[nq];  // array to contain the quantiles (thresholds)
	// quantiles stuff (bkg only)
	const Int_t nqq = 5;
	Double_t xqq[nqq]; // probability content
	xqq[0] = (1. - .950)/2.;
	xqq[1] = (1. - .683)/2.;
	xqq[2] = 0.5;
	xqq[3] = .5 + 0.683/2.;
	xqq[4] = .5 + .950/2.;
	Double_t yqq[nqq]; // quantiles of the bkg only hp

	// set of points to scan
	Double_t TrueSignal[] = {1,10,100,1000};
	const Double_t TrueBkg = 1000.;
	
	// set variables to plot
	TGraph gr_cut;
	int i = 0;
	TH1D hist_teststat;
	TH1D hist_solofondo;
	for(auto ss:TrueSignal)
	{
		if(debug)cout << "[FindWhereToCutLogL] DEBUG: ss = " << ss << endl;
		hist_teststat = BuildTestStatisticsHisto(ss,TrueBkg,ss,Ntoys,debug);
		hist_solofondo = BuildTestStatisticsHisto(0.,TrueBkg,ss,Ntoys,debug); // solo fondo!!
		hist_teststat.GetQuantiles(nq,yq,xq);
		hist_solofondo.GetQuantiles(nqq,yqq,xqq);
		gr_cut.SetPoint(i,ss,yq[0]);
		cout << "#" << ss << "\t" << yq[0] << endl; // signal - threshold
		cout << "##" << ss << "\t"; // signal - pval @ -2sigma -1sigma median +1sigma +2sigma
		for(auto qq:yqq)cout << qq << "\t";
		cout << endl;
		i++;
	}
	return gr_cut;
}


