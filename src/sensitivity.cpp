#include "./coverage.cpp" // include many nice functions we can re-use right now

// returns -2log(L(par,nuis^^)/L(par^,nuis^)) where nuis^^ maximize L @ fixed par, par^,nuis^ minimize L as a function of both
// TestSignal == par
// NOTICE: min(-2log(...)) = 0. always
int gGetLogLValue; // global variable holding status of fit: 0/1 = fail / ok If == 0 -> function returns < 0.
Double_t GetLogLValue(double ExpSignal,double ExpBkg,double TestSignal,bool debug = false,double xmin = -10.,double xmax = 10.,int nbin = 200)
{
	
	TH1D h_data = GenerateDataFromExpectation(ExpSignal,ExpBkg,xmin,xmax,nbin);
	
	// define function to fit data -> WARNING! need to normalize with binwidth
	// The correct function to perform a fit is <binwidth>*f(x) (with Integral(f dx) = 1 over fit range)
	double binwidth = (xmax - xmin)/(double)nbin;
	TF1 f_fit("f_fit",Form("[0]/(%lf)+[1]*(%lf)*TMath::Gaus(x,[2],[3],true)",(double)nbin,binwidth),xmin,xmax);
	f_fit.SetParNames("Fbkg","Fsig","Fmean","Fsigma");
	f_fit.SetParameter("Fbkg",ExpBkg); // set starting point
	f_fit.SetParLimits(0,0.,1e6);
	f_fit.SetParameter("Fsig",ExpSignal); // set starting point
	f_fit.SetParLimits(1,0.,1e6);
	f_fit.FixParameter(2,0.);
	f_fit.FixParameter(3,1.);
	
	
	// fit data and retrieve results
	double fcn_min,Fs,Fb,Fs_err,Fb_err;
	TFitResultPtr fitresult = h_data.Fit(&f_fit,"QSRL0"); // Q = quiet S = get res R = range of TH1 L = likelihood 0 = do not plot fit result
	Int_t fit_status_code; // this is 0 if the fit is ok, else see documentation https://root.cern.ch/doc/master/classTH1.html (Fit Status)
	fit_status_code = (Int_t) fitresult;
	if(fit_status_code != 0)
	{
		cerr << "[GetLogLValue] ERROR in main fit! ExpSignal = " << ExpSignal << " ExpBkg = " << ExpBkg << endl;
		gGetLogLValue = 0;
		return -1.;
	}else{
		gGetLogLValue = 1;
	}
	fcn_min = fitresult->MinFcnValue(); // get -logL value @ minimum
	Fb = fitresult->Parameter(0);
	Fb_err = fitresult->ParError(0);
	Fs = fitresult->Parameter(1);
	Fs_err = fitresult->ParError(1);

	
	
	if(debug)
	{
		cout << "[GetLogLValue] DEBUG: SIGNAL(fit): " << Fs << " +- " << Fs_err << endl;
		cout << "[GetLogLValue] DEBUG: BACKGROUND(fit): " << Fb << " +- " << Fb_err << endl;
		cout << "[GetLogLValue] DEBUG: FCN(fit): " << fcn_min << endl;
	}
	
	// repeat fit for the alternate value == TestSignal
	// fix the signal to the alternate_value
	f_fit.FixParameter(1,TestSignal);
	fitresult = h_data.Fit(&f_fit,"QSRLN"); // Q = quiet S = get res R = range of TH1 L = likelihood N = do not screw my plots
	fit_status_code = (Int_t) fitresult;
	if(fit_status_code != 0)cerr << "[GetLogLValue] ERROR in \"profile\" fit! ExpSignal = " << ExpSignal << " ExpBkg = " << ExpBkg << " alternate_signal = " << TestSignal << endl;
	return 2.*(fitresult->MinFcnValue() - fcn_min); // return -2logL(par,nuis^^) + 2logL(par^,nuis^)
}



// Given the TRUE values ExpSignal / ExpBkg
// Given the TEST value for TestSignal
// -> Calls GetLogLValue to get the profile likelihood @ TestSignal
// -> Builds histo of profileL(TestSignal) repeating the procedure NumToyMC times
TH1D BuildTestStatisticsHisto(double ExpSignal,double ExpBkg,double TestSignal,const int NumToyMC,bool debug = false)
{
	vector<Double_t> test_stat;
	Double_t tmp;
	Double_t test_max = 0.;
	int i,j;
	for(i=0;i<NumToyMC;i++)
	{
		if(debug && i%100 == 0)cout << "[BuildTestStatisticsHisto] DEBUG: cycle "<< i << "/" << NumToyMC << endl;
		tmp = GetLogLValue(ExpSignal,ExpBkg,TestSignal,false);
		if(gGetLogLValue==1)
		{
			if(tmp > test_max)test_max=tmp; // compute maximum value of test stat
			test_stat.push_back(tmp);
		}else{
			if(debug)cerr << "[BuildTestStatisticsHisto] WARNING: GetLogLValue main fit failed." << endl;
		}
	}
	TH1D hist_teststat("","",(int)(1.1*test_max/0.01),0.,1.1*test_max); // no-name histo prevents memory leak (?)
	for(i=0;i<test_stat.size();i++)
	{
		hist_teststat.Fill(test_stat[i]);
	}
	
	return hist_teststat;
}

// integrate from xmin to end of histo = return sum of bin contents
double GetIntegralTH1D(TH1D h,double xmin)
{
	TAxis *axis = h.GetXaxis();
	int bmin = axis->FindBin(xmin);
	double integral = h.Integral(bmin,h.GetNbinsX());
	integral -= h.GetBinContent(bmin)*(xmin-axis->GetBinLowEdge(bmin))/axis->GetBinWidth(bmin);
	return integral;
}
// integrate the whole histo excluded under/overflow
double GetIntegralTH1D(TH1D h)
{
	return h.Integral(1,h.GetNbinsX());
}
// ratio of integrals (xmin --> end) / (begin --> end)
double GetPvalueTH1D(TH1D h,double xmin)
{
	return GetIntegralTH1D(h,xmin)/GetIntegralTH1D(h);
}


// Given ExpSignal, ExpBkg, TestSignal writes to stdout
// #B = ExpBkg
// #SB  ExpSignal -2logL_threshold == 90% quantile of pdf(test-stat | s+b)
// #-B  ExpSignal -2sigma -1sigma median +1sigma +2sigma --> p-values of pdf(test-stat | s+b) @ given quantiles of pdf(test-stat | b)
void PrintThresholdAndSensitivity(const double ExpSignal,const double ExpBkg,const double TestSignal,const int Ntoys,const bool debug = false, const bool save = false)
{
	cout << "#B \t" << ExpBkg << endl;
	
	TH1D hist_teststat;
	TH1D hist_solofondo;

	// s + b hypothesis -> threshold computation
	// - - - - - - - - -  - - - - - - - - - - - - - - - -  - - - - -- - - - -
	// quantiles stuff
	const Int_t nq = 1; // # of quantiles to compute
	Double_t xq[nq];  // position where to compute the quantiles in [0,1]
	xq[0] = .9;
	Double_t yq[nq];  // array to contain the quantiles (thresholds)
	hist_teststat = BuildTestStatisticsHisto(ExpSignal,ExpBkg,ExpSignal,Ntoys,debug);
	hist_teststat.GetQuantiles(nq,yq,xq);
	cout << "#SB \t" << ExpSignal << "\t" << yq[0] << endl; // signal - threshold

	// b only hypothesis -> sensitivity computation
	// - - - - - - - - -  - - - - - - - - - - - - - - - -  - - - - -- - - - -
	// quantiles stuff (bkg only)
	const Int_t nqq = 5;
	Double_t xqq[nqq]; // probability content
	xqq[0] = (1. - .950)/2.;
	xqq[1] = (1. - .683)/2.;
	xqq[2] = 0.5;
	xqq[3] = .5 + 0.683/2.;
	xqq[4] = .5 + .950/2.;
	Double_t yqq[nqq]; // quantiles of the bkg only hp
	hist_solofondo = BuildTestStatisticsHisto(0.,ExpBkg,ExpSignal,Ntoys,debug); // only bkg!!
	hist_solofondo.GetQuantiles(nqq,yqq,xqq);
	cout << "#-B \t" << ExpSignal << "\t"; // signal - pvalue @ -2sigma -1sigma median +1sigma +2sigma
	for(auto qq:yqq)cout << GetPvalueTH1D(hist_teststat,qq) << "\t";
	cout << endl;
	
	// save histo in order to plot later
	TFile *fp;
	if(save)
	{
		fp = new TFile(Form("s%d.root",(int)ExpSignal),"new");
		if(fp->IsOpen())
		{
			hist_solofondo.Write("b");
			hist_teststat.Write("sb");
			fp->Close();
			delete fp;
		}else{
			cerr << "[PrintThresholdAndSensitivity] ERROR: could not open ROOT file." << endl;
		}
	}
}

