// random generator
TRandom3 rangen;

// returns an histogram ("h_data") with EXACTLY Nsignal gaussian entries mu = 0, sigma = 1
// and EXACTLY Nbkg background entries
// you can tune the range and the binning with xmin, xmax, nbin
TH1D GenerateData(int Nsignal,int Nbkg,double xmin = -10.,double xmax = 10.,int nbin = 200)
{
	
	// define the gaussian shape of the signal
	TF1 f_signal("f_signal","TMath::Gaus(x,[0],[1],true)",xmin,xmax);
	f_signal.SetParNames("mean","sigma");
	f_signal.SetParameter("mean",0.);
	f_signal.SetParameter("sigma",1.);
	
	// histogram to store data
	TH1D h_data("h_data","Data histogram",nbin,xmin,xmax);
	
	// simulate signal + background
	h_data.FillRandom("f_signal",Nsignal); // inject signal
	double tmp;
	for(int i=0;i<Nbkg;i++) // inject background
	{
		tmp = rangen.Rndm();
		tmp = xmin + tmp*TMath::Abs(xmax-xmin);
		h_data.Fill(tmp);
	}
	return h_data;
}

// The same as GenerateData but Nsignal and Nbkg are extracted Poisson from the expected values
TH1D GenerateDataFromExpectation(double ExpSignal,double ExpBkg,double xmin = -10.,double xmax = 10.,int nbin = 200)
{
	int trueS,trueB;
	trueS = rangen.Poisson(ExpSignal);
	trueB = rangen.Poisson(ExpBkg);
	if(trueS == 0 && trueB == 0)cout << "[GenerateDataFromExpectation] WARNING: S(true) = 0 and B(true) = 0 -> empty TH1D." << endl;
	return GenerateData(trueS,trueB,xmin,xmax,nbin);
}

// Generates data with "GenerateDataFromExpectation"
// scans -logL around the minimum for the Nsignal parameter (+- 3 sigma_fit) and returns TGraph of -2logL (profile-likelihood-ratio)
// COMMENT: By def. the prof. likelihood ratio is -2log{L(Ns,nuisance^^)/L(Ns^,nuisance^)} = -2logL(Ns,nuisance^^) + cost --> drop the cost.
// Fills Nsignal_best_fit with the signal value that minimizes -2logL
TGraph ScanLogL(double ExpSignal,double ExpBkg,Double_t& Nsignal_best_fit,bool produceplot = false,bool debug = false,double xmin = -10.,double xmax = 10.,int nbin = 200)
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
	
	TCanvas *c1;
	if(produceplot)
	{
		c1 = new TCanvas();// otherwise once it goes out of scope I am lost
					    //c1->Divide(2); // divide in 2 pads horizontally
	}
	
	// fit data and retrieve results
	double fcn_min,Fs,Fb,Fs_err,Fb_err;
	TFitResultPtr fitresult = h_data.Fit(&f_fit,"QSRL"); // Q = quiet S = get res R = range of TH1 L = likelihood
	Int_t fit_status_code; // this is 0 if the fit is ok, else see documentation https://root.cern.ch/doc/master/classTH1.html (Fit Status)
	fit_status_code = (Int_t) fitresult;
	if(fit_status_code != 0)cout << "[ScanLogL] ERROR in main fit! ExpSignal = " << ExpSignal << " ExpBkg = " << ExpBkg << endl;
	fcn_min = fitresult->MinFcnValue(); // get -logL value @ minimum
	Fb = fitresult->Parameter(0);
	Fb_err = fitresult->ParError(0);
	Fs = fitresult->Parameter(1);
	Fs_err = fitresult->ParError(1);
	if(produceplot){
		//c1->cd(1);
		h_data.DrawClone();
	}

	
	if(debug)
	{
		cout << "SIGNAL: " << Fs << " +- " << Fs_err << endl;
		cout << "BACKGROUND: " << Fb << " +- " << Fb_err << endl;
		cout << "FCN: " << fcn_min << endl;
	}
	// write result of the best fit value for the signal component
	Nsignal_best_fit = Fs;
	// scan the likelihood for the signal
	vector <double> alternate_signal;
	vector <double> alternate_fcn;
	double step = TMath::Min(Fs_err/20.,sqrt(TMath::Abs(Fs))/20. ); // define step for scanning likelihood
	double alternate_signal_tmp = Fs - 10.*Fs_err; // define starting point for scanning likelihood
	while(alternate_signal_tmp < Fs + 10.*Fs_err)
	{
		// fix the signal to the alternate_value
		f_fit.FixParameter(1,alternate_signal_tmp);
		// perform fit again
		fitresult = h_data.Fit(&f_fit,"QSRLN"); // Q = quiet S = get res R = range of TH1 L = likelihood N = do not screw my plots
		if(fit_status_code != 0)cout << "[ScanLogL] ERROR in \"profile\" fit! ExpSignal = " << ExpSignal << " ExpBkg = " << ExpBkg << " alternate_signal = " << alternate_signal_tmp << endl;
		alternate_fcn.push_back(fitresult->MinFcnValue()); // store FCN= -logL value @ new minimum
		alternate_signal.push_back(alternate_signal_tmp);  // store signal used to compute FCN = -logL
		alternate_signal_tmp += step;
	}
	
	// make TGraph of the results
	TGraph graph;
	int i;
	for(i=0;i<alternate_signal.size();i++)
	{
		graph.SetPoint(i,alternate_signal[i],2.*alternate_fcn[i]);
	}
	graph.SetMarkerStyle(kFullCircle);
	TCanvas *c2;
	if(produceplot)
	{
		// dunno why cannot manage to put those 2 plots in the same canvas
		c2 = new TCanvas();
		//c1->cd(2);
		graph.DrawClone("AP");
	}
	
	return graph;
}

// Computes a scan of the -2logL with the input parameters
// Gets the confidence interval @ the specified confidence level
// returns true if the ExpSignal is inside confidence interval, false otherwise
bool CheckCoverage(double CL,double ExpSignal,double ExpBkg,bool ProducePlot=false,bool VerboseDebug=false,double xmin = -10.,double xmax = 10.,int nbin = 200)
{
	// consistency check
	if(CL <= 0. || CL > 1.)
	{
		cout << "[CheckCoverage] ERROR: CL not in 0. <-> 1. Setting to .9" << endl;
		CL = .9;
	}
	Double_t best_fit_signal; // will be filled by ScanLogL
	TGraph profile_likelihood_ratio = ScanLogL(ExpSignal,ExpBkg,best_fit_signal,ProducePlot,false,xmin,xmax,nbin);
	//if(ProducePlot)profile_likelihood_ratio.DrawClone("AP");
	int npts = profile_likelihood_ratio.GetN(); // store here the # of points in the graph
	Double_t ll_min = profile_likelihood_ratio.Eval(best_fit_signal);
	Double_t ll_threshold = TMath::ChisquareQuantile(CL,1);
	if(VerboseDebug)
	{
		cout << "[CheckCoverage] DEBUG: ll_threshold = "<< ll_threshold << endl;
		cout << "[CheckCoverage] DEBUG: ll_min = "<< ll_threshold << endl;
		cout << "[CheckCoverage] DEBUG: par_min = " << profile_likelihood_ratio.GetX()[0] << endl;
		cout << "[CheckCoverage] DEBUG: par_max = " << profile_likelihood_ratio.GetX()[npts-1] << endl;

	}
	// ^^^^ WARNING the convergence (Wilks' theorem) is assumed! Not true for few events
	
	double diff;
	double diff_old;
	double sig_min = profile_likelihood_ratio.GetX()[2];
	double sig_max = profile_likelihood_ratio.GetX()[npts-1];
	double sig_step = (sig_max - sig_min)/((double)npts);
	if(VerboseDebug)cout << "[CheckCoverage] DEBUG: par_step = " << sig_step << endl;
	double sig_now = sig_min+sig_step;
	int flagmin = 0;
	int flagmax = 0;
	while(sig_now < sig_max)
	{
		if(VerboseDebug)cout << "[CheckCoverage] DEBUG: sig_now / diff(now-min) = " << sig_now << " / " << diff << endl;
		diff_old=profile_likelihood_ratio.Eval(sig_now-sig_step)-ll_min;
		diff=profile_likelihood_ratio.Eval(sig_now)-ll_min;
		if(diff_old >= ll_threshold && diff <= ll_threshold)
		{
			// update minimum of interval
			// this should happen only once
			flagmin++;
			sig_min = sig_now;
		}
		if(diff_old <= ll_threshold && diff >= ll_threshold)
		{
			// update maximum of interval
			// this should happen only once
			flagmax++;
			sig_max = sig_now;
		}
		sig_now += sig_step;
	}
	// consistency checks
	if(flagmin == 0 && flagmax == 0)
	{
		cout << "[CheckCoverage] WARNING: -2logL threshold never crossed. Setting interval = extrema of TGraph." << endl;
		cout << "[CheckCoverage] CL = " << CL << " ExpSignal  = " << ExpSignal << " ExpBkg = " << ExpBkg << endl;
	}else if(flagmin == 0 && flagmax == 1)
	{
		cout << "[CheckCoverage] WARNING: -2logL threshold (left) not crossed. Setting interval minimum = extremum of TGraph." << endl;
	}else if(flagmin == 1 && flagmax == 0)
	{
		cout << "[CheckCoverage] WARNING: -2logL threshold (right) not crossed. Setting interval maximum = extremum of TGraph." << endl;
	}else if(flagmin != 1 && flagmax != 1)
	{
		cout << "[CheckCoverage] ERROR: -2logL threshold multiple (pathologic) crossing!" << endl;
	}
	if(sig_min >= sig_max)cout << "[CheckCoverage] ERROR: min >= max !! Interval is WRONG!" << endl;
	
	
	if(VerboseDebug)cout << "[CheckCoverage] DEBUG: interval is "<< sig_min << " <-> " << sig_max << endl;
	if(sig_min <= ExpSignal && ExpSignal <= sig_max)
	{
		return true;
	}else{
		return false;
	}
}

// input number of bkg events
void exe(int Nb,bool debug)
{

	//TCanvas *c2 = new TCanvas();
	
	int nhit;
	int ntot = 1e4;
	TGraphErrors gr_coverage;
	double p;
	int i,jj;
	jj=0;
	for(int Ns:{10,100,500,1000})
	{
		if(debug)cout << "[EXE] DEBUG: Iteration Ns = " << Ns << endl;
		nhit = 0;
		for(i=0;i<ntot;i++)
		{
			if( CheckCoverage(0.9,Ns,Nb) )nhit++;
			if((i+1)%100==0 && debug)cout << "[EXE] DEBUG: Sub-iteration i / ntot = " << i+1 << " / " << ntot << endl;
		}
		p = nhit/(double)ntot;
		gr_coverage.SetPoint(jj,(Double_t)Ns,p);
		gr_coverage.SetPointError(jj,0.,sqrt(p*(1.-p)/(double)ntot));
		jj++;
	}
	gr_coverage.DrawClone("AP");
	//grafico.DrawClone("AP*");
}

// repeats the CheckCoverage for the specified
// CL -> confidence level
// Ns -> number of true signal events
// Nb -> number of true background events
// Ntry -> number of repetitions
// prints out some numbers: Ns, coverage
// on another line: err-coverage (1/sqrt(Ntry) * sqrt(coverage*(1-coverage)))
void ComputeCoverage(double CL,double Ns,double Nb,int Ntry)
{
	int Nhit = 0;
	int i;
	for(i=0;i<Ntry;i++)if( CheckCoverage(CL,Ns,Nb) )Nhit++;
	double c = Nhit/(double)Ntry; // coverage
	cout << Ns << "," << c << endl;
	cout << "Coverage error: " << sqrt(c*(1.-c))/sqrt(Ntry) << endl;
}

void bug()
{
	cout << "Looking for a bug in 23rd iteration of CheckCoverage(0.9,10,100)" << endl;
	for(int i = 0;i<18;i++){
		CheckCoverage(0.9,10,100);
		cout << "[bug()] i = " << i << endl;
	}
	CheckCoverage(0.9,10,100,true,true);
	// 1st point of likelihood scan screwed up
	// why error is so big? --> interval of TGraph so big? (this can affect the step)
}
