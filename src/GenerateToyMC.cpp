// THIS MACRO STUDIES THE RESOLUTION / CORRELATION OF 3 PARAMETERS IN A FIT: # BKG # SIGNAL1 #SIGNAL2

/*!
 * some global variables to define fixed values of the model
 * please call global variables g_<name>
 */
const Double_t g_xMin = 0.;
const Double_t g_xMax = 200.;
const Int_t g_xNbin = 1000;
const Double_t g_xBinWidth = TMath::Abs(g_xMax-g_xMin)/(double)g_xNbin;
// The Truth
const Double_t g_N1true = 1e6;
const Double_t g_N2true = 1e5;
const Double_t g_N3true = 1e5;

TRandom3 g_r3;

// generates a ROOT file where it stores Ntoys fit results for the 3 components of the model + outcome flag of the fit
// each dataset (toy) is generated from the true model (exp bkg + gaus 1 + gaus 2)
// the fit is performed with mean and/or sigma evenually modified if the multiplier (argument) is set (!= 1.)
void GenerateToyMC(TString output_filename,const int Ntoys,const double gauss_mean_multiplier = 1.,const double gauss_sigma_multiplier = 1.)
{
	// TRUE model == used for generation
	// background model
	TF1 pdf1_true("pdf1_true",Form("(%lf)*(TMath::Exp(-x/200.)/(200.*(%lf)))",g_xBinWidth,(TMath::E()-1.)/TMath::E()),g_xMin,g_xMax);
	// signal model
	TF1 pdf2_true("pdf2_true",Form("(%lf)*TMath::Gaus(x,50,30,true)",g_xBinWidth),g_xMin,g_xMax);
	TF1 pdf3_true("pdf3_true",Form("(%lf)*TMath::Gaus(x,60,10,true)",g_xBinWidth),g_xMin,g_xMax);
	for(auto f : {pdf1_true,pdf2_true,pdf3_true} )f.SetNpx(1e2*g_xNbin);//set 100 pts for every bin
	
	
	// FIT model == used for fitting
	// background model
	TF1 pdf1("pdf1",Form("(%lf)*(TMath::Exp(-x/200.)/(200.*(%lf)))",g_xBinWidth,(TMath::E()-1.)/TMath::E()),g_xMin,g_xMax);
	// signal model
	TF1 pdf2("pdf2",Form("(%lf)*TMath::Gaus(x,50*(%lf),30*(%lf),true)",g_xBinWidth*1.05019,gauss_mean_multiplier,gauss_sigma_multiplier),g_xMin,g_xMax);
	TF1 pdf3("pdf3",Form("(%lf)*TMath::Gaus(x,60*(%lf),10*(%lf),true)",g_xBinWidth,gauss_mean_multiplier,gauss_sigma_multiplier),g_xMin,g_xMax);
	for(auto f : {pdf1,pdf2,pdf3} )f.SetNpx(1e2*g_xNbin);//set 100 pts for every bin
	// implement lambda function that sees all variables by reference [&]
	// ah.. if you have a TF1 funz --> funz(.3) == funz.Eval(.3) did you know that? I didn't...
	auto lalalambda = [&](double *x,double *p){return p[0]*pdf1(x) + p[1]*pdf2(x) + p[2]*pdf3(x);};
	const int NumberOfParameters = 3;
	TF1 model("model",lalalambda,g_xMin,g_xMax,NumberOfParameters);
	model.SetParNames("N1","N2","N3");
	model.SetParameter("N1",g_N1true);
	model.SetParameter("N2",g_N2true);
	model.SetParameter("N3",g_N3true);

	// histogram for the generated data
	TH1D h_data("","Generated dataset",g_xNbin,g_xMin,g_xMax);
	TFitResultPtr frp; // fit result pointer
	Int_t fsc; // fit status code: 0 <-> fit converged https://root.cern.ch/doc/master/classTH1.html (Fit Status)
	Double_t fN1; // fit result for N1
	Double_t fN2; //  '    '     '  N2
	Double_t fN3; //  '    '     '  N3
	
	// TTree to store fit results
	TFile fileroot(output_filename,"recreate");
	TTree t3("tree","");
	t3.Branch("fsc",&fsc); // create branch and link it to the specified variable
	t3.Branch("fN1",&fN1);
	t3.Branch("fN2",&fN2);
	t3.Branch("fN3",&fN3);
	for(int i=0;i<Ntoys;i++)
	{
		// fill histo with random numbers according to Poisson(The Truth)
		h_data.FillRandom("pdf1", g_r3.Poisson(g_N1true) );
		h_data.FillRandom("pdf2", g_r3.Poisson(g_N2true) );
		h_data.FillRandom("pdf3", g_r3.Poisson(g_N3true) );
		// fit the data with the model
		frp = h_data.Fit("model","IQLSR0"); // Integral in bin normalized with binwidth, Quiet printout, Likelihood, return reSult, Range of TF1, 0 do not draw
									// * in the V = verbose printout I find error definition = 0.5
									// * I interpret this as the threshold for -logL instead of -2logL..
		fsc = (Int_t) frp; if(fsc != 0)cerr << "[--] ERROR: fit did not converge." << endl;
		fN1 = frp->Parameter(0);
		fN2 = frp->Parameter(1);
		fN3 = frp->Parameter(2);
		if((i+1)%1000==0)
		{
			cerr << "[GenerateToyMC] INFO: " << i+1 << " / " << Ntoys << endl;
			t3.Write(); // update rootfile on disk
		}
		t3.Fill();
		h_data.Reset(); // Reset this histogram: contents, errors, etc.
	}

	fileroot.Close();
	
}

// returns 1 if everything went well
int ProduceCorrelationPlots(const TString input_filename)
{
	TFile f_input(input_filename,"read");
	if(f_input.IsZombie())return 0;
	// read file == get TTree
	TTree *p_t3;
	f_input.GetObject("tree",p_t3);
	// to make scatter plots first use the "Draw" method of TTree (that make TGraphs)
	// fsc == 0 checks the fit has converged properly
	p_t3->Draw("fN1:fN2>>h_12","fsc==0"); // Y-axis : X-axis
	p_t3->Draw("fN1:fN3>>h_13","fsc==0");
	p_t3->Draw("fN2:fN3>>h_23","fsc==0");
	// make 1-D distributions
	p_t3->Draw("fN1>>h_11","fsc==0");
	p_t3->Draw("fN2>>h_22","fsc==0");
	p_t3->Draw("fN3>>h_33","fsc==0");
	// get the histograms just created from the current directory
	TH2D *ph_12 = (TH2D*) gDirectory->Get("h_12");
	TH2D *ph_13 = (TH2D*) gDirectory->Get("h_13");
	TH2D *ph_23 = (TH2D*) gDirectory->Get("h_23");
	TH1D *ph_11 = (TH1D*) gDirectory->Get("h_11");
	TH1D *ph_22 = (TH1D*) gDirectory->Get("h_22");
	TH1D *ph_33 = (TH1D*) gDirectory->Get("h_33");
	// set axis names
	ph_12->GetYaxis()->SetTitle("N1 (fit)");
	ph_12->GetXaxis()->SetTitle("N2 (fit)");
	ph_13->GetYaxis()->SetTitle("N1 (fit)");
	ph_13->GetXaxis()->SetTitle("N3 (fit)");
	ph_23->GetYaxis()->SetTitle("N2 (fit)");
	ph_23->GetXaxis()->SetTitle("N3 (fit)");
	ph_11->GetXaxis()->SetTitle("N1 (fit)");
	ph_22->GetXaxis()->SetTitle("N2 (fit)");
	ph_33->GetXaxis()->SetTitle("N3 (fit)");
	// set style
	gStyle->SetOptFit(1);
	gStyle->SetOptStat(0);
	// fit histo (diagonal)
	ph_11->Fit("gaus","L");
	ph_22->Fit("gaus","L");
	ph_33->Fit("gaus","L");
	
	
	TCanvas *can = new TCanvas();
	can->Divide(3,3);
	can->cd(2);
	ph_12->DrawCopy("COLZ");
	can->cd(3);
	ph_13->DrawCopy("COLZ");
	can->cd(6);
	ph_23->DrawCopy("COLZ");
	can->cd(1);
	ph_11->DrawCopy();
	can->cd(5);
	ph_22->DrawCopy();
	can->cd(9);
	ph_33->DrawCopy();


	return 1;
}

// This computes the normalization of the gaussian due to the finite range in which it is defined
void ComputeIntegral()
{
	TF1 pdf2("pdf2",Form("(%lf)*TMath::Gaus(x,50*(%lf),30*(%lf),true)",g_xBinWidth,1.0,1.0),g_xMin,g_xMax);
	pdf2.SetNpx(1e5);
	TF1 pdf3("pdf3",Form("(%lf)*TMath::Gaus(x,60*(%lf),10*(%lf),true)",g_xBinWidth,1.0,1.0),g_xMin,g_xMax);
	pdf3.SetNpx(1e5);

	Double_t errore2;
	Double_t integrale2 = pdf2.IntegralOneDim(g_xMin,g_xMax,1e-9,1e-9,errore2);
	cout << "PDF2" << endl;
	cout << "binwidth: \t" << g_xBinWidth << endl;
	cout << "integrale: \t" << integrale2 << endl;
	cout << "errore: \t" << errore2 << endl;
	cout << "normaliz: \t" << g_xBinWidth/integrale2 << endl;
	cout << "- - - - - - - - - - - - - - - - - - - -" << endl;

	Double_t errore3;
	Double_t integrale3 = pdf3.IntegralOneDim(g_xMin,g_xMax,1e-9,1e-9,errore2);
	cout << "PDF3" << endl;
	cout << "binwidth: \t" << g_xBinWidth << endl;
	cout << "integrale: \t" << integrale3 << endl;
	cout << "errore: \t" << errore3 << endl;
	cout << "normaliz: \t" << g_xBinWidth/integrale3 << endl;
	cout << "- - - - - - - - - - - - - - - - - - - -" << endl;
}
