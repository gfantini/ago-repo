/*!
 * some global variables to define fixed values of the model
 * please call global variables g_<name>
 */
const Double_t g_xMin = 0.;
const Double_t g_xMax = 200.;
const Int_t g_xNbin = 1000;
const Double_t g_xBinWidth = TMath::Abs(g_xMax-g_xMin)/(double)g_xNbin;


void correbias(Double_t N1_true,Double_t N2_true,Double_t N3_true)
{
	// background model
	TF1 pdf1("pdf1",Form("%lf*(TMath::Exp(-x/200.)/(200.*%lf))",g_xBinWidth,(TMath::E()-1.)/TMath::E()),g_xMin,g_xMax);
	pdf1.SetNpx(1e2*g_xNbin);

	//	pdf1.DrawCopy();
	
	// signal model
	TF1 pdf2("pdf2",Form("%lf*(TMath::Gaus(x,50,30,true)",g_xBinWidth),g_xMin,g_xMax);
	TF1 pdf3("pdf3",Form(),g_xMin,g_xMax);
	
}
