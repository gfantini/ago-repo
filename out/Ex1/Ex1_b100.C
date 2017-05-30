void Ex1_b100()
{
//=========Macro generated from canvas: c1/c1
//=========  (Tue May 30 10:02:29 2017) by ROOT version6.08/04
   TCanvas *c1 = new TCanvas("c1", "c1",9,67,1116,528);
   c1->Range(-1.843808,0.8796212,3.583973,0.9656344);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetLogx();
   c1->SetFrameBorderMode(0);
   c1->SetFrameBorderMode(0);
   
   Double_t _fx1001[6] = {
   0.1,
   5,
   10,
   100,
   500,
   1000};
   Double_t _fy1001[6] = {
   0.950015,
   0.9183,
   0.8926,
   0.9028,
   0.903,
   0.9012};
   Double_t _fex1001[6] = {
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t _fey1001[6] = {
   0.00218318,
   0.002739071,
   0.003096211,
   0.002962299,
   0.002959578,
   0.002983933};
   TGraphErrors *gre = new TGraphErrors(6,_fx1001,_fy1001,_fex1001,_fey1001);
   gre->SetName("");
   gre->SetTitle("");
   gre->SetFillColor(1);
   gre->SetMarkerStyle(24);
   
   TH1F *Graph_Graph_Graph10011001 = new TH1F("Graph_Graph_Graph10011001","",100,0.05,1099.5);
   Graph_Graph_Graph10011001->SetMinimum(0.8882225);
   Graph_Graph_Graph10011001->SetMaximum(0.9570331);
   Graph_Graph_Graph10011001->SetDirectory(0);
   Graph_Graph_Graph10011001->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   Graph_Graph_Graph10011001->SetLineColor(ci);
   Graph_Graph_Graph10011001->GetXaxis()->SetLabelFont(42);
   Graph_Graph_Graph10011001->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph10011001->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph10011001->GetXaxis()->SetTitleFont(42);
   Graph_Graph_Graph10011001->GetYaxis()->SetLabelFont(42);
   Graph_Graph_Graph10011001->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph10011001->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph10011001->GetYaxis()->SetTitleFont(42);
   Graph_Graph_Graph10011001->GetZaxis()->SetLabelFont(42);
   Graph_Graph_Graph10011001->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph10011001->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph10011001->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph_Graph10011001);
   
   gre->Draw("ap");
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
