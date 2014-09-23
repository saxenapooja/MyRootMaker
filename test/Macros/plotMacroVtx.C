void plotMacro(TString refHist,  TString refitHist, float xmin_, float xmax_, float ymin_, float ymax_, TString xtitle_, TString ytitle_, TString outFileName_, TString inFileName_)
{

  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  gStyle->SetOptStat(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleH(0.07);
  gStyle->SetTitleFontSize(0.1);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleOffset(1.3,"y");
  
  //TFile *f1 = TFile::Open("GluGluToHToTauTau_M-125_MC_v4_anal.root", "READONLY");
  //TFile *f1 = TFile::Open("VBF_HToTauTau_M-125_MC_v4_anal.root", "READONLY");
  TFile *f1 = TFile::Open(inFileName_, "READONLY");
  
  TH1F* hPVorig = (TH1F*)f1->Get(refHist);
  hPVorig->GetXaxis()->SetRangeUser(xmin_, xmax_);
  //hPVorig->GetYaxis()->SetRangeUser(ymin_, ymax_);
  hPVorig->SetLineColor(kRed);
  hPVorig->GetXaxis()->SetTitle(xtitle_);
  hPVorig->GetYaxis()->SetTitle(ytitle_);

  TH1F* hPVrefit = (TH1F*)f1->Get(refitHist);
  hPVrefit->SetLineColor(kBlue);
  
  hPVorig->Draw();
  hPVrefit->Draw("same");

  TLegend* leg = new TLegend(0.6,0.75,0.8,0.88,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.03);
  leg->SetHeader("#splitline{CMS Simulation}{ #sqrt{s}=8 TeV}");

  leg->AddEntry(hPVorig,"Original PV");
  leg->AddEntry(hPVrefit,"Re-fitted PV");
  
  leg->Draw();
  
  c1->SaveAs("plots/"+outFileName_+".png");
}

void plotMacroAll()
{
  /*
  plotMacro("VtxXRes",  "ReFitVtxXRes", 0.0, 0.04, 0., 1.0, "Resolution X-position", "Arbitrary Units", "GGFH_VtxResolution_X");
  plotMacro("VtxYRes",  "ReFitVtxYRes", 0.0, 0.04, 0., 1.0, "Resolution Y-position", "Arbitrary Units", "GGFH_VtxResolution_Y");
  plotMacro("VtxZRes",  "ReFitVtxZRes", -0.02, 0.02, 0., 1.0, "Resolution Z-position", "Arbitrary Units", "GGFH_VtxResolution_Z");
  */
  /*
  plotMacro("VtxXRes",  "ReFitVtxXRes", 0.0, 0.04, 0., 1.0, "Resolution X-position", "Arbitrary Units", "VBFH_VtxResolution_X");
  plotMacro("VtxYRes",  "ReFitVtxYRes", 0.0, 0.04, 0., 1.0, "Resolution Y-position", "Arbitrary Units", "VBFH_VtxResolution_Y");
  plotMacro("VtxZRes",  "ReFitVtxZRes", -0.02, 0.02, 0., 1.0, "Resolution Z-position", "Arbitrary Units", "VBFH_VtxResolution_Z");
  */

  plotMacro("VtxXErr",  "ReFitVtxXErr", 0.0, 0.02, 0., 1.0, "Error X-position", "Arbitrary Units", "VBFH_VtxError_X", "VBF_HToTauTau_M-125_MC_v4_noBS_anal.root");
  plotMacro("VtxYErr",  "ReFitVtxYErr", 0.0, 0.02, 0., 1.0, "Error Y-position", "Arbitrary Units", "VBFH_VtxError_Y", "VBF_HToTauTau_M-125_MC_v4_noBS_anal.root");
  plotMacro("VtxZErr",  "ReFitVtxZErr", 0.0, 0.02, 0., 1.0, "Error Z-position", "Arbitrary Units", "VBFH_VtxError_Z", "VBF_HToTauTau_M-125_MC_v4_noBS_anal.root");

  plotMacro("VtxXErr",  "ReFitVtxXErr", 0.0, 0.02, 0., 1.0, "Error X-position", "Arbitrary Units", "GGFH_VtxError_X", "GluGluToHToTauTau_M-125_MC_v4_noBS_anal.root");
  plotMacro("VtxYErr",  "ReFitVtxYErr", 0.0, 0.02, 0., 1.0, "Error Y-position", "Arbitrary Units", "GGFH_VtxError_Y", "GluGluToHToTauTau_M-125_MC_v4_noBS_anal.root");
  plotMacro("VtxZErr",  "ReFitVtxZErr", 0.0, 0.02, 0., 1.0, "Error Z-position", "Arbitrary Units", "GGFH_VtxError_Z", "GluGluToHToTauTau_M-125_MC_v4_noBS_anal.root");


  plotMacro("DeltaPVwrtGenX",  "DeltaRfPVwrtGenX", -0.01, 0.01, 0., 1.0, "Gen X-position - Reco X-position", "Arbitrary Units", "VBFH_PosDiffWrtGen_X", "VBF_HToTauTau_M-125_MC_v4_noBS_anal.root");
  plotMacro("DeltaPVwrtGenY",  "DeltaRfPVwrtGenY", -0.01, 0.01, 0., 1.0, "Gen Y-position - Reco Y-position", "Arbitrary Units", "VBFH_PosDiffWrtGen_Y", "VBF_HToTauTau_M-125_MC_v4_noBS_anal.root");
  plotMacro("DeltaPVwrtGenZ",  "DeltaRfPVwrtGenZ", -0.01, 0.01, 0., 1.0, "Gen Z-position - Reco Z-position", "Arbitrary Units", "VBFH_PosDiffWrtGen_Z", "VBF_HToTauTau_M-125_MC_v4_noBS_anal.root");

  plotMacro("DeltaPVwrtGenX",  "DeltaRfPVwrtGenX", -0.01, 0.01, 0., 1.0, "Gen X-position - Reco X-position", "Arbitrary Units", "GGFH_PosDiffWrtGen_X", "GluGluToHToTauTau_M-125_MC_v4_noBS_anal.root");
  plotMacro("DeltaPVwrtGenY",  "DeltaRfPVwrtGenY", -0.01, 0.01, 0., 1.0, "Gen Y-position - Reco Y-position", "Arbitrary Units", "GGFH_PosDiffWrtGen_Y", "GluGluToHToTauTau_M-125_MC_v4_noBS_anal.root");
  plotMacro("DeltaPVwrtGenZ",  "DeltaRfPVwrtGenZ", -0.01, 0.01, 0., 1.0, "Gen Z-position - Reco Z-position", "Arbitrary Units", "GGFH_PosDiffWrtGen_Z", "GluGluToHToTauTau_M-125_MC_v4_noBS_anal.root");

  plotMacro("ResPVwrtGenX",  "ResRfPVwrtGenX", -5.0, 5.0, 0., 1.0, "(Gen X - Reco X) / Reco X Error", "Arbitrary Units", "VBFH_PosResWrtGen_X", "VBF_HToTauTau_M-125_MC_v4_noBS_anal.root");
  plotMacro("ResPVwrtGenY",  "ResRfPVwrtGenY", -5.0, 5.0, 0., 1.0, "(Gen Y - Reco Y) / Reco Y Error", "Arbitrary Units", "VBFH_PosResWrtGen_Y", "VBF_HToTauTau_M-125_MC_v4_noBS_anal.root");
  plotMacro("ResPVwrtGenZ",  "ResRfPVwrtGenZ", -5.0, 5.0, 0., 1.0, "(Gen Z - Reco Z) / Reco Z Error", "Arbitrary Units", "VBFH_PosResWrtGen_Z", "VBF_HToTauTau_M-125_MC_v4_noBS_anal.root");

  plotMacro("ResPVwrtGenX",  "ResRfPVwrtGenX", -5.0, 5.0, 0., 1.0, "(Gen X - Reco X) / Reco X Error", "Arbitrary Units", "GGFH_PosResWrtGen_X", "GluGluToHToTauTau_M-125_MC_v4_noBS_anal.root");
  plotMacro("ResPVwrtGenY",  "ResRfPVwrtGenY", -5.0, 5.0, 0., 1.0, "(Gen Y - Reco Y) / Reco Y Error", "Arbitrary Units", "GGFH_PosResWrtGen_Y", "GluGluToHToTauTau_M-125_MC_v4_noBS_anal.root");
  plotMacro("ResPVwrtGenZ",  "ResRfPVwrtGenZ", -5.0, 5.0, 0., 1.0, "(Gen Z - Reco Z) / Reco Z Error", "Arbitrary Units", "GGFH_PosResWrtGen_Z", "GluGluToHToTauTau_M-125_MC_v4_noBS_anal.root");

}

void plotMacroSample(TString histname_, float xmin_, float xmax_, float ymin_, float ymax_, TString xtitle_, TString ytitle_, TString outFileName_)
{

  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  gStyle->SetOptStat(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleH(0.07);
  gStyle->SetTitleFontSize(0.1);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleOffset(1.3,"y");

  TFile *f1 = TFile::Open("VBF_HToTauTau_M-125_MC_v4_noBS_anal.root", "READONLY");
  TH1F* hPVvbf = (TH1F*)f1->Get(histname_);
  hPVvbf->GetXaxis()->SetRangeUser(xmin_, xmax_);
  //hPVvbf->GetYaxis()->SetRangeUser(ymin_, ymax_);
  hPVvbf->SetLineColor(kRed);
  hPVvbf->GetXaxis()->SetTitle(xtitle_);
  hPVvbf->GetYaxis()->SetTitle(ytitle_);

  TFile *f2 = TFile::Open("GluGluToHToTauTau_M-125_MC_v4_noBS_anal.root", "READONLY");
  TH1F* hPVggf = (TH1F*)f2->Get(histname_);
  hPVggf->SetLineColor(kBlue);

  hPVvbf->DrawNormalized();
  hPVggf->DrawNormalized("same");

  TLegend* leg = new TLegend(0.6,0.75,0.8,0.88,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.03);
  leg->SetHeader("#splitline{CMS Simulation}{ #sqrt{s}=8 TeV}");

  leg->AddEntry(hPVvbf,"VBF H#rightarrow#tau#tau");
  leg->AddEntry(hPVggf,"GGF H#rightarrow#tau#tau");

  leg->Draw();

  c1->SaveAs("plots/"+outFileName_+".png");
}

void plotMacroSampleAll()
{
  plotMacroSample("DeltaPVX", -0.02, 0.02, 0., 1.0, "Difference X-position", "Arbitrary Units", "Position_Difference_X");
  plotMacroSample("DeltaPVY", -0.02, 0.02, 0., 1.0, "Difference Y-position", "Arbitrary Units", "Position_Difference_Y");
  plotMacroSample("DeltaPVZ", -0.02, 0.02, 0., 1.0, "Difference Z-position", "Arbitrary Units", "Position_Difference_Z");
}
  
void plotMacroPhi(TString histname_, float xmin_, float xmax_, float ymin_, float ymax_, TString xtitle_, TString ytitle_, TString outFileName_)
{

  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  gStyle->SetOptStat(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleH(0.07);
  gStyle->SetTitleFontSize(0.1);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleOffset(1.3,"y");

  TFile *f1 = TFile::Open("GluGluToHToTauTau_M-125_MC_v5_vtxWithBS_anal_pt10.root", "READONLY");
  TH1F* hggf = (TH1F*)f1->Get(histname_);
  hggf->GetXaxis()->SetRangeUser(xmin_, xmax_);
  //hggf->GetYaxis()->SetRangeUser(ymin_, ymax_);
  hggf->SetMaximum(hggf->GetMaximum()*2.0);
  hggf->SetLineColor(kRed);
  hggf->GetXaxis()->SetTitle(xtitle_);
  hggf->GetYaxis()->SetTitle(ytitle_);

  TFile *f2 = TFile::Open("SUSYGluGluToHToTauTau_M-120_MC_v5_vtxWithBS_anal_pt10.root", "READONLY");
  TH1F* hSusyggf = (TH1F*)f2->Get(histname_);
  hSusyggf->SetLineColor(kBlue);

  TFile *f3 = TFile::Open("DYJetsToLL_M-50_MC_v5_vtxWithBS_anal_pt10.root", "READONLY");
  TH1F* hdy = (TH1F*)f3->Get(histname_);
  hdy->SetLineColor(kGreen);

  hggf->DrawNormalized();
  hSusyggf->DrawNormalized("same");
  hdy->DrawNormalized("same");

  TLegend* leg = new TLegend(0.6,0.7,0.8,0.88,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.03);
  leg->SetHeader("#splitline{CMS Simulation}{ #sqrt{s}=8 TeV}");

  leg->AddEntry(hggf,"SM GGF H#rightarrow#tau#tau");
  leg->AddEntry(hSusyggf,"MSSM GGF H#rightarrow#tau#tau");
  leg->AddEntry(hdy,"Z#rightarrow#tau#tau");

  leg->Draw();

  c1->SaveAs("plots/"+outFileName_+".png");
}

void plotMacroPhiAll()
{
  plotMacroPhi("CPPhiStar", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_pt10");
  plotMacroPhi("CPPhiStar_opv", 0., 3.2, 0., 0.15, "#phi^{*} (Radian)", "Arbitrary Units", "cp_phi_star_opv_pt10");
  plotMacroPhi("CPPhiLab", 0., 3.2, 0., 0.15, "#phi^{lab} (Radian)", "Arbitrary Units", "cp_phi_lab_pt10");
}
