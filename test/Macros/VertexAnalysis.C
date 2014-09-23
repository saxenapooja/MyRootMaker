//--------------------------------------------------------------------------------------------------
// VertexAnalysis
//
// Macro to study Vertex parameters for Higgs CP-studies
// Authors: A. Nayak
//--------------------------------------------------------------------------------------------------
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "TString.h"
#include "TPad.h"
#include <TFile.h>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TAxis.h>
#include <TGraphAsymmErrors.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TMath.h>
#include <TROOT.h>
#include <TStyle.h>
#include "TCut.h"
#include "Math/PtEtaPhiE4D.h"
#include "Math/LorentzVector.h"
#include "Math/Point3D.h"
#include "Math/GenVector/CoordinateSystemTags.h"
#include "Math/GenVector/Boost.h"
#include "Math/Vector3D.h"

#include <string>
#include <map>
#include <iostream>
#include <iomanip>

#define DEBUG false
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LV;
typedef ROOT::Math::XYZPointD Point3D;
typedef ROOT::Math::XYZVectorD PV;

using namespace std;

TH1F * DrawOverflow(TH1F *h)
{
  // This function paint the histogram h with an extra bin for overflows
  UInt_t nx    = h->GetNbinsX()+1;
  Double_t *xbins= new Double_t[nx+1];
  for (UInt_t i=0;i<nx;i++)
    xbins[i]=h->GetBinLowEdge(i+1);
  xbins[nx]=xbins[nx-1]+h->GetBinWidth(nx);
  char *tempName= new char[strlen(h->GetName())+10];
  sprintf(tempName,"%swtOverFlow",h->GetName());
  // Book a temporary histogram having ab extra bin for overflows
  TH1F *htmp = new TH1F(tempName, h->GetTitle(), nx, xbins);
  // Reset the axis labels
  htmp->SetXTitle(h->GetXaxis()->GetTitle());
  htmp->SetYTitle(h->GetYaxis()->GetTitle());
  // Fill the new hitogram including the extra bin for overflows
  for (UInt_t i=1; i<=nx; i++){
    htmp->Fill(htmp->GetBinCenter(i), h->GetBinContent(i));
    htmp->SetBinError(i, h->GetBinError(i));
  }

  // Fill the underflows
  htmp->Fill(h->GetBinLowEdge(1)-1, h->GetBinContent(0));
  htmp->SetBinError(0, h->GetBinError(0));
  // Restore the number of entries
  htmp->SetEntries(h->GetEntries());

  // FillStyle and color
  htmp->SetLineWidth(h->GetLineWidth());
  htmp->SetLineColor(h->GetLineColor());
  htmp->SetFillStyle(h->GetFillStyle());
  htmp->SetFillColor(h->GetFillColor());

  return htmp;
}

std::map<std::string, TH1F*>  createHistograms(){

  std::map<std::string, TH1F*> histos_;

  histos_["DeltaPVX"] = new TH1F("DeltaPVX", "", 100, -0.1, 0.1);  
  histos_["DeltaPVY"] = new TH1F("DeltaPVY", "", 100, -0.1, 0.1);
  histos_["DeltaPVZ"] = new TH1F("DeltaPVZ", "", 100, -0.1, 0.1);

  histos_["DeltaPVwrtGenX"] = new TH1F("DeltaPVwrtGenX", "", 100, -0.1, 0.1);
  histos_["DeltaPVwrtGenY"] = new TH1F("DeltaPVwrtGenY", "", 100, -0.1, 0.1);
  histos_["DeltaPVwrtGenZ"] = new TH1F("DeltaPVwrtGenZ", "", 100, -0.1, 0.1);

  histos_["DeltaRfPVwrtGenX"] = new TH1F("DeltaRfPVwrtGenX", "", 100, -0.1, 0.1);
  histos_["DeltaRfPVwrtGenY"] = new TH1F("DeltaRfPVwrtGenY", "", 100, -0.1, 0.1);
  histos_["DeltaRfPVwrtGenZ"] = new TH1F("DeltaRfPVwrtGenZ", "", 100, -0.1, 0.1);

  histos_["ResPVwrtGenX"] = new TH1F("ResPVwrtGenX", "", 100, -10.0, 10.0);
  histos_["ResPVwrtGenY"] = new TH1F("ResPVwrtGenY", "", 100, -10.0, 10.0);
  histos_["ResPVwrtGenZ"] = new TH1F("ResPVwrtGenZ", "", 100, -10.0, 10.0);

  histos_["ResRfPVwrtGenX"] = new TH1F("ResRfPVwrtGenX", "", 100, -10.0, 10.0);
  histos_["ResRfPVwrtGenY"] = new TH1F("ResRfPVwrtGenY", "", 100, -10.0, 10.0);
  histos_["ResRfPVwrtGenZ"] = new TH1F("ResRfPVwrtGenZ", "", 100, -10.0, 10.0);

  histos_["VtxXRes"] = new TH1F("VtxXRes",  "", 100, -0.1, 0.1);
  histos_["VtxYRes"] = new TH1F("VtxYRes",  "", 100, -0.1, 0.1);
  histos_["VtxZRes"] = new TH1F("VtxZRes",  "", 100, -0.1, 0.1);

  histos_["ReFitVtxXRes"] = new TH1F("ReFitVtxXRes",  "", 100, -0.1, 0.1);
  histos_["ReFitVtxYRes"] = new TH1F("ReFitVtxYRes",  "", 100, -0.1, 0.1);
  histos_["ReFitVtxZRes"] = new TH1F("ReFitVtxZRes",  "", 100, -0.1, 0.1);
  
  histos_["VtxXErr"] = new TH1F("VtxXErr",  "", 100, -0.05, 0.05);
  histos_["VtxYErr"] = new TH1F("VtxYErr",  "", 100, -0.05, 0.05);
  histos_["VtxZErr"] = new TH1F("VtxZErr",  "", 100, -0.05, 0.05);

  histos_["ReFitVtxXErr"] = new TH1F("ReFitVtxXErr",  "", 100, -0.05, 0.05);
  histos_["ReFitVtxYErr"] = new TH1F("ReFitVtxYErr",  "", 100, -0.05, 0.05);
  histos_["ReFitVtxZErr"] = new TH1F("ReFitVtxZErr",  "", 100, -0.05, 0.05);

  histos_["CPPhiStar"] = new TH1F("CPPhiStar",  "", 20, 0.0, 4.0);
  histos_["CPPhiLab"] = new TH1F("CPPhiLab",  "", 20, 0.0, 4.0);

  histos_["CPPhiStar_opv"] = new TH1F("CPPhiStar_opv",  "", 20, 0.0, 4.0);
  histos_["CPPhiStar_gen"] = new TH1F("CPPhiStar_gen",  "", 20, 0.0, 4.0);
  
  return histos_;
}
  
void VertexAnalysis(TString inputFilePath_, TString outFileName_)
{

  TChain *chain = new TChain("vertexAnalyzer/tree");
  chain->Add(inputFilePath_+"/*.root");

  TFile *outFile_ = TFile::Open(outFileName_+".root", "RECREATE");
  outFile_->SetCompressionLevel( 9 );

  //std::map<std::string, TH1F*> histos_;
  std::map<std::string, TH1F*> histos_ = createHistograms();
  
  unsigned long run_,event_,lumi_;
  //int index_;
  
  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > *diTauLegsP4_ 
    = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > (); 
  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > *diTauLegsLchP4_
    = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > ();
  std::vector< ROOT::Math::XYZPointD > *diTauLegsPCA_ 
    = new std::vector< ROOT::Math::XYZPointD > ();
  std::vector< ROOT::Math::XYZPointD > *diTauLegsPCAOPV_
    = new std::vector< ROOT::Math::XYZPointD > ();
  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* genTausP4_
    = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > ();
  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* genVP4_
    = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > ();
  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* genTauPSonsP4_
    = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > ();
  std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >* genTauNSonsP4_
    = new std::vector< ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > ();

  std::vector< int > *genTausPid_ = new std::vector< int >();
  std::vector< int > *genTausCharge_ = new std::vector< int >();
  std::vector< int > *genTausStatus_ = new std::vector< int >();
  std::vector< int > *genTauPSonsPid_ = new std::vector< int >();
  std::vector< int > *genTauPSonsCharge_ = new std::vector< int >();
  std::vector< int > *genTauPSonsStatus_ = new std::vector< int >();
  std::vector< int > *genTauNSonsPid_ = new std::vector< int >();
  std::vector< int > *genTauNSonsCharge_ = new std::vector< int >();
  std::vector< int > *genTauNSonsStatus_ = new std::vector< int >();

  int genVPid_;
  int tightestHPSDB3HWPLeg1_, tightestHPSDB3HWPLeg2_;
  float hpsDB3HLeg1_, hpsDB3HLeg2_, diTauCharge_, chargeLeg1_;
  int decayModeLeg1_, decayModeLeg2_;

  float VtxX_, VtxY_, VtxZ_, VtxXErr_, VtxYErr_, VtxZErr_;
  float genVtxX_, genVtxY_, genVtxZ_;
  
  float ReFitVtxX_, ReFitVtxY_, ReFitVtxZ_, ReFitVtxXErr_, ReFitVtxYErr_, ReFitVtxZErr_;
  int   ReFitVtxRho_;
  float ReFitVtxNdof_;

  float dxyLeg1_, dxyLeg2_, dxyLeg1Err_, dxyLeg2Err_;
  float dzLeg1_, dzLeg2_, dzLeg1Err_, dzLeg2Err_;
  float dxyOPVLeg1_, dxyOPVLeg2_, dxyOPVLeg1Err_, dxyOPVLeg2Err_;
  float dzOPVLeg1_, dzOPVLeg2_, dzOPVLeg1Err_, dzOPVLeg2Err_;


  //Set branches
  chain->SetBranchStatus("run", 1);
  chain->SetBranchStatus("event", 1);
  chain->SetBranchStatus("lumi", 1);
  chain->SetBranchStatus("diTauLegsP4", 1);
  chain->SetBranchStatus("diTauLegsLchP4", 1);
  chain->SetBranchStatus("diTauLegsPCA", 1);
  chain->SetBranchStatus("diTauLegsPCAOPV", 1);
  chain->SetBranchStatus("genTausP4", 1);
  chain->SetBranchStatus("genVP4", 1);
  chain->SetBranchStatus("genVPid", 1);
  chain->SetBranchStatus("genTausPid", 1);
  chain->SetBranchStatus("genTausCharge", 1);
  chain->SetBranchStatus("genTausStatus", 1);
  chain->SetBranchStatus("genTauPSonsP4", 1);
  chain->SetBranchStatus("genTauNSonsP4", 1);
  chain->SetBranchStatus("genTauPSonsPid", 1);
  chain->SetBranchStatus("genTauPSonsCharge", 1);
  chain->SetBranchStatus("genTauPSonsStatus", 1);
  chain->SetBranchStatus("genTauNSonsPid", 1);
  chain->SetBranchStatus("genTauNSonsCharge", 1);
  chain->SetBranchStatus("genTauNSonsStatus", 1);
  chain->SetBranchStatus("tightestHPSDB3HWPLeg1", 1);
  chain->SetBranchStatus("tightestHPSDB3HWPLeg2", 1);
  chain->SetBranchStatus("hpsDB3HLeg1", 1);
  chain->SetBranchStatus("hpsDB3HLeg2", 1);
  chain->SetBranchStatus("diTauCharge", 1);
  chain->SetBranchStatus("chargeLeg1", 1);
  chain->SetBranchStatus("decayModeLeg1", 1);
  chain->SetBranchStatus("decayModeLeg2", 1);
  chain->SetBranchStatus("genVtxX", 1);
  chain->SetBranchStatus("genVtxY", 1);
  chain->SetBranchStatus("genVtxZ", 1);
  chain->SetBranchStatus("VtxX", 1);
  chain->SetBranchStatus("VtxY", 1);
  chain->SetBranchStatus("VtxZ", 1);
  chain->SetBranchStatus("VtxXErr", 1);
  chain->SetBranchStatus("VtxYErr", 1);
  chain->SetBranchStatus("VtxZErr", 1);
  chain->SetBranchStatus("ReFitVtxX", 1);
  chain->SetBranchStatus("ReFitVtxY", 1);
  chain->SetBranchStatus("ReFitVtxZ", 1);
  chain->SetBranchStatus("ReFitVtxXErr", 1);
  chain->SetBranchStatus("ReFitVtxYErr", 1);
  chain->SetBranchStatus("ReFitVtxZErr", 1);
  chain->SetBranchStatus("ReFitVtxRho", 1);
  chain->SetBranchStatus("ReFitVtxNdof", 1);
  chain->SetBranchStatus("dxyLeg1", 1);
  chain->SetBranchStatus("dxyLeg2", 1);
  chain->SetBranchStatus("dxyLeg1Err", 1);
  chain->SetBranchStatus("dxyLeg2Err", 1);
  chain->SetBranchStatus("dzLeg1", 1);
  chain->SetBranchStatus("dzLeg2", 1);
  chain->SetBranchStatus("dzLeg1Err", 1);
  chain->SetBranchStatus("dzLeg2Err", 1);
  chain->SetBranchStatus("dxyOPVLeg1", 1);
  chain->SetBranchStatus("dxyOPVLeg2", 1);
  chain->SetBranchStatus("dxyOPVLeg1Err", 1);
  chain->SetBranchStatus("dxyOPVLeg2Err", 1);
  chain->SetBranchStatus("dzOPVLeg1", 1);
  chain->SetBranchStatus("dzOPVLeg2", 1);
  chain->SetBranchStatus("dzOPVLeg1Err", 1);
  chain->SetBranchStatus("dzOPVLeg2Err", 1);

  chain->SetBranchAddress("run", &run_);
  chain->SetBranchAddress("event", &event_);
  chain->SetBranchAddress("lumi", &lumi_);
  chain->SetBranchAddress("diTauLegsP4", &diTauLegsP4_);
  chain->SetBranchAddress("diTauLegsLchP4", &diTauLegsLchP4_);
  chain->SetBranchAddress("diTauLegsPCA", &diTauLegsPCA_);
  chain->SetBranchAddress("diTauLegsPCAOPV", &diTauLegsPCAOPV_);
  chain->SetBranchAddress("genTausP4", &genTausP4_);
  chain->SetBranchAddress("genVP4", &genVP4_);
  chain->SetBranchAddress("genVPid", &genVPid_);
  chain->SetBranchAddress("genTausPid", &genTausPid_);
  chain->SetBranchAddress("genTausCharge", &genTausCharge_);
  chain->SetBranchAddress("genTausStatus", &genTausStatus_);
  chain->SetBranchAddress("genTauPSonsP4", &genTauPSonsP4_);
  chain->SetBranchAddress("genTauNSonsP4", &genTauNSonsP4_);
  chain->SetBranchAddress("genTauPSonsPid", &genTauPSonsPid_);
  chain->SetBranchAddress("genTauPSonsCharge", &genTauPSonsCharge_);
  chain->SetBranchAddress("genTauPSonsStatus", &genTauPSonsStatus_);
  chain->SetBranchAddress("genTauNSonsPid", &genTauNSonsPid_);
  chain->SetBranchAddress("genTauNSonsCharge", &genTauNSonsCharge_);
  chain->SetBranchAddress("genTauNSonsStatus", &genTauNSonsStatus_);
  chain->SetBranchAddress("tightestHPSDB3HWPLeg1", &tightestHPSDB3HWPLeg1_);
  chain->SetBranchAddress("tightestHPSDB3HWPLeg2", &tightestHPSDB3HWPLeg2_);
  chain->SetBranchAddress("hpsDB3HLeg1", &hpsDB3HLeg1_);
  chain->SetBranchAddress("hpsDB3HLeg2", &hpsDB3HLeg2_);
  chain->SetBranchAddress("diTauCharge", &diTauCharge_);
  chain->SetBranchAddress("chargeLeg1", &chargeLeg1_);
  chain->SetBranchAddress("decayModeLeg1", &decayModeLeg1_);
  chain->SetBranchAddress("decayModeLeg2", &decayModeLeg2_);
  chain->SetBranchAddress("genVtxX", &genVtxX_);
  chain->SetBranchAddress("genVtxY", &genVtxY_);
  chain->SetBranchAddress("genVtxZ", &genVtxZ_);
  chain->SetBranchAddress("VtxX", &VtxX_);
  chain->SetBranchAddress("VtxY", &VtxY_);
  chain->SetBranchAddress("VtxZ", &VtxZ_);
  chain->SetBranchAddress("VtxXErr", &VtxXErr_);
  chain->SetBranchAddress("VtxYErr", &VtxYErr_);
  chain->SetBranchAddress("VtxZErr", &VtxZErr_);
  chain->SetBranchAddress("ReFitVtxX", &ReFitVtxX_);
  chain->SetBranchAddress("ReFitVtxY", &ReFitVtxY_);
  chain->SetBranchAddress("ReFitVtxZ", &ReFitVtxZ_);
  chain->SetBranchAddress("ReFitVtxXErr", &ReFitVtxXErr_);
  chain->SetBranchAddress("ReFitVtxYErr", &ReFitVtxYErr_);
  chain->SetBranchAddress("ReFitVtxZErr", &ReFitVtxZErr_);
  chain->SetBranchAddress("ReFitVtxRho", &ReFitVtxRho_);
  chain->SetBranchAddress("ReFitVtxNdof", &ReFitVtxNdof_);
  chain->SetBranchAddress("dxyLeg1", &dxyLeg1_);
  chain->SetBranchAddress("dxyLeg2", &dxyLeg2_);
  chain->SetBranchAddress("dxyLeg1Err", &dxyLeg1Err_);
  chain->SetBranchAddress("dxyLeg2Err", &dxyLeg2Err_);
  chain->SetBranchAddress("dzLeg1", &dzLeg1_);
  chain->SetBranchAddress("dzLeg2", &dzLeg2_);
  chain->SetBranchAddress("dzLeg1Err", &dzLeg1Err_);
  chain->SetBranchAddress("dzLeg2Err", &dzLeg2Err_);
  chain->SetBranchAddress("dxyOPVLeg1", &dxyOPVLeg1_);
  chain->SetBranchAddress("dxyOPVLeg2", &dxyOPVLeg2_);
  chain->SetBranchAddress("dxyOPVLeg1Err", &dxyOPVLeg1Err_);
  chain->SetBranchAddress("dxyOPVLeg2Err", &dxyOPVLeg2Err_);
  chain->SetBranchAddress("dzOPVLeg1", &dzOPVLeg1_);
  chain->SetBranchAddress("dzOPVLeg2", &dzOPVLeg2_);
  chain->SetBranchAddress("dzOPVLeg1Err", &dzOPVLeg1Err_);
  chain->SetBranchAddress("dzOPVLeg2Err", &dzOPVLeg2Err_);

  
  ///////////////////////
  // LOOP OVER ENTRIES //
  ///////////////////////

  int nEntries    = chain->GetEntries() ;
  unsigned int lastRun_, lastLumi_, lastEvent_=0;
  bool foundDiTauPair_ = false;
  for (int n = 0; n < nEntries ; n++) {
    
    chain->GetEntry(n);
    if(n%1000==0) std::cout << n <<"/"<<nEntries<< std::endl;

    //compute decay plane in generator level
    if(event_ != lastEvent_){ //take one entry per event
	if(genVPid_ == 25){
	  //rest frame of Higgs
	  ROOT::Math::Boost boost_to_rf_h((*genVP4_)[0].BoostToCM());
	  
	  LV gentauP_lab(0, 0, 0, 0); LV gentauN_lab(0, 0, 0, 0);
	  if(genTausP4_->size() == 2 && genTausCharge_->size() == 2){
	    if((*genTausCharge_)[0] > 0 && (*genTausCharge_)[1] < 0){
	      gentauP_lab = (*genTausP4_)[0];
	      gentauN_lab = (*genTausP4_)[1];
	    }
	    else if((*genTausCharge_)[0] < 0 && (*genTausCharge_)[1] > 0){
	      gentauP_lab = (*genTausP4_)[1];
              gentauN_lab = (*genTausP4_)[0];
            } 
	  }
	  
	  LV gentauP_rf = gentauP_lab.pt()>0 ? boost_to_rf_h(gentauP_lab) : LV(0, 0, 0, 0);
	  LV gentauN_rf = gentauN_lab.pt()>0 ? boost_to_rf_h(gentauN_lab) : LV(0, 0, 0, 0);

	  LV genpionP_lab(0, 0, 0, 0); LV genpionN_lab(0, 0, 0, 0);
	  if(genTauPSonsP4_->size()>0 && genTauPSonsP4_->size()<3 && genTauPSonsPid_->size() > 0){ //take only one prong
	    for(size_t i = 0; i<genTauPSonsP4_->size(); i++){
	      if(abs((*genTauPSonsPid_)[i]) == 211)genpionP_lab = (*genTauPSonsP4_)[i];
	    }
	  }
	  if(genTauNSonsP4_->size()>0 && genTauNSonsP4_->size()<3 && genTauNSonsPid_->size() > 0){ //take only one prong
            for(size_t i = 0; i<genTauNSonsP4_->size(); i++){
              if(abs((*genTauNSonsPid_)[i]) == 211)genpionN_lab = (*genTauNSonsP4_)[i];
            }
          }

	  LV genpionP_rf = genpionP_lab.pt()>0 ? boost_to_rf_h(genpionP_lab) : LV(0, 0, 0, 0);
          LV genpionN_rf = genpionN_lab.pt()>0 ? boost_to_rf_h(genpionN_lab) : LV(0, 0, 0, 0);

	  if(gentauP_rf.pt()>0 && gentauN_rf.pt()>0 && genpionP_rf.pt()>0 && genpionN_rf.pt()>0){
	    //Get only the position component of the vector
	    PV gentauP_rf_3d = gentauP_rf.Vect();
	    PV gentauN_rf_3d = gentauN_rf.Vect();
	    PV genpionP_rf_3d = genpionP_rf.Vect();
            PV genpionN_rf_3d = genpionN_rf.Vect();

	    //Get the normal to the decay plane
	    PV dplane_normal_1 = gentauP_rf_3d.Cross(genpionP_rf_3d);
	    PV dplane_normal_2 = genpionN_rf_3d.Cross(gentauN_rf_3d); //gentauN_rf_3d.Cross(genpionN_rf_3d);
	    PV dplane_normal_1_u = dplane_normal_1/sqrt(dplane_normal_1.mag2());
	    PV dplane_normal_2_u = dplane_normal_2/sqrt(dplane_normal_2.mag2());

	    //Get CP angle in ZMF
	    double phi_star_gen = TMath::ACos(dplane_normal_1_u.Dot(dplane_normal_2_u));

	    histos_["CPPhiStar_gen"]->Fill(phi_star_gen);
	  }
	}
    }//end of Generator studies

    //Take only one pair per events
    //tau pairs are ordered in 1st OS, the sum pT
    if(foundDiTauPair_ && run_ == lastRun_ &&
       lumi_ == lastLumi_ && event_ == lastEvent_)
      continue;

    foundDiTauPair_ = false;
    lastRun_ = run_;
    lastLumi_ = lumi_; 
    lastEvent_ = event_;
    
    bool passDecayMode = false;
    if(decayModeLeg1_ == 0 && decayModeLeg2_ == 0) passDecayMode = true;

    //Apply Medium DeltaBeta3Hits isolation
    bool passIsolation = false;
    if(tightestHPSDB3HWPLeg1_ >= 1 && tightestHPSDB3HWPLeg2_ >= 1)
      passIsolation = true;

    //apply kinematic cuts
    double ptL1     = (*diTauLegsP4_)[0].Pt();
    double ptL2     = (*diTauLegsP4_)[1].Pt();
    double etaL1    = (*diTauLegsP4_)[0].Eta();
    double etaL2    = (*diTauLegsP4_)[1].Eta();
    bool passKineCuts = false;
    if(ptL1 > 30 && ptL2 > 30 && TMath::Abs(etaL1) < 2.1 && TMath::Abs(etaL2) < 2.1)
      passKineCuts = true;

    bool passLeadChHadCut = false;
    if((*diTauLegsLchP4_)[0].Pt() > 10.0 && (*diTauLegsLchP4_)[1].Pt() > 10.0)
      passLeadChHadCut = true;

    //Apply a dZ and dxy cut wrt original vertex
    //Information not available in the current ntuple

    bool passDiTauCharge = false;
    if(diTauCharge_ == 0) passDiTauCharge = true;

    //use this pair only
    if(passDecayMode && passIsolation)
      foundDiTauPair_ = true;

    if(!passDecayMode) continue;
    if(!passIsolation) continue;
    if(!passKineCuts) continue;
    if(!passLeadChHadCut) continue;
    if(!passDiTauCharge) continue;

    //Fill Histograms

    //Vertex Resolution
    histos_["VtxXRes"]->Fill(VtxXErr_/VtxX_);
    histos_["VtxYRes"]->Fill(VtxYErr_/VtxY_);
    histos_["VtxZRes"]->Fill(VtxZErr_/VtxZ_);

    histos_["ReFitVtxXRes"]->Fill(ReFitVtxXErr_/ReFitVtxX_);
    histos_["ReFitVtxYRes"]->Fill(ReFitVtxYErr_/ReFitVtxY_);
    histos_["ReFitVtxZRes"]->Fill(ReFitVtxZErr_/ReFitVtxZ_);
    
    histos_["VtxXErr"]->Fill(VtxXErr_);
    histos_["VtxYErr"]->Fill(VtxYErr_);
    histos_["VtxZErr"]->Fill(VtxZErr_);

    histos_["ReFitVtxXErr"]->Fill(ReFitVtxXErr_);
    histos_["ReFitVtxYErr"]->Fill(ReFitVtxYErr_);
    histos_["ReFitVtxZErr"]->Fill(ReFitVtxZErr_);
    
    histos_["DeltaPVX"]->Fill(VtxX_ - ReFitVtxX_);
    histos_["DeltaPVY"]->Fill(VtxY_ - ReFitVtxY_);
    histos_["DeltaPVZ"]->Fill(VtxZ_ - ReFitVtxZ_);

    histos_["DeltaPVwrtGenX"]->Fill(genVtxX_ - VtxX_);
    histos_["DeltaPVwrtGenY"]->Fill(genVtxY_ - VtxY_);
    histos_["DeltaPVwrtGenZ"]->Fill(genVtxZ_ - VtxZ_);

    histos_["DeltaRfPVwrtGenX"]->Fill(genVtxX_ - ReFitVtxX_);
    histos_["DeltaRfPVwrtGenY"]->Fill(genVtxY_ - ReFitVtxY_);
    histos_["DeltaRfPVwrtGenZ"]->Fill(genVtxZ_ - ReFitVtxZ_);

    histos_["ResPVwrtGenX"]->Fill((genVtxX_ - VtxX_)/VtxXErr_);
    histos_["ResPVwrtGenY"]->Fill((genVtxY_ - VtxY_)/VtxYErr_);
    histos_["ResPVwrtGenZ"]->Fill((genVtxZ_ - VtxZ_)/VtxZErr_);

    histos_["ResRfPVwrtGenX"]->Fill((genVtxX_ - ReFitVtxX_)/ReFitVtxXErr_);
    histos_["ResRfPVwrtGenY"]->Fill((genVtxY_ - ReFitVtxY_)/ReFitVtxYErr_);
    histos_["ResRfPVwrtGenZ"]->Fill((genVtxZ_ - ReFitVtxZ_)/ReFitVtxZErr_);

    //std::cout<<" diTau pT"<<(*diTauLegsP4_)[0].Pt()+(*diTauLegsP4_)[1].Pt()<<std::endl;
    
    //Make Pi+Pi- Zero-momentum-frame
    LV pion1_lab = ((*diTauLegsLchP4_)[0]);
    LV pion2_lab = ((*diTauLegsLchP4_)[1]);
    LV pion_pair_lab = pion1_lab + pion2_lab;
    
    ROOT::Math::Boost boost_to_rf(pion_pair_lab.BoostToCM());
    //ROOT::Math::Boost boost_to_lab(boost_to_rf.Inverse());

    LV pion1_rf = boost_to_rf(pion1_lab);
    LV pion2_rf = boost_to_rf(pion2_lab);

    //Define the normalized impact parameter vectors
    double ipv1_lab_mag = std::sqrt(pow((*diTauLegsPCA_)[0].x() - ReFitVtxX_, 2) + pow((*diTauLegsPCA_)[0].y() - ReFitVtxY_, 2) +
				    pow((*diTauLegsPCA_)[0].z() - ReFitVtxZ_, 2));
    double ipv2_lab_mag = std::sqrt(pow((*diTauLegsPCA_)[1].x() - ReFitVtxX_, 2) + pow((*diTauLegsPCA_)[1].y() - ReFitVtxY_, 2) +
                                    pow((*diTauLegsPCA_)[1].z() - ReFitVtxZ_, 2));
    LV ipv1_lab(((*diTauLegsPCA_)[0].x() - ReFitVtxX_)/ipv1_lab_mag, ((*diTauLegsPCA_)[0].y() - ReFitVtxY_)/ipv1_lab_mag, 
		((*diTauLegsPCA_)[0].z() - ReFitVtxZ_)/ipv1_lab_mag, 0);
    LV ipv2_lab(((*diTauLegsPCA_)[1].x() - ReFitVtxX_)/ipv2_lab_mag, ((*diTauLegsPCA_)[1].y() - ReFitVtxY_)/ipv2_lab_mag, 
		((*diTauLegsPCA_)[1].z() - ReFitVtxZ_)/ipv2_lab_mag, 0);
    
    //boost to ZMF
    LV ipv1_rf = boost_to_rf(ipv1_lab);
    LV ipv2_rf = boost_to_rf(ipv2_lab);
    
    //Get only the position component of the vector
    PV ipv1_rf_3d = ipv1_rf.Vect();
    PV ipv2_rf_3d = ipv2_rf.Vect();

    PV pion1_rf_3d = pion1_rf.Vect();
    PV pion2_rf_3d = pion2_rf.Vect();

    //Get the unit (normalized) vector along pion momentum 
    PV pion1_rf_3d_u = pion1_rf_3d/TMath::Sqrt(pion1_rf_3d.mag2());
    PV pion2_rf_3d_u = pion2_rf_3d/TMath::Sqrt(pion2_rf_3d.mag2());

    //Get the longitudinal component of IP vector parallel to pion momenta
    PV ipv1_rf_3d_l = ipv1_rf_3d.Dot(pion1_rf_3d_u)*pion1_rf_3d_u;
    PV ipv2_rf_3d_l = ipv2_rf_3d.Dot(pion2_rf_3d_u)*pion2_rf_3d_u;

    //Get IP vector normal to pion momenta
    PV ipv1_rf_3d_t = ipv1_rf_3d - ipv1_rf_3d_l;
    PV ipv2_rf_3d_t = ipv2_rf_3d - ipv2_rf_3d_l;

    //Get normalized normal IP vector
    PV ipv1_rf_3d_t_u = ipv1_rf_3d_t/TMath::Sqrt(ipv1_rf_3d_t.mag2());
    PV ipv2_rf_3d_t_u = ipv2_rf_3d_t/TMath::Sqrt(ipv2_rf_3d_t.mag2());

    //Get CP angle in ZMF
    double phi_star = TMath::ACos(ipv1_rf_3d_t_u.Dot(ipv2_rf_3d_t_u));

    histos_["CPPhiStar"]->Fill(phi_star);
    
    //Get CP angle in lab frame
    PV ipv1_lab_3d(((*diTauLegsPCA_)[0].x() - ReFitVtxX_)/ipv1_lab_mag, ((*diTauLegsPCA_)[0].y() - ReFitVtxY_)/ipv1_lab_mag,
		   ((*diTauLegsPCA_)[0].z() - ReFitVtxZ_)/ipv1_lab_mag);
    PV ipv2_lab_3d(((*diTauLegsPCA_)[1].x() - ReFitVtxX_)/ipv2_lab_mag, ((*diTauLegsPCA_)[1].y() - ReFitVtxY_)/ipv2_lab_mag,
		   ((*diTauLegsPCA_)[1].z() - ReFitVtxZ_)/ipv2_lab_mag);

    double phi_lab = TMath::ACos(ipv1_lab_3d.Dot(ipv2_lab_3d));
    histos_["CPPhiLab"]->Fill(phi_lab);

    //Check this with original vertex
    //Define the normalized impact parameter vectors
    double ipv1_lab_mag_opv = std::sqrt(pow((*diTauLegsPCAOPV_)[0].x() - VtxX_, 2) + pow((*diTauLegsPCAOPV_)[0].y() - VtxY_, 2) +
					pow((*diTauLegsPCAOPV_)[0].z() - VtxZ_, 2));
    double ipv2_lab_mag_opv = std::sqrt(pow((*diTauLegsPCAOPV_)[1].x() - VtxX_, 2) + pow((*diTauLegsPCAOPV_)[1].y() - VtxY_, 2) +
					pow((*diTauLegsPCAOPV_)[1].z() - VtxZ_, 2));
    LV ipv1_lab_opv(((*diTauLegsPCAOPV_)[0].x() - VtxX_)/ipv1_lab_mag, ((*diTauLegsPCAOPV_)[0].y() - VtxY_)/ipv1_lab_mag,
		    ((*diTauLegsPCAOPV_)[0].z() - VtxZ_)/ipv1_lab_mag, 0);
    LV ipv2_lab_opv(((*diTauLegsPCAOPV_)[1].x() - VtxX_)/ipv2_lab_mag, ((*diTauLegsPCAOPV_)[1].y() - VtxY_)/ipv2_lab_mag,
		    ((*diTauLegsPCAOPV_)[1].z() - VtxZ_)/ipv2_lab_mag, 0);

    //boost to ZMF
    LV ipv1_opv_rf = boost_to_rf(ipv1_lab_opv);
    LV ipv2_opv_rf = boost_to_rf(ipv2_lab_opv);

    //Get only the position component of the vector
    PV ipv1_opv_rf_3d = ipv1_opv_rf.Vect();
    PV ipv2_opv_rf_3d = ipv2_opv_rf.Vect();

    //Get the longitudinal component of IP vector parallel to pion momenta
    PV ipv1_opv_rf_3d_l = ipv1_opv_rf_3d.Dot(pion1_rf_3d_u)*pion1_rf_3d_u;
    PV ipv2_opv_rf_3d_l = ipv2_opv_rf_3d.Dot(pion2_rf_3d_u)*pion2_rf_3d_u;

    //Get IP vector normal to pion momenta
    PV ipv1_opv_rf_3d_t = ipv1_opv_rf_3d - ipv1_opv_rf_3d_l;
    PV ipv2_opv_rf_3d_t = ipv2_opv_rf_3d - ipv2_opv_rf_3d_l;

    //Get normalized normal IP vector
    PV ipv1_opv_rf_3d_t_u = ipv1_opv_rf_3d_t/TMath::Sqrt(ipv1_opv_rf_3d_t.mag2());
    PV ipv2_opv_rf_3d_t_u = ipv2_opv_rf_3d_t/TMath::Sqrt(ipv2_opv_rf_3d_t.mag2());

    //Get CP angle in ZMF
    double phi_star_opv = TMath::ACos(ipv1_opv_rf_3d_t_u.Dot(ipv2_opv_rf_3d_t_u));

    histos_["CPPhiStar_opv"]->Fill(phi_star_opv);
  }

  outFile_->Write();
  outFile_->Close();
  
  delete diTauLegsPCA_; delete diTauLegsPCAOPV_; delete diTauLegsLchP4_; delete diTauLegsP4_;
}

void VertexAnalysisAll()
{
  
  VertexAnalysis("/nfs/dust/cms/user/anayak/CMS/Ntuple_VertexStudy/GluGluToHToTauTau_M-125_MC_v7_vtxWithBS/", "GluGluToHToTauTau_M-125_MC_v7_vtxWithBS_anal_pt10");
  VertexAnalysis("/nfs/dust/cms/user/anayak/CMS/Ntuple_VertexStudy/VBF_HToTauTau_M-125_MC_v7_vtxWithBS/", "VBF_HToTauTau_M-125_MC_v7_vtxWithBS_anal_pt10");
  //VertexAnalysis("/nfs/dust/cms/user/anayak/CMS/Ntuple_VertexStudy/SUSYGluGluToHToTauTau_M-120_MC_v7_vtxWithBS/", "SUSYGluGluToHToTauTau_M-120_MC_v7_vtxWithBS_anal_pt10");
  //VertexAnalysis("/nfs/dust/cms/user/anayak/CMS/Ntuple_VertexStudy/DYJetsToLL_M-50_MC_v7_vtxWithBS/", "DYJetsToLL_M-50_MC_v7_vtxWithBS_anal_pt10");
} 
