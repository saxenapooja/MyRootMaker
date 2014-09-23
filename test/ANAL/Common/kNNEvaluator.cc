
#include <cassert>
#include <fstream>

#include <TFile.h>

#include "Common/kNNEvaluator.h"

SinglekNNEvaluator::SinglekNNEvaluator(const char* weightsFile)
{
  // If the MVA couldn't be booked than there is no trained MVA yet. Don't bail though --
  // we'll keep running and create trees to be used to train the MVA. In a second run we
  // can then apply the MVA.
  reader = new TMVA::Reader("Silent");

  // Register the variables with the reader
  reader->AddVariable("LeptonPt", &leptonPt_);
  reader->AddVariable("JetPt", &jetPt_);
  reader->AddVariable("NJets", &nJets20_);
 
  try
  {
    std::cout << "Loading kNN training from " << weightsFile << "..." << std::endl;

    std::ifstream stream(weightsFile);
    if(!stream)
      throw std::runtime_error("Weight file does not exist");

    // Note BookMVA crashes the program if it fails...
    reader->BookMVA("kNN", weightsFile);
  }
  catch(const std::exception& ex)
  {
    delete reader;
    reader = NULL;

    throw ex;
  }
}

double SinglekNNEvaluator::mvaValue(double leptonPt, double jetPt, unsigned int nJets20) const
{
  leptonPt_ = leptonPt;
  jetPt_ = jetPt;
  nJets20_ = nJets20;
  const double value = reader->EvaluateMVA("kNN");
  return value;
}

std::string kNNTrainingEra(const char* data_era)
{
  if(strcmp(data_era, "2011") == 0 || strcmp(data_era, "fall11") == 0 || strcmp(data_era, "fall-1jet") == 0 || strcmp(data_era, "fall11-2jets") == 0 || strcmp(data_era, "fall11-3jets") == 0 || strcmp(data_era, "fall11-4jets") == 0)
    return "2011";
  else if(strcmp(data_era, "2012") == 0 || strcmp(data_era, "summer12") == 0 || strcmp(data_era, "fall-1jet") == 0 || strcmp(data_era, "summer12-2jets") == 0 || strcmp(data_era, "summer12-3jets") == 0 || strcmp(data_era, "summer12-4jets") == 0)
    return "2012";
  else 
    throw std::runtime_error("No such era: " + std::string(data_era));
}

std::string weightsSuffix(bool wz_zz)
{
  if(wz_zz) return "";
  else return "-NoWZZZ";
}

kNNEvaluator::kNNEvaluator(const char* data_era, const char* version, const char* iso, const char* nkNN, const char* sample, const char* region):
  leading_muon(("/scratch/hh/dust/naf/cms/user/aburgmei/tautau-fullyhadronic/WH_Analysis/WH_Analysis/kNNTraining/weights-" + std::string(sample) + "-" + std::string(version) + "-" + kNNTrainingEra(data_era) + "-lead-" + std::string(region) + "-" + std::string(iso) + "-" + std::string(nkNN) + "/kNN_kNN.weights.xml").c_str()),
  subleading_muon(("/scratch/hh/dust/naf/cms/user/aburgmei/tautau-fullyhadronic/WH_Analysis/WH_Analysis/kNNTraining/weights-" + std::string(sample) + "-" + std::string(version) + "-" + kNNTrainingEra(data_era) + "-sub-" + std::string(region) + "-" + std::string(iso) + "-" + std::string(nkNN) + "/kNN_kNN.weights.xml").c_str())
{
}

double kNNEvaluator::getLeadingMuon(double lepPt, double jetPt, unsigned int nJets20) const
{
  return leading_muon.mvaValue(lepPt, jetPt, nJets20);
}

double kNNEvaluator::getSubleadingMuon(double lepPt, double jetPt, unsigned int nJets20) const
{
  return subleading_muon.mvaValue(lepPt, jetPt, nJets20);
}

#if 0
FRkNNEvaluator::FRkNNEvaluator(const char* data_era, const char* version, const char* iso, const char* variable, bool wzzz)
{
	file = new TFile(("/scratch/hh/dust/naf/cms/user/aburgmei/tautau-fullyhadronic/WH_Analysis/WH_Analysis/kNNTraining/" + kNNTrainingEra(data_era) + ".root").c_str());
	qcd_leading_muon = dynamic_cast<TF1*>(file->Get((std::string(version) + "/" + std::string(iso) + "/" + std::string(variable) + "/ExpFit_QCDlead" + weightsSuffix(wzzz)).c_str()));
	qcd_subleading_muon = dynamic_cast<TF1*>(file->Get((std::string(version) + "/" + std::string(iso) + "/" + std::string(variable) + "/ExpFit_QCDsub" + weightsSuffix(wzzz)).c_str()));
	wjets_leading_muon = dynamic_cast<TF1*>(file->Get((std::string(version) + "/" + std::string(iso) + "/" + std::string(variable) + "/ExpFit_WJetslead" + weightsSuffix(wzzz)).c_str()));
	wjets_subleading_muon = dynamic_cast<TF1*>(file->Get((std::string(version) + "/" + std::string(iso) + "/" + std::string(variable) + "/ExpFit_WJetssub" + weightsSuffix(wzzz)).c_str()));

	if(strcmp(variable, "Pt") == 0) variable_ = LEP_PT;
	else if(strcmp(variable, "JetPt") == 0) variable_ = JET_PT;
	else if(strcmp(variable, "JetPt_000_020") == 0) variable_ = JET_PT_000_020;
	else if(strcmp(variable, "JetPt_020_040") == 0) variable_ = JET_PT_020_040;
	else if(strcmp(variable, "JetPt_040_inf") == 0) variable_ = JET_PT_040_INF;
	else throw std::runtime_error("Unknown variable: " + std::string(variable));

	if((variable_ != JET_PT_000_020 && !qcd_leading_muon) || !qcd_subleading_muon || (variable_ != JET_PT_000_020 && !wjets_leading_muon) || !wjets_subleading_muon)
		throw std::runtime_error("Could not load kNN fake rates for variable" + std::string(variable));
}

double FRkNNEvaluator::getQCDLeadingMuon(double lepPt, double jetPt, unsigned int nJets20) const
{
  switch(variable_)
  {
  case LEP_PT: return qcd_leading_muon->Eval(lepPt);
  case JET_PT: return qcd_leading_muon->Eval(jetPt);
  case JET_PT_000_020: assert(lepPt >= 0. && lepPt <= 20.); return qcd_leading_muon->Eval(jetPt);
  case JET_PT_020_040: assert(lepPt >= 20. && lepPt <= 40.); return qcd_leading_muon->Eval(jetPt);
  case JET_PT_040_INF: assert(lepPt >= 40); return qcd_leading_muon->Eval(jetPt);
  default: assert(false); return 0.0;
  }
}

double FRkNNEvaluator::getQCDSubleadingMuon(double lepPt, double jetPt, unsigned int nJets20) const
{
  switch(variable_)
  {
  case LEP_PT: return qcd_subleading_muon->Eval(lepPt);
  case JET_PT: return qcd_subleading_muon->Eval(jetPt);
  case JET_PT_000_020: assert(lepPt >= 0. && lepPt <= 20.); return qcd_subleading_muon->Eval(jetPt);
  case JET_PT_020_040: assert(lepPt >= 20. && lepPt <= 40.); return qcd_subleading_muon->Eval(jetPt);
  case JET_PT_040_INF: assert(lepPt >= 40); return qcd_subleading_muon->Eval(jetPt);
  default: assert(false); return 0.0;
  }
}

double FRkNNEvaluator::getWJetsLeadingMuon(double lepPt, double jetPt, unsigned int nJets20) const
{
  switch(variable_)
  {
  case LEP_PT: return wjets_leading_muon->Eval(lepPt);
  case JET_PT: return wjets_leading_muon->Eval(jetPt);
  case JET_PT_000_020: assert(lepPt >= 0. && lepPt <= 20.); return wjets_leading_muon->Eval(jetPt);
  case JET_PT_020_040: assert(lepPt >= 20. && lepPt <= 40.); return wjets_leading_muon->Eval(jetPt);
  case JET_PT_040_INF: assert(lepPt >= 40); return wjets_leading_muon->Eval(jetPt);
  default: assert(false); return 0.0;
  }
}

double FRkNNEvaluator::getWJetsSubleadingMuon(double lepPt, double jetPt, unsigned int nJets20) const
{
  switch(variable_)
  {
  case LEP_PT: return wjets_subleading_muon->Eval(lepPt);
  case JET_PT: return wjets_subleading_muon->Eval(jetPt);
  case JET_PT_000_020: assert(lepPt >= 0. && lepPt <= 20.); return wjets_subleading_muon->Eval(jetPt);
  case JET_PT_020_040: assert(lepPt >= 20. && lepPt <= 40.); return wjets_subleading_muon->Eval(jetPt);
  case JET_PT_040_INF: assert(lepPt >= 40); return wjets_subleading_muon->Eval(jetPt);
  default: assert(false); return 0.0;
  }
}

FRZkNNEvaluator::FRZkNNEvaluator(const char* data_era, const char* version, const char* iso, bool wzzz)
{
	const char* variable = "JetPt";

	file = new TFile(("/scratch/hh/dust/naf/cms/user/aburgmei/tautau-fullyhadronic/WH_Analysis/WH_Analysis/kNNTraining/" + kNNTrainingEra(data_era) + ".root").c_str());
	zjets_leading_muon = dynamic_cast<TF1*>(file->Get((std::string(version) + "/" + std::string(iso) + "/" + std::string(variable) + "/ExpFit_ZJetslead" + weightsSuffix(wzzz)).c_str()));
	zjets_subleading_muon = dynamic_cast<TF1*>(file->Get((std::string(version) + "/" + std::string(iso) + "/" + std::string(variable) + "/ExpFit_ZJetssub" + weightsSuffix(wzzz)).c_str()));

	if(!zjets_leading_muon || !zjets_subleading_muon)
		throw std::runtime_error("Could not load kNN fake rates for variable" + std::string(variable));
}

double FRZkNNEvaluator::getQCDLeadingMuon(double lepPt, double jetPt, unsigned int nJets20) const
{
  return zjets_leading_muon->Eval(jetPt);
}

double FRZkNNEvaluator::getQCDSubleadingMuon(double lepPt, double jetPt, unsigned int nJets20) const
{
  return zjets_subleading_muon->Eval(jetPt);
}

double FRZkNNEvaluator::getWJetsLeadingMuon(double lepPt, double jetPt, unsigned int nJets20) const
{
  return zjets_leading_muon->Eval(jetPt);
}

double FRZkNNEvaluator::getWJetsSubleadingMuon(double lepPt, double jetPt, unsigned int nJets20) const
{
  return zjets_subleading_muon->Eval(jetPt);
}

FR2DkNNEvaluator::FR2DkNNEvaluator(const char* data_era, const char* version, const char* iso, bool wzzz):
  fr000_020(data_era, version, iso, "JetPt_000_020", wzzz),
  fr020_040(data_era, version, iso, "JetPt_020_040", wzzz),
  fr040_inf(data_era, version, iso, "JetPt_040_inf", wzzz)
{
}

double FR2DkNNEvaluator::getQCDLeadingMuon(double lepPt, double jetPt, unsigned int nJets20) const
{
  if(lepPt >= 0. && lepPt <= 20.) return fr000_020.getQCDLeadingMuon(lepPt, jetPt, nJets20);
  else if(lepPt >= 20. && lepPt <= 40.) return fr020_040.getQCDLeadingMuon(lepPt, jetPt, nJets20);
  else return fr040_inf.getQCDLeadingMuon(lepPt, jetPt, nJets20);
}

double FR2DkNNEvaluator::getQCDSubleadingMuon(double lepPt, double jetPt, unsigned int nJets20) const
{
  if(lepPt >= 0. && lepPt <= 20.) return fr000_020.getQCDSubleadingMuon(lepPt, jetPt, nJets20);
  else if(lepPt >= 20. && lepPt <= 40.) return fr020_040.getQCDSubleadingMuon(lepPt, jetPt, nJets20);
  else return fr040_inf.getQCDSubleadingMuon(lepPt, jetPt, nJets20);
}

double FR2DkNNEvaluator::getWJetsLeadingMuon(double lepPt, double jetPt, unsigned int nJets20) const
{
  if(lepPt >= 0. && lepPt <= 20.) return fr000_020.getWJetsLeadingMuon(lepPt, jetPt, nJets20);
  else if(lepPt >= 20. && lepPt <= 40.) return fr020_040.getWJetsLeadingMuon(lepPt, jetPt, nJets20);
  else return fr040_inf.getWJetsLeadingMuon(lepPt, jetPt, nJets20);
}

double FR2DkNNEvaluator::getWJetsSubleadingMuon(double lepPt, double jetPt, unsigned int nJets20) const
{
  if(lepPt >= 0. && lepPt <= 20.) return fr000_020.getWJetsSubleadingMuon(lepPt, jetPt, nJets20);
  else if(lepPt >= 20. && lepPt <= 40.) return fr020_040.getWJetsSubleadingMuon(lepPt, jetPt, nJets20);
  else return fr040_inf.getWJetsSubleadingMuon(lepPt, jetPt, nJets20);
}

MaurokNNEvaluator::MaurokNNEvaluator(const char* data_era)
{
	if(kNNTrainingEra(data_era) != "2012")
		throw std::runtime_error("Mauro kNN training only available for 2012 data");

	TMVA::Reader** readers[] = { &qcd_leading_muon, &qcd_subleading_muon, &wjets_leading_muon, &wjets_subleading_muon };
	const char* files[] = { 
		"/scratch/hh/dust/naf/cms/user/aburgmei/tautau-fullyhadronic/WH_Analysis/WH_Analysis/kNNTraining/mauro/mm_qcd_pt10_h2taucuts_muonInfo_k100.kNN.weights.xml",
		"/scratch/hh/dust/naf/cms/user/aburgmei/tautau-fullyhadronic/WH_Analysis/WH_Analysis/kNNTraining/mauro/mm_qcd_pt10_h2taucuts020_muonInfo_k100.kNN.weights.xml",
		"/scratch/hh/dust/naf/cms/user/aburgmei/tautau-fullyhadronic/WH_Analysis/WH_Analysis/kNNTraining/mauro/mm_wjets_pt10_h2taucuts_muonInfo_k100.kNN.weights.xml",
		"/scratch/hh/dust/naf/cms/user/aburgmei/tautau-fullyhadronic/WH_Analysis/WH_Analysis/kNNTraining/mauro/mm_wjets_pt10_h2taucuts020_muonInfo_k100.kNN.weights.xml"
	};

	for(unsigned int i = 0; i < 4; ++i)
	{
		TMVA::Reader** reader = readers[i];
		*reader = new TMVA::Reader("Silent");

		// Register the variables with the reader
		(*reader)->AddVariable("muonJetPt", &jetPt_);
		(*reader)->AddVariable("muonPt", &leptonPt_);
		(*reader)->AddVariable("numJets20", &nJets20_);
 
		(*reader)->BookMVA("kNN", files[i]);
		std::cout << "Loaded Mauro kNN training from " << files[i] << std::endl;
	}
}

double MaurokNNEvaluator::getQCDLeadingMuon(double lepPt, double jetPt, unsigned int nJets20) const
{
	jetPt_ = jetPt;
	leptonPt_ = lepPt;
	nJets20_ = nJets20;
	const double value = qcd_leading_muon->EvaluateMVA("kNN");
	return value;
}

double MaurokNNEvaluator::getQCDSubleadingMuon(double lepPt, double jetPt, unsigned int nJets20) const
{
	jetPt_ = jetPt;
	leptonPt_ = lepPt;
	nJets20_ = nJets20;
	const double value = qcd_subleading_muon->EvaluateMVA("kNN");
	return value;
}

double MaurokNNEvaluator::getWJetsLeadingMuon(double lepPt, double jetPt, unsigned int nJets20) const
{
	jetPt_ = jetPt;
	leptonPt_ = lepPt;
	nJets20_ = nJets20;
	const double value = wjets_leading_muon->EvaluateMVA("kNN");
	return value;
}

double MaurokNNEvaluator::getWJetsSubleadingMuon(double lepPt, double jetPt, unsigned int nJets20) const
{
	jetPt_ = jetPt;
	leptonPt_ = lepPt;
	nJets20_ = nJets20;
	const double value = wjets_subleading_muon->EvaluateMVA("kNN");
	return value;
}
#endif
