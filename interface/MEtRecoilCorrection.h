#ifndef _HTAUTAU_METRECOILCORRECTION_H_
#define _HTAUTAU_METRECOILCORRECTION_H_

#include "MyRootMaker/MyRootMaker/interface/RecoilCorrector.hh"
#include "TLorentzVector.h"
#include "Math/VectorUtil.h"
//#include "Math/LorentzVector.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

class MEtRecoilCorrection
{
 public:
  MEtRecoilCorrection(string sampleName_);
  ~MEtRecoilCorrection(){ delete corrector_;};

  TLorentzVector RecoilCorrectedMET(TLorentzVector pfMet_, LorentzVector lepLeg1_, LorentzVector lepLeg2_, LorentzVector genBoson_, string sampleName_, int nJets_);

 private:
  RecoilCorrector *corrector_;
};


MEtRecoilCorrection::MEtRecoilCorrection(string sampleName_)
{
  if(sampleName_.find("ZJets") !=std::string::npos) {
    corrector_ = new RecoilCorrector("/nfs/dust/cms/user/pooja/scratch/plot-macro/RecoilCorrector_v7/recoilfits/recoilfit_zmm53X_2012_njet.root");
    corrector_->addMCFile("/nfs/dust/cms/user/pooja/scratch/plot-macro/RecoilCorrector_v7/recoilfits/recoilfit_zmm53X_2012_njet.root");
    corrector_->addDataFile("/nfs/dust/cms/user/pooja/scratch/plot-macro/RecoilCorrector_v7/recoilfits/recoilfit_datamm53X_2012_njet.root");
  }
  else if(sampleName_.find("WJets")!=std::string::npos){
    corrector_ = new RecoilCorrector("/nfs/dust/cms/user/pooja/scratch/plot-macro/RecoilCorrector_v7/recoilfits/recoilfit_wjets53X_20pv_njet.root");
    corrector_->addMCFile("/nfs/dust/cms/user/pooja/scratch/plot-macro/RecoilCorrector_v7/recoilfits/recoilfit_zmm53X_2012_njet.root");
    corrector_->addDataFile("/nfs/dust/cms/user/pooja/scratch/plot-macro/RecoilCorrector_v7/recoilfits/recoilfit_datamm53X_2012_njet.root");
  }
  else if(sampleName_.find("Higgs")!=std::string::npos){
    corrector_ = new RecoilCorrector("/nfs/dust/cms/user/pooja/scratch/plot-macro/RecoilCorrector_v7/recoilfits/recoilfit_higgs53X_20pv_njet.root");
    corrector_->addMCFile("/nfs/dust/cms/user/pooja/scratch/plot-macro/RecoilCorrector_v7/recoilfits/recoilfit_zmm53X_2012_njet.root");
    corrector_->addDataFile("/nfs/dust/cms/user/pooja/scratch/plot-macro/RecoilCorrector_v7/recoilfits/recoilfit_datamm53X_2012_njet.root");
  }
  else
    corrector_ = 0;
}


TLorentzVector MEtRecoilCorrection::RecoilCorrectedMET(TLorentzVector pfMet_, LorentzVector lepLeg1_, LorentzVector lepLeg2_, LorentzVector genBoson_, string sampleName_, int nJets_)
{
  double newPfMetPt_ = pfMet_.Pt(); 
  double newPfMetPhi_ = pfMet_.Phi();
  
  double u1 = 0; double u2 = 0;
  LorentzVector diLepton_(0, 0, 0, 0);
  
  // case -I   
  if( (sampleName_.find("Zjets")!=string::npos) || (sampleName_.find("Higgs")!=string::npos) ) 
    diLepton_ = lepLeg1_ + lepLeg2_;
  
  
  // case -II   
  if(sampleName_.find("WJets")!=string::npos)
    diLepton_ = lepLeg1_;
  
  
  // TypeI Correction
  TLorentzVector corMET_;
  corrector_->CorrectType1(newPfMetPt_,newPfMetPhi_, genBoson_.pt() ,genBoson_.phi(), diLepton_.Pt(), diLepton_.Phi(), u1, u2 , 0 , 0, TMath::Min(nJets_,2) );
  corMET_.SetPtEtaPhiM(newPfMetPt_, 0, newPfMetPhi_, 0);
 
  return corMET_;
}

#endif // _HTAUTAU

