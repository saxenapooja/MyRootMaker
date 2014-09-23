#include "Common/Cutflow.h"
#include "Common/CutFlowTree.h"
#include <Math/VectorUtil.h>

static double calcMt(double px1, double py1, double px2, double py2)
{
        double pt1 = TMath::Sqrt(px1*px1 + py1*py1);
        double pt2 = TMath::Sqrt(px2*px2 + py2*py2);

        double p1Dotp2 = px1*px2 + py1*py2;
        double cosAlpha = p1Dotp2/(pt1*pt2);

        return TMath::Sqrt(2*pt1*pt2*(1 - cosAlpha));
}

CutFlowTree::Variables::Variables(const Muon& muon, const TauCand& OStau, const TauCand& SStau, const TLorentzVector& met, unsigned int njets30, int regions, bool case2):
  regions(regions), case2(case2)
{
  init(muon, OStau, SStau, met, njets30);
}

CutFlowTree::Variables::Variables(const Electron& electron, const TauCand& OStau, const TauCand& SStau, const TLorentzVector& met, unsigned int njets30, int regions, bool case2):
  regions(regions), case2(case2)
{
  init(electron, OStau, SStau, met, njets30);
}

void CutFlowTree::Variables::init(const TLorentzVector& lepton, const TauCand& OStau, const TauCand& SStau, const TLorentzVector& met, unsigned int njets30)
{
  mvaVars.leadPt = TMath::Max(OStau.Pt(), SStau.Pt());
  mvaVars.subPt = TMath::Min(OStau.Pt(), SStau.Pt());
  mvaVars.OSPt = OStau.Pt();
  mvaVars.OSEta = OStau.Eta();
  mvaVars.SSPt = SStau.Pt();
  mvaVars.SSEta = SStau.Eta();
  //mvaVars.Lt = lepton.Pt() + OStau.Pt() + SStau.Pt();
  mvaVars.Met = met.Et();
  mvaVars.deltaRDiTau = ROOT::Math::VectorUtil::DeltaR(OStau.p4, SStau.p4);
  mvaVars.deltaPhiDiTau = ROOT::Math::VectorUtil::DeltaPhi(OStau.p4, SStau.p4);
  mvaVars.deltaEtaDiTau = OStau.Eta() - SStau.Eta();
  //mvaVars.etaDiTau = (OStau.p4 + SStau.p4).Eta();
  //mvaVars.etaLepTauTau = (lepton + OStau.p4 + SStau.p4).Eta();
  //mvaVars.deltaRLepOSTau = ROOT::Math::VectorUtil::DeltaR(lepton, OStau.p4);
  //mvaVars.deltaRLepSSTau = ROOT::Math::VectorUtil::DeltaR(lepton, SStau.p4);
  mvaVars.ptRatio = (OStau.p4 + SStau.p4).Pt() / (OStau.Pt() + SStau.Pt());

  const TVector3 zDir(0,0,1);
  const TVector3 lepDir( lepton.Px(), lepton.Py(), lepton.Pz());
  const TVector3 tau1Dir( OStau.p4.Px(), OStau.p4.Py(), OStau.p4.Pz());
  const TVector3 tau2Dir( SStau.p4.Px(), SStau.p4.Py(), SStau.p4.Pz());
  const TVector3 ditauDir= tau1Dir.Cross(tau2Dir);
  const TVector3 zlep = zDir.Cross(lepDir);
  mvaVars.angle = zlep.Angle(ditauDir);

  //mvaVars.nJets30 = njets30;
  //mvaVars.etaLeadTau = OStau.Pt() > SStau.Pt() ? OStau.Eta() : SStau.Eta();
  //mvaVars.etaSubTau = OStau.Pt() < SStau.Pt() ? OStau.Eta() : SStau.Eta();
  mvaVars.transMassLepMET = calcMt(lepton.Px(), lepton.Py(), met.Px(), met.Py());
  //mvaVars.transMassOSTauMET = calcMt(OStau.p4.Px(), OStau.p4.Py(), met.Px(), met.Py());
  //mvaVars.transMassSSTauMET = calcMt(SStau.p4.Px(), SStau.p4.Py(), met.Px(), met.Py());
  //mvaVars.deltaPhiLepMET = ROOT::Math::VectorUtil::DeltaPhi(lepton, met);

  //mvaVars.dxyOSTau = OStau.Dxy;
  //mvaVars.dxySSTau = SStau.Dxy;
  //mvaVars.dzOSTau = OStau.Dz;
  //mvaVars.dzSSTau = SStau.Dz;
}

CutFlowTree::CutFlowTree(const char* name, const char* desc)
{
  TDirectory* prevDirectory = gDirectory;
  
  gDirectory->mkdir(name);
  gDirectory->cd(name);

  mvaTree = new TTree("signal", desc);
  mvaTreeTau12Fake = new TTree("tau12fake", (std::string(desc) + " of both #tau fakes").c_str());
  mvaTreeTauSSFake = new TTree("tauSSfake", (std::string(desc) + " of same sign #tau fakes").c_str());

  mvaTree->Branch("LeadPt", &mvaVars.mvaVars.leadPt, "LeadPt/F");
  mvaTree->Branch("SubPt", &mvaVars.mvaVars.subPt, "SubPt/F");
  mvaTree->Branch("OSPt", &mvaVars.mvaVars.OSPt, "OSPt/F");
  mvaTree->Branch("SSPt", &mvaVars.mvaVars.SSPt, "SSPt/F");
  //mvaTree->Branch("Lt", &mvaVars.mvaVars.Lt, "Lt/F");
  mvaTree->Branch("Met", &mvaVars.mvaVars.Met, "MEt/F");
  mvaTree->Branch("DeltaRDiTau", &mvaVars.mvaVars.deltaRDiTau, "DeltaRDiTau/F");
  mvaTree->Branch("DeltaPhiDiTau", &mvaVars.mvaVars.deltaPhiDiTau, "DeltaPhiDiTau/F");
  mvaTree->Branch("DeltaEtaDiTau", &mvaVars.mvaVars.deltaEtaDiTau, "DeltaEtaDiTau/F");
  //mvaTree->Branch("EtaDiTau", &mvaVars.mvaVars.etaDiTau, "EtaDiTau/F");
  //mvaTree->Branch("EtaLepTauTau", &mvaVars.mvaVars.etaLepTauTau, "EtaLepTauTau/F");
  //mvaTree->Branch("DeltaRLepOSTau", &mvaVars.mvaVars.deltaRLepOSTau, "DeltaRLepOSTau/F");
  //mvaTree->Branch("DeltaRLepSSTau", &mvaVars.mvaVars.deltaRLepSSTau, "DeltaRLepSSTau/F");
  mvaTree->Branch("PtRatio", &mvaVars.mvaVars.ptRatio, "PtRatio/F");
  mvaTree->Branch("Angle", &mvaVars.mvaVars.angle, "Angle/F");
  //mvaTree->Branch("NJets30", &mvaVars.mvaVars.nJets30, "NJets30/i");
  //mvaTree->Branch("EtaLeadTau", &mvaVars.mvaVars.etaLeadTau, "EtaLeadTau/F");
  //mvaTree->Branch("EtaSubTau", &mvaVars.mvaVars.etaSubTau, "EtaSubTau/F");
  mvaTree->Branch("TransMassLepMET", &mvaVars.mvaVars.transMassLepMET, "TransMassLepMET/F");
  //mvaTree->Branch("TransMassOSTauMET", &mvaVars.mvaVars.transMassOSTauMET, "TransMassOSTauMET/F");
  //mvaTree->Branch("TransMassSSTauMET", &mvaVars.mvaVars.transMassSSTauMET, "TransMassSSTauMET/F");
  //mvaTree->Branch("DeltaPhiLepMET", &mvaVars.mvaVars.deltaPhiLepMET, "DeltaPhiLepMET/F");
  //mvaTree->Branch("DxyOSTau", &mvaVars.mvaVars.dxyOSTau, "DxyOSTau/F");
  //mvaTree->Branch("DxySSTau", &mvaVars.mvaVars.dxySSTau, "DxySSTau/F");
  //mvaTree->Branch("DzOSTau", &mvaVars.mvaVars.dzOSTau, "DzOSTau/F");
  //mvaTree->Branch("DzSSTau", &mvaVars.mvaVars.dzSSTau, "DzSSTau/F");
  mvaTree->Branch("Weight", &mvaVars.weight, "Weight/F");
  mvaTree->Branch("WWeight", &mvaVars.Wweight, "WWeight/F");
  mvaTree->Branch("ZWeight", &mvaVars.Zweight, "ZWeight/F");

  mvaTreeTau12Fake->Branch("LeadPt", &mvaVarsTau12Fake.mvaVars.leadPt, "LeadPt/F");
  mvaTreeTau12Fake->Branch("SubPt", &mvaVarsTau12Fake.mvaVars.subPt, "SubPt/F");
  mvaTreeTau12Fake->Branch("OSPt", &mvaVarsTau12Fake.mvaVars.OSPt, "OSPt/F");
  mvaTreeTau12Fake->Branch("SSPt", &mvaVarsTau12Fake.mvaVars.SSPt, "SSPt/F");
  //mvaTreeTau12Fake->Branch("Lt", &mvaVarsTau12Fake.mvaVars.Lt, "Lt/F");
  mvaTreeTau12Fake->Branch("Met", &mvaVarsTau12Fake.mvaVars.Met, "MEt/F");
  mvaTreeTau12Fake->Branch("DeltaRDiTau", &mvaVarsTau12Fake.mvaVars.deltaRDiTau, "DeltaRDiTau/F");
  mvaTreeTau12Fake->Branch("DeltaPhiDiTau", &mvaVarsTau12Fake.mvaVars.deltaPhiDiTau, "DeltaPhiDiTau/F");
  mvaTreeTau12Fake->Branch("DeltaEtaDiTau", &mvaVarsTau12Fake.mvaVars.deltaEtaDiTau, "DeltaEtaDiTau/F");
  //mvaTreeTau12Fake->Branch("EtaDiTau", &mvaVarsTau12Fake.mvaVars.etaDiTau, "EtaDiTau/F");
  //mvaTreeTau12Fake->Branch("EtaLepTauTau", &mvaVarsTau12Fake.mvaVars.etaLepTauTau, "EtaLepTauTau/F");
  //mvaTreeTau12Fake->Branch("DeltaRLepOSTau", &mvaVarsTau12Fake.mvaVars.deltaRLepOSTau, "DeltaRLepOSTau/F");
  //mvaTreeTau12Fake->Branch("DeltaRLepSSTau", &mvaVarsTau12Fake.mvaVars.deltaRLepSSTau, "DeltaRLepSSTau/F");
  mvaTreeTau12Fake->Branch("PtRatio", &mvaVarsTau12Fake.mvaVars.ptRatio, "PtRatio/F");
  mvaTreeTau12Fake->Branch("Angle", &mvaVarsTau12Fake.mvaVars.angle, "Angle/F");
  //mvaTreeTau12Fake->Branch("NJets30", &mvaVarsTau12Fake.mvaVars.nJets30, "NJets30/i");
  //mvaTreeTau12Fake->Branch("EtaLeadTau", &mvaVarsTau12Fake.mvaVars.etaLeadTau, "EtaLeadTau/F");
  //mvaTreeTau12Fake->Branch("EtaSubTau", &mvaVarsTau12Fake.mvaVars.etaSubTau, "EtaSubTau/F");
  mvaTreeTau12Fake->Branch("TransMassLepMET", &mvaVarsTau12Fake.mvaVars.transMassLepMET, "TransMassLepMET/F");
  //mvaTreeTau12Fake->Branch("TransMassOSTauMET", &mvaVarsTau12Fake.mvaVars.transMassOSTauMET, "TransMassOSTauMET/F");
  //mvaTreeTau12Fake->Branch("TransMassSSTauMET", &mvaVarsTau12Fake.mvaVars.transMassSSTauMET, "TransMassSSTauMET/F");
  //mvaTreeTau12Fake->Branch("DeltaPhiLepMET", &mvaVarsTau12Fake.mvaVars.deltaPhiLepMET, "DeltaPhiLepMET/F");
  //mvaTreeTau12Fake->Branch("DxyOSTau", &mvaVarsTau12Fake.mvaVars.dxyOSTau, "DxyOSTau/F");
  //mvaTreeTau12Fake->Branch("DxySSTau", &mvaVarsTau12Fake.mvaVars.dxySSTau, "DxySSTau/F");
  //mvaTreeTau12Fake->Branch("DzOSTau", &mvaVarsTau12Fake.mvaVars.dzOSTau, "DzOSTau/F");
  //mvaTreeTau12Fake->Branch("DzSSTau", &mvaVarsTau12Fake.mvaVars.dzSSTau, "DzSSTau/F");
  mvaTreeTau12Fake->Branch("Weight", &mvaVarsTau12Fake.weight, "Weight/F");
  mvaTreeTau12Fake->Branch("WWeight", &mvaVarsTau12Fake.Wweight, "WWeight/F");
  mvaTreeTau12Fake->Branch("ZWeight", &mvaVarsTau12Fake.Zweight, "ZWeight/F");

  mvaTreeTauSSFake->Branch("LeadPt", &mvaVarsTauSSFake.mvaVars.leadPt, "LeadPt/F");
  mvaTreeTauSSFake->Branch("SubPt", &mvaVarsTauSSFake.mvaVars.subPt, "SubPt/F");
  mvaTreeTauSSFake->Branch("OSPt", &mvaVarsTauSSFake.mvaVars.OSPt, "OSPt/F");
  mvaTreeTauSSFake->Branch("SSPt", &mvaVarsTauSSFake.mvaVars.SSPt, "SSPt/F");
  //mvaTreeTauSSFake->Branch("Lt", &mvaVarsTauSSFake.mvaVars.Lt, "Lt/F");
  mvaTreeTauSSFake->Branch("Met", &mvaVarsTauSSFake.mvaVars.Met, "MEt/F");
  mvaTreeTauSSFake->Branch("DeltaRDiTau", &mvaVarsTauSSFake.mvaVars.deltaRDiTau, "DeltaRDiTau/F");
  mvaTreeTauSSFake->Branch("DeltaPhiDiTau", &mvaVarsTauSSFake.mvaVars.deltaPhiDiTau, "DeltaPhiDiTau/F");
  mvaTreeTauSSFake->Branch("DeltaEtaDiTau", &mvaVarsTauSSFake.mvaVars.deltaEtaDiTau, "DeltaEtaDiTau/F");
  //mvaTreeTauSSFake->Branch("EtaDiTau", &mvaVarsTauSSFake.mvaVars.etaDiTau, "EtaDiTau/F");
  //mvaTreeTauSSFake->Branch("EtaLepTauTau", &mvaVarsTauSSFake.mvaVars.etaLepTauTau, "EtaLepTauTau/F");
  //mvaTreeTauSSFake->Branch("DeltaRLepOSTau", &mvaVarsTauSSFake.mvaVars.deltaRLepOSTau, "DeltaRLepOSTau/F");
  //mvaTreeTauSSFake->Branch("DeltaRLepSSTau", &mvaVarsTauSSFake.mvaVars.deltaRLepSSTau, "DeltaRLepSSTau/F");
  mvaTreeTauSSFake->Branch("PtRatio", &mvaVarsTauSSFake.mvaVars.ptRatio, "PtRatio/F");
  mvaTreeTauSSFake->Branch("Angle", &mvaVarsTauSSFake.mvaVars.angle, "Angle/F");
  //mvaTreeTauSSFake->Branch("NJets30", &mvaVarsTauSSFake.mvaVars.nJets30, "NJets30/i");
  //mvaTreeTauSSFake->Branch("EtaLeadTau", &mvaVarsTauSSFake.mvaVars.etaLeadTau, "EtaLeadTau/F");
  //mvaTreeTauSSFake->Branch("EtaSubTau", &mvaVarsTauSSFake.mvaVars.etaSubTau, "EtaSubTau/F");
  mvaTreeTauSSFake->Branch("TransMassLepMET", &mvaVarsTauSSFake.mvaVars.transMassLepMET, "TransMassLepMET/F");
  //mvaTreeTauSSFake->Branch("TransMassOSTauMET", &mvaVarsTauSSFake.mvaVars.transMassOSTauMET, "TransMassOSTauMET/F");
  //mvaTreeTauSSFake->Branch("TransMassSSTauMET", &mvaVarsTauSSFake.mvaVars.transMassSSTauMET, "TransMassSSTauMET/F");
  //mvaTreeTauSSFake->Branch("DeltaPhiLepMET", &mvaVarsTauSSFake.mvaVars.deltaPhiLepMET, "DeltaPhiLepMET/F");
  //mvaTreeTauSSFake->Branch("DxyOSTau", &mvaVarsTauSSFake.mvaVars.dxyOSTau, "DxyOSTau/F");
  //mvaTreeTauSSFake->Branch("DxySSTau", &mvaVarsTauSSFake.mvaVars.dxySSTau, "DxySSTau/F");
  //mvaTreeTauSSFake->Branch("DzOSTau", &mvaVarsTauSSFake.mvaVars.dzOSTau, "DzOSTau/F");
  //mvaTreeTauSSFake->Branch("DzSSTau", &mvaVarsTauSSFake.mvaVars.dzSSTau, "DzSSTau/F");
  mvaTreeTauSSFake->Branch("Weight", &mvaVarsTauSSFake.weight, "Weight/F");
  mvaTreeTauSSFake->Branch("WWeight", &mvaVarsTauSSFake.Wweight, "WWeight/F");
  mvaTreeTauSSFake->Branch("ZWeight", &mvaVarsTauSSFake.Zweight, "ZWeight/F");

  prevDirectory->cd();
}

void CutFlowTree::fill(const Variables& vars, double weight)
{
  const int regions = vars.regions;
  const bool case2 = vars.case2;

//   const double wosw = fakeRateFuncs.TauOS->getPtEta(FakeRateFunctions::WJETS_FUNC, vars.mvaVars.OSPt, vars.mvaVars.OSEta, case2);
//   const double wssw = fakeRateFuncs.TauSS->getPtEta(FakeRateFunctions::WJETS_FUNC, vars.mvaVars.SSPt, vars.mvaVars.SSEta, case2);
//   const double wosz = fakeRateFuncs.TauOS->getPtEta(FakeRateFunctions::ZJETS_FUNC, vars.mvaVars.OSPt, vars.mvaVars.OSEta, case2);
//   const double wssz = fakeRateFuncs.TauSS->getPtEta(FakeRateFunctions::ZJETS_FUNC, vars.mvaVars.SSPt, vars.mvaVars.SSEta, case2);

  if(regions & 1) //CutflowRegions::SIGNAL)
  {
    mvaVars.mvaVars = vars.mvaVars;
    mvaVars.weight = weight;
    mvaVars.Wweight = 1.;
    mvaVars.Zweight = 1.;
    mvaTree->Fill();
  }

  if(regions & 2) //CutflowRegions::BACKGROUND)
  {
    mvaVarsTauSSFake.mvaVars = vars.mvaVars;
    mvaVarsTauSSFake.weight = weight;
//     mvaVarsTauSSFake.Wweight = wssw/(1.0 - wssw);
//     mvaVarsTauSSFake.Zweight = wssz/(1.0 - wssz);
    mvaTreeTauSSFake->Fill();
  }

  if(regions & 8) //CutflowRegions::CONTROL_BACKGROUND)
  {
    mvaVarsTau12Fake.mvaVars = vars.mvaVars;
    mvaVarsTau12Fake.weight = weight;
//     mvaVarsTau12Fake.Wweight = wosw/(1.0 - wosw) * wssw/(1.0 - wssw);
//     mvaVarsTau12Fake.Zweight = wosz/(1.0 - wosz) * wssz/(1.0 - wssz);
    mvaTreeTau12Fake->Fill();
  }
}
