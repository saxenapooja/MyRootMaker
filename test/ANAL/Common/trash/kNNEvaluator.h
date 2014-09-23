
#ifndef _WHANALYSIS_MVAEVALUATOR_H_
#define _WHANALYSIS_MVAEVALUATOR_H_

#include <TMVA/Reader.h>
#include <TF1.h>

class SinglekNNEvaluator
{
public:
  SinglekNNEvaluator(const char* weightsFile);

  double mvaValue(double lepPt, double jetPt, unsigned int nJets20) const;

protected:
  TMVA::Reader* reader; 

  mutable float leptonPt_;
  mutable float jetPt_;
  mutable float nJets20_;
};

class Basic_kNNEvaluator
{
public:
  Basic_kNNEvaluator() {}
  virtual ~Basic_kNNEvaluator() {}

  virtual double getLeadingMuon(double lepPt, double jetPt, unsigned int nJets20) const = 0;
  virtual double getSubleadingMuon(double lepPt, double jetPt, unsigned int nJets20) const = 0;
};

class kNNEvaluator: public Basic_kNNEvaluator
{
public:
  kNNEvaluator(const char* data_era, const char* version, const char* iso, const char* nkNN, const char* sample, const char* region);

  virtual double getLeadingMuon(double lepPt, double jetPt, unsigned int nJets20) const;
  virtual double getSubleadingMuon(double lepPt, double jetPt, unsigned int nJets20) const;
protected:
  SinglekNNEvaluator leading_muon;
  SinglekNNEvaluator subleading_muon;
};

#if 0
class FRkNNEvaluator: public Basic_kNNEvaluator
{
public:
  FRkNNEvaluator(const char* data_era, const char* version, const char* iso, const char* variable, bool wz_zz);

  virtual double getLeadingMuon(double lepPt, double jetPt, unsigned int nJets20) const;
  virtual double getSubleadingMuon(double lepPt, double jetPt, unsigned int nJets20) const;
protected:
  TFile* file;
  TF1* leading_muon;
  TF1* subleading_muon;

  enum Variable { LEP_PT, JET_PT, JET_PT_000_020, JET_PT_020_040, JET_PT_040_INF };
  Variable variable_;
};

class FRZkNNEvaluator: public Basic_kNNEvaluator
{
public:
  FRZkNNEvaluator(const char* data_era, const char* version, const char* iso, bool wz_zz);

  virtual double getQCDLeadingMuon(double lepPt, double jetPt, unsigned int nJets20) const;
  virtual double getQCDSubleadingMuon(double lepPt, double jetPt, unsigned int nJets20) const;
  virtual double getWJetsLeadingMuon(double lepPt, double jetPt, unsigned int nJets20) const;
  virtual double getWJetsSubleadingMuon(double lepPt, double jetPt, unsigned int nJets20) const;
protected:
  TFile* file;
  TF1* zjets_leading_muon;
  TF1* zjets_subleading_muon;
};

class FR2DkNNEvaluator: public Basic_kNNEvaluator
{
public:
  FR2DkNNEvaluator(const char* data_era, const char* version, const char* iso, bool wz_zz);

  virtual double getQCDLeadingMuon(double lepPt, double jetPt, unsigned int nJets20) const;
  virtual double getQCDSubleadingMuon(double lepPt, double jetPt, unsigned int nJets20) const;
  virtual double getWJetsLeadingMuon(double lepPt, double jetPt, unsigned int nJets20) const;
  virtual double getWJetsSubleadingMuon(double lepPt, double jetPt, unsigned int nJets20) const;
protected:
  FRkNNEvaluator fr000_020;
  FRkNNEvaluator fr020_040;
  FRkNNEvaluator fr040_inf;
};

class MaurokNNEvaluator: public Basic_kNNEvaluator
{
public:
  MaurokNNEvaluator(const char* data_era);

  virtual double getQCDLeadingMuon(double lepPt, double jetPt, unsigned int nJets20) const;
  virtual double getQCDSubleadingMuon(double lepPt, double jetPt, unsigned int nJets20) const;
  virtual double getWJetsLeadingMuon(double lepPt, double jetPt, unsigned int nJets20) const;
  virtual double getWJetsSubleadingMuon(double lepPt, double jetPt, unsigned int nJets20) const;
protected:

  TMVA::Reader* qcd_leading_muon; 
  TMVA::Reader* qcd_subleading_muon; 
  TMVA::Reader* wjets_leading_muon; 
  TMVA::Reader* wjets_subleading_muon; 

  mutable float leptonPt_;
  mutable float jetPt_;
  mutable float nJets20_;
};
#endif

#endif // _WHANALYSIS_MVAEVALUATOR_H_
