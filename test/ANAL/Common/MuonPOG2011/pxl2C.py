#!/usr/bin/env python
import pickle
f = open('MuonEfficiencies2011_42X.pkl', 'r')
dict = pickle.load(f)

def ptOutput(ptOut):
	return """
if(pt >= 10. && pt < 20.) return %g;
else if(pt >= 20. && pt < 30.) return %g;
else if(pt >= 30. && pt < 40.) return %g;
else if(pt >= 40. && pt < 50.) return %g;
else if(pt >= 50. && pt < 60.) return %g;
else if(pt >= 60. && pt < 80.) return %g;
else if(pt >= 80.) return %g;
else assert(false);""" % (
		ptOut['10_20']['data/mc']['efficiency_ratio'],
		ptOut['20_30']['data/mc']['efficiency_ratio'],
		ptOut['30_40']['data/mc']['efficiency_ratio'],
		ptOut['40_50']['data/mc']['efficiency_ratio'],
		ptOut['50_60']['data/mc']['efficiency_ratio'],
		ptOut['60_80']['data/mc']['efficiency_ratio'],
		ptOut['80_250']['data/mc']['efficiency_ratio']
	)

def etaOutput(etaOut):
	return """
if(fabs(eta) < 1.2)
{
	%s
}
else
{
	%s
}""" % (
	ptOutput(etaOut['pt_abseta<1.2']),
	ptOutput(etaOut['pt_abseta>1.2']),
	)

def funcOutput(dict, funcName):
	etaOut = dict[funcName]

	return """
static double %s(double pt, double eta)
{
	%s
}""" % (funcName, etaOutput(etaOut))

efficiencies = ['PFTIGHT_nL8_2011A', 'PFTIGHT_nL8_2011B', 'combRelPFISO12_2011A', 'combRelPFISO12_2011B']
out_text = "\n\n".join([funcOutput(dict, x) for x in efficiencies])
open('efficiencyFunctionsIDIso.C', 'w').write(out_text)
