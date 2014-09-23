
static double PFTIGHT_nL8_2011A(double pt, double eta)
{
	
if(fabs(eta) < 1.2)
{
	
if(pt >= 10. && pt < 20.) return 0.98822;
else if(pt >= 20. && pt < 30.) return 0.993888;
else if(pt >= 30. && pt < 40.) return 0.993331;
else if(pt >= 40. && pt < 50.) return 0.993305;
else if(pt >= 50. && pt < 60.) return 0.994317;
else if(pt >= 60. && pt < 80.) return 0.986547;
else if(pt >= 80.) return 0.987065;
else assert(false);
}
else
{
	
if(pt >= 10. && pt < 20.) return 0.995575;
else if(pt >= 20. && pt < 30.) return 0.994378;
else if(pt >= 30. && pt < 40.) return 0.991713;
else if(pt >= 40. && pt < 50.) return 0.991514;
else if(pt >= 50. && pt < 60.) return 0.990772;
else if(pt >= 60. && pt < 80.) return 0.989951;
else if(pt >= 80.) return 0.993597;
else assert(false);
}
}


static double PFTIGHT_nL8_2011B(double pt, double eta)
{
	
if(fabs(eta) < 1.2)
{
	
if(pt >= 10. && pt < 20.) return 0.985259;
else if(pt >= 20. && pt < 30.) return 0.99125;
else if(pt >= 30. && pt < 40.) return 0.987533;
else if(pt >= 40. && pt < 50.) return 0.988953;
else if(pt >= 50. && pt < 60.) return 0.992382;
else if(pt >= 60. && pt < 80.) return 0.983605;
else if(pt >= 80.) return 0.979509;
else assert(false);
}
else
{
	
if(pt >= 10. && pt < 20.) return 0.968447;
else if(pt >= 20. && pt < 30.) return 0.96397;
else if(pt >= 30. && pt < 40.) return 0.962589;
else if(pt >= 40. && pt < 50.) return 0.966384;
else if(pt >= 50. && pt < 60.) return 0.965024;
else if(pt >= 60. && pt < 80.) return 0.966253;
else if(pt >= 80.) return 0.975424;
else assert(false);
}
}


static double combRelPFISO12_2011A(double pt, double eta)
{
	
if(fabs(eta) < 1.2)
{
	
if(pt >= 10. && pt < 20.) return 0.955401;
else if(pt >= 20. && pt < 30.) return 0.981241;
else if(pt >= 30. && pt < 40.) return 0.988022;
else if(pt >= 40. && pt < 50.) return 0.992949;
else if(pt >= 50. && pt < 60.) return 0.993477;
else if(pt >= 60. && pt < 80.) return 0.994019;
else if(pt >= 80.) return 0.993763;
else assert(false);
}
else
{
	
if(pt >= 10. && pt < 20.) return 0.985755;
else if(pt >= 20. && pt < 30.) return 0.990614;
else if(pt >= 30. && pt < 40.) return 0.999142;
else if(pt >= 40. && pt < 50.) return 0.999988;
else if(pt >= 50. && pt < 60.) return 0.999754;
else if(pt >= 60. && pt < 80.) return 0.99811;
else if(pt >= 80.) return 0.99369;
else assert(false);
}
}


static double combRelPFISO12_2011B(double pt, double eta)
{
	
if(fabs(eta) < 1.2)
{
	
if(pt >= 10. && pt < 20.) return 0.905921;
else if(pt >= 20. && pt < 30.) return 0.940926;
else if(pt >= 30. && pt < 40.) return 0.968547;
else if(pt >= 40. && pt < 50.) return 0.980621;
else if(pt >= 50. && pt < 60.) return 0.986195;
else if(pt >= 60. && pt < 80.) return 0.988684;
else if(pt >= 80.) return 0.994216;
else assert(false);
}
else
{
	
if(pt >= 10. && pt < 20.) return 0.942515;
else if(pt >= 20. && pt < 30.) return 0.947486;
else if(pt >= 30. && pt < 40.) return 0.964939;
else if(pt >= 40. && pt < 50.) return 0.98663;
else if(pt >= 50. && pt < 60.) return 0.991718;
else if(pt >= 60. && pt < 80.) return 0.993393;
else if(pt >= 80.) return 1.00566;
else assert(false);
}
}
