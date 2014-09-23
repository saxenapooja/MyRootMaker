//////////////FR functions(pt) 1D
static Double_t fake_rate_tau_pt_20(Double_t pt) {
    if (pt<0) return 1.0;
    if (pt<=100){
        //1D
        Double_t mpv=-4.585;
        Double_t  sigma=3.378;
        Double_t costante=0.002719;
        
        Double_t f=((TMath::Landau(pt,mpv,sigma))+costante);
        
        
        return f;}
    else{
        return fake_rate_tau_pt_20(100);}
}


//////////////FR functions(pt) 3 regions in eta
static Double_t fake_rate_tau_pt_eta(Double_t pt, Double_t eta) {
    if (pt<0) return 1.0;
    if (pt<=100){
        //eta1
        if (TMath::Abs(eta) < 0.8){
            Double_t mpv=2.677;
            Double_t sigma=1.91;
            Double_t costante=0.003176;
        
        Double_t f1=((TMath::Landau(pt,mpv,sigma))+costante);
        
            return f1;}
        //eta2
        if ((TMath::Abs(eta) > 0.8) && (TMath::Abs(eta) < 1.6)){
            Double_t mpv=0.8542;
            Double_t  sigma=2.905;
            Double_t costante=0.004203;
            
            Double_t f2=((TMath::Landau(pt,mpv,sigma))+costante);
            
            return f2;}
        //eta3
        if ((TMath::Abs(eta) > 1.6) && (TMath::Abs(eta) < 2.3)){
            Double_t mpv=6.706;
            Double_t  sigma=2.216;
            Double_t costante=0.007026;
            
            Double_t f3=((TMath::Landau(pt,mpv,sigma))+costante);
            
            return f3;}
        
    }
    else{
        return fake_rate_tau_pt_eta(100,eta);}
}

//////////////FR function(pt,eta) 10
static Double_t fake_rate_tau(Double_t pt, Double_t eta) {
    //if (pt<0) return 1.0;
    if (pt<=100 && (TMath::Abs(eta) < 2.3)){
        
        Double_t p_par[18] = {2.9794,1.15446,0.256635,-0.304231,0.225319,1.52836,-0.243571,0.453055,0.0622961,-0.131668,0.00430408,-0.000362758,0.000812893,0.000604475,0.000987013,0.000170702,-0.000138295,-6.86732e-05};
        ///parametri!!!!
        Double_t p0 = p_par[0]+p_par[1]*eta+p_par[2]*TMath::Power(eta,2)+p_par[3]*TMath::Power(eta,3)+p_par[4]*TMath::Power(eta,4);
        Double_t p1 = p_par[5]+p_par[6]*eta+p_par[7]*TMath::Power(eta,2)+p_par[8]*TMath::Power(eta,3)+p_par[9]*TMath::Power(eta,4);
        Double_t p2 = p_par[10]+p_par[11]*eta+p_par[12]*TMath::Power(eta,2)+p_par[13]*TMath::Power(eta,3)+p_par[14]*TMath::Power(eta,4)+p_par[15]*TMath::Power(eta,5)+p_par[16]*TMath::Power(eta,6)+p_par[17]*TMath::Power(eta,7);
        Double_t result = TMath::Landau(pt,p0,p1)+p2;
        return result;
    }
    else{
        if (pt>100){
            return fake_rate_tau(100,eta);}
        
        if (pt>100){
            return fake_rate_tau(100,eta);}
        
    }
}
