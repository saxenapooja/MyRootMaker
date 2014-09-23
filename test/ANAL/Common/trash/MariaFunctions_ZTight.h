//////////////FR functions(pt) 1D
static Double_t fake_rate_tau_pt_20(Double_t pt) {
    if (pt<0) return 1.0;
    if (pt<=100){
        //1D
        Double_t mpv=-16.7;
        Double_t  sigma=7.405;
        Double_t costante=0.001794;
        
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
            Double_t mpv=-5.336;
            Double_t sigma=4.384;
            Double_t costante=0.005082;
        
        Double_t f1=((TMath::Landau(pt,mpv,sigma))+costante);
        
            return f1;}
        //eta2
        if ((TMath::Abs(eta) > 0.8) && (TMath::Abs(eta) < 1.6)){
            Double_t mpv=-13.71;
            Double_t  sigma=7.583;
            Double_t costante=0.0008783;
            
            Double_t f2=((TMath::Landau(pt,mpv,sigma))+costante);
            
            return f2;}
        //eta3
        if ((TMath::Abs(eta) > 1.6) && (TMath::Abs(eta) < 2.3)){
            Double_t mpv=-10.54;
            Double_t  sigma=6.317;
            Double_t costante=0.009449;
            
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
        
        Double_t p_par[14] = {0.411829,-0.330438,3.35952,0.0806,-0.344626,3.1413,0.109505,-0.309343,-0.0359541,0.00541517,-0.000719796,0.00968051,0.000338601,-0.00135203};
        ///parametri!!!!
        Double_t p0 = p_par[0]+p_par[1]*eta+p_par[2]*TMath::Power(eta,2)+p_par[3]*TMath::Power(eta,3)+p_par[4]*TMath::Power(eta,4);
        Double_t p1 = p_par[5]+p_par[6]*eta+p_par[7]*TMath::Power(eta,2)+p_par[8]*TMath::Power(eta,3);
        Double_t p2 = p_par[9]+p_par[10]*eta+p_par[11]*TMath::Power(eta,2)+p_par[12]*TMath::Power(eta,3)+p_par[13]*TMath::Power(eta,4);
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
