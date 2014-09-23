//////////////FR functions(pt) 1D
static Double_t fake_rate_tau_pt_20(Double_t pt) {
    if (pt<0) return 1.0;
    if (pt<=100){
        //1D
        Double_t mpv=-12.47;
        Double_t  sigma=7.108;
        Double_t costante=0.003016;
        
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
            Double_t mpv=-3.077;
            Double_t sigma=4.431;
            Double_t costante=0.005742;
        
        Double_t f1=((TMath::Landau(pt,mpv,sigma))+costante);
        
            return f1;}
        //eta2
        if ((TMath::Abs(eta) > 0.8) && (TMath::Abs(eta) < 1.6)){
            Double_t mpv=-14.16;
            Double_t  sigma=8.5;
            Double_t costante=-0.0006049;
            
            Double_t f2=((TMath::Landau(pt,mpv,sigma))+costante);
            
            return f2;}
        //eta3
        if ((TMath::Abs(eta) > 1.6) && (TMath::Abs(eta) < 2.3)){
            Double_t mpv=-8.79;
            Double_t  sigma=6.853;
            Double_t costante=0.00795;
            
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
        
        Double_t p_par[15] = {1.85101,-0.346399,2.41707,0.0679917,-0.168245,3.06144,0.102339,0.0513099,-0.0319107,-0.0789404,0.00762862,-0.000676062,0.00774476,0.000376936,-0.000894677};
        ///parametri!!!!
        Double_t p0 = p_par[0]+p_par[1]*eta+p_par[2]*TMath::Power(eta,2)+p_par[3]*TMath::Power(eta,3)+p_par[4]*TMath::Power(eta,4);
        Double_t p1 = p_par[5]+p_par[6]*eta+p_par[7]*TMath::Power(eta,2)+p_par[8]*TMath::Power(eta,3)+p_par[9]*TMath::Power(eta,4);
        Double_t p2 = p_par[10]+p_par[11]*eta+p_par[12]*TMath::Power(eta,2)+p_par[13]*TMath::Power(eta,3)+p_par[14]*TMath::Power(eta,4);
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
