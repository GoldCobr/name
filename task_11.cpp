double hist_1[100]; double hist_2[100];double error_hist_1[100]; double error_hist_2[100]; 
double x[100];
double chi2;

double func(double x, double *par){
    return par[0]/(sqrt(2*M_PI)*par[2]) * TMath::Gaus(x, par[1], par[2]);
}
double func_1(double* x, double *par){
    return par[0]/(sqrt(2*M_PI)*par[2]) * TMath::Gaus(x[0], par[1], par[2]) + par[3];
}

void fcn(int &npar, double *gin, double &f, double *par, int iflag)
{
   const int nbins = 100;
   int i;
 
//calculate chisquare
    double chisq = 0;
    chi2 = 0;
    for (int i = 0; i < nbins; i++) {
        double delta1 = TMath::Poisson(hist_1[i], func(x[i], par) + par[3]);
        double delta2 = TMath::Poisson(hist_2[i], par[3]);
        chisq += -2. * log(delta1) - 2. * log(delta2);
    }
    chi2 = chisq;
    f = chisq;
}

void read_fill(std::string name_file, TH1D* hist){
    std::ifstream new_file(name_file);
    double_t x;
    if(new_file.is_open()){
        while(new_file >> x){
            hist->Fill(x);
        }
    }
    else
        std::cout << "Didn't read lol" << std::endl;
    new_file.close();
}

void task_11() {
    TH1D* hist1 = new TH1D("h1", "data from data_1; units;enteries", 100, 500, 600);
    TH1D* hist2 = new TH1D("h2", "data from data_2; units;enteries", 100, 500, 600);

    std::string hist1_data = "data_1.txt";
    std::string hist2_data = "data_2.txt";
    
    read_fill(hist1_data, hist1);
    read_fill(hist2_data, hist2);
    auto *canv = new TCanvas("c", "", 1200, 900);
    canv->Divide(1,2);
    canv->cd(1);
    hist1->Draw("e");

    for (int i = 0; i < 100; i++) {
        hist_1[i] = hist1->GetBinContent(i);
        hist_2[i] = hist2->GetBinContent(i);
    }
    
    for (int i = 0; i < 100; i++) x[i] = i + 500;
    
    //The errors values
    // for (int i = 0; i < 100; i++) {
    //     if (hist_1[i] == 0) error_hist_1[i] = sqrt(3.09);
    //     else {
    //         error_hist_1[i] = sqrt(hist_1[i]);
    //     }
    // }

    // for (int i = 0; i < 100; i++) {
    //     if (hist_2[i] == 0){
    //         error_hist_2[i] = sqrt(3.09);
    //     }
    //     else{
    //         error_hist_2[i] = sqrt(hist_2[i]);
    //     }
    // }
    
    TMinuit *gMinuit = new TMinuit(4);
    gMinuit->SetFCN(fcn);

    double arglist[10];
    int ierflg = 0;

    arglist[0] = 1;
    gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);

    // Set starting values and step sizes for parameters
    static double vstart[4] = {0.1, 550., 10., 10.};
    static double step[4] = {0.1, 0.1 , 0.1 , 0.1};
    gMinuit->mnparm(0, "ampl", vstart[0], step[0], 0, 0, ierflg);
    gMinuit->mnparm(1, "mean", vstart[1], step[1], 0, 0, ierflg);
    gMinuit->mnparm(2, "sigma", vstart[2], step[2], 0, 0, ierflg);
    gMinuit->mnparm(3, "const", vstart[3], step[3], 0, 0, ierflg);

    arglist[0] = 500;
    arglist[1] = 0.1;
    gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);

    double amin, edm, errdef;
    int nvpar, nparx, icstat;
    gMinuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);
    gMinuit->mnprin(3, amin);
    
    double ampl, mean, sigma, Const;
    double ampl_err, mean_err, sigma_err, Const_err;
    gMinuit->GetParameter(0, ampl, ampl_err);
    gMinuit->GetParameter(1, mean, mean_err);
    gMinuit->GetParameter(2, sigma, sigma_err);
    gMinuit->GetParameter(3, Const, Const_err);

    auto fitfunc1 = new TF1("f1", "func_1", 500, 600,4);
    fitfunc1->SetParameter(0, ampl);
    fitfunc1->SetParameter(1, mean);
    fitfunc1->SetParameter(2, sigma);
    fitfunc1->SetParameter(3, Const);
    fitfunc1->Draw("Same");

    auto fitfunc2 = new TF1("f2", "[3]", 500, 600);
    fitfunc2->SetParameter(3, Const);
    canv->cd(2);
    hist2->Draw("e");
    fitfunc2->Draw("Same");

    std::cout << "chi2: " << chi2 <<std::endl;
    std::cout << "number of events = " << ampl << " ampl_err: "<< ampl_err<<std::endl;

    TFile* ofile = new TFile("ofile_task_11.root", "RECREATE");
    hist1->Write();
    hist2->Write();
    fitfunc1->Write();
    fitfunc2->Write();
    ofile->Close();
}