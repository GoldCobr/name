double fitfunc(double* x, double* par){
    //return par[0]*exp(-x[0]) + par[1] * exp(-0.5*pow((x[0] - par[3]), 2)/(pow(par[2],2)));
    return par[0]*exp(-x[0]) + par[1]*TMath::Landau(x[0], par[2], par[3], false);
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

void task_12(){    
    TH1D* hist1 = new TH1D("h1", "data from data3; units;enteries", 100, 0, 10);
    std::string hist1_data = "data3.txt";
    
    read_fill(hist1_data, hist1);

    TF1* fitting = new TF1("fitting", "fitfunc", 0, 10, 4);

    fitting->FixParameter(0, 21);    
    fitting->FixParameter(1, 93);
    fitting->FixParameter(2, 4.8);
    fitting->FixParameter(3, 0.87);

    auto c = new TCanvas("c", "c");

    hist1->Draw();
    hist1->Fit(fitting, "R");
    c->Draw();

}