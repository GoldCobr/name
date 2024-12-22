double fitfunc(double* x, double* par){
    return 22*exp(par[0]*x[0]) + par[1] * exp(-0.5*pow((x[0] - par[3]), 2)/(pow(par[2],2))) + 4;
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

    fitting->FixParameter(0, -2);    
    fitting->FixParameter(1, 14);
    fitting->FixParameter(2, 1);
    fitting->FixParameter(3, 5.4);


    auto c = new TCanvas("c", "c");

    hist1->Draw();
    hist1->Fit(fitting, "R");
    c->Draw();

}