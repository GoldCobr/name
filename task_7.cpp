#include <TTree.h>
#include <TFile.h>

double fit_func(double *x, double *par){
    return par[0]*pow(x[0],par[1])*TMath::Exp(par[2]*x[0]);
}

void task_7() {
    //читаем дерево
    TFile ifile("m3pimc.root", "read");
    TTree* h10 = (TTree*)ifile.Get("h10");
    TTree* MyTree = new TTree("MyTree", "tree with data");

    //создаем дерево и файл сохранения
    TFile* ofile = new TFile("newroot.root", "RECREATE");
    TTree* filtered_tree = h10->CopyTree("Isrfilter == 1 && chi2_3p < 30");

    filtered_tree->SetBranchStatus("*", 0);
    for (auto item: {"nph","eph","thetaph","phiph"}) { filtered_tree->SetBranchStatus(item, 1);}

    MyTree = filtered_tree->CloneTree();

    TH1F* hist_eph = new TH1F("hist_eph", "eph;Mev;Events", 50, 0., 9.);
    
    auto canvas = new TCanvas("task_7", "task_7");
    hist_eph->SetXTitle("E_{#gamma}, MeV");
    hist_eph->SetYTitle("Num");
    canvas->SetLogy();

    TF1 *fitting = new TF1("fitting","fit_func", 0, 9, 3);
    MyTree->Fit("fitting","eph");

    MyTree->GetHistogram()->GetXaxis()->SetTitle("E_{#gamma}, GeV");
    MyTree->GetHistogram()->GetYaxis()->SetTitle("entries");
    MyTree->GetHistogram()->SetTitle("ph enegy");
    gPad->SetLogy();
    
    //MyTree->Draw("eph");
    //fitting->Draw("Same");
    
    ofile->cd();
    MyTree->Write();
    //canvas->Write();
    double max = MyTree->GetMaximum("eph");
    double min = MyTree->GetMinimum("eph");
    std::cout <<"maximum: "<< max << " minimum: " << min <<endl;

    ofile->Close();
    ifile.Close();

}