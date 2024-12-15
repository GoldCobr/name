#include <TTree.h>
#include <TFile.h>

void task_9(){
    TFile file("newroot.root", "read");
    TFile outfile("ofile_task8.root","read");
    TTree *new_tree =(TTree *)file.Get("h10");
    //TTree* main_tree = (TTree*)ifile2.Get("ofile_task8");
    
    int nph; float eph[45]; float thetaph[45]; float phiph[45]; 
    std::vector<float> thetaph_graph(18,0);
    std::vector<float> phiph_graph(36,0);
        
    new_tree->SetBranchAddress("nph",&nph);
    new_tree->SetBranchAddress("eph",eph);
    new_tree->SetBranchAddress("thetaph",thetaph);
    new_tree->SetBranchAddress("phiph",phiph);

    for (int i = 0; i < new_tree->GetEntries(); i++) {
        new_tree->GetEntry(i);
        for (int j = 0; j < nph - 1; j++) {
            for (int k = 1; k < nph; k++) {
                TLorentzVector vector1(eph[j], 0, 0, eph[j]);
                vector1.SetTheta(thetaph[j]);
                vector1.SetPhi(phiph[j]);
                TLorentzVector vector2(eph[k], 0, 0, eph[k]);
                vector2.SetTheta(thetaph[k]);
                vector2.SetPhi(phiph[k]);
                TLorentzVector vector_sum= vector1+vector2;

                TVector3 v1(sin(thetaph[j])*cos(phiph[j]), sin(thetaph[j])*sin(phiph[j]), cos(phiph[j]));
                TVector3 v2(sin(thetaph[k])*cos(phiph[k]), sin(thetaph[k])*sin(phiph[k]), cos(phiph[k]));
                auto ang = v1.Angle(v2);

                double m_inv = sqrt(2. * eph[j] * eph[k] * (1. - cos(ang)));
                if (0.1 <= m_inv && m_inv <= 0.2){
                    int th = vector_sum.Theta()*18/TMath::Pi();
                    int ph = (vector_sum.Phi()+TMath::Pi())*18/TMath::Pi();
                    thetaph_graph[th]=thetaph_graph[th]+1;
                    phiph_graph[ph]=phiph_graph[ph]+1;
                }
            }
        }
    }
    auto c1 = new TCanvas("c1","c1",600,500);

    

    auto graph1 = new TGraphErrors();
    graph1->SetTitle("dependence on #theta");
    for(int i=0; i<18; i++){
        graph1->SetPoint(i, i*10, thetaph_graph[i]);
        graph1->SetPointError(i,0,sqrt(thetaph_graph[i]));
    }
    graph1->SetMarkerColor(kRed);
    graph1->SetMarkerStyle(21);

    auto graph2 = new TGraphErrors();
    graph2->SetTitle("dependence on #phi");
    for(int i = 0; i<36; i++ ){
        graph2->SetPoint(i, i*10, phiph_graph[i]);
        graph2->SetPointError(i,0,sqrt(phiph_graph[i]));
    }
    graph2->SetMarkerColor(kBlue);
    graph2->SetMarkerStyle(21);

    auto my_graph = new TMultiGraph;
    my_graph->SetTitle("dependence of the number of candidates for #pi^{0} on the angle and polar azimuthal angle;angle [degree];candidates");
    my_graph->Add(graph1);
    my_graph->Add(graph2);
    my_graph->Draw("ap");

    c1->BuildLegend(0.48,0.8,0.9,0.9);
    TDirectory* subdir = outfile.mkdir("subdir","subdir title");
    subdir->cd();
    my_graph->Write();
    subdir->pwd();
}