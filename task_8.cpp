#include <TTree.h>
#include <TFile.h>

void task_8() {
    TFile ifile("newroot.root", "read");
    TTree* MyTree = (TTree*)ifile.Get("h10");

    int nph; float eph[45]; float thetaph[45]; float phiph[45]; int p0_candidate; std::vector<double> inv_mass_vec;

    MyTree->SetBranchAddress("nph",&nph);
    MyTree->SetBranchAddress("eph",eph);
    MyTree->SetBranchAddress("thetaph",thetaph);
    MyTree->SetBranchAddress("phiph",phiph);

    std::vector<int> count_candidates;
    TH1D* h_inv_mass = new TH1D("h_inv_mass","Inv Mass; Mass, [Gev];enteries", 40, 0., 0.3);
    TH1D* h_angle = new TH1D("h_angle", "Angle between photons;angle, [rad];enteries", 50, 0, TMath::Pi()+0.2);

    //graph_theta = ROOT.TGraph();
    //graph_phi = ROOT.TGraph();

    for (int i = 0; i < MyTree->GetEntries(); i++) {
        MyTree->GetEntry(i);
        p0_candidate = 0;
        inv_mass_vec.erase(inv_mass_vec.begin(), inv_mass_vec.end());
        for (int j = 0; j < nph - 1; j++) {
            for (int k = j; k < nph - 1; k++) {
                //расчет угла через tvector3
                TVector3 v1(sin(thetaph[j])*cos(phiph[j]), sin(thetaph[j])*sin(phiph[j]), cos(phiph[j]));
                TVector3 v2(sin(thetaph[k])*cos(phiph[k]), sin(thetaph[k])*sin(phiph[k]), cos(phiph[k]));
                auto ang = v1.Angle(v2);
                double m_inv = sqrt(2. * eph[j] * eph[k] * (1. - cos(ang)));
                h_angle->Fill(ang);
                if (0.1 <= m_inv && m_inv <= 0.2){
                    p0_candidate++;
                    inv_mass_vec.push_back(m_inv);
                }
            }
        }
        count_candidates.push_back(p0_candidate);
        if (p0_candidate == 2){
            for (double val : inv_mass_vec){
                h_inv_mass->Fill(val);
            }
        } 
    }
    TFile* ofile = new TFile("ofile_task8.root", "RECREATE");
    h_inv_mass->Write();
    h_angle->Write();
    ofile->Close();
}