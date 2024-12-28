int number = 1000000;

double energy = 1020;
double Pi_mass = 139.57;
double K_mass = 497.61;

double radius_detector = 30;
double len_detector = 50;

double min_energy = 40;
double C_TAU = 2.69;

void task_6()
{

    TRandom3 *rnd = new TRandom3();
    rnd->SetSeed(10);

    TH1D *Ks_Theta = new TH1D("", "", 100, 0, TMath::Pi());
    Ks_Theta->SetTitle("distirbution of angle #theta_{Ks} lab;#theta_{Ks}, [rad]");
    
    TH1D *Ks_Phi   = new TH1D("", "", 100, 0, TMath::Pi());
    Ks_Phi->SetTitle("distirbution of angle #phi_{Ks} lab;#phi_{Ks}, [rad]");

    TH1D *Pi_Theta = new TH1D("", "", 100, 0, TMath::Pi());
    Pi_Theta->SetTitle("distribution of angle #theta_{#pi} lab;#theta_{Pi_{#pm}}, [rad]");
    
    TH1D *Pi_Phi = new TH1D("", "", 100, 0, TMath::Pi());
    Pi_Phi->SetTitle("distribution of angle #phi_{#pi} lab;#phi_{Pi_{#pm}}, [rad]");
    
    TH1D *Ks_leng = new TH1D("", "", 100, 0, 20);
    Ks_leng->SetTitle("distirbution of len Ks;length, [cm]");

    TF1 *KsThetaDist = new TF1("KsTheta", "sin(x)*sin(x)*sin(x)", 0, TMath::Pi());

    for (int i = 0; i < number; i++) {
        
        double Ks_theta = KsThetaDist->GetRandom();
        double Ks_phi   = 2 * TMath::Pi() * rnd->Rndm() - TMath::Pi();

        double Ks_cosTheta = TMath::Cos(Ks_theta);
        double Ks_sinTheta = TMath::Sin(Ks_theta);
        double Ks_cosPhi = TMath::Cos(Ks_phi);
        double Ks_sinPhi = TMath::Sin(Ks_phi);
            
        double Ks_energy = energy/2;
        double Ks_momentum = TMath::Sqrt(pow(Ks_energy,2) - pow(K_mass, 2));

        TLorentzVector Ks(0, 0, 0, K_mass);

        double KsBx = Ks_momentum*Ks_sinTheta*Ks_cosPhi / Ks_energy;
        double KsBy = Ks_momentum*Ks_sinTheta*Ks_sinPhi / Ks_energy;
        double KsBz = Ks_momentum*Ks_cosTheta/ Ks_energy;

        double length = - TMath::Log(rnd->Rndm()) * C_TAU * Ks_momentum/Ks_energy; 
        TVector3 Decay(length*Ks_sinTheta*Ks_cosPhi, length*Ks_sinTheta*Ks_sinPhi,length*Ks_cosTheta);

        if (TMath::Abs(Decay.Z()) < len_detector && TMath::Abs(Decay.X()) < radius_detector) {
            double Pi_cosTheta = 2 * rnd->Rndm() - 1;
            double Pi_sinTheta = TMath::Sqrt(1 - pow(Pi_cosTheta,2));

            double Pi_phi   = 2 * TMath::Pi() * rnd->Rndm() - TMath::Pi();
            
            double Pi_cosPhi = TMath::Cos(Pi_phi);
            double Pi_sinPhi = TMath::Sin(Pi_phi);
            
            double Pi_energy = (K_mass - Pi_mass) * rnd->Rndm() + Pi_mass;
            double Pi_Momentum = TMath::Sqrt(Pi_energy*Pi_energy - pow(Pi_mass, 2));

            double Pi_Px = Pi_Momentum*Pi_sinTheta*Pi_cosPhi;
            double Pi_Py = Pi_Momentum*Pi_sinTheta*Pi_sinPhi;
            double Pi_Pz = Pi_Momentum*Pi_cosTheta;

            TLorentzVector Pi_Plus(Pi_Px, Pi_Py, Pi_Pz, Pi_energy);
            TLorentzVector Pi_Minus(-Pi_Px, -Pi_Py, -Pi_Pz, Pi_energy);

            TLorentzRotation T;
            T.Boost(KsBx, KsBy, KsBz);

            Ks = T.VectorMultiplication(Ks);
            Pi_Plus = T.VectorMultiplication(Pi_Plus);
            Pi_Minus = T.VectorMultiplication(Pi_Minus);

            if (Pi_Plus.P() < min_energy || Pi_Minus.P() < min_energy) continue;
           
            Ks_Theta->Fill(Ks_theta);
            Ks_Phi->Fill(Ks_phi);

            Pi_Theta->Fill(Pi_Plus.Theta());
            Pi_Theta->Fill(Pi_Minus.Theta());

            Pi_Phi->Fill(Pi_Plus.Phi());
            Pi_Phi->Fill(Pi_Minus.Phi());

            Ks_leng->Fill(length);

        }
    }

    TCanvas *c = new TCanvas("", "");
    c->Divide(3, 2);

    c->cd(1);
    Ks_Theta->Draw();
    
    c->cd(2);
    Pi_Theta->Draw();

    c->cd(3);
    Ks_leng->Draw();

    c->cd(4);
    Ks_Phi->Draw();    
    
    c->cd(5);
    Pi_Phi->Draw();

}
