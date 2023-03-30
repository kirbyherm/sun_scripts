// Kirby Hermansen 
// 2022-10-19
// Purpose: To create a histogram of SuN's TAS from geant sim

#include "TFile.h"
#include "TChain.h"
#include "TH1F.h"
#ifndef UTILS_H
#define UTILS_H
#include "utils.h"
#endif

#include <iostream>
#include <string>
#include <vector>
#include <fstream>

void tas_sim(int thread=0){
    extern const ConfigParams p;
    std::vector<std::string> templates = p.Get_Levels();
    std::vector<std::string> filelist_samp = templates;
    std::string isotope = p.Get_Isotope();
    int tas_lower_limit = p.Get_TAS_Window()[0];
    int tas_upper_limit = p.Get_TAS_Window()[1];
    int bin_width = 1;
    for (int j = 0; j<filelist_samp.size(); j++) {
        std::string energy=filelist_samp[j];
        
        int energy_int = std::stoi(energy);
        TFile * f1;
        f1 = TFile::Open(Form("output_%s_%s.root",isotope.c_str(),energy.c_str()));
        TTree *t;
        f1->GetObject("t",t); 
        TH1F * h1;
        std::string cuts = p.Get_simCuts();
        t->Draw(Form("eneAll>>TAS(%d,%d,%d)",int((tas_upper_limit-tas_lower_limit)/bin_width),tas_lower_limit,tas_upper_limit),Form("%s",cuts.c_str()),"");
        h1 = ((TH1F* ) (gDirectory->Get("TAS"))); 
        if( thread <1){
            std::ofstream myFile(Form("tas_new/output_%s.csv",energy.c_str()));
            for ( int i = 1; i<h1->GetNbinsX()+1; i++){
                myFile << i << "," << h1->GetBinContent(i) << "\n";
            }
            myFile.close();
        }else{
            std::cout << h1->GetNbinsX() << ", " << h1->Integral() << std::endl;
        std::cout << thread << std::endl;
            std::ofstream myFile(Form("%d/tas_new/output_%s.csv",thread,energy.c_str()));
            for ( int i = 1; i<h1->GetNbinsX()+1; i++){
                myFile << i << "," << h1->GetBinContent(i) << "\n";
            }
            myFile.close();
        }
    }
    return;
}

int main(int argc, char **argv){
    int thread = 0;
    if( argc > 1){
        std::cout << thread << std::endl;
        thread = std::stoi(argv[1]);
        std::cout << thread << std::endl;
    }
    tas_sim(thread);
    return 0;
}
 

