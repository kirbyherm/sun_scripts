// Kirby Hermansen 
// 2022-10-19
// Purpose: To create a histogram of SuN's multiplicity from geant sim 

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

void mult_sim(int thread=0){
    extern const ConfigParams p;
    std::vector<std::string> templates = p.Get_Levels();
    std::vector<std::string> filelist_samp = templates;
    std::string cuts = p.Get_simCuts();
    std::string isotope = p.Get_Isotope();
    int tas_lower_limit = p.Get_simTAS_Window()[0];
    int tas_upper_limit = p.Get_simTAS_Window()[1];
    for (int j = 0; j<filelist_samp.size(); j++) {
        std::string energy=filelist_samp[j];
        
        int energy_int = std::stoi(energy);
        TFile * f1;
        f1 = TFile::Open(Form("output_%s_%s.root",isotope.c_str(),energy.c_str()));
        TTree *t;
        f1->GetObject("t",t); 
        TH1F * h1;
        t->Draw("multi>>mult(10,0,10)",Form("%s",cuts.c_str()),"");
        h1 = ((TH1F* ) (gDirectory->Get("mult"))); 
        if( thread <1){
            std::ofstream myFile(Form("mult_new/output_%s.csv",energy.c_str()));
            for ( int i = 1; i<h1->GetNbinsX()+1; i++){
                myFile << i << "," << h1->GetBinContent(i) << "\n";
            }
            myFile.close();
        }else{
            std::ofstream myFile(Form("%d/mult_new/output_%s.csv",thread,energy.c_str()));
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
    mult_sim(thread);
    return 0;
}
 

