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
//    std::vector<std::string> filelist_samp = {"0","364","1100","1500","1600","2000","2200","2500","2600","2700","2800","2900","3000","3150","3300","3450","3600","3750","3900","4100","4300","4500","4700","4900","5100","5300","5500","5700","5900","6100","6300","6500","6700","6900","7100","7400","7700","8000","8300","8600","8900","9200","9500","9800","10100","10400","10800","11200","11600","12000","12400"};
    //std::vector<std::string> filelist_samp = {"0","364"};
    extern const ConfigParams p;
    std::vector<std::string> templates = p.Get_Levels();
    std::vector<std::string> filelist_samp = templates;
    //std::vector<std::string> filelist_samp = {"0","364","1100","1500_m1","1500_m2"};
    //std::vector<std::string> filelist_samp = {"0","364","1100","1500","1600"};
    //std::vector<std::string> filelist_samp = {"0","364","1100","1500","1600","2000_m1","2000_m2"};
    //std::vector<std::string> filelist_samp = {"0","364","1100","1500","1600","2000","2200_m1","2200_m2"};
    //std::vector<std::string> filelist_samp = {"0","364","1100","1500","1600","2000","2200","2500_m1","2500_m2"};
    //std::vector<std::string> filelist_samp = {"0","364","1100","1500","1600","2000","2200","2500","2600"};
    std::string isotope = p.Get_Isotope();
    int tas_lower_limit = p.Get_simTAS_Window()[0];
    int tas_upper_limit = p.Get_simTAS_Window()[1];
    int bin_width = 1;
    for (int j = 0; j<filelist_samp.size(); j++) {
        std::string energy=filelist_samp[j];
        
        int energy_int = std::stoi(energy);
        TFile * f1;
        f1 = TFile::Open(Form("output_%s_%s.root",isotope.c_str(),energy.c_str()));
        TTree *t;
        f1->GetObject("t",t); 
        TH1F * h1;
        //t->Draw("eneAll>>TAS(10,0,10)","eneAll>140&&eneAll<220","");
        //t->Draw("eneAll>>TAS(10000,1,10001)","","");
//        t->Draw("eneAll>>TAS(13000,1,13001)","eneDSSD>0","");
        t->Draw(Form("eneAll>>TAS(%d,%d,%d)",int((tas_upper_limit-tas_lower_limit)/bin_width),tas_lower_limit,tas_upper_limit),Form("eneAll>=%d&&eneAll<%d&&eneDSSD>0",tas_lower_limit, tas_upper_limit),"");
        //t->Draw("eneAll>>TAS(13000,1,13001)","eneAll>900&&eneAll<1300&&eneDSSD>0","");
        //t->Draw("eneAll>>TAS(13000,1,13001)","eneAll>1200&&eneAll<1700&&eneDSSD>0","");
        //t->Draw("eneAll>>TAS(13000,1,13001)","eneAll>1550&&eneAll<1760&&eneDSSD>0","");
        //t->Draw("eneAll>>TAS(13000,1,13001)","eneAll>1800&&eneAll<2050&&eneDSSD>0","");
        //t->Draw("eneAll>>TAS(13000,1,13001)","eneAll>2040&&eneAll<2320&&eneDSSD>0","");
        //t->Draw("eneAll>>TAS(13000,1,13001)","eneAll>2200&&eneAll<2500&&eneDSSD>0","");
        //t->Draw("eneAll>>TAS(13000,1,13001)","eneAll>2500&&eneAll<2700&&eneDSSD>0","");
        //t->Draw("eneAll>>TAS(13000,1,13001)","eneAll>2650&&eneAll<10001&&eneDSSD>0","");
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
 

