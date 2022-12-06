

#include "TFile.h"
#include "TChain.h"
#include "TH1F.h"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>

void mult_sim(int thread=0){
    std::vector<std::string> filelist_samp = {"0","61","113","174","744","1200","1240","1690","1731","1750","1754","2000","2036","2178","1700500","2289","2475","2476","2500","2600","2700","2800","2900","3000","3150","3300","3450","3600","3750","3900","4100","4300","4500","4700","4900","5100","5300","5500","5700","5900","6100","6300","6500","6700","6900","7100","7400","7700","8000","8300","8600","8900","9200","9500","9800","10100","10400"};//,"Test"};
    //std::vector<std::string> filelist_samp = {"0","364"};
    //std::vector<std::string> filelist_samp = {"0","364","1100_m1","1100_m2"};
    //std::vector<std::string> filelist_samp = {"0","364","1100","1500_m1","1500_m2"};
    //std::vector<std::string> filelist_samp = {"0","364","1100","1500","1600"};
    //std::vector<std::string> filelist_samp = {"0","364","1100","1500","1600","2000_m1","2000_m2"};
    //std::vector<std::string> filelist_samp = {"0","364","1100","1500","1600","2000","2200_m1","2200_m2"};
    //std::vector<std::string> filelist_samp = {"0","364","1100","1500","1600","2000","2200","2500_m1","2500_m2"};
    //std::vector<std::string> filelist_samp = {"0","364","1100","1500","1600","2000","2200","2500","2600"};
    for (int j = 0; j<filelist_samp.size(); j++) {
        std::string energy=filelist_samp[j];
        
        int energy_int = std::stoi(energy);
        TFile * f1;
        if (energy_int < 2476 || energy_int == 1700500)
        f1 = TFile::Open(Form("output_57Ti_%s.root",energy.c_str()));
        else
        f1 = TFile::Open(Form("Enew%s.root",energy.c_str()));
        TTree *t;
        f1->GetObject("t",t); 
        TH1F * h1;
        t->Draw("multi>>mult(10,0,10)","eneDSSD>0","");
        //t->Draw("multi>>mult(10,0,10)","eneAll>320&&eneAll<415&&eneDSSD>0","");
        //t->Draw("multi>>mult(10,0,10)","eneAll>900&&eneAll<1300&&eneDSSD>0","");
        //t->Draw("multi>>mult(10,0,10)","eneAll>1200&&eneAll<1700&&eneDSSD>0","");
        //t->Draw("multi>>mult(10,0,10)","eneAll>1550&&eneAll<1760&&eneDSSD>0","");
        //t->Draw("multi>>mult(10,0,10)","eneAll>1800&&eneAll<2050&&eneDSSD>0","");
        //t->Draw("multi>>mult(10,0,10)","eneAll>2040&&eneAll<2320&&eneDSSD>0","");
        //t->Draw("multi>>mult(10,0,10)","eneAll>2200&&eneAll<2500&&eneDSSD>0","");
        //t->Draw("multi>>mult(10,0,10)","eneAll>2500&&eneAll<2700&&eneDSSD>0","");
        //t->Draw("multi>>mult(10,0,10)","eneAll>2650&&eneAll<10001&&eneDSSD>0","");
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
 

