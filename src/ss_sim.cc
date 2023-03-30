#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#ifndef UTILS_H
#define UTILS_H
#include "utils.h"
#endif

extern const ConfigParams p;

void write_to_csv_ss(int thread=-1, std::string energy="2475", TH1F* h1 = new TH1F()){

        if( thread <1){
            std::ofstream myFile(Form("ss_new/Enew%s.csv",energy.c_str()));
            for ( int i = 1; i<h1->GetNbinsX()+1; i++){
                myFile << i << "," << h1->GetBinContent(i) << "\n";
            }
            myFile.close();
        }else{
            std::ofstream myFile(Form("%d/ss_new/Enew%s.csv",thread,energy.c_str()));
            for ( int i = 1; i<h1->GetNbinsX()+1; i++){
                myFile << i << "," << h1->GetBinContent(i) << "\n";
            }
            myFile.close();
        }

return;
}
 
TH1F * sumSeg(std::string filename, std::string histname){
    TFile *f0 = TFile::Open(filename.c_str());
    TTree* t;
    f0->GetObject("t",t); 
    std::vector<std::string> side = {"B","T"}; 
    std::vector<TH1F*> co60_h;
    std::vector<TH1F*> sim_h;
    int tas_lower_limit = p.Get_TAS_Window()[0];
    int tas_upper_limit = p.Get_TAS_Window()[1];
    int ss_lower_limit = p.Get_SS_Window()[0];
    int ss_upper_limit = p.Get_SS_Window()[1];
    int bin_width = 1;
    TH1F * sumseg_sim = new TH1F(Form("E%s",histname.c_str()),"simulated sum of seg",int((ss_upper_limit-ss_lower_limit)/bin_width),ss_lower_limit,ss_upper_limit);
//    TH1F * sumseg_co60 = new TH1F("sumSegCo60","Co60 sum of seg",4000,1,4001);
    double total_co60 = 0;
    double total_sim = 0;
    for(int i =0; i<side.size(); i++){
        for(int j=1; j<5; j++){
            if (p.Use_Middle_Segments_Only() && (j==1 || j==4))
                continue;
    std::string cuts = p.Get_simCuts();
    t->Draw(Form("ene%s%d>>simSeg_%s_%s%d(%d,%d,%d)",side[i].c_str(),j,histname.c_str(),side[i].c_str(),j,int((ss_upper_limit-ss_lower_limit)/bin_width),ss_lower_limit,ss_upper_limit),Form("%s",cuts.c_str()),"");
    TH1F *simSeg = ((TH1F* ) (gDirectory->Get(Form("simSeg_%s_%s%d",histname.c_str(),side[i].c_str(),j)))); 
    simSeg->Draw();
    simSeg->SetLineColor(2);
    sim_h.push_back(simSeg);
    sumseg_sim->Add(sim_h.back());
        }
    }
    std::cout<<histname.c_str()<<" "<<sumseg_sim->Integral(450,560)/sumseg_sim->Integral()<<std::endl;
//    for(int i =0; i<side.size(); i++){
//        for(int j=1; j<5; j++){
//    TCanvas *c1 = new TCanvas();
//    sim_h[i*4+j-1]->Draw("hist");
//        }
//    }
    TCanvas *c2 = new TCanvas();
    sumseg_sim->SetLineColor(2);
    sumseg_sim->Draw("hist");
    sumseg_sim->GetXaxis()->SetRangeUser(1,2500);

    return sumseg_sim;
}

void ss_sim(int thread=0){


    std::vector<std::string> templates = p.Get_Levels();
    std::string isotope = p.Get_Isotope();
    TFile *f = new TFile("RAINIER_ss.root", "RECREATE");
    std::vector<std::string> filelist_samp = templates;//"174",,"2036","2475"};

    for (int i = 0; i<filelist_samp.size(); i++) {
        std::vector<std::string> mult_samp = {"m1","m2","m3"};
        mult_samp = {""};
        for (int j = 0; j<mult_samp.size(); j++) {
         
        std::string filename = "output_"+isotope+"_"+filelist_samp[i]+".root";
        TH1F * sumseg_x = new TH1F();
        sumseg_x = sumSeg(filename, filelist_samp[i]);
        write_to_csv_ss(thread, filelist_samp[i]+mult_samp[j], sumseg_x);
        f->cd();
        sumseg_x->Write();
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
    ss_sim(thread);
    return 0;
}

