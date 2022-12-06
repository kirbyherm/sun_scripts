#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include <iostream>
#include <string>
#include <vector>
#include <fstream>

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
//    TFile *f1 = TFile::Open("new_bgs_segs.root");
//    TFile *f0 = TFile::Open("output_60Co_center.root");
//    TFile *f0 = TFile::Open("output_60Co_+1.2z.root");
//    TFile *f1 = TFile::Open("bgs_segs_3005.root");
    TTree* t;
    f0->GetObject("t",t); 
    std::vector<std::string> side = {"B","T"}; 
    std::vector<TH1F*> co60_h;
    std::vector<TH1F*> sim_h;
    TH1F * sumseg_sim = new TH1F(Form("E%s",histname.c_str()),"simulated sum of seg",10000,1,10001);
//    TH1F * sumseg_co60 = new TH1F("sumSegCo60","Co60 sum of seg",4000,1,4001);
    double total_co60 = 0;
    double total_sim = 0;
    for(int i =0; i<side.size(); i++){
        for(int j=1; j<5; j++){
//    t->Draw(Form("ene%s%d>>simSeg_%s_%s%d(12000,1,12001)",side[i].c_str(),j,histname.c_str(),side[i].c_str(),j),"","");
    t->Draw(Form("ene%s%d>>simSeg_%s_%s%d(10000,1,10001)",side[i].c_str(),j,histname.c_str(),side[i].c_str(),j),"eneDSSD>0","");
//    t->Draw(Form("ene%s%d>>simSeg_%s_%s%d(12000,1,12001)",side[i].c_str(),j,histname.c_str(),side[i].c_str(),j),"eneAll>140&&eneAll<220&&eneDSSD>0","");
//    t->Draw(Form("ene%s%d>>simSeg_%s_%s%d(12000,1,12001)",side[i].c_str(),j,histname.c_str(),side[i].c_str(),j),"eneAll>1650&&eneAll<1850&&eneDSSD>0","");
//    t->Draw(Form("ene%s%d>>simSeg_%s_%s%d(12000,1,12001)",side[i].c_str(),j,histname.c_str(),side[i].c_str(),j),"eneAll>1915&&eneAll<2115&&eneDSSD>0","");
//    t->Draw(Form("ene%s%d>>simSeg_%s_%s%d(12000,1,12001)",side[i].c_str(),j,histname.c_str(),side[i].c_str(),j),"eneAll>2300&&eneAll<2600&&eneDSSD>0","");
    TH1F *simSeg = ((TH1F* ) (gDirectory->Get(Form("simSeg_%s_%s%d",histname.c_str(),side[i].c_str(),j)))); 
    simSeg->Draw();
//    TH1F *co60Seg;
//    co60Seg = (TH1F*)(f1->Get(Form("co60bgs_E_SuN%s%d",side[i].c_str(),j)));    
    simSeg->SetLineColor(2);
    sim_h.push_back(simSeg);
//    co60_h.push_back(co60Seg);
    if (j==2 || j==3)
//    std::cout<<co60Seg->Integral(100,4000)<<std::endl;
//    std::cout<<co60Seg->Integral()/simSeg->Integral()<<std::endl;
//    total_co60 += co60Seg->Integral();
//    total_sim += simSeg->Integral();
    {
        sumseg_sim->Add(sim_h.back());
//        sumseg_co60->Add(co60_h.back());
    }
        }
    }
    std::cout<<histname.c_str()<<" "<<sumseg_sim->Integral(450,560)/sumseg_sim->Integral()<<std::endl;
//    std::cout<<total_co60/total_sim<<std::endl;
//    std::cout<<co60_h[1]->Integral()/co60_h[2]->Integral()<<std::endl;
//    std::cout<<co60_h[5]->Integral()/co60_h[6]->Integral()<<std::endl;
//    std::cout<<sim_h[1]->Integral()/sim_h[2]->Integral()<<std::endl;
//    std::cout<<sim_h[5]->Integral()/sim_h[6]->Integral()<<std::endl;
    for(int i =0; i<side.size(); i++){
        for(int j=1; j<5; j++){
    TCanvas *c1 = new TCanvas();
//    sim_h[i*4+j-1]->Scale(total_co60/total_sim);
//    sim_h[i*4+j-1]->GetYaxis()->SetRangeUser(0,std::max(sim_h[i*4+j-1]->GetMaximum(),co60_h[i*4+j-1]->GetMaximum())+10.0);
    sim_h[i*4+j-1]->Draw("hist");
//    co60_h[i*4+j-1]->Draw("hist same");
//    save_to_png(c1,Form("E_SuN%s%d",side[i].c_str(),j));
        }
    }
    TCanvas *c2 = new TCanvas();
//    sumseg_sim->Scale(total_co60/total_sim);
    sumseg_sim->SetLineColor(2);
    sumseg_sim->Draw("hist");
    sumseg_sim->GetXaxis()->SetRangeUser(1,2500);
//    sumseg_co60->Draw("same hist");
//    std::string savefilename = "sumSeg" + filename.substr(3,filename.size());
//    TLine *l1 = new TLine(560,0,560,20000);
//    l1->Draw("same");
//    TLine *l2 = new TLine(450,0,450,20000);
//    l2->Draw("same");
//    save_to_png(c2,savefilename.c_str());

    return sumseg_sim;
}

void ss_sim(int thread=0){


    std::vector<std::string> filelist = {"2476","2500","2600","2700","2800","2900","3000","3150","3300","3450","3600","3750","3900","4100","4300","4500","4700","4900","5100","5300","5500","5700","5900","6100","6300","6500","6700","6900","7100","7400","7700","8000","8300","8600","8900","9200","9500","9800","10100","10400"};//,"Test"};
//    std::vector<std::string> filelist = {"0","2476","2500","2600","2700","2800","2900","3000","3150","3300","3450","3650","3850","4050","4250","4450","4650","4850","5050","5250","5500","5700","5900","6100","6300","6500","6700","6900","7100","8000","8300","8600","8900","9200","9500","9800","10100","10400"};//,"Test"};
////    std::vector<std::string> filelist = {"Test1000","Test1400","Test2000","Test5000","Test8000"};
    TFile *f = new TFile("RAINIER_ss.root", "RECREATE");
    for (int i = 0; i<filelist.size(); i++) {

        std::string filename = "Enew"+filelist[i]+".root";
        TH1F * sumseg_x = new TH1F();
        sumseg_x = sumSeg(filename, filelist[i]);
//        sumseg_x->Rebin(4);
        write_to_csv_ss(thread, filelist[i], sumseg_x);
        f->cd();
        sumseg_x->Write();
    }
//    std::vector<std::string> filelist_samp = {"0","61","113","174","1690","1731","1754","2036","2178","2289","2475"};
    std::vector<std::string> filelist_samp = {"0","61","113","174","744","1200","1240","1690","1731","1750","1754","2000","2036","2178","1700500","2289","2475"};
//    std::vector<std::string> filelist_samp = {"0","174_m1","174_m2","1731_m1","1731_m2","1731_m3","1754_m1","1754_m2","1754_m3","2036_m1","2036_m2","2036_m3","2036_m4","2036_m5","2475_m1","2475_m2","2475_m3"};
//    std::vector<std::string> filelist_samp = {"174_m1","174_m2"};
//    std::vector<std::string> filelist_samp = {"1731","1754"};//"174",,"2036","2475"};
//    std::vector<std::string> filelist_samp = {"2036"};//"174",,"2036","2475"};
//    std::vector<std::string> filelist_samp = {"2475"};//"174",,"2036","2475"};
    for (int i = 0; i<filelist_samp.size(); i++) {
//        std::vector<std::string> mult_samp = {"m1","m2"};
        std::vector<std::string> mult_samp = {"m1","m2","m3"};
        mult_samp = {""};
        for (int j = 0; j<mult_samp.size(); j++) {
         
//        std::string filename = "../output_57Ti_"+filelist_samp[i]+"_" + mult_samp[j]+".root";
        std::string filename = "output_57Ti_"+filelist_samp[i]+".root";
        TH1F * sumseg_x = new TH1F();
        sumseg_x = sumSeg(filename, filelist_samp[i]);
//        sumseg_x->Rebin(4);
//        int j = 0;
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

