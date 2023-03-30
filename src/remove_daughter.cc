#include <vector>
#ifndef UTILS_H
#define UTILS_H
#include "utils.h"
#endif

#include "TH2F.h"

extern const ConfigParams p;
int tas_bin_width = p.Get_TAS_Bin_Width();
int tas_lower_limit = p.Get_TAS_Window()[0];
int tas_upper_limit = p.Get_TAS_Window()[1];
int ss_bin_width = p.Get_SS_Bin_Width();
int ss_lower_limit = p.Get_SS_Window()[0];
int ss_upper_limit = p.Get_SS_Window()[1];
double scale_factor = p.Get_Scale_Factor();
std::string isotope = p.Get_Isotope();
std::string isotope_daughter = p.Get_Isotope_Daughter();

void remove_daughter_tas_ss(){

    TFile *f0 = new TFile("output_tas_ss_parent.root","READ");
    TFile *f1 = new TFile("output_tas_ss_daughter.root","READ");

    TH2F* parent_bgs = ((TH2F*) f0->Get(Form("%s_tas_ss_bgs",isotope.c_str()))); 
    TH2F* daughter_bgs = ((TH2F*) f1->Get(Form("%s_tas_ss_bgs",isotope_daughter.c_str()))); 
    
    TFile *f2 = new TFile(Form("%s_tas_ss_final.root",isotope.c_str()),"RECREATE");
    parent_bgs->Write();
    daughter_bgs->Write();
    daughter_bgs->Scale(scale_factor);
    parent_bgs->Add(daughter_bgs, -1);

    TH2F * parent_out = new TH2F(Form("%s_tas_ss_final",isotope.c_str()),Form("%s_tas_ss_final",isotope.c_str()),int((ss_upper_limit-ss_lower_limit)/ss_bin_width),ss_lower_limit,ss_upper_limit,int((tas_upper_limit-tas_lower_limit)/tas_bin_width),tas_lower_limit,tas_upper_limit);    
    parent_out->Add(parent_bgs);
    for (int i = 1; i < parent_out->GetNbinsX()+1; i++){
        for (int j = 1; j < parent_out->GetNbinsY()+1; j++){
//            parent_out->SetBinContent(i,j, parent_bgs->GetBinContent(i,j));
//            parent_out->SetBinError(i,j, parent_bgs->GetBinError(i,j));
//            if (parent_out->GetBinError(i,j) < 1)
//            parent_out->SetBinError(i,j, 1);
//            parent_out->SetBinContent(i,j, parent_bgs->GetBinContent(i,j));
//            parent_out->SetBinError(i,j, parent_bgs->GetBinError(i,j));
//            if (parent_out->GetBinError(i,j) < 1)
//            parent_out->SetBinError(i,j, 1);
        }
    }
    parent_out->Write();
    f2->Close();


    return;
}

void remove_daughter_ss(){

    TFile *f0 = new TFile("output_ss_parent.root","READ");
    TFile *f1 = new TFile("output_ss_daughter.root","READ");

    TH1F* parent_bgs = ((TH1F*) f0->Get(Form("%s_ss_bgs",isotope.c_str()))); 
    TH1F* daughter_bgs = ((TH1F*) f1->Get(Form("%s_ss_bgs",isotope_daughter.c_str()))); 
    
    daughter_bgs->SetLineColor(kRed);
    TFile *f2 = new TFile(Form("%s_ss_final.root",isotope.c_str()),"RECREATE");
    parent_bgs->Write();
    daughter_bgs->Write();
//    add_bg_error(parent_bgs, daughter_bgs, scale_factor);
    daughter_bgs->Scale(scale_factor);
    parent_bgs->Add(daughter_bgs, -1);

    TH1F * parent_out = new TH1F(Form("%s_ss_final",isotope.c_str()),Form("%s_ss_final",isotope.c_str()),int((ss_upper_limit-ss_lower_limit)/ss_bin_width),ss_lower_limit,ss_upper_limit);    
    parent_out->Add(parent_bgs);
    TH1F* parent_error = ((TH1F*) f0->Get(Form("%s_ss_bgs",isotope.c_str()))); 
    TH1F* daughter_error = ((TH1F*) f1->Get(Form("%s_ss_bgs",isotope_daughter.c_str()))); 
    daughter_error->Scale(scale_factor);
    for (int i = 1; i < parent_out->GetNbinsX()+1; i++){
        parent_out->SetBinContent(i, parent_bgs->GetBinContent(i));
        parent_out->SetBinError(i, parent_bgs->GetBinError(i));
//        std::cout << parent_error->GetBinContent(i) << ", " << daughter_error->GetBinContent(i) << ", " << parent_out->GetBinContent(i)<< std::endl;
//        parent_out->SetBinError(i, (sqrt(pow(parent_error->GetBinError(i)/std::max(parent_error->GetBinContent(i),1.0),2)+pow(daughter_bgs->GetBinError(i)/std::max(daughter_bgs->GetBinContent(i),1.0),2)))*(std::max(parent_out->GetBinContent(i),1.0)));
        if (parent_out->GetBinError(i) < 1)
        parent_out->SetBinError(i, 1);
    }
    parent_out->SetLineColor(kGreen);
    parent_out->Write();
    f2->Close();


    return;
}
void remove_daughter_mult(){

    TFile *f0 = new TFile("output_mult_parent.root","READ");
    TFile *f1 = new TFile("output_mult_daughter.root","READ");

    TH1F* parent_bgs = ((TH1F*) f0->Get(Form("%s_mult_bgs",isotope.c_str()))); 
    TH1F* daughter_bgs = ((TH1F*) f1->Get(Form("%s_mult_bgs",isotope_daughter.c_str()))); 
    
    daughter_bgs->SetLineColor(kRed);
    TFile *f2 = new TFile(Form("%s_mult_final.root",isotope.c_str()),"RECREATE");
    parent_bgs->Write();
    daughter_bgs->Write();
//    add_bg_error(parent_bgs, daughter_bgs, scale_factor);
    daughter_bgs->Scale(scale_factor);
    parent_bgs->Add(daughter_bgs, -1);

    TH1F * parent_out = new TH1F(Form("%s_mult_final",isotope.c_str()),Form("%s_mult_final",isotope.c_str()),20,0,20);    
    parent_out->Add(parent_bgs);
    TH1F* parent_error = ((TH1F*) f0->Get(Form("%s_mult_bgs",isotope.c_str()))); 
    for (int i = 1; i < parent_out->GetNbinsX()+1; i++){
        parent_out->SetBinContent(i, parent_bgs->GetBinContent(i));
        parent_out->SetBinError(i, parent_bgs->GetBinError(i));
//        parent_out->SetBinError(i, (sqrt(pow(parent_error->GetBinError(i)/std::max(parent_error->GetBinContent(i),1.0),2)+pow(daughter_bgs->GetBinError(i)/std::max(daughter_bgs->GetBinContent(i),1.0),2)))*(std::max(parent_out->GetBinContent(i),1.0)));
        if (parent_out->GetBinError(i) < 1)
        parent_out->SetBinError(i, 1);
    }
    parent_out->SetLineColor(kGreen);
    parent_out->Write();
    f2->Close();


    return;
    
}
void remove_daughter_tas(){

    TFile *f0 = new TFile("output_tas_parent.root","READ");
    TFile *f1 = new TFile("output_tas_daughter.root","READ");

    TH1F* parent_bgs = ((TH1F*) f0->Get(Form("%s_tas_bgs",isotope.c_str()))); 
    TH1F* daughter_bgs = ((TH1F*) f1->Get(Form("%s_tas_bgs",isotope_daughter.c_str()))); 
    
    daughter_bgs->SetLineColor(kRed);
    TFile *f2 = new TFile(Form("%s_tas_final.root",isotope.c_str()),"RECREATE");
    parent_bgs->Write();
    daughter_bgs->Write();
//    add_bg_error(parent_bgs, daughter_bgs, scale_factor);
    daughter_bgs->Scale(scale_factor);
    parent_bgs->Add(daughter_bgs, -1);

    TH1F * parent_out = new TH1F(Form("%s_tas_final",isotope.c_str()),Form("%s_tas_final",isotope.c_str()),int((tas_upper_limit-tas_lower_limit)/tas_bin_width),tas_lower_limit,tas_upper_limit);    
    parent_out->Add(parent_bgs);
    TH1F* parent_error = ((TH1F*) f0->Get(Form("%s_tas_bgs",isotope.c_str()))); 
    for (int i = 1; i < parent_out->GetNbinsX()+1; i++){
        parent_out->SetBinContent(i, parent_bgs->GetBinContent(i));
        parent_out->SetBinError(i, parent_bgs->GetBinError(i));
//        std::cout << parent_error->GetBinContent(i) << ", " << daughter_bgs->GetBinContent(i) << std::endl;
//        parent_out->SetBinError(i, (sqrt(pow(parent_error->GetBinError(i)/std::max(parent_error->GetBinContent(i),parent_error->GetBinError(i)),2)+pow(daughter_bgs->GetBinError(i)/std::max(daughter_bgs->GetBinContent(i),daughter_bgs->GetBinError(i)),2)))*(std::max(parent_out->GetBinContent(i),1.0)));
        if (parent_out->GetBinError(i) < 1)
        parent_out->SetBinError(i, 1);
    }
    parent_out->SetLineColor(kGreen);
    parent_out->Write();
    f2->Close();


    return;
}

void remove_daughter() {

    remove_daughter_ss();
    remove_daughter_tas();
    remove_daughter_mult();
    remove_daughter_tas_ss();

    return;
}

int main() {
    remove_daughter();
    return 0;
}
