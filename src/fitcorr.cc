#include "TFile.h"
#include "TH1F.h"
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include "TString.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TMinuit.h"

#include "utils.h"


using namespace std;

double parentN0;
double daughterl = TMath::Log(2)/347.0;
double gdaughterl = TMath::Log(2)/21000.1;
double ndaughterl = TMath::Log(2)/216.0;
double ngdaughterl = TMath::Log(2)/32000.0;
double daughterhl = p.Get_Daughter_HL();
double gdaughterhl = p.Get_Granddaughter_HL();
double ndaughterhl = p.Get_Neutron_Daughter_HL();
double ngdaughterhl = p.Get_Neutron_Granddaughter_HL();
const int bin_width = p.Get_Decay_Bin_Width();
int dssd_window = p.Get_DSSD_Window();
const int bins = 2000/bin_width;
std::vector<Double_t> errorc(bins);
std::vector<Double_t> stuffs(bins);
std::vector<Double_t> targlist(10);
std::vector<Double_t> binno(bins);
std::vector<Double_t> xlim = {0, 2000};



//void save_to_png( TCanvas *c1, std::string name ){
//
//    TImage *img = TImage::Create();
//    img->FromPad(c1);
//    img->WriteImage(Form("%s.png",name.c_str()));
//    return;
//}
//
//


//bin = 1;
//static const bins;
//Functions

//Double_t bateman(double t,Double_t *targlist);

//void Cint::G__CallFunc::SetArgArray(long* p, int narg=-1);

Double_t bateman(double x, Double_t *par)
   {
   
      Double_t C1, dC1, dC2, ndC1, ndC2, N;
      Double_t dC1n, dndC1, dC2n, dndC2;
      Double_t gdC1, gdC2, gdC3;
   
      //Bateman equation for parent
      C1 = par[1];                                       /* Initial parent activity */
   
      //Bateman equations for daughter
      dC1 = daughterl*par[1]/((daughterl-par[0]));
      dC2 = daughterl*par[1]/((par[0]-daughterl));
   
      //Bateman equations for gdaughter
      gdC1 = gdaughterl*par[1]*daughterl/((daughterl-par[0])*(gdaughterl-par[0]));
      gdC2 = gdaughterl*par[1]*daughterl/((par[0]-daughterl)*(gdaughterl-daughterl));
      gdC3 = gdaughterl*par[1]*daughterl/((gdaughterl-par[0])*(gdaughterl-daughterl));
   
      //Bateman equation for neutron daughter
      ndC1 = ndaughterl*par[1]/((ndaughterl-par[0]));
      ndC2 = ndaughterl*par[1]/((par[0]-ndaughterl));
   
      N = C1*(TMath::Exp(-par[0]*x)) + (1-par[3])*(dC1*TMath::Exp(-par[0]*x) + dC2*TMath::Exp(-daughterl*x)) + (par[3])*(ndC1*TMath::Exp(-par[0]*x) + ndC2*TMath::Exp(-ndaughterl*x))+ /* + par[2]*TMath::Exp(-par[3]*x[0])*/ (1-par[3])*(gdC1*TMath::Exp(-par[0]*x) + gdC2*TMath::Exp(-daughterl*x) + gdC3*TMath::Exp(-gdaughterl*x))+ par[2];
//      N = C1*(TMath::Exp(-par[0]*x))/* + (1-par[3])*(dC1*TMath::Exp(-par[0]*x) + dC2*TMath::Exp(-daughterl*x)) + (par[3])*(ndC1*TMath::Exp(-par[0]*x) + ndC2*TMath::Exp(-ndaughterl*x))+  + par[2]*TMath::Exp(-par[3]*x[0]) (1-par[3])*(gdC1*TMath::Exp(-par[0]*x) + gdC2*TMath::Exp(-daughterl*x) + gdC3*TMath::Exp(-gdaughterl*x))*/+ par[2];
   
     /* par[0] = parent half-life
        par[1] = initial parent activity
        par[2] = background
        par[3] = Pn */
   
     Double_t value = N;
     return value;
   
   }
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)   /* This calculates Xi square */
{

   Int_t i;

   //calculate chisquare
   Double_t chisq = 0;
   Double_t delta;
   for (i=0;i<bins;i++) {
     binno[i] = i;
     if(errorc[i] >0) 
       {
         delta = (stuffs[i]-bateman(binno[i],par))/errorc[i] /*(TMath::Sqrt(func(i+(bin/2),par)))*/ ;    
         //cout << delta << endl;
          //cout << bateman(binno,par) << endl;
       }
       else delta = 0;
     chisq += delta*delta;
     }
   f = chisq;
  //f = bateman->GetChisquare();
}

Double_t parent_f(Double_t *x, Double_t *par){
   Double_t parentval = par[0]*TMath::Exp(-par[1]*x[0]);
   /*par[0] : Parent initial activity
     par[1] : Parent half-life */
   return parentval;
}



Double_t daughter_f(Double_t *x, Double_t *par){
   Double_t daughterval = par[0]*TMath::Exp(-par[2]*x[0]) + par[1]*TMath::Exp(-daughterl*x[0]); 
   /*par[0] : dC1
     par[1] : dC2
     par[2] : parent half-life*/
   return daughterval;  
}

Double_t gdaughter_f(Double_t *x, Double_t *par){
   Double_t gdaughterval = par[0]*TMath::Exp(-par[3]*x[0]) + par[1]*TMath::Exp(-daughterl*x[0]) + par[2]*TMath::Exp(-gdaughterl*x[0]); 
   /*par[0] : gdC1
     par[1] : gdC2
     par[2] : gdC3
     par[3] : parent half-life*/
   return gdaughterval;  
}
Double_t ndaughter_f(Double_t *x, Double_t *par){
   Double_t ndaughterval = par[0]*TMath::Exp(-par[2]*x[0]) + par[1]*TMath::Exp(-ndaughterl*x[0]); 
   /*par[0] : dC1
     par[1] : dC2
     par[2] : parent half-life*/
   return ndaughterval;  
}

Double_t background_f(Double_t *x, Double_t *par){
   Double_t bg = par[0];
   return bg;
}

   
//Define total bateman function
Double_t total_f(Double_t *x, Double_t *par)
{

   Double_t C1, dC1, dC2, ndC1, ndC2, N;
   Double_t dC1n, dndC1, dC2n, dndC2;
   Double_t gdC1, gdC2, gdC3;

   //Bateman equation for parent
   C1 = par[1];                                       /* Initial parent activity */

   //Bateman equations for daughter
   dC1 = daughterl*par[1]/((daughterl-par[0]));
   dC2 = daughterl*par[1]/((par[0]-daughterl));

   //Bateman equations for gdaughter
   gdC1 = gdaughterl*par[1]*daughterl/((daughterl-par[0])*(gdaughterl-par[0]));
   gdC2 = gdaughterl*par[1]*daughterl/((par[0]-daughterl)*(gdaughterl-daughterl));
   gdC3 = gdaughterl*par[1]*daughterl/((gdaughterl-par[0])*(gdaughterl-daughterl));
   
   //Bateman equation for neutron daughter
   ndC1 = ndaughterl*par[1]/((ndaughterl-par[0]));
   ndC2 = ndaughterl*par[1]/((par[0]-ndaughterl));

   N = C1*(TMath::Exp(-par[0]*x[0])) + (1-par[3])*(dC1*TMath::Exp(-par[0]*x[0]) + dC2*TMath::Exp(-daughterl*x[0])) + (par[3])*(ndC1*TMath::Exp(-par[0]*x[0]) + ndC2*TMath::Exp(-ndaughterl*x[0]))+(1-par[3])*(gdC1*TMath::Exp(-par[0]*x[0]) + gdC2*TMath::Exp(-daughterl*x[0]) + gdC3*TMath::Exp(-gdaughterl*x[0]))+ par[2];
//   N = C1*(TMath::Exp(-par[0]*x[0])) /*+ (1-par[3])*(dC1*TMath::Exp(-par[0]*x[0]) + dC2*TMath::Exp(-daughterl*x[0])) + (par[3])*(ndC1*TMath::Exp(-par[0]*x[0]) + ndC2*TMath::Exp(-ndaughterl*x[0]))+(1-par[3])*(gdC1*TMath::Exp(-par[0]*x[0]) + gdC2*TMath::Exp(-daughterl*x[0]) + gdC3*TMath::Exp(-gdaughterl*x[0]))*/+ par[2];

  /* par[0] = parent half-life
     par[1] = initial parent activity
     par[2] = background
     par[3] = Pn */

  Double_t value = N;
  return value;

}

//void fitcorr(const char* inputfile){
Double_t single_fitcorr(std::string isotope="sc54", int time_window = 1000, double pn_value = 1.0){
    TFile *fOut = new TFile("output_fitcorr.root","RECREATE");
   //Get inputfile
    std::string str_time_window = Form("%d",time_window);
    TChain* t = new TChain("t");
    TChain* t_bg = new TChain("t");
//    TChain* corr_flags = new TChain("correlation_flags");
    int total_implants = 0;
    int good_implants = 0;
    int good_decays = 0;
    int bad_decays_t = 0;
    int bad_decays_x = 0;
    std::string pixels = "5";
    std::string scint = "3650";
//    char pixels{"1","3","5","7","9","11","13","15"};
//    std::cout << pixels[0] << std::endl;
    TH1F *h1 = new TH1F("h1","h1",50,0,50);
//    for (int pixel=0; pixel<8; pixel++){
//    for (int i=2154; i<2169; i++)
    std::string dssd_window_str = Form("%d",dssd_window);
    for (int i=2054; i<2169; i++)
    {
//        if( i == 2087)
//            continue;
        std::string run_str = Form("%d",i);
        std::string filepath = "/mnt/analysis/e17028/correlated_rootfiles/"+str_time_window+"ms/true/Run" + run_str + "_SuN_" + str_time_window + "ms_" + isotope + "_" + dssd_window_str + "px_bound.root";
        std::string filepath_bg = "/mnt/analysis/e17028/correlated_rootfiles/"+str_time_window+"ms/bg/BWRun" + run_str + "_SuN_" + str_time_window + "ms_" + isotope + "_" + dssd_window_str + "px_bound.root";
//        std::string filepath = "/mnt/analysis/e17028/correlated_rootfiles/"+str_time_window+"ms/true/Run" + run_str + "_SuN_" + str_time_window + "ms_" + isotope + "_9px_bound.root";
//        std::string filepath = "./Run" + run_str + "_SuN_" + str_time_window + "ms_" + isotope + ".root";
        const char * c = filepath.c_str();
//        std::cout << filepath << std::endl;
        t->Add(c);
        TFile *fIn = TFile::Open(c);
        TH1F* corr_flags = ((TH1F*) gDirectory->Get("correlation_flags"));
        h1->Add(corr_flags);
        total_implants+=corr_flags->Integral(5,5);
        good_implants+=corr_flags->Integral(13,13);
//        good_implants+=corr_flags->Integral(11,11);
        good_decays+=corr_flags->Integral(29,29);
        bad_decays_t+=corr_flags->Integral(25,25);
        bad_decays_x+=corr_flags->Integral(21,21);
        const char * c_bg = filepath_bg.c_str();
//        std::cout << filepath_bg << std::endl;
        t_bg->Add(c_bg);
        TFile *fIn_bg = TFile::Open(c_bg);
    } 
//   std::string input_temp = filepath;
//   char* inputfile = new char[input_temp.size()+1];
//   std::copy(input_temp.begin(), input_temp.end(), inputfile);
//   inputfile[input_temp.size()] = '\0';
//   std::cout<<input_temp<<std::endl;
//   char* outputfile = new char[output_temp.size()+1];
//   std::copy(output_temp.begin(), output_temp.end(), outputfile);
//   outputfile[output_temp.size()] = '\0';
//   char * inputfile = input_temp.c_str();
//   TFile *fIn = TFile::Open(inputfile);

//   TChain* t = new TChain("t");
//   t->Add(inputfile);
   TCanvas *c0 = new TCanvas();
    std::string seg_cut = "(";
    for (int i_seg=0; i_seg<8; i_seg++){
//        std::cout << seg_cut.size() << std::endl;
        if (seg_cut.size() > 1 && i_seg < 8)
            seg_cut += "||";
        { 
            seg_cut += Form("(SuN_Seg_high_gain_cal[%d]>1525&&SuN_Seg_high_gain_cal[%d]<1650)",i_seg,i_seg);
        }
    }
    seg_cut += ")";

   
   //Get histogram
   //ft->Draw("correlation_implant_decay_time>>hcorr(1000,0,1000000000)","correlation_implant_decay_time > 0 && SuN_TAS_high_gain_cal > 388 && SuN_TAS_high_gain_cal < 465");
//   std::string filter = "SuN_TAS_high_gain_cal>0&&SuN_TAS_high_gain_cal<10001&&correlation_implant_decay_time>1250000&&correlation_implant_decay_time<2000000000";
//   t->Draw(Form("correlation_implant_decay_time*0.000001>>hcorr(%d,1,%d)",(int)(bins),(int)(xlim[1]*1)),Form("%s",filter.c_str()),"");
//   t->Draw(Form("correlation_implant_decay_time*0.000001>>hcorr(%d,%d,%d)",(int)(bins),(int)(xlim[0]),(int)(xlim[1]*1+1)),Form("correlation_implant_decay_time > 0 && correlation_implant_decay_time < %d && E_veto_raw < 25 && %s",(int)(xlim[1]*1e6), seg_cut.c_str()),"");
   t->Draw(Form("correlation_implant_decay_time*0.000001>>hcorr(%d,%d,%d)",(int)(bins),(int)(xlim[0]),(int)(xlim[1]*1+1)),Form("correlation_implant_decay_time > 0 && correlation_implant_decay_time < %d && E_veto_raw < 25",(int)(xlim[1]*1e6)),"");
//   t_bg->Draw(Form("correlation_implant_decay_time*0.000001>>hcorr_bg(%d,%d,%d)",(int)(bins),(int)(xlim[0]),(int)(xlim[1]*1+1)),Form("correlation_implant_decay_time > 0 && correlation_implant_decay_time < %d && E_veto_raw <25 && %s",(int)(xlim[1]*1e6), seg_cut.c_str()),"");
   t_bg->Draw(Form("correlation_implant_decay_time*0.000001>>hcorr_bg(%d,%d,%d)",(int)(bins),(int)(xlim[0]),(int)(xlim[1]*1+1)),Form("correlation_implant_decay_time > 0 && correlation_implant_decay_time < %d && E_veto_raw <25",(int)(xlim[1]*1e6)),"");
   TH1F* histo = ((TH1F*) gDirectory->Get("hcorr"));
   TH1F* histo_bg = ((TH1F*) gDirectory->Get("hcorr_bg"));
//   printf("%s",t_bg->Print());
//   histo->GetYaxis()->SetRangeUser(0,260);
//   histo_bg->GetYaxis()->SetRangeUser(0,260);
   //t->Draw("correlation_implant_decay_time>>hcorr(1000,0,1000000000)","td_NERO_raw>0");
//   TH1F *histo = (TH1F*)hcorr->Clone("histo");
   histo->Draw();
   histo_bg->SetLineColor(kRed);
   histo_bg->Draw("same");
   printf("%d\n", histo->GetBinContent(900));
   save_to_png(c0, "raw");
//   histo->Add(histo,histo_bg,1,-1);
   printf("%d\n", histo->GetBinContent(900));
   printf("%d\n", histo->FindLastBinAbove(0,1));
   //Fill the array from the histogram
   for(int i = 0; i < histo->FindLastBinAbove(0,1); i++){
     //cout << "i = " << i << " ";
     //cout << histo->GetBinContent(i);
     stuffs[i] = histo->GetBinContent(i);
//     cout << " Stuffs is " << stuffs[i];
     errorc[i] = histo->GetBinError(i);
//     cout << " and Errorc is " << errorc[i] << endl;
   }
   //Run Minuit
   TMinuit *gMinuit = new TMinuit(4); //TMinuit(number of params) 
   gMinuit->SetFCN(fcn);         //fcn is the function to minimise, in this case chisq   
   
   Double_t arglist[10];
   Int_t ierflg = 0;
   arglist[0] = 1;
   gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
   
   cout << "aftermnexcm" << endl;
  /* param[0]: lambda_parent in ms
     param[1]: parent initial activity 
     param[2]: Background initial activity
     param[3]: Pn of parent */
   
//   int Bin1 = histo->GetBinContent(1);
   double Bin1 = histo->GetMaximum();
//   int BinBG = histo->GetBinContent(xlim[1]-10);
   double BinBG = histo->Integral(bins-10,bins)/10.0;
//   BinBG = 1;
//   BinBG = TMath::Max(BinBG,(double)0);
   printf("Bin1: %d, BinBG: %d\n", Bin1, BinBG);
//   int BinBG = 30;
   Double_t vstart[4] = {TMath::Log(2)/100,Bin1,BinBG,0.10};   //start param
   Double_t step[4] = {.00001,1,0.01,.01};        // Step size used to minimize the parameters
  
  /* Set the fit variables to their initial states from above */

  /*gMinuit->mnparm(par no., parname, starting guess, step size, 
    limit low, limit upper,command (ierflg=0 if mnparm is successful) */
  
  //Parent half-life, param[0]
  gMinuit->mnparm(0, "parent half-life", vstart[0], step[0],TMath::Log(2)/(1000),TMath::Log(2)/(1.0),ierflg);
//  gMinuit->mnparm(0, "parent half-life", vstart[0], step[0],0,10000,ierflg);
   cout <<"after parent" << endl;
  //Initial parent activity, param[1]
   gMinuit->mnparm(1, "N0" , vstart[1], step[1], Bin1*0.0 ,Bin1*4.5,ierflg);
  
  //Background, param[2]
   gMinuit->mnparm(2, "Cbackground", vstart[2], step[2], BinBG*0.0, BinBG*1.9, ierflg); 
//   gMinuit->mnparm(2, "Cbackground", vstart[2], step[2], BinBG*0.0, BinBG*0.1, ierflg); 

  //Pn, param[3]
    if( pn_value < 0.4) {
    if(pn_value >= 0.02)
   gMinuit->mnparm(3, "Pn", vstart[3], step[3],pn_value-0.02,pn_value+0.02,ierflg); 
    else
   gMinuit->mnparm(3, "Pn", vstart[3], step[3],0.00,pn_value+0.02,ierflg); 
    }
    else {
   gMinuit->mnparm(3, "Pn", vstart[3], step[3],0,0.1,ierflg); 
    }
//   gMinuit->mnparm(3, "Pn", vstart[3], step[3],0,pn_value,ierflg); 
   
   cout << "Managed to Set Params" << endl; 
  
  // Now ready for minimization step
   arglist[0] = 10000;
   arglist[1] = 0.01;
   cout << "fin" << endl;

   //Run MIGRAD
   gMinuit->mnexcm("MIGRAD", arglist,2,ierflg); 
   
   Double_t amin, edm, errdef;
   Int_t nvpar, nparx,icstat;
   gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
   gMinuit->mnprin(3,amin);


   //Get the minimised parameters
   Double_t ParentT[2];
   gMinuit->GetParameter(0,ParentT[0],ParentT[1]);
   cout << "Set the parent" << endl;
   cout << ParentT[0] << ", " << ParentT[1] << endl;
 
   Double_t N0[2];
   gMinuit->GetParameter(1,N0[0],N0[1]);
   cout << "Set the init activity" << endl;
   cout << Bin1 << endl;
   
   Double_t Background[2];
   gMinuit->GetParameter(2,Background[0],Background[1]);
   cout << "Set the background" << endl;
   
   Double_t Pn[2];
   gMinuit->GetParameter(3,Pn[0],Pn[1]);
   cout << "Set the Pn" << endl;
  
   
   //Get the co-efficients
   double dC1 = (1-Pn[0])*daughterl*N0[0]/((daughterl-ParentT[0]));
   double dC2 = (1-Pn[0])*daughterl*N0[0]/((ParentT[0]-daughterl));

   double gdC1 = (1-Pn[0])*daughterl*N0[0]*gdaughterl/((daughterl-ParentT[0])*(gdaughterl-ParentT[0]));     
   double gdC2 = (1-Pn[0])*daughterl*N0[0]*gdaughterl/((ParentT[0]-daughterl)*(gdaughterl-daughterl));  
   double gdC3 = (1-Pn[0])*daughterl*N0[0]*gdaughterl/((gdaughterl-ParentT[0])*(gdaughterl-daughterl)); 

   double ndC1 = Pn[0]*ndaughterl*N0[0]/(ndaughterl-ParentT[0]);
   double ndC2 = Pn[0]*ndaughterl*N0[0]/(ParentT[0]-ndaughterl);


   TCanvas *c1 = new TCanvas();

   // Make the functions
//   TF1 *dparent = new TF1("dparent", parent,xlim[0],xlim[1],2);
   TF1 *dparent = new TF1("dparent",parent_f,xlim[0],bins,2);
   dparent->SetParameter(0,N0[0]);
   dparent->SetParameter(1,ParentT[0]);
   dparent->SetParNames("N0","ParentT");
   dparent->SetLineColor(kBlack);
   cout << dparent->GetParameter(0) << ", " << dparent->GetParameter(1) << endl;
   
//   TF1 *daughter = new TF1("daughter", daughter,xlim[0],xlim[1],3);
   TF1 *daughter = new TF1("daughter", daughter_f,xlim[0],bins,3);
   daughter->SetParameter(0,dC1);
   daughter->SetParameter(1,dC2);
   daughter->SetParameter(2,ParentT[0]);
   daughter->SetLineColor(kGreen);
   cout << daughter->GetParameter(0) << ", " << daughter->GetParameter(1) << ", " << daughter->GetParameter(2) << endl;
   
//   TF1 *daughter = new TF1("daughter", daughter,xlim[0],xlim[1],3);
   TF1 *gdaughter = new TF1("gdaughter", gdaughter_f,xlim[0],bins,4);
   gdaughter->SetParameter(0,gdC1);
   gdaughter->SetParameter(1,gdC2);
   gdaughter->SetParameter(2,gdC3);
   gdaughter->SetParameter(3,ParentT[0]);
   gdaughter->SetLineColor(kYellow);
   cout << gdaughter->GetParameter(0) << ", " << gdaughter->GetParameter(1) << ", " << gdaughter->GetParameter(2) << endl;
   
//   TF1 *ndaughter = new TF1("ndaughter", ndaughter,xlim[0],xlim[1],3);
   TF1 *ndaughter = new TF1("ndaughter", ndaughter_f,xlim[0],bins,3);
   ndaughter->SetParameter(0,ndC1);
   ndaughter->SetParameter(1,ndC2);
   ndaughter->SetParameter(2,ParentT[0]);
   ndaughter->SetLineColor(kCyan);
   cout << ndaughter->GetParameter(0) << ", " << ndaughter->GetParameter(1) << ", " << ndaughter->GetParameter(2) << endl;

//   TF1 *background = new TF1("background",background,xlim[0],xlim[1],1);
   TF1 *background = new TF1("background",background_f,xlim[0],bins,1);
   background->SetParameter(0,Background[0]);
   background->SetLineColor(kRed);
   cout << background->GetParameter(0) << endl;

//   TF1 *totalrate = new TF1("totalrate",total,xlim[0],xlim[1],4);
   TF1 *totalrate = new TF1("totalrate",total_f,xlim[0],bins,4);
   totalrate->SetParameter(0,ParentT[0]);
//   totalrate->SetParameter(0,TMath::Log(2)*2/51.1);
   totalrate->SetParameter(1,N0[0]);
   totalrate->SetParameter(2,Background[0]);
   totalrate->SetParameter(3,Pn[0]);
   totalrate->SetLineColor(kBlue);
   cout << totalrate->GetParameter(0) << ", " << totalrate->GetParameter(1) << ", " << totalrate->GetParameter(2) << ", " << totalrate->GetParameter(3) << endl;
   dparent->Draw();
   daughter->Draw("same");
   gdaughter->Draw("same");
   ndaughter->Draw("same");
   background->Draw("same");
   totalrate->Draw("same");
   
//   TH1F *decaycurve = new TH1F("decaycurve","decaycurve",bins,xlim[0],xlim[1]);
   TH1F *decaycurve = new TH1F("decaycurve","decaycurve",bins,xlim[0],bins);
   std::cout << dparent->Integral(0.0,xlim[1]) << ", " << daughter->Integral(0.0,xlim[1])<< ", " << ndaughter->Integral(0.0,xlim[1])<< std::endl;
   for (int i=0; i<bins; i++){
     int content = histo->GetBinContent(i);
     int error = histo->GetBinError(i);
//     printf("%d\n", content);
     decaycurve->SetBinContent(i,content);
     decaycurve->SetBinError(i,error);
//     if(i<bins-1){
//     std::cout << dparent->Integral(0.0,xlim[1])-dparent->Integral(0,i) << ", " << (dparent->Integral(0,i))*0.92-daughter->Integral(0.0,i) << ", " << (dparent->Integral(0,i))*0.08-ndaughter->Integral(0.0,i) << std::endl;
//     }
   }
   decaycurve->Draw("same");
   dparent->SetTitle(Form("presumed-%s, calcT_1/2 = %0.fms", isotope.c_str(), TMath::Log(2)/ParentT[0]*histo->GetBinWidth(2)));
//   dparent->SetTitleX(0.5);
//   dparent->SetTitleAlign(23);
   parentN0 = N0[0]/ParentT[0];//*1000*histo->GetBinWidth(0)/1e9;
//   c1->SetLogy();
   dparent->GetYaxis()->SetRangeUser(1., decaycurve->GetMaximum()+ 10.0);
   dparent->SetTitle("");
   dparent->GetYaxis()->SetTitle(Form("Decays/%dms",bin_width));
   dparent->GetXaxis()->SetTitle(Form("Decay Time (%dms bins)",bin_width));
   TLegend *leg = new TLegend(0.5,0.6,0.8,0.8);
   leg->AddEntry(totalrate,"Total");
   leg->AddEntry(dparent, p.Get_Isotope().c_str());// %0.1fms",TMath::Log(2)/(ParentT[0])*(histo->GetBinWidth(2))));
   leg->AddEntry(background, "Background");
   leg->AddEntry(daughter, p.Get_Isotope_Daughter().c_str());
   leg->AddEntry(gdaughter, p.Get_Isotope_Granddaughter().c_str());
   leg->AddEntry(ndaughter,p.Get_Isotope_Neutron_Daughter().c_str());
   decaycurve->SetStats(0);
   leg->SetBorderSize(0);
   leg->Draw("Same");
   save_to_png( c1, Form("%s-%s-dt",isotope.c_str(),pixels.c_str()) );
   cout << "done"<< " "<< parentN0 << " bin_width "<< histo->GetBinWidth(2)<<endl;
   cout << "half-life " << TMath::Log(2)/(ParentT[0])*(histo->GetBinWidth(2)) << endl;
   cout << totalrate->Integral(0.0,xlim[1])<< " " << decaycurve->Integral(0,xlim[1])<< endl;
   cout << daughter->Integral(0.0,xlim[1])<< " " << ndaughter->Integral(0,xlim[1])<< endl;
   cout << "decays: " << dparent->Integral(0.0,xlim[1]) << " background: " <<  background->Integral(0,xlim[1])<< endl;
   std::cout<<total_implants<<" "<<good_implants<< " " << dparent->Integral(0.0,xlim[1])<<std::endl;
   std::cout<<good_decays<< " " << bad_decays_t << " " << bad_decays_x <<std::endl;
//   int parent_time = 300/bin_width;
//   int daughter_time_lo = 600/bin_width;
//   int daughter_time_hi = 1200/bin_width;
   int parent_time_lo = p.Get_Parent_Decay_Window()[0]/bin_width;
   int parent_time_hi = p.Get_Parent_Decay_Window()[1]/bin_width;
   int daughter_time_lo = p.Get_Daughter_Decay_Window()[0]/bin_width;
   int daughter_time_hi = p.Get_Daughter_Decay_Window()[1]/bin_width;
   cout << dparent->Integral(parent_time_lo,parent_time_hi) << ", " << daughter->Integral(parent_time_lo,parent_time_hi) << ", " << ndaughter->Integral(parent_time_lo,parent_time_hi)  <<  ", " << gdaughter->Integral(parent_time_lo,parent_time_hi)  <<  ", " << dparent->Integral(daughter_time_lo,daughter_time_hi) << ", " << daughter->Integral(daughter_time_lo,daughter_time_hi)  << ", " << ndaughter->Integral(daughter_time_lo,daughter_time_hi) << ", " << gdaughter->Integral(daughter_time_lo, daughter_time_hi)  <<  endl;
   cout << dparent->Integral(parent_time_lo,parent_time_hi)/dparent->Integral(daughter_time_lo,daughter_time_hi) << ", " << daughter->Integral(parent_time_lo,parent_time_hi)/daughter->Integral(daughter_time_lo,daughter_time_hi) <<  ", " << gdaughter->Integral(parent_time_lo,parent_time_hi)/gdaughter->Integral(daughter_time_lo, daughter_time_hi)  <<  endl;
   cout << background->Integral(parent_time_lo,parent_time_hi) << ", " << background->Integral(daughter_time_lo,daughter_time_hi) << endl;
   cout << totalrate->Integral(parent_time_lo,parent_time_hi) << ", " << totalrate->Integral(daughter_time_lo,daughter_time_hi) << endl;
//   }
   fOut->cd();
   decaycurve->Write();
   fOut->Close();
   return TMath::Log(2)/ParentT[0]*histo->GetBinWidth(2);
}


Double_t fitcorr_bgs(double pn_value = p.Get_Pn() ){
//    daughters = std::make_tuple ("ti56", 190.0, 120.0); 
//    std::cout<<std::get<0>(daughters)<< " " << std::get<1>(daughters) <<std::endl;   
    std::string isotope = p.Get_Isotope();
    int time_window = 2000;
//    std::vector<double> daughters = halflives(isotope) ;
    daughterl = TMath::Log(2) / (daughterhl / xlim[1]*bins);
    gdaughterl = TMath::Log(2) / (gdaughterhl / xlim[1]*bins);
    ndaughterl = TMath::Log(2) / (ndaughterhl / xlim[1]*bins);
//    ngdaughterl = TMath::Log(2) / (ngdaughterhl / xlim[1]*bins);
    Double_t calc_hl = single_fitcorr( isotope, time_window, pn_value );
    return calc_hl;
}

int main(){
    fitcorr_bgs();
    return 0;
}
