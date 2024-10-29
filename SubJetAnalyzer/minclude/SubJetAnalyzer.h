#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

#include <string>

#include <iostream>
#include <iomanip>

#include <vector>
#include "math.h"
#include <numeric>
#include <string>
#include <fstream>
#include <algorithm>
#include <iterator>

#include <TStyle.h>
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TCutG.h"
#include "TRandom3.h"
#include "TGraphErrors.h" 
#include "coordinateTools.h" 
#include "1d2d_constants.h" 



#ifndef SubJetAnalyzer_h
#define SubJetAnalyzer_h

using TMath::ATan;
using TMath::Exp;
using namespace std;

class SubJetAnalyzer{
public:
   TFile          *fFile;
   TTree          *fChain;
   //long int        TotalEntry;
   std::vector<std::string> SDfileList;
   std::string      store_dir;
   
   //SDtrackTree
   std::vector< float > *jetPt;
   std::vector< float > *jetEta;
   std::vector< float > *jetPhi;
   std::vector< int   > *jetChargedMultiplicity;
   std::vector< std::vector<int>  > *has_parents;
   std::vector< std::vector<float> > *jetZg;
   std::vector< std::vector<float> > *jetRg;
   std::vector< std::vector<float> > *subjet1_pt;
   std::vector< std::vector<float> > *subjet1_eta;
   std::vector< std::vector<float> > *subjet1_phi;
   std::vector< std::vector<float> > *subjet1_m;
   std::vector< std::vector<float> > *subjet2_pt;
   std::vector< std::vector<float> > *subjet2_eta;
   std::vector< std::vector<float> > *subjet2_phi;
   std::vector< std::vector<float> > *subjet2_m;
   std::vector< std::vector<float> > *dau_pt;
   std::vector< std::vector<int> >   *dau_chg;
   std::vector< std::vector<float> > *dau_eta;
   std::vector< std::vector<float> > *dau_phi;
   std::vector< std::vector<int> >   *dau_pid;
   std::vector< std::vector< std::vector<int> > * > dau_Flag;   // if flag=0 , groomed
                                                                // if flag=1 , this particle belongs to subjet1
                                                                // if flga=2 , this particle belongs to subjet2   

   TBranch        *b_jetPt;
   TBranch        *b_jetEta;
   TBranch        *b_jetPhi;
   TBranch        *b_jetChargedMultiplicity;
   TBranch        *b_has_parents;
   TBranch        *b_jetZg;
   TBranch        *b_jetRg;
   TBranch        *b_subjet1_pt;
   TBranch        *b_subjet1_eta;
   TBranch        *b_subjet1_phi;
   TBranch        *b_subjet1_m;
   TBranch        *b_subjet2_pt;
   TBranch        *b_subjet2_eta;
   TBranch        *b_subjet2_phi;
   TBranch        *b_subjet2_m;
   TBranch        *b_dau_pt;
   TBranch        *b_dau_chg;
   TBranch        *b_dau_eta;
   TBranch        *b_dau_phi;
   TBranch        *b_dau_pid;


   std::vector<TBranch *> b_dau_Flag ;

   
   std::vector<double> Classifier_lo;
   std::vector<double> Classifier_hi;

   int mode;
   int beta_index;
   double R;


   SubJetAnalyzer(std::vector<std::string> _SDfileList,std::string _store_dir);//,int _mode,int _beta_index,double _R);
   virtual ~SubJetAnalyzer();
   virtual void     gen_background(TH2D *hraw,TH2D *hbkg,long int NENT);
   virtual void     LoadData(TTree *tree);
   virtual void     CalculateCorrelation();// calculate 2D correlation function
   virtual void     DrawFlow();// calcualte 1D correlation function 
   virtual void     Draw2DCorrelation();// calculate Fourier cofficient  
   virtual double   Cal_SubStructure_classifier(int ijet);



};


SubJetAnalyzer::SubJetAnalyzer(std::vector<std::string> _SDfileList,std::string _store_dir )//,int _mode=0,int _beta_index=0,double _R=0.8)
   :fChain(0),SDfileList(_SDfileList),store_dir(_store_dir),mode(0),beta_index(2),R(0.8)
{
   
}







SubJetAnalyzer::~SubJetAnalyzer()
{

}


void SubJetAnalyzer::gen_background(TH2D *hraw,TH2D *hbkg,long int NENT){
    
        int backMult =10;
        long int XENT = ((1+floor(sqrt(1+(4*2*backMult*NENT))))/2) ;////////////////////
        vector<double> A_ETA_Cor(XENT,0.);
        vector<double> A_PHI_Cor(XENT,0.);
        for(int x = 0; x<XENT; x++){
                gRandom->SetSeed(0);
                double WEta1_Cor, WPhi1_Cor;//making the pseudoparticles
                hraw->GetRandom2(WEta1_Cor, WPhi1_Cor);
                A_ETA_Cor[x] = WEta1_Cor;
                A_PHI_Cor[x] = WPhi1_Cor;
        }
        for(long int i = 0; i < (XENT-1); i++){
            for(long int j = (i+1); j < XENT; j++){

                double WdeltaEta_Cor = (A_ETA_Cor[i]-A_ETA_Cor[j]);
                double WdeltaPhi_Cor = (TMath::ACos(TMath::Cos(A_PHI_Cor[i]-A_PHI_Cor[j])));
                hbkg->Fill(WdeltaEta_Cor, WdeltaPhi_Cor,   1);//./XENT);
                hbkg->Fill(-WdeltaEta_Cor, WdeltaPhi_Cor,  1);//./XENT);
                hbkg->Fill(WdeltaEta_Cor, -WdeltaPhi_Cor,  1);//./XENT);
                hbkg->Fill(-WdeltaEta_Cor, -WdeltaPhi_Cor, 1);//./XENT);
                hbkg->Fill(WdeltaEta_Cor, 2*TMath::Pi() - WdeltaPhi_Cor, 1);//./XENT);
                hbkg->Fill(-WdeltaEta_Cor,2*TMath::Pi() - WdeltaPhi_Cor, 1);//./XENT);

            }
        }
}



void SubJetAnalyzer::LoadData(TTree *tree){
   
   b_dau_Flag.clear();
   
   for(int ibeta=0;ibeta<beta_bin;ibeta++){
      b_dau_Flag.push_back(0);
   }
   

   jetPt=0;
   jetEta=0;
   jetPhi=0;
   jetChargedMultiplicity=0;
   has_parents=0;
   jetZg=0;
   jetRg=0;
   subjet1_pt=0;
   subjet1_eta=0;
   subjet1_phi=0;
   subjet1_m=0;
   subjet2_pt=0;
   subjet2_eta=0;
   subjet2_phi=0;
   subjet2_m=0;
   dau_pt=0;
   dau_chg=0;
   dau_eta=0;
   dau_phi=0;
   dau_pid=0;
   dau_Flag.clear();
   for(int i=0;i<beta_bin;i++){
      dau_Flag.push_back(0);
   }
   
   if (!tree) {
         std::cout<<"Input tree does not exist "<<std::endl;
         return;
   }
   
   fChain = tree;
   //TotalEntry=fChain->GetEntriesFast();
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("jetPt"                 ,&jetPt                 ,&b_jetPt                 );
   fChain->SetBranchAddress("jetEta"                ,&jetEta                ,&b_jetEta                );
   fChain->SetBranchAddress("jetPhi"                ,&jetPhi                ,&b_jetPhi                );
   fChain->SetBranchAddress("jetChargedMultiplicity",&jetChargedMultiplicity,&b_jetChargedMultiplicity);
   fChain->SetBranchAddress("has_parents"           ,&has_parents           ,&b_has_parents           );
   fChain->SetBranchAddress("jetZg"                 ,&jetZg                 ,&b_jetZg                 );
   fChain->SetBranchAddress("jetRg"                 ,&jetRg                 ,&b_jetRg                 );
   fChain->SetBranchAddress("subjet1_pt"            ,&subjet1_pt            ,&b_subjet1_pt            );
   fChain->SetBranchAddress("subjet1_eta"           ,&subjet1_eta           ,&b_subjet1_eta           );
   fChain->SetBranchAddress("subjet1_phi"           ,&subjet1_phi           ,&b_subjet1_phi           );
   fChain->SetBranchAddress("subjet1_m"             ,&subjet1_m             ,&b_subjet1_m             );
   fChain->SetBranchAddress("subjet2_pt"            ,&subjet2_pt            ,&b_subjet2_pt            );
   fChain->SetBranchAddress("subjet2_eta"           ,&subjet2_eta           ,&b_subjet2_eta           );
   fChain->SetBranchAddress("subjet2_phi"           ,&subjet2_phi           ,&b_subjet2_phi           );
   fChain->SetBranchAddress("subjet2_m"             ,&subjet2_m             ,&b_subjet2_m             );
   fChain->SetBranchAddress("dau_pt"                ,&dau_pt                ,&b_dau_pt                );
   fChain->SetBranchAddress("dau_chg"               ,&dau_chg               ,&b_dau_chg               );
   fChain->SetBranchAddress("dau_eta"               ,&dau_eta               ,&b_dau_eta               );
   fChain->SetBranchAddress("dau_phi"               ,&dau_phi               ,&b_dau_phi               );
   fChain->SetBranchAddress("dau_pid"               ,&dau_pid               ,&b_dau_pid               );

   
   for(int i=0;i<beta_bin;i++){
      fChain->SetBranchAddress( Form("dau_flag_beta_%d",(int)beta_SD[i] )         ,& (dau_Flag[i] )         ,&(b_dau_Flag[i]) );
   }

   


}



double SubJetAnalyzer::Cal_SubStructure_classifier(int ijet){
   if( mode==0 ) return (*jetZg)[ijet][beta_index] ;
   if( mode==1 ) return (*jetRg)[ijet][beta_index] ;
   if( mode==2 ) return (*jetZg)[ijet][beta_index] * ( TMath::Power( (*jetRg)[ijet][beta_index] / (R) , beta_SD[beta_index]     )  ) ;
   else return 0;
}





void SubJetAnalyzer::CalculateCorrelation(){
   //Define Classifier
   int Classifier=1;
   if(mode ==0){//cut Zg
      Classifier=sdZg_bin;
      Classifier_lo.clear();
      Classifier_hi.clear();
      for(int itmp=0;itmp<sdZg_bin;itmp++){ Classifier_lo.push_back( sdZg_lo[itmp]  ); Classifier_hi.push_back( sdZg_hi[itmp]  );}
   }
   if(mode ==1){//cut Rg  
      Classifier=sdRg_bin;
      Classifier_lo.clear();
      Classifier_hi.clear();
      for(int itmp=0;itmp<sdRg_bin;itmp++){ Classifier_lo.push_back( sdRg_lo[itmp]  ); Classifier_hi.push_back( sdRg_hi[itmp]  );}
   }
   if(mode ==2){//cut Zg * (Thetag)^{beta}  
      Classifier=sdZgTg_bin;
      Classifier_lo.clear();
      Classifier_hi.clear();
      for(int itmp=0;itmp<sdZgTg_bin;itmp++){ Classifier_lo.push_back( sdZgTg_lo[itmp]  ); Classifier_hi.push_back( sdZgTg_hi[itmp]  );}
   }   
   if(mode>2 || mode<0) { std::cout<<"wrong mode!"<<std::endl; return ;}
   //End Define Classifier
   const int Classifier_bin=Classifier;
   TH1::SetDefaultSumw2(kTRUE);
   TH2::SetDefaultSumw2(kTRUE);
   double track_eta_lim = 5.0;
   //Initializing Histograms
   TH2D* hJet_Pass     = new TH2D("hJet_Pass"  ,"hJet_Pass"  , trackbin,bin0,trackbin,Classifier_bin,bin0,Classifier_bin);

   TH2D* hEPDrawCor[trackbin][Classifier_bin][ptbin];
   TH2D* hSignalShiftedCor[trackbin][Classifier_bin][ptbin];
   TH2D* hBckrndShiftedCor[trackbin][Classifier_bin][ptbin];
  


   


   TH2D* hSignal_groomed_groomed[trackbin][Classifier_bin][ptbin];      // use groomed particles to construct particle 
   TH2D* hSignal_nongroomed_nongroomed[trackbin][Classifier_bin][ptbin];// use non-groomed particles to construct particle

   TH2D* hBkg_groomed_groomed[trackbin][Classifier_bin][ptbin];         
   TH2D* hBkg_nongroomed_nongroomed[trackbin][Classifier_bin][ptbin];

   TH2D* hEPDrawCor_groomed[trackbin][Classifier_bin][ptbin];
   TH2D* hEPDrawCor_nongroomed[trackbin][Classifier_bin][ptbin];
   
   TH1D *h_Nch_dis[trackbin][Classifier_bin];
   for(int i=0;i<trackbin;i++){
      for(int j=0;j<Classifier_bin;j++){
         h_Nch_dis[i][j]=new TH1D(Form("h_Nch_dis_Nch_%d_%d_Class_%d_%d",trackbinbounds[i],trackbinboundsUpper[i],(int)(100*Classifier_lo[j]),(int)(100*Classifier_hi[j])),Form("h_Nch_dis_Nch_%d_%d_Class_%d_%d",trackbinbounds[i],trackbinboundsUpper[i],(int)(100*Classifier_lo[j]),(int)(100*Classifier_hi[j]) ),120,0,120);
      }
   }
                                                                                

   long int Pairs[trackbin][Classifier_bin][ptbin]             = {0};
   long int g_gPairs[trackbin][Classifier_bin][ptbin]          = {0};
   long int ng_ngPairs[trackbin][Classifier_bin][ptbin]        = {0};

   for(int i = 0; i < trackbin; i++){
      for(int j = 0; j < Classifier_bin; j++){
         for(int k=0; k < ptbin ;k++){
                 Pairs[i][j][k]=0;
                 g_gPairs[i][j][k]=0;
                 ng_ngPairs[i][j][k]=0;
         }
      }
   }
   

   for(int wtrk = 1; wtrk < trackbin +1 ;wtrk++){
        for(int wcla = 1; wcla<Classifier_bin+1; wcla++){
            for(int wppt = 1; wppt<ptbin+1; wppt++){
               hBckrndShiftedCor[wtrk-1][wcla-1][wppt-1]                         = new TH2D(Form("hBckrndS_Cor_trk_%d_ppt_%d_Class_%d",wtrk,wppt,wcla) ,Form("hBckrndS_trk_%d_ppt_%d_Class_%d",wtrk,wppt,wcla) ,41,-(20*EtaBW)-(0.5*EtaBW),(20*EtaBW)+(0.5*EtaBW),33,-(8*PhiBW)-0.5*PhiBW,(24*PhiBW)+0.5*PhiBW);
               hSignalShiftedCor[wtrk-1][wcla-1][wppt-1]                         = new TH2D(Form("hSignalS_Cor_trk_%d_ppt_%d_Class_%d",wtrk,wppt,wcla) ,Form("hSignalS_trk_%d_ppt_%d_Class_%d",wtrk,wppt,wcla) ,41,-(20*EtaBW)-(0.5*EtaBW),(20*EtaBW)+(0.5*EtaBW),33,-(8*PhiBW)-0.5*PhiBW,(24*PhiBW)+0.5*PhiBW);
               hEPDrawCor[wtrk-1][wcla-1][wppt-1]                                = new TH2D(Form("hEPDraw_Cor_trk_%d_ppt_%d_Class_%d",wtrk,wppt,wcla) ,Form( "hEPDraw_trk_%d_ppt_%d_Class_%d",wtrk,wppt,wcla) , EPD_xb   , EPD_xlo, EPD_xhi , EPD_yb      , EPD_ylo    , EPD_yhi);
               
               hSignal_groomed_groomed[wtrk-1][wcla-1][wppt-1]                   =new TH2D(Form("hSignalS_G_G_Cor_trk_%d_ppt_%d_Class_%d",wtrk,wppt,wcla) ,Form("hSignalS_G_G_trk_%d_ppt_%d_Class_%d",wtrk,wppt,wcla) ,41,-(20*EtaBW)-(0.5*EtaBW),(20*EtaBW)+(0.5*EtaBW),33,-(8*PhiBW)-0.5*PhiBW,(24*PhiBW)+0.5*PhiBW);
               hSignal_nongroomed_nongroomed[wtrk-1][wcla-1][wppt-1]             =new TH2D(Form("hSignalS_NG_NG_Cor_trk_%d_ppt_%d_Class_%d",wtrk,wppt,wcla) ,Form("hSignalS_NG_NG_trk_%d_ppt_%d_Class_%d",wtrk,wppt,wcla) ,41,-(20*EtaBW)-(0.5*EtaBW),(20*EtaBW)+(0.5*EtaBW),33,-(8*PhiBW)-0.5*PhiBW,(24*PhiBW)+0.5*PhiBW);
               hBkg_groomed_groomed[wtrk-1][wcla-1][wppt-1]                      =new TH2D(Form("hBckrndS_G_G_Cor_trk_%d_ppt_%d_Class_%d",wtrk,wppt,wcla) ,Form("hBckrndS_G_G_trk_%d_ppt_%d_Class_%d",wtrk,wppt,wcla) ,41,-(20*EtaBW)-(0.5*EtaBW),(20*EtaBW)+(0.5*EtaBW),33,-(8*PhiBW)-0.5*PhiBW,(24*PhiBW)+0.5*PhiBW);
               hBkg_nongroomed_nongroomed[wtrk-1][wcla-1][wppt-1]                =new TH2D(Form("hBckrndS_NG_NG_Cor_trk_%d_ppt_%d_Class_%d",wtrk,wppt,wcla) ,Form("hBckrndS_NG_NG_trk_%d_ppt_%d_Class_%d",wtrk,wppt,wcla) ,41,-(20*EtaBW)-(0.5*EtaBW),(20*EtaBW)+(0.5*EtaBW),33,-(8*PhiBW)-0.5*PhiBW,(24*PhiBW)+0.5*PhiBW);
               hEPDrawCor_groomed[wtrk-1][wcla-1][wppt-1]                        =new TH2D(Form("hEPDraw_G_Cor_trk_%d_ppt_%d_Class_%d",wtrk,wppt,wcla) ,Form( "hEPDraw_G_trk_%d_ppt_%d_Class_%d",wtrk,wppt,wcla) , EPD_xb   , EPD_xlo, EPD_xhi , EPD_yb      , EPD_ylo    , EPD_yhi);
               hEPDrawCor_nongroomed[wtrk-1][wcla-1][wppt-1]                     =new TH2D(Form("hEPDraw_NG_Cor_trk_%d_ppt_%d_Class_%d",wtrk,wppt,wcla) ,Form( "hEPDraw_NG_trk_%d_ppt_%d_Class_%d",wtrk,wppt,wcla) , EPD_xb   , EPD_xlo, EPD_xhi , EPD_yb      , EPD_ylo    , EPD_yhi);
            }
        }

   }
    //END Initializing Histograms 


   //print mode information 
   if(mode==0){std::cout<<"Cut on Zg"<<std::endl;}
   if(mode==1){std::cout<<"Cut on Rg"<<std::endl;}
   if(mode==2){std::cout<<"Cut on Zg*(Thetag)^(beta)"<<std::endl;}

   for(int itmp=0;itmp<Classifier_bin;itmp++){
      std::cout<<"CLASS "<<itmp<<" : "<<"["<<Classifier_lo[itmp]<<" , "<<Classifier_hi[itmp]<<"]";
   }

   

   std::cout << "Start generating 2D correlation histograms" << std::endl;
   //*******************ENTER FILE LOOP*******************************
   for(int ifile=0;ifile<SDfileList.size();ifile++){
      std::cout<<std::endl;
      std::cout << "Generating 2D correlation histograms" << std::endl;
      fFile=0;
      fFile = TFile::Open(SDfileList.at(ifile).c_str(),"READ");
      if(!fFile){
         std::cout<<SDfileList.at(ifile)<<" can not open , try next file"<<std::endl;
         continue;
      }
      TTree *tree = (TTree*)fFile->Get("SDtrackTree");
      if(!tree){
         std::cout<<SDfileList.at(ifile)<<" can not get tree,try next file"<<std::endl;
         continue;
      }
      LoadData(tree);

      std::cout << "File " << ifile+1 << " out of " << SDfileList.size() << std::endl;
      Long64_t nentries = fChain->GetEntriesFast();
      std::cout<<"Total Entries is:"<<nentries<<std::endl;
      std::cout<<std::endl;





    
      //*************ENTER EVENT LOOP*******************
      for (int ievent = 0; ievent<nentries ; ievent ++){
         
         if( ievent%10000==0 ) {std::cout<<"ievent="<<ievent<<std::endl; } 
         fChain->GetEntry(ievent);
         int jetCounter =(*jetPt).size();// the number of jet in ievent
         if( jetCounter == 0 ) continue;// There is no jet in this event 
         //*************ENTER JET LOOP*******************
         for(int ijet=0 ; ijet < jetCounter ; ijet++){
            long int NNtrk = ( (*dau_pt)[ijet] ).size();// number of daughter particles in ijet
             
            //if( (*has_parents)[ijet][beta_index]==-1  ) {continue; N_noparents++;  }//count the number of jets that do not have subjets, with given beta
            if( fabs(  (*jetEta) [ijet] ) > jetEtaCut      ) continue;//****************!!!!!****************
            if(        (*jetPt ) [ijet]   < jetPtCut_Jet   ) continue;//********************!!!!!******************************
            double  jet_HLT_weight                           = 1.0;


            
            std::vector<int> ClassifierBool(Classifier_bin,0);
            std::vector<int> tkBool(trackbin,0);
            double  n_ChargeMult_DCA_labPt_Eta_exclusion =0;
            
            //loop all daughter particle in ijet
            for(int A_trk=0; A_trk < NNtrk ; A_trk++ ){
               if(      (*dau_chg)[ijet][A_trk]    == 0   )              continue ; // not charged
               if(fabs( (*dau_pt )[ijet][A_trk] )  <  0.3 )              continue;//lab pt
               if(fabs( (*dau_eta)[ijet][A_trk] )  >  2.4 )              continue;//lab eta
               n_ChargeMult_DCA_labPt_Eta_exclusion += 1;
         
            }//end loop all daughter particle in ijet


           
            double SubStructure_classifier=Cal_SubStructure_classifier(ijet);// this relies on mode 
                


                

            for(int i = 0; i < trackbin; i++){
               for(int j =0; j < Classifier_bin ; j++ ){
                  if( n_ChargeMult_DCA_labPt_Eta_exclusion >= trackbinbounds[i] && n_ChargeMult_DCA_labPt_Eta_exclusion < trackbinboundsUpper[i] && SubStructure_classifier > Classifier_lo[j] && SubStructure_classifier < Classifier_hi[j] ){
                     tkBool[i] = 1;
                     ClassifierBool[j]=1;
                     hJet_Pass->Fill(i,j);
                     h_Nch_dis[i][j]->Fill(n_ChargeMult_DCA_labPt_Eta_exclusion);
                  }
               }
            }
               
            




            //int Ntrig[trackbin][Classifier_bin][ptbin]             = {0};
            double  NtrigCorrected[trackbin][Classifier_bin][ptbin] = {0};            
            double  Ntrig_g_Corrected[trackbin][Classifier_bin][ptbin] = {0};
            double  Ntrig_ng_Corrected[trackbin][Classifier_bin][ptbin] = {0};


            
            for(int i = 0; i < trackbin; i++){
                  for(int j = 0; j < Classifier_bin; j++){
                     for(int k=0; k < ptbin ;k++){
                        NtrigCorrected[i][j][k] = 0;
                        Ntrig_g_Corrected[i][j][k]=0; 
                        Ntrig_ng_Corrected[i][j][k] =0;

                  }
               }
            }
           
                
            vector< vector<int> > A_ptBool;
            A_ptBool.resize(NNtrk);
            for(int temp=0;temp<NNtrk;temp++){
               A_ptBool[temp].resize(ptbin);
            }


            //start loop daughter particles :A_trk
            for(int A_trk = 0 ; A_trk < NNtrk; A_trk++ ){

               if((* (dau_chg) )[ijet][A_trk] == 0)     continue;
               if((* (dau_pt)  )[ijet][A_trk] < 0.3)      continue;
               if(fabs((* (dau_eta) )[ijet][A_trk]) > 2.4) continue;


               double jet_dau_pt = ptWRTJet((double)(*jetPt)[ijet], (double)(* jetEta )[ijet], (double)(* jetPhi )[ijet], (double)(* dau_pt )[ijet][A_trk], (double)(* dau_eta )[ijet][A_trk], (double)(* dau_phi )[ijet][A_trk]);

               for(int i = 0; i < ptbin; i++){
                  if(jet_dau_pt > ptbinbounds_lo[i] && jet_dau_pt < ptbinbounds_hi[i]){
                     A_ptBool[A_trk][i] = 1; // defined within the loop of jets, find the pt range of  A_trk
                  }
               }

                  
               double Atrk_weight=1.0; 
               for(int i = 0; i < trackbin; i++){
                  for(int j = 0; j < Classifier_bin; j++){
                     for(int k=0; k < ptbin ;k++){

                           if(tkBool[i] +ClassifierBool[j] +A_ptBool[A_trk][k] == 3){
                              
                              NtrigCorrected[i][j][k] += (1.0/Atrk_weight);
                              //Ntrig[i][j][k] += 1; // defined within the loop of jets,count the number of the passed daughter particle in ijet
                              
                              if(   (*dau_Flag[beta_index])[ijet][A_trk] == 0 ) {Ntrig_g_Corrected[i][j][k]+= (1.0/Atrk_weight);}
                              if(   (*dau_Flag[beta_index])[ijet][A_trk]      ) {Ntrig_ng_Corrected[i][j][k]+= (1.0/Atrk_weight);}

                           //if(    (  Ntrig_g_Corrected[i][j][k]+Ntrig_ng_Corrected[i][j][k] ) !=  NtrigCorrected[i][j][k]      ){std::cout<<Ntrig_g_Corrected[i][j][k]<<" "<< Ntrig_ng_Corrected[i][j][k] << " " << NtrigCorrected[i][j][k]<<std::endl;}
                              
                           }
                    
                        }//end loop i
                  }//end loop j
               }//end loop k

               
              


            }//end loop daughter particles:A_trk
             
            //only for double check 
            for(int i = 0; i < trackbin; i++){
               for(int j = 0; j < Classifier_bin; j++){
                  for(int k=0; k < ptbin ;k++){
                     
                     if(    (  (int)Ntrig_g_Corrected[i][j][k]+(int)Ntrig_ng_Corrected[i][j][k] ) != (int) NtrigCorrected[i][j][k]      ){
                        
                        std::cout<<Ntrig_g_Corrected[i][j][k]<<" "<< Ntrig_ng_Corrected[i][j][k] << " " << NtrigCorrected[i][j][k]<<std::endl;
                     } 

                  }
               }
            }



            //**********START A_trk LOOP**************
            for(int A_trk = 0;A_trk < NNtrk ; A_trk++){
               if( (* dau_chg )[ijet][A_trk] == 0  ) continue;
               if( (* dau_pt  )[ijet][A_trk] < 0.3 ) continue;
               if(fabs( (*dau_eta)[ijet][A_trk] ) > 2.4) continue;


                    
               //calculate dauhter particles' pt , eta , phi in jet frame 
               double jet_dau_pt    =  ptWRTJet((double)(*jetPt)[ijet], (double)(*jetEta)[ijet]     , (double)(*jetPhi)[ijet]       , (double)(* dau_pt  )[ijet][A_trk], (double)(*dau_eta )[ijet][A_trk], (double)(* dau_phi )[ijet][A_trk]);
               double jet_dau_eta   = etaWRTJet((double)(*jetPt)[ijet], (double)(*jetEta)[ijet]    , (double)(*jetPhi)[ijet]       , (double)(* dau_pt  )[ijet][A_trk], (double)(*dau_eta )[ijet][A_trk],  (double)(* dau_phi )[ijet][A_trk]);
               double jet_dau_phi   = phiWRTJet((double)(*jetPt)[ijet], (double)(*jetEta)[ijet]     , (double)(*jetPhi)[ijet]       , (double)(* dau_pt   )[ijet][A_trk], (double)(*dau_eta)[ijet][A_trk], (double)(* dau_phi )[ijet][A_trk]);


                    



               double jet_dau_theta = 2*ATan(Exp(-(jet_dau_eta)));
               if(jet_dau_eta > track_eta_lim) continue;
               //getting et pt efficiency for A track in beam frame
               double Atrk_weight = 1.0;//(hReco2D[thisEffTable]->GetBinContent(hReco2D[thisEffTable]->FindBin( (*dau_pt)[ijet][A_trk] , (*dau_eta)[ijet][A_trk] )));
                    


                
                
               for(int i = 0; i < trackbin; i++){
                  for(int j = 0; j < Classifier_bin; j++){
                     for(int k=0; k <ptbin;k++){

                        if(tkBool[i] +ClassifierBool[j] +A_ptBool[A_trk][k] == 3){
                                
                              if( !NtrigCorrected[i][j][k] ){
                                 std::cout<<"NtrigCorrected[i][j][k]"<<std::endl;
                                 continue;
                              }
                              hEPDrawCor[i][j][k]->Fill(jet_dau_eta, jet_dau_phi, 1.0*jet_HLT_weight/(Atrk_weight * NtrigCorrected[i][j][k] ));//NtrigCorrected is definded within the loop  of jets,count the number of the passed daughter particle in ijet
                              if( (*dau_Flag[beta_index])[ijet][A_trk]==0 ){
                                 hEPDrawCor_groomed[i][j][k]->Fill(jet_dau_eta, jet_dau_phi, 1.0*jet_HLT_weight/(Atrk_weight * Ntrig_g_Corrected[i][j][k] ));
                              }
                              if( (*dau_Flag[beta_index])[ijet][A_trk]    ){
                                 hEPDrawCor_nongroomed[i][j][k]->Fill(jet_dau_eta, jet_dau_phi, 1.0*jet_HLT_weight/(Atrk_weight * Ntrig_ng_Corrected[i][j][k] ));
                              }
                        }
                     }
                  }
               }



               if(A_trk == NNtrk - 1) continue;
                    //***************START T_trk LOOP**********************
                    for(long int T_trk = A_trk+1 ; T_trk< NNtrk ; T_trk ++ ){
                        if(((*dau_chg))[ijet][T_trk] == 0) continue;
                        if(((*dau_pt))[ijet][T_trk] < 0.3) continue;
                        if(fabs(((*dau_eta))[ijet][T_trk]) > 2.4) continue;


                    //calculate pt,eta,phi of daughter particle in jet frame 
                        double T_jet_dau_pt    =  ptWRTJet( (double)(*jetPt)[ijet], (double)(*jetEta)[ijet]     , (double)(*jetPhi)[ijet]       , (double)(*dau_pt)[ijet][T_trk], (double)(*dau_eta)[ijet][T_trk], (double)(*dau_phi)[ijet][T_trk]);
                        double T_jet_dau_eta   = etaWRTJet( (double)(*jetPt)[ijet], (double)(*jetEta)[ijet]     , (double)(*jetPhi)[ijet]       , (double)(*dau_pt)[ijet][T_trk], (double)(*dau_eta)[ijet][T_trk], (double)(*dau_phi)[ijet][T_trk]);
                        double T_jet_dau_phi   = phiWRTJet( (double)(*jetPt)[ijet], (double)(*jetEta)[ijet]     , (double)(*jetPhi)[ijet]       , (double)(*dau_pt)[ijet][T_trk], (double)(*dau_eta)[ijet][T_trk], (double)(*dau_phi)[ijet][T_trk]);

                        if(T_jet_dau_eta > track_eta_lim) continue;
                        double deltaEta = (jet_dau_eta - T_jet_dau_eta);
                        double deltaPhi = (TMath::ACos(TMath::Cos(jet_dau_phi - T_jet_dau_phi)));
                        double Ttrk_weight =1.0; //(hReco2D[thisEffTable]->GetBinContent(hReco2D[thisEffTable]->FindBin( (*dau_pt)[ijet][T_trk] , (*dau_eta)[ijet][T_trk] )));


                            for(        int i = 0; i < trackbin; i++){
                                for(    int j = 0; j < Classifier_bin;  j++){
                                    for( int k=0; k<ptbin ; k++ ){

                                        if(tkBool[i] +ClassifierBool[j] +A_ptBool[A_trk][k] + A_ptBool[T_trk][k] == 4){

                                            
                                          Pairs[i][j][k]+=1;// the number of pairs including both groomed-groomed paris and nongroomed-nongroomed pairs

                                          if(Atrk_weight == 0){/* cout  << "Atrk_weight = 0!" << endl;
                                             */
                                             continue;}
                                          if(Ttrk_weight == 0){/* cout << "Ttrk_weight = 0!" << endl;
                                            */
                                             continue;}
                                          if(NtrigCorrected[i][j][k] == 0) { cout << "NtrigCorrected[i][j] = 0!" << endl; continue;}

                                           

                                       
                                          hSignalShiftedCor[i][j][k]->Fill(deltaEta, deltaPhi,                  1.0*jet_HLT_weight/(Atrk_weight * Ttrk_weight * NtrigCorrected[i][j][k] ));
                                          hSignalShiftedCor[i][j][k]->Fill(-deltaEta, deltaPhi,                 1.0*jet_HLT_weight/(Atrk_weight * Ttrk_weight * NtrigCorrected[i][j][k] ));
                                          hSignalShiftedCor[i][j][k]->Fill(deltaEta, -deltaPhi,                 1.0*jet_HLT_weight/(Atrk_weight * Ttrk_weight * NtrigCorrected[i][j][k] ));
                                          hSignalShiftedCor[i][j][k]->Fill(-deltaEta, -deltaPhi,                1.0*jet_HLT_weight/(Atrk_weight * Ttrk_weight * NtrigCorrected[i][j][k] ));
                                          hSignalShiftedCor[i][j][k]->Fill( deltaEta,2*TMath::Pi() - deltaPhi,  1.0*jet_HLT_weight/(Atrk_weight * Ttrk_weight * NtrigCorrected[i][j][k] ));
                                          hSignalShiftedCor[i][j][k]->Fill(-deltaEta,2*TMath::Pi() - deltaPhi,  1.0*jet_HLT_weight/(Atrk_weight * Ttrk_weight * NtrigCorrected[i][j][k] ));
                                          

                                          // groomed_groomed pairs 
                                          if(   (*dau_Flag[beta_index])[ijet][T_trk]==0  && (*dau_Flag[beta_index])[ijet][A_trk]==0      ){
                                            hSignal_groomed_groomed[i][j][k]->Fill(deltaEta, deltaPhi,                  1.0*jet_HLT_weight/(Atrk_weight * Ttrk_weight * Ntrig_g_Corrected[i][j][k] ));
                                            hSignal_groomed_groomed[i][j][k]->Fill(-deltaEta, deltaPhi,                 1.0*jet_HLT_weight/(Atrk_weight * Ttrk_weight * Ntrig_g_Corrected[i][j][k] ));
                                            hSignal_groomed_groomed[i][j][k]->Fill(deltaEta, -deltaPhi,                 1.0*jet_HLT_weight/(Atrk_weight * Ttrk_weight * Ntrig_g_Corrected[i][j][k] ));
                                            hSignal_groomed_groomed[i][j][k]->Fill(-deltaEta, -deltaPhi,                1.0*jet_HLT_weight/(Atrk_weight * Ttrk_weight * Ntrig_g_Corrected[i][j][k] ));
                                            hSignal_groomed_groomed[i][j][k]->Fill( deltaEta,2*TMath::Pi() - deltaPhi,  1.0*jet_HLT_weight/(Atrk_weight * Ttrk_weight * Ntrig_g_Corrected[i][j][k] ));
                                            hSignal_groomed_groomed[i][j][k]->Fill(-deltaEta,2*TMath::Pi() - deltaPhi,  1.0*jet_HLT_weight/(Atrk_weight * Ttrk_weight * Ntrig_g_Corrected[i][j][k] ));
                                            g_gPairs[i][j][k]++;

                                          }

                                          //nongroomed-nongroomed pairs 
                                          if(   (*dau_Flag[beta_index] )[ijet][T_trk]  &&  (*dau_Flag[beta_index] )[ijet][A_trk]          ){
                                            hSignal_nongroomed_nongroomed[i][j][k]->Fill(deltaEta, deltaPhi,                  1.0*jet_HLT_weight/(Atrk_weight * Ttrk_weight * Ntrig_ng_Corrected[i][j][k] ));
                                            hSignal_nongroomed_nongroomed[i][j][k]->Fill(-deltaEta, deltaPhi,                 1.0*jet_HLT_weight/(Atrk_weight * Ttrk_weight * Ntrig_ng_Corrected[i][j][k] ));
                                            hSignal_nongroomed_nongroomed[i][j][k]->Fill(deltaEta, -deltaPhi,                 1.0*jet_HLT_weight/(Atrk_weight * Ttrk_weight * Ntrig_ng_Corrected[i][j][k] ));
                                            hSignal_nongroomed_nongroomed[i][j][k]->Fill(-deltaEta, -deltaPhi,                1.0*jet_HLT_weight/(Atrk_weight * Ttrk_weight * Ntrig_ng_Corrected[i][j][k] ));
                                            hSignal_nongroomed_nongroomed[i][j][k]->Fill( deltaEta,2*TMath::Pi() - deltaPhi,  1.0*jet_HLT_weight/(Atrk_weight * Ttrk_weight * Ntrig_ng_Corrected[i][j][k] ));
                                            hSignal_nongroomed_nongroomed[i][j][k]->Fill(-deltaEta,2*TMath::Pi() - deltaPhi,  1.0*jet_HLT_weight/(Atrk_weight * Ttrk_weight * Ntrig_ng_Corrected[i][j][k] ));
                                            ng_ngPairs[i][j][k]++;

                                          }




                                       }
                                    }
                                }
                            }





                    }//***************END  T_Trk LOOP**********************
                }//**********END A_trk LOOP 
         }//END JET LOOP
        }//END EVENT LOOP 
        fFile->Close();
    }//END FILE LOOP


        std::cout<<"Start generating background"<<std::endl;    
        for(int wtrk = 1; wtrk < trackbin+1; wtrk++){
            std::cout<< wtrk<<"/"<<trackbin<<std::endl;
            for(int wcla =1;wcla < Classifier_bin+1;wcla++){
                for(int wppt = 1; wppt < ptbin+1; wppt++){
                    
                    gen_background(hEPDrawCor[wtrk-1][wcla-1][wppt-1],hBckrndShiftedCor[wtrk-1][wcla-1][wppt-1],Pairs[wtrk-1][wcla-1][wppt-1] );
                }
            }
        }

        std::cout<<"Start generating groommed  background"<<std::endl;    
        for(int wtrk = 1; wtrk < trackbin+1; wtrk++){
            std::cout<< wtrk<<"/"<<trackbin<<std::endl;
            for(int wcla =1;wcla < Classifier_bin+1;wcla++){
                for(int wppt = 1; wppt < ptbin+1; wppt++){
                    
                    gen_background(hEPDrawCor_groomed[wtrk-1][wcla-1][wppt-1],hBkg_groomed_groomed[wtrk-1][wcla-1][wppt-1],g_gPairs[wtrk-1][wcla-1][wppt-1] );
                }
            }
        }


        std::cout<<"Start  generating non_groomed  background"<<std::endl;    
        for(int wtrk = 1; wtrk < trackbin+1; wtrk++){
            std::cout<< wtrk<<"/"<<trackbin<<std::endl;
            for(int wcla =1;wcla < Classifier_bin+1;wcla++){
                for(int wppt = 1; wppt < ptbin+1; wppt++){
                    
                    gen_background(hEPDrawCor_nongroomed[wtrk-1][wcla-1][wppt-1],hBkg_nongroomed_nongroomed[wtrk-1][wcla-1][wppt-1],ng_ngPairs[wtrk-1][wcla-1][wppt-1] );
                }
            }
        }

        
        std::string fout_dir=store_dir+Form("/hpp_beta_%d.root",(int)beta_SD[beta_index]);



        TFile *fout=new TFile(fout_dir.c_str(),"RECREATE");
        hJet_Pass->Write();
        for(int wtrk =1; wtrk <trackbin+1; wtrk++){
            for(int wcla =1; wcla<Classifier_bin+1; wcla++){
                for(int wppt =1; wppt <ptbin+1; wppt++){
                    hSignalShiftedCor[wtrk-1][wcla-1][wppt-1]->Write( Form("hSigS_Cor_%d_to_%d_%d_to_%d_%d_to_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1] ,(int)(10*ptbinbounds_lo[wppt-1]),(int)(10*ptbinbounds_hi[wppt-1]),(int)(100*Classifier_lo[wcla-1]),(int)(100*Classifier_hi[wcla-1]) ) );
                    hBckrndShiftedCor[wtrk-1][wcla-1][wppt-1]->Write( Form("hBckS_Cor_%d_to_%d_%d_to_%d_%d_to_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1] ,(int)(10*ptbinbounds_lo[wppt-1]),(int)(10*ptbinbounds_hi[wppt-1]),(int)(100*Classifier_lo[wcla-1]),(int)(100*Classifier_hi[wcla-1]) ) );
                    hEPDrawCor[wtrk-1][wcla-1][wppt-1]->Write( Form("hEPD_Cor_%d_to_%d_%d_to_%d_%d_to_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1] ,(int)(10*ptbinbounds_lo[wppt-1]),(int)(10*ptbinbounds_hi[wppt-1]),(int)(100*Classifier_lo[wcla-1]),(int)(100*Classifier_hi[wcla-1]) ) );

                    hEPDrawCor_groomed[wtrk-1][wcla-1][wppt-1]->Write( Form("hEPD_Cor_g_%d_to_%d_%d_to_%d_%d_to_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1] ,(int)(10*ptbinbounds_lo[wppt-1]),(int)(10*ptbinbounds_hi[wppt-1]),(int)(100*Classifier_lo[wcla-1]),(int)(100*Classifier_hi[wcla-1]) )  );
                    hEPDrawCor_nongroomed[wtrk-1][wcla-1][wppt-1]->Write( Form("hEPD_Cor_ng_%d_to_%d_%d_to_%d_%d_to_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1] ,(int)(10*ptbinbounds_lo[wppt-1]),(int)(10*ptbinbounds_hi[wppt-1]),(int)(100*Classifier_lo[wcla-1]),(int)(100*Classifier_hi[wcla-1]) ) );
                    hSignal_groomed_groomed[wtrk-1][wcla-1][wppt-1]->Write( Form("hSigS_gg_Cor_%d_to_%d_%d_to_%d_%d_to_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1] ,(int)(10*ptbinbounds_lo[wppt-1]),(int)(10*ptbinbounds_hi[wppt-1]),(int)(100*Classifier_lo[wcla-1]),(int)(100*Classifier_hi[wcla-1]) ));
                    hSignal_nongroomed_nongroomed[wtrk-1][wcla-1][wppt-1]->Write( Form("hSigS_ngng_Cor_%d_to_%d_%d_to_%d_%d_to_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1] ,(int)(10*ptbinbounds_lo[wppt-1]),(int)(10*ptbinbounds_hi[wppt-1]),(int)(100*Classifier_lo[wcla-1]),(int)(100*Classifier_hi[wcla-1]) ) );
                    hBkg_nongroomed_nongroomed[wtrk-1][wcla-1][wppt-1]->Write( Form("hBckS_ngng_Cor_%d_to_%d_%d_to_%d_%d_to_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1] ,(int)(10*ptbinbounds_lo[wppt-1]),(int)(10*ptbinbounds_hi[wppt-1]),(int)(100*Classifier_lo[wcla-1]),(int)(100*Classifier_hi[wcla-1]) ) );
                    hBkg_groomed_groomed[wtrk-1][wcla-1][wppt-1]->Write( Form("hBckS_gg_Cor_%d_to_%d_%d_to_%d_%d_to_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1] ,(int)(10*ptbinbounds_lo[wppt-1]),(int)(10*ptbinbounds_hi[wppt-1]),(int)(100*Classifier_lo[wcla-1]),(int)(100*Classifier_hi[wcla-1]) ) );

                   

                  
                }
            }
        }

        for(int i=0;i<trackbin;i++){
            for(int j=0;j<Classifier_bin;j++){
               h_Nch_dis[i][j]->Write();
            }
        }
        
        
        
     
    
        fout->Close();





}

























#endif
