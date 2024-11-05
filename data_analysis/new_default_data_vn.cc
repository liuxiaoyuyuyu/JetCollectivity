#define MyClass_cxx

#include "include/MyTrim.h"
#include "include/coordinateTools.h"
//#include "include/high_vn.h"




bool F_eventpass(std::vector< float > *jetPt, int jetnumber, float jetPtCut){

    if(jetPt->size() < 1)           return false;
    if(jetnumber > 200)              return false;
    if((*jetPt)[0] < jetPtCut) return false;
    return true;
}

bool F_jetpass(std::vector< float > * jetEta, std::vector< float > * jetPt, int     ijet, float   jetPtCut){
    if(fabs( (*jetEta)[ijet])   >jetEtaCut)   return false;
    if((*jetPt)[ijet]           <jetPtCut)    return false;
    return true;
}

bool isSubstring(string s1, string s2)
{
    int M = s1.length();
    int N = s2.length();

    /* A loop to slide pat[] one by one */
    for (int i = 0; i <= N - M; i++) {
        int j;

        /* For current index i, check for
         *  pattern match */
        for (j = 0; j < M; j++)
            if (s2[i + j] != s1[j])
                break;

        if (j == M)
            return 1;
    }

    return 0;
}


void MyClass::GenBackGround(long  Npairs,TH2D *hraw,TH2D *hbkg){
    int backMult =10;

    long  NENT =  Npairs;
    long  XENT =  ((1+floor(sqrt(1+(4*2*backMult*NENT))))/2) ;

    //float A_ETA_Cor[XENT] = {0};
    //float A_PHI_Cor[XENT] = {0};
    std::vector<float> A_ETA_Cor(XENT);
    std::vector<float> A_PHI_Cor(XENT);

    for(int x = 0; x<XENT; x++){
        gRandom->SetSeed(0);

        double WEta1_Cor, WPhi1_Cor;//making the pseudoparticles
        //hEPDrawCor[wtrk-1][wppt-1][wpPU-1]->GetRandom2(WEta1_Cor, WPhi1_Cor);
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


void MyClass::Loop(int job, std::string fList){

    TH1::SetDefaultSumw2(kTRUE);
    TH2::SetDefaultSumw2(kTRUE);

    bool smear_bool = 0;

    const int i40 = 40;
    const int i470 = 470;
    const int i3500 = 3500;
    const int i50 = 50;
    const int i500 = 500;
    const int i1500 = 1500;
    const int i150 = 150;
    const int i100 = 100;
    const int PUbin = 1;
    double track_eta_lim = 5.0;
    //Initializing Histograms
    //TH2D* hPairs        = new TH2D("hPairs","hPairs",  trackbin, bin0,trackbin, ptbin, bin0, ptbin);
    TH2D *hPairs[beta_bin][SDclassifer_bin];
    for(int ibeta=0;ibeta<beta_bin;ibeta++){
        for(int iclass=0;iclass<SDclassifer_bin;iclass++){
            hPairs[ibeta][iclass]= new TH2D( Form("hPairs_beta_%d_SDclass_%d_%d",(int)beta_SD[ibeta],(int)(100*SDclassifer_low[iclass]),(int)(100*SDclassifer_high[iclass]) ),Form("hPairs_beta_%d_SDclass_%d_%d",(int)beta_SD[ibeta],(int)(100*SDclassifer_low[iclass]),(int)(100*SDclassifer_high[iclass]) ),  trackbin, bin0,trackbin, ptbin, bin0, ptbin);
        }
    }

    TH2D* hNtrig        = new TH2D("hNtrig","hNtrig",  trackbin, bin0,trackbin, ptbin, bin0, ptbin);
    TH2D* hNtrigCorrected        = new TH2D("hNtrigCorrected","hNtrigCorrected",  trackbin, bin0,trackbin, ptbin, bin0, ptbin);

    TH1D* hEvent_Pass   = new TH1D("hEvent_Pass","hEvent_Pass", trackbin,bin0,trackbin);
    TH1D* hJet_Pass     = new TH1D("hJet_Pass"  ,"hJet_Pass"  , trackbin,bin0,trackbin);
    TH1D* hJet_Pass550     = new TH1D("hJet_Pass550"  ,"hJet_Pass550"  , trackbin,bin0,trackbin);
    TH1D* hJet_Pass550_hltCor     = new TH1D("hJet_Pass550_hltCor"  ,"hJet_Pass550_hltCor"  , trackbin,bin0,trackbin);

    TH2D* h2Jet_Pass[beta_bin]; 
    for(int ibeta=0;ibeta<beta_bin;ibeta++){
       h2Jet_Pass[ibeta]= new TH2D( Form( "h2Jet_Pass_beta_%d",(int)beta_SD[ibeta] )  ,Form("h2Jet_Pass_beta_%d",(int)beta_SD[ibeta]) , trackbin,bin0,trackbin,SDclassifer_bin,bin0,SDclassifer_bin);
    }


    TH1D* hJetPt = new TH1D("hJetPt"  ,"hJetPt" ,i150,i100,i3500);
    TH1D* hJetPt_wo_ptcut = new TH1D("hJetPt_wo_ptcut"  ,"hJetPt_wo_ptcut" ,i150,i100,i3500);
    TH1D* hBinDist_cor_single = new TH1D("hBinDist_cor_single","hBinDist_cor_single",bin360,bin0,bin120);
    TH1D* hBinDist_unc_single = new TH1D("hBinDist_unc_single","hBinDist_unc_single",bin120,bin0,bin120);

    TH2D* h_lab_JetMult_pT=new TH2D("lab_JetMult_pT","lab_JetMult_pT",300,100,3500,120,0.5,120.5);
    TH2D* h_lab_JetMult_phi=new TH2D("lab_JetMult_phi","lab_JetMult_phi",30,-TMath::Pi(),TMath::Pi(),120,0.5,120.5);
    TH2D* h_lab_JetMult_eta=new TH2D("lab_JetMult_eta","lab_JetMult_eta",34,-1.7,1.7,120,0.5,120.5);

    TH1D* h_jet_jT=new TH1D("jet_jT","jet_jT",200,0,20);
    TH1D* h_jet_etastar=new TH1D("jet_etastar","jet_etastar",100,0,10);

    TH2D* hEPDrawCor[beta_bin][SDclassifer_bin][trackbin][ptbin][PUbin];
    TH2D* hSignalShiftedCor[beta_bin][SDclassifer_bin][trackbin][ptbin][PUbin];
    TH2D* hBckrndShiftedCor[beta_bin][SDclassifer_bin][trackbin][ptbin][PUbin];

    TH1D* hTotalWeight                = new TH1D("hTotalWeight","hTotalWeight"                          , bin500,-bin2,bin3);
    TH1D* hAvg_Atrk_Weight            = new TH1D("hAvg_Atrk_Weight","hAvg_Atrk_Weight"                  , bin500,-bin2,bin3);
    TH1D* hAvg_NtrigCorrected_Bump    = new TH1D("hAvg_NtrigCorrected_Bump","hAvg_NtrigCorrected_Bump"  , bin22 ,-bin2,bin20);

    TH1D* hBinDist_cor[beta_bin][SDclassifer_bin][trackbin];
    TH1D* hBinDist_unc[beta_bin][SDclassifer_bin][trackbin];

    for(int wtrk = 1; wtrk<trackbin+1; wtrk++){
        //hBinDist_cor[wtrk-1]    = new TH1D(Form("hBinDist_cor_%d",wtrk),Form("hBinDist_cor_%d",wtrk), bin360, bin0, bin120);
        //hBinDist_unc[wtrk-1]    = new TH1D(Form("hBinDist_unc_%d",wtrk),Form("hBinDist_unc_%d",wtrk), bin120, bin0, bin120);
        for(int ibeta=1; ibeta<beta_bin+1;ibeta++){
            for(int iclass=1;iclass<SDclassifer_bin+1;iclass++){
                hBinDist_cor[ibeta-1][iclass-1][wtrk-1]=new TH1D(Form("hBinDist_cor_%d_%d_%d",ibeta,iclass,wtrk),Form("hBinDist_cor_%d_%d_%d",ibeta,iclass,wtrk),bin360,bin0,bin120) ;
                hBinDist_unc[ibeta-1][iclass-1][wtrk-1]=new TH1D(Form("hBinDist_unc_%d_%d_%d",ibeta,iclass,wtrk),Form("hBinDist_unc_%d_%d_%d",ibeta,iclass,wtrk),bin120,bin0,bin120) ;
                 for(int wppt = 1; wppt<ptbin+1; wppt++){
                     for(int wpPU = 1; wpPU<PUbin+1; wpPU++){
                        hBckrndShiftedCor[ibeta-1][iclass-1][wtrk-1][wppt-1][wpPU-1]   = new TH2D(Form("hBckrndS_Cor_Nch_%d_%d_beta_%d_SDClass_%d_%d_pt_%d_%d_PU_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1],(int)beta_SD[ibeta-1],(int)(100*SDclassifer_low[iclass-1]),(int)(100*SDclassifer_high[iclass-1]),(int)(10*ptbinbounds_lo[wppt-1]),(int)(10*ptbinbounds_hi[wppt-1]),wpPU ) ,Form("hBckrndS_Cor_Nch_%d_%d_beta_%d_SDClass_%d_%d_pt_%d_%d_PU_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1],(int)beta_SD[ibeta-1],(int)(100*SDclassifer_low[iclass-1]),(int)(100*SDclassifer_high[iclass-1]),(int)(10*ptbinbounds_lo[wppt-1]),(int)(10*ptbinbounds_hi[wppt-1]),wpPU ) ,41,-(20*EtaBW)-(0.5*EtaBW),(20*EtaBW)+(0.5*EtaBW),33,-(8*PhiBW)-0.5*PhiBW,(24*PhiBW)+0.5*PhiBW);
                        hSignalShiftedCor[ibeta-1][iclass-1][wtrk-1][wppt-1][wpPU-1]   = new TH2D(Form("hSignalS_Cor_Nch_%d_%d_beta_%d_SDClass_%d_%d_pt_%d_%d_PU_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1],(int)beta_SD[ibeta-1],(int)(100*SDclassifer_low[iclass-1]),(int)(100*SDclassifer_high[iclass-1]),(int)(10*ptbinbounds_lo[wppt-1]),(int)(10*ptbinbounds_hi[wppt-1]),wpPU ) ,Form("hSignalS_Cor_Nch_%d_%d_beta_%d_SDClass_%d_%d_pt_%d_%d_PU_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1],(int)beta_SD[ibeta-1],(int)(100*SDclassifer_low[iclass-1]),(int)(100*SDclassifer_high[iclass-1]),(int)(10*ptbinbounds_lo[wppt-1]),(int)(10*ptbinbounds_hi[wppt-1]),wpPU ) ,41,-(20*EtaBW)-(0.5*EtaBW),(20*EtaBW)+(0.5*EtaBW),33,-(8*PhiBW)-0.5*PhiBW,(24*PhiBW)+0.5*PhiBW);
                               hEPDrawCor[ibeta-1][iclass-1][wtrk-1][wppt-1][wpPU-1]   = new TH2D( Form("hEPDraw_Cor_Nch_%d_%d_beta_%d_SDClass_%d_%d_pt_%d_%d_PU_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1],(int)beta_SD[ibeta-1],(int)(100*SDclassifer_low[iclass-1]),(int)(100*SDclassifer_high[iclass-1]),(int)(10*ptbinbounds_lo[wppt-1]),(int)(10*ptbinbounds_hi[wppt-1]),wpPU ) ,Form("hEPDraw_Cor_Nch_%d_%d_beta_%d_SDClass_%d_%d_pt_%d_%d_PU_%d" ,trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1],(int)beta_SD[ibeta-1],(int)(100*SDclassifer_low[iclass-1]),(int)(100*SDclassifer_high[iclass-1]),(int)(10*ptbinbounds_lo[wppt-1]),(int)(10*ptbinbounds_hi[wppt-1]),wpPU )  , EPD_xb   , EPD_xlo, EPD_xhi , EPD_yb      , EPD_ylo    , EPD_yhi);
                    }
                }
            }
        }
   }

    TH2D* hjet_eta_phi_before_veto=new TH2D(Form("jet_eta_phi_before_veto"),Form("jet_eta_phi_before_veto"),82,-5.191,5.191,72,-TMath::Pi(),TMath::Pi());
    
    TH2D* hjet_eta_phi_after_veto=new TH2D(Form("jet_eta_phi_after_veto"),Form("jet_eta_phi_after_veto"),82,-5.191,5.191,72,-TMath::Pi(),TMath::Pi());


    //Jet SoftDrop Subtructure

    //for all Nch
    TH2D *hNch_Zg[beta_bin];                                                                                              
    TH2D *hNch_Rg[beta_bin];                                                                                                
    TH2D *hNch_ZgTgB[beta_bin];                                                                                            
    TH2D *hZg_Rg[beta_bin];

    TH1D *hZg_dis[beta_bin];
    TH1D *hRg_dis[beta_bin];
    TH1D *hZgTgB_dis[beta_bin];
    TH1D *hSDJetMass_dis[beta_bin];

    //for different Nch
    TH1D *hZgTgB_dis_Nch[beta_bin][trackbin];                                                                                       
    TH1D *hZg_dis_Nch[beta_bin][trackbin];                                                                                     
    TH1D *hRg_dis_Nch[beta_bin][trackbin];  
    TH1D *hSDJetMass_dis_Nch[beta_bin][trackbin]; 

    TH1D *hStoredSDJetMassNch=new TH1D("hStoredSDJetMassNch","hStoredSDJetMassNch",300,0,300); 
    TH1D *hCalSDJetMassNch=new TH1D("hCalSDJetMassNch","hCalSDJetMassNch",300,0,300); 

    for(int ibeta=0;ibeta<beta_bin;ibeta++){




        hNch_ZgTgB[ibeta]=new TH2D(Form("Nch_VS_ZgThetagBeta_beta_%d",(int)beta_SD[ibeta]),Form("Nch_VS_ZgThetagBeta_beta_%d",(int)beta_SD[ibeta]), 120 , 0 , 120 , 105 , -0.05 , 1.0 );
        hNch_Rg[ibeta]=new TH2D(Form("Nch_VS_R_g_beta_%d",(int)beta_SD[ibeta]),Form("Nch_VS_R_g_beta_%d",(int)beta_SD[ibeta])                    , 120 , 0 , 120 , 50  , -0.05 , 1   );
        hNch_Zg[ibeta]=new TH2D(Form("Nch_VS_Z_g_beta_%d",(int)beta_SD[ibeta]),Form("Nch_VS_Z_g_beta_%d",(int)beta_SD[ibeta])                    , 120 , 0 , 120 , 50  , -0.05 , 0.6 );
        hZg_Rg[ibeta]=new TH2D( Form("hZg_Vs_Rg_beta_%d",(int)beta_SD[ibeta]) ,Form("hZg_Vs_Rg_beta_%d",(int)beta_SD[ibeta])          ,50 ,   -0.05   ,0.6    ,50   ,-0.05   ,1   ); 
        
        hZg_dis[ibeta]=new TH1D(Form("hZg_dis_beta_%d",  (int)beta_SD[ibeta]) ,Form("hZg_dis_beta_%d",  (int)beta_SD[ibeta])          ,50 ,   -0.05   ,0.55   );
        hRg_dis[ibeta]=new TH1D(Form("hRg_dis_beta_%d",  (int)beta_SD[ibeta]) ,Form("hRg_dis_beta_%d",  (int)beta_SD[ibeta])          ,50 ,   -0.05   ,0.55   );
        hZgTgB_dis[ibeta]=new TH1D(Form("hZgTgB_dis_beta_%d",  (int)beta_SD[ibeta]) ,Form("hZgTgB_dis_beta_%d",  (int)beta_SD[ibeta]) ,100,   0       ,1.0    );
        hSDJetMass_dis[ibeta]=new TH1D(Form("hSDJetMass_dis_beta_%d",(int)beta_SD[ibeta]),Form("hSDJetMass_dis_beta_%d",(int)beta_SD[ibeta]), 4000,0. ,4000   );

        for(int iNch=0;iNch<trackbin;iNch++){
            hZgTgB_dis_Nch[ibeta][iNch]=new TH1D(Form("ZgThetagBeta_dis_Nch_%d_%d_beta_%d",trackbinbounds[iNch],trackbinboundsUpper[iNch],(int)beta_SD[ibeta]),Form("ZgThetagBeta_dis_Nch_%d_%d_beta_%d",trackbinbounds[iNch],trackbinboundsUpper[iNch],(int)beta_SD[ibeta]),100,0,1.0);
            hZg_dis_Nch[ibeta][iNch]=new TH1D(Form("ZgDis_Nch_%d_%d_beta_%d",trackbinbounds[iNch],trackbinboundsUpper[iNch],(int)beta_SD[ibeta]),Form("ZgDis_Nch_%d_%d_beta_%d",trackbinbounds[iNch],trackbinboundsUpper[iNch],(int)beta_SD[ibeta]),50,-0.05,0.55);
            hRg_dis_Nch[ibeta][iNch]=new TH1D(Form("RgDis_Nch_%d_%d_beta_%d",trackbinbounds[iNch],trackbinboundsUpper[iNch],(int)beta_SD[ibeta]),Form("RgDis_Nch_%d_%d_beta_%d",trackbinbounds[iNch],trackbinboundsUpper[iNch],(int)beta_SD[ibeta]),50,-0.05,0.8);
            hSDJetMass_dis_Nch[ibeta][iNch]=new TH1D(Form("hSDJetMass_dis_beta_%d_Nch_%d_%d",(int)beta_SD[ibeta],trackbinbounds[iNch],trackbinboundsUpper[iNch]),Form("hSDJetMass_dis_beta_%d_Nch_%d_%d",(int)beta_SD[ibeta],trackbinbounds[iNch],trackbinboundsUpper[iNch]),4000,0,4000);

        }
    }



    //NEW THING                        
    //NEW THING                        
    //NEW THING                         
    /*
    TH2D* hReco2D[fileList.size()];
    TH2D* hGen2D[fileList.size()];
    TH1D* hdid500;
    TH1D* hdid400;

    TFile *f_jet_HLT_lookup = new TFile("did400500_v2_all.root");
    hdid500 = (TH1D*)f_jet_HLT_lookup->Get("d500")->Clone("did500");
    hdid400 = (TH1D*)f_jet_HLT_lookup->Get("d400")->Clone("did400");
    hdid500->Divide(hdid400);
    */
    std::cout<<"Init SoftDrop"<<std::endl;
    InitSoftDrop();

    //************************************************************************
    fastjet::JetDefinition jet_def(fastjet::cambridge_algorithm ,fastjet::JetDefinition::max_allowable_R);
    //************************************************************************
    std::cout << "Starting event loop" << std::endl;
    std::cout << "Total Number of Files in this Job: " << fileList.size() << std::endl;

    //*******************************************ENTER FILE LOOP**********************************************************
    for(int f = 0; f<fileList.size(); f++){
        int f_from_file = f;
        fFile = TFile::Open(fileList.at(f).c_str(),"read");
        TTree *tree = (TTree*)fFile->Get("analyzerOffline/trackTree");
        Init(tree);

        std::cout << "File " << f+1 << " out of " << fileList.size() << std::endl;
        Long64_t nbytes = 0, nb = 0;
        Long64_t nentries = fChain->GetEntriesFast();
        cout<<"Total Entries is:"<<endl;
        cout<< nentries <<endl;

        /*
        vector<string> era_vec;
        vector<string> matched_cor_table_vec;
        era_vec.push_back("2018D_"); matched_cor_table_vec.push_back("18_v2");
        era_vec.push_back("2018C_"); matched_cor_table_vec.push_back("18_v2");
        era_vec.push_back("2018B_"); matched_cor_table_vec.push_back("18_v2");
        era_vec.push_back("2018A_"); matched_cor_table_vec.push_back("18_v2");
        era_vec.push_back("2017B_"); matched_cor_table_vec.push_back("17_v2");
        era_vec.push_back("2017C_"); matched_cor_table_vec.push_back("17_v2");
        era_vec.push_back("2017D_"); matched_cor_table_vec.push_back("17_v2");
        era_vec.push_back("2017E_"); matched_cor_table_vec.push_back("17_v2");
        era_vec.push_back("2017F_"); matched_cor_table_vec.push_back("17_v2");
        era_vec.push_back("2016Bv2_"); matched_cor_table_vec.push_back("16apv_v2");
        era_vec.push_back("2016C_"); matched_cor_table_vec.push_back("16apv_v2");
        era_vec.push_back("2016D_"); matched_cor_table_vec.push_back("16apv_v2");
        era_vec.push_back("2016E_"); matched_cor_table_vec.push_back("16apv_v2");
        era_vec.push_back("2016F-HIPM_"); matched_cor_table_vec.push_back("16apv_v2");
        era_vec.push_back("2016F_"); matched_cor_table_vec.push_back("16_v2");
        era_vec.push_back("2016G_"); matched_cor_table_vec.push_back("16_v2");
        era_vec.push_back("2016H_"); matched_cor_table_vec.push_back("16_v2");
        int i_keep=999;
        for(int i=0; i<era_vec.size(); i++){
            if(isSubstring(era_vec[i],fileList.at(f).c_str())){
                i_keep=i;
                continue;
            }
        }

        string filename = "~/StorageArea/"+matched_cor_table_vec[i_keep]+".root";
        TFile *f_pt_eta_DCA_lookup = new TFile(filename.c_str());

        hReco2D[f] = (TH2D*)f_pt_eta_DCA_lookup->Get("h2_MAIN_DCA_Dau_Reco_Pt_Eta_Lab_All")->Clone(Form("h2_DCA_Dau_Reco_Pt_Eta_Lab_All_%d",f));
        hGen2D[f] = (TH2D*)f_pt_eta_DCA_lookup->Get("h2_Dau_Gen_Pt_Eta_Lab_All")->Clone(Form("h2_Dau_Gen_Pt_Eta_Lab_All_%d",f));
        hReco2D[f]->Divide(hGen2D[f]);
        int thisEffTable =f_from_file;
        */

        //TFile* jet_veto_file[2];
        //jet_veto_file[0]=new TFile("~/StorageArea/Summer22_23Sep2023_RunCD_v1.root","read");
        //jet_veto_file[1]=new TFile("~/StorageArea/Summer22EE_23Sep2023_RunEFG_v1.root","read");
        //TH2D* jet_veto_map=(TH2D*)jet_veto_file[0]->Get("jetvetomap");

        
        //********************************************ENTERING EVENT LOOP**************************************************
        for (Long64_t ievent=0; ievent <nentries; ievent ++){
            Long64_t jevent = LoadTree(ievent);
            nb = fChain->GetEntry(ievent);   nbytes += nb;

            //if(!F_eventpass(jetPt, jetN, jetPtCut_Event)){
            //Use jetPtCut_Dist=100 for QA plots.
            if(!F_eventpass(jetPt, jetN, jetPtCut_Dist)){
                continue;
            }
            
            hEvent_Pass->Fill(1);

            int jetCounter = jetPt->size();
            if(jetCounter == 0) continue;

            /*
            //apply jet veto map
            //jet veto, reference: https://cms-jerc.web.cern.ch/Recommendations/#run-3
            //"The safest procedure would be to veto events if ANY jet with a loose selection lies in the veto regions."
            bool jetvetoBool=0;
            for(int ijet=0; ijet < jetCounter; ijet++){
                hjet_eta_phi_before_veto->Fill((*jetEta)[ijet],(*jetPhi)[ijet]);
                
                //if( (*jetPt)[ijet] > 15){
                //    if(jet_veto_map->GetBinContent(jet_veto_map->FindBin((*jetEta)[ijet],(*jetPhi)[ijet]))>0){
                //        jetvetoBool=1;
                //        break;
                //    } 
                //}

            }

            if(jetvetoBool) continue;
            */

            //**********************************************ENTERING JET LOOP*****************************************
            for(int ijet=0; ijet < jetCounter; ijet++){

                //std::cout<<"*****************"<<std::endl;
                //std::cout<<"stored softjet mass="<<(*jetSoftDropMass)[ijet]<<std::endl;
                hStoredSDJetMassNch->Fill( (*jetSoftDropMass)[ijet] );

                //the vector that will store all daughter particles 
                std::vector< fastjet::PseudoJet > particles; 
                //fille all the daughter particles;
                for(int idau=0;idau<(*dau_pt)[ijet].size();idau++){
                    TLorentzVector v;
                    int current_dau_pid=(*dau_pid)[ijet][idau];
                    //v.SetPtEtaPhiM((*dau_pt)[ijet][idau],(*dau_eta)[ijet][idau],(*dau_phi)[ijet][idau],particleData.m0(current_dau_pid));
                    v.SetPtEtaPhiM((*dau_pt)[ijet][idau],(*dau_eta)[ijet][idau],(*dau_phi)[ijet][idau],(*dau_mass)[ijet][idau] );
                    fastjet::PseudoJet particle( v.Px(),v.Py(),v.Pz(),v.E() );
                    particles.push_back(particle);
                }

                //define CA alogrithm
                fastjet::ClusterSequence cs(particles,jet_def);
                vector<fastjet::PseudoJet> jets = cs.inclusive_jets();//only one jet shoule be resonable 

                if(jets.size()>1) std::cout<<"waring with CA "<<std::endl;

                std::vector<double > Zgs;
                std::vector<double > Rgs;
                std::vector<double > ZgTgBs;
                std::vector<double > SDJetMass; 
                for(int ibeta=0;ibeta<beta_bin;ibeta++){
                    fastjet::PseudoJet sd_jet = (*(sd[ibeta]) ) ( jets[0] );

                    if( !sd_jet.has_structure_of<fastjet::contrib::SoftDrop>() ){
                        std::cout<<"fail in SoftDrop"<<std::endl;
                        Zgs.push_back(-1);
                        Rgs.push_back(-1);
                        ZgTgBs.push_back(-1);
                        SDJetMass.push_back(-1);
                        continue;
                    } 

                    fastjet::PseudoJet parent1;
                    fastjet::PseudoJet parent2;

                    if(!sd_jet.has_parents(parent1,parent2) ){
                        Zgs.push_back(-1);
                        Rgs.push_back(-1);
                        ZgTgBs.push_back(-1);
                        SDJetMass.push_back(sd_jet.m());
                        continue;
                    }

                    Zgs.push_back(    sd_jet.structure_of<fastjet::contrib::SoftDrop>().symmetry()  );
                    Rgs.push_back(    sd_jet.structure_of<fastjet::contrib::SoftDrop>().delta_R()   );
                    ZgTgBs.push_back( Zgs[ibeta]*TMath::Power(Rgs[ibeta]/R,beta_SD[ibeta]) );
                    SDJetMass.push_back(sd_jet.m());
                    //std::cout<<"calculated softdrop jet mass:"<<std::endl;
                    //std::cout<<"beta="<<beta_SD[ibeta]<<" "<<" softdrop jet mass="<<sd_jet.m()<<std::endl;
                    if(beta_SD[ibeta]==0) hCalSDJetMassNch->Fill( sd_jet.m() );

                }
                //std::cout<<"**********************"<<std::endl;
                //fill substructure histogram:do this after indentify Nch and cut on jet pt




                hjet_eta_phi_after_veto->Fill((*jetEta)[ijet],(*jetPhi)[ijet]);
                
                long int NNtrk = (dau_pt->at(ijet)).size();
                gRandom->SetSeed(0);
                double eta_smear;
                eta_smear=0;
                if( fabs(((*jetEta)[ijet])+eta_smear) > jetEtaCut ) continue;//*****************************CUT On JEt ETA*************************
                //if( (*jetPt)[ijet] < jetPtCut_Jet   ) continue;
                //hJetPt->Fill((*jetPt)[ijet]);
                // filling distrivutions within track bins
                // ALSO VERY IMPORTANLTY changing the tkBool to 1 for this particular jet. This will be usefull later wen I create conditons for filling other historgams.

                double jet_HLT_weight = 1.0;
                /*
                if((*jetPt)[ijet] < 880 && (*jetPt)[ijet] > 550){
                    jet_HLT_weight = 1.0/(hdid500->GetBinContent(hdid500->FindBin((*jetPt)[ijet]) ) );
                }
                */
                int tkBool[trackbin] = {0};
                //int sdBool[beta_bin][SDclassifer_bin]={0};

                std::vector<std::vector<int>> sdBool(beta_bin,std::vector<int>(SDclassifer_bin)) ;



                int n_ChargeMult_DCA_labPt_Eta_exclusion =0;
                double n_ChargeMult_DCA_labPt_Eta_exclusion_Cor =0;

                for(int  A_trk=0; A_trk < NNtrk; A_trk++ ){
                    if((*dau_chg)[ijet][A_trk] == 0) continue;//charge
                    if(fabs((*dau_pt)[ijet][A_trk])  < 0.3)     continue;//lab pt
                    if(fabs((*dau_eta)[ijet][A_trk]) > 2.4)     continue;//lab eta

                    double dauptnow = (*dau_pt)[ijet][A_trk]; double pterr_R_pt =0; double dcaXY = 0; double dcaZ = 0;
                    pterr_R_pt= ( (*dau_ptError)[ijet][A_trk] ) / ( (*dau_pt)[ijet][A_trk] );
                    dcaXY = (*dau_XYDCAsig)[ijet][A_trk];
                    dcaZ = (*dau_ZDCAsig)[ijet][A_trk];

                    if(pterr_R_pt   > 0.1   && dauptnow > 0.5)  continue;
                    if(fabs(dcaXY)  > 3     && dauptnow > 0.5)  continue;
                    if(fabs(dcaZ)   > 3     && dauptnow > 0.5)  continue;
                    n_ChargeMult_DCA_labPt_Eta_exclusion += 1;
                    //double nUnc_weight = (hReco2D[thisEffTable]->GetBinContent(hReco2D[thisEffTable]->FindBin( (*dau_pt)[ijet][A_trk] , (*dau_eta)[ijet][A_trk] )));
                    //n_ChargeMult_DCA_labPt_Eta_exclusion_Cor += (1.0/nUnc_weight);
                }


                


                h_lab_JetMult_pT->Fill( (*jetPt)[ijet],n_ChargeMult_DCA_labPt_Eta_exclusion);
                h_lab_JetMult_phi->Fill((*jetPhi)[ijet],n_ChargeMult_DCA_labPt_Eta_exclusion);
                h_lab_JetMult_eta->Fill((*jetEta)[ijet],n_ChargeMult_DCA_labPt_Eta_exclusion);

                hJetPt_wo_ptcut->Fill((*jetPt)[ijet]);    
                if( (*jetPt)[ijet] < jetPtCut_Jet   ) continue;//****************************Cut ON Jet PT**********************
                hJetPt->Fill((*jetPt)[ijet]);







                //********************************************start fill SoftDrop subtructure histgrom**************************************************
                for(int ibeta=0;ibeta<beta_bin;ibeta++){
                    hNch_Zg[ibeta]->Fill(n_ChargeMult_DCA_labPt_Eta_exclusion    , Zgs[ibeta]                                                );                                                                                              
                    hNch_Rg[ibeta]->Fill(n_ChargeMult_DCA_labPt_Eta_exclusion    , Rgs[ibeta]                                                );                                                                                                
                    hNch_ZgTgB[ibeta]->Fill(n_ChargeMult_DCA_labPt_Eta_exclusion , Zgs[ibeta]*TMath::Power(Rgs[ibeta]/R,beta_SD[ibeta])      );                                                                                            
                    hZg_Rg[ibeta]->Fill(   Zgs[ibeta]                            , Rgs[ibeta]                                                );

                    hZg_dis[ibeta]->Fill(    Zgs[ibeta]           );
                    hRg_dis[ibeta]->Fill(    Rgs[ibeta]           );
                    hZgTgB_dis[ibeta]->Fill( ZgTgBs[ibeta]        );
                    hSDJetMass_dis[ibeta]->Fill( SDJetMass[ibeta] );
   
    

                    for(int iNch=0;iNch<trackbin;iNch++){
                        if( n_ChargeMult_DCA_labPt_Eta_exclusion >= trackbinbounds[iNch] && n_ChargeMult_DCA_labPt_Eta_exclusion < trackbinboundsUpper[iNch]){
                            hZgTgB_dis_Nch[ibeta][iNch]->Fill(  ZgTgBs[ibeta]           );
                            hZg_dis_Nch[ibeta][iNch]->Fill(    Zgs[ibeta]               );
                            hRg_dis_Nch[ibeta][iNch]->Fill(    Rgs[ibeta]               );
                            hSDJetMass_dis_Nch[ibeta][iNch]->Fill(     SDJetMass[ibeta] );

                            

                        }


                    }
                }
                //********************************************end    fill SoftDrop subtructure histgrom**************************************************

                std::vector<double> current_jetClassifier;
                if(mode==1) current_jetClassifier=Zgs;
                else if(mode==2) current_jetClassifier=Rgs;
                else if(mode==3) current_jetClassifier=ZgTgBs;
                else {
                    std::cout<<"Please select correct cutmode"<<std::endl;
                    return ;
                }

                //hBinDist_cor_single            ->Fill(n_ChargeMult_DCA_labPt_Eta_exclusion_Cor, 1.0*jet_HLT_weight);
                hBinDist_unc_single            ->Fill(n_ChargeMult_DCA_labPt_Eta_exclusion, 1.0*jet_HLT_weight);
                for(int i = 0; i < trackbin; i++){
                    if( n_ChargeMult_DCA_labPt_Eta_exclusion >= trackbinbounds[i] && n_ChargeMult_DCA_labPt_Eta_exclusion < trackbinboundsUpper[i]){
                        tkBool[i] = 1;
                        hJet_Pass                   ->Fill(i);
                        if((*jetPt)[ijet] >= jetPtCut_Event) hJet_Pass550                           ->Fill(i);
                        if((*jetPt)[ijet] >= jetPtCut_Event) hJet_Pass550_hltCor                   ->Fill(i,1.0*jet_HLT_weight);

                        //hBinDist_cor[i]             ->Fill(n_ChargeMult_DCA_labPt_Eta_exclusion_Cor,1.0*jet_HLT_weight);
                        //hBinDist_unc[i]             ->Fill(n_ChargeMult_DCA_labPt_Eta_exclusion,1.0*jet_HLT_weight);
                    }
                }
                for(int ibeta=0;ibeta<beta_bin;ibeta++){
                    for(int iclass=0;iclass<SDclassifer_bin;iclass++){
                        if( (current_jetClassifier)[ibeta] >= SDclassifer_low[iclass] && (current_jetClassifier)[ibeta] < SDclassifer_high[iclass]  ){
                            //h2Jet_Pass[ibeta]->Fill()
                            sdBool[ibeta][iclass]=1;
                        }
                    }   
                }


                for(int ibeta=0;ibeta<beta_bin;ibeta++){
                    for(int iclass=0;iclass<SDclassifer_bin;iclass++){
                        for(int iNch=0;iNch<trackbin;iNch++){
                            if(n_ChargeMult_DCA_labPt_Eta_exclusion >= trackbinbounds[iNch] && n_ChargeMult_DCA_labPt_Eta_exclusion < trackbinboundsUpper[iNch] && (current_jetClassifier)[ibeta] >= SDclassifer_low[iclass] && (current_jetClassifier)[ibeta] < SDclassifer_high[iclass]  ){
                                h2Jet_Pass[ibeta]->Fill(iNch,iclass);
                                hBinDist_cor[ibeta][iclass][iNch]             ->Fill(n_ChargeMult_DCA_labPt_Eta_exclusion_Cor,1.0*jet_HLT_weight);
                                hBinDist_unc[ibeta][iclass][iNch]            ->Fill(n_ChargeMult_DCA_labPt_Eta_exclusion,1.0*jet_HLT_weight);
                            }


                        }
                    }
                }
                   


                int Ntrig[trackbin][ptbin] = {0}; 
                double NtrigCorrected[trackbin][ptbin] = {0};
                int A_ptBool[NNtrk][ptbin] = {0};

                //int Ntrig[beta_bin][SDclassifer_bin][trackbin][ptbin]={0};                There is no nedd to expand Ntrig and NtrigCorreced,since no matter what the SDclassifier of a jet is, this jet always has the same Ntrig for given pt range.In addition, track should also be redunant
                //double hNtrigCorrected[beta_bin][SDclassifer_bin][trackbin][ptbin]={0};


                // VERY IMPORTANT calculating daughter pt wrt to jet axis.
                // So this needs to be 2d vector, for pt bin and for daughter index.
                // In each case I create a true falsse for that daughter falling in to the specific pt bin.
                // Note this is NOT jet x daughter. It's pt bin x daughtr
                //****************************************START A_trk LOOP**********************************
                for(int  A_trk=0; A_trk < NNtrk; A_trk++ ){
                    if((*dau_chg)[ijet][A_trk] == 0) continue;
                    if((*dau_pt)[ijet][A_trk] < 0.3) continue;
                    if(fabs((*dau_eta)[ijet][A_trk]) > 2.4) continue;

                    double dauptnow = (*dau_pt)[ijet][A_trk];
                    double pterr_R_pt =0;
                    double dcaXY = 0;
                    double dcaZ = 0;
                    pterr_R_pt= ( (*dau_ptError)[ijet][A_trk] ) / ( (*dau_pt)[ijet][A_trk] );
                    dcaXY = (*dau_XYDCAsig)[ijet][A_trk];
                    dcaZ = (*dau_ZDCAsig)[ijet][A_trk];
                    if(pterr_R_pt   > 0.1   && dauptnow > 0.5) continue;
                    if(fabs(dcaXY)  > 3     && dauptnow > 0.5) continue;
                    if(fabs(dcaZ)   > 3     && dauptnow > 0.5) continue;

                    double jet_dau_pt = ptWRTJet((double)(*jetPt)[ijet], (double)(*jetEta)[ijet], (double)(*jetPhi)[ijet], (double)(*dau_pt)[ijet][A_trk], (double)(*dau_eta)[ijet][A_trk], (double)(*dau_phi)[ijet][A_trk]);

                    for(int i = 0; i < ptbin; i++){
                        if(jet_dau_pt >= ptbinbounds_lo[i] && jet_dau_pt < ptbinbounds_hi[i]){
                            A_ptBool[A_trk][i] = 1;
                        }
                    }

                    //NEW THING                        
                    //NEW THING                        
                    //NEW THING                        

                    //getting et pt efficiency for A track in beam frame
                    //double Atrk_weight = (hReco2D[thisEffTable]->GetBinContent(hReco2D[thisEffTable]->FindBin( (*dau_pt)[ijet][A_trk] , (*dau_eta)[ijet][A_trk] )));
                    double Atrk_weight=1.0; 

                    //in the case where the total jet mult and the individual daughter pt is acceptable for this track bin and pt bin, we increase the Ntrig count.
                    for(int i = 0; i < trackbin; i++){
                        for(int j = 0; j < ptbin; j++){
                            if(tkBool[i] + A_ptBool[A_trk][j] == 2){
                                NtrigCorrected[i][j] += (1.0/Atrk_weight);
                                Ntrig[i][j] += 1;
                            }
                        }
                    }




                }//here ends the boolean array creations.
                //********************************************END A_trk LOOP*********************************************
                //continuation of main loops. Here is where the 2D Corr plots are created using the above booleans and 
                for(int i = 0; i < trackbin; i++){
                    for(int j = 0; j < ptbin; j++){
                        hNtrig->Fill(i,j,Ntrig[i][j]);
                        //hNtrigCorrected->Fill(i,j, NtrigCorrected[i][j]);
                        //hAvg_NtrigCorrected_Bump->Fill(NtrigCorrected[i][j] - Ntrig[i][j]);
                    }
                }
                //******************************************START A_trk LOOP*****************************************************
                for(int  A_trk=0; A_trk < NNtrk; A_trk++ ){
                    //this should be redundant if it passes the bools above? i guess it helps skip daughters faster. maybe i can reindex and run through the daughters quickly by aranging all the charged dauhghter sat the front.
                    if((*dau_chg)[ijet][A_trk] == 0) continue;
                    if((*dau_pt)[ijet][A_trk] < 0.3) continue;
                    if(fabs((*dau_eta)[ijet][A_trk]) > 2.4) continue;
                    //if((*highPurity)[ijet][A_trk] == 0) continue;

                    double dauptnow = (*dau_pt)[ijet][A_trk];
                    double pterr_R_pt =0;
                    double dcaXY = 0;
                    double dcaZ = 0;
                    pterr_R_pt= ( (*dau_ptError)[ijet][A_trk] ) / ( (*dau_pt)[ijet][A_trk] );
                    dcaXY = (*dau_XYDCAsig)[ijet][A_trk];
                    dcaZ = (*dau_ZDCAsig)[ijet][A_trk];
                    if(pterr_R_pt   > 0.1   && dauptnow > 0.5) continue;
                    if(fabs(dcaXY)  > 3     && dauptnow > 0.5) continue;
                    if(fabs(dcaZ)   > 3     && dauptnow > 0.5) continue;

                    gRandom->SetSeed(0);
                    double phi_smear;
                        phi_smear=0;

                    double jet_dau_pt    =  ptWRTJet((double)(*jetPt)[ijet], (double)(*jetEta)[ijet] +eta_smear    , (double)(*jetPhi)[ijet] +phi_smear      , (double)(*dau_pt)[ijet][A_trk], (double)(*dau_eta)[ijet][A_trk], (double)(*dau_phi)[ijet][A_trk]);
                    double jet_dau_eta   = etaWRTJet((double)(*jetPt)[ijet], (double)(*jetEta)[ijet] +eta_smear    , (double)(*jetPhi)[ijet] +phi_smear      , (double)(*dau_pt)[ijet][A_trk], (double)(*dau_eta)[ijet][A_trk], (double)(*dau_phi)[ijet][A_trk]);
                    double jet_dau_phi   = phiWRTJet((double)(*jetPt)[ijet], (double)(*jetEta)[ijet] +eta_smear    , (double)(*jetPhi)[ijet] +phi_smear      , (double)(*dau_pt)[ijet][A_trk], (double)(*dau_eta)[ijet][A_trk], (double)(*dau_phi)[ijet][A_trk]);

                    h_jet_jT->Fill(jet_dau_pt);
                    h_jet_etastar->Fill(jet_dau_eta);

                    double jet_dau_theta = 2*ATan(Exp(-(jet_dau_eta)));
                    if(jet_dau_eta > track_eta_lim) continue;
                    //getting et pt efficiency for A track in beam frame
                    double Atrk_weight = 1.0;//(hReco2D[thisEffTable]->GetBinContent(hReco2D[thisEffTable]->FindBin( (*dau_pt)[ijet][A_trk] , (*dau_eta)[ijet][A_trk] )));
                    hAvg_Atrk_Weight->Fill(1.0/(Atrk_weight));

                    /*for(int i = 0; i < trackbin; i++){
                        for(int j = 0; j < ptbin; j++){
                            if(tkBool[i] + A_ptBool[A_trk][j] == 2){
                                int k_PU=0;
                                hEPDrawCor[i][j][k_PU]->Fill(jet_dau_eta, jet_dau_phi, 1.0*jet_HLT_weight/(Atrk_weight * NtrigCorrected[i][j] ));
                            }
                        }
                    }*/

                    //*****************************************START FILL The hEPDraw***********************************
                    for(int ibeta=0;ibeta<beta_bin;ibeta++){
                        for(int iclass=0;iclass<SDclassifer_bin;iclass++){
                            for(int iNch=0;iNch<trackbin;iNch++){
                                for(int ipt=0;ipt<ptbin;ipt++){
                                    if( (tkBool[iNch] + A_ptBool[ipt][A_trk] + sdBool[ibeta][iclass])==3){
                                        int i_PU=0;
                                        hEPDrawCor[ibeta][iclass][iNch][ipt][i_PU]->Fill(jet_dau_eta, jet_dau_phi, 1.0*jet_HLT_weight/(Atrk_weight * NtrigCorrected[iNch][ipt] ));
                                    }
                                }
                            }
                        }
                    }
                    //*****************************************END  FILL The hEPDraw***********************************

                    if(A_trk == NNtrk - 1) continue;
                    //******************************************START T_trk**************************************************
                    for(long int T_trk=A_trk+1; T_trk< NNtrk; T_trk++ ){

                        if((*dau_chg)[ijet][T_trk] == 0) continue;
                        if((*dau_pt)[ijet][T_trk] < 0.3) continue;
                        if(fabs((*dau_eta)[ijet][T_trk]) > 2.4) continue;

                        double T_dauptnow = (*dau_pt)[ijet][T_trk];
                        double T_pterr_R_pt =0;
                        double T_dcaXY = 0;
                        double T_dcaZ = 0;
                        T_pterr_R_pt= ( (*dau_ptError)[ijet][T_trk] ) / ( (*dau_pt)[ijet][T_trk] );
                        T_dcaXY = (*dau_XYDCAsig)[ijet][T_trk];
                        T_dcaZ = (*dau_ZDCAsig)[ijet][T_trk];
                        if(T_pterr_R_pt   > 0.1   && T_dauptnow > 0.5) continue;
                        if(fabs(T_dcaXY)  > 3     && T_dauptnow > 0.5) continue;
                        if(fabs(T_dcaZ)   > 3     && T_dauptnow > 0.5) continue;

                        double T_jet_dau_pt    =  ptWRTJet( (double)(*jetPt)[ijet], (double)(*jetEta)[ijet] +eta_smear    , (double)(*jetPhi)[ijet] +phi_smear      , (double)(*dau_pt)[ijet][T_trk], (double)(*dau_eta)[ijet][T_trk], (double)(*dau_phi)[ijet][T_trk]);
                        double T_jet_dau_eta   = etaWRTJet( (double)(*jetPt)[ijet], (double)(*jetEta)[ijet] +eta_smear    , (double)(*jetPhi)[ijet] +phi_smear      , (double)(*dau_pt)[ijet][T_trk], (double)(*dau_eta)[ijet][T_trk], (double)(*dau_phi)[ijet][T_trk]);
                        double T_jet_dau_phi   = phiWRTJet( (double)(*jetPt)[ijet], (double)(*jetEta)[ijet] +eta_smear    , (double)(*jetPhi)[ijet] +phi_smear      , (double)(*dau_pt)[ijet][T_trk], (double)(*dau_eta)[ijet][T_trk], (double)(*dau_phi)[ijet][T_trk]);


                        if(T_jet_dau_eta > track_eta_lim) continue;

                        double deltaEta = (jet_dau_eta - T_jet_dau_eta);
                        double deltaPhi = (TMath::ACos(TMath::Cos(jet_dau_phi - T_jet_dau_phi)));
                        //NEW THING                        
                        //NEW THING                        
                        //NEW THING                        
                        //getting et pt efficiency for T track in beam frame
                        double Ttrk_weight =1.0; //(hReco2D[thisEffTable]->GetBinContent(hReco2D[thisEffTable]->FindBin( (*dau_pt)[ijet][T_trk] , (*dau_eta)[ijet][T_trk] )));
                        
                        

                       
                        //**************************************************START FILL Signal**************************************************
                        for(int ibeta=0;ibeta<beta_bin;ibeta++){
                            for(int iclass=0;iclass<SDclassifer_bin;iclass++){
                                for(int iNch=0;iNch<trackbin;iNch++){
                                    for(int ipt=0;ipt<ptbin;ipt++){
                                        if(tkBool[iNch]+A_ptBool[A_trk][ipt]+A_ptBool[T_trk][ipt]+sdBool[ibeta][iclass]==4 ){
                                            //cout << "1.0/(Atrk_weight*Ttrk_weight*NtrigCorrected[i][j]) " << 1.0/(Atrk_weight * Ttrk_weight * NtrigCorrected[i][j] ) << " and components: " << Atrk_weight << " " << Ttrk_weight << " " << NtrigCorrected[i][j] << endl;

                                            hPairs[ibeta][iclass]->Fill(iNch,ipt);
                                            if(Atrk_weight == 0){/* cout  << "Atrk_weight = 0!" << endl;
                                            */
                                                continue;}
                                            if(Ttrk_weight == 0){/* cout << "Ttrk_weight = 0!" << endl;
                                            */
                                                continue;}
                                            if(NtrigCorrected[iNch][ipt] == 0){ cout << "NtrigCorrected[i][j] = 0!" << endl; continue;}

                                            hTotalWeight->Fill(1.0 /(Atrk_weight * Ttrk_weight * NtrigCorrected[iNch][ipt] ));

                                            int i_PU=0;
                                            hSignalShiftedCor[ibeta][iclass][iNch][ipt][i_PU]->Fill(deltaEta, deltaPhi,                  1.0*jet_HLT_weight/(Atrk_weight * Ttrk_weight * NtrigCorrected[iNch][ipt] ));
                                            hSignalShiftedCor[ibeta][iclass][iNch][ipt][i_PU]->Fill(-deltaEta, deltaPhi,                 1.0*jet_HLT_weight/(Atrk_weight * Ttrk_weight * NtrigCorrected[iNch][ipt] ));
                                            hSignalShiftedCor[ibeta][iclass][iNch][ipt][i_PU]->Fill(deltaEta, -deltaPhi,                 1.0*jet_HLT_weight/(Atrk_weight * Ttrk_weight * NtrigCorrected[iNch][ipt] ));
                                            hSignalShiftedCor[ibeta][iclass][iNch][ipt][i_PU]->Fill(-deltaEta, -deltaPhi,                1.0*jet_HLT_weight/(Atrk_weight * Ttrk_weight * NtrigCorrected[iNch][ipt] ));
                                            hSignalShiftedCor[ibeta][iclass][iNch][ipt][i_PU]->Fill( deltaEta,2*TMath::Pi() - deltaPhi,  1.0*jet_HLT_weight/(Atrk_weight * Ttrk_weight * NtrigCorrected[iNch][ipt] ));
                                            hSignalShiftedCor[ibeta][iclass][iNch][ipt][i_PU]->Fill(-deltaEta,2*TMath::Pi() - deltaPhi,  1.0*jet_HLT_weight/(Atrk_weight * Ttrk_weight * NtrigCorrected[iNch][ipt] ));


                                        }
                                    }
                                }
                            }    
                        }    
                        //*************************************************END FILL Signal*********************************************************
                    }//****************************************END T_trk LOOP***********************************************
                }//*****************************************END A_trk LOOP*************************************************
            }//*********************************************END JET LOOP***********************************************

        }//***************************************************END EVENT LOOP*****************************************
        
        fFile->Close();
    }//***********************************************END FILE LOOP*****************************************
    //f in filelist end 

    
    /*for(int wtrk = 1; wtrk < trackbin+1; wtrk++){
        std::cout << wtrk << "/" << trackbin << std::endl;
        for(int wppt = 1; wppt < ptbin+1; wppt++){ 
            std::cout << wppt << "/" << ptbin << std::endl;
            for(int wpPU = 1; wpPU < PUbin+1; wpPU++){
                std::cout << wpPU << "/" << PUbin << std::endl;
                //Nent is the number of pairs in the signal which we will try to 10x
                //Xent is the number of pseudoparticles requried such that when we build the pairs nCp = Xent CHOOSE 2 will give 
                //us 10 times as many pairs as we have in the signal histogrm.

                //long int NENT =  hPairsPU[wpPU-1]->GetBinContent(wtrk, wppt);
               
            }
        }
    }*/


    for(int ibeta=0;ibeta<beta_bin;ibeta++){
        for(int iclass=0;iclass<SDclassifer_bin;iclass++){
            for(int iNch=0;iNch<trackbin;iNch++){
                for(int ipt=0;ipt<ptbin;ipt++){
                    int i_PU=0;
                    GenBackGround((int)hPairs[ibeta][iclass]->GetBinContent(iNch+1,ipt+1),hEPDrawCor[ibeta][iclass][iNch][ipt][i_PU],hBckrndShiftedCor[ibeta][iclass][iNch][ipt][i_PU]);
                }
            }
        }
    }

    //string subList = fList.substr(fList.size() - 3);
    //TFile* fS_tempA = new TFile(Form("/eos/cms/store/group/phys_heavyions/xiaoyul/Run3_2022_root_out/allNch/jetveto_CD/job_D_%s.root",subList.c_str()), "recreate");
    TFile *fout=new TFile("test.root","RECREATE");
    for(int wtrk = 1; wtrk<trackbin+1; wtrk++){
        //hBinDist_cor[wtrk-1]->Write();    
        //hBinDist_unc[wtrk-1]->Write();   
        for(int ibeta=1; ibeta<beta_bin+1;ibeta++){
            for(int iclass=1;iclass<SDclassifer_bin+1;iclass++){
                hBinDist_cor[ibeta-1][iclass-1][wtrk-1]->Write();    
                hBinDist_unc[ibeta-1][iclass-1][wtrk-1]->Write();  
                 for(int wppt = 1; wppt<ptbin+1; wppt++){
                     for(int wpPU = 1; wpPU<PUbin+1; wpPU++){
                        hBckrndShiftedCor[ibeta-1][iclass-1][wtrk-1][wppt-1][wpPU-1]->Write();   
                        hSignalShiftedCor[ibeta-1][iclass-1][wtrk-1][wppt-1][wpPU-1]->Write();   
                               hEPDrawCor[ibeta-1][iclass-1][wtrk-1][wppt-1][wpPU-1]->Write();   
                    }
                }
            }
        }
   }

   


    //Write SoftDrop Subtructure

    for(int ibeta=0;ibeta<beta_bin;ibeta++){
        h2Jet_Pass[ibeta]->Write();
        hNch_Zg[ibeta]->Write();                                                                                              
        hNch_Rg[ibeta]->Write();                                                                                                
        hNch_ZgTgB[ibeta]->Write();                                                                                            
        hZg_Rg[ibeta]->Write();

        hZg_dis[ibeta]->Write();
        hRg_dis[ibeta]->Write();
        hZgTgB_dis[ibeta]->Write();
        hSDJetMass_dis[ibeta]->Write();

        for(int iNch=0;iNch<trackbin;iNch++){
            hZgTgB_dis_Nch[ibeta][iNch]->Write();                                                                                       
            hZg_dis_Nch[ibeta][iNch]->Write();                                                                                     
            hRg_dis_Nch[ibeta][iNch]->Write();
            hSDJetMass_dis_Nch[ibeta][iNch]->Write();    
        }
    }


    //hPairs->Write();
    for(int ibeta=0;ibeta<beta_bin;ibeta++){
        for(int iclass=0;iclass<SDclassifer_bin;iclass++){
            hPairs[ibeta][iclass]->Write();
        }
    }
    hCalSDJetMassNch->Write();
    hStoredSDJetMassNch->Write();
    hNtrig->Write();
    hNtrigCorrected->Write();
    hTotalWeight->Write();
    hAvg_Atrk_Weight->Write();
    hAvg_NtrigCorrected_Bump->Write();
    hJetPt->Write();
    hJetPt_wo_ptcut->Write(); 
    hEvent_Pass->Write();
    hJet_Pass->Write();
    hJet_Pass550->Write();
    hJet_Pass550_hltCor->Write();
    hBinDist_cor_single->Write();
    hBinDist_unc_single->Write();
    
    h_jet_jT->Write();
    h_jet_etastar->Write();
    h_lab_JetMult_pT->Write();
    h_lab_JetMult_phi->Write();
    h_lab_JetMult_eta->Write();

    hjet_eta_phi_before_veto->Write();
    hjet_eta_phi_after_veto->Write();

    //fS_tempA->Close();
    fout->Close();
}


//Code enters execution here
int main(int argc, const char* argv[])
{   
    if(argc != 4)
    {
        std::cout << "Usage: Z_mumu_Channel <fileList> <jobNumber> <nJobs>" << std::endl;
        return 1;
    }


    //read input parameters
    std::string fList = argv[1];
    std::string buffer;
    std::vector<std::string> listOfFiles;
    std::ifstream inFile(fList.data());

    int job = (int)std::atoi(argv[2]);
    int nJobs = (int)std::atoi(argv[3]);


    //read the file list and spit it into a vector of strings based on how the parallelization is to be done
    //each vector is a separate subset of the fileList based on the job number
    if(!inFile.is_open())
    {
        std::cout << "Error opening jet file. Exiting." <<std::endl;
        return 1;
    }
    else
    {
        int line = 0;
        while(true)
        {
            inFile >> buffer;
            if(inFile.eof()) break;
            if( line%nJobs == job) listOfFiles.push_back(buffer);
            line++;
        }
    }

/*
    int job=1;

    std::string fList="MATCH_68.root";
    std::vector<std::string> listOfFiles;
    listOfFiles.push_back(fList);
*/
    //create the MyClass Object
    MyClass m = MyClass(listOfFiles);
    m.InitCutmode(3);
    m.Loop(job, fList);

    return 0;
}

