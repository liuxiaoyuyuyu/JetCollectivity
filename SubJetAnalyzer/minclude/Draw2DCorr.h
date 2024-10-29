#include "1d2d_constants.h"
#include "SubJetAnalyzer.h"
 
#ifndef Draw2DCorr_C
#define Draw2DCorr_C

void SubJetAnalyzer::Draw2DCorrelation(){
    
    std::cout<<"Start Draw2DCorrelation"<<std::endl;
    const int Classifier_bin=Classifier_lo.size(); 
    TH2D* hSignal[trackbin][Classifier_bin][ptbin]; 
    TH2D* hBkg[trackbin][Classifier_bin][ptbin];
    TH2D* h2DCorr[trackbin][Classifier_bin][ptbin];

    TH2D* hSignal_gg[trackbin][Classifier_bin][ptbin]; 
    TH2D* hBkg_gg[trackbin][Classifier_bin][ptbin];
    TH2D* h2DCorr_gg[trackbin][Classifier_bin][ptbin];


    TH2D* hSignal_ngng[trackbin][Classifier_bin][ptbin]; 
    TH2D* hBkg_ngng[trackbin][Classifier_bin][ptbin];
    TH2D* h2DCorr_ngng[trackbin][Classifier_bin][ptbin];

    //signal file 
    std::string fin_dir=store_dir+Form("/hpp_beta_%d.root",(int)beta_SD[beta_index]);
    TFile* fin=new TFile(fin_dir.c_str(),"READ");
    
    for(int i=0;i<trackbin;i++){
        for(int j=0;j<Classifier_bin;j++){
            for(int k=0;k<ptbin;k++){
                hSignal[i][j][k]=(TH2D *)fin->Get(Form("hSigS_Cor_%d_to_%d_%d_to_%d_%d_to_%d",trackbinbounds[i],trackbinboundsUpper[i] ,(int)(10*ptbinbounds_lo[k]),(int)(10*ptbinbounds_hi[k]),(int)(100*Classifier_lo[j]),(int)(100*Classifier_hi[j]) ) ); 
                hBkg[i][j][k]=(TH2D *)fin->Get( Form("hBckS_Cor_%d_to_%d_%d_to_%d_%d_to_%d",trackbinbounds[i],trackbinboundsUpper[i] ,(int)(10*ptbinbounds_lo[k]),(int)(10*ptbinbounds_hi[k]),(int)(100*Classifier_lo[j]),(int)(100*Classifier_hi[j]) ) );

                hSignal_ngng[i][j][k]=(TH2D *)fin->Get( Form("hSigS_ngng_Cor_%d_to_%d_%d_to_%d_%d_to_%d",trackbinbounds[i],trackbinboundsUpper[i] ,(int)(10*ptbinbounds_lo[k]),(int)(10*ptbinbounds_hi[k]),(int)(100*Classifier_lo[j]),(int)(100*Classifier_hi[j]) ) );
                hBkg_ngng[i][j][k]=(TH2D *)fin->Get( Form("hBckS_ngng_Cor_%d_to_%d_%d_to_%d_%d_to_%d",trackbinbounds[i],trackbinboundsUpper[i] ,(int)(10*ptbinbounds_lo[k]),(int)(10*ptbinbounds_hi[k]),(int)(100*Classifier_lo[j]),(int)(100*Classifier_hi[j]) ) );

                hSignal_gg[i][j][k]=(TH2D *)fin->Get( Form("hSigS_gg_Cor_%d_to_%d_%d_to_%d_%d_to_%d",trackbinbounds[i],trackbinboundsUpper[i] ,(int)(10*ptbinbounds_lo[k]),(int)(10*ptbinbounds_hi[k]),(int)(100*Classifier_lo[j]),(int)(100*Classifier_hi[j]) ) );
                hBkg_gg[i][j][k]=(TH2D *)fin->Get( Form("hBckS_gg_Cor_%d_to_%d_%d_to_%d_%d_to_%d",trackbinbounds[i],trackbinboundsUpper[i] ,(int)(10*ptbinbounds_lo[k]),(int)(10*ptbinbounds_hi[k]),(int)(100*Classifier_lo[j]),(int)(100*Classifier_hi[j]) ) );

            }
        }
    }
    std::string fout_dir=store_dir+Form("/Corr2D_beta_%d.root",(int)beta_SD[beta_index]);
    TFile *fout=new TFile(fout_dir.c_str(),"RECREATE");
    
    TH2D* hJetPass = (TH2D*)fin->Get("hJet_Pass");
    for(int i=0;i<trackbin;i++){
        for(int j=0;j<Classifier_bin;j++){
            for(int k=0;k<ptbin;k++){
                h2DCorr[i][j][k]=(TH2D*)hSignal[i][j][k]->Clone(Form("h2DCorr_Cor_%d_to_%d_%d_to_%d_%d_to_%d",trackbinbounds[i],trackbinboundsUpper[i] ,(int)(10*ptbinbounds_lo[k]),(int)(10*ptbinbounds_hi[k]),(int)(100*Classifier_lo[j]),(int)(100*Classifier_hi[j]) ));
                //h2DCorr[i][j][k]->Reset();
                double x_bw=hSignal[i][j][k]->GetXaxis()->GetBinWidth(1);
                double y_bw=hSignal[i][j][k]->GetYaxis()->GetBinWidth(1);
                double tmp_njet=hJetPass->GetBinContent(i+1,j+1);
                double tmp_b00=hBkg[i][j][k]->GetBinContent(hBkg[i][j][k]->FindBin(0,0));
                h2DCorr[i][j][k]->Divide(hBkg[i][j][k]);
                h2DCorr[i][j][k]->Scale( tmp_b00/ (tmp_njet * x_bw * y_bw) );
                h2DCorr[i][j][k]->Write();
            }
        }
    }

    for(int i=0;i<trackbin;i++){
        for(int j=0;j<Classifier_bin;j++){
            for(int k=0;k<ptbin;k++){
                h2DCorr_gg[i][j][k]=(TH2D*)hSignal_gg[i][j][k]->Clone(Form("h2DCorr_gg_Cor_%d_to_%d_%d_to_%d_%d_to_%d",trackbinbounds[i],trackbinboundsUpper[i] ,(int)(10*ptbinbounds_lo[k]),(int)(10*ptbinbounds_hi[k]),(int)(100*Classifier_lo[j]),(int)(100*Classifier_hi[j]) ));
                //h2DCorr[i][j][k]->Reset();
                double x_bw=hSignal_gg[i][j][k]->GetXaxis()->GetBinWidth(1);
                double y_bw=hSignal_gg[i][j][k]->GetYaxis()->GetBinWidth(1);
                double tmp_njet=hJetPass->GetBinContent(i+1,j+1);
                double tmp_b00=hBkg_gg[i][j][k]->GetBinContent(hBkg_gg[i][j][k]->FindBin(0,0));
                h2DCorr_gg[i][j][k]->Divide(hBkg_gg[i][j][k]);
                h2DCorr_gg[i][j][k]->Scale( tmp_b00/ (tmp_njet * x_bw * y_bw) );
                h2DCorr_gg[i][j][k]->Write();
            }
        }
    }

    for(int i=0;i<trackbin;i++){
        for(int j=0;j<Classifier_bin;j++){
            for(int k=0;k<ptbin;k++){
                h2DCorr_ngng[i][j][k]=(TH2D*)hSignal_ngng[i][j][k]->Clone(Form("h2DCorr_ngng_Cor_%d_to_%d_%d_to_%d_%d_to_%d",trackbinbounds[i],trackbinboundsUpper[i] ,(int)(10*ptbinbounds_lo[k]),(int)(10*ptbinbounds_hi[k]),(int)(100*Classifier_lo[j]),(int)(100*Classifier_hi[j]) ));
                //h2DCorr[i][j][k]->Reset();
                double x_bw=hSignal_ngng[i][j][k]->GetXaxis()->GetBinWidth(1);
                double y_bw=hSignal_ngng[i][j][k]->GetYaxis()->GetBinWidth(1);
                double tmp_njet=hJetPass->GetBinContent(i+1,j+1);
                double tmp_b00=hBkg_ngng[i][j][k]->GetBinContent(hBkg_ngng[i][j][k]->FindBin(0,0));
                h2DCorr_ngng[i][j][k]->Divide(hBkg_ngng[i][j][k]);
                h2DCorr_ngng[i][j][k]->Scale( tmp_b00/ (tmp_njet * x_bw * y_bw) );
                h2DCorr_ngng[i][j][k]->Write();
            }
        }
    }

    fin->Close();
    fout->Close();
    delete fin;
    delete fout;

}

#endif