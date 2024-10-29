#include "1d2d_constants.h"
#include "SubJetAnalyzer.h"

#ifndef DrawFlow_C
#define DrawFlow_C

void SubJetAnalyzer::DrawFlow(){
    std::cout<<std::endl;
    std::cout<<"Starting DrawFolow"<<std::endl;

    const int Classifier_bin=Classifier_lo.size(); 

    TH2D* hSignal[trackbin][Classifier_bin][ptbin];
    TH2D* hBkg[trackbin][Classifier_bin][ptbin];
    TH1D* h1DFlow[trackbin][Classifier_bin][ptbin];

    TH2D* hSignal_gg[trackbin][Classifier_bin][ptbin];
    TH2D* hBkg_gg[trackbin][Classifier_bin][ptbin];
    TH1D* h1DFlow_gg[trackbin][Classifier_bin][ptbin];

    TH2D* hSignal_ngng[trackbin][Classifier_bin][ptbin];
    TH2D* hBkg_ngng[trackbin][Classifier_bin][ptbin];
    TH1D* h1DFlow_ngng[trackbin][Classifier_bin][ptbin];
    
    std::string fin_dir=store_dir+Form("/hpp_beta_%d.root",(int)beta_SD[beta_index]);
    
    TFile *fin=TFile::Open( fin_dir.c_str() ,"READ" );
    TH2D *hJetPass=(TH2D*)fin->Get("hJet_Pass");


    int YPlo=28;
    int YPhi=34;
    double v2[Classifier_bin][ptbin][trackbin]={0.};
    double v2Error[Classifier_bin][ptbin][trackbin]={0.};

    double v2_gg[Classifier_bin][ptbin][trackbin]={0.};
    double v2Error_gg[Classifier_bin][ptbin][trackbin]={0.};

    double v2_ngng[Classifier_bin][ptbin][trackbin]={0.};
    double v2Error_ngng[Classifier_bin][ptbin][trackbin]={0.};


    double trackbinmid[trackbin]={0.};
    double trackbinErr[trackbin]={0.};


    for(int i=0;i<trackbin;i++){
        trackbinmid[i]=(trackbinbounds[i]+trackbinboundsUpper[i])/2.;

    }
    
    std::string function = "[0]/(TMath::Pi()*2)*(1+2*([1]*TMath::Cos(x)+[2]*TMath::Cos(2*x)+[3]*TMath::Cos(3*x)+[4]*TMath::Cos(4*x)+[5]*TMath::Cos(5*x)))";


    for(int i=0; i<trackbin;i++){
        for(int j=0;j<Classifier_bin;j++){
            for(int k=0;k<ptbin;k++){
                

                hSignal[i][j][k]=(TH2D *)fin->Get( Form("hSigS_Cor_%d_to_%d_%d_to_%d_%d_to_%d",trackbinbounds[i],trackbinboundsUpper[i] ,(int)(10*ptbinbounds_lo[k]),(int)(10*ptbinbounds_hi[k]),(int)(100*Classifier_lo[j]),(int)(100*Classifier_hi[j]) ) );
                hBkg[i][j][k]=(TH2D *)fin->Get( Form("hBckS_Cor_%d_to_%d_%d_to_%d_%d_to_%d",trackbinbounds[i],trackbinboundsUpper[i] ,(int)(10*ptbinbounds_lo[k]),(int)(10*ptbinbounds_hi[k]),(int)(100*Classifier_lo[j]),(int)(100*Classifier_hi[j]) ) );

                double eta_bw=hSignal[i][j][k]->GetXaxis()->GetBinWidth(1); //x:delta eta; y:delta phi
                double phi_bw=hSignal[i][j][k]->GetYaxis()->GetBinWidth(1);
                hSignal[i][j][k]->Scale( 1.0/(hJetPass->GetBinContent(i+1,j+1)) );
                hSignal[i][j][k]->Scale(1./(phi_bw*eta_bw));

                TH1D *histfit1 = (TH1D*) hSignal[i][j][k]->ProjectionY("",YPlo,YPhi)->Clone();
                TH1D *histfit2 = (TH1D*) hBkg[i][j][k]->ProjectionY("",YPlo,YPhi)->Clone();
                histfit1->Divide(histfit2);
                histfit1->Scale(hBkg[i][j][k]->GetMaximum());
                h1DFlow[i][j][k]=(TH1D*)histfit1->Clone(); 
                TF1 func1("deltaPhi1", function.c_str(), -0.5*TMath::Pi(), 1.5*TMath::Pi());
                func1.SetParameter(0, histfit1->GetMaximum());
                func1.SetParameter(1, 0.1);
                func1.SetParameter(2, 0.1);
                func1.SetParameter(3, 0.1);
                func1.SetParameter(4, 0.1);
                func1.SetParameter(5, 0.1);
                h1DFlow[i][j][k]->Fit(&func1, "m E q");

                v2[j][k][i]=func1.GetParameter(2);
                v2Error[j][k][i]=(func1.GetParError(2))*TMath::Sqrt(2.);

            }
        }
    }

    for(int i=0; i<trackbin;i++){
        for(int j=0;j<Classifier_bin;j++){
            for(int k=0;k<ptbin;k++){
                

                hSignal_gg[i][j][k]=(TH2D *)fin->Get( Form("hSigS_gg_Cor_%d_to_%d_%d_to_%d_%d_to_%d",trackbinbounds[i],trackbinboundsUpper[i] ,(int)(10*ptbinbounds_lo[k]),(int)(10*ptbinbounds_hi[k]),(int)(100*Classifier_lo[j]),(int)(100*Classifier_hi[j]) ) );
                hBkg_gg[i][j][k]=(TH2D *)fin->Get( Form("hBckS_gg_Cor_%d_to_%d_%d_to_%d_%d_to_%d",trackbinbounds[i],trackbinboundsUpper[i] ,(int)(10*ptbinbounds_lo[k]),(int)(10*ptbinbounds_hi[k]),(int)(100*Classifier_lo[j]),(int)(100*Classifier_hi[j]) ) );

                double eta_bw=hSignal_gg[i][j][k]->GetXaxis()->GetBinWidth(1); //x:delta eta; y:delta phi
                double phi_bw=hSignal_gg[i][j][k]->GetYaxis()->GetBinWidth(1);
                hSignal_gg[i][j][k]->Scale( 1.0/(hJetPass->GetBinContent(i+1,j+1)) );
                hSignal_gg[i][j][k]->Scale(1./(phi_bw*eta_bw));

                TH1D *histfit1 = (TH1D*) hSignal_gg[i][j][k]->ProjectionY("",YPlo,YPhi)->Clone();
                TH1D *histfit2 = (TH1D*) hBkg_gg[i][j][k]->ProjectionY("",YPlo,YPhi)->Clone();
                histfit1->Divide(histfit2);
                histfit1->Scale(hBkg_gg[i][j][k]->GetMaximum());
                h1DFlow_gg[i][j][k]=(TH1D*)histfit1->Clone(); 
                TF1 func1("deltaPhi1", function.c_str(), -0.5*TMath::Pi(), 1.5*TMath::Pi());
                func1.SetParameter(0, histfit1->GetMaximum());
                func1.SetParameter(1, 0.1);
                func1.SetParameter(2, 0.1);
                func1.SetParameter(3, 0.1);
                func1.SetParameter(4, 0.1);
                func1.SetParameter(5, 0.1);
                h1DFlow_gg[i][j][k]->Fit(&func1, "m E q");

                v2_gg[j][k][i]=func1.GetParameter(2);
                v2Error_gg[j][k][i]=(func1.GetParError(2))*TMath::Sqrt(2.);

            }
        }
    }


    for(int i=0; i<trackbin;i++){
        for(int j=0;j<Classifier_bin;j++){
            for(int k=0;k<ptbin;k++){
                

                hSignal_ngng[i][j][k]=(TH2D *)fin->Get( Form("hSigS_ngng_Cor_%d_to_%d_%d_to_%d_%d_to_%d",trackbinbounds[i],trackbinboundsUpper[i] ,(int)(10*ptbinbounds_lo[k]),(int)(10*ptbinbounds_hi[k]),(int)(100*Classifier_lo[j]),(int)(100*Classifier_hi[j]) ) );
                hBkg_ngng[i][j][k]=(TH2D *)fin->Get( Form("hBckS_ngng_Cor_%d_to_%d_%d_to_%d_%d_to_%d",trackbinbounds[i],trackbinboundsUpper[i] ,(int)(10*ptbinbounds_lo[k]),(int)(10*ptbinbounds_hi[k]),(int)(100*Classifier_lo[j]),(int)(100*Classifier_hi[j]) ) );

                double eta_bw=hSignal_ngng[i][j][k]->GetXaxis()->GetBinWidth(1); //x:delta eta; y:delta phi
                double phi_bw=hSignal_ngng[i][j][k]->GetYaxis()->GetBinWidth(1);
                hSignal_ngng[i][j][k]->Scale( 1.0/(hJetPass->GetBinContent(i+1,j+1)) );
                hSignal_ngng[i][j][k]->Scale(1./(phi_bw*eta_bw));

                TH1D *histfit1 = (TH1D*) hSignal_ngng[i][j][k]->ProjectionY("",YPlo,YPhi)->Clone();
                TH1D *histfit2 = (TH1D*) hBkg_ngng[i][j][k]->ProjectionY("",YPlo,YPhi)->Clone();
                histfit1->Divide(histfit2);
                histfit1->Scale(hBkg_ngng[i][j][k]->GetMaximum());
                h1DFlow_ngng[i][j][k]=(TH1D*)histfit1->Clone(); 
                TF1 func1("deltaPhi1", function.c_str(), -0.5*TMath::Pi(), 1.5*TMath::Pi());
                func1.SetParameter(0, histfit1->GetMaximum());
                func1.SetParameter(1, 0.1);
                func1.SetParameter(2, 0.1);
                func1.SetParameter(3, 0.1);
                func1.SetParameter(4, 0.1);
                func1.SetParameter(5, 0.1);
                h1DFlow_ngng[i][j][k]->Fit(&func1, "m E q");

                v2_ngng[j][k][i]=func1.GetParameter(2);
                v2Error_ngng[j][k][i]=(func1.GetParError(2))*TMath::Sqrt(2.);

            }
        }
    }





    std::string fout_dir=store_dir+Form("/flow_beta_%d.root",(int)beta_SD[beta_index]);
    TFile *fout=new TFile( fout_dir.c_str() , "RECREATE");
    TGraphErrors *gr[Classifier_bin][ptbin];
    TGraphErrors *gr_gg[Classifier_bin][ptbin];
    TGraphErrors *gr_ngng[Classifier_bin][ptbin];
    for(int j=0;j<Classifier_bin;j++){
        for(int k=0;k<ptbin;k++){
            gr[j][k]=new TGraphErrors(trackbin,trackbinmid,v2[j][k],trackbinErr,v2Error[j][k]);
            gr[j][k]->Write(Form("v2_Nch_Class_%d_%d_pt_%d_%d",(int)(Classifier_lo[j]*100),(int)(Classifier_hi[j]*100),(int)(ptbinbounds_lo[k]*10),(int)(ptbinbounds_hi[k]*10) ));
        }
    }


     for(int j=0;j<Classifier_bin;j++){
        for(int k=0;k<ptbin;k++){
            gr_gg[j][k]=new TGraphErrors(trackbin,trackbinmid,v2_gg[j][k],trackbinErr,v2Error_gg[j][k]);
            gr_gg[j][k]->Write(Form("v2_Nch_gg_Class_%d_%d_pt_%d_%d",(int)(Classifier_lo[j]*100),(int)(Classifier_hi[j]*100),(int)(ptbinbounds_lo[k]*10),(int)(ptbinbounds_hi[k]*10) ));
        }
    }


     for(int j=0;j<Classifier_bin;j++){
        for(int k=0;k<ptbin;k++){
            gr_ngng[j][k]=new TGraphErrors(trackbin,trackbinmid,v2_ngng[j][k],trackbinErr,v2Error_ngng[j][k]);
            gr_ngng[j][k]->Write(Form("v2_Nch_ngng_Class_%d_%d_pt_%d_%d",(int)(Classifier_lo[j]*100),(int)(Classifier_hi[j]*100),(int)(ptbinbounds_lo[k]*10),(int)(ptbinbounds_hi[k]*10) ));
        }
    }





    for(int i=0;i<trackbin;i++){
        for(int j=0;j<Classifier_bin;j++){
            for(int k=0;k<ptbin;k++){
                h1DFlow[i][j][k]->Write(Form("1d_corr_track_%d_%d_Class_%d_%d_pt_%d_%d",trackbinbounds[i],trackbinboundsUpper[i],(int)(Classifier_lo[j]*100),(int)(Classifier_hi[j]*100),(int)(ptbinbounds_lo[k]*10),(int)(ptbinbounds_hi[k]*10) ));
                h1DFlow_gg[i][j][k]->Write(Form("1d_corr_gg_track_%d_%d_Class_%d_%d_pt_%d_%d",trackbinbounds[i],trackbinboundsUpper[i],(int)(Classifier_lo[j]*100),(int)(Classifier_hi[j]*100),(int)(ptbinbounds_lo[k]*10),(int)(ptbinbounds_hi[k]*10) ));
                h1DFlow_ngng[i][j][k]->Write(Form("1d_corr_ngng_track_%d_%d_Class_%d_%d_pt_%d_%d",trackbinbounds[i],trackbinboundsUpper[i],(int)(Classifier_lo[j]*100),(int)(Classifier_hi[j]*100),(int)(ptbinbounds_lo[k]*10),(int)(ptbinbounds_hi[k]*10) ));
            }
        }
    }

    fin->Close();
    fout->Close();
    delete fin;
    delete fout;

}



#endif
