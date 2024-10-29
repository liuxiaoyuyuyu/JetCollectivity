

#include "1d2d_constants.h"
#include "coordinateTools.h" 






class SDSubstructureWorker{
public:
   TFile          *fFile;
   TTree          *fChain;
   long int            TotalEntry;
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
   //TBranch        *b_dau_DropFlag[beta_bin];


   std::vector<TBranch *> b_dau_Flag ;

   



   SDSubstructureWorker(std::vector<std::string>_SDfileList,std::string _store_dir);
   virtual ~SDSubstructureWorker();
   virtual void     LoadData(TTree *tree);
   virtual void     GetMyEntry(int entry);
   virtual void     SDpainter();




};



SDSubstructureWorker::SDSubstructureWorker(std::vector<std::string> _SDfileList,std::string _store_dir)
	:fChain(0),SDfileList(_SDfileList),store_dir(_store_dir) 
{
	
}




SDSubstructureWorker::~SDSubstructureWorker()
{

}



void SDSubstructureWorker::LoadData(TTree *tree){
   
   b_dau_Flag.clear();
   
   for(int ibeta=0;ibeta<beta_bin;ibeta++){
      b_dau_Flag.push_back(0);
   }
   

   jetPt=0;
   jetEta=0;
   jetPhi;
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
   TotalEntry=fChain->GetEntriesFast();
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("jetPt"                 ,&jetPt                 ,&b_jetPt                 );
   fChain->SetBranchAddress("jetEta"                ,&jetEta                 ,&b_jetEta                );
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

void SDSubstructureWorker::GetMyEntry(int entry){
	if (!fChain) return 0;
	fChain->GetEntry(entry);
}


void SDSubstructureWorker::SDpainter(){
	double R=0.8;
	int waring_Flag=0;
	//DEFINE Historgam here 
   
	//IN THE LAB FRAME, pt ,phi, eta of jet distribution
    TH2D* h_lab_JetMult_pT=new TH2D("lab_Jet_Nch_VS_pT","lab_Jet_Nch_VS_pT",300,100,6000,120,0.5,120.5);                     //
    TH2D* h_lab_JetMult_phi=new TH2D("lab_Jet_Nch_VS_phi","lab_Jet_Nch_VS_phi",30,-TMath::Pi(),TMath::Pi(),120,0.5,120.5);   //
    TH2D* h_lab_JetMult_eta=new TH2D("lab_Jet_Nch_VS_eta","lab_Jet_Nch_VS_eta",34,-1.7,1.7,120,0.5,120.5);                   //
    //IN THE JET FRAME ,pt etastar of jet daughter particles
    TH1D* h_jet_jT[trackbin];       																						 //
    TH1D* h_jet_etastar[trackbin];																							 //

    for(int i=0;i<trackbin;i++){
    	h_jet_jT[i]=new TH1D(Form("h_jet_jT_dis_Nch_%d_%d",trackbinbounds[i],trackbinboundsUpper[i]),Form("h_jet_jT_dis_Nch_%d_%d",trackbinbounds[i],trackbinboundsUpper[i]),200,0,20);
    	h_jet_etastar[i]=new TH1D(Form("h_jet_etastar_dis_Nch_%d_%d",trackbinbounds[i],trackbinboundsUpper[i]),Form("h_jet_etastar_dis_Nch_%d_%d",trackbinbounds[i],trackbinboundsUpper[i]),100,0,10);
    }

    TH1D* h_jet_groomed_jT[beta_bin][trackbin];                                                                               //
    TH1D* h_jet_non_groomed_jT[beta_bin][trackbin];   																		  //
    TH1D* h_jet_groomed_etastar[beta_bin][trackbin];																		  //
    TH1D* h_jet_non_groomed_etastar[beta_bin][trackbin];																	  //

    for(int i=0;i<beta_bin;i++){
    	for(int j=0;j<trackbin;j++){
    		h_jet_groomed_jT[i][j]   = new TH1D(Form("h_jet_groomed_jT_dis_Nch_%d_%d_beta_%d",trackbinbounds[j],trackbinboundsUpper[j],(int)beta_SD[i]),Form("h_jet_groomed_jT_dis_Nch_%d_%d_beta_%d",trackbinbounds[j],trackbinboundsUpper[j],(int)beta_SD[i]),200,0,20); 
    		h_jet_non_groomed_jT[i][j] = new TH1D(Form("h_jet_non_groomed_jT_dis_Nch_%d_%d_beta_%d",trackbinbounds[j],trackbinboundsUpper[j],(int)beta_SD[i]),Form("h_jet_non_groomed_jT_dis_Nch_%d_%d_beta_%d",trackbinbounds[j],trackbinboundsUpper[j],(int)beta_SD[i]),200,0,20); 
    		h_jet_groomed_etastar[i][j] =new TH1D(Form("h_jet_groomed_etastar_dis_Nch_%d_%d_beta_%d",trackbinbounds[j],trackbinboundsUpper[j],(int)beta_SD[i]),Form("h_jet_groomed_etastar_dis_Nch_%d_%d_beta_%d",trackbinbounds[j],trackbinboundsUpper[j],(int)beta_SD[i]),100,0,10); 
    		h_jet_non_groomed_etastar[i][j]=new TH1D(Form("h_jet_non_groomed_etastar_dis_Nch_%d_%d_beta_%d",trackbinbounds[j],trackbinboundsUpper[j],(int)beta_SD[i]),Form("h_non_jet_groomed_etastar_dis_Nch_%d_%d_beta_%d",trackbinbounds[j],trackbinboundsUpper[j],(int)beta_SD[i]),100,0,10); 


    	}
    }

    //TH2D* h_jet_jt_Nch =new TH2D("h_jet_jt_Nch","h_jet_jt_Nch",200,0,20,120,0,120 ) ;
    //TH2D* h_jet_etastar_Nch=new TH2D("h_jet_etastar_Nch","h_jet_etastar_Nch",100,0,10,120,0,120);

    //TH2D* h_jet_groomed_jt_Nch[beta_bin];
    //TH2D* h_jet_grommed_etastar_Nch[beta_bin];


    //TH2D* h_jet_non_groomed_jt_Nch[beta_bin];
    //TH2D* h_jet_non_grommed_etastar_Nch[beta_bin];

   TH2D *hNch_Zg[beta_bin];      																							//
	TH2D *hNch_Rg[beta_bin];																								//
	TH2D *hNch_ZgRg[beta_bin];																								//
	TH2D *hZg_Rg[beta_bin];																									//
	TH2D *hlogZg_logRg[beta_bin];																							//

	TH2D *hZg_Rg_Nch[beta_bin][trackbin];																					//
	TH2D *hlogZg_logRg_Nch[beta_bin][trackbin];																				//
	TH1D *h_Ngroomed[beta_bin][trackbin];																					//
	TH1D *h_ratio_groomed[beta_bin][trackbin];																				//
	//this actually should be Zg * (Rg/R)^{\beta}
	TH1D *hZgRg_dis[beta_bin][trackbin];																						//
	TH1D *hZg_dis[beta_bin][trackbin];																						//
	TH1D *hRg_dis[beta_bin][trackbin];																						//

   TH2D *hZg_Ngroomed[beta_bin];
   TH2D *hZg_Rgroomed[beta_bin];
   TH2D *hRg_Ngroomed[beta_bin];
   TH2D *hRg_Rgroomed[beta_bin];
   TH2D *hZgTg_Ngroomed[beta_bin];
   TH2D *hZgTg_Rgroomed[beta_bin];


   TH2D *hZg_Ngroomed_Nch[beta_bin][trackbin];
   TH2D *hZg_Rgroomed_Nch[beta_bin][trackbin];
   TH2D *hRg_Ngroomed_Nch[beta_bin][trackbin];
   TH2D *hRg_Rgroomed_Nch[beta_bin][trackbin];
   TH2D *hZgTg_Ngroomed_Nch[beta_bin][trackbin];
   TH2D *hZgTg_Rgroomed_Nch[beta_bin][trackbin];

   TH1D *hdelta_phi_star_j1j2[beta_bin][trackbin];// In jet frame , the delte phi* beween two subjets.

	
	for(int i=0;i<beta_bin;i++){
		//h_jet_groomed_jt_Nch[i]=new TH2D(Form("h_jet_groomed_jt_VS_Nch_beta_%d",(int)beta_SD[i] ),Form("h_jet_groomed_jt_VS_Nch_beta_%d",(int)beta_SD[i] ),200,0,20,120,0,120);
		//h_jet_grommed_etastar_Nch[i]=new TH2D(Form("h_jet_groomed_etastar_VS_Nch_beta_%d",(int)beta_SD[i] ),Form("h_jet_groomed_etastar_VS_Nch_beta_%d",(int)beta_SD[i] ),100,0,10,120,0,120);

		//h_jet_non_groomed_jt_Nch[i]=new TH2D(Form("h_jet_non_groomed_jt_VS_Nch_beta_%d",(int)beta_SD[i] ),Form("h_jet_groomed_jt_VS_Nch_beta_%d",(int)beta_SD[i] ),200,0,20,120,0,120);
		//h_jet_non_grommed_etastar_Nch[i]=new TH2D(Form("h_jet_non_groomed_etastar_VS_Nch_beta_%d",(int)beta_SD[i] ),Form("h_jet_non_groomed_etastar_VS_Nch_beta_%d",(int)beta_SD[i] ),100,0,10,120,0,120);

		


		hNch_ZgRg[i]=new TH2D(Form("Nch_VS_ZgThetagBeta_beta_%d",(int)beta_SD[i]),Form("Nch_VS_ZgThetagBeta_beta_%d",(int)beta_SD[i]), 120 , 0 , 120 , 105 , -0.05 , 1.0 );
		hNch_Rg[i]=new TH2D(Form("Nch_VS_R_g_beta_%d",(int)beta_SD[i]),Form("Nch_VS_R_g_beta_%d",(int)beta_SD[i]), 120 , 0 , 120 , 50 , -0.05 , 1   );
		hNch_Zg[i]=new TH2D(Form("Nch_VS_Z_g_beta_%d",(int)beta_SD[i]),Form("Nch_VS_Z_g_beta_%d",(int)beta_SD[i]), 120 , 0 , 120 , 50 , -0.05 , 0.6 );
		hZg_Rg[i]=new TH2D(Form("Z_g_VS_R_g_beta_%d",(int)beta_SD[i]),Form("Z_g_VS_R_g_beta_%d", (int)beta_SD[i]), 50 , -0.05 , 0.6 , 50 , -0.05 , 1);
		hlogZg_logRg[i]=new TH2D(Form("logZ_g_VS_logR_g_beta_%d",(int)beta_SD[i]),Form("logZ_g_VS_logR_g_beta_%d", (int)beta_SD[i]), 50 , -3 , -0.3 , 50 , -3 , 0);


      hZg_Ngroomed[i]=new TH2D(Form("Zg_Vs_NumberGroomed_beta_%d",(int)beta_SD[i]),Form("Zg_Vs_NumberGroomed_beta_%d",(int)beta_SD[i]),65,-0.05,0.6,120,0,120);
      hZg_Rgroomed[i]=new TH2D(Form("Zg_Vs_RatioGroomed_beta_%d",(int)beta_SD[i]),Form("Zg_Vs_RatioGroomed_beta_%d",(int)beta_SD[i]),65,-0.05,0.6,110,-0.05,1.05);
      hRg_Ngroomed[i]=new TH2D(Form("Rg_Vs_NumberGroomed_beta_%d",(int)beta_SD[i]),Form("Rg_Vs_NumberGroomed_beta_%d",(int)beta_SD[i]),110,-0.05,1.05,120,0,120);
      hRg_Rgroomed[i]=new TH2D(Form("Rg_Vs_RatioGroomed_beta_%d",(int)beta_SD[i]),Form("Rg_Vs_RatioGroomed_beta_%d",(int)beta_SD[i]),110,-0.05,1.05,110,-0.05,1.05);

      hZgTg_Ngroomed[i]=new TH2D(Form("ZgTgBeta_Vs_NumberGroomed_beta_%d",(int)beta_SD[i]),Form("ZgTgBeta_Vs_NumberGroomed_beta_%d",(int)beta_SD[i]),110,-0.05,1.05,120,0,120);
      hZgTg_Rgroomed[i]=new TH2D(Form("ZgTgBeta_Vs_RatioGroomed_beta_%d",(int)beta_SD[i]),Form("ZgTgBeta_Vs_RatioGroomed_beta_%d",(int)beta_SD[i]),110,-0.05,1.05,110,-0.05,1.05);



		for(int j=0;j<trackbin;j++){
			hZg_Rg_Nch[i][j]=new TH2D( Form("Zg_VS_Rg_Nch_%d_%d_beta_%d",trackbinbounds[j],trackbinboundsUpper[j],(int)beta_SD[i]) ,Form("Zg_VS_Rg_Nch_%d_%d_beta_%d",trackbinbounds[j],trackbinboundsUpper[j],(int)beta_SD[i]) ,50 , -0.05 , 0.6 , 50 , -0.05 , 1 );
			hlogZg_logRg_Nch[i][j]=new TH2D( Form("logZg_VS_logRg_Nch_%d_%d_beta_%d",trackbinbounds[j],trackbinboundsUpper[j],(int)beta_SD[i]) ,Form("logZg_VS_logRg_Nch_%d_%d_beta_%d",trackbinbounds[j],trackbinboundsUpper[j],(int)beta_SD[i]) , 50 , -3 , -0.3 , 50 , -3 , 0 );

			h_Ngroomed[i][j]=new TH1D( Form("Groomed_number_dis_Nch_%d_%d_beta_%d",trackbinbounds[j],trackbinboundsUpper[j],(int)beta_SD[i]), Form("Groomed_number_dis_Nch_%d_%d_beta_%d",trackbinbounds[j],trackbinboundsUpper[j],(int)beta_SD[i]) , 120,0,120 );
			h_ratio_groomed[i][j]=new TH1D( Form("Groomed_ratio_dis_Nch_%d_%d_beta_%d",trackbinbounds[j],trackbinboundsUpper[j],(int)beta_SD[i]), Form("Groomed_ratio_dis_Nch_%d_%d_beta_%d",trackbinbounds[j],trackbinboundsUpper[j],(int)beta_SD[i]) , 60,0,1.2);
			hZgRg_dis[i][j]=new TH1D(Form("ZgThetagBeta_dis_Nch_%d_%d_beta_%d",trackbinbounds[j],trackbinboundsUpper[j],(int)beta_SD[i]),Form("ZgThetagBeta_dis_Nch_%d_%d_beta_%d",trackbinbounds[j],trackbinboundsUpper[j],(int)beta_SD[i]),100,0,1.0);
			
			hZg_dis[i][j]=new TH1D(Form("ZgDis_Nch_%d_%d_beta_%d",trackbinbounds[j],trackbinboundsUpper[j],(int)beta_SD[i]),Form("ZgDis_Nch_%d_%d_beta_%d",trackbinbounds[j],trackbinboundsUpper[j],(int)beta_SD[i]),50,-0.05,0.55);
			hRg_dis[i][j]=new TH1D(Form("RgDis_Nch_%d_%d_beta_%d",trackbinbounds[j],trackbinboundsUpper[j],(int)beta_SD[i]),Form("RgDis_Nch_%d_%d_beta_%d",trackbinbounds[j],trackbinboundsUpper[j],(int)beta_SD[i]),50,-0.05,0.8);
			

         hZg_Ngroomed_Nch[i][j]=new TH2D(Form("Zg_Vs_NumberGroomed_beta_%d_Nch_%d_%d",(int)beta_SD[i],trackbinbounds[j],trackbinboundsUpper[j] ),Form("Zg_Vs_NumberGroomed_beta_%d_Nch_%d_%d",(int)beta_SD[i],trackbinbounds[j],trackbinboundsUpper[j]),65,-0.05,0.6,120,0,120);
         hZg_Rgroomed_Nch[i][j]=new TH2D(Form("Zg_Vs_RatioGroomed_beta_%d_Nch_%d_%d",(int)beta_SD[i], trackbinbounds[j],trackbinboundsUpper[j] ),Form("Zg_Vs_RatioGroomed_beta_%d_Nch_%d_%d",(int)beta_SD[i],trackbinbounds[j],trackbinboundsUpper[j]),65,-0.05,0.6,110,-0.05,1.05);
         hRg_Ngroomed_Nch[i][j]=new TH2D(Form("Rg_Vs_NumberGroomed_beta_%d_Nch_%d_%d",(int)beta_SD[i],trackbinbounds[j],trackbinboundsUpper[j]),Form("Rg_Vs_NumberGroomed_beta_%d_Nch_%d_%d",(int)beta_SD[i],trackbinbounds[j],trackbinboundsUpper[j]),110,-0.05,1.05,120,0,120);
         hRg_Rgroomed_Nch[i][j]=new TH2D(Form("Rg_Vs_RatioGroomed_beta_%d_Nch_%d_%d",(int)beta_SD[i], trackbinbounds[j],trackbinboundsUpper[j]),Form("Rg_Vs_RatioGroomed_beta_%d_Nch_%d_%d",(int)beta_SD[i],trackbinbounds[j],trackbinboundsUpper[j]),110,-0.05,1.05,110,-0.05,1.05);

         hZgTg_Ngroomed_Nch[i][j]=new TH2D(Form("ZgTgBeta_Vs_NumberGroomed_beta_%d_Nch_%d_%d",(int)beta_SD[i],trackbinbounds[j],trackbinboundsUpper[j]),Form("ZgTgBeta_Vs_NumberGroomed_beta_%d_Nch_%d_%d",(int)beta_SD[i],trackbinbounds[j],trackbinboundsUpper[j]),110,-0.05,1.05,120,0,120);
         hZgTg_Rgroomed_Nch[i][j]=new TH2D(Form("ZgTgBeta_Vs_RatioGroomed_beta_%d_Nch_%d_%d",(int)beta_SD[i],trackbinbounds[j],trackbinboundsUpper[j]),Form("ZgTgBeta_Vs_RatioGroomed_beta_%d_Nch_%d_%d",(int)beta_SD[i],trackbinbounds[j],trackbinboundsUpper[j]),110,-0.05,1.05,110,-0.05,1.05);

         hdelta_phi_star_j1j2[i][j]=new TH1D(Form("delta_phi_star_between_subjet1_subjet2_beta_%d_Nch_%d_%d",(int)beta_SD[i],trackbinbounds[j],trackbinboundsUpper[j]) ,Form("delta_phi_star_between_subjet1_subjet2_beta_%d_Nch_%d_%d",(int)beta_SD[i],trackbinbounds[j],trackbinboundsUpper[j]),410,-0.1,4.00 );


		}
	}

	//END DEFINE Histogram
	

	std::cout<<"SD painter starts working"<<std::endl;
	std::cout<<"There are "<<SDfileList.size() << " in total"<<std::endl;


   JetDefinition jet_def(cambridge_algorithm ,JetDefinition::max_allowable_R );

	//*********************ENTER LOOP OF FIlE*****************************
	for(int ifile=0;ifile<SDfileList.size();ifile++){
		std::cout<<"Current file :" <<ifile+1 <<"/"<<SDfileList.size()<<std::endl;
		fFile=TFile::Open(SDfileList[ifile].c_str(),"READ");
		if(!fFile){
			std::cout<<"Can not open Current file "<<std::endl;
			continue;
		}
		TTree *tree=(TTree*)fFile->Get("SDtrackTree");
		LoadData(tree);
		std::cout<<"Total number of events of current file :"<<TotalEntry<<std::endl;
		//*******************************ENTER LOOP OF EVENT***********************
		for(int ievent=0;ievent<TotalEntry;ievent++){
			GetMyEntry(ievent);
			if(ievent%10000==0)std::cout<<"ievent="<<ievent<<std::endl;
         
			//*****************************ENTER LOOP OF JET***********************
			for(int ijet=0;ijet<(*jetPt).size();ijet++){

				long int NNtrk = ( (dau_pt)->at(ijet) ).size();// the number of daughter particles in this jet
	  			if( fabs( (*jetEta)[ijet] ) > jetEtaCut ) continue;
	  			if( fabs( (*jetPt)[ijet]  )< jetPtCut_Jet   ) continue;

    			int 	  Nch     =0;
    			int     Ngroomed[beta_bin]={0};
    			int     NnonG   [beta_bin]={0};

    			for(int A_trk=0; A_trk < NNtrk ; A_trk++ ){
    				   if( (* dau_chg )[ijet][A_trk] == 0 ) 	continue ; // not charged
            		if(fabs((* dau_pt )[ijet][A_trk])  < 0.3)     continue;//lab pt
                	if(fabs( (* dau_eta ) [ijet][A_trk]) > 2.4)     continue;//lab eta
                	Nch+= 1;
            }
            

            h_lab_JetMult_pT->Fill( (*jetPt)[ijet] , Nch );                   
    			h_lab_JetMult_phi->Fill((*jetPhi)[ijet], Nch );                   
 			   h_lab_JetMult_eta->Fill((*jetEta)[ijet], Nch );


            
           


 			    //**********************ENTER trackbin LOOP****************************************************
 			    for(int itrack=0;itrack<trackbin;itrack++){
 			    	if(  (Nch>=trackbinbounds[itrack]) &&  (Nch<trackbinboundsUpper[itrack]) ){
                 
 			    		//***********************ENTER inntrk LOOP**********************************************
 			    		for(int inntrk=0;inntrk<NNtrk;inntrk++){
 			    			   if( (* dau_chg )[ijet][inntrk] == 0 ) 	continue ; // not charged
            				if(fabs((* dau_pt )[ijet][inntrk])  < 0.3)     continue;//lab pt
                			if(fabs( (* dau_eta ) [ijet][inntrk]) > 2.4)     continue;//lab eta


                			double jet_dau_pt    =  ptWRTJet((double)((*jetPt))[ijet], (double)((*jetEta))[ijet]     , (double)((*jetPhi))[ijet]       , (double)(* (dau_pt)  )[ijet][inntrk], (double)((*dau_eta) )[ijet][inntrk], (double)(* (dau_phi) )[ijet][inntrk]);
                    		double jet_dau_eta   = etaWRTJet((double)((*jetPt))[ijet], (double)((*jetEta))[ijet]     , (double)((*jetPhi))[ijet]       , (double)(* (dau_pt)  )[ijet][inntrk], (double)((*dau_eta) )[ijet][inntrk], (double)(* (dau_phi)  )[ijet][inntrk]);
                    		double jet_dau_phi   = phiWRTJet((double)((*jetPt))[ijet], (double)((*jetEta))[ijet]     , (double)((*jetPhi))[ijet]       , (double)(* (dau_pt)   )[ijet][inntrk], (double)((*dau_eta))[ijet][inntrk], (double)(* (dau_phi)  )[ijet][inntrk]);

                    		h_jet_jT[itrack]->Fill(jet_dau_pt);
   							h_jet_etastar[itrack]->Fill(jet_dau_eta);


   							 //**************************ENTER ibeta LOOP **************************************
   							 for(int ibeta=0;ibeta<beta_bin;ibeta++){
                           
                           
   							 	if(  (*dau_Flag[ibeta])[ijet][inntrk] == 0  ) {
                              
   							 		h_jet_groomed_jT[ibeta][itrack]->Fill(jet_dau_pt);
   							 		h_jet_groomed_etastar[ibeta][itrack]->Fill(jet_dau_eta);
   							 		Ngroomed[ibeta]++;
   							 		
   							 	}

   							 	if(  ( (*dau_Flag[ibeta])[ijet][inntrk] ==1 ) || ( (*dau_Flag[ibeta])[ijet][inntrk] ==2 )  ) {
   							 		h_jet_non_groomed_jT[ibeta][itrack]->Fill(jet_dau_pt);
   							 		h_jet_non_groomed_etastar[ibeta][itrack]->Fill(jet_dau_eta);
   							 		NnonG[ibeta]++;
   							 	}


   							 }//*******************************END ibeta LOOP **********************************
 			    		}//*************************************END inntrk LOOP*********************************
 			    	}
 			    }//***********************END trackbin LOOP*****************************************************


 			    //*****only for check******
 			   for(int ibeta=0;ibeta<beta_bin;ibeta++){
              
               
 			    	if( (Ngroomed[ibeta] + NnonG[ibeta] )!=Nch ) {
                  std::cout<<"Warning"<<std::endl; waring_Flag++;
                  std::cout<<"Nch:"<<Nch<<std::endl;
                  std::cout<<"Ngroomed:"<<Ngroomed[ibeta]<<std::endl;
                  std::cout<<"NnonG:"<<NnonG[ibeta]<<std::endl;

               }
 			   }



  				//***********************************ENTER ibeta LOOP*******************************************
 			    for(int ibeta=0;ibeta<beta_bin;ibeta++){

 			    	hNch_Zg[ibeta]->Fill( Nch ,(*jetZg)[ijet][ibeta] );
					hNch_Rg[ibeta]->Fill( Nch ,(*jetRg)[ijet][ibeta] );
					hNch_ZgRg[ibeta]->Fill( Nch,(*jetZg)[ijet][ibeta] *  TMath::Power( (*jetRg)[ijet][ibeta]/R, beta_SD[ibeta]  )    );
					hZg_Rg[ibeta]->Fill( (*jetZg)[ijet][ibeta],(*jetRg)[ijet][ibeta] );
					hlogZg_logRg[ibeta]->Fill( TMath::Log( (*jetZg)[ijet][ibeta] ) / TMath::Log(10.)   , TMath::Log( (*jetRg)[ijet][ibeta] )/ TMath::Log(10.)   );

               hZg_Ngroomed[ibeta]->Fill( (*jetZg)[ijet][ibeta], Ngroomed[ibeta]  );
               hZg_Rgroomed[ibeta]->Fill( (*jetZg)[ijet][ibeta], (Ngroomed[ibeta]*1.0) / (Nch*1.0) );
               hRg_Ngroomed[ibeta]->Fill( (*jetRg)[ijet][ibeta], Ngroomed[ibeta] );
               hRg_Rgroomed[ibeta]->Fill( (*jetRg)[ijet][ibeta], (Ngroomed[ibeta]*1.0) / (Nch*1.0) );

               hZgTg_Ngroomed[ibeta]->Fill( (*jetZg)[ijet][ibeta] *  TMath::Power( (*jetRg)[ijet][ibeta]/R, beta_SD[ibeta]  ) ,   Ngroomed[ibeta]                        );
               hZgTg_Rgroomed[ibeta]->Fill( (*jetZg)[ijet][ibeta] *  TMath::Power( (*jetRg)[ijet][ibeta]/R, beta_SD[ibeta]  ) ,  (Ngroomed[ibeta]*1.0) / (Nch*1.0)  );

               //calculate the phi_star of subjet1 and subjet2 in jet frame.
               double jet1_phi_star   = phiWRTJet( (double)((*jetPt))[ijet], (double)((*jetEta))[ijet]     , (double)((*jetPhi))[ijet]       , (double)(*subjet1_pt   )[ijet][ibeta], (double)( *subjet1_eta )[ijet][ibeta], (double)(*subjet1_phi)[ijet][ibeta]);
               double jet2_phi_satr   = phiWRTJet( (double)((*jetPt))[ijet], (double)((*jetEta))[ijet]     , (double)((*jetPhi))[ijet]       , (double)(*subjet2_pt   )[ijet][ibeta], (double)( *subjet2_eta )[ijet][ibeta], (double)(*subjet2_phi)[ijet][ibeta]);
               double delta_phi_star_sj1_sj2= TMath::ACos(TMath::Cos(jet1_phi_star - jet2_phi_satr));




					//********************************ENTER itrack LOOP****************************************
					for(int itrack=0;itrack<trackbin;itrack++){
						if(  (Nch>=trackbinbounds[itrack]) &&  (Nch<trackbinboundsUpper[itrack]) ){
							hZg_Rg_Nch[ibeta][itrack]->Fill( (*jetZg)[ijet][ibeta],(*jetRg)[ijet][ibeta] );
						   hlogZg_logRg_Nch[ibeta][itrack]->Fill( TMath::Log( (*jetZg)[ijet][ibeta] ) / TMath::Log(10.)   , TMath::Log( (*jetRg)[ijet][ibeta] )/ TMath::Log(10.) );
								

							h_Ngroomed[ibeta][itrack]->Fill(Ngroomed[ibeta]);
							h_ratio_groomed[ibeta][itrack]->Fill( (Ngroomed[ibeta]*1.0) / (Nch*1.0) );

							//this actually should be Zg * (Rg/R)^{\beta}
							hZgRg_dis[ibeta][itrack]->Fill( (*jetZg)[ijet][ibeta] *  TMath::Power( (*jetRg)[ijet][ibeta]/R, beta_SD[ibeta]  )  );
							hZg_dis[ibeta][itrack]->Fill( (*jetZg)[ijet][ibeta] );
							hRg_dis[ibeta][itrack]->Fill( (*jetRg)[ijet][ibeta] );



                     hZg_Ngroomed_Nch[ibeta][itrack]->Fill( (*jetZg)[ijet][ibeta], Ngroomed[ibeta]  );
                     hZg_Rgroomed_Nch[ibeta][itrack]->Fill( (*jetZg)[ijet][ibeta], (Ngroomed[ibeta]*1.0) / (Nch*1.0) );
                     hRg_Ngroomed_Nch[ibeta][itrack]->Fill( (*jetRg)[ijet][ibeta], Ngroomed[ibeta] );
                     hRg_Rgroomed_Nch[ibeta][itrack]->Fill( (*jetRg)[ijet][ibeta], (Ngroomed[ibeta]*1.0) / (Nch*1.0) );

                     hZgTg_Ngroomed_Nch[ibeta][itrack]->Fill( (*jetZg)[ijet][ibeta] *  TMath::Power( (*jetRg)[ijet][ibeta]/R, beta_SD[ibeta]  ) ,   Ngroomed[ibeta]                         );
                     hZgTg_Rgroomed_Nch[ibeta][itrack]->Fill( (*jetZg)[ijet][ibeta] *  TMath::Power( (*jetRg)[ijet][ibeta]/R, beta_SD[ibeta]  ) ,  (Ngroomed[ibeta]*1.0) / (Nch*1.0)  );

                     hdelta_phi_star_j1j2[ibeta][itrack]->Fill( delta_phi_star_sj1_sj2 );

							
						}
					}//*******************************END itrack LOOP******************************************
 			    }//***********************************END ibeta LOOP********************************************





			}//*******************************ENTER LOOP OF JET*******************
		}//*******************************END LOOP OF EVENT*************************
		fFile->Close();	
	}//********************************END LOOP OF FILE***************************************


    TFile *fout=new TFile(store_dir.c_str(),"RECREATE");
	 h_lab_JetMult_pT->Write();                     //
    h_lab_JetMult_phi->Write();   //
    h_lab_JetMult_eta->Write();                   //
    //IN THE JET FRAME ,pt etastar of jet daughter particles
   																							 //

    for(int i=0;i<trackbin;i++){
    	h_jet_jT[i]->Write();
    	h_jet_etastar[i]->Write();
    }

    
    for(int i=0;i<beta_bin;i++){
    	for(int j=0;j<trackbin;j++){
    		h_jet_groomed_jT[i][j]   ->Write();
    		h_jet_non_groomed_jT[i][j] ->Write();
    		h_jet_groomed_etastar[i][j] ->Write();
    		h_jet_non_groomed_etastar[i][j]->Write();


    	}
    }

   

	
	for(int i=0;i<beta_bin;i++){
		


		hNch_Zg[i]->Write();
		hNch_Rg[i]->Write();
		hNch_ZgRg[i]->Write();
		hZg_Rg[i]->Write();
		hlogZg_logRg[i]->Write();


      hZg_Ngroomed[i]->Write();
      hZg_Rgroomed[i]->Write();
      hRg_Ngroomed[i]->Write();
      hRg_Rgroomed[i]->Write();

      hZgTg_Ngroomed[i]->Write();
      hZgTg_Rgroomed[i]->Write();



		for(int j=0;j<trackbin;j++){
			hZg_Rg_Nch[i][j]->Write();
			hlogZg_logRg_Nch[i][j]->Write();

			h_Ngroomed[i][j]->Write();
			h_ratio_groomed[i][j]->Write();
			hZgRg_dis[i][j]->Write();
			
			hZg_dis[i][j]->Write();
			hRg_dis[i][j]->Write();

         hZg_Ngroomed_Nch[i][j]->Write();
         hZg_Rgroomed_Nch[i][j]->Write();
         hRg_Ngroomed_Nch[i][j]->Write();
         hRg_Rgroomed_Nch[i][j]->Write();

         hZgTg_Ngroomed_Nch[i][j]->Write();
         hZgTg_Rgroomed_Nch[i][j]->Write();

         hdelta_phi_star_j1j2[i][j]->Write();

		}	
	}

	fout->Close();




	if(waring_Flag==0) std::cout<<"GOOD LUCK!"<<std::endl;

}








