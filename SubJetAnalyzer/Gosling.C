//#include "minclude/SDfactory.h"
//#include "minclude/SubstructureWorker.h"

#include "minclude/SubJetAnalyzer.h"
#include "minclude/DrawFlow.h"
#include "minclude/Draw2DCorr.h"








int main(){
	std::vector<std::string> dir_in;
	std::vector<std::string> dir_out;
	
	

	//***************Pythia8 inclusive**************************** 
	//std::string store="/home/fl52/ppJetAnalysis/Pythia_CP5_pp/new_result/Substructure_inclusive.root";
	//dir_in.push_back( "/home/fl52/ppJetAnalysis/Pythia_CP5_pp/Pythia_CP5_pp_rootfile/pp_highMultGen_CP5_inclusive.root" );
	//dir_out.push_back("/home/fl52/ppJetAnalysis/Pythia_CP5_pp/new_SDrootfile/SDpp_highMultGen_CP5_inclusive.root");
	std::string store="/home/fl52/ppJetAnalysis/Pythia_CP5_pp/new_result";
	//for(int i=0;i<18;i++){
	//	//dir_in.push_back( Form("/home/fl52/ppJetAnalysis/Pythia_CP5_pp/Pythia_CP5_pp_rootfile/pp_highMultGen_CP5_inclusive_%d.root",i+1) );
		//dir_out.push_back( Form("/home/fl52/ppJetAnalysis/Pythia_CP5_pp/new_SDrootfile/SDpp_highMultGen_CP5_inclusive_%d.root",i+1) );
	//}
	//************************************************************


	//*****************Pythia8 high Nch***************************
	//std::string store="/home/fl52/ppJetAnalysis/Pythia_CP5_pp/new_result/Substructure_highNch.root";
	for(int i=0;i<132;i++){
	//	dir_in.push_back( Form("/home/fl52/ppJetAnalysis/Pythia_CP5_pp/Pythia_CP5_pp_rootfile/pp_highMultGen_nChGT60_%d.root ",i+1) );
		dir_out.push_back( Form("/home/fl52/ppJetAnalysis/Pythia_CP5_pp/new_SDrootfile/SDpp_highMultGen_nChGT60_%d.root ",i+1) );
	}
	//************************************************************

	//*****************parton_cascade***************************
	//std::string store="/home/fl52/ppJetAnalysis/pp_parton_cascade/new_result/Substructure_cascade.root";
	/*for(int i=62;i<111;i++){
		dir_in.push_back( Form("/home/fl52/ppJetAnalysis/pp_parton_cascade/rootfile/pp_parton_cascade_%d.root",i) );
		dir_out.push_back( Form("/home/fl52/ppJetAnalysis/pp_parton_cascade/new_SDrootfile/SDpp_parton_cascade_%d.root ",i) );
	}
	//************************************************************

	for(int i=0;i<10;i++){
		dir_in.push_back( Form("/home/fl52/parton_cascade_high_Nch/temp/parton_cascadeNch60_%d.root",i+1) );
		dir_out.push_back( Form("/home/fl52/ppJetAnalysis/pp_parton_cascade/new_SDrootfile/SDparton_cascadeNch60_%d.root ",i+1) );
	}*/

	//SDfactory sdF(dir_in,dir_out);
	//sdF.SDtreePlanter();
	//SDSubstructureWorker worker(dir_out,store);
	//worker.SDpainter();
	SubJetAnalyzer sjAna(dir_out,store);
	sjAna.CalculateCorrelation();
	sjAna.DrawFlow();
	sjAna.Draw2DCorrelation();
	return 0;

}