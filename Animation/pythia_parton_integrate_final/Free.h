void PartonShower::FreeStreaming()
{
    //string random_str = string(argc[1]);
    int Spatial_mode = 1; // 0: use the accumulant method; 1: uses the free-streaming to formation time; 
   
    // Output the parton information
    string output_filename;
    output_filename = "parton_info.dat";
    ofstream output_parton(output_filename.c_str());  //-----zhao

    if( ! output_parton.is_open() ) {
        cout << "cannot open output file:"<< endl
             << output_filename << endl;
        return ;
    }




    // Number of events to generate
    int numEvent = 0;
    int Nparton = 0;
    for (int j = 0; j < (*event).size(); ++j) {
            
            
            if ((*event)[j].isFinal() && (*event)[j].isParton()) {
                Nparton++;
                
            }
    }

    output_parton << "# pdgid  px  py  pz  energy  x  y  z  t" << std::endl;
    output_parton << "# event id " << numEvent  <<  ", Number_of_parton  "  << std::endl;
    output_parton << Nparton << std::endl;
    for (int j = 0; j < (*event).size(); ++j) {
        Particle& particle = (*event)[j];
        
        if ( particle.isFinal() && particle.isParton() ) {
            double position[3] = {0.0}; // spatial information of partons (x, y, z)
            position[0] = 0.0;
            position[1] = 0.0;
            position[2] = 0.0;
            //......formation time + ......
            //int ishower = 0;
            double p0[4] = {0.0}; // particle's four momentum
            double p4[4] = {0.0}; // mother's four momentum
		    double qt, time_step;
		    double timeplus = 0.0;
            int IDmom1, IDmom2;
            int timebreaker = 0;
            int IDmom0 = j;
            //... start to calculate the formation and the spatial information of partons ...
            while (timebreaker == 0) {
                int IDiii = IDmom0;
                if (abs( (*event)[IDiii].status() )==23 || abs( (*event)[IDiii].status() )==21 ||
                    abs( (*event)[IDiii].status() )==12 ) timebreaker=1;			
                IDmom1 = (*event)[IDiii].mother1();
                IDmom2 = (*event)[IDiii].mother2();
                if (IDmom1==IDmom2 && IDmom1==0) timebreaker=1;
                if (IDmom1==IDmom2 && IDmom1>0) IDmom0=IDmom1;
                if (IDmom1>0 && IDmom2==0) {
                    IDmom0=IDmom1;
                    double IDdaughter1 = (*event)[IDmom0].daughter1();
                    double IDdaughter2 = (*event)[IDmom0].daughter2();
                    if (IDdaughter1 != IDdaughter2 && IDdaughter1>0 && IDdaughter2>0) {
                        p4[0] = (*event)[IDdaughter1].e()+(*event)[IDdaughter2].e();
                        p4[1] = (*event)[IDdaughter1].px()+(*event)[IDdaughter2].px();
                        p4[2] = (*event)[IDdaughter1].py()+(*event)[IDdaughter2].py();
                        p4[3] = (*event)[IDdaughter1].pz()+(*event)[IDdaughter2].pz();
                        double x_split = (*event)[IDiii].e()/p4[0];
                        if (x_split>1) x_split=1.0/x_split;
                        p0[0] = (*event)[IDiii].e();
                        p0[1] = (*event)[IDiii].px();
                        p0[2] = (*event)[IDiii].py();
                        p0[3] = (*event)[IDiii].pz(); 
                        //double pt_daughter = sqrt(pow(p0[1],2)+pow(p0[2],2));
                        //double pt_mother = sqrt(pow(p4[1],2)+pow(p4[2],2));			  

                        rotate(p4[1],p4[2],p4[3],p0,1);
                        qt = sqrt(pow(p0[1],2)+pow(p0[2],2));
                        rotate(p4[1],p4[2],p4[3],p0,-1);
                        double kt_daughter=qt;
                        if (x_split<0.5) {
                            if(kt_daughter > 0.0001) {
                                time_step = 2.0*p4[0]*x_split*(1-x_split)/pow(kt_daughter,2)*0.197;
                                timeplus = timeplus + time_step;
                                position[0] = position[0] + time_step * p4[1]/p4[0];
                                position[1] = position[1] + time_step * p4[2]/p4[0];
                                position[2] = position[2] + time_step * p4[3]/p4[0];
                            }
                        } else {
                            if(kt_daughter > 0.0001) {
                                time_step = 1/p4[0];
                                timeplus = timeplus + time_step;
                                position[0] = position[0] + time_step * p4[1]/p4[0];
                                position[1] = position[1] + time_step * p4[2]/p4[0];
                                position[2] = position[2] + time_step * p4[3]/p4[0];
                            }
                        }
                    }
                }//if(IDmom1>0 && IDmom2==0)

                if (IDmom1<IDmom2 && IDmom1>0 && IDmom2>0) {
                    if((*event)[IDmom1].e()>(*event)[IDmom2].e()) {
                        IDmom0=IDmom1;			
                    }
                    if((*event)[IDmom1].e()<=(*event)[IDmom2].e()) {
                        IDmom0=IDmom2;			
                    }
                    double IDdaughter1=(*event)[IDmom0].daughter1();
                    double IDdaughter2=(*event)[IDmom0].daughter2();
                    if (IDdaughter1 != IDdaughter2 && IDdaughter1>0 && IDdaughter2>0) {
                        p4[0]=(*event)[IDdaughter1].e()+(*event)[IDdaughter2].e();
                        p4[1]=(*event)[IDdaughter1].px()+(*event)[IDdaughter2].px();
                        p4[2]=(*event)[IDdaughter1].py()+(*event)[IDdaughter2].py();
                        p4[3]=(*event)[IDdaughter1].pz()+(*event)[IDdaughter2].pz();
                        double x_split = (*event)[IDiii].e()/p4[0];
                        if(x_split>1) x_split=1.0/x_split;

                        p0[0]=(*event)[IDiii].e();
                        p0[1]=(*event)[IDiii].px();
                        p0[2]=(*event)[IDiii].py();
                        p0[3]=(*event)[IDiii].pz(); 
                        rotate(p4[1],p4[2],p4[3],p0,1);
                        qt=sqrt(pow(p0[1],2)+pow(p0[2],2));
                        rotate(p4[1],p4[2],p4[3],p0,-1);
                        double kt_daughter=qt;

                        if (x_split<0.5) {
			                if (kt_daughter > 0.0001) {
                                time_step = 2.0*p4[0]*x_split*(1-x_split)/pow(kt_daughter,2)*0.197;
                                timeplus = timeplus + time_step;
                                position[0] = position[0] + time_step * p4[1]/p4[0];
                                position[1] = position[1] + time_step * p4[2]/p4[0];
                                position[2] = position[2] + time_step * p4[3]/p4[0];
                                    
                            }
                        } else {
                            if(kt_daughter > 0.0001) {
                                time_step = 1/p4[0];
                                timeplus = timeplus + time_step;
                                position[0] = position[0] + time_step * p4[1]/p4[0];
                                position[1] = position[1] + time_step * p4[2]/p4[0];
                                position[2] = position[2] + time_step * p4[3]/p4[0];
                            }
                        }
                    }
			
		            if ( (IDdaughter1 > 0 && IDdaughter2 == 0) || (IDdaughter1 > 0 && IDdaughter2 == IDdaughter1)) {
                        p4[0] = (*event)[IDdaughter1].e();
                        p4[1] = (*event)[IDdaughter1].px();
                        p4[2] = (*event)[IDdaughter1].py();
                        p4[3] = (*event)[IDdaughter1].pz();
                        double x_split = (*event)[IDiii].e()/p4[0];
                        if (x_split>1) x_split=1.0/x_split; // revise the mother and daughter
                        p0[0] = (*event)[IDiii].e();
                        p0[1] = (*event)[IDiii].px();
                        p0[2] = (*event)[IDiii].py();
                        p0[3] = (*event)[IDiii].pz(); 
                        //double pt_daughter=sqrt(pow(p0[1],2)+pow(p0[2],2));
                        //double pt_mother=sqrt(pow(p4[1],2)+pow(p4[2],2));			 
                        rotate(p4[1],p4[2],p4[3],p0,1); // rotate into the 
                        qt = sqrt(pow(p0[1],2)+pow(p0[2],2));
                        rotate(p4[1],p4[2],p4[3],p0,-1);
                        double kt_daughter = qt;
                        //double Q2 = 1.0/(x_split*(1-x_split)/pow(kt_daughter,2));
                        if(x_split<0.5){
			                if (kt_daughter > 0.0001) {						
                                time_step = 2.0*p4[0]*x_split*(1-x_split)/pow(kt_daughter,2)*0.197;
                                timeplus = timeplus + time_step;
                                position[0] = position[0] + time_step * p4[1]/p4[0];
                                position[1] = position[1] + time_step * p4[2]/p4[0];
                                position[2] = position[2] + time_step * p4[3]/p4[0];					
			                }
                        } else {
                            if(kt_daughter > 0.0001) {
                                time_step = 1/p4[0];
                                timeplus = timeplus + time_step;
                                position[0] = position[0] + time_step * p4[1]/p4[0];
                                position[1] = position[1] + time_step * p4[2]/p4[0];
                                position[2] = position[2] + time_step * p4[3]/p4[0];
                            }
                        }
                    }
                        
                }//if(IDmom1<IDmom2 && IDmom1>0 && IDmom2>0)

                if (IDmom1>IDmom2 && IDmom1>0 && IDmom2>0) {
                    if ((*event)[IDmom1].e()>(*event)[IDmom2].e()) {
                        IDmom0=IDmom1;			
                    }
                    if ((*event)[IDmom1].e()<=(*event)[IDmom2].e()) {
                        IDmom0=IDmom2;			
                    }			
                    double IDdaughter1=(*event)[IDmom0].daughter1();
                    double IDdaughter2=(*event)[IDmom0].daughter2();
                    if (IDdaughter1 != IDdaughter2 && IDdaughter1>0 && IDdaughter2>0) {
                        p4[0]=(*event)[IDdaughter1].e()+ (*event)[IDdaughter2].e();
                        p4[1]=(*event)[IDdaughter1].px()+(*event)[IDdaughter2].px();
                        p4[2]=(*event)[IDdaughter1].py()+(*event)[IDdaughter2].py();
                        p4[3]=(*event)[IDdaughter1].pz()+(*event)[IDdaughter2].pz();
                        double x_split=(*event)[IDiii].e()/p4[0];
                        if(x_split>1) x_split=1.0/x_split;
                        p0[0]=(*event)[IDiii].e();
                        p0[1]=(*event)[IDiii].px();
                        p0[2]=(*event)[IDiii].py();
                        p0[3]=(*event)[IDiii].pz(); 
                        rotate(p4[1],p4[2],p4[3],p0,1);
                        qt=sqrt(pow(p0[1],2)+pow(p0[2],2));
                        rotate(p4[1],p4[2],p4[3],p0,-1);
                        double kt_daughter=qt;
                        if(x_split<0.5){
			                if (kt_daughter > 0.0001) {						
                                time_step = 2.0*p4[0]*x_split*(1-x_split)/pow(kt_daughter,2)*0.197;
                                timeplus = timeplus + time_step;
                                position[0] = position[0] + time_step * p4[1]/p4[0];
                                position[1] = position[1] + time_step * p4[2]/p4[0];
                                position[2] = position[2] + time_step * p4[3]/p4[0];					
			                }
                        } else {
                            if(kt_daughter > 0.0001) {
                                time_step = 1/p4[0];
                                timeplus = timeplus + time_step;
                                position[0] = position[0] + time_step * p4[1]/p4[0];
                                position[1] = position[1] + time_step * p4[2]/p4[0];
                                position[2] = position[2] + time_step * p4[3]/p4[0];
                            }
                        }
                    }
                        
                    if ( (IDdaughter1 > 0 && IDdaughter2 == 0) || (IDdaughter1 > 0 && IDdaughter2 == IDdaughter1)) {
                        p4[0] = (*event)[IDdaughter1].e();
                        p4[1] = (*event)[IDdaughter1].px();
                        p4[2] = (*event)[IDdaughter1].py();
                        p4[3] = (*event)[IDdaughter1].pz();
                        double x_split = (*event)[IDiii].e()/p4[0];
                        if (x_split>1) x_split=1.0/x_split; // revise the mother and daughter
                        p0[0] = (*event)[IDiii].e();
                        p0[1] = (*event)[IDiii].px();
                        p0[2] = (*event)[IDiii].py();
                        p0[3] = (*event)[IDiii].pz(); 
                        //double pt_daughter=sqrt(pow(p0[1],2)+pow(p0[2],2));
                        //double pt_mother=sqrt(pow(p4[1],2)+pow(p4[2],2));			 
                        rotate(p4[1],p4[2],p4[3],p0,1); // rotate into the 
                        qt = sqrt(pow(p0[1],2)+pow(p0[2],2));
                        rotate(p4[1],p4[2],p4[3],p0,-1);
                        double kt_daughter = qt;
                        //double Q2 = 1.0/(x_split*(1-x_split)/pow(kt_daughter,2));
                        if(x_split<0.5){
			                if (kt_daughter > 0.0001) {						
                                time_step = 2.0*p4[0]*x_split*(1-x_split)/pow(kt_daughter,2)*0.197;
                                timeplus = timeplus + time_step;
                                position[0] = position[0] + time_step * p4[1]/p4[0];
                                position[1] = position[1] + time_step * p4[2]/p4[0];
                                position[2] = position[2] + time_step * p4[3]/p4[0];					
			                }
                        } else {
                            if(kt_daughter > 0.0001) {
                                time_step = 1/p4[0];
                                timeplus = timeplus + time_step;
                                position[0] = position[0] + time_step * p4[1]/p4[0];
                                position[1] = position[1] + time_step * p4[2]/p4[0];
                                position[2] = position[2] + time_step * p4[3]/p4[0];
                            }
                        }
                    }
                        
                }//if(IDmom1>IDmom2 && IDmom1>0 && IDmom2>0)
                    
                if (IDmom1==IDmom2 && IDmom1>0) {
                    IDmom0=IDmom1;
                    int IDdaughter1 = (*event)[IDmom0].daughter1();
                    int IDdaughter2 = (*event)[IDmom0].daughter2();
                        
                    if (IDdaughter1 != IDdaughter2 && IDdaughter1>0 && IDdaughter2>0) {
                        p4[0] = (*event)[IDdaughter1].e()+(*event)[IDdaughter2].e();
                        p4[1] = (*event)[IDdaughter1].px()+(*event)[IDdaughter2].px();
                        p4[2] = (*event)[IDdaughter1].py()+(*event)[IDdaughter2].py();
                        p4[3] = (*event)[IDdaughter1].pz()+(*event)[IDdaughter2].pz();
                        double x_split = (*event)[IDiii].e()/p4[0];
                        if (x_split>1) x_split=1.0/x_split; // revise the mother and daughter
                        p0[0] = (*event)[IDiii].e();
                        p0[1] = (*event)[IDiii].px();
                        p0[2] = (*event)[IDiii].py();
                        p0[3] = (*event)[IDiii].pz(); 
                        rotate(p4[1],p4[2],p4[3],p0,1); // rotate into the 
                        qt = sqrt(pow(p0[1],2)+pow(p0[2],2));
                        rotate(p4[1],p4[2],p4[3],p0,-1);
                        double kt_daughter = qt;
                        if(x_split<0.5){
			                if (kt_daughter > 0.0001) {						
                                time_step = 2.0*p4[0]*x_split*(1-x_split)/pow(kt_daughter,2)*0.197;
                                timeplus = timeplus + time_step;
                                position[0] = position[0] + time_step * p4[1]/p4[0];
                                position[1] = position[1] + time_step * p4[2]/p4[0];
                                position[2] = position[2] + time_step * p4[3]/p4[0];					
			                }
                        } else {
                            if(kt_daughter > 0.0001) {
                                time_step = 1/p4[0];
                                timeplus = timeplus + time_step;
                                position[0] = position[0] + time_step * p4[1]/p4[0];
                                position[1] = position[1] + time_step * p4[2]/p4[0];
                                position[2] = position[2] + time_step * p4[3]/p4[0];
                            }
                        }
                    }
                        
                    if ( (IDdaughter1 > 0 && IDdaughter2 == 0) || (IDdaughter1 > 0 && IDdaughter2 == IDdaughter1)) {
                        p4[0] = (*event)[IDdaughter1].e();
                        p4[1] = (*event)[IDdaughter1].px();
                        p4[2] = (*event)[IDdaughter1].py();
                        p4[3] = (*event)[IDdaughter1].pz();
                        double x_split = (*event)[IDiii].e()/p4[0];
                        if (x_split>1) x_split=1.0/x_split; // revise the mother and daughter
                        p0[0] = (*event)[IDiii].e();
                        p0[1] = (*event)[IDiii].px();
                        p0[2] = (*event)[IDiii].py();
                        p0[3] = (*event)[IDiii].pz(); 
                        rotate(p4[1],p4[2],p4[3],p0,1); // rotate into the 
                        qt = sqrt(pow(p0[1],2)+pow(p0[2],2));
                        rotate(p4[1],p4[2],p4[3],p0,-1);
                        double kt_daughter = qt;
                        //if (kt_daughter <= 0.0001) kt_daughter = 0.0001;
                        if(x_split<0.5){
			                 if (kt_daughter > 0.0001) {						
                                time_step = 2.0*p4[0]*x_split*(1-x_split)/pow(kt_daughter,2)*0.197;
                                timeplus = timeplus + time_step;
                                position[0] = position[0] + time_step * p4[1]/p4[0];
                                position[1] = position[1] + time_step * p4[2]/p4[0];
                                position[2] = position[2] + time_step * p4[3]/p4[0];					
			                }
                        } else {
                            if(kt_daughter > 0.0001) {
                                time_step = 1/p4[0];
                                timeplus = timeplus + time_step;
                                position[0] = position[0] + time_step * p4[1]/p4[0];
                                position[1] = position[1] + time_step * p4[2]/p4[0];
                                position[2] = position[2] + time_step * p4[3]/p4[0];
                            }
                        }
                    }
                }
                    
            } //while(timebreaker == 0)

            // Extract relevant information
            int pdgId = particle.id();
            /*
            int status = particle.status();
            int mother1 = particle.mother1();
            int mother2 = particle.mother2();
            int daughter1 = particle.daughter1();
            int daughter2 = particle.daughter2();
            */
            double px = particle.px();
            double py = particle.py();
            double pz = particle.pz();
            double energy = particle.e();
            if (Spatial_mode == 1) {
                position[0] = 0.0 + timeplus * px / energy;
                position[1] = 0.0 + timeplus * py / energy;
                position[2] = 0.0 + timeplus * pz / energy;
            }
            output_parton << pdgId << "  " 
                            << px << "  " << py << "  " << pz << "  " << energy <<"  "
                            << position[0] << "  "<< position[1] << "  " << position[2] << "  " << timeplus << "  "
                            << particle.col() << "  " << particle.acol()
                            << std::endl;
            
        }//******if (is parton and is final)******
    }//*******loop of j************* 
    output_parton.close();
    
}










