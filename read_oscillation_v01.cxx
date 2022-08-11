#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include "stdlib.h"
using namespace std;

#include<vector>
#include<map>
#include<set>

#include "WCPLEEANA/TOsc.h"

#include "WCPLEEANA/Configure_Osc.h"

#include "TApplication.h"

//#include <chrono> // timer
//auto time_start = chrono::high_resolution_clock::now();
//auto time_stop = chrono::high_resolution_clock::now();
//auto time_duration = chrono::duration_cast<chrono::seconds>(time_stop - time_start);
//cout<<endl<<" ---> check time duration "<<time_duration.count()<<endl<<endl;
//milliseconds, minutes

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// MAIN //////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
  TString roostr = "";

  cout<<endl<<" ---> A story ..."<<endl<<endl;

  int ifile = 1;
  double scaleF_POT_BNB  = 1;
  double scaleF_POT_NuMI = 1;
  int display = 0;

  int it14 = 0;
  int idm2 = 0;
  int it24 = 0;
  
  for(int i=1; i<argc; i++) {
    if( strcmp(argv[i],"-f")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>ifile ) ) { cerr<<" ---> Error ifile !"<<endl; exit(1); }
    }
    if( strcmp(argv[i],"-pbnb")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>scaleF_POT_BNB ) ) { cerr<<" ---> Error scaleF_POT_BNB !"<<endl; exit(1); }
    }
    if( strcmp(argv[i],"-pnumi")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>scaleF_POT_NuMI ) ) { cerr<<" ---> Error scaleF_POT_NuMI !"<<endl; exit(1); }
    }    
    if( strcmp(argv[i],"-d")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>display ) ) { cerr<<" ---> Error display !"<<endl; exit(1); }
    }
    if( strcmp(argv[i],"-it14")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>it14 ) ) { cerr<<" ---> Error it14 !"<<endl; exit(1); }
    }
    if( strcmp(argv[i],"-idm2")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>idm2 ) ) { cerr<<" ---> Error idm2 !"<<endl; exit(1); }
    }
    if( strcmp(argv[i],"-it24")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>it24 ) ) { cerr<<" ---> Error it24 !"<<endl; exit(1); }
    }    
  }

  ///////////////////////////////////////////////////////////
  
  if( !display ) {
    gROOT->SetBatch( 1 );
  }
  
  TApplication theApp("theApp",&argc,argv);
  
  /////////////////////////////////////////////////////////// Draw style

  gStyle->SetOptStat(0);
  //gStyle->SetPalette(kBird);

  double snWidth = 2;

  // use medium bold lines and thick markers
  gStyle->SetLineWidth(snWidth);
  gStyle->SetFrameLineWidth(snWidth);
  gStyle->SetHistLineWidth(snWidth);
  gStyle->SetFuncWidth(snWidth);
  gStyle->SetGridWidth(snWidth);
  gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.2);
  gStyle->SetEndErrorSize(4);
  gStyle->SetEndErrorSize(0);

  ///////////////////////////////////////////////////////////

  TOsc *osc_test = new TOsc();

  ///////

  osc_test->scaleF_POT_BNB  = scaleF_POT_BNB;
  osc_test->scaleF_POT_NuMI = scaleF_POT_NuMI;
     
  ///////

  osc_test->flag_syst_dirt   = Configure_Osc::flag_syst_dirt;
  osc_test->flag_syst_mcstat = Configure_Osc::flag_syst_mcstat;
  osc_test->flag_syst_flux   = Configure_Osc::flag_syst_flux;
  osc_test->flag_syst_geant  = Configure_Osc::flag_syst_geant;
  osc_test->flag_syst_Xs     = Configure_Osc::flag_syst_Xs;
  osc_test->flag_syst_det    = Configure_Osc::flag_syst_det;
  
  ///////

  osc_test->flag_NuMI_nueCC_from_intnue       = Configure_Osc::flag_NuMI_nueCC_from_intnue;
  osc_test->flag_NuMI_nueCC_from_overlaynumu  = Configure_Osc::flag_NuMI_nueCC_from_overlaynumu;
  osc_test->flag_NuMI_nueCC_from_appnue       = Configure_Osc::flag_NuMI_nueCC_from_appnue;
  osc_test->flag_NuMI_nueCC_from_appnumu      = Configure_Osc::flag_NuMI_nueCC_from_appnumu;
  osc_test->flag_NuMI_nueCC_from_overlaynueNC = Configure_Osc::flag_NuMI_nueCC_from_overlaynueNC;
  osc_test->flag_NuMI_nueCC_from_overlaynumuNC= Configure_Osc::flag_NuMI_nueCC_from_overlaynumuNC;
  
  osc_test->flag_NuMI_numuCC_from_overlaynumu   = Configure_Osc::flag_NuMI_numuCC_from_overlaynumu;
  osc_test->flag_NuMI_numuCC_from_overlaynue    = Configure_Osc::flag_NuMI_numuCC_from_overlaynue;
  osc_test->flag_NuMI_numuCC_from_appnue        = Configure_Osc::flag_NuMI_numuCC_from_appnue;
  osc_test->flag_NuMI_numuCC_from_appnumu       = Configure_Osc::flag_NuMI_numuCC_from_appnumu;
  osc_test->flag_NuMI_numuCC_from_overlaynumuNC = Configure_Osc::flag_NuMI_numuCC_from_overlaynumuNC;
  osc_test->flag_NuMI_numuCC_from_overlaynueNC  = Configure_Osc::flag_NuMI_numuCC_from_overlaynueNC;  
  
  osc_test->flag_NuMI_CCpi0_from_overlaynumu  = Configure_Osc::flag_NuMI_CCpi0_from_overlaynumu;
  osc_test->flag_NuMI_CCpi0_from_appnue       = Configure_Osc::flag_NuMI_CCpi0_from_appnue;
  osc_test->flag_NuMI_CCpi0_from_overlaynumuNC= Configure_Osc::flag_NuMI_CCpi0_from_overlaynumuNC;
  osc_test->flag_NuMI_CCpi0_from_overlaynueNC = Configure_Osc::flag_NuMI_CCpi0_from_overlaynueNC;
  
  osc_test->flag_NuMI_NCpi0_from_overlaynumu  = Configure_Osc::flag_NuMI_NCpi0_from_overlaynumu;
  osc_test->flag_NuMI_NCpi0_from_appnue       = Configure_Osc::flag_NuMI_NCpi0_from_appnue;
  osc_test->flag_NuMI_NCpi0_from_overlaynumuNC= Configure_Osc::flag_NuMI_NCpi0_from_overlaynumuNC;
  osc_test->flag_NuMI_NCpi0_from_overlaynueNC = Configure_Osc::flag_NuMI_NCpi0_from_overlaynueNC;


  ///////
  
  osc_test->flag_BNB_nueCC_from_intnue       = Configure_Osc::flag_BNB_nueCC_from_intnue;
  osc_test->flag_BNB_nueCC_from_overlaynumu  = Configure_Osc::flag_BNB_nueCC_from_overlaynumu;
  osc_test->flag_BNB_nueCC_from_appnue       = Configure_Osc::flag_BNB_nueCC_from_appnue;
  osc_test->flag_BNB_nueCC_from_appnumu      = Configure_Osc::flag_BNB_nueCC_from_appnumu;
  osc_test->flag_BNB_nueCC_from_overlaynueNC = Configure_Osc::flag_BNB_nueCC_from_overlaynueNC;
  osc_test->flag_BNB_nueCC_from_overlaynumuNC= Configure_Osc::flag_BNB_nueCC_from_overlaynumuNC;
  
  osc_test->flag_BNB_numuCC_from_overlaynumu   = Configure_Osc::flag_BNB_numuCC_from_overlaynumu;
  osc_test->flag_BNB_numuCC_from_overlaynue    = Configure_Osc::flag_BNB_numuCC_from_overlaynue;
  osc_test->flag_BNB_numuCC_from_appnue        = Configure_Osc::flag_BNB_numuCC_from_appnue;
  osc_test->flag_BNB_numuCC_from_appnumu       = Configure_Osc::flag_BNB_numuCC_from_appnumu;
  osc_test->flag_BNB_numuCC_from_overlaynumuNC = Configure_Osc::flag_BNB_numuCC_from_overlaynumuNC;
  osc_test->flag_BNB_numuCC_from_overlaynueNC  = Configure_Osc::flag_BNB_numuCC_from_overlaynueNC;  
  
  osc_test->flag_BNB_CCpi0_from_overlaynumu  = Configure_Osc::flag_BNB_CCpi0_from_overlaynumu;
  osc_test->flag_BNB_CCpi0_from_appnue       = Configure_Osc::flag_BNB_CCpi0_from_appnue;
  osc_test->flag_BNB_CCpi0_from_overlaynumuNC= Configure_Osc::flag_BNB_CCpi0_from_overlaynumuNC;
  osc_test->flag_BNB_CCpi0_from_overlaynueNC = Configure_Osc::flag_BNB_CCpi0_from_overlaynueNC;
  
  osc_test->flag_BNB_NCpi0_from_overlaynumu  = Configure_Osc::flag_BNB_NCpi0_from_overlaynumu;
  osc_test->flag_BNB_NCpi0_from_appnue       = Configure_Osc::flag_BNB_NCpi0_from_appnue;
  osc_test->flag_BNB_NCpi0_from_overlaynumuNC= Configure_Osc::flag_BNB_NCpi0_from_overlaynumuNC;
  osc_test->flag_BNB_NCpi0_from_overlaynueNC = Configure_Osc::flag_BNB_NCpi0_from_overlaynueNC;
  
  /////// set only one time
  
  osc_test->Set_default_cv_cov(Configure_Osc::default_cv_file,
			       Configure_Osc::default_dirtadd_file,
			       Configure_Osc::default_mcstat_file,
			       Configure_Osc::default_fluxXs_dir,
			       Configure_Osc::default_detector_dir);
  
  osc_test->Set_oscillation_base();
  
  /////// Set_oscillation_pars(double val_dm2_41, double val_sin2_2theta_14, double val_sin2_theta_24, double val_sin2_theta_34)
  
  double val_dm2_41         = 7.3;
  double val_sin2_2theta_14 = 0.36;
  double val_sin2_theta_24  = 0;
  double val_sin2_theta_34  = 0;

  /// standard order
  val_dm2_41         = 7.3;
  val_sin2_2theta_14 = 0.236;
  val_sin2_theta_24  = 0.005;  
  osc_test->Set_oscillation_pars(val_dm2_41, val_sin2_2theta_14, val_sin2_theta_24, val_sin2_theta_34);  
  osc_test->Apply_oscillation();
  osc_test->Set_apply_POT();// meas, CV, COV: all ready
  //osc_test->Set_meas2fitdata();
  osc_test->Set_asimov2fitdata();
  
  ///////
  //osc_test->Plot_user();
  //osc_test->Minimization_OscPars_FullCov(7.3, 0.1, 0.0048, 0, "str_flag_fixpar");
  
  /////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////// Test speed and memory
  /////////////////////////////////////////////////////////// Test speed and memory

  if( 1 ) {
    cout<<endl<<" TEST "<<endl<<endl;

    //auto time_start = chrono::high_resolution_clock::now();
    //auto time_stop = chrono::high_resolution_clock::now();
    //auto time_duration = chrono::duration_cast<chrono::seconds>(time_stop - time_start);
    //cout<<endl<<" ---> check time duration "<<time_duration.count()<<endl<<endl;
    //milliseconds, minutes

    val_dm2_41         = 2;
    val_sin2_2theta_14 = 0.3;
    val_sin2_theta_24  = 0.05;  
    osc_test->Set_oscillation_pars(val_dm2_41, val_sin2_2theta_14, val_sin2_theta_24, val_sin2_theta_34);  
    osc_test->Apply_oscillation();
    osc_test->Set_apply_POT();// meas, CV, COV: all ready  
    osc_test->Set_meas2fitdata();

    for(int idx=1; idx<=100; idx++) {
      
      auto time_start = chrono::high_resolution_clock::now();
      for(int jdx=1; jdx<=100; jdx++) {
	//cout<<jdx<<endl;
	
	double pars_4v[4] = {2, 0.3, 0.05, 0};// dm2, sin^2_theta_14, sin^2_theta_24
	double chi2_4v = osc_test->FCN( pars_4v );
      }
      auto time_stop = chrono::high_resolution_clock::now();
      auto time_duration = chrono::duration_cast<chrono::milliseconds>(time_stop - time_start);
      cout<<" ---> time spent/100 calculation: "<<time_duration.count()<<endl;
    }
    
  }
   
  /////////////////////////////////////////////////////////// Profiling at each grid of (dm2, sin2_2t14), plus solution of T14, aaa
  /////////////////////////////////////////////////////////// Profiling at each grid of (dm2, sin2_2t14), plus solution of T14, aaa

  if( 0 ) {
    cout<<endl<<" ---> Profiling at each grid"<<endl;
 
    int bins_theta = 18;
    int bins_dm2   = 54;  
    TH2D *h2_space = new TH2D("h2_space_whole", "h2_space_whole", bins_theta, -1.2, 0, bins_dm2, -0.7, 2);

    ///////
    val_dm2_41         = 0;
    val_sin2_2theta_14 = 0.1;
    val_sin2_theta_24  = 0.01;  
    osc_test->Set_oscillation_pars(val_dm2_41, val_sin2_2theta_14, val_sin2_theta_24, val_sin2_theta_34);  
    osc_test->Apply_oscillation();
    osc_test->Set_apply_POT();// meas, CV, COV: all ready  
    //osc_test->Set_meas2fitdata();
    osc_test->Set_asimov2fitdata();
    
    double pars_4v[4] = {0};

    ///////
    
    int    min_status             = 10;
    int    usr_idm2               = 0;
    int    usr_it14               = 0;
    double min_chi2               = 0;
    double min_dm2_41_val         = 0;
    double min_sin2_2theta_14_val = 0;
    double min_sin2_theta_24_val  = 0;
    double min_sin2_theta_34_val  = 0;
    double min_dm2_41_err         = 0;
    double min_sin2_2theta_14_err = 0;
    double min_sin2_theta_24_err  = 0;
    double min_sin2_theta_34_err  = 0;
    
    roostr = TString::Format("sub_fit_%02d_%02d.root", idm2, it14);
    TFile *subroofile = new TFile(roostr, "recreate");
    TTree *tree = new TTree("tree", "tree");
    
    tree->Branch( "min_status",             &min_status,             "min_status/I" );
    tree->Branch( "usr_idm2",               &usr_idm2,               "usr_idm2/I" );
    tree->Branch( "usr_it14",               &usr_it14,               "usr_it14/I" );    
    tree->Branch( "min_chi2",               &min_chi2,               "min_chi2/D" );
    tree->Branch( "min_dm2_41_val",         &min_dm2_41_val,         "min_dm2_41_val/D" );
    tree->Branch( "min_sin2_2theta_14_val", &min_sin2_2theta_14_val, "min_sin2_2theta_14_val/D" );
    tree->Branch( "min_sin2_theta_24_val",  &min_sin2_theta_24_val,  "min_sin2_theta_24_val/D" );
    tree->Branch( "min_sin2_theta_34_val",  &min_sin2_theta_34_val,  "min_sin2_theta_34_val/D" );
    tree->Branch( "min_dm2_41_err",         &min_dm2_41_err,         "min_dm2_41_err/D" );
    tree->Branch( "min_sin2_2theta_14_err", &min_sin2_2theta_14_err, "min_sin2_2theta_14_err/D" );
    tree->Branch( "min_sin2_theta_24_err",  &min_sin2_theta_24_err,  "min_sin2_theta_24_err/D" );
    tree->Branch( "min_sin2_theta_34_err",  &min_sin2_theta_34_err,  "min_sin2_theta_34_err/D" );
    
    ///////

    double array_t24[300] = {0};
    for(int ii=0; ii<=100; ii++) {
      array_t24[ii] = ii*0.0001;
    }
    for(int ii=1; ii<200; ii++) {
      array_t24[ii+100] = ii*0.005;
    }
    
    for(int jdm2=idm2; jdm2<=idm2; jdm2++) {      
      for(int jt14=it14; jt14<=it14; jt14++) {
	cout<<" ------> jdm2 "<<jdm2<<" jt14 "<<jt14<<endl;

	double obj_chi2 = 1e10;
	double obj_dm2 = pow(10, h2_space->GetYaxis()->GetBinCenter(jdm2));
	double obj_t14 = pow(10, h2_space->GetXaxis()->GetBinCenter(jt14));
	double obj_t24 = 0;

	
	double y = obj_t14;// sin2_2t14
	double x = ( 1 + sqrt(1-y) )/2;
	double plus_sin2_T14 = x;	      
	obj_t14 = plus_sin2_T14;

	
	for(int idx=0; idx<300; idx++) {
	  if(idx%10==0) cout<<"       ---> processing t24idx "<<idx<<endl;
	  
	  //double val_t24 = idx * 0.001;
	  double val_t24 = array_t24[idx];
	  double val_chi2 = 0;
	  
	  if( val_t24<=1 ) {	    
	    pars_4v[0] = obj_dm2;
	    pars_4v[1] = obj_t14;	
	    pars_4v[2] = val_t24;
	    val_chi2 = osc_test->FCN( pars_4v );
	    
	    if( val_chi2 < obj_chi2 ) {
	      obj_chi2 = val_chi2;
	      obj_t24 = val_t24;
	    }
	  }	  
	  
	}// for(int idx=0; idx<100; jdx++)

	if( obj_t24==0 ) obj_t24 = 1e-8;
	cout<<endl<<" ---> objtesta "<<obj_t24<<endl;
	
	osc_test->Minimization_OscPars_FullCov(obj_dm2, obj_t14, obj_t24, 0, "dm2_t14");

	/////////// patch
	/////////// patch
	
	if( osc_test->minimization_status!=0 ) {
	  if( obj_t24<1e-7 ) {
	    cout<<endl<<" ---> pactch_AA_minimization "<<endl<<endl;
	    osc_test->Minimization_OscPars_FullCov(obj_dm2, obj_t14, 1e-7, 0, "dm2_t14");
	  }
	}
	if( osc_test->minimization_status!=0 ) {
	  if( obj_t24<1e-7 ) {
	    cout<<endl<<" ---> pactch_AB_minimization "<<endl<<endl;
	    osc_test->Minimization_OscPars_FullCov(obj_dm2, obj_t14, 1e-6, 0, "dm2_t14");
	  }
	}
	if( osc_test->minimization_status!=0 ) {
	  if( obj_t24<1e-7 ) {
	    cout<<endl<<" ---> pactch_AC_minimization "<<endl<<endl;
	    osc_test->Minimization_OscPars_FullCov(obj_dm2, obj_t14, 1e-5, 0, "dm2_t14");
	  }
	}


	if( osc_test->minimization_status!=0 ) {
	  cout<<endl<<" ---> pactch_BA_minimization "<<endl<<endl;
	  osc_test->Minimization_OscPars_FullCov(obj_dm2, obj_t14, osc_test->minimization_sin2_theta_24_val, 0, "dm2_t14");
	}
	if( osc_test->minimization_status!=0 ) {
	  cout<<endl<<" ---> pactch_BB_minimization "<<endl<<endl;
	  osc_test->Minimization_OscPars_FullCov(obj_dm2, obj_t14, osc_test->minimization_sin2_theta_24_val * 1.1, 0, "dm2_t14");
	}
	
	///////////
	///////////
	
	min_status            = osc_test->minimization_status;
	usr_idm2              = idm2;
	usr_it14              = it14;
	min_chi2              = osc_test->minimization_chi2;
	min_dm2_41_val        = osc_test->minimization_dm2_41_val;
	min_sin2_2theta_14_val= osc_test->minimization_sin2_2theta_14_val;
	min_sin2_theta_24_val = osc_test->minimization_sin2_theta_24_val;
	if(min_sin2_theta_24_val<1e-5) min_sin2_theta_24_val = 0;
	min_sin2_theta_34_val = osc_test->minimization_sin2_theta_34_val;
	min_dm2_41_err        = osc_test->minimization_dm2_41_err;
	min_sin2_2theta_14_err= osc_test->minimization_sin2_2theta_14_err;
	min_sin2_theta_24_err = osc_test->minimization_sin2_theta_24_err;
	min_sin2_theta_34_err = osc_test->minimization_sin2_theta_34_err;
	tree->Fill();
	
      }// for(int jt14=1; jt14<=bins_theta; jt14++)
    }// for(int jdm2=1; jdm2<=bins_dm2; jdm2++)

    tree->Write();
    subroofile->Close();
  }

             
  /////////////////////////////////////////////////////////// frequentist CLs: 3+1 analysis, dm2 vs. sin2_2Tee fixing t24 at the "best-fit of each grid"
  /////////////////////////////////////////////////////////// frequentist CLs: 3+1 analysis, dm2 vs. sin2_2Tee fixing t24 at the "best-fit of each grid"
  
  if( 0 ) {
    
    cout<<endl;
    cout<<" ---> frequentist CLs"<<endl;
   
    const int bins_theta = 18;
    const int bins_dm2   = 54;        
    TH2D *h2_space = new TH2D("h2_space_whole", "h2_space_whole", bins_theta, -1.2, 0, bins_dm2, -0.7, 2);
  
    ///////
    ///////
    
    double dm2_val = pow(10, h2_space->GetYaxis()->GetBinCenter(idm2));
    double t14_val = pow(10, h2_space->GetXaxis()->GetBinCenter(it14));
    double t24_val = 0;
    double obj_t14 = 0;
    
    double array_t24[bins_dm2+1][bins_theta+1] = {{0}};
    double array_obj_t14[bins_dm2+1][bins_theta+1] = {{0}};
    double array_chi2[bins_dm2+1][bins_theta+1] = {{0}};
    
    if( 1 ) {
      TFile *userfile = new TFile("./solution_plus.root", "read");
      TTree *tree_plus = (TTree*)userfile->Get("tree");

      // Declaration of leaf types
      Int_t           min_status;
      Int_t           usr_idm2;
      Int_t           usr_it14;
      Double_t        min_chi2;
      Double_t        min_dm2_41_val;
      Double_t        min_sin2_2theta_14_val;
      Double_t        min_sin2_theta_24_val;

      // List of branches
      TBranch        *b_min_status;   //!
      TBranch        *b_usr_idm2;   //!
      TBranch        *b_usr_it14;   //!
      TBranch        *b_min_chi2;   //!
      TBranch        *b_min_dm2_41_val;   //!
      TBranch        *b_min_sin2_2theta_14_val;   //!
      TBranch        *b_min_sin2_theta_24_val;   //!

      // Set branch addresses and branch pointers 
      tree_plus->SetBranchAddress("min_status", &min_status, &b_min_status);
      tree_plus->SetBranchAddress("usr_idm2", &usr_idm2, &b_usr_idm2);
      tree_plus->SetBranchAddress("usr_it14", &usr_it14, &b_usr_it14);
      tree_plus->SetBranchAddress("min_chi2", &min_chi2, &b_min_chi2);
      tree_plus->SetBranchAddress("min_dm2_41_val", &min_dm2_41_val, &b_min_dm2_41_val);
      tree_plus->SetBranchAddress("min_sin2_2theta_14_val", &min_sin2_2theta_14_val, &b_min_sin2_2theta_14_val);
      tree_plus->SetBranchAddress("min_sin2_theta_24_val", &min_sin2_theta_24_val, &b_min_sin2_theta_24_val);
      
      //
      int entries = tree_plus->GetEntries();
      bool flag_tree_plus = false;
      for(int ientry=0; ientry<entries; ientry++) {
	tree_plus->GetEntry(ientry);
	if( usr_idm2==idm2 && usr_it14==it14) {	  

	  double plus_chi2      = min_chi2;
	  double plus_sin2_T14  = min_sin2_2theta_14_val;
	  double plus_sin2_T24  = min_sin2_theta_24_val;

	  t24_val = plus_sin2_T24;
	  obj_t14 = plus_sin2_T14;
	  
	  flag_tree_plus = true;
	  break;
	}
      }

      for(int ientry=0; ientry<entries; ientry++) {
	tree_plus->GetEntry(ientry);


	double plus_chi2      = min_chi2;
	double plus_sin2_T14  = min_sin2_2theta_14_val;
	double plus_sin2_T24  = min_sin2_theta_24_val;
	
	array_t24[usr_idm2][usr_it14] = min_sin2_theta_24_val;
	array_obj_t14[usr_idm2][usr_it14] = plus_sin2_T14;
	array_chi2[usr_idm2][usr_it14] = min_chi2;
      }   
      
      delete tree_plus;
      delete userfile;

      if( !flag_tree_plus ) {
	cout<<endl;
	cout<<" ERROR: do not find the grid minimum of t24 at idm2/it14: "<<idm2<<"\t"<<it14<<endl;
	cout<<endl;
      }
    }

    //cout<<endl<<" ---> gridmin t24: "<<idm2<<"\t"<<it14<<"\t"<<t24_val<<endl<<endl;
    //return 0;
    
    double pars_4v[4] = {dm2_val, obj_t14, t24_val, 0};
    double pars_3v[4] = {0};

    ///////////////////////////////////////////////////////////
    
    if( 1 ) {
      cout<<" ---> Asimov/Pseudo CLs "<<endl;

      vector<int> vec_dm2;
      vector<int> vec_t14;      
      vector<double> vec_chi2_4v;
      vector<double> vec_chi2_3v;
      
      roostr = TString::Format("sub_sm_CLs_%06d.root", ifile);
      TFile *subroofile = new TFile(roostr, "recreate");
      TTree *tree = new TTree("tree", "tree");
      
      tree->Branch( "vec_dm2",  &vec_dm2 );
      tree->Branch( "vec_t14",  &vec_t14 );
      tree->Branch( "vec_chi2_4v",  &vec_chi2_4v );
      tree->Branch( "vec_chi2_3v",  &vec_chi2_3v );

      //////// 3v
      osc_test->Set_oscillation_pars(0, 0, 0, 0);  
      osc_test->Apply_oscillation();
      osc_test->Set_apply_POT();// meas, CV, COV: all ready
      
      //osc_test->Set_meas2fitdata();
      osc_test->Set_asimov2fitdata();      
      //osc_test->Set_toy_variations( 1 );
      //osc_test->Set_toy2fitdata( 1 );
      
      for(int jdm2=1; jdm2<=bins_dm2; jdm2++) {
	cout<<" ---> processing jdm2 "<<jdm2<<endl;
	
	for(int jt14=1; jt14<=bins_theta; jt14++) {
	  //cout<<" ---> processing tt "<<jt14<<endl;
	  
	  double dm2 = pow(10, h2_space->GetYaxis()->GetBinCenter(jdm2));
	  double t14 = array_obj_t14[jdm2][jt14];
	  double t24 = array_t24[jdm2][jt14];
	  pars_4v[0] = dm2;
	  pars_4v[1] = t14;
	  pars_4v[2] = t24;

	  double chi2_4v = osc_test->FCN( pars_4v );
	  double chi2_3v = osc_test->FCN( pars_3v );
	  
	  if((t24!=0) && 0) {
	    cout<<" nozero: "<<jdm2<<"\t"<<t14<<"\t t24: "<<t24<<"\t min_chi2: "
		<<array_chi2[jdm2][jt14]<<"\t"<<chi2_4v<<endl;
	  }
	  
	  vec_dm2.push_back(jdm2);
	  vec_t14.push_back(jt14);
	  vec_chi2_4v.push_back( chi2_4v );
	  vec_chi2_3v.push_back( chi2_3v );
	}
      }
      
      tree->Fill();

      
      tree->Write();
      subroofile->Close();
	
      return 0;
    }
    
    /////////////////////////////////////////////////////

    
    int ntoys = 2000;

    
    vector<double> vec_chi2_4v_on_4vPseudo;
    vector<double> vec_chi2_3v_on_4vPseudo;      
    vector<double> vec_chi2_4v_on_3vPseudo;
    vector<double> vec_chi2_3v_on_3vPseudo;
      
    roostr = TString::Format("sub_frequentist_CLs_dm2_%02d_tt_%02d.root", idm2, it14);
    TFile *subroofile = new TFile(roostr, "recreate");
    TTree *tree = new TTree("tree", "tree");
    
    tree->Branch( "idm2",          &idm2,     "idm2/I" );
    tree->Branch( "it14",          &it14,     "it14/I" );
    tree->Branch( "vec_chi2_4v_on_4vPseudo",  &vec_chi2_4v_on_4vPseudo );
    tree->Branch( "vec_chi2_3v_on_4vPseudo",  &vec_chi2_3v_on_4vPseudo );
    tree->Branch( "vec_chi2_4v_on_3vPseudo",  &vec_chi2_4v_on_3vPseudo );
    tree->Branch( "vec_chi2_3v_on_3vPseudo",  &vec_chi2_3v_on_3vPseudo );

    //////// 4v
    osc_test->Set_oscillation_pars(pars_4v[0], pars_4v[1], pars_4v[2], pars_4v[3]);  
    osc_test->Apply_oscillation();
    osc_test->Set_apply_POT();// meas, CV, COV: all ready
    TMatrixD matrix_4v_COV  = osc_test->matrix_eff_newworld_abs_syst_total;
    TMatrixD matrix_4v_pred = osc_test->matrix_eff_newworld_pred;
    
    //////// 3v
    osc_test->Set_oscillation_pars(0, 0, 0, 0);  
    osc_test->Apply_oscillation();
    osc_test->Set_apply_POT();// meas, CV, COV: all ready
    TMatrixD matrix_3v_COV  = osc_test->matrix_eff_newworld_abs_syst_total;
    TMatrixD matrix_3v_pred = osc_test->matrix_eff_newworld_pred;
  
    ////////
    
    for(int itoy=1; itoy<=ntoys; itoy++) {

      cout<<" ---> processing itoy "<<itoy<<endl;
      
      /////////////////////// 4v pseudo
    
      // osc_test->Set_oscillation_pars(dm2_val, t14_val, 0, 0);  
      // osc_test->Apply_oscillation();
      // osc_test->Set_apply_POT();// meas, CV, COV: all ready
      osc_test->matrix_eff_newworld_abs_syst_total = matrix_4v_COV;
      osc_test->matrix_eff_newworld_pred = matrix_4v_pred;
      osc_test->Set_toy_variations( 1 );
      osc_test->Set_toy2fitdata( 1 );
      double chi2_4v_on_4vPseudo = osc_test->FCN( pars_4v );
      double chi2_3v_on_4vPseudo = osc_test->FCN( pars_3v );

      vec_chi2_4v_on_4vPseudo.push_back(chi2_4v_on_4vPseudo);
      vec_chi2_3v_on_4vPseudo.push_back(chi2_3v_on_4vPseudo);
      
      /////////////////////// 3v pseudo
    
      // osc_test->Set_oscillation_pars(0, 0, 0, 0);  
      // osc_test->Apply_oscillation();
      // osc_test->Set_apply_POT();// meas, CV, COV: all ready
      osc_test->matrix_eff_newworld_abs_syst_total = matrix_3v_COV;
      osc_test->matrix_eff_newworld_pred = matrix_3v_pred;
      osc_test->Set_toy_variations( 1 );
      osc_test->Set_toy2fitdata( 1 );      
      double chi2_4v_on_3vPseudo = osc_test->FCN( pars_4v );
      double chi2_3v_on_3vPseudo = osc_test->FCN( pars_3v );

      vec_chi2_4v_on_3vPseudo.push_back(chi2_4v_on_3vPseudo);
      vec_chi2_3v_on_3vPseudo.push_back(chi2_3v_on_3vPseudo);

      ////////////////////////

      //cout<<"  ---> "<<chi2_4v_on_3vPseudo - chi2_3v_on_3vPseudo <<endl;
      //cout<<"  ---> check1 "<<pars_4v[0]<<"\t"<<pars_4v[1]<<"\t"<<pars_4v[2]<<"\t"<<pars_4v[3]<<endl;
    }
    
    tree->Fill();
    tree->Write();
    subroofile->Close();

    delete h2_space;
    
  }
  


  /////////////////////////////////////////////////////////// frequentist CLs: 3+1 analysis, dm2 vs. sin2_2Tue fixing t24 at 0.005
  /////////////////////////////////////////////////////////// frequentist CLs: 3+1 analysis, dm2 vs. sin2_2Tue fixing t24 at 0.005
  
  if( 0 ) {
    
    cout<<endl;
    cout<<" ---> frequentist CLs"<<endl;
    
    const int bins_theta = 60;
    const int bins_dm2   = 60;
        
    /////// X: sin22t14, 1e-4 -> 1   ---> "log10()" ---> -4 -> 0
    /////// Y: m41^2,    1e-1 -> 1e2  ---> "log10()" ---> -1 -> 2
    TH2D *h2_space = new TH2D("h2_space_whole", "h2_space_whole", bins_theta, -4, 0, bins_dm2, -1, 2);

    double dm2_val = pow(10, h2_space->GetYaxis()->GetBinCenter(idm2));
    double t14_val = pow(10, h2_space->GetXaxis()->GetBinCenter(it14));
    double t24_val = 0;
    double obj_t14 = 0;
    
    double array_t24[bins_dm2+1][bins_theta+1] = {{0}};
    double array_obj_t14[bins_dm2+1][bins_theta+1] = {{0}};
    double array_chi2[bins_dm2+1][bins_theta+1] = {{0}};
    
    /////////////////
    
    double plus_sin2_2Tue = t14_val;
    double plus_sin2_T24  = 0.005;
    double plus_sin2_2Tee = plus_sin2_2Tue/plus_sin2_T24;
    double y = plus_sin2_2Tee;
    double x = ( 1 + sqrt(1-y) )/2;
    double plus_sin2_T14 = x;	      
    obj_t14 = plus_sin2_T14;

    cout<<" ---> check "<<t14_val<<"\t"<<4*obj_t14*(1-obj_t14)*0.005<<endl;

    //return 0;
    
    //////////////////

    //cout<<endl<<" ---> gridmin t24: "<<idm2<<"\t"<<it14<<"\t"<<t24_val<<endl<<endl;
    //return 0;
    
    double pars_4v[4] = {dm2_val, obj_t14, 0.005, 0};
    double pars_3v[4] = {0};

    ///////////////////////////////////////////////////////////
    
    if( 0 ) {
      cout<<" ---> Asimov/Pseudo CLs "<<endl;

      vector<int> vec_dm2;
      vector<int> vec_t14;      
      vector<double> vec_chi2_4v;
      vector<double> vec_chi2_3v;
      
      roostr = TString::Format("sub_sm_CLs_%06d.root", ifile);
      TFile *subroofile = new TFile(roostr, "recreate");
      TTree *tree = new TTree("tree", "tree");
      
      tree->Branch( "vec_dm2",  &vec_dm2 );
      tree->Branch( "vec_t14",  &vec_t14 );
      tree->Branch( "vec_chi2_4v",  &vec_chi2_4v );
      tree->Branch( "vec_chi2_3v",  &vec_chi2_3v );

      //////// 3v
      osc_test->Set_oscillation_pars(0, 0, 0, 0);  
      osc_test->Apply_oscillation();
      osc_test->Set_apply_POT();// meas, CV, COV: all ready
      
      //osc_test->Set_meas2fitdata();
      //osc_test->Set_asimov2fitdata();      
      osc_test->Set_toy_variations( 1 );
      osc_test->Set_toy2fitdata( 1 );

      for(int jdm2=1; jdm2<=bins_dm2; jdm2++) {
	cout<<" ---> processing jdm2 "<<jdm2<<endl;
	
	for(int jt14=1; jt14<=bins_theta; jt14++) {
	  //cout<<" ---> processing tt "<<jt14<<endl;
	  
	  double dm2 = pow(10, h2_space->GetYaxis()->GetBinCenter(jdm2));
	  double t14 = pow(10, h2_space->GetXaxis()->GetBinCenter(jt14));
	  double t24 = 0.005;

	  double plus_sin2_2Tue = t14;
	  double plus_sin2_T24  = t24;
	  double plus_sin2_2Tee = plus_sin2_2Tue/plus_sin2_T24;
	  double y = plus_sin2_2Tee;
	  double x = ( 1 + sqrt(1-y) )/2;
	  double plus_sin2_T14 = x;	      
	  double usr_t14 = plus_sin2_T14;

	  //cout<<" ---> check "<<jt14<<"\t"<<t14<<"\t"<<4*usr_t14*(1-usr_t14)*t24<<endl;
	  
	  if( jt14<26 ) {
	    pars_4v[0] = dm2;
	    pars_4v[1] = usr_t14;
	    pars_4v[2] = t24;
	    
	    double chi2_4v = osc_test->FCN( pars_4v );
	    double chi2_3v = osc_test->FCN( pars_3v );	  
	    
	    vec_dm2.push_back(jdm2);
	    vec_t14.push_back(jt14);
	    vec_chi2_4v.push_back( chi2_4v );
	    vec_chi2_3v.push_back( chi2_3v );
	  }
	  else {
	    vec_dm2.push_back(jdm2);
	    vec_t14.push_back(jt14);
	    vec_chi2_4v.push_back( 0 );
	    vec_chi2_3v.push_back( 0 );
	  }
	}
      }
      
      tree->Fill();

      
      tree->Write();
      subroofile->Close();
	
      return 0;
    }
    
    /////////////////////////////////////////////////////

    
    int ntoys = 1000;

    
    vector<double> vec_chi2_4v_on_4vPseudo;
    vector<double> vec_chi2_3v_on_4vPseudo;      
    vector<double> vec_chi2_4v_on_3vPseudo;
    vector<double> vec_chi2_3v_on_3vPseudo;
      
    roostr = TString::Format("sub_frequentist_CLs_dm2_%02d_tt_%02d.root", idm2, it14);
    TFile *subroofile = new TFile(roostr, "recreate");
    TTree *tree = new TTree("tree", "tree");
    
    tree->Branch( "idm2",          &idm2,     "idm2/I" );
    tree->Branch( "it14",          &it14,     "it14/I" );
    tree->Branch( "vec_chi2_4v_on_4vPseudo",  &vec_chi2_4v_on_4vPseudo );
    tree->Branch( "vec_chi2_3v_on_4vPseudo",  &vec_chi2_3v_on_4vPseudo );
    tree->Branch( "vec_chi2_4v_on_3vPseudo",  &vec_chi2_4v_on_3vPseudo );
    tree->Branch( "vec_chi2_3v_on_3vPseudo",  &vec_chi2_3v_on_3vPseudo );

    //////// 4v
    osc_test->Set_oscillation_pars(pars_4v[0], pars_4v[1], pars_4v[2], pars_4v[3]);  
    osc_test->Apply_oscillation();
    osc_test->Set_apply_POT();// meas, CV, COV: all ready
    TMatrixD matrix_4v_COV  = osc_test->matrix_eff_newworld_abs_syst_total;
    TMatrixD matrix_4v_pred = osc_test->matrix_eff_newworld_pred;
    
    //////// 3v
    osc_test->Set_oscillation_pars(0, 0, 0, 0);  
    osc_test->Apply_oscillation();
    osc_test->Set_apply_POT();// meas, CV, COV: all ready
    TMatrixD matrix_3v_COV  = osc_test->matrix_eff_newworld_abs_syst_total;
    TMatrixD matrix_3v_pred = osc_test->matrix_eff_newworld_pred;
  
    ////////
    
    for(int itoy=1; itoy<=ntoys; itoy++) {

      cout<<" ---> processing itoy "<<itoy<<endl;

      if( it14<26 ) {      
	/////////////////////// 4v pseudo
    
	// osc_test->Set_oscillation_pars(dm2_val, t14_val, 0, 0);  
	// osc_test->Apply_oscillation();
	// osc_test->Set_apply_POT();// meas, CV, COV: all ready
	osc_test->matrix_eff_newworld_abs_syst_total = matrix_4v_COV;
	osc_test->matrix_eff_newworld_pred = matrix_4v_pred;
	osc_test->Set_toy_variations( 1 );
	osc_test->Set_toy2fitdata( 1 );
	double chi2_4v_on_4vPseudo = osc_test->FCN( pars_4v );
	double chi2_3v_on_4vPseudo = osc_test->FCN( pars_3v );

	vec_chi2_4v_on_4vPseudo.push_back(chi2_4v_on_4vPseudo);
	vec_chi2_3v_on_4vPseudo.push_back(chi2_3v_on_4vPseudo);
      
	/////////////////////// 3v pseudo
    
	// osc_test->Set_oscillation_pars(0, 0, 0, 0);  
	// osc_test->Apply_oscillation();
	// osc_test->Set_apply_POT();// meas, CV, COV: all ready
	osc_test->matrix_eff_newworld_abs_syst_total = matrix_3v_COV;
	osc_test->matrix_eff_newworld_pred = matrix_3v_pred;
	osc_test->Set_toy_variations( 1 );
	osc_test->Set_toy2fitdata( 1 );      
	double chi2_4v_on_3vPseudo = osc_test->FCN( pars_4v );
	double chi2_3v_on_3vPseudo = osc_test->FCN( pars_3v );

	vec_chi2_4v_on_3vPseudo.push_back(chi2_4v_on_3vPseudo);
	vec_chi2_3v_on_3vPseudo.push_back(chi2_3v_on_3vPseudo);
      }
      else {
	vec_chi2_4v_on_4vPseudo.push_back(0);
	vec_chi2_3v_on_4vPseudo.push_back(0);
	vec_chi2_4v_on_3vPseudo.push_back(0);
	vec_chi2_3v_on_3vPseudo.push_back(0);	
      }
      
    }
    
    tree->Fill();
    tree->Write();
    subroofile->Close();

    delete h2_space;
    
  }

  ///////////////////////////////////////////////////////////

  cout<<endl;
  cout<<" ------------------------------ check at the final step ------------------------------"<<endl;
  cout<<" ------------------------------ check at the final step ------------------------------"<<endl;

  cout<<endl;
  cout<<TString::Format(" ---> display(-d) %d, ifile(-f) %d, scaleF_POT_BNB(-pbnb) %6.4f, scaleF_POT_NuMI(-pnumi) %6.4f, theta14(-it14) %d, dm2(-idm2) %d, theta24(-it24) %d",
			display, ifile, osc_test->scaleF_POT_BNB, osc_test->scaleF_POT_NuMI, it14, idm2, it24)<<endl;
  
  cout<<endl;
  cout<<TString::Format(" ---> flag_syst_dirt    %d", osc_test->flag_syst_dirt)<<endl;
  cout<<TString::Format(" ---> flag_syst_mcstat  %d", osc_test->flag_syst_mcstat)<<endl;
  cout<<TString::Format(" ---> flag_syst_flux    %d", osc_test->flag_syst_flux)<<endl;
  cout<<TString::Format(" ---> flag_syst_geant   %d", osc_test->flag_syst_geant)<<endl;
  cout<<TString::Format(" ---> flag_syst_Xs      %d", osc_test->flag_syst_Xs)<<endl;
  cout<<TString::Format(" ---> flag_syst_det     %d", osc_test->flag_syst_det)<<endl;

  cout<<endl;
  cout<<" ---> Finished sucessfully"<<endl;
  
  cout<<endl;
  if( display ) {
    cout<<" Enter Ctrl+c to end the program"<<endl;
    cout<<" Enter Ctrl+c to end the program"<<endl;
    cout<<endl;
    theApp.Run();
  }
  
  return 0;
}
