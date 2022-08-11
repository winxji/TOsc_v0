void read_obj()
{
  TString roostr = "";

  TString map_detectorfile_str[11];
  map_detectorfile_str[1] = "cov_LYDown.root";
  map_detectorfile_str[2] = "cov_LYRayleigh.root";
  map_detectorfile_str[3] = "cov_Recomb2.root";
  map_detectorfile_str[4] = "cov_SCE.root";
  map_detectorfile_str[5] = "cov_WMdEdx.root";
  map_detectorfile_str[6] = "cov_WMThetaXZ.root";
  map_detectorfile_str[7] = "cov_WMThetaYZ.root";
  map_detectorfile_str[8] = "cov_WMX.root";
  map_detectorfile_str[9] = "cov_WMYZ.root";
  map_detectorfile_str[10]= "cov_LYatt.root";

  for(int idx=1; idx<=10; idx++) {
    if(idx==5) continue;

    TFile *file_input = new TFile("../DetVar/"+map_detectorfile_str[idx], "read");

    roostr = TString::Format("frac_cov_det_mat_%d", idx);
    TMatrixD *matrix_temp = (TMatrixD*)file_input->Get(roostr);

    int rows = matrix_temp->GetNrows();
    cout<<" ---> rows "<<rows<<endl;

    TMatrixD matrix_user(rows, rows);
    matrix_user = (*matrix_temp);
    
    // BNB and NuMI: 26 bins per channel
    //
    // 1-7   overlay
    // 8-14  EXT
    // 15-21 nueapp
    //
    // 22-28 overlay
    // 29-35 EXT
    // 36-42 nueapp
    //

    int baseAA = 0;
    int baseBB = 0;

    //////////// BNB and BNB
    ////////////
    
    baseAA = 0;
    baseBB = 26*14;

    for(int idx=1; idx<=26*7; idx++) {
      for(int jdx=1; jdx<=26*7; jdx++) {

	int i_AA = idx-1 + baseAA;
	int j_AA = jdx-1 + baseAA;	
	double cov_ij_AA = (*matrix_temp)(i_AA, j_AA);

	int i_BB = 0;
	int j_BB = 0;

	//////
	i_BB = idx-1 + baseBB;
	j_BB = jdx-1 + baseBB;
	matrix_user(i_BB, j_BB) = cov_ij_AA;
	
	//////
	i_BB = idx-1 + baseAA;
	j_BB = jdx-1 + baseBB;
	matrix_user(i_BB, j_BB) = cov_ij_AA;
		
	//////
	i_BB = idx-1 + baseBB;
	j_BB = jdx-1 + baseAA;
	matrix_user(i_BB, j_BB) = cov_ij_AA;
      }
    }

    //////////// NuMI and NuMI
    ////////////
    
    baseAA = 26*21;
    baseBB = 26*21 + 26*14;

    for(int idx=1; idx<=26*7; idx++) {
      for(int jdx=1; jdx<=26*7; jdx++) {

	int i_AA = idx-1 + baseAA;
	int j_AA = jdx-1 + baseAA;	
	double cov_ij_AA = (*matrix_temp)(i_AA, j_AA);

	int i_BB = 0;
	int j_BB = 0;

	//////
	i_BB = idx-1 + baseBB;
	j_BB = jdx-1 + baseBB;
	matrix_user(i_BB, j_BB) = cov_ij_AA;
	
	//////
	i_BB = idx-1 + baseAA;
	j_BB = jdx-1 + baseBB;
	matrix_user(i_BB, j_BB) = cov_ij_AA;
		
	//////
	i_BB = idx-1 + baseBB;
	j_BB = jdx-1 + baseAA;
	matrix_user(i_BB, j_BB) = cov_ij_AA;
      }
    }

    //////////// BNB and NuMI
    ////////////

    for(int idx=1; idx<=26*7; idx++) {
      for(int jdx=1; jdx<=26*7; jdx++) {

	int i_AA = idx-1;
	int j_AA = jdx-1 + 26*21;	
	double cov_ij_AA = (*matrix_temp)(i_AA, j_AA);

	int i_BB = 0;
	int j_BB = 0;

	//////////
	i_BB = idx-1;
	j_BB = jdx-1 + 26*21 + 26*14;
	matrix_user(i_BB, j_BB) = cov_ij_AA;
	matrix_user(j_BB, i_BB) = cov_ij_AA;

	//////////
	i_BB = idx-1 + 26*14;
	j_BB = jdx-1 + 26*21 + 26*14;
	matrix_user(i_BB, j_BB) = cov_ij_AA;
	matrix_user(j_BB, i_BB) = cov_ij_AA;
	
	//////////
	i_BB = idx-1 + 26*14;
	j_BB = jdx-1 + 26*21;
	matrix_user(i_BB, j_BB) = cov_ij_AA;
	matrix_user(j_BB, i_BB) = cov_ij_AA;
	
      }
    }
    
    
    //////////////////////////////
    
    delete file_input;

    TFile *outfile = new TFile(map_detectorfile_str[idx], "recreate");
    roostr = TString::Format("frac_cov_det_mat_%d", idx);
    matrix_user.Write(roostr);
    outfile->Close();
    
  }


  
}
