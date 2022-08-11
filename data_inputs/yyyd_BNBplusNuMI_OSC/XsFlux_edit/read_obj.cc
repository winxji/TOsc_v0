void read_obj()
{
  TString roostr = "";

  int rows = 26*42;

  TMatrixD matrix_user[18];
 
  
  for(int idx=1; idx<=17; idx++) {
    cout<<" ---> "<<idx<<endl;
    
    matrix_user[idx].Clear(); matrix_user[idx].ResizeTo(rows, rows);
    
    roostr = TString::Format("../XsFlux/cov_%d.root", idx);
    TFile *roofile = new TFile(roostr, "read");

    TMatrixD *matrix_cov = (TMatrixD*)roofile->Get( TString::Format("frac_cov_xf_mat_%d", idx) );
    //cout<<matrix_cov->GetNrows()<<endl;
    
    for(int i=0; i<rows; i++) {
      for(int j=0; j<rows; j++) {
	double cov = (*matrix_cov)(i,j);

	bool ff_bnb = false;
	if(i<26*21 && j<26*21) ff_bnb = true;

	bool ff_numi = false;
	if(i>=26*21 && j>=26*21) ff_numi = true;

	///////////
	
	if( ff_bnb || ff_numi ) {
	  int itest = 0;
	}
	else {
	  cov = 0;
	}
	
	///////////
	
	if( ff_numi ) {
	  if( idx!=3) cov = 0;
	}

	///////////

	if( idx>=14 ) cov = (*matrix_cov)(i,j);
	
	///////////
	
	(matrix_user[idx])(i,j) = cov;
      }
    }
    
    delete roofile;
    
  }

  for(int ifile=1; ifile<=17; ifile++) {
    roostr = TString::Format("cov_%d.root", ifile);
    TFile *roofile = new TFile(roostr, "recreate");
    matrix_user[ifile].Write( TString::Format("frac_cov_xf_mat_%d", ifile) );
    delete roofile;
  }

}
