void check()
{
  TString roostr = "";

  TFile *file = new TFile("merge.root", "read");
  for(int idx=1; idx<=1000; idx++) {
    roostr = TString::Format("histo_%d", idx);
    TH1D *h1d_temp = (TH1D*)file->Get(roostr);
    if( h1d_temp==NULL ) break;
    cout<<idx<<"\t"<<h1d_temp->Integral()<<endl;
    delete h1d_temp;
  }

}
