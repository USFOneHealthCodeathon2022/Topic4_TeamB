data {
  int N_unk_dis;
}
parameters {
  vector<lower=-1,upper=1>[N_unk_dis] disease_state;
}
model {

}
