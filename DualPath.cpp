//#define ARMA_NO_DEBUG
#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS
#include<armadillo>
#include <time.h>
using namespace arma;

int main() {
  // Users need to specify two things:
  // 1/ The text file that contains the data
  arma::Mat<double> yv;
  yv.load("data.txt");
  // 2/ Number of lambdas along the solution path at which we solve the GFL problem
  int len=100;
  // The whole solution path for the GFL will be saved in "GFL_dualpath.txt"
  
  
  double eps = 0.0001;
  int checkk=0;
  int n = yv.n_rows;
  int p = yv.n_cols;
  int n1=n-1;
  int ite=10;
  vec theta = zeros(n);
  int i,l,j,k,g,it,la;
  int len=n;
  vec inside;  
  mat D(p*n1,p*n,fill::zeros);
  for(i = 0;i<p*n1;i=i+p){
    D(span(i,i+p-1),span(i,i+p-1))=eye(p,p);
    D(span(i,i+p-1),span(i+p,p+i+p-1))=-eye(p,p);
  }
  int Dnr = D.n_rows;
  vec my = sum(yv,1)/n;
  vec y = vectorise(yv);
  mat tempy(p*n,p,fill::zeros);
  for(i=0;i<p*n;i=i+p){
    tempy(span(i,i+p-1),span(0,p-1))=eye(p,p);
  }
  vec avg = tempy * my;
  mat tD = D.t();
  mat uni = eye(p*n,p*n);
  mat gtD = solve(D.t(),uni,solve_opts::fast);
  
  vec umax = gtD * (y-avg);
  vec umax2 = umax % umax;
  mat umaxm = reshape(umax2,p,n);
  rowvec rs = sqrt(sum(umaxm));
  uvec activeu = find(zeros<vec>(1)==1);
  
  double lambmax = max(rs);
  mat lambseq = linspace<vec>(0, lambmax,len);
  mat thetas(n*p,len,fill::zeros); 
  mat us(n1*p,len+1,fill::zeros); 
  
  mat precals = D/2;
  vec first = precals * y;
  vec u(Dnr+p,fill::zeros);
  uvec tempid;
  vec tempactiveu2;
  vec tempactiveu ; 
  uvec dti;
  vec newu =u;
  
  uvec lfinalind = find(linspace<vec>(0, p-1, p)==linspace<vec>(0, p-1, p));
  uvec finalind = lfinalind + p;
  uvec ufinalind = finalind+p;
  uvec tfinalind=finalind+p*(n-2);
  finalind=join_cols(finalind,tfinalind);
  
  
  for(i=p;i<Dnr;i=i+p){
    finalind = join_cols(finalind,lfinalind);
    finalind = join_cols(finalind,ufinalind);
    lfinalind = lfinalind+p;
    ufinalind = ufinalind+p;
    
  }
  vec minside;
  vec diff;
  vec dif1;
  uvec tempid11;
  vec diff11;
  uvec fdiff;
  uvec latesubsort;
  vec latediff;
  uvec ifd;
  int un=0;
  vec dif2;
  uvec lala;
  vec lightu;
  uvec combine;
  mat mdif1;
  uvec umdif1;
  rowvec tsum;
  vec temu;
  for(i=0;i<len;i++){
    vec sil = ones(2)*i+1;
    sil(1)=10;
    int ms = min(sil);
    int check=0;
    double lamb = lambseq(i,0);
    if(lamb==0){ theta = y;  u.subvec(0,Dnr-1)=zeros(umax.n_elem);}
    else if(lamb == lambmax){ theta = avg; u.subvec(0,Dnr-1)=umax;} 
    else {
      u.subvec(0,Dnr-1) = sign(us.unsafe_col(i)) * lamb;
      if(activeu.n_rows !=0){u.rows(activeu) = us.unsafe_col(i)(activeu)* 2 - us.unsafe_col(i-1)(activeu) ;}
      newu =u;
      tempactiveu = ones<vec>(Dnr) * (-3);
      int cou = 0;
      while(check==0){
        temu = u;
        if(cou>0){ms=2;}
        cou++;
        for(it = 0;it<ms; it++){
          for(j=0;j<Dnr;j=j+p){
            inside = first.subvec(j,j+p-1)-((u(finalind.subvec(2*j,2*j+p-1)) + u(finalind.subvec(2*j+p,2*j+2*p-1)))* (-0.5) );
            double n2 = sqrt(as_scalar(inside.t()*inside));
            if(n2>lamb){u.subvec(j,j+p-1) = lamb*inside/n2;}
            else{u.subvec(j,j+p-1) = inside; 
            }
          }
        }
        dif1 = u-temu;
        newu = u;
        
        if(max(abs(dif1))>eps/1000  ){
          int checkf=0;
          mdif1 = reshape(dif1,p,n);
          tsum=sum(abs(mdif1));
          umdif1 = find(tsum>eps*p/1000);
          diff11 = zeros(u.n_elem);
          while(checkf==0){
            for(g =0;g<umdif1.n_elem;g++){
              j=umdif1(g)*p;
              inside = first.subvec(j,j+p-1)-((u(finalind.subvec(2*j,2*j+p-1)) + u(finalind.subvec(2*j+p,2*j+2*p-1)))* (-0.5) );
              double n2 = sqrt(as_scalar(inside.t()*inside));
              if(n2>lamb){newu.subvec(j,j+p-1) = lamb*inside/n2;}
              else{newu.subvec(j,j+p-1) = inside; 
              }
              diff11.subvec(j,j+p-1) = newu.subvec(j,j+p-1)-u.subvec(j,j+p-1);
              u.subvec(j,j+p-1) = newu.subvec(j,j+p-1);
            }
            if((max(abs(diff11))<eps)){checkf = 1;}
          }
          if(umdif1.n_elem != n-1){
            tempid11 = find(tsum <= eps*p/1000 );
            diff = zeros(u.n_elem);
            for(g =0;g<tempid11.n_elem-1;g++){
              j=tempid11(g)*p;
              inside = first.subvec(j,j+p-1)-((u(finalind.subvec(2*j,2*j+p-1)) + u(finalind.subvec(2*j+p,2*j+2*p-1)))* (-0.5) );
              double n2 = sqrt(as_scalar(inside.t()*inside));
              if(n2>lamb){newu.subvec(j,j+p-1) = lamb*inside/n2;}
              else{newu.subvec(j,j+p-1) = inside; 
              }
              diff.subvec(j,j+p-1) = newu.subvec(j,j+p-1)-u.subvec(j,j+p-1);
              u.subvec(j,j+p-1) = newu.subvec(j,j+p-1);
            }
            if((max(abs(diff))<eps)){check = 1;}
          } else {check = 1;}
          
        } else {check = 1;}
        
      }
    }
    if((lamb>0) && (lamb<lambmax)){ theta = y-tD*u.subvec(0,Dnr-1);}
    us.unsafe_col(i+1) = u.subvec(0,Dnr-1);
    thetas.unsafe_col(i) = theta;
    uvec idu = find((round(abs(us.unsafe_col(i+1)) * 1000) / 1000) < (round(lamb * 1000) / 1000));
    activeu = idu;
  }
  thetas.save("GFL_dualpath.txt", arma_ascii);
 }
