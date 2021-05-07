/*
 * funaux.h
 *
 *  Created on: 29 août 2016
 *      Author: anthony
 */

#ifndef C_FUNAUX_H_
#define C_FUNAUX_H_

#include "armadillo"
#include <string>
#include "lorenz95.h"

using namespace arma;
using namespace std;

template<class T>
T H(const T& X, string name){
	T Y(X);
	if(name.compare("id")!=0){
		cout << "Unknown observation operator" << endl;
	}
	return Y;
}

vec moy(mat E){
	return mean(E,1);
}
mat ano(mat E){
	int d = E.n_cols;
	vec un(d,fill::ones);
	return (E-moy(E)*un.t())*pow(d-1,-0.5);
}
mat ens(vec x, mat A){
	int d = A.n_cols;
	vec un(d,fill::ones);
	return x*un.t()+A*pow(d-1,0.5);
}
mat SR(mat A, int m){
	int Ns(A.n_rows);
	int n(min(Ns,m));
	mat U,V;vec s;
	svd(U,s,V,A);
	mat S(A.n_rows,m,fill::zeros);
	S.submat(0,0,n-1,n-1) = diagmat(s.subvec(0,n-1));
	vec a(m,fill::ones);
	a = a/norm(a);
	a(m-1) -= 1;
	mat V_(m,m,fill::eye);
	V_ = V_ - 2*a*a.t()/pow(norm(a),2);
	mat Q(m,m,fill::eye);
	mat Q_(m-1,m-1);mat R_(m-1,m-1);
	qr(Q_,R_,mat(m-1,m-1,fill::randn));
	Q.submat(0,0,m-2,m-2) = Q_;
	return U*S*Q*V_;
}

void gen_Xq(const string str_Q, mat& Xq){
	vector<string> tokens;
	string token;
	istringstream tokenStream(str_Q);
	while (getline(tokenStream, token, '-'))
	{
		tokens.push_back(token);
	}
	if (tokens[0].compare("I")==0){
		double q; istringstream(tokens[1]) >> q;
		Xq = pow(q,0.5)*mat(Xq.n_rows,Xq.n_cols,fill::eye);
	}else if (tokens[0].compare("diffusion")==0){
		mat Q_ (Xq.n_rows,Xq.n_rows);
		double q; istringstream(tokens[1]) >> q;
		unsigned int i,j; double k,d;
		for (i=0;i<Xq.n_rows;i++){
			for (j=0;j<Xq.n_rows;j++){
				if (i==j){k=0.1;}
				d = min(abs(int(i-j)),40-abs(int(i-j)));
				d = pow(d,2);
				Q_(i,j)= q*( exp(-d/30.0) + k );
			}
		}
		vec s(Xq.n_cols);
		if (Xq.n_cols<Xq.n_rows){
			eigs_sym(s,Xq,sp_mat(Q_),Xq.n_cols);
		}else{
			eig_sym(s,Xq,Q_);
		}
		Xq = Xq*diagmat(pow(s,0.5));
	}
}
mat ortho(const int& N, const bool& RR){
		mat U(N,N,fill::eye), R(N,N);
		if (RR){
			qr(U,R,mat(N,N,fill::randn));
		}
		vec v = U.col(0) - pow(N,-0.5)*vec(N,fill::ones);
		U = (mat(N,N,fill::eye) - 2*pow(norm(v),-2)*v*v.t())*U;
		return U.submat(0,1,N-1,N-1);
}
class dU{

	public:
	dU(){
		;
	}
	map<int,mat> d;
	mat get(int N){
		map<int,mat>::iterator d_it(d.find(N));
		mat U;
		if (d_it==d.end()){
			U = ortho(N,false);
			d[N] = U;
		}else{
			U = d_it->second;
		}
		return U;
	}
};

void print(string name, map<string, double> ocfg){
	cout << endl;
	cout << "Outputs of:  " << name << endl;
	for(map<string,double>::iterator it = ocfg.begin(); it != ocfg.end(); ++it) {
	  cout << it->first << " = " << it->second << endl;
	}
	cout << endl;
}
#endif /* C_FUNAUX_H_ */
