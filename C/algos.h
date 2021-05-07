#ifndef C_ALGOS_H_
#define C_ALGOS_H_

#include <iostream>
#include <string>
#include <map>
#include <math.h>
#include <time.h>

#include "armadillo"
#include "lorenz95.h"
#include "funaux.h"

using namespace arma;
using namespace std;

map<string,double> IEnKS_bundle(map<string, string> icfg){
	// Initialisations /////////////////////////////////////////////////////////
	unsigned int jmax, Ns, Nmb, Nspin, Ndt, L, S, Na;
	int Ptol;
	double sigmaB, sigmaR, sigmaQ, Ppas, inflation, burn;
	string str_M, str_H;
	istringstream(icfg["Ns"]) >> Ns;
	istringstream(icfg["Nmb"]) >> Nmb;
	istringstream(icfg["Nspin"]) >> Nspin;
	istringstream(icfg["Ndt"]) >> Ndt;
	istringstream(icfg["L"]) >> L;
	istringstream(icfg["S"]) >> S;
	istringstream(icfg["Na"]) >> Na;
	istringstream(icfg["sigmaB"]) >> sigmaB;
	istringstream(icfg["sigmaR"]) >> sigmaR;
	istringstream(icfg["sigmaQ"]) >> sigmaQ;
	istringstream(icfg["Ppas"]) >> Ppas;
	istringstream(icfg["inflation"]) >> inflation;
	istringstream(icfg["jmax"]) >> jmax;
	istringstream(icfg["Ptol"]) >> Ptol;
	istringstream(icfg["burn"]) >> burn;
	str_M = icfg["str_M"];
	str_H = icfg["str_H"];
	unsigned int i(0), i_(0), j(0), js(0), l, n_mod(0), NM(Nmb+1), K(L-S+1);
	double errS(0), errF(0), crit, pas(pow(10,Ppas));
	mat F, Eb(Ns,NM), Ea(Ns,NM), Xb(Ns,Nmb), U(ortho(NM,false)), V,
		C(Nmb,Nmb), E(Ns,NM), Urand, I(mat(Nmb,Nmb,fill::eye));
	vec x0, f, xt ,xt0(3.0*vec(Ns,fill::ones)+vec(Ns,fill::randn)), xtS(Ns), xtL(Ns),
		y(Ns*S), xb(Ns), v, w(Nmb), dw(Nmb), g(Nmb);
	clock_t t, t_tot;
	double dt_mod(0), dt_svd(0);
	span l_,s_;
	M(xt0,Nspin,str_M);// Spin-up
	////////////////////////////////////////////////////////////////////////////

	// Assimilation ////////////////////////////////////////////////////////////
	t_tot=clock();
	try{
	while (i<Na){
		cout << "\r" << i*100.0/Na << " %";

		// Génération des observations
		xt = xt0;
		for(l=0;l<=L;l++){
			if (i==0 and l==0){
				Eb = sigmaB*mat(Ns,NM,fill::randn);
				Eb.each_col() += xt;
			}
			if (l>=K){
				l_ = span((l-K)*Ns,(l-K+1)*Ns-1);
				y(l_) = H(xt,str_H) + sigmaR*vec(Ns,fill::randn);
			}
			if(l==L){
				xtL = xt;
			}
			M(xt,Ndt,str_M);
			xt += sigmaQ*vec(Ns,fill::randn);
			if (l+1==S){
				xtS = xt;
			}
		}
		// Extension de la dimension
		xb = mean(Eb,1);
		Xb = inflation*pow(Nmb,-0.5)*Eb*U;
		w.zeros();
		j = 0; crit = 1.;
		do{
			g = w; C = mat(Nmb,Nmb,fill::eye);
			x0 = xb+Xb*w;
			E = pas*Xb*U.t(); E.each_col() += x0;
			for(l=0;l<=L;l++){
				if (l>0){
					t = clock();
					M(E,Ndt,str_M);
					if(i>=burn*Na){
						dt_mod += double(clock()-t);
						n_mod += Ndt*E.n_cols;
					}
				}
				if(l>=K){
					V = H(E,str_H);
					f = mean(V,1);
					F = V*U/pas;
					V = F.t()*pow(sigmaR,-2);
					l_ = span((l-K)*Ns,(l-K+1)*Ns-1);
					g -= V*(y(l_)-f) ;
					C += V*F;
				}
			}
			if (C.is_finite() and g.is_finite()){
				dw = solve(C,g);
			}else{
				cout << "infinite derivatives" << endl;
				throw 1;
			}
			w -= dw;
			crit = norm(dw)*pow(Nmb,-0.5);
			j++;
		}while( (crit > pow(10,Ptol)) and (j<jmax));

		if (i>=burn*Na){
			errS += norm(xt0-x0);
			errF += norm(xtL-mean(E,1));
			js += j-1;
			i_++;
		}

		// Génération de l'ensemble analysé
		t = clock();
		eig_sym(v,V,C);
		if(i>burn*Na){dt_svd += double(clock()-t);}
		// Propagation
		Ea = Xb*pow(Nmb,0.5)*V*diagmat(pow(v,-0.5))*ortho(NM,true).t();
		Ea.each_col() += xb+Xb*w;
		for (l=1;l<=S;l++){
			t = clock();
			M(Ea,Ndt,str_M);
			if(i>=burn*Na){
				dt_mod += double(clock()-t);
				n_mod += Ndt*Ea.n_cols;
			}
		}

		// Marginalisation
		Eb = Ea;
		xt0 = xtS;
		i += S;
	}}catch(...){
		errF = -1.0;
		errS = -1.0;
		i_ = 1;
		cout << "Echec de l'assimilation: " << icfg["name"] << endl;
		cout << "i= " << i << endl;
		// throw 1;
//		cout << "dbg= " << dbg << endl;
	}
	////////////////////////////////////////////////////////////////////////////
	map<string,double> ocfg;
	ocfg["t_tot"] = double(clock()-t_tot)/CLOCKS_PER_SEC;
	ocfg["t_mod"] = dt_mod/CLOCKS_PER_SEC;
	ocfg["t_svd"] = dt_svd/CLOCKS_PER_SEC;
	ocfg["n_mod"] = n_mod/(S*i_);
	ocfg["j"] = double(js)/i_;
	ocfg["errF"] = errF/(i_*pow(Ns,0.5));
	ocfg["errS"] = errS/(i_*pow(Ns,0.5));

	print(icfg["name"],ocfg);

	return ocfg;
}

map<string,double> IEnKF_Q_bundle(map<string, string> icfg){
	int Nmb, Nmq, Ns, Na, No, Ndt, Nspin, Ptol, jmax;
	string str_M, str_H;
	double sigmaB, sigmaR, sigmaQ, inflation, Ppas, burn;
	istringstream(icfg["Ppas"]) >> Ppas;
	istringstream(icfg["Nmb"]) >> Nmb;
	istringstream(icfg["Nmq"]) >> Nmq;
	istringstream(icfg["Ns"]) >> Ns;
	istringstream(icfg["Na"]) >> Na;
	istringstream(icfg["No"]) >> No;
	istringstream(icfg["Ndt"]) >> Ndt;
	istringstream(icfg["Nspin"]) >> Nspin;
	istringstream(icfg["sigmaB"]) >> sigmaB;
	istringstream(icfg["sigmaR"]) >> sigmaR;
	istringstream(icfg["sigmaQ"]) >> sigmaQ;
	istringstream(icfg["Ptol"]) >> Ptol;
	istringstream(icfg["jmax"]) >> jmax;
	istringstream(icfg["inflation"]) >> inflation;
	istringstream(icfg["burn"]) >> burn;
	str_M = icfg["str_M"];
	str_H = icfg["str_H"];

	int i(0), i_(0), j, NMb(Nmb+1), NMq(Nmq+1), NM(NMb+NMq);
	mat I(NM,NM,fill::eye);
	vec x1t = 3.0*vec(Ns,fill::ones)+vec(Ns,fill::randn);
	M(x1t,Nspin,str_M);
	mat E1a(sigmaB*mat(Ns,NMb,fill::randn));
	E1a.each_col() += x1t;
	vec x2t,x1a,w(NM),dw,y2,x1,x2,s,g;
	mat A1a,HA2,HA2q,U,V,E1,E2,HA,A2,A,C;
	double crit, tol(pow(10,Ptol)),errF(0),errS(0),js(0), n_mod(0), pas(pow(10,Ppas)),
			dt_svd(0), dt_mod(0);
	clock_t t_tot,t;
	//Construction de A2q tq A2q*A2q.T = q**2*I
	//et A2q*1=0 avec une matrice de householder
	span l_, s_;
	mat A2q(Ns,NMq,fill::zeros);
	l_ = span(0,Ns-1); s_ = span(0,Nmq-1);
	A2q(l_,s_) = sigmaQ*mat(Ns,Nmq,fill::eye);
	vec v(NMq,fill::ones);
	v = v/norm(v);
	v(NMq-1) -= 1;
	A2q = A2q*(mat(NMq,NMq,fill::eye)-2*v*v.t()/pow(norm(v),2));

	t_tot = clock();
	try{
	while(i<Na){
		x1a = moy(E1a);
		A1a = ano(E1a);
		w.zeros();
		j = 0; crit = 1.;

		//Génération d'une obs
		x2t = x1t;
		M(x2t,Ndt,str_M);
		x2t += sigmaQ*vec(Ns,fill::randn);
		y2 = H(x2t,str_H) + sigmaB*vec(No,fill::randn);

		//Analysis
		do{
			x1 = x1a + A1a*w.subvec(0,NMb-1);
			E1 = ens(x1,pas*A1a);
			E2 = E1;
			t = clock();
			M(E2,Ndt,str_M);
			if(i>=burn*Na){
				dt_mod += double(clock()-t);
				n_mod += Ndt*E2.n_cols;
			}
			HA2 = ano(mat(H(E2,str_H)))/pas;
			HA2q = ano(mat(H(ens(moy(E2),pas*A2q),str_H)))/pas;
			HA = join_rows(HA2,HA2q);
			x2 = moy(E2) + A2q*w.subvec(NMb,NM-1);
			g = w - HA.t()*(y2-H(x2,str_H))*pow(sigmaR,-2);
			C = I + HA.t()*HA*pow(sigmaR,-2);
			if (C.is_finite() and g.is_finite()){
				dw = -solve(C,g);
			}else{
				cout << "infinite derivatives" << endl;
				throw 1;
			}
			crit = norm(dw)*pow(NM,-0.5);
			w += dw;
			j++;
		}while( (crit > tol) and (j<jmax));
		t = clock();
		eig_sym(s,V,C);
		if(i>=burn*Na){dt_svd += double(clock()-t);}
		A2 = ano(E2)/pas;
		A = join_rows(A2,A2q)*V*diagmat(pow(s,-0.5))*V.t();
		t = clock();
		A2 = inflation*SR(A,NMb);
		if(i>=burn*Na){dt_svd += double(clock()-t);}
		E2 = ens(x2,A2);

		//Errors
		if(i>=burn*Na){
			errS += norm(x1-x1t);
			errF += norm(x2-x2t);
			js += j-1;
			i_ ++;
		}

		//Propagation
		E1a = E2;
		x1t = x2t;

		i += 1;

		cout << i*100.0/(Na-1) << " % \r";

	}}catch(...){
		errF = -1.0;
		errS = -1.0;
		js = -1.0;
		i_ = 1;
		cout << "Echec de l'assimilation: " << icfg["name"] << endl;
	}
	map<string,double> ocfg;
	ocfg["t_tot"] = double(clock()-t_tot)/CLOCKS_PER_SEC;
	ocfg["t_mod"] = dt_mod/CLOCKS_PER_SEC;
	ocfg["t_svd"] = dt_svd/CLOCKS_PER_SEC;
	ocfg["n_mod"] = n_mod/i_;
	ocfg["j"] = double(js)/i_;
	ocfg["errF"] = errF/(i_*pow(Ns,0.5));
	ocfg["errS"] = errS/(i_*pow(Ns,0.5));

	print(icfg["name"],ocfg);

	return ocfg;
}

map<string,double> IEnKF_Q_transform(map<string, string> icfg){
	int Nmb, Nmq, Ns, Na, No, Ndt, Nspin, Ptol, jmax;
	string str_Q, str_M, str_H;
	double sigmaB, sigmaR, sigmaQ, inflation, burn;
	istringstream(icfg["Nmb"]) >> Nmb;
	istringstream(icfg["Nmq"]) >> Nmq;
	istringstream(icfg["Ns"]) >> Ns;
	istringstream(icfg["Na"]) >> Na;
	istringstream(icfg["No"]) >> No;
	istringstream(icfg["Ndt"]) >> Ndt;
	istringstream(icfg["Nspin"]) >> Nspin;
	istringstream(icfg["sigmaB"]) >> sigmaB;
	istringstream(icfg["sigmaR"]) >> sigmaR;
	istringstream(icfg["sigmaQ"]) >> sigmaQ;
	istringstream(icfg["Ptol"]) >> Ptol;
	istringstream(icfg["jmax"]) >> jmax;
	istringstream(icfg["inflation"]) >> inflation;
	istringstream(icfg["burn"]) >> burn;
	str_M = icfg["str_M"];
	str_H = icfg["str_H"];
	int n_mod(0), NMb(Nmb+1), NMq(Nmq+1), NM(NMb+NMq), i(0), i_(0), j(0);
	mat I(NM,NM,fill::eye);
	vec zeros(NM,fill::zeros);
	vec x1t = 3.0*vec(Ns,fill::ones)+vec(Ns,fill::randn);
	M(x1t,Nspin,str_M);
	mat E1a = ens(x1t,sigmaB*mat(Ns,NMb,fill::randn));
	vec x2t,x1a,w,dw,y2,x1,x2,s,g;
	mat A1a,HA2,HA2q,D,U,V,T,iT,E1,E2,HA,A2,A,C;

	double dt_svd(0), dt_mod(0), crit(1), tol(pow(10,Ptol)),errF(0),errS(0),js(0),prop(0);
	clock_t t_tot;
	span l_, s_;

	//Construction de A2q tq A2q*A2q.T = q**2*I
	//et A2q*1=0 avec une matrice de householder
	mat A2q(Ns,NMq,fill::zeros);
	l_  = span(0,Ns-1); s_ = span(0,Nmq-1);
	A2q(l_,s_) = sigmaQ*mat(Ns,Nmq,fill::eye);
	vec v(NMq,fill::ones);
	v = v/norm(v);
	v(NMq-1) -= 1;
	A2q = A2q*(mat(NMq,NMq,fill::eye)-2*v*v.t()/pow(norm(v),2));

	t_tot = clock();
	try{
	while(i<Na){
		x1a = moy(E1a);
		A1a = ano(E1a);
		D = I;
		w = zeros;
		j = 0; crit = 1.;

		//Génération d'une obs
		x2t = x1t;
		M(x2t,Ndt,str_M);
		x2t += sigmaQ*vec(Ns,fill::randn);
		y2 = H(x2t,str_H) + sigmaB*vec(No,fill::randn);

		//Analysis
		do{
			x1 = x1a + A1a*w.subvec(0,NMb-1);
			eig_sym(s,V,D.submat(0,0,NMb-1,NMb-1));
			T = V*diagmat(pow(s,0.5))*V.t();
			iT = V*diagmat(pow(s,-0.5))*V.t();
			E1 = ens(x1,A1a*T);
			E2 = E1;
			M(E2,Ndt,str_M); n_mod += E2.n_cols;
			if(i>=burn*Na){
				prop += Ndt*E2.n_cols;
			}
			HA2 = ano(mat(H(E2,str_H)))*iT;
			HA2q = ano(mat(H(ens(moy(E2),A2q),str_H)));
			HA = join_rows(HA2,HA2q);
			x2 = moy(E2) + A2q*w.subvec(NMb,NM-1);
			g = w - HA.t()*(y2-H(x2,str_H))*pow(sigmaR,-2);
			C = I + HA.t()*HA*pow(sigmaR,-2);
			if (C.is_finite() and g.is_finite()){
				D = inv(C);
				dw = -solve(C,g);
			}else{
				cout << "infinite derivatives" << endl;
				throw 1;
			}
			crit = norm(dw)*pow(NM,-0.5);
			w += dw;
			j++;
		}while( (crit > tol) and (j<jmax));
		if (D.is_finite()){
			eig_sym(s,V,D);
		}else{
			cout << "infinite derivatives" << endl;
			throw 1;
		}
		A2 = ano(E2)*iT;
		A = join_rows(A2,A2q)*V*diagmat(pow(s,0.5))*V.t();
		A2 = inflation*SR(A,NMb);
		E2 = ens(x2,A2);

		//Errors
		if (i>=burn*Na){
			errS += norm(x1-x1t);
			errF += norm(x2-x2t);
			js += j-1;
			i_ ++;
		}

		//Propagation
		E1a = E2;
		x1t = x2t;

		i += 1;
		cout << i*100.0/(Na-1) << " % \r";

	}}catch(...){
		errF = -1.0;
		errS = -1.0;
		js = -1.0;
		i_ = 1;
		cout << "Echec de l'assimilation: " << icfg["name"] << endl;
	}
	map<string,double> ocfg;
	ocfg["t_tot"] = double(clock()-t_tot)/CLOCKS_PER_SEC;
	ocfg["t_mod"] = dt_mod/CLOCKS_PER_SEC;
	ocfg["t_svd"] = dt_svd/CLOCKS_PER_SEC;
	ocfg["n_mod"] = n_mod/i_;
	ocfg["j"] = double(js)/i_;
	ocfg["errF"] = errF/(i_*pow(Ns,0.5));
	ocfg["errS"] = errS/(i_*pow(Ns,0.5));

	print(icfg["name"],ocfg);

	return ocfg;
}

map<string,double> IEnKS_Q(map<string,string> icfg){
	// Initialisations /////////////////////////////////////////////////////////
	int jmax, Ns, Nmb, Nmq, Nspin, Ndt, L, S, G, Na;
	int Ptol, Ppas, dbg(0), Nm0(0);
	double sigmaB,  sigmaR, sigmaQ, inflation, burn;
	string str_M, str_H, str_S, str_lin, str_stateUp;
	istringstream(icfg["Ns"]) >> Ns;
	istringstream(icfg["Nmb"]) >> Nmb;
	istringstream(icfg["Nmq"]) >> Nmq;
	istringstream(icfg["Nspin"]) >> Nspin;
	istringstream(icfg["Ndt"]) >> Ndt;
	istringstream(icfg["L"]) >> L;
	istringstream(icfg["S"]) >> S;
	istringstream(icfg["G"]) >> G;
	istringstream(icfg["Na"]) >> Na;
	istringstream(icfg["sigmaB"]) >> sigmaB;
	istringstream(icfg["sigmaR"]) >> sigmaR;
	istringstream(icfg["sigmaQ"]) >> sigmaQ;
	istringstream(icfg["Ppas"]) >> Ppas;
	istringstream(icfg["inflation"]) >> inflation;
	istringstream(icfg["jmax"]) >> jmax;
	istringstream(icfg["Ptol"]) >> Ptol;
	istringstream(icfg["burn"]) >> burn;
	str_M = icfg["str_M"];
	str_H = icfg["str_H"];
	str_lin = icfg["str_lin"];
	str_stateUp = icfg["str_stateUp"];
	//////////////////////////////////////////////////////////////////////////////

	// Definitions ///////////////////////////////////////////////////////////////
	int i(0), i_(0), j, js(0), l, ncv(0), n_mod(0),
		K(L-S+1), D(max(G,K-1)), Nm(Nmb+(L-D)*Nmq), Nm_(Nm+S*Nmq), NM(Nm+1),
		NM_(Nm_+1);
	vector<string> tokens;
	string token;
	istringstream tokenStream(icfg["str_S"]);
	while (getline(tokenStream, token, '_'))
	{
		tokens.push_back(token);
	}
	str_S = tokens[0];
	if (tokens.size()>1){
		istringstream(tokens[1]) >> Nm0;
	}else if(str_S.compare("full")==0){
		Nm0 = Nm;
	}

	double errS(0), errF(0), stop, pas(pow(10,Ppas)), tol(pow(10,Ptol)), dt_svd(0), dt_mod(0);
	vec v, va((L+S+1)*Ns), dw, f, g, w(Nm), vb, xt0(3.0*vec(Ns,fill::ones)+vec(Ns,fill::randn)), y(S*Ns), xtS, xtL, xt;
	M(xt0,Nspin,str_M);// Spin-up
	mat Q, R, Ea, Va((L+S+1)*Ns,Nm_), F, V, C, W, iW, Vb, E((L+1)*Ns,NM,fill::zeros), Vq, P(Nm,Nm,fill::zeros), I(Nm,Nm,fill::eye);
	for (l=0;l<Nm;l++){
		P(l,Nm-l-1) = 1;
	}
	Vq = sigmaQ*mat(Ns,Nmq,fill::eye);
	dU U;
	clock_t t_tot;
	span l_,s_;
	////////////////////////////////////////////////////////////////////////////

	// Assimilation ////////////////////////////////////////////////////////////
	t_tot = clock();
	try{
	while (i<Na){
		cout << "\r" << i*100.0/Na << " %";

		// Génération des observations /////////////////////////////////////////
		xt = xt0;
		for(l=0;l<=L;l++){
			if (i==0 and l<=G){
				l_ = span(l*Ns,(l+1)*Ns-1); s_ = span::all;
				E(l_,s_) = sigmaB*mat(Ns,NM,fill::randn);
				E(l_,s_).each_col() += xt;
			}
			if (i==0 and l>G){
				l_ = span(l*Ns,(l+1)*Ns-1); s_ = span::all;
				E(l_,s_) = sigmaQ*mat(Ns,NM,fill::randn);
			}
			if (i==0 and l==L){
				vb = mean(E,1);
				Vb = pow(Nm,-0.5)*E*U.get(NM);
			}
			if (l>=K){
				l_ = span((l-K)*Ns,(l+1-K)*Ns-1);
				y(l_) = H(xt,str_H) + sigmaR*vec(Ns,fill::randn);
			}
			if (l==L){
				xtL = xt;
			}
			M(xt,Ndt,str_M);
			xt += sigmaQ*vec(Ns,fill::randn);
			if (l+1==S){
				xtS = xt;
			}
		}
		////////////////////////////////////////////////////////////////////////

		// Analyse /////////////////////////////////////////////////////////////
		w.zeros();
		W = pas*I;
		iW = I/pas;
		j = 0;
		stop = 1;
		do{
			g = w; C = I;
			E = pow(Nm,0.5)*Vb*W*U.get(NM).t();
			E.each_col() += (vb+Vb*w);dbg++;
			for(l=0;l<=L;l++){
				l_ = span(l*Ns,(l+1)*Ns-1);s_ = span::all;
				if (l>G){
					if (str_S.compare("sp")==0){Nm0=Nm_-(L+S-max(D,l-1))*Nmq;}
					E(l_,s_) += M_surrogate(
							E(span((l-1)*Ns,l*Ns-1),s_),
							Nm0, Ndt,
							str_S, str_M,
							U
							);
					if(i>=burn*Na){n_mod += Ndt*(Nm0+1);}
				}
				if (l>=K){
					V = H(E(l_,s_),str_H);
					f = mean(V,1);
					F = pow(Nm,-0.5)*V*U.get(NM)*iW;
					V = pow(sigmaR,-2)*F.t();
					g += V*(f-y(span((l-K)*Ns,(l-K+1)*Ns-1)));
					C += V*F;dbg++;
				}
			}
			if (C.is_finite() and g.is_finite()){
				chol(V, P*C*P,"lower");
				if(str_lin.compare("bundle")==0){
					solve(dw,trimatl(V),-P*g);
					solve(dw,trimatu(V.t()),dw);
					dw = P*dw;
				}else if (str_lin.compare("transform")==0){
					solve(W,trimatl(V),I);
					W = P*W.t()*P;
					iW = P*V.t()*P;
					dw = -W*W.t()*g;dbg++;
				}
				w += dw;
				stop = pow(Nm,-0.5)*norm(dw);
			}else{
				cout << "infinite derivatives" << endl;
				throw 1;
			}
			j++;
		}while( (stop > tol) and (j<jmax));
		if(j >= jmax){
			ncv ++;
		}
		////////////////////////////////////////////////////////////////////////

		// Erreurs /////////////////////////////////////////////////////////////
		if (i>=burn*Na){
			l_ = span(0,Ns-1); s_ = span::all;
			errS += norm(xt0-mean(E(l_,s_),1));
			l_ = span(L*Ns,(L+1)*Ns-1);
			errF += norm(xtL-mean(E(l_,s_),1));
			js += j-1;
			i_++;
		}
		////////////////////////////////////////////////////////////////////////

		// Propagation /////////////////////////////////////////////////////////
		if (str_lin.compare("bundle")==0){
			solve(W,trimatl(V),I);
			W = P*W.t()*P;
		}
		va.zeros(); Va.zeros();
		l_ = span(0,(L+1)*Ns-1); s_ = span(0,Nm-1);
		va(l_) = vb + Vb*w;
		Va(l_,s_) = Vb*W;
		for (l=1;l<=S;l++){
			l_ = span((L+l)*Ns,(L+l+1)*Ns-1); s_ = span(Nm+(l-1)*Nmq,Nm+l*Nmq-1);
			Va(l_,s_) = Vq;
		}
		Ea = pow(Nm_,0.5)*Va*U.get(NM_).t();
		Ea.each_col() += va;
		for (l=G+1;l<=G+S;l++){
			if (str_stateUp.compare("nlin")==0 or l>L){
				l_ = span(l*Ns,(l+1)*Ns-1), s_ = span::all;
				if (str_S.compare("sp")==0){Nm0=Nm_-(L+S-max(D,l-1))*Nmq;}
				Ea(l_,s_) +=  M_surrogate(
						Ea(span((l-1)*Ns,l*Ns-1),s_),
						Nm0, Ndt, str_S, str_M, U);dbg++;
				if(i>=burn*Na){
					n_mod += Ndt*(Nm0+1);
				}
			}else if (str_stateUp.compare("lin")==0 and l<=L){
				l_ = span(l*Ns,(l+1)*Ns-1), s_ = span::all;
				Ea(l_,s_) = pow(double(Nm_)/Nm,0.5)*E(l_,s_)*U.get(NM)*iW*W*mat(Nm,Nm_,fill::eye)*U.get(NM_).t();dbg++;
				Ea(l_,s_).each_col() += mean(E(l_,s_),1);dbg++;
			}
		}
		l_ = span(S*Ns,(L+S+1)*Ns-1), s_ = span::all;
		vb = mean(Ea(l_,s_),1);
		V = pow(Nm_,-0.5)*Ea(l_,s_)*U.get(NM_);dbg++;
		l_ = span(0,(D+1)*Ns-1);
		if (Nmb<(D+1)*Ns){
			eigs_sym(v,V,sp_mat(V(l_,s_)*V(l_,s_).t()),Nmb);dbg++;
		}else{
			eig_sym(v,V,V(l_,s_)*V(l_,s_).t());
		}
		qr(Q,R,mat(Nmb,Nmb,fill::randn));
		Vb.zeros();
		l_ = span(0,(D+1)*Ns-1); s_ = span(0,Nmb-1);
		Vb(l_,s_) = inflation*V*diagmat(pow(v,0.5))*Q;dbg++;
		for (l=D+1;l<=L;l++){
			l_ = span(l*Ns,(l+1)*Ns-1); s_ = span(Nmb+(l-D-1)*Nmq,Nmb+(l-D)*Nmq-1);
			Vb(l_,s_) = Vq;dbg++;
		}

		xt0 = xtS;
		i += S;
		////////////////////////////////////////////////////////////////////////
	}}catch(...){
		errF = -1.0;
		errS = -1.0;
		i_ = 1;
		cout << "Echec de l'assimilation: " << icfg["name"] << endl;
		cout << "i= " << i << endl;
		cout << "j= " << j << endl;
		cout << "l= " << l << endl;
		cout << "dbg= " << dbg << endl;
	}

	// Post-traitement /////////////////////////////////////////////////////////
	map<string,double> ocfg;
		ocfg["t_tot"] = double(clock()-t_tot)/CLOCKS_PER_SEC;
		ocfg["t_mod"] = dt_mod/CLOCKS_PER_SEC;
		ocfg["t_svd"] = dt_svd/CLOCKS_PER_SEC;
		ocfg["n_mod"] = n_mod/(S*i_);
		ocfg["j"] = double(js)/i_;
		ocfg["errF"] = errF/(i_*pow(Ns,0.5));
		ocfg["errS"] = errS/(i_*pow(Ns,0.5));

		print(icfg["name"],ocfg);

		return ocfg;}

map<string,double> algo_selector(map<string, string> icfg){
	map<string,double> ocfg;
	if (icfg["algo"].compare("IEnKS_bundle")==0){
		ocfg = IEnKS_bundle(icfg);
	}else if(icfg["algo"].compare("IEnKF_Q_bundle")==0){
		ocfg = IEnKF_Q_bundle(icfg);
	}else if(icfg["algo"].compare("IEnKF_Q_transform")==0){
		ocfg = IEnKF_Q_transform(icfg);
	}else if(icfg["algo"].compare("IEnKS_Q")==0){
		ocfg = IEnKS_Q(icfg);
	}else{
		cout << "Unknown algorithm: " << icfg["algo"] << endl;
	}
	return ocfg;

}
#endif /* C_ALGOS_H_ */
