#include <Rcpp.h>

// Nelson/Siegel loss function for bonds 
RcppExport SEXP objfct_ns_bonds_Cpp(SEXP beta, SEXP m, SEXP cf, SEXP w, SEXP p){

	Rcpp::NumericVector betac(beta);
	Rcpp::NumericMatrix mc(m);
	Rcpp::NumericMatrix cfc(cf);
	Rcpp::NumericVector wc(w);
	Rcpp::NumericVector pc(p);
	Rcpp::NumericVector phat(pc.size());

	int i = 0;
	int j = 0;
	double s;
	double mse = 0;
	
	for (j = 0; j<mc.cols(); j++) {
		i = 0;
		while (i<mc.rows() && mc(i,j)>0) {
			s =  (betac[0] + betac[1] * ((1 - exp(-mc(i,j)/betac[3]))/(mc(i,j)/betac[3])) +
				  betac[2] * (((1 - exp(-mc(i,j)/betac[3]))/(mc(i,j)/betac[3])) - exp(-mc(i,j)/betac[3])))/100;
		
			phat[j] += cfc(i,j)*exp(-s*mc(i,j));
			i++;
		}
		mse += pow(pc[j] - phat[j],2)*wc[j];
	}

	return Rcpp::wrap(mse);
}

// Nelson/Siegel grid loss function for bonds 
RcppExport SEXP objfct_ns_bonds_gridCpp(SEXP beta, SEXP tau, SEXP m, SEXP cf, SEXP w, SEXP p){

	Rcpp::NumericVector betac(beta);
	Rcpp::NumericVector tauc(tau);
	Rcpp::NumericMatrix mc(m);
	Rcpp::NumericMatrix cfc(cf);
	Rcpp::NumericVector wc(w);
	Rcpp::NumericVector pc(p);
	Rcpp::NumericVector phat(pc.size());

	int i = 0;
	int j = 0;
	double s;
	double mse = 0;
	
	for (j = 0; j<mc.cols(); j++) {
		i = 0;
		while (i<mc.rows() && mc(i,j)>0) {
			s =  (betac[0] + betac[1] * ((1 - exp(-mc(i,j)/tauc[0]))/(mc(i,j)/tauc[0])) +
				  betac[2] * (((1 - exp(-mc(i,j)/tauc[0]))/(mc(i,j)/tauc[0])) - exp(-mc(i,j)/tauc[0])))/100;
		
			phat[j] += cfc(i,j)*exp(-s*mc(i,j));
			i++;
		}
		mse += pow(pc[j] - phat[j],2)*wc[j];
	}

	return Rcpp::wrap(mse);
}

// Svensson loss function for bonds 
RcppExport SEXP objfct_sv_bonds_Cpp(SEXP beta, SEXP m, SEXP cf, SEXP w, SEXP p){

	Rcpp::NumericVector betac(beta);
	Rcpp::NumericMatrix mc(m);
	Rcpp::NumericMatrix cfc(cf);
	Rcpp::NumericVector wc(w);
	Rcpp::NumericVector pc(p);
	Rcpp::NumericVector phat(pc.size());

	int i = 0;
	int j = 0;
	double s;
	double mse = 0;
	
	for (j = 0; j<mc.cols(); j++) {
		i = 0;
		while (i<mc.rows() && mc(i,j)>0) {
			s =  (betac[0] + betac[1] * ((1 - exp(-mc(i,j)/betac[3]))/(mc(i,j)/betac[3])) +
				  betac[2] * (((1 - exp(-mc(i,j)/betac[3]))/(mc(i,j)/betac[3])) - exp(-mc(i,j)/betac[3])) +
				  betac[4] * (((1 - exp(-mc(i,j)/betac[5]))/(mc(i,j)/betac[5])) - exp(-mc(i,j)/betac[5])))/100;
		
			phat[j] += cfc(i,j)*exp(-s*mc(i,j));
			i++;
		}
		mse += pow(pc[j] - phat[j],2)*wc[j];
	}

	return Rcpp::wrap(mse);
}

// Svensson grid loss function for bonds
RcppExport SEXP objfct_sv_bonds_gridCpp(SEXP beta, SEXP tau, SEXP m, SEXP cf, SEXP w, SEXP p){

	Rcpp::NumericVector betac(beta);
	Rcpp::NumericVector tauc(tau);
	Rcpp::NumericMatrix mc(m);
	Rcpp::NumericMatrix cfc(cf);
	Rcpp::NumericVector wc(w);
	Rcpp::NumericVector pc(p);
	Rcpp::NumericVector phat(pc.size());

	int i = 0;
	int j = 0;
	double s;
	double mse = 0;
	
	for (j = 0; j<mc.cols(); j++) {
		i = 0;
		while (i<mc.rows() && mc(i,j)>0) {
			s =  (betac[0] + betac[1] * ((1 - exp(-mc(i,j)/tauc[0]))/(mc(i,j)/tauc[0])) +
				  betac[2] * (((1 - exp(-mc(i,j)/tauc[0]))/(mc(i,j)/tauc[0])) - exp(-mc(i,j)/tauc[0])) +
				  betac[3] * (((1 - exp(-mc(i,j)/tauc[1]))/(mc(i,j)/tauc[1])) - exp(-mc(i,j)/tauc[1])))/100;
		
			phat[j] += cfc(i,j)*exp(-s*mc(i,j));
			i++;
		}
		mse += pow(pc[j] - phat[j],2)*wc[j];
	}

	return Rcpp::wrap(mse);
}

// Adjusted Svensson loss function for bonds
RcppExport SEXP objfct_asv_bonds_Cpp(SEXP beta, SEXP m, SEXP cf, SEXP w, SEXP p){

	Rcpp::NumericVector betac(beta);
	Rcpp::NumericMatrix mc(m);
	Rcpp::NumericMatrix cfc(cf);
	Rcpp::NumericVector wc(w);
	Rcpp::NumericVector pc(p);
	Rcpp::NumericVector phat(pc.size());

	int i = 0;
	int j = 0;
	double s;
	double mse = 0;
	
	for (j = 0; j<mc.cols(); j++) {
		i = 0;
		while (i<mc.rows() && mc(i,j)>0) {
			s =  (betac[0] + betac[1] * ((1 - exp(-mc(i,j)/betac[3]))/(mc(i,j)/betac[3])) +
				  betac[2] * (((1 - exp(-mc(i,j)/betac[3]))/(mc(i,j)/betac[3])) - exp(-mc(i,j)/betac[3])) +
				  betac[4] * (((1 - exp(-mc(i,j)/betac[5]))/(mc(i,j)/betac[5])) - exp(-2*mc(i,j)/betac[5])))/100;
		
			phat[j] += cfc(i,j)*exp(-s*mc(i,j));
			i++;
		}
		mse += pow(pc[j] - phat[j],2)*wc[j];
	}

	return Rcpp::wrap(mse);
}

// Adjusted Svensson grid loss function for bonds
RcppExport SEXP objfct_asv_bonds_gridCpp(SEXP beta, SEXP tau, SEXP m, SEXP cf, SEXP w, SEXP p){

	Rcpp::NumericVector betac(beta);
	Rcpp::NumericVector tauc(tau);
	Rcpp::NumericMatrix mc(m);
	Rcpp::NumericMatrix cfc(cf);
	Rcpp::NumericVector wc(w);
	Rcpp::NumericVector pc(p);
	Rcpp::NumericVector phat(pc.size());

	int i = 0;
	int j = 0;
	double s;
	double mse = 0;
	
	for (j = 0; j<mc.cols(); j++) {
		i = 0;
		while (i<mc.rows() && mc(i,j)>0) {
			s =  (betac[0] + betac[1] * ((1 - exp(-mc(i,j)/tauc[0]))/(mc(i,j)/tauc[0])) +
				  betac[2] * (((1 - exp(-mc(i,j)/tauc[0]))/(mc(i,j)/tauc[0])) - exp(-mc(i,j)/tauc[0])) +
				  betac[3] * (((1 - exp(-mc(i,j)/tauc[1]))/(mc(i,j)/tauc[1])) - exp(-2*mc(i,j)/tauc[1])))/100;
		
			phat[j] += cfc(i,j)*exp(-s*mc(i,j));
			i++;
		}
		mse += pow(pc[j] - phat[j],2)*wc[j];
	}

	return Rcpp::wrap(mse);
}

// Gradient of Svensson grid loss function for bonds
RcppExport SEXP grad_sv_bonds_gridCpp(SEXP beta, SEXP tau, SEXP m, SEXP cf, SEXP w, SEXP p){
	Rcpp::NumericVector betac(beta);
	Rcpp::NumericVector tauc(tau);
	Rcpp::NumericMatrix mc(m);
	Rcpp::NumericMatrix cfc(cf);
	Rcpp::NumericVector wc(w);
	Rcpp::NumericVector pc(p);
	
	const int N = mc.rows();
	const int M = mc.cols();
	
	Rcpp::NumericMatrix emt1(N, M);
	Rcpp::NumericMatrix emt2(N, M);
	Rcpp::NumericMatrix t1emt1(N, M);
	Rcpp::NumericMatrix emt1tm(N, M);
	Rcpp::NumericMatrix emt2tm(N, M);
	Rcpp::NumericMatrix acf(N, M);
	Rcpp::NumericMatrix dm(N, M);
	Rcpp::NumericVector csacf(M);
	Rcpp::NumericVector csdm(M);
	Rcpp::NumericVector csdt1emt1(M);
	Rcpp::NumericVector csdmem1tm(M);
	Rcpp::NumericVector csdmem2tm(M);
	Rcpp::NumericVector b(M);
	Rcpp::NumericVector gbeta(4);
	
	int i = 0;
	int j = 0;
		
	for (j = 0; j<M; j++) {
		i = 0;
		while (i<N && mc(i,j)>0) {
			emt1(i,j) = exp(-mc(i,j)/tauc[0]);
			emt2(i,j) = exp(-mc(i,j)/tauc[1]);
			t1emt1(i,j) = tauc[0]*(1 - emt1(i,j));
			emt1tm(i,j) = - emt1(i,j) + t1emt1(i,j)/mc(i,j);
			emt2tm(i,j) = - emt2(i,j) + tauc[1]*(1 - emt2(i,j))/mc(i,j);					   
			acf(i,j) = cfc(i,j)*exp((-betac[0] - betac[2]*emt1tm(i,j) - betac[3]*emt2tm(i,j) - (betac[1]*t1emt1(i,j))/mc(i,j))*mc(i,j)/100);
			dm(i,j) = acf(i,j)/100*mc(i,j);
			csacf[j] += acf(i,j);
			csdm[j] += dm(i,j) ;
			csdt1emt1[j] += acf(i,j)/100*t1emt1(i,j);
			csdmem1tm[j] += dm(i,j)*emt1tm(i,j);
			csdmem2tm[j] += dm(i,j)*emt2tm(i,j);
			i++;
		}
		b[j] = -2*wc[j]*(pc[j] - csacf[j]);
		gbeta[0] += b[j]*(-csdm[j]);
		gbeta[1] += b[j]*(-csdt1emt1[j]);
		gbeta[2] += b[j]*(-csdmem1tm[j]);
		gbeta[3] += b[j]*(-csdmem2tm[j]);
	}
	
	return gbeta;
}


// Gradient of Svensson grid loss function for yields
RcppExport SEXP grad_sv_gridCpp(SEXP beta, SEXP tau, SEXP m, SEXP y){
	Rcpp::NumericVector betac(beta);
	Rcpp::NumericVector tauc(tau);
	Rcpp::NumericVector mc(m);
	Rcpp::NumericVector yc(y);
	Rcpp::NumericVector yhat(yc.size());
	Rcpp::NumericVector gbeta(4);
	
	int j = 0;
	
	for (j = 0; j<mc.size(); j++) {
		
		gbeta[0] += -2*(-betac[0] - betac[2]*(-exp(-mc[j]/tauc[0]) + (tauc[0]*(1 - exp(-mc[j]/tauc[0])))/mc[j]) - 
						betac[3]*(-exp(-mc[j]/tauc[1]) + (tauc[1]*(1 - exp(-mc[j]/tauc[1])))/mc[j]) - 
						(betac[1]*tauc[0]*(1 - exp(-mc[j]/tauc[0])))/mc[j] + yc[j]);
		
		gbeta[1] +=  (-2*tauc[0]*(1 - exp(-mc[j]/tauc[0]))*(-betac[0] - betac[2]*(-exp(-mc[j]/tauc[0]) + (tauc[0]*(1 - exp(-mc[j]/tauc[0])))/mc[j]) - 
					 betac[3]*(-exp(-mc[j]/tauc[1]) + (tauc[1]*(1 - exp(-mc[j]/tauc[1])))/mc[j]) - 
															(betac[1]*tauc[0]*(1 - exp(-mc[j]/tauc[0])))/mc[j] + yc[j]))/mc[j];
		
		gbeta[2] += -2*(-exp(-mc[j]/tauc[0]) + (tauc[0]*(1 - exp(-mc[j]/tauc[0])))/mc[j])*
		(-betac[0] - betac[2]*(-exp(-mc[j]/tauc[0]) + (tauc[0]*(1 - exp(-mc[j]/tauc[0])))/mc[j]) - 
		 betac[3]*(-exp(-mc[j]/tauc[1]) + (tauc[1]*(1 - exp(-mc[j]/tauc[1])))/mc[j]) - 
		 (betac[1]*tauc[0]*(1 - exp(-mc[j]/tauc[0])))/mc[j] + yc[j]);
		
		
		gbeta[3] += -2*(-exp(-mc[j]/tauc[1]) + (tauc[1]*(1 - exp(-mc[j]/tauc[1])))/mc[j])*
		(-betac[0] - betac[2]*(-exp(-mc[j]/tauc[0]) + (tauc[0]*(1 - exp(-mc[j]/tauc[0])))/mc[j]) - 
		 betac[3]*(-exp(-mc[j]/tauc[1]) + (tauc[1]*(1 - exp(-mc[j]/tauc[1])))/mc[j]) - 
		 (betac[1]*tauc[0]*(1 - exp(-mc[j]/tauc[0])))/mc[j] + yc[j]);
	}
	
	return gbeta;
}
