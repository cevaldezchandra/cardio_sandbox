import copy
import pandas as pd
import numpy as np
import scipy as sp

class EncPredictor:
	def __init__(self, chi2s, ld_scores, n_idvs, posterior_betas=None):
		self.chi2s = chi2s # np array
		self.ld_scores = ld_scores # np array
		self.n_idvs = n_idvs
		self.posterior_betas = posterior_betas
		self.est_heritability = None
		self.beta_hats = None
		self.predictions = None

	def estimate_ldsc_heritability(self):
		# Input: n_idvs = number of individuals
		# Output: Estimates heritiability using LD-score
		n_snps = self.chi2s.shape[0]
		ldsc_est_heritability = ((self.chi2s.mean() - 1) * n_snps) / (self.ld_scores.mean() * self.n_idvs)
		self.est_heritability = ldsc_est_heritability
		return self.est_heritability

	@staticmethod
	def static_predict(genos, betas):
		prediction = np.dot(genos, betas)
		return prediction

	def predict(self, genos):
		# Input: genos = genotypes
		# Output: Predictions using betas from our method
		self.predictions = np.dot(genos, self.posterior_betas)
		return self.predictions

	def predict_using_prs(self, genos, sig=0.05, scale=100):
		# Define a significance threshold 
		threshold = (sig/len(self.beta_hats))*scale
		self.posterior_betas = self.beta_hats
		# Set non-significant betas to 0
		self.posterior_betas[abs(self.posterior_betas)<threshold]=0
		self.predict(genos)

	def predict_using_inf(self, n_idvs, genos, z_scores=False, ld_radius=100, D=None):
		# Input: n_idvs = number of individuals, genos = genotypes, ld_radius = LD block radius
		# Output: Prediction using the inf method 
		self.posterior_betas = self.compute_posterior_betas_inf(self.beta_hats, self.est_heritability, n_idvs, ld_radius, D)
		self.predict(genos)
		if z_scores:
			scale_factor = np.sqrt(np.var(self.predictions) / self.est_heritability)
			self.predictions = self.predictions / scale_factor

		return self.predictions

	def predict_using_gibbs(self, n_idvs, genos, z_scores=False, ld_radius=100, burn_in=10, iters=100, D=None):
		# Input: n_idvs = number of individuals, genos = genotypes, ld_radius = LD block radius, burn_in = burn in iterations, iters = number of iterations
		# Output: Prediction using the mcmc based method 
		self.posterior_betas, start_betas = self.compute_posterior_betas_gibbs(None, self.est_heritability, self.beta_hats, ld_radius, n_idvs, z_scores=z_scores, burn_in=burn_in, iters=iters, p_causal=0.1, zero_jump_prob=0.05, D=D)
		self.predict(genos)
		if z_scores:
			scale_factor = np.sqrt(np.var(self.predictions) / self.est_heritability)
			self.predictions = self.predictions / scale_factor
			
		return self.predictions 

	def compute_local_posterior_betas_inf(self, local_D, est_heritability, n_snps, n_idvs, local_betas):
		# Input: local_D = local LD matrix, est_heritability = estimated heritability, n_snps = number of snps, n_idvs = number of individuals, local_betas = local betas
		# Output: LD adjusted posterior betas
		diag_matrix = (n_snps / (n_idvs * est_heritability)) * np.identity(local_D.shape[0])
		local_D_weighted_inv = np.linalg.inv(diag_matrix + local_D)
		local_posterior_betas = np.dot(local_D_weighted_inv, local_betas)
		return local_posterior_betas

	def compute_posterior_betas_inf(self, beta_hats, est_heritability, ld_radius, n_idvs, D=None):
		# This function is the same as compute_local_posterior_betas_inf but we are iterating over all SNPs using a window
		# Input: See above for compute_local_posterior_betas_inf
		# Output: Same as compute_local_posterior_betas_inf
		n_snps = len(beta_hats)
		posterior_betas = np.zeros(n_snps)
		for s in xrange(0, n_snps):
			start_snp_idx = max(0, s - ld_radius)
			stop_snp_idx = min(n_snps, s + ld_radius + 1)
			local_betas = beta_hats[start_snp_idx:stop_snp_idx]
			local_D = D[start_snp_idx:stop_snp_idx, start_snp_idx:stop_snp_idx]	# TODO: compute local_D or get from DB
			local_posterior_betas = self.compute_local_posterior_betas_inf(local_D, est_heritability, n_snps, n_idvs, local_betas)
			posterior_betas[start_snp_idx:stop_snp_idx] = local_posterior_betas
		self.posterior_betas = posterior_betas
		return self.posterior_betas

	def compute_local_posterior_betas_gibbs(self, current_betas, current_posterior_means, 
											local_D, est_heritability, n_snps, n_idvs, local_betas, 
											global_snp_idx, ld_radius, p_causal, alpha, rand_ps, rand_norms,
											mp, hdmp, hdmpn, hdmp_hdmpn, c_const, d_const):
		# Input: current_betas = current updated betas, current_posterior_means = current beta means, local_D = local LD matrix, est_heritability = estimated heritability
		# n_snps = number of SNPs, n_idvs = number of individuals, local_betas = betas in current window, global_snp_idx = index of SNP in entire block,
		# ld_radius = LD block size divided by 2, p_causal = percent causal, for the rest of the definitions see function compute_posterior_betas_gibbs
		# Output: Posterior betas and the means
		focal_snp_idx = min(ld_radius, global_snp_idx)
		local_betas_fixed = np.copy(self.beta_hats, order='C')
		local_betas_var = np.copy(local_betas, order='C')
		local_betas_var[focal_snp_idx] = 0
		res_beta_hat_i = local_betas_fixed[focal_snp_idx] - np.dot(local_D[focal_snp_idx, :], local_betas_var)
		res_beta_sq = res_beta_hat_i ** 2
		postp = 1
		dbetasq_const = d_const * np.exp(-res_beta_sq * n_idvs / 2.0)
		if np.isreal(dbetasq_const):
			numerator = c_const * np.exp(-res_beta_sq / (2.0 * hdmpn))
			if numerator == 0:
				postp = 0
			elif np.isreal(numerator):
				postp = numerator / (numerator + dbetasq_const)
				assert np.isreal(postp)
		current_posterior_means[global_snp_idx] = hdmp_hdmpn * postp * res_beta_hat_i
		proposed_beta = 0
		if rand_ps[global_snp_idx] < (postp * alpha):
			proposed_beta = rand_norms[global_snp_idx] + hdmp_hdmpn * res_beta_hat_i
		current_betas[global_snp_idx] = proposed_beta
		return (current_betas, current_posterior_means)

	def compute_posterior_betas_gibbs(self, start_betas, est_heritability, beta_hats, ld_radius, n_idvs, z_scores=False, burn_in=10, iters=100, p_causal=0.1, zero_jump_prob=0.05, D=None):
		# This function computes the posterior betas and the means using an MCMC approach
		# See compute_local_posterior_betas_gibbs and predict_using_gibbs for input definitions
		# Output: Posterior betas and the means
		n_snps = len(beta_hats)
		if start_betas is None:
			start_betas = self.compute_posterior_betas_inf(beta_hats, est_heritability, ld_radius, n_idvs, D)
		current_betas = np.copy(start_betas)
		current_posterior_means = np.zeros(n_snps)
		average_betas = np.zeros(n_snps)
		mp = p_causal * n_snps
		hdmp = 1
		if not z_scores:
			hdmp = est_heritability / mp # Note that this should be 1 for z-scores, how should we force this to happen
		hdmpn = hdmp + 1.0 / n_idvs
		hdmp_hdmpn = hdmp / hdmpn
		c_const = p_causal / np.sqrt(hdmpn)
		d_const = (1.0 - p_causal) / (np.sqrt(1.0 / n_idvs))
		for i in xrange(0, iters+burn_in):
			sigma_g = max(0.00001, np.sum(current_betas ** 2))
			alpha = min(1.0 - zero_jump_prob, 1.0 / sigma_g, (est_heritability + 1.0 / np.sqrt(n_idvs)) / sigma_g)
			rand_ps = np.random.random(n_snps)
			# rand_norms = sp.stats.norm.rvs(0, (hdmp_hdmpn) * (1 / n_idvs), size=n_snps)
			rand_norms = np.random.normal(loc=0.0, scale=(hdmp_hdmpn) * (1.0 / n_idvs), size=n_snps)
			#print "Iter: ", i
			for s in xrange(0, n_snps):
				start_snp_idx = max(0, s - ld_radius)
				stop_snp_idx = min(n_snps, s + ld_radius + 1)
				local_betas = current_betas[start_snp_idx:stop_snp_idx]
				local_D = D[start_snp_idx:stop_snp_idx, start_snp_idx:stop_snp_idx]	# TODO: compute local_D or get from DB
				current_betas, current_posterior_means = self.compute_local_posterior_betas_gibbs(current_betas, current_posterior_means, 
																									local_D, est_heritability, n_snps, n_idvs, local_betas, 
																									s, ld_radius, p_causal, alpha, rand_ps, rand_norms,
																									mp, hdmp, hdmpn, hdmp_hdmpn, c_const, d_const)
			if i > burn_in:
				average_betas = average_betas + current_posterior_means
		average_betas = average_betas / (iters - burn_in)
		self.posterior_betas = average_betas
		return (self.posterior_betas, start_betas)

if __name__ == '__main__':
	main()

