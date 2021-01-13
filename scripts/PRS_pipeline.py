import os
import argparse
import pandas as pd
import numpy as np
import sys

sys.path.append("/Users/crystalvaldez/Dropbox (Encompass Bioscience)/cv_sandbox/scripts")

import vcf_reader_prs as vcp
import enc_predictor as erp
from sklearn import linear_model
from sklearn.preprocessing import Imputer


def main():
	parser = argparse.ArgumentParser(description='')
	parser.add_argument('-v', '--input-vcf', action='store', dest='vcf', default='/extra/adni/output/adni_clumped.vcf', help='')
	parser.add_argument('-b', '--betas_file', action='store', dest='betas_file', default='/extra/adni/data/IGAP_stage_1.txt', help='')
	parser.add_argument('-a', '--pca-file', action='store', dest='pca_file', default='/extra/adni/data/adni.eigenvec', help='')
	parser.add_argument('-o', '--output-file', action='store', dest='output_file', default='test.txt', help='')

	args = parser.parse_args()

	# load the files that we need
	betas = load_betas(args.betas_file)

	# do predictions
	predictions = predict_on_vcf(args.vcf, args.pca_file, betas)
        predictions.to_csv(args.output_file, index=True, sep='\t', na_rep='NA')

def predict_on_vcf(vcf, pca_file, betas):
	complement_map = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
	vcf_reader = vcp.genoVCFreader(vcf)
	predictor = erp.EncPredictor(None, None, None)
	predictions = pd.DataFrame(0, index=remove_family_id_from_sample_id(vcf_reader.get_sample_list()), columns=['predictions'])
	rsids_to_flip = []
	try:
		while True:
			# generate genotypes by matching effect allele with VCF
			rsid, chrom, pos, ref, alt, genos = vcf_reader.next()
			snp_gwas_data = betas[ betas['RSID']==rsid ]
			snp_beta = snp_gwas_data['Beta']
			snp_beta_se = snp_gwas_data['SE']
			gwas_alt = snp_gwas_data['Effect_allele']
			gwas_ref = snp_gwas_data['Non_Effect_allele']
			if gwas_alt.values[0] == ref:
				rsids_to_flip.append(rsid)
			elif gwas_alt.values[0] == complement_map[ref]:
				rsids_to_flip.append(rsid)
			if len(rsids_to_flip) > 0:
				new_record = vcp.flip_alleles([[rsid, chrom, pos, ref, alt, genos]], rsids_to_flip)
				rsid, chrom, pos, ref, alt, genos = new_record[0]

			# fill in missing genotypes with 0, they will have no effect on the prediction
			genos = np.reshape(genos, (genos.shape[0], 1))

			# normalize genotypes using PCA
			pcs = load_pcs(pca_file)
			z_genos, est_freqs = normalize_genos(genos, pcs, predictions)
			filled_genos = np.nan_to_num(z_genos)

			# do the prediction
			zscore = np.array(snp_beta.values[0] / snp_beta_se.values[0])
			snp_predictions = predictor.static_predict(filled_genos, zscore)
                        #print z_genos[0:10], zscore, snp_predictions[0:10], est_freqs[0:10], genos[0:10]
                        #raw_input()
			predictions['predictions'] += snp_predictions[:, 0]
	except StopIteration:
		pass

	return predictions


def remove_family_id_from_sample_id(sample_ids):
	new_sample_ids = []
	for sid in sample_ids:
		data = sid.split('_')
		new_sid = '_'.join(data[1:])
		new_sample_ids.append(new_sid)
	return new_sample_ids


def normalize_genos(genos, pcs, predictions):
	# fit the PCs
	imp = Imputer(missing_values='NaN', strategy='mean', axis=0)
	imp_genos = imp.fit_transform(genos)
	ols = linear_model.LinearRegression()

        #print pcs.ix[predictions.index].values
	ols.fit(pcs.ix[predictions.index].values, imp_genos)

		# predict and compute the frequencies
	geno_preds = ols.predict(pcs.ix[predictions.index].values)
	est_freqs = geno_preds / 2

	# fill in the outliers
	epsilon = 0.001
	est_freqs[ est_freqs < 0 ] = epsilon
	est_freqs[ est_freqs > 1 ] = 1 - epsilon

	# normalize the genotypes
	z_genos = (genos - 2*est_freqs) / np.sqrt(2*est_freqs*(1-est_freqs))

	return z_genos, est_freqs


def load_betas(betas_file):
	betas = pd.read_csv(betas_file, sep='\t')
	return betas


def load_pcs(pca_file):
	pcs = pd.read_csv(pca_file, sep=' ', header=None)
	pcs_index = pcs.iloc[: , 1]
	pcs = pcs.iloc[:, 2:12]
	pcs.index = map(str, pcs_index)
	pcs.columns = range(0, pcs.shape[1])
	return pcs


if __name__ == '__main__':
	main()
