import numpy as np
import pandas as pd

def load_datasets(df):
	"""Computes the separate datasets based on a single input dataframe, df.
	"""
	# ========== Raw abundance ==========
	d_means = df.groupby('Group').mean().reset_index()
	d_sd = df.groupby("Group").agg(np.std, ddof=0).reset_index().drop(["Sample"], axis=1)

	# transpose
	d_means2 = d_means.set_index('Group').T
	d_sd2 = d_sd.set_index('Group').T
	d_means2.columns = ["ctrl_mean", "trtmt_mean"]
	d_sd2.columns = ["ctrl_sd", "trtmt_sd"]

	# merge
	dz = pd.concat([d_means2, d_sd2], axis=1).sort_values(by="ctrl_mean")

	# ========== log abundance ==========
	colnames_ls = list(df.columns)[2:df.shape[1]]
	min_nonzero = np.sort(np.unique(df[colnames_ls].values.flatten()))[1]/2
	d_t = df.replace(0, min_nonzero)
	d_log = df[["Sample", "Group"]]

	for cnm in colnames_ls:
	    new_colname = "log_"+cnm
	    d_log[new_colname] = np.log2(d_t[cnm])

	dl_means = d_log.groupby('Group').mean().reset_index()
	dl_sd = d_log.groupby("Group").agg(np.std, ddof=0).reset_index().drop(["Sample"], axis=1)
	dl_sd = dl_sd
	# transpose
	dl_means2 = dl_means.set_index('Group').T
	dl_sd2 = dl_sd.set_index('Group').T

	dl_means2.columns = ["ctrl_mean", "trtmt_mean"]
	dl_sd2.columns = ["ctrl_sd", "trtmt_sd"]
	dz2 = pd.concat([dl_means2, dl_sd2], axis=1).sort_values(by="ctrl_mean")

	# ========== log and median-normalized ==========
	row_medians_ls = []
	for i in range(d_log.shape[0]):
	    row_medians_ls.append(np.median(list(d_log.iloc[i])[2:d_log.shape[1]]))

	contents = []
	for i in range(len(row_medians_ls)):
	    row = np.array(d_log.iloc[i])
	    row[2:len(row)] -= row_medians_ls[i]
	    contents.append(row)
	d_logm = pd.DataFrame(data=contents, columns=["Sample", "Group"] + ["logm_"+x for x in colnames_ls])

	dlm_means = d_logm.groupby('Group').mean().reset_index()
	dlm_sd = d_logm.groupby("Group").agg(np.std, ddof=0).reset_index().drop(["Sample"], axis=1)
	# transpose
	dlm_means2 = dlm_means.set_index('Group').T
	dlm_sd2 = dlm_sd.set_index('Group').T

	dlm_means2.columns = ["ctrl_mean", "trtmt_mean"]
	dlm_sd2.columns = ["ctrl_sd", "trtmt_sd"]
	dz3 = pd.concat([dlm_means2, dlm_sd2], axis=1).sort_values(by="ctrl_mean")

	return(d_log, d_logm, dz, dz2, dz3)