import pandas as pd

pheno_df = pd.read_csv(pheno_path)
diag_df = pd.read_csv(diag_path)


# exclude:
# ILD - J84
# Heart failure: I50
ild = pd.read_csv(ild_path)
heart_failure = pd.read_csv(heart_failure_path)

exclude_ids = set(ild['Subject_Id']).union(set(heart_failure['Subject_Id']))
pheno_df = pheno_df[~pheno_df['Subject_Id'].isin(exclude_ids)].copy()


################# get exacerbation #################
# ICD9 491.21, 493.22, ICD10 J44.0, J44.1
exac_df_all = diag_df[
    (((diag_df['Code_Type'] == 'ICD10') & (diag_df['Code'].isin(['J440', 'J441']))) |
    ((diag_df['Code_Type'] == 'ICD9') & (diag_df['Code'].isin(['49121', '49322']))))
].copy()

exac_df_all['Date'] = pd.to_datetime(exac_df_all["Date"])
# first exacerbation
exac_df_all.loc[:, 'first_date'] = exac_df_all.groupby('Subject_Id')['Date'].transform('min')

# keep records that are one year from the first exacerbation event
exac_df = exac_df_all[exac_df_all['Date'] <= (exac_df_all['first_date'] + pd.DateOffset(years=1))]

# multiple exacerbation events within a 15 days period is considered one event
exac_df = exac_df.sort_values(by=['Subject_Id', 'Date'])
exac_df['exacerbation'] = 0

# to keep track of the last exacerbation date
last_exacerbation_date = {}

def make_exacerbation(row):
    person_id = row['Subject_Id']
    curr_date = row['Date']

    if person_id not in last_exacerbation_date:
        # first exacerbation for this person
        last_exacerbation_date[person_id] = curr_date
        return 1
    
    elif (curr_date - last_exacerbation_date[person_id]).days > 15:
        # if the current date is more than 15 days after the last counted exacerbation
        last_exacerbation_date[person_id] = curr_date
        return 1
    
    else:
        return 0
    
exac_df['exacerbation'] = exac_df.apply(make_exacerbation, axis=1)

# get total exacerbations for each person
exac_counts = exac_df.groupby('Subject_Id')['exacerbation'].sum().reset_index()
exac_counts = exac_counts.rename(columns={'exacerbation': 'total_exacerbations'})


################# get COPD #################
# ICD9 491.xx, 492.xx, 496, ICD10 codes J41.x, J42, J43.x, J44.x
copd_df_all = diag_df[
    (((diag_df['Code_Type'] == 'ICD9') & (diag_df['Code'].fillna("").str.startswith(("491", "492", "496")))) |
    ((diag_df['Code_Type'] == 'ICD10') & (diag_df['Code'].fillna("").str.startswith(("J41", "J42", "J43", "J44")))))
].copy()

# filter for  ≥1 inpatient code or ≥3 outpatient code
copd_grouped = copd_df_all.groupby(['Subject_Id', 'Inpatient_Outpatient']).size().unstack(fill_value=0)

copd_cases = copd_grouped[
    (copd_grouped.get('Inpatient', 0) >= 1) |
    (copd_grouped.get('Outpatient', 0) >= 3)
].copy()

copd_cases['copd'] = 1
copd_cases = copd_cases.reset_index()
copd_df = copd_cases[['Subject_Id', 'copd']]

# merge with pheno
pheno_df = pheno_df.merge(exac_counts, on='Subject_Id', how='left')
pheno_df = pheno_df.merge(copd_df, on='Subject_Id', how='left')

pheno_df['total_exacerbations'].fillna(0, inplace=True)
pheno_df['copd'].fillna(0, inplace=True)
pheno_df['total_exacerbations'] = pheno_df['total_exacerbations'].astype('int64')
pheno_df['copd'] = pheno_df['copd'].astype('int64')

# people with exacerbation are considered copd cases
pheno_df.loc[pheno_df['total_exacerbations'] > 0, 'copd'] = 1
