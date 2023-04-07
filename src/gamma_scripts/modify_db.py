import pandas as pd

df = pd.read_hdf('fit_3d_db_30s_2000.h5')
df['tas1'] = df['tas']
df['branch_sum1'] = df['branch_sum']
df = df.drop(['branch_sum','tas'], axis=1)

df['tas'] = df['branch_sum1']
df['branch_sum'] = df['tas1']
df_new = df.drop(['branch_sum1','tas1'], axis=1)
df_new.to_hdf('fit_3d_db_30s_2000.h5',key='df')

