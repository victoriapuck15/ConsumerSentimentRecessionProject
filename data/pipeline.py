import pandas as pd

data = pd.read_csv('DP_LIVE_08112023214649475.csv')

data = data[["LOCATION", "TIME", "Value"]]

data = data.pivot(index='TIME', columns='LOCATION', values='Value')
data.columns = [col[:3] for col in data.columns]
data.index.name = 'Dates'
data.index = pd.to_datetime(data.index)
data = data[data.index.year >= 2017]
data = data[data.index.month.isin([1, 4, 7, 10])]
data.index = data.index.map(lambda x: x.replace(day=31) if x.month in [1, 7, 10] else x.replace(day=30))
data.index = data.index.strftime('%m/%d/%y')

data.to_excel("data/cci_final2.xlsx")