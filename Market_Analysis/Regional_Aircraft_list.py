import pandas as pd
aircraft = pd.read_csv('aircrafts.txt',names = ["Name","ICAO","IATA","Capacity","Country"],delimiter=';',na_values="-").dropna()
del aircraft['Country']
aircraft = aircraft[aircraft['Capacity'] <= 100]
aircraft = aircraft[aircraft['Capacity'] > 20]
aircraft_list = aircraft["IATA"].tolist()

