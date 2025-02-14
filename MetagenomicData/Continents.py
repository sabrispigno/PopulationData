import pandas as pd
# 1. Crea un dizionario che mappa ogni codice paese al suo continente
iso_to_continent = {
    "ITA": "Europa",
    "USA": "Nord America",
    "GBR": "Europa",
    "DEU": "Europa",
    "SWE": "Europa",
    "FJI": "Oceania",
    "CHN": "Asia",
    "SGP": "Asia",
    "IDN": "Asia",
    "MYS": "Asia",
    "BRN": "Asia",
    "KAZ": "Asia",  # (Transcontinentale, ma spesso associato all’Asia)
    "BGD": "Asia",
    "CAN": "Nord America",
    "DNK": "Europa",
    "LUX": "Europa",
    "FRA": "Europa",
    "TZA": "Africa",
    "CMR": "Africa",
    "NLD": "Europa",
    "NOR": "Europa",
    "SVK": "Europa",
    "HUN": "Europa",
    "EST": "Europa",
    "ISL": "Europa",
    "FIN": "Europa",
    "PHL": "Asia",
    "ESP": "Europa",
    "AUT": "Europa",
    "IND": "Asia",
    "MDG": "Africa",
    "ETH": "Africa",
    "GHA": "Africa",
    "SLV": "Nord America",  # El Salvador è in America Centrale, ma statisticamente Nord America
    "ISR": "Asia",
    "JPN": "Asia",
    "KOR": "Asia",
    "RUS": "Europa/Asia",  # Transcontinentale
    "PER": "Sud America",
    "MNG": "Asia",
    "LBR": "Africa",
    # ... (aggiungi qui tutti i codici che ti servono)
}

def add_continent_column(df, country_col="country_code", new_col="continent"):
    """
    Aggiunge la colonna `new_col` con il continente corrispondente ai codici
    ISO nella colonna `country_col` del DataFrame df.
    """
    df[new_col] = df[country_col].map(iso_to_continent)
    return df

